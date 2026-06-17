//! Gray flux-limited radiation diffusion (ADR-0006), built test-first.
//!
//! The 1D rad-hydro takes the radiation step as an **operator-split implicit diffusion solve on
//! the same staggered mesh** (ADR-0022). This module holds the pieces that solve in their own
//! right, before any coupling to the gas:
//!
//! - [`thomas_solve`] — the tridiagonal (Thomas) solver every implicit 1D diffusion step reduces
//!   to. A pure helper with a closed form, so it gets a direct unit test (per `CLAUDE.md`).
//! - [`flux_limiter`] — the **Levermore–Pomraning** limiter `λ(R)` that ties the two regimes of
//!   ADR-0006 together: `λ → 1/3` (Fickian diffusion, `τ ≫ 1`) and `λ → 1/R` (free-streaming,
//!   `τ ≪ 1`, capping the flux at `cE` so radiation never moves faster than light).
//!
//! The coupled gray-FLD substep (Planck-mean emission/absorption + Rosseland-mean flux-limited
//! diffusion) is assembled on top of these in B3b.

/// Solve the tridiagonal system `M x = rhs`, where row `i` of `M` is
/// `sub[i]·x[i−1] + diag[i]·x[i] + sup[i]·x[i+1]`. `sub[0]` and `sup[n−1]` lie outside the matrix
/// and are ignored. This is the Thomas algorithm — Gaussian elimination specialized to a
/// tridiagonal, `O(n)` and allocation-light. No pivoting: valid for the diagonally dominant,
/// positive-definite systems that implicit diffusion produces (always the case here, since every
/// off-diagonal is `≤ 0` and the diagonal exceeds their magnitude).
///
/// # Panics
/// Panics unless all four slices have the same length `n ≥ 1`.
#[must_use]
pub fn thomas_solve(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    assert!(n >= 1, "empty tridiagonal system");
    assert!(
        sub.len() == n && sup.len() == n && rhs.len() == n,
        "tridiagonal slices must share one length"
    );
    // Forward sweep: eliminate the sub-diagonal, carrying modified super-diagonal `c` and rhs `d`.
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];
    c[0] = sup[0] / diag[0];
    d[0] = rhs[0] / diag[0];
    for i in 1..n {
        let denom = diag[i] - sub[i] * c[i - 1];
        c[i] = sup[i] / denom;
        d[i] = (rhs[i] - sub[i] * d[i - 1]) / denom;
    }
    // Back substitution.
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d[i] - c[i] * x[i + 1];
    }
    x
}

/// Levermore–Pomraning flux limiter `λ(R) = (1/R)(coth R − 1/R)`, where the dimensionless
/// `R = |∇E| / (κ_R ρ E)` measures how steep the radiation gradient is on the scale of a mean free
/// path. The flux-limited diffusion coefficient is `D = c λ(R) / (κ_R ρ)`, and the flux magnitude
/// `|F| = D|∇E| = c λ R E`.
///
/// Asymptotics (the whole point, ADR-0006):
/// - `R → 0` (optically thick / shallow gradient): `λ → 1/3`, so `D → c/(3 κ_R ρ)` — ordinary
///   **Fickian** radiation diffusion.
/// - `R → ∞` (optically thin / steep gradient): `λ → 1/R`, so `|F| → c E` — **free-streaming**,
///   the flux saturated at light speed, never superluminal.
///
/// `λ` is smooth and monotone decreasing on `[1/3, 0)`; the small-`R` series avoids the `0/0` in
/// the closed form.
#[must_use]
pub fn flux_limiter(r: f64) -> f64 {
    if r < 1e-3 {
        // coth R − 1/R = R/3 − R³/45 + …  ⇒  λ = 1/3 − R²/45 + …
        (r * r).mul_add(-1.0 / 45.0, 1.0 / 3.0)
    } else {
        // Stable for all R ≥ 1e-3: coth R = 1/tanh R → 1 and 1/R → 0 as R → ∞, giving λ → 1/R.
        (1.0 / r.tanh() - 1.0 / r) / r
    }
}

#[cfg(test)]
mod tests {
    use super::{flux_limiter, thomas_solve};
    use approx::assert_relative_eq;

    /// Thomas solve against a hand-built system whose answer is known: the 1D Poisson matrix
    /// `tridiag(−1, 2, −1)` applied to `x = [1,2,3,4]` gives a known rhs; recover `x`.
    #[test]
    fn thomas_solves_known_system() {
        let x_true = [1.0, 2.0, 3.0, 4.0];
        let sub = [0.0, -1.0, -1.0, -1.0];
        let diag = [2.0, 2.0, 2.0, 2.0];
        let sup = [-1.0, -1.0, -1.0, 0.0];
        // rhs = M·x_true.
        let n = 4;
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            rhs[i] = diag[i] * x_true[i];
            if i > 0 {
                rhs[i] += sub[i] * x_true[i - 1];
            }
            if i + 1 < n {
                rhs[i] += sup[i] * x_true[i + 1];
            }
        }
        let x = thomas_solve(&sub, &diag, &sup, &rhs);
        for (got, want) in x.iter().zip(x_true.iter()) {
            assert_relative_eq!(got, want, max_relative = 1e-13);
        }
    }

    /// A 1×1 system is just division.
    #[test]
    fn thomas_solves_scalar() {
        let x = thomas_solve(&[0.0], &[4.0], &[0.0], &[12.0]);
        assert_relative_eq!(x[0], 3.0, max_relative = 1e-15);
    }

    /// Limiter asymptotics: Fickian `1/3` as `R → 0`, free-streaming `1/R` as `R → ∞`.
    #[test]
    fn limiter_has_correct_asymptotes() {
        assert_relative_eq!(flux_limiter(1e-8), 1.0 / 3.0, max_relative = 1e-12);
        assert_relative_eq!(flux_limiter(0.0), 1.0 / 3.0, max_relative = 1e-12);
        for &r in &[1e3, 1e5, 1e7] {
            assert_relative_eq!(flux_limiter(r), 1.0 / r, max_relative = 1e-3);
        }
    }

    /// `λ` is continuous across the series/closed-form switch and monotone decreasing in between.
    #[test]
    fn limiter_is_continuous_and_monotone() {
        // Continuity at the cutoff: series just below vs closed form just above.
        let below = flux_limiter(1e-3 - 1e-9);
        let just_above = (1.0 / 1e-3_f64.tanh() - 1.0 / 1e-3) / 1e-3;
        assert_relative_eq!(below, just_above, max_relative = 1e-6);
        // Monotone decreasing over a wide range.
        let rs = [1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e4];
        let mut prev = flux_limiter(0.0);
        for &r in &rs {
            let v = flux_limiter(r);
            assert!(v < prev, "λ not decreasing at R={r}: {v} ≥ {prev}");
            assert!(v > 0.0 && v <= 1.0 / 3.0, "λ out of (0, 1/3] at R={r}: {v}");
            prev = v;
        }
    }

    // --- erfc half-space: the implicit linear-diffusion assembly built on `thomas_solve` ---

    /// `erf` via Abramowitz & Stegun 7.1.26 (|error| ≤ 1.5e-7), enough for a ~1e-3 acceptance.
    fn erf(x: f64) -> f64 {
        let t = 1.0 / (0.327_591_1_f64.mul_add(x.abs(), 1.0));
        let poly = t
            * (0.254_829_592
                + t * (-0.284_496_736
                    + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
        let y = 1.0 - poly * (-x * x).exp();
        if x < 0.0 { -y } else { y }
    }

    fn erfc(x: f64) -> f64 {
        1.0 - erf(x)
    }

    /// Diffusion into a cold half-space with a fixed surface value has the closed-form solution
    /// `E(x, t) = E₀ · erfc(x / (2√(Dt)))`. We assemble one **backward-Euler** step of linear
    /// diffusion (constant `D`, Dirichlet `E₀` at `x = 0`, zero-flux far end) as a tridiagonal
    /// system and solve it with [`thomas_solve`], marching to `t`; the result must match `erfc` to
    /// the scheme's discretization error. This verifies the implicit-diffusion machinery end to
    /// end before any radiation coupling.
    #[test]
    fn implicit_diffusion_matches_erfc_half_space() {
        let n = 400usize;
        let l = 1.0;
        let dx = l / n as f64;
        let diff = 1.0; // D
        let dt = 1e-5;
        let t_final = 0.01; // √(Dt) = 0.1 ≪ L, so the far boundary stays cold
        let e0 = 1.0;
        let r = diff * dt / (dx * dx);

        let mut e = vec![0.0; n];
        let sub = {
            let mut s = vec![-r; n];
            s[0] = 0.0;
            s
        };
        let sup = {
            let mut s = vec![-r; n];
            s[n - 1] = 0.0;
            s
        };
        let mut diag = vec![1.0 + 2.0 * r; n];
        diag[0] = 1.0 + 3.0 * r; // Dirichlet half-cell to the boundary at x=0
        diag[n - 1] = 1.0 + r; // zero-flux far end

        let steps = (t_final / dt).round() as usize;
        for _ in 0..steps {
            let mut rhs = e.clone();
            rhs[0] += 2.0 * r * e0; // Dirichlet boundary contribution
            e = thomas_solve(&sub, &diag, &sup, &rhs);
        }

        let mut linf = 0.0_f64;
        for (j, &ej) in e.iter().enumerate() {
            let x = (j as f64 + 0.5) * dx;
            let exact = e0 * erfc(x / (2.0 * (diff * t_final).sqrt()));
            linf = linf.max((ej - exact).abs());
        }
        assert!(linf < 5e-3, "erfc half-space L∞ error = {linf:e}");
    }
}
