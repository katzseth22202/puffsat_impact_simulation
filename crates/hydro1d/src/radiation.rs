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
//! The coupled gray-FLD substep ([`fld_substep`]) is assembled on top of these (B3b): Planck-mean
//! emission/absorption `c χ_P (a T⁴ − E)` exchanging energy with the matter, and Rosseland-mean
//! flux-limited diffusion transporting it, solved **linearized-implicit** so the stiff
//! emission/absorption and diffusion are both stable at a hydro-sized timestep. Radiation pressure
//! and radiation work on the gas are deferred (revisited at the 16 km/s anchor, B5); this step
//! moves radiation *energy* and sets the matter temperature that drives opacity.

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

/// Physical constants for the gray radiation step (consistent units; SI for the production run,
/// dimensionless for the verification tests).
#[derive(Debug, Clone, Copy)]
pub struct RadConstants {
    /// Speed of light `c`.
    pub c: f64,
    /// Radiation constant `a` (so the equilibrium radiation energy density is `a T⁴`).
    pub a: f64,
}

/// Which flux limiter the diffusion step uses (ADR-0006).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Limiter {
    /// Levermore–Pomraning `λ(R)` — the production default. Reduces to Fickian `1/3` when optically
    /// thick and caps the flux at `cE` when thin (see [`flux_limiter`]).
    LevermorePomraning,
    /// Forced Fickian `λ ≡ 1/3` — ordinary (un-limited) radiation diffusion. Used to verify the
    /// solver against pure-diffusion benchmarks (Su–Olson Marshak), where the limiter would
    /// otherwise pull the solution off the analytic diffusion answer near steep fronts.
    Fick,
}

/// Radiation boundary condition at a domain end.
#[derive(Debug, Clone, Copy)]
pub enum RadBc {
    /// Zero-flux (symmetry / reflecting) boundary.
    Reflecting,
    /// Fixed radiation energy density at the boundary face (half a cell beyond the edge center).
    Dirichlet(f64),
    /// Incident-flux (**Marshak**) boundary: an external source of radiation energy density `e_inc`
    /// drives the surface through the net inward flux `F(0) = (c/2)(e_inc − E_surface)` — the
    /// diffusion-limit incident-current condition `cE/4 + F/2 = c·e_inc/4`. This is the radiation
    /// source-surface used by the Su–Olson Marshak benchmark (and the wall absorber of B4).
    Marshak(f64),
    /// A Marshak surface with its conductance scaled by `transmission ∈ (0, 1]` — an intervening gray
    /// absorbing layer (the Rung E vapor curtain, ADR-0014) of transmission `1/(1+τ)` in series with
    /// the surface, throttling the net flux to `transmission·(c/2)(e_inc − E_surface)`. The
    /// `(1−transmission)` fraction is *retained in the radiation field* near the surface (it does not
    /// leave through the throttled boundary), so the FLD couples it back into the gas self-consistently
    /// — the energy-conserving "shielded wall". `transmission = 1` is exactly [`RadBc::Marshak`].
    MarshakAttenuated { e_inc: f64, transmission: f64 },
}

/// The per-cell material state the radiation step reads (a frozen 1D Lagrangian mesh). All slices
/// have length `N` except `center_spacing`, the `N−1` center-to-center distances of the interior
/// faces. Opacities are **per length** (`χ = κ ρ`); `cv_vol` is the heat capacity per unit volume
/// `ρ ∂e/∂T`.
#[derive(Debug, Clone, Copy)]
pub struct Medium<'a> {
    /// Cell widths.
    pub dx: &'a [f64],
    /// Center-to-center spacings of the `N−1` interior faces.
    pub center_spacing: &'a [f64],
    /// Matter temperature `T`.
    pub temp: &'a [f64],
    /// Heat capacity per unit volume `ρ ∂e/∂T`.
    pub cv_vol: &'a [f64],
    /// Planck-mean absorption coefficient `χ_P = κ_P ρ` (emission/absorption source; ADR-0006).
    pub chi_planck: &'a [f64],
    /// Rosseland-mean coefficient `χ_R = κ_R ρ` (flux-limited diffusion; ADR-0006).
    pub chi_ross: &'a [f64],
    /// Optional volumetric radiation source `S` (energy / volume / time), e.g. the Su–Olson drive.
    pub source: Option<&'a [f64]>,
}

/// Diffusion coefficient `D = c λ / χ_R` at a face. With [`Limiter::LevermorePomraning`], `λ = λ(R)`
/// with `R = |∇E| / (χ_R E)` evaluated from the (lagged) radiation field; with [`Limiter::Fick`],
/// `λ ≡ 1/3` (ordinary diffusion). A transparent face (`χ_R ≤ 0`) carries no diffusive coupling.
fn face_diffusion(limiter: Limiter, c: f64, chi_r: f64, e_lo: f64, e_hi: f64, spacing: f64) -> f64 {
    if chi_r <= 0.0 {
        return 0.0;
    }
    let lambda = match limiter {
        Limiter::Fick => 1.0 / 3.0,
        Limiter::LevermorePomraning => {
            let e_face = 0.5 * (e_lo + e_hi);
            let r = if e_face > 0.0 {
                (e_hi - e_lo).abs() / (spacing * chi_r * e_face)
            } else {
                0.0
            };
            flux_limiter(r)
        }
    };
    c * lambda / chi_r
}

/// One **linearized-implicit** gray flux-limited-diffusion substep over `dt` on the frozen mesh
/// `medium`. Updates the radiation energy density `e_rad` in place and returns the per-cell matter
/// **internal-energy change** `Δe` (energy / volume) for the caller to deposit in the gas and
/// re-invert the EOS for the new temperature.
///
/// # Method
///
/// Per cell, linearize the emission `a T⁴ ≈ a T_n⁴ + β δT` (`β = 4 a T_n³`) and eliminate the
/// local matter response analytically, which folds a factor `f = C_v / (C_v + dt c χ_P β)` into the
/// exchange coefficient `k = f c χ_P`. The radiation update is then the **linear** system
///
/// ```text
/// E_j/dt − ∇·(D ∇E)_j + k_j E_j = E_j^n/dt + k_j a T_j⁴ + S_j
/// ```
///
/// — tridiagonal in `E^{n+1}` (Rosseland-mean flux-limited `D` at faces, lagged in the limiter),
/// solved with [`thomas_solve`]. The matter then takes `Δe_j = dt k_j (E^{n+1}_j − a T_j⁴)`, which
/// is *exactly* the energy the radiation lost to the local coupling, so matter + radiation energy is
/// conserved to round-off (diffusion is conservative by construction; the only sinks are the
/// boundary faces and `S`).
///
/// Returning **energy** rather than `δT = Δe / C_v` keeps the step robust where the heat capacity is
/// small or steeply varying (e.g. the `C_v = α T³` Su–Olson cold front, or ionization in the
/// production water table): the caller advances `e` and inverts the EOS for `T`, instead of dividing
/// a finite absorbed energy by a vanishing `C_v`.
///
/// # Panics
/// Panics unless the `medium` slice lengths are consistent (`N` cells, `N−1` interior faces) and
/// match `e_rad`.
#[must_use]
pub fn fld_substep(
    medium: &Medium<'_>,
    e_rad: &mut [f64],
    bc_left: RadBc,
    bc_right: RadBc,
    dt: f64,
    consts: RadConstants,
    limiter: Limiter,
) -> Vec<f64> {
    let n = medium.dx.len();
    assert!(n >= 1, "empty medium");
    assert!(
        medium.temp.len() == n
            && medium.cv_vol.len() == n
            && medium.chi_planck.len() == n
            && medium.chi_ross.len() == n
            && e_rad.len() == n
            && medium.center_spacing.len() == n - 1,
        "inconsistent medium/e_rad lengths"
    );
    let RadConstants { c, a } = consts;

    let mut sub = vec![0.0; n];
    let mut diag = vec![0.0; n];
    let mut sup = vec![0.0; n];
    let mut rhs = vec![0.0; n];

    // Local emission/absorption coupling, with the matter back-reaction folded into `k`.
    let mut beta = vec![0.0; n];
    let mut k = vec![0.0; n];
    for j in 0..n {
        beta[j] = 4.0 * a * medium.temp[j].powi(3);
        let denom = medium.cv_vol[j] + dt * c * medium.chi_planck[j] * beta[j];
        let f = if denom > 0.0 {
            medium.cv_vol[j] / denom
        } else {
            0.0
        };
        k[j] = f * c * medium.chi_planck[j];
        diag[j] = 1.0 / dt + k[j];
        rhs[j] = e_rad[j] / dt + k[j] * a * medium.temp[j].powi(4);
        if let Some(s) = medium.source {
            rhs[j] += s[j];
        }
    }

    // Interior-face flux-limited diffusion (conservative: one face feeds both its cells).
    for i in 0..n - 1 {
        let chi_r = 0.5 * (medium.chi_ross[i] + medium.chi_ross[i + 1]);
        let d = face_diffusion(
            limiter,
            c,
            chi_r,
            e_rad[i],
            e_rad[i + 1],
            medium.center_spacing[i],
        );
        let flux = d / medium.center_spacing[i];
        let w_lo = flux / medium.dx[i];
        let w_hi = flux / medium.dx[i + 1];
        diag[i] += w_lo;
        sup[i] -= w_lo;
        diag[i + 1] += w_hi;
        sub[i + 1] -= w_hi;
    }

    // Boundary faces: Dirichlet contributes a diffusive half-cell flux; Marshak contributes the
    // incident-current surface flux F = (c/2)(e_inc − E_edge); Reflecting contributes nothing.
    let last = n - 1;
    let bctx = BoundaryCtx { limiter, c, medium };
    add_boundary(&bctx, bc_left, 0, e_rad[0], &mut diag, &mut rhs);
    add_boundary(&bctx, bc_right, last, e_rad[last], &mut diag, &mut rhs);

    let e_new = thomas_solve(&sub, &diag, &sup, &rhs);

    // Matter response: exactly the energy the radiation lost to the local coupling. Returned as an
    // internal-energy change `Δe` (not `δT`); `k` already carries the implicit factor `f`, so this
    // is conservative against the radiation update and stays finite as `C_v → 0`.
    let mut delta_e = vec![0.0; n];
    for j in 0..n {
        delta_e[j] = dt * k[j] * (e_new[j] - a * medium.temp[j].powi(4));
    }
    e_rad.copy_from_slice(&e_new);
    delta_e
}

/// Inputs shared by both ends when folding a boundary face into the tridiagonal system: the flux
/// limiter, the speed of light, and the medium. Bundled so [`add_boundary`] stays a tidy helper.
struct BoundaryCtx<'a, 'm> {
    limiter: Limiter,
    c: f64,
    medium: &'a Medium<'m>,
}

/// Fold one end's boundary condition into the tridiagonal `diag`/`rhs` for its `edge` cell, given
/// the lagged radiation energy density `e_edge` there. Both ends share the same algebra:
/// - [`RadBc::Reflecting`]: zero flux, no contribution.
/// - [`RadBc::Dirichlet`]: a diffusive flux across the half-cell to a fixed `e_b` at the face (the
///   limiter sees `|∇E|`, so the `(e_b, e_edge)` order does not matter).
/// - [`RadBc::Marshak`]: the incident-current surface flux `F = (c/2)(e_inc − E_edge)`, implicit in
///   the edge `E` — a surface conductance `(c/2)/dx`, independent of the interior diffusion `D`.
fn add_boundary(
    ctx: &BoundaryCtx<'_, '_>,
    bc: RadBc,
    edge: usize,
    e_edge: f64,
    diag: &mut [f64],
    rhs: &mut [f64],
) {
    let medium = ctx.medium;
    match bc {
        RadBc::Reflecting => {}
        RadBc::Dirichlet(e_b) => {
            let dist = 0.5 * medium.dx[edge];
            let d = face_diffusion(ctx.limiter, ctx.c, medium.chi_ross[edge], e_b, e_edge, dist);
            let w = d / (dist * medium.dx[edge]);
            diag[edge] += w;
            rhs[edge] += w * e_b;
        }
        RadBc::Marshak(e_inc) => {
            let w = 0.5 * ctx.c / medium.dx[edge];
            diag[edge] += w;
            rhs[edge] += w * e_inc;
        }
        RadBc::MarshakAttenuated {
            e_inc,
            transmission,
        } => {
            let w = transmission * 0.5 * ctx.c / medium.dx[edge];
            diag[edge] += w;
            rhs[edge] += w * e_inc;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Limiter, Medium, RadBc, RadConstants, face_diffusion, fld_substep, flux_limiter,
        thomas_solve,
    };
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

        // SAFE: t_final, dt > 0 ⇒ the rounded ratio is a small non-negative integer (well below
        // usize::MAX), so neither truncation nor sign loss can occur.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
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

    // --- B3b: coupled gray-FLD substep ---

    /// Newton-solve `c_v T + a T⁴ = total` for the matter–radiation equilibrium temperature.
    fn equilibrium_temperature(cv: f64, a: f64, total: f64) -> f64 {
        let mut t = 1.0;
        for _ in 0..100 {
            let g = cv.mul_add(t, a * t.powi(4)) - total;
            let gp = cv + 4.0 * a * t.powi(3);
            t -= g / gp;
        }
        t
    }

    /// A box of matter and radiation out of balance relaxes to the state where `E = aT⁴`, while
    /// total (matter + radiation) energy is conserved to round-off at **every** step. The final
    /// temperature must match the energy-conserving equilibrium computed independently.
    #[test]
    fn equilibration_relaxes_to_energy_conserving_state() {
        let n = 4;
        let consts = RadConstants { c: 1.0, a: 1.0 };
        let dx = vec![1.0; n];
        let center_spacing = vec![1.0; n - 1];
        let cv_vol = vec![1.0; n];
        let chi_planck = vec![1.0; n];
        let chi_ross = vec![1.0; n];
        let mut temp = vec![1.0; n]; // aT⁴ = 1
        let mut e_rad = vec![5.0; n]; // radiation hotter than matter
        let energy = |t: &[f64], e: &[f64]| -> f64 {
            (0..n).map(|j| dx[j] * (cv_vol[j] * t[j] + e[j])).sum()
        };
        let total0 = energy(&temp, &e_rad);

        let dt = 0.05;
        for _ in 0..1000 {
            let medium = Medium {
                dx: &dx,
                center_spacing: &center_spacing,
                temp: &temp,
                cv_vol: &cv_vol,
                chi_planck: &chi_planck,
                chi_ross: &chi_ross,
                source: None,
            };
            let delta_e = fld_substep(
                &medium,
                &mut e_rad,
                RadBc::Reflecting,
                RadBc::Reflecting,
                dt,
                consts,
                Limiter::LevermorePomraning,
            );
            // Constant C_v = 1, so the (linear) EOS inversion is T += Δe / C_v = Δe.
            for (t, de) in temp.iter_mut().zip(delta_e.iter()) {
                *t += de;
            }
            assert_relative_eq!(energy(&temp, &e_rad), total0, max_relative = 1e-10);
        }

        let per_vol = total0 / dx.iter().sum::<f64>();
        let t_eq = equilibrium_temperature(1.0, 1.0, per_vol);
        let e_eq = t_eq.powi(4); // a = 1
        for j in 0..n {
            assert_relative_eq!(temp[j], t_eq, max_relative = 1e-6);
            assert_relative_eq!(e_rad[j], e_eq, max_relative = 1e-6);
        }
    }

    /// With a spatially non-uniform radiation field and reflecting ends, the flux-limited diffusion
    /// term redistributes energy while the coupling exchanges it with matter — and the total stays
    /// conserved to round-off (the diffusion is conservative) at every step, with the field
    /// spreading out (max − min shrinking) as it diffuses.
    #[test]
    fn diffusion_plus_coupling_conserves_energy() {
        let n = 12;
        let consts = RadConstants { c: 1.0, a: 1.0 };
        let dx = vec![1.0; n];
        let center_spacing = vec![1.0; n - 1];
        let cv_vol = vec![1.0; n];
        let chi_planck = vec![0.2; n];
        let chi_ross = vec![0.5; n];
        let mut temp = vec![1.0; n];
        // A radiation bump in the middle.
        let mut e_rad: Vec<f64> = (0..n)
            .map(|j| if (4..8).contains(&j) { 8.0 } else { 1.0 })
            .collect();
        let energy = |t: &[f64], e: &[f64]| -> f64 {
            (0..n).map(|j| dx[j] * (cv_vol[j] * t[j] + e[j])).sum()
        };
        let total0 = energy(&temp, &e_rad);
        let spread0 = e_rad.iter().copied().fold(f64::MIN, f64::max)
            - e_rad.iter().copied().fold(f64::MAX, f64::min);

        let dt = 0.05;
        for _ in 0..200 {
            let medium = Medium {
                dx: &dx,
                center_spacing: &center_spacing,
                temp: &temp,
                cv_vol: &cv_vol,
                chi_planck: &chi_planck,
                chi_ross: &chi_ross,
                source: None,
            };
            let delta_e = fld_substep(
                &medium,
                &mut e_rad,
                RadBc::Reflecting,
                RadBc::Reflecting,
                dt,
                consts,
                Limiter::LevermorePomraning,
            );
            // Constant C_v = 1, so the (linear) EOS inversion is T += Δe / C_v = Δe.
            for (t, de) in temp.iter_mut().zip(delta_e.iter()) {
                *t += de;
            }
            assert_relative_eq!(energy(&temp, &e_rad), total0, max_relative = 1e-10);
        }

        let spread = e_rad.iter().copied().fold(f64::MIN, f64::max)
            - e_rad.iter().copied().fold(f64::MAX, f64::min);
        assert!(
            spread < spread0,
            "diffusion did not smooth the field: {spread} vs {spread0}"
        );
    }

    // --- B3c-1: free-streaming cap (Levermore–Pomraning limiter) ---

    /// The radiative flux a face carries in the substep is `F = D·(E_lo − E_hi)/spacing` with the
    /// flux-limited `D = c λ(R)/χ_R` ([`face_diffusion`]). Because `λ(R)·R = coth R − 1/R ≤ 1`, this
    /// gives `|F| = c·λ(R)·R·E_face ≤ c·E_face` for **every** gradient — radiation never streams
    /// faster than light — and it **saturates** at `c·E_face` as the gradient steepens (the
    /// optically-thin / `R → ∞` free-streaming limit, ADR-0006). This is the per-face guarantee that
    /// keeps the assembled scheme from transporting radiation superluminally in thin gas.
    #[test]
    fn free_streaming_caps_flux_at_ce() {
        let c = 3.0; // arbitrary light speed; the cap is |F| ≤ c·E_face for any c
        let e_lo = 1.0;
        let e_hi = 0.0;
        let e_face = 0.5 * (e_lo + e_hi);
        let cap = c * e_face;

        // Across thin→thick media and a fixed jump, the flux never exceeds the c·E cap.
        for &chi_r in &[1e-4, 1e-2, 1.0, 10.0, 100.0] {
            for &spacing in &[1e-2, 0.1, 1.0, 10.0] {
                let d = face_diffusion(Limiter::LevermorePomraning, c, chi_r, e_lo, e_hi, spacing);
                let flux = (d * (e_lo - e_hi) / spacing).abs();
                assert!(
                    flux <= cap * (1.0 + 1e-12),
                    "superluminal flux: |F|={flux} > c·E={cap} (χ_R={chi_r}, spacing={spacing})"
                );
            }
        }

        // Free-streaming saturation: as the gradient steepens (R → ∞), |F| → c·E_face.
        let chi_r = 1e-6; // very thin
        let spacing = 1.0;
        let d = face_diffusion(Limiter::LevermorePomraning, c, chi_r, e_lo, e_hi, spacing);
        let flux = (d * (e_lo - e_hi) / spacing).abs();
        assert_relative_eq!(flux, cap, max_relative = 1e-4);
    }

    // --- B4b: cold black absorber wall (ADR-0005, radiative loss channel 1a) ---

    /// Total radiation energy `Σ E_j dx_j` (matter decoupled, so this is the only radiation bucket).
    fn rad_energy(e_rad: &[f64], dx: &[f64]) -> f64 {
        e_rad.iter().zip(dx).map(|(&e, &d)| e * d).sum()
    }

    /// The cold black absorber wall is just `RadBc::Marshak(0.0)`: with no incident radiation the
    /// Marshak current is the pure outflow `F = −(c/2)E_surface` — radiation streams into the wall
    /// and nothing comes back (ADR-0005: the radiative loss *is* the flux that reaches the wall).
    ///
    /// With matter decoupled (`χ_P = 0`, no emission/absorption) and conservative interior diffusion,
    /// the only sink is the wall, so the discrete energy balance is **exact**:
    /// `E_total⁺ = E_total⁻ − dt·(c/2)·E_wall⁺`. We verify that accounting every step, plus that the
    /// energy decreases monotonically, stays non-negative, and a meaningful fraction is actually
    /// drained (so the check is not passing on ~zero flux).
    #[test]
    fn cold_black_absorber_drains_energy_at_half_c_e() {
        let n = 100;
        let (c, a) = (2.0, 1.0);
        let consts = RadConstants { c, a };
        let dx = vec![0.01; n];
        let center_spacing = vec![0.01; n - 1];
        let temp = vec![0.0; n]; // T = 0 and χ_P = 0 ⇒ matter fully decoupled (aT⁴ = 0, k = 0)
        let cv_vol = vec![1.0; n];
        let chi_planck = vec![0.0; n];
        let chi_ross = vec![1.0; n];
        let mut e_rad = vec![1.0; n]; // uniform radiation; the left wall drains it from step one
        let dt = 1e-3;

        let initial = rad_energy(&e_rad, &dx);
        for _ in 0..1000 {
            let medium = Medium {
                dx: &dx,
                center_spacing: &center_spacing,
                temp: &temp,
                cv_vol: &cv_vol,
                chi_planck: &chi_planck,
                chi_ross: &chi_ross,
                source: None,
            };
            let before = rad_energy(&e_rad, &dx);
            let _ = fld_substep(
                &medium,
                &mut e_rad,
                RadBc::Marshak(0.0), // cold black absorber at x = 0
                RadBc::Reflecting,   // closed far end
                dt,
                consts,
                Limiter::Fick,
            );
            let after = rad_energy(&e_rad, &dx);

            // Exact flux accounting: the drop equals the implicit wall outflow dt·(c/2)·E_wall⁺.
            let wall_loss = dt * 0.5 * c * e_rad[0];
            assert_relative_eq!(
                before - after,
                wall_loss,
                max_relative = 1e-9,
                epsilon = 1e-14
            );
            assert!(after <= before, "energy must not grow through an absorber");
            assert!(
                e_rad.iter().all(|&e| e >= -1e-12),
                "radiation energy stays non-negative"
            );
        }
        let drained = (initial - rad_energy(&e_rad, &dx)) / initial;
        assert!(
            drained > 0.3,
            "absorber drained only {drained:.3} of the energy"
        );
    }

    /// Contrast / control: with both ends `Reflecting` and matter decoupled, the same diffusion is
    /// a closed system — total radiation energy is conserved to round-off. This isolates the loss in
    /// the absorber test above as the wall, not a diffusion artifact, and checks the reflecting BC.
    #[test]
    fn reflecting_walls_conserve_radiation_energy() {
        let n = 100;
        let consts = RadConstants { c: 2.0, a: 1.0 };
        let dx = vec![0.01; n];
        let center_spacing = vec![0.01; n - 1];
        let temp = vec![0.0; n];
        let cv_vol = vec![1.0; n];
        let chi_planck = vec![0.0; n];
        let chi_ross = vec![1.0; n];
        // A non-uniform blob so diffusion actually moves energy around between the closed walls.
        let mut e_rad: Vec<f64> = (0..n)
            .map(|j| {
                let x = (j as f64 + 0.5) / n as f64;
                (-(x - 0.5) * (x - 0.5) / 0.01).exp()
            })
            .collect();

        let total0 = rad_energy(&e_rad, &dx);
        for _ in 0..400 {
            let medium = Medium {
                dx: &dx,
                center_spacing: &center_spacing,
                temp: &temp,
                cv_vol: &cv_vol,
                chi_planck: &chi_planck,
                chi_ross: &chi_ross,
                source: None,
            };
            let _ = fld_substep(
                &medium,
                &mut e_rad,
                RadBc::Reflecting,
                RadBc::Reflecting,
                1e-3,
                consts,
                Limiter::Fick,
            );
            assert_relative_eq!(rad_energy(&e_rad, &dx), total0, max_relative = 1e-10);
        }
    }
}
