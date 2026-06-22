//! Transient 1D heat conduction in the semi-infinite conducting solid behind the rigid wall
//! (ADR-0005, Rung-B B4). The wall's **conductive-to-plate** loss (loss channel 2, ADR-0016) is the
//! heat flux that crosses the gas/solid interface into this solid. We model the solid as a 1D
//! **backward-Euler implicit** heat-conduction mesh on the same Thomas solver the radiation step
//! uses ([`crate::radiation::thomas_solve`]); a fixed deep far boundary held at the initial
//! temperature stands in for the semi-infinite extent — valid while the thermal penetration depth
//! `√(αt)` stays shallow compared to the mesh depth (the same `√(Dt) ≪ L` guard as the radiation
//! erfc test). The meshed solid (vs an effusivity convolution) was chosen so the same machinery
//! later carries the two-layer ablator-on-SiC stack and the Phase-2 ablating wall.
//!
//! Built test-first against the textbook **step-temperature** solution: a solid initially at `T₀`
//! whose surface is suddenly held at `Tₛ` relaxes to `T(x,t) = T₀ + (Tₛ−T₀)·erfc(x / 2√(αt))`, and
//! the interface heat flux is `q(t) = e (Tₛ−T₀) / √(πt)` with **effusivity** `e = √(kρc) = k/√α` —
//! the single material parameter ADR-0005 says sets the conductive loss.

use crate::radiation::thomas_solve;

/// A semi-infinite conducting solid discretized as a uniform 1D cell-centered mesh, with the
/// gas/solid interface at `x = 0` and `x` increasing into the solid. Node `j` sits at the cell
/// center `x = (j + ½)·dx`. Holds the thermal diffusivity `α = k/(ρc)` and the conductivity `k`;
/// together these fix the effusivity `e = k/√α = √(kρc)`.
#[derive(Debug, Clone)]
pub struct Solid {
    /// Cell-centered node temperatures.
    temp: Vec<f64>,
    /// Initial (and deep far-boundary) temperature `T₀`.
    t_init: f64,
    /// Cell width.
    dx: f64,
    /// Thermal diffusivity `α = k/(ρc)`.
    alpha: f64,
    /// Conductivity `k` (sets the interface flux and, with `α`, the effusivity).
    k: f64,
}

impl Solid {
    /// A solid of `cells` cells spanning depth `depth` (from the interface), uniformly at `t_init`,
    /// with thermal diffusivity `alpha` and conductivity `k`.
    ///
    /// # Panics
    /// Panics unless `cells ≥ 1` and `depth`, `alpha`, `k` are all positive.
    #[must_use]
    pub fn new(cells: usize, depth: f64, t_init: f64, alpha: f64, k: f64) -> Self {
        assert!(cells >= 1, "solid needs at least one cell");
        assert!(
            depth > 0.0 && alpha > 0.0 && k > 0.0,
            "depth, alpha, k must be positive"
        );
        Self {
            temp: vec![t_init; cells],
            t_init,
            dx: depth / cells as f64,
            alpha,
            k,
        }
    }

    /// Effusivity `e = √(kρc) = k/√α` — the material parameter that sets the conductive loss.
    #[must_use]
    pub fn effusivity(&self) -> f64 {
        self.k / self.alpha.sqrt()
    }

    /// Number of mesh cells.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.temp.len()
    }

    /// Cell-center depth `x = (j + ½)·dx` of node `j` (distance into the solid from the interface).
    #[must_use]
    pub fn center(&self, j: usize) -> f64 {
        (j as f64 + 0.5) * self.dx
    }

    /// Temperature of node `j`.
    #[must_use]
    pub fn temperature(&self, j: usize) -> f64 {
        self.temp[j]
    }

    /// Advance one backward-Euler step of size `dt` with the interface (`x = 0`) held at
    /// `surface_temp` (Dirichlet, half a cell from node 0) and the deep far boundary held at the
    /// initial temperature (the semi-infinite truncation). Returns the interface heat flux **into**
    /// the solid, `q = k·(Tₛ − T₀cell)/(dx/2)` — the per-step conductive loss the gas sees (positive
    /// when heat flows from the hot gas into the cooler solid).
    pub fn step_surface_temp(&mut self, surface_temp: f64, dt: f64) -> f64 {
        let n = self.temp.len();
        let r = self.alpha * dt / (self.dx * self.dx);

        // Backward-Euler heat equation: −r·T[j−1] + (1+2r)·T[j] − r·T[j+1] = Tⁿ[j]. Both ends are
        // Dirichlet at the half-cell face (conductance 2r): the interface at `surface_temp`, the
        // deep far end at the initial temperature (the semi-infinite truncation).
        let mut sub = vec![-r; n];
        let mut sup = vec![-r; n];
        let mut diag = vec![2.0f64.mul_add(r, 1.0); n];
        let mut rhs = self.temp.clone();

        sub[0] = 0.0;
        diag[0] = 3.0f64.mul_add(r, 1.0);
        rhs[0] += 2.0 * r * surface_temp;

        sup[n - 1] = 0.0;
        diag[n - 1] = 3.0f64.mul_add(r, 1.0);
        rhs[n - 1] += 2.0 * r * self.t_init;

        self.temp = thomas_solve(&sub, &diag, &sup, &rhs);

        // Interface flux into the solid from the surface half-cell gradient: q = k·(Tₛ − T[0])/(dx/2).
        2.0 * self.k * (surface_temp - self.temp[0]) / self.dx
    }
}

#[cfg(test)]
mod tests {
    use super::Solid;
    use std::f64::consts::PI;

    /// Abramowitz–Stegun 7.1.26 — `erf` to ~1e-7, plenty for an analytic reference.
    fn erf(x: f64) -> f64 {
        let t = 1.0 / 0.327_591_1f64.mul_add(x.abs(), 1.0);
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

    /// Step-temperature semi-infinite solid (the B4a acceptance): with the surface held at `Tₛ` and
    /// the solid initially at `T₀`, the profile must match `T₀ + (Tₛ−T₀)·erfc(x / 2√(αt))` and the
    /// interface flux must match `e·(Tₛ−T₀)/√(πt)` with `e = √(kρc) = k/√α`. This verifies both the
    /// implicit conduction machinery and the effusivity wiring end to end.
    #[test]
    fn step_temperature_matches_erfc_and_effusivity_flux() {
        let (cells, depth) = (400usize, 1.0);
        let (alpha, k) = (0.5, 2.0); // ρc = k/α = 4, so e = k/√α = √(kρc) = √8 = 2.828…
        let (t0, ts) = (0.0, 1.0);
        let dt = 1e-5_f64;
        let t_final = 0.01_f64; // √(αt) ≈ 0.07 ≪ depth ⇒ the far boundary stays at T₀

        let mut solid = Solid::new(cells, depth, t0, alpha, k);
        assert!(
            (solid.effusivity() - 8.0_f64.sqrt()).abs() < 1e-12,
            "effusivity wiring: e = k/√α = √8"
        );

        // SAFE: t_final, dt > 0 ⇒ a small non-negative integer; no truncation or sign loss.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (t_final / dt).round() as usize;
        let mut flux = 0.0;
        for _ in 0..steps {
            flux = solid.step_surface_temp(ts, dt);
        }

        // Profile: erfc, affine-shifted by T₀ and scaled by (Tₛ−T₀).
        let mut linf = 0.0_f64;
        for j in 0..solid.cells() {
            let xi = solid.center(j) / (2.0 * (alpha * t_final).sqrt());
            let exact = (ts - t0).mul_add(erfc(xi), t0);
            linf = linf.max((solid.temperature(j) - exact).abs());
        }
        assert!(linf < 5e-3, "step-temp profile L∞ = {linf:e}");

        // Interface flux: effusivity·ΔT/√(πt). The half-cell gradient is first-order, but erfc is
        // nearly linear at the surface, so a few-percent tolerance is honest at this resolution.
        let q_exact = solid.effusivity() * (ts - t0) / (PI * t_final).sqrt();
        let rel = (flux - q_exact).abs() / q_exact;
        assert!(
            rel < 0.03,
            "interface flux {flux:e} vs analytic {q_exact:e} (rel {rel:e})"
        );
    }
}
