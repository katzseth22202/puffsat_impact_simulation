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

    /// One **coupled** backward-Euler conduction step over the gas *and* the solid solved as a single
    /// tridiagonal system (B-flux, ADR-0005). This is the fix for the inviscid kernel's missing
    /// gas-side thermal resistance: instead of pinning the interface to the bulk gas temperature (the
    /// over-draining `step_surface_temp` coupling), the gas gets its own conductivity `k_gas`, so the
    /// **interface temperature emerges from flux continuity** between the two conducting media.
    ///
    /// The unknowns are ordered `[gas_{n−1} … gas_0 | solid_0 … solid_{m−1}]` (gas cell 0 at the
    /// `x = 0` interface, increasing index toward the trailing free surface; solid `x` increasing into
    /// the plate), so physically adjacent cells are adjacent in the matrix. Every face uses the
    /// **series resistance of its two half-cells**, `G = 1 / (½w_L/k_L + ½w_R/k_R)`, which handles the
    /// non-uniform Lagrangian gas mesh, the variable `k_gas`, *and* the gas|solid material jump in one
    /// formula (it reduces to `k/dx` for the uniform solid interior, matching `step_surface_temp`).
    /// The gas far end is adiabatic (the vacuum/free surface); the solid deep end is held Dirichlet at
    /// `t_init` (the semi-infinite truncation, valid while `√(αt) ≪ depth`).
    ///
    /// Updates the solid temperatures in place and returns the per-gas-cell temperature change plus
    /// the new-time interface flux `q = G_iface·(T_gas0 − T_solid0)` — the gas's conductive loss
    /// (channel 2) over the step is `q·dt` (positive when heat flows from the hot gas into the plate).
    ///
    /// # Panics
    /// Panics unless every gas slice has the same length `≥ 1` and `dt > 0`.
    pub fn step_coupled(&mut self, gas: &GasConductionState<'_>, dt: f64) -> CoupledStep {
        let n_gas = gas.temp.len();
        assert!(n_gas >= 1, "need at least one gas cell");
        assert!(
            gas.dx.len() == n_gas && gas.cv_vol.len() == n_gas && gas.k_gas.len() == n_gas,
            "gas conduction arrays must all have the same length"
        );
        assert!(dt > 0.0, "dt must be positive");

        let n_solid = self.temp.len();
        let n = n_gas + n_solid;

        // Per global-cell width, conductivity, capacity/dt, and old temperature. Global index
        // k < n_gas is gas cell j = n_gas−1−k (gas reversed); k ≥ n_gas is solid cell k−n_gas.
        let rho_c_solid = self.k / self.alpha; // ρc of the solid = k/α
        let mut width = vec![0.0; n];
        let mut kk = vec![0.0; n];
        let mut cap = vec![0.0; n]; // C_i / dt
        let mut told = vec![0.0; n];
        for k in 0..n_gas {
            let j = n_gas - 1 - k;
            width[k] = gas.dx[j];
            kk[k] = gas.k_gas[j];
            cap[k] = gas.cv_vol[j] * gas.dx[j] / dt;
            told[k] = gas.temp[j];
        }
        for jj in 0..n_solid {
            let k = n_gas + jj;
            width[k] = self.dx;
            kk[k] = self.k;
            cap[k] = rho_c_solid * self.dx / dt;
            told[k] = self.temp[jj];
        }

        // Face conductances between consecutive global cells (series resistance of the two halves).
        let g_face: Vec<f64> = (0..n - 1)
            .map(|k| 1.0 / (0.5 * width[k] / kk[k] + 0.5 * width[k + 1] / kk[k + 1]))
            .collect();
        // Deep Dirichlet end: half a solid cell to the far boundary held at t_init.
        let g_deep = 2.0 * self.k / self.dx;

        // Assemble: −G_L·u_{i−1} + (C_i/dt + G_L + G_R)·u_i − G_R·u_{i+1} = (C_i/dt)·u_i^old.
        // Gas far end (i = 0) is adiabatic (no left face); solid deep end adds the Dirichlet face.
        let mut sub = vec![0.0; n];
        let mut sup = vec![0.0; n];
        let mut diag = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for i in 0..n {
            let gl = if i > 0 { g_face[i - 1] } else { 0.0 };
            let gr = if i < n - 1 { g_face[i] } else { 0.0 };
            sub[i] = -gl;
            sup[i] = -gr;
            diag[i] = cap[i] + gl + gr;
            rhs[i] = cap[i] * told[i];
        }
        diag[n - 1] += g_deep;
        rhs[n - 1] += g_deep * self.t_init;

        let u = thomas_solve(&sub, &diag, &sup, &rhs);

        self.temp.copy_from_slice(&u[n_gas..]);
        let gas_dtemp: Vec<f64> = (0..n_gas).map(|j| u[n_gas - 1 - j] - gas.temp[j]).collect();

        // New-time interface flux: gas_0 is global n_gas−1, solid_0 is global n_gas.
        let interface_flux = g_face[n_gas - 1] * (u[n_gas - 1] - u[n_gas]);
        CoupledStep {
            gas_dtemp,
            interface_flux,
        }
    }
}

/// The gas-side inputs to one [`Solid::step_coupled`] conduction step, borrowed from the current gas
/// state. All slices are per gas cell in kernel order (cell 0 at the `x = 0` wall interface): `dx`
/// the cell width, `temp` the temperature `T(ρ, e)`, `cv_vol = ρ c_v` the volumetric heat capacity,
/// and `k_gas` the thermal conductivity (ADR-0005).
#[derive(Debug, Clone, Copy)]
pub struct GasConductionState<'a> {
    /// Per-cell width.
    pub dx: &'a [f64],
    /// Per-cell temperature.
    pub temp: &'a [f64],
    /// Per-cell volumetric heat capacity `ρ c_v`.
    pub cv_vol: &'a [f64],
    /// Per-cell gas thermal conductivity `k_gas`.
    pub k_gas: &'a [f64],
}

/// Outcome of one [`Solid::step_coupled`] coupled conduction step.
#[derive(Debug, Clone)]
pub struct CoupledStep {
    /// Per-gas-cell temperature change `Tⁿ⁺¹ − Tⁿ`, in the same (kernel) order as the input gas cells.
    pub gas_dtemp: Vec<f64>,
    /// New-time interface heat flux into the solid `q = G_iface·(T_gas0 − T_solid0)` [W/m²]; the gas's
    /// conductive loss (channel 2) over the step is `q·dt`.
    pub interface_flux: f64,
}

#[cfg(test)]
mod tests {
    use super::{GasConductionState, Solid};
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

    /// Run a *uniform* conducting gas slab (constant `k_g`, `cv_vol_g`) against `solid` for `steps`
    /// coupled steps of size `dt`, applying the returned gas ΔT each step. Returns the final gas
    /// temperatures, the gas cell width, and the cumulative interface impulse `Σ q·dt` (the gas's
    /// total conductive loss). The solid is mutated in place (read it via `solid.temperature`).
    #[allow(clippy::too_many_arguments)] // a test fixture spelling out the contact-problem parameters
    fn run_contact(
        solid: &mut Solid,
        n_gas: usize,
        depth_g: f64,
        t_gas: f64,
        cv_vol_g: f64,
        k_g: f64,
        dt: f64,
        steps: usize,
    ) -> (Vec<f64>, f64, f64) {
        let dx_g = depth_g / n_gas as f64;
        let dx = vec![dx_g; n_gas];
        let cv = vec![cv_vol_g; n_gas];
        let kg = vec![k_g; n_gas];
        let mut temp = vec![t_gas; n_gas];
        let mut impulse = 0.0;
        for _ in 0..steps {
            let gas = GasConductionState {
                dx: &dx,
                temp: &temp,
                cv_vol: &cv,
                k_gas: &kg,
            };
            let step = solid.step_coupled(&gas, dt);
            for (t, d) in temp.iter_mut().zip(step.gas_dtemp.iter()) {
                *t += d;
            }
            impulse += step.interface_flux * dt;
        }
        (temp, dx_g, impulse)
    }

    /// `√(αt)` for the contact-problem parameters, shared by the oracle tests.
    const ALPHA_S: f64 = 0.5;
    const K_S: f64 = 10.0;
    const ALPHA_G: f64 = 0.5; // same α as the solid, so the effusivity ratio is the conductivity ratio
    const K_G: f64 = 1.0; // ⇒ e_s/e_g = 10, the cold high-effusivity plate the gas meets

    /// L∞ of the two-media contact solution (gas + solid sides) against the analytic erf profiles
    /// about the effusivity-weighted interface `T_i`. `t_final` must keep `√(αt) ≪ depth`.
    fn contact_linf(n: usize, dt: f64, t_final: f64) -> f64 {
        let (t_gas, t_solid, depth) = (1.0, 0.0, 1.0);
        let cv_vol_g = K_G / ALPHA_G; // ρc_v of the gas
        let e_g = K_G / ALPHA_G.sqrt();
        let e_s = K_S / ALPHA_S.sqrt();
        let t_i = (e_g * t_gas + e_s * t_solid) / (e_g + e_s);

        // SAFE: t_final, dt > 0 ⇒ a small non-negative integer; no truncation or sign loss.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (t_final / dt).round() as usize;
        let mut solid = Solid::new(n, depth, t_solid, ALPHA_S, K_S);
        let (gas_temp, dx_g, _) =
            run_contact(&mut solid, n, depth, t_gas, cv_vol_g, K_G, dt, steps);

        let mut linf = 0.0_f64;
        for (j, &tj) in gas_temp.iter().enumerate() {
            let x = (j as f64 + 0.5) * dx_g;
            let exact = (t_gas - t_i).mul_add(erf(x / (2.0 * (ALPHA_G * t_final).sqrt())), t_i);
            linf = linf.max((tj - exact).abs());
        }
        for j in 0..solid.cells() {
            let x = solid.center(j);
            let exact = (t_solid - t_i).mul_add(erf(x / (2.0 * (ALPHA_S * t_final).sqrt())), t_i);
            linf = linf.max((solid.temperature(j) - exact).abs());
        }
        linf
    }

    /// Two-semi-infinite-media contact (the B-flux acceptance): a uniform conducting gas against a
    /// 10×-more-effusive solid. The interface jumps to the effusivity-weighted
    /// `T_i = (e_g·T_g + e_s·T_s)/(e_g + e_s)` and each side relaxes as an erf profile about it —
    /// verifying flux continuity and the emergent interface temperature end to end.
    #[test]
    fn coupled_two_media_contact_matches_effusivity_interface_and_erf() {
        let (t_gas, t_solid, depth) = (1.0, 0.0, 1.0);
        let (dt, t_final) = (1e-5_f64, 0.01_f64); // √(αt) ≈ 0.0707 ≪ depth ⇒ semi-infinite both sides
        let cv_vol_g = K_G / ALPHA_G;
        let e_g = K_G / ALPHA_G.sqrt();
        let e_s = K_S / ALPHA_S.sqrt();
        let t_i = (e_g * t_gas + e_s * t_solid) / (e_g + e_s);

        assert!(
            contact_linf(400, dt, t_final) < 1e-2,
            "contact profile L∞ too large"
        );

        // The effusivity weighting explicitly: extrapolate each side's first two cells to x = 0 → T_i.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (t_final / dt).round() as usize;
        let mut solid = Solid::new(400, depth, t_solid, ALPHA_S, K_S);
        let (gas_temp, _, impulse) =
            run_contact(&mut solid, 400, depth, t_gas, cv_vol_g, K_G, dt, steps);
        let gas_iface = 1.5f64.mul_add(gas_temp[0], -0.5 * gas_temp[1]);
        let solid_iface = 1.5f64.mul_add(solid.temperature(0), -0.5 * solid.temperature(1));
        assert!(
            (gas_iface - t_i).abs() < 5e-3 && (solid_iface - t_i).abs() < 5e-3,
            "interface T_i={t_i}: gas {gas_iface}, solid {solid_iface}"
        );
        assert!(
            impulse > 0.0,
            "heat must flow from the hot gas into the cold plate"
        );
    }

    /// The contact-solution error shrinks under joint refinement (`dt ∝ dx²`, so both the 1st-order
    /// backward-Euler and 2nd-order spatial errors scale as `dx²`) — the order-of-accuracy gate.
    #[test]
    fn coupled_contact_converges_under_refinement() {
        let coarse = contact_linf(200, 4e-5, 0.01);
        let fine = contact_linf(400, 1e-5, 0.01);
        assert!(
            fine < coarse,
            "error grew under refinement: {coarse:e} -> {fine:e}"
        );
        let ratio = coarse / fine;
        assert!(
            ratio > 2.5,
            "convergence ratio {ratio:.2} below the expected ~2nd order"
        );
    }

    /// Energy closure: with the gas far end adiabatic, the gas's total energy loss equals the
    /// cumulative interface impulse `Σ q·dt` exactly (backward-Euler conservation), and matches the
    /// solid's energy gain while the deep end stays cold (shallow penetration).
    #[test]
    fn coupled_conduction_conserves_energy() {
        let (t_gas, t_solid, depth) = (1.0, 0.0, 1.0);
        let (dt, n) = (1e-5_f64, 400usize);
        let cv_vol_g = K_G / ALPHA_G;
        let mut solid = Solid::new(n, depth, t_solid, ALPHA_S, K_S);
        let (gas_temp, dx_g, impulse) =
            run_contact(&mut solid, n, depth, t_gas, cv_vol_g, K_G, dt, 500);

        let gas_lost: f64 = gas_temp
            .iter()
            .map(|&t| cv_vol_g * dx_g * (t_gas - t))
            .sum();
        let rho_c_solid = K_S / ALPHA_S;
        let dx_s = depth / n as f64;
        let solid_gained: f64 = (0..solid.cells())
            .map(|j| rho_c_solid * dx_s * (solid.temperature(j) - t_solid))
            .sum();

        // Gas only exchanges across the interface (far end adiabatic) ⇒ exact to solver round-off.
        assert!(
            (gas_lost - impulse).abs() / impulse < 1e-9,
            "gas loss {gas_lost:e} vs interface impulse {impulse:e}"
        );
        // Solid gain ≈ gas loss while the deep Dirichlet end leaks negligibly.
        assert!(
            (solid_gained - gas_lost).abs() / gas_lost < 1e-2,
            "solid gain {solid_gained:e} vs gas loss {gas_lost:e}"
        );
    }

    /// The gas conductivity is the knob the whole rung adds: a very conductive gas presents **no
    /// internal thermal resistance** (it stays isothermal — the over-drain regime, where the plate is
    /// driven by the bulk gas temperature with nothing to throttle it), while a very insulating gas
    /// chokes the interface flux to ≈ 0 and leaves the plate cold. The physical `k_gas` lives between.
    #[test]
    fn coupled_conduction_gas_conductivity_limits_bracket_the_interface() {
        let (t_gas, t_solid, depth) = (1.0, 0.0, 1.0);
        let (dt, n, steps) = (1e-5_f64, 200usize, 200usize);
        let cv_vol_g = 1.0; // fixed ρc_v; vary k_gas to swing the gas effusivity

        // Huge k_gas: the gas develops no internal gradient (infinite gas conductance ⇒ isothermal),
        // so there is no gas-side resistance — the plate sees the bulk gas temperature directly.
        let mut solid_hi = Solid::new(n, depth, t_solid, ALPHA_S, K_S);
        let (gas_hi, _, impulse_hi) =
            run_contact(&mut solid_hi, n, depth, t_gas, cv_vol_g, 1e6, dt, steps);
        let spread = gas_hi.iter().copied().fold(f64::MIN, f64::max)
            - gas_hi.iter().copied().fold(f64::MAX, f64::min);
        assert!(
            spread < 1e-3,
            "huge k_gas: gas should be isothermal (no internal resistance), spread {spread:e}"
        );

        // Tiny k_gas: the gas insulates, so almost no heat reaches the plate and it stays cold.
        let mut solid_lo = Solid::new(n, depth, t_solid, ALPHA_S, K_S);
        let (_, _, impulse_lo) =
            run_contact(&mut solid_lo, n, depth, t_gas, cv_vol_g, 1e-6, dt, steps);
        assert!(
            impulse_lo < 1e-3 * impulse_hi && solid_lo.temperature(0) < t_solid + 1e-3,
            "tiny k_gas: impulse {impulse_lo:e} (vs {impulse_hi:e}) and plate must stay cold"
        );
    }
}
