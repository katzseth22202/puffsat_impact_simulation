//! Staggered-grid Lagrangian hydrodynamics with von Neumann–Richtmyer artificial viscosity
//! (ADR-0022), 1D planar, ideal-gas EOS `p = (γ−1) ρ e`.
//!
//! # Discretization
//!
//! The mesh is **staggered**: `N+1` nodes carry position `x_i` and velocity `u_i`; the `N`
//! cells between them carry a (Lagrangian-conserved) mass `m_j`, specific internal energy
//! `e_j`, and derived density `ρ_j = m_j / (x_{j+1} − x_j)` and pressure `p_j`. Shocks are
//! captured by an artificial-viscosity pressure `q_j`, active only in compression — no Riemann
//! solver appears anywhere in the kernel (ADR-0022).
//!
//! # Time integration
//!
//! One step is **velocity Verlet** (kick–drift–kick), 2nd-order in time, so that together with
//! the 2nd-order-in-space staggered differencing the scheme converges at rate 2 in smooth flow
//! (the convergence test). With acceleration `a = −∂(p+q)/∂m` at the nodes:
//! 1. **half-kick** `u ← u + ½dt·a(tⁿ)`;
//! 2. **drift** `x ← x + dt·u`; recompute `ρ`;
//! 3. **energy** update from `de = −(p̄ + q) dV` with `p̄` time-centered, solved implicitly for
//!    the ideal gas (stable, energy-conserving), giving `pⁿ⁺¹`;
//! 4. **half-kick** `u ← u + ½dt·a(tⁿ⁺¹)`.
//!
//! Endpoints are held fixed (`u = 0`); that is exact for the Sod tube and for a standing
//! acoustic wave in a rigid-walled tube (the convergence test). Reflecting/vacuum wall
//! boundaries for the momentum-limit tests (ADR-0001) are a later increment.

use crate::Primitive;

/// CFL number for the explicit timestep.
const CFL: f64 = 0.4;

/// Artificial-viscosity coefficients. The **quadratic** term damps strong shocks and is
/// `O(Δx²)` in smooth flow (it preserves 2nd-order accuracy); the **linear** term suppresses
/// post-shock oscillations but is `O(Δx)` (it degrades accuracy to 1st order), so smooth
/// order-of-accuracy tests use [`Viscosity::QUADRATIC_ONLY`].
#[derive(Debug, Clone, Copy)]
pub struct Viscosity {
    /// Quadratic coefficient `c_q`.
    pub quadratic: f64,
    /// Linear coefficient `c_l`.
    pub linear: f64,
}

impl Viscosity {
    /// Standard von Neumann–Richtmyer coefficients (quadratic + linear), for shock problems.
    pub const VON_NEUMANN_RICHTMYER: Self = Self {
        quadratic: 2.0,
        linear: 0.5,
    };
    /// Quadratic only — for smooth flow where the linear term would cap the convergence rate.
    pub const QUADRATIC_ONLY: Self = Self {
        quadratic: 2.0,
        linear: 0.0,
    };
}

/// A 1D Lagrangian gas column on a staggered mesh.
#[derive(Debug, Clone)]
pub struct Tube {
    gamma: f64,
    viscosity: Viscosity,
    /// Node positions, length `N + 1`.
    x: Vec<f64>,
    /// Node velocities, length `N + 1`.
    u: Vec<f64>,
    /// Cell masses (conserved), length `N`.
    mass: Vec<f64>,
    /// Cell specific internal energy, length `N`.
    energy: Vec<f64>,
}

impl Tube {
    /// Build a tube from cell-centered primitive initial conditions on the node grid `x`
    /// (length `cells + 1`). All `cells` slices share that node grid.
    ///
    /// # Panics
    /// Panics if `x.len() != rho.len() + 1` (one more node than cells).
    #[must_use]
    pub fn new(
        x: Vec<f64>,
        rho: &[f64],
        vel: &[f64],
        pressure: &[f64],
        gamma: f64,
        viscosity: Viscosity,
    ) -> Self {
        let cells = rho.len();
        assert_eq!(x.len(), cells + 1, "need one more node than cells");
        let mass: Vec<f64> = (0..cells).map(|j| rho[j] * (x[j + 1] - x[j])).collect();
        let energy: Vec<f64> = (0..cells)
            .map(|j| pressure[j] / ((gamma - 1.0) * rho[j]))
            .collect();
        // Node velocities: average of adjacent cell velocities; ends take their neighbor.
        let nodes = cells + 1;
        let mut u = vec![0.0; nodes];
        for (i, ui) in u.iter_mut().enumerate() {
            let left = vel[i.saturating_sub(1)];
            let right = vel[i.min(cells - 1)];
            *ui = 0.5 * (left + right);
        }
        Self {
            gamma,
            viscosity,
            x,
            u,
            mass,
            energy,
        }
    }

    /// The standard Sod shock tube on `x ∈ [0, 1]` with `cells` cells: a diaphragm at `x = 0.5`
    /// separating `(ρ,u,p) = (1, 0, 1)` on the left from `(0.125, 0, 0.1)` on the right.
    #[must_use]
    pub fn sod(cells: usize, gamma: f64) -> Self {
        let dx = 1.0 / cells as f64;
        let x: Vec<f64> = (0..=cells).map(|i| i as f64 * dx).collect();
        let mut rho = vec![0.0; cells];
        let mut pressure = vec![0.0; cells];
        for j in 0..cells {
            let center = (j as f64 + 0.5) * dx;
            if center < 0.5 {
                rho[j] = 1.0;
                pressure[j] = 1.0;
            } else {
                rho[j] = 0.125;
                pressure[j] = 0.1;
            }
        }
        let vel = vec![0.0; cells];
        Self::new(
            x,
            &rho,
            &vel,
            &pressure,
            gamma,
            Viscosity::VON_NEUMANN_RICHTMYER,
        )
    }

    /// Number of cells.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.mass.len()
    }

    /// Density of cell `j`, `ρ_j = m_j / (x_{j+1} − x_j)`.
    #[must_use]
    pub fn density(&self, j: usize) -> f64 {
        self.mass[j] / (self.x[j + 1] - self.x[j])
    }

    /// Pressure of cell `j` from the ideal-gas EOS, `p = (γ−1) ρ e`.
    #[must_use]
    pub fn pressure(&self, j: usize) -> f64 {
        (self.gamma - 1.0) * self.density(j) * self.energy[j]
    }

    /// Cell-centered velocity (average of the two bounding node velocities).
    #[must_use]
    pub fn velocity(&self, j: usize) -> f64 {
        0.5 * (self.u[j] + self.u[j + 1])
    }

    /// Current (Lagrangian-moved) center of cell `j`.
    #[must_use]
    pub fn center(&self, j: usize) -> f64 {
        0.5 * (self.x[j] + self.x[j + 1])
    }

    /// Current width of cell `j`, `x_{j+1} − x_j` (a quadrature weight for cell-centered fields).
    #[must_use]
    pub fn width(&self, j: usize) -> f64 {
        self.x[j + 1] - self.x[j]
    }

    /// The cell-centered primitive state of cell `j`.
    #[must_use]
    pub fn primitive(&self, j: usize) -> Primitive {
        Primitive::new(self.density(j), self.velocity(j), self.pressure(j))
    }

    /// Ideal-gas sound speed in cell `j`.
    fn sound_speed(&self, j: usize) -> f64 {
        (self.gamma * self.pressure(j) / self.density(j)).sqrt()
    }

    /// Artificial-viscosity pressure of cell `j`: quadratic + linear, active only under
    /// compression (`Δu = u_{j+1} − u_j < 0`), else zero.
    fn artificial_viscosity(&self, j: usize) -> f64 {
        let du = self.u[j + 1] - self.u[j];
        if du < 0.0 {
            let rho = self.density(j);
            rho * (self.viscosity.quadratic * du * du
                - self.viscosity.linear * self.sound_speed(j) * du)
        } else {
            0.0
        }
    }

    /// Nodal accelerations `a_i = −(P_j − P_{j−1}) / m̄_i` from the total pressure `P = p + q`,
    /// with node mass `m̄_i = ½(m_{j−1} + m_j)`. Endpoints are fixed, so their acceleration is 0.
    fn node_accelerations(&self) -> Vec<f64> {
        let cells = self.cells();
        let total_p: Vec<f64> = (0..cells)
            .map(|j| self.pressure(j) + self.artificial_viscosity(j))
            .collect();
        let mut accel = vec![0.0; self.x.len()];
        for i in 1..cells {
            let node_mass = 0.5 * (self.mass[i - 1] + self.mass[i]);
            accel[i] = -(total_p[i] - total_p[i - 1]) / node_mass;
        }
        accel
    }

    /// CFL-limited timestep over all cells, `dt = CFL · min_j Δx_j / c_j`.
    fn stable_dt(&self) -> f64 {
        let dt = (0..self.cells())
            .map(|j| (self.x[j + 1] - self.x[j]) / self.sound_speed(j))
            .fold(f64::INFINITY, f64::min);
        CFL * dt
    }

    /// Advance one step of size `dt` with velocity Verlet (kick–drift–kick).
    fn step(&mut self, dt: f64) {
        // 1. Half-kick to uⁿ⁺¹ᐟ²; endpoints have zero acceleration so stay fixed.
        let accel = self.node_accelerations();
        for (ui, ai) in self.u.iter_mut().zip(accel.iter()) {
            *ui += 0.5 * dt * ai;
        }

        // 2. Drift the mesh, remembering the time-n specific volume and pressure.
        let v_old: Vec<f64> = (0..self.cells()).map(|j| 1.0 / self.density(j)).collect();
        let p_old: Vec<f64> = (0..self.cells()).map(|j| self.pressure(j)).collect();
        for (xi, ui) in self.x.iter_mut().zip(self.u.iter()) {
            *xi += dt * ui;
        }

        // 3. Implicit ideal-gas energy update, `de = −(p̄ + q) dV`, p̄ = ½(p_old + p_new):
        //    e_new = [e_old − (½p_old + q)(V_new − V_old)] / [1 + ½(γ−1)(V_new − V_old)/V_new].
        for j in 0..self.cells() {
            let v_new = 1.0 / self.density(j);
            let dv = v_new - v_old[j];
            let q = self.artificial_viscosity(j);
            let numer = self.energy[j] - (0.5 * p_old[j] + q) * dv;
            let denom = 1.0 + 0.5 * (self.gamma - 1.0) * dv / v_new;
            self.energy[j] = numer / denom;
        }

        // 4. Half-kick to uⁿ⁺¹ using the updated (time-n+1) pressures.
        let accel = self.node_accelerations();
        for (ui, ai) in self.u.iter_mut().zip(accel.iter()) {
            *ui += 0.5 * dt * ai;
        }
    }

    /// Advance the solution to `t_end`, choosing a CFL-limited step each time and clipping the
    /// final step to land exactly on `t_end`.
    pub fn run_to(&mut self, t_end: f64) {
        let mut t = 0.0;
        while t < t_end {
            let dt = self.stable_dt().min(t_end - t);
            self.step(dt);
            t += dt;
        }
    }
}
