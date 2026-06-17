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
//! Each end carries a [`Boundary`]: a rigid [`Boundary::Wall`] (node held at `u = 0` —
//! reflecting; exact for the Sod tube and the rigid-walled standing wave of the convergence
//! test) or a [`Boundary::Free`] vacuum surface (driven outward by the interior pressure with
//! `p = 0` outside). The slug-into-wall bounce ([`Tube::slug`], [`Tube::run_bounce`]) pairs a
//! reflecting wall with a trailing free surface to measure the restitution `e_eff` (ADR-0001).

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

/// Boundary condition at an end node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Rigid reflecting wall: the end node is held fixed (`u = 0`).
    Wall,
    /// Free surface / vacuum: the end node is accelerated by the interior pressure with `p = 0`
    /// outside, so the gas re-expands outward into vacuum.
    Free,
}

/// Result of a slug-into-wall bounce (ADR-0001): the wall impulse and the restitution it implies.
#[derive(Debug, Clone, Copy)]
pub struct BounceResult {
    /// Time-integrated wall force `J_wall = ∫ P_wall dt`.
    pub wall_impulse: f64,
    /// Incident axial momentum magnitude `p_in` (the slug's initial momentum).
    pub incident_momentum: f64,
    /// Gas momentum still in flight when the run stopped (the rebound, signed away from wall).
    pub residual_momentum: f64,
    /// Effective restitution `e_eff = J_wall / p_in − 1` (ADR-0001).
    pub e_eff: f64,
    /// Peak wall force seen during the bounce (the tail guard stops at `10⁻³` of this).
    pub peak_wall_force: f64,
}

/// A 1D Lagrangian gas column on a staggered mesh.
#[derive(Debug, Clone)]
pub struct Tube {
    gamma: f64,
    viscosity: Viscosity,
    left: Boundary,
    right: Boundary,
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
            let v_left = vel[i.saturating_sub(1)];
            let v_right = vel[i.min(cells - 1)];
            *ui = 0.5 * (v_left + v_right);
        }
        let mut tube = Self {
            gamma,
            viscosity,
            left: Boundary::Wall,
            right: Boundary::Wall,
            x,
            u,
            mass,
            energy,
        };
        tube.enforce_wall_velocities();
        tube
    }

    /// Pin the velocity of any node that sits against a rigid [`Boundary::Wall`] to zero.
    fn enforce_wall_velocities(&mut self) {
        if self.left == Boundary::Wall {
            self.u[0] = 0.0;
        }
        if self.right == Boundary::Wall {
            let last = self.u.len() - 1;
            self.u[last] = 0.0;
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

    /// A finite cold gas slug on `x ∈ [0, 1]` moving toward a rigid wall at `x = 0`, with a free
    /// (vacuum) surface trailing at `x = 1` — the momentum-limit bounce harness (ADR-0001).
    ///
    /// Normalized to `ρ₀ = 1`, `v = 1`, so the incident Mach number `M = v / c₀` is set purely by
    /// the (cold) pressure `p₀ = ρ₀ v² / (γ M²)`. Lowering `M` warms the slug toward the elastic
    /// (acoustic) limit; raising it cools the slug toward the strong-shock ceiling.
    #[must_use]
    pub fn slug(cells: usize, mach: f64, gamma: f64) -> Self {
        let dx = 1.0 / cells as f64;
        let x: Vec<f64> = (0..=cells).map(|i| i as f64 * dx).collect();
        let p0 = 1.0 / (gamma * mach * mach); // c₀ = v/M = 1/M with ρ₀ = 1, v = 1
        let rho = vec![1.0; cells];
        let pressure = vec![p0; cells];
        let vel = vec![-1.0; cells]; // moving toward the wall at x = 0
        let mut tube = Self::new(
            x,
            &rho,
            &vel,
            &pressure,
            gamma,
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        tube.right = Boundary::Free;
        tube.enforce_wall_velocities();
        tube
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

    /// Ideal-gas sound speed in cell `j`, `c = sqrt(γ p / ρ)`. A vacuum/near-vacuum cell
    /// (`p ≤ 0` or `ρ ≤ 0`, reachable at the free surface) has no acoustic signal, so `c = 0`.
    fn sound_speed(&self, j: usize) -> f64 {
        let p = self.pressure(j);
        let rho = self.density(j);
        if p > 0.0 && rho > 0.0 {
            (self.gamma * p / rho).sqrt()
        } else {
            0.0
        }
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
    /// with node mass `m̄_i = ½(m_{j−1} + m_j)`. A [`Boundary::Wall`] end stays at zero
    /// acceleration (held fixed); a [`Boundary::Free`] end is driven by the one interior cell
    /// against vacuum (`P = 0` outside), so it accelerates outward.
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
        if self.left == Boundary::Free {
            accel[0] = -total_p[0] / (0.5 * self.mass[0]);
        }
        if self.right == Boundary::Free {
            let last = self.x.len() - 1;
            accel[last] = total_p[cells - 1] / (0.5 * self.mass[cells - 1]);
        }
        accel
    }

    /// CFL-limited timestep, `dt = CFL · min_j Δx_j / (c_j + |Δu_j|)`. The signal speed is the
    /// sound speed **plus** the cell's compression rate `|Δu_j| = |u_{j+1} − u_j|`: in a
    /// Lagrangian frame a cell can be crushed by the relative node motion as well as traversed by
    /// sound, and at supersonic (high-Mach) inflow the `|Δu|` term is what stops a node from
    /// overrunning a full cell width in one step and tangling the mesh.
    fn stable_dt(&self) -> f64 {
        let dt = (0..self.cells())
            .map(|j| {
                let signal = self.sound_speed(j) + (self.u[j + 1] - self.u[j]).abs();
                self.width(j) / signal
            })
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
        // The CFL above keeps the per-step compression modest, so the denominator stays well
        // away from its zero at `dV/V_new = −2/(γ−1)`; the `max(0)` floor is a positivity
        // safety net for strong expansion into vacuum, where the gas can cool past zero internal
        // energy numerically. It is never exercised by the smooth/shock interior tests.
        for j in 0..self.cells() {
            let v_new = 1.0 / self.density(j);
            let dv = v_new - v_old[j];
            let q = self.artificial_viscosity(j);
            let numer = self.energy[j] - (0.5 * p_old[j] + q) * dv;
            let denom = 1.0 + 0.5 * (self.gamma - 1.0) * dv / v_new;
            self.energy[j] = (numer / denom).max(0.0);
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

    /// Total axial momentum `Σ_i m̄_i u_i` carried by the nodes (boundary nodes own a half-cell).
    fn total_momentum(&self) -> f64 {
        let cells = self.cells();
        self.u
            .iter()
            .enumerate()
            .map(|(i, &ui)| {
                let m_left = if i == 0 { 0.0 } else { self.mass[i - 1] };
                let m_right = if i == cells { 0.0 } else { self.mass[i] };
                0.5 * (m_left + m_right) * ui
            })
            .sum()
    }

    /// Force the gas exerts on the rigid wall at `x = 0`: the total pressure `p + q` of cell 0.
    fn wall_force(&self) -> f64 {
        self.pressure(0) + self.artificial_viscosity(0)
    }

    /// Fire the slug at the wall and integrate the bounce until the wall force decays to `10⁻³`
    /// of its peak (ADR-0001's tail guard) or a safety step cap is hit, returning the wall
    /// impulse and the restitution it implies.
    ///
    /// The wall impulse is accumulated with the **trapezoidal** rule, which matches the
    /// velocity-Verlet momentum update exactly — so the conservation identity
    /// `J_wall == p_final − p_initial` holds to round-off (the elastic bookkeeping check).
    pub fn run_bounce(&mut self) -> BounceResult {
        let p_initial = self.total_momentum();
        let incident = p_initial.abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = self.wall_force();
        let max_steps = 400 * self.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            if force_old < 0.5 * peak {
                past_peak = true;
            }
            if past_peak && force_old < 1e-3 * peak {
                break;
            }
            let dt = self.stable_dt();
            self.step(dt);
            let force_new = self.wall_force();
            wall_impulse += 0.5 * dt * (force_old + force_new);
            force_old = force_new;
        }

        let residual = self.total_momentum();
        BounceResult {
            wall_impulse,
            incident_momentum: incident,
            residual_momentum: residual,
            e_eff: wall_impulse / incident - 1.0,
            peak_wall_force: peak,
        }
    }
}
