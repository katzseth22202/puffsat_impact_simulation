//! Staggered-grid Lagrangian hydrodynamics with von Neumann–Richtmyer artificial viscosity
//! (ADR-0022), 1D planar. The equation of state is pluggable (ADR-0022): the kernel is generic
//! over [`Eos`], with rung A's [`IdealGas`] (`p = (γ−1) ρ e`) and rung B's tabulated
//! [`TableEos`] both satisfying the same bare `p(ρ, e)` interface.
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
//!    the (possibly tabulated) EOS by Newton iteration, giving `pⁿ⁺¹`;
//! 4. **half-kick** `u ← u + ½dt·a(tⁿ⁺¹)`.
//!
//! Each end carries a [`Boundary`]: a rigid [`Boundary::Wall`] (node held at `u = 0` —
//! reflecting; exact for the Sod tube and the rigid-walled standing wave of the convergence
//! test) or a [`Boundary::Free`] vacuum surface (driven outward by the interior pressure with
//! `p = 0` outside). The slug-into-wall bounce ([`Tube::slug`], [`Tube::run_bounce`]) pairs a
//! reflecting wall with a trailing free surface to measure the restitution `e_eff` (ADR-0001).

use crate::Primitive;
use crate::conduction::{GasConductionState, Solid};
use crate::eos::{Eos, IdealGas, TableEos};
use crate::radiation::{Limiter, Medium, RadBc, RadConstants, fld_substep};

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

/// A 1D Lagrangian gas column on a staggered mesh, carrying its equation of state `E`.
#[derive(Debug, Clone)]
pub struct Tube<E: Eos> {
    eos: E,
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

impl Tube<IdealGas> {
    /// Build an **ideal-gas** tube from cell-centered primitive initial conditions on the node
    /// grid `x` (length `cells + 1`). Convenience wrapper over [`Tube::with_eos`] for rung A's
    /// `p = (γ−1) ρ e`; all `cells` slices share that node grid.
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
        Self::with_eos(x, rho, vel, pressure, IdealGas::new(gamma), viscosity)
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
}

impl<E: Eos> Tube<E> {
    /// Build a tube with an arbitrary [`Eos`] from cell-centered primitive initial conditions on
    /// the node grid `x` (length `cells + 1`). The initial `e` is seeded from the initial `p` via
    /// [`Eos::energy_from_pressure`].
    ///
    /// # Panics
    /// Panics if `x.len() != rho.len() + 1` (one more node than cells).
    #[must_use]
    pub fn with_eos(
        x: Vec<f64>,
        rho: &[f64],
        vel: &[f64],
        pressure: &[f64],
        eos: E,
        viscosity: Viscosity,
    ) -> Self {
        let cells = rho.len();
        assert_eq!(x.len(), cells + 1, "need one more node than cells");
        let mass: Vec<f64> = (0..cells).map(|j| rho[j] * (x[j + 1] - x[j])).collect();
        let energy: Vec<f64> = (0..cells)
            .map(|j| eos.energy_from_pressure(rho[j], pressure[j]))
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
            eos,
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

    /// Pressure of cell `j` from the EOS, `p(ρ_j, e_j)`. A fully cooled cell (`e ≤ 0`, reachable
    /// at the free surface) exerts no pressure.
    #[must_use]
    pub fn pressure(&self, j: usize) -> f64 {
        let e = self.energy[j];
        if e > 0.0 {
            self.eos.pressure(self.density(j), e)
        } else {
            0.0
        }
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

    /// EOS sound speed in cell `j`, `c_s(ρ_j, e_j)`. A vacuum/near-vacuum cell (`e ≤ 0` or
    /// `ρ ≤ 0`, reachable at the free surface) has no acoustic signal, so `c = 0`.
    fn sound_speed(&self, j: usize) -> f64 {
        let rho = self.density(j);
        let e = self.energy[j];
        if e > 0.0 && rho > 0.0 {
            self.eos.sound_speed(rho, e)
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

    /// Implicit time-centered energy update for one cell: solve
    /// `e_new = e_old − (½(p_old + p(ρ_new, e_new)) + q)·dV` for `e_new` by Newton iteration on
    /// `g(e) = e − e_old + (½(p_old + p(ρ_new, e)) + q)·dV`, with `g'(e) = 1 + ½ dV ∂p/∂e`. For an
    /// ideal gas `g` is linear, so the first step is exact (reproducing rung A); for a tabulated
    /// EOS a handful of steps converge. The result is floored at 0 — the positivity safety net
    /// for strong expansion into vacuum (never exercised by the smooth/shock interior tests).
    fn update_energy(eos: &E, rho_new: f64, e_old: f64, p_old: f64, q: f64, dv: f64) -> f64 {
        let mut e = e_old;
        for _ in 0..100 {
            let p = eos.pressure(rho_new, e);
            let g = e - e_old + (0.5 * (p_old + p) + q) * dv;
            let gp = 1.0 + 0.5 * dv * eos.dp_de(rho_new, e);
            if gp <= 0.0 {
                break; // CFL keeps the step away from the EOS-dependent singularity; bail safely.
            }
            let step = g / gp;
            e -= step;
            if step.abs() <= 1e-13 * (e.abs() + 1e-300) {
                break;
            }
        }
        e.max(0.0)
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

        // 3. Implicit time-centered energy update `de = −(p̄ + q) dV`, p̄ = ½(p_old + p_new),
        //    solved per cell for the (possibly tabulated) EOS.
        for j in 0..self.cells() {
            let rho_new = self.density(j);
            let v_new = 1.0 / rho_new;
            let dv = v_new - v_old[j];
            let q = self.artificial_viscosity(j);
            self.energy[j] =
                Self::update_energy(&self.eos, rho_new, self.energy[j], p_old[j], q, dv);
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

    /// Total (Lagrangian-conserved) gas mass `Σ_j m_j`. Conserved by the hydro step; a wall-sticking
    /// condensation sink (Rung C) is the only operator that removes it, so the drop equals the stuck
    /// mass (the mass-sink closure check exercises this).
    #[cfg(test)]
    fn total_mass(&self) -> f64 {
        self.mass.iter().sum()
    }

    /// Force the gas exerts on the rigid wall at `x = 0`: the total pressure `p + q` of cell 0.
    fn wall_force(&self) -> f64 {
        self.pressure(0) + self.artificial_viscosity(0)
    }

    /// Fire the slug at the wall and integrate the bounce until the wall force decays to `10⁻³`
    /// of its peak (ADR-0001's tail guard) or a safety step cap is hit, returning the wall
    /// impulse and the restitution it implies.
    ///
    /// The wall impulse is accumulated with the **trapezoidal** rule, mirroring the
    /// velocity-Verlet momentum update — so the conservation identity
    /// `J_wall == p_final − p_initial` holds to the scheme's `O(Δx)` consistency (≈4e-4 at 400
    /// cells, not round-off: `J_wall` samples the wall-pressure history independently of the
    /// momentum bookkeeping). This is the elastic bookkeeping check (ADR-0001).
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

    /// Fire the slug at an **idealized absorbing wall** — ADR-0001's dead-stick (`f → 0.5`) limit.
    /// Integrate the wall impulse only up to **stagnation**: the instant the net gas momentum
    /// first reaches zero, the gas has been brought fully to rest. Suppressing the rebound there
    /// is the perfect-momentum-sink idealization, so by conservation `J_wall = p_in` and
    /// `e_eff → 0`. This is the degenerate, fully-absorbing limit of the real wall (ADR-0005): a
    /// lossless gas cannot physically stick (it would re-expand to the bounce ceiling), so rung A
    /// *imposes* the stop rather than modeling a loss channel; later rungs replace it with a
    /// realistic stick/condensation model.
    ///
    /// While gas is compressed against the wall the wall force is `≥ 0`, so the total momentum
    /// rises monotonically from `−p_in` and crosses zero exactly once. That crossing is bracketed
    /// and the final step linearly interpolated (`dp/dt = F_wall`), so the reported impulse is not
    /// biased by a full-step overshoot past stagnation. The residual `O(Δx)` error is the same
    /// impulse-vs-momentum consistency as [`Self::run_bounce`].
    pub fn run_stick_bounce(&mut self) -> BounceResult {
        let p_initial = self.total_momentum();
        let incident = p_initial.abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut p_old = p_initial;
        let mut force_old = self.wall_force();
        let max_steps = 400 * self.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            let dt = self.stable_dt();
            self.step(dt);
            let force_new = self.wall_force();
            let p_new = self.total_momentum();

            if p_new >= 0.0 {
                // Stagnation lies within this step. Momentum rises linearly across the step
                // (`dp/dt = F_wall`), so it crosses zero at fraction θ = −p_old/(p_new − p_old);
                // integrate the wall force trapezoidally only over [0, θ·dt].
                let theta = (-p_old / (p_new - p_old)).clamp(0.0, 1.0);
                let force_cross = force_old + theta * (force_new - force_old);
                wall_impulse += 0.5 * theta * dt * (force_old + force_cross);
                return BounceResult {
                    wall_impulse,
                    incident_momentum: incident,
                    residual_momentum: 0.0, // gas held at rest; rebound suppressed
                    e_eff: wall_impulse / incident - 1.0,
                    peak_wall_force: peak,
                };
            }

            wall_impulse += 0.5 * dt * (force_old + force_new);
            p_old = p_new;
            force_old = force_new;
        }

        // Never stagnated within the step cap (a bug, or a Mach so low the slug barely slows):
        // report what we have so the test surfaces it.
        BounceResult {
            wall_impulse,
            incident_momentum: incident,
            residual_momentum: self.total_momentum(),
            e_eff: wall_impulse / incident - 1.0,
            peak_wall_force: peak,
        }
    }
}

/// Per-cell inputs for one radiation substep, built from the gas state of a [`Tube<TableEos>`] (B5).
/// Owns the `Vec`s so a borrowed [`Medium`] can be constructed from it for a single `fld_substep`.
#[derive(Debug)]
struct RadFields {
    dx: Vec<f64>,
    center_spacing: Vec<f64>,
    temp: Vec<f64>,
    cv_vol: Vec<f64>,
    chi_planck: Vec<f64>,
    chi_ross: Vec<f64>,
}

impl RadFields {
    /// Borrow as a [`Medium`] for one radiation substep (no external volumetric source).
    fn medium(&self) -> Medium<'_> {
        Medium {
            dx: &self.dx,
            center_spacing: &self.center_spacing,
            temp: &self.temp,
            cv_vol: &self.cv_vol,
            chi_planck: &self.chi_planck,
            chi_ross: &self.chi_ross,
            source: None,
        }
    }
}

/// Owned per-cell gas-side inputs for one coupled conduction substep (B-flux), built from the gas
/// state of a `Tube<TableEos>` so a borrowed [`GasConductionState`] can be handed to
/// [`Solid::step_coupled`]. Carries the specific `c_v` as well, for the post-solve energy update
/// `Δe = c_v·ΔT` that mirrors the operator's per-cell heat capacity `C = ρ c_v dx`.
#[derive(Debug)]
struct GasConductionFields {
    dx: Vec<f64>,
    temp: Vec<f64>,
    cv_vol: Vec<f64>,
    k_gas: Vec<f64>,
    cv: Vec<f64>,
}

impl GasConductionFields {
    /// Borrow as a [`GasConductionState`] for one [`Solid::step_coupled`] call.
    fn state(&self) -> GasConductionState<'_> {
        GasConductionState {
            dx: &self.dx,
            temp: &self.temp,
            cv_vol: &self.cv_vol,
            k_gas: &self.k_gas,
        }
    }
}

/// One gas-side conduction substep (B-flux, ADR-0005), shared by [`CoupledBounce`] and
/// [`CondensingBounce`]: cool the gas into the wall `solid` via the combined gas+solid operator
/// ([`Solid::step_coupled`]), updating each gas cell's energy by `Δe = c_v·ΔT`, and return the
/// channel-2 interface loss `q·dt`. A no-op returning `0` when there is no wall, or when the table
/// carries no `k_gas` (no gas transport data — the high-v plasma table; its transport is deferred).
fn conduction_into_wall(tube: &mut Tube<TableEos>, wall: &mut Option<Solid>, dt: f64) -> f64 {
    let Some(fields) = tube.gas_conduction_fields() else {
        return 0.0;
    };
    let Some(solid) = wall.as_mut() else {
        return 0.0;
    };
    let step = solid.step_coupled(&fields.state(), dt);
    for (j, &dtemp) in step.gas_dtemp.iter().enumerate() {
        tube.energy[j] = fields.cv[j].mul_add(dtemp, tube.energy[j]).max(0.0);
    }
    step.interface_flux * dt
}

impl Tube<TableEos> {
    /// A cold gas slug fired at a rigid wall (`x = 0`) with a trailing free (vacuum) surface
    /// (`x = length`), in **SI** units — the production analogue of [`Tube::slug`] (ADR-0001) used
    /// by the `e_eff(ρ)` sweep (`crates/sweep`).
    ///
    /// `cells` cells uniformly fill `x ∈ [0, length]` at density `rho_impact`, all coasting toward
    /// the wall at `u = −v`. The cold cloud's pressure (and hence specific internal energy) is seeded
    /// from the table at `(rho_impact, t0)`: `p₀ = p(ρ, T₀)`, and [`Tube::with_eos`] then recovers the
    /// consistent `e₀ = e(ρ, T₀)` by inverting that pressure. The incident Mach number `v / c_s(ρ, T₀)`
    /// is therefore set entirely by how cold `T₀` is, exactly as in the normalized [`Tube::slug`].
    #[must_use]
    pub fn slug_si(
        cells: usize,
        rho_impact: f64,
        v: f64,
        length: f64,
        t0: f64,
        eos: TableEos,
        viscosity: Viscosity,
    ) -> Self {
        let x: Vec<f64> = (0..=cells)
            .map(|i| i as f64 / cells as f64 * length)
            .collect();
        let p0 = eos.table().pressure(rho_impact, t0);
        let rho = vec![rho_impact; cells];
        let pressure = vec![p0; cells];
        let vel = vec![-v; cells];
        let mut tube = Self::with_eos(x, &rho, &vel, &pressure, eos, viscosity);
        tube.right = Boundary::Free;
        tube.enforce_wall_velocities();
        tube
    }

    /// Build the radiation medium from the current gas state: per cell the temperature `T(ρ, e)`,
    /// the volumetric heat capacity `ρ c_v`, and the per-length opacities `χ = κ ρ` — Planck for
    /// emission/absorption, Rosseland for the flux-limited diffusion (ADR-0006). The mesh geometry
    /// (`dx`, the `N−1` center-to-center spacings) is read from the current Lagrangian node
    /// positions, so it tracks the moving mesh each step.
    fn radiation_fields(&self) -> RadFields {
        let n = self.cells();
        let table = self.eos.table();
        let dx: Vec<f64> = (0..n).map(|j| self.width(j)).collect();
        let center_spacing: Vec<f64> = (0..n - 1)
            .map(|i| self.center(i + 1) - self.center(i))
            .collect();
        let mut temp = vec![0.0; n];
        let mut cv_vol = vec![0.0; n];
        let mut chi_planck = vec![0.0; n];
        let mut chi_ross = vec![0.0; n];
        for j in 0..n {
            let rho = self.density(j);
            let t = self.eos.temperature(rho, self.energy[j]);
            temp[j] = t;
            cv_vol[j] = rho * table.cv(rho, t);
            chi_planck[j] = rho * table.kappa_planck(rho, t);
            chi_ross[j] = rho * table.kappa_rosseland(rho, t);
        }
        RadFields {
            dx,
            center_spacing,
            temp,
            cv_vol,
            chi_planck,
            chi_ross,
        }
    }

    /// Build the gas-side conduction inputs from the current gas state (B-flux): per cell the width,
    /// temperature `T(ρ, e)`, volumetric heat capacity `ρ c_v`, thermal conductivity `k_gas(ρ, e)`,
    /// and specific `c_v`. Returns `None` when the table carries no `k_gas` (no gas transport data —
    /// e.g. the high-v plasma table, whose transport is the deferred B-flux sibling), in which case
    /// the gas gets no conduction. `k_gas` is a table-wide property, so a single cell-0 probe decides.
    fn gas_conduction_fields(&self) -> Option<GasConductionFields> {
        self.eos.k_gas(self.density(0), self.energy[0])?;
        let n = self.cells();
        let table = self.eos.table();
        let mut dx = vec![0.0; n];
        let mut temp = vec![0.0; n];
        let mut cv_vol = vec![0.0; n];
        let mut k_gas = vec![0.0; n];
        let mut cv = vec![0.0; n];
        for j in 0..n {
            let rho = self.density(j);
            let t = self.eos.temperature(rho, self.energy[j]);
            let cvj = table.cv(rho, t);
            dx[j] = self.width(j);
            temp[j] = t;
            cv[j] = cvj;
            cv_vol[j] = rho * cvj;
            k_gas[j] = table
                .k_gas(rho, t)
                .expect("k_gas present (probed at cell 0)");
        }
        Some(GasConductionFields {
            dx,
            temp,
            cv_vol,
            k_gas,
            cv,
        })
    }
}

/// Outcome of a coupled bounce: the [`BounceResult`] plus the three energy loss channels
/// (per unit wall area) the rigid wall splits the deficit into (ADR-0016).
#[derive(Debug, Clone, Copy)]
pub struct CoupledBounceResult {
    /// The usual restitution/impulse bookkeeping.
    pub bounce: BounceResult,
    /// Channel 1a — radiation absorbed at the wall (`x = 0`), `∫ (c/2)(E₀ − e_inc) dt`.
    pub loss_radiative_wall: f64,
    /// Channel 1b — radiation escaping to space at the far (re-expansion) end.
    pub loss_escape_space: f64,
    /// Channel 2 — heat conducted into the wall solid, `∫ q_wall dt`.
    pub loss_conductive: f64,
}

/// A radiation + conduction coupled slug bounce (B5b). Holds the gas [`Tube<TableEos>`], the
/// per-cell radiation energy density `e_rad`, an optional wall conducting [`Solid`], the physical
/// constants, and the accumulated loss channels. One step is **Lie-split** at the hydro `dt`: the
/// hydro update, then one implicit gray-FLD substep (radiation transport + matter exchange), then
/// the conductive wall loss — each operator reading the state the previous one left (ADR-0006,
/// ADR-0005). Radiation work/pressure on the gas is deferred (ADR-0006); this moves radiation
/// *energy* and the loss it carries off through the wall and to space.
#[derive(Debug)]
pub struct CoupledBounce {
    tube: Tube<TableEos>,
    e_rad: Vec<f64>,
    wall: Option<Solid>,
    consts: RadConstants,
    limiter: Limiter,
    /// Radiation BC at the wall (`x = 0`); the cold black absorber is `Marshak(0)` (ADR-0005).
    bc_wall: RadBc,
    /// Radiation BC at the far (re-expansion) end; escape to space is `Marshak(0)`.
    bc_space: RadBc,
    loss_radiative_wall: f64,
    loss_escape_space: f64,
    loss_conductive: f64,
}

impl CoupledBounce {
    /// Wrap a tabulated-EOS `tube` (and an optional wall `solid`) for a coupled bounce. Radiation
    /// starts in local equilibrium `e_rad = a T⁴`; the wall and far end default to cold black
    /// absorbers (`Marshak(0)`) — loss channels 1a and 1b. Pass `wall = None` to disable conduction.
    #[must_use]
    pub fn new(
        tube: Tube<TableEos>,
        wall: Option<Solid>,
        consts: RadConstants,
        limiter: Limiter,
    ) -> Self {
        let e_rad: Vec<f64> = (0..tube.cells())
            .map(|j| {
                let rho = tube.density(j);
                let t = tube.eos.temperature(rho, tube.energy[j]);
                consts.a * t.powi(4)
            })
            .collect();
        Self {
            tube,
            e_rad,
            wall,
            consts,
            limiter,
            bc_wall: RadBc::Marshak(0.0),
            bc_space: RadBc::Marshak(0.0),
            loss_radiative_wall: 0.0,
            loss_escape_space: 0.0,
            loss_conductive: 0.0,
        }
    }

    /// Total radiation field energy `Σ e_rad_j dx_j` (per unit wall area).
    #[must_use]
    pub fn radiation_energy(&self) -> f64 {
        (0..self.tube.cells())
            .map(|j| self.e_rad[j] * self.tube.width(j))
            .sum()
    }

    /// Total matter internal energy `Σ m_j e_j` (per unit wall area).
    #[must_use]
    pub fn matter_internal_energy(&self) -> f64 {
        (0..self.tube.cells())
            .map(|j| self.tube.mass[j] * self.tube.energy[j])
            .sum()
    }

    /// One implicit gray-FLD substep: transport radiation, exchange energy with the matter, and tally
    /// the radiation leaving through any absorbing boundary (B4b accounting, post-solve `e_rad`).
    fn radiation_substep(&mut self, dt: f64) {
        let fields = self.tube.radiation_fields();
        let delta_e = fld_substep(
            &fields.medium(),
            &mut self.e_rad,
            self.bc_wall,
            self.bc_space,
            dt,
            self.consts,
            self.limiter,
        );
        // Deposit the matter's share of the exchange: Δe is energy/volume, so Δe/ρ is specific.
        for (j, &de) in delta_e.iter().enumerate() {
            let rho = self.tube.density(j);
            self.tube.energy[j] = (self.tube.energy[j] + de / rho).max(0.0);
        }
        // Net outflow through an absorbing end: (c/2)(E_edge − e_inc) per area per time.
        let c = self.consts.c;
        if let RadBc::Marshak(e_inc) = self.bc_wall {
            self.loss_radiative_wall += dt * 0.5 * c * (self.e_rad[0] - e_inc);
        }
        if let RadBc::Marshak(e_inc) = self.bc_space {
            let last = self.tube.cells() - 1;
            self.loss_escape_space += dt * 0.5 * c * (self.e_rad[last] - e_inc);
        }
    }

    /// Conductive wall loss (channel 2), via the **gas-side conduction operator** (B-flux, ADR-0005).
    /// Rather than pinning the interface to the near-wall gas temperature and draining the single wall
    /// cell — the inviscid kernel's missing gas-side resistance, which over-drained the thin cell and
    /// collapsed the bounce at high Mach — this solves the gas *and* the wall solid as one conduction
    /// system ([`Solid::step_coupled`]), so the interface temperature emerges from flux continuity.
    /// The per-cell gas energy is updated by `Δe = c_v·ΔT`, and the interface flux is tallied as the
    /// channel-2 loss. No-ops when there is no wall, or when the table carries no `k_gas` (the high-v
    /// plasma table — its transport is the deferred B-flux sibling, so high-v conduction stays off).
    fn conduction_substep(&mut self, dt: f64) {
        self.loss_conductive += conduction_into_wall(&mut self.tube, &mut self.wall, dt);
    }

    /// One coupled step: hydro, then radiation, then conduction (Lie split at the hydro `dt`).
    fn coupled_step(&mut self, dt: f64) {
        self.tube.step(dt);
        self.radiation_substep(dt);
        self.conduction_substep(dt);
    }

    /// Fire the coupled slug at the wall, integrating to the same `10⁻³`-of-peak tail guard as
    /// [`Tube::run_bounce`], and return the restitution plus the loss-channel decomposition.
    pub fn run(&mut self) -> CoupledBounceResult {
        let incident = self.tube.total_momentum().abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = self.tube.wall_force();
        let max_steps = 400 * self.tube.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            if force_old < 0.5 * peak {
                past_peak = true;
            }
            if past_peak && force_old < 1e-3 * peak {
                break;
            }
            let dt = self.tube.stable_dt();
            self.coupled_step(dt);
            let force_new = self.tube.wall_force();
            wall_impulse += 0.5 * dt * (force_old + force_new);
            force_old = force_new;
        }

        CoupledBounceResult {
            bounce: BounceResult {
                wall_impulse,
                incident_momentum: incident,
                residual_momentum: self.tube.total_momentum(),
                e_eff: wall_impulse / incident - 1.0,
                peak_wall_force: peak,
            },
            loss_radiative_wall: self.loss_radiative_wall,
            loss_escape_space: self.loss_escape_space,
            loss_conductive: self.loss_conductive,
        }
    }
}

/// Outcome of a condensing bounce (Rung C / B-flux): the restitution plus the loss channels.
#[derive(Debug, Clone, Copy)]
pub struct CondensingBounceResult {
    /// The usual restitution/impulse bookkeeping.
    pub bounce: BounceResult,
    /// Channel 3 — energy carried off by condensate that stuck to the wall (per unit wall area;
    /// includes the latent heat, which is already inside `e`).
    pub loss_condensation: f64,
    /// Channel 2 — heat conducted into the wall solid (per unit wall area). `0` without a wall (the
    /// adiabatic upper bound); positive once a conducting wall cools the near-wall gas (B-flux).
    pub loss_conductive: f64,
    /// Condensate mass deposited at the wall (per unit area) — the gas's lost mass, for the
    /// mass-sink closure check.
    pub stuck_mass: f64,
}

/// A low-v condensing slug bounce (Rung C / B-flux): each step is the hydro update, then an optional
/// gas-side **conduction** substep into a cold wall (channel 2), then the irreversible
/// **wall-deposition mass sink** (ADR-0004 channel 3). The other condensation channel, **bulk
/// vapor-pressure collapse**, needs no code here: it lives in the two-phase EOS (`p → p_sat(T)`,
/// latent heat folded into `e`), so it acts through the ordinary pressure the hydro already sees.
///
/// Conduction and deposition are **sequenced, not independent** (ADR-0004 amendment): wall deposition
/// only fires once the cold plate cools the near-wall gas below `T_sat`, which raises `liquid_frac` in
/// the wall cell. With no wall the run is the *adiabatic upper bound* (the original Rung C); adding a
/// wall (B-flux) activates the dormant deposition that makes 3.2 km/s the condensation-dominated worst
/// case the design anticipates.
///
/// The sink removes the condensate that *newly forms* in the wall cell each step (so the total stuck
/// mass tracks the cumulative condensation, not the step count — dt-convergent), carrying off its
/// mass, momentum, and energy. Removing rebounding mass is what lowers `e_eff` below the lossless
/// ceiling; the dead-stick floor `e_eff → 0` is recovered when the compression impulse alone
/// (`= p_in`) is delivered and no gas rebounds.
#[derive(Debug)]
pub struct CondensingBounce {
    tube: Tube<TableEos>,
    /// Sticking coefficient `α ∈ [0, 1]` (baseline 1 — the pessimistic equilibrium bound, ADR-0004).
    alpha: f64,
    /// Optional cold conducting wall (B-flux). `None` is the adiabatic upper bound (radiation is
    /// negligible at 3.2 km/s, design §3, so it stays off either way).
    wall: Option<Solid>,
    loss_condensation: f64,
    loss_conductive: f64,
    stuck_mass: f64,
    /// Liquid mass in the wall cell at the end of the previous step — the baseline against which
    /// newly-condensed (hence newly-sticking) mass is measured.
    m_liq_prev: f64,
}

impl CondensingBounce {
    /// Wrap a tabulated-EOS `tube` for an *adiabatic* condensing bounce (no wall) with wall sticking
    /// coefficient `alpha` — the original Rung C upper bound.
    #[must_use]
    pub fn new(tube: Tube<TableEos>, alpha: f64) -> Self {
        Self::new_with_wall(tube, alpha, None)
    }

    /// Wrap a tabulated-EOS `tube` for a condensing bounce with sticking coefficient `alpha` and an
    /// optional cold conducting `wall` (B-flux). The wall cools the near-wall gas below `T_sat`, which
    /// is what *drives* the wall-deposition sink (the deposition channel is conduction-gated,
    /// ADR-0004). Conduction engages only if the table also carries `k_gas`.
    #[must_use]
    pub fn new_with_wall(tube: Tube<TableEos>, alpha: f64, wall: Option<Solid>) -> Self {
        let m_liq_prev = Self::wall_liquid_mass(&tube);
        Self {
            tube,
            alpha,
            wall,
            loss_condensation: 0.0,
            loss_conductive: 0.0,
            stuck_mass: 0.0,
            m_liq_prev,
        }
    }

    /// Gas-side conduction into the cold wall (channel 2, B-flux): the same operator
    /// [`CoupledBounce`] uses. This is what cools the near-wall gas below `T_sat` and thereby
    /// *activates* the deposition sink. No-op without a wall or `k_gas`.
    fn conduction_substep(&mut self, dt: f64) {
        self.loss_conductive += conduction_into_wall(&mut self.tube, &mut self.wall, dt);
    }

    /// Condensed (liquid) mass currently in the wall cell, `liquid_frac(ρ₀, e₀) · m₀`.
    fn wall_liquid_mass(tube: &Tube<TableEos>) -> f64 {
        let e0 = tube.energy[0];
        if e0 <= 0.0 {
            return 0.0;
        }
        tube.eos.liquid_fraction(tube.density(0), e0) * tube.mass[0]
    }

    /// Wall-deposition sink (channel 3): the liquid that has newly condensed in the wall cell since
    /// the previous step sticks irreversibly with coefficient `α`, removing its mass (a Lagrangian
    /// mass sink — its momentum leaves with it) and energy from the gas.
    fn condensation_substep(&mut self) {
        let e0 = self.tube.energy[0];
        if e0 <= 0.0 {
            return;
        }
        let m_liq = Self::wall_liquid_mass(&self.tube);
        let newly = m_liq - self.m_liq_prev; // liquid that condensed since the last step
        if newly > 0.0 && self.alpha > 0.0 {
            let stuck = (self.alpha * newly).min(self.tube.mass[0]);
            self.loss_condensation += stuck * e0; // energy leaving with the condensate (incl. latent)
            self.stuck_mass += stuck;
            self.tube.mass[0] -= stuck; // mass sink; the condensate's momentum leaves with the mass
        }
        // Re-baseline to the liquid still present (intensive `liquid_frac` is unchanged this step).
        self.m_liq_prev = Self::wall_liquid_mass(&self.tube);
    }

    /// Fire the condensing slug at the wall, integrating to the same `10⁻³`-of-peak tail guard as
    /// [`Tube::run_bounce`], and return the restitution plus the condensation loss.
    pub fn run(&mut self) -> CondensingBounceResult {
        let incident = self.tube.total_momentum().abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = self.tube.wall_force();
        let max_steps = 400 * self.tube.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            if force_old < 0.5 * peak {
                past_peak = true;
            }
            if past_peak && force_old < 1e-3 * peak {
                break;
            }
            let dt = self.tube.stable_dt();
            self.tube.step(dt);
            // Conduction first: the cold wall cools the near-wall gas below T_sat, raising
            // liquid_frac, which the deposition sink then sticks (the two channels are sequenced).
            self.conduction_substep(dt);
            self.condensation_substep();
            let force_new = self.tube.wall_force();
            wall_impulse += 0.5 * dt * (force_old + force_new);
            force_old = force_new;
        }

        CondensingBounceResult {
            bounce: BounceResult {
                wall_impulse,
                incident_momentum: incident,
                residual_momentum: self.tube.total_momentum(),
                e_eff: wall_impulse / incident - 1.0,
                peak_wall_force: peak,
            },
            loss_condensation: self.loss_condensation,
            loss_conductive: self.loss_conductive,
            stuck_mass: self.stuck_mass,
        }
    }
}

/// Parameters of the quasi-steady **ablating wall** (Rung E, ADR-0014). The incoming wall flux boils
/// off ablator at the surface energy balance `ṁ = q_in / Q*`; the vapor is injected as a cold mass
/// source at the wall. `Q*` is the effective heat of ablation (silicone ~2–10 MJ/kg literature);
/// `t_vapor` is the (cold) temperature the injected vapor enters at — conservatively the cloud's own
/// `T₀`, so the vapor adds little enthalpy and the recovery is not flattered.
#[derive(Debug, Clone, Copy)]
pub struct Ablation {
    /// Effective heat of ablation `Q*` [J/kg]. Larger `Q*` ⇒ less ablation; `Q* → ∞` ⇒ rigid wall.
    q_star: f64,
    /// Injected-vapor temperature [K] — the cold vapor's specific energy is read from the table here.
    t_vapor: f64,
}

impl Ablation {
    /// An ablation model with heat of ablation `q_star` [J/kg] and cold-vapor injection temperature
    /// `t_vapor` [K] (E2 blowing / E3 vapor-shielding parameters are added by later slices).
    ///
    /// # Panics
    /// Panics unless `q_star > 0` and `t_vapor > 0`.
    #[must_use]
    pub fn new(q_star: f64, t_vapor: f64) -> Self {
        assert!(
            q_star > 0.0 && t_vapor > 0.0,
            "Q* and t_vapor must be positive"
        );
        Self { q_star, t_vapor }
    }
}

/// Outcome of an ablating bounce (Rung E): the restitution/loss decomposition of [`CoupledBounce`]
/// plus the ablation bookkeeping (per unit wall area).
#[derive(Debug, Clone, Copy)]
pub struct AblatingBounceResult {
    /// The usual restitution/impulse bookkeeping.
    pub bounce: BounceResult,
    /// Channel 1a — radiation absorbed at the wall.
    pub loss_radiative_wall: f64,
    /// Channel 1b — radiation escaping to space at the re-expansion end.
    pub loss_escape_space: f64,
    /// Channel 2 — heat conducted into the wall solid (`0` when `wall = None`).
    pub loss_conductive: f64,
    /// Total vapor mass boiled off and injected at the wall, per unit area [kg/m²].
    pub ablated_mass: f64,
    /// Energy consumed by ablation `Σ Q*·ṁ·dt` [J/m²] — by the quasi-steady balance this equals the
    /// incoming wall flux it converts (`loss_radiative_wall + loss_conductive`).
    pub loss_ablation: f64,
}

/// A radiation + conduction slug bounce against a **quasi-steady ablating wall** (Rung E, ADR-0014):
/// [`CoupledBounce`] plus an ablation **mass source** — the inverse of [`CondensingBounce`]'s wall
/// sink. Each step is Lie-split `hydro → radiation → conduction → **ablation**`: the incoming wall
/// flux `q_in = q_rad_wall + q_cond` boils off `ṁ = q_in/Q*` of ablator (the surface energy balance),
/// injected as cold vapor into the wall cell. The rigid [`CoupledBounce`] is the conservative floor;
/// the ablating wall is the best-estimate refinement (ADR-0013), and `Q* → ∞` recovers the floor
/// exactly. Runs with `wall = None` as its realistic high-v config (the high-v table carries no
/// `k_gas`, so conduction — hence blowing — is off; the recovery is shielding + mass injection).
#[derive(Debug)]
pub struct AblatingBounce {
    tube: Tube<TableEos>,
    e_rad: Vec<f64>,
    wall: Option<Solid>,
    consts: RadConstants,
    limiter: Limiter,
    bc_wall: RadBc,
    bc_space: RadBc,
    ablation: Ablation,
    loss_radiative_wall: f64,
    loss_escape_space: f64,
    loss_conductive: f64,
    ablated_mass: f64,
    loss_ablation: f64,
}

impl AblatingBounce {
    /// Wrap a tabulated-EOS `tube` (and optional wall `solid`) for an ablating bounce with the given
    /// `ablation` model. Radiation starts in local equilibrium `e_rad = a T⁴`; the wall and far end
    /// are cold black absorbers (`Marshak(0)`), as in [`CoupledBounce::new`].
    #[must_use]
    pub fn new(
        tube: Tube<TableEos>,
        wall: Option<Solid>,
        consts: RadConstants,
        limiter: Limiter,
        ablation: Ablation,
    ) -> Self {
        let e_rad: Vec<f64> = (0..tube.cells())
            .map(|j| {
                let rho = tube.density(j);
                let t = tube.eos.temperature(rho, tube.energy[j]);
                consts.a * t.powi(4)
            })
            .collect();
        Self {
            tube,
            e_rad,
            wall,
            consts,
            limiter,
            bc_wall: RadBc::Marshak(0.0),
            bc_space: RadBc::Marshak(0.0),
            ablation,
            loss_radiative_wall: 0.0,
            loss_escape_space: 0.0,
            loss_conductive: 0.0,
            ablated_mass: 0.0,
            loss_ablation: 0.0,
        }
    }

    /// One implicit gray-FLD substep (as [`CoupledBounce::radiation_substep`]), returning the
    /// **radiative wall flux absorbed this step** `dt·(c/2)(E₀ − e_inc)` — the radiative part of the
    /// `q_in` that drives ablation.
    fn radiation_substep(&mut self, dt: f64) -> f64 {
        let fields = self.tube.radiation_fields();
        let delta_e = fld_substep(
            &fields.medium(),
            &mut self.e_rad,
            self.bc_wall,
            self.bc_space,
            dt,
            self.consts,
            self.limiter,
        );
        for (j, &de) in delta_e.iter().enumerate() {
            let rho = self.tube.density(j);
            self.tube.energy[j] = (self.tube.energy[j] + de / rho).max(0.0);
        }
        let c = self.consts.c;
        let mut wall_flux = 0.0;
        if let RadBc::Marshak(e_inc) = self.bc_wall {
            wall_flux = dt * 0.5 * c * (self.e_rad[0] - e_inc);
            self.loss_radiative_wall += wall_flux;
        }
        if let RadBc::Marshak(e_inc) = self.bc_space {
            let last = self.tube.cells() - 1;
            self.loss_escape_space += dt * 0.5 * c * (self.e_rad[last] - e_inc);
        }
        wall_flux
    }

    /// Conductive wall loss (channel 2), via the shared gas-side operator. Returns the flux conducted
    /// into the wall this step — the conductive part of `q_in`. `0` when `wall = None` or no `k_gas`.
    fn conduction_substep(&mut self, dt: f64) -> f64 {
        let q = conduction_into_wall(&mut self.tube, &mut self.wall, dt);
        self.loss_conductive += q;
        q
    }

    /// Ablation substep: the quasi-steady surface energy balance. The incoming wall flux `q_in`
    /// (energy/area this step) boils off `ṁ·dt = q_in/Q*` of ablator, injected as **cold vapor** into
    /// the wall cell — the inverse of [`CondensingBounce`]'s wall sink. The vapor enters at the wall
    /// node (`u = 0`), adding no axial momentum directly; it acts on `e_eff` through the near-wall
    /// pressure it raises (E1) and, in later slices, by shielding radiation (E3). Mass is added to
    /// cell 0 and its specific energy mass-blended toward the cold vapor's.
    fn ablation_substep(&mut self, q_in: f64) {
        if q_in <= 0.0 {
            return;
        }
        let dm = q_in / self.ablation.q_star;
        if dm <= 0.0 {
            return;
        }
        let rho_wall = self.tube.density(0);
        let e_vapor = self
            .tube
            .eos
            .table()
            .energy(rho_wall, self.ablation.t_vapor);
        let m_old = self.tube.mass[0];
        let e_old = self.tube.energy[0];
        self.tube.mass[0] = m_old + dm;
        self.tube.energy[0] = dm.mul_add(e_vapor, m_old * e_old) / (m_old + dm);
        self.ablated_mass += dm;
        self.loss_ablation += q_in; // = dm · Q*
    }

    /// One ablating step: hydro, radiation, conduction, then ablation (Lie split at the hydro `dt`).
    fn ablating_step(&mut self, dt: f64) {
        self.tube.step(dt);
        let q_rad = self.radiation_substep(dt);
        let q_cond = self.conduction_substep(dt);
        self.ablation_substep((q_rad + q_cond).max(0.0));
    }

    /// Fire the ablating slug at the wall, integrating to the same `10⁻³`-of-peak tail guard as
    /// [`CoupledBounce::run`], and return the restitution, loss decomposition, and ablation tally.
    pub fn run(&mut self) -> AblatingBounceResult {
        let incident = self.tube.total_momentum().abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = self.tube.wall_force();
        let max_steps = 400 * self.tube.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            if force_old < 0.5 * peak {
                past_peak = true;
            }
            if past_peak && force_old < 1e-3 * peak {
                break;
            }
            let dt = self.tube.stable_dt();
            self.ablating_step(dt);
            let force_new = self.tube.wall_force();
            wall_impulse += 0.5 * dt * (force_old + force_new);
            force_old = force_new;
        }

        AblatingBounceResult {
            bounce: BounceResult {
                wall_impulse,
                incident_momentum: incident,
                residual_momentum: self.tube.total_momentum(),
                e_eff: wall_impulse / incident - 1.0,
                peak_wall_force: peak,
            },
            loss_radiative_wall: self.loss_radiative_wall,
            loss_escape_space: self.loss_escape_space,
            loss_conductive: self.loss_conductive,
            ablated_mass: self.ablated_mass,
            loss_ablation: self.loss_ablation,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{CondensingBounce, CoupledBounce, Tube, Viscosity};
    use crate::conduction::Solid;
    use crate::eos::TableEos;
    use crate::radiation::{Limiter, RadBc, RadConstants};
    use approx::assert_relative_eq;
    use tables::Table;

    const GAMMA: f64 = 1.4;

    /// An ideal-gas EOS table (`e = T`, so `c_v = 1`) with opacity power laws scaled by the
    /// coefficients `(kr, kp)`: `κ_R = kr·ρ²·T^-3.5`, `κ_P = kp·ρ·T^-2`. Tiny coefficients make
    /// the gas effectively **transparent** (`κ → 0`: no emission/absorption, the radiation-off
    /// regression); `O(1)` coefficients give a strongly coupled, lossy gas. Power laws in `(ρ, T)`,
    /// so the table's log-log interpolation is exact (the `χ = κ ρ` builder is checkable).
    fn gas_table(kr: f64, kp: f64) -> TableEos {
        let n = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64)) // 0.01 … 10
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64)) // 0.05 … 200
            .collect();
        let (mut p, mut e, mut cs, mut kr_v, mut kp_v, mut kg_v) = (
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(t); // e = T ⇒ c_v = 1
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
                kr_v.push(kr * r.powf(2.0) * t.powf(-3.5)); // κ_R(ρ,T)
                kp_v.push(kp * r * t.powf(-2.0)); // κ_P(ρ,T)
                // A synthetic gas conductivity for the B-flux conduction operator: small enough that
                // the gas effusivity √(k_gas·ρ·c_v) sits well below the test solids' effusivity, so
                // the gas-side resistance is real and the cold plate cannot over-drain the wall cell.
                kg_v.push(0.05);
            }
        }
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": {
                "p": p, "e": e, "c_s": cs,
                "kappa_rosseland": kr_v, "kappa_planck": kp_v, "k_gas": kg_v,
            },
        });
        TableEos::new(Table::from_json(&json.to_string()).unwrap())
    }

    /// A cold gas slug fired at a wall (like [`Tube::slug`]) but carrying a tabulated EOS, for the
    /// coupled-bounce tests. A thin wrapper over the production [`Tube::slug_si`] in normalized units
    /// (`ρ₀ = 1`, `v = 1`, `L = 1`): for the `e = T` synthetic tables a cold `T₀ = 1/(γ(γ−1)M²)` gives
    /// `p₀ = (γ−1)ρ₀T₀ = 1/(γM²)`, i.e. incident Mach number `M = v/c₀`.
    fn slug_with_table(cells: usize, mach: f64, eos: TableEos) -> Tube<TableEos> {
        let t0 = 1.0 / (GAMMA * (GAMMA - 1.0) * mach * mach);
        Tube::slug_si(
            cells,
            1.0,
            1.0,
            1.0,
            t0,
            eos,
            Viscosity::VON_NEUMANN_RICHTMYER,
        )
    }

    /// An ideal-gas EOS table (`e = T`, like [`gas_table`]) carrying a synthetic `liquid_frac` field
    /// that rises gently with compression (`liquid_frac = clamp((ρ − 2)/20, 0, 1)`), so a slug
    /// compressing at the wall condenses and the wall-sticking sink (Rung C) engages. `liquid_frac`
    /// is an independent field exercising the *sink mechanism*; the real two-phase physics lives in
    /// the low-v table (C1).
    fn condensing_table() -> TableEos {
        let n: usize = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64)) // 0.01 … 10
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64)) // 0.05 … 200
            .collect();
        let (mut p, mut e, mut cs, mut kr, mut kp, mut lf) = (
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(t);
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
                kr.push(1e-10); // transparent: radiation off at low-v
                kp.push(1e-10);
                lf.push(((r - 2.0) / 20.0).clamp(0.0, 1.0)); // condenses as the gas compresses
            }
        }
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": { "p": p, "e": e, "c_s": cs,
                        "kappa_rosseland": kr, "kappa_planck": kp, "liquid_frac": lf },
        });
        TableEos::new(Table::from_json(&json.to_string()).unwrap())
    }

    /// A synthetic two-phase table (`e = T`, `c_v = 1`) whose `liquid_frac` requires **both**
    /// compression *and* cooling: `clamp((ρ−2)/8)·clamp((T_thr−T)/T_thr)` with `T_thr = 3`. So hot
    /// dense gas is vapor (`liquid_frac ≈ 0`) and only condenses once cooled below the threshold — the
    /// **conduction-gated** deposition of the ADR-0004 amendment. Carries `k_gas` so the B-flux
    /// conduction operator engages.
    fn cooling_condensing_table() -> TableEos {
        let n: usize = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64))
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64))
            .collect();
        let t_thr = 3.0;
        let mut p = Vec::new();
        let mut e = Vec::new();
        let mut cs = Vec::new();
        let mut kr = Vec::new();
        let mut kp = Vec::new();
        let mut lf = Vec::new();
        let mut kg = Vec::new();
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(t);
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
                kr.push(1e-10);
                kp.push(1e-10);
                let comp = ((r - 2.0) / 8.0).clamp(0.0, 1.0); // condenses only once compressed
                let cool = ((t_thr - t) / t_thr).clamp(0.0, 1.0); // …and only once cooled
                lf.push(comp * cool);
                kg.push(0.05);
            }
        }
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": { "p": p, "e": e, "c_s": cs, "kappa_rosseland": kr,
                        "kappa_planck": kp, "liquid_frac": lf, "k_gas": kg },
        });
        TableEos::new(Table::from_json(&json.to_string()).unwrap())
    }

    /// With `α = 0` the wall-sticking sink is inert, so the condensing bounce reproduces the plain
    /// [`Tube::run_bounce`] exactly and no mass is removed.
    #[allow(clippy::float_cmp)] // α=0 removes nothing → exactly zero stuck mass
    #[test]
    fn condensing_alpha_zero_matches_plain_bounce() {
        let tube = slug_with_table(150, 4.0, condensing_table());
        let plain = tube.clone().run_bounce().e_eff;
        let cond = CondensingBounce::new(tube, 0.0).run();
        assert_relative_eq!(cond.bounce.e_eff, plain, max_relative = 1e-12);
        assert_eq!(cond.stuck_mass, 0.0);
    }

    /// `α = 1`: the wall sink removes condensate (mass-sink closure: the gas's lost mass equals the
    /// tallied stuck mass), which lowers `e_eff` below the lossless ceiling while staying physical.
    #[test]
    fn condensing_mass_sink_closes_and_lowers_e_eff() {
        let tube = slug_with_table(150, 4.0, condensing_table());
        let m_init = tube.total_mass();
        let e_plain = tube.clone().run_bounce().e_eff;
        let mut cb = CondensingBounce::new(tube, 1.0);
        let r = cb.run();

        // The sink is the only mass remover: the tube's lost mass == the tallied stuck mass.
        assert_relative_eq!(
            m_init - cb.tube.total_mass(),
            r.stuck_mass,
            max_relative = 1e-10
        );
        assert!(r.stuck_mass > 0.0, "no condensate stuck");
        assert!(r.loss_condensation > 0.0);
        assert!(
            r.bounce.e_eff < e_plain,
            "α=1 did not lower e_eff: {} vs {e_plain}",
            r.bounce.e_eff
        );
        assert!(
            r.bounce.e_eff > 0.0,
            "e_eff went non-physical: {}",
            r.bounce.e_eff
        );
    }

    /// The condensation sink is dt-convergent (it tracks newly-condensed mass, not step count), so
    /// `e_eff(α=1)` is grid-converged under refinement.
    #[test]
    fn condensing_e_eff_is_grid_convergent() {
        let coarse = CondensingBounce::new(slug_with_table(100, 4.0, condensing_table()), 1.0)
            .run()
            .bounce
            .e_eff;
        let fine = CondensingBounce::new(slug_with_table(200, 4.0, condensing_table()), 1.0)
            .run()
            .bounce
            .e_eff;
        assert!(
            (coarse - fine).abs() < 2e-2,
            "condensing e_eff not grid-convergent: {coarse} vs {fine}"
        );
    }

    /// **Wall deposition is conduction-gated (B-flux gate, ADR-0004 amendment).** A hot, dense gas is
    /// vapor (`liquid_frac ≈ 0`) until something cools it. With a cold conducting wall the gas-side
    /// conduction operator cools the wall cell below the condensation threshold, raising `liquid_frac`
    /// so the sink deposits mass; with no wall (the adiabatic upper bound) the cell stays hot and
    /// nothing deposits. This is the sequencing the rung exists to demonstrate: conduction (channel 2)
    /// *drives* deposition (channel 3). Driven by the conduction + condensation substeps directly (no
    /// hydro) to isolate the gating.
    #[allow(clippy::float_cmp)] // adiabatic path deposits *exactly* nothing → 0.0
    #[test]
    fn conduction_gates_wall_deposition() {
        let make = |wall: Option<Solid>| {
            let cells = 40;
            let x: Vec<f64> = (0..=cells).map(|i| i as f64 / cells as f64).collect();
            let rho = vec![3.0; cells]; // dense (ρ > 2): the compression factor is on
            let vel = vec![0.0; cells];
            let pressure = vec![(GAMMA - 1.0) * 3.0 * 5.0; cells]; // T = 5: hot ⇒ vapor (liquid_frac ≈ 0)
            let tube = Tube::with_eos(
                x,
                &rho,
                &vel,
                &pressure,
                cooling_condensing_table(),
                Viscosity::VON_NEUMANN_RICHTMYER,
            );
            CondensingBounce::new_with_wall(tube, 1.0, wall)
        };

        // Cold conducting wall: cooling the wall cell across the threshold drives deposition.
        let mut lossy = make(Some(Solid::new(200, 1.0, 0.05, 0.5, 20.0)));
        for _ in 0..400 {
            lossy.conduction_substep(1e-3);
            lossy.condensation_substep();
        }
        assert!(
            lossy.loss_conductive > 0.0,
            "conduction must remove heat from the gas"
        );
        assert!(
            lossy.stuck_mass > 0.0,
            "cooling the wall cell below T_sat must drive wall deposition"
        );

        // Adiabatic (no wall): the cell stays hot, liquid_frac ≈ 0, nothing deposits.
        let mut adiabatic = make(None);
        for _ in 0..400 {
            adiabatic.conduction_substep(1e-3);
            adiabatic.condensation_substep();
        }
        assert_eq!(
            adiabatic.stuck_mass, 0.0,
            "no conduction ⇒ no cooling ⇒ no deposition (the dormant channel)"
        );
    }

    /// DIAGNOSTIC (ignored): does the real low-v water slug ever condense at the wall (or anywhere)
    /// during an *adiabatic* bounce? Run with `cargo test -p hydro1d -- --ignored --nocapture
    /// diag_lowv`. Tracks the peak liquid fraction reached at the wall cell and across all cells.
    #[test]
    #[ignore = "diagnostic; needs data/tables/water_lowv.json"]
    fn diag_lowv_condensation() {
        use tables::Table;
        let table = Table::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/tables/water_lowv.json"
        ))
        .unwrap();
        let eos = TableEos::new(table);
        for &t0 in &[350.0_f64, 450.0, 600.0] {
            let mut tube = Tube::slug_si(
                200,
                0.32,
                3200.0,
                1.0,
                t0,
                eos.clone(),
                Viscosity::VON_NEUMANN_RICHTMYER,
            );
            let (mut max_wall_lf, mut max_any_lf) = (0.0_f64, 0.0_f64);
            let (mut peak, mut past, mut fold) = (0.0_f64, false, tube.wall_force());
            for _ in 0..(400 * tube.cells() + 10_000) {
                peak = peak.max(fold);
                if fold < 0.5 * peak {
                    past = true;
                }
                if past && fold < 1e-3 * peak {
                    break;
                }
                let dt = tube.stable_dt();
                tube.step(dt);
                max_wall_lf =
                    max_wall_lf.max(tube.eos.liquid_fraction(tube.density(0), tube.energy[0]));
                for j in 0..tube.cells() {
                    let lf = tube.eos.liquid_fraction(tube.density(j), tube.energy[j]);
                    max_any_lf = max_any_lf.max(lf);
                }
                fold = tube.wall_force();
            }
            println!(
                "t0={t0}: peak wall liquid_frac={max_wall_lf:.4}, peak any-cell={max_any_lf:.4}"
            );
        }
    }

    /// The radiation-medium builder (B5a) reads the gas state correctly: `T = e` (here `c_v = 1`),
    /// the volumetric heat capacity `ρ c_v = ρ`, and the per-length opacities `χ = κ ρ` from the
    /// table; plus the moving-mesh geometry (`dx`, `N−1` center spacings).
    #[test]
    fn radiation_fields_builds_temperature_cv_and_chi() {
        let cells = 4;
        let x: Vec<f64> = (0..=cells).map(|i| i as f64 / cells as f64).collect();
        let rho = vec![2.0; cells];
        let vel = vec![0.0; cells];
        let pressure = vec![8.0; cells]; // T = e = p/((γ−1)ρ) = 8/(0.4·2) = 10, inside the grid
        let tube = Tube::with_eos(
            x,
            &rho,
            &vel,
            &pressure,
            gas_table(0.7, 0.3),
            Viscosity::VON_NEUMANN_RICHTMYER,
        );

        let rf = tube.radiation_fields();
        let (r, t) = (2.0_f64, 10.0_f64);
        for j in 0..cells {
            assert_relative_eq!(rf.temp[j], t, max_relative = 1e-9);
            assert_relative_eq!(rf.cv_vol[j], r * 1.0, max_relative = 1e-6); // c_v = ∂e/∂T = 1
            assert_relative_eq!(
                rf.chi_planck[j],
                r * 0.3 * r * t.powf(-2.0),
                max_relative = 1e-9
            );
            assert_relative_eq!(
                rf.chi_ross[j],
                r * 0.7 * r.powf(2.0) * t.powf(-3.5),
                max_relative = 1e-9
            );
        }
        assert_eq!(rf.dx.len(), cells);
        assert_eq!(rf.center_spacing.len(), cells - 1);
        assert_relative_eq!(rf.dx[0], 1.0 / cells as f64, max_relative = 1e-12);
        assert_relative_eq!(
            rf.center_spacing[0],
            1.0 / cells as f64,
            max_relative = 1e-12
        );
    }

    /// **Radiation-off regression (B5b gate).** Matter exchanges energy with the field only through
    /// the Planck (emission/absorption) opacity `κ_P`, so killing `κ_P → 0` decouples the gas while
    /// a normal Rosseland `κ_R` keeps the diffusion solve well-conditioned (the radiation still
    /// diffuses and leaks, it just never touches the matter). With no conducting wall either, the
    /// Lie-split coupled bounce must then reproduce the pure-hydro [`Tube::run_bounce`] restitution
    /// to round-off — proving the operator split injects nothing spurious into the hydro path.
    #[test]
    fn radiation_off_matches_pure_hydro_bounce() {
        let (cells, mach) = (200, 5.0);
        let consts = RadConstants { c: 1.0, a: 1e-3 };
        let tube = slug_with_table(cells, mach, gas_table(0.7, 1e-12));
        let reference = tube.clone().run_bounce();
        let coupled = CoupledBounce::new(tube, None, consts, Limiter::LevermorePomraning).run();
        assert_relative_eq!(coupled.bounce.e_eff, reference.e_eff, max_relative = 1e-6);
    }

    /// **Energy balance (B5b gate).** In a closed static box (both walls, reflecting radiation
    /// ends, no conduction) with the radiation pushed out of equilibrium, matter + radiation
    /// energy is conserved to machine precision: the FLD substep exchanges exactly the energy the
    /// radiation gains/loses to the local coupling, and the uniform box never moves (zero `PdV`
    /// work), so internal energy changes only through that conservative exchange.
    #[test]
    fn coupling_conserves_matter_plus_radiation_energy() {
        let cells = 10;
        let x: Vec<f64> = (0..=cells).map(|i| i as f64 / cells as f64).collect();
        let rho = vec![1.0; cells];
        let vel = vec![0.0; cells];
        let pressure = vec![4.0; cells]; // (γ−1)·ρ·T = 0.4·1·10 = 4 ⇒ T = e = 10
        let tube = Tube::with_eos(
            x,
            &rho,
            &vel,
            &pressure,
            gas_table(1.0, 1.0),
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        let consts = RadConstants { c: 1.0, a: 1e-3 };
        let mut cb = CoupledBounce::new(tube, None, consts, Limiter::LevermorePomraning);
        // Closed box: no radiation escapes either end.
        cb.bc_wall = RadBc::Reflecting;
        cb.bc_space = RadBc::Reflecting;
        // Push radiation below equilibrium so matter ⇄ radiation actually exchanges energy.
        for e in &mut cb.e_rad {
            *e *= 0.5;
        }
        let total0 = cb.matter_internal_energy() + cb.radiation_energy();
        for _ in 0..50 {
            cb.coupled_step(0.01);
        }
        let total1 = cb.matter_internal_energy() + cb.radiation_energy();
        assert_relative_eq!(total1, total0, max_relative = 1e-10);
    }

    /// **Loss direction + non-negative channels (B5b gate).** Turning on radiative coupling (loss
    /// to the wall and to space) and a conducting wall bleeds energy the lossless bounce keeps, so
    /// the rebound — hence `e_eff` — drops. Each loss channel is a non-negative drain, and the
    /// total deficit is strictly positive.
    #[test]
    fn losses_lower_e_eff_and_split_into_nonneg_channels() {
        let (cells, mach) = (200, 5.0);
        let consts = RadConstants { c: 1.0, a: 1e-3 };
        let lossless = CoupledBounce::new(
            slug_with_table(cells, mach, gas_table(0.7, 1e-12)), // κ_P→0 decouples matter
            None,
            consts,
            Limiter::LevermorePomraning,
        )
        .run();
        let wall = Solid::new(400, 1.0, 0.0, 0.5, 2.0);
        let lossy = CoupledBounce::new(
            slug_with_table(cells, mach, gas_table(1.0, 1.0)),
            Some(wall),
            consts,
            Limiter::LevermorePomraning,
        )
        .run();

        // Losses bleed internal energy ⇒ a weaker rebound ⇒ lower restitution.
        assert!(
            lossy.bounce.e_eff < lossless.bounce.e_eff,
            "lossy e_eff {} should be below lossless {}",
            lossy.bounce.e_eff,
            lossless.bounce.e_eff
        );
        // Each channel is a non-negative drain, and at least one fired.
        assert!(lossy.loss_radiative_wall >= 0.0);
        assert!(lossy.loss_escape_space >= 0.0);
        assert!(lossy.loss_conductive >= 0.0);
        assert!(lossy.loss_radiative_wall + lossy.loss_escape_space + lossy.loss_conductive > 0.0);
    }

    /// **No over-drain (B-flux gate).** A single coupled conduction substep against a cold, highly
    /// effusive wall removes only a bounded share of the thin wall cell's energy: the gas-side
    /// resistance throttles the interface flux. The old gas-temperature-Dirichlet coupling drained the
    /// cell to zero in one step for exactly this regime (the B5d-1 `e_eff = −0.97` collapse) — even
    /// with a deliberately large `dt`, the operator leaves the wall cell hot and positive.
    #[test]
    fn coupled_conduction_does_not_over_drain_the_wall_cell() {
        let cells = 50;
        let x: Vec<f64> = (0..=cells).map(|i| i as f64 / cells as f64).collect();
        let rho = vec![1.0; cells];
        let vel = vec![0.0; cells];
        let pressure = vec![(GAMMA - 1.0) * 1.0 * 100.0; cells]; // T = e = 100 (hot, inside the grid)
        let tube = Tube::with_eos(
            x,
            &rho,
            &vel,
            &pressure,
            gas_table(0.7, 1e-12),
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        let consts = RadConstants { c: 1.0, a: 1e-3 };
        let wall = Solid::new(400, 1.0, 0.0, 0.5, 50.0); // cold (T=0), effusivity 50/√0.5 ≫ the gas's
        let mut cb = CoupledBounce::new(tube, Some(wall), consts, Limiter::LevermorePomraning);

        let e0 = cb.tube.energy[0];
        cb.conduction_substep(1.0); // a deliberately large dt — the over-drain regime
        assert!(
            cb.tube.energy[0] > 0.0 && cb.tube.energy[0] < e0,
            "wall cell should cool but never drain to zero: {} (was {e0})",
            cb.tube.energy[0]
        );
        assert!(
            cb.loss_conductive > 0.0 && cb.loss_conductive.is_finite(),
            "conduction must bleed a finite, positive amount: {}",
            cb.loss_conductive
        );
    }
}
