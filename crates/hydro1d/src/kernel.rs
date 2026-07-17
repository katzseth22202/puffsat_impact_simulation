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
    ///
    /// Includes the artificial viscosity `q`, whose first-impact spike `≈ c_q·ρv²` dominates
    /// this peak (it tracks `c_q`, not the physics). Use [`Self::peak_wall_pressure`] for the
    /// physical load on the plate.
    pub peak_wall_force: f64,
    /// Peak **physical** wall pressure: max over the bounce of the wall cell's EOS pressure
    /// `p(0)`, excluding artificial viscosity. Converges under refinement to the reflected-shock
    /// stagnation pressure `≈ (γ_eff+1)/2 · ρv²` (≈1.1·ρv² for the water EOS at 11–16 km/s) —
    /// the correct facesheet survivability load (ADR-0010 correction).
    pub peak_wall_pressure: f64,
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

/// Shared wall-impulse integration loop for a bounce (ADR-0001): peak-detect the wall force, then
/// integrate `J_wall` trapezoidally until it decays to `10⁻³` of its peak. [`Tube::run_bounce`],
/// [`CoupledBounce::run`], [`CondensingBounce::run`], and [`AblatingBounce::run`] each wrap
/// [`Self::run_bounce_loop`] — they differ only in what one `step` does (bare hydro vs.
/// radiation/conduction/condensation/ablation substeps), not in how the impulse is accumulated or
/// the run is terminated. [`Tube::run_stick_bounce`] and [`Tube::run_bounce_frozen_rebound`] keep
/// their own loops: their termination is a stagnation crossing, not this tail guard.
trait BounceStepper {
    /// Axial force the gas exerts on the wall this instant (includes artificial viscosity).
    fn wall_force(&self) -> f64;
    /// Physical wall pressure this instant (excludes artificial viscosity).
    fn wall_pressure(&self) -> f64;
    /// CFL-stable timestep for the current state.
    fn stable_dt(&self) -> f64;
    /// Total gas momentum (signed).
    fn total_momentum(&self) -> f64;
    /// Cell count, for the `max_steps` safety cap.
    fn cells(&self) -> usize;
    /// Advance the state by one step of `dt` — the part that varies per bounce variant.
    fn step(&mut self, dt: f64);

    /// Fire the slug at the wall and integrate to the `10⁻³`-of-peak tail guard, returning the
    /// wall impulse and the restitution it implies.
    fn run_bounce_loop(&mut self) -> BounceResult {
        let incident = self.total_momentum().abs();
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut peak_pressure: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = self.wall_force();
        let max_steps = 400 * self.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            peak_pressure = peak_pressure.max(self.wall_pressure());
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

        BounceResult {
            wall_impulse,
            incident_momentum: incident,
            residual_momentum: self.total_momentum(),
            e_eff: wall_impulse / incident - 1.0,
            peak_wall_force: peak,
            peak_wall_pressure: peak_pressure,
        }
    }
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

    /// Physical wall pressure: the EOS pressure `p` of cell 0 alone, excluding the artificial
    /// viscosity `q`. The impulse must integrate `p + q` (AV carries real momentum flux through
    /// the smeared shock), but the *peak* of `p + q` is a numerical artifact `≈ c_q·ρv²`; the
    /// peak of `p` is the physical stagnation load.
    fn wall_pressure(&self) -> f64 {
        self.pressure(0)
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
        self.run_bounce_loop()
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
        let mut peak_pressure: f64 = 0.0;
        let mut p_old = p_initial;
        let mut force_old = self.wall_force();
        let max_steps = 400 * self.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            peak_pressure = peak_pressure.max(self.wall_pressure());
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
                    peak_wall_pressure: peak_pressure,
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
            peak_wall_pressure: peak_pressure,
        }
    }
}

impl<E: Eos> BounceStepper for Tube<E> {
    fn wall_force(&self) -> f64 {
        Tube::wall_force(self)
    }
    fn wall_pressure(&self) -> f64 {
        Tube::wall_pressure(self)
    }
    fn stable_dt(&self) -> f64 {
        Tube::stable_dt(self)
    }
    fn total_momentum(&self) -> f64 {
        Tube::total_momentum(self)
    }
    fn cells(&self) -> usize {
        Tube::cells(self)
    }
    fn step(&mut self, dt: f64) {
        Tube::step(self, dt);
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
    conduction_into_wall_scaled(tube, wall, dt, 1.0)
}

/// As [`conduction_into_wall`], but with the conducted heat attenuated by the **blowing factor**
/// `phi ∈ (0, 1]` (Rung E, E2): the injected vapor curtain intercepts a fraction `(1 − φ)` of the
/// heat the gas would conduct to the plate and convects it back into the gas, so the gas net-cools by
/// `φ·ΔT` and the plate receives `φ` of the flux. `phi = 1` is the unblown conduction the rigid-wall
/// bounces use (an exact no-op of the attenuation).
///
/// Attenuating the *delivered heat* — rather than the gas-side conductance — is the physically right
/// lever: over a full bounce the integrated conductive loss is rate-insensitive (the near-wall gas
/// dumps its enthalpy into the deep solid sink whatever the conductance), so only intercepting the
/// heat actually reduces channel 2.
fn conduction_into_wall_scaled(
    tube: &mut Tube<TableEos>,
    wall: &mut Option<Solid>,
    dt: f64,
    phi: f64,
) -> f64 {
    let Some(fields) = tube.gas_conduction_fields() else {
        return 0.0;
    };
    let Some(solid) = wall.as_mut() else {
        return 0.0;
    };
    let step = solid.step_coupled(&fields.state(), dt);
    for (j, &dtemp) in step.gas_dtemp.iter().enumerate() {
        tube.energy[j] = (phi * fields.cv[j]).mul_add(dtemp, tube.energy[j]).max(0.0);
    }
    phi * step.interface_flux * dt
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

/// Outcome of a sudden-freeze bounce ([`Tube::run_bounce_frozen_rebound`]): the usual
/// [`BounceResult`] plus the freeze-instant diagnostics of the frozen-recombination bounding run.
#[derive(Debug, Clone, Copy)]
pub struct FreezeBounceResult {
    /// The usual restitution/impulse bookkeeping over the *whole* bounce (both EOS phases).
    pub bounce: BounceResult,
    /// Mass-weighted mean density at the freeze instant (global momentum zero) [kg/m³];
    /// `0` if the gas never turned around within the step cap.
    pub rho_star: f64,
    /// Mass-weighted mean temperature at the freeze instant [K]; `0` if never turned around.
    pub t_star: f64,
    /// Total internal-energy jump `Σ m_j (e_new − e_old)` across the EOS-swap re-seed (per unit
    /// wall area) — the splice-consistency diagnostic. Exactly the per-cell mismatch between the
    /// frozen table (one composition for the whole slug) and each cell's own state at the freeze
    /// instant; `0` when no swap happened.
    pub swap_energy_jump: f64,
}

impl Tube<TableEos> {
    /// Mass-weighted mean `(ρ, T)` over the slug — the freeze reference state recorded at
    /// turnaround (the probe output the frozen-composition table is generated from).
    fn mass_weighted_state(&self) -> (f64, f64) {
        let mut m_tot = 0.0;
        let mut rho_mean = 0.0;
        let mut t_mean = 0.0;
        for j in 0..self.cells() {
            let m = self.mass[j];
            let rho = self.density(j);
            m_tot += m;
            rho_mean += m * rho;
            t_mean += m * self.eos.temperature(rho, self.energy[j]);
        }
        (rho_mean / m_tot, t_mean / m_tot)
    }

    /// Swap the tube's EOS for `new_eos`, re-seeding every cell's specific internal energy at its
    /// **current temperature**: `e_j → e_new(ρ_j, T_j)` with `T_j` inverted from the old EOS. This
    /// is the sudden-freeze splice: temperature (and, to the accuracy of the shared freeze
    /// composition, pressure) is continuous across the swap, and any chemical energy the new
    /// table carries as a constant offset is inert thereafter. Cells whose temperature inversion
    /// clamps at a grid edge (vacuum-cooled free-surface tail, including `e = 0` floored cells)
    /// carry their out-of-table energy *deficit/excess* across the tables' zero-point shift at
    /// that edge — they must not keep their raw energy, or a new table with a large chemical
    /// offset would leave them far below its energy floor, where the pressure response is flat
    /// (temperature clamped) and neighbouring cells crush them without limit (a `dt → 0`
    /// collapse). Returns the total energy jump `Σ m_j (e_new − e_old)` as the consistency
    /// diagnostic.
    fn swap_eos_reseed_temperature(&mut self, new_eos: TableEos) -> f64 {
        let t_grid = self.eos.table().t_grid();
        let (t_lo, t_hi) = (t_grid[0], t_grid[t_grid.len() - 1]);
        let mut jump = 0.0;
        for j in 0..self.cells() {
            let rho = self.density(j);
            let e_old = self.energy[j];
            let t = self.eos.temperature(rho, e_old);
            let e_new = if t <= t_lo {
                new_eos.table().energy(rho, t_lo) - (self.eos.table().energy(rho, t_lo) - e_old)
            } else if t >= t_hi {
                new_eos.table().energy(rho, t_hi) + (e_old - self.eos.table().energy(rho, t_hi))
            } else {
                new_eos.table().energy(rho, t)
            };
            jump += self.mass[j] * (e_new - e_old);
            self.energy[j] = e_new;
        }
        self.eos = new_eos;
        jump
    }

    /// [`Tube::run_bounce`] with a **sudden freeze at turnaround** — the frozen-recombination
    /// bounding run (equilibrium in, frozen out; the classic nozzle-flow sudden-freeze
    /// approximation applied at the instant of global stagnation, where the chemical store is
    /// maximal).
    ///
    /// Compression runs on the tube's own (equilibrium) EOS. At the first instant the total gas
    /// momentum crosses zero, the mass-weighted `(ρ*, T*)` is recorded and — if `frozen` is
    /// `Some` — the EOS is swapped for the frozen-composition table via the temperature-continuous
    /// re-seed, so the rebound returns **no** chemical (dissociation/ionization) energy. With
    /// `frozen = None` this is a plain bounce that additionally reports `(ρ*, T*)` (the probe
    /// mode that feeds the frozen-table generator).
    pub fn run_bounce_frozen_rebound(&mut self, frozen: Option<TableEos>) -> FreezeBounceResult {
        let p_initial = self.total_momentum();
        let incident = p_initial.abs();
        let mut frozen = frozen;
        let mut wall_impulse = 0.0;
        let mut peak: f64 = 0.0;
        let mut peak_pressure: f64 = 0.0;
        let mut past_peak = false;
        let mut turned_around = false;
        let mut rho_star = 0.0;
        let mut t_star = 0.0;
        let mut swap_energy_jump = 0.0;
        let mut force_old = self.wall_force();
        let mut p_old = p_initial;
        let max_steps = 400 * self.cells() + 10_000;

        for _ in 0..max_steps {
            peak = peak.max(force_old);
            peak_pressure = peak_pressure.max(self.wall_pressure());
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

            let p_new = self.total_momentum();
            if !turned_around && p_old < 0.0 && p_new >= 0.0 {
                turned_around = true;
                (rho_star, t_star) = self.mass_weighted_state();
                if let Some(new_eos) = frozen.take() {
                    swap_energy_jump = self.swap_eos_reseed_temperature(new_eos);
                    // The wall force is (mildly) discontinuous across the swap; re-read it so the
                    // next trapezoid uses the frozen-EOS pressure.
                    force_old = self.wall_force();
                }
            }
            p_old = p_new;
        }

        let residual = self.total_momentum();
        FreezeBounceResult {
            bounce: BounceResult {
                wall_impulse,
                incident_momentum: incident,
                residual_momentum: residual,
                e_eff: wall_impulse / incident - 1.0,
                peak_wall_force: peak,
                peak_wall_pressure: peak_pressure,
            },
            rho_star,
            t_star,
            swap_energy_jump,
        }
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
        let bounce = self.run_bounce_loop();
        CoupledBounceResult {
            bounce,
            loss_radiative_wall: self.loss_radiative_wall,
            loss_escape_space: self.loss_escape_space,
            loss_conductive: self.loss_conductive,
        }
    }
}

impl BounceStepper for CoupledBounce {
    fn wall_force(&self) -> f64 {
        self.tube.wall_force()
    }
    fn wall_pressure(&self) -> f64 {
        self.tube.wall_pressure()
    }
    fn stable_dt(&self) -> f64 {
        self.tube.stable_dt()
    }
    fn total_momentum(&self) -> f64 {
        self.tube.total_momentum()
    }
    fn cells(&self) -> usize {
        self.tube.cells()
    }
    fn step(&mut self, dt: f64) {
        self.coupled_step(dt);
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

    /// One condensing step: hydro, then conduction (channel 2), then wall deposition (channel 3).
    /// Conduction first: the cold wall cools the near-wall gas below `T_sat`, raising
    /// `liquid_frac`, which the deposition sink then sticks (the two channels are sequenced).
    fn condensing_step(&mut self, dt: f64) {
        self.tube.step(dt);
        self.conduction_substep(dt);
        self.condensation_substep();
    }

    /// Fire the condensing slug at the wall, integrating to the same `10⁻³`-of-peak tail guard as
    /// [`Tube::run_bounce`], and return the restitution plus the condensation loss.
    pub fn run(&mut self) -> CondensingBounceResult {
        let bounce = self.run_bounce_loop();
        CondensingBounceResult {
            bounce,
            loss_condensation: self.loss_condensation,
            loss_conductive: self.loss_conductive,
            stuck_mass: self.stuck_mass,
        }
    }
}

impl BounceStepper for CondensingBounce {
    fn wall_force(&self) -> f64 {
        self.tube.wall_force()
    }
    fn wall_pressure(&self) -> f64 {
        self.tube.wall_pressure()
    }
    fn stable_dt(&self) -> f64 {
        self.tube.stable_dt()
    }
    fn total_momentum(&self) -> f64 {
        self.tube.total_momentum()
    }
    fn cells(&self) -> usize {
        self.tube.cells()
    }
    fn step(&mut self, dt: f64) {
        self.condensing_step(dt);
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
    /// Blowing-correction coefficient (E2; `0` = off). The injected vapor thickens the boundary
    /// layer, attenuating the conducted heat delivered to the plate by `φ = 1/(1 + B)` with the
    /// dimensionless blowing rate `B = blowing · ablated_mass / m_cloud`. Cuts the conductive wall
    /// loss where conduction is active; **null at the high-v anchors** (no `k_gas`, conduction off).
    blowing: f64,
    /// Vapor gray opacity `κ_vapor` [m²/kg] for shielding (E3; `0` = off). The ablated vapor forms a
    /// near-wall absorbing curtain of optical depth `τ_vapor = κ_vapor · ablated_mass` and transmission
    /// `1/(1 + τ_vapor)`; it throttles the wall's radiative conductance ([`RadBc::MarshakAttenuated`]),
    /// retaining the intercepted radiation in the near-wall field (which couples back into the gas)
    /// rather than letting it drain to the cold plate. This cuts the radiative wall loss (channel 1a)
    /// and raises `e_eff` — the **live recovery lever at the transitional dip** (ADR-0012/0014).
    kappa_vapor: f64,
}

impl Ablation {
    /// An ablation model with heat of ablation `q_star` [J/kg] and cold-vapor injection temperature
    /// `t_vapor` [K]. Blowing (E2) is off by default; set it with [`Ablation::with_blowing`].
    ///
    /// # Panics
    /// Panics unless `q_star > 0` and `t_vapor > 0`.
    #[must_use]
    pub fn new(q_star: f64, t_vapor: f64) -> Self {
        assert!(
            q_star > 0.0 && t_vapor > 0.0,
            "Q* and t_vapor must be positive"
        );
        Self {
            q_star,
            t_vapor,
            blowing: 0.0,
            kappa_vapor: 0.0,
        }
    }

    /// Set the vapor gray opacity `κ_vapor` [m²/kg] for shielding (E3). `> 0` builds a near-wall
    /// absorbing layer (`τ_vapor = κ_vapor · ablated_mass`) that intercepts incoming radiation before
    /// the plate; `0` recovers the bare-wall radiative loss.
    ///
    /// # Panics
    /// Panics if `kappa_vapor < 0`.
    #[must_use]
    pub fn with_vapor_opacity(mut self, kappa_vapor: f64) -> Self {
        assert!(kappa_vapor >= 0.0, "vapor opacity must be non-negative");
        self.kappa_vapor = kappa_vapor;
        self
    }

    /// Set the blowing-correction coefficient (E2). `> 0` makes the vapor curtain thicken as ablation
    /// proceeds, monotonically cutting the conductive wall loss; `0` recovers the unblown conduction.
    ///
    /// # Panics
    /// Panics if `coeff < 0`.
    #[must_use]
    pub fn with_blowing(mut self, coeff: f64) -> Self {
        assert!(coeff >= 0.0, "blowing coefficient must be non-negative");
        self.blowing = coeff;
        self
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
    /// Initial cloud mass per area `ρ·L` — the reference for the dimensionless blowing rate (E2).
    m_ref: f64,
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
        let m_ref: f64 = tube.mass.iter().sum();
        Self {
            tube,
            e_rad,
            wall,
            consts,
            limiter,
            bc_wall: RadBc::Marshak(0.0),
            bc_space: RadBc::Marshak(0.0),
            ablation,
            m_ref,
            loss_radiative_wall: 0.0,
            loss_escape_space: 0.0,
            loss_conductive: 0.0,
            ablated_mass: 0.0,
            loss_ablation: 0.0,
        }
    }

    /// One implicit gray-FLD substep (as [`CoupledBounce::radiation_substep`]), returning the
    /// **radiative wall flux that reaches the plate this step** — the radiative part of the `q_in`
    /// that drives ablation.
    ///
    /// **Vapor shielding (E3):** the ablated vapor forms a near-wall absorbing curtain of optical
    /// depth `τ_vapor = κ_vapor · ablated_mass` (lagged on the ablation so far, no within-step
    /// circularity) and transmission `trans = 1/(1 + τ_vapor)`. Rather than draining the full Marshak
    /// flux and re-injecting the intercepted part (which over-pressurizes the thin wall cell and
    /// destabilizes the bounce), the curtain is folded into the boundary condition itself
    /// ([`RadBc::MarshakAttenuated`]): the wall's radiative conductance is scaled by `trans`, so the
    /// `(1 − trans)` it intercepts is *retained in the near-wall radiation field* by the implicit solve
    /// and couples back into the gas self-consistently (raising the near-wall pressure → bounce). This
    /// cuts the radiative wall loss (channel 1a) and raises `e_eff`, energy-conservingly and stably.
    fn radiation_substep(&mut self, dt: f64) -> f64 {
        let e_inc = match self.bc_wall {
            RadBc::Marshak(e) => e,
            _ => 0.0,
        };
        let tau_v = self.ablation.kappa_vapor * self.ablated_mass;
        let transmission = 1.0 / (1.0 + tau_v);
        let bc_wall = RadBc::MarshakAttenuated {
            e_inc,
            transmission,
        };

        let fields = self.tube.radiation_fields();
        let delta_e = fld_substep(
            &fields.medium(),
            &mut self.e_rad,
            bc_wall,
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
        // The true wall loss is the attenuated flux through the throttled surface.
        let wall_flux = dt * transmission * 0.5 * c * (self.e_rad[0] - e_inc);
        self.loss_radiative_wall += wall_flux;
        if let RadBc::Marshak(e_inc) = self.bc_space {
            let last = self.tube.cells() - 1;
            self.loss_escape_space += dt * 0.5 * c * (self.e_rad[last] - e_inc);
        }
        wall_flux
    }

    /// The blowing-reduction factor `φ = 1/(1 + B)` (E2, ADR-0014), with the dimensionless blowing
    /// rate `B = blowing · ablated_mass / m_cloud` lagged on the ablation accumulated so far (so there
    /// is no within-step circularity). `φ = 1` when blowing is off or before any ablation.
    fn blowing_factor(&self) -> f64 {
        let b = self.ablation.blowing * self.ablated_mass / self.m_ref;
        1.0 / (1.0 + b)
    }

    /// Conductive wall loss (channel 2), via the shared gas-side operator with the **blowing
    /// correction** (E2): the conducted heat is attenuated by `φ` so the injected vapor curtain
    /// intercepts the flux into the plate. Returns the (blown) flux conducted this step — the
    /// conductive part of `q_in`. `0` when `wall = None` or the table carries no `k_gas` (the high-v
    /// anchors, where blowing is therefore null).
    fn conduction_substep(&mut self, dt: f64) -> f64 {
        let phi = self.blowing_factor();
        let q = conduction_into_wall_scaled(&mut self.tube, &mut self.wall, dt, phi);
        self.loss_conductive += q;
        q
    }

    /// Ablation substep: the quasi-steady surface energy balance. The incoming wall flux `q_in`
    /// (energy/area this step) boils off `ṁ·dt = q_in/Q*` of ablator, injected as **cold vapor** into
    /// the near-wall layer — the inverse of [`CondensingBounce`]'s wall sink. The vapor enters at the
    /// wall node (`u = 0`), adding no axial momentum directly; it acts on `e_eff` through the near-wall
    /// pressure it raises (E1) and by shielding radiation (E3).
    ///
    /// The injected mass is **spread mass-weighted over the first `K = cells/20` cells** rather than
    /// dumped into cell 0. Concentrating every step's `dm` in the single wall cell quintuples its mass
    /// over the bounce, spiking its density and collapsing the CFL step at the dense corner (the run
    /// then starves on `max_steps`). Spreading over the thin near-wall layer keeps the same total
    /// injected mass and ablation enthalpy — so E1's mass/energy closure is unchanged — while keeping
    /// each cell's density bounded.
    fn ablation_substep(&mut self, q_in: f64) {
        if q_in <= 0.0 {
            return;
        }
        let dm = q_in / self.ablation.q_star;
        if dm <= 0.0 {
            return;
        }
        let k = (self.tube.cells() / 20).max(1);
        let m_layer: f64 = self.tube.mass[..k].iter().sum();
        if m_layer <= 0.0 {
            return;
        }
        for j in 0..k {
            let dm_j = dm * self.tube.mass[j] / m_layer;
            let rho_j = self.tube.density(j);
            let e_vapor = self.tube.eos.table().energy(rho_j, self.ablation.t_vapor);
            let m_old = self.tube.mass[j];
            let e_old = self.tube.energy[j];
            self.tube.mass[j] = m_old + dm_j;
            self.tube.energy[j] = dm_j.mul_add(e_vapor, m_old * e_old) / (m_old + dm_j);
        }
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
        let bounce = self.run_bounce_loop();
        AblatingBounceResult {
            bounce,
            loss_radiative_wall: self.loss_radiative_wall,
            loss_escape_space: self.loss_escape_space,
            loss_conductive: self.loss_conductive,
            ablated_mass: self.ablated_mass,
            loss_ablation: self.loss_ablation,
        }
    }
}

impl BounceStepper for AblatingBounce {
    fn wall_force(&self) -> f64 {
        self.tube.wall_force()
    }
    fn wall_pressure(&self) -> f64 {
        self.tube.wall_pressure()
    }
    fn stable_dt(&self) -> f64 {
        self.tube.stable_dt()
    }
    fn total_momentum(&self) -> f64 {
        self.tube.total_momentum()
    }
    fn cells(&self) -> usize {
        self.tube.cells()
    }
    fn step(&mut self, dt: f64) {
        self.ablating_step(dt);
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

    /// DIAGNOSTIC (ignored): trace the wall cell through the coupled bounce that used to collapse
    /// (the 2026-07-16 resolution-onset finding — rho = 0.849, 16 km/s, 600 cells collapsed while
    /// 300 was healthy). Prints the wall cell's `(e, T, rho, p, E_rad, aT⁴)` trajectory and counts
    /// energy-floor hits. This trace *disproved* the "radiative over-drain to the energy floor"
    /// hypothesis (zero floor hits; the exchange tracks LTE cleanly) and pinned the real cause:
    /// the radiatively-cooled wall cell was compressed past the old table's ρ = 20 kg/m³ ceiling,
    /// where the clamped `p(ρ)` no longer arrests the Lagrangian compression (ρ ran to ~2.8e4,
    /// cell width → 0, dt → 0, run stalled mid-infall). With the extended table (ceiling
    /// 1000 kg/m³) the compression self-arrests near ρ ≈ 27 and the bounce completes on the
    /// plateau (`e_eff ≈ 0.638`, zero floor hits). Run with
    /// `cargo test -p hydro1d --release -- --ignored --nocapture diag_radiative`.
    #[test]
    #[ignore = "diagnostic; needs data/tables/water.json"]
    fn diag_radiative_collapse_wall_cell_trace() {
        use crate::eos::Eos as _;
        let table = Table::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/tables/water.json"
        ))
        .unwrap();
        let consts = RadConstants {
            c: 2.997_924_58e8,
            a: 7.565_733e-16,
        };
        let tube = Tube::slug_si(
            600,
            0.849,
            16_000.0,
            1.5,
            400.0,
            TableEos::new(table),
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        let mut cb = CoupledBounce::new(tube, None, consts, Limiter::LevermorePomraning);
        let incident = cb.tube.total_momentum().abs();
        let mut peak: f64 = 0.0;
        let mut past_peak = false;
        let mut force_old = cb.tube.wall_force();
        let mut wall_impulse = 0.0;
        let mut floor_hits = 0usize;
        let mut t_sim = 0.0;
        let max_steps = 400 * cb.tube.cells() + 10_000;
        for step in 0..max_steps {
            peak = peak.max(force_old);
            if force_old < 0.5 * peak {
                past_peak = true;
            }
            if past_peak && force_old < 1e-3 * peak {
                println!("terminated by tail guard at step {step}");
                break;
            }
            let dt = cb.tube.stable_dt();
            t_sim += dt;
            cb.coupled_step(dt);
            if cb.tube.energy[0] <= 0.0 {
                floor_hits += 1;
            }
            let force_new = cb.tube.wall_force();
            wall_impulse += 0.5 * dt * (force_old + force_new);
            force_old = force_new;
            if step % 2_000 == 0 || (step < 200 && step % 20 == 0) {
                let e0 = cb.tube.energy[0];
                let rho0 = cb.tube.density(0);
                let t0 = cb.tube.eos.temperature(rho0, e0.max(0.0));
                let at4 = consts.a * t0.powi(4);
                println!(
                    "step {step:>6} t={t_sim:.3e}: wall e={e0:.3e} T={t0:.0} rho={rho0:.2} \
                     p={:.3e} E_rad={:.3e} aT4={at4:.3e} mom/inc={:+.3} floor_hits={floor_hits}",
                    cb.tube.pressure(0),
                    cb.e_rad[0],
                    cb.tube.total_momentum() / incident,
                );
            }
        }
        println!(
            "final: e_eff={:.4} floor_hits={floor_hits} (of wall cell)",
            wall_impulse / incident - 1.0
        );
    }

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

    // ---- Sudden-freeze splice (frozen-recombination bounding run) ------------------------------

    /// An ideal-gas table with a tunable heat capacity and an inert energy offset:
    /// `e = cv·T + off`, `p = (γ−1)ρT`, `c_s = √(γ(γ−1)T)`, on `n` log-spaced nodes per axis.
    /// With `off = 0` every field is a power law in `(ρ, T)`, so the log-log interpolation is
    /// exact; a nonzero `off` breaks that exactness, which is why the offset test uses a dense
    /// grid. `cv` rescales the thermal energy pool at fixed pressure — the splice direction knob.
    fn gas_table_cv_offset(cv: f64, off: f64, n: usize) -> TableEos {
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64)) // 0.01 … 10
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64)) // 0.05 … 200
            .collect();
        let (mut p, mut e, mut cs) = (Vec::new(), Vec::new(), Vec::new());
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(cv * t + off);
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
            }
        }
        let one = vec![1e-10; n * n];
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": {
                "p": p, "e": e, "c_s": cs,
                "kappa_rosseland": one, "kappa_planck": one,
            },
        });
        TableEos::new(Table::from_json(&json.to_string()).unwrap())
    }

    /// With `frozen = None` the sudden-freeze runner is a plain [`Tube::run_bounce`] that
    /// additionally reports a physical turnaround state: same restitution, `ρ*` above the initial
    /// density (the slug is shock-compressed at stagnation), `T*` above the cold inflow.
    #[allow(clippy::float_cmp)] // exact: no swap happens, so the jump is the literal 0.0
    #[test]
    fn frozen_rebound_probe_matches_plain_bounce() {
        let table = gas_table_cv_offset(1.0, 0.0, 8);
        let mach = 5.0;
        let t0 = 1.0 / (GAMMA * (GAMMA - 1.0) * mach * mach);
        let plain = slug_with_table(200, mach, table.clone()).run_bounce();
        let probe = slug_with_table(200, mach, table).run_bounce_frozen_rebound(None);

        assert_relative_eq!(probe.bounce.e_eff, plain.e_eff, max_relative = 1e-12);
        assert_relative_eq!(
            probe.bounce.wall_impulse,
            plain.wall_impulse,
            max_relative = 1e-12
        );
        assert!(
            probe.rho_star > 1.0,
            "turnaround slug should be compressed above ρ₀ = 1: {}",
            probe.rho_star
        );
        assert!(
            probe.t_star > t0,
            "turnaround slug should be shock-heated above T₀ = {t0}: {}",
            probe.t_star
        );
        assert_eq!(probe.swap_energy_jump, 0.0);
    }

    /// Swapping in a clone of the *same* table at turnaround is a no-op: the temperature re-seed
    /// round-trips `e → T → e` exactly (the inversion is the analytic inverse of the forward
    /// interpolation), so the restitution is unchanged and the energy jump vanishes.
    #[test]
    fn frozen_rebound_swap_same_table_is_noop() {
        let table = gas_table_cv_offset(1.0, 0.0, 8);
        let plain = slug_with_table(200, 5.0, table.clone()).run_bounce();
        let swapped =
            slug_with_table(200, 5.0, table.clone()).run_bounce_frozen_rebound(Some(table));

        assert_relative_eq!(swapped.bounce.e_eff, plain.e_eff, max_relative = 1e-9);
        let e_scale = swapped.bounce.incident_momentum; // O(ρLv): a per-area energy/momentum scale
        assert!(
            swapped.swap_energy_jump.abs() < 1e-9 * e_scale,
            "same-table re-seed should not move energy: {}",
            swapped.swap_energy_jump
        );
    }

    /// A constant energy offset in the frozen table (`e = T + C`) is **inert**: it shifts the
    /// energy zero point but exchanges nothing with the thermal pool, so the bounce is unchanged.
    /// This is the trust the splice puts in the frozen chemical energy `e_chem` — locked, carried,
    /// never returned. Dense grid: the offset breaks the log-log power-law exactness, so agreement
    /// is to interpolation error, not round-off.
    #[test]
    fn frozen_rebound_inert_offset_does_not_change_the_bounce() {
        let n = 192;
        let plain = slug_with_table(200, 5.0, gas_table_cv_offset(1.0, 0.0, n)).run_bounce();
        let offset = slug_with_table(200, 5.0, gas_table_cv_offset(1.0, 0.0, n))
            .run_bounce_frozen_rebound(Some(gas_table_cv_offset(1.0, 20.0, n)));

        assert_relative_eq!(offset.bounce.e_eff, plain.e_eff, max_relative = 2e-3);
        // The re-seed jump is exactly the added offset: C per unit mass, total mass = ρ₀L = 1.
        assert_relative_eq!(offset.swap_energy_jump, 20.0, max_relative = 2e-2);
    }

    /// Splice direction: swapping to a table with a **larger** thermal pool at the same `(p, T)`
    /// (`e = 2T`) strengthens the rebound, a **smaller** pool (`e = T/2`) weakens it — the
    /// mechanics the frozen-composition bound relies on (the frozen EOS holds less *returnable*
    /// energy than the equilibrium one).
    #[test]
    fn frozen_rebound_thermal_pool_sets_the_rebound_direction() {
        let e_plain = slug_with_table(200, 5.0, gas_table_cv_offset(1.0, 0.0, 8))
            .run_bounce()
            .e_eff;
        let e_rich = slug_with_table(200, 5.0, gas_table_cv_offset(1.0, 0.0, 8))
            .run_bounce_frozen_rebound(Some(gas_table_cv_offset(2.0, 0.0, 8)))
            .bounce
            .e_eff;
        let e_poor = slug_with_table(200, 5.0, gas_table_cv_offset(1.0, 0.0, 8))
            .run_bounce_frozen_rebound(Some(gas_table_cv_offset(0.5, 0.0, 8)))
            .bounce
            .e_eff;

        assert!(
            e_rich > e_plain && e_plain > e_poor,
            "rebound should order with the swapped-in thermal pool: rich {e_rich} vs plain \
             {e_plain} vs poor {e_poor}"
        );
    }
}
