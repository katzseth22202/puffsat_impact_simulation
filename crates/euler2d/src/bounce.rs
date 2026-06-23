//! Flat-plate slug bounce and the `eta_capture` extraction (ADR-0003, D2).
//!
//! A finite cold gas slug (radius `r_foot`, length `L`, speed `v` toward the plate at `z = 0`) in
//! near-vacuum — the 2D axisymmetric analogue of the 1D `hydro1d` slug (ADR-0001). The plate
//! reflects for `r ≤ r_plate` and lets gas past its edge escape (§7). The axial wall impulse
//!
//! ```text
//! J_wall = ∫_t ∫_{r ≤ r_plate} p(r,t) dA dt   (n̂·ẑ = 1, flat plate),
//! ```
//!
//! is integrated until the plate force decays to `10⁻³` of its global peak (ADR-0001; the
//! "stay-below-after-the-peak" rule also covers a concave plate's focused secondary peak, which a
//! flat plate does not produce).
//!
//! `eta_capture` is the **lossless 2D/1D wall-impulse ratio** (ADR-0003). We form it from two runs
//! of *this* kernel so the common gas-dynamic re-expansion (and the scheme error) cancel:
//! - **free** — the finite cloud on the finite plate, gas free to spread and escape;
//! - **confined** — the same column as a plane wave (cloud fills the radius, reflecting `r`-walls),
//!   the perfectly-collimated 1D limit.
//!
//! Each run's `J_wall / p_in` is a restitution-like number (`→ 1 + e_eff` for the plane wave);
//! `eta_capture = (J_wall/p_in)_free / (J_wall/p_in)_confined` is then pure geometry, `= 1` for the
//! confined case by construction and `< 1` once the rebound spreads radially or overshoots the
//! plate edge.

use crate::kernel::{Bc, Grid2D};
use crate::state::Prim;

/// A slug-bounce configuration (normalized `ρ₀ = 1`, `v = 1`; incident Mach `M` sets the cold
/// pressure `p₀ = 1/(γM²)`, as in the 1D `hydro1d` slug).
#[derive(Debug, Clone, Copy)]
pub struct SlugConfig {
    /// Adiabatic index γ (the calibrated effective-γ, ADR-0008).
    pub gamma: f64,
    /// Incident Mach number `M = v/c₀`.
    pub mach: f64,
    /// Cloud footprint radius `r_foot`.
    pub r_foot: f64,
    /// Cloud axial length `L`.
    pub length: f64,
    /// Plate radius `r_plate` (reflecting for `r ≤ r_plate`).
    pub r_plate: f64,
    /// Domain size and resolution.
    pub r_max: f64,
    pub z_max: f64,
    pub nr: usize,
    pub nz: usize,
    /// Confined (plane-wave / 1D-limit) run: reflecting outer `r`-wall, cloud fills the radius.
    pub confined: bool,
}

/// The result of one slug bounce.
#[derive(Debug, Clone, Copy)]
pub struct Bounce2D {
    /// Time-integrated axial plate impulse `J_wall`.
    pub wall_impulse: f64,
    /// Incident axial momentum `p_in = |Σ ρu_z dV|` at `t = 0`.
    pub incident_momentum: f64,
    /// Peak axial plate force during the bounce.
    pub peak_force: f64,
    /// Number of time steps taken.
    pub steps: usize,
}

impl Bounce2D {
    /// `J_wall / p_in` — the restitution-like throughput (`→ 1 + e_eff` for the plane wave).
    #[must_use]
    pub fn restitution_ratio(&self) -> f64 {
        self.wall_impulse / self.incident_momentum
    }
}

/// `eta_capture = (J_wall/p_in)_free / (J_wall/p_in)_confined` — the lossless 2D/1D ratio (ADR-0003),
/// pure capture geometry with the common re-expansion divided out.
#[must_use]
pub fn eta_capture(free: &Bounce2D, confined: &Bounce2D) -> f64 {
    free.restitution_ratio() / confined.restitution_ratio()
}

/// Run one flat-plate slug bounce and return its wall impulse / incident momentum / peak force.
#[must_use]
pub fn run_slug_bounce(cfg: &SlugConfig) -> Bounce2D {
    let dz = cfg.z_max / cfg.nz as f64;
    let dr = cfg.r_max / cfg.nr as f64;
    let mut g = Grid2D::new(cfg.nz, cfg.nr, dz, dr, cfg.gamma);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect; // axis
    g.bc_rhi = if cfg.confined {
        Bc::Reflect // plane-wave confinement: no sideways escape
    } else {
        Bc::Transmissive // gas spreading past the domain edge leaves
    };
    g.bc_zlo = Bc::Reflect; // the plate (edge governed by plate_radius)
    g.bc_zhi = Bc::Transmissive; // rebounding gas leaves the far end
    g.set_plate_radius(Some(cfg.r_plate));

    // Cold slug against the plate (z ∈ [0, L]) moving in at v = 1, in near-vacuum ambient.
    let v = 1.0;
    let rho0 = 1.0;
    let p0 = 1.0 / (cfg.gamma * cfg.mach * cfg.mach);
    let rho_amb = 1.0e-3;
    let p_amb = p0 * 1.0e-3;
    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        if z < cfg.length && r < cfg.r_foot {
            Prim::new(rho0, -v, 0.0, p0)
        } else {
            Prim::new(rho_amb, 0.0, 0.0, p_amb)
        }
    });

    let incident_momentum = g.axial_momentum().abs();

    // Integrate the plate impulse (trapezoid) until the force decays past the cutoff after its peak.
    let mut wall_impulse = 0.0;
    let mut peak = 0.0_f64;
    let mut past_peak = false;
    let mut force_old = g.plate_force();
    let max_steps = 400 * cfg.nz + 50_000;
    let mut steps = 0;
    while steps < max_steps {
        let dt = g.stable_dt();
        g.step(dt);
        let force_new = g.plate_force();
        wall_impulse += 0.5 * dt * (force_old + force_new);
        peak = peak.max(force_new);
        if force_new < 0.999 * peak {
            past_peak = true;
        }
        steps += 1;
        if past_peak && force_new < 1.0e-3 * peak {
            break;
        }
        force_old = force_new;
    }

    Bounce2D {
        wall_impulse,
        incident_momentum,
        peak_force: peak,
        steps,
    }
}
