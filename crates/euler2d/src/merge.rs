//! The merged fireball's free expansion — the initial-value problem behind N1 and N7.
//!
//! # What this simulates, and why it is the right question
//!
//! `eq:reflection_baseline` asks what a nozzle that only *reflected* would return. Its derivation
//! spreads `cos(theta)` uniformly over `[-1, 1]` — an **isotropic** thermal remainder, `alpha = 1/3`
//! — while `sec:minimum_nozzle` says in the same breath that the stagnation region "squirts
//! radially into the field". N1 exists because those cannot both be true and the baseline goes as
//! `sqrt(2 alpha/pi)`, so the disagreement is worth a factor of two.
//!
//! The quantity the baseline needs is the velocity distribution the plume would reach **with no
//! nozzle at all**: the merged fireball expanding freely until its pressure has done its work. That
//! is a clean initial-value problem — no inflow boundary, no bore wall, no field, no radiation — and
//! it is what this module sets up.
//!
//! # The analytic bracket this has to land inside
//!
//! Two limits bound it and they are three times apart in the baseline, which is why an argument
//! could not settle it:
//!
//! - **Shape-following.** A body expanding self-similarly keeps its own aspect ratio, and for the
//!   flown 23 m × 3.02 m column that gives `alpha = 0.91` — strongly *prolate*, baseline 0.76.
//! - **Gradient-following.** Free expansion into vacuum is driven by pressure gradients, and a long
//!   thin body's steepest gradient is across its **short** axis, so it accelerates fastest
//!   radially — *oblate*, `alpha < 1/3`, baseline below 0.46 and falling.
//!
//! The paper assumes 1/3, which sits between them. Shape-prolate does not imply velocity-prolate,
//! and that is exactly the step no closed form takes for you.
//!
//! # What is initialised, and what that assumes
//!
//! The **merged** state, not the collision that produces it: 238 kg filling the column uniformly at
//! the bag density, carrying the pulse energy split into bulk drift (`f_d`) and internal energy
//! (the rest). That isolates N1's question — *given the fireball the paper describes, what does free
//! expansion do to it?* — from the separate question of whether the snowplow actually produces a
//! uniform fireball, which needs the rod resolved and is the follow-on run.

use crate::kernel::{Bc, Grid2D};
use crate::state::Prim;

/// Geometry and energy budget of one merged pulse.
#[derive(Debug, Clone, Copy)]
pub struct MergeConfig {
    /// Effective γ for the calibrated 2D track (ADR-0008).
    pub gamma: f64,
    /// Projectile mass [kg].
    pub m_projectile: f64,
    /// Slug mass [kg].
    pub m_slug: f64,
    /// Closing speed [m/s].
    pub closing_speed: f64,
    /// Column length [m] the merged mass occupies.
    pub column_length: f64,
    /// Column radius [m].
    pub column_radius: f64,
    /// Domain half-extents: the fireball must stay inside for the moments to be volume-complete.
    pub z_max: f64,
    pub r_max: f64,
    pub nz: usize,
    pub nr: usize,
    /// Ambient density as a fraction of the fireball's — a stand-in for vacuum, which a Godunov
    /// scheme cannot take literally. Checked for irrelevance by `test_ambient_does_not_set_alpha`.
    pub ambient_fraction: f64,
    /// Fill a **sphere** of `column_radius` instead of the cylinder — the control case.
    ///
    /// The paper's bag was originally a 5.4 m sphere and became a 23 m column, and the reason was
    /// the launch fairing rather than the plume: "length is the cheap dimension for a rocket to
    /// carry and diameter is the expensive one". Running the same mass, density and energy in a
    /// sphere isolates whether the anisotropy is caused by that packing decision.
    pub spherical: bool,
}

impl MergeConfig {
    /// The flown Jupiter-only pulse: 25 kg into 213 kg at `k = 8.5`, in the 23 m × 3.02 m column.
    #[must_use]
    pub fn flown(closing_speed: f64, nz: usize, nr: usize) -> Self {
        Self {
            gamma: 5.0 / 3.0,
            m_projectile: 25.0,
            m_slug: 213.0,
            closing_speed,
            column_length: 23.0,
            column_radius: 3.02,
            z_max: 140.0,
            r_max: 40.0,
            nz,
            nr,
            ambient_fraction: 1e-6,
            spherical: false,
        }
    }

    /// The control: the same 238 kg at the same density in the **5.4 m spherical bag**
    /// the paper started with, before the launch envelope stretched it into a column.
    #[must_use]
    pub fn spherical_bag(closing_speed: f64, nz: usize, nr: usize) -> Self {
        let mut cfg = Self::flown(closing_speed, nz, nr);
        // Radius chosen so the sphere holds the column's volume: same mass, same
        // density, same energy. Only the shape differs.
        let volume =
            std::f64::consts::PI * cfg.column_radius * cfg.column_radius * cfg.column_length;
        cfg.column_radius = (3.0 * volume / (4.0 * std::f64::consts::PI)).cbrt();
        cfg.spherical = true;
        cfg
    }

    /// Merged mass `m_p + m_s` [kg].
    #[must_use]
    pub fn merged_mass(&self) -> f64 {
        self.m_projectile + self.m_slug
    }

    /// Slug ratio `k = m_s/m_p`.
    #[must_use]
    pub fn slug_ratio(&self) -> f64 {
        self.m_slug / self.m_projectile
    }

    /// Pulse energy `½ mu w²` on the reduced mass — 62.9 GJ at the flown 75 km/s.
    #[must_use]
    pub fn pulse_energy(&self) -> f64 {
        let mu = self.m_projectile * self.m_slug / self.merged_mass();
        0.5 * mu * self.closing_speed * self.closing_speed
    }

    /// `f_d = 1/(1+k)` — the clean-merge drift share the paper substitutes into the baseline.
    ///
    /// N7 asks whether the realised value matches; here it is an *input*, because this run starts
    /// from an already-merged uniform body. A run that resolves the snowplow would produce it.
    #[must_use]
    pub fn drift_fraction(&self) -> f64 {
        1.0 / (1.0 + self.slug_ratio())
    }

    /// Centre-of-mass speed of the merged body, `w/(1+k)`.
    #[must_use]
    pub fn drift_speed(&self) -> f64 {
        self.closing_speed * self.drift_fraction()
    }

    /// Volume the merged mass occupies [m³] — the cylinder, or the sphere in the control case.
    #[must_use]
    pub fn volume(&self) -> f64 {
        if self.spherical {
            4.0 / 3.0 * std::f64::consts::PI * self.column_radius.powi(3)
        } else {
            std::f64::consts::PI * self.column_radius * self.column_radius * self.column_length
        }
    }

    /// Uniform density of the merged body [kg/m³] — 0.361 at the flown geometry.
    ///
    /// Not the paper's 0.323: `BAG_RHO` is the 213 kg slug over 660 m³, before the projectile is
    /// added. The *merged* 238 kg in the same volume is denser by `(1+k)/k`.
    #[must_use]
    pub fn column_density(&self) -> f64 {
        self.merged_mass() / self.volume()
    }

    /// Specific internal energy [J/kg]: the pulse energy less the bulk drift, over the merged mass.
    #[must_use]
    pub fn specific_internal(&self) -> f64 {
        let internal = self.pulse_energy() * (1.0 - self.drift_fraction());
        internal / self.merged_mass()
    }

    /// Initial pressure of the column, `(γ−1) ρ e_int`.
    #[must_use]
    pub fn column_pressure(&self) -> f64 {
        (self.gamma - 1.0) * self.column_density() * self.specific_internal()
    }
}

/// Build the initial grid: a uniform drifting column in a near-vacuum, free on every boundary.
///
/// The column is centred axially so it can expand both ways without meeting a boundary first, and
/// the axis carries the reflecting condition axisymmetry requires.
#[must_use]
pub fn init_merge_grid(cfg: &MergeConfig) -> Grid2D {
    let dz = cfg.z_max / (cfg.nz as f64);
    let dr = cfg.r_max / (cfg.nr as f64);
    let mut g = Grid2D::new(cfg.nz, cfg.nr, dz, dr, cfg.gamma);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Transmissive;
    g.bc_zlo = Bc::Transmissive;
    g.bc_zhi = Bc::Transmissive;

    let rho = cfg.column_density();
    let p = cfg.column_pressure();
    let uz = cfg.drift_speed();
    let z_lo = 0.5 * (cfg.z_max - cfg.column_length);
    let z_hi = z_lo + cfg.column_length;
    let amb_rho = rho * cfg.ambient_fraction;
    let amb_p = p * cfg.ambient_fraction;

    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        let inside = if cfg.spherical {
            let dz_c = z - 0.5 * cfg.z_max;
            (dz_c * dz_c + r * r).sqrt() <= cfg.column_radius
        } else {
            z >= z_lo && z <= z_hi && r <= cfg.column_radius
        };
        if inside {
            Prim::new(rho, uz, 0.0, p)
        } else {
            Prim::new(amb_rho, 0.0, 0.0, amb_p)
        }
    });
    g
}

/// Sound-crossing time of the column's **short** axis [s] — the timescale the expansion runs on.
///
/// The radius, not the length: free expansion is set by the steepest gradient, and for a long thin
/// column that is across the radius. Running for a few of these is what converts internal energy to
/// directed motion, which is the state `eq:reflection_baseline` is about.
#[must_use]
pub fn radial_sound_crossing(cfg: &MergeConfig) -> f64 {
    let c = (cfg.gamma * cfg.column_pressure() / cfg.column_density()).sqrt();
    cfg.column_radius / c
}
