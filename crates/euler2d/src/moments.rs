//! Velocity-distribution moments of a fireball — the diagnostic N1 and N7 are asked in terms of.
//!
//! The companion repo's `docs/nozzle_asks_for_impact_sim.md` asks two things of the merged plume
//! that the existing aggregate diagnostics (`plate_force`, `axial_momentum`) cannot answer:
//!
//! - **N1** — `alpha = <v_z²>/<v²>`, mass-weighted over the **drift-subtracted** remainder.
//!   `eq:reflection_baseline` assumes the remainder is isotropic (`alpha = 1/3`), while
//!   `sec:minimum_nozzle` says in the same breath that the stagnation region "squirts radially into
//!   the field". Those cannot both be right, and the baseline goes as `sqrt(2 alpha/pi)`.
//! - **N7** — `f_d`, the share of pulse energy carried as bulk drift, against the clean-merge
//!   `1/(1+k)` the paper substitutes. A snowplow deposits momentum differently from a uniform
//!   merge, so the realised value is not obviously the formula's.
//!
//! # What "thermal remainder" means for a fluid code, and why it needs the free expansion
//!
//! In a Euler code the resolved velocity field is *bulk* motion; the random thermal motion lives in
//! the internal energy as a scalar pressure and is isotropic **by construction**. So measuring
//! `alpha` on a freshly merged fireball, while most of its energy is still internal, would report
//! the anisotropy of the bulk flow and miss the part `eq:reflection_baseline` is actually about.
//!
//! The quantity the baseline needs is the velocity distribution the plume ends up with once its
//! pressure has done its work — i.e. **after free expansion into vacuum**, where internal energy has
//! converted to directed kinetic energy. That is the no-nozzle limit the reflection baseline assumes,
//! and it is what these moments should be evaluated on. Evaluating them early is a diagnostic of the
//! run, not an answer to N1.
//!
//! # Conventions
//!
//! Mass-weighted throughout, with the axisymmetric cell volume `2 pi r dr dz` when the grid is
//! cylindrical (the common `2 pi` cancels out of every ratio here, but is kept for legibility).
//!
//! In axisymmetric flow without swirl, a cell's velocity is `(u_r cos phi, u_r sin phi, u_z)` in 3D,
//! so `|v|² = u_z² + u_r²` and the two transverse components together carry `u_r²`. A Hubble-like
//! isotropic expansion therefore returns `alpha = 1/3` exactly, which is this module's acceptance
//! test rather than an approximation of one.

use crate::kernel::Grid2D;

/// Mass-weighted velocity moments of the fluid on a grid.
///
/// Every field is in the grid's own frame; `alpha` and `f_d` are the two the asks want.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Moments {
    /// Total mass `Σ ρ dV` (per radian × 2π when axisymmetric).
    pub mass: f64,
    /// Centre-of-mass axial velocity — the bulk drift the merge produces.
    pub v_cm_z: f64,
    /// Mass-weighted mean **radial expansion speed** — *not* a drift.
    ///
    /// In axisymmetric coordinates `u_r >= 0` everywhere for an expanding cloud, so this sum does
    /// not cancel. The true 3D transverse momentum is nevertheless **zero by symmetry**: the
    /// outward motion at azimuth `phi` is cancelled by that at `phi + pi`. Folding this into the
    /// drift kinetic energy was a real bug here — it made a freely expanding fireball look like it
    /// was carrying six times its own bulk momentum.
    pub v_expansion_r: f64,
    /// Total kinetic energy in the grid frame, `Σ ½ m |u|²`.
    pub kinetic_total: f64,
    /// Kinetic energy of the bulk drift alone, `½ M |v_cm|²`.
    pub kinetic_drift: f64,
    /// Kinetic energy about the centre of mass, `kinetic_total − kinetic_drift`.
    pub kinetic_thermal: f64,
    /// Internal energy `Σ p/(γ−1) dV` — energy not yet converted to motion.
    ///
    /// **The gate on whether `alpha` means anything yet.** While this is a large share of the
    /// total, the expansion is unfinished and `alpha` describes the bulk flow rather than the
    /// distribution `eq:reflection_baseline` is about.
    pub internal: f64,
    /// `f_d` — drift share of the *kinetic* energy, `kinetic_drift/kinetic_total`.
    pub drift_fraction_kinetic: f64,
    /// `f_d` — drift share of total (kinetic + internal) energy. This is the paper's definition,
    /// and it only equals the kinetic one once the expansion has finished.
    pub drift_fraction_total: f64,
    /// `alpha = <v_z²>/<v²>` about the centre of mass — the N1 quantity, drift removed.
    pub alpha: f64,
    /// The same ratio without removing the drift, reported beside it because the ask is ambiguous
    /// on which it wants and they differ by exactly the `f_d` that N7 is about.
    pub alpha_lab: f64,
}

impl Moments {
    /// Share of the total energy still sitting as pressure rather than motion.
    ///
    /// `alpha` is only the answer to N1 once this is small; see the module docs.
    #[must_use]
    pub fn unconverted_fraction(&self) -> f64 {
        let total = self.kinetic_total + self.internal;
        if total <= 0.0 {
            return 0.0;
        }
        self.internal / total
    }
}

/// Mass-weighted velocity moments over every fluid cell of a grid.
///
/// Solid cells under an immersed plate are skipped: their contents are ghost fill, not fluid.
// The axial/radial second-moment accumulators are intrinsically similarly named — they are the
// two components of one velocity moment — so the lint is silenced for this routine, matching
// how `kernel::sweep` handles the same pairing.
#[allow(clippy::similar_names)]
#[must_use]
pub fn moments(grid: &Grid2D) -> Moments {
    let gamma = grid.gamma();
    let mut mass = 0.0;
    let mut mom_z = 0.0;
    let mut mom_r = 0.0;
    let mut ke = 0.0;
    let mut internal = 0.0;
    // Second pass needs the centre of mass, so accumulate the raw sums first.
    let mut sum_uz2 = 0.0;
    let mut sum_ur2 = 0.0;
    let mut sum_uz = 0.0;

    for iz in 0..grid.nz() {
        for ir in 0..grid.nr() {
            if grid.is_solid(iz, ir) {
                continue;
            }
            let w = grid.prim(iz, ir);
            let dv = grid.cell_volume(ir);
            let m = w.rho * dv;
            mass += m;
            mom_z += m * w.uz;
            mom_r += m * w.ur;
            ke += 0.5 * m * (w.uz * w.uz + w.ur * w.ur);
            internal += w.p / (gamma - 1.0) * dv;
            sum_uz2 += m * w.uz * w.uz;
            sum_ur2 += m * w.ur * w.ur;
            sum_uz += m * w.uz;
        }
    }

    if mass <= 0.0 {
        return Moments {
            mass: 0.0,
            v_cm_z: 0.0,
            v_expansion_r: 0.0,
            kinetic_total: 0.0,
            kinetic_drift: 0.0,
            kinetic_thermal: 0.0,
            internal: 0.0,
            drift_fraction_kinetic: 0.0,
            drift_fraction_total: 0.0,
            alpha: 0.0,
            alpha_lab: 0.0,
        };
    }

    let v_cm_z = mom_z / mass;
    let v_expansion_r = mom_r / mass;
    // Axial only: axisymmetry leaves no net transverse momentum to carry drift energy.
    let kinetic_drift = 0.5 * mass * v_cm_z * v_cm_z;
    let kinetic_thermal = (ke - kinetic_drift).max(0.0);

    // <w_z²> about the centre of mass: Σ m (u_z − v_cm)² = Σ m u_z² − 2 v_cm Σ m u_z + M v_cm².
    let sum_wz2 = (sum_uz2 - 2.0 * v_cm_z * sum_uz + mass * v_cm_z * v_cm_z).max(0.0);
    // The transverse pair carries u_r² together; the drift is axial, so it does not enter here.
    let sum_wr2 = sum_ur2;

    let denom_cm = sum_wz2 + sum_wr2;
    let denom_lab = sum_uz2 + sum_ur2;
    let total_energy = ke + internal;

    Moments {
        mass,
        v_cm_z,
        v_expansion_r,
        kinetic_total: ke,
        kinetic_drift,
        kinetic_thermal,
        internal,
        drift_fraction_kinetic: if ke > 0.0 { kinetic_drift / ke } else { 0.0 },
        drift_fraction_total: if total_energy > 0.0 {
            kinetic_drift / total_energy
        } else {
            0.0
        },
        alpha: if denom_cm > 0.0 {
            sum_wz2 / denom_cm
        } else {
            0.0
        },
        alpha_lab: if denom_lab > 0.0 {
            sum_uz2 / denom_lab
        } else {
            0.0
        },
    }
}

/// Histogram of `v_z/|v|` about the centre of mass, mass-weighted into `bins` over `[−1, 1]`.
///
/// N1 asks for this "if the distribution is far from a spheroid", where a single second moment
/// stops being a fair summary. Bins are returned normalised to sum to 1; an isotropic distribution
/// is flat, because `cos(theta)` is uniform on `[−1, 1]` for an isotropic velocity field.
///
/// # Panics
/// If `bins` is zero.
#[must_use]
pub fn cosine_histogram(grid: &Grid2D, bins: usize) -> Vec<f64> {
    assert!(bins > 0, "need at least one bin");
    let m = moments(grid);
    let mut hist = vec![0.0; bins];
    let mut total = 0.0;
    for iz in 0..grid.nz() {
        for ir in 0..grid.nr() {
            if grid.is_solid(iz, ir) {
                continue;
            }
            let w = grid.prim(iz, ir);
            let wz = w.uz - m.v_cm_z;
            let speed = (wz * wz + w.ur * w.ur).sqrt();
            if speed <= 0.0 {
                continue;
            }
            let mass = w.rho * grid.cell_volume(ir);
            // cos(theta) = w_z/|w|, mapped from [−1, 1] onto [0, bins).
            let frac = (wz / speed + 1.0) * 0.5 * (bins as f64);
            // SAFE: `frac` is finite and clamped into [0, bins−1] before the cast.
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let idx = (frac.floor().clamp(0.0, (bins - 1) as f64)) as usize;
            hist[idx] += mass;
            total += mass;
        }
    }
    if total > 0.0 {
        for h in &mut hist {
            *h /= total;
        }
    }
    hist
}
