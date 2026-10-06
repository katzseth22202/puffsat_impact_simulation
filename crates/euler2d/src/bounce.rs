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
use crate::plate::PlateProfile;
use crate::state::Prim;

/// The plate the slug bounces off.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlateShape {
    /// Grid-aligned flat plate at `z = 0` (the verified D2 path; the confined denominator and the
    /// original flat-plate result use it).
    FlatGridAligned,
    /// Immersed shallow-concave dish of depth-to-diameter ratio `d/D` (ADR-0021); `d/D = 0` is a
    /// raised flat plate imposed through the *same* immersed boundary (ADR-0023). Running both the
    /// flat baseline (`d/D = 0`) and the concave plates through the IBM keeps the curvature gain
    /// free of any grid-alignment artifact.
    Dish {
        /// Depth-to-diameter ratio `d/D` (`D = 2·r_plate`).
        d_over_d: f64,
    },
    /// The dish with a wall rising from its rim: a skirt (`flare = 0`) or a flared bell
    /// (ADR-0055). Zero `skirt_height` runs exactly as [`Self::Dish`].
    Cup {
        /// Floor depth-to-diameter ratio `d/D`.
        d_over_d: f64,
        /// Wall height above the rim.
        skirt_height: f64,
        /// Outward flare `dr/dz` of the wall (0 is a straight cylinder).
        flare: f64,
        /// Wall thickness in `r` (two or more cells).
        thickness: f64,
    },
}

/// The immersed profile for a dished or cupped plate, `None` for the grid-aligned flat plate.
/// The floor is raised four cells off the domain bottom so solid always underlies it.
fn immersed_profile(cfg: &SlugConfig) -> Option<PlateProfile> {
    let dz = cfg.z_max / cfg.nz as f64;
    let z0 = 4.0 * dz;
    match cfg.shape {
        PlateShape::FlatGridAligned => None,
        PlateShape::Dish { d_over_d } => Some(PlateProfile::Dish {
            r_plate: cfg.r_plate,
            z0,
            depth: d_over_d * 2.0 * cfg.r_plate,
        }),
        PlateShape::Cup {
            d_over_d,
            skirt_height,
            flare,
            thickness,
        } => Some(PlateProfile::Cup {
            r_plate: cfg.r_plate,
            z0,
            depth: d_over_d * 2.0 * cfg.r_plate,
            skirt_height,
            flare,
            thickness,
        }),
    }
}

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
    /// The plate shape (flat or shallow-concave dish).
    pub shape: PlateShape,
    /// Edge taper width as a fraction of `r_foot` (pulse-shape axis, design §13): the cloud density
    /// ramps linearly from full at `r = (1 − taper_frac)·r_foot` to zero at `r = r_foot`, with the
    /// core density raised by [`taper_mass_factor`] so the total mass equals the top-hat's
    /// (equal-mass contract, ADR-0028). `0.0` is the *exact* top-hat — the identity every
    /// pre-§13 result was produced under.
    pub taper_frac: f64,
    /// Radial divergence (pulse-shape axis, design §13): the cloud initializes with
    /// `v_r = alpha_div·v·r/r_foot` — an expanding puffed cloud at impact. Axially inert (adds no
    /// axial momentum at `t = 0`). `0.0` is the purely-axial identity.
    pub alpha_div: f64,
}

/// Relative density weight of the tapered edge profile at radius `r`: `1` in the core
/// (`r ≤ (1 − taper_frac)·r_foot`), a linear ramp to `0` at `r = r_foot`, `0` outside. At
/// `taper_frac = 0` this is exactly the top-hat indicator (design §13).
#[must_use]
pub fn taper_weight(r: f64, r_foot: f64, taper_frac: f64) -> f64 {
    if r >= r_foot {
        return 0.0;
    }
    if taper_frac <= 0.0 {
        return 1.0;
    }
    let r1 = r_foot * (1.0 - taper_frac);
    if r <= r1 {
        1.0
    } else {
        (r_foot - r) / (r_foot - r1)
    }
}

/// Core-density multiplier restoring the top-hat mass under the edge taper (the equal-mass
/// contract, design §13/ADR-0028): `m_tophat / m_taper` for the linear ramp, from the closed form
/// `∫₀^1 w(r)·2r dr = r₁² + [1/3 − r₁² + (2/3)r₁³]/taper_frac` with `r₁ = 1 − taper_frac`
/// (scale-invariant in `r_foot`). Exactly `1` at `taper_frac = 0`.
#[must_use]
pub fn taper_mass_factor(taper_frac: f64) -> f64 {
    if taper_frac <= 0.0 {
        return 1.0;
    }
    let r1 = 1.0 - taper_frac;
    let band = (1.0 / 3.0 - r1 * r1 + 2.0 / 3.0 * r1 * r1 * r1) / taper_frac;
    1.0 / (r1 * r1 + band)
}

/// Mass-weighted areal-density statistics of the tapered profile, as factors on the top-hat `Σ`
/// (the §13 Σ bookkeeping): `Σ(r) = Σ_tophat · s · w(r)` with `s = taper_mass_factor`, so the
/// factors are profile-only and scale-invariant. All are `1` at zero taper.
#[derive(Debug, Clone, Copy)]
pub struct TaperSigmaStats {
    /// Mass-weighted mean `Σ` factor `s²·∫w²·2r dr` — where the headline `e_eff` is evaluated.
    pub mean: f64,
    /// Core `Σ` factor `s` — the profile's Σ-max bound point.
    pub hi: f64,
    /// `Σ` factor at the radius enclosing 90% of the pulse mass — the representative Σ-low bound
    /// point (the extreme ramp tail below 10% of the mass cannot move the mass-weighted `e_eff`
    /// materially; ADR-0028 records the bound as a named halt condition, not a silent average).
    pub lo90: f64,
}

/// Compute the [`TaperSigmaStats`] of the linear edge-ramp profile at `taper_frac` (design §13).
/// The mean uses the closed form `∫₀^1 w²·2r dr = r₁² + [F(1) − F(r₁)]/taper_frac²` with
/// `F(r) = r² − (4/3)r³ + (1/2)r⁴`; the 90%-mass radius is found by bisection on the cumulative
/// mass (monotone), then mapped through `w(r)`.
#[must_use]
pub fn taper_sigma_stats(taper_frac: f64) -> TaperSigmaStats {
    if taper_frac <= 0.0 {
        return TaperSigmaStats {
            mean: 1.0,
            hi: 1.0,
            lo90: 1.0,
        };
    }
    let s = taper_mass_factor(taper_frac);
    let r1 = 1.0 - taper_frac;
    let f = |r: f64| r * r - 4.0 / 3.0 * r * r * r + 0.5 * r * r * r * r;
    let w2_int = r1 * r1 + (f(1.0) - f(r1)) / (taper_frac * taper_frac);
    // Cumulative mass fraction M(r) = s·∫₀^r w·2r dr (M(1) = 1 by the equal-mass construction).
    let cum = |r: f64| {
        if r <= r1 {
            s * r * r
        } else {
            let g = |x: f64| x * x - 2.0 / 3.0 * x * x * x;
            s * (r1 * r1 + (g(r) - g(r1)) / taper_frac)
        }
    };
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if cum(mid) < 0.9 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let r90 = 0.5 * (lo + hi);
    TaperSigmaStats {
        mean: s * s * w2_int,
        hi: s,
        lo90: s * taper_weight(r90, 1.0, taper_frac),
    }
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
    /// Peak *local* facesheet pressure (max surface-cell pressure over time) — the survivability
    /// concentration that a concave plate focuses above the flat stagnation value (Rung S).
    pub peak_local_pressure: f64,
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

/// Build and initialize the slug-bounce grid from `cfg` — the plate, the boundary conditions, and
/// the cold slug (with the §13 edge taper and radial divergence; both `0` reproduce the exact
/// top-hat every pre-§13 result was produced under). Public so the identity/bookkeeping tests can
/// inspect the initial state directly (ADR-0028).
#[must_use]
pub fn init_slug_grid(cfg: &SlugConfig) -> Grid2D {
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
    g.bc_zhi = Bc::Transmissive; // rebounding gas leaves the far end

    // Configure the plate and find the height the slug starts above. A grid-aligned flat plate sits
    // at z = 0; an immersed dish is raised a few cells off the floor so a solid layer always
    // underlies it, and the slug starts above the whole dish (its rim).
    let z_floor = match cfg.shape {
        PlateShape::FlatGridAligned => {
            g.bc_zlo = Bc::Reflect; // the plate (edge governed by plate_radius)
            g.set_plate_radius(Some(cfg.r_plate));
            0.0
        }
        PlateShape::Dish { .. } | PlateShape::Cup { .. } => {
            g.bc_zlo = Bc::Transmissive; // the immersed surface, not z = 0, is the wall
            let profile = immersed_profile(cfg).expect("a dished or cupped plate");
            g.set_plate_profile(Some(profile));
            profile.z_surface(cfg.r_plate) // the rim
        }
    };

    // Cold slug just above the plate (z ∈ [z_floor, z_floor + L]) moving in at v = 1, in
    // near-vacuum. The tapered density profile carries the equal-mass core multiplier; a tapered
    // cell whose density would fall below the ambient floor is ambient (the extreme edge tail —
    // mass error second-order, and the `f` assembly normalizes by the *measured* p_in anyway).
    let v = 1.0;
    let rho0 = 1.0;
    let p0 = 1.0 / (cfg.gamma * cfg.mach * cfg.mach);
    let rho_amb = 1.0e-3;
    let p_amb = p0 * 1.0e-3;
    let core_rho = rho0 * taper_mass_factor(cfg.taper_frac);
    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        let rho = core_rho * taper_weight(r, cfg.r_foot, cfg.taper_frac);
        if z >= z_floor && z < z_floor + cfg.length && rho > rho_amb {
            let v_r = cfg.alpha_div * v * r / cfg.r_foot;
            Prim::new(rho, -v, v_r, p0)
        } else {
            Prim::new(rho_amb, 0.0, 0.0, p_amb)
        }
    });
    g
}

/// Run one flat-plate slug bounce and return its wall impulse / incident momentum / peak force.
#[must_use]
pub fn run_slug_bounce(cfg: &SlugConfig) -> Bounce2D {
    run_bounce(init_slug_grid(cfg), cfg.nz)
}

/// A spray cloud resting on the plate under the slug (ADR-0055): `depth` deep over the slug's
/// footprint and radial profile, holding `mass_ratio` times the slug's mass per unit area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SprayCloud {
    /// Cloud depth above the plate (same length units as the slug).
    pub depth: f64,
    /// Spray mass per slug mass, the injection ratio `k`.
    pub mass_ratio: f64,
    /// Gap between the plate and the bottom of the cloud, left as ambient gas. `0` puts the cloud
    /// on the plate.
    pub standoff: f64,
}

/// The slug grid with a resting spray cloud between it and the plate. Over the local floor height
/// `z_f(r)` (0 for the flat plate, the dish surface otherwise) the cloud fills
/// `z ∈ [z_f + s, z_f + s + depth)`, `s` the standoff, at density `k·ρ₀·L/depth` and the slug's
/// radial profile, at rest and at the slug's cold pressure; the gap below it is ambient. The slug
/// sits flat on top of the cloud at the rim, `z ∈ [z_rim + s + depth, z_rim + s + depth + L)`.
#[must_use]
pub fn init_spray_grid(cfg: &SlugConfig, spray: SprayCloud) -> Grid2D {
    init_spray_grid_with_ambient(cfg, spray, AMBIENT_DEFAULT)
}

/// The ambient gas's density and pressure as fractions of the slug's: the vacuum stand-in.
pub const AMBIENT_DEFAULT: f64 = 1.0e-3;

/// [`init_spray_grid`] with the ambient at `ambient` times the slug's density and cold pressure.
/// The default is a negligible share of a spray cup's mass, but not of a domain much wider than the
/// slug (one necklace sphere on its plate share, the plug study), where it must be thinner.
#[must_use]
pub fn init_spray_grid_with_ambient(cfg: &SlugConfig, spray: SprayCloud, ambient: f64) -> Grid2D {
    init_spray(cfg, spray, ambient, None)
}

/// A trailing pearl (ADR-0055): a second PuffSat of the lead's architecture (same footprint, length,
/// radial profile, speed and temperature) flying in formation behind it, holding `mass_ratio` times
/// the lead's mass. Its gas meets the rebound off the plate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrailingPearl {
    /// Pearl mass per lead mass.
    pub mass_ratio: f64,
    /// Ambient gap between the lead's tail and the pearl's front (same length units as the slug).
    pub gap: f64,
}

/// [`init_spray_grid`] with a trailing pearl: over the lead, after `gap` of ambient, the pearl
/// fills `z ∈ [top + L + gap, top + 2L + gap)` at `mass_ratio` times the lead's density and cold
/// pressure (so the same temperature and Mach number), moving at the lead's speed.
#[must_use]
pub fn init_spray_grid_with_pearl(
    cfg: &SlugConfig,
    spray: SprayCloud,
    pearl: TrailingPearl,
) -> Grid2D {
    init_spray(cfg, spray, AMBIENT_DEFAULT, Some(pearl))
}

fn init_spray(
    cfg: &SlugConfig,
    spray: SprayCloud,
    ambient: f64,
    pearl: Option<TrailingPearl>,
) -> Grid2D {
    let profile = immersed_profile(cfg);
    let floor = |r: f64| profile.map_or(0.0, |p| p.z_surface(r));
    let z_rim = floor(cfg.r_plate);
    let mut g = init_slug_grid(&SlugConfig {
        length: 0.0,
        ..*cfg
    });
    let dz = cfg.z_max / cfg.nz as f64;
    let dr = cfg.r_max / cfg.nr as f64;
    let v = 1.0;
    let p0 = 1.0 / (cfg.gamma * cfg.mach * cfg.mach);
    let rho_amb = ambient;
    let p_amb = p0 * ambient;
    let core_rho = taper_mass_factor(cfg.taper_frac);
    let spray_rho = spray.mass_ratio * cfg.length / spray.depth;
    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        let weight = core_rho * taper_weight(r, cfg.r_foot, cfg.taper_frac);
        let bottom = floor(r) + spray.standoff;
        let top = z_rim + spray.standoff + spray.depth;
        if z >= bottom && z < bottom + spray.depth && spray_rho * weight > rho_amb {
            Prim::new(spray_rho * weight, 0.0, 0.0, p0)
        } else if z >= top && z < top + cfg.length && weight > rho_amb {
            let v_r = cfg.alpha_div * v * r / cfg.r_foot;
            Prim::new(weight, -v, v_r, p0)
        } else if let Some(pearl) = pearl.filter(|q| {
            let front = top + cfg.length + q.gap;
            z >= front && z < front + cfg.length && q.mass_ratio * weight > rho_amb
        }) {
            let v_r = cfg.alpha_div * v * r / cfg.r_foot;
            Prim::new(pearl.mass_ratio * weight, -v, v_r, p0 * pearl.mass_ratio)
        } else {
            Prim::new(rho_amb, 0.0, 0.0, p_amb)
        }
    });
    g
}

/// Run one spray-plate bounce: the slug arrives through the resting cloud onto the flat plate.
/// `incident_momentum` is the slug's, since the cloud starts at rest.
#[must_use]
pub fn run_spray_bounce(cfg: &SlugConfig, spray: SprayCloud) -> Bounce2D {
    run_bounce(init_spray_grid(cfg, spray), cfg.nz)
}

/// [`run_spray_bounce`] with a trailing pearl. `incident_momentum` is the lead's and the pearl's.
#[must_use]
pub fn run_spray_bounce_with_pearl(
    cfg: &SlugConfig,
    spray: SprayCloud,
    pearl: TrailingPearl,
) -> Bounce2D {
    run_bounce(init_spray_grid_with_pearl(cfg, spray, pearl), cfg.nz)
}

/// [`run_spray_bounce`] with the ambient at `ambient` times the slug's density and pressure.
#[must_use]
pub fn run_spray_bounce_with_ambient(
    cfg: &SlugConfig,
    spray: SprayCloud,
    ambient: f64,
) -> Bounce2D {
    run_bounce(init_spray_grid_with_ambient(cfg, spray, ambient), cfg.nz)
}

/// The pressure on the cup wall's inner face through a bounce (ADR-0055): `pressure[i][j]` is row
/// `z[j]`'s face pressure at `time[i]`, sampled at the start and after every step.
#[derive(Debug, Clone, PartialEq)]
pub struct SkirtHistory {
    /// Row heights (cell centres from the rim to the wall's lip).
    pub z: Vec<f64>,
    /// Sample times.
    pub time: Vec<f64>,
    /// Face pressure per sample, one entry per row.
    pub pressure: Vec<Vec<f64>>,
    /// The dish floor's risers (the bowl's band): heights, ring radii, and pressure per sample.
    pub band_z: Vec<f64>,
    pub band_r: Vec<f64>,
    pub band_pressure: Vec<Vec<f64>>,
}

/// [`run_spray_bounce`] that also records the skirt's face pressure every step. Recording only
/// reads the grid, so the bounce is identical to the plain run.
#[must_use]
pub fn run_spray_bounce_with_skirt_history(
    cfg: &SlugConfig,
    spray: SprayCloud,
) -> (Bounce2D, SkirtHistory) {
    let g = init_spray_grid(cfg, spray);
    let z: Vec<f64> = g.skirt_face_pressure().iter().map(|&(z, _)| z).collect();
    let risers = g.floor_riser_pressure();
    let mut history = SkirtHistory {
        z,
        time: Vec::new(),
        pressure: Vec::new(),
        band_z: risers.iter().map(|&(z, _, _)| z).collect(),
        band_r: risers.iter().map(|&(_, r, _)| r).collect(),
        band_pressure: Vec::new(),
    };
    let bounce = run_bounce_observed(g, cfg.nz, |g, t| {
        history.time.push(t);
        history
            .pressure
            .push(g.skirt_face_pressure().iter().map(|&(_, p)| p).collect());
        history.band_pressure.push(
            g.floor_riser_pressure()
                .iter()
                .map(|&(_, _, p)| p)
                .collect(),
        );
    });
    (bounce, history)
}

/// Integrate the plate impulse on an initialized grid until the force has decayed (ADR-0001).
fn run_bounce(g: Grid2D, nz: usize) -> Bounce2D {
    run_bounce_observed(g, nz, |_, _| {})
}

/// [`run_bounce`] calling `observe(grid, t)` at the start and after every step.
fn run_bounce_observed(
    mut g: Grid2D,
    nz: usize,
    mut observe: impl FnMut(&Grid2D, f64),
) -> Bounce2D {
    let incident_momentum = g.axial_momentum().abs();
    let mut t = 0.0;
    observe(&g, t);

    // Integrate the plate impulse (trapezoid) until the force has stayed below the 10⁻³-of-peak
    // cutoff for a sustained window — not merely crossed it once. A concave plate refocuses the
    // rebound into a secondary peak (§7); requiring a sustained quiet tail keeps the integration
    // from truncating in the trough before that secondary impulse arrives. The window resets
    // whenever the force climbs back above the cutoff (i.e. a secondary peak), so it captures it.
    let mut wall_impulse = 0.0;
    let mut peak = 0.0_f64;
    let mut peak_local = 0.0_f64;
    let mut past_peak = false;
    let mut steps_below = 0_usize;
    let window = 40;
    let mut force_old = g.plate_force();
    let max_steps = 400 * nz + 50_000;
    let mut steps = 0;
    while steps < max_steps {
        let dt = g.stable_dt();
        g.step(dt);
        t += dt;
        observe(&g, t);
        let force_new = g.plate_force();
        wall_impulse += 0.5 * dt * (force_old + force_new);
        peak = peak.max(force_new);
        peak_local = peak_local.max(g.max_plate_pressure());
        if force_new < 0.999 * peak {
            past_peak = true;
        }
        if force_new < 1.0e-3 * peak {
            steps_below += 1;
        } else {
            steps_below = 0;
        }
        steps += 1;
        if past_peak && steps_below >= window {
            break;
        }
        force_old = force_new;
    }

    Bounce2D {
        wall_impulse,
        incident_momentum,
        peak_force: peak,
        peak_local_pressure: peak_local,
        steps,
    }
}

#[cfg(test)]
mod tests {
    use super::{
        PlateShape, SlugConfig, init_slug_grid, run_slug_bounce, taper_mass_factor,
        taper_sigma_stats, taper_weight,
    };

    /// A small dished free-bounce config the shape tests perturb (taper/divergence default 0).
    fn base_cfg() -> SlugConfig {
        SlugConfig {
            gamma: 1.4,
            mach: 5.0,
            r_foot: 1.0,
            length: 0.6,
            r_plate: 2.0,
            r_max: 2.8,
            z_max: 3.1,
            nr: 56,
            nz: 40,
            confined: false,
            shape: PlateShape::Dish { d_over_d: 0.10 },
            taper_frac: 0.0,
            alpha_div: 0.0,
        }
    }

    /// Identity regression (ADR-0028): `taper_frac = 0`, `alpha_div = 0` initializes the *exact*
    /// top-hat — unit density and purely axial `v = 1` inside `{z-window, r < r_foot}`, ambient
    /// outside, zero radial momentum — so every pre-§13 result is reproduced unchanged.
    // Exact float `==` is intentional: the identity claim is bit-exactness, not closeness.
    #[allow(clippy::float_cmp)]
    #[test]
    fn zero_taper_zero_divergence_is_the_exact_top_hat() {
        let cfg = base_cfg();
        let g = init_slug_grid(&cfg);
        let dz = cfg.z_max / cfg.nz as f64;
        let dr = cfg.r_max / cfg.nr as f64;
        // The dish raises the slug floor to z0 + depth (see init_slug_grid).
        let z_floor = 4.0 * dz + 0.10 * 2.0 * cfg.r_plate;
        let mut cloud_cells = 0;
        for iz in 0..cfg.nz {
            for ir in 0..cfg.nr {
                let z = (iz as f64 + 0.5) * dz;
                let r = (ir as f64 + 0.5) * dr;
                let p = g.prim(iz, ir);
                if z >= z_floor && z < z_floor + cfg.length && r < cfg.r_foot {
                    assert_eq!(p.rho, 1.0, "cloud density must be exactly 1 at ({iz},{ir})");
                    assert_eq!(p.uz, -1.0);
                    assert_eq!(p.ur, 0.0);
                    cloud_cells += 1;
                } else {
                    assert_eq!(p.rho, 1.0e-3, "ambient at ({iz},{ir})");
                    assert_eq!(p.uz, 0.0);
                }
            }
        }
        assert!(cloud_cells > 0, "no cloud cells initialized");
        assert_eq!(g.radial_momentum(), 0.0);
    }

    /// The taper mass factor's closed form matches a numeric integral of the ramp profile, and is
    /// exactly 1 at zero taper (the identity the §13 axis is one-sided from).
    #[allow(clippy::float_cmp)] // the zero-taper identity is exact by construction
    #[test]
    fn taper_mass_factor_matches_numeric_integral() {
        assert_eq!(taper_mass_factor(0.0), 1.0);
        for taper in [0.1, 0.2, 0.3] {
            // ∫₀^1 w(r)·2r dr by midpoint rule; the factor is its reciprocal.
            let n = 200_000_usize;
            let mass: f64 = (0..n)
                .map(|i| {
                    let r = (i as f64 + 0.5) / n as f64;
                    taper_weight(r, 1.0, taper) * 2.0 * r / n as f64
                })
                .sum();
            let factor = taper_mass_factor(taper);
            assert!(
                (factor - 1.0 / mass).abs() < 1e-6,
                "taper={taper}: closed form {factor} vs numeric {}",
                1.0 / mass
            );
            assert!(factor > 1.0, "equal mass needs a raised core: {factor}");
        }
    }

    /// Equal-mass bookkeeping (design §13): the tapered cloud initializes the same incident axial
    /// momentum (= mass, since `u_z = −1` uniformly) as the top-hat, to discretization tolerance —
    /// and the discrepancy shrinks under grid refinement.
    #[test]
    fn tapered_cloud_conserves_incident_momentum() {
        let p_in = |taper: f64, scale: usize| {
            let cfg = SlugConfig {
                taper_frac: taper,
                nr: 56 * scale,
                nz: 40 * scale,
                ..base_cfg()
            };
            init_slug_grid(&cfg).axial_momentum().abs()
        };
        let rel = |scale: usize| (p_in(0.3, scale) - p_in(0.0, scale)).abs() / p_in(0.0, scale);
        let coarse = rel(1);
        let fine = rel(2);
        assert!(coarse < 0.03, "coarse-grid mass error {coarse} > 3%");
        assert!(
            fine < coarse,
            "mass error must shrink under refinement: {coarse} -> {fine}"
        );
    }

    /// Radial divergence is axially inert (design §13): `alpha_div` leaves the incident axial
    /// momentum bit-identical, and the initialized radial momentum matches the analytic
    /// `Σ ρ·(α·v·r/r_foot)·r dr dz` over the cloud cells.
    #[test]
    fn divergence_is_axially_inert_and_matches_analytic_radial_momentum() {
        let cfg0 = base_cfg();
        let cfg_div = SlugConfig {
            alpha_div: 0.1,
            ..base_cfg()
        };
        let g0 = init_slug_grid(&cfg0);
        let g = init_slug_grid(&cfg_div);
        assert_eq!(
            g.axial_momentum().to_bits(),
            g0.axial_momentum().to_bits(),
            "divergence must not touch the axial momentum"
        );
        // Discrete analytic expectation over the same cell centers (2π dropped, as in the kernel).
        let dz = cfg_div.z_max / cfg_div.nz as f64;
        let dr = cfg_div.r_max / cfg_div.nr as f64;
        let z_floor = 4.0 * dz + 0.10 * 2.0 * cfg_div.r_plate;
        let mut expected = 0.0;
        for iz in 0..cfg_div.nz {
            for ir in 0..cfg_div.nr {
                let z = (iz as f64 + 0.5) * dz;
                let r = (ir as f64 + 0.5) * dr;
                if z >= z_floor && z < z_floor + cfg_div.length && r < cfg_div.r_foot {
                    expected += 1.0 * (0.1 * r / cfg_div.r_foot) * r * dr * dz;
                }
            }
        }
        let measured = g.radial_momentum();
        assert!(
            (measured - expected).abs() <= 1e-12 * expected,
            "radial momentum {measured} vs analytic {expected}"
        );
    }

    /// The taper Σ statistics' closed forms match numeric integration of the ramp profile, are all
    /// exactly 1 at zero taper, and are ordered `1 < mean < hi` with `lo90 < hi` (a raised core
    /// and a mass-weighted mean above the top-hat; `lo90` only drops below 1 once the taper is
    /// wide enough that the 90%-mass radius sits deep in the ramp — checked at 0.3).
    #[allow(clippy::float_cmp)] // the zero-taper identity is exact by construction
    #[test]
    fn taper_sigma_stats_match_numeric_integration() {
        let id = taper_sigma_stats(0.0);
        assert_eq!((id.mean, id.hi, id.lo90), (1.0, 1.0, 1.0));
        for taper in [0.1, 0.3] {
            let st = taper_sigma_stats(taper);
            let s = taper_mass_factor(taper);
            // Numeric mass-weighted mean Σ factor: ∫ (s·w)²·2r dr (the s·w mass weight integrates
            // to 1 by the equal-mass construction).
            let n = 200_000_usize;
            let (mut mean, mut cum, mut r90) = (0.0, 0.0, f64::NAN);
            for i in 0..n {
                let r = (i as f64 + 0.5) / n as f64;
                let sw = s * taper_weight(r, 1.0, taper);
                mean += sw * sw * 2.0 * r / n as f64;
                cum += sw * 2.0 * r / n as f64;
                if cum >= 0.9 && r90.is_nan() {
                    r90 = r;
                }
            }
            assert!(
                (st.mean - mean).abs() < 1e-4,
                "mean {} vs numeric {mean}",
                st.mean
            );
            assert_eq!(st.hi, s);
            let lo_num = s * taper_weight(r90, 1.0, taper);
            assert!(
                (st.lo90 - lo_num).abs() < 1e-3,
                "lo90 {} vs numeric {lo_num}",
                st.lo90
            );
            assert!(1.0 < st.mean && st.mean < st.hi && st.lo90 < st.hi);
        }
        assert!(
            taper_sigma_stats(0.3).lo90 < 1.0,
            "a 30% taper's 90%-mass radius sits deep in the diluted ramp"
        );
    }

    /// Taper plausibility (design §13 exit criterion): concentrating the mass toward the axis must
    /// not *hurt* capture — the dished free restitution with a 30% edge taper stays at or above the
    /// top-hat value (small tolerance for grid noise), and within a sane band of it.
    #[test]
    fn taper_direction_is_plausible_on_the_dished_bounce() {
        let ratio = |taper: f64| {
            run_slug_bounce(&SlugConfig {
                taper_frac: taper,
                ..base_cfg()
            })
            .restitution_ratio()
        };
        let top_hat = ratio(0.0);
        let tapered = ratio(0.3);
        assert!(
            tapered > top_hat - 0.02,
            "taper should not lose capture: top-hat {top_hat:.4}, tapered {tapered:.4}"
        );
        assert!(
            (tapered - top_hat).abs() < 0.15,
            "taper is a *slight* shape change: top-hat {top_hat:.4}, tapered {tapered:.4}"
        );
    }

    /// Very-strong-shock dished bounce stays bounded and terminates (M=40 rim-corner regression).
    ///
    /// Before the dish rim's vertical side face entered the immersed-boundary geometry, solid
    /// cells near the rim were mirrored across the faraway *top* surface; the resulting spurious
    /// radial fluxes at the rim corner acted as an energy source — bounded at M ≲ 20 but
    /// self-exciting at M = 40 (`J/p_in` → 12 at production resolution, → 10⁵² at this coarse
    /// one). The physical restitution here is ~1.5; steps must stay well below the driver's cap
    /// (the pathology also manifests as a never-decaying phantom force that runs the integration
    /// to `max_steps`).
    #[test]
    fn strong_shock_dished_bounce_terminates_with_physical_restitution() {
        let b = run_slug_bounce(&SlugConfig {
            gamma: 1.4,
            mach: 40.0,
            r_foot: 1.0,
            length: 0.6,
            r_plate: 2.0,
            r_max: 2.8,
            z_max: 3.1,
            nr: 56,
            nz: 40,
            confined: false,
            shape: PlateShape::Dish { d_over_d: 0.10 },
            taper_frac: 0.0,
            alpha_div: 0.0,
        });
        let ratio = b.restitution_ratio();
        assert!(
            ratio > 1.0 && ratio < 2.5,
            "unphysical restitution {ratio} (rim-corner blow-up?)"
        );
        assert!(
            b.steps < 400 * 40 + 50_000,
            "force never decayed: ran to the step cap ({} steps)",
            b.steps
        );
    }
}
