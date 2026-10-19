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
fn run_bounce_observed(g: Grid2D, nz: usize, observe: impl FnMut(&Grid2D, f64)) -> Bounce2D {
    run_bounce_observed_final(g, nz, observe).0
}

/// [`run_bounce_observed`], also handing back the grid as it stood when the integration stopped.
fn run_bounce_observed_final(
    mut g: Grid2D,
    nz: usize,
    mut observe: impl FnMut(&Grid2D, f64),
) -> (Bounce2D, Grid2D) {
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

    (
        Bounce2D {
            wall_impulse,
            incident_momentum,
            peak_force: peak,
            peak_local_pressure: peak_local,
            steps,
        },
        g,
    )
}

/// A point-source fireball for the X-ray collision study's charge (Result 10): a projectile has
/// threaded the bore and stopped in a polyethylene charge astern of the plate. The merged mass is
/// started as a uniform hot sphere at rest in its own frame, drifting **away** from the plate.
///
/// Units: the fireball's mass is 1 (with the common `2π` dropped, as in the impulse) and its
/// internal energy is 1/2, so its eventual rms expansion speed is about 1 and impulses read
/// directly in units of `M u_rms`. Lengths are physical (metres).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointFireball {
    /// Height of the sphere's centre above the plate's vertex (the focus, for a paraboloid).
    pub height: f64,
    /// Initial sphere radius (several cells; the far-field answer should not depend on it).
    pub radius: f64,
    /// Drift speed away from the plate, in units of `u_rms`: `w/u_rms` with `w = v/(1+k)`.
    pub drift: f64,
    /// Ambient density as a fraction of the sphere's (a vacuum stand-in; keep its mass ≪ 1).
    pub ambient: f64,
    /// Open bore radius at the vertex: its impulse is tallied separately and reported as lost.
    pub hole_radius: f64,
}

/// The result of one point-fireball run.
#[derive(Debug, Clone, Copy)]
pub struct FireballBounce {
    /// Axial plate impulse over the whole plate, in units of `M u_rms`.
    pub wall_impulse: f64,
    /// The part of `wall_impulse` from inside the bore radius (lost through the real hole).
    pub hole_impulse: f64,
    /// Peak axial plate force.
    pub peak_force: f64,
    /// Peak local surface pressure, in units of `M u_rms^2 / m^3` (the `2π` restored).
    pub peak_local_pressure: f64,
    /// The ambient's mass as a share of the fireball's.
    pub ambient_mass_share: f64,
    pub steps: usize,
    /// Peak local pressure on the plate outside the open bore, same units.
    pub peak_pressure_outside_bore: f64,
    /// Where and when that peak came: its radius [m] and time (in the run's units).
    pub peak_radius: f64,
    pub peak_time: f64,
    /// Gas accounted for at the end (out through the outer faces plus still on the grid), as a
    /// share of the fireball's; short of 1 by what went down the bore or out the bottom.
    pub exit_mass_share: f64,
    /// Direction factor of that gas, `P_z / Σ m|u|`: 1 if every parcel leaves along the axis.
    pub exit_direction: f64,
    /// Speed-uniformity factor, `(Σ m|u|)^2 / (M Σ m u^2)`: 1 if every parcel has one speed.
    pub exit_uniformity: f64,
    /// What lands on the rings of [`COLLAR_BANDS`] (the bore collar's footprint).
    pub bands: [BandLoad; N_BANDS],
    /// Axial impulse on the plate within `(i + 1) PROFILE_DR` of the axis, in the units of
    /// `wall_impulse`: the fine radial profile the bands are too coarse to give.
    pub profile: [f64; N_PROFILE],
}

/// Number of radii in [`FireballBounce::profile`].
pub const N_PROFILE: usize = 64;
/// Spacing of [`FireballBounce::profile`]'s radii [m].
pub const PROFILE_DR: f64 = 0.25;

/// Number of radial bands in [`COLLAR_BANDS`].
pub const N_BANDS: usize = 6;

/// Radial bands `[r_lo, r_hi)` [m] tallied for the bore collar, beginning at the bore's edge and
/// running out over the 20 m plate, so the collar's edge is where a ring's load meets the plate's.
pub const COLLAR_BANDS: [(f64, f64); N_BANDS] = [
    (0.0, 0.6),
    (0.6, 1.0),
    (1.0, 2.0),
    (2.0, 3.5),
    (3.5, 6.0),
    (6.0, 10.0),
];

/// The load a ring of the plate takes over a pulse, per unit ring area and per unit (physical, `2π`
/// restored) fireball mass.
/// The arriving flow is read [`SAMPLE_CELLS`] cells above the surface, so the wall's own stagnated
/// layer does not hide it; only gas moving toward the plate counts.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BandLoad {
    /// Inner and outer radius [m] (the inner one raised to the bore's edge).
    pub r_lo: f64,
    pub r_hi: f64,
    /// Axial impulse on the ring per area.
    pub impulse: f64,
    /// Mass arriving per area, and its momentum-weighted speed `Σ m u_n^2 / Σ m u_n`.
    pub mass_in: f64,
    pub speed_in: f64,
    /// Total energy (kinetic plus enthalpy) carried in per area.
    pub energy_in: f64,
    /// Peak surface pressure on the ring.
    pub peak_pressure: f64,
}

/// Cells above the plate's surface where a band's arriving flow is read.
const SAMPLE_CELLS: usize = 3;

/// The surface fluid cell of column `ir` (the first fluid cell above a solid one), if any.
fn surface_cell(g: &Grid2D, ir: usize) -> Option<usize> {
    (1..g.nz()).find(|&iz| !g.is_solid(iz, ir) && g.is_solid(iz - 1, ir))
}

/// One step's contribution to each band's arriving mass, momentum and energy (`2π` dropped), and
/// its current peak surface pressure.
fn band_step(
    g: &Grid2D,
    bands: &[(f64, f64); N_BANDS],
    dt: f64,
) -> [(f64, f64, f64, f64); N_BANDS] {
    let dr = g.dr();
    let gamma = g.gamma();
    let mut out = [(0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64); N_BANDS];
    for ir in 0..g.nr() {
        let r = (ir as f64 + 0.5) * dr;
        let Some(b) = bands.iter().position(|&(lo, hi)| r >= lo && r < hi) else {
            continue;
        };
        let Some(iz_s) = surface_cell(g, ir) else {
            continue;
        };
        out[b].3 = out[b].3.max(g.prim(iz_s, ir).p);
        let iz = iz_s + SAMPLE_CELLS;
        if iz >= g.nz() || g.is_solid(iz, ir) {
            continue;
        }
        let w = g.prim(iz, ir);
        let un = -w.uz;
        if un > 0.0 {
            let area = r * dr;
            let dm = w.rho * un * area * dt;
            let h = 0.5 * w.rho * (w.uz * w.uz + w.ur * w.ur) + gamma / (gamma - 1.0) * w.p;
            out[b].0 += dm;
            out[b].1 += dm * un;
            out[b].2 += h * un * area * dt;
        }
    }
    out
}

/// Mass, axial momentum, `Σ m|u|` and `Σ m u^2` carried out through the far-`z` and outer-`r`
/// faces in `dt` (`2π` dropped), from the outermost cells' states.
fn outflow_moments(g: &Grid2D, dt: f64) -> (f64, f64, f64, f64) {
    let (nz, nr) = (g.nz(), g.nr());
    let (dz, dr) = (g.dz(), g.dr());
    let (mut m, mut pz, mut ms, mut mu2) = (0.0, 0.0, 0.0, 0.0);
    let mut add = |w: crate::state::Prim, un: f64, area: f64| {
        if un > 0.0 {
            let dm = w.rho * un * area * dt;
            let speed = w.uz.hypot(w.ur);
            m += dm;
            pz += dm * w.uz;
            ms += dm * speed;
            mu2 += dm * speed * speed;
        }
    };
    for ir in 0..nr {
        let w = g.prim(nz - 1, ir);
        add(w, w.uz, (ir as f64 + 0.5) * dr * dr);
    }
    for iz in 0..nz {
        if g.is_solid(iz, nr - 1) {
            continue;
        }
        let w = g.prim(iz, nr - 1);
        add(w, w.ur, nr as f64 * dr * dz);
    }
    (m, pz, ms, mu2)
}

/// Mass, axial momentum, `Σ m|u|` and `Σ m u^2` of the fluid on the grid (`2π` dropped).
fn exit_moments(g: &Grid2D) -> (f64, f64, f64, f64) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let (mut m, mut pz, mut ms, mut mu2) = (0.0, 0.0, 0.0, 0.0);
    for iz in 0..g.nz() {
        for ir in 0..g.nr() {
            if g.is_solid(iz, ir) {
                continue;
            }
            let w = g.prim(iz, ir);
            let dm = w.rho * g.cell_volume(ir) / two_pi;
            let speed = w.uz.hypot(w.ur);
            m += dm;
            pz += dm * w.uz;
            ms += dm * speed;
            mu2 += dm * speed * speed;
        }
    }
    (m, pz, ms, mu2)
}

/// Integrate an initialized fireball grid onto its plate: the impulse, the bore's share, the peak
/// pressure and the departing gas's moments. Mass and impulse are per unit fireball mass.
fn run_fireball_grid(
    g: Grid2D,
    nz: usize,
    hole_radius: f64,
    fireball_mass: f64,
    ambient_mass_share: f64,
) -> FireballBounce {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut hole_impulse = 0.0;
    let mut hole_old = 0.0;
    let mut t_old = 0.0;
    let mut peak_outside = 0.0_f64;
    let mut peak_at = (0.0, 0.0);
    let mut out = (0.0, 0.0, 0.0, 0.0);
    let bands = COLLAR_BANDS.map(|(lo, hi)| (lo.max(hole_radius), hi));
    let mut band_force_old = [0.0; N_BANDS];
    let mut band_impulse = [0.0; N_BANDS];
    let mut band_flow = [(0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64); N_BANDS];
    let mut profile = [0.0; N_PROFILE];
    let mut profile_old = [0.0; N_PROFILE];
    let (r, g) = run_bounce_observed_final(g, nz, |grid, t| {
        let dt = t - t_old;
        for (i, (acc, old)) in profile.iter_mut().zip(profile_old.iter_mut()).enumerate() {
            let f = grid.plate_force_within((i + 1) as f64 * PROFILE_DR);
            *acc += 0.5 * dt * (*old + f);
            *old = f;
        }
        let step = band_step(grid, &bands, dt);
        for (b, &(lo, hi)) in bands.iter().enumerate() {
            let f = grid.plate_force_within(hi) - grid.plate_force_within(lo);
            band_impulse[b] += 0.5 * dt * (band_force_old[b] + f);
            band_force_old[b] = f;
            band_flow[b].0 += step[b].0;
            band_flow[b].1 += step[b].1;
            band_flow[b].2 += step[b].2;
            band_flow[b].3 = band_flow[b].3.max(step[b].3);
        }
        let hole_new = grid.plate_force_within(hole_radius);
        hole_impulse += 0.5 * (t - t_old) * (hole_old + hole_new);
        hole_old = hole_new;
        let (p_out, r_out) = grid.max_plate_pressure_beyond_at(hole_radius);
        if p_out > peak_outside {
            peak_outside = p_out;
            peak_at = (r_out, t);
        }
        let f = outflow_moments(grid, t - t_old);
        out = (out.0 + f.0, out.1 + f.1, out.2 + f.2, out.3 + f.3);
        t_old = t;
    });
    let grid_moments = exit_moments(&g);
    let (m, pz, ms, mu2) = (
        out.0 + grid_moments.0,
        out.1 + grid_moments.1,
        out.2 + grid_moments.2,
        out.3 + grid_moments.3,
    );
    FireballBounce {
        wall_impulse: r.wall_impulse / fireball_mass,
        hole_impulse: hole_impulse / fireball_mass,
        peak_force: r.peak_force / fireball_mass,
        peak_local_pressure: r.peak_local_pressure / two_pi / fireball_mass,
        ambient_mass_share,
        steps: r.steps,
        peak_pressure_outside_bore: peak_outside / two_pi / fireball_mass,
        peak_radius: peak_at.0,
        peak_time: peak_at.1,
        exit_mass_share: m / fireball_mass,
        exit_direction: if ms > 0.0 { pz / ms } else { 0.0 },
        exit_uniformity: if mu2 > 0.0 { ms * ms / (m * mu2) } else { 0.0 },
        bands: std::array::from_fn(|b| {
            let (lo, hi) = bands[b];
            let area = 0.5 * (hi * hi - lo * lo);
            let (dm, dp, de, pk) = band_flow[b];
            BandLoad {
                r_lo: lo,
                r_hi: hi,
                impulse: band_impulse[b] / area / (two_pi * fireball_mass),
                mass_in: dm / area / (two_pi * fireball_mass),
                speed_in: if dm > 0.0 { dp / dm } else { 0.0 },
                energy_in: de / area / (two_pi * fireball_mass),
                peak_pressure: pk / (two_pi * fireball_mass),
            }
        }),
        profile: profile.map(|j| j / fireball_mass),
    }
}

/// Run a [`PointFireball`] onto the dished or cupped plate of `cfg` (`cfg.shape` must be a
/// `Dish` or `Cup`; `gamma`, `r_plate` and the domain are read from `cfg`, the slug fields are not).
///
/// # Panics
///
/// If `cfg.shape` is the grid-aligned flat plate, which has no immersed profile.
#[must_use]
pub fn run_fireball_bounce(cfg: &SlugConfig, fb: PointFireball) -> FireballBounce {
    let dz = cfg.z_max / cfg.nz as f64;
    let dr = cfg.r_max / cfg.nr as f64;
    let mut g = Grid2D::new(cfg.nz, cfg.nr, dz, dr, cfg.gamma);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Transmissive;
    g.bc_zhi = Bc::Transmissive;
    g.bc_zlo = Bc::Transmissive;
    let profile = immersed_profile(cfg).expect("a dished or cupped plate");
    g.set_plate_profile(Some(profile));
    let z_c = profile.z_surface(0.0) + fb.height;

    let inside = |iz: usize, ir: usize| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        (z - z_c).hypot(r) < fb.radius
    };
    // Volumes with the common 2π dropped, so mass and impulse share one convention.
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut v_ball = 0.0;
    let mut v_domain = 0.0;
    for iz in 0..cfg.nz {
        for ir in 0..cfg.nr {
            let v = g.cell_volume(ir) / two_pi;
            v_domain += v;
            if inside(iz, ir) {
                v_ball += v;
            }
        }
    }
    let rho = 1.0 / v_ball;
    let p = (cfg.gamma - 1.0) * 0.5 / v_ball;
    let rho_amb = fb.ambient * rho;
    let p_amb = rho_amb * 1.0e-6;
    g.init(|iz, ir| {
        if inside(iz, ir) {
            Prim::new(rho, fb.drift, 0.0, p)
        } else {
            Prim::new(rho_amb, 0.0, 0.0, p_amb)
        }
    });

    run_fireball_grid(
        g,
        cfg.nz,
        fb.hole_radius,
        1.0,
        rho_amb * (v_domain - v_ball),
    )
}

/// The real collision behind a [`PointFireball`] (X-ray study, Result 10 mixing check): a disk
/// that has just threaded the bore, moving away from the plate at speed 1, about to strike a
/// cylindrical polyethylene charge of `k` times its mass. With `toward_plate` it is Result 11's
/// overtake instead: the disk arrives from beyond the charge, heading for the plate. Both are γ-law gas; only masses, areal
/// densities and the closing speed matter to the snowplow, so the disk is a resolvable slab of its
/// true areal density. Impulse is per unit **merged** mass; the caller reads it against that
/// mass's eventual `u_rms = sqrt(2 Q / M)`, as for the point fireball.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiskCharge {
    /// Disk and charge radius (they share the bore's footprint).
    pub radius: f64,
    /// Height of the charge's centre above the plate's vertex.
    pub height: f64,
    /// Charge length along the axis (solid polyethylene is short, foam long).
    pub charge_length: f64,
    /// Charge-to-disk mass ratio `k`.
    pub k: f64,
    /// The disk slab's thickness on the grid (several cells).
    pub disk_thickness: f64,
    /// Gap between the disk and the charge at the start.
    pub gap: f64,
    /// Ambient density as a fraction of the charge's.
    pub ambient: f64,
    /// Open bore radius at the vertex.
    pub hole_radius: f64,
    /// The overtake case (X-ray study, Result 11): the disk arrives from beyond the charge,
    /// moving *toward* the plate, so the merged mass drifts onto it. `false` is Result 10's disk
    /// that has come through the bore and strikes the charge's plate-side face.
    pub toward_plate: bool,
}

/// Internal energy on the grid, `sum p/(gamma - 1) dV` with the common `2π` dropped (the
/// convention the disk's unit mass is set in).
fn cold_heat(g: &Grid2D, gamma: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    (0..g.nz())
        .flat_map(|iz| (0..g.nr()).map(move |ir| (iz, ir)))
        .map(|(iz, ir)| g.prim(iz, ir).p / (gamma - 1.0) * g.cell_volume(ir) / two_pi)
        .sum()
}

/// Run a [`DiskCharge`] collision onto the dished or cupped plate of `cfg`.
///
/// # Panics
///
/// If `cfg.shape` is the grid-aligned flat plate, which has no immersed profile.
#[must_use]
pub fn run_disk_charge_bounce(cfg: &SlugConfig, dc: DiskCharge) -> FireballBounce {
    let dz = cfg.z_max / cfg.nz as f64;
    let dr = cfg.r_max / cfg.nr as f64;
    let mut g = Grid2D::new(cfg.nz, cfg.nr, dz, dr, cfg.gamma);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Transmissive;
    g.bc_zhi = Bc::Transmissive;
    g.bc_zlo = Bc::Transmissive;
    let profile = immersed_profile(cfg).expect("a dished or cupped plate");
    g.set_plate_profile(Some(profile));
    let z_vertex = profile.z_surface(0.0);
    let charge_lo = z_vertex + dc.height - 0.5 * dc.charge_length;
    let charge_hi = charge_lo + dc.charge_length;
    let (disk_lo, disk_hi, disk_speed) = if dc.toward_plate {
        let lo = charge_hi + dc.gap;
        (lo, lo + dc.disk_thickness, -1.0)
    } else {
        let hi = charge_lo - dc.gap;
        (hi - dc.disk_thickness, hi, 1.0)
    };
    let region = |iz: usize, ir: usize| -> u8 {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        if r >= dc.radius {
            0
        } else if z >= charge_lo && z < charge_hi {
            2
        } else if z >= disk_lo && z < disk_hi {
            1
        } else {
            0
        }
    };
    let two_pi = 2.0 * std::f64::consts::PI;
    let (mut v_disk, mut v_charge, mut v_domain) = (0.0, 0.0, 0.0);
    for iz in 0..cfg.nz {
        for ir in 0..cfg.nr {
            let v = g.cell_volume(ir) / two_pi;
            v_domain += v;
            match region(iz, ir) {
                1 => v_disk += v,
                2 => v_charge += v,
                _ => {}
            }
        }
    }
    // Disk mass 1 at speed 1, charge mass k at rest, both cold.
    let rho_disk = 1.0 / v_disk;
    let rho_charge = dc.k / v_charge;
    let p_cold = rho_disk / (cfg.gamma * cfg.mach * cfg.mach);
    let rho_amb = dc.ambient * rho_charge;
    g.init(|iz, ir| match region(iz, ir) {
        1 => Prim::new(rho_disk, disk_speed, 0.0, p_cold),
        // The charge starts at the disk's cold *pressure*. At the disk's temperature it would hold
        // k times the disk's specific heat: at Mach 20 and k = 30 that is ~28% of the collision's
        // own heat, and ~90% at k = 100 (the defect corrected 2026-10-15).
        // A long foam column holds ~50 disk volumes, so its cold pressure is scaled down to keep
        // its stored heat under a quarter of the disk's own (itself ~1% of the disk's energy).
        2 => Prim::new(
            rho_charge,
            0.0,
            0.0,
            p_cold * (0.25 * v_disk / v_charge).min(1.0),
        ),
        _ => Prim::new(rho_amb, 0.0, 0.0, p_cold * 1.0e-6),
    });
    assert!(
        cold_heat(&g, cfg.gamma) < 0.02 * 0.5,
        "the disk and charge must start cold against the disk's kinetic energy of 1/2"
    );
    let merged = 1.0 + dc.k;
    run_fireball_grid(
        g,
        cfg.nz,
        dc.hole_radius,
        merged,
        rho_amb * (v_domain - v_disk - v_charge) / merged,
    )
}

#[cfg(test)]
mod tests {
    use super::{
        PlateShape, PointFireball, SlugConfig, init_slug_grid, run_fireball_bounce,
        run_slug_bounce, taper_mass_factor, taper_sigma_stats, taper_weight,
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

    /// The X-ray study's point fireball (Result 10): a fireball receding faster from the plate
    /// delivers less impulse, the bore's share is a small part of it, and the ambient is light.
    #[test]
    fn point_fireball_recedes_and_the_bore_is_a_small_share() {
        let cfg = SlugConfig {
            r_plate: 2.0,
            r_max: 2.4,
            z_max: 2.4,
            nr: 48,
            nz: 48,
            shape: PlateShape::Dish { d_over_d: 0.30 },
            ..base_cfg()
        };
        let fireball = |drift: f64| PointFireball {
            height: 2.0 / (8.0 * 0.30),
            radius: 0.25,
            drift,
            ambient: 1.0e-8,
            hole_radius: 0.1,
        };
        let still = run_fireball_bounce(&cfg, fireball(0.0));
        let fleeing = run_fireball_bounce(&cfg, fireball(0.4));
        assert!(still.wall_impulse > 0.0);
        assert!(fleeing.wall_impulse < still.wall_impulse);
        assert!(still.hole_impulse > 0.0 && still.hole_impulse < 0.1 * still.wall_impulse);
        assert!(still.ambient_mass_share < 1.0e-3);
    }
}
