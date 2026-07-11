//! Rung-B production sweep (B5d-1): the high-v half of the paper's `f(v)` — `e_eff(ρ_impact)` at
//! 16 km/s — from the coupled radiation+conduction slug bounce ([`CoupledBounce`]).
//!
//! Each `ρ_impact` builds a cold water slug ([`Tube::slug_si`]) coasting into a rigid wall, fires
//! the coupled bounce, and records the restitution `e_eff` plus the radiative loss channels (1a
//! absorbed at the wall, 1b escaping to space). The grid is swept in parallel with rayon
//! (ADR-0002); results are written as JSONL — one JSON object per line, appendable/crash-resilient
//! (ADR-0019) — for the Python frontier extraction (B5d-2).
//!
//! At 16 km/s the stagnated gas is `τ ≫ 1` (design §3), so `e_eff` is EOS/gas-dynamics-dominated
//! and opacity-insensitive; the opacity in the loaded table is an interim bracket (B5c-2), and the
//! B5d-3 sensitivity sweep demonstrates `e_eff` does not move with it.
//!
//! **Conductive channel (2): operator landed (B-flux), high-v transport still gated.** The inviscid
//! Lagrangian gas has no gas-side thermal resistance, so a cold-wall semi-infinite [`Solid`] of very
//! high effusivity once extracted, in its very first step, more heat than the thin wall gas cell
//! contains — that over-drain zeroed the wall cell and collapsed the bounce (the original deferral,
//! 2026-06-22). The B-flux gas-side conduction operator ([`Solid::step_coupled`]) now fixes the
//! *mechanism*: the gas carries its own conductivity `k_gas`, so the interface flux is finite and the
//! over-drain cannot recur. But it engages only where the table provides `k_gas`, and the **high-v
//! plasma table has none** — plasma transport (Spitzer-like conductivity) is the deferred B-flux
//! sibling alongside the real opacity table. `e_eff` is loss-insensitive (0.63 with no conduction vs
//! 0.64 lossless at M≈30), so this `e_eff(ρ)` pass still runs with `wall = None` and reports
//! `loss_conductive = 0` pending that high-v transport data. (The low-v `CoolProp` table *does* carry
//! `k_gas`, so the low-v path activates the operator.)

use std::fs;
use std::io::Write as _;
use std::path::Path;

use euler2d::bounce::{PlateShape, SlugConfig, eta_capture, run_slug_bounce};
use hydro1d::conduction::Solid;
use hydro1d::eos::TableEos;
use hydro1d::kernel::{AblatingBounce, Ablation, CondensingBounce, CoupledBounce, Tube, Viscosity};
use hydro1d::radiation::{Limiter, RadConstants};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use tables::Table;

const TABLE_PATH: &str = "data/tables/water.json";
const RESULT_PATH: &str = "data/results/sweep.jsonl";
const TABLE_PATH_LOWV: &str = "data/tables/water_lowv.json";
const RESULT_PATH_LOWV: &str = "data/results/sweep_lowv.jsonl";
const RESULT_PATH_TRANS_EOS: &str = "data/results/sweep_transitional_eos.jsonl";
const RESULT_PATH_TRANS_RAD: &str = "data/results/sweep_transitional_rad.jsonl";
const RESULT_PATH_GEOMETRY: &str = "data/results/sweep_geometry.jsonl";
const RESULT_PATH_ABLATING: &str = "data/results/sweep_ablating.jsonl";
const RESULT_PATH_FROZEN_PROBE: &str = "data/results/frozen_probe.jsonl";
const RESULT_PATH_FROZEN: &str = "data/results/sweep_frozen.jsonl";
const TABLE_DIR_FROZEN: &str = "data/tables/frozen";
const TABLE_PATH_JUPITER: &str = "data/tables/water_jupiter.json";
const RESULT_PATH_JUPITER: &str = "data/results/sweep_jupiter.jsonl";
const RESULT_PATH_GEOMETRY_M40: &str = "data/results/sweep_geometry_m40.jsonl";

/// The production `ρ_impact` grid [kg/m³] (design §3-4): the impact-density axis of `e_eff(ρ)`.
const RHO_GRID: [f64; 4] = [0.16, 0.32, 0.48, 0.64];

/// The transitional-anchor velocity grid [m/s] (ADR-0012): dense across the ~5–9 km/s partial-
/// ionization window where `e_eff` may dip, plus context points up to the 16 km/s anchor. Below
/// ~5 km/s the high-v `eos_water` dissociation chemistry degrades (the low-v `CoolProp` package takes
/// over), so the sweep starts there; it reuses the same `water.json` table the 16 km/s pass loads.
const V_GRID: [f64; 8] = [
    5_000.0, 6_000.0, 7_000.0, 8_000.0, 9_000.0, 11_000.0, 13_000.0, 16_000.0,
];

/// Fixed configuration for the 16 km/s `e_eff(ρ)` pass (cited; the footprint/Σ and 3.2/8 km/s
/// anchors are deferred to their own rungs — see the plan's "Out of scope").
#[derive(Debug, Clone, Copy)]
struct Config {
    /// Impact speed [m/s].
    v: f64,
    /// Cold cloud temperature [K]; with `c_s(ρ, T₀)` it sets the (very high) incident Mach number.
    t0: f64,
    /// Initial slug column length [m] — a nominal column giving `τ ≫ 1`.
    length: f64,
    /// Gas cells.
    gas_cells: usize,
    /// Radiation constants (SI): `c` [m/s] and `a = 4σ/c` [J/m³/K⁴].
    consts: RadConstants,
    /// Flux limiter (production default Levermore–Pomraning).
    limiter: Limiter,
}

impl Config {
    /// The production configuration (16 km/s, cold 400 K cloud, ~300 gas cells). The conductive
    /// wall is deferred (see the module docstring), so the bounce runs with `wall = None`.
    fn production() -> Self {
        Self {
            v: 16_000.0,
            t0: 400.0,
            length: 1.0,
            gas_cells: 300,
            consts: RadConstants {
                c: 2.997_924_58e8,
                a: 7.565_733e-16,
            },
            limiter: Limiter::LevermorePomraning,
        }
    }
}

/// One JSONL output row: the swept input, the restitution result, and the ADR-0016 loss
/// decomposition (per unit wall area).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct Record {
    /// Impact density `ρ_impact` [kg/m³].
    rho_impact: f64,
    /// Impact speed `v` [m/s].
    v: f64,
    /// Effective restitution `e_eff = J_wall/p_in − 1` (ADR-0001): 0 = stick, 1 = elastic.
    e_eff: f64,
    /// Peak wall force during the bounce (`p + q`; peak dominated by the AV spike `≈ c_q·ρv²`).
    peak_wall_force: f64,
    /// Peak **physical** wall pressure (EOS `p` only, AV excluded) — the facesheet survivability
    /// load, `≈ (γ_eff+1)/2 · ρv²` (ADR-0010 correction). Defaults to 0 when reading pre-fix
    /// JSONL so stale data is detectable downstream.
    #[serde(default)]
    peak_wall_pressure: f64,
    /// Incident axial momentum `p_in`.
    incident_momentum: f64,
    /// Gas momentum still in flight at stop (the rebound).
    residual_momentum: f64,
    /// Time-integrated wall force `J_wall`.
    wall_impulse: f64,
    /// Channel 1a — radiation absorbed at the wall.
    loss_radiative_wall: f64,
    /// Channel 1b — radiation escaping to space at the re-expansion end.
    loss_escape_space: f64,
    /// Channel 2 — heat conducted into the wall solid. **Structurally 0 this pass**: the conductive
    /// channel is deferred to its own rung (see the module docstring). The field is retained so the
    /// JSONL schema is stable for when a gas-side boundary-layer conduction model lands.
    loss_conductive: f64,
    /// Channel 3 — energy carried off by condensate that stuck to the wall (Rung C, low-v). `0` for
    /// the high-v plasma pass (no condensation); the dominant loss at 3.2 km/s (ADR-0004).
    loss_condensation: f64,
}

/// Run one coupled bounce at `rho_impact` with `table` and `cfg`; return its JSONL row.
fn run_one(rho_impact: f64, table: &Table, cfg: &Config) -> Record {
    let eos = TableEos::new(table.clone());
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    // `wall = None`: the conductive channel is deferred to its own rung (module docstring).
    let result = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
    Record {
        rho_impact,
        v: cfg.v,
        e_eff: result.bounce.e_eff,
        peak_wall_force: result.bounce.peak_wall_force,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        incident_momentum: result.bounce.incident_momentum,
        residual_momentum: result.bounce.residual_momentum,
        wall_impulse: result.bounce.wall_impulse,
        loss_radiative_wall: result.loss_radiative_wall,
        loss_escape_space: result.loss_escape_space,
        loss_conductive: result.loss_conductive,
        loss_condensation: 0.0, // no condensation in the high-v plasma
    }
}

/// Sweep `rho_grid` in parallel (rayon), preserving the input order. Each run is independent, so the
/// result is deterministic regardless of the parallel schedule.
fn run_sweep(rho_grid: &[f64], table: &Table, cfg: &Config) -> Vec<Record> {
    rho_grid
        .par_iter()
        .map(|&rho| run_one(rho, table, cfg))
        .collect()
}

/// A cold conducting wall behind the rigid plate (B-flux, ADR-0005): a semi-infinite [`Solid`] of
/// `cells` cells over `depth` (m), initially at `t_init` (K), with diffusivity `alpha` (m²/s) and
/// conductivity `k` (W/m/K). The effusivity `√(kρc) = k/√α` sets the conductive loss and, with the
/// cold `t_init`, drives the near-wall cooling that condenses (and deposits) the gas.
#[derive(Debug, Clone, Copy)]
struct WallParams {
    cells: usize,
    depth: f64,
    t_init: f64,
    alpha: f64,
    k: f64,
}

/// Fixed configuration for the 3.2 km/s low-v anchor (Rung C / B-flux): the cold cloud is warm enough
/// that every swept `ρ` starts as single-phase vapor (incident Mach ≈ 6); radiation is off (design
/// §3). With a cold conducting `wall` (B-flux) the near-wall gas is cooled below `T_sat`, activating
/// the wall-deposition sink (channels 2 + 3); `wall = None` is the adiabatic upper bound.
#[derive(Debug, Clone, Copy)]
struct LowvConfig {
    v: f64,
    t0: f64,
    length: f64,
    gas_cells: usize,
    /// Wall sticking coefficient (baseline 1 — the pessimistic equilibrium bound, ADR-0004).
    alpha: f64,
    /// Cold conducting wall (B-flux); `None` for the adiabatic upper bound.
    wall: Option<WallParams>,
}

impl LowvConfig {
    fn anchor() -> Self {
        Self {
            v: 3200.0,
            t0: 450.0,
            length: 1.0,
            gas_cells: 300,
            alpha: 1.0,
            // SiC-like plate (the design baseline): k ≈ 120 W/m/K, ρc ≈ 2.24e6 J/m³/K ⇒
            // α ≈ 5.4e-5 m²/s, effusivity ≈ 1.6e4 ≫ water vapor's; cold at 300 K. Depth (5 mm) and
            // cell count keep the thermal penetration √(αt) ≪ depth over the ~µs bounce.
            wall: Some(WallParams {
                cells: 200,
                depth: 5.0e-3,
                t_init: 300.0,
                alpha: 5.4e-5,
                k: 120.0,
            }),
        }
    }
}

/// Run one condensing bounce at `rho_impact` (Rung C / B-flux low-v); return its JSONL row. Radiation
/// is off (design §3); with a wall the conductive (channel 2) and condensation (channel 3) losses are
/// both populated, otherwise only condensation (the adiabatic upper bound).
fn run_one_lowv(rho_impact: f64, table: &Table, cfg: &LowvConfig) -> Record {
    let eos = TableEos::new(table.clone());
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let wall = cfg
        .wall
        .map(|w| Solid::new(w.cells, w.depth, w.t_init, w.alpha, w.k));
    let result = CondensingBounce::new_with_wall(tube, cfg.alpha, wall).run();
    Record {
        rho_impact,
        v: cfg.v,
        e_eff: result.bounce.e_eff,
        peak_wall_force: result.bounce.peak_wall_force,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        incident_momentum: result.bounce.incident_momentum,
        residual_momentum: result.bounce.residual_momentum,
        wall_impulse: result.bounce.wall_impulse,
        loss_radiative_wall: 0.0, // radiation off at 3.2 km/s (design §3)
        loss_escape_space: 0.0,
        loss_conductive: result.loss_conductive, // B-flux: gas-side conduction into the cold plate
        loss_condensation: result.loss_condensation,
    }
}

/// Sweep the low-v anchor in parallel (rayon), preserving input order.
fn run_sweep_lowv(rho_grid: &[f64], table: &Table, cfg: &LowvConfig) -> Vec<Record> {
    rho_grid
        .par_iter()
        .map(|&rho| run_one_lowv(rho, table, cfg))
        .collect()
}

/// Run one **EOS-only** bounce at `rho_impact` (transitional anchor, ADR-0012): pure gas dynamics with
/// the `eos_water` EOS via [`Tube::run_bounce`] — *no* radiation or conduction. This isolates the
/// opacity-independent dissociation/ionization specific-heat feature (the part of the transitional
/// `e_eff(v)` that is computable without the deferred real opacity table). All loss channels are 0.
fn run_one_eos(rho_impact: f64, table: &Table, cfg: &Config) -> Record {
    let eos = TableEos::new(table.clone());
    let mut tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let bounce = tube.run_bounce();
    Record {
        rho_impact,
        v: cfg.v,
        e_eff: bounce.e_eff,
        peak_wall_force: bounce.peak_wall_force,
        peak_wall_pressure: bounce.peak_wall_pressure,
        incident_momentum: bounce.incident_momentum,
        residual_momentum: bounce.residual_momentum,
        wall_impulse: bounce.wall_impulse,
        loss_radiative_wall: 0.0, // radiation off: the EOS-only feature
        loss_escape_space: 0.0,
        loss_conductive: 0.0,
        loss_condensation: 0.0,
    }
}

/// Sweep the transitional `v × ρ` grid in parallel (ADR-0012), returning **two** row sets in input
/// order: the EOS-only curve ([`run_one_eos`]) and the radiation-on comparison curve ([`run_one`],
/// interim opacity, `wall = None` — the 16 km/s pass's config). The gap between them is the
/// radiative-uncertainty band (pending the real opacity table). `base` supplies everything but `v`,
/// which is swept from `v_grid`.
fn run_sweep_transitional(
    v_grid: &[f64],
    rho_grid: &[f64],
    table: &Table,
    base: &Config,
) -> (Vec<Record>, Vec<Record>) {
    let grid: Vec<(f64, f64)> = v_grid
        .iter()
        .flat_map(|&v| rho_grid.iter().map(move |&rho| (v, rho)))
        .collect();
    grid.par_iter()
        .map(|&(v, rho)| {
            let cfg = Config { v, ..*base };
            (run_one_eos(rho, table, &cfg), run_one(rho, table, &cfg))
        })
        .collect::<Vec<_>>()
        .into_iter()
        .unzip()
}

// ---- Frozen-recombination bounding sweep (audit finding 3) --------------------------------------
//
// The equilibrium EOS returns dissociation/ionization energy during the rebound; if recombination
// *freezes* (three-body rates collapse as the rebounding gas rarefies) that energy stays locked and
// `e_eff` drops below the equilibrium value — the one approximation the audit flagged as stacked
// optimistically. Two EOS-only bounding runs per transitional case bracket the freeze timing:
//
// - **frozen-rebound** (freeze *after* the plate): equilibrium in, composition frozen at the
//   mass-weighted turnaround state for the rebound (`Tube::run_bounce_frozen_rebound`, the
//   sudden-freeze splice) — the pessimistic bound.
// - **frozen-throughout** (freeze *before* the plate): pure molecular H2O, no chemical sink at all
//   (`data/tables/frozen/h2o.json`) — the optimistic bound.
//
// Two-stage pipeline: `--frozen-probe` records each case's turnaround state `(ρ*, T*)`; the Python
// table generator (`puffsat.tables --frozen-from-probe`) freezes the equilibrium composition there
// and emits one table per case into `data/tables/frozen/`; `--frozen` then runs the three curves.

/// One `--frozen-probe` row: the case axes plus the mass-weighted turnaround state the
/// frozen-composition table is generated at (and the equilibrium EOS-only `e_eff` for reference).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct FrozenProbeRecord {
    v: f64,
    rho_impact: f64,
    /// Mass-weighted mean density at global turnaround [kg/m³].
    rho_star: f64,
    /// Mass-weighted mean temperature at global turnaround [K].
    t_star: f64,
    /// Equilibrium EOS-only restitution (same run).
    e_eff_eq: f64,
}

/// One `--frozen` row: the equilibrium EOS-only curve and its two freeze-timing brackets.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
struct FrozenRecord {
    v: f64,
    rho_impact: f64,
    /// Equilibrium EOS-only restitution (chemistry returns its energy) — the study's curve.
    e_eff_eq: f64,
    /// Sudden-freeze-at-turnaround restitution (freeze *after* the plate) — pessimistic bound.
    e_eff_frozen_rebound: f64,
    /// Pure-H2O no-chemistry restitution (freeze *before* the plate) — optimistic bound.
    e_eff_frozen_all: f64,
    /// Freeze state the per-case table was generated at (echoed from this run's turnaround).
    rho_star: f64,
    t_star: f64,
    /// EOS-swap re-seed energy jump as a fraction of the incident kinetic energy — the splice
    /// consistency diagnostic (small ⇒ the single freeze composition represents the slug well).
    swap_energy_jump_frac: f64,
}

/// The per-case frozen-table path — the naming contract shared with
/// `puffsat.tables.frozen_table_name`.
fn frozen_table_path(v: f64, rho_impact: f64) -> String {
    // SAFE: v is a positive velocity grid value ≤ 16000, exactly representable and in range.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let v_int = v.round() as u64;
    format!("{TABLE_DIR_FROZEN}/v{v_int:05}_rho{rho_impact:.2}.json")
}

/// Probe one EOS-only bounce for its turnaround state (no swap).
fn run_one_frozen_probe(
    v: f64,
    rho_impact: f64,
    table: &Table,
    base: &Config,
) -> FrozenProbeRecord {
    let cfg = Config { v, ..*base };
    let eos = TableEos::new(table.clone());
    let mut tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let r = tube.run_bounce_frozen_rebound(None);
    FrozenProbeRecord {
        v,
        rho_impact,
        rho_star: r.rho_star,
        t_star: r.t_star,
        e_eff_eq: r.bounce.e_eff,
    }
}

/// Run one frozen-recombination case: the equilibrium EOS-only bounce, the sudden-freeze splice
/// (per-case frozen table), and the pure-H2O no-chemistry bounce.
fn run_one_frozen(
    v: f64,
    rho_impact: f64,
    table: &Table,
    frozen_tbl: &Table,
    h2o_tbl: &Table,
    base: &Config,
) -> FrozenRecord {
    let cfg = Config { v, ..*base };
    let slug = |eos_tbl: &Table| {
        Tube::slug_si(
            cfg.gas_cells,
            rho_impact,
            cfg.v,
            cfg.length,
            cfg.t0,
            TableEos::new(eos_tbl.clone()),
            Viscosity::VON_NEUMANN_RICHTMYER,
        )
    };

    let eq = slug(table).run_bounce();
    let frozen = slug(table).run_bounce_frozen_rebound(Some(TableEos::new(frozen_tbl.clone())));
    let all = slug(h2o_tbl).run_bounce();

    // Incident kinetic energy per unit wall area: ½ p_in v (p_in = ρLv).
    let ke_in = 0.5 * frozen.bounce.incident_momentum * cfg.v;
    FrozenRecord {
        v,
        rho_impact,
        e_eff_eq: eq.e_eff,
        e_eff_frozen_rebound: frozen.bounce.e_eff,
        e_eff_frozen_all: all.e_eff,
        rho_star: frozen.rho_star,
        t_star: frozen.t_star,
        swap_energy_jump_frac: frozen.swap_energy_jump / ke_in,
    }
}

/// Sweep the transitional `v × ρ` grid in parallel for the frozen-probe pass (input order).
fn run_sweep_frozen_probe(
    v_grid: &[f64],
    rho_grid: &[f64],
    table: &Table,
    base: &Config,
) -> Vec<FrozenProbeRecord> {
    let grid: Vec<(f64, f64)> = v_grid
        .iter()
        .flat_map(|&v| rho_grid.iter().map(move |&rho| (v, rho)))
        .collect();
    grid.par_iter()
        .map(|&(v, rho)| run_one_frozen_probe(v, rho, table, base))
        .collect()
}

// ---- Geometry sweep (Rung D follow-on): eta_capture(curvature × L/D × r_foot/R) -----------------
//
// The 2D `eta_capture` track (ADR-0003), driven by the radiation-free axisymmetric Euler kernel
// (`euler2d`, ADR-0023). `eta_capture` is scale-invariant (Euler + geometry, no intrinsic length),
// so the footprint is fixed at `r_foot = 1` WLOG and the plate radius follows from `r_foot/R`; the
// cloud length follows from `L/D = L/(2·r_foot)`. Each case is sized to its own cloud (so a fixed
// cell count gives comparable resolution) and forms `eta_capture = (J_wall/p_in)_free,dish /
// (J_wall/p_in)_confined,planewave`, the same-kernel ratio that cancels the common re-expansion.

/// Plate curvature axis: flat (`d/D = 0`) + the two shallow-concave depths (ADR-0021).
const GEO_D_OVER_D: [f64; 3] = [0.0, 0.10, 0.15];
/// Cloud aspect `L/D` (disk → roughly unity); long cylinders are survivability-driven and out of
/// this first sweep's scope.
const GEO_L_OVER_D: [f64; 3] = [0.3, 0.6, 1.0];
/// Footprint coverage `r_foot/R` (design §sweep 0.3–1.0; `R` fixed, the shared knob).
const GEO_RFOOT_OVER_R: [f64; 3] = [0.3, 0.5, 0.7];
/// Incident-Mach anchors. `eta_capture` is geometry-dominated and only weakly Mach-dependent, so two
/// anchors bracket that dependence; D7 pairs them with the 1D `e_eff` velocity anchors. The physical
/// incident Mach of the production cloud is M ≈ 21 at the 11 km/s dip and M ≈ 32 at 16 km/s
/// (T₀ = 400–450 K water vapor, c_s ≈ 500 m/s), so the anchors sit at the strong-shock end rather
/// than the marginal M = 5 the first pass used (2026-07 audit).
const GEO_MACH: [f64; 2] = [10.0, 20.0];

/// Fixed resolution for the geometry sweep (cells per case; the domain is sized per case so this is
/// a comparable resolution across cases). Coarse for wall-time; the `diag_*` tables refine.
#[derive(Debug, Clone, Copy)]
struct GeoConfig {
    gamma: f64,
    nr: usize,
    nz: usize,
}

impl GeoConfig {
    fn production() -> Self {
        // 112×80: the 2026-07 audit showed 56×40 is not converged for the deep-dish/tight-footprint
        // corner (eta 1.024 → 0.989 at 2×, 0.993 at 3×); 2× is within ~0.5% of 3×.
        Self {
            gamma: 1.4,
            nr: 112,
            nz: 80,
        }
    }
}

/// One geometry-sweep row: the case (curvature, shape, footprint, Mach) and its `eta_capture` with
/// the two restitution ratios it is formed from.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct GeoRecord {
    /// Plate depth-to-diameter ratio `d/D` (0 = flat).
    d_over_d: f64,
    /// Cloud aspect `L/D`.
    l_over_d: f64,
    /// Footprint coverage `r_foot/R`.
    r_foot_over_r: f64,
    /// Incident Mach number.
    mach: f64,
    /// `eta_capture = (J_wall/p_in)_free / (J_wall/p_in)_confined` (ADR-0003).
    eta_capture: f64,
    /// `J_wall/p_in` of the free (finite-cloud, dished) run.
    restitution_free: f64,
    /// `J_wall/p_in` of the confined (plane-wave) denominator.
    restitution_confined: f64,
    /// Peak axial plate force of the free run (the survivability proxy, reported not gated here).
    peak_force: f64,
    /// Peak *local* facesheet pressure of the free run — the concave focusing concentration the
    /// survivability frontier divides by its flat (`d/D = 0`) counterpart (Rung S, ADR-0010/0021).
    peak_local_pressure: f64,
}

/// Run one `eta_capture` case: a free dished bounce over its confined plane-wave denominator. The
/// domain is sized to the cloud (`r_foot = 1` WLOG); `r_plate = r_foot/(r_foot/R)`, `L = (L/D)·2`.
fn run_eta_case(
    d_over_d: f64,
    l_over_d: f64,
    r_foot_over_r: f64,
    mach: f64,
    cfg: &GeoConfig,
) -> GeoRecord {
    let r_foot = 1.0;
    let r_plate = r_foot / r_foot_over_r;
    let length = l_over_d * 2.0 * r_foot;
    let r_max = r_plate * 1.4; // room past the rim for gas to escape (§7)
    let depth = d_over_d * 2.0 * r_plate;
    let z_max = depth + 2.0 * length + 1.5; // dish + cloud + rebound headroom
    let free = run_slug_bounce(&SlugConfig {
        gamma: cfg.gamma,
        mach,
        r_foot,
        length,
        r_plate,
        r_max,
        z_max,
        nr: cfg.nr,
        nz: cfg.nz,
        confined: false,
        shape: PlateShape::Dish { d_over_d },
    });
    // The plane-wave denominator: same column (L, ρ, Mach), cloud fills the radius, flat plate.
    let confined = run_slug_bounce(&SlugConfig {
        gamma: cfg.gamma,
        mach,
        r_foot: r_max,
        length,
        r_plate: r_max,
        r_max,
        z_max,
        nr: 8,
        nz: cfg.nz,
        confined: true,
        shape: PlateShape::FlatGridAligned,
    });
    GeoRecord {
        d_over_d,
        l_over_d,
        r_foot_over_r,
        mach,
        eta_capture: eta_capture(&free, &confined),
        restitution_free: free.restitution_ratio(),
        restitution_confined: confined.restitution_ratio(),
        peak_force: free.peak_force,
        peak_local_pressure: free.peak_local_pressure,
    }
}

/// Sweep the full (curvature × `L/D` × `r_foot/R` × Mach) grid in parallel (rayon, ADR-0002).
fn run_geometry_sweep(cfg: &GeoConfig) -> Vec<GeoRecord> {
    let cases: Vec<(f64, f64, f64, f64)> = GEO_MACH
        .iter()
        .flat_map(|&m| {
            GEO_RFOOT_OVER_R.iter().flat_map(move |&rf| {
                GEO_L_OVER_D
                    .iter()
                    .flat_map(move |&ld| GEO_D_OVER_D.iter().map(move |&dd| (dd, ld, rf, m)))
            })
        })
        .collect();
    cases
        .par_iter()
        .map(|&(dd, ld, rf, m)| run_eta_case(dd, ld, rf, m, cfg))
        .collect()
}

// ---- Ablating-wall recovery sweep (Rung E, ADR-0014) --------------------------------------------
//
// The Phase-2 best-estimate refinement above the rigid-wall floor (ADR-0013): a thin vapor layer
// boils off the plate and (E3) shields incoming radiation before it reaches the cold wall, recovering
// the radiative wall loss as near-wall pressure → bounce. Each case runs the rigid `CoupledBounce`
// (the conservative floor) and the `AblatingBounce` (shielding + mass injection) at the same config,
// so `recovery = e_eff_ablating − e_eff_rigid` is the pure ablating gain. The interim Kramers opacity
// is structurally `τ ≫ 1` at the anchors (ADR-0012), where radiation is trapped and the shield has
// little to do; the opacity-scale knob ([`Table::with_opacity_scale`]) scales it down to manufacture
// the wall-reaching `τ ≲ 1` regime, so the recovery is reported as a **τ-bracket** over the scale.
// Runs `wall = None` (the high-v table carries no `k_gas`; conduction — hence blowing — is off).

/// The two velocity anchors [m/s]: the transitional dip (~11 km/s, the worst-case `e_eff`, ADR-0012)
/// and the 16 km/s high-v anchor (the `f = 0.8` recovery-lever decision, ADR-0009).
const ABL_V: [f64; 2] = [11_000.0, 16_000.0];
/// Effective heat-of-ablation axis [J/kg] (Q*, 2–10 MJ/kg; ADR-0014). Larger Q* ⇒ less mass boils off
/// per unit incident flux ⇒ a thinner shield.
const ABL_Q_STAR: [f64; 3] = [2.0e6, 5.0e6, 10.0e6];
/// Opacity-scale axis spanning the τ regime (ADR-0012): the interim Kramers opacity (1×, `τ ≫ 1`)
/// scaled down toward the wall-reaching `τ ≲ 1` window where shielding recovers the most, plus a 10×
/// over-trapped point to bound the high-τ end.
const ABL_OPACITY_SCALE: [f64; 4] = [0.01, 0.1, 1.0, 10.0];
/// Vapor curtain gray opacity `κ_vapor` [m²/kg] (E3): the ablation-product absorbing layer's optical
/// depth is `τ_v = κ_vapor · ablated_mass`. Calibrated (the `diag_ablating_magnitudes` probe) so a
/// quasi-steady ablated mass gives an `O(1)` curtain at the anchors — the regime where the shield is
/// a meaningful but not saturating correction.
const ABL_KAPPA_VAPOR: f64 = 2.0e2;
/// Gas cells for the ablating sweep: `e_eff` is an integral momentum ratio and converges fast, so a
/// coarser-than-production mesh keeps the (v × ρ × scale × Q*) grid tractable.
const ABL_GAS_CELLS: usize = 200;

/// One ablating-sweep row: the case axes plus the rigid floor, the ablating restitution, and the
/// ablation bookkeeping (per unit wall area). `recovery = e_eff_ablating − e_eff_rigid`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct AblatingRecord {
    /// Impact speed [m/s].
    v: f64,
    /// Impact density `ρ_impact` [kg/m³].
    rho_impact: f64,
    /// Opacity scale applied to the table ([`Table::with_opacity_scale`]); 1 = interim Kramers.
    opacity_scale: f64,
    /// Effective heat of ablation `Q*` [J/kg].
    q_star: f64,
    /// Vapor curtain gray opacity `κ_vapor` [m²/kg].
    kappa_vapor: f64,
    /// Rigid coupled-bounce `e_eff` at this `(v, ρ, opacity_scale)` — the conservative floor.
    e_eff_rigid: f64,
    /// Ablating-wall `e_eff` (shielding + mass injection) — the best estimate.
    e_eff_ablating: f64,
    /// `e_eff_ablating − e_eff_rigid` — the pure ablating recovery.
    recovery: f64,
    /// Vapor mass boiled off and injected at the wall, per unit area [kg/m²].
    ablated_mass: f64,
    /// `ablated_mass / (ρ_impact · length)` — the quasi-steady fraction of the cloud (ADR-0014).
    ablated_fraction: f64,
    /// Channel 1a (shielded) — radiation that reaches the plate.
    loss_radiative_wall: f64,
    /// Channel 1b — radiation escaping to space.
    loss_escape_space: f64,
    /// Energy spent ablating `Σ Q*·ṁ·dt` [J/m²].
    loss_ablation: f64,
    /// Peak wall force of the ablating bounce.
    peak_wall_force: f64,
}

/// Build the cold water slug for the ablating sweep at `(rho_impact, v)` on `table` (already
/// opacity-scaled), with `cfg`'s mesh, length, and `t0`.
fn ablating_slug(rho_impact: f64, cfg: &Config, table: &Table) -> Tube<TableEos> {
    Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        TableEos::new(table.clone()),
        Viscosity::VON_NEUMANN_RICHTMYER,
    )
}

/// Run the rigid floor and the three `Q*` ablating bounces at one `(v, ρ, opacity_scale)` case;
/// return the three rows (one per `Q*`). The rigid floor is computed once and shared across `Q*`.
fn run_ablating_case(
    v: f64,
    rho_impact: f64,
    scale: f64,
    base: &Config,
    base_tbl: &Table,
) -> Vec<AblatingRecord> {
    let table = base_tbl.with_opacity_scale(scale);
    let cfg = Config { v, ..*base };
    let rigid = CoupledBounce::new(
        ablating_slug(rho_impact, &cfg, &table),
        None,
        cfg.consts,
        cfg.limiter,
    )
    .run();
    let cloud_mass = rho_impact * cfg.length;
    ABL_Q_STAR
        .iter()
        .map(|&q_star| {
            let ablation = Ablation::new(q_star, cfg.t0).with_vapor_opacity(ABL_KAPPA_VAPOR);
            let abl = AblatingBounce::new(
                ablating_slug(rho_impact, &cfg, &table),
                None,
                cfg.consts,
                cfg.limiter,
                ablation,
            )
            .run();
            AblatingRecord {
                v,
                rho_impact,
                opacity_scale: scale,
                q_star,
                kappa_vapor: ABL_KAPPA_VAPOR,
                e_eff_rigid: rigid.bounce.e_eff,
                e_eff_ablating: abl.bounce.e_eff,
                recovery: abl.bounce.e_eff - rigid.bounce.e_eff,
                ablated_mass: abl.ablated_mass,
                ablated_fraction: abl.ablated_mass / cloud_mass,
                loss_radiative_wall: abl.loss_radiative_wall,
                loss_escape_space: abl.loss_escape_space,
                loss_ablation: abl.loss_ablation,
                peak_wall_force: abl.bounce.peak_wall_force,
            }
        })
        .collect()
}

/// The real high-v EOS/opacity table (`data/tables/water.json`), loaded relative to the workspace
/// root the sweep binary is launched from (the Makefile runs it there).
fn base_table() -> Table {
    Table::load(TABLE_PATH).expect("load data/tables/water.json (run `make tables` first)")
}

/// Sweep the (v × ρ × opacity-scale) grid in parallel (rayon, ADR-0002); each case emits one row per
/// `Q*`. Input order is preserved (v outer, then ρ, then scale). The base table is loaded once and
/// opacity-scaled per case in-process.
fn run_ablating_sweep(base: &Config, base_tbl: &Table) -> Vec<AblatingRecord> {
    let cases: Vec<(f64, f64, f64)> = ABL_V
        .iter()
        .flat_map(|&v| {
            RHO_GRID
                .iter()
                .flat_map(move |&rho| ABL_OPACITY_SCALE.iter().map(move |&scale| (v, rho, scale)))
        })
        .collect();
    cases
        .par_iter()
        .flat_map(|&(v, rho, scale)| run_ablating_case(v, rho, scale, base, base_tbl))
        .collect()
}

// ---- Jupiter-retrograde 69 km/s scenario sweep (special scenario, 2026-07) ----------------------
//
// A 100 kg pulse at 69 km/s (Jupiter-retrograde encounter, "Sorry No ISRU" launch-capability
// study) on a large plate (<= 100 t). The survivability ceiling `rho <= P_limit/(c_stag v²)`
// lands at ~0.07 kg/m³ — dilute enough that the stagnated slab sits near tau ~ 1 instead of the
// design's tau >> 1, so `e_eff` is *opacity-sensitive* here and the sweep brackets it over the
// table's kappa scale. The extended-grid table (multi-stage O Saha ladder, T to 1.2e6 K) keeps
// the ~1.5e5 K stagnation state on-grid. Slug length is swept too: the survivable cloud is
// physically ~10 m long, and the radiative loss integrates over the (longer) bounce time.

/// Impact speed [m/s] of the Jupiter-retrograde scenario.
const JUP_V: f64 = 69_000.0;
/// Impact-density grid [kg/m³] spanning the 400 MPa survivability ceiling (~0.07) both ways;
/// 0.01 probes the dilute optically-thin corner where the coarse grid's e_eff was still rising.
const JUP_RHO: [f64; 6] = [0.01, 0.02, 0.04, 0.07, 0.11, 0.16];
/// Slug column lengths [m]: log-spaced from the 1 m production convention to the realistic
/// stretched cloud (~10 m), filling the ridge between the coarse grid's two endpoints.
const JUP_LENGTH: [f64; 5] = [1.0, 2.0, 4.0, 8.0, 12.0];
/// Opacity-scale bracket (tau ~ 1 regime: e_eff genuinely moves with opacity here); log-spaced,
/// keeping the 0.1/1/10 anchors the frontier's kappa bracket slices on.
const JUP_OPACITY_SCALE: [f64; 5] = [0.1, 0.3, 1.0, 3.0, 10.0];

/// One Jupiter-scenario row: the swept `(rho, length, opacity_scale)` case at 69 km/s with the
/// restitution, the physical peak wall pressure (plate sizing), and the radiative loss split.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct JupiterRecord {
    v: f64,
    rho_impact: f64,
    /// Slug column length [m] (production convention is 1 m; the survivable cloud is ~10 m).
    length: f64,
    /// Opacity scale applied to the table (tau ~ 1 bracket).
    opacity_scale: f64,
    e_eff: f64,
    /// Peak physical wall pressure (EOS `p`, AV excluded) — the survivability load (ADR-0010).
    peak_wall_pressure: f64,
    incident_momentum: f64,
    wall_impulse: f64,
    loss_radiative_wall: f64,
    loss_escape_space: f64,
}

/// Run one 69 km/s coupled bounce at `(rho, length)` on an opacity-scaled table.
fn run_one_jupiter(rho_impact: f64, length: f64, scale: f64, base_tbl: &Table) -> JupiterRecord {
    let table = base_tbl.with_opacity_scale(scale);
    let cfg = Config {
        v: JUP_V,
        length,
        ..Config::production()
    };
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        TableEos::new(table),
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let result = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
    JupiterRecord {
        v: cfg.v,
        rho_impact,
        length,
        opacity_scale: scale,
        e_eff: result.bounce.e_eff,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        incident_momentum: result.bounce.incident_momentum,
        wall_impulse: result.bounce.wall_impulse,
        loss_radiative_wall: result.loss_radiative_wall,
        loss_escape_space: result.loss_escape_space,
    }
}

/// Sweep the Jupiter (rho × length × opacity-scale) grid in parallel (rayon), input order.
fn run_jupiter_sweep(base_tbl: &Table) -> Vec<JupiterRecord> {
    let cases: Vec<(f64, f64, f64)> = JUP_RHO
        .iter()
        .flat_map(|&rho| {
            JUP_LENGTH
                .iter()
                .flat_map(move |&len| JUP_OPACITY_SCALE.iter().map(move |&s| (rho, len, s)))
        })
        .collect();
    cases
        .par_iter()
        .map(|&(rho, len, s)| run_one_jupiter(rho, len, s, base_tbl))
        .collect()
}

/// Write `records` as JSONL (one object per line, ADR-0019), creating the parent dir and replacing
/// the file. Generic over the row type so the transitional (`Record`) and geometry (`GeoRecord`)
/// sweeps share it.
fn write_rows<T: Serialize>(path: &str, records: &[T]) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent)?;
    }
    let mut out = fs::File::create(path)?;
    for r in records {
        writeln!(out, "{}", serde_json::to_string(r)?)?;
    }
    Ok(())
}

/// `--frozen-probe`: run the transitional grid EOS-only and record each case's mass-weighted
/// turnaround state; the Python table generator freezes the composition there.
fn frozen_probe_mode() -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH)?;
    let rows = run_sweep_frozen_probe(&V_GRID, &RHO_GRID, &table, &Config::production());
    write_rows(RESULT_PATH_FROZEN_PROBE, &rows)?;
    for r in &rows {
        println!(
            "rust: v={:.0} rho={:.2} -> turnaround rho*={:.3} T*={:.0} K (e_eff_eq={:.4})",
            r.v, r.rho_impact, r.rho_star, r.t_star, r.e_eff_eq,
        );
    }
    println!(
        "rust: wrote {} probe rows -> {RESULT_PATH_FROZEN_PROBE}",
        rows.len()
    );
    Ok(())
}

/// `--frozen`: the three-curve frozen-recombination bounding sweep (equilibrium vs
/// sudden-freeze-at-turnaround vs pure-H2O-no-chemistry). Per-case frozen tables are loaded
/// serially up front (cheap next to the bounces themselves).
fn frozen_sweep_mode() -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH)?;
    let h2o_tbl = Table::load(format!("{TABLE_DIR_FROZEN}/h2o.json"))?;
    let base = Config::production();
    let mut cases = Vec::new();
    for &v in &V_GRID {
        for &rho in &RHO_GRID {
            cases.push((v, rho, Table::load(frozen_table_path(v, rho))?));
        }
    }
    let rows: Vec<FrozenRecord> = cases
        .par_iter()
        .map(|(v, rho, frozen_tbl)| run_one_frozen(*v, *rho, &table, frozen_tbl, &h2o_tbl, &base))
        .collect();
    write_rows(RESULT_PATH_FROZEN, &rows)?;
    for r in &rows {
        println!(
            "rust: v={:.0} rho={:.2} -> e_eff eq={:.4} frozen-rebound={:.4} frozen-all={:.4} (jump {:.2e})",
            r.v,
            r.rho_impact,
            r.e_eff_eq,
            r.e_eff_frozen_rebound,
            r.e_eff_frozen_all,
            r.swap_energy_jump_frac,
        );
    }
    println!(
        "rust: wrote {} frozen rows -> {RESULT_PATH_FROZEN}",
        rows.len()
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // `--lowv` selects the 3.2 km/s condensing anchor (Rung C); `--transitional` the ADR-0012 velocity
    // sweep (two files); otherwise the 16 km/s high-v pass. Optional positional `[table] [result]`
    // override the high-v/low-v defaults (the B5d-3 opacity scan sweeps each scaled table into its own
    // JSONL).
    let args: Vec<String> = std::env::args().skip(1).collect();

    // Transitional anchor (ADR-0012): sweep V_GRID × RHO_GRID with the high-v table, emitting the
    // EOS-only curve and the radiation-on comparison curve into two fixed files.
    if args.iter().any(|a| a == "--transitional") {
        let table = Table::load(TABLE_PATH)?;
        let (eos, rad) = run_sweep_transitional(&V_GRID, &RHO_GRID, &table, &Config::production());
        write_rows(RESULT_PATH_TRANS_EOS, &eos)?;
        write_rows(RESULT_PATH_TRANS_RAD, &rad)?;
        for (e, r) in eos.iter().zip(rad.iter()) {
            println!(
                "rust: v={:.0} rho={:.2} -> e_eff_eos={:.4} e_eff_rad={:.4} (rad band {:.4})",
                e.v,
                e.rho_impact,
                e.e_eff,
                r.e_eff,
                e.e_eff - r.e_eff,
            );
        }
        println!(
            "rust: wrote {} eos + {} rad rows -> {RESULT_PATH_TRANS_EOS} , {RESULT_PATH_TRANS_RAD}",
            eos.len(),
            rad.len(),
        );
        return Ok(());
    }

    // Frozen-recombination probe: record each transitional case's turnaround state, from which the
    // Python side generates the per-case frozen-composition tables.
    if args.iter().any(|a| a == "--frozen-probe") {
        return frozen_probe_mode();
    }

    // Frozen-recombination bounding sweep: equilibrium vs sudden-freeze-at-turnaround vs
    // pure-H2O-no-chemistry, per transitional case. Needs the per-case frozen tables (make
    // tables-frozen).
    if args.iter().any(|a| a == "--frozen") {
        return frozen_sweep_mode();
    }

    // Ablating-wall recovery sweep (Rung E, ADR-0014): the rigid floor vs the shielding+injection
    // ablating wall over (v × ρ × opacity-scale × Q*), reporting recovery as a τ-bracket.
    if args.iter().any(|a| a == "--ablating") {
        let base_tbl = base_table();
        let cfg = Config {
            gas_cells: ABL_GAS_CELLS,
            ..Config::production()
        };
        let rows = run_ablating_sweep(&cfg, &base_tbl);
        write_rows(RESULT_PATH_ABLATING, &rows)?;
        for r in &rows {
            println!(
                "rust: v={:.0} rho={:.2} scale={:.2} Q*={:.1e} -> e_eff rigid={:.4} abl={:.4} (recovery {:+.4}, ablated {:.2}%)",
                r.v,
                r.rho_impact,
                r.opacity_scale,
                r.q_star,
                r.e_eff_rigid,
                r.e_eff_ablating,
                r.recovery,
                100.0 * r.ablated_fraction,
            );
        }
        println!(
            "rust: wrote {} ablating rows -> {RESULT_PATH_ABLATING}",
            rows.len()
        );
        return Ok(());
    }

    // Jupiter-retrograde 69 km/s scenario sweep: e_eff(rho × length × opacity-scale) with the
    // extended-grid table (multi-stage O ladder). Needs `make tables-jupiter` first.
    if args.iter().any(|a| a == "--jupiter") {
        let table = Table::load(TABLE_PATH_JUPITER)?;
        let rows = run_jupiter_sweep(&table);
        write_rows(RESULT_PATH_JUPITER, &rows)?;
        for r in &rows {
            println!(
                "rust: rho={:.3} L={:>4.1} scale={:>5.2} -> e_eff={:.4} peak_p={:.3e} (1a={:.3e} 1b={:.3e})",
                r.rho_impact,
                r.length,
                r.opacity_scale,
                r.e_eff,
                r.peak_wall_pressure,
                r.loss_radiative_wall,
                r.loss_escape_space,
            );
        }
        println!(
            "rust: wrote {} jupiter rows -> {RESULT_PATH_JUPITER}",
            rows.len()
        );
        return Ok(());
    }

    // High-Mach spot check for the Jupiter scenario: the full geometry grid at M = 40 (the
    // strong-shock plateau check past the production anchors 10/20).
    if args.iter().any(|a| a == "--geometry-m40") {
        let cfg = GeoConfig::production();
        let cases: Vec<(f64, f64, f64)> = GEO_RFOOT_OVER_R
            .iter()
            .flat_map(|&rf| {
                GEO_L_OVER_D
                    .iter()
                    .flat_map(move |&ld| GEO_D_OVER_D.iter().map(move |&dd| (dd, ld, rf)))
            })
            .collect();
        let rows: Vec<GeoRecord> = cases
            .par_iter()
            .map(|&(dd, ld, rf)| run_eta_case(dd, ld, rf, 40.0, &cfg))
            .collect();
        write_rows(RESULT_PATH_GEOMETRY_M40, &rows)?;
        for r in &rows {
            println!(
                "rust: d/D={:.2} L/D={:.2} r_foot/R={:.2} M=40 -> eta_capture={:.4}",
                r.d_over_d, r.l_over_d, r.r_foot_over_r, r.eta_capture,
            );
        }
        println!(
            "rust: wrote {} M=40 geometry rows -> {RESULT_PATH_GEOMETRY_M40}",
            rows.len()
        );
        return Ok(());
    }

    // Geometry sweep (Rung D follow-on): eta_capture(curvature × L/D × r_foot/R) from the euler2d
    // kernel; no EOS/opacity table needed (radiation-free, effective-γ, ADR-0008).
    if args.iter().any(|a| a == "--geometry") {
        let rows = run_geometry_sweep(&GeoConfig::production());
        write_rows(RESULT_PATH_GEOMETRY, &rows)?;
        for r in &rows {
            println!(
                "rust: d/D={:.2} L/D={:.2} r_foot/R={:.2} M={:.0} -> eta_capture={:.4} (free {:.4} / confined {:.4}, peak F={:.3e})",
                r.d_over_d,
                r.l_over_d,
                r.r_foot_over_r,
                r.mach,
                r.eta_capture,
                r.restitution_free,
                r.restitution_confined,
                r.peak_force,
            );
        }
        println!(
            "rust: wrote {} geometry rows -> {RESULT_PATH_GEOMETRY}",
            rows.len()
        );
        return Ok(());
    }

    let lowv = args.iter().any(|a| a == "--lowv");
    let positional: Vec<&String> = args.iter().filter(|a| !a.starts_with("--")).collect();
    let (def_table, def_result) = if lowv {
        (TABLE_PATH_LOWV, RESULT_PATH_LOWV)
    } else {
        (TABLE_PATH, RESULT_PATH)
    };
    let table_path = positional.first().map_or(def_table, |s| s.as_str());
    let result_path = positional.get(1).map_or(def_result, |s| s.as_str());

    let table = Table::load(table_path)?;
    let records = if lowv {
        run_sweep_lowv(&RHO_GRID, &table, &LowvConfig::anchor())
    } else {
        run_sweep(&RHO_GRID, &table, &Config::production())
    };

    if let Some(parent) = Path::new(result_path).parent() {
        fs::create_dir_all(parent)?;
    }
    // A fresh sweep replaces the file; one JSON object per line (ADR-0019).
    let mut out = fs::File::create(result_path)?;
    for r in &records {
        writeln!(out, "{}", serde_json::to_string(r)?)?;
        println!(
            "rust: rho={:.3} -> e_eff={:.4}  (peak F={:.3e}, losses 1a={:.3e} 1b={:.3e} 2={:.3e} 3={:.3e})",
            r.rho_impact,
            r.e_eff,
            r.peak_wall_force,
            r.loss_radiative_wall,
            r.loss_escape_space,
            r.loss_conductive,
            r.loss_condensation,
        );
    }
    println!("rust: wrote {} rows -> {result_path}", records.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        ABL_KAPPA_VAPOR, ABL_Q_STAR, AblatingRecord, Config, GeoConfig, GeoRecord, LowvConfig,
        Record, frozen_table_path, run_ablating_case, run_eta_case, run_one, run_one_frozen,
        run_sweep, run_sweep_frozen_probe, run_sweep_lowv, run_sweep_transitional,
    };
    use hydro1d::radiation::{Limiter, RadConstants};
    use tables::Table;

    const GAMMA: f64 = 1.4;

    /// A tiny ideal-gas table (`e = T`, so `c_v = 1`) with moderate Kramers opacities — power laws
    /// in `(ρ, T)`, so the loader's log-log interpolation is exact. Stands in for the real water
    /// table so the sweep machinery (`slug_si` → `CoupledBounce::run` → JSONL row) is testable fast.
    fn tiny_ideal_table() -> Table {
        let n: usize = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64)) // 0.01 … 10
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64)) // 0.05 … 200
            .collect();
        let (mut p, mut e, mut cs, mut kr, mut kp) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(t);
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
                kr.push(0.7 * r.powf(2.0) * t.powf(-3.5));
                kp.push(0.3 * r * t.powf(-2.0));
            }
        }
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": { "p": p, "e": e, "c_s": cs, "kappa_rosseland": kr, "kappa_planck": kp },
        });
        Table::from_json(&json.to_string()).unwrap()
    }

    /// A small, fast configuration (normalized units, like the kernel's coupled-bounce gate tests):
    /// `v = 1`, cold `T₀` giving Mach ≈ 5, a coarse gas mesh, and **weak** radiation so the
    /// radiative channels are exercised (positive) yet remain a perturbation on the lossless ceiling
    /// — keeping `e_eff` inside the physical `(0, 1)` band. Like production, this runs `wall = None`
    /// (the conductive channel is deferred), so `loss_conductive` is structurally 0.
    fn tiny_config() -> Config {
        let mach = 5.0;
        Config {
            v: 1.0,
            t0: 1.0 / (GAMMA * (GAMMA - 1.0) * mach * mach), // p₀ = 1/(γM²); c₀ = v/M
            length: 1.0,
            gas_cells: 40,
            consts: RadConstants { c: 1.0, a: 1e-5 },
            limiter: Limiter::Fick,
        }
    }

    /// The driver produces one well-formed row per ρ, in input order, with a physical restitution
    /// `0 < e_eff < 1` and non-negative loss channels — the B5d-1 schema/invariant gate.
    // Exact float `==` is intentional here: rho/v/loss_conductive flow through verbatim (input
    // passthrough, exact-zero deferral), not via arithmetic, so equality is the correct assertion.
    #[allow(clippy::float_cmp)]
    #[test]
    fn sweep_rows_are_well_formed_and_physical() {
        let table = tiny_ideal_table();
        let cfg = tiny_config();
        let rho_grid = [1.0, 2.0];
        let records = run_sweep(&rho_grid, &table, &cfg);

        assert_eq!(records.len(), rho_grid.len());
        for (rec, &rho) in records.iter().zip(rho_grid.iter()) {
            assert_eq!(rec.rho_impact, rho); // order preserved
            assert_eq!(rec.v, cfg.v);
            assert!(rec.e_eff.is_finite());
            assert!(
                rec.e_eff > 0.0 && rec.e_eff < 1.0,
                "e_eff out of (0,1): {}",
                rec.e_eff
            );
            assert!(rec.incident_momentum > 0.0);
            assert!(rec.wall_impulse > 0.0);
            assert!(rec.peak_wall_force > 0.0);
            assert!(
                rec.peak_wall_pressure > 0.0 && rec.peak_wall_pressure < rec.peak_wall_force,
                "physical pressure peak must be positive and below the p+q peak"
            );
            assert!(rec.loss_radiative_wall >= 0.0);
            assert!(rec.loss_escape_space >= 0.0);
            assert_eq!(rec.loss_conductive, 0.0); // conductive channel deferred (wall = None)
            assert_eq!(rec.loss_condensation, 0.0); // no condensation in the high-v plasma
        }
    }

    /// One row round-trips through the JSONL schema (serialize → parse → equal): the Python reader
    /// (B5d-2) sees exactly the fields written.
    // Exact float `==` is intentional: the inputs round-trip verbatim through the decimal text.
    #[allow(clippy::float_cmp)]
    #[test]
    fn record_jsonl_roundtrips() {
        let table = tiny_ideal_table();
        let cfg = tiny_config();
        let rec = run_one(1.0, &table, &cfg);
        let line = serde_json::to_string(&rec).unwrap();
        let back: Record = serde_json::from_str(&line).unwrap();
        // Inputs and the headline result survive the JSONL hop exactly; the tiny loss values may
        // differ by a ULP through the decimal text, which the Python reader does not depend on.
        assert_eq!(back.rho_impact, rec.rho_impact);
        assert_eq!(back.v, rec.v);
        assert_eq!(back.e_eff.to_bits(), rec.e_eff.to_bits());
        // The schema carries the swept inputs and the three named loss channels.
        for key in [
            "rho_impact",
            "v",
            "e_eff",
            "peak_wall_force",
            "peak_wall_pressure",
            "loss_radiative_wall",
            "loss_escape_space",
            "loss_conductive",
            "loss_condensation",
        ] {
            assert!(line.contains(key), "missing field {key}");
        }
    }

    /// An `e = T` ideal-gas table with a `liquid_frac` ramp rising with compression, so the low-v
    /// condensing sweep actually condenses (mirrors the kernel's `condensing_table`).
    fn tiny_condensing_table() -> Table {
        let n: usize = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64))
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64))
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
                kr.push(1e-10);
                kp.push(1e-10);
                lf.push(((r - 2.0) / 20.0).clamp(0.0, 1.0));
            }
        }
        let json = serde_json::json!({
            "rho_grid": rho_grid, "T_grid": t_grid, "shape": [n, n],
            "fields": { "p": p, "e": e, "c_s": cs,
                        "kappa_rosseland": kr, "kappa_planck": kp, "liquid_frac": lf },
        });
        Table::from_json(&json.to_string()).unwrap()
    }

    /// The low-v (condensing) sweep produces well-formed rows: `0 < e_eff < 1`, only the condensation
    /// channel populated (radiation/conduction off), and it is exercised (`> 0`) where the gas
    /// condenses.
    #[allow(clippy::float_cmp)] // exact-zero deferral for the off channels; verbatim inputs
    #[test]
    fn lowv_sweep_rows_are_well_formed_and_condense() {
        let table = tiny_condensing_table();
        let mach = 5.0;
        let cfg = LowvConfig {
            v: 1.0,
            t0: 1.0 / (GAMMA * (GAMMA - 1.0) * mach * mach), // M ≈ 5
            length: 1.0,
            gas_cells: 40,
            alpha: 1.0,
            wall: None, // adiabatic plumbing check; conduction is exercised by the kernel tests
        };
        let rho_grid = [1.0, 2.0];
        let records = run_sweep_lowv(&rho_grid, &table, &cfg);

        assert_eq!(records.len(), rho_grid.len());
        for (rec, &rho) in records.iter().zip(rho_grid.iter()) {
            assert_eq!(rec.rho_impact, rho);
            assert!(
                rec.e_eff > 0.0 && rec.e_eff < 1.0,
                "e_eff out of (0,1): {}",
                rec.e_eff
            );
            // Only the condensation channel is active at low-v.
            assert_eq!(rec.loss_radiative_wall, 0.0);
            assert_eq!(rec.loss_escape_space, 0.0);
            assert_eq!(rec.loss_conductive, 0.0);
            assert!(rec.loss_condensation > 0.0, "no condensation loss recorded");
        }
    }

    /// The transitional sweep (ADR-0012) returns two row sets — EOS-only and radiation-on — over the
    /// full `v × ρ` grid in input order (`v` outer, `ρ` inner, matching `run_sweep_transitional`).
    /// Both are well-formed (`0 < e_eff < 1`); `v`/`ρ` pass through verbatim; the EOS-only curve is a
    /// lossless upper bound (its `e_eff ≥` the radiation-on point, which carries the radiative band);
    /// the EOS rows report zero loss while the rad rows' channels are non-negative.
    // Exact float `==` is intentional: `v`/`ρ` flow through verbatim and the EOS losses are exact 0.
    #[allow(clippy::float_cmp)]
    #[test]
    fn transitional_sweep_rows_well_formed() {
        let table = tiny_ideal_table();
        let base = tiny_config();
        // Normalized v grid (Mach = v/c_s(t0) ranges ≈ 4–6, all supersonic) × two densities.
        let v_grid = [0.8, 1.0, 1.2];
        let rho_grid = [0.5, 1.0];
        let (eos, rad) = run_sweep_transitional(&v_grid, &rho_grid, &table, &base);

        assert_eq!(eos.len(), v_grid.len() * rho_grid.len());
        assert_eq!(rad.len(), eos.len());

        for (idx, (re, rr)) in eos.iter().zip(rad.iter()).enumerate() {
            // Grid is `v` outer, `ρ` inner (see `run_sweep_transitional`).
            let v = v_grid[idx / rho_grid.len()];
            let rho = rho_grid[idx % rho_grid.len()];
            // The two curves are the same (v, ρ) point, differing only in physics.
            assert_eq!(re.v, v);
            assert_eq!(re.rho_impact, rho);
            assert_eq!(rr.v, v);
            assert_eq!(rr.rho_impact, rho);

            for rec in [re, rr] {
                assert!(
                    rec.e_eff > 0.0 && rec.e_eff < 1.0,
                    "e_eff out of (0,1): {}",
                    rec.e_eff
                );
                assert!(rec.incident_momentum > 0.0);
                assert!(rec.wall_impulse > 0.0);
            }

            // The two curves are the same (v, ρ) point under different physics, so they agree to
            // within this harness's radiative band. That band is tiny here because `tiny_config`
            // uses deliberately weak radiation (`a = 1e-5`) to keep `e_eff` in (0, 1); at this
            // coupling the FLD's internal transport can nudge `e_eff` either way, so we assert
            // closeness, not a strict ordering. The production-scale upper bound (EOS 0.64 ≥ rad
            // 0.63 at 16 km/s) is an end-to-end check; the kernel's losses-lower-e_eff gate pins
            // the loss sign at meaningful coupling.
            assert!(
                (re.e_eff - rr.e_eff).abs() < 1e-3,
                "EOS-only e_eff {} and radiation-on {} diverge beyond the weak-radiation band",
                re.e_eff,
                rr.e_eff
            );
            // EOS-only carries no losses; the radiation-on channels are non-negative.
            assert_eq!(re.loss_radiative_wall, 0.0);
            assert_eq!(re.loss_escape_space, 0.0);
            assert_eq!(re.loss_conductive, 0.0);
            assert_eq!(re.loss_condensation, 0.0);
            assert!(rr.loss_radiative_wall >= 0.0);
            assert!(rr.loss_escape_space >= 0.0);
        }

        // `v` actually varies across the grid (not a constant column).
        assert!(eos.iter().any(|r| r.v != eos[0].v));
    }

    /// The per-case frozen-table path matches the Python generator's naming contract
    /// (`puffsat.tables.frozen_table_name`): zero-padded integer velocity, two-decimal density.
    #[test]
    fn frozen_table_path_matches_python_contract() {
        assert_eq!(
            frozen_table_path(5_000.0, 0.16),
            "data/tables/frozen/v05000_rho0.16.json"
        );
        assert_eq!(
            frozen_table_path(16_000.0, 0.64),
            "data/tables/frozen/v16000_rho0.64.json"
        );
    }

    /// The frozen-probe sweep covers the grid in input order with physical turnaround states
    /// (compressed above ρ_impact, heated above T₀), and `run_one_frozen` degenerates correctly
    /// when every role is played by the *same* table: the three curves coincide and the splice
    /// diagnostic is ~0.
    #[allow(clippy::float_cmp)] // verbatim pass-through of the case axes
    #[test]
    fn frozen_sweep_rows_well_formed_and_degenerate_correctly() {
        let table = tiny_ideal_table();
        let base = tiny_config();
        let v_grid = [1.0, 1.2];
        let rho_grid = [1.0, 2.0];

        let probe = run_sweep_frozen_probe(&v_grid, &rho_grid, &table, &base);
        assert_eq!(probe.len(), v_grid.len() * rho_grid.len());
        for (idx, r) in probe.iter().enumerate() {
            assert_eq!(r.v, v_grid[idx / rho_grid.len()]);
            assert_eq!(r.rho_impact, rho_grid[idx % rho_grid.len()]);
            assert!(
                r.rho_star > r.rho_impact,
                "turnaround should be compressed: rho*={} vs rho={}",
                r.rho_star,
                r.rho_impact
            );
            assert!(r.t_star > base.t0);
            assert!(r.e_eff_eq > 0.0 && r.e_eff_eq < 1.0);
        }

        let rec = run_one_frozen(1.0, 1.0, &table, &table, &table, &base);
        assert!(
            (rec.e_eff_frozen_rebound - rec.e_eff_eq).abs() < 1e-9,
            "same-table splice should be a no-op: {} vs {}",
            rec.e_eff_frozen_rebound,
            rec.e_eff_eq
        );
        assert_eq!(rec.e_eff_frozen_all, rec.e_eff_eq);
        assert!(rec.swap_energy_jump_frac.abs() < 1e-9);
    }

    /// DIAGNOSTIC (ignored): load the real water table and isolate which loss channel breaks the
    /// high-Mach bounce. Run with `cargo test -p sweep -- --ignored --nocapture diag`.
    #[test]
    #[ignore = "diagnostic; needs data/tables/water.json"]
    fn diag_high_mach_channels() {
        use hydro1d::conduction::Solid;
        use hydro1d::eos::TableEos;
        use hydro1d::kernel::{CoupledBounce, Tube, Viscosity};
        let table = Table::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/tables/water.json"
        ))
        .unwrap();
        let cells = 150usize;
        let consts = RadConstants {
            c: 2.997_924_58e8,
            a: 7.565_733e-16,
        };
        let mk = || {
            Tube::slug_si(
                cells,
                0.32,
                16_000.0,
                1.0,
                400.0,
                TableEos::new(table.clone()),
                Viscosity::VON_NEUMANN_RICHTMYER,
            )
        };
        let wall = || Solid::new(400, 2.0e-3, 300.0, 120.0 / (3210.0 * 700.0), 120.0);
        for (name, w, a) in [
            ("full (wall+rad)", true, consts.a),
            ("no conduction  ", false, consts.a),
            ("lossless-ish   ", false, 1e-30),
        ] {
            let wall_opt = if w { Some(wall()) } else { None };
            let c2 = RadConstants { c: consts.c, a };
            let r = CoupledBounce::new(mk(), wall_opt, c2, Limiter::LevermorePomraning).run();
            println!(
                "{name}: e_eff={:+.4} resid/inc={:+.3} peak={:.2e} loss2={:.2e}",
                r.bounce.e_eff,
                r.bounce.residual_momentum / r.bounce.incident_momentum,
                r.bounce.peak_wall_force,
                r.loss_conductive,
            );
        }
    }

    /// The parallel sweep is deterministic: each run is independent, so repeated sweeps agree
    /// bit-for-bit on `e_eff`.
    #[test]
    fn sweep_is_deterministic() {
        let table = tiny_ideal_table();
        let cfg = tiny_config();
        let rho_grid = [1.0, 2.0];
        let a = run_sweep(&rho_grid, &table, &cfg);
        let b = run_sweep(&rho_grid, &table, &cfg);
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.e_eff.to_bits(), y.e_eff.to_bits());
        }
    }

    /// A tiny geometry config — coarse enough that `eta_capture` cases run fast in a unit test.
    fn tiny_geo() -> GeoConfig {
        GeoConfig {
            gamma: GAMMA,
            nr: 24,
            nz: 18,
        }
    }

    /// A geometry case yields a well-formed positive `eta_capture` from two restitution ratios, and
    /// the case parameters pass through into the row. The upper bound is generous (not 1): a concave
    /// plate can *over-collimate* — bend the rebound to be more axial than the flat plane-wave limit
    /// — so `eta_capture` legitimately exceeds 1 for the best short-disk/deep-dish cases.
    #[test]
    fn geometry_case_is_well_formed() {
        let r = run_eta_case(0.10, 0.6, 0.5, 5.0, &tiny_geo());
        assert!((r.d_over_d - 0.10).abs() < 1e-12 && (r.l_over_d - 0.6).abs() < 1e-12);
        assert!((r.r_foot_over_r - 0.5).abs() < 1e-12 && (r.mach - 5.0).abs() < 1e-12);
        assert!(
            r.eta_capture > 0.0 && r.eta_capture < 1.5,
            "eta_capture out of range: {}",
            r.eta_capture
        );
        assert!(r.restitution_free > 0.0 && r.restitution_confined > 0.0);
    }

    /// The curvature gain survives into the sweep driver: a shallow-concave plate captures more
    /// axial momentum than the flat plate at the same cloud — `eta_capture(0.15) > eta_capture(0)`.
    #[test]
    fn geometry_concave_beats_flat() {
        let cfg = tiny_geo();
        let flat = run_eta_case(0.0, 0.6, 0.5, 5.0, &cfg).eta_capture;
        let concave = run_eta_case(0.15, 0.6, 0.5, 5.0, &cfg).eta_capture;
        assert!(
            concave > flat,
            "expected concave to beat flat: flat {flat:.4}, concave {concave:.4}"
        );
    }

    /// A `GeoRecord` round-trips through the JSONL boundary (ADR-0019): the case parameters (round
    /// numbers) exactly, the computed `eta_capture` and ratios to round-off (`serde_json`'s default
    /// float parse can drift a ULP — the boundary is plaintext consumed as f64 by Python, so that is
    /// the contract, not bit-equality).
    #[test]
    #[allow(clippy::float_cmp)] // exact compares on the verbatim round-number case inputs
    fn geo_record_jsonl_roundtrip() {
        let r = run_eta_case(0.10, 0.6, 0.5, 5.0, &tiny_geo());
        let line = serde_json::to_string(&r).unwrap();
        let back: GeoRecord = serde_json::from_str(&line).unwrap();
        assert_eq!(back.d_over_d, r.d_over_d);
        assert_eq!(back.l_over_d, r.l_over_d);
        assert_eq!(back.r_foot_over_r, r.r_foot_over_r);
        assert_eq!(back.mach, r.mach);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-12 * a.abs().max(1.0);
        assert!(close(back.eta_capture, r.eta_capture));
        assert!(close(back.restitution_free, r.restitution_free));
        assert!(close(back.restitution_confined, r.restitution_confined));
        assert!(close(back.peak_force, r.peak_force));
        assert!(close(back.peak_local_pressure, r.peak_local_pressure));
    }

    /// One ablating case emits a well-formed row per `Q*`: the rigid floor and the ablating
    /// restitution are both physical `(0, 1)`, the recovery is their difference, the ablated mass is a
    /// non-negative small fraction, and the case axes (v, ρ, scale, κ_vapor) pass through. The rigid
    /// floor is shared across the `Q*` rows. Uses a tiny radiating table + coarse config for speed.
    #[allow(clippy::float_cmp)] // verbatim input passthrough (v/ρ/scale/Q*/κ_vapor), not arithmetic
    #[test]
    fn ablating_sweep_rows_well_formed() {
        let base_tbl = tiny_ideal_table();
        let cfg = tiny_config(); // v = 1, M ≈ 5, 40 cells, weak radiation
        let rows = run_ablating_case(cfg.v, 1.0, 1.0, &cfg, &base_tbl);

        assert_eq!(rows.len(), ABL_Q_STAR.len());
        for (rec, &q_star) in rows.iter().zip(ABL_Q_STAR.iter()) {
            assert_eq!(rec.v, cfg.v);
            assert_eq!(rec.rho_impact, 1.0);
            assert_eq!(rec.opacity_scale, 1.0);
            assert_eq!(rec.q_star, q_star);
            assert_eq!(rec.kappa_vapor, ABL_KAPPA_VAPOR);
            assert!(
                rec.e_eff_rigid > 0.0 && rec.e_eff_rigid < 1.0,
                "rigid e_eff out of (0,1): {}",
                rec.e_eff_rigid
            );
            assert!(
                rec.e_eff_ablating > 0.0 && rec.e_eff_ablating < 1.0,
                "ablating e_eff out of (0,1): {}",
                rec.e_eff_ablating
            );
            assert!((rec.recovery - (rec.e_eff_ablating - rec.e_eff_rigid)).abs() < 1e-12);
            assert!(rec.ablated_mass >= 0.0 && rec.ablated_fraction >= 0.0);
            assert!(rec.loss_radiative_wall >= 0.0 && rec.loss_escape_space >= 0.0);
        }
        // The rigid floor is the same physics at every Q* (shared across the rows).
        assert_eq!(rows[0].e_eff_rigid, rows[1].e_eff_rigid);
        assert_eq!(rows[1].e_eff_rigid, rows[2].e_eff_rigid);
    }

    /// An `AblatingRecord` round-trips through the JSONL boundary (ADR-0019): the Python `--axis
    /// ablating` reader sees exactly the fields written.
    #[allow(clippy::float_cmp)] // round-number inputs survive the decimal text verbatim
    #[test]
    fn ablating_record_jsonl_roundtrip() {
        let rec = AblatingRecord {
            v: 16_000.0,
            rho_impact: 0.32,
            opacity_scale: 0.1,
            q_star: 5.0e6,
            kappa_vapor: ABL_KAPPA_VAPOR,
            e_eff_rigid: 0.55,
            e_eff_ablating: 0.61,
            recovery: 0.06,
            ablated_mass: 4.0e-3,
            ablated_fraction: 0.0125,
            loss_radiative_wall: 1.0e6,
            loss_escape_space: 1.0e5,
            loss_ablation: 2.0e4,
            peak_wall_force: 3.0e8,
        };
        let line = serde_json::to_string(&rec).unwrap();
        let back: AblatingRecord = serde_json::from_str(&line).unwrap();
        assert_eq!(back, rec);
        for key in [
            "v",
            "rho_impact",
            "opacity_scale",
            "q_star",
            "kappa_vapor",
            "e_eff_rigid",
            "e_eff_ablating",
            "recovery",
            "ablated_mass",
            "ablated_fraction",
        ] {
            assert!(line.contains(key), "missing field {key}");
        }
    }

    /// DIAGNOSTIC (ignored): calibrate `κ_vapor` and the opacity-scale grid against the real water
    /// table — print the ablated mass/fraction, the (shielded) radiative wall loss, and the recovery
    /// at the two anchors over a κ_vapor × Q* × scale grid. Run with
    /// `cargo test -p sweep -- --ignored --nocapture diag_ablating`.
    #[test]
    #[ignore = "diagnostic; needs data/tables/water.json"]
    fn diag_ablating_magnitudes() {
        use hydro1d::kernel::{AblatingBounce, Ablation, CoupledBounce, Tube, Viscosity};
        use tables::Table;
        let base = Table::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/tables/water.json"
        ))
        .unwrap();
        let consts = super::Config::production().consts;
        let (t0, length, cells) = (400.0, 1.0, super::ABL_GAS_CELLS);
        let rho = 0.64;
        for v in super::ABL_V {
            for scale in [0.1, 1.0] {
                let table = base.with_opacity_scale(scale);
                let mk = || {
                    Tube::slug_si(
                        cells,
                        rho,
                        v,
                        length,
                        t0,
                        super::TableEos::new(table.clone()),
                        Viscosity::VON_NEUMANN_RICHTMYER,
                    )
                };
                let rigid =
                    CoupledBounce::new(mk(), None, consts, super::Limiter::LevermorePomraning)
                        .run()
                        .bounce
                        .e_eff;
                for kappa_vapor in [0.0, super::ABL_KAPPA_VAPOR, 800.0] {
                    for q_star in super::ABL_Q_STAR {
                        let abl = AblatingBounce::new(
                            mk(),
                            None,
                            consts,
                            super::Limiter::LevermorePomraning,
                            Ablation::new(q_star, t0).with_vapor_opacity(kappa_vapor),
                        )
                        .run();
                        println!(
                            "v={v:.0} scale={scale:.2} kv={kappa_vapor:.0} Q*={q_star:.0e}: rigid={rigid:+.4} abl={:+.4} rec={:+.4} m_abl={:.3e} ({:.2}%) loss1a={:.3e} tau_v={:.3}",
                            abl.bounce.e_eff,
                            abl.bounce.e_eff - rigid,
                            abl.ablated_mass,
                            100.0 * abl.ablated_mass / (rho * length),
                            abl.loss_radiative_wall,
                            kappa_vapor * abl.ablated_mass,
                        );
                    }
                }
            }
        }
    }
}
