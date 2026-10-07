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

mod progress;

use std::collections::HashMap;
use std::fs;
use std::io::Write as _;
use std::path::Path;

use euler2d::bounce::{
    PlateShape, SlugConfig, SprayCloud, TrailingPearl, eta_capture, init_spray_grid,
    run_slug_bounce, run_spray_bounce, run_spray_bounce_with_ambient, run_spray_bounce_with_pearl,
    run_spray_bounce_with_skirt_history, taper_sigma_stats,
};
use hydro1d::conduction::Solid;
use hydro1d::eos::{Eos as _, TableEos};
use hydro1d::kernel::{
    AblatingBounce, Ablation, Boundary, CondensingBounce, CoupledBounce, LagrangianState,
    TransportAudit, Tube, Viscosity,
};
use hydro1d::radiation::{Limiter, RadConstants};
use serde::{Deserialize, Serialize};
use tables::Table;

use crate::progress::par_map_with_progress;

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
// The 69 km/s freeze-timing bracket (ADR-0026 instrument at the realistic-cloud anchor).
const RESULT_PATH_FROZEN_PROBE_JUPITER: &str = "data/results/frozen_probe_jupiter.jsonl";
const RESULT_PATH_FROZEN_JUPITER: &str = "data/results/sweep_frozen_jupiter.jsonl";
const TABLE_DIR_FROZEN_JUPITER: &str = "data/tables/frozen_jupiter";
const RESULT_PATH_GEOMETRY_M40: &str = "data/results/sweep_geometry_m40.jsonl";
// Design §7's remaining `r_foot/R` nodes (0.8/0.9/1.0), in their own file: five frontier CSVs
// consume `sweep_geometry_m40.jsonl` and must not shift. Same `GeoRecord` schema, so the analysis
// concatenates the two.
const RESULT_PATH_GEOMETRY_WIDE: &str = "data/results/sweep_geometry_wide.jsonl";
const RESULT_PATH_GEOMETRY_WIDE_HEADROOM: &str = "data/results/sweep_geometry_wide_headroom.jsonl";
// The heavy-plate 16–28 km/s special scenario (design §12.1, ADR-0027): reuses the Jupiter
// extended-grid table; its own result files + freeze-bracket directory.
const RESULT_PATH_HEAVYPLATE: &str = "data/results/sweep_heavyplate.jsonl";
const RESULT_PATH_FROZEN_PROBE_HEAVYPLATE: &str = "data/results/frozen_probe_heavyplate.jsonl";
// Turnaround states at the *diagnostic* velocities, in a separate file on purpose. The frozen-table
// builder emits one Saha table per row of the probe it reads (`build_frozen_tables_from_probe`), so
// folding these into the Q4 probe above would build ~36 expensive tables that nothing consumes.
// These rows feed only the LTE check (Q5) and the opacity bracket (Q18), which need a turnaround
// (rho*, T*) at every state they report on — without them the new 45-63 km/s tau-check rows would be
// silently dropped for want of a composition.
const RESULT_PATH_PROBE_HEAVYPLATE_DIAG: &str = "data/results/frozen_probe_heavyplate_diag.jsonl";
const RESULT_PATH_FROZEN_HEAVYPLATE: &str = "data/results/sweep_frozen_heavyplate.jsonl";
const TABLE_DIR_FROZEN_HEAVYPLATE: &str = "data/tables/frozen_heavyplate";

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
    par_map_with_progress("baseline", rho_grid, |&rho| run_one(rho, table, cfg))
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
    par_map_with_progress("lowv", rho_grid, |&rho| run_one_lowv(rho, table, cfg))
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
    par_map_with_progress("transitional", &grid, |&(v, rho)| {
        let cfg = Config { v, ..*base };
        (run_one_eos(rho, table, &cfg), run_one(rho, table, &cfg))
    })
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
fn frozen_table_path(dir: &str, v: f64, rho_impact: f64) -> String {
    // SAFE: v is a positive velocity grid value ≤ 69000, exactly representable and in range.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let v_int = v.round() as u64;
    format!("{dir}/v{v_int:05}_rho{rho_impact:.2}.json")
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
    par_map_with_progress("frozen-probe", &grid, |&(v, rho)| {
        run_one_frozen_probe(v, rho, table, base)
    })
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
/// The rest of design §7's specified `r_foot/R` range. The original grid stopped at 0.7, which is
/// where `contour::R_FOOT_BOX` took its upper edge from — so the box edge was the *sweep's* limit,
/// not a delivery constraint. It binds: the 25 kg / 5 m core plate needs `r_foot/R` 0.87 at 45 km/s
/// and 1.00 at 55 km/s to hold the eta-optimal `L/D = 0.3` on the survivable contour, and clamping
/// it at 0.7 forces a cloud stretch that costs ~0.05 in `f`. Past `r_foot/R = 1` the cloud is wider
/// than the plate and the spill is unambiguous, so 1.0 is the physical end of the useful range.
const GEO_RFOOT_WIDE: [f64; 3] = [0.8, 0.9, 1.0];
/// Radial domain margin past the plate rim, as a multiple of `r_plate` (§7: room for gas to escape).
/// Held at the original value so the extended grid is directly comparable to the published rows.
const GEO_RIM_HEADROOM: f64 = 1.4;
/// Headroom multipliers for the domain-convergence check. `r_max = GEO_RIM_HEADROOM * r_plate`
/// measures escape room from the *plate* rim, so as `r_foot/R -> 1` the room past the *cloud* edge
/// shrinks from `3.7 r_foot` (at 0.3) to `0.4 r_foot` (at 1.0). A too-tight domain would confine the
/// rebound and inflate `eta_capture` exactly where the new rows are, so it is measured, not assumed.
const GEO_RIM_HEADROOM_CHECK: [f64; 3] = [1.4, 2.0, 2.8];
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
    /// The same quantity from the **confined** (plane-wave) run.
    ///
    /// Survivability is classified as `c_stag rho v^2 * focusing`, whose first factor is the 1D
    /// plane-wave stagnation pressure and whose second is the 2D ratio `P_local(concave) /
    /// P_local(flat)`. That product is the true concave peak only if a *flat* 2D plate sees the
    /// plane-wave load to begin with. Recording the plane-wave peak next to the free one makes that
    /// premise measurable rather than assumed, and because both come from this kernel the scheme
    /// error is common-mode and divides out — the same construction `eta_capture` uses (ADR-0003).
    peak_local_pressure_confined: f64,
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
    run_eta_case_with_headroom(
        d_over_d,
        l_over_d,
        r_foot_over_r,
        mach,
        cfg,
        GEO_RIM_HEADROOM,
    )
}

/// [`run_eta_case`] with the radial domain margin exposed, for the domain-convergence check.
fn run_eta_case_with_headroom(
    d_over_d: f64,
    l_over_d: f64,
    r_foot_over_r: f64,
    mach: f64,
    cfg: &GeoConfig,
    rim_headroom: f64,
) -> GeoRecord {
    let r_foot = 1.0;
    let r_plate = r_foot / r_foot_over_r;
    let length = l_over_d * 2.0 * r_foot;
    let r_max = r_plate * rim_headroom; // room past the rim for gas to escape (§7)
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
        taper_frac: 0.0,
        alpha_div: 0.0,
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
        taper_frac: 0.0,
        alpha_div: 0.0,
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
        peak_local_pressure_confined: confined.peak_local_pressure,
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
    par_map_with_progress("geometry", &cases, |&(dd, ld, rf, m)| {
        run_eta_case(dd, ld, rf, m, cfg)
    })
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
    par_map_with_progress("ablating", &cases, |&(v, rho, scale)| {
        run_ablating_case(v, rho, scale, base, base_tbl)
    })
    .into_iter()
    .flatten()
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
/// Slug length [m] the freeze-timing bracket is anchored at: the longest, most-realistic cloud
/// (`JUP_LENGTH`'s top), where the coupled headline `e_eff` is quoted — so the EOS-only frozen
/// delta translates directly onto it.
const JUP_FROZEN_LENGTH: f64 = 12.0;
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
    /// Whether the bounce reached its tail guard (`BounceResult::converged`, Q6). Same contract
    /// as the heavy-plate record: `false` means no result, not a low one.
    converged: bool,
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
        converged: result.bounce.converged,
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
    par_map_with_progress("jupiter", &cases, |&(rho, len, s)| {
        run_one_jupiter(rho, len, s, base_tbl)
    })
}

// ---- Spray plate, argon arm: radiation share of a premixed pulse (ADR-0055) ---------------------
//
// The necklace charge's target state: the PuffSat's water and `k = 10` of spray already merged. In
// the plate frame the mixture closes at `w/(1+k)` carrying the merge heat `w² k / (2(1+k)²)`, less
// what reaching the table's energy reference costs (vaporizing the spray and, for argon, atomizing
// the water). One material per run, so this is the premixed case only.

const SPRAY_K: f64 = 10.0;
/// Footprint area [m²]: half the radius of the 20 m tapered plate.
const SPRAY_FOOTPRINT: f64 = std::f64::consts::PI * 25.0;
/// (closing speed [m/s], PuffSat mass per pulse [kg]) for the 100 t plate's 8 MN·s pulse.
const SPRAY_SPEEDS: [(f64, f64); 2] = [(45_580.0, 47.0), (65_130.0, 33.0)];
/// Spray-cloud depth [m] the merged pulse occupies.
const SPRAY_DEPTH: [f64; 5] = [0.5, 1.0, 2.0, 4.0, 8.0];
const RESULT_PATH_SPRAY: &str = "data/results/water_plate/spray_radiation.jsonl";

#[derive(Debug, Clone, Copy)]
struct SprayMaterial {
    name: &'static str,
    table: &'static str,
    /// Energy [J/kg of mixture] paid between the arriving gas and the table's reference.
    charge: f64,
}

const SPRAY_MATERIALS: [SprayMaterial; 2] = [
    SprayMaterial {
        name: "argon",
        table: "data/tables/spray_k10.json",
        // Water atomization (50.94 MJ/kg, never returned in the atomic table) on 1/11 of the mass,
        // and vaporizing liquid argon at 87 K (~0.116 MJ/kg to the 0 K atom reference) on 10/11.
        charge: 50.94e6 / 11.0 + 0.116e6 * 10.0 / 11.0,
    },
    SprayMaterial {
        name: "water",
        table: "data/tables/water_jupiter.json",
        // Liquid water near 300 K sits ~1.89 MJ/kg below the table's 0 K vapor reference.
        charge: 1.89e6 * 10.0 / 11.0,
    },
];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayRecord {
    material: String,
    w: f64,
    depth: f64,
    /// PuffSat mass per footprint area [kg/m²].
    sigma_puffsat: f64,
    rho0: f64,
    t0: f64,
    e0: f64,
    /// Wall impulse per arriving PuffSat momentum, `J / (m w)`.
    beta: f64,
    /// `beta` over the overtake ceiling `1 + sqrt(1 + k)`.
    ceiling_share: f64,
    /// Radiation absorbed by the plate face over the pulse's kinetic energy `½ m w²`.
    face_radiation_share: f64,
    /// Radiation escaping to space over the same.
    escape_share: f64,
    peak_wall_pressure: f64,
    converged: bool,
}

/// The merged pulse as a hot slab closing on the plate: `(tube, config, sigma_puffsat, rho0, t0, e0)`.
fn spray_slab(
    mat: SprayMaterial,
    table: &Table,
    w: f64,
    m: f64,
    depth: f64,
) -> (Tube<TableEos>, Config, f64, f64, f64, f64) {
    let sigma_p = m / SPRAY_FOOTPRINT;
    let rho0 = sigma_p * (1.0 + SPRAY_K) / depth;
    let v0 = w / (1.0 + SPRAY_K);
    let e0 = 0.5 * w * w * SPRAY_K / (1.0 + SPRAY_K).powi(2) - mat.charge;
    let eos = TableEos::new(table.clone());
    let t0 = eos.temperature(rho0, e0);
    let cfg = Config {
        v: v0,
        length: depth,
        ..Config::production()
    };
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho0,
        v0,
        depth,
        t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    (tube, cfg, sigma_p, rho0, t0, e0)
}

fn run_one_spray(mat: SprayMaterial, table: &Table, w: f64, m: f64, depth: f64) -> SprayRecord {
    let (tube, cfg, sigma_p, rho0, t0, e0) = spray_slab(mat, table, w, m, depth);
    let result = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
    let ke = 0.5 * sigma_p * w * w;
    let beta = result.bounce.wall_impulse / result.bounce.incident_momentum;
    SprayRecord {
        material: mat.name.to_string(),
        w,
        depth,
        sigma_puffsat: sigma_p,
        rho0,
        t0,
        e0,
        beta,
        ceiling_share: beta / (1.0 + (1.0 + SPRAY_K).sqrt()),
        face_radiation_share: result.loss_radiative_wall / ke,
        escape_share: result.loss_escape_space / ke,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        converged: result.bounce.converged,
    }
}

fn cmd_spray(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases = Vec::new();
    for mat in SPRAY_MATERIALS {
        let table = Table::load(mat.table)?;
        for (w, m) in SPRAY_SPEEDS {
            for depth in SPRAY_DEPTH {
                cases.push((mat, table.clone(), w, m, depth));
            }
        }
    }
    let rows = par_map_with_progress("spray", &cases, |(mat, table, w, m, depth)| {
        run_one_spray(*mat, table, *w, *m, *depth)
    });
    if let Some(dir) = Path::new(RESULT_PATH_SPRAY).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_SPRAY, "spray", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} L={:>3.1} T0={:>7.0} -> beta={:.3} ({:.3} of ceiling) \
             face_rad={:.2e} escape={:.2e} peak_p={:.3e} converged={}",
            r.material,
            r.w,
            r.depth,
            r.t0,
            r.beta,
            r.ceiling_share,
            r.face_radiation_share,
            r.escape_share,
            r.peak_wall_pressure,
            r.converged,
        );
    })
}

// ---- Spray plate, step 1b: the same pulse on an ablating, vapor-shielded face (ADR-0055) --------
//
// The film boils at `ṁ = q_in / Q*` and its vapor forms a curtain of optical depth
// `κ_vapor · ablated mass` in front of the steel (Rung E, ADR-0014). The vapor is carried with the
// cloud's own EOS, so only its shielding and its mass are represented, not its chemistry.

/// Spray-cloud depths [m] at the working point found in step 1.
const SPRAY_ABL_DEPTH: [f64; 2] = [1.0, 2.0];
/// Effective heat of ablation [J/kg], ADR-0014's 2-10 MJ/kg band.
const SPRAY_ABL_Q_STAR: [f64; 3] = [2.0e6, 5.0e6, 10.0e6];
/// Vapor gray opacity [m²/kg]: 0 unshielded; 200 the Rung E calibration; 1500 the child repo's
/// 20%-carbon film at Bond & Bergstrom's 7.5 m²/g (550 nm); 1e4 an EUV photoabsorption estimate.
const SPRAY_ABL_KAPPA: [f64; 4] = [0.0, 200.0, 1500.0, 1.0e4];
const RESULT_PATH_SPRAY_ABL: &str = "data/results/water_plate/spray_ablating.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayAblatingRecord {
    material: String,
    w: f64,
    depth: f64,
    q_star: f64,
    kappa_vapor: f64,
    beta: f64,
    ceiling_share: f64,
    /// Radiation that still reaches the steel, over the pulse's kinetic energy.
    face_radiation_share: f64,
    escape_share: f64,
    /// Film boiled off per pulse over the whole footprint [kg].
    film_per_pulse: f64,
    /// Film per pulse over the PuffSat mass per pulse.
    film_over_puffsat: f64,
    peak_wall_pressure: f64,
    converged: bool,
}

fn run_one_spray_ablating(
    mat: SprayMaterial,
    table: &Table,
    (w, m): (f64, f64),
    depth: f64,
    q_star: f64,
    kappa_vapor: f64,
) -> SprayAblatingRecord {
    let (tube, cfg, sigma_p, ..) = spray_slab(mat, table, w, m, depth);
    let ablation = Ablation::new(q_star, cfg.t0).with_vapor_opacity(kappa_vapor);
    let result = AblatingBounce::new(tube, None, cfg.consts, cfg.limiter, ablation).run();
    let ke = 0.5 * sigma_p * w * w;
    let beta = result.bounce.wall_impulse / result.bounce.incident_momentum;
    let film = result.ablated_mass * SPRAY_FOOTPRINT;
    SprayAblatingRecord {
        material: mat.name.to_string(),
        w,
        depth,
        q_star,
        kappa_vapor,
        beta,
        ceiling_share: beta / (1.0 + (1.0 + SPRAY_K).sqrt()),
        face_radiation_share: result.loss_radiative_wall / ke,
        escape_share: result.loss_escape_space / ke,
        film_per_pulse: film,
        film_over_puffsat: film / m,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        converged: result.bounce.converged,
    }
}

fn cmd_spray_ablating(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases = Vec::new();
    for mat in SPRAY_MATERIALS {
        let table = Table::load(mat.table)?;
        for speed in SPRAY_SPEEDS {
            for depth in SPRAY_ABL_DEPTH {
                for q_star in SPRAY_ABL_Q_STAR {
                    for kappa in SPRAY_ABL_KAPPA {
                        cases.push((mat, table.clone(), speed, depth, q_star, kappa));
                    }
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-ablating",
        &cases,
        |(mat, table, speed, depth, q, kappa)| {
            run_one_spray_ablating(*mat, table, *speed, *depth, *q, *kappa)
        },
    );
    if let Some(dir) = Path::new(RESULT_PATH_SPRAY_ABL).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_SPRAY_ABL, "spray-ablating", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} L={:.0} Q*={:.0e} kv={:>6.0} -> share={:.3} face={:.2e} \
             esc={:.2e} film={:>7.1} kg ({:.2} x m) converged={}",
            r.material,
            r.w,
            r.depth,
            r.q_star,
            r.kappa_vapor,
            r.ceiling_share,
            r.face_radiation_share,
            r.escape_share,
            r.film_per_pulse,
            r.film_over_puffsat,
            r.converged,
        );
    })
}

// ---- Spray plate, step 2a: mixing quality with two materials and radiation (ADR-0055) -----------
//
// The PuffSat's water gas and the cold spray start as separate layers, ordered from the plate
// outward: spray, PuffSat, spray, PuffSat, ... with `n` pairs. `n = 1` is a stratified pulse; larger
// `n` stands in for a pulse that has partly penetrated the cloud before the layers interact.
// Lagrangian cells never exchange mass, so finer interleaving is the only mixing represented.
// The water arm's spray starts as vapor at 400 K, which flatters water by ~2% of the merge heat.

/// Number of spray/PuffSat layer pairs.
const SPRAY_MIX_PAIRS: [usize; 3] = [1, 4, 16];
/// Cells across the whole column.
const SPRAY_MIX_CELLS: usize = 384;
const RESULT_PATH_SPRAY_MIX: &str = "data/results/water_plate/spray_mixing.jsonl";

/// The face treatment: `None` is a bare cold black face; `Some(κ)` is the Rung E ablating face at
/// Q* = 5 MJ/kg with vapor opacity κ (5e3 the sourced value, 100 a curtain that has lost its soot).
const SPRAY_MIX_FACES: [Option<f64>; 3] = [None, Some(5.0e3), Some(100.0)];
const SPRAY_MIX_Q_STAR: f64 = 5.0e6;

#[derive(Debug, Clone, Copy)]
struct SprayArm {
    name: &'static str,
    spray_table: &'static str,
}

const SPRAY_ARMS: [SprayArm; 2] = [
    SprayArm {
        name: "argon",
        spray_table: "data/tables/argon.json",
    },
    SprayArm {
        name: "water",
        spray_table: "data/tables/water_jupiter.json",
    },
];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayMixRecord {
    material: String,
    w: f64,
    /// Spray-cloud depth [m].
    depth: f64,
    /// Arriving PuffSat pulse length [m].
    pulse_length: f64,
    pairs: usize,
    /// Vapor opacity of the ablating face [m²/kg]; absent for the bare face.
    kappa_vapor: Option<f64>,
    beta: f64,
    ceiling_share: f64,
    face_radiation_share: f64,
    escape_share: f64,
    film_per_pulse: f64,
    peak_wall_pressure: f64,
    converged: bool,
}

/// The layered column: spray (at rest) and PuffSat water gas (closing at `w`), each `depth` deep in
/// total, split into `pairs` alternating layers starting with spray at the plate.
fn spray_layers(
    water: &Table,
    spray: &Table,
    (w, m): (f64, f64),
    (depth_spray, depth_pulse): (f64, f64),
    pairs: usize,
    cells: usize,
) -> (Tube<TableEos>, f64) {
    spray_layers_with_gap(
        water,
        spray,
        (w, m),
        (depth_spray, depth_pulse),
        pairs,
        cells,
        0.0,
    )
}

/// Density of the near-empty standoff between plate and spray cloud [kg/m³] (2x the table floor).
const SPRAY_GAP_RHO: f64 = 2.0e-4;
/// Cells across the standoff, when there is one.
const SPRAY_GAP_CELLS: usize = 16;

/// [`spray_layers`] with a standoff of `gap` metres of near-empty spray gas between the plate and
/// the cloud. `gap = 0` builds exactly the column [`spray_layers`] always built.
fn spray_layers_with_gap(
    water: &Table,
    spray: &Table,
    speed: (f64, f64),
    depths: (f64, f64),
    pairs: usize,
    cells: usize,
    gap: f64,
) -> (Tube<TableEos>, f64) {
    spray_layers_with_gap_k(water, spray, speed, depths, pairs, cells, gap, SPRAY_K)
}

/// [`spray_layers_with_gap`] at an explicit injection ratio `k`.
#[allow(clippy::too_many_arguments)]
fn spray_layers_with_gap_k(
    water: &Table,
    spray: &Table,
    (w, m): (f64, f64),
    (depth_spray, depth_pulse): (f64, f64),
    pairs: usize,
    cells: usize,
    gap: f64,
    k: f64,
) -> (Tube<TableEos>, f64) {
    let sigma_p = m / SPRAY_FOOTPRINT;
    let layers = 2 * pairs;
    let per_layer = cells / layers;
    let (dx_s, dx_p) = (
        depth_spray / (pairs * per_layer) as f64,
        depth_pulse / (pairs * per_layer) as f64,
    );
    let t0 = Config::production().t0;
    let (rho_s, rho_p) = (sigma_p * k / depth_spray, sigma_p / depth_pulse);
    let mut positions = vec![0.0];
    let (mut mass, mut vel, mut energy, mut index) = (vec![], vec![], vec![], vec![]);
    if gap > 0.0 {
        let dx = gap / SPRAY_GAP_CELLS as f64;
        for _ in 0..SPRAY_GAP_CELLS {
            positions.push(positions[positions.len() - 1] + dx);
            mass.push(SPRAY_GAP_RHO * dx);
            vel.push(0.0);
            energy.push(spray.energy(SPRAY_GAP_RHO, t0));
            index.push(1);
        }
    }
    for layer in 0..layers {
        let (r, u, tbl, mat, dx) = if layer % 2 == 0 {
            (rho_s, 0.0, spray, 1, dx_s)
        } else {
            (rho_p, -w, water, 0, dx_p)
        };
        for _ in 0..per_layer {
            positions.push(positions[positions.len() - 1] + dx);
            mass.push(r * dx);
            vel.push(u);
            energy.push(tbl.energy(r, t0));
            index.push(mat);
        }
    }
    // Interior nodes carry the mass-weighted velocity of their two half-cells, so the column's
    // momentum is exactly the PuffSat's. The plate node is at rest; the outer node is free.
    let n = mass.len();
    let mut node_velocities = vec![0.0; n + 1];
    for i in 1..n {
        node_velocities[i] =
            (mass[i - 1] * vel[i - 1] + mass[i] * vel[i]) / (mass[i - 1] + mass[i]);
    }
    node_velocities[n] = vel[n - 1];
    let eos = TableEos::layered(
        vec![TableEos::new(water.clone()), TableEos::new(spray.clone())],
        index,
    );
    let state = LagrangianState {
        positions,
        node_velocities,
        cell_masses: mass,
        specific_energies: energy,
    };
    let tube = Tube::from_lagrangian(
        state,
        eos,
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    (tube, sigma_p)
}

fn run_one_spray_mix(
    arm: SprayArm,
    tables: (&Table, &Table),
    speed: (f64, f64),
    depths: (f64, f64),
    pairs: usize,
    face: Option<f64>,
    cells: usize,
) -> SprayMixRecord {
    let face = face.map(|kappa| (kappa, SPRAY_MIX_Q_STAR));
    run_one_spray_mix_q(arm, tables, speed, depths, pairs, face, cells)
}

/// [`run_one_spray_mix`] with the face given as `(κ_vapor, Q*)`, for films other than the 5 MJ/kg
/// Rung E default.
fn run_one_spray_mix_q(
    arm: SprayArm,
    (water, spray): (&Table, &Table),
    speed: (f64, f64),
    depths: (f64, f64),
    pairs: usize,
    face: Option<(f64, f64)>,
    cells: usize,
) -> SprayMixRecord {
    let w = speed.0;
    let (tube, sigma_p) = spray_layers(water, spray, speed, depths, pairs, cells);
    let cfg = Config::production();
    let (bounce, rad_wall, escape, film) = match face {
        None => {
            let r = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
            (r.bounce, r.loss_radiative_wall, r.loss_escape_space, 0.0)
        }
        Some((kappa, q_star)) => {
            let ablation = Ablation::new(q_star, cfg.t0).with_vapor_opacity(kappa);
            let r = AblatingBounce::new(tube, None, cfg.consts, cfg.limiter, ablation).run();
            (
                r.bounce,
                r.loss_radiative_wall,
                r.loss_escape_space,
                r.ablated_mass * SPRAY_FOOTPRINT,
            )
        }
    };
    let ke = 0.5 * sigma_p * w * w;
    let beta = bounce.wall_impulse / (sigma_p * w);
    SprayMixRecord {
        material: arm.name.to_string(),
        w,
        depth: depths.0,
        pulse_length: depths.1,
        pairs,
        kappa_vapor: face.map(|(kappa, _)| kappa),
        beta,
        ceiling_share: beta / (1.0 + (1.0 + SPRAY_K).sqrt()),
        face_radiation_share: rad_wall / ke,
        escape_share: escape / ke,
        film_per_pulse: film,
        peak_wall_pressure: bounce.peak_wall_pressure,
        converged: bounce.converged,
    }
}

fn cmd_spray_mixing(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for depth in SPRAY_ABL_DEPTH {
                for pairs in SPRAY_MIX_PAIRS {
                    for face in SPRAY_MIX_FACES {
                        cases.push((arm, spray.clone(), speed, depth, pairs, face));
                    }
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-mixing",
        &cases,
        |(arm, spray, speed, depth, pairs, face)| {
            run_one_spray_mix(
                *arm,
                (&water, spray),
                *speed,
                (*depth, *depth),
                *pairs,
                *face,
                SPRAY_MIX_CELLS,
            )
        },
    );
    if let Some(dir) = Path::new(RESULT_PATH_SPRAY_MIX).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_SPRAY_MIX, "spray-mixing", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} L={:.0} n={:>2} kv={:>6} -> share={:.3} face={:.2e} \
             esc={:.2e} film={:>6.1} kg peak_p={:.2e} converged={}",
            r.material,
            r.w,
            r.depth,
            r.pairs,
            r.kappa_vapor
                .map_or("bare".to_string(), |k| format!("{k:.0}")),
            r.ceiling_share,
            r.face_radiation_share,
            r.escape_share,
            r.film_per_pulse,
            r.peak_wall_pressure,
            r.converged,
        );
    })
}

const RESULT_PATH_SPRAY_MIX_CONV: &str = "data/results/water_plate/spray_mixing_convergence.jsonl";

/// Grid convergence of step 2a on its extreme cases: stratified and 16 pairs, both fluids, at the
/// cold speed, 1 m deep, bare face, at 192 / 384 / 768 / 1536 cells.
fn cmd_spray_mixing_convergence(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for pairs in [1, 16] {
            for cells in [192, 384, 768, 1536] {
                cases.push((arm, spray.clone(), pairs, cells));
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-mixing-convergence",
        &cases,
        |(arm, spray, pairs, cells)| {
            let mut r = run_one_spray_mix(
                *arm,
                (&water, spray),
                SPRAY_SPEEDS[0],
                (1.0, 1.0),
                *pairs,
                None,
                *cells,
            );
            r.depth = *cells as f64; // reused as the cell count in this file only
            r
        },
    );
    emit_scenario(
        RESULT_PATH_SPRAY_MIX_CONV,
        "spray-mixing-convergence",
        &rows,
        |r| {
            println!(
                "rust: {:>5} n={:>2} cells={:>5.0} -> share={:.4} face={:.3e} peak_p={:.3e}",
                r.material,
                r.pairs,
                r.depth,
                r.ceiling_share,
                r.face_radiation_share,
                r.peak_wall_pressure,
            );
        },
    )
}

// ---- Spray plate, step 2a': keeping an unmerged pulse off the plate (ADR-0055) -----------------
//
// The worst case of step 2a, stratified (one pair), on a bare face, with the spray-cloud depth and
// the arriving pulse length swept independently. The column masses are fixed by `k`, so a deeper
// cloud is more dilute and puts more distance between the contact and the plate, and a longer
// pulse is more dilute and arrives over a longer time.

const SPRAY_LEVER_DEPTHS: [f64; 3] = [1.0, 4.0, 16.0];
const SPRAY_LEVER_CELLS: usize = 768;
const RESULT_PATH_SPRAY_LEVERS: &str = "data/results/water_plate/spray_levers.jsonl";

fn cmd_spray_levers(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for depth_spray in SPRAY_LEVER_DEPTHS {
                for depth_pulse in SPRAY_LEVER_DEPTHS {
                    cases.push((arm, spray.clone(), speed, (depth_spray, depth_pulse)));
                }
            }
        }
    }
    let rows = par_map_with_progress("spray-levers", &cases, |(arm, spray, speed, depths)| {
        run_one_spray_mix(
            *arm,
            (&water, spray),
            *speed,
            *depths,
            1,
            None,
            SPRAY_LEVER_CELLS,
        )
    });
    emit_scenario(RESULT_PATH_SPRAY_LEVERS, "spray-levers", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} cloud={:>4.0} pulse={:>4.0} -> share={:.3} face={:.2e} \
             esc={:.2e} peak_p={:.2e} converged={}",
            r.material,
            r.w,
            r.depth,
            r.pulse_length,
            r.ceiling_share,
            r.face_radiation_share,
            r.escape_share,
            r.peak_wall_pressure,
            r.converged,
        );
    })
}

const RESULT_PATH_SPRAY_LEVERS_SHIELDED: &str =
    "data/results/water_plate/spray_levers_shielded.jsonl";

/// The candidate design point (16 m cloud, 1 m and 4 m pulses) on the ablating faces, for its film
/// cost per pulse.
fn cmd_spray_levers_shielded(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for pulse in [1.0, 4.0] {
                for face in [Some(5.0e3), Some(100.0)] {
                    cases.push((arm, spray.clone(), speed, pulse, face));
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-levers-shielded",
        &cases,
        |(arm, spray, speed, pulse, face)| {
            run_one_spray_mix(
                *arm,
                (&water, spray),
                *speed,
                (16.0, *pulse),
                1,
                *face,
                SPRAY_LEVER_CELLS,
            )
        },
    );
    emit_scenario(
        RESULT_PATH_SPRAY_LEVERS_SHIELDED,
        "spray-levers-shielded",
        &rows,
        |r| {
            println!(
                "rust: {:>5} w={:>5.0} pulse={:>2.0} kv={:>5} -> share={:.3} face={:.2e} film={:>6.1} kg peak_p={:.2e}",
                r.material,
                r.w,
                r.pulse_length,
                r.kappa_vapor
                    .map_or("bare".to_string(), |k| format!("{k:.0}")),
                r.ceiling_share,
                r.face_radiation_share,
                r.film_per_pulse,
                r.peak_wall_pressure,
            );
        },
    )
}

const RESULT_PATH_SPRAY_LEVERS_CONV: &str =
    "data/results/water_plate/spray_levers_convergence.jsonl";

/// Grid convergence of step 2a' at its candidate design point: a 16 m cloud, 1 m and 4 m pulses,
/// both fluids, both speeds, at 384 / 768 / 1536 cells.
fn cmd_spray_levers_convergence(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for pulse in [1.0, 4.0] {
                for cells in [384, 768, 1536] {
                    cases.push((arm, spray.clone(), speed, pulse, cells));
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-levers-convergence",
        &cases,
        |(arm, spray, speed, pulse, cells)| {
            let mut r = run_one_spray_mix(
                *arm,
                (&water, spray),
                *speed,
                (16.0, *pulse),
                1,
                None,
                *cells,
            );
            r.pairs = *cells; // reused as the cell count in this file only
            r
        },
    );
    emit_scenario(
        RESULT_PATH_SPRAY_LEVERS_CONV,
        "spray-levers-convergence",
        &rows,
        |r| {
            println!(
                "rust: {:>5} w={:>5.0} pulse={:>2.0} cells={:>4} -> share={:.4} face={:.3e} peak_p={:.3e}",
                r.material,
                r.w,
                r.pulse_length,
                r.pairs,
                r.ceiling_share,
                r.face_radiation_share,
                r.peak_wall_pressure,
            );
        },
    )
}

// ---- Spray plate, step 2b: rim spill and lateral relief in 2-D (ADR-0055) ----------------------
//
// The ADR-0003/0008 geometry method: the resting spray cloud and the arriving slug as one
// effective-γ gas, run free (spilling past the 10 m plate) and confined (plane wave). Their ratios
// give `eta_capture` (impulse kept) and the lateral relief of the peak facesheet pressure. Both
// multiply the 1-D real-physics column of step 2a' at the same depths. Lengths are in units of the
// 4 m pulse: footprint radius 1.25 (5 m), plate radius 2.5 (10 m), domain radius 3.5.

/// Spray-cloud depths [m].
const SPRAY_2D_DEPTHS_M: [usize; 4] = [4, 8, 12, 16];
const SPRAY_2D_PULSE_M: f64 = 4.0;
const SPRAY_2D_GAMMA: f64 = 1.4;
const SPRAY_2D_MACH: f64 = 20.0;
/// Radial cells across the 3.5-unit domain at production resolution; axial cells per unit length.
const SPRAY_2D_NR: usize = 112;
const SPRAY_2D_NZ_PER_UNIT: usize = 32;
const RESULT_PATH_SPRAY_2D: &str = "data/results/water_plate/spray_2d.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct Spray2dRecord {
    depth: f64,
    nr: usize,
    nz: usize,
    eta_capture: f64,
    /// Peak facesheet pressure, free over confined: the lateral relief factor.
    relief: f64,
    ratio_free: f64,
    ratio_confined: f64,
}

fn spray_2d_case(depth_m: usize, refine: usize) -> Spray2dRecord {
    let unit = SPRAY_2D_PULSE_M;
    let depth = depth_m as f64 / unit;
    // The column (depth + 1 pulse length) plus 2 units of headroom for the rebound; depth_m is a
    // multiple of 4, so this is an integer number of units.
    let nz = refine * SPRAY_2D_NZ_PER_UNIT * (depth_m / 4 + 3);
    let z_max = (depth_m / 4 + 3) as f64;
    let spray = SprayCloud {
        depth,
        mass_ratio: SPRAY_K,
        standoff: 0.0,
    };
    let base = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate: 10.0 / unit,
        r_max: 14.0 / unit,
        z_max,
        nr: refine * SPRAY_2D_NR,
        nz,
        confined: false,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let free = run_spray_bounce(&base, spray);
    let confined = run_spray_bounce(
        &SlugConfig {
            r_foot: base.r_max,
            r_plate: base.r_max,
            nr: 8,
            confined: true,
            ..base
        },
        spray,
    );
    Spray2dRecord {
        depth: depth_m as f64,
        nr: base.nr,
        nz,
        eta_capture: eta_capture(&free, &confined),
        relief: free.peak_local_pressure / confined.peak_local_pressure,
        ratio_free: free.restitution_ratio(),
        ratio_confined: confined.restitution_ratio(),
    }
}

fn cmd_spray_2d(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases: Vec<(usize, usize)> = SPRAY_2D_DEPTHS_M.iter().map(|&d| (d, 1)).collect();
    cases.extend([(4, 2), (16, 2)]);
    let geometry =
        par_map_with_progress("spray-2d", &cases, |&(d, refine)| spray_2d_case(d, refine));
    emit_scenario(RESULT_PATH_SPRAY_2D, "spray-2d", &geometry, |r| {
        println!(
            "rust: 2-D cloud={:>2.0} m nr={:>3} nz={:>3} -> eta_capture={:.4} relief={:.3} \
             (free {:.4}, confined {:.4})",
            r.depth, r.nr, r.nz, r.eta_capture, r.relief, r.ratio_free, r.ratio_confined,
        );
    })?;

    // The 1-D real-physics column at the same depths, bare and on the sourced shielded face.
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases_1d = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for &d in &SPRAY_2D_DEPTHS_M {
                for face in [None, Some(5.0e3)] {
                    cases_1d.push((arm, spray.clone(), speed, d as f64, face));
                }
            }
        }
    }
    let rows_1d = par_map_with_progress(
        "spray-2d-columns",
        &cases_1d,
        |(arm, spray, speed, depth, face)| {
            run_one_spray_mix(
                *arm,
                (&water, spray),
                *speed,
                (*depth, SPRAY_2D_PULSE_M),
                1,
                *face,
                SPRAY_LEVER_CELLS,
            )
        },
    );
    emit_scenario(
        "data/results/water_plate/spray_2d_columns.jsonl",
        "spray-2d-columns",
        &rows_1d,
        |r| {
            let g = geometry
                .iter()
                .find(|g| g.depth == r.depth && g.nr == SPRAY_2D_NR)
                .expect("a production 2-D row at every depth");
            println!(
                "rust: {:>5} w={:>5.0} cloud={:>2.0} kv={:>5} -> 1-D share={:.3} peak={:.2e} \
                 film={:>5.1} kg | 2-D share={:.3} eta_jet={:.3} peak={:.2e}",
                r.material,
                r.w,
                r.depth,
                r.kappa_vapor
                    .map_or("bare".to_string(), |k| format!("{k:.0}")),
                r.ceiling_share,
                r.peak_wall_pressure,
                r.film_per_pulse,
                r.ceiling_share * g.eta_capture,
                (r.ceiling_share * g.eta_capture * (1.0 + (1.0 + SPRAY_K).sqrt()) - 1.0)
                    / (1.0 + SPRAY_K).sqrt(),
                r.peak_wall_pressure * g.relief,
            );
        },
    )
}

// ---- Spray plate, step 2c: merge just before the plate (ADR-0055) ------------------------------
//
// The spray cloud stands off the plate by a near-empty gap, so the PuffSat merges inside the cloud
// and the merged gas crosses the gap before it presses on the plate. Stratified pulse, bare face,
// 768 cells. Each row also reports how long the face stays above each candidate allowable, from
// the wall-pressure history.

// Gaps > 0 fail numerically here: the near-empty gap cells are crushed far off the table and
// the time step collapses. A gap needs the 2-D Eulerian kernel; this sweep keeps gap = 0 for the
// spike widths.
const SPRAY_STANDOFF_GAPS: [f64; 1] = [0.0];
const SPRAY_STANDOFF_DEPTHS: [f64; 3] = [1.0, 2.0, 4.0];
const SPRAY_STANDOFF_PULSE: f64 = 4.0;
/// Candidate face allowables [Pa] for the time-above-limit columns.
const SPRAY_ALLOWABLES: [f64; 4] = [0.4e9, 0.7e9, 1.0e9, 1.5e9];
const RESULT_PATH_SPRAY_STANDOFF: &str = "data/results/water_plate/spray_standoff.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayStandoffRecord {
    material: String,
    w: f64,
    gap: f64,
    depth: f64,
    ceiling_share: f64,
    face_radiation_share: f64,
    escape_share: f64,
    peak_wall_pressure: f64,
    /// Seconds the face pressure exceeds each of [`SPRAY_ALLOWABLES`].
    time_above: [f64; 4],
    converged: bool,
}

/// Total time a sampled `(t, p)` history spends above `level` (trapezoid over sample intervals).
fn time_above(history: &[(f64, f64)], level: f64) -> f64 {
    history
        .windows(2)
        .filter(|w| 0.5 * (w[0].1 + w[1].1) > level)
        .map(|w| w[1].0 - w[0].0)
        .sum()
}

fn run_one_spray_standoff(
    arm: SprayArm,
    (water, spray): (&Table, &Table),
    speed: (f64, f64),
    depth: f64,
    gap: f64,
) -> SprayStandoffRecord {
    let w = speed.0;
    let (tube, sigma_p) = spray_layers_with_gap(
        water,
        spray,
        speed,
        (depth, SPRAY_STANDOFF_PULSE),
        1,
        SPRAY_LEVER_CELLS,
        gap,
    );
    let cfg = Config::production();
    let (r, history) = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run_with_history();
    let ke = 0.5 * sigma_p * w * w;
    let beta = r.bounce.wall_impulse / (sigma_p * w);
    SprayStandoffRecord {
        material: arm.name.to_string(),
        w,
        gap,
        depth,
        ceiling_share: beta / (1.0 + (1.0 + SPRAY_K).sqrt()),
        face_radiation_share: r.loss_radiative_wall / ke,
        escape_share: r.loss_escape_space / ke,
        peak_wall_pressure: r.bounce.peak_wall_pressure,
        time_above: SPRAY_ALLOWABLES.map(|level| time_above(&history, level)),
        converged: r.bounce.converged,
    }
}

fn cmd_spray_standoff(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for depth in SPRAY_STANDOFF_DEPTHS {
                for gap in SPRAY_STANDOFF_GAPS {
                    cases.push((arm, spray.clone(), speed, depth, gap));
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-standoff",
        &cases,
        |(arm, spray, speed, depth, gap)| {
            run_one_spray_standoff(*arm, (&water, spray), *speed, *depth, *gap)
        },
    );
    if let Some(dir) = Path::new(RESULT_PATH_SPRAY_STANDOFF).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_SPRAY_STANDOFF, "spray-standoff", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} cloud={:.0} gap={:>3.1} -> share={:.3} face={:.2e} \
             peak={:.2e} t>0.4/0.7/1.0/1.5GPa={:.0}/{:.0}/{:.0}/{:.0} us converged={}",
            r.material,
            r.w,
            r.depth,
            r.gap,
            r.ceiling_share,
            r.face_radiation_share,
            r.peak_wall_pressure,
            r.time_above[0] * 1e6,
            r.time_above[1] * 1e6,
            r.time_above[2] * 1e6,
            r.time_above[3] * 1e6,
            r.converged,
        );
    })
}

// ---- Walled nozzle: the blast inside the 20 m^3 near-term chamber (wall spike and rise) -----
//
// A 1-D spherical Euler run (`euler2d::sphere`) of the rod's energy released at the chamber's
// centre. It records the wall-pressure history: the reflected spike, the rise to the
// quasi-static pressure, and the reverberations. Fills: a dry uniform charge, or the parent's
// liquid share as an equilibrium heavy-gas shell near the wall (drops of 0.5-1 mm shatter and
// follow the gas within microseconds). Delivery: one charge, two halves T_b/2 apart, or the
// 1/4-1/2-1/4 train. The gas is γ-law, so pressures are read against the run's own
// quasi-static pressure and scaled to the solved 496 / 818 bar in Python
// (`puffsat.walled_nozzle.chamber_blast`).

const BLAST_RADIUS: f64 = 1.684;
const BLAST_CELLS: usize = 1684;
const BLAST_ENERGY: f64 = 0.5 * 2.5 * 75_000.0 * 75_000.0;
const BLAST_ROD_PLUG_KG: f64 = 6.9;
const BLAST_DEPOSIT_RADIUS: f64 = 0.15;
const BLAST_FILL_PA: f64 = 5.0e5;
const BLAST_T_END: f64 = 6.0e-3;
const BLAST_BIN: f64 = 0.5e-6;
/// Half the breathing period of the parent's layups (1.10-1.11 ms, `chamber_fatigue.py`).
const BLAST_HALF_PERIOD: f64 = 0.553e-3;
const RESULT_PATH_CHAMBER_BLAST: &str = "data/results/walled_nozzle/near_term/chamber_blast.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct ChamberBlastRecord {
    chamber: String,
    gamma: f64,
    fill: String,
    liquid_fraction: f64,
    shell_m: f64,
    graded: bool,
    delivery: String,
    peak_pa: f64,
    first_arrival_s: f64,
    /// Bin-max wall pressure every `BLAST_BIN`.
    time_s: Vec<f64>,
    pressure_pa: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
struct BlastCase {
    chamber: &'static str,
    charge_kg: f64,
    gamma: f64,
    liquid_fraction: f64,
    shell_m: f64,
    graded: bool,
    delivery: &'static str,
}

fn chamber_blast_one(c: BlastCase) -> ChamberBlastRecord {
    let volume = 4.0 / 3.0 * std::f64::consts::PI * BLAST_RADIUS.powi(3);
    let gas_rho = c.charge_kg * (1.0 - c.liquid_fraction) / volume;
    let inner = BLAST_RADIUS - c.shell_m;
    let shell_volume = volume - 4.0 / 3.0 * std::f64::consts::PI * inner.max(0.0).powi(3);
    let liquid_kg = c.charge_kg * c.liquid_fraction;
    // Graded: liquid density rises linearly from zero at the shell's inner edge to the wall.
    let graded_norm = if c.graded && c.shell_m > 0.0 {
        let n = 400;
        let dr = c.shell_m / f64::from(n);
        (0..n)
            .map(|k| {
                let r = inner + (f64::from(k) + 0.5) * dr;
                4.0 * std::f64::consts::PI * r * r * dr * (r - inner) / c.shell_m
            })
            .sum::<f64>()
    } else {
        1.0
    };
    let fill = move |r: f64| {
        let liquid = if c.liquid_fraction > 0.0 && r >= inner {
            if c.graded {
                liquid_kg * (r - inner) / c.shell_m / graded_norm
            } else {
                liquid_kg / shell_volume
            }
        } else {
            0.0
        };
        (gas_rho + liquid, 0.0, BLAST_FILL_PA)
    };
    let mut sphere = euler2d::sphere::Sphere1D::new(BLAST_CELLS, BLAST_RADIUS, c.gamma, fill);
    let pieces: &[(f64, f64)] = match c.delivery {
        "split 2" => &[(0.5, 0.0), (0.5, BLAST_HALF_PERIOD)],
        "train 3" => &[
            (0.25, 0.0),
            (0.5, BLAST_HALF_PERIOD),
            (0.25, 2.0 * BLAST_HALF_PERIOD),
        ],
        _ => &[(1.0, 0.0)],
    };
    let mut bins: Vec<f64> = Vec::new();
    let mut first_arrival = f64::NAN;
    let mut t = 0.0;
    for (k, &(share, at)) in pieces.iter().enumerate() {
        if at > t {
            t = sphere.run_to(t, at, |time, p| {
                record_bin(&mut bins, &mut first_arrival, time, p)
            });
        }
        sphere.deposit(
            BLAST_ENERGY * share,
            BLAST_ROD_PLUG_KG * share,
            BLAST_DEPOSIT_RADIUS,
        );
        let next = pieces.get(k + 1).map_or(BLAST_T_END, |p| p.1);
        t = sphere.run_to(t, next, |time, p| {
            record_bin(&mut bins, &mut first_arrival, time, p)
        });
    }
    let time_s = (0..bins.len())
        .map(|i| (i as f64 + 0.5) * BLAST_BIN)
        .collect();
    ChamberBlastRecord {
        chamber: c.chamber.to_string(),
        gamma: c.gamma,
        fill: if c.liquid_fraction > 0.0 {
            "liquid shell"
        } else {
            "dry"
        }
        .to_string(),
        liquid_fraction: c.liquid_fraction,
        shell_m: c.shell_m,
        graded: c.graded,
        delivery: c.delivery.to_string(),
        peak_pa: bins.iter().copied().fold(0.0, f64::max),
        first_arrival_s: first_arrival,
        time_s,
        pressure_pa: bins,
    }
}

/// Keep the maximum wall pressure in each `BLAST_BIN`, and the first time it doubles the fill.
fn record_bin(bins: &mut Vec<f64>, first: &mut f64, time: f64, p: f64) {
    // SAFE: time is non-negative and bounded by BLAST_T_END, so the bin index is small.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let i = (time / BLAST_BIN) as usize;
    if bins.len() <= i {
        bins.resize(i + 1, 0.0);
    }
    bins[i] = bins[i].max(p);
    if first.is_nan() && p > 2.0 * BLAST_FILL_PA {
        *first = time;
    }
}

fn cmd_chamber_blast(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases = Vec::new();
    for (chamber, charge_kg) in [("methane 7000 K", 68.763), ("hydrogen 5500 K", 59.913)] {
        for gamma in [1.2, 1.4] {
            let fills: [(f64, f64, bool); 6] = [
                (0.0, 0.0, false),
                (0.4, 0.4, false),
                (0.4, 0.8, false),
                (0.7, 0.4, false),
                (0.7, 0.8, false),
                (0.7, 0.8, true),
            ];
            for (liquid_fraction, shell_m, graded) in fills {
                for delivery in ["single", "split 2", "train 3"] {
                    cases.push(BlastCase {
                        chamber,
                        charge_kg,
                        gamma,
                        liquid_fraction,
                        shell_m,
                        graded,
                        delivery,
                    });
                }
            }
        }
    }
    let rows = par_map_with_progress("chamber-blast", &cases, |c| chamber_blast_one(*c));
    if let Some(dir) = Path::new(RESULT_PATH_CHAMBER_BLAST).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_CHAMBER_BLAST, "chamber-blast", &rows, |r| {
        println!(
            "rust: {:>15} g={:.1} {:>12} liq={:.1} shell={:.1} graded={} {:>8} -> peak={:.2e} Pa \
             arrival={:.0} us",
            r.chamber,
            r.gamma,
            r.fill,
            r.liquid_fraction,
            r.shell_m,
            r.graded,
            r.delivery,
            r.peak_pa,
            r.first_arrival_s * 1e6,
        );
    })
}

// ---- Walled nozzle: the blast in a rocket-shaped 20 m^3 chamber (2-D axisymmetric) ----------
//
// The sphere is the worst case: the blast reaches the whole wall at once and refocuses at the
// centre. A rocket chamber spreads the arrival over its wall. The question is whether its
// convergence focuses the shock instead. Four 20 m^3 contours (`euler2d::vessel`), each with the
// rod's energy released one port-face radius inside the port. The wall pressure is recorded along
// the contour and the port face.

const VESSEL_THROAT_R: f64 = 0.2178; // sqrt(0.149 m^2 / pi)
const VESSEL_VOLUME: f64 = 20.0;
/// Port opening radius [m]: the 15 cm port of the parent's `fig:rod_port_plug`.
const VESSEL_PORT_R: f64 = 0.075;
const VESSEL_T_END: f64 = 4.0e-3;
const VESSEL_BIN: f64 = 2.0e-6;
const RESULT_PATH_VESSEL_BLAST: &str = "data/results/walled_nozzle/near_term/vessel_blast.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct VesselBlastRecord {
    shape: String,
    deposition: String,
    chamber: String,
    gamma: f64,
    cell_m: f64,
    r_c: f64,
    z_c: f64,
    z_hi: f64,
    curve: f64,
    #[serde(default)]
    h_head: f64,
    /// Chamber volume [m^3], bulge included.
    #[serde(default)]
    volume_m3: f64,
    deposit_z: f64,
    /// Station positions along the wall: side wall rows, then the port face.
    station_z: Vec<f64>,
    station_r: Vec<f64>,
    time_s: Vec<f64>,
    /// `pressure_pa[bin][station]`, the bin maximum.
    pressure_pa: Vec<Vec<f64>>,
    /// Mass [kg] and convected axial momentum [N s] that left through the throat.
    #[serde(default)]
    throat_mass_out: f64,
    #[serde(default)]
    throat_momentum_out: f64,
}

/// Solve the free length (cylinder end `z_c`, or `z_hi` with no cylinder) for 20 m^3.
/// [`vessel_for_volume`] at another volume and throat radius.
fn vessel_sized(
    r_c: f64,
    nose: Option<f64>,
    curve: f64,
    head: f64,
    volume: f64,
    r_t: f64,
) -> euler2d::vessel::Vessel {
    let make = |len: f64| {
        let (z_c, z_hi) = match nose {
            Some(n) => (head + len, head + len + n),
            None => (0.0, len),
        };
        euler2d::vessel::Vessel {
            r_c,
            z_c,
            z_hi,
            r_t,
            curve,
            h_head: head,
            r_port: if head > 0.0 { VESSEL_PORT_R } else { 0.0 },
            bulge: euler2d::vessel::Bulge::NONE,
        }
    };
    let (mut lo, mut hi) = (0.0, 40.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if make(mid).volume() < volume {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    make(0.5 * (lo + hi))
}

fn vessel_for_volume(
    r_c: f64,
    nose: Option<f64>,
    curve: f64,
    head: f64,
) -> euler2d::vessel::Vessel {
    let make = |len: f64| {
        let (z_c, z_hi) = match nose {
            Some(n) => (head + len, head + len + n),
            None => (0.0, len),
        };
        euler2d::vessel::Vessel {
            r_c,
            z_c,
            z_hi,
            r_t: VESSEL_THROAT_R,
            curve,
            h_head: head,
            r_port: if head > 0.0 { VESSEL_PORT_R } else { 0.0 },
            bulge: euler2d::vessel::Bulge::NONE,
        }
    };
    let (mut lo, mut hi) = (0.0, 20.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if make(mid).volume() < VESSEL_VOLUME {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    make(0.5 * (lo + hi))
}

#[allow(clippy::too_many_lines)]
/// One piece of the rod's energy, or of cold added mass: a sphere or spherical shell
/// (`z0 == z1`, between `inner` and `radius`) or an axial cylinder or disk, carrying heat [J],
/// mass [kg] and an axial velocity [m/s].
#[derive(Debug, Clone, Copy)]
struct DepositSeg {
    z0: f64,
    z1: f64,
    radius: f64,
    inner: f64,
    heat: f64,
    mass: f64,
    uz: f64,
}

impl DepositSeg {
    fn contains(&self, z: f64, r: f64) -> bool {
        if self.z1 - self.z0 < 1e-9 {
            let d = (z - self.z0).hypot(r);
            d >= self.inner && d < self.radius
        } else {
            r < self.radius && z >= self.z0 && z <= self.z1
        }
    }
}

#[derive(Debug, Clone)]
struct VesselCase {
    shape: &'static str,
    vessel: euler2d::vessel::Vessel,
    chamber: &'static str,
    charge_kg: f64,
    gamma: f64,
    cell: f64,
    t_end: f64,
    deposition: &'static str,
    segs: Vec<DepositSeg>,
    /// A rod entering through the port: `(radius, density, speed, duration)`. The port cells
    /// inside `radius` are held at the rod's state while it passes.
    inflow: Option<(f64, f64, f64, f64)>,
    /// A dense layer of the charge: `fraction` of it in the band, the rest uniform elsewhere, all
    /// at one temperature (so the band's pressure is higher). Held by membranes that burst
    /// instantly when the blast arrives, or standing in for a fine liquid mist (small droplets
    /// follow the gas within millimetres, so a dense gas is their equilibrium limit).
    fill_band: Option<FillBand>,
}

/// Where a dense band of the charge sits: `z0 <= z < z1` and `r0 <= r < r1`.
#[derive(Debug, Clone, Copy)]
struct FillBand {
    z0: f64,
    z1: f64,
    r0: f64,
    r1: f64,
    fraction: f64,
}

impl FillBand {
    /// A layer across the whole chamber section.
    fn slab(z0: f64, z1: f64, fraction: f64) -> Self {
        Self {
            z0,
            z1,
            r0: 0.0,
            r1: f64::INFINITY,
            fraction,
        }
    }

    fn contains(&self, z: f64, r: f64) -> bool {
        z >= self.z0 && z < self.z1 && r >= self.r0 && r < self.r1
    }
}

/// The point blast one port-face radius inside the port, all heat (the original setup).
fn point_deposit(vessel: &euler2d::vessel::Vessel) -> Vec<DepositSeg> {
    let z = vessel.r_c.min(0.5 * vessel.z_hi);
    vec![DepositSeg {
        z0: z,
        z1: z,
        radius: BLAST_DEPOSIT_RADIUS,
        inner: 0.0,
        heat: BLAST_ENERGY,
        mass: BLAST_ROD_PLUG_KG,
        uz: 0.0,
    }]
}

#[allow(clippy::too_many_lines)]
fn vessel_blast_one(c: &VesselCase) -> VesselBlastRecord {
    use euler2d::kernel::{Bc, Grid2D};
    use euler2d::state::Prim;
    let (vessel, cell, gamma) = (c.vessel, c.cell, c.gamma);
    // SAFE: chamber lengths and the cell size are positive and small, so the counts fit a usize.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let (nz, nr) = (
        (vessel.z_hi / cell).ceil() as usize,
        (vessel.r_max() / cell).ceil() as usize + 3,
    );
    let mut g = Grid2D::new(nz, nr, cell, cell, gamma);
    g.set_axisymmetric(true);
    g.bc_zlo = Bc::Reflect;
    g.bc_zhi = Bc::Transmissive;
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect;
    g.set_vessel(Some(vessel));
    g.set_parallel(true); // bit-identical to serial (euler2d vessel test)
    let centre = |iz: usize, ir: usize| ((iz as f64 + 0.5) * cell, (ir as f64 + 0.5) * cell);
    // Fill: uniform, or a dense band held by membranes. Volumes are summed over the fluid cells so
    // the charge is exact on the grid.
    let (v_band, v_fluid) = (0..nz)
        .flat_map(|iz| (0..nr).map(move |ir| (iz, ir)))
        .filter(|&(iz, ir)| {
            let (z, r) = centre(iz, ir);
            !vessel.is_solid(z, r)
        })
        .fold((0.0, 0.0), |(b, f), (iz, ir)| {
            let (z, r) = centre(iz, ir);
            let v = 2.0 * std::f64::consts::PI * (ir as f64 + 0.5) * cell * cell * cell;
            let in_band = c.fill_band.is_some_and(|band| band.contains(z, r));
            (if in_band { b + v } else { b }, f + v)
        });
    let uniform = c.charge_kg / v_fluid;
    let fill_rho = |z: f64, r: f64| match c.fill_band {
        Some(band) if band.contains(z, r) => band.fraction * c.charge_kg / v_band,
        Some(band) => (1.0 - band.fraction) * c.charge_kg / (v_fluid - v_band),
        None => uniform,
    };
    // Each segment is spread over the cells whose centres it contains, normalised by their actual
    // volume, so heat, mass and momentum are exact whatever the grid.
    let volumes: Vec<f64> = c
        .segs
        .iter()
        .map(|seg| {
            (0..nz)
                .flat_map(|iz| (0..nr).map(move |ir| (iz, ir)))
                .filter(|&(iz, ir)| {
                    let (z, r) = centre(iz, ir);
                    seg.contains(z, r)
                })
                .map(|(_, ir)| 2.0 * std::f64::consts::PI * (ir as f64 + 0.5) * cell * cell * cell)
                .sum()
        })
        .collect();
    g.init(|iz, ir| {
        let (z, r) = centre(iz, ir);
        let base = fill_rho(z, r);
        let (mut rho, mut mom, mut heat) = (base, 0.0, 0.0);
        for (seg, &v) in c.segs.iter().zip(&volumes) {
            if v > 0.0 && seg.contains(z, r) {
                rho += seg.mass / v;
                mom += seg.mass * seg.uz / v;
                heat += seg.heat / v;
            }
        }
        Prim::new(
            rho,
            mom / rho,
            0.0,
            BLAST_FILL_PA * base / uniform + (gamma - 1.0) * heat,
        )
    });
    let stations = g.vessel_wall_pressures();
    let keep: Vec<usize> = (0..stations.len()).step_by(2).collect();
    let mut bins: Vec<Vec<f64>> = Vec::new();
    let mut t = 0.0;
    let (mut throat_mass, mut throat_momentum) = (0.0, 0.0);
    let t_end = c.t_end;
    while t < t_end {
        if let Some((radius, rho, speed, until)) = c.inflow {
            if t < until {
                let rod = euler2d::state::Prim::new(rho, speed, 0.0, BLAST_FILL_PA);
                for ir in 0..nr {
                    if (ir as f64 + 0.5) * cell < radius {
                        g.set_prim(0, ir, rod);
                    }
                }
            }
        }
        let dt = g.stable_dt().min(t_end - t);
        g.step(dt);
        t += dt;
        // SAFE: t is non-negative and bounded by t_end.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let b = (t / VESSEL_BIN) as usize;
        if bins.len() <= b {
            bins.resize(b + 1, vec![0.0; keep.len()]);
        }
        let walls = g.vessel_wall_pressures();
        for (j, &k) in keep.iter().enumerate() {
            bins[b][j] = bins[b][j].max(walls[k].2);
        }
        // Outflow through the throat plane (the last row, inside the throat radius).
        let r_t = vessel.r_wall(vessel.z_hi);
        for ir in 0..nr {
            let r = (ir as f64 + 0.5) * cell;
            if r >= r_t {
                break;
            }
            let w = g.prim(nz - 1, ir);
            if w.uz > 0.0 {
                let area = 2.0 * std::f64::consts::PI * r * cell;
                throat_mass += w.rho * w.uz * area * dt;
                throat_momentum += w.rho * w.uz * w.uz * area * dt;
            }
        }
    }
    VesselBlastRecord {
        shape: c.shape.to_string(),
        deposition: c.deposition.to_string(),
        chamber: c.chamber.to_string(),
        gamma,
        cell_m: cell,
        r_c: vessel.r_c,
        z_c: vessel.z_c,
        z_hi: vessel.z_hi,
        curve: vessel.curve,
        h_head: vessel.h_head,
        volume_m3: vessel.volume(),
        deposit_z: c.segs.first().map_or(0.0, |d| d.z0),
        station_z: keep.iter().map(|&k| stations[k].0).collect(),
        station_r: keep.iter().map(|&k| stations[k].1).collect(),
        time_s: (0..bins.len())
            .map(|i| (i as f64 + 0.5) * VESSEL_BIN)
            .collect(),
        pressure_pa: bins,
        throat_mass_out: throat_mass,
        throat_momentum_out: throat_momentum,
    }
}

/// The deposition geometries (parent `sec:material_chamber_plug`), all at the rod's full energy:
/// the 0.98 m foam plug with its front face one port-face radius inside the port, 64% of the
/// energy released there as heat and 36% left as the merged 6.9 kg moving inward at 27 km/s, and
/// ~15% of the rod eroded crossing the gas before the plug.
fn deposition_cases(v: euler2d::vessel::Vessel) -> Vec<(&'static str, Vec<DepositSeg>)> {
    let e = BLAST_ENERGY;
    let m = BLAST_ROD_PLUG_KG;
    let front = v.r_c;
    let plug = (front, front + 0.98);
    let long = (0.5, v.z_c + 0.5 * (v.z_hi - v.z_c));
    let slug_v = 27_000.0;
    let slug_ke = 0.5 * m * slug_v * slug_v;
    let line = |z: (f64, f64), heat: f64, mass: f64, uz: f64| DepositSeg {
        z0: z.0,
        z1: z.1,
        radius: 0.10,
        inner: 0.0,
        heat,
        mass,
        uz,
    };
    // Cold polyethylene added near the plug as a sacrificial tamper: a 2 cm disk of 0.4 m radius
    // 10 cm on the port side of the plug's front face, or a shell 0.30-0.35 m around that face.
    let disk = |kg: f64| DepositSeg {
        z0: front - 0.12,
        z1: front - 0.10,
        radius: 0.40,
        inner: 0.0,
        heat: 0.0,
        mass: kg,
        uz: 0.0,
    };
    let shell = |kg: f64| DepositSeg {
        z0: front,
        z1: front,
        radius: 0.35,
        inner: 0.30,
        heat: 0.0,
        mass: kg,
        uz: 0.0,
    };
    let tamped = |extra: DepositSeg| vec![line(plug, e - slug_ke, m, slug_v), extra];
    vec![
        ("point", point_deposit(&v)),
        ("plug line, heat", vec![line(plug, e, m, 0.0)]),
        ("plug line + slug", vec![line(plug, e - slug_ke, m, slug_v)]),
        (
            "track + plug + slug",
            vec![
                line((0.2, front), 0.15 * e, 0.0, 0.0),
                line(plug, 0.85 * e - slug_ke, m, slug_v),
            ],
        ),
        ("long column, heat", vec![line(long, e, m, 0.0)]),
        (
            "long column + slug",
            vec![line(long, e - slug_ke, m, slug_v)],
        ),
        ("plug + slug + 2 kg port disk", tamped(disk(2.0))),
        ("plug + slug + 5 kg port disk", tamped(disk(5.0))),
        ("plug + slug + 3 kg shell", tamped(shell(3.0))),
        ("plug + slug + 10 kg shell", tamped(shell(10.0))),
    ]
}

/// A long plug of total axis mass `total_kg` (the paper's 4.4 kg plug included) over `length`
/// from `front`, with uniform or back-weighted density (3x denser at the back). The rod accretes
/// it slice by slice in sticky collisions: with `p = m_rod w` fixed, each slice releases
/// `p^2/2 (1/m_before - 1/m_after)` as heat where it sits. The final merged mass moves on at
/// `p / (m_rod + total)` from the column's last slice. Energy and momentum are exact; the hot
/// accreted mass is placed where it was struck rather than carried with the moving body.
fn long_plug_segments(
    front: f64,
    length: f64,
    total_kg: f64,
    back_weighted: bool,
) -> Vec<DepositSeg> {
    long_plug_segments_with(front, length, total_kg, back_weighted, false)
}

/// [`long_plug_segments`], optionally with a **hot remnant**: the heat of the collisions in the
/// back half of the column travels with the merged body instead of staying where it was struck,
/// a bound on how much the remnant spreads before it reaches the nose and throat.
fn long_plug_segments_with(
    front: f64,
    length: f64,
    total_kg: f64,
    back_weighted: bool,
    hot_remnant: bool,
) -> Vec<DepositSeg> {
    const SLICES: u32 = 20;
    let rod = 2.5;
    let p = rod * 75_000.0;
    let weight = |k: u32| {
        let x = (f64::from(k) + 0.5) / f64::from(SLICES);
        if back_weighted { 0.5 + x } else { 1.0 }
    };
    let norm: f64 = (0..SLICES).map(weight).sum();
    let dz = length / f64::from(SLICES);
    let mut moving = rod;
    let mut segs = Vec::new();
    let mut carried = 0.0;
    for k in 0..SLICES {
        let dm = total_kg * weight(k) / norm;
        let mut heat = 0.5 * p * p * (1.0 / moving - 1.0 / (moving + dm));
        moving += dm;
        if hot_remnant && 2 * k >= SLICES {
            carried += heat;
            heat = 0.0;
        }
        let z0 = front + f64::from(k) * dz;
        let last = k + 1 == SLICES;
        segs.push(DepositSeg {
            z0,
            z1: z0 + dz,
            radius: 0.10,
            inner: 0.0,
            heat: if last { heat + carried } else { heat },
            // The last slice carries the whole merged body; the others only their heat.
            mass: if last { rod + total_kg } else { 0.0 },
            uz: if last { p / (rod + total_kg) } else { 0.0 },
        });
    }
    segs
}

/// The rod as material, not as energy: 2.5 kg of polyethylene (950 kg/m^3), 3.5 cm radius and
/// 0.7 m long, its tip just inside the port, moving in along the axis at 75 km/s.
fn rod_segment() -> DepositSeg {
    DepositSeg {
        z0: 0.05,
        z1: 0.75,
        radius: 0.035,
        inner: 0.0,
        heat: 0.0,
        mass: 2.5,
        uz: 75_000.0,
    }
}

/// A plug column of `total_kg` at rest over `[front, front + length]`, radius set by its mean
/// density, in 20 slices, uniform or back-weighted (3x denser at the back).
fn plug_column(
    front: f64,
    length: f64,
    total_kg: f64,
    density: f64,
    back: bool,
) -> Vec<DepositSeg> {
    const SLICES: u32 = 20;
    let radius = (total_kg / (density * length * std::f64::consts::PI)).sqrt();
    let weight = |k: u32| {
        let x = (f64::from(k) + 0.5) / f64::from(SLICES);
        if back { 0.5 + x } else { 1.0 }
    };
    let norm: f64 = (0..SLICES).map(weight).sum();
    let dz = length / f64::from(SLICES);
    (0..SLICES)
        .map(|k| DepositSeg {
            z0: front + f64::from(k) * dz,
            z1: front + f64::from(k + 1) * dz,
            radius,
            inner: 0.0,
            heat: 0.0,
            mass: total_kg * weight(k) / norm,
            uz: 0.0,
        })
        .collect()
}

fn cmd_vessel_blast(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let shapes = [
        (
            "cylinder + 45 deg cone",
            vessel_for_volume(1.4, Some(1.4 - VESSEL_THROAT_R), 0.0, 0.0),
        ),
        (
            "cylinder + elliptical nose",
            vessel_for_volume(1.4, Some(1.4), 1.0, 0.0),
        ),
        ("curved cone", vessel_for_volume(1.8, None, 1.0, 0.0)),
        ("straight cone", vessel_for_volume(1.8, None, 0.0, 0.0)),
    ];
    let domed = (
        "domed head + elliptical nose",
        vessel_for_volume(1.4, Some(1.4), 1.0, 0.7),
    );
    let gases = [("methane 7000 K", 68.763), ("hydrogen 5500 K", 59.913)];
    let case = |shape: &'static str,
                vessel: euler2d::vessel::Vessel,
                gas: (&'static str, f64),
                gamma: f64,
                cell: f64,
                t_end: f64,
                deposition: &'static str,
                segs: Vec<DepositSeg>| VesselCase {
        shape,
        vessel,
        chamber: gas.0,
        charge_kg: gas.1,
        gamma,
        cell,
        t_end,
        deposition,
        segs,
        inflow: None,
        fill_band: None,
    };
    let mut cases = Vec::new();
    let path;
    if args.iter().any(|a| a == "--shape-40b") {
        // The bulge's plateau cut the hammer as 1/R^2, but its steep closing shoulder (0.8 m
        // over 0.5 m, ~68 deg) faced the downstream-skewed blast and took 2.3-3 GPa. Close the
        // bulge over 2-3 m instead (max wall angle ~32 / ~23 deg). On the 40 m^3 layout; the
        // volume grows. 0.5 cm, uniform fill.
        use euler2d::vessel::Bulge;
        path = "data/results/walled_nozzle/near_term/vessel_blast_shape_40b.jsonl";
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let r_t = (40.0 / (VESSEL_VOLUME / 0.149) / std::f64::consts::PI).sqrt();
        let cone_len = (1.4 - r_t) / 12.0_f64.to_radians().tan();
        let long = vessel_sized(1.4, Some(cone_len), 0.0, 0.7, 40.0, r_t);
        let bulge = |dr: f64, z1: f64, ramp_out: f64| Bulge {
            dr,
            z0: 1.4,
            z1,
            ramp: 0.5,
            ramp_out,
        };
        for (label, b, length) in [
            (
                "bulge r 2.2 m z 1.4-2.8, 2 m taper, plug 1.5 m",
                bulge(0.8, 2.8, 2.0),
                1.5,
            ),
            (
                "bulge r 2.2 m z 1.4-2.8, 3 m taper, plug 1.5 m",
                bulge(0.8, 2.8, 3.0),
                1.5,
            ),
            (
                "bulge r 2.2 m z 1.4-3.6, 2 m taper, plug 2.5 m",
                bulge(0.8, 3.6, 2.0),
                2.5,
            ),
            (
                "bulge r 1.8 m z 1.4-2.8, 2 m taper, plug 1.5 m",
                bulge(0.4, 2.8, 2.0),
                1.5,
            ),
        ] {
            let mut segs = vec![rod_segment()];
            segs.extend(plug_column(0.8, length, 10.0, 500.0, true));
            cases.push(case(
                "shape 40 m^3 tapered bulge",
                long.with_bulge(b),
                fill,
                1.2,
                0.005,
                1.2e-3,
                label,
                segs,
            ));
        }
    } else if args.iter().any(|a| a == "--shape-40") {
        // The 40 m^3 long cone's hammer beside the plug follows a strong line blast,
        // p ~ (E/L)/R^2, independent of the gas density (the cushion runs). So stand the wall off
        // locally (a bulge at the plug on the 40 m^3 layout, adding volume -- at fixed volume it
        // shortens the cylinder until the cone starts beside the plug), widen
        // the whole chamber, or spread the stopped energy over a longer plug. Uniform fill, rod
        // and 10 kg plug as material, 0.5 cm (1 cm smears the rod).
        use euler2d::vessel::Bulge;
        path = "data/results/walled_nozzle/near_term/vessel_blast_shape_40.jsonl";
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let r_t = (40.0 / (VESSEL_VOLUME / 0.149) / std::f64::consts::PI).sqrt();
        let cone = |r_c: f64, deg: f64| (r_c - r_t) / deg.to_radians().tan();
        let bulge = |dr: f64, z0: f64, z1: f64| Bulge {
            dr,
            z0,
            z1,
            ramp: 0.5,
            ramp_out: 0.5,
        };
        let long = vessel_sized(1.4, Some(cone(1.4, 12.0)), 0.0, 0.7, 40.0, r_t);
        let plug = |front: f64, length: f64| plug_column(front, length, 10.0, 500.0, true);
        let shapes: [(&'static str, euler2d::vessel::Vessel, f64); 6] = [
            (
                "40 m^3 cone, bulge r 1.8 m at z 1.4-2.8",
                long.with_bulge(bulge(0.4, 1.4, 2.8)),
                1.5,
            ),
            (
                "40 m^3 cone, bulge r 2.2 m at z 1.4-2.8",
                long.with_bulge(bulge(0.8, 1.4, 2.8)),
                1.5,
            ),
            (
                "40 m^3 wide r 1.7 m, 20 deg cone",
                vessel_sized(1.7, Some(cone(1.7, 20.0)), 0.0, 0.85, 40.0, r_t),
                1.5,
            ),
            (
                "40 m^3 cone, plug 0.8-3.3 m (2.5 m long)",
                vessel_sized(1.4, Some(cone(1.4, 12.0)), 0.0, 0.7, 40.0, r_t),
                2.5,
            ),
            (
                "40 m^3 cone, bulge r 2.2 m at z 1.4-3.6, plug 2.5 m",
                long.with_bulge(bulge(0.8, 1.4, 3.6)),
                2.5,
            ),
            (
                "40 m^3 cone, bulge r 1.8 m at z 1.4-3.6, plug 2.5 m",
                long.with_bulge(bulge(0.4, 1.4, 3.6)),
                2.5,
            ),
        ];
        for (label, v, length) in shapes {
            let mut segs = vec![rod_segment()];
            segs.extend(plug(0.8, length));
            cases.push(case(
                "shape 40 m^3",
                v,
                fill,
                1.2,
                0.005,
                1.2e-3,
                label,
                segs,
            ));
        }
    } else if args.iter().any(|a| a == "--cushion-40") {
        // The 40 m^3 long cone's remaining hot spot is the cylinder beside the plug, where the jet
        // stops and blasts sideways. Put most of the charge (as mist or membrane-held gas) in a
        // ring around the plug, z 0.5-3.0 m: either filling plug-to-wall, or as a layer against
        // the wall. Rod and 10 kg plug (0.8-2.3 m) as material, 0.5 cm.
        path = "data/results/walled_nozzle/near_term/vessel_blast_cushion_40.jsonl";
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let r_t = (40.0 / (VESSEL_VOLUME / 0.149) / std::f64::consts::PI).sqrt();
        let cone_len = (1.4 - r_t) / 12.0_f64.to_radians().tan();
        let v40 = vessel_sized(1.4, Some(cone_len), 0.0, 0.7, 40.0, r_t);
        for (r0, frac, label) in [
            (
                0.15,
                0.5,
                "40 m^3 cone, 50% of charge in ring plug-to-wall, z 0.5-3 m",
            ),
            (
                0.15,
                0.8,
                "40 m^3 cone, 80% of charge in ring plug-to-wall, z 0.5-3 m",
            ),
            (
                0.9,
                0.5,
                "40 m^3 cone, 50% of charge in wall layer r>0.9, z 0.5-3 m",
            ),
            (
                0.9,
                0.8,
                "40 m^3 cone, 80% of charge in wall layer r>0.9, z 0.5-3 m",
            ),
        ] {
            let mut segs = vec![rod_segment()];
            segs.extend(plug_column(0.8, 1.5, 10.0, 500.0, true));
            let mut c = case(
                "long cone 40 m^3",
                v40,
                fill,
                1.2,
                0.005,
                1.2e-3,
                label,
                segs,
            );
            c.fill_band = Some(FillBand {
                z0: 0.5,
                z1: 3.0,
                r0,
                r1: f64::INFINITY,
                fraction: frac,
            });
            cases.push(c);
        }
    } else if args.iter().any(|a| a == "--membrane") {
        // A dense gas layer in front of the nose, held by thin membranes (and a throat membrane),
        // so the jet ploughs into most of the charge before it reaches the wall. Rod and 10 kg plug
        // (0.8-2.3 m) as material, 0.5 cm.
        path = "data/results/walled_nozzle/near_term/vessel_blast_membrane.jsonl";
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let r_t = (40.0 / (VESSEL_VOLUME / 0.149) / std::f64::consts::PI).sqrt();
        let cone_len = (1.4 - r_t) / 12.0_f64.to_radians().tan();
        let v40 = vessel_sized(1.4, Some(cone_len), 0.0, 0.7, 40.0, r_t);
        let rod_plug = || {
            let mut segs = vec![rod_segment()];
            segs.extend(plug_column(0.8, 1.5, 10.0, 500.0, true));
            segs
        };
        for (frac, label20, label40) in [
            (
                0.5,
                "20 m^3, 50% of charge in the last metre",
                "40 m^3 cone, 50% of charge in 1 m before the cone",
            ),
            (
                0.8,
                "20 m^3, 80% of charge in the last metre",
                "40 m^3 cone, 80% of charge in 1 m before the cone",
            ),
        ] {
            let mut c = case(
                domed.0,
                domed.1,
                fill,
                1.2,
                0.005,
                1.2e-3,
                label20,
                rod_plug(),
            );
            c.fill_band = Some(FillBand::slab(domed.1.z_hi - 1.0, domed.1.z_hi, frac));
            cases.push(c);
            let mut c = case(
                "long cone 40 m^3",
                v40,
                fill,
                1.2,
                0.005,
                1.2e-3,
                label40,
                rod_plug(),
            );
            c.fill_band = Some(FillBand::slab(v40.z_c - 1.0, v40.z_c, frac));
            cases.push(c);
        }
    } else if args.iter().any(|a| a == "--nose-study-40") {
        // The long 12-degree cone at 40 m^3 (throat grown to keep V/A* = 134 m), plug in place
        // and with the rod entering through the port.
        path = "data/results/walled_nozzle/near_term/vessel_blast_nose_study_40.jsonl";
        let r_t = (40.0 / (VESSEL_VOLUME / 0.149) / std::f64::consts::PI).sqrt();
        let cone_len = (1.4 - r_t) / 12.0_f64.to_radians().tan();
        let v40 = vessel_sized(1.4, Some(cone_len), 0.0, 0.7, 40.0, r_t);
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let plug = |front: f64| plug_column(front, 1.5, 10.0, 500.0, true);
        let mut segs = vec![rod_segment()];
        segs.extend(plug(0.8));
        cases.push(case(
            "long cone 40 m^3",
            v40,
            fill,
            1.2,
            0.005,
            1.2e-3,
            "40 m^3 long cone, 10 kg 0.8-2.3 m",
            segs,
        ));
        let mut c = case(
            "long cone 40 m^3",
            v40,
            fill,
            1.2,
            0.005,
            1.2e-3,
            "40 m^3 long cone, rod enters, 10 kg 0.4-1.9 m",
            plug(0.4),
        );
        c.inflow = Some((0.035, 950.0, 75_000.0, 0.7 / 75_000.0));
        cases.push(c);
    } else if args.iter().any(|a| a == "--nose-study") {
        // Give the jet room and a glancing nose: convergence of the best case at 0.25 cm, a long
        // 12-degree conical nose, a narrower longer chamber, and the rod entering through the port
        // so the plug can sit closer to it.
        path = "data/results/walled_nozzle/near_term/vessel_blast_nose_study.jsonl";
        let cone_len = (1.4 - VESSEL_THROAT_R) / 12.0_f64.to_radians().tan();
        let long_cone = vessel_for_volume(1.4, Some(cone_len), 0.0, 0.7);
        let narrow = vessel_for_volume(1.2, Some(1.2), 1.0, 0.6);
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        let plug = |front: f64| plug_column(front, 1.5, 10.0, 500.0, true);
        let rod_in = Some((0.035, 950.0, 75_000.0, 0.7 / 75_000.0));
        let mut with_rod = |shape: &'static str,
                            v: euler2d::vessel::Vessel,
                            cell: f64,
                            t_end: f64,
                            label: &'static str,
                            front: f64| {
            let mut segs = vec![rod_segment()];
            segs.extend(plug(front));
            cases.push(case(shape, v, fill, 1.2, cell, t_end, label, segs));
        };
        with_rod(
            domed.0,
            domed.1,
            0.0025,
            0.8e-3,
            "elliptical, 10 kg 0.8-2.3 m @ 0.25 cm",
            0.8,
        );
        with_rod(
            "long 12 deg cone",
            long_cone,
            0.005,
            1.2e-3,
            "long cone, 10 kg 0.8-2.3 m",
            0.8,
        );
        with_rod(
            "narrow r 1.2 m",
            narrow,
            0.005,
            1.2e-3,
            "narrow, 10 kg 0.8-2.3 m",
            0.8,
        );
        for (shape, v, label) in [
            (domed.0, domed.1, "elliptical, rod enters, 10 kg 0.4-1.9 m"),
            (
                "long 12 deg cone",
                long_cone,
                "long cone, rod enters, 10 kg 0.4-1.9 m",
            ),
        ] {
            let mut c = case(shape, v, fill, 1.2, 0.005, 1.2e-3, label, plug(0.4));
            c.inflow = rod_in;
            cases.push(c);
        }
    } else if args.iter().any(|a| a == "--plug-position") {
        // Give the hot remnant room to spread: shorter 10 kg plugs that finish earlier. The rod
        // starts inside the chamber (0.05-0.75 m), so a plug cannot start before ~0.8 m without
        // modelling the rod's entry through the port. Rod and plug as material.
        path = "data/results/walled_nozzle/near_term/vessel_blast_plug_position.jsonl";
        for (front, length, label) in [
            (0.8, 1.5, "rod -> 10 kg plug, 0.8-2.3 m, back"),
            (0.8, 1.0, "rod -> 10 kg plug, 0.8-1.8 m, back"),
            (0.8, 0.75, "rod -> 10 kg plug, 0.8-1.55 m, back"),
            (1.0, 1.0, "rod -> 10 kg plug, 1.0-2.0 m, back"),
            (1.0, 0.75, "rod -> 10 kg plug, 1.0-1.75 m, back"),
        ] {
            let mut segs = vec![rod_segment()];
            segs.extend(plug_column(front, length, 10.0, 500.0, true));
            let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
            cases.push(case(domed.0, domed.1, fill, 1.2, 0.01, 2.0e-3, label, segs));
        }
    } else if args.iter().any(|a| a == "--rod-plug") {
        // The rod and plug as material: no deposition model. 10 kg plugs of frozen methane in a
        // polyethylene container (~500 kg/m^3), and the paper's 4.4 kg foam plug as reference.
        let fine = args.iter().any(|a| a == "--fine");
        path = if fine {
            "data/results/walled_nozzle/near_term/vessel_blast_rod_plug_fine.jsonl"
        } else {
            "data/results/walled_nozzle/near_term/vessel_blast_rod_plug.jsonl"
        };
        let cell = if fine { 0.005 } else { 0.01 };
        for (total, front, length, density, back, label) in [
            (
                4.4,
                1.4,
                0.98,
                475.0,
                false,
                "rod -> paper plug 4.4 kg, 1.4-2.38 m",
            ),
            (
                10.0,
                1.0,
                2.5,
                500.0,
                true,
                "rod -> 10 kg plug, 1.0-3.5 m, back",
            ),
            (
                10.0,
                1.4,
                1.75,
                500.0,
                true,
                "rod -> 10 kg plug, 1.4-3.15 m, back",
            ),
            (
                10.0,
                0.8,
                1.5,
                500.0,
                true,
                "rod -> 10 kg plug, 0.8-2.3 m, back",
            ),
        ] {
            let mut segs = vec![rod_segment()];
            segs.extend(plug_column(front, length, total, density, back));
            let fill = (gases[0].0, gases[0].1 - (total - 4.4));
            cases.push(case(domed.0, domed.1, fill, 1.2, cell, 2.0e-3, label, segs));
        }
    } else if args.iter().any(|a| a == "--plug-clear-of-nose") {
        // Plugs that end in the cylinder, short of the nose (z_c = 2.48 m), bounded pessimistically
        // with the hot remnant: 0.8 -> 2.3 m, back-weighted.
        path = "data/results/walled_nozzle/near_term/vessel_blast_plug_clear.jsonl";
        for (total, label) in [
            (10.0, "10 kg, 0.8-2.3 m, back, hot remnant"),
            (20.0, "20 kg, 0.8-2.3 m, back, hot remnant"),
            (40.0, "40 kg, 0.8-2.3 m, back, hot remnant"),
        ] {
            let fill = (gases[0].0, gases[0].1 - (total - 4.4));
            cases.push(case(
                domed.0,
                domed.1,
                fill,
                1.2,
                0.01,
                2.0e-3,
                label,
                long_plug_segments_with(0.8, 1.5, total, true, true),
            ));
        }
        for (total, label) in [
            (10.0, "10 kg, 0.8-2.3 m, back, cold remnant"),
            (40.0, "40 kg, 0.8-2.3 m, back, cold remnant"),
        ] {
            let fill = (gases[0].0, gases[0].1 - (total - 4.4));
            cases.push(case(
                domed.0,
                domed.1,
                fill,
                1.2,
                0.01,
                2.0e-3,
                label,
                long_plug_segments_with(0.8, 1.5, total, true, false),
            ));
        }
    } else if args.iter().any(|a| a == "--long-plug-hot") {
        path = "data/results/walled_nozzle/near_term/vessel_blast_long_plug_hot.jsonl";
        for (total, length, label) in [
            (4.4, 1.75, "4.4 kg, 1.75 m, back, hot remnant"),
            (4.4, 2.5, "4.4 kg, 2.5 m, back, hot remnant"),
            (10.0, 1.75, "10 kg, 1.75 m, back, hot remnant"),
            (10.0, 2.5, "10 kg, 2.5 m, back, hot remnant"),
        ] {
            let front = domed.1.r_c.min(3.5 - length);
            let fill = (gases[0].0, gases[0].1 - (total - 4.4));
            cases.push(case(
                domed.0,
                domed.1,
                fill,
                1.2,
                0.01,
                2.0e-3,
                label,
                long_plug_segments_with(front, length, total, true, true),
            ));
        }
    } else if args.iter().any(|a| a == "--long-plug-light") {
        // The paper's 4.4 kg plug, only lengthened: does length alone remove the hammer?
        path = "data/results/walled_nozzle/near_term/vessel_blast_long_plug_light.jsonl";
        for (length, back, label) in [
            (1.75, false, "4.4 kg, 1.75 m"),
            (2.5, false, "4.4 kg, 2.5 m"),
            (1.75, true, "4.4 kg, 1.75 m, back"),
            (2.5, true, "4.4 kg, 2.5 m, back"),
        ] {
            let front = domed.1.r_c.min(3.5 - length);
            cases.push(case(
                domed.0,
                domed.1,
                gases[0],
                1.2,
                0.01,
                2.0e-3,
                label,
                long_plug_segments(front, length, 4.4, back),
            ));
        }
    } else if args.iter().any(|a| a == "--long-plug-fine") {
        // Grid convergence of the best long plug: 10 kg over 2.5 m, back-weighted, 2/1/0.5 cm.
        path = "data/results/walled_nozzle/near_term/vessel_blast_long_plug_fine.jsonl";
        let front = domed.1.r_c.min(3.5 - 2.5);
        let fill = (gases[0].0, gases[0].1 - (10.0 - 4.4));
        for (cell, label) in [
            (0.02, "10 kg, 2.5 m, back @ 2 cm"),
            (0.01, "10 kg, 2.5 m, back @ 1 cm"),
            (0.005, "10 kg, 2.5 m, back @ 0.5 cm"),
        ] {
            cases.push(case(
                domed.0,
                domed.1,
                fill,
                1.2,
                cell,
                2.0e-3,
                label,
                long_plug_segments(front, 2.5, 10.0, true),
            ));
        }
    } else if args.iter().any(|a| a == "--long-plug") {
        path = "data/results/walled_nozzle/near_term/vessel_blast_long_plug.jsonl";
        let labels = [
            [
                ["10 kg, 1 m", "10 kg, 1.75 m", "10 kg, 2.5 m"],
                ["20 kg, 1 m", "20 kg, 1.75 m", "20 kg, 2.5 m"],
                ["40 kg, 1 m", "40 kg, 1.75 m", "40 kg, 2.5 m"],
            ],
            [
                [
                    "10 kg, 1 m, back",
                    "10 kg, 1.75 m, back",
                    "10 kg, 2.5 m, back",
                ],
                [
                    "20 kg, 1 m, back",
                    "20 kg, 1.75 m, back",
                    "20 kg, 2.5 m, back",
                ],
                [
                    "40 kg, 1 m, back",
                    "40 kg, 1.75 m, back",
                    "40 kg, 2.5 m, back",
                ],
            ],
        ];
        for (b, back) in [false, true].into_iter().enumerate() {
            for (i, total) in [10.0_f64, 20.0, 40.0].into_iter().enumerate() {
                for (j, length) in [1.0, 1.75, 2.5].into_iter().enumerate() {
                    let front = domed.1.r_c.min(3.5 - length);
                    // The mass beyond the paper's 4.4 kg plug is taken from the charge.
                    let remaining = gases[0].1 - (total - 4.4);
                    assert!(
                        remaining > 0.0,
                        "a {total} kg plug exceeds the methane charge"
                    );
                    let fill = (gases[0].0, remaining);
                    cases.push(case(
                        domed.0,
                        domed.1,
                        fill,
                        1.2,
                        0.01,
                        2.0e-3,
                        labels[b][i][j],
                        long_plug_segments(front, length, total, back),
                    ));
                }
            }
        }
    } else if args.iter().any(|a| a == "--deposition") {
        path = "data/results/walled_nozzle/near_term/vessel_blast_deposition.jsonl";
        cases = deposition_cases(domed.1)
            .into_iter()
            .map(|(label, segs)| case(domed.0, domed.1, gases[0], 1.2, 0.01, 2.0e-3, label, segs))
            .collect();
    } else if args.iter().any(|a| a == "--domed") {
        path = "data/results/walled_nozzle/near_term/vessel_blast_domed.jsonl";
        for gas in gases {
            cases.push(case(
                domed.0,
                domed.1,
                gas,
                1.2,
                0.02,
                VESSEL_T_END,
                "point",
                point_deposit(&domed.1),
            ));
        }
        for cell in [0.01, 0.005] {
            cases.push(case(
                domed.0,
                domed.1,
                gases[0],
                1.2,
                cell,
                1.5e-3,
                "point",
                point_deposit(&domed.1),
            ));
        }
    } else if args.iter().any(|a| a == "--fine") {
        path = "data/results/walled_nozzle/near_term/vessel_blast_fine.jsonl";
        let v = shapes[1].1;
        cases.push(case(
            shapes[1].0,
            v,
            gases[0],
            1.2,
            0.01,
            VESSEL_T_END,
            "point",
            point_deposit(&v),
        ));
    } else {
        path = RESULT_PATH_VESSEL_BLAST;
        for (shape, v) in shapes {
            for gas in gases {
                for gamma in [1.2, 1.4] {
                    cases.push(case(
                        shape,
                        v,
                        gas,
                        gamma,
                        0.02,
                        VESSEL_T_END,
                        "point",
                        point_deposit(&v),
                    ));
                }
            }
        }
    }
    let rows = par_map_with_progress("vessel-blast", &cases, vessel_blast_one);
    if let Some(dir) = Path::new(path).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(path, "vessel-blast", &rows, |r| {
        let peak = r.pressure_pa.iter().flatten().copied().fold(0.0, f64::max);
        println!(
            "rust: {:>27} {:>22} {:>15} g={:.1} cell={:.3} -> peak={:.2e} Pa",
            r.shape, r.deposition, r.chamber, r.gamma, r.cell_m, peak
        );
    })
}

// ---- Spray plate: the face-pressure history at the 150 t design point (fatigue and spall) ----
//
// The stratified (unmerged) pulse on a bare face at k = 8.52 and the 12 MN·s pulse's PuffSat
// masses, 4 m cloud and 4 m pulse, gap 0. This is the sharpest load the face sees. The whole
// `(t, p)` history is written so the steel's stress-wave and fatigue analysis
// (`puffsat.water_plate.plate_fatigue`) can drive the plate with it.

/// (closing speed [m/s], PuffSat mass per pulse [kg]) for 12 MN·s at η_jet 0.6 (β = 2.85).
const FACE_HISTORY_SPEEDS: [(f64, f64); 2] = [(45_580.0, 92.0), (65_130.0, 65.0)];
const FACE_HISTORY_K: f64 = 8.52;
const FACE_HISTORY_DEPTH: f64 = 4.0;
const RESULT_PATH_FACE_HISTORY: &str = "data/results/water_plate/spray_face_history.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct FaceHistoryRecord {
    material: String,
    w: f64,
    puffsat_kg: f64,
    k: f64,
    depth: f64,
    peak_wall_pressure: f64,
    /// Wall impulse per unit area [Pa s].
    wall_impulse: f64,
    time_s: Vec<f64>,
    pressure_pa: Vec<f64>,
    converged: bool,
}

fn cmd_spray_face_history(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in FACE_HISTORY_SPEEDS {
            cases.push((arm, spray.clone(), speed));
        }
    }
    let rows = par_map_with_progress("spray-face-history", &cases, |(arm, spray, speed)| {
        let (tube, _) = spray_layers_with_gap_k(
            &water,
            spray,
            *speed,
            (FACE_HISTORY_DEPTH, FACE_HISTORY_DEPTH),
            1,
            SPRAY_LEVER_CELLS,
            0.0,
            FACE_HISTORY_K,
        );
        let cfg = Config::production();
        let (r, history) =
            CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run_with_history();
        FaceHistoryRecord {
            material: arm.name.to_string(),
            w: speed.0,
            puffsat_kg: speed.1,
            k: FACE_HISTORY_K,
            depth: FACE_HISTORY_DEPTH,
            peak_wall_pressure: r.bounce.peak_wall_pressure,
            wall_impulse: r.bounce.wall_impulse,
            time_s: history.iter().map(|h| h.0).collect(),
            pressure_pa: history.iter().map(|h| h.1).collect(),
            converged: r.bounce.converged,
        }
    });
    if let Some(dir) = Path::new(RESULT_PATH_FACE_HISTORY).parent() {
        fs::create_dir_all(dir)?;
    }
    emit_scenario(RESULT_PATH_FACE_HISTORY, "spray-face-history", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} m={:.0} kg -> peak={:.2e} Pa samples={} converged={}",
            r.material,
            r.w,
            r.puffsat_kg,
            r.peak_wall_pressure,
            r.time_s.len(),
            r.converged,
        );
    })
}

// ---- Spray plate, step 2c: merge before the plate, in 2-D (ADR-0055) --------------------------
//
// The 4 m cloud lifted off the 10 m plate by a gap the Eulerian kernel carries as ambient gas, so
// the PuffSat merges in the cloud and the merged gas crosses the gap. Each gap runs free and
// confined. Both are compared to the confined no-gap run, the 2-D twin of the valid 1-D column.
// The ratios scale step 2a's real-physics no-gap numbers (effective-γ geometry, as in step 2b).

/// Gaps [m] between plate and cloud.
const SPRAY_2D_GAPS_M: [usize; 5] = [0, 1, 2, 4, 8];
const SPRAY_2D_GAP_CLOUD_M: usize = 4;
const RESULT_PATH_SPRAY_2D_GAP: &str = "data/results/water_plate/spray_2d_standoff.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct Spray2dGapRecord {
    gap: f64,
    nr: usize,
    nz: usize,
    /// Free run at this gap, over the confined no-gap run: impulse and peak facesheet pressure.
    impulse_vs_1d: f64,
    peak_vs_1d: f64,
    /// Confined run at this gap, over the confined no-gap run: what the gap does without spill.
    confined_impulse_vs_1d: f64,
    confined_peak_vs_1d: f64,
}

fn spray_2d_gap_case(gap_m: usize, refine: usize) -> Spray2dGapRecord {
    let unit = SPRAY_2D_PULSE_M;
    // Units of the 4 m pulse; gaps and depth are multiples of 1 m, so quarters of a unit.
    let column_quarters = SPRAY_2D_GAP_CLOUD_M + gap_m + 4 + 8;
    let z_max = column_quarters as f64 / 4.0;
    let nz = refine * SPRAY_2D_NZ_PER_UNIT * column_quarters / 4;
    let config = |gap: usize, confined: bool| {
        let base = SlugConfig {
            gamma: SPRAY_2D_GAMMA,
            mach: SPRAY_2D_MACH,
            r_foot: 5.0 / unit,
            length: 1.0,
            r_plate: 10.0 / unit,
            r_max: 14.0 / unit,
            z_max,
            nr: refine * SPRAY_2D_NR,
            nz,
            confined: false,
            shape: PlateShape::FlatGridAligned,
            taper_frac: 0.0,
            alpha_div: 0.0,
        };
        let cfg = if confined {
            SlugConfig {
                r_foot: base.r_max,
                r_plate: base.r_max,
                nr: 8,
                confined: true,
                ..base
            }
        } else {
            base
        };
        let spray = SprayCloud {
            depth: SPRAY_2D_GAP_CLOUD_M as f64 / unit,
            mass_ratio: SPRAY_K,
            standoff: gap as f64 / unit,
        };
        run_spray_bounce(&cfg, spray)
    };
    let reference = config(0, true);
    let free = config(gap_m, false);
    let confined = config(gap_m, true);
    Spray2dGapRecord {
        gap: gap_m as f64,
        nr: refine * SPRAY_2D_NR,
        nz,
        impulse_vs_1d: free.restitution_ratio() / reference.restitution_ratio(),
        peak_vs_1d: free.peak_local_pressure / reference.peak_local_pressure,
        confined_impulse_vs_1d: confined.restitution_ratio() / reference.restitution_ratio(),
        confined_peak_vs_1d: confined.peak_local_pressure / reference.peak_local_pressure,
    }
}

fn cmd_spray_2d_standoff(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases: Vec<(usize, usize)> = SPRAY_2D_GAPS_M.iter().map(|&g| (g, 1)).collect();
    cases.push((2, 2));
    let rows = par_map_with_progress("spray-2d-standoff", &cases, |&(g, refine)| {
        spray_2d_gap_case(g, refine)
    });
    // The valid 1-D real-physics no-gap columns: 4 m cloud, 4 m pulse (step 2a').
    let columns: Vec<SprayMixRecord> = fs::read_to_string(RESULT_PATH_SPRAY_LEVERS)?
        .lines()
        .map(serde_json::from_str::<SprayMixRecord>)
        .collect::<Result<_, _>>()?;
    emit_scenario(RESULT_PATH_SPRAY_2D_GAP, "spray-2d-standoff", &rows, |r| {
        println!(
            "rust: gap={:>3.0} m nr={:>3} -> impulse x{:.3} peak x{:.3} (confined: x{:.3}, x{:.3})",
            r.gap,
            r.nr,
            r.impulse_vs_1d,
            r.peak_vs_1d,
            r.confined_impulse_vs_1d,
            r.confined_peak_vs_1d,
        );
        if r.nr == SPRAY_2D_NR {
            for c in columns.iter().filter(|c| {
                c.depth == SPRAY_2D_GAP_CLOUD_M as f64 && c.pulse_length == SPRAY_2D_PULSE_M
            }) {
                let share = c.ceiling_share * r.impulse_vs_1d;
                println!(
                    "rust:   {:>5} w={:>5.0}: share {:.3} (eta_jet {:.3}), peak {:.2e} Pa",
                    c.material,
                    c.w,
                    share,
                    (share * (1.0 + (1.0 + SPRAY_K).sqrt()) - 1.0) / (1.0 + SPRAY_K).sqrt(),
                    c.peak_wall_pressure * r.peak_vs_1d,
                );
            }
        }
    })
}

// ---- Spray plate, step 2d: containing the spill with a dish, skirt or flared wall (ADR-0055) ----
//
// The 4 m cloud, 1 m off the floor, on the 10 m plate shaped as a dish, a cup with a straight
// skirt, or a cup with a flared wall. Each runs free and is compared to the confined flat no-gap
// run (the 2-D twin of the 1-D column); the ratios scale step 2a's real-physics 4 m numbers.
// Lengths in units of the 4 m pulse. The domain is 4.5 units wide so a flared lip fits; the cell
// size is the step 2b production size, 0.03125.

const SPRAY_CUP_DR: f64 = 0.03125;
const SPRAY_CUP_NR: usize = 144; // 4.5 units
const SPRAY_CUP_NZ: usize = 192; // 6 units
/// Wall thickness in the 2-D model: four cells, so ghost cells always have solid behind them.
const SPRAY_CUP_WALL: f64 = 4.0 * SPRAY_CUP_DR;
const RESULT_PATH_SPRAY_CUP: &str = "data/results/water_plate/spray_2d_cup.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayCupRecord {
    shape: String,
    refine: usize,
    impulse_vs_1d: f64,
    peak_vs_1d: f64,
    /// Steel added by the wall per centimetre of wall thickness [t], at 7850 kg/m^3.
    wall_tonnes_per_cm: f64,
}

/// `(label, d/D, wall height [m], flare)`; a zero height is the plain dish.
type CupShape = (&'static str, f64, f64, f64);

const SPRAY_CUP_SHAPES: [CupShape; 12] = [
    ("open plate", 0.0, 0.0, 0.0),
    ("dish 0.10", 0.10, 0.0, 0.0),
    ("dish 0.15", 0.15, 0.0, 0.0),
    ("skirt 2 m", 0.0, 2.0, 0.0),
    ("skirt 4 m", 0.0, 4.0, 0.0),
    ("skirt 8 m", 0.0, 8.0, 0.0),
    ("flare 0.3, 4 m", 0.0, 4.0, 0.3),
    ("flare 0.6, 4 m", 0.0, 4.0, 0.6),
    ("dish 0.10 + skirt 4 m", 0.10, 4.0, 0.0),
    ("dish 0.10 + flare 0.3, 4 m", 0.10, 4.0, 0.3),
    ("dish 0.10 + skirt 8 m", 0.10, 8.0, 0.0),
    ("dish 0.10 + skirt 12 m", 0.10, 12.0, 0.0),
];

fn spray_cup_case(shape: CupShape, refine: usize) -> SprayCupRecord {
    let (label, d_over_d, height_m, flare) = shape;
    let unit = SPRAY_2D_PULSE_M;
    let r_plate = 10.0 / unit;
    let plate_shape = if height_m > 0.0 {
        PlateShape::Cup {
            d_over_d,
            skirt_height: height_m / unit,
            flare,
            thickness: SPRAY_CUP_WALL,
        }
    } else {
        PlateShape::Dish { d_over_d }
    };
    let base = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate,
        r_max: SPRAY_CUP_NR as f64 * SPRAY_CUP_DR,
        z_max: SPRAY_CUP_NZ as f64 * SPRAY_CUP_DR,
        nr: refine * SPRAY_CUP_NR,
        nz: refine * SPRAY_CUP_NZ,
        confined: false,
        shape: plate_shape,
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let spray = |standoff: f64| SprayCloud {
        depth: 1.0,
        mass_ratio: SPRAY_K,
        standoff,
    };
    let free = run_spray_bounce(&base, spray(1.0 / unit));
    let reference = run_spray_bounce(
        &SlugConfig {
            r_foot: base.r_max,
            r_plate: base.r_max,
            nr: 8,
            confined: true,
            shape: PlateShape::FlatGridAligned,
            ..base
        },
        spray(0.0),
    );
    let r1 = 10.0;
    let r2 = r1 + flare * height_m;
    let wall_area = std::f64::consts::PI * (r1 + r2) * height_m * (1.0 + flare * flare).sqrt();
    SprayCupRecord {
        shape: label.to_string(),
        refine,
        impulse_vs_1d: free.restitution_ratio() / reference.restitution_ratio(),
        peak_vs_1d: free.peak_local_pressure / reference.peak_local_pressure,
        wall_tonnes_per_cm: wall_area * 0.01 * 7850.0 / 1000.0,
    }
}

fn cmd_spray_2d_cup(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut cases: Vec<(CupShape, usize)> = SPRAY_CUP_SHAPES.iter().map(|&s| (s, 1)).collect();
    cases.push((SPRAY_CUP_SHAPES[0], 2));
    cases.push((SPRAY_CUP_SHAPES[8], 2));
    let rows = par_map_with_progress("spray-2d-cup", &cases, |&(shape, refine)| {
        spray_cup_case(shape, refine)
    });
    let columns: Vec<SprayMixRecord> = fs::read_to_string(RESULT_PATH_SPRAY_LEVERS)?
        .lines()
        .map(serde_json::from_str::<SprayMixRecord>)
        .collect::<Result<_, _>>()?;
    emit_scenario(RESULT_PATH_SPRAY_CUP, "spray-2d-cup", &rows, |r| {
        println!(
            "rust: {:<28} x{} -> impulse x{:.3} peak x{:.3} wall {:>5.1} t/cm",
            r.shape, r.refine, r.impulse_vs_1d, r.peak_vs_1d, r.wall_tonnes_per_cm,
        );
        if r.refine == 1 {
            for c in columns
                .iter()
                .filter(|c| c.depth == 4.0 && c.pulse_length == 4.0)
            {
                let share = c.ceiling_share * r.impulse_vs_1d;
                println!(
                    "rust:     {:>5} w={:>5.0}: share {:.3} eta_jet {:.3} peak {:.2e} Pa",
                    c.material,
                    c.w,
                    share,
                    (share * (1.0 + (1.0 + SPRAY_K).sqrt()) - 1.0) / (1.0 + SPRAY_K).sqrt(),
                    c.peak_wall_pressure * r.peak_vs_1d,
                );
            }
        }
    })
}

// ---- Spray plate: oil against pitch on the shielded face (ADR-0055, Q40) ------------------------
//
// The 4 m cloud, 4 m pulse column on the Rung E ablating face, with each film's effective heat of
// ablation: carbon-loaded oil at 0.9 / 1.6 / 3.1 MJ/kg (docs/spray_plate_film_properties.md) and
// pitch at 59.3 MJ/kg (the wall solver's pitch removal heat), at the sourced vapor opacity and
// the soot-lost bound.

const SPRAY_FILMS: [(&str, f64); 4] = [
    ("oil low", 0.9e6),
    ("oil", 1.6e6),
    ("oil high", 3.1e6),
    ("pitch", 59.3e6),
];
const RESULT_PATH_SPRAY_FILM: &str = "data/results/water_plate/spray_film.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayFilmRecord {
    film: String,
    q_star: f64,
    run: SprayMixRecord,
}

fn cmd_spray_film(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let water = Table::load("data/tables/water_jupiter.json")?;
    let mut cases = Vec::new();
    for arm in SPRAY_ARMS {
        let spray = Table::load(arm.spray_table)?;
        for speed in SPRAY_SPEEDS {
            for film in SPRAY_FILMS {
                for kappa in [5.0e3, 100.0] {
                    cases.push((arm, spray.clone(), speed, film, kappa));
                }
            }
        }
    }
    let rows = par_map_with_progress("spray-film", &cases, |(arm, spray, speed, film, kappa)| {
        SprayFilmRecord {
            film: film.0.to_string(),
            q_star: film.1,
            run: run_one_spray_mix_q(
                *arm,
                (&water, spray),
                *speed,
                (4.0, 4.0),
                1,
                Some((*kappa, film.1)),
                SPRAY_LEVER_CELLS,
            ),
        }
    });
    emit_scenario(RESULT_PATH_SPRAY_FILM, "spray-film", &rows, |r| {
        println!(
            "rust: {:>5} w={:>5.0} {:<8} kv={:>5.0} -> film {:>6.1} kg/pulse share {:.3} face {:.2e}",
            r.run.material,
            r.run.w,
            r.film,
            r.run.kappa_vapor.unwrap_or(0.0),
            r.run.film_per_pulse,
            r.run.ceiling_share,
            r.run.face_radiation_share,
        );
    })
}

// ---- Spray plate: tuning the injection ratio k on the cup (ADR-0055, Q39) -----------------------
//
// For each k: the 1-D real-physics column at both mixing bounds (stratified layers, and the
// premixed slab on that k's mixture table), bare and on the pitch-shielded face, times the 2-D
// containment of the working cup (dish 0.10 + 4 m skirt, cloud 1 m off the floor) at that k's
// spray density. The 4 m cloud, 4 m pulse and the PuffSat masses per pulse are held fixed.

const SPRAY_K_VALUES: [f64; 5] = [4.0, 6.0, 8.5, 10.0, 14.0];
const PITCH_Q_STAR: f64 = 59.3e6;
const RESULT_PATH_SPRAY_K: &str = "data/results/water_plate/spray_k.jsonl";

fn spray_k_table(k: f64) -> String {
    if k == 10.0 {
        "data/tables/spray_k10.json".to_string()
    } else {
        format!("data/tables/spray_k{k}.json")
    }
}

/// `(tube, sigma_puffsat)`: one spray layer (at rest) on the plate, the PuffSat layer above it.
fn k_stratified(water: &Table, spray: &Table, (w, m): (f64, f64), k: f64) -> (Tube<TableEos>, f64) {
    let sigma_p = m / SPRAY_FOOTPRINT;
    let per_layer = SPRAY_LEVER_CELLS / 2;
    let dx = 4.0 / per_layer as f64;
    let t0 = Config::production().t0;
    let (rho_s, rho_p) = (sigma_p * k / 4.0, sigma_p / 4.0);
    let mut positions = vec![0.0];
    let (mut mass, mut vel, mut energy, mut index) = (vec![], vec![], vec![], vec![]);
    for (r, u, tbl, mat) in [(rho_s, 0.0, spray, 1), (rho_p, -w, water, 0)] {
        for _ in 0..per_layer {
            positions.push(positions[positions.len() - 1] + dx);
            mass.push(r * dx);
            vel.push(u);
            energy.push(tbl.energy(r, t0));
            index.push(mat);
        }
    }
    let n = mass.len();
    let mut node_velocities = vec![0.0; n + 1];
    for i in 1..n {
        node_velocities[i] =
            (mass[i - 1] * vel[i - 1] + mass[i] * vel[i]) / (mass[i - 1] + mass[i]);
    }
    node_velocities[n] = vel[n - 1];
    let eos = TableEos::layered(
        vec![TableEos::new(water.clone()), TableEos::new(spray.clone())],
        index,
    );
    let state = LagrangianState {
        positions,
        node_velocities,
        cell_masses: mass,
        specific_energies: energy,
    };
    let tube = Tube::from_lagrangian(
        state,
        eos,
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    (tube, sigma_p)
}

/// `(tube, sigma_puffsat)`: the merged slab, 4 m deep, closing at `w/(1+k)` with the merge heat
/// less what reaching the table's reference costs (water atomization for argon's atomic table).
fn k_premixed(table: &Table, argon: bool, (w, m): (f64, f64), k: f64) -> (Tube<TableEos>, f64) {
    let sigma_p = m / SPRAY_FOOTPRINT;
    let rho0 = sigma_p * (1.0 + k) / 4.0;
    let v0 = w / (1.0 + k);
    let charge = if argon {
        50.94e6 / (1.0 + k) + 0.116e6 * k / (1.0 + k)
    } else {
        1.89e6 * k / (1.0 + k)
    };
    let e0 = 0.5 * w * w * k / (1.0 + k).powi(2) - charge;
    let eos = TableEos::new(table.clone());
    let t0 = eos.temperature(rho0, e0);
    let cfg = Config::production();
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho0,
        v0,
        4.0,
        t0,
        eos,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    (tube, sigma_p)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayKRecord {
    k: f64,
    material: String,
    w: f64,
    /// "stratified" or "premixed".
    bound: String,
    /// 1-D share of the overtake ceiling, bare face.
    share_1d: f64,
    /// 2-D cup containment at this k (free cup over confined flat no-gap).
    containment: f64,
    relief: f64,
    /// Delivered share and η_jet on the cup.
    share: f64,
    eta_jet: f64,
    /// Impulse per PuffSat momentum `J/(m w)`.
    beta: f64,
    /// Impulse per kilogram consumed (PuffSat + spray + film), as a fraction of `w`.
    impulse_per_consumed: f64,
    peak_pressure: f64,
    film_per_pulse: f64,
    face_radiation_share: f64,
}

/// `(containment, relief)` of the working cup at injection ratio `k`.
fn k_containment(k: f64) -> (f64, f64) {
    let unit = SPRAY_2D_PULSE_M;
    let base = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate: 10.0 / unit,
        r_max: SPRAY_CUP_NR as f64 * SPRAY_CUP_DR,
        z_max: SPRAY_CUP_NZ as f64 * SPRAY_CUP_DR,
        nr: SPRAY_CUP_NR,
        nz: SPRAY_CUP_NZ,
        confined: false,
        shape: PlateShape::Cup {
            d_over_d: 0.10,
            skirt_height: 4.0 / unit,
            flare: 0.0,
            thickness: SPRAY_CUP_WALL,
        },
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let spray = |standoff: f64| SprayCloud {
        depth: 1.0,
        mass_ratio: k,
        standoff,
    };
    let free = run_spray_bounce(&base, spray(1.0 / unit));
    let reference = run_spray_bounce(
        &SlugConfig {
            r_foot: base.r_max,
            r_plate: base.r_max,
            nr: 8,
            confined: true,
            shape: PlateShape::FlatGridAligned,
            ..base
        },
        spray(0.0),
    );
    (
        free.restitution_ratio() / reference.restitution_ratio(),
        free.peak_local_pressure / reference.peak_local_pressure,
    )
}

fn cmd_spray_k(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let containment = par_map_with_progress("spray-k-2d", &SPRAY_K_VALUES, |&k| k_containment(k));
    let water = Table::load("data/tables/water_jupiter.json")?;
    let argon = Table::load("data/tables/argon.json")?;
    let mut cases = Vec::new();
    for (ik, &k) in SPRAY_K_VALUES.iter().enumerate() {
        let mixture = Table::load(spray_k_table(k))?;
        for material in ["argon", "water"] {
            for speed in SPRAY_SPEEDS {
                for bound in ["stratified", "premixed"] {
                    cases.push((ik, k, material, speed, bound, mixture.clone()));
                }
            }
        }
    }
    let rows = par_map_with_progress(
        "spray-k",
        &cases,
        |(ik, k, material, speed, bound, mixture)| {
            let is_argon = *material == "argon";
            let spray = if is_argon { &argon } else { &water };
            let build = || {
                if *bound == "stratified" {
                    k_stratified(&water, spray, *speed, *k)
                } else if is_argon {
                    k_premixed(mixture, true, *speed, *k)
                } else {
                    k_premixed(&water, false, *speed, *k)
                }
            };
            let cfg = Config::production();
            let (tube, sigma_p) = build();
            let bare = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
            let (tube, _) = build();
            let ablation = Ablation::new(PITCH_Q_STAR, cfg.t0).with_vapor_opacity(5.0e3);
            let shielded = AblatingBounce::new(tube, None, cfg.consts, cfg.limiter, ablation).run();
            let (w, m) = *speed;
            let ceiling = 1.0 + (1.0 + k).sqrt();
            let share_1d = bare.bounce.wall_impulse / (sigma_p * w) / ceiling;
            let (eta2d, relief) = containment[*ik];
            let share = share_1d * eta2d;
            let beta = share * ceiling;
            let film = shielded.ablated_mass * SPRAY_FOOTPRINT;
            SprayKRecord {
                k: *k,
                material: (*material).to_string(),
                w,
                bound: (*bound).to_string(),
                share_1d,
                containment: eta2d,
                relief,
                share,
                eta_jet: (beta - 1.0) / (1.0 + k).sqrt(),
                beta,
                impulse_per_consumed: beta / (1.0 + k + film / m),
                peak_pressure: bare.bounce.peak_wall_pressure * relief,
                film_per_pulse: film,
                face_radiation_share: bare.loss_radiative_wall / (0.5 * sigma_p * w * w),
            }
        },
    );
    emit_scenario(RESULT_PATH_SPRAY_K, "spray-k", &rows, |r| {
        println!(
            "rust: k={:>4.1} {:>5} w={:>5.0} {:<10} -> eta_jet {:.3} beta {:.2} J/consumed {:.3} w \
             peak {:.2e} film {:>4.1} kg (1-D {:.3} x cup {:.3})",
            r.k,
            r.material,
            r.w,
            r.bound,
            r.eta_jet,
            r.beta,
            r.impulse_per_consumed,
            r.peak_pressure,
            r.film_per_pulse,
            r.share_1d,
            r.containment,
        );
    })
}

// ---- Plug study step 2: one compact sphere into one plug, in the cup, in 2-D (ADR-0055, Q49) -----
//
// One necklace sphere (10 cm radius, 15 cm long, ~4.7 kg of water) buries itself in a plug of the
// same radius holding k = 10 times its mass (1.58 m of solid polyethylene, density ~0.95x the
// sphere's), on the axis. The plate is one sphere's share of the 20 m plate (3.2 m radius); its
// skirt stands in for the neighbouring fireballs. Lengths in metres, the sphere's speed 1, its
// density 1, effective γ 1.4. Reported as J/(m w) directly, beside the spray cup on the same method.

const PLUG_K: f64 = 10.0;
const PLUG_R: f64 = 0.1;
const PLUG_SPHERE_L: f64 = 0.15;
const PLUG_LENGTH: f64 = 1.58;
const PLUG_PLATE_R: f64 = 3.16;
const PLUG_DR: f64 = 0.02;
/// The ambient (vacuum stand-in) as a fraction of the sphere's density. The default 10⁻³ would weigh
/// ~5 times the plug on this domain, tamp the fireball and hold the plate above the quiet cutoff, so
/// the run integrated to its step cap (β 52, 12 times the ceiling); at 10⁻⁶ it is ~5% of the
/// sphere (`crates/euler2d/tests/plug_ambient.rs`).
const PLUG_AMBIENT: f64 = 1.0e-6;
const RESULT_PATH_PLUG_2D: &str = "data/results/water_plate/plug_2d.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct Plug2dRecord {
    case: String,
    /// `J / (m w)`: plate impulse per incoming momentum.
    beta: f64,
    /// `beta` over the overtake ceiling `1 + √(1+k)`.
    ceiling_share: f64,
    peak_local_pressure: f64,
    nr: usize,
    nz: usize,
}

/// One plug case on cells `PLUG_DR / refine` wide (`refine` 1 is the production grid).
fn plug_case(label: &str, standoff: f64, skirt: f64, refine: usize) -> Plug2dRecord {
    // Domain: the plate plus a 0.8 m margin; height covers the plug, the sphere and the skirt.
    let r_max = PLUG_PLATE_R + 0.84;
    let nr = 200 * refine;
    let dr = PLUG_DR / refine as f64;
    let top = (standoff + PLUG_LENGTH + PLUG_SPHERE_L + 2.0).max(skirt + 2.0);
    // Round the height up to a whole number of cells without a float-to-int cast.
    let mut nz = 0usize;
    while (nz as f64) * dr < top {
        nz += 1;
    }
    let shape = if skirt > 0.0 {
        PlateShape::Cup {
            d_over_d: 0.10,
            skirt_height: skirt,
            flare: 0.0,
            thickness: 4.0 * PLUG_DR,
        }
    } else {
        PlateShape::Dish { d_over_d: 0.10 }
    };
    let cfg = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: PLUG_R,
        length: PLUG_SPHERE_L,
        r_plate: PLUG_PLATE_R,
        r_max,
        z_max: nz as f64 * dr,
        nr,
        nz,
        confined: false,
        shape,
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let r = run_spray_bounce_with_ambient(
        &cfg,
        SprayCloud {
            depth: PLUG_LENGTH,
            mass_ratio: PLUG_K,
            standoff,
        },
        PLUG_AMBIENT,
    );
    let beta = r.restitution_ratio();
    Plug2dRecord {
        case: label.to_string(),
        beta,
        ceiling_share: beta / (1.0 + (1.0 + PLUG_K).sqrt()),
        peak_local_pressure: r.peak_local_pressure,
        nr,
        nz,
    }
}

/// The working spray cup (dish 0.10 + 4 m skirt, 4 m cloud 1 m off the floor) on the same 2-D
/// method, as `J/(m w)` directly.
fn spray_cup_beta() -> Plug2dRecord {
    let unit = SPRAY_2D_PULSE_M;
    let cfg = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate: 10.0 / unit,
        r_max: SPRAY_CUP_NR as f64 * SPRAY_CUP_DR,
        z_max: SPRAY_CUP_NZ as f64 * SPRAY_CUP_DR,
        nr: SPRAY_CUP_NR,
        nz: SPRAY_CUP_NZ,
        confined: false,
        shape: PlateShape::Cup {
            d_over_d: 0.10,
            skirt_height: 4.0 / unit,
            flare: 0.0,
            thickness: SPRAY_CUP_WALL,
        },
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let r = run_spray_bounce(
        &cfg,
        SprayCloud {
            depth: 1.0,
            mass_ratio: SPRAY_K,
            standoff: 1.0 / unit,
        },
    );
    let beta = r.restitution_ratio();
    Plug2dRecord {
        case: "spray cup (reference)".to_string(),
        beta,
        ceiling_share: beta / (1.0 + (1.0 + SPRAY_K).sqrt()),
        peak_local_pressure: r.peak_local_pressure,
        nr: SPRAY_CUP_NR,
        nz: SPRAY_CUP_NZ,
    }
}

/// `--plug-2d` runs every case and the spray reference; `--plug-2d N` runs only case `N` (0-5, or
/// 6 for the spray reference) and writes `plug_2d_N.jsonl`, so a slow case cannot cost the others.
/// `--plug-2d N R` runs case `N` on an `R`-times finer grid into `plug_2d_N_xR.jsonl`.
fn cmd_plug_2d(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let only: Option<usize> = args
        .iter()
        .skip_while(|a| *a != "--plug-2d")
        .nth(1)
        .and_then(|a| a.parse().ok());
    let refine: usize = args
        .iter()
        .skip_while(|a| *a != "--plug-2d")
        .nth(2)
        .and_then(|a| a.parse().ok())
        .unwrap_or(1);
    let cases: [(&str, f64, f64); 6] = [
        ("plug 1 m off, dish only", 1.0, 0.0),
        ("plug 1 m off, skirt 4 m", 1.0, 4.0),
        ("plug 2 m off, skirt 4 m", 2.0, 4.0),
        ("plug 1 m off, skirt 8 m", 1.0, 8.0),
        ("plug 2 m off, skirt 8 m", 2.0, 8.0),
        ("plug 4 m off, skirt 8 m", 4.0, 8.0),
    ];
    let (rows, path) = match only {
        Some(i) if i < cases.len() => {
            let (label, standoff, skirt) = cases[i];
            let path = if refine > 1 {
                format!("data/results/water_plate/plug_2d_{i}_x{refine}.jsonl")
            } else {
                format!("data/results/water_plate/plug_2d_{i}.jsonl")
            };
            (vec![plug_case(label, standoff, skirt, refine)], path)
        }
        Some(_) => (
            vec![spray_cup_beta()],
            "data/results/water_plate/plug_2d_spray.jsonl".to_string(),
        ),
        None => {
            let mut rows = par_map_with_progress("plug-2d", &cases, |&(label, standoff, skirt)| {
                plug_case(label, standoff, skirt, 1)
            });
            rows.push(spray_cup_beta());
            (rows, RESULT_PATH_PLUG_2D.to_string())
        }
    };
    emit_scenario(&path, "plug-2d", &rows, |r| {
        println!(
            "rust: {:<26} -> beta {:.3} ({:.3} of ceiling) peak {:.3e} [{}x{}]",
            r.case, r.beta, r.ceiling_share, r.peak_local_pressure, r.nr, r.nz,
        );
    })?;
    println!(
        "rust: ballistic plug fireball, all toward-plate gas caught: beta 2.739 (0.635 of ceiling)"
    );
    Ok(())
}

// ---- Spray plate: the skirt's outward load (ADR-0055, Q56) --------------------------------------
//
// The working spray (4 m cloud, 1 m off the floor, k = 10) in cups with a straight skirt, recording
// the gas pressure on the skirt's inner face row by row through the bounce. `skirt_hoop.py` scales
// the histories to the real pulse and sizes the liner and fibre wrap.

const RESULT_PATH_SPRAY_SKIRT: &str = "data/results/water_plate/spray_2d_skirt.jsonl";

/// `(label, d/D, skirt height [m])`.
const SKIRT_LOAD_CASES: [(&str, f64, f64); 3] = [
    ("dish 0.10 + skirt 4 m", 0.10, 4.0),
    ("dish 0.10 + skirt 8 m", 0.10, 8.0),
    ("flat + skirt 4 m", 0.0, 4.0),
];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayskirtRecord {
    case: String,
    d_over_d: f64,
    skirt_m: f64,
    /// Length unit [m] (the pulse length) and the plate radius in that unit.
    unit_m: f64,
    r_plate: f64,
    /// Floor impulse, `∫ Σ p r dr dt` with the common 2π dropped (as `Bounce2D`).
    wall_impulse: f64,
    incident_momentum: f64,
    /// Row heights, sample times and face pressures (`pressure[i][j]` at `time[i]`, row `z[j]`).
    z: Vec<f64>,
    time: Vec<f64>,
    pressure: Vec<Vec<f64>>,
}

fn spray_skirt_case(case: (&str, f64, f64)) -> SprayskirtRecord {
    let (label, d_over_d, skirt_m) = case;
    let unit = SPRAY_2D_PULSE_M;
    let cfg = SlugConfig {
        shape: PlateShape::Cup {
            d_over_d,
            skirt_height: skirt_m / unit,
            flare: 0.0,
            thickness: SPRAY_CUP_WALL,
        },
        ..pearl_cup_cfg(SPRAY_CUP_NZ)
    };
    let (r, history) = run_spray_bounce_with_skirt_history(&cfg, pearl_cloud(SPRAY_K));
    SprayskirtRecord {
        case: label.to_string(),
        d_over_d,
        skirt_m,
        unit_m: unit,
        r_plate: cfg.r_plate,
        wall_impulse: r.wall_impulse,
        incident_momentum: r.incident_momentum,
        z: history.z,
        time: history.time,
        pressure: history.pressure,
    }
}

/// `--spray-2d-skirt`: the skirt's face-pressure histories for each case.
fn cmd_spray_2d_skirt(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let rows = par_map_with_progress("spray-2d-skirt", &SKIRT_LOAD_CASES, |&c| {
        spray_skirt_case(c)
    });
    emit_scenario(RESULT_PATH_SPRAY_SKIRT, "spray-2d-skirt", &rows, |r| {
        println!(
            "rust: {:<24} -> beta {:.3}, {} wall rows, {} samples",
            r.case,
            r.wall_impulse / r.incident_momentum,
            r.z.len(),
            r.time.len(),
        );
    })
}

// ---- Spray plate: the 150 t plate's shape (ADR-0055, Q62-Q63) ---------------------------------
//
// Dish depth against skirt height at the paper's k = 8.52: each shape's impulse against the 1-D
// confined column, and the face-pressure histories on its skirt and on its dish's risers, which
// `plate_shape.py` sizes into liner and wrap. Histories keep every fourth step to bound the file.

const SHAPE_K: f64 = 8.52;
const SHAPE_DEPTHS: [f64; 5] = [0.10, 0.15, 0.20, 0.25, 0.30];
const SHAPE_SKIRTS_M: [f64; 3] = [0.0, 2.0, 4.0];
const SHAPE_DECIMATE: usize = 4;
const RESULT_PATH_SPRAY_SHAPE: &str = "data/results/water_plate/spray_2d_shape.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayShapeRecord {
    d_over_d: f64,
    skirt_m: f64,
    k: f64,
    unit_m: f64,
    r_plate: f64,
    impulse_vs_1d: f64,
    beta: f64,
    wall_impulse: f64,
    incident_momentum: f64,
    skirt_z: Vec<f64>,
    band_z: Vec<f64>,
    band_r: Vec<f64>,
    time: Vec<f64>,
    skirt_pressure: Vec<Vec<f64>>,
    band_pressure: Vec<Vec<f64>>,
}

fn spray_shape_case(d_over_d: f64, skirt_m: f64) -> SprayShapeRecord {
    let unit = SPRAY_2D_PULSE_M;
    let r_plate = 10.0 / unit;
    // Headroom for the dish's depth on top of the default domain, in whole cells.
    let depth = d_over_d * 2.0 * r_plate;
    let mut extra = 0usize;
    while (extra as f64) * SPRAY_CUP_DR < depth {
        extra += 1;
    }
    let nz = SPRAY_CUP_NZ + extra;
    let shape = if skirt_m > 0.0 {
        PlateShape::Cup {
            d_over_d,
            skirt_height: skirt_m / unit,
            flare: 0.0,
            thickness: SPRAY_CUP_WALL,
        }
    } else {
        PlateShape::Dish { d_over_d }
    };
    let base = SlugConfig {
        shape,
        ..pearl_cup_cfg(nz)
    };
    let (free, h) = run_spray_bounce_with_skirt_history(&base, pearl_cloud(SHAPE_K));
    let reference = run_spray_bounce(
        &SlugConfig {
            r_foot: base.r_max,
            r_plate: base.r_max,
            nr: 8,
            confined: true,
            shape: PlateShape::FlatGridAligned,
            ..base
        },
        SprayCloud {
            standoff: 0.0,
            ..pearl_cloud(SHAPE_K)
        },
    );
    let keep = |i: &usize| i % SHAPE_DECIMATE == 0;
    SprayShapeRecord {
        d_over_d,
        skirt_m,
        k: SHAPE_K,
        unit_m: unit,
        r_plate,
        impulse_vs_1d: free.restitution_ratio() / reference.restitution_ratio(),
        beta: free.restitution_ratio(),
        wall_impulse: free.wall_impulse,
        incident_momentum: free.incident_momentum,
        skirt_z: h.z,
        band_z: h.band_z,
        band_r: h.band_r,
        time: (0..h.time.len()).filter(keep).map(|i| h.time[i]).collect(),
        skirt_pressure: (0..h.pressure.len())
            .filter(keep)
            .map(|i| h.pressure[i].clone())
            .collect(),
        band_pressure: (0..h.band_pressure.len())
            .filter(keep)
            .map(|i| h.band_pressure[i].clone())
            .collect(),
    }
}

/// `--spray-2d-shape`: every dish depth and skirt height at k = 8.52, with wall histories.
fn cmd_spray_2d_shape(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let cases: Vec<(f64, f64)> = SHAPE_DEPTHS
        .iter()
        .flat_map(|&d| SHAPE_SKIRTS_M.iter().map(move |&h| (d, h)))
        .collect();
    let rows = par_map_with_progress("spray-2d-shape", &cases, |&(d, h)| spray_shape_case(d, h));
    emit_scenario(RESULT_PATH_SPRAY_SHAPE, "spray-2d-shape", &rows, |r| {
        println!(
            "rust: dish {:.2} + skirt {:>2.0} m -> impulse x{:.3} of 1-D, beta {:.3}, {} skirt rows, \
             {} band rows",
            r.d_over_d,
            r.skirt_m,
            r.impulse_vs_1d,
            r.beta,
            r.skirt_z.len(),
            r.band_z.len(),
        );
    })
}

// ---- Spray plate: a trailing pearl in formation behind the lead (ADR-0055, Q54) -----------------
//
// The working cup (dish 0.10 + 4 m skirt, 4 m cloud 1 m off the floor, k = 10 on the lead) with a
// second PuffSat of the same architecture trailing the lead by `gap` (tail to front), holding a
// `share` of the lead's mass. Each case is scored on the total incoming mass: the spray is
// `k_eff = k / (1 + share)` times it, and the fair baseline is one PuffSat carrying that total mass
// into the same cloud (the same cup run at `k_eff`). Lengths in units of the 4 m pulse.

const PEARL_SHARES: [f64; 4] = [0.05, 0.10, 0.20, 0.30];
/// Gaps from the lead's tail to the pearl's front [m]; centre to centre is 4 m more.
const PEARL_GAPS_M: [f64; 5] = [5.0, 10.0, 15.0, 25.0, 40.0];
const RESULT_PATH_SPRAY_PEARL: &str = "data/results/water_plate/spray_2d_pearl.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayPearlRecord {
    share: f64,
    gap_m: f64,
    /// Plate impulse per total incoming momentum (lead + pearl).
    beta: f64,
    /// `(beta - 1) / √(1 + k_eff)`, `k_eff = k / (1 + share)`.
    eta_jet: f64,
    /// The single-PuffSat baseline at the same total mass.
    eta_jet_single: f64,
    peak_local_pressure: f64,
    nz: usize,
}

/// The working cup's configuration with `nz` cells of height.
fn pearl_cup_cfg(nz: usize) -> SlugConfig {
    let unit = SPRAY_2D_PULSE_M;
    SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate: 10.0 / unit,
        r_max: SPRAY_CUP_NR as f64 * SPRAY_CUP_DR,
        z_max: nz as f64 * SPRAY_CUP_DR,
        nr: SPRAY_CUP_NR,
        nz,
        confined: false,
        shape: PlateShape::Cup {
            d_over_d: 0.10,
            skirt_height: 4.0 / unit,
            flare: 0.0,
            thickness: SPRAY_CUP_WALL,
        },
        taper_frac: 0.0,
        alpha_div: 0.0,
    }
}

fn pearl_cloud(mass_ratio: f64) -> SprayCloud {
    SprayCloud {
        depth: 1.0,
        mass_ratio,
        standoff: 1.0 / SPRAY_2D_PULSE_M,
    }
}

/// `eta_jet` of the single-PuffSat cup at spray ratio `k`.
fn pearl_single_eta(k: f64) -> f64 {
    let r = run_spray_bounce(&pearl_cup_cfg(SPRAY_CUP_NZ), pearl_cloud(k));
    (r.restitution_ratio() - 1.0) / (1.0 + k).sqrt()
}

fn spray_pearl_case(share: f64, gap_m: f64, eta_jet_single: f64) -> SprayPearlRecord {
    let unit = SPRAY_2D_PULSE_M;
    let gap = gap_m / unit;
    // Height: standoff + cloud + lead + gap + pearl, plus the default domain's 3.75 units of
    // headroom above the lead, rounded up to whole cells without a float-to-int cast.
    let top = 1.0 / unit + 1.0 + 1.0 + gap + 1.0 + 3.75;
    let mut nz = SPRAY_CUP_NZ;
    while (nz as f64) * SPRAY_CUP_DR < top {
        nz += 1;
    }
    let r = run_spray_bounce_with_pearl(
        &pearl_cup_cfg(nz),
        pearl_cloud(SPRAY_K),
        TrailingPearl {
            mass_ratio: share,
            gap,
        },
    );
    let beta = r.restitution_ratio();
    let k_eff = SPRAY_K / (1.0 + share);
    SprayPearlRecord {
        share,
        gap_m,
        beta,
        eta_jet: (beta - 1.0) / (1.0 + k_eff).sqrt(),
        eta_jet_single,
        peak_local_pressure: r.peak_local_pressure,
        nz,
    }
}

/// `--spray-2d-pearl`: every share and gap on the working cup, each beside its single-PuffSat
/// baseline at the same total mass.
fn cmd_spray_2d_pearl(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let singles = par_map_with_progress("pearl-baseline", &PEARL_SHARES, |&share| {
        pearl_single_eta(SPRAY_K / (1.0 + share))
    });
    let cases: Vec<(f64, f64, f64)> = PEARL_SHARES
        .iter()
        .zip(&singles)
        .flat_map(|(&share, &single)| PEARL_GAPS_M.iter().map(move |&gap| (share, gap, single)))
        .collect();
    let rows = par_map_with_progress("spray-2d-pearl", &cases, |&(share, gap, single)| {
        spray_pearl_case(share, gap, single)
    });
    emit_scenario(RESULT_PATH_SPRAY_PEARL, "spray-2d-pearl", &rows, |r| {
        println!(
            "rust: pearl {:>4.0}% gap {:>4.0} m -> beta {:.3} eta_jet {:.3} vs single {:.3} \
             ({:+.3}) peak {:.3e} [nz {}]",
            100.0 * r.share,
            r.gap_m,
            r.beta,
            r.eta_jet,
            r.eta_jet_single,
            r.eta_jet - r.eta_jet_single,
            r.peak_local_pressure,
            r.nz,
        );
    })
}

// ---- Spray plate: an inward lip on the cup's skirt (ADR-0055, Q50) -----------------------------
//
// The working cup (dish 0.10, 4 m skirt, 4 m cloud 1 m off the floor) with its skirt leaning inward
// (negative flare), so its mouth narrows into a wide throat. Besides the impulse and the floor's
// peak pressure, this records the peak pressure the skirt and lip feel and how long the gas stays
// inside the cup (time for the mass inside to fall to 1/e of its peak). Impulse is scaled onto the
// 1-D real-physics 4 m column, as in step 2d.

const SPRAY_LIP_FLARES: [f64; 5] = [0.0, -0.2, -0.4, -0.6, -1.0];
const RESULT_PATH_SPRAY_LIP: &str = "data/results/water_plate/spray_2d_lip.jsonl";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct SprayLipRecord {
    flare: f64,
    /// Mouth diameter [m].
    mouth: f64,
    impulse_vs_1d: f64,
    floor_peak_vs_1d: f64,
    /// Peak pressure on the skirt and lip, over the confined reference's floor peak.
    wall_peak_vs_1d: f64,
    /// Time the gas stays in the cup, in units of the pulse's crossing time `4 m / w`.
    residence: f64,
}

fn spray_lip_case(flare: f64) -> SprayLipRecord {
    let unit = SPRAY_2D_PULSE_M;
    let base = SlugConfig {
        gamma: SPRAY_2D_GAMMA,
        mach: SPRAY_2D_MACH,
        r_foot: 5.0 / unit,
        length: 1.0,
        r_plate: 10.0 / unit,
        r_max: SPRAY_CUP_NR as f64 * SPRAY_CUP_DR,
        z_max: SPRAY_CUP_NZ as f64 * SPRAY_CUP_DR,
        nr: SPRAY_CUP_NR,
        nz: SPRAY_CUP_NZ,
        confined: false,
        shape: PlateShape::Cup {
            d_over_d: 0.10,
            skirt_height: 4.0 / unit,
            flare,
            thickness: SPRAY_CUP_WALL,
        },
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let spray = |standoff: f64| SprayCloud {
        depth: 1.0,
        mass_ratio: SPRAY_K,
        standoff,
    };
    let reference = run_spray_bounce(
        &SlugConfig {
            r_foot: base.r_max,
            r_plate: base.r_max,
            nr: 8,
            confined: true,
            shape: PlateShape::FlatGridAligned,
            ..base
        },
        spray(0.0),
    );

    // The free run, stepped here so the wall pressure and the cup's contents can be read.
    let mut g = init_spray_grid(&base, spray(1.0 / unit));
    let incident = g.axial_momentum().abs();
    let (dz, dr) = (base.z_max / base.nz as f64, base.r_max / base.nr as f64);
    let z_rim = 4.0 * dz + 0.10 * 2.0 * base.r_plate;
    let z_top = z_rim + 4.0 / unit;
    let inside = |iz: usize, ir: usize| {
        let (z, r) = ((iz as f64 + 0.5) * dz, (ir as f64 + 0.5) * dr);
        let r_in = base.r_plate + flare * (z - z_rim).max(0.0);
        z < z_top && r < r_in.min(base.r_plate)
    };
    let cup_mass = |g: &euler2d::kernel::Grid2D| -> f64 {
        let mut m = 0.0;
        for iz in 0..g.nz() {
            for ir in 0..g.nr() {
                if inside(iz, ir) && !g.is_solid(iz, ir) {
                    m += g.prim(iz, ir).rho * g.cell_volume(ir);
                }
            }
        }
        m
    };
    let wall_pressure = |g: &euler2d::kernel::Grid2D| -> f64 {
        let mut p: f64 = 0.0;
        for iz in 0..g.nz() {
            for ir in 1..g.nr() - 1 {
                let r = (ir as f64 + 0.5) * dr;
                let beside_wall = g.is_solid(iz, ir + 1) || g.is_solid(iz, ir - 1);
                if r > 0.5 * base.r_plate && !g.is_solid(iz, ir) && beside_wall {
                    p = p.max(g.prim(iz, ir).p);
                }
            }
        }
        p
    };
    let (mut t, mut impulse, mut peak_force, mut floor_peak, mut wall_peak) =
        (0.0, 0.0, 0.0_f64, 0.0_f64, 0.0_f64);
    let (mut mass_peak, mut t_peak, mut t_drop) = (0.0_f64, 0.0, f64::NAN);
    let mut force_old = g.plate_force();
    // Arm the tail guard only once the pulse has arrived: under an inward overhang the ambient
    // gas pushing up on the lip makes the starting net force slightly negative, and an unarmed
    // guard would stop the run at once.
    let ambient = force_old.abs().max(1e-30);
    let (mut past_peak, mut quiet) = (false, 0usize);
    for _ in 0..(400 * base.nz + 50_000) {
        let dt = g.stable_dt();
        g.step(dt);
        t += dt;
        let force = g.plate_force();
        impulse += 0.5 * dt * (force_old + force);
        force_old = force;
        peak_force = peak_force.max(force);
        floor_peak = floor_peak.max(g.max_plate_pressure());
        wall_peak = wall_peak.max(wall_pressure(&g));
        let m = cup_mass(&g);
        if m > mass_peak {
            mass_peak = m;
            t_peak = t;
        } else if t_drop.is_nan() && m < mass_peak / std::f64::consts::E {
            t_drop = t;
        }
        let armed = peak_force > 1.0e3 * ambient;
        if armed && force < 0.999 * peak_force {
            past_peak = true;
        }
        quiet = if armed && force < 1.0e-3 * peak_force {
            quiet + 1
        } else {
            0
        };
        if past_peak && quiet >= 40 {
            break;
        }
    }
    SprayLipRecord {
        flare,
        mouth: 2.0 * (10.0 + flare * 4.0),
        impulse_vs_1d: impulse / incident / reference.restitution_ratio(),
        floor_peak_vs_1d: floor_peak / reference.peak_local_pressure,
        wall_peak_vs_1d: wall_peak / reference.peak_local_pressure,
        residence: t_drop - t_peak,
    }
}

fn cmd_spray_2d_lip(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let rows = par_map_with_progress("spray-2d-lip", &SPRAY_LIP_FLARES, |&f| spray_lip_case(f));
    let columns: Vec<SprayMixRecord> = fs::read_to_string(RESULT_PATH_SPRAY_LEVERS)?
        .lines()
        .map(serde_json::from_str::<SprayMixRecord>)
        .collect::<Result<_, _>>()?;
    emit_scenario(RESULT_PATH_SPRAY_LIP, "spray-2d-lip", &rows, |r| {
        println!(
            "rust: flare {:+.1} (mouth {:>4.1} m) -> impulse x{:.3} floor peak x{:.3} wall peak \
             x{:.3} residence {:.2} crossings",
            r.flare, r.mouth, r.impulse_vs_1d, r.floor_peak_vs_1d, r.wall_peak_vs_1d, r.residence,
        );
        for c in columns
            .iter()
            .filter(|c| c.depth == 4.0 && c.pulse_length == 4.0)
        {
            let share = c.ceiling_share * r.impulse_vs_1d;
            println!(
                "rust:     {:>5} w={:>5.0}: eta_jet {:.3}",
                c.material,
                c.w,
                (share * (1.0 + (1.0 + SPRAY_K).sqrt()) - 1.0) / (1.0 + SPRAY_K).sqrt(),
            );
        }
    })
}

// ---- Heavy-plate 16–28 km/s scenario sweep (special scenario, design §12.1 / ADR-0027) ----------
//
// A 100 kg pulse on a tripled 30 m-diameter (`R = 15 m`) pusher plate of mass `≤ 40 t`, swept
// 16–28 km/s at 0.5 km/s (25 anchors; 16 km/s overlaps the core study anchor as a consistency
// check). It reuses the extended-grid Jupiter table (real TOPS/OPLIB opacity, `T` to 1.2×10⁶ K),
// which already covers the 16–28 km/s stagnation window. Because `R` grows 3× while `m` grows 4×,
// the Σ contract dilutes the column below the core baseline (`ρ ≈ 0.01–0.6`): survivability-easy on
// the facesheet, but flirting with `τ ~ 1` at the wide-footprint corner — hence the opacity τ-check.
//
// `e_eff` is the equilibrium coupled bounce (`wall = None`) on the real-opacity table across all 25
// velocities, at a fixed representative stretched-cloud length; the `L` axis is spot-checked at the
// bracket anchors (`τ ≫ 1` makes `e_eff` `L`-insensitive). `eta_capture` is reused from the M = 40
// geometry sweep (scale-invariant), and the whole-plate structural go/no-go is the closed-form
// companion (ADR-0027, `puffsat.structure`), decoupled from `f(v)`.

/// The **two-leg** velocity grid (Q17), 16–63 km/s in [`HEAVY_N_V`] anchors, built by
/// [`heavy_v_grid`] from these endpoints.
///
/// Two legs because the two halves of the range are known to different degrees. Over 16–43 the
/// existing study already showed `e_eff` varying smoothly, so 3 km/s is reconnaissance. Over
/// 45–63 — the extension's own range — the curve has never been sampled and the `τ` structure is
/// unknown, so it gets 1 km/s. The 43 → 45 seam is a 2 km/s join: finer than the coarse leg,
/// coarser than the fine one, so nothing is lost across it.
const HEAVY_V_COARSE_LO: f64 = 16_000.0;
const HEAVY_V_COARSE_STEP: f64 = 3_000.0;
const HEAVY_N_V_COARSE: usize = 10;
const HEAVY_V_FINE_LO: f64 = 45_000.0;
const HEAVY_V_FINE_STEP: f64 = 1_000.0;
const HEAVY_N_V_FINE: usize = 19;
const HEAVY_N_V: usize = HEAVY_N_V_COARSE + HEAVY_N_V_FINE;
/// The freeze-timing / `L`-sensitivity bracket anchors [m/s] (design §12.1: 16 / 22 / 28 km/s).
/// All three must land on the coarse leg — Q4 runs the full frozen-table pipeline at exactly these.
const HEAVY_V_ANCHORS: [f64; 3] = [16_000.0, 22_000.0, 28_000.0];
/// Impact-density grid [kg/m³] (Q11), 0.01–0.58: the dilute wide-footprint corner up to the dense
/// tight-disk corner the Σ contract reaches at `m = 100 kg`, `R = 15 m`.
///
/// Deliberately non-uniform. It is **dense through 0.06–0.20**, the steep on-contour band where
/// `ρ_ceiling(v)` actually lives and where Q15 interpolates the curve, and **coarse across the flat
/// top**, where `e_eff` barely moves and extra points would only cost run time. The three points
/// below 0.06 exist to support Q15's validation points, not to resolve a region the contour spends
/// time in.
const HEAVY_RHO: [f64; 12] = [
    0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.28, 0.40, 0.58,
];
/// Fixed representative stretched-cloud length [m] the headline `e_eff(v, ρ)` grid runs at (design
/// §12.1: `L ≈ 8–10 m`). The Python frontier reads this slice as its headline.
const HEAVY_LENGTH: f64 = 10.0;
/// `L`-sensitivity spot-check lengths [m]; the headline `HEAVY_LENGTH` = 10 m is the reference from
/// the main grid, so these bracket it below and above, giving each spot-checked state a three-point
/// `L` curve.
const HEAVY_L_SPOT: [f64; 2] = [6.0, 14.0];
/// Velocities [m/s] the diagnostic rows (`L`-sensitivity and the opacity τ-check) run at.
///
/// **Deliberately wider than the old single 28 km/s point.** Both diagnostics were previously
/// justified by `τ ≫ 1` holding across the 16–28 band, and both of those justifications expire in
/// the extension: `L` sets when the slab goes transparent, and the τ-check exists precisely to
/// detect that. 45/55/63 km/s are where the claims are most likely to fail, so they are where the
/// checks now run. Distinct from [`HEAVY_V_ANCHORS`], which stays at three points because Q4's
/// freeze-timing pipeline behind it is expensive.
const HEAVY_DIAG_V: [f64; 6] = [16_000.0, 22_000.0, 28_000.0, 45_000.0, 55_000.0, 63_000.0];
/// Density subset [kg/m³] for the diagnostic rows — a spot-check spans the range rather than
/// resolving it, and the full 12-point grid at every diagnostic velocity would triple the sweep for
/// no extra signal. Every entry must be in [`HEAVY_RHO`], so the opacity bracket can pair these rows
/// with the `κ = 1` headline row and the probe's turnaround state at the same `(v, ρ)`.
const HEAVY_SPOT_RHO: [f64; 4] = [0.01, 0.04, 0.10, 0.28];
/// Opacity scales for the τ-check (scale 1.0 is the headline, already in the main grid). If `e_eff`
/// does not move across these, `τ ≫ 1` is confirmed there and the headline is robust.
const HEAVY_TAU_SCALES: [f64; 4] = [0.1, 0.3, 3.0, 10.0];

/// One heavy-plate row: the swept `(v, ρ, length, opacity_scale)` case with the restitution, the
/// physical peak wall pressure (facesheet sizing), the incident/impulse pair (whole-plate structural
/// companion, ADR-0027), and the radiative loss split. Same schema as [`JupiterRecord`] (a sibling
/// special scenario), kept a distinct type so the arm stays self-contained.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct HeavyPlateRecord {
    v: f64,
    rho_impact: f64,
    /// Slug column length [m] (headline `HEAVY_LENGTH`; the spot-check rows vary it).
    length: f64,
    /// Opacity scale applied to the table (τ-check bracket; headline 1.0).
    opacity_scale: f64,
    e_eff: f64,
    /// Peak physical wall pressure (EOS `p`, AV excluded) — the survivability load (ADR-0010).
    peak_wall_pressure: f64,
    /// Incident momentum per unit wall area `ρLv` [Pa·s] (ADR-0027 areal-impulse input).
    incident_momentum: f64,
    /// Delivered wall impulse per unit area `≈ (1 + e_eff)·ρLv` [Pa·s] (ADR-0027 areal impulse).
    wall_impulse: f64,
    loss_radiative_wall: f64,
    loss_escape_space: f64,
    /// Whether the bounce reached its tail guard rather than exhausting its step budget
    /// (`BounceResult::converged`, Q6). A `false` row stopped mid-infall and its `e_eff` is an
    /// artifact, not a low restitution — the analysis must reject it, not average it in.
    converged: bool,
}

/// Run one heavy-plate coupled bounce at `(v, rho, length)` on an opacity-scaled table.
fn run_one_heavyplate(
    v: f64,
    rho_impact: f64,
    length: f64,
    scale: f64,
    base_tbl: &Table,
) -> HeavyPlateRecord {
    let table = base_tbl.with_opacity_scale(scale);
    let cfg = Config {
        v,
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
    HeavyPlateRecord {
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
        converged: result.bounce.converged,
    }
}

/// The 25-anchor heavy-plate velocity grid [m/s]: 16–28 km/s at 0.5 km/s.
fn heavy_v_grid() -> Vec<f64> {
    let coarse = (0..HEAVY_N_V_COARSE).map(|i| HEAVY_V_COARSE_LO + HEAVY_V_COARSE_STEP * i as f64);
    let fine = (0..HEAVY_N_V_FINE).map(|i| HEAVY_V_FINE_LO + HEAVY_V_FINE_STEP * i as f64);
    let mut grid = Vec::with_capacity(HEAVY_N_V);
    grid.extend(coarse.chain(fine));
    grid
}

/// Sweep the heavy-plate grid in parallel (rayon), input order: the headline `v(25) × ρ(7)` at fixed
/// `L` and `κ = 1`, plus the `L`-sensitivity spot rows at the anchors and the opacity τ-check rows at
/// the 28 km/s top.
/// The heavy-plate case list `(v, ρ, length, opacity_scale)`: the headline grid plus the two
/// diagnostic families. Factored out so its shape is testable without running any physics — a
/// duplicated or missing case is invisible in the output but corrupts every average taken over it.
fn heavyplate_cases() -> Vec<(f64, f64, f64, f64)> {
    let mut cases: Vec<(f64, f64, f64, f64)> = Vec::new();
    // Headline grid: all 29 velocities × 12 densities at the representative length, real opacity.
    for &v in &heavy_v_grid() {
        for &rho in &HEAVY_RHO {
            cases.push((v, rho, HEAVY_LENGTH, 1.0));
        }
    }
    // L-sensitivity spot-check (the headline L = HEAVY_LENGTH comes from the grid above, so each
    // spot-checked state ends up with a three-point L curve).
    for &v in &HEAVY_DIAG_V {
        for &rho in &HEAVY_SPOT_RHO {
            for &l in &HEAVY_L_SPOT {
                cases.push((v, rho, l, 1.0));
            }
        }
    }
    // Opacity τ-check (κ = 1 comes from the main grid).
    for &v in &HEAVY_DIAG_V {
        for &rho in &HEAVY_SPOT_RHO {
            for &s in &HEAVY_TAU_SCALES {
                cases.push((v, rho, HEAVY_LENGTH, s));
            }
        }
    }
    cases
}

fn run_heavyplate_sweep(base_tbl: &Table) -> Vec<HeavyPlateRecord> {
    let cases = heavyplate_cases();
    par_map_with_progress("heavyplate", &cases, |&(v, rho, l, s)| {
        run_one_heavyplate(v, rho, l, s, base_tbl)
    })
}

// ---- Sₙ transport check of the FLD escape channel (Q9/Q21, ADR-0012) ----------------------------
//
// FLD is exact where `τ ≫ 1` and degrades as the cloud thins. Extending the study to 63 km/s pushes
// the re-expansion toward `τ ~ 1`, so before trusting `e_eff` up there the escape-to-space channel
// gets audited against an independent radiation model: a gray discrete-ordinates solve on the same
// states (`hydro1d::transport`). The audit is one-way — it observes the FLD run without changing
// it — so what it produces is a **bias estimate on a loss channel**, not a corrected `e_eff`.
//
// The gate (design §12.1 step 5): a bias above 10% on a channel that carries meaningful energy
// escalates to a coupled M1 solver. Below that, the bias is recorded as the radiation-model error
// bar and the FLD result stands. `escape_share_of_ke` is reported alongside precisely so a large
// relative bias on a negligible channel is not mistaken for a large error on `f`.

const RESULT_PATH_TRANSPORT: &str = "data/results/sweep_transport_check.jsonl";
const RESULT_PATH_TRANSPORT_RES: &str = "data/results/sweep_transport_resolution.jsonl";

/// Velocities [m/s] the audit runs at: the heavy-plate band's three bracket anchors, the four new
/// anchors the extension adds, and the 69 km/s Jupiter anchor that brackets it from above.
const TRANSPORT_V: [f64; 8] = [
    16_000.0, 22_000.0, 28_000.0, 35_000.0, 45_000.0, 55_000.0, 63_000.0, 69_000.0,
];
/// Impact densities [kg/m³]: the three most dilute points of the heavy-plate grid — where `τ` is
/// smallest and FLD is weakest — plus one dense point as the `τ ≫ 1` control.
const TRANSPORT_RHO: [f64; 4] = [0.01, 0.02, 0.04, 0.3];
/// Angular resolution `S₃₂` (16 ordinates per hemisphere). The escape-flux convergence test in
/// `hydro1d::transport` puts S64 within 5·10⁻³ even for the limb-brightened `τ = 0.1` case, so S32
/// is comfortable for the `τ ≳ 0.5` states here.
const TRANSPORT_HALF_ORDINATES: usize = 16;
/// Keep one sampled comparison per this many substeps, bounding the trace on a ~10⁴-step run.
const TRANSPORT_SAMPLE_STRIDE: usize = 50;

/// One transport-audit row: what FLD said the escape channel carried, what Sₙ says, and how much
/// the channel matters.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct TransportCheckRecord {
    v: f64,
    rho_impact: f64,
    /// Gas cell count. Fixed at the production value for the main grid; the resolution study is
    /// the axis that varies it, because the space-facing cell's optical depth -- and therefore
    /// FLD's free-surface error -- is set by how finely the mesh resolves the photosphere.
    gas_cells: usize,
    length: f64,
    e_eff: f64,
    /// FLD's own escape integral (kernel channel 1b).
    loss_escape_space: f64,
    /// Sₙ escape integral on the Rosseland mean — the closure-comparison tally.
    transport_escape_rosseland: f64,
    /// Sₙ escape integral on the Planck mean — the mean-selection tally.
    transport_escape_planck: f64,
    /// `(Sₙ_rosseland − FLD)/FLD`. **The gate reads this.**
    relative_bias: f64,
    /// `(Sₙ_planck − FLD)/FLD`; its spread against `relative_bias` is the mean-selection band.
    relative_bias_planck: f64,
    /// Worst instantaneous `|bias|` over the run, and the slab `τ` where it happened.
    worst_relative_difference: f64,
    worst_optical_depth: f64,
    /// Escape-weighted mean slab optical depth: where the escaping energy actually came from.
    flux_weighted_optical_depth: f64,
    /// Escape-weighted optical depth of the space-facing cell alone. FLD's Marshak boundary is only
    /// meaningful while this is well under 1; above that its error is the mesh, not the closure.
    flux_weighted_surface_cell_depth: f64,
    /// Escape-weighted **Planck** optical depth of the space-facing cell — the number that says
    /// whether FLD's Marshak boundary is being applied to a cell it can legitimately describe.
    flux_weighted_surface_cell_depth_planck: f64,
    /// Escape-weighted `FLD flux / σT⁴` of the space-facing *cell*. Above 1 is expected wherever
    /// that cell is optically thin — the escaping radiation was emitted deeper and hotter — so this
    /// measures how far below the last cell the photosphere sits, not a physical violation.
    flux_weighted_fld_over_blackbody: f64,
    /// Escape loss as a fraction of the slug's incident kinetic energy. A bias on a channel worth
    /// 0.1% of the energy budget cannot move `f`, however large the ratio looks.
    escape_share_of_ke: f64,
    converged: bool,
}

/// Run one audited coupled bounce and reduce it to a [`TransportCheckRecord`].
fn run_one_transport_check(v: f64, rho_impact: f64, base_tbl: &Table) -> TransportCheckRecord {
    run_one_transport_check_at(v, rho_impact, Config::production().gas_cells, base_tbl)
}

/// The same audited bounce at an explicit cell count, for the resolution study.
fn run_one_transport_check_at(
    v: f64,
    rho_impact: f64,
    gas_cells: usize,
    base_tbl: &Table,
) -> TransportCheckRecord {
    let cfg = Config {
        v,
        length: HEAVY_LENGTH,
        gas_cells,
        ..Config::production()
    };
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        TableEos::new(base_tbl.clone()),
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let mut bounce = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).with_transport_audit(
        TransportAudit::new(TRANSPORT_HALF_ORDINATES, TRANSPORT_SAMPLE_STRIDE),
    );
    let result = bounce.run();
    let audit = bounce
        .transport_audit()
        .expect("audit was attached before the run");
    // Incident KE per unit wall area: p = ρLv and KE = ½ρLv², so KE = p·v/2.
    let incident_ke = 0.5 * result.bounce.incident_momentum * cfg.v;
    TransportCheckRecord {
        v: cfg.v,
        rho_impact,
        gas_cells: cfg.gas_cells,
        length: cfg.length,
        e_eff: result.bounce.e_eff,
        loss_escape_space: result.loss_escape_space,
        transport_escape_rosseland: audit.transport_escape_rosseland,
        transport_escape_planck: audit.transport_escape_planck,
        relative_bias: audit.relative_bias(),
        relative_bias_planck: audit.relative_bias_planck(),
        worst_relative_difference: audit.worst_relative_difference,
        worst_optical_depth: audit.worst_optical_depth,
        flux_weighted_optical_depth: audit.flux_weighted_optical_depth,
        flux_weighted_surface_cell_depth: audit.flux_weighted_surface_cell_depth,
        flux_weighted_surface_cell_depth_planck: audit.flux_weighted_surface_cell_depth_planck,
        flux_weighted_fld_over_blackbody: audit.flux_weighted_fld_over_blackbody,
        escape_share_of_ke: result.loss_escape_space / incident_ke,
        converged: result.bounce.converged,
    }
}

// The resolution gate (2026-08-18). Before building a fix for FLD's free surface, two things have
// to be measured, because both change what the fix is worth:
//
//   1. **Is the diagnosis right?** If the escape bias does not shrink as the mesh refines, the
//      mechanism is not resolution and the fix would be aimed at the wrong target.
//   2. **Does `e_eff` care?** The bias reaches the deliverable only through the energy books, and
//      `e_eff` is an integral momentum ratio that converges fast. If `e_eff` is already
//      mesh-converged while the escape bias is not, the escape error is largely not propagating
//      into `f`, and the first-order `delta_f` estimate is an overestimate.
//
// Refining is *not* proposed as the fix: this is a Lagrangian mesh, the space-facing cell's optical
// depth spans 0.2-48 across the grid, and convergence is first order -- no static grading reaches
// all of it. This is a measurement.

/// Cell counts the resolution gate sweeps. Production is 300.
const TRANSPORT_RES_CELLS: [usize; 4] = [150, 300, 600, 1200];
/// Representative `(v, ρ)` states for the gate: the dilute high-v corner that carries the most
/// escape energy, the mid-band anchor, the dense case whose surface cell is deepest in the Planck
/// mean, and the published 28 km/s heavy-plate top.
const TRANSPORT_RES_STATES: [(f64, f64); 4] = [
    (69_000.0, 0.02),
    (45_000.0, 0.04),
    (28_000.0, 0.01),
    (69_000.0, 0.30),
];

const RESULT_PATH_MESH: &str = "data/results/sweep_mesh_convergence.jsonl";

/// One mesh-convergence row: `e_eff` alone, at an explicit cell count.
///
/// Deliberately **unaudited**. The Sₙ audit is a bit-identical observer (proved by
/// `transport_audit_changes_nothing_and_tallies_the_escape_channel`), so dropping it leaves `e_eff`
/// unchanged while removing an Sₙ solve per substep — which is what makes 4800-cell meshes
/// affordable. The question here is whether the *deliverable* has converged, not what the audit
/// says about it.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct MeshConvergenceRecord {
    v: f64,
    rho_impact: f64,
    length: f64,
    gas_cells: usize,
    e_eff: f64,
    peak_wall_pressure: f64,
    incident_momentum: f64,
    wall_impulse: f64,
    loss_radiative_wall: f64,
    loss_escape_space: f64,
    converged: bool,
}

/// Cell ladder for the deep mesh-convergence study. Production is 300; each rung doubles, and cost
/// grows roughly as cells² (the CFL step count scales with the mesh too).
const MESH_CELLS: [usize; 5] = [300, 600, 1200, 2400, 4800];

/// Run one unaudited coupled bounce at an explicit cell count.
fn run_one_mesh_case(
    v: f64,
    rho_impact: f64,
    gas_cells: usize,
    base_tbl: &Table,
) -> MeshConvergenceRecord {
    let cfg = Config {
        v,
        length: HEAVY_LENGTH,
        gas_cells,
        ..Config::production()
    };
    let tube = Tube::slug_si(
        cfg.gas_cells,
        rho_impact,
        cfg.v,
        cfg.length,
        cfg.t0,
        TableEos::new(base_tbl.clone()),
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let result = CoupledBounce::new(tube, None, cfg.consts, cfg.limiter).run();
    MeshConvergenceRecord {
        v: cfg.v,
        rho_impact,
        length: cfg.length,
        gas_cells,
        e_eff: result.bounce.e_eff,
        peak_wall_pressure: result.bounce.peak_wall_pressure,
        incident_momentum: result.bounce.incident_momentum,
        wall_impulse: result.bounce.wall_impulse,
        loss_radiative_wall: result.loss_radiative_wall,
        loss_escape_space: result.loss_escape_space,
        converged: result.bounce.converged,
    }
}

/// Sweep the deep mesh-convergence ladder over the gate's four states.
fn run_mesh_convergence_sweep(base_tbl: &Table) -> Vec<MeshConvergenceRecord> {
    let cases: Vec<(f64, f64, usize)> = TRANSPORT_RES_STATES
        .iter()
        .flat_map(|&(v, rho)| MESH_CELLS.iter().map(move |&n| (v, rho, n)))
        .collect();
    par_map_with_progress("mesh-convergence", &cases, |&(v, rho, n)| {
        run_one_mesh_case(v, rho, n, base_tbl)
    })
}

/// Sweep the resolution gate: the same states at four mesh resolutions.
fn run_transport_resolution_sweep(base_tbl: &Table) -> Vec<TransportCheckRecord> {
    let cases: Vec<(f64, f64, usize)> = TRANSPORT_RES_STATES
        .iter()
        .flat_map(|&(v, rho)| TRANSPORT_RES_CELLS.iter().map(move |&n| (v, rho, n)))
        .collect();
    par_map_with_progress("transport-resolution", &cases, |&(v, rho, n)| {
        run_one_transport_check_at(v, rho, n, base_tbl)
    })
}

/// Sweep the transport audit over `(v, ρ)` in parallel (rayon), input order.
fn run_transport_check_sweep(base_tbl: &Table) -> Vec<TransportCheckRecord> {
    let cases: Vec<(f64, f64)> = TRANSPORT_V
        .iter()
        .flat_map(|&v| TRANSPORT_RHO.iter().map(move |&rho| (v, rho)))
        .collect();
    par_map_with_progress("transport-check", &cases, |&(v, rho)| {
        run_one_transport_check(v, rho, base_tbl)
    })
}

// ---- Pulse-shape sensitivity sweep (design §13, ADR-0028) ---------------------------------------
//
// Raw `f(shape)` at the FIXED best-survivable baseline design (`d/D = 0.1` headline with a flat
// cross-check, `L/D = 0.3`, `r_foot/R = 0.5`, M = 20): the pulse varies over the assumed shape box
// while the plate, domain, and grid stay frozen — ADR-0028 decision 3 deliberately inverts the
// geometry sweep's size-the-domain-to-the-cloud convention, so IBM stair-step jitter cannot land
// at the `Δf ~ 0.005–0.02` differences being resolved. Two output files:
//
// - 2D (`sweep_shape_geometry.jsonl`): `eta_capture` per shape sample on the fixed grid, plus
//   refined-resolution repeats of a noise-floor subset (`σ_noise`; a flagged cliff must survive
//   these before it is called physical).
// - 1D (`sweep_shape_sigma.jsonl`): fresh equilibrium `e_eff` at each sample's physical `(ρ, L)`
//   via the Σ contract at fixed pulse mass and plate radius — never interpolated from frontier
//   data (ADR-0028 decision 2). The taper axis runs at the mass-weighted mean Σ with Σ-hi/lo90
//   bound rows (the named halt condition); divergence leaves Σ untouched (no 1D rows — the
//   Python assembly reuses the nominal `e_eff`).

const RESULT_PATH_SHAPE_2D: &str = "data/results/sweep_shape_geometry.jsonl";
const RESULT_PATH_SHAPE_1D: &str = "data/results/sweep_shape_sigma.jsonl";

/// The fixed plate curvatures: the shallow-concave headline and the flat cross-check (design §13).
const SHAPE_D_OVER_D: [f64; 2] = [0.10, 0.0];
/// The geometry track's strong-shock Mach anchor (the nominal design point was pinned at M = 20).
const SHAPE_MACH: f64 = 20.0;
/// The nominal pulse shape (the best-survivable baseline, design §13).
const SHAPE_NOM_RFOOT_OVER_R: f64 = 0.5;
const SHAPE_NOM_L_OVER_D: f64 = 0.3;
/// Two-sided box sampling for the `r_foot/R` and `L/D` axes: ±5/10/20% around nominal.
const SHAPE_REL: [f64; 7] = [0.80, 0.90, 0.95, 1.0, 1.05, 1.10, 1.20];
/// One-sided edge-taper axis (ramp width as a fraction of `r_foot`; 0 = the top-hat nominal).
const SHAPE_TAPER: [f64; 4] = [0.0, 0.10, 0.20, 0.30];
/// One-sided radial-divergence axis (`α` in `v_r = α·v·r/r_foot`).
const SHAPE_ALPHA: [f64; 4] = [0.0, 0.033, 0.067, 0.10];
/// The frozen plate/domain (ADR-0028 decision 3): nominal footprint `r_foot = 1` (`r_plate = 2`,
/// as in the geometry sweep's `r_foot/R = 0.5` case), domain sized once to the largest box sample
/// (`r_foot ≤ 1.2`, `L ≤ 0.72` above the 0.4-deep dish).
const SHAPE_R_PLATE: f64 = 2.0;
const SHAPE_R_MAX: f64 = 2.8;
const SHAPE_Z_MAX: f64 = 3.4;
/// `(resolution_scale, nr, nz)`: the base grid (the converged geometry-sweep cell size on the
/// fixed domain) plus the 1.5×/2× refinements the noise-floor subset re-runs (≥ 3 repeats,
/// ADR-0028).
const SHAPE_RES: [(f64, usize, usize); 3] = [(1.0, 112, 88), (1.5, 168, 132), (2.0, 224, 176)];
/// Physical Σ-contract anchors for the 1D arm (design §2: 25 kg pulse on the R = 5 m plate) at
/// the two quoted velocity anchors (the ~11 km/s dip and 16 km/s).
const SHAPE_PULSE_MASS_KG: f64 = 25.0;
const SHAPE_PLATE_RADIUS_M: f64 = 5.0;
const SHAPE_V: [f64; 2] = [11_000.0, 16_000.0];

/// One pulse-shape sample: the axis it perturbs (the Python `S_x` extraction key) and the four
/// §13 shape coordinates. One-factor-at-a-time — `S` is a per-axis derivative, not a response
/// surface.
#[derive(Debug, Clone, Copy)]
struct ShapeSample {
    axis: &'static str,
    r_foot_over_r: f64,
    l_over_d: f64,
    taper_frac: f64,
    alpha_div: f64,
}

/// The nominal sample (the fixed design point; `axis = "nominal"`).
fn shape_nominal() -> ShapeSample {
    ShapeSample {
        axis: "nominal",
        r_foot_over_r: SHAPE_NOM_RFOOT_OVER_R,
        l_over_d: SHAPE_NOM_L_OVER_D,
        taper_frac: 0.0,
        alpha_div: 0.0,
    }
}

/// The shape-box samples: the nominal once, then each axis perturbed alone.
fn shape_samples() -> Vec<ShapeSample> {
    let mut v = vec![shape_nominal()];
    for &rel in &SHAPE_REL {
        if (rel - 1.0).abs() > 1e-12 {
            v.push(ShapeSample {
                axis: "r_foot_over_r",
                r_foot_over_r: SHAPE_NOM_RFOOT_OVER_R * rel,
                ..shape_nominal()
            });
            v.push(ShapeSample {
                axis: "l_over_d",
                l_over_d: SHAPE_NOM_L_OVER_D * rel,
                ..shape_nominal()
            });
        }
    }
    for &taper in &SHAPE_TAPER[1..] {
        v.push(ShapeSample {
            axis: "taper_frac",
            taper_frac: taper,
            ..shape_nominal()
        });
    }
    for &alpha in &SHAPE_ALPHA[1..] {
        v.push(ShapeSample {
            axis: "alpha_div",
            alpha_div: alpha,
            ..shape_nominal()
        });
    }
    v
}

/// The noise-floor subset re-run at the refined resolutions: the nominal, the two samples that
/// stress the grid hardest (the widest footprint — most rim cells — and the widest taper — the
/// shallowest density gradient across cells), plus the two samples the first base-resolution
/// pass flagged as cliff candidates (`L/D` −10%, taper 20%) — the §13 rule is that a flag must
/// survive refinement before it is called physical, so the flagged samples need refined repeats.
fn shape_noise_subset() -> Vec<ShapeSample> {
    vec![
        shape_nominal(),
        ShapeSample {
            axis: "r_foot_over_r",
            r_foot_over_r: SHAPE_NOM_RFOOT_OVER_R * 1.20,
            ..shape_nominal()
        },
        ShapeSample {
            axis: "taper_frac",
            taper_frac: SHAPE_TAPER[3],
            ..shape_nominal()
        },
        ShapeSample {
            axis: "l_over_d",
            l_over_d: SHAPE_NOM_L_OVER_D * 0.90,
            ..shape_nominal()
        },
        ShapeSample {
            axis: "taper_frac",
            taper_frac: SHAPE_TAPER[2],
            ..shape_nominal()
        },
    ]
}

/// One 2D shape-sweep row: the sample, the plate/resolution it ran on, and the `eta_capture`
/// pieces (with the *measured* initialized `p_in` — the §13 normalization, never the analytic
/// value).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct ShapeRecord {
    axis: String,
    d_over_d: f64,
    r_foot_over_r: f64,
    l_over_d: f64,
    taper_frac: f64,
    alpha_div: f64,
    mach: f64,
    resolution_scale: f64,
    eta_capture: f64,
    restitution_free: f64,
    restitution_confined: f64,
    /// Measured initialized `p_in` of the free run.
    incident_momentum: f64,
    /// Peak local facesheet pressure of the free run (the survivability focusing input, Rung S).
    peak_local_pressure: f64,
}

/// Cloud geometry of a sample on the frozen grid: `r_foot = (r_foot/R)·r_plate`,
/// `L = 2·(L/D)·r_foot`.
fn shape_cloud(sample: &ShapeSample) -> (f64, f64) {
    let r_foot = sample.r_foot_over_r * SHAPE_R_PLATE;
    (r_foot, 2.0 * sample.l_over_d * r_foot)
}

/// Run one free shape bounce on the frozen plate/domain at resolution `(nr, nz)`.
fn run_shape_free(
    sample: &ShapeSample,
    d_over_d: f64,
    nr: usize,
    nz: usize,
) -> euler2d::bounce::Bounce2D {
    let (r_foot, length) = shape_cloud(sample);
    run_slug_bounce(&SlugConfig {
        gamma: 1.4,
        mach: SHAPE_MACH,
        r_foot,
        length,
        r_plate: SHAPE_R_PLATE,
        r_max: SHAPE_R_MAX,
        z_max: SHAPE_Z_MAX,
        nr,
        nz,
        confined: false,
        shape: PlateShape::Dish { d_over_d },
        taper_frac: sample.taper_frac,
        alpha_div: sample.alpha_div,
    })
}

/// The plane-wave confined denominator for a column of `length` — always the top-hat 1D-limit
/// reference (taper/divergence are free-run shape perturbations, not 1D-limit properties).
fn run_shape_confined(length: f64, nz: usize) -> euler2d::bounce::Bounce2D {
    run_slug_bounce(&SlugConfig {
        gamma: 1.4,
        mach: SHAPE_MACH,
        r_foot: SHAPE_R_MAX,
        length,
        r_plate: SHAPE_R_MAX,
        r_max: SHAPE_R_MAX,
        z_max: SHAPE_Z_MAX,
        nr: 8,
        nz,
        confined: true,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    })
}

/// The `(sample, d/D, resolution index)` 2D case list: the full box × both plates at the base
/// resolution, plus the noise-floor subset × both plates at the refinements.
fn shape_2d_cases() -> Vec<(ShapeSample, f64, usize)> {
    let mut cases = Vec::new();
    for s in shape_samples() {
        for &dd in &SHAPE_D_OVER_D {
            cases.push((s, dd, 0));
        }
    }
    for s in shape_noise_subset() {
        for &dd in &SHAPE_D_OVER_D {
            for res in [1, 2] {
                cases.push((s, dd, res));
            }
        }
    }
    cases
}

/// Run the full 2D shape sweep in parallel (rayon): confined denominators once per unique
/// `(length, resolution)` (the `r_foot/R` and `L/D` axes share the length set `L = 0.6·rel`),
/// then every free case against its denominator.
fn run_shape_2d_sweep() -> Vec<ShapeRecord> {
    let cases = shape_2d_cases();
    let mut keys: Vec<(u64, usize)> = cases
        .iter()
        .map(|(s, _, res)| (shape_cloud(s).1.to_bits(), *res))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    let confined: HashMap<(u64, usize), f64> =
        par_map_with_progress("shape-2d confined", &keys, |&(len_bits, res)| {
            let b = run_shape_confined(f64::from_bits(len_bits), SHAPE_RES[res].2);
            ((len_bits, res), b.restitution_ratio())
        })
        .into_iter()
        .collect();
    par_map_with_progress("shape-2d free", &cases, |&(s, dd, res)| {
        let (scale, nr, nz) = SHAPE_RES[res];
        let free = run_shape_free(&s, dd, nr, nz);
        let rc = confined[&(shape_cloud(&s).1.to_bits(), res)];
        ShapeRecord {
            axis: s.axis.to_string(),
            d_over_d: dd,
            r_foot_over_r: s.r_foot_over_r,
            l_over_d: s.l_over_d,
            taper_frac: s.taper_frac,
            alpha_div: s.alpha_div,
            mach: SHAPE_MACH,
            resolution_scale: scale,
            eta_capture: free.restitution_ratio() / rc,
            restitution_free: free.restitution_ratio(),
            restitution_confined: rc,
            incident_momentum: free.incident_momentum,
            peak_local_pressure: free.peak_local_pressure,
        }
    })
}

/// One 1D Σ-contract row of the shape study: the shape sample it serves, the Σ role (`"sample"`
/// for the headline runs, `"taper_mean"`/`"sigma_hi"`/`"sigma_lo90"` for the taper bookkeeping),
/// and the equilibrium coupled-bounce result at the physical `(ρ, L)`.
///
/// Each row carries a **two-resolution validity protocol** (2026-07-16 finding, root-caused and
/// fixed 2026-07-17): with the old ρ ≤ 20 kg/m³ table the coupled runs had a resolution-onset
/// radiative collapse — the radiatively-cooled wall cell was compressed past the table's density
/// ceiling, where the *clamped* `p(ρ)` no longer arrests the Lagrangian compression, so the cell
/// width → 0, `dt` → 0, and the run stalled mid-infall with an unphysical `e_eff`. The onset `dx`
/// coarsened as `ρv²` rose (ρ = 0.64 collapsed at 1200 cells, ρ ≈ 1.2–1.7 by 200–300). The fix is
/// the node-preserving ρ-grid extension to 1000 kg/m³ in `tables.py` (plus the Newton-iterated
/// backward-Euler exchange hardening in `hydro1d::radiation`); `diag_shape_high_sigma` verifies
/// the plateau now extends past every measured old onset. The protocol is retained as a cheap
/// regression guard: the row records the 300-cell headline, a 150-cell coarse cross-check, and
/// the EOS-only reference; the Python assembly accepts the headline only when the two coupled
/// runs agree, falls back to the coarse value otherwise, and flags the row.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
struct Shape1DRecord {
    v: f64,
    axis: String,
    sigma_role: String,
    r_foot_over_r: f64,
    l_over_d: f64,
    taper_frac: f64,
    /// Areal density `Σ = ρ·L` [kg/m²] the row runs at (the taper roles scale it off nominal).
    sigma: f64,
    rho_impact: f64,
    length: f64,
    /// Coupled `e_eff` at the production 300-cell convention (the headline when valid).
    e_eff: f64,
    /// Coupled `e_eff` at 150 cells — the stable-window cross-check/fallback.
    e_eff_coarse: f64,
    /// EOS-only `e_eff` (no radiation) at 300 cells — the smoothness reference the validity
    /// check anchors to (the gas dynamics is verified smooth across the whole Σ box).
    e_eff_eos: f64,
    peak_wall_pressure: f64,
    /// Peak physical wall pressure of the 150-cell coupled run (used when the fine run is
    /// rejected).
    peak_wall_pressure_coarse: f64,
    incident_momentum: f64,
    wall_impulse: f64,
}

/// The Σ contract at the fixed 25 kg / `R = 5 m` design, matching `puffsat.analysis.
/// impact_density` (the cross-language contract): `ρ = m/(2π·(L/D)·r_foot³)` with
/// `r_foot = (r_foot/R)·R`, and `L = 2·(L/D)·r_foot`.
fn shape_rho_length(r_foot_over_r: f64, l_over_d: f64) -> (f64, f64) {
    let r_foot = r_foot_over_r * SHAPE_PLATE_RADIUS_M;
    let length = 2.0 * l_over_d * r_foot;
    let rho = SHAPE_PULSE_MASS_KG / (2.0 * std::f64::consts::PI * l_over_d * r_foot.powi(3));
    (rho, length)
}

/// The 1D case list `(v, sample, sigma_role, Σ factor)`: every footprint/aspect sample fresh at
/// its own `(ρ, L)`; each taper sample at its mass-weighted mean Σ, with the Σ-hi/lo90 bound rows
/// at the widest taper (the ADR-0028 halt-condition check); divergence emits no 1D rows.
fn shape_1d_cases() -> Vec<(f64, ShapeSample, &'static str, f64)> {
    let mut cases = Vec::new();
    for &v in &SHAPE_V {
        for s in shape_samples() {
            match s.axis {
                "alpha_div" => {}
                "taper_frac" => {
                    let st = taper_sigma_stats(s.taper_frac);
                    cases.push((v, s, "taper_mean", st.mean));
                    if (s.taper_frac - SHAPE_TAPER[3]).abs() < 1e-12 {
                        cases.push((v, s, "sigma_hi", st.hi));
                        cases.push((v, s, "sigma_lo90", st.lo90));
                    }
                }
                _ => cases.push((v, s, "sample", 1.0)),
            }
        }
    }
    cases
}

// The §13 three-point frozen-chemistry spot-check at the dip anchor: Σ-nominal and the footprint
// box endpoints, through the standard ADR-0026 two-stage frozen pipeline on their own directory.
const RESULT_PATH_FROZEN_PROBE_SHAPE: &str = "data/results/frozen_probe_shape.jsonl";
const RESULT_PATH_FROZEN_SHAPE: &str = "data/results/sweep_frozen_shape.jsonl";
const TABLE_DIR_FROZEN_SHAPE: &str = "data/tables/frozen_shape";

/// The three dip-anchor spot-check densities (design §13): the Σ box endpoints and Σ-nominal as
/// Σ-equivalent densities at the *fixed nominal length* (`Σ = ρ·L` is the 1D knob at `τ ≫ 1`, so
/// fixing `L` and moving `ρ` probes the same axis the footprint box moves — and keeps the frozen
/// pipeline's fixed-`length` `Config` contract). Descending Σ: tightest footprint, nominal,
/// widest.
fn shape_frozen_rho_grid() -> [f64; 3] {
    let (_, l_nom) = shape_rho_length(SHAPE_NOM_RFOOT_OVER_R, SHAPE_NOM_L_OVER_D);
    let sigma = |rf_over_r: f64| {
        let r_foot = rf_over_r * SHAPE_PLATE_RADIUS_M;
        SHAPE_PULSE_MASS_KG / (std::f64::consts::PI * r_foot * r_foot)
    };
    [
        sigma(SHAPE_NOM_RFOOT_OVER_R * SHAPE_REL[0]) / l_nom,
        sigma(SHAPE_NOM_RFOOT_OVER_R) / l_nom,
        sigma(SHAPE_NOM_RFOOT_OVER_R * SHAPE_REL[SHAPE_REL.len() - 1]) / l_nom,
    ]
}

/// The 150-cell coarse cross-check resolution (deeper inside the radiative-collapse stable
/// window; see the [`Shape1DRecord`] docs).
const SHAPE_1D_COARSE_CELLS: usize = 150;

/// Run one 1D shape row: the equilibrium coupled bounce (production config, `wall = None`) at the
/// sample's Σ-contract `(ρ, L)`, with the taper roles scaling `ρ` at fixed `L` (Σ moves through
/// the density, exactly as the radius-dependent profile does) — at both validity-protocol
/// resolutions, plus the EOS-only smoothness reference.
fn run_shape_1d_case(
    v: f64,
    sample: &ShapeSample,
    role: &str,
    sigma_factor: f64,
    table: &Table,
) -> Shape1DRecord {
    let (rho0, length) = shape_rho_length(sample.r_foot_over_r, sample.l_over_d);
    let rho = rho0 * sigma_factor;
    let cfg = Config {
        v,
        length,
        ..Config::production()
    };
    let slug = |cells: usize| {
        Tube::slug_si(
            cells,
            rho,
            cfg.v,
            cfg.length,
            cfg.t0,
            TableEos::new(table.clone()),
            Viscosity::VON_NEUMANN_RICHTMYER,
        )
    };
    let fine = CoupledBounce::new(slug(cfg.gas_cells), None, cfg.consts, cfg.limiter).run();
    let coarse =
        CoupledBounce::new(slug(SHAPE_1D_COARSE_CELLS), None, cfg.consts, cfg.limiter).run();
    let eos_only = slug(cfg.gas_cells).run_bounce();
    Shape1DRecord {
        v,
        axis: sample.axis.to_string(),
        sigma_role: role.to_string(),
        r_foot_over_r: sample.r_foot_over_r,
        l_over_d: sample.l_over_d,
        taper_frac: sample.taper_frac,
        sigma: rho * length,
        rho_impact: rho,
        length,
        e_eff: fine.bounce.e_eff,
        e_eff_coarse: coarse.bounce.e_eff,
        e_eff_eos: eos_only.e_eff,
        peak_wall_pressure: fine.bounce.peak_wall_pressure,
        peak_wall_pressure_coarse: coarse.bounce.peak_wall_pressure,
        incident_momentum: fine.bounce.incident_momentum,
        wall_impulse: fine.bounce.wall_impulse,
    }
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

/// Write `rows` to `result_path` (ADR-0019 JSONL), echo each row with `echo`, then print the
/// "wrote N rows" trailer — the write + per-row echo + summary shape every scenario command
/// shares; only `echo`'s row format differs, since each scenario's JSONL schema differs.
fn emit_scenario<T: Serialize>(
    result_path: &str,
    label: &str,
    rows: &[T],
    mut echo: impl FnMut(&T),
) -> Result<(), Box<dyn std::error::Error>> {
    write_rows(result_path, rows)?;
    for r in rows {
        echo(r);
    }
    println!("rust: wrote {} {label} rows -> {result_path}", rows.len());
    Ok(())
}

/// `--frozen-probe`: run `v_grid × rho_grid` EOS-only on `table_path` and record each case's
/// mass-weighted turnaround state; the Python table generator freezes the composition there.
/// Parameterized so the transitional (ADR-0026) and 69 km/s Jupiter passes share it.
fn frozen_probe_mode(
    v_grid: &[f64],
    rho_grid: &[f64],
    table_path: &str,
    base: &Config,
    result_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(table_path)?;
    let rows = run_sweep_frozen_probe(v_grid, rho_grid, &table, base);
    write_rows(result_path, &rows)?;
    for r in &rows {
        println!(
            "rust: v={:.0} rho={:.2} -> turnaround rho*={:.3} T*={:.0} K (e_eff_eq={:.4})",
            r.v, r.rho_impact, r.rho_star, r.t_star, r.e_eff_eq,
        );
    }
    println!("rust: wrote {} probe rows -> {result_path}", rows.len());
    Ok(())
}

/// `--frozen`: the three-curve frozen-recombination bounding sweep (equilibrium vs
/// sudden-freeze-at-turnaround vs pure-H2O-no-chemistry). Per-case frozen tables (in `table_dir`)
/// are loaded serially up front (cheap next to the bounces themselves). Parameterized so the
/// transitional and 69 km/s Jupiter passes share it.
fn frozen_sweep_mode(
    v_grid: &[f64],
    rho_grid: &[f64],
    table_path: &str,
    table_dir: &str,
    base: &Config,
    result_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(table_path)?;
    let h2o_tbl = Table::load(format!("{table_dir}/h2o.json"))?;
    let mut cases = Vec::new();
    for &v in v_grid {
        for &rho in rho_grid {
            cases.push((v, rho, Table::load(frozen_table_path(table_dir, v, rho))?));
        }
    }
    let rows: Vec<FrozenRecord> =
        par_map_with_progress("frozen", &cases, |(v, rho, frozen_tbl)| {
            run_one_frozen(*v, *rho, &table, frozen_tbl, &h2o_tbl, base)
        });
    write_rows(result_path, &rows)?;
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
    println!("rust: wrote {} frozen rows -> {result_path}", rows.len());
    Ok(())
}

// Transitional anchor (ADR-0012): sweep V_GRID × RHO_GRID with the high-v table, emitting the
// EOS-only curve and the radiation-on comparison curve into two fixed files.
fn cmd_transitional(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
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
    Ok(())
}

// Jupiter 69 km/s freeze-timing bracket (ADR-0026 instrument at the L=12 m realistic-cloud
// anchor, on the extended-grid table). Needs `make tables-jupiter` / `tables-frozen-jupiter`
// first. (`--frozen-probe-jupiter` / `--frozen-jupiter`.)
fn cmd_frozen_probe_jupiter(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = Config {
        v: JUP_V,
        length: JUP_FROZEN_LENGTH,
        ..Config::production()
    };
    frozen_probe_mode(
        &[JUP_V],
        &JUP_RHO,
        TABLE_PATH_JUPITER,
        &base,
        RESULT_PATH_FROZEN_PROBE_JUPITER,
    )
}

fn cmd_frozen_jupiter(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = Config {
        v: JUP_V,
        length: JUP_FROZEN_LENGTH,
        ..Config::production()
    };
    frozen_sweep_mode(
        &[JUP_V],
        &JUP_RHO,
        TABLE_PATH_JUPITER,
        TABLE_DIR_FROZEN_JUPITER,
        &base,
        RESULT_PATH_FROZEN_JUPITER,
    )
}

// Heavy-plate 16–28 km/s freeze-timing bracket (ADR-0026 instrument at the 16 / 22 / 28 km/s
// anchors, on the reused Jupiter extended-grid table). Needs `make tables-jupiter` /
// `tables-frozen-heavyplate` first. (`base.v` is overridden per grid velocity inside the mode.)
fn cmd_probe_heavyplate_diag(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = Config {
        v: HEAVY_DIAG_V[0],
        length: HEAVY_LENGTH,
        ..Config::production()
    };
    frozen_probe_mode(
        &HEAVY_DIAG_V,
        &HEAVY_SPOT_RHO,
        TABLE_PATH_JUPITER,
        &base,
        RESULT_PATH_PROBE_HEAVYPLATE_DIAG,
    )
}

fn cmd_frozen_probe_heavyplate(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = Config {
        v: HEAVY_V_ANCHORS[0],
        length: HEAVY_LENGTH,
        ..Config::production()
    };
    frozen_probe_mode(
        &HEAVY_V_ANCHORS,
        &HEAVY_RHO,
        TABLE_PATH_JUPITER,
        &base,
        RESULT_PATH_FROZEN_PROBE_HEAVYPLATE,
    )
}

fn cmd_frozen_heavyplate(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base = Config {
        v: HEAVY_V_ANCHORS[0],
        length: HEAVY_LENGTH,
        ..Config::production()
    };
    frozen_sweep_mode(
        &HEAVY_V_ANCHORS,
        &HEAVY_RHO,
        TABLE_PATH_JUPITER,
        TABLE_DIR_FROZEN_HEAVYPLATE,
        &base,
        RESULT_PATH_FROZEN_HEAVYPLATE,
    )
}

// Pulse-shape frozen spot-check (design §13): the three dip-anchor Σ points through the standard
// two-stage frozen pipeline (probe -> Python tables -> three-curve run), at the nominal
// Σ-contract length. Needs `make tables` / `tables-frozen-shape` first.
fn cmd_frozen_probe_shape(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let (_, l_nom) = shape_rho_length(SHAPE_NOM_RFOOT_OVER_R, SHAPE_NOM_L_OVER_D);
    let base = Config {
        v: SHAPE_V[0],
        length: l_nom,
        ..Config::production()
    };
    frozen_probe_mode(
        &[SHAPE_V[0]],
        &shape_frozen_rho_grid(),
        TABLE_PATH,
        &base,
        RESULT_PATH_FROZEN_PROBE_SHAPE,
    )
}

fn cmd_frozen_shape(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let (_, l_nom) = shape_rho_length(SHAPE_NOM_RFOOT_OVER_R, SHAPE_NOM_L_OVER_D);
    let base = Config {
        v: SHAPE_V[0],
        length: l_nom,
        ..Config::production()
    };
    frozen_sweep_mode(
        &[SHAPE_V[0]],
        &shape_frozen_rho_grid(),
        TABLE_PATH,
        TABLE_DIR_FROZEN_SHAPE,
        &base,
        RESULT_PATH_FROZEN_SHAPE,
    )
}

// Frozen-recombination probe: record each transitional case's turnaround state, from which the
// Python side generates the per-case frozen-composition tables.
fn cmd_frozen_probe(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    frozen_probe_mode(
        &V_GRID,
        &RHO_GRID,
        TABLE_PATH,
        &Config::production(),
        RESULT_PATH_FROZEN_PROBE,
    )
}

// Frozen-recombination bounding sweep: equilibrium vs sudden-freeze-at-turnaround vs
// pure-H2O-no-chemistry, per transitional case. Needs the per-case frozen tables (make
// tables-frozen).
fn cmd_frozen(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    frozen_sweep_mode(
        &V_GRID,
        &RHO_GRID,
        TABLE_PATH,
        TABLE_DIR_FROZEN,
        &Config::production(),
        RESULT_PATH_FROZEN,
    )
}

// Ablating-wall recovery sweep (Rung E, ADR-0014): the rigid floor vs the shielding+injection
// ablating wall over (v × ρ × opacity-scale × Q*), reporting recovery as a τ-bracket.
fn cmd_ablating(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let base_tbl = base_table();
    let cfg = Config {
        gas_cells: ABL_GAS_CELLS,
        ..Config::production()
    };
    let rows = run_ablating_sweep(&cfg, &base_tbl);
    emit_scenario(RESULT_PATH_ABLATING, "ablating", &rows, |r| {
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
    })
}

// Jupiter-retrograde 69 km/s scenario sweep: e_eff(rho × length × opacity-scale) with the
// extended-grid table (multi-stage O ladder). Needs `make tables-jupiter` first.
fn cmd_jupiter(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH_JUPITER)?;
    let rows = run_jupiter_sweep(&table);
    emit_scenario(RESULT_PATH_JUPITER, "jupiter", &rows, |r| {
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
    })
}

// Heavy-plate 16–28 km/s scenario sweep: the headline e_eff(v × ρ) grid + the L-sensitivity spot
// rows + the opacity τ-check, all on the reused Jupiter extended-grid table. Needs
// `make tables-jupiter` first.
fn cmd_transport_check(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH_JUPITER)?;
    let rows = run_transport_check_sweep(&table);
    emit_scenario(RESULT_PATH_TRANSPORT, "transport-check", &rows, |r| {
        println!(
            "rust: v={:>5.0} rho={:.3} -> tau_w={:.3e} bias_R={:+.2}% bias_P={:+.2}% \
             dtauR={:.2e} dtauP={:.2e} fld/BB={:.2e} escape={:.3e} ({:.4}% of KE) e_eff={:.4}",
            r.v,
            r.rho_impact,
            r.flux_weighted_optical_depth,
            100.0 * r.relative_bias,
            100.0 * r.relative_bias_planck,
            r.flux_weighted_surface_cell_depth,
            r.flux_weighted_surface_cell_depth_planck,
            r.flux_weighted_fld_over_blackbody,
            r.loss_escape_space,
            100.0 * r.escape_share_of_ke,
            r.e_eff,
        );
    })
}

fn cmd_mesh_convergence(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH_JUPITER)?;
    let rows = run_mesh_convergence_sweep(&table);
    emit_scenario(RESULT_PATH_MESH, "mesh-convergence", &rows, |r| {
        println!(
            "rust: v={:>5.0} rho={:.3} cells={:>5} -> e_eff={:.6} peak_p={:.4e} escape={:.4e} conv={}",
            r.v,
            r.rho_impact,
            r.gas_cells,
            r.e_eff,
            r.peak_wall_pressure,
            r.loss_escape_space,
            r.converged,
        );
    })
}

fn cmd_transport_resolution(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH_JUPITER)?;
    let rows = run_transport_resolution_sweep(&table);
    emit_scenario(
        RESULT_PATH_TRANSPORT_RES,
        "transport-resolution",
        &rows,
        |r| {
            println!(
                "rust: v={:>5.0} rho={:.3} cells={:>5} -> dtauP={:.3e} bias_R={:+.2}% bias_P={:+.2}% \
             e_eff={:.6} escape={:.4}% of KE",
                r.v,
                r.rho_impact,
                r.gas_cells,
                r.flux_weighted_surface_cell_depth_planck,
                100.0 * r.relative_bias,
                100.0 * r.relative_bias_planck,
                r.e_eff,
                100.0 * r.escape_share_of_ke,
            );
        },
    )
}

fn cmd_heavyplate(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH_JUPITER)?;
    let rows = run_heavyplate_sweep(&table);
    emit_scenario(RESULT_PATH_HEAVYPLATE, "heavy-plate", &rows, |r| {
        println!(
            "rust: v={:>5.0} rho={:.3} L={:>4.1} scale={:>5.2} -> e_eff={:.4} peak_p={:.3e} J={:.3e}",
            r.v,
            r.rho_impact,
            r.length,
            r.opacity_scale,
            r.e_eff,
            r.peak_wall_pressure,
            r.wall_impulse,
        );
    })
}

// Pulse-shape sensitivity sweep (design §13, ADR-0028): the fixed-grid 2D box + the fresh
// Σ-contract 1D runs, into two JSONL files. Needs `make tables` first (the 1D arm). `--1d-only`
// skips the (slow, unchanged) 2D box — a dev shortcut for 1D-schema re-runs.
fn cmd_shape(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let table = Table::load(TABLE_PATH)?;
    if !args.iter().any(|a| a == "--1d-only") {
        let rows2d = run_shape_2d_sweep();
        emit_scenario(RESULT_PATH_SHAPE_2D, "2D", &rows2d, |r| {
            println!(
                "rust: [{}] d/D={:.2} rf/R={:.2} L/D={:.3} taper={:.2} alpha={:.3} res={:.1} -> eta={:.4}",
                r.axis,
                r.d_over_d,
                r.r_foot_over_r,
                r.l_over_d,
                r.taper_frac,
                r.alpha_div,
                r.resolution_scale,
                r.eta_capture,
            );
        })?;
    }
    let cases = shape_1d_cases();
    let rows1d: Vec<Shape1DRecord> =
        par_map_with_progress("shape-1d", &cases, |(v, s, role, fac)| {
            run_shape_1d_case(*v, s, role, *fac, &table)
        });
    emit_scenario(RESULT_PATH_SHAPE_1D, "1D", &rows1d, |r| {
        println!(
            "rust: [{}/{}] v={:.0} Sigma={:.3} rho={:.3} L={:.2} -> e_eff={:.4} (coarse {:.4}, eos {:.4}) peak_p={:.3e}",
            r.axis,
            r.sigma_role,
            r.v,
            r.sigma,
            r.rho_impact,
            r.length,
            r.e_eff,
            r.e_eff_coarse,
            r.e_eff_eos,
            r.peak_wall_pressure,
        );
    })
}

/// One domain-convergence row: `eta_capture` at a given radial margin, at the corner where the
/// margin is tightest (`r_foot/R = 1.0`, cloud edge flush with the plate rim).
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
struct HeadroomRecord {
    d_over_d: f64,
    l_over_d: f64,
    r_foot_over_r: f64,
    mach: f64,
    /// `r_max / r_plate`.
    rim_headroom: f64,
    eta_capture: f64,
}

// Design §7 specified `r_foot/R` over 0.3–1.0; the implemented grid stopped at 0.7 and
// `contour::R_FOOT_BOX` inherited that as its upper edge. This completes the specified range at the
// same M = 40 anchor, and measures the domain-convergence risk that comes with it.
fn cmd_geometry_wide(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let cfg = GeoConfig::production();
    let cases: Vec<(f64, f64, f64)> = GEO_RFOOT_WIDE
        .iter()
        .flat_map(|&rf| {
            GEO_L_OVER_D
                .iter()
                .flat_map(move |&ld| GEO_D_OVER_D.iter().map(move |&dd| (dd, ld, rf)))
        })
        .collect();
    let rows: Vec<GeoRecord> = par_map_with_progress("geometry-wide", &cases, |&(dd, ld, rf)| {
        run_eta_case(dd, ld, rf, 40.0, &cfg)
    });
    emit_scenario(
        RESULT_PATH_GEOMETRY_WIDE,
        "wide-footprint geometry",
        &rows,
        |r| {
            println!(
                "rust: d/D={:.2} L/D={:.2} r_foot/R={:.2} M=40 -> eta_capture={:.4}",
                r.d_over_d, r.l_over_d, r.r_foot_over_r, r.eta_capture,
            );
        },
    )?;

    // Domain convergence at the tightest corner: only 0.4 r_foot of escape room past the cloud edge
    // at the default margin, against 3.7 at r_foot/R = 0.3. Run the headline curvature across L/D.
    let checks: Vec<(f64, f64, f64)> = GEO_RFOOT_WIDE
        .iter()
        .flat_map(|&rf| {
            GEO_L_OVER_D
                .iter()
                .flat_map(move |&ld| GEO_RIM_HEADROOM_CHECK.iter().map(move |&h| (rf, ld, h)))
        })
        .collect();
    let head: Vec<HeadroomRecord> =
        par_map_with_progress("geometry-wide-headroom", &checks, |&(rf, ld, h)| {
            let r = run_eta_case_with_headroom(0.10, ld, rf, 40.0, &cfg, h);
            HeadroomRecord {
                d_over_d: r.d_over_d,
                l_over_d: r.l_over_d,
                r_foot_over_r: r.r_foot_over_r,
                mach: r.mach,
                rim_headroom: h,
                eta_capture: r.eta_capture,
            }
        });
    emit_scenario(
        RESULT_PATH_GEOMETRY_WIDE_HEADROOM,
        "wide-footprint domain check",
        &head,
        |r| {
            println!(
                "rust: L/D={:.2} r_foot/R={:.2} r_max/r_plate={:.1} -> eta_capture={:.4}",
                r.l_over_d, r.r_foot_over_r, r.rim_headroom, r.eta_capture,
            );
        },
    )
}

// High-Mach spot check for the Jupiter scenario: the full geometry grid at M = 40 (the
// strong-shock plateau check past the production anchors 10/20).
fn cmd_geometry_m40(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let cfg = GeoConfig::production();
    let cases: Vec<(f64, f64, f64)> = GEO_RFOOT_OVER_R
        .iter()
        .flat_map(|&rf| {
            GEO_L_OVER_D
                .iter()
                .flat_map(move |&ld| GEO_D_OVER_D.iter().map(move |&dd| (dd, ld, rf)))
        })
        .collect();
    let rows: Vec<GeoRecord> = par_map_with_progress("geometry-m40", &cases, |&(dd, ld, rf)| {
        run_eta_case(dd, ld, rf, 40.0, &cfg)
    });
    emit_scenario(RESULT_PATH_GEOMETRY_M40, "M=40 geometry", &rows, |r| {
        println!(
            "rust: d/D={:.2} L/D={:.2} r_foot/R={:.2} M=40 -> eta_capture={:.4}",
            r.d_over_d, r.l_over_d, r.r_foot_over_r, r.eta_capture,
        );
    })
}

// Geometry sweep (Rung D follow-on): eta_capture(curvature × L/D × r_foot/R) from the euler2d
// kernel; no EOS/opacity table needed (radiation-free, effective-γ, ADR-0008).
fn cmd_geometry(_args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let rows = run_geometry_sweep(&GeoConfig::production());
    emit_scenario(RESULT_PATH_GEOMETRY, "geometry", &rows, |r| {
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
    })
}

// The default scenario: `--lowv` selects the 3.2 km/s condensing anchor (Rung C); otherwise the
// 16 km/s high-v pass. Optional positional `[table] [result]` override the high-v/low-v defaults
// (the B5d-3 opacity scan sweeps each scaled table into its own JSONL).
fn cmd_baseline(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
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

    let label = if lowv { "low-v" } else { "baseline" };
    emit_scenario(result_path, label, &records, |r| {
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
    })
}

/// A scenario command: parse whatever it needs from `args`, run its sweep, and write its JSONL.
type CmdFn = fn(&[String]) -> Result<(), Box<dyn std::error::Error>>;

/// Flag → handler. `main` looks up the first matching flag and dispatches to it (the exact-string
/// flags are mutually exclusive in practice, so table order does not affect behavior); anything
/// unmatched falls through to [`cmd_baseline`] (the `--lowv` / positional-override default).
const SCENARIOS: &[(&str, CmdFn)] = &[
    ("--transitional", cmd_transitional),
    ("--frozen-probe-jupiter", cmd_frozen_probe_jupiter),
    ("--frozen-jupiter", cmd_frozen_jupiter),
    ("--probe-heavyplate-diag", cmd_probe_heavyplate_diag),
    ("--frozen-probe-heavyplate", cmd_frozen_probe_heavyplate),
    ("--frozen-heavyplate", cmd_frozen_heavyplate),
    ("--frozen-probe-shape", cmd_frozen_probe_shape),
    ("--frozen-shape", cmd_frozen_shape),
    ("--frozen-probe", cmd_frozen_probe),
    ("--frozen", cmd_frozen),
    ("--ablating", cmd_ablating),
    ("--mesh-convergence", cmd_mesh_convergence),
    ("--transport-resolution", cmd_transport_resolution),
    ("--transport-check", cmd_transport_check),
    ("--spray-2d-lip", cmd_spray_2d_lip),
    ("--plug-2d", cmd_plug_2d),
    ("--spray-2d-pearl", cmd_spray_2d_pearl),
    ("--spray-2d-skirt", cmd_spray_2d_skirt),
    ("--spray-2d-shape", cmd_spray_2d_shape),
    ("--spray-k", cmd_spray_k),
    ("--spray-film", cmd_spray_film),
    ("--spray-2d-cup", cmd_spray_2d_cup),
    ("--spray-2d-standoff", cmd_spray_2d_standoff),
    ("--vessel-blast", cmd_vessel_blast),
    ("--chamber-blast", cmd_chamber_blast),
    ("--spray-face-history", cmd_spray_face_history),
    ("--spray-standoff", cmd_spray_standoff),
    ("--spray-2d", cmd_spray_2d),
    ("--spray-levers-shielded", cmd_spray_levers_shielded),
    ("--spray-levers-convergence", cmd_spray_levers_convergence),
    ("--spray-levers", cmd_spray_levers),
    ("--spray-mixing-convergence", cmd_spray_mixing_convergence),
    ("--spray-mixing", cmd_spray_mixing),
    ("--spray-ablating", cmd_spray_ablating),
    ("--spray", cmd_spray),
    ("--jupiter", cmd_jupiter),
    ("--heavyplate", cmd_heavyplate),
    ("--shape", cmd_shape),
    ("--geometry-m40", cmd_geometry_m40),
    ("--geometry-wide", cmd_geometry_wide),
    ("--geometry", cmd_geometry),
];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    for &(flag, handler) in SCENARIOS {
        if args.iter().any(|a| a == flag) {
            return handler(&args);
        }
    }
    cmd_baseline(&args)
}

#[cfg(test)]
mod tests {
    use super::{
        ABL_KAPPA_VAPOR, ABL_Q_STAR, AblatingRecord, Config, GeoConfig, GeoRecord, HEAVY_DIAG_V,
        HEAVY_L_SPOT, HEAVY_LENGTH, HEAVY_N_V, HEAVY_RHO, HEAVY_SPOT_RHO, HEAVY_TAU_SCALES,
        HEAVY_V_ANCHORS, HeavyPlateRecord, LowvConfig, Record, SHAPE_ALPHA, SHAPE_NOM_L_OVER_D,
        SHAPE_NOM_RFOOT_OVER_R, SHAPE_REL, SHAPE_TAPER, SHAPE_V, Shape1DRecord, ShapeRecord,
        TABLE_DIR_FROZEN, TABLE_DIR_FROZEN_HEAVYPLATE, TABLE_DIR_FROZEN_JUPITER, frozen_table_path,
        heavy_v_grid, heavyplate_cases, run_ablating_case, run_eta_case, run_one, run_one_frozen,
        run_sweep, run_sweep_frozen_probe, run_sweep_lowv, run_sweep_transitional, shape_1d_cases,
        shape_2d_cases, shape_frozen_rho_grid, shape_nominal, shape_rho_length, shape_samples,
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
            frozen_table_path(TABLE_DIR_FROZEN, 5_000.0, 0.16),
            "data/tables/frozen/v05000_rho0.16.json"
        );
        assert_eq!(
            frozen_table_path(TABLE_DIR_FROZEN, 16_000.0, 0.64),
            "data/tables/frozen/v16000_rho0.64.json"
        );
        // The 69 km/s Jupiter freeze bracket shares the contract on its own directory.
        assert_eq!(
            frozen_table_path(TABLE_DIR_FROZEN_JUPITER, 69_000.0, 0.07),
            "data/tables/frozen_jupiter/v69000_rho0.07.json"
        );
        // The heavy-plate freeze bracket, likewise, on its own directory (per-anchor velocities).
        assert_eq!(
            frozen_table_path(TABLE_DIR_FROZEN_HEAVYPLATE, 22_000.0, 0.30),
            "data/tables/frozen_heavyplate/v22000_rho0.30.json"
        );
    }

    /// The heavy-plate velocity grid is the **two-leg 29-anchor 16–63 km/s** grid (Q17): a coarse
    /// 3 km/s reconnaissance leg over 16–43, where the existing 16–28 study already showed `e_eff`
    /// varying smoothly, and a fine 1 km/s leg over 45–63, the extension's own range, where the
    /// curve has never been sampled and the `τ` structure is unknown.
    ///
    /// The 43 → 45 gap is deliberate, not an off-by-one: it is the seam between the two legs, and
    /// spacing there is 2 km/s — finer than the coarse leg, coarser than the fine one — so no
    /// resolution is lost across the join.
    ///
    /// The three freeze-bracket anchors (16 / 22 / 28 km/s) must land on the grid, since Q4 runs
    /// the full probe → frozen-table → three-curve pipeline at exactly those velocities.
    #[test]
    fn heavy_v_grid_is_29_anchors_over_two_legs() {
        let grid = heavy_v_grid();
        assert_eq!(grid.len(), HEAVY_N_V);
        assert_eq!(grid.len(), 29);
        assert!((grid[0] - 16_000.0).abs() < 1e-9);
        assert!((grid[grid.len() - 1] - 63_000.0).abs() < 1e-9);

        let coarse: Vec<f64> = grid.iter().copied().filter(|&v| v <= 43_000.0).collect();
        let fine: Vec<f64> = grid.iter().copied().filter(|&v| v >= 45_000.0).collect();
        assert_eq!(coarse.len(), 10, "16–43 km/s at 3 km/s is 10 anchors");
        assert_eq!(fine.len(), 19, "45–63 km/s at 1 km/s is 19 anchors");

        for pair in coarse.windows(2) {
            assert!(
                (pair[1] - pair[0] - 3_000.0).abs() < 1e-9,
                "coarse leg must step 3 km/s: {} -> {}",
                pair[0],
                pair[1]
            );
        }
        for pair in fine.windows(2) {
            assert!(
                (pair[1] - pair[0] - 1_000.0).abs() < 1e-9,
                "fine leg must step 1 km/s: {} -> {}",
                pair[0],
                pair[1]
            );
        }
        // Strictly increasing overall, and the seam is the 2 km/s join.
        assert!(grid.windows(2).all(|w| w[1] > w[0]));
        assert!((fine[0] - coarse[coarse.len() - 1] - 2_000.0).abs() < 1e-9);

        for &a in &HEAVY_V_ANCHORS {
            assert!(
                grid.iter().any(|&v| (v - a).abs() < 1e-9),
                "bracket anchor {a} not on the swept grid"
            );
        }
    }

    /// The heavy-plate case list: the headline grid plus both diagnostic families, with **no
    /// duplicates**.
    ///
    /// The duplicate check is the substantive one. A case listed twice runs twice, writes two
    /// identical rows, and then silently double-weights that state in every average the Python
    /// frontier takes — a corruption that produces confident, wrong numbers and shows up nowhere in
    /// the output. It is easy to reintroduce, because the diagnostic families deliberately overlap
    /// the headline grid in `(v, ρ)` and are kept distinct only by `L` or `opacity_scale`.
    #[test]
    fn heavyplate_case_list_covers_headline_and_diagnostics_without_duplicates() {
        let cases = heavyplate_cases();
        let headline = HEAVY_N_V * HEAVY_RHO.len();
        let l_spot = HEAVY_DIAG_V.len() * HEAVY_SPOT_RHO.len() * HEAVY_L_SPOT.len();
        let tau = HEAVY_DIAG_V.len() * HEAVY_SPOT_RHO.len() * HEAVY_TAU_SCALES.len();
        assert_eq!(cases.len(), headline + l_spot + tau);
        assert_eq!(cases.len(), 348 + 48 + 96);

        let mut keys: Vec<(u64, u64, u64, u64)> = cases
            .iter()
            .map(|&(v, r, l, s)| (v.to_bits(), r.to_bits(), l.to_bits(), s.to_bits()))
            .collect();
        keys.sort_unstable();
        let before = keys.len();
        keys.dedup();
        assert_eq!(before, keys.len(), "duplicate case in the heavy-plate list");

        // Every headline (v, ρ) is present exactly once at the reference L and κ = 1.
        for &v in &heavy_v_grid() {
            for &rho in &HEAVY_RHO {
                let n = cases
                    .iter()
                    .filter(|&&(cv, cr, cl, cs)| {
                        cv == v && cr == rho && cl == HEAVY_LENGTH && cs == 1.0
                    })
                    .count();
                assert_eq!(n, 1, "headline ({v}, {rho}) appears {n} times");
            }
        }

        // The diagnostic densities must be drawn from the headline grid, or the opacity bracket
        // cannot pair a τ-check row with its κ = 1 row and its probe turnaround state.
        for &rho in &HEAVY_SPOT_RHO {
            assert!(
                HEAVY_RHO.iter().any(|&r| (r - rho).abs() < 1e-12),
                "spot density {rho} is not on the headline grid"
            );
        }
        // Likewise the diagnostic velocities, and they must reach the top of the new range — the
        // whole point of widening them is to test the τ ≫ 1 claim where it is likeliest to fail.
        let grid = heavy_v_grid();
        for &v in &HEAVY_DIAG_V {
            assert!(
                grid.iter().any(|&g| (g - v).abs() < 1e-9),
                "diagnostic velocity {v} is not on the swept grid"
            );
        }
        assert!(
            HEAVY_DIAG_V.iter().any(|&v| v >= 63_000.0),
            "diagnostics must reach the top of the extended range"
        );
        // κ = 1 is never duplicated into the τ family; it comes from the headline grid.
        assert!(!HEAVY_TAU_SCALES.iter().any(|&s| (s - 1.0).abs() < 1e-12));
    }

    /// The Q11 density grid: 12 points, ascending, dense through the steep on-contour band
    /// 0.06–0.20 and coarse across the flat top, with only three dilute points below 0.06 —
    /// those exist to support Q15's validation points, not to resolve a region the contour spends
    /// time in.
    #[test]
    fn heavy_rho_grid_is_twelve_points_dense_on_contour() {
        assert_eq!(HEAVY_RHO.len(), 12);
        assert!(HEAVY_RHO.windows(2).all(|w| w[1] > w[0]), "must ascend");
        assert!((HEAVY_RHO[0] - 0.01).abs() < 1e-12);
        assert!((HEAVY_RHO[HEAVY_RHO.len() - 1] - 0.58).abs() < 1e-12);

        let on_contour = HEAVY_RHO
            .iter()
            .filter(|&&r| (0.06..=0.20).contains(&r))
            .count();
        assert_eq!(
            on_contour, 6,
            "steep on-contour band 0.06–0.20 gets 6 points"
        );
        assert_eq!(
            HEAVY_RHO.iter().filter(|&&r| r < 0.06).count(),
            3,
            "three dilute points, for Q15 validation only"
        );

        // Spacing widens monotonically above the on-contour band — the flat top is deliberately
        // coarse, and a regression that "helpfully" evened out the grid would cost run time in the
        // region the answer does not depend on.
        let top: Vec<f64> = HEAVY_RHO.iter().copied().filter(|&r| r >= 0.15).collect();
        let gaps: Vec<f64> = top.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(
            gaps.windows(2).all(|g| g[1] > g[0]),
            "flat-top spacing should widen: {gaps:?}"
        );
    }

    /// A `HeavyPlateRecord` round-trips through the JSONL boundary (ADR-0019): the Python heavyplate
    /// reader sees exactly the fields written — the case axes, the restitution, the facesheet peak,
    /// the areal-impulse pair (ADR-0027 structural companion), and the radiative loss split.
    #[allow(clippy::float_cmp)] // round-number inputs survive the decimal text verbatim
    #[test]
    fn heavyplate_record_jsonl_roundtrip() {
        let rec = HeavyPlateRecord {
            v: 22_000.0,
            rho_impact: 0.08,
            length: 10.0,
            opacity_scale: 1.0,
            e_eff: 0.74,
            peak_wall_pressure: 4.6e7,
            incident_momentum: 1.76e4,
            wall_impulse: 3.06e4,
            loss_radiative_wall: 1.0e6,
            loss_escape_space: 1.0e5,
            converged: true,
        };
        let line = serde_json::to_string(&rec).unwrap();
        let back: HeavyPlateRecord = serde_json::from_str(&line).unwrap();
        assert_eq!(back, rec);
        for key in [
            "v",
            "rho_impact",
            "length",
            "opacity_scale",
            "e_eff",
            "peak_wall_pressure",
            "incident_momentum",
            "wall_impulse",
            "loss_radiative_wall",
            "loss_escape_space",
            "converged",
        ] {
            assert!(line.contains(key), "missing field {key}");
        }
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

    /// The shape-box sample list (design §13): exactly one nominal, one-factor-at-a-time on every
    /// other sample (S is a per-axis derivative), and the expected per-axis counts (6 two-sided
    /// footprint/aspect samples each, 3 one-sided taper/divergence samples each).
    #[test]
    fn shape_samples_are_one_factor_at_a_time() {
        let samples = shape_samples();
        assert_eq!(samples.len(), 19);
        let nom = shape_nominal();
        let count = |axis: &str| samples.iter().filter(|s| s.axis == axis).count();
        assert_eq!(count("nominal"), 1);
        assert_eq!(count("r_foot_over_r"), 6);
        assert_eq!(count("l_over_d"), 6);
        assert_eq!(count("taper_frac"), 3);
        assert_eq!(count("alpha_div"), 3);
        for s in &samples {
            let moved = [
                (s.r_foot_over_r - nom.r_foot_over_r).abs() > 1e-12,
                (s.l_over_d - nom.l_over_d).abs() > 1e-12,
                (s.taper_frac - nom.taper_frac).abs() > 1e-12,
                (s.alpha_div - nom.alpha_div).abs() > 1e-12,
            ]
            .iter()
            .filter(|&&m| m)
            .count();
            assert_eq!(
                moved,
                usize::from(s.axis != "nominal"),
                "sample [{}] must perturb exactly its own axis",
                s.axis
            );
        }
    }

    /// The Rust Σ contract matches `puffsat.analysis.impact_density` (the cross-language contract
    /// the 1D arm's densities come from): `ρ = m/(2π·(L/D)·r_foot³)`, `L = 2·(L/D)·r_foot` at
    /// m = 25 kg, R = 5 m.
    #[test]
    fn shape_sigma_contract_matches_python_impact_density() {
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12 * b.abs();
        let (rho, l) = shape_rho_length(0.5, 0.3);
        assert!(close(rho, 0.848_826_363_156_775_2), "nominal rho {rho}");
        assert!(close(l, 1.5), "nominal length {l}");
        let (rho, l) = shape_rho_length(0.4, 0.3);
        assert!(close(rho, 1.657_863_990_540_576_5), "endpoint rho {rho}");
        assert!(close(l, 1.2), "endpoint length {l}");
    }

    /// The 2D case list: the full 19-sample box × both plates at base resolution, plus the
    /// 5-sample noise-floor/flag-refinement subset × both plates × the two refinements (≥ 3
    /// repeats, ADR-0028) — and only 7 unique column lengths (the footprint and aspect axes
    /// share `L = 0.6·rel`).
    #[test]
    fn shape_2d_case_list_covers_box_and_noise_subset() {
        let cases = shape_2d_cases();
        assert_eq!(cases.len(), 19 * 2 + 5 * 2 * 2);
        assert_eq!(cases.iter().filter(|(_, _, res)| *res == 0).count(), 38);
        let mut lengths: Vec<u64> = cases
            .iter()
            .map(|(s, _, _)| super::shape_cloud(s).1.to_bits())
            .collect();
        lengths.sort_unstable();
        lengths.dedup();
        assert_eq!(
            lengths.len(),
            7,
            "footprint and aspect axes share the length set"
        );
    }

    /// The 1D case list (per anchor): fresh Σ runs for the 13 footprint/aspect samples (incl.
    /// nominal), a mean-Σ run per taper sample, the Σ-hi/lo90 bound rows at the widest taper only,
    /// and *no* divergence rows (Σ untouched).
    #[test]
    fn shape_1d_case_list_covers_sigma_bookkeeping() {
        let cases = shape_1d_cases();
        assert_eq!(cases.len(), SHAPE_V.len() * (13 + 3 + 2));
        assert!(cases.iter().all(|(_, s, _, _)| s.axis != "alpha_div"));
        for &v in &SHAPE_V {
            #[allow(clippy::float_cmp)] // v passes through the case builder verbatim
            let at_v: Vec<_> = cases.iter().filter(|(cv, ..)| *cv == v).collect();
            assert_eq!(
                at_v.iter()
                    .filter(|(_, _, role, _)| *role == "taper_mean")
                    .count(),
                SHAPE_TAPER.len() - 1
            );
            for role in ["sigma_hi", "sigma_lo90"] {
                let bound: Vec<_> = at_v.iter().filter(|(_, _, r, _)| *r == role).collect();
                assert_eq!(bound.len(), 1, "{role} only at the widest taper");
                assert!((bound[0].1.taper_frac - SHAPE_TAPER[3]).abs() < 1e-12);
            }
        }
        // The Σ factors bracket the mean for the bound rows.
        let factor = |role: &str| {
            cases
                .iter()
                .find(|(_, _, r, _)| *r == role)
                .map(|(_, _, _, f)| *f)
                .unwrap()
        };
        assert!(factor("sigma_lo90") < factor("sigma_hi"));
    }

    /// A `ShapeRecord` and a `Shape1DRecord` round-trip through the JSONL boundary (ADR-0019):
    /// the Python shape reader sees exactly the fields written.
    #[allow(clippy::float_cmp)] // round-number inputs survive the decimal text verbatim
    #[test]
    fn shape_records_jsonl_roundtrip() {
        let rec2d = ShapeRecord {
            axis: "r_foot_over_r".to_string(),
            d_over_d: 0.10,
            r_foot_over_r: 0.55,
            l_over_d: 0.3,
            taper_frac: 0.0,
            alpha_div: 0.0,
            mach: 20.0,
            resolution_scale: 1.0,
            eta_capture: 0.98,
            restitution_free: 1.42,
            restitution_confined: 1.45,
            incident_momentum: 0.61,
            peak_local_pressure: 620.0,
        };
        let line = serde_json::to_string(&rec2d).unwrap();
        let back: ShapeRecord = serde_json::from_str(&line).unwrap();
        assert_eq!(back, rec2d);
        for key in [
            "axis",
            "resolution_scale",
            "eta_capture",
            "incident_momentum",
        ] {
            assert!(line.contains(key), "missing field {key}");
        }
        let rec1d = Shape1DRecord {
            v: 11_000.0,
            axis: "taper_frac".to_string(),
            sigma_role: "taper_mean".to_string(),
            r_foot_over_r: 0.5,
            l_over_d: 0.3,
            taper_frac: 0.3,
            sigma: 1.54,
            rho_impact: 1.027,
            length: 1.5,
            e_eff: 0.57,
            e_eff_coarse: 0.569,
            e_eff_eos: 0.576,
            peak_wall_pressure: 1.2e8,
            peak_wall_pressure_coarse: 1.19e8,
            incident_momentum: 1.7e4,
            wall_impulse: 2.7e4,
        };
        let line = serde_json::to_string(&rec1d).unwrap();
        let back: Shape1DRecord = serde_json::from_str(&line).unwrap();
        assert_eq!(back, rec1d);
        for key in [
            "sigma_role",
            "sigma",
            "rho_impact",
            "e_eff",
            "e_eff_coarse",
            "e_eff_eos",
        ] {
            assert!(line.contains(key), "missing field {key}");
        }
    }

    /// The frozen spot-check densities: three distinct Σ-equivalent points in descending order
    /// (tightest footprint, nominal, widest), the middle one the nominal Σ-contract density, and
    /// all distinct under the two-decimal frozen-table naming contract.
    #[test]
    fn shape_frozen_grid_is_three_distinct_descending_sigma_points() {
        let grid = shape_frozen_rho_grid();
        assert!(grid[0] > grid[1] && grid[1] > grid[2]);
        let (rho_nom, _) = shape_rho_length(SHAPE_NOM_RFOOT_OVER_R, SHAPE_NOM_L_OVER_D);
        assert!((grid[1] - rho_nom).abs() < 1e-12);
        let names: Vec<String> = grid.iter().map(|r| format!("{r:.2}")).collect();
        assert_eq!(names.len(), 3);
        assert!(names[0] != names[1] && names[1] != names[2]);
        // The endpoints follow Σ ∝ 1/r_foot²: (0.5/0.4)² and (0.5/0.6)² of nominal.
        assert!((grid[0] / grid[1] - 1.5625).abs() < 1e-9);
        assert!((grid[2] / grid[1] - (0.5f64 / 0.6).powi(2)).abs() < 1e-9);
        // Guard the one-sided axes' sampling extents the case builders rely on.
        assert!((SHAPE_REL[0] - 0.8).abs() < 1e-12 && (SHAPE_REL[6] - 1.2).abs() < 1e-12);
        assert!((SHAPE_ALPHA[3] - 0.1).abs() < 1e-12 && (SHAPE_NOM_L_OVER_D - 0.3).abs() < 1e-12);
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

    /// DIAGNOSTIC (ignored): verify the **radiative collapse fix** across the measured onset
    /// ladder at 16 km/s (2026-07-16 finding; see the [`Shape1DRecord`] docs). With the old
    /// ρ ≤ 20 kg/m³ table the coupled solve collapsed (run stalls mid-infall, unphysical `e_eff`)
    /// past a critical refinement that coarsened with `ρ`: ρ = 0.64 at 1200 cells, 0.849/0.99 at
    /// 600, 1.164 at 300, 1.658 at 200 — the radiatively-cooled wall cell compressed past the
    /// table ceiling, where the clamped `p(ρ)` no longer arrests the Lagrangian compression.
    /// With the node-preservingly extended table (ceiling 1000 kg/m³) the compression
    /// self-arrests, so every case below must complete on the same flat plateau as its coarse
    /// siblings (< ~0.002 spread per ρ; EOS-only reference `e_eff` 0.638–0.647 throughout).
    /// The `(ρ, L)` pairs are the shape box's Σ contract: `ρ = 0.849·rel⁻³`, `L = 1.5·rel`.
    /// Run with `cargo test -p sweep --release -- --ignored --nocapture diag_shape`.
    #[test]
    #[ignore = "diagnostic; needs data/tables/water.json"]
    fn diag_shape_high_sigma() {
        use hydro1d::eos::TableEos;
        use hydro1d::kernel::{CoupledBounce, Tube, Viscosity};
        let table = Table::load(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/tables/water.json"
        ))
        .unwrap();
        let consts = RadConstants {
            c: 2.997_924_58e8,
            a: 7.565_733e-16,
        };
        // (rho, L, cell ladder): each ladder spans the old collapse onset and one refinement past
        // it (2x), plus a coarse in-plateau anchor.
        let cases: [(f64, f64, &[usize]); 5] = [
            (0.638, 1.65, &[300, 600, 1200]),  // old onset 1200
            (0.849, 1.50, &[300, 600, 1200]),  // old onset 600
            (0.990, 1.425, &[300, 600, 1200]), // old onset 600
            (1.164, 1.35, &[150, 300, 600]),   // old onset 300
            (1.658, 1.20, &[100, 200, 400]),   // old onset 200
        ];
        for (rho, length, ladder) in cases {
            for &cells in ladder {
                let mk = || {
                    Tube::slug_si(
                        cells,
                        rho,
                        16_000.0,
                        length,
                        400.0,
                        TableEos::new(table.clone()),
                        Viscosity::VON_NEUMANN_RICHTMYER,
                    )
                };
                let eos_only = mk().run_bounce();
                let coupled =
                    CoupledBounce::new(mk(), None, consts, Limiter::LevermorePomraning).run();
                println!(
                    "rho={rho:.3} L={length:.3} cells={cells}: eos-only e_eff={:+.4} | coupled \
                     e_eff={:+.4} resid/inc={:+.3} loss1a={:.3e} loss1b={:.3e}",
                    eos_only.e_eff,
                    coupled.bounce.e_eff,
                    coupled.bounce.residual_momentum / coupled.bounce.incident_momentum,
                    coupled.loss_radiative_wall,
                    coupled.loss_escape_space,
                );
            }
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
