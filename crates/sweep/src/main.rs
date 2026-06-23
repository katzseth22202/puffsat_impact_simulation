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
//! `loss_conductive = 0` pending that high-v transport data. (The low-v CoolProp table *does* carry
//! `k_gas`, so the low-v path activates the operator.)

use std::fs;
use std::io::Write as _;
use std::path::Path;

use hydro1d::conduction::Solid;
use hydro1d::eos::TableEos;
use hydro1d::kernel::{CondensingBounce, CoupledBounce, Tube, Viscosity};
use hydro1d::radiation::{Limiter, RadConstants};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use tables::Table;

const TABLE_PATH: &str = "data/tables/water.json";
const RESULT_PATH: &str = "data/results/sweep.jsonl";
const TABLE_PATH_LOWV: &str = "data/tables/water_lowv.json";
const RESULT_PATH_LOWV: &str = "data/results/sweep_lowv.jsonl";

/// The production `ρ_impact` grid [kg/m³] (design §3-4): the impact-density axis of `e_eff(ρ)`.
const RHO_GRID: [f64; 4] = [0.16, 0.32, 0.48, 0.64];

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
    /// Peak wall force during the bounce.
    peak_wall_force: f64,
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // `--lowv` selects the 3.2 km/s condensing anchor (Rung C); otherwise the 16 km/s high-v pass.
    // Optional positional `[table_path] [result_path]` override the defaults (the B5d-3 opacity scan
    // uses them to sweep each opacity-scaled table into its own JSONL).
    let args: Vec<String> = std::env::args().skip(1).collect();
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
    use super::{Config, LowvConfig, Record, run_one, run_sweep, run_sweep_lowv};
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
}
