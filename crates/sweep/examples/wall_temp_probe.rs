//! Diagnostic: the peak **plate-facing** (wall-cell) gas temperature over the production coupled
//! bounce, tabulated across the full v × ρ grid. This is the reflected-shock stagnation temperature
//! the facesheet actually sees — the ablation/thermal-load companion to `peak_wall_pressure`, and
//! notably hotter than the mass-weighted mean-freeze temperature the frozen probe logs (the wall
//! cell is doubly shocked). Mirrors `run_one` in the sweep binary exactly (same `water.json` table,
//! `t0`, `length`, cell count, radiation constants, limiter, and `wall = None`), so the numbers are
//! the production bounce, not a reduced model. Run from the workspace root:
//!
//!     cargo run --release -p sweep --example wall_temp_probe

use hydro1d::eos::TableEos;
use hydro1d::kernel::{CoupledBounce, Tube, Viscosity};
use hydro1d::radiation::{Limiter, RadConstants};
use tables::Table;

const RHO_GRID: [f64; 4] = [0.16, 0.32, 0.48, 0.64];
const V_GRID: [f64; 8] = [
    5_000.0, 6_000.0, 7_000.0, 8_000.0, 9_000.0, 11_000.0, 13_000.0, 16_000.0,
];
const T0: f64 = 400.0;
const LENGTH: f64 = 1.0;
const GAS_CELLS: usize = 300;
const CONSTS: RadConstants = RadConstants {
    c: 2.997_924_58e8,
    a: 7.565_733e-16,
};

fn main() {
    let table = Table::load("data/tables/water.json")
        .expect("load data/tables/water.json (run `make tables` first)");

    println!("     v[m/s]   rho    e_eff   p_peak[GPa]   T_wall_peak[K]");
    for &v in &V_GRID {
        for &rho in &RHO_GRID {
            let eos = TableEos::new(table.clone());
            let tube = Tube::slug_si(
                GAS_CELLS,
                rho,
                v,
                LENGTH,
                T0,
                eos,
                Viscosity::VON_NEUMANN_RICHTMYER,
            );
            let r = CoupledBounce::new(tube, None, CONSTS, Limiter::LevermorePomraning).run();
            println!(
                "{:>10.0} {:>5.2} {:>8.4} {:>12.4} {:>14.0}",
                v,
                rho,
                r.bounce.e_eff,
                r.bounce.peak_wall_pressure / 1e9,
                r.bounce.peak_wall_temperature,
            );
        }
    }
}
