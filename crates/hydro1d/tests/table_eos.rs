//! B2 regression: the **tabulated** EOS path reproduces rung A's analytic ideal gas.
//!
//! A table that *encodes* `p = (γ−1) ρ e` (with heat capacity `c_v = 1`, so `e = T`,
//! `p = (γ−1) ρ T`, `c_s = √(γ(γ−1)T)`) is fed to the kernel through [`TableEos`], and the Sod
//! shock tube is run through it. Because the table's log-log interpolation is *exact* for a
//! power-law ideal gas and the `e → T` inversion is bisected to machine precision, the tabulated
//! solver must
//!
//! 1. match the analytic-EOS run cell-for-cell (the EOS abstraction + inversion add no error), and
//! 2. independently clear rung A's Sod acceptance bar (L1 density error `< 0.0015`).

use hydro1d::Primitive;
use hydro1d::eos::{IdealGas, TableEos};
use hydro1d::kernel::{Tube, Viscosity};
use hydro1d::riemann::solve;
use tables::Table;

const GAMMA: f64 = 1.4;
const CELLS: usize = 400;
const T_END: f64 = 0.2;

/// A table encoding the ideal gas on log-spaced grids spanning the Sod state space with margin.
/// Power laws in `(ρ, T)`, so log-log bilinear interpolation reproduces them exactly.
fn ideal_gas_table() -> Table {
    let n = 8;
    let rho_grid: Vec<f64> = (0..n)
        .map(|i| 1e-3 * 1e5_f64.powf(i as f64 / (n - 1) as f64)) // 1e-3 … 100
        .collect();
    let t_grid: Vec<f64> = (0..n)
        .map(|j| 1e-3 * 1e6_f64.powf(j as f64 / (n - 1) as f64)) // 1e-3 … 1000
        .collect();
    let mut p = Vec::new();
    let mut e = Vec::new();
    let mut cs = Vec::new();
    for &r in &rho_grid {
        for &t in &t_grid {
            p.push((GAMMA - 1.0) * r * t);
            e.push(t);
            cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
        }
    }
    let one = vec![1.0; n * n];
    let json = serde_json::json!({
        "rho_grid": rho_grid,
        "T_grid": t_grid,
        "shape": [n, n],
        "fields": { "p": p, "e": e, "c_s": cs, "kappa_rosseland": one, "kappa_planck": one },
    });
    Table::from_json(&json.to_string()).unwrap()
}

/// The standard Sod initial conditions on `x ∈ [0, 1]` (diaphragm at 0.5).
fn sod_ic() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let dx = 1.0 / CELLS as f64;
    let x: Vec<f64> = (0..=CELLS).map(|i| i as f64 * dx).collect();
    let mut rho = vec![0.0; CELLS];
    let mut pressure = vec![0.0; CELLS];
    for j in 0..CELLS {
        if (j as f64 + 0.5) * dx < 0.5 {
            rho[j] = 1.0;
            pressure[j] = 1.0;
        } else {
            rho[j] = 0.125;
            pressure[j] = 0.1;
        }
    }
    (x, rho, vec![0.0; CELLS], pressure)
}

#[test]
fn table_eos_solves_sod_like_ideal_gas() {
    let (x, rho, vel, p) = sod_ic();
    let visc = Viscosity::VON_NEUMANN_RICHTMYER;
    let mut ideal = Tube::with_eos(x.clone(), &rho, &vel, &p, IdealGas::new(GAMMA), visc);
    let mut table = Tube::with_eos(x, &rho, &vel, &p, TableEos::new(ideal_gas_table()), visc);
    ideal.run_to(T_END);
    table.run_to(T_END);

    // (a) the tabulated path matches the analytic-EOS path cell-for-cell.
    let max_diff = (0..CELLS)
        .map(|j| (table.density(j) - ideal.density(j)).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_diff < 1e-6,
        "table vs ideal-gas max |Δρ| = {max_diff:e}"
    );

    // (b) the tabulated path independently clears rung A's Sod acceptance bar.
    let exact = solve(
        Primitive::new(1.0, 0.0, 1.0),
        Primitive::new(0.125, 0.0, 0.1),
        GAMMA,
    );
    let l1: f64 = (0..CELLS)
        .map(|j| {
            let xi = (table.center(j) - 0.5) / T_END;
            (table.density(j) - exact.sample(xi).rho).abs()
        })
        .sum::<f64>()
        / CELLS as f64;
    assert!(l1 < 0.0015, "table-EOS Sod L1 density error = {l1}");
}
