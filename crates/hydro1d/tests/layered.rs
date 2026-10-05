//! A layered `TableEos`: each Lagrangian cell carries its own material (ADR-0055 step 2a).
//!
//! The spray plate's mixing runs put `PuffSat` water and argon in neighbouring cells. Each material
//! has its own table and its own energy reference, so every per-cell query must reach that cell's
//! table. The seams, agreed before these tests were written:
//!
//!  1. the layered EOS answers per cell exactly as the standalone material does;
//!  2. a tube seeded from per-cell pressures reports those pressures back;
//!  3. radiation reads each cell's own material: an isothermal A|B|A slab loses radiation exactly
//!     as an all-A slab does, because nothing passes between cells at one temperature;
//!  4. (the existing suite) single-material runs are unchanged.

use hydro1d::eos::{Eos, TableEos};
use hydro1d::kernel::{CoupledBounce, Tube, Viscosity};
use hydro1d::radiation::{Limiter, RadConstants};
use tables::Table;

const CONSTS: RadConstants = RadConstants { c: 3.0, a: 1.0 };
const KAPPA: f64 = 50.0;

/// Ideal monatomic gas of relative molar mass `mu` with an energy offset `e0`:
/// `p = ρT/μ`, `e = 1.5 T/μ + e0`. Power laws in (ρ, T), so log-log interpolation is exact.
fn ideal(mu: f64, e0: f64) -> TableEos {
    let n = 9;
    let rho_grid: Vec<f64> = (0..n)
        .map(|i| 1e-2 * 1e4_f64.powf(i as f64 / (n - 1) as f64))
        .collect();
    let t_grid: Vec<f64> = (0..n)
        .map(|j| 1e-2 * 1e4_f64.powf(j as f64 / (n - 1) as f64))
        .collect();
    let (mut p, mut e, mut cs) = (Vec::new(), Vec::new(), Vec::new());
    for &r in &rho_grid {
        for &t in &t_grid {
            p.push(r * t / mu);
            e.push(1.5 * t / mu + e0);
            cs.push((5.0 / 3.0 * t / mu).sqrt());
        }
    }
    let kap = vec![KAPPA; n * n];
    let json = serde_json::json!({
        "rho_grid": rho_grid,
        "T_grid": t_grid,
        "shape": [n, n],
        "fields": { "p": p, "e": e, "c_s": cs, "kappa_rosseland": kap, "kappa_planck": kap },
    });
    TableEos::new(Table::from_json(&json.to_string()).unwrap())
}

fn light() -> TableEos {
    ideal(1.0, 0.0)
}

/// Four times heavier, with an energy offset, so reading it through the light table goes wrong.
fn heavy() -> TableEos {
    ideal(4.0, 2.0)
}

/// An isothermal, pressure-matched stationary slab of `cells` cells on [0, 1]. `index[j]` picks
/// the material (0 light, 1 heavy); density scales with molar mass so pressure is uniform.
fn isothermal_slab(index: Vec<usize>, t0: f64) -> Tube<TableEos> {
    let cells = index.len();
    let x: Vec<f64> = (0..=cells).map(|i| i as f64 / cells as f64).collect();
    let rho: Vec<f64> = index
        .iter()
        .map(|&m| if m == 0 { 1.0 } else { 4.0 })
        .collect();
    let p = vec![t0; cells];
    let vel = vec![0.0; cells];
    let eos = TableEos::layered(vec![light(), heavy()], index);
    Tube::with_eos(x, &rho, &vel, &p, eos, Viscosity::VON_NEUMANN_RICHTMYER)
}

#[test]
fn layered_eos_answers_each_cell_as_its_own_material() {
    let eos = TableEos::layered(vec![light(), heavy()], vec![0, 1, 1, 0]);
    let (rho, e) = (2.0, 5.0);
    for (cell, standalone) in [(0, light()), (1, heavy()), (2, heavy()), (3, light())] {
        assert_eq!(
            eos.for_cell(cell).pressure(rho, e),
            standalone.pressure(rho, e)
        );
        assert_eq!(
            eos.for_cell(cell).temperature(rho, e),
            standalone.temperature(rho, e)
        );
    }
}

#[test]
fn layered_tube_reports_the_pressures_it_was_seeded_with() {
    let tube = isothermal_slab(vec![0, 0, 1, 1, 1, 0], 3.0);
    for j in 0..tube.cells() {
        let p = tube.pressure(j);
        assert!((p - 3.0).abs() < 1e-9 * 3.0, "cell {j}: p = {p}, want 3");
    }
}

#[test]
fn isothermal_layered_slab_radiates_as_its_outer_material_does() {
    let t0 = 2.0;
    let (cells, edge) = (60, 10);
    let layered_index: Vec<usize> = (0..cells)
        .map(|j| usize::from(j >= edge && j < cells - edge))
        .collect();
    let mut layered = CoupledBounce::new(
        isothermal_slab(layered_index, t0),
        None,
        CONSTS,
        Limiter::LevermorePomraning,
    );
    let mut uniform = CoupledBounce::new(
        isothermal_slab(vec![0; cells], t0),
        None,
        CONSTS,
        Limiter::LevermorePomraning,
    );
    // Seeded in equilibrium: e_rad = a T0^4 in every cell.
    let expected = CONSTS.a * t0.powi(4);
    assert!((layered.radiation_energy() - expected).abs() < 1e-9 * expected);

    let duration = 0.05;
    layered.advance_for(duration);
    uniform.advance_for(duration);
    let (l, u) = (layered.radiation_energy(), uniform.radiation_energy());
    assert!(u < expected, "the walls must leak some radiation");
    assert!(
        (l - u).abs() < 1e-6 * u,
        "layered {l} vs uniform {u}: interior cells exchanged radiation at one temperature"
    );
}
