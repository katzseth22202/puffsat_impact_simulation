//! Order-of-accuracy acceptance (ADR-0023): on a *smooth* problem the MUSCL-Hancock scheme must
//! converge at its formal (≈ 2nd-order) rate under grid refinement — the strongest single
//! correctness signal for the reconstruction (CLAUDE.md). Shock-capturing is tested separately
//! (Sod, Sedov, Noh).
//!
//! Test problem: a smooth entropy wave `ρ = 1 + 0.2 sin(2πz)` on uniform velocity `u_z = 1` and
//! uniform pressure (an exact contact wave that simply advects), periodic in `z`, uniform in `r`.
//! The exact solution at `t` is the initial profile shifted by `u_z t`. The L1 density error is
//! measured at three resolutions and the observed order read from successive halvings.

use euler2d::kernel::{Bc, Grid2D};
use euler2d::state::Prim;
use std::f64::consts::TAU;

const GAMMA: f64 = 1.4;
const T_END: f64 = 0.1;

fn rho0(z: f64) -> f64 {
    1.0 + 0.2 * (TAU * z).sin()
}

/// L1 density error after advecting the smooth entropy wave to `T_END` at `nz` axial cells.
fn l1_density_error(nz: usize) -> f64 {
    let nr = 4;
    let dz = 1.0 / nz as f64;
    let dr = 1.0 / nr as f64;
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.bc_zlo = Bc::Periodic;
    g.bc_zhi = Bc::Periodic;
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect;
    g.init(|iz, _ir| {
        let z = (iz as f64 + 0.5) * dz;
        Prim::new(rho0(z), 1.0, 0.0, 1.0)
    });
    g.run_to(T_END);

    let ir = nr / 2;
    let mut l1 = 0.0;
    for iz in 0..nz {
        let z = (iz as f64 + 0.5) * dz;
        let z_src = (z - T_END).rem_euclid(1.0); // exact: advected by u_z = 1, periodic on [0,1]
        l1 += (g.prim(iz, ir).rho - rho0(z_src)).abs() * dz;
    }
    l1
}

#[test]
fn second_order_on_smooth_advection() {
    let e1 = l1_density_error(32);
    let e2 = l1_density_error(64);
    let e3 = l1_density_error(128);
    let o1 = (e1 / e2).log2();
    let o2 = (e2 / e3).log2();
    // van Leer MUSCL-Hancock is formally 2nd order; the TVD limiter clips the smooth extrema, so
    // the L1 rate sits a touch below 2. Require clearly-better-than-first-order convergence.
    assert!(
        o1 > 1.6 && o2 > 1.6,
        "convergence orders {o1:.2}, {o2:.2} (L1 = {e1:.2e}, {e2:.2e}, {e3:.2e})"
    );
}
