//! Sod shock-tube acceptance test (ADR-0023): the 2D HLLC Godunov kernel must reproduce the exact
//! Riemann solution along the sweep direction, embedded in 2D.
//!
//! The tube varies only in `z`, is uniform in `r`, and has reflecting `r`-walls — so the r-sweep is
//! a no-op and every `r`-row must match the exact 1D self-similar solution. This pins the
//! z-direction HLLC flux + the conservative update + the transverse (r) machinery's no-op behavior.
//! The exact solution is `hydro1d`'s test-only exact Riemann solver (the same oracle the 1D Sod
//! test uses), so the two kernels are checked against one analytic answer. The *order* of accuracy
//! is tested separately on a smooth problem (D1).

use euler2d::kernel::{Bc, Grid2D};
use euler2d::state::Prim;
use hydro1d::Primitive;
use hydro1d::riemann::solve;

const GAMMA: f64 = 1.4;
const NZ: usize = 400;
const NR: usize = 4;
const T_END: f64 = 0.2;

/// Mean absolute density error of a fixed `r`-row against the exact self-similar Sod solution
/// (same metric as the 1D `hydro1d` Sod test, for a like-for-like comparison).
fn sod_mean_abs_density_error() -> f64 {
    let dz = 1.0 / NZ as f64;
    let dr = 1.0 / NR as f64;
    let mut g = Grid2D::new(NZ, NR, dz, dr, GAMMA);
    // Uniform in r with reflecting r-walls (the r-sweep does nothing); transmissive z-ends (waves
    // do not reach them by t = 0.2).
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect;
    g.init(|iz, _ir| {
        let z = (iz as f64 + 0.5) * dz;
        if z < 0.5 {
            Prim::new(1.0, 0.0, 0.0, 1.0)
        } else {
            Prim::new(0.125, 0.0, 0.0, 0.1)
        }
    });
    g.run_to(T_END);

    let exact = solve(
        Primitive::new(1.0, 0.0, 1.0),
        Primitive::new(0.125, 0.0, 0.1),
        GAMMA,
    );

    let ir = NR / 2;
    let mut sum = 0.0;
    for iz in 0..NZ {
        let z = (iz as f64 + 0.5) * dz;
        let xi = (z - 0.5) / T_END;
        sum += (g.prim(iz, ir).rho - exact.sample(xi).rho).abs();
    }
    sum / NZ as f64
}

#[test]
fn sod_density_matches_exact_in_l1() {
    let err = sod_mean_abs_density_error();
    // Measured 0.00177 at 400 cells with the second-order MUSCL-Hancock scheme (D1) — close to the
    // 1D AV kernel's 0.00129, and a 4× improvement over the first-order Godunov's 0.00738. The
    // bound guards against scheme regressions (a real bug gives an O(1) miss).
    assert!(err < 2.5e-3, "Sod mean-abs density error = {err}");
}

/// Every `r`-row is identical: the embedded-1D problem stays uniform in `r` (no spurious transverse
/// transport from the r-sweep).
#[test]
fn rows_stay_uniform_in_r() {
    let dz = 1.0 / NZ as f64;
    let dr = 1.0 / NR as f64;
    let mut g = Grid2D::new(NZ, NR, dz, dr, GAMMA);
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect;
    g.init(|iz, _ir| {
        let z = (iz as f64 + 0.5) * dz;
        if z < 0.5 {
            Prim::new(1.0, 0.0, 0.0, 1.0)
        } else {
            Prim::new(0.125, 0.0, 0.0, 0.1)
        }
    });
    g.run_to(T_END);
    for iz in 0..NZ {
        let ref_rho = g.prim(iz, 0).rho;
        for ir in 1..NR {
            assert!(
                (g.prim(iz, ir).rho - ref_rho).abs() < 1e-12,
                "row not uniform at iz={iz}, ir={ir}"
            );
        }
    }
}
