//! Sod shock-tube acceptance test (ADR-0022): the staggered AV-Lagrangian kernel must
//! reproduce the exact Riemann solution to tolerance. This is rung A's shock-capturing
//! correctness check — the *order* of accuracy is tested separately on a smooth problem.

use hydro1d::Primitive;
use hydro1d::kernel::Tube;
use hydro1d::riemann::solve;

const GAMMA: f64 = 1.4;
const CELLS: usize = 400;
const T_END: f64 = 0.2;

/// Mean absolute density error of the kernel against the exact self-similar solution, sampled
/// at each (Lagrangian-moved) cell center.
fn sod_l1_density_error() -> f64 {
    let mut tube = Tube::sod(CELLS, GAMMA);
    tube.run_to(T_END);

    let exact = solve(
        Primitive::new(1.0, 0.0, 1.0),
        Primitive::new(0.125, 0.0, 0.1),
        GAMMA,
    );

    let mut sum = 0.0;
    for j in 0..tube.cells() {
        let xi = (tube.center(j) - 0.5) / T_END;
        sum += (tube.density(j) - exact.sample(xi).rho).abs();
    }
    sum / tube.cells() as f64
}

#[test]
fn sod_density_matches_exact_in_l1() {
    // Measured 0.00129 at 400 cells; this bound guards against regressions in the scheme.
    let l1 = sod_l1_density_error();
    assert!(l1 < 0.0015, "Sod L1 density error = {l1}");
}
