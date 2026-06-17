//! Order-of-accuracy test (ADR-0022): the kernel must converge at its formal **rate 2** in
//! smooth flow. The problem is a standing acoustic wave in a rigid-walled tube (`u = 0` at both
//! ends — the kernel's fixed-endpoint BC), evolved one full period. Artificial viscosity is
//! quadratic-only so its `O(Δx²)` term cannot mask the rate (the linear term would be `O(Δx)`
//! and cap convergence at first order).
//!
//! # Why self-convergence rather than a linear-acoustics reference
//!
//! Comparing to the *linear* solution has no usable amplitude window: this scheme's error
//! constant is small, so the discretization error is squeezed between the f64 round-off floor
//! (at small amplitude) and the `O(ε²)` nonlinearity floor (at large amplitude). Instead we use
//! **Richardson self-convergence** — compare the kernel to itself at `N, 2N, 4N, 8N`. All
//! solutions approach the same (nonlinear) continuum solution at rate 2, so the nonlinearity
//! cancels in the differences and the amplitude can be large enough to bury round-off.
//!
//! The grid-independent functional we converge is the fundamental-mode amplitude
//! `a = (2/L) ∫₀ᴸ (ρ − ρ₀) cos(kx) dx` (a cell-centered midpoint quadrature), which avoids any
//! cross-grid mesh-alignment issue. With `a(N) = a_∞ + C·Δx² + …`, successive differences shrink
//! by 4× per refinement, i.e. `log₂(Δaₙ / Δaₙ₊₁) → 2`.

use hydro1d::kernel::{Tube, Viscosity};
use std::f64::consts::PI;

const GAMMA: f64 = 1.4;
const RHO0: f64 = 1.0;
const P0: f64 = 1.0;
const L: f64 = 1.0;
const EPS: f64 = 1e-2; // large amplitude: round-off is buried; nonlinearity cancels in self-diff

/// Fundamental-mode density amplitude after one wave period at resolution `cells`.
fn modal_amplitude(cells: usize) -> f64 {
    let c0 = (GAMMA * P0 / RHO0).sqrt();
    let k = PI / L;
    let period = 2.0 * PI / (c0 * k);

    let dx = L / cells as f64;
    let x: Vec<f64> = (0..=cells).map(|i| i as f64 * dx).collect();
    let vel = vec![0.0; cells];
    let mut rho = vec![0.0; cells];
    let mut pressure = vec![0.0; cells];
    for j in 0..cells {
        let xc = (j as f64 + 0.5) * dx;
        let dp = EPS * P0 * (k * xc).cos();
        pressure[j] = P0 + dp;
        rho[j] = RHO0 + dp / (c0 * c0);
    }

    let mut tube = Tube::new(x, &rho, &vel, &pressure, GAMMA, Viscosity::QUADRATIC_ONLY);
    tube.run_to(period);

    let integral: f64 = (0..tube.cells())
        .map(|j| (tube.density(j) - RHO0) * (k * tube.center(j)).cos() * tube.width(j))
        .sum();
    (2.0 / L) * integral
}

#[test]
fn acoustic_wave_converges_at_second_order() {
    let a: Vec<f64> = [32usize, 64, 128, 256, 512]
        .iter()
        .map(|&n| modal_amplitude(n))
        .collect();
    // Successive differences and the rate at which they shrink (→ 2 for a 2nd-order scheme).
    let diffs: Vec<f64> = a.windows(2).map(|w| (w[0] - w[1]).abs()).collect();
    let rates: Vec<f64> = diffs.windows(2).map(|w| (w[0] / w[1]).log2()).collect();

    // Measured rates [2.018, 2.004, 2.001] tightening to 2 as the grid refines; this band
    // guards the 2nd-order property against a regression to first order (which reads ≈1).
    assert!(
        rates.iter().all(|&r| (1.9..=2.1).contains(&r)),
        "expected second-order convergence; diffs = {diffs:?}, rates = {rates:?}"
    );
}
