//! Cylindrical Noh implosion acceptance (ADR-0023): a strong converging shock with a closed-form
//! exact solution that directly exercises the axisymmetric radial update.
//!
//! Setup (Noh 1987, cylindrical ν = 2, γ = 5/3): cold gas (`p ≈ 0`) of unit density flowing inward
//! at `u_r = −1`, uniform in `z`. A shock forms at the axis and moves out at speed `(γ−1)/2 = 1/3`.
//! Exact post-shock state (`r < t/3`): `u = 0`, `ρ = ((γ+1)/(γ−1))² = 16`. Pre-shock the convergence
//! amplifies density as `ρ = 1 + t/r`. We check the post-shock density in a band that avoids both
//! the axis (the classic wall-heating anomaly) and the shock front, plus the shock location.

use euler2d::kernel::{Bc, Grid2D};
use euler2d::state::Prim;

const GAMMA: f64 = 5.0 / 3.0;

#[test]
fn cylindrical_noh_post_shock_density_and_shock_radius() {
    let nr = 200;
    let nz = 4;
    let dr = 1.0 / nr as f64;
    let dz = dr;
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect; // axis
    g.bc_rhi = Bc::Transmissive; // far field keeps feeding the inflow (shock never reaches it)
    g.bc_zlo = Bc::Reflect; // uniform in z (no-op)
    g.bc_zhi = Bc::Reflect;
    g.init(|_iz, _ir| Prim::new(1.0, 0.0, -1.0, 1e-6));

    let t_end = 0.6; // shock at r = t/3 = 0.2
    g.run_to(t_end);

    let iz = nz / 2;
    let r_c = |ir: usize| (ir as f64 + 0.5) * dr;

    // Post-shock density in a band behind the shock (r ∈ [0.05, 0.15]), away from the axis.
    let band: Vec<f64> = (0..nr)
        .filter(|&ir| (0.05..0.15).contains(&r_c(ir)))
        .map(|ir| g.prim(iz, ir).rho)
        .collect();
    let mean = band.iter().sum::<f64>() / band.len() as f64;
    assert!(
        (mean - 16.0).abs() / 16.0 < 0.08,
        "cylindrical Noh post-shock density {mean:.2} (exact 16)"
    );

    // Shock radius: the outermost cell whose density exceeds the midpoint of pre/post (≈ 10) marks
    // the front; it should sit near r = t/3 = 0.2.
    let r_shock = (0..nr)
        .filter(|&ir| g.prim(iz, ir).rho > 10.0)
        .map(r_c)
        .fold(0.0_f64, f64::max);
    assert!(
        (r_shock - 0.2).abs() < 0.02,
        "cylindrical Noh shock radius {r_shock:.3} (exact 0.2)"
    );
}
