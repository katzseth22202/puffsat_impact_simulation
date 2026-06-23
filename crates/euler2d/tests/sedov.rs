//! Cylindrical Sedov–Taylor blast acceptance (ADR-0023): a self-similar diverging blast that tests
//! the axisymmetric radial update on an expanding strong shock.
//!
//! Energy `E_L` (per unit length) is deposited in a few cells near the axis (uniform in `z`), the
//! ambient cold. The cylindrical (ν = 2) self-similar shock radius is `R_s(t) = ξ₀ (E_L t²/ρ₀)^{1/4}`
//! with `ξ₀ ≈ 1.0` for γ = 1.4. We check (a) the self-similar power law `R_s ∝ t^{1/2}` between two
//! times (the tight, constant-free check) and (b) the absolute radius against the similarity law
//! (a looser bound, since the blast has finite initial size and the captured shock smears inward).

use euler2d::kernel::{Bc, Grid2D};
use euler2d::state::Prim;
use std::f64::consts::TAU;

const GAMMA: f64 = 1.4;
const E_L: f64 = 1.0; // deposited energy per unit length
const R_BLAST: f64 = 0.03; // initial energy-deposition radius (≪ shock radius at the test times)

/// Radius of the shock front (the density maximum along `r`) on the mid-`z` row.
fn shock_radius(g: &Grid2D, dr: f64) -> f64 {
    let iz = g.nz() / 2;
    let mut best = (0.0, 0.0_f64); // (radius, density)
    for ir in 0..g.nr() {
        let rho = g.prim(iz, ir).rho;
        if rho > best.1 {
            best = ((ir as f64 + 0.5) * dr, rho);
        }
    }
    best.0
}

fn sedov_grid(nr: usize) -> (Grid2D, f64) {
    let nz = 4;
    let dr = 1.0 / nr as f64;
    let mut g = Grid2D::new(nz, nr, dr, dr, GAMMA);
    g.set_axisymmetric(true);
    g.bc_rlo = Bc::Reflect; // axis
    g.bc_rhi = Bc::Transmissive;
    g.bc_zlo = Bc::Reflect; // uniform in z
    g.bc_zhi = Bc::Reflect;

    // Distribute E_L over the near-axis cells: internal energy per unit length in annular cell i is
    // 2π r_c dr · p/(γ−1), so set p_in to hit the target sum.
    let r_c = |ir: usize| (ir as f64 + 0.5) * dr;
    let vol_sum: f64 = (0..nr)
        .map(r_c)
        .filter(|&rc| rc < R_BLAST)
        .map(|rc| TAU * rc * dr)
        .sum();
    let p_in = E_L * (GAMMA - 1.0) / vol_sum;
    g.init(|_iz, ir| {
        let p = if r_c(ir) < R_BLAST { p_in } else { 1e-5 };
        Prim::new(1.0, 0.0, 0.0, p)
    });
    (g, dr)
}

/// Analytic cylindrical self-similar shock radius (ξ₀ ≈ 1.0 for γ = 1.4).
fn sedov_radius(t: f64) -> f64 {
    (E_L * t * t / 1.0).powf(0.25)
}

#[test]
fn cylindrical_sedov_self_similar_radius() {
    let (mut g, dr) = sedov_grid(200);
    let t1 = 0.16;
    let t2 = 0.36;
    // `run_to` advances the state by its argument of simulated time, so step by the increment.
    g.run_to(t1);
    let r1 = shock_radius(&g, dr);
    g.run_to(t2 - t1);
    let r2 = shock_radius(&g, dr);

    // (a) Self-similar power law R_s ∝ t^{1/2}: the ratio is constant-free and should be tight.
    let observed = r2 / r1;
    let expected = (t2 / t1).powf(0.5);
    assert!(
        (observed - expected).abs() / expected < 0.05,
        "Sedov R_s ratio {observed:.3} vs t^1/2 law {expected:.3} (r1={r1:.3}, r2={r2:.3})"
    );

    // (b) Absolute radius against the similarity law (looser: finite blast size + shock smearing).
    for (t, r) in [(t1, r1), (t2, r2)] {
        let analytic = sedov_radius(t);
        assert!(
            (r - analytic).abs() / analytic < 0.15,
            "Sedov R_s({t}) = {r:.3}, similarity law {analytic:.3}"
        );
    }
}
