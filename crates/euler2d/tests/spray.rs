//! The spray plate in 2-D (ADR-0055 step 2b): a resting spray cloud on the plate with the PuffSat's
//! slug arriving above it, as one effective-γ gas (the ADR-0003/0008 geometry method).
//!
//! The seams, agreed before these tests were written:
//!
//!  1. bookkeeping: the spray holds `k` times the slug's mass and is at rest, so the incident
//!     momentum is the slug's alone;
//!  2. cross-code: the confined (plane-wave) run matches the independent 1-D Lagrangian kernel's
//!     layered bounce at the same γ, densities and depths;
//!  3. geometry: free never beats confined, and a wide, short column captures more than a tall,
//!     narrow one;
//!  4. (the existing suite) the slug-only paths are unchanged.

use euler2d::bounce::{
    Bounce2D, PlateShape, SlugConfig, SprayCloud, eta_capture, init_spray_grid, run_spray_bounce,
};
use hydro1d::eos::IdealGas;
use hydro1d::kernel::{Boundary, LagrangianState, Tube, Viscosity};

const GAMMA: f64 = 1.4;
const MACH: f64 = 5.0;
const K: f64 = 10.0;
const DEPTH: f64 = 2.0;
const SPRAY: SprayCloud = SprayCloud {
    depth: DEPTH,
    mass_ratio: K,
};

fn cfg(r_foot: f64, r_plate: f64, r_max: f64, nr: usize, nz: usize, confined: bool) -> SlugConfig {
    SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot,
        length: 1.0,
        r_plate,
        r_max,
        z_max: 8.0,
        nr,
        nz,
        confined,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    }
}

fn confined(nz: usize) -> SlugConfig {
    cfg(5.0, 5.0, 5.0, 8, nz, true)
}

#[test]
fn spray_holds_k_times_the_slug_mass_at_rest() {
    let g = init_spray_grid(&confined(120), SPRAY);
    let (mut spray, mut slug, mut slug_momentum) = (0.0_f64, 0.0_f64, 0.0_f64);
    for iz in 0..g.nz() {
        for ir in 0..g.nr() {
            let p = g.prim(iz, ir);
            if p.rho < 1.0e-2 {
                continue; // ambient
            }
            let m = p.rho * g.cell_volume(ir);
            if p.uz == 0.0 {
                spray += m;
            } else {
                slug += m;
                slug_momentum += m * p.uz;
            }
        }
    }
    assert!(
        (spray / slug - K).abs() < 1e-9 * K,
        "spray/slug = {}",
        spray / slug
    );
    // `axial_momentum` drops the common 2π (see its docs); `cell_volume` includes it.
    let p_in = 2.0 * std::f64::consts::PI * g.axial_momentum();
    assert!((p_in - slug_momentum).abs() < 1e-9 * slug_momentum.abs());
    assert!(
        (slug_momentum + slug).abs() < 1e-9 * slug,
        "slug moves at v = 1"
    );
}

/// The same layered column in the 1-D Lagrangian kernel: spray (at rest, density `k L / depth`)
/// against the wall, the slug (density 1, speed 1) outside it, one cold pressure throughout.
fn layered_1d_ratio(cells_per_m: usize) -> f64 {
    let p0 = 1.0 / (GAMMA * MACH * MACH);
    let rho_s = K * 1.0 / DEPTH;
    let zones = [(DEPTH, rho_s, 0.0), (1.0, 1.0, -1.0)];
    let mut positions = vec![0.0];
    let (mut mass, mut vel, mut energy) = (vec![], vec![], vec![]);
    for (len, rho, u) in zones {
        let n = (len * cells_per_m as f64).round() as usize;
        let dx = len / n as f64;
        for _ in 0..n {
            positions.push(positions[positions.len() - 1] + dx);
            mass.push(rho * dx);
            vel.push(u);
            energy.push(p0 / ((GAMMA - 1.0) * rho));
        }
    }
    let n = mass.len();
    let mut node_velocities = vec![0.0; n + 1];
    for i in 1..n {
        node_velocities[i] =
            (mass[i - 1] * vel[i - 1] + mass[i] * vel[i]) / (mass[i - 1] + mass[i]);
    }
    node_velocities[n] = vel[n - 1];
    let state = LagrangianState {
        positions,
        node_velocities,
        cell_masses: mass,
        specific_energies: energy,
    };
    let mut tube = Tube::from_lagrangian(
        state,
        IdealGas::new(GAMMA),
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let r = tube.run_bounce();
    r.wall_impulse / r.incident_momentum
}

#[test]
fn confined_spray_bounce_matches_the_1d_layered_bounce() {
    let ratio_2d = run_spray_bounce(&confined(120), SPRAY).restitution_ratio();
    let ratio_1d = layered_1d_ratio(200);
    let rel = (ratio_2d - ratio_1d).abs() / ratio_1d;
    assert!(
        rel < 0.08,
        "confined 2D J/p_in = {ratio_2d:.4}, 1D = {ratio_1d:.4} (rel {rel:.3})"
    );
}

#[test]
fn free_never_beats_confined_and_a_wide_column_captures_more() {
    let denom: Bounce2D = run_spray_bounce(&confined(80), SPRAY);
    let tall = run_spray_bounce(&cfg(1.0, 2.0, 3.0, 24, 80, false), SPRAY);
    let wide = run_spray_bounce(&cfg(4.0, 4.5, 5.0, 40, 80, false), SPRAY);
    let (eta_tall, eta_wide) = (eta_capture(&tall, &denom), eta_capture(&wide, &denom));
    assert!(
        eta_tall <= 1.01 && eta_wide <= 1.01,
        "{eta_tall}, {eta_wide}"
    );
    assert!(eta_wide > eta_tall, "wide {eta_wide} vs tall {eta_tall}");
}
