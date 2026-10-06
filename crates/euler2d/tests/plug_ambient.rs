//! The ambient gas must stand in for vacuum (ADR-0055, plug study). The default ambient, 10⁻³ of
//! the slug's density, is a negligible share of a spray cup's mass. Against one 10 cm necklace
//! sphere on its 3.2 m plate share it weighs ~5 times the plug and ~50 times the sphere, so it
//! tamps the fireball and keeps a standing load on the plate above the 10⁻³-of-peak quiet cutoff;
//! plug case 0 integrated that load to the step cap and reported 12 times the overtake ceiling.
//!
//! The seams, agreed before these tests were written:
//!
//!  1. the default ambient runs exactly as before (bit-identical);
//!  2. a narrow plug in thin ambient stays under the overtake ceiling and ends on a quiet tail.

use euler2d::bounce::{
    AMBIENT_DEFAULT, PlateShape, SlugConfig, SprayCloud, run_spray_bounce,
    run_spray_bounce_with_ambient,
};

const K: f64 = 10.0;

/// One sphere (5 cells in radius) over a plug of `K` times its mass, standing off a wide dish.
fn narrow_plug() -> (SlugConfig, SprayCloud) {
    let dr = 0.02;
    let nz = 83;
    let cfg = SlugConfig {
        gamma: 1.4,
        mach: 20.0,
        r_foot: 0.1,
        length: 0.15,
        r_plate: 1.2,
        r_max: 1.5,
        z_max: nz as f64 * dr,
        nr: 75,
        nz,
        confined: false,
        shape: PlateShape::Dish { d_over_d: 0.10 },
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let spray = SprayCloud {
        depth: 0.6,
        mass_ratio: K,
        standoff: 0.3,
    };
    (cfg, spray)
}

#[test]
fn the_default_ambient_runs_exactly_as_before() {
    let (mut cfg, spray) = narrow_plug();
    // A short, coarse copy: bit-identity needs no resolution.
    cfg.nr = 30;
    cfg.nz = 34;
    cfg.z_max = 34.0 * 0.05;
    let a = run_spray_bounce(&cfg, spray);
    let b = run_spray_bounce_with_ambient(&cfg, spray, AMBIENT_DEFAULT);
    assert_eq!(a.wall_impulse.to_bits(), b.wall_impulse.to_bits());
    assert_eq!(a.steps, b.steps);
}

#[test]
fn a_narrow_plug_in_thin_ambient_stays_under_the_ceiling_and_ends_quietly() {
    let (cfg, spray) = narrow_plug();
    let r = run_spray_bounce_with_ambient(&cfg, spray, 1.0e-6);
    let beta = r.restitution_ratio();
    let ceiling = 1.0 + (1.0 + K).sqrt();
    assert!(
        beta <= ceiling,
        "beta {beta} above the overtake ceiling {ceiling}"
    );
    assert!(beta > 0.0, "beta {beta}: the plate must be pushed");
    let cap = 400 * cfg.nz + 50_000;
    assert!(
        r.steps < cap,
        "ran to the step cap ({cap}): the tail never went quiet"
    );
}
