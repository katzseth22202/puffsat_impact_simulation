//! The trailing pearl (ADR-0055): a second, lighter PuffSat of the same architecture flying in
//! formation a set distance behind the lead, so its gas meets the rebound off the spray plate.
//!
//! The seams, agreed before these tests were written:
//!
//!  1. a pearl of zero mass runs exactly as no pearl (bit-identical);
//!  2. bookkeeping: the pearl holds its share of the lead's mass and moves at the lead's speed, so
//!     the incident momentum grows by that share;
//!  3. causality: the plate force is the no-pearl force until the pearl's gas can arrive;
//!  4. (the existing suite) every other path is unchanged.

use euler2d::bounce::{
    PlateShape, SlugConfig, SprayCloud, TrailingPearl, init_spray_grid, init_spray_grid_with_pearl,
    run_spray_bounce, run_spray_bounce_with_pearl,
};

const GAMMA: f64 = 1.4;
const MACH: f64 = 5.0;
const SPRAY: SprayCloud = SprayCloud {
    depth: 0.5,
    mass_ratio: 10.0,
    standoff: 0.0,
};

/// A confined column (plane-wave) with room above the lead for the pearl.
fn column(nz: usize) -> SlugConfig {
    SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot: 2.0,
        length: 1.0,
        r_plate: 2.0,
        r_max: 2.0,
        z_max: 8.0,
        nr: 6,
        nz,
        confined: true,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    }
}

/// A free dished case, so the zero-mass identity also covers the immersed path.
fn free_dish() -> SlugConfig {
    SlugConfig {
        r_foot: 1.0,
        r_plate: 2.0,
        r_max: 3.0,
        nr: 24,
        nz: 80,
        confined: false,
        shape: PlateShape::Dish { d_over_d: 0.1 },
        ..column(80)
    }
}

#[test]
fn a_pearl_of_zero_mass_runs_exactly_as_no_pearl() {
    let pearl = TrailingPearl {
        mass_ratio: 0.0,
        gap: 2.0,
    };
    for cfg in [column(80), free_dish()] {
        let a = run_spray_bounce(&cfg, SPRAY);
        let b = run_spray_bounce_with_pearl(&cfg, SPRAY, pearl);
        assert_eq!(a.wall_impulse.to_bits(), b.wall_impulse.to_bits());
        assert_eq!(a.incident_momentum.to_bits(), b.incident_momentum.to_bits());
        assert_eq!(a.steps, b.steps);
    }
}

#[test]
fn the_pearl_holds_its_share_of_the_lead_mass_at_the_lead_speed() {
    let share = 0.2;
    let g = init_spray_grid_with_pearl(
        &column(160),
        SPRAY,
        TrailingPearl {
            mass_ratio: share,
            gap: 2.0,
        },
    );
    // Lead top: depth 0.5 + length 1.0 = 1.5; pearl from 1.5 + gap 2.0 = 3.5 to 4.5.
    let (mut lead, mut pearl, mut pearl_momentum) = (0.0_f64, 0.0_f64, 0.0_f64);
    for iz in 0..g.nz() {
        for ir in 0..g.nr() {
            let p = g.prim(iz, ir);
            if p.uz == 0.0 {
                continue; // spray and ambient are at rest
            }
            let m = p.rho * g.cell_volume(ir);
            if (iz as f64 + 0.5) * (8.0 / 160.0) > 3.0 {
                pearl += m;
                pearl_momentum += m * p.uz;
            } else {
                lead += m;
            }
        }
    }
    assert!(
        (pearl / lead - share).abs() < 1e-9,
        "pearl/lead = {}",
        pearl / lead
    );
    assert!(
        (pearl_momentum + pearl).abs() < 1e-9 * pearl,
        "the pearl moves at v = 1"
    );
    let alone = init_spray_grid(&column(160), SPRAY).axial_momentum();
    let with = g.axial_momentum();
    assert!(
        (with / alone - (1.0 + share)).abs() < 1e-9,
        "incident momentum grew by {} not {share}",
        with / alone - 1.0
    );
}

/// The lead's shock crosses the dense cloud at ~0.2 v and loads the plate near t ≈ 2.3, so the pearl
/// trails by 6 (its front at `z = 7.5`). No gas moves toward the plate faster than the cold slug's
/// speed plus its free-expansion front, `1 + 2 c0 / (γ - 1)`. Until half that crossing time the
/// plate force is the no-pearl force, though the lead has by then loaded the plate.
#[test]
fn the_plate_force_is_unchanged_until_the_pearl_can_arrive() {
    let pearl = TrailingPearl {
        mass_ratio: 0.2,
        gap: 6.0,
    };
    let cfg = SlugConfig {
        z_max: 10.0,
        ..column(200)
    };
    let mut a = init_spray_grid(&cfg, SPRAY);
    let mut b = init_spray_grid_with_pearl(&cfg, SPRAY, pearl);
    let start = a.plate_force();
    let speed = 1.0 + 2.0 * (1.0 / MACH) / (GAMMA - 1.0);
    let t = 0.5 * 7.5 / speed;
    // Step both with the same dt: the pearl's bow shock in the ambient sets a smaller stable step,
    // and the comparison is of causality, not of two time-step sequences.
    let mut now = 0.0;
    while now < t {
        let dt = a.stable_dt().min(b.stable_dt()).min(t - now);
        a.step(dt);
        b.step(dt);
        now += dt;
    }
    let (fa, fb) = (a.plate_force(), b.plate_force());
    assert!(
        fa > 2.0 * start,
        "the lead has not loaded the plate yet: {start:e} -> {fa:e}"
    );
    assert!(
        (fa - fb).abs() < 1e-9 * fa,
        "the pearl moved the plate force before it could arrive: {fa:e} vs {fb:e}"
    );
}
