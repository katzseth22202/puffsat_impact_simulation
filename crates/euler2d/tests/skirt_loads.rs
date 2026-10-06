//! The skirt's outward load (ADR-0055, Q56): the gas pressure on the cup wall's inner face, row by
//! row up the wall, recorded through the bounce so the hoop can be sized.
//!
//! The seams, agreed before these tests were written:
//!
//!  1. a cup full of still gas at pressure `p` reads exactly `p` on every wall row, and has one row
//!     per cell between the rim and the wall top;
//!  2. a plate with no wall reports no rows;
//!  3. recording is passive: the bounce with its skirt history is bit-identical to the plain one.

use euler2d::bounce::{
    PlateShape, SlugConfig, SprayCloud, init_slug_grid, run_spray_bounce,
    run_spray_bounce_with_skirt_history,
};
use euler2d::state::Prim;

/// A flat-floored cup on 0.125 cells: the floor sits four cells up (z = 0.5) and the 1.0 wall
/// rises to z = 1.5, so the wall rows are the eight cells centred at 0.5625 ... 1.4375. The wall's
/// inner face is at r = 2.0, two cells thick.
fn cup(skirt_height: f64) -> SlugConfig {
    SlugConfig {
        gamma: 1.4,
        mach: 5.0,
        r_foot: 1.0,
        length: 1.0,
        r_plate: 2.0,
        r_max: 3.0,
        z_max: 4.0,
        nr: 24,
        nz: 32,
        confined: false,
        shape: PlateShape::Cup {
            d_over_d: 0.0,
            skirt_height,
            flare: 0.0,
            thickness: 0.25,
        },
        taper_frac: 0.0,
        alpha_div: 0.0,
    }
}

#[test]
fn still_gas_reads_its_pressure_on_every_wall_row() {
    let p = 0.37;
    let mut g = init_slug_grid(&cup(1.0));
    g.init(|_, _| Prim::new(1.0, 0.0, 0.0, p));
    let rows = g.skirt_face_pressure();
    let heights: Vec<f64> = (4..12).map(|iz| (f64::from(iz) + 0.5) * 0.125).collect();
    assert_eq!(rows.len(), heights.len(), "rows {rows:?}");
    for ((z, q), want) in rows.iter().zip(&heights) {
        assert!((z - want).abs() < 1e-12, "row at {z}, want {want}");
        assert_eq!(q.to_bits(), p.to_bits(), "row at {z} reads {q}");
    }
}

#[test]
fn a_plate_without_a_wall_has_no_rows() {
    for shape in [
        PlateShape::Dish { d_over_d: 0.1 },
        PlateShape::Cup {
            d_over_d: 0.1,
            skirt_height: 0.0,
            flare: 0.0,
            thickness: 0.25,
        },
    ] {
        let g = init_slug_grid(&SlugConfig { shape, ..cup(1.0) });
        assert!(g.skirt_face_pressure().is_empty());
    }
}

#[test]
fn recording_the_skirt_leaves_the_bounce_unchanged() {
    let cfg = SlugConfig {
        r_max: 3.0,
        z_max: 8.0,
        nz: 64,
        ..cup(1.0)
    };
    let spray = SprayCloud {
        depth: 1.0,
        mass_ratio: 10.0,
        standoff: 0.25,
    };
    let plain = run_spray_bounce(&cfg, spray);
    let (with, history) = run_spray_bounce_with_skirt_history(&cfg, spray);
    assert_eq!(plain.wall_impulse.to_bits(), with.wall_impulse.to_bits());
    assert_eq!(plain.peak_force.to_bits(), with.peak_force.to_bits());
    assert_eq!(plain.steps, with.steps);
    assert_eq!(
        history.time.len(),
        with.steps + 1,
        "one sample per step plus the start"
    );
    assert!(history.time.windows(2).all(|w| w[1] > w[0]));
    assert!(
        history
            .pressure
            .iter()
            .all(|row| row.len() == history.z.len())
    );
    assert_eq!(history.z.len(), 8);
}

// ---- The bowl's steep band (Q63) ----------------------------------------------------------------
//
// The seams, agreed before these tests were written:
//
//  1. still gas at `p` reads `p` on every riser of the dish floor, and the risers' radial
//     projected area `Σ r dz` matches the closed form `(2/3) r_plate · depth` of the paraboloid
//     `z = z0 + depth (r / r_plate)²`;
//  2. a flat floor has no risers;
//  3. recording is passive.

/// A bare dish on 0.05 cells: floor on the axis at z0 = 0.2 (four cells), rim 1.2 higher.
fn bowl(d_over_d: f64) -> SlugConfig {
    SlugConfig {
        r_max: 4.0,
        z_max: 6.0,
        nr: 80,
        nz: 120,
        shape: PlateShape::Dish { d_over_d },
        ..cup(0.0)
    }
}

/// The risers' projected area against the closed form over the same rows, at `refine` times the
/// base resolution. Each riser sits at the first solid cell's inner face, within half a cell of the
/// surface `r_s(z) = r_plate √((z − z0) / depth)`, so the sum carries no bias; the rows end where
/// the last floor cell inside the rim does.
fn riser_area_error(refine: usize) -> f64 {
    let p = 0.37;
    let cfg = SlugConfig {
        nr: 80 * refine,
        nz: 120 * refine,
        ..bowl(0.3) // r_plate 2.0, depth 1.2, z0 four base cells = 0.2
    };
    let mut g = init_slug_grid(&cfg);
    g.init(|_, _| Prim::new(1.0, 0.0, 0.0, p));
    let risers = g.floor_riser_pressure();
    assert!(risers.len() > 20, "risers {}", risers.len());
    for &(_, _, q) in &risers {
        assert_eq!(q.to_bits(), p.to_bits());
    }
    let dz = 6.0 / (120 * refine) as f64;
    let z0 = 4.0 * dz;
    let projected: f64 = risers.iter().map(|&(_, r, _)| r * dz).sum();
    let top = risers.last().map_or(z0, |&(z, _, _)| z + 0.5 * dz);
    let bottom = risers.first().map_or(z0, |&(z, _, _)| z - 0.5 * dz);
    let integral = |z: f64| 2.0 / 3.0 * 2.0 * (z - z0).max(0.0).powf(1.5) / 1.2_f64.sqrt();
    let exact = integral(top) - integral(bottom);
    (projected - exact).abs() / exact
}

/// Within a percent on the base grid, and closer on a finer one.
#[test]
fn still_gas_reads_its_pressure_on_every_riser_of_the_bowl() {
    let coarse = riser_area_error(1);
    let fine = riser_area_error(2);
    assert!(coarse < 0.01, "coarse relative error {coarse}");
    assert!(fine < coarse, "error {coarse} -> {fine} under refinement");
}

#[test]
fn a_flat_floor_has_no_risers() {
    let mut g = init_slug_grid(&bowl(0.0));
    g.init(|_, _| Prim::new(1.0, 0.0, 0.0, 1.0));
    assert!(g.floor_riser_pressure().is_empty());
}

#[test]
fn recording_the_bowl_leaves_the_bounce_unchanged() {
    let cfg = SlugConfig {
        z_max: 8.0,
        nz: 160,
        ..bowl(0.2)
    };
    let spray = SprayCloud {
        depth: 1.0,
        mass_ratio: 8.52,
        standoff: 0.25,
    };
    let plain = run_spray_bounce(&cfg, spray);
    let (with, history) = run_spray_bounce_with_skirt_history(&cfg, spray);
    assert_eq!(plain.wall_impulse.to_bits(), with.wall_impulse.to_bits());
    assert_eq!(plain.steps, with.steps);
    assert!(history.z.is_empty(), "a bare dish has no skirt rows");
    assert!(!history.band_r.is_empty());
    assert_eq!(history.band_r.len(), history.band_z.len());
    assert!(
        history
            .band_pressure
            .iter()
            .all(|row| row.len() == history.band_r.len())
    );
    assert_eq!(history.band_pressure.len(), history.time.len());
}
