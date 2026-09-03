//! Acceptance tests for the fireball velocity moments (Rung 4, N1 and N7).
//!
//! The exit criterion is an analytic one and it is known before the diagnostic exists: **a
//! self-similar isotropic expansion must return `alpha = 1/3`**. That is the value
//! `eq:reflection_baseline` assumes, so a diagnostic that cannot reproduce it cannot be trusted to
//! report a departure from it — which is the whole of N1.
//!
//! The oblate and prolate cases bracket it on both sides, because a diagnostic that returned 1/3
//! for everything would also pass the isotropic test.

use euler2d::kernel::Grid2D;
use euler2d::moments::{cosine_histogram, moments};
use euler2d::state::Prim;

const GAMMA: f64 = 5.0 / 3.0;

/// A Hubble-like expansion `u = (k_z z, k_r r)` inside an ellipsoid, on an axisymmetric mesh.
///
/// `k_z == k_r` is a self-similar isotropic expansion; unequal rates make it oblate or prolate
/// without changing anything else.
fn hubble_grid(nz: usize, nr: usize, k_z: f64, k_r: f64, radius: f64) -> Grid2D {
    let dz = 2.0 * radius / (nz as f64);
    let dr = radius / (nr as f64);
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.set_axisymmetric(true);
    let z0 = radius; // centre the sphere at z = radius so the mesh spans [0, 2R]
    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz - z0;
        let r = (ir as f64 + 0.5) * dr;
        if (z * z + r * r).sqrt() <= radius {
            // Cold: all the energy is in the resolved velocity field, which is what alpha measures.
            Prim::new(1.0, k_z * z, k_r * r, 1e-12)
        } else {
            Prim::new(1e-9, 0.0, 0.0, 1e-12)
        }
    });
    g
}

#[test]
fn isotropic_expansion_returns_one_third() {
    // THE acceptance test. For a uniform sphere <z²> = <x²> = <y²> = R²/5, so the transverse pair
    // carries twice the axial share and alpha = 1/3 exactly, independent of the expansion rate.
    let g = hubble_grid(120, 60, 1.0, 1.0, 1.0);
    let m = moments(&g);
    assert!(
        (m.alpha - 1.0 / 3.0).abs() < 0.01,
        "isotropic expansion gave alpha = {}, expected 1/3",
        m.alpha
    );
}

#[test]
fn isotropic_result_is_independent_of_expansion_rate() {
    // alpha is a shape, not a speed: doubling the rate must not move it.
    let slow = moments(&hubble_grid(120, 60, 1.0, 1.0, 1.0)).alpha;
    let fast = moments(&hubble_grid(120, 60, 7.5, 7.5, 1.0)).alpha;
    assert!((slow - fast).abs() < 1e-9, "{slow} vs {fast}");
}

#[test]
fn a_pancake_reports_below_one_third_and_a_cigar_above() {
    // The sensitivity N1 is about. Squashing the axial expansion must drive alpha down, which is
    // the "squirts radially into the field" case; stretching it must drive alpha up.
    let pancake = moments(&hubble_grid(120, 60, 0.3, 1.0, 1.0)).alpha;
    let cigar = moments(&hubble_grid(120, 60, 3.0, 1.0, 1.0)).alpha;
    assert!(pancake < 0.33, "pancake gave alpha = {pancake}");
    assert!(cigar > 0.34, "cigar gave alpha = {cigar}");
    assert!(pancake < cigar);
}

#[test]
fn refining_the_mesh_converges_on_one_third() {
    // The isotropic answer is exact, so the only error is the staircased sphere boundary; it must
    // shrink under refinement rather than sit at whatever the coarse mesh happened to give.
    let coarse = (moments(&hubble_grid(60, 30, 1.0, 1.0, 1.0)).alpha - 1.0 / 3.0).abs();
    let fine = (moments(&hubble_grid(160, 80, 1.0, 1.0, 1.0)).alpha - 1.0 / 3.0).abs();
    assert!(fine < coarse, "coarse {coarse}, fine {fine}");
}

#[test]
fn a_uniformly_translating_body_is_all_drift_and_no_thermal_spread() {
    // f_d = 1 and the drift-subtracted remainder is empty. This is the `sqrt(f_d)` limit of
    // eq:reflection_baseline, where a mirror returns the drift and nothing else.
    let mut g = Grid2D::new(40, 20, 0.05, 0.05, GAMMA);
    g.set_axisymmetric(true);
    g.init(|_, _| Prim::new(1.0, 800.0, 0.0, 1e-12));
    let m = moments(&g);
    assert!((m.v_cm_z - 800.0).abs() < 1e-9);
    assert!((m.drift_fraction_kinetic - 1.0).abs() < 1e-9);
    assert!(m.kinetic_thermal / m.kinetic_total < 1e-12);
}

#[test]
fn drift_is_removed_before_alpha_is_taken() {
    // The ask is ambiguous between the drift-subtracted and drift-included readings, and they
    // differ by exactly the f_d N7 is about. Adding a bulk boost must leave `alpha` alone and move
    // `alpha_lab`, or the two are not the quantities they claim to be.
    let still = moments(&hubble_grid(120, 60, 1.0, 1.0, 1.0));
    let dz = 2.0 / 120.0;
    let dr = 1.0 / 60.0;
    let mut boosted = Grid2D::new(120, 60, dz, dr, GAMMA);
    boosted.set_axisymmetric(true);
    boosted.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz - 1.0;
        let r = (ir as f64 + 0.5) * dr;
        if (z * z + r * r).sqrt() <= 1.0 {
            Prim::new(1.0, z + 5.0, r, 1e-12)
        } else {
            Prim::new(1e-9, 0.0, 0.0, 1e-12)
        }
    });
    let m = moments(&boosted);
    assert!(
        (m.alpha - still.alpha).abs() < 1e-6,
        "drift changed alpha: {} vs {}",
        m.alpha,
        still.alpha
    );
    assert!(
        m.alpha_lab > m.alpha + 0.2,
        "alpha_lab should absorb the boost: {} vs {}",
        m.alpha_lab,
        m.alpha
    );
}

#[test]
fn internal_energy_is_reported_so_an_unfinished_expansion_cannot_be_mistaken_for_an_answer() {
    // alpha only means what N1 wants once the pressure has done its work. A hot, still gas is the
    // degenerate case: no motion at all, so any alpha reported would be meaningless — and
    // `unconverted_fraction` must say so.
    let mut g = Grid2D::new(30, 15, 0.05, 0.05, GAMMA);
    g.set_axisymmetric(true);
    g.init(|_, _| Prim::new(1.0, 0.0, 0.0, 1.0e6));
    let m = moments(&g);
    assert!(m.internal > 0.0);
    assert!(m.unconverted_fraction() > 0.999);
}

#[test]
fn an_isotropic_distribution_has_a_flat_cosine_histogram() {
    // For an isotropic velocity field cos(theta) is uniform on [−1, 1], so every bin holds the
    // same mass. This is the check N1 asks for when a single second moment is not a fair summary.
    let g = hubble_grid(160, 80, 1.0, 1.0, 1.0);
    let hist = cosine_histogram(&g, 8);
    let expected = 1.0 / 8.0;
    for (i, h) in hist.iter().enumerate() {
        assert!(
            (h - expected).abs() < 0.02,
            "bin {i} held {h}, expected about {expected}"
        );
    }
    let total: f64 = hist.iter().sum();
    assert!((total - 1.0).abs() < 1e-9);
}

#[test]
fn a_pancake_piles_mass_into_the_transverse_bins() {
    // The histogram must show the shape the moment only summarises: a radially-favoured expansion
    // puts its mass near cos(theta) = 0 rather than spreading evenly.
    let hist = cosine_histogram(&hubble_grid(160, 80, 0.2, 1.0, 1.0), 8);
    let middle = hist[3] + hist[4];
    let ends = hist[0] + hist[7];
    assert!(middle > 2.0 * ends, "middle {middle}, ends {ends}");
}
