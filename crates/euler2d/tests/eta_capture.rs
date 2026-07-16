//! `eta_capture` acceptance (ADR-0003): the lossless 2D/1D wall-impulse ratio, flat (D2) and
//! shallow-concave (D5).
//!
//! `eta_capture = (J_wall/p_in)_free / (J_wall/p_in)_confined` (see [`euler2d::bounce`]). The
//! confined (plane-wave) run is the perfectly-collimated 1D limit and is cross-checked against the
//! independent 1D `hydro1d` kernel; the free finite-cloud run gives `eta_capture < 1` (radial
//! rebound + edge escape), rising toward 1 as the cloud widens (less radial relief). A shallow
//! **concave** plate (ADR-0021) bends the outward rebound back toward the axis, raising the captured
//! axial momentum — `eta_capture` rises with the dish depth. The committed tests use coarse grids
//! for speed; the headline numbers are the ignored `diag_*` tables.

use euler2d::bounce::{Bounce2D, PlateShape, SlugConfig, eta_capture, run_slug_bounce};
use hydro1d::kernel::Tube;

const GAMMA: f64 = 1.4;
const MACH: f64 = 5.0;

/// The confined (plane-wave / 1D-limit) bounce: cloud fills the radius, reflecting outer wall, flat
/// grid-aligned plate. Few radial cells (uniform in r), so it is cheap.
fn confined(nz: usize) -> Bounce2D {
    run_slug_bounce(&SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot: 5.0,
        length: 1.0,
        r_plate: 5.0,
        r_max: 5.0,
        z_max: 3.0,
        nr: 8,
        nz,
        confined: true,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    })
}

/// A free finite-cloud bounce on a large grid-aligned flat plate (the D2 footprint study).
fn free(r_foot: f64, nr: usize, nz: usize) -> Bounce2D {
    run_slug_bounce(&SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot,
        length: 1.0,
        r_plate: 4.0,
        r_max: 5.0,
        z_max: 3.0,
        nr,
        nz,
        confined: false,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    })
}

/// A free finite-cloud bounce on a plate of the given shape (`d/D = 0` is the immersed flat
/// baseline; `> 0` is shallow concave). `r_plate = 2`, `r_max = 3` so the rebound spreads across the
/// tilted dish before escaping the rim.
fn free_shape(shape: PlateShape, r_foot: f64, nr: usize, nz: usize) -> Bounce2D {
    run_slug_bounce(&SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot,
        length: 1.0,
        r_plate: 2.0,
        r_max: 3.0,
        z_max: 3.0,
        nr,
        nz,
        confined: false,
        shape,
        taper_frac: 0.0,
        alpha_div: 0.0,
    })
}

/// The confined plane-wave bounce reproduces the independent 1D `hydro1d` kernel's `1 + e_eff` — an
/// Eulerian-Godunov vs Lagrangian-AV cross-check that the 2D bounce physics is right.
#[test]
fn confined_matches_1d_bounce() {
    let ratio_2d = confined(60).restitution_ratio();
    let e_eff_1d = {
        let mut tube = Tube::slug(400, MACH, GAMMA);
        tube.run_bounce().e_eff
    };
    let ratio_1d = 1.0 + e_eff_1d;
    let rel = (ratio_2d - ratio_1d).abs() / ratio_1d;
    assert!(
        rel < 0.08,
        "confined 2D 1+e_eff = {ratio_2d:.4}, 1D = {ratio_1d:.4} (rel {rel:.3})"
    );
}

/// `eta_capture` is a genuine capture fraction in `(0, 1)`, and a wider cloud re-collimates more
/// (slower radial relief) so it captures more axial momentum — `eta_capture` rises with the
/// footprint, toward the 1D ceiling of 1.
#[test]
fn flat_plate_eta_capture_is_a_fraction_rising_with_footprint() {
    let denom = confined(28);
    let eta_narrow = eta_capture(&free(0.6, 32, 28), &denom);
    let eta_wide = eta_capture(&free(2.0, 32, 28), &denom);
    assert!(
        eta_narrow > 0.0 && eta_wide < 1.0,
        "eta_capture out of (0,1): narrow {eta_narrow:.3}, wide {eta_wide:.3}"
    );
    assert!(
        eta_wide > eta_narrow,
        "expected eta to rise with footprint: narrow {eta_narrow:.3}, wide {eta_wide:.3}"
    );
}

/// The immersed flat plate (`Dish` with `d/D = 0`, a raised flat wall through the IBM) reproduces the
/// verified grid-aligned flat-plate `eta_capture` for the same cloud — the consistency gate that
/// makes the curvature gain (concave vs immersed-flat, both via the IBM) trustworthy, since the
/// IBM's boundary error does not cancel against the grid-aligned confined denominator.
#[test]
fn immersed_flat_matches_grid_aligned_flat() {
    let denom = confined(32);
    let r_foot = 1.0;
    let (nr, nz) = (40, 32);
    let eta_ibm = eta_capture(
        &free_shape(PlateShape::Dish { d_over_d: 0.0 }, r_foot, nr, nz),
        &denom,
    );
    // Same cloud/plate, grid-aligned flat (r_plate = 2, r_max = 3 to match `free_shape`).
    let eta_grid = eta_capture(
        &run_slug_bounce(&SlugConfig {
            gamma: GAMMA,
            mach: MACH,
            r_foot,
            length: 1.0,
            r_plate: 2.0,
            r_max: 3.0,
            z_max: 3.0,
            nr,
            nz,
            confined: false,
            shape: PlateShape::FlatGridAligned,
            taper_frac: 0.0,
            alpha_div: 0.0,
        }),
        &denom,
    );
    let rel = (eta_ibm - eta_grid).abs() / eta_grid;
    assert!(
        rel < 0.10,
        "IBM-flat eta {eta_ibm:.4} vs grid-aligned flat {eta_grid:.4} (rel {rel:.3})"
    );
}

/// A shallow-concave plate re-collimates the rebound toward the axis, so `eta_capture` **rises with
/// the dish depth** `d/D` — the recovery lever (ADR-0021). All three runs go through the same
/// immersed boundary (`d/D = 0, 0.10, 0.15`), so the gain is pure curvature.
#[test]
fn concave_eta_capture_rises_with_depth() {
    let denom = confined(32);
    let r_foot = 1.0;
    let (nr, nz) = (48, 32);
    let eta_flat = eta_capture(
        &free_shape(PlateShape::Dish { d_over_d: 0.0 }, r_foot, nr, nz),
        &denom,
    );
    let eta_mild = eta_capture(
        &free_shape(PlateShape::Dish { d_over_d: 0.10 }, r_foot, nr, nz),
        &denom,
    );
    let eta_deep = eta_capture(
        &free_shape(PlateShape::Dish { d_over_d: 0.15 }, r_foot, nr, nz),
        &denom,
    );
    assert!(
        eta_mild > eta_flat && eta_deep > eta_mild,
        "expected eta to rise with depth: flat {eta_flat:.4}, 0.10 {eta_mild:.4}, 0.15 {eta_deep:.4}"
    );
}

/// A concave plate focuses the rebound, so its peak *local* facesheet pressure exceeds the flat
/// plate's near-uniform stagnation pressure at the same cloud/M — the survivability concentration
/// factor (Rung S, ADR-0010/0021), the quantitative cousin of why the deep dish is foreclosed.
#[test]
fn concave_focuses_local_peak_above_flat() {
    let r_foot = 1.0;
    let (nr, nz) = (48, 32);
    let flat = free_shape(PlateShape::Dish { d_over_d: 0.0 }, r_foot, nr, nz);
    let deep = free_shape(PlateShape::Dish { d_over_d: 0.15 }, r_foot, nr, nz);
    let focus = deep.peak_local_pressure / flat.peak_local_pressure;
    assert!(
        focus > 1.0 && focus < 8.0,
        "concave should focus the local peak above flat: flat {:.3e}, deep {:.3e}, factor {focus:.3}",
        flat.peak_local_pressure,
        deep.peak_local_pressure,
    );
}

/// The flat-plate peak local pressure is positive and the same order as the r-weighted mean implied
/// by the integrated plate force (`peak_force / (r_plate²/2)`) — the local and integrated loads agree.
#[test]
fn flat_local_peak_ties_to_integrated_force() {
    let r_plate = 2.0;
    let flat = free_shape(PlateShape::Dish { d_over_d: 0.0 }, 1.0, 48, 32);
    let mean = flat.peak_force / (0.5 * r_plate * r_plate);
    assert!(flat.peak_local_pressure > 0.0);
    assert!(
        flat.peak_local_pressure > 0.2 * mean && flat.peak_local_pressure < 8.0 * mean,
        "flat local peak {:.3e} should track the r-weighted mean {:.3e}",
        flat.peak_local_pressure,
        mean,
    );
}

/// DIAGNOSTIC (ignored): the flat-plate `eta_capture` vs footprint table at a finer resolution, plus
/// the confined-vs-1D cross-check. Run with `cargo test -p euler2d --test eta_capture -- --ignored
/// --nocapture diag`.
#[test]
#[ignore = "headline table; slow (fine grid)"]
fn diag_eta_capture_table() {
    let denom = confined(60);
    let e_eff_1d = {
        let mut tube = Tube::slug(400, MACH, GAMMA);
        tube.run_bounce().e_eff
    };
    eprintln!(
        "confined 1+e_eff = {:.4}  (1D {:.4})",
        denom.restitution_ratio(),
        1.0 + e_eff_1d
    );
    for r_foot in [0.5, 1.0, 1.5, 2.0] {
        let b = free(r_foot, 100, 60);
        eprintln!(
            "r_foot/L={r_foot:.1}: eta_capture = {:.4}  (peak F={:.3e})",
            eta_capture(&b, &denom),
            b.peak_force,
        );
    }
}

/// DIAGNOSTIC (ignored): the shallow-concave `eta_capture` vs depth table at finer resolution — the
/// headline curvature gain. Run with `cargo test -p euler2d --test eta_capture -- --ignored
/// --nocapture diag_concave`.
#[test]
#[ignore = "headline table; slow (fine grid)"]
fn diag_concave_eta_capture_table() {
    let denom = confined(60);
    let (nr, nz) = (80, 48);
    for r_foot in [1.0, 1.5] {
        eprintln!("--- r_foot = {r_foot:.1} ---");
        for d_over_d in [0.0, 0.10, 0.15] {
            let b = free_shape(PlateShape::Dish { d_over_d }, r_foot, nr, nz);
            eprintln!(
                "  d/D={d_over_d:.2}: eta_capture = {:.4}  (peak F={:.3e}, steps={})",
                eta_capture(&b, &denom),
                b.peak_force,
                b.steps,
            );
        }
    }
}
