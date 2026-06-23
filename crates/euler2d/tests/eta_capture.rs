//! Flat-plate `eta_capture` acceptance (ADR-0003, D2): the lossless 2D/1D wall-impulse ratio.
//!
//! `eta_capture = (J_wall/p_in)_free / (J_wall/p_in)_confined` (see [`euler2d::bounce`]). The
//! confined (plane-wave) run is the perfectly-collimated 1D limit and is cross-checked against the
//! independent 1D `hydro1d` kernel; the free finite-cloud run gives `eta_capture < 1` (radial
//! rebound + edge escape), rising toward 1 as the cloud widens (less radial relief). The committed
//! tests use coarse grids for speed; the headline number is the ignored `diag_eta_capture_table`.

use euler2d::bounce::{Bounce2D, SlugConfig, eta_capture, run_slug_bounce};
use hydro1d::kernel::Tube;

const GAMMA: f64 = 1.4;
const MACH: f64 = 5.0;

/// The confined (plane-wave / 1D-limit) bounce: cloud fills the radius, reflecting outer wall. Few
/// radial cells (uniform in r), so it is cheap.
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
    })
}

/// A free finite-cloud bounce on a large plate (catches the spread, so `eta_capture` measures the
/// rebound-axiality floor).
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
