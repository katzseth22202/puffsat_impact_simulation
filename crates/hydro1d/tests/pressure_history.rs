//! The wall-pressure history of a bounce (ADR-0055 step 2a, spike width for the face allowable).
//!
//! Seam, agreed before this test was written: an ideal-gas slab (ρ₀ = 1, u = 1 toward the wall,
//! free rear face) holds the wall at the reflected-shock pressure `p₁` for exactly
//! `T = t_m (1 + D/c₁)`. `D` is the reflected shock's speed, `t_m = 1/(D + 1 + c₀)` is when the rear
//! face's expansion head meets it, and `c₁` is the sound speed of the stagnated gas the head then
//! crosses. Both `p₁` and `T` are closed forms, independent of the kernel.

use hydro1d::kernel::Tube;

const GAMMA: f64 = 1.4;
const MACH: f64 = 5.0;

/// `(p₁, T)` for the slug: reflected-shock pressure and plateau duration.
fn plateau() -> (f64, f64) {
    let u = 1.0;
    let c0 = 1.0 / MACH;
    let p0 = 1.0 / (GAMMA * MACH * MACH);
    let d = -(3.0 - GAMMA) / 4.0 * u + (((GAMMA + 1.0) / 4.0 * u).powi(2) + c0 * c0).sqrt();
    let p1 = p0 + u * (d + u);
    let rho1 = (d + u) / d;
    let c1 = (GAMMA * p1 / rho1).sqrt();
    let t_m = 1.0 / (d + u + c0);
    (p1, t_m * (1.0 + d / c1))
}

#[test]
fn wall_holds_the_reflected_shock_pressure_for_the_closed_form_time() {
    let (p1, t_plateau) = plateau();
    let mut tube = Tube::slug(800, MACH, GAMMA);
    let (_, history) = tube.run_bounce_with_history();

    // The level: the median wall pressure over the middle of the predicted plateau.
    let mut mid: Vec<f64> = history
        .iter()
        .filter(|(t, _)| *t > 0.3 * t_plateau && *t < 0.7 * t_plateau)
        .map(|&(_, p)| p)
        .collect();
    mid.sort_by(f64::total_cmp);
    let level = mid[mid.len() / 2];
    assert!(
        (level - p1).abs() < 0.02 * p1,
        "plateau {level:.5} vs reflected-shock p1 {p1:.5}"
    );

    // The duration: from first reaching to last holding 95% of p1.
    let above: Vec<f64> = history
        .iter()
        .filter(|(_, p)| *p >= 0.95 * p1)
        .map(|&(t, _)| t)
        .collect();
    let width = above[above.len() - 1] - above[0];
    assert!(
        (width - t_plateau).abs() < 0.05 * t_plateau,
        "plateau lasts {width:.4} vs closed form {t_plateau:.4}"
    );
}
