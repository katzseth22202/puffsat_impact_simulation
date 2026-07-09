//! Momentum-limit smoke tests (ADR-0001): a finite cold gas slug bounces off a rigid wall and
//! re-expands into vacuum. On the elastic (reflecting) wall, two checks of different character:
//!
//! - a value-independent **momentum-conservation** check (the wall impulse integrated from the
//!   wall-pressure history matches the gas momentum change), and
//! - the lossless **bounce ceiling** `e_eff < 1` (set by γ) — the *true* upper bound on a finite
//!   slug's restitution.
//!
//! Note the slug's `e_eff` cannot reach the elastic limit `e_eff → 1`: that is the
//! *small-amplitude* (acoustic) limit, demonstrated by the convergence standing wave, not
//! something a finite high-amplitude slug can produce (a low-Mach slug is pressure-dominated, not
//! a gentle projectile). See ADR-0001.

use hydro1d::kernel::{BounceResult, Tube};

const GAMMA: f64 = 5.0 / 3.0; // monatomic; the bounce default (ADR-0001)

fn bounce(cells: usize, mach: f64) -> BounceResult {
    Tube::slug(cells, mach, GAMMA).run_bounce()
}

/// Conservation: the wall impulse `J_wall` (integrated from the wall-pressure history) equals the
/// gas momentum change `p_final − p_initial` to the scheme's consistency tolerance — an `O(Δx)`
/// error that refines to zero (≈4e-4 at 400 cells), not round-off, because `J_wall` is sampled
/// independently of the momentum bookkeeping.
#[test]
fn elastic_bounce_conserves_momentum() {
    for &mach in &[3.0, 5.0, 8.0] {
        let r = bounce(400, mach);
        let rel = (r.wall_impulse - (r.residual_momentum + r.incident_momentum)).abs()
            / r.incident_momentum;
        assert!(
            rel < 1e-3,
            "M={mach}: momentum-conservation rel err = {rel:.2e}"
        );
    }
}

/// The lossless bounce ceiling sits strictly between dead-stick (`e_eff = 0`) and elastic
/// (`e_eff = 1`), and at high Mach approaches a γ-set, Mach-independent strong-shock value
/// (≈0.81 for γ = 5/3) — the "true bounce ceiling" reported by rung A (ADR-0001).
#[test]
fn bounce_ceiling_is_physical_and_mach_independent() {
    let e8 = bounce(400, 8.0).e_eff;
    let e12 = bounce(400, 12.0).e_eff;
    assert!(
        (0.5..1.0).contains(&e8),
        "ceiling outside the physical band (0.5, 1): {e8}"
    );
    assert!(
        (e8 - e12).abs() < 1e-2,
        "ceiling not Mach-converged: {e8} vs {e12}"
    );
}

/// The ceiling is resolution-converged (1st-order, shock-limited): doubling the grid shifts it by
/// well under 1%.
#[test]
fn bounce_ceiling_is_resolution_converged() {
    let e200 = bounce(200, 8.0).e_eff;
    let e400 = bounce(400, 8.0).e_eff;
    assert!(
        (e200 - e400).abs() < 5e-3,
        "ceiling not resolution-converged: {e200} vs {e400}"
    );
}

/// Dead-stick limit (`f → 0.5`): the idealized absorbing wall stops at stagnation, bringing the
/// gas to rest with no rebound, so by momentum conservation `J_wall = p_in` and `e_eff → 0`
/// (ADR-0001). The residual is the scheme's `O(Δx)` impulse-vs-momentum consistency, independent
/// of Mach (the stop is imposed, not gas-dynamic).
#[test]
fn dead_stick_absorbs_all_momentum() {
    for &mach in &[3.0, 5.0, 8.0] {
        let e = Tube::slug(400, mach, GAMMA).run_stick_bounce().e_eff;
        assert!(
            e.abs() < 1e-2,
            "M={mach}: dead-stick e_eff = {e:.2e}, want ≈ 0"
        );
    }
}

/// The **physical** peak wall pressure excludes the artificial-viscosity spike (ADR-0010
/// correction). For a γ = 1.4 slug at M = 5, the reflected shock that stagnates the gas gives
/// `p_peak ≈ (γ+1)/2 · ρ₀v² ≈ 1.2·ρ₀v²` (slightly above, at finite Mach) — while the reported
/// `peak_wall_force` (`p + q`) is dominated by the first-impact AV spike `≈ c_q·ρ₀v² = 2·ρ₀v²`
/// under the production VNR viscosity. The two must not be conflated.
#[test]
fn peak_wall_pressure_is_physical_not_av_spike() {
    let gamma = 1.4;
    // Tube::slug normalizes ρ₀ = 1, v = 1, so ρ₀v² = 1.
    let r = Tube::slug(400, 5.0, gamma).run_bounce();
    assert!(
        (1.0..1.5).contains(&r.peak_wall_pressure),
        "peak p(0)/ρv² = {:.3}, want ≈ (γ+1)/2 ≈ 1.2 (reflected shock)",
        r.peak_wall_pressure
    );
    assert!(
        r.peak_wall_force > 1.8,
        "peak (p+q)/ρv² = {:.3}, want ≳ c_q = 2 (AV first-impact spike)",
        r.peak_wall_force
    );
    assert!(
        r.peak_wall_pressure < r.peak_wall_force,
        "physical pressure peak must sit below the p+q peak"
    );
}
