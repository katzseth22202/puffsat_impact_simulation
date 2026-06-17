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
