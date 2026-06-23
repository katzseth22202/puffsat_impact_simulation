//! Acceptance tests for the ghost-cell immersed-boundary reflecting wall (ADR-0023 amendment, D4).
//!
//! The curved/inclined plate is imposed by mirroring fluid into the solid cells across the **true**
//! local surface normal. Two kinematic gates pin that the true normal (not a grid-aligned staircase)
//! is used — exactly the property `eta_capture` (a rebound-angle measurement) depends on:
//!
//! - **free-slip tangency** — uniform flow *parallel* to an inclined wall stays uniform (a staircase
//!   would inject spurious normal velocity at every step);
//! - **normal-only wall impulse** — a gas blob bounced off an inclined wall changes its total
//!   momentum *along the wall normal* (a free-slip wall exerts no tangential force); a grid-aligned
//!   reflect would push it straight up (`∥ ẑ`) instead.
//!
//! Plus a `plate_force` unit check: the axial wall force is the surface pressure on the *projected*
//! annulus regardless of tilt.

use euler2d::kernel::Grid2D;
use euler2d::plate::PlateProfile;
use euler2d::state::Prim;

const GAMMA: f64 = 1.4;

/// Sum `(mz, mr)` over fluid (non-solid) cells — direction is all we need, so the uniform planar
/// cell volume `dz·dr` is dropped.
fn total_momentum(g: &Grid2D) -> (f64, f64) {
    let mut mz = 0.0;
    let mut mr = 0.0;
    for iz in 0..g.nz() {
        for ir in 0..g.nr() {
            if !g.is_solid(iz, ir) {
                mz += g.cons(iz, ir).mz;
                mr += g.cons(iz, ir).mr;
            }
        }
    }
    (mz, mr)
}

/// Total fluid mass (drops the uniform cell volume) — the escape gate for the impulse test.
fn fluid_mass(g: &Grid2D) -> f64 {
    let mut m = 0.0;
    for iz in 0..g.nz() {
        for ir in 0..g.nr() {
            if !g.is_solid(iz, ir) {
                m += g.cons(iz, ir).rho;
            }
        }
    }
    m
}

/// **D4a — free-slip tangency.** A uniform stream directed *along* an inclined flat wall is a steady
/// state of a true-normal reflecting boundary: the mirror is the identity for tangent flow, so the
/// field must stay uniform. A staircase boundary would deflect it and ripple the pressure.
#[test]
fn free_slip_tangency_keeps_a_parallel_stream_uniform() {
    let (nz, nr) = (40, 20);
    let (dz, dr) = (0.05, 0.05);
    let slope = 0.5;
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.bc_zlo = euler2d::kernel::Bc::Transmissive; // the immersed surface is the wall, not z=0
    g.bc_zhi = euler2d::kernel::Bc::Transmissive;
    g.bc_rlo = euler2d::kernel::Bc::Transmissive;
    g.bc_rhi = euler2d::kernel::Bc::Transmissive;
    g.set_plate_profile(Some(PlateProfile::InclinedPlane { z0: 0.2, slope }));

    // Velocity along the wall tangent t̂ = (slope, 1)/√(1+slope²); speed 1.
    let inv = 1.0 / (1.0 + slope * slope).sqrt();
    let (uz0, ur0) = (slope * inv, inv);
    g.init(|_, _| Prim::new(1.0, uz0, ur0, 1.0));

    for _ in 0..50 {
        let dt = g.stable_dt();
        g.step(dt);
    }

    // Interior fluid cells must be unchanged (away from the open boundaries where outflow noise can
    // accumulate).
    for iz in 4..nz - 4 {
        for ir in 2..nr - 2 {
            if g.is_solid(iz, ir) {
                continue;
            }
            let w = g.prim(iz, ir);
            assert!((w.rho - 1.0).abs() < 1e-6, "rho moved: {}", w.rho);
            assert!((w.p - 1.0).abs() < 1e-6, "p moved: {}", w.p);
            assert!(
                (w.uz - uz0).abs() < 1e-6 && (w.ur - ur0).abs() < 1e-6,
                "velocity moved"
            );
        }
    }
}

/// **D4b — normal-only wall impulse (the anti-staircase gate).** A gas blob dropped straight down
/// onto an inclined flat wall can only receive a *normal* impulse (a free-slip wall exerts no
/// tangential force), so the change in total fluid momentum is `∥ n̂`. A grid-aligned reflect would
/// give `∥ ẑ`, whose tangential fraction is `slope` — here 0.5 — so the < 0.15 gate cleanly
/// separates the true-normal mirror from a staircase.
#[test]
fn wall_impulse_is_along_the_true_normal() {
    let (nz, nr) = (50, 40);
    let (dz, dr) = (0.05, 0.05); // domain 2.5 × 2.0
    let slope = 0.5;
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.bc_zlo = euler2d::kernel::Bc::Transmissive;
    g.bc_zhi = euler2d::kernel::Bc::Transmissive;
    g.bc_rlo = euler2d::kernel::Bc::Transmissive;
    g.bc_rhi = euler2d::kernel::Bc::Transmissive;
    g.set_plate_profile(Some(PlateProfile::InclinedPlane { z0: 0.2, slope }));

    // A finite blob, away from both r-boundaries, moving down in near-vacuum.
    let v = 3.0;
    g.init(|iz, ir| {
        let z = (iz as f64 + 0.5) * dz;
        let r = (ir as f64 + 0.5) * dr;
        if (1.0..1.5).contains(&z) && (0.6..1.1).contains(&r) {
            Prim::new(1.0, -v, 0.0, 1.0)
        } else {
            Prim::new(1.0e-3, 0.0, 0.0, 1.0e-3)
        }
    });

    let p0 = total_momentum(&g);
    let m0 = fluid_mass(&g);

    // Integrate through the compression phase and stop once the gas peels back off the wall (force
    // past its peak and falling) — the normal impulse is delivered by then, before the gas slides
    // downhill far enough to spill off the wall edge.
    let mut peak = 0.0_f64;
    let mut past_peak = false;
    for _ in 0..400 {
        let dt = g.stable_dt();
        g.step(dt);
        let f = g.plate_force();
        peak = peak.max(f);
        if f < 0.999 * peak {
            past_peak = true;
        }
        if past_peak && f < 0.6 * peak {
            break;
        }
    }

    // Escape gate: the bounce must still be contained (else boundary outflow would corrupt Δp).
    let m1 = fluid_mass(&g);
    assert!(
        (m1 - m0).abs() / m0 < 0.03,
        "mass escaped ({m0:.4} → {m1:.4}); widen the domain"
    );

    let p1 = total_momentum(&g);
    let (dpz, dpr) = (p1.0 - p0.0, p1.1 - p0.1);
    let inv = 1.0 / (1.0 + slope * slope).sqrt();
    let (nz_, nr_) = (inv, -slope * inv); // n̂
    let (tz_, tr_) = (slope * inv, inv); // t̂
    let dp_n = dpz * nz_ + dpr * nr_;
    let dp_t = dpz * tz_ + dpr * tr_;
    assert!(dp_n.abs() > 1e-3, "no measurable wall impulse");
    let tangential_fraction = dp_t.abs() / dp_n.abs();
    assert!(
        tangential_fraction < 0.15,
        "wall impulse not normal: |Δp·t̂|/|Δp·n̂| = {tangential_fraction:.3} (staircase ≈ {slope})"
    );
}

/// **`plate_force` projection-cancellation.** The axial force on a tilted surface is the surface
/// pressure integrated over the *projected* annulus — the `√(1+s²)` of the arc length cancels the
/// `1/√(1+s²)` of `n̂·ẑ`. So for a uniform pressure the dish's axial force equals `p·Σ_{r≤r_plate}
/// r·dr`, independent of depth; this also checks that `plate_force` reads the surface (first fluid)
/// cell, not `iz = 0` (which is solid under the raised dish).
#[test]
fn plate_force_uses_surface_cell_and_cancels_tilt() {
    let (nz, nr) = (30, 20);
    let (dz, dr) = (0.1, 0.1);
    let r_plate = 1.5;
    let p_uniform = 2.0;
    let mut g = Grid2D::new(nz, nr, dz, dr, GAMMA);
    g.set_plate_profile(Some(PlateProfile::Dish {
        r_plate,
        z0: 0.2,
        depth: 0.3,
    }));
    g.init(|_, _| Prim::new(1.0, 0.0, 0.0, p_uniform));

    let expected: f64 = (0..nr)
        .map(|ir| (ir as f64 + 0.5) * dr)
        .filter(|&r| r <= r_plate)
        .map(|r| p_uniform * r * dr)
        .sum();
    assert!(
        (g.plate_force() - expected).abs() < 1e-12,
        "plate_force {} vs expected {expected}",
        g.plate_force()
    );
}
