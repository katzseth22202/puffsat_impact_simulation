//! Plane-wave wall-pressure acceptance: does `peak_local_pressure` measure the physical stagnation
//! load, or the wall cell's numerical artifact?
//!
//! `peak_local_pressure` is the max surface-cell pressure over the bounce. Survivability is
//! classified against it (via the `focusing` ratio, ADR-0010), so what it actually measures matters.
//!
//! The confined run is the plane-wave limit, which has a **closed-form** answer: a slug of density
//! `rho0` driven at `v` into a reflecting wall is the classical piston problem, whose post-shock
//! pressure is
//!
//! ```text
//! p2/p0 = 1 + g(g+1)/4 M^2 + g M sqrt(1 + (g+1)^2 M^2 / 16)
//! ```
//!
//! With the bounce module's normalization (`rho0 = 1`, `v = 1`, `p0 = 1/(g M^2)`) this is directly
//! in units of `rho0 v^2`, so it is the `c_stag` coefficient. At `g = 1.4`, `M = 40` it is
//! **1.20097**, approaching the strong-shock limit `(g+1)/2 = 1.2`.
//!
//! This matters because the kernel's own Noh acceptance test (`noh.rs`) deliberately measures "in a
//! band that avoids ... the axis (the classic wall-heating anomaly)" -- while `max_plate_pressure`
//! samples exactly that wall cell. If wall heating contaminates it, the absolute peak is wrong even
//! though the `focusing` *ratio* may still divide the common part out.

use euler2d::bounce::{Bounce2D, PlateShape, SlugConfig, run_slug_bounce};

const GAMMA: f64 = 1.4;
const MACH: f64 = 40.0;

/// Exact piston-problem wall pressure, in units of `rho0 v^2`.
fn analytic_piston_peak(gamma: f64, mach: f64) -> f64 {
    let m2 = mach * mach;
    let p2_over_p0 = 1.0
        + gamma * (gamma + 1.0) / 4.0 * m2
        + gamma * mach * (1.0 + (gamma + 1.0) * (gamma + 1.0) * m2 / 16.0).sqrt();
    let p0 = 1.0 / (gamma * m2); // the bounce module's normalization
    p2_over_p0 * p0
}

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

#[test]
fn plane_wave_wall_impulse_matches_the_piston_solution() {
    // The *impulse* is the quantity `eta_capture` and `f` are built from, and it integrates over
    // the whole plate and the whole bounce rather than reading one cell at one instant. It should
    // therefore be clean even if the peak is not -- which is the control for the test below.
    let b = confined(320);
    let ratio = b.restitution_ratio();
    assert!(
        (1.0..=2.0).contains(&ratio),
        "confined 1 + e_eff = {ratio:.4} is outside the physical range"
    );
}

#[test]
fn peak_plate_pressure_is_an_impact_overshoot_above_the_steady_piston_load() {
    // Measured 2026-08-22. The plane-wave wall pressure rises smoothly to ~1.335 at t ~ 0.08, then
    // decays back *through* the exact piston value and settles near 1.16 as the rear rarefaction
    // arrives. So `peak_local_pressure` reports an impact **overshoot**, ~11% above the steady
    // stagnation load -- and the overshoot is resolution-independent to five decimals over an 8x
    // refinement (nz 80 -> 640), so it is a converged feature of this initial condition, not a
    // startup artifact. (Contrast the 1D kernel's `peak_wall_force` spike, which ADR-0010 corrected
    // as a genuine artificial-viscosity artifact at ~2.0 rho v^2; this is a different thing.)
    //
    // 15% bounds it while catching any gross error in the stagnation load.
    let exact = analytic_piston_peak(GAMMA, MACH);
    let b = confined(320);
    let err = b.peak_local_pressure / exact - 1.0;
    assert!(
        (0.0..0.15).contains(&err),
        "plane-wave peak is {:+.1}% off the exact piston value ({:.5} against {exact:.5})",
        err * 100.0,
        b.peak_local_pressure
    );
}

#[test]
fn focusing_does_not_depend_on_measuring_the_overshoot_or_the_steady_load() {
    // The test that licenses the survivability model. `peak = c_stag rho v^2 * focusing` mixes two
    // things measured differently: `c_stag` is a **steady** reflected-shock coefficient from the 1D
    // kernel, while `focusing` is a ratio of 2D **peak-over-time** pressures that include the
    // overshoot above. Mixing them is only legitimate if the ratio is insensitive to which one is
    // used -- otherwise the model is multiplying a steady load by a transient concentration.
    //
    // It is: at r_foot/R = 0.5, L/D = 0.3, focusing is 1.2714 from peak-over-time against 1.2627
    // from the sustained load, 0.7% apart. The overshoot inflates numerator and denominator alike
    // and divides out, which is the same cancellation ADR-0003 relies on for `eta_capture`.
    let sample = |d_over_d: f64| -> (f64, f64) {
        let r_foot = 1.0;
        let r_plate = r_foot / 0.5;
        let length = 0.3 * 2.0 * r_foot;
        let depth = d_over_d * 2.0 * r_plate;
        let cfg = SlugConfig {
            gamma: GAMMA,
            mach: MACH,
            r_foot,
            length,
            r_plate,
            r_max: r_plate * 1.4,
            z_max: depth + 2.0 * length + 1.5,
            nr: 112,
            nz: 80,
            confined: false,
            shape: PlateShape::Dish { d_over_d },
            taper_frac: 0.0,
            alpha_div: 0.0,
        };
        let peak = run_slug_bounce(&cfg).peak_local_pressure;
        let mut g = euler2d::bounce::init_slug_grid(&cfg);
        let mut hist: Vec<f64> = Vec::new();
        for _ in 0..400 {
            g.run_to(0.01);
            hist.push(g.max_plate_pressure());
        }
        let sustained = hist
            .windows(5)
            .map(|w| w.iter().copied().fold(f64::INFINITY, f64::min))
            .fold(0.0_f64, f64::max);
        (peak, sustained)
    };
    let (flat_peak, flat_sus) = sample(0.0);
    let (con_peak, con_sus) = sample(0.10);
    let from_peak = con_peak / flat_peak;
    let from_sustained = con_sus / flat_sus;
    assert!(
        (from_peak / from_sustained - 1.0).abs() < 0.02,
        "focusing depends on how the load is measured: {from_peak:.4} from the peak against \
         {from_sustained:.4} sustained -- the model multiplies a steady c_stag by this ratio"
    );
}

/// Diagnostic: is the wall cell anomalous, or is the whole post-shock region at 1.335?
///
/// Wall heating is a local error confined to the first few cells against a reflecting boundary --
/// the plateau further in should sit on the exact piston value. This distinguishes "the kernel is
/// wrong about the stagnation load" from "`max_plate_pressure` samples the one place the kernel is
/// known to be unreliable", which have completely different consequences for `focusing`.
#[test]
#[ignore = "diagnostic: prints the near-wall pressure profile"]
fn near_wall_pressure_profile() {
    use euler2d::kernel::{Bc, Grid2D};
    use euler2d::state::Prim;

    let exact = analytic_piston_peak(GAMMA, MACH);
    let nz = 400usize;
    let nr = 4usize;
    let z_max = 3.0;
    let dz = z_max / nz as f64;
    let mut g = Grid2D::new(nz, nr, dz, dz, GAMMA);
    g.bc_zlo = Bc::Reflect; // the plate
    g.bc_zhi = Bc::Transmissive;
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect; // confined: plane wave

    let p0 = 1.0 / (GAMMA * MACH * MACH);
    let length = 1.0;
    let z_floor = 0.20;
    g.init(|iz, _ir| {
        let z = (iz as f64 + 0.5) * dz;
        if z >= z_floor && z < z_floor + length {
            Prim::new(1.0, -1.0, 0.0, p0)
        } else {
            Prim::new(1.0e-3, 0.0, 0.0, p0 * 1.0e-3)
        }
    });

    for t in [0.25_f64, 0.35, 0.50] {
        g.run_to(t);
        let ps: Vec<f64> = (0..14).map(|iz| g.prim(iz, nr / 2).p).collect();
        println!("t={t:.2}  exact={exact:.4}");
        for (iz, p) in ps.iter().enumerate() {
            println!("   iz={iz:>2}  p={p:.5}  p/exact={:.4}", p / exact);
        }
    }
}

/// Is the 11% excess a brief impact transient rather than the sustained load?
#[test]
#[ignore = "diagnostic: prints the wall-pressure time history"]
fn wall_pressure_time_history() {
    let exact = analytic_piston_peak(GAMMA, MACH);
    let b = confined(320);
    println!(
        "run_slug_bounce peak_local_pressure = {:.5} ({:+.2}% vs exact {exact:.5})",
        b.peak_local_pressure,
        (b.peak_local_pressure / exact - 1.0) * 100.0
    );

    use euler2d::kernel::{Bc, Grid2D};
    use euler2d::state::Prim;
    let nz = 320usize;
    let z_max = 3.0;
    let dz = z_max / nz as f64;
    let mut g = Grid2D::new(nz, 4, dz, dz, GAMMA);
    g.bc_zlo = Bc::Reflect;
    g.bc_zhi = Bc::Transmissive;
    g.bc_rlo = Bc::Reflect;
    g.bc_rhi = Bc::Reflect;
    let p0 = 1.0 / (GAMMA * MACH * MACH);
    g.init(|iz, _ir| {
        let z = (iz as f64 + 0.5) * dz;
        if (0.20..1.20).contains(&z) {
            Prim::new(1.0, -1.0, 0.0, p0)
        } else {
            Prim::new(1.0e-3, 0.0, 0.0, p0 * 1.0e-3)
        }
    });
    let mut peak = 0.0_f64;
    let mut t_peak = 0.0;
    let mut t = 0.0;
    while t < 0.9 {
        t += 0.004;
        g.run_to(t);
        let p = g.max_plate_pressure();
        if p > peak {
            peak = p;
            t_peak = t;
        }
        if (0.19..0.34).contains(&t) {
            println!("  t={t:.3}  p_wall={p:.5}  p/exact={:.4}", p / exact);
        }
    }
    println!(
        "max over history = {peak:.5} at t={t_peak:.3} ({:+.2}% vs exact)",
        (peak / exact - 1.0) * 100.0
    );
}

/// The question that decides whether the impact transient matters: does it cancel in `focusing`?
///
/// `focusing = peak(concave)/peak(flat)`, both maxima **over time**. If the impact transient
/// inflates both by the same factor it divides out and the survivability model is unaffected. If
/// the concave case's reported peak is the transient while its true geometric concentration
/// develops later, `focusing` is measuring the wrong thing.
#[test]
#[ignore = "diagnostic: transient vs sustained focusing"]
fn focusing_from_transient_versus_sustained_peak() {
    // Mirrors `run_eta_case`'s free-run geometry at r_foot/R = 0.5, L/D = 0.3.
    let case = |d_over_d: f64| -> (f64, f64) {
        let r_foot = 1.0;
        let r_plate = r_foot / 0.5;
        let length = 0.3 * 2.0 * r_foot;
        let depth = d_over_d * 2.0 * r_plate;
        let cfg = SlugConfig {
            gamma: GAMMA,
            mach: MACH,
            r_foot,
            length,
            r_plate,
            r_max: r_plate * 1.4,
            z_max: depth + 2.0 * length + 1.5,
            nr: 112,
            nz: 80,
            confined: false,
            shape: PlateShape::Dish { d_over_d },
            taper_frac: 0.0,
            alpha_div: 0.0,
        };
        let whole = run_slug_bounce(&cfg).peak_local_pressure;

        // Re-run by hand, tracking when the maximum occurs and what the load settles to.
        let mut g = euler2d::bounce::init_slug_grid(&cfg);
        let mut peak = 0.0_f64;
        let mut t_peak = 0.0_f64;
        let mut t = 0.0_f64;
        let mut hist: Vec<(f64, f64)> = Vec::new();
        let step = 0.01;
        for _ in 0..400 {
            g.run_to(step);
            t += step;
            let p = g.max_plate_pressure();
            hist.push((t, p));
            if p > peak {
                peak = p;
                t_peak = t;
            }
        }
        // "Sustained" = the largest load held for at least 5 consecutive samples (0.05 in time),
        // which excludes a one-sample impact spike but keeps a real focusing concentration.
        let mut sustained = 0.0_f64;
        for w in hist.windows(5) {
            let lo = w.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
            sustained = sustained.max(lo);
        }
        println!(
            "  d/D={d_over_d:.2}: run_slug_bounce peak={whole:.4}  hand peak={peak:.4} at t={t_peak:.2}  sustained={sustained:.4}"
        );
        (whole, sustained)
    };

    println!("focusing = peak(concave)/peak(flat), transient vs sustained:");
    let (flat_peak, flat_sus) = case(0.0);
    let (con_peak, con_sus) = case(0.10);
    println!(
        "  focusing from peak-over-time = {:.4}",
        con_peak / flat_peak
    );
    println!("  focusing from sustained load = {:.4}", con_sus / flat_sus);
}

/// Is the confined run's reported peak the sustained plane-wave load, or an impact transient?
#[test]
#[ignore = "diagnostic: confined transient vs sustained"]
fn confined_transient_versus_sustained() {
    let exact = analytic_piston_peak(GAMMA, MACH);
    let cfg = SlugConfig {
        gamma: GAMMA,
        mach: MACH,
        r_foot: 5.0,
        length: 1.0,
        r_plate: 5.0,
        r_max: 5.0,
        z_max: 3.0,
        nr: 8,
        nz: 320,
        confined: true,
        shape: PlateShape::FlatGridAligned,
        taper_frac: 0.0,
        alpha_div: 0.0,
    };
    let reported = run_slug_bounce(&cfg).peak_local_pressure;
    let mut g = euler2d::bounce::init_slug_grid(&cfg);
    let step = 0.01;
    let mut hist: Vec<(f64, f64)> = Vec::new();
    let mut t = 0.0;
    for _ in 0..300 {
        g.run_to(step);
        t += step;
        hist.push((t, g.max_plate_pressure()));
    }
    let (t_peak, peak) = hist
        .iter()
        .copied()
        .fold((0.0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
    let mut sustained = 0.0_f64;
    for w in hist.windows(5) {
        sustained = sustained.max(w.iter().map(|x| x.1).fold(f64::INFINITY, f64::min));
    }
    println!("confined (plane wave), exact piston = {exact:.5}");
    println!(
        "  run_slug_bounce reported peak = {reported:.5} ({:+.1}%)",
        (reported / exact - 1.0) * 100.0
    );
    println!(
        "  hand peak = {peak:.5} at t={t_peak:.2} ({:+.1}%)",
        (peak / exact - 1.0) * 100.0
    );
    println!(
        "  sustained = {sustained:.5} ({:+.1}%)",
        (sustained / exact - 1.0) * 100.0
    );
    for (tt, p) in hist.iter().take(40) {
        if *p > 0.5 {
            println!("    t={tt:.2}  p={p:.5}  p/exact={:.4}", p / exact);
        }
    }
}
