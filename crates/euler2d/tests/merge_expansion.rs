//! N1's answer: the anisotropy a freely expanding merged fireball actually reaches.
//!
//! `eq:reflection_baseline` assumes the thermal remainder is isotropic (`alpha = 1/3`).
//! `sec:minimum_nozzle` says the stagnation region "squirts radially into the field". The two
//! analytic limits — shape-following (`alpha = 0.91` for the flown column) and gradient-following
//! (`alpha < 1/3`) — bracket the answer three-fold in the baseline, so only a run settles it.
//!
//! # The measurement window, which is narrow and asserted rather than assumed
//!
//! Two clocks run against each other. The expansion has to finish (internal energy converted to
//! motion) before `alpha` means what the baseline needs; and the fireball has to stay inside the
//! domain, because a transmissive boundary silently deletes the fastest material — which is
//! exactly the material carrying the anisotropy. Those windows overlap at `t/t_c` ≈ 1.5–5.0,
//! and every test here is taken inside it.
//!
//! An earlier version of this file measured at `t/t_c` = 4, where 93% of the mass had already
//! left, and its mass check passed only because a staircased initial cylinder was over-massed by
//! the same amount it then lost. Both errors are pinned below.

use euler2d::merge::{MergeConfig, init_merge_grid, radial_sound_crossing};
use euler2d::moments::{Moments, cosine_histogram, moments};

/// The measurement point: expansion essentially complete, mass still fully inside.
const MEASURE_AT: f64 = 4.0;

fn expand(cfg: &MergeConfig, crossings: f64) -> Moments {
    let mut g = init_merge_grid(cfg);
    g.run_to(crossings * radial_sound_crossing(cfg));
    moments(&g)
}

#[test]
fn the_initial_state_carries_the_energy_budget_it_was_given() {
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    // 0.361, not the paper's 0.323: `BAG_RHO` is the 213 kg slug over 660 m³, before the
    // projectile is added. The *merged* 238 kg in the same volume is denser by (1+k)/k.
    assert!(
        (cfg.column_density() - 0.361).abs() < 0.005,
        "{}",
        cfg.column_density()
    );
    assert!(
        (cfg.pulse_energy() - 62.9e9).abs() / 62.9e9 < 0.01,
        "{}",
        cfg.pulse_energy()
    );
    // The flown masses give k = 213/25 = 8.52, not the 8.5 the paper states -- a 0.2%
    // inconsistency in the paper's own numbers, and harmless, but it is the masses that
    // define the physics here so they are what this pins.
    assert!(
        (cfg.slug_ratio() - 8.52).abs() < 1e-12,
        "{}",
        cfg.slug_ratio()
    );
    assert!((cfg.drift_fraction() - 1.0 / 9.52).abs() < 1e-12);

    let m = moments(&init_merge_grid(&cfg));
    assert!((m.mass - 238.0).abs() / 238.0 < 0.01, "mass {}", m.mass);
    assert!(
        m.unconverted_fraction() > 0.8,
        "{}",
        m.unconverted_fraction()
    );
    assert!((m.v_cm_z - cfg.drift_speed()).abs() / cfg.drift_speed() < 0.01);
}

#[test]
fn axial_momentum_and_mass_are_conserved_through_the_measurement_window() {
    // The validity gate. Compared against the run's OWN initial mass, not the nominal 238 kg —
    // comparing to the nominal let a staircase error and a boundary loss cancel each other.
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    let m0 = moments(&init_merge_grid(&cfg));
    let m1 = expand(&cfg, MEASURE_AT);
    assert!(
        (m1.mass - m0.mass).abs() / m0.mass < 1e-3,
        "{} -> {}",
        m0.mass,
        m1.mass
    );
    let p0 = m0.mass * m0.v_cm_z;
    let p1 = m1.mass * m1.v_cm_z;
    assert!((p1 - p0).abs() / p0 < 1e-3, "axial momentum {p0} -> {p1}");
}

#[test]
fn the_transverse_mean_is_an_expansion_speed_and_not_a_drift() {
    // In (z, r) coordinates `u_r >= 0` everywhere for an expanding cloud, so its mass-weighted mean
    // is large and positive — but the true 3D transverse momentum is zero by symmetry. Folding it
    // into the drift energy made a freely expanding fireball look like it carried six times its own
    // bulk momentum, and produced f_d = 0.78 against an input of 0.105.
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    let m = expand(&cfg, MEASURE_AT);
    assert!(
        m.v_expansion_r > 10_000.0,
        "expansion should be vigorous: {}",
        m.v_expansion_r
    );
    assert!(
        (m.drift_fraction_kinetic - cfg.drift_fraction()).abs() < 0.05,
        "drift share {} should track the input f_d {}",
        m.drift_fraction_kinetic,
        cfg.drift_fraction()
    );
}

#[test]
fn the_expansion_converts_pressure_into_motion() {
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    let early = expand(&cfg, 0.4);
    let measured = expand(&cfg, MEASURE_AT);
    assert!(measured.unconverted_fraction() < early.unconverted_fraction());
    assert!(
        measured.unconverted_fraction() < 0.15,
        "expansion unfinished at the measurement point: {}",
        measured.unconverted_fraction()
    );
}

#[test]
fn the_fireball_is_a_pancake_not_the_assumed_sphere() {
    // **The finding.** Free expansion lands on the gradient-following limit, not the
    // shape-following one and not the paper's isotropic assumption.
    let m = expand(&MergeConfig::flown(75_000.0, 700, 200), MEASURE_AT);
    assert!(
        m.alpha < 0.15,
        "alpha = {} — expected a strong pancake",
        m.alpha
    );
    assert!(
        m.alpha > 0.05,
        "alpha = {} — implausibly flat, suspect the run",
        m.alpha
    );
}

#[test]
fn alpha_is_stationary_across_the_measurement_window() {
    // If it were still moving, the reported number would be a timestamp rather than an answer.
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    let a = expand(&cfg, 1.5).alpha;
    let b = expand(&cfg, 3.0).alpha;
    let c = expand(&cfg, 5.0).alpha;
    // Inside the window alpha is flat to well under the difference from the assumed 1/3.
    // The third decimal is not settled -- it moves a little with the timestep sequence -- so
    // this pins the band the finding actually rests on rather than a spurious precision.
    for x in [a, b, c] {
        assert!(
            (0.07..0.12).contains(&x),
            "alpha left its band: {a} {b} {c}"
        );
    }
    assert!(
        (b - a).abs() < 0.01,
        "still moving inside the window: {a} {b}"
    );
}

#[test]
fn alpha_is_converged_under_mesh_refinement() {
    let coarse = expand(&MergeConfig::flown(75_000.0, 350, 100), MEASURE_AT).alpha;
    let medium = expand(&MergeConfig::flown(75_000.0, 525, 150), MEASURE_AT).alpha;
    let fine = expand(&MergeConfig::flown(75_000.0, 700, 200), MEASURE_AT).alpha;
    let d1 = (medium - coarse).abs();
    let d2 = (fine - medium).abs();
    assert!(
        d2 <= d1.max(0.01),
        "not converging: {coarse} -> {medium} -> {fine}"
    );
}

#[test]
fn the_ambient_stand_in_for_vacuum_does_not_set_the_answer() {
    let mut thin = MergeConfig::flown(75_000.0, 700, 200);
    thin.ambient_fraction = 1e-7;
    let mut thick = MergeConfig::flown(75_000.0, 700, 200);
    thick.ambient_fraction = 1e-5;
    let a = expand(&thin, MEASURE_AT).alpha;
    let b = expand(&thick, MEASURE_AT).alpha;
    assert!((a - b).abs() < 0.02, "ambient set alpha: {a} vs {b}");
}

#[test]
fn alpha_does_not_depend_on_closing_speed() {
    // `alpha` is a shape, and the geometry is identical across the burn — only the energy scale
    // differs, and `radial_sound_crossing` normalises that out. Reporting four speeds as if they
    // were independent measurements would be presenting one number four times.
    let hot = expand(&MergeConfig::flown(75_000.0, 700, 200), MEASURE_AT).alpha;
    let cold = expand(&MergeConfig::flown(45_580.0, 700, 200), MEASURE_AT).alpha;
    // Agreement to ~8 significant figures; the residue is float accumulation over
    // different timestep sequences, not physics.
    assert!((hot - cold).abs() < 1e-6, "{hot} vs {cold}");
}

/// The deliverable. `cargo test --release -p euler2d --test merge_expansion -- --nocapture report`
#[test]
fn report() {
    let cfg = MergeConfig::flown(75_000.0, 700, 200);
    let tc = radial_sound_crossing(&cfg);
    let m0 = moments(&init_merge_grid(&cfg));
    println!("\n== N1: alpha of the freely expanding merged fireball ==");
    println!(
        "column {} m x {} m, {:.1} kg at {:.4} kg/m^3, radial sound crossing {:.3} ms",
        cfg.column_length,
        cfg.column_radius,
        m0.mass,
        cfg.column_density(),
        tc * 1e3
    );
    println!(
        "brackets: shape-following 0.906 | isotropic (assumed) 0.333 | gradient-following < 0.333\n"
    );
    println!(
        "{:>6} {:>9} {:>8} {:>9} {:>8} {:>8}",
        "t/tc", "mass kg", "alpha", "baseline", "unconv", "f_d"
    );
    for i in [1, 2, 3, 4, 5] {
        let t = f64::from(i) * 0.3;
        let m = expand(&cfg, t);
        println!(
            "{:6.2} {:9.2} {:8.4} {:9.4} {:8.4} {:8.4}",
            t,
            m.mass,
            m.alpha,
            (2.0 * m.alpha / std::f64::consts::PI).sqrt(),
            m.unconverted_fraction(),
            m.drift_fraction_kinetic
        );
    }

    let m = expand(&cfg, MEASURE_AT);
    let baseline = (2.0 * m.alpha / std::f64::consts::PI).sqrt();
    println!(
        "\nat t/tc = {MEASURE_AT}: alpha = {:.4}, baseline = {:.4}",
        m.alpha, baseline
    );
    println!(
        "against the assumed alpha = 1/3, baseline 0.4607 -> a {:.0}% reduction",
        100.0 * (1.0 - baseline / 0.4607)
    );

    let sphere = expand(&MergeConfig::spherical_bag(75_000.0, 700, 200), MEASURE_AT);
    println!(
        "\nCONTROL -- same 238 kg, same density, same energy, in the 5.4 m SPHERE the paper\n\
         started with before the launch fairing stretched it into a column:"
    );
    println!(
        "  sphere: alpha = {:.4}, baseline = {:.4}",
        sphere.alpha,
        (2.0 * sphere.alpha / std::f64::consts::PI).sqrt()
    );
    println!(
        "  column: alpha = {:.4}, baseline = {:.4}",
        m.alpha, baseline
    );
    println!(
        "  -> the shape is what does it, by a factor of {:.1}x in alpha.",
        sphere.alpha / m.alpha
    );

    let mut g = init_merge_grid(&cfg);
    g.run_to(MEASURE_AT * tc);
    println!("\ncos(theta) = v_z/|v| histogram for the COLUMN (flat = isotropic, 0.100/bin):");
    for (i, h) in cosine_histogram(&g, 10).iter().enumerate() {
        let lo = -1.0 + 0.2 * f64::from(i32::try_from(i).unwrap_or(0));
        let bar = "#".repeat(((h * 200.0) as usize).min(120));
        println!("  [{lo:+.1}, {:+.1})  {h:.4}  {bar}", lo + 0.2);
    }
}

#[test]
fn a_spherical_bag_of_the_same_mass_and_energy_is_nearly_isotropic() {
    // **The control, and the point of the whole finding.** Same 238 kg, same density, same energy
    // — only the shape differs. A sphere expands isotropically by symmetry, so it must return
    // `alpha` near 1/3, and the column must not.
    //
    // The paper's bag *was* a 5.4 m sphere and became a 23 m column because "length is the cheap
    // dimension for a rocket to carry and diameter is the expensive one". That is a launch-fairing
    // decision, and it is what breaks `eq:reflection_baseline`'s isotropy assumption.
    let sphere = expand(&MergeConfig::spherical_bag(75_000.0, 700, 200), MEASURE_AT);
    let column = expand(&MergeConfig::flown(75_000.0, 700, 200), MEASURE_AT);
    assert!(
        (sphere.alpha - 1.0 / 3.0).abs() < 0.05,
        "the sphere should be isotropic: alpha = {}",
        sphere.alpha
    );
    assert!(
        sphere.alpha > 3.0 * column.alpha,
        "shape must be what separates them: sphere {} vs column {}",
        sphere.alpha,
        column.alpha
    );
}
