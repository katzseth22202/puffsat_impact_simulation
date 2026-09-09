//! Real-water acceptance fixtures are generated independently by the Python
//! partition functions and reversible source, never by the Rust implementation.

use serde::Deserialize;
use water_plate::parcel::{Chemistry, Context, Knot, Settings, State, integrate};
use water_plate::{ResidualCurve, ResidualKnot};

#[derive(Deserialize)]
struct Check {
    temp: f64,
    rho: f64,
    y: [f64; 6],
    u: [f64; 6],
    cv: [f64; 6],
    log_reference: [f64; 6],
    source: [f64; 6],
}
#[derive(Deserialize)]
struct Relaxation {
    rho: f64,
    y0: [f64; 6],
    temperature0: f64,
    spectators: f64,
    equilibrium_y: [f64; 6],
    equilibrium_temperature: f64,
    duration: f64,
}
#[derive(Deserialize)]
struct Fixture {
    chemistry: Chemistry,
    checks: Vec<Check>,
    relaxation: Relaxation,
}

fn fixture() -> Fixture {
    serde_json::from_str(include_str!("fixtures/neutral_parcel.json")).unwrap()
}
fn residual() -> ResidualCurve {
    ResidualCurve {
        knots: [-20.0, 20.0]
            .map(|log_rho| ResidualKnot {
                log_rho,
                values: [0.0; 3],
                slopes: [0.0; 3],
            })
            .to_vec(),
    }
}

#[test]
fn rust_matches_independent_partition_energies_and_reversible_sources() {
    let f = fixture();
    for check in f.checks {
        let reference = f.chemistry.log_reference(check.temp);
        let (u, cv) = f.chemistry.species_thermo(check.temp);
        let (source, _, _) = f.chemistry.source(check.y, check.rho, check.temp);
        let scale = check.source.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for i in 0..6 {
            assert!((reference[i] - check.log_reference[i]).abs() < 1e-11);
            assert!((u[i] - check.u[i]).abs() < 1e-10 * u[i].abs().max(1.0));
            assert!((cv[i] - check.cv[i]).abs() < 1e-12);
            assert!((source[i] - check.source[i]).abs() < 1e-10 * scale);
        }
        let h = source[0] + 2.0 * source[2] + source[3] + 2.0 * source[5];
        let o = source[1] + source[3] + 2.0 * source[4] + source[5];
        assert!(h.abs().max(o.abs()) < 1e-14 * scale);
    }
}

#[test]
fn analytic_source_jacobian_includes_colliders_and_temperature() {
    let f = fixture();
    for check in f.checks {
        let (_, jac, ft) = f.chemistry.source(check.y, check.rho, check.temp);
        for j in 0..7 {
            let mut lo = check.y;
            let mut hi = check.y;
            let mut tl = check.temp;
            let mut th = check.temp;
            let dx = if j == 6 {
                check.temp * 1e-5
            } else {
                check.y[j] * 1e-5
            };
            if j == 6 {
                tl -= dx;
                th += dx;
            } else {
                lo[j] -= dx;
                hi[j] += dx;
            }
            let a = f.chemistry.source(lo, check.rho, tl).0;
            let b = f.chemistry.source(hi, check.rho, th).0;
            let measured: [f64; 6] = std::array::from_fn(|i| (b[i] - a[i]) / (2.0 * dx));
            let scale = measured.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
            for i in 0..6 {
                let expected = if j == 6 { ft[i] } else { jac[i][j] };
                assert!((measured[i] - expected).abs() < 2e-6 * scale);
            }
        }
    }
}

#[test]
fn finite_water_perturbation_returns_to_known_fixed_energy_equilibrium() {
    let f = fixture();
    let r = f.relaxation;
    let residual = residual();
    let y = r.y0;
    let c = Context {
        chem: &f.chemistry,
        residual: &residual,
        spectators: r.spectators,
        atoms: [
            y[0] + 2.0 * y[2] + y[3] + 2.0 * y[5],
            y[1] + y[3] + 2.0 * y[4] + y[5],
        ],
        energy_offset: 0.0,
        pressure_r_offset: 0.0,
    };
    let state = State {
        molecules: [y[2], y[3], y[4], y[5]],
        temperature: r.temperature0,
    };
    let knots = [
        Knot {
            time_s: 0.0,
            rho: r.rho,
            forcing_j_kg: 0.0,
        },
        Knot {
            time_s: r.duration,
            rho: r.rho,
            forcing_j_kg: 0.0,
        },
    ];
    let out = integrate(
        &c,
        state,
        &knots,
        Settings {
            relative_tolerance: 1e-5,
            absolute_population_tolerance: 1e-11,
            maximum_step_s: r.duration / 20.0,
            parent_forcing: false,
        },
    );
    assert_eq!(out.status, "completed");
    let last = out.samples.last().unwrap();
    assert!(last.energy_error_j_kg.abs() < 0.1);
    assert!((last.state.temperature - r.equilibrium_temperature).abs() < 1e-3);
    for (actual, expected) in c.populations(last.state).iter().zip(r.equilibrium_y) {
        assert!((actual - expected).abs() < 1e-7);
        assert!(*actual > 0.0);
    }
    assert!(last.net_bond_return_j_kg > 0.0);
}

#[test]
fn matched_equilibrium_does_not_spontaneously_release_bond_energy() {
    let f = fixture();
    let r = f.relaxation;
    let residual = residual();
    let y = r.equilibrium_y;
    let c = Context {
        chem: &f.chemistry,
        residual: &residual,
        spectators: r.spectators,
        atoms: [
            y[0] + 2.0 * y[2] + y[3] + 2.0 * y[5],
            y[1] + y[3] + 2.0 * y[4] + y[5],
        ],
        energy_offset: 0.0,
        pressure_r_offset: 0.0,
    };
    let state = State {
        molecules: [y[2], y[3], y[4], y[5]],
        temperature: r.equilibrium_temperature,
    };
    let knots = [
        Knot {
            time_s: 0.0,
            rho: r.rho,
            forcing_j_kg: 0.0,
        },
        Knot {
            time_s: r.duration,
            rho: r.rho,
            forcing_j_kg: 0.0,
        },
    ];
    let settings = Settings {
        relative_tolerance: 1e-5,
        absolute_population_tolerance: 1e-11,
        maximum_step_s: r.duration / 20.0,
        parent_forcing: false,
    };
    let out = integrate(&c, state, &knots, settings);
    assert_eq!(out.status, "completed");
    assert!(out.samples.last().unwrap().net_bond_return_j_kg.abs() < 0.01);
    let bad = integrate(
        &c,
        state,
        &knots,
        Settings {
            absolute_population_tolerance: f64::NAN,
            ..settings
        },
    );
    assert_ne!(bad.status, "completed");
    assert_eq!(bad.accepted_steps, 0);
}

#[test]
fn reaction_validation_rejects_negative_rates_and_broken_element_bookkeeping() {
    let mut chem = fixture().chemistry;
    assert!(chem.validate().is_ok());
    chem.reactions[0].a = -1.0;
    assert!(chem.validate().is_err());
    chem.reactions[0].a = 1.0;
    chem.reactions[0].products[0] += 1;
    assert!(chem.validate().is_err());
}
