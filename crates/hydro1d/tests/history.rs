//! Conservative trajectory access for the water-layer reference: force impulse
//! must equal nodal momentum change, and a phase energy offset must be inert.

use approx::assert_relative_eq;
use hydro1d::eos::{Eos, IdealGas};
use hydro1d::kernel::{Boundary, LagrangianState, Tube, Viscosity};

#[derive(Debug, Clone)]
struct ShiftedGas {
    shift: f64,
}

impl Eos for ShiftedGas {
    fn pressure(&self, rho: f64, e: f64) -> f64 {
        0.4 * rho * (e - self.shift)
    }
    fn sound_speed(&self, _rho: f64, e: f64) -> f64 {
        (1.4 * 0.4 * (e - self.shift)).sqrt()
    }
    fn energy_from_pressure(&self, rho: f64, p: f64) -> f64 {
        p / (0.4 * rho) + self.shift
    }
    fn dp_de(&self, rho: f64, _e: f64) -> f64 {
        0.4 * rho
    }
    fn temperature(&self, _rho: f64, e: f64) -> f64 {
        e - self.shift
    }
}

fn initial(shift: f64) -> LagrangianState {
    let n = 40;
    LagrangianState {
        positions: (0..=n).map(|i| f64::from(i) / f64::from(n)).collect(),
        node_velocities: vec![0.0; 41],
        cell_masses: vec![1.0 / 40.0; 40],
        specific_energies: vec![2e6 + shift; 40],
    }
}

#[test]
fn integrated_boundary_forces_close_the_momentum_ledger() {
    let mut tube = Tube::from_lagrangian(
        initial(0.0),
        IdealGas::new(1.4),
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let start = tube.total_momentum();
    let mut impulse = 0.0;
    for _ in 0..100 {
        let step = tube.advance(1e-5);
        impulse += step.left_wall_impulse - step.right_wall_impulse;
    }
    assert_relative_eq!(tube.total_mass(), 1.0, epsilon = 1e-12);
    assert_relative_eq!(tube.total_momentum() - start, impulse, max_relative = 1e-12);
    assert!((tube.total_energy() / 2e6 - 1.0).abs() < 0.002);
}

#[test]
fn a_common_energy_shift_does_not_change_the_expansion() {
    let run = |shift| {
        let mut tube = Tube::from_lagrangian(
            initial(shift),
            ShiftedGas { shift },
            Boundary::Wall,
            Boundary::Free,
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        tube.run_to(1e-4);
        tube
    };
    let plain = run(0.0);
    let shifted = run(3e6);
    for j in 0..plain.cells() {
        assert_relative_eq!(plain.density(j), shifted.density(j), max_relative = 1e-10);
        assert_relative_eq!(plain.pressure(j), shifted.pressure(j), max_relative = 1e-10);
    }
    assert_relative_eq!(
        plain.total_energy() + 3e6,
        shifted.total_energy(),
        max_relative = 1e-12
    );
}

#[test]
fn smaller_cfl_steps_reduce_shock_energy_drift() {
    let drift = |fraction| {
        let mut state = initial(0.0);
        state.node_velocities[1..].fill(-5000.0);
        let mut tube = Tube::from_lagrangian(
            state,
            IdealGas::new(1.4),
            Boundary::Wall,
            Boundary::Free,
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        let energy = tube.total_energy();
        let mut time = 0.0;
        let mut worst: f64 = 0.0;
        while time < 0.001 {
            time += tube.advance_with_cfl_fraction(0.001 - time, fraction).dt;
            worst = worst.max((tube.total_energy() / energy - 1.0).abs());
        }
        worst
    };
    let coarse = drift(1.0);
    let medium = drift(0.5);
    let fine = drift(0.25);
    assert!(medium < 0.7 * coarse, "{coarse} -> {medium}");
    assert!(fine < 0.7 * medium, "{medium} -> {fine}");
    assert!(fine < 0.005, "{fine}");
}
