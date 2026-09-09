//! Matched same-state equilibrium and parcel-frozen restart experiment.

use hydro1d::eos::{Eos, TableEos};
use hydro1d::kernel::{Boundary, LagrangianState, Tube, Viscosity};
use serde::Deserialize;
use serde_json::json;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;
use water_plate::{FrozenCell, FrozenGas, ParcelEos, ResidualCurve};

#[derive(Debug, Deserialize)]
struct State {
    positions: Vec<f64>,
    node_velocities: Vec<f64>,
    cell_masses: Vec<f64>,
    specific_energies: Vec<f64>,
}
impl State {
    fn copy(&self) -> LagrangianState {
        LagrangianState {
            positions: self.positions.clone(),
            node_velocities: self.node_velocities.clone(),
            cell_masses: self.cell_masses.clone(),
            specific_energies: self.specific_energies.clone(),
        }
    }
}
#[derive(Debug, Deserialize)]
struct Parcel {
    gas: FrozenGas,
    ionizing: Option<water_plate::selective::IonizingGas>,
    temperature_k: f64,
    pressure_pa: f64,
    expansion_rate_s: f64,
    source_pressure_pa: f64,
    source_energy_j_kg: f64,
}
#[derive(Debug, Deserialize)]
struct Restart {
    name: String,
    start_time_s: f64,
    initial_wall_impulse_ns: f64,
    area_m2: f64,
    storage_shift: f64,
    incoming_momentum_ns: f64,
    timestep_fraction: f64,
    state: State,
    parcels: Vec<Parcel>,
    output_times_s: Vec<f64>,
}
#[derive(Debug, Deserialize)]
struct Inputs {
    residual_curve: ResidualCurve,
    restarts: Vec<Restart>,
}

fn row(out: &mut impl Write, value: &serde_json::Value) -> Result<(), Box<dyn Error>> {
    serde_json::to_writer(&mut *out, value)?;
    writeln!(out)?;
    Ok(())
}

fn domain_error(tube: &Tube<ParcelEos>, eos: &ParcelEos, table: &TableEos) -> Option<String> {
    let state = tube.lagrangian_state();
    for (j, &energy) in state.specific_energies.iter().enumerate() {
        let rho = tube.density(j);
        if !rho.is_finite()
            || rho < table.table().rho_grid()[0]
            || rho > *table.table().rho_grid().last().unwrap()
        {
            return Some(format!("density domain: cell {j}, rho={rho}"));
        }
        let (lo, hi) = match eos.for_cell(j) {
            ParcelEos::Equilibrium(t) => (t.table().energy(rho, 273.2), t.table().energy(rho, 2e6)),
            ParcelEos::Frozen(t) => (
                t.at_temperature(rho, 1000.0).energy,
                t.at_temperature(rho, 2e6).energy,
            ),
            ParcelEos::Indexed(_) => unreachable!(),
        };
        if !energy.is_finite() || energy < lo * (1.0 - 1e-9) || energy > hi * (1.0 + 1e-9) {
            return Some(format!(
                "energy domain: cell {j}, e={energy}, limits=[{lo},{hi}]"
            ));
        }
        if let ParcelEos::Frozen(t) = eos.for_cell(j) {
            let s = t.at_temperature(rho, t.invert_energy(rho, energy));
            let cs2 = s.pressure_r + s.pressure_t * (s.pressure / (rho * rho) - s.energy_r) / s.cv;
            if !(s.pressure > 0.0 && s.cv > 0.0 && cs2 > 0.0 && cs2.is_finite()) {
                return Some(format!("unstable frozen state: cell {j}, {s:?}"));
            }
        }
    }
    None
}

fn run(
    input: &Restart,
    table: &Arc<TableEos>,
    curve: &Arc<ResidualCurve>,
    frozen: bool,
    molecular: bool,
    out: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    assert_eq!(input.parcels.len(), input.state.cell_masses.len());
    assert!((input.output_times_s[0] - input.start_time_s).abs() < 1e-15);
    assert!(input.output_times_s.windows(2).all(|w| w[0] < w[1]));
    let mut materials = Vec::new();
    let mut abs_energy_offset = 0.0;
    let mut max_pressure_correction: f64 = 0.0;
    let mut max_source_pressure_error: f64 = 0.0;
    let mut max_source_energy_error: f64 = 0.0;
    let mut max_start_pressure_jump: f64 = 0.0;
    let mut max_start_temperature_jump: f64 = 0.0;
    let mut bond = 0.0;
    let mut ion = 0.0;
    let mut compressing_mass = 0.0;
    for (j, p) in input.parcels.iter().enumerate() {
        let mass = input.state.cell_masses[j] * input.area_m2;
        let rho =
            input.state.cell_masses[j] / (input.state.positions[j + 1] - input.state.positions[j]);
        let e = input.state.specific_energies[j];
        bond += mass * p.gas.bond_store;
        ion += mass * p.gas.ion_store;
        if p.expansion_rate_s < 0.0 {
            compressing_mass += mass;
        }
        if frozen {
            let mut cell = FrozenCell {
                gas: p.gas.clone(),
                ionizing: if molecular {
                    Some(p.ionizing.clone().ok_or("missing selective coefficients")?)
                } else {
                    None
                },
                residual: Arc::clone(curve),
                storage_shift: input.storage_shift,
                energy_offset: 0.0,
                pressure_r_offset: 0.0,
            };
            let original = cell.at_temperature(rho, p.temperature_k);
            // These compare the independent Rust potential/thermal evaluator
            // with Python's source EOS before matching the tabulated parent.
            max_source_pressure_error = max_source_pressure_error
                .max((original.pressure / p.source_pressure_pa - 1.0).abs());
            max_source_energy_error =
                max_source_energy_error.max((original.energy / p.source_energy_j_kg - 1.0).abs());
            cell.energy_offset = e - original.energy;
            cell.pressure_r_offset = (p.pressure_pa - original.pressure) / (rho * p.temperature_k);
            abs_energy_offset += mass * cell.energy_offset.abs();
            max_pressure_correction =
                max_pressure_correction.max((original.pressure / p.pressure_pa - 1.0).abs());
            let eos = ParcelEos::Frozen(cell);
            max_start_pressure_jump =
                max_start_pressure_jump.max((eos.pressure(rho, e) / p.pressure_pa - 1.0).abs());
            max_start_temperature_jump = max_start_temperature_jump
                .max((eos.temperature(rho, e) / p.temperature_k - 1.0).abs());
            materials.push(eos);
        } else {
            materials.push(ParcelEos::Equilibrium(Arc::clone(table)));
        }
        // Also reject a mismatched parent table in the equilibrium restart.
        // Matching only the frozen branch could otherwise hide a wrong input table.
        let material = materials.last().unwrap();
        max_start_pressure_jump =
            max_start_pressure_jump.max((material.pressure(rho, e) / p.pressure_pa - 1.0).abs());
        max_start_temperature_jump = max_start_temperature_jump
            .max((material.temperature(rho, e) / p.temperature_k - 1.0).abs());
    }
    let eos = ParcelEos::Indexed(materials);
    let mut tube = Tube::from_lagrangian(
        input.state.copy(),
        eos.clone(),
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let initial_energy = tube.total_energy();
    let physical_energy =
        (initial_energy - input.storage_shift * tube.total_mass()) * input.area_m2;
    if max_start_pressure_jump > 1e-8
        || max_start_temperature_jump > 1e-8
        || max_source_pressure_error > 1e-3
        || max_source_energy_error > 1e-3
    {
        return Err(format!("splice or cross-language EOS acceptance failed: p_jump={max_start_pressure_jump}, T_jump={max_start_temperature_jump}, source_p={max_source_pressure_error}, source_e={max_source_energy_error}").into());
    }
    let branch = if frozen && molecular {
        "molecular_frozen"
    } else if frozen {
        "parcel_frozen"
    } else {
        "equilibrium"
    };
    let mut time = input.start_time_s;
    let mut impulse = input.initial_wall_impulse_ns;
    let mut max_energy_error: f64 = 0.0;
    let mut steps = 0;
    let mut status = String::from("completed_time_window");
    'outputs: for &target in &input.output_times_s {
        loop {
            if let Some(error) = domain_error(&tube, &eos, table) {
                status = format!("invalid_domain: {error}");
                break 'outputs;
            }
            if time >= target {
                break;
            }
            if steps >= 1_000_000 {
                status = String::from("step_budget_exhausted");
                break 'outputs;
            }
            let step = tube.advance_with_cfl_fraction(target - time, input.timestep_fraction);
            impulse += step.left_wall_impulse * input.area_m2;
            time += step.dt;
            steps += 1;
            max_energy_error =
                max_energy_error.max((tube.total_energy() - initial_energy).abs() * input.area_m2);
        }
        let state = tube.lagrangian_state();
        let min_temp = (0..tube.cells())
            .map(|j| {
                eos.for_cell(j)
                    .temperature(tube.density(j), state.specific_energies[j])
            })
            .fold(f64::INFINITY, f64::min);
        row(
            out,
            &json!({"record":"history","name":input.name,"branch":branch,"time_s":time,
            "wall_impulse_ns":impulse,"wall_pressure_pa":tube.pressure(0),"minimum_temperature_k":min_temp,
            "energy_error_j":(tube.total_energy()-initial_energy)*input.area_m2}),
        )?;
    }
    row(
        out,
        &json!({"record":"summary","name":input.name,"branch":branch,"status":status,"time_s":time,"steps":steps,
        "start_time_s":input.start_time_s,"initial_wall_impulse_ns":input.initial_wall_impulse_ns,"wall_impulse_ns":impulse,
        "momentum_error_ns":impulse-input.incoming_momentum_ns-tube.total_momentum()*input.area_m2,
        "max_energy_error_j":max_energy_error,"restart_physical_energy_j":physical_energy,
        "abs_matching_energy_offsets_j":abs_energy_offset,"max_matching_pressure_correction_fraction":max_pressure_correction,
        "max_start_pressure_jump_fraction":max_start_pressure_jump,"max_start_temperature_jump_fraction":max_start_temperature_jump,
        "max_source_pressure_error_fraction":max_source_pressure_error,"max_source_energy_error_fraction":max_source_energy_error,
        "initial_bond_store_j":bond,"initial_ion_store_j":ion,"initial_compressing_mass_kg":compressing_mass}),
    )?;
    eprintln!(
        "{} {branch}: {status}; J={impulse:.6e}, t={time:.6e}, abs match E={abs_energy_offset:.3e}, max match p={max_pressure_correction:.3e}",
        input.name
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if !(args.len() == 4 || (args.len() == 5 && args[4] == "molecular")) {
        return Err("usage: freeze TABLE.json INPUTS.json OUTPUT.jsonl [molecular]".into());
    }
    let molecular = args.len() == 5;
    let table = Arc::new(TableEos::new(tables::Table::load(&args[1])?));
    let input: Inputs = serde_json::from_reader(File::open(&args[2])?)?;
    let curve = Arc::new(input.residual_curve);
    let mut out = BufWriter::new(File::create(&args[3])?);
    row(
        &mut out,
        &json!({"record":"metadata","input_path":args[2],"table_path":args[1],
        "model":"same-state constrained composition; residual exponent 2; T>=1000 K",
        "molecular_frozen_ionization_active":molecular}),
    )?;
    for restart in &input.restarts {
        for frozen in [false, true] {
            run(restart, &table, &curve, frozen, molecular, &mut out)?;
            out.flush()?;
        }
    }
    Ok(())
}
