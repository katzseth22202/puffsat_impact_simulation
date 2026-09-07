//! Planar prepositioned-water reference, using the existing Lagrangian kernel.
//! x points away from the plate (opposite the study's prograde z). No injection
//! boundary, lateral escape, radiation, conduction, or finite-rate chemistry.

use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

use hydro1d::eos::{Eos, TableEos};
use hydro1d::kernel::{Boundary, LagrangianState, Tube, Viscosity};
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Debug, Deserialize)]
struct Config {
    name: String,
    area_m2: f64,
    projectile_kg: f64,
    water_kg: f64,
    speed_m_s: f64,
    pulse_length_m: f64,
    layer_length_m: f64,
    projectile_energy_j_kg: f64,
    liquid_energy_j_kg: f64,
    energy_shift_j_kg: f64,
    projectile_cells: usize,
    water_cells: usize,
    #[serde(default = "full_cfl")]
    timestep_fraction: f64,
    output_times_s: Vec<f64>,
}

fn full_cfl() -> f64 {
    1.0
}

#[derive(Debug)]
struct Initial {
    state: LagrangianState,
    initial_wall_impulse: f64,
    input_energy: f64,
    incoming_momentum: f64,
}

fn initialize(config: &Config) -> Initial {
    assert!(config.area_m2 > 0.0 && config.speed_m_s > 0.0);
    let mut x = vec![0.0];
    let mut masses = Vec::new();
    let mut energy = Vec::new();
    let mut velocity = Vec::new();
    for (n, mass, length, e, v) in [
        (
            config.water_cells,
            config.water_kg,
            config.layer_length_m,
            config.liquid_energy_j_kg,
            0.0,
        ),
        (
            config.projectile_cells,
            config.projectile_kg,
            config.pulse_length_m,
            config.projectile_energy_j_kg,
            -config.speed_m_s,
        ),
    ] {
        if n == 0 {
            assert_eq!(mass, 0.0);
            continue;
        }
        assert!(mass > 0.0 && length > 0.0);
        for _ in 0..n {
            x.push(x.last().copied().unwrap() + length / n as f64);
            masses.push(mass / (config.area_m2 * n as f64));
            energy.push(e + config.energy_shift_j_kg);
            velocity.push(v);
        }
    }
    let n = masses.len();
    assert!(n > 0);
    let input_energy: f64 = (0..n)
        .map(|j| masses[j] * (energy[j] + 0.5 * velocity[j].powi(2)))
        .sum();
    let incoming_momentum = config.projectile_kg * config.speed_m_s / config.area_m2;
    let mut u = vec![0.0; n + 1];
    for i in 1..n {
        u[i] = (masses[i - 1] * velocity[i - 1] + masses[i] * velocity[i])
            / (masses[i - 1] + masses[i]);
    }
    u[n] = velocity[n - 1];
    // Each cell puts half its mass at each node. Mass-weighted projection
    // preserves momentum internally; convert its KE deficit locally to heat.
    // Pinning the wall node has a separate, explicit initial impulse credit.
    let initial_wall_impulse = -0.5 * masses[0] * velocity[0];
    for j in 0..n {
        energy[j] += 0.25 * ((velocity[j] - u[j]).powi(2) + (velocity[j] - u[j + 1]).powi(2));
    }
    Initial {
        state: LagrangianState {
            positions: x,
            node_velocities: u,
            cell_masses: masses,
            specific_energies: energy,
        },
        initial_wall_impulse,
        input_energy,
        incoming_momentum,
    }
}

#[derive(Debug, Serialize)]
struct Cell {
    x_m: f64,
    width_m: f64,
    mass_kg: f64,
    rho_kg_m3: f64,
    temperature_k: f64,
    energy_j_kg: f64,
    pressure_pa: f64,
    velocity_m_s: f64,
    /// Positive for expansion: -d ln(rho)/dt = (u_right-u_left)/width.
    expansion_rate_s: f64,
}

fn write_row(out: &mut impl Write, value: &serde_json::Value) -> Result<(), Box<dyn Error>> {
    serde_json::to_writer(&mut *out, value)?;
    writeln!(out)?;
    Ok(())
}

fn snapshot(tube: &Tube<TableEos>, eos: &TableEos, config: &Config) -> Vec<Cell> {
    let state = tube.lagrangian_state();
    (0..tube.cells())
        .map(|j| {
            let rho = tube.density(j);
            let e = state.specific_energies[j];
            Cell {
                x_m: tube.center(j),
                width_m: tube.width(j),
                mass_kg: state.cell_masses[j] * config.area_m2,
                rho_kg_m3: rho,
                temperature_k: eos.temperature(rho, e),
                energy_j_kg: e - config.energy_shift_j_kg,
                pressure_pa: tube.pressure(j),
                velocity_m_s: tube.velocity(j),
                expansion_rate_s: (state.node_velocities[j + 1] - state.node_velocities[j])
                    / tube.width(j),
            }
        })
        .collect()
}

fn domain_error(tube: &Tube<TableEos>, eos: &TableEos) -> Option<String> {
    let table = eos.table();
    let state = tube.lagrangian_state();
    let (rlo, rhi) = (table.rho_grid()[0], *table.rho_grid().last().unwrap());
    let (tlo, thi) = (table.t_grid()[0], *table.t_grid().last().unwrap());
    for (j, &energy) in state.specific_energies.iter().enumerate() {
        let rho = tube.density(j);
        if !rho.is_finite() || rho < rlo || rho > rhi {
            return Some(format!("cell {j}: density {rho} outside [{rlo},{rhi}]"));
        }
        let (elo, ehi) = (table.energy(rho, tlo), table.energy(rho, thi));
        if !energy.is_finite() || energy < elo * (1.0 - 1e-9) || energy > ehi * (1.0 + 1e-9) {
            return Some(format!("cell {j}: energy {energy} outside [{elo},{ehi}]"));
        }
    }
    None
}

/// A uniform stationary layer cannot change its wall pressure before the head
/// of the vacuum rarefaction arrives. This needs only the undisturbed state,
/// not a metastable-liquid continuation of the cold free surface into ice.
fn undisturbed_wall_window(
    eos: &impl Eos,
    rho: f64,
    energy: f64,
    length: f64,
    end_time: f64,
) -> Option<(f64, f64)> {
    let arrival = length / eos.sound_speed(rho, energy);
    (end_time < arrival).then(|| (eos.pressure(rho, energy), arrival))
}

fn run(config: &Config, eos: &TableEos, out: &mut impl Write) -> Result<(), Box<dyn Error>> {
    assert!(!config.output_times_s.is_empty());
    assert!(config.output_times_s[0] == 0.0);
    assert!(config.output_times_s.windows(2).all(|t| t[0] < t[1]));
    let initial = initialize(config);
    let mut tube = Tube::from_lagrangian(
        initial.state,
        eos.clone(),
        Boundary::Wall,
        Boundary::Free,
        Viscosity::VON_NEUMANN_RICHTMYER,
    );
    let end_time = *config.output_times_s.last().unwrap();
    if config.projectile_kg == 0.0 && domain_error(&tube, eos).is_none() {
        if let Some((pressure, arrival)) = undisturbed_wall_window(
            eos,
            tube.density(0),
            config.liquid_energy_j_kg + config.energy_shift_j_kg,
            config.layer_length_m,
            end_time,
        ) {
            write_row(
                out,
                &json!({
                    "record": "summary", "name": config.name,
                    "status": "completed_analytic_wall_window", "time_s": end_time, "steps": 0,
                    "wall_impulse_ns": pressure * config.area_m2 * end_time,
                    "initial_projection_impulse_ns": 0.0,
                    "initial_rarefaction_arrival_s": arrival,
                    "control_method": "undisturbed wall until rarefaction head arrival; no spatial history",
                    "area_m2": config.area_m2, "layer_length_m": config.layer_length_m,
                    "projectile_cells": 0, "water_cells": config.water_cells,
                    "timestep_fraction": config.timestep_fraction,
                    "peak_wall_pressure_pa": pressure, "final_wall_pressure_pa": pressure,
                    "input_energy_j": config.water_kg * config.liquid_energy_j_kg,
                }),
            )?;
            eprintln!(
                "{}: completed_analytic_wall_window; arrival={arrival:.5e} s, J={:.6e} N s",
                config.name,
                pressure * config.area_m2 * end_time
            );
            return Ok(());
        }
    }
    let mut impulse = initial.initial_wall_impulse;
    let mut time = 0.0;
    let mut steps = 0;
    let mut peak_pressure: f64 = 0.0;
    let mut max_energy_error: f64 = 0.0;
    let mut status = String::from("completed_time_window");
    'outputs: for &target in &config.output_times_s {
        loop {
            if let Some(error) = domain_error(&tube, eos) {
                status = format!("invalid_table_domain: {error}");
                break 'outputs;
            }
            if time >= target {
                break;
            }
            if steps >= 1_000_000 {
                status = String::from("step_budget_exhausted");
                break 'outputs;
            }
            let step = tube.advance_with_cfl_fraction(target - time, config.timestep_fraction);
            if step.dt < 1e-14 {
                status = String::from("timestep_stalled");
                break 'outputs;
            }
            impulse += step.left_wall_impulse;
            time += step.dt;
            steps += 1;
            peak_pressure = peak_pressure.max(tube.pressure(0));
            max_energy_error =
                max_energy_error.max((tube.total_energy() - initial.input_energy).abs());
        }
        write_row(
            out,
            &json!({
                "record": "snapshot", "name": config.name, "time_s": time,
                "wall_impulse_ns": impulse * config.area_m2,
                "gas_momentum_x_ns": tube.total_momentum() * config.area_m2,
                "energy_error_j": (tube.total_energy() - initial.input_energy) * config.area_m2,
                "water_cells": config.water_cells, "cells": snapshot(&tube, eos, config),
                "lagrangian_state": {
                    "positions": tube.lagrangian_state().positions,
                    "node_velocities": tube.lagrangian_state().node_velocities,
                    "cell_masses": tube.lagrangian_state().cell_masses,
                    "specific_energies": tube.lagrangian_state().specific_energies,
                },
            }),
        )?;
    }
    write_row(
        out,
        &json!({
            "record": "summary", "name": config.name, "status": status,
            "time_s": time, "steps": steps,
            "projectile_cells": config.projectile_cells, "water_cells": config.water_cells,
            "timestep_fraction": config.timestep_fraction,
            "layer_length_m": config.layer_length_m, "area_m2": config.area_m2,
            "wall_impulse_ns": impulse * config.area_m2,
            "initial_projection_impulse_ns": initial.initial_wall_impulse * config.area_m2,
            "momentum_error_ns": (impulse - initial.incoming_momentum - tube.total_momentum()) * config.area_m2,
            "max_energy_error_j": max_energy_error * config.area_m2,
            "input_energy_j": (initial.input_energy - config.energy_shift_j_kg * tube.total_mass()) * config.area_m2,
            "peak_wall_pressure_pa": peak_pressure,
            "final_wall_pressure_pa": tube.pressure(0),
        }),
    )?;
    eprintln!(
        "{}: {status}; t={time:.5e}, steps={steps}, J={:.6e} N s",
        config.name,
        impulse * config.area_m2
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: water_plate TABLE.json CONFIGS.json OUTPUT.jsonl".into());
    }
    let eos = TableEos::new(tables::Table::load(&args[1])?);
    let configs: Vec<Config> = serde_json::from_reader(File::open(&args[2])?)?;
    let mut out = BufWriter::new(File::create(&args[3])?);
    let raw_table: serde_json::Value = serde_json::from_reader(File::open(&args[1])?)?;
    let raw_configs: serde_json::Value = serde_json::from_reader(File::open(&args[2])?)?;
    write_row(
        &mut out,
        &json!({
            "record": "metadata", "table_path": args[1], "config_path": args[2],
            "table_provenance": raw_table["provenance"], "configs": raw_configs,
        }),
    )?;
    for config in &configs {
        run(config, &eos, &mut out)?;
        out.flush()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use hydro1d::eos::IdealGas;

    #[test]
    fn cold_control_requires_the_rarefaction_not_to_have_arrived() {
        let eos = IdealGas::new(1.4);
        let arrival = 1.0 / eos.sound_speed(5.0, 2e6);
        let (p, measured_arrival) =
            undisturbed_wall_window(&eos, 5.0, 2e6, 1.0, arrival / 2.0).unwrap();
        assert_relative_eq!(p, 4e6, max_relative = 1e-12);
        assert_relative_eq!(measured_arrival, arrival, max_relative = 1e-12);
        assert!(undisturbed_wall_window(&eos, 5.0, 2e6, 1.0, arrival).is_none());
    }

    #[test]
    fn heterogeneous_initial_projection_conserves_energy_and_accounts_for_wall_impulse() {
        for (water, projectile) in [(212.5, 25.0), (0.0, 25.0), (212.5, 0.0)] {
            let config = Config {
                name: String::new(),
                area_m2: 40.0,
                water_kg: water,
                projectile_kg: projectile,
                speed_m_s: 45580.0,
                layer_length_m: 0.1,
                pulse_length_m: 4.558,
                liquid_energy_j_kg: -1.914e6,
                projectile_energy_j_kg: 0.56e6,
                energy_shift_j_kg: 3e6,
                water_cells: if water > 0.0 { 170 } else { 0 },
                timestep_fraction: 1.0,
                projectile_cells: if projectile > 0.0 { 20 } else { 0 },
                output_times_s: vec![0.0],
            };
            let initial = initialize(&config);
            let tube = Tube::from_lagrangian(
                initial.state,
                IdealGas::new(1.4),
                Boundary::Wall,
                Boundary::Free,
                Viscosity::VON_NEUMANN_RICHTMYER,
            );
            assert_relative_eq!(
                tube.total_energy(),
                initial.input_energy,
                max_relative = 1e-12
            );
            assert_relative_eq!(
                tube.total_momentum() + initial.incoming_momentum,
                initial.initial_wall_impulse,
                epsilon = 1e-9
            );
            assert_relative_eq!(
                tube.total_mass() * config.area_m2,
                water + projectile,
                max_relative = 1e-12
            );
        }
    }
}
