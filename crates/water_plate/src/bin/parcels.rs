//! Conditional finite-rate chemistry along prescribed, not recomputed, flow.

use serde::Deserialize;
use serde_json::json;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use water_plate::ResidualCurve;
use water_plate::parcel::{Chemistry, Context, Knot, Settings, State, integrate};

#[derive(Debug, Deserialize)]
struct Parcel {
    index: usize,
    mass_kg: f64,
    initial_y: [f64; 6],
    initial_temperature: f64,
    initial_energy: f64,
    initial_pressure: f64,
    initial_bond: f64,
    spectators: f64,
    history: Vec<Knot>,
    parent_returns: serde_json::Value,
    initial_ionized_nuclei_fraction: f64,
}

#[derive(Debug, Deserialize)]
struct Variant {
    name: String,
    chemistry: Chemistry,
}

#[derive(Debug, Deserialize)]
struct Input {
    residual_curve: ResidualCurve,
    parcels: Vec<Parcel>,
    settings: Settings,
    variants: Vec<Variant>,
    forcing_modes: Vec<bool>,
    report_times: Vec<f64>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: parcels INPUT.json OUTPUT.jsonl".into());
    }
    let input: Input = serde_json::from_reader(BufReader::new(File::open(&args[1])?))?;
    let mut out = BufWriter::new(File::create(&args[2])?);
    serde_json::to_writer(
        &mut out,
        &json!({"record":"metadata","input_path":args[1],
        "settings":input.settings,"kinetics_validated":false}),
    )?;
    writeln!(out)?;
    for variant in &input.variants {
        for &forcing in &input.forcing_modes {
            let mut completed = 0;
            for parcel in &input.parcels {
                let y = parcel.initial_y;
                let state = State {
                    molecules: [y[2], y[3], y[4], y[5]],
                    temperature: parcel.initial_temperature,
                };
                let mut c = Context {
                    chem: &variant.chemistry,
                    residual: &input.residual_curve,
                    atoms: [
                        y[0] + 2.0 * y[2] + y[3] + 2.0 * y[5],
                        y[1] + y[3] + 2.0 * y[4] + y[5],
                    ],
                    spectators: parcel.spectators,
                    energy_offset: 0.0,
                    pressure_r_offset: 0.0,
                };
                let rho = parcel.history[0].rho;
                let (energy, pressure, _) = c.thermo(state, rho)?;
                c.energy_offset = parcel.initial_energy - energy;
                c.pressure_r_offset =
                    (parcel.initial_pressure - pressure) / (rho * state.temperature);
                let settings = Settings {
                    parent_forcing: forcing,
                    ..input.settings
                };
                let mut result = integrate(&c, state, &parcel.history, settings);
                completed += usize::from(result.status == "completed");
                let max_energy_error = result
                    .samples
                    .iter()
                    .map(|s| s.energy_error_j_kg.abs())
                    .fold(0.0_f64, f64::max);
                let min_temperature = result
                    .samples
                    .iter()
                    .map(|s| s.state.temperature)
                    .fold(f64::INFINITY, f64::min);
                result.samples.retain(|s| {
                    input
                        .report_times
                        .iter()
                        .any(|t| (s.time_s - t).abs() < 1e-12)
                });
                serde_json::to_writer(
                    &mut out,
                    &json!({"record":"parcel","variant":variant.name,
                    "parent_forcing":forcing,"index":parcel.index,"mass_kg":parcel.mass_kg,
                    "initial_bond_j_kg":parcel.initial_bond,"initial_energy_j_kg":parcel.initial_energy,
                    "initial_ionized_nuclei_fraction":parcel.initial_ionized_nuclei_fraction,
                    "energy_offset_j_kg":c.energy_offset,"pressure_r_offset":c.pressure_r_offset,
                    "max_energy_error_j_kg":max_energy_error,"minimum_sampled_temperature_k":min_temperature,
                    "parent_returns":parcel.parent_returns,"trajectory":result}),
                )?;
                writeln!(out)?;
            }
            out.flush()?;
            eprintln!(
                "{} parent_forcing={forcing}: {completed}/{} complete",
                variant.name,
                input.parcels.len()
            );
        }
    }
    Ok(())
}
