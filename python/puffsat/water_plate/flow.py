"""Inputs and offline analysis for the planar prepositioned-water reference.

This step puts cold water in a layer touching a rigid wall, then impacts it with
the specified pulse. It is a lossless 1D reference, not the distributed injector
or a finite-plate performance calculation. The liquid energy source is explicit;
there is no prescribed initial mixed-gas density or temperature.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import cast

import numpy as np

from puffsat import eos_water, recombination
from puffsat.water_plate import flow_eos
from puffsat.water_plate.profiles import IncomingPulse
from puffsat.water_plate.thermodynamics import OUTPUT_DIR, energy_reference


def write_inputs(
    cells: int,
    layer: float,
    samples: int = 64,
    timestep_fraction: float = 0.25,
    coupling_samples: int = 0,
) -> Path:
    """Matched combined/unsprayed/water-only controls at k=8.5 and 45.58 km/s."""
    if cells <= 0 or cells % 2 or not math.isfinite(layer) or layer <= 0 or samples < 2:
        raise ValueError("need a positive even projectile cell count and positive layer depth")
    if not 0 < timestep_fraction <= 1:
        raise ValueError("timestep fraction must be in (0,1]")
    if coupling_samples < 0 or coupling_samples == 1:
        raise ValueError("coupling samples must be zero or at least 2")
    pulse = IncomingPulse(3.75, 100e-6)
    area = math.pi * pulse.radius_m**2
    times = np.unique(
        np.concatenate(
            (
                [0.0],
                np.geomspace(1e-7, 0.001, samples),
                np.linspace(100e-6, 250e-6, coupling_samples),
            )
        )
    ).tolist()
    initial_energy = flow_eos.pressure_energy(pulse.density_kg_m3, 400.0)[1] - flow_eos.ENERGY_SHIFT
    configs = []
    for name, projectile, water in (
        ("combined", 25.0, 212.5),
        ("unsprayed", 25.0, 0.0),
        ("water_only", 0.0, 212.5),
    ):
        configs.append(
            {
                "name": f"{name}_h{layer:g}_n{cells}",
                "area_m2": area,
                "projectile_kg": projectile,
                "water_kg": water,
                "speed_m_s": pulse.speed_m_s,
                "pulse_length_m": pulse.length_m,
                "layer_length_m": layer,
                "projectile_energy_j_kg": initial_energy,
                "liquid_energy_j_kg": energy_reference().liquid_energy_j_kg,
                "energy_shift_j_kg": flow_eos.ENERGY_SHIFT,
                "projectile_cells": cells if projectile else 0,
                "water_cells": 17 * cells // 2 if water else 0,
                "timestep_fraction": timestep_fraction,
                "output_times_s": times,
            }
        )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"flow_inputs_h{layer:g}_n{cells}.json"
    path.write_text(json.dumps(configs) + "\n")
    print(f"Wrote {path}; initial cold-water bulk density={212.5 / (area * layer):.6g} kg/m3")
    return path


CsvValue = str | int | float | bool


def _write(path: Path, rows: list[dict[str, CsvValue]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@dataclass(frozen=True)
class FlowCell:
    """Typed boundary for a Rust-generated material cell snapshot."""

    x_m: float
    width_m: float
    mass_kg: float
    rho_kg_m3: float
    temperature_k: float
    energy_j_kg: float
    pressure_pa: float
    velocity_m_s: float
    # NaN marks older files without this measurement; never substitute a clock.
    expansion_rate_s: float = math.nan


@dataclass(frozen=True)
class Snapshot:
    """A fixed-mass Lagrangian history at a sampled time."""

    name: str
    time_s: float
    impulse_ns: float
    water_cells: int
    cells: list[FlowCell]


def read_output(path: Path) -> tuple[list[Snapshot], list[dict[str, CsvValue]]]:
    """Decode the locally generated JSONL at a typed boundary."""
    snapshots = []
    summaries = []
    with path.open() as stream:
        for line in stream:
            raw = json.loads(line)
            if raw["record"] == "metadata":
                continue
            if raw["record"] == "summary":
                summaries.append(cast(dict[str, CsvValue], raw))
            elif raw["record"] == "snapshot":
                snapshots.append(
                    Snapshot(
                        str(raw["name"]),
                        float(raw["time_s"]),
                        float(raw["wall_impulse_ns"]),
                        int(raw["water_cells"]),
                        [FlowCell(**c) for c in raw["cells"]],
                    )
                )
            else:
                raise ValueError(f"unknown flow record: {raw['record']}")
    return snapshots, summaries


@lru_cache(maxsize=32768)
def cell_stores(cell: FlowCell) -> tuple[float, float]:
    """Gas-EOS bond/ion proxies. Below 1000 K water is taken as molecular.

    Above 1000 K these decompose the ideal-gas part of the extended potential,
    not the dense-fluid residual. Dense/hot closure mass is reported separately.
    """
    rho, temp = cell.rho_kg_m3, cell.temperature_k
    if temp <= flow_eos.COLD_MAX_T:
        return 0.0, 0.0
    comp = eos_water.composition(rho, temp)
    bond = eos_water.bond_energy_held(rho, temp)
    ion = (
        comp.n_hp * eos_water.IP_H
        + sum(n * energy for n, energy in zip(comp.n_o_ions, eos_water.E_O_CUM, strict=True))
    ) / rho
    return bond, ion


def history_rows(snapshots: list[Snapshot]) -> list[dict[str, CsvValue]]:
    """Report stores and drop from each parcel's sampled prior peak versus impulse.

    This does not identify a causal chemistry contribution to wall impulse.
    Sampled peaks need temporal convergence; they are not assumed exact maxima.
    """
    peaks: dict[str, list[float]] = {}
    ion_peaks: dict[str, list[float]] = {}
    previous: dict[str, Snapshot] = {}
    rows: list[dict[str, CsvValue]] = []
    for snap in snapshots:
        prior = peaks.setdefault(snap.name, [0.0] * len(snap.cells))
        prior_ion = ion_peaks.setdefault(snap.name, [0.0] * len(snap.cells))
        bond_total = ion_total = drop = water_peak = dense_hot = closure_mass = 0.0
        ion_drop = water_bond = water_ion = expanding_bond = expanding_ion = 0.0
        slow_bond = slow_ion = fast_bond = fast_ion = 0.0
        max_water_rho = max_water_temp = 0.0
        high_pressure_mass = interval_bond = interval_ion = bond_x = ion_x = 0.0
        before = previous.get(snap.name)
        for j, cell in enumerate(snap.cells):
            bond, ion = cell_stores(cell)
            prior[j] = max(prior[j], bond)
            prior_ion[j] = max(prior_ion[j], ion)
            bond_total += cell.mass_kg * bond
            ion_total += cell.mass_kg * ion
            drop += cell.mass_kg * (prior[j] - bond)
            ion_drop += cell.mass_kg * (prior_ion[j] - ion)
            if before is not None:
                old = before.cells[j]
                old_bond, old_ion = cell_stores(old)
                db = cell.mass_kg * max(0.0, old_bond - bond)
                di = cell.mass_kg * max(0.0, old_ion - ion)
                interval_bond += db
                interval_ion += di
                # Midpoint location is a sampling approximation, not a resolved reaction front.
                bond_x += db * (old.x_m + cell.x_m) / 2
                ion_x += di * (old.x_m + cell.x_m) / 2
            clocks = reaction_damkohler(cell)
            if clocks is not None:
                da_bond, da_ion = clocks
                expanding_bond += cell.mass_kg * bond
                expanding_ion += cell.mass_kg * ion
                slow_bond += cell.mass_kg * bond * (da_bond < 0.1)
                slow_ion += cell.mass_kg * ion * (da_ion < 0.1)
                fast_bond += cell.mass_kg * bond * (da_bond > 10)
                fast_ion += cell.mass_kg * ion * (da_ion > 10)
            if cell.temperature_k >= 5000 and cell.rho_kg_m3 >= 10:
                dense_hot += cell.mass_kg
            if cell.temperature_k > flow_eos.COLD_MAX_T and cell.rho_kg_m3 >= 100:
                closure_mass += cell.mass_kg
            if cell.pressure_pa > 1e9:
                high_pressure_mass += cell.mass_kg
            if j < snap.water_cells:
                water_bond += cell.mass_kg * bond
                water_ion += cell.mass_kg * ion
                water_peak += cell.mass_kg * prior[j]
                max_water_rho = max(max_water_rho, cell.rho_kg_m3)
                max_water_temp = max(max_water_temp, cell.temperature_k)
        rows.append(
            {
                "name": snap.name,
                "time_s": snap.time_s,
                "wall_impulse_ns": snap.impulse_ns,
                "bond_store_j": bond_total,
                "ion_store_j": ion_total,
                "bond_drop_from_sampled_parcel_peaks_j": drop,
                "ion_drop_from_sampled_parcel_peaks_j": ion_drop,
                "bond_interval_store_decrease_j": interval_bond,
                "ion_interval_store_decrease_j": interval_ion,
                "bond_decrease_midpoint_x_m": bond_x / interval_bond if interval_bond > 0 else "",
                "ion_decrease_midpoint_x_m": ion_x / interval_ion if interval_ion > 0 else "",
                "water_bond_store_j": water_bond,
                "water_ion_store_j": water_ion,
                "water_sampled_peak_bond_store_j": water_peak,
                "expanding_gas_bond_store_j": expanding_bond,
                "expanding_gas_ion_store_j": expanding_ion,
                "expanding_bond_store_da_oh_lt0p1_j": slow_bond,
                "expanding_bond_store_da_oh_gt10_j": fast_bond,
                "expanding_ion_store_da_coll_lt0p1_j": slow_ion,
                "expanding_ion_store_da_coll_gt10_j": fast_ion,
                "hot_dense_mass_kg_Tge5000_rhoge10": dense_hot,
                "dense_closure_mass_kg_Tgt1000_rhoge100": closure_mass,
                "high_pressure_mass_kg_pgt1GPa": high_pressure_mass,
                "max_water_density_kg_m3": max_water_rho,
                "max_water_temperature_k": max_water_temp,
            }
        )
        previous[snap.name] = snap
    return rows


def reaction_damkohler(cell: FlowCell) -> tuple[float, float] | None:
    """OH-partner molecular and collisional-ion proxies on the actual local clock.

    Only expanding gas above the real-fluid transition is screened. Missing
    velocity gradients, compression and stationary cells have no expansion Da.
    These remain rate proxies, not a finite-rate chemical network or heat credit.
    """
    rate = cell.expansion_rate_s
    if not math.isfinite(rate) or rate <= 0 or cell.temperature_k < 1800:
        return None
    comp = eos_water.composition(cell.rho_kg_m3, cell.temperature_k)
    tau_bond = recombination.atom_recombination_time(
        cell.temperature_k, comp.n_neutral_heavy, comp.n_oh
    )
    ion_rate = recombination.three_body_coefficient(cell.temperature_k) * comp.n_e**2
    return 1 / (rate * tau_bond), ion_rate / rate


def analyze(path: Path) -> None:
    """Keep invalid-domain results explicit; calculate stores on actual cell histories."""
    snapshots, summaries = read_output(path)
    _write(path.with_name(path.stem + "_summary.csv"), summaries)
    rows = history_rows(snapshots)
    _write(path.with_name(path.stem + "_history.csv"), rows)
    gain = matched_gain(summaries)
    if gain is not None:
        print(f"Provisional planar finite-window gain G = {gain:.8g} N s")
    else:
        print("Invalid, unconserved, unmatched or incomplete control: no matched gain reported.")
    for row in summaries:
        print(row)


def matched_gain(summaries: list[dict[str, CsvValue]]) -> float | None:
    """Require matched time/geometry, valid domains and <=0.5% energy drift.

    This is a diagnostic gate, not a claim of resolution convergence. The input
    energy's magnitude normalizes the cold-only control (its energy zero is
    negative); collision controls are dominated by the incoming kinetic energy.
    """
    by_case = {str(row["name"]).split("_h")[0]: row for row in summaries}
    if len(summaries) != 3 or set(by_case) != {"combined", "unsprayed", "water_only"}:
        return None
    for row in summaries:
        if (
            str(row["name"]).startswith("water_only_h")
            and row["status"] == "completed_analytic_wall_window"
        ):
            if float(row["time_s"]) >= float(row["initial_rarefaction_arrival_s"]):
                return None
            continue
        if row["status"] != "completed_time_window":
            return None
        energy_scale = abs(float(row["input_energy_j"]))
        if float(row["max_energy_error_j"]) > 0.005 * energy_scale:
            return None
        if abs(float(row["momentum_error_ns"])) > 1e-8 * max(
            1.0, abs(float(row["wall_impulse_ns"]))
        ):
            return None
    for field in ("time_s", "area_m2", "layer_length_m"):
        if any(
            not math.isclose(float(row[field]), float(summaries[0][field]), rel_tol=1e-12)
            for row in summaries[1:]
        ):
            return None
    return (
        float(by_case["combined"]["wall_impulse_ns"])
        - float(by_case["unsprayed"]["wall_impulse_ns"])
        - float(by_case["water_only"]["wall_impulse_ns"])
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=20)
    parser.add_argument("--layer", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--timestep-fraction", type=float, default=0.25)
    parser.add_argument(
        "--coupling-samples",
        type=int,
        default=0,
        help="extra samples across 100-250 microseconds at the 1 m anchor",
    )
    parser.add_argument("--analyze", type=Path)
    args = parser.parse_args()
    if args.analyze:
        analyze(args.analyze)
    else:
        write_inputs(
            args.cells, args.layer, args.samples, args.timestep_fraction, args.coupling_samples
        )


if __name__ == "__main__":
    main()
