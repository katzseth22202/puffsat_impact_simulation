"""Prepare energy-preserving parcel-frozen restarts of the stratified reference.

All chemical reactions are disabled after the chosen switch; this is a
counterfactual closure test, not a physical freeze-out model or an additive
measurement of molecular recombination work. The equilibrium control restarts
the identical staggered state. Dense-fluid physics remains provisional.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate import flow_eos


@dataclass(frozen=True)
class Vibration:
    theta_k: float
    amplitude_j_kg: float


@dataclass(frozen=True)
class FrozenGas:
    gas_constant: float
    linear_cv: float
    chemical_energy: float
    vibrations: list[Vibration]
    bond_store: float
    ion_store: float

    def thermal(self, temp: float) -> float:
        return (
            self.chemical_energy
            + self.linear_cv * temp
            + sum(v.amplitude_j_kg / math.expm1(v.theta_k / temp) for v in self.vibrations)
        )


def gas_parameters(rho: float, temp: float) -> FrozenGas:
    y = ew.frozen_composition(rho, temp)
    factor = ew.K_B / ew.M_H2O
    mono = y.y_h + y.y_o + y.y_hp + sum(y.y_o_ions) + y.y_e
    diatomic = sum(y.y_intermediates)
    modes = [Vibration(t, factor * y.y_h2o * t) for t in ew._THETA_VIB]
    modes.extend(
        Vibration(d.theta_vib, factor * fraction * d.theta_vib)
        for d, fraction in zip(ew.INTERMEDIATES, y.y_intermediates, strict=True)
    )
    bond = ew.bond_energy_held(rho, temp)
    ion = (
        y.y_hp * ew.IP_H + sum(a * b for a, b in zip(y.y_o_ions, ew.E_O_CUM, strict=True))
    ) / ew.M_H2O
    return FrozenGas(
        factor * (mono + y.y_h2o + diatomic),
        factor * (1.5 * mono + 3 * y.y_h2o + 2.5 * diatomic),
        bond + ion,
        modes,
        bond,
        ion,
    )


def residual_curve(points: int = 512) -> dict[str, object]:
    """Hermite data for one density-dependent potential; no separate p/e interpolation."""
    knots = []
    for rho in np.geomspace(1e-5, 2000, points):
        r = float(rho)
        ref = flow_eos.residual_reference(r)
        knots.append(
            {
                "log_rho": math.log(r),
                "values": [ref.energy, ref.entropy, ref.cv],
                "slopes": [
                    (ref.pressure - 1000 * ref.pressure_t) / r,
                    -ref.pressure_t / r,
                    -1000 * ref.pressure_tt / r,
                ],
            }
        )
    return {"knots": knots}


def sample_times(start: float, end: float) -> list[float]:
    """Resolve the loading tail without introducing samples before the restart."""
    if not (math.isfinite(start) and math.isfinite(end) and 0 < start < end):
        raise ValueError("require finite 0 < switch < end")
    transition_end = max(start, min(end, 250e-6))
    early = np.linspace(start, transition_end, 101)
    late = np.geomspace(transition_end, end, 40) if transition_end < end else np.array([end])
    return [float(t) for t in np.unique(np.clip(np.concatenate((early, late)), start, end))]


def prepare(parent: Path, times: list[float], output: Path, end: float = 0.001) -> None:
    # Any stays confined to this locally produced JSON boundary. Calculations
    # receive explicit floats or typed gas parameters, never an untyped EOS API.
    selected: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {}
    with parent.open() as stream:
        for line in stream:
            raw = json.loads(line)
            if raw["record"] == "metadata":
                metadata = raw
            elif (
                raw["record"] == "snapshot"
                and not raw["name"].startswith("water_only")
                and any(
                    math.isclose(float(raw["time_s"]), t, rel_tol=0, abs_tol=1e-12) for t in times
                )
            ):
                selected.append(raw)
    if len(selected) != 2 * len(times):
        raise ValueError("need exact combined and unsprayed checkpoints at every switch time")
    provenance = metadata["table_provenance"]
    if provenance["dense_hot_cv_residual_power"] != 2:
        raise ValueError("the frozen continuation currently implements residual exponent 2 only")
    shift = float(provenance["energy_storage_shift_j_kg"])
    configs = {str(c["name"]): c for c in metadata["configs"]}
    restarts = []
    for snap in selected:
        if "lagrangian_state" not in snap:
            raise ValueError(
                "checkpoint must include exact staggered arrays, not reconstructed nodes"
            )
        start = float(snap["time_s"])
        if not start < end:
            raise ValueError("restart end must be later than switch time")
        states = []
        for cell in snap["cells"]:
            rho, temp = float(cell["rho_kg_m3"]), float(cell["temperature_k"])
            if temp < 1000:
                raise ValueError(
                    "all selected parcels must be above the dense-fluid matching temperature"
                )
            gas = gas_parameters(rho, temp)
            p, e = flow_eos.pressure_energy(rho, temp)
            states.append(
                {
                    "gas": asdict(gas),
                    "temperature_k": temp,
                    "pressure_pa": float(cell["pressure_pa"]),
                    "source_pressure_pa": p,
                    "source_energy_j_kg": e,
                    "expansion_rate_s": float(cell["expansion_rate_s"]),
                }
            )
        samples = sample_times(start, end)
        c = configs[snap["name"]]
        restarts.append(
            {
                "name": f"{snap['name']}_freeze_{start * 1e6:g}us",
                "start_time_s": start,
                "initial_wall_impulse_ns": float(snap["wall_impulse_ns"]),
                "area_m2": float(c["area_m2"]),
                "storage_shift": shift,
                "water_cells": int(snap["water_cells"]),
                "incoming_momentum_ns": float(c["projectile_kg"]) * float(c["speed_m_s"]),
                "timestep_fraction": float(c["timestep_fraction"]),
                "state": snap["lagrangian_state"],
                "parcels": states,
                "output_times_s": samples,
            }
        )
    payload = {
        "parent_path": str(parent),
        "parent_metadata": metadata,
        "residual_curve": residual_curve(),
        "restarts": restarts,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload) + "\n")
    print(f"Wrote {output}: {len(restarts)} matched restarts")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path)
    parser.add_argument("--times-us", nargs="+", type=float, default=[153.5, 193.0])
    parser.add_argument("--end", type=float, default=0.001)
    parser.add_argument(
        "--output", type=Path, default=Path("data/results/water_plate/freeze_inputs.json")
    )
    args = parser.parse_args()
    prepare(args.parent, [t * 1e-6 for t in args.times_us], args.output, args.end)


if __name__ == "__main__":
    main()
