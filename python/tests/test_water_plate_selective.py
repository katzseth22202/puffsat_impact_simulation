"""Constrained chemistry must preserve inventories, its splice and thermodynamics."""

import json
import math
from pathlib import Path

import pytest

from puffsat import eos_water as ew
from puffsat.water_plate.freeze import gas_parameters
from puffsat.water_plate.selective import IonizingGas, ionizing_parameters
from puffsat.water_plate.selective_report import matched_inputs


@pytest.mark.parametrize(
    "rho,temp",
    [
        (0.01, 1e5),
        (1.0, 16000.0),
        (300.0, 8000.0),
        (1e-5, 1000.0),
        (2000.0, 1000.0),
        (1e-5, 2e6),
        (2000.0, 2e6),
    ],
)
def test_selective_equilibrium_matches_parent_at_switch(rho: float, temp: float) -> None:
    gas = ionizing_parameters(rho, temp)
    assert gas.pressure_energy(rho, temp) == pytest.approx(ew.pressure_energy(rho, temp), rel=1e-9)
    y = ew.frozen_composition(rho, temp)
    rows = gas.populations(rho, temp)
    assert rows[0] == pytest.approx([y.y_h, y.y_hp], rel=1e-8, abs=1e-14)
    assert rows[1] == pytest.approx([y.y_o, *y.y_o_ions], rel=1e-8, abs=1e-14)


@pytest.mark.parametrize(
    "rho,temp", [(1e-5, 1000.0), (2000.0, 3000.0), (0.1, 16000.0), (10.0, 2e6)]
)
def test_selective_inventory_charge_and_mass_action(rho: float, temp: float) -> None:
    gas = ionizing_parameters(1.0, 16000.0)
    rows = gas.populations(rho, temp)
    electrons = sum(z * p for row in rows for z, p in enumerate(row))
    for a, row in zip(gas.ladders, rows, strict=True):
        assert sum(row) == pytest.approx(a.abundance, rel=1e-13)
        for z in range(1, len(row)):
            if row[z] > 1e-200 and row[z - 1] > 1e-200:
                lhs = math.log(row[z] / row[z - 1]) + math.log(electrons * rho / ew.M_H2O)
                rhs = a.log_prefactors[z] - a.log_prefactors[z - 1] + 1.5 * math.log(temp)
                rhs -= (a.energies[z] - a.energies[z - 1]) / (ew.K_B * temp)
                assert lhs == pytest.approx(rhs, abs=1e-10)
    assert gas.base.chemical_energy == gas.base.bond_store


@pytest.mark.parametrize("rho,temp", [(0.01, 5000.0), (1.0, 16000.0), (100.0, 1e5)])
def test_selective_maxwell_and_positive_cv(rho: float, temp: float) -> None:
    gas = ionizing_parameters(1.0, 16000.0)
    dr, dt = rho * 1e-4, temp * 1e-4
    p, _ = gas.pressure_energy(rho, temp)
    pr, er = (
        (b - a) / (2 * dr)
        for a, b in zip(
            gas.pressure_energy(rho - dr, temp), gas.pressure_energy(rho + dr, temp), strict=True
        )
    )
    pt, cv = (
        (b - a) / (2 * dt)
        for a, b in zip(
            gas.pressure_energy(rho, temp - dt), gas.pressure_energy(rho, temp + dt), strict=True
        )
    )
    assert abs(rho * rho * er - p + temp * pt) / p < 1e-6
    assert cv > 0
    assert pr + pt * (p / rho**2 - er) / cv > 0


def test_no_available_atomic_inventory_reduces_to_fixed_composition() -> None:
    base = gas_parameters(1.0, 1000.0)
    gas = IonizingGas(base, [])
    assert gas.populations(1.0, 5000.0) == []
    assert gas.pressure_energy(1.0, 5000.0) == (base.gas_constant * 5000, base.thermal(5000))


@pytest.mark.parametrize("rho,temp", [(0.0, 1000.0), (1.0, 999.0), (1.0, math.inf)])
def test_invalid_selective_state_is_rejected(rho: float, temp: float) -> None:
    gas = ionizing_parameters(1.0, 16000.0)
    with pytest.raises(ValueError, match="requires"):
        gas.pressure_energy(rho, temp)


def test_ordered_report_checks_actual_restart_arrays(tmp_path: Path) -> None:
    payload = {
        "residual_curve": {"knots": [1, 2]},
        "restarts": [
            {
                "name": "test",
                "state": {"positions": [0, 1]},
                "parcels": [{"gas": {"gas_constant": 1.0}}],
                "start_time_s": 1,
                "initial_wall_impulse_ns": 2,
                "area_m2": 3,
                "storage_shift": 4,
                "timestep_fraction": 0.25,
                "output_times_s": [1, 2],
            }
        ],
    }
    left, right = tmp_path / "left.json", tmp_path / "right.json"
    left.write_text(json.dumps(payload))
    right.write_text(json.dumps(payload))
    a: dict[str, object] = {
        "input_path": str(left),
        "parent_path": "parent",
        "parent_metadata": {},
        "cold_wall_impulse_rate_ns_per_s": 1,
        "cold_head_arrival_s": 3,
    }
    b = {**a, "input_path": str(right)}
    assert matched_inputs(a, b)
    # Same filenames and summary values must not disguise different conserved arrays.
    changed = json.loads(right.read_text())
    changed["restarts"][0]["state"]["positions"][1] = 1.01
    right.write_text(json.dumps(changed))
    assert not matched_inputs(a, b)
    changed = json.loads(left.read_text())
    changed["restarts"][0]["parcels"][0]["gas"]["gas_constant"] = 1.01
    right.write_text(json.dumps(changed))
    assert not matched_inputs(a, b)
