"""Cold-path contracts for nonlinear parcel source export and net-return reports."""

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from puffsat import eos_water as ew
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate.parcel import chemistry, fixture, prepare
from puffsat.water_plate.parcel_report import summarize


def test_rust_fixture_is_reproducible_from_independent_python_thermochemistry() -> None:
    saved = json.loads(Path("crates/water_plate/tests/fixtures/neutral_parcel.json").read_text())
    generated = json.loads(json.dumps(fixture()))
    assert saved == generated


def test_li_replaces_all_three_hydrogen_channels_and_preserves_units() -> None:
    model = chemistry("li_h2")
    rates = model["reactions"]
    assert isinstance(rates, list)
    assert len(rates) == 8
    h2 = [r for r in rates if r["reactants"] == (0, 0, 1, 0, 0, 0)]
    assert len(h2) == 1
    rate = h2[0]
    assert rate["products"] == (2, 0, 0, 0, 0, 0)
    assert rate["a"] == pytest.approx(4.577e19 * 1e-6 / ew.N_A)
    assert rate["colliders"] == (1, 1, 2.5, 1, 1, 12)
    assert not any(r["reactants"] == (2, 0, 0, 0, 0, 0) for r in rates)


def test_group_scales_do_not_scale_exchanges_or_duplicate_steam_channels() -> None:
    base = chemistry()["reactions"]
    slow = chemistry("association_1e6_slower")["reactions"]
    removed = chemistry("without_direct_and_steam_h2")["reactions"]
    assert isinstance(base, list) and isinstance(slow, list) and isinstance(removed, list)
    for r, a, b, c in zip(kn.REACTIONS, base, slow, removed, strict=True):
        assert b["a"] / a["a"] == pytest.approx(1e-6 if r.colliders is not None else 1)
        assert c["a"] / a["a"] == pytest.approx(0 if r.source_number in (14, 15) else 1)
        nu = np.array(b["products"]) - np.array(b["reactants"])
        np.testing.assert_array_equal(kn.ATOM_COUNTS @ nu, [0, 0])


def record(mass: float, net: float | None, error: float = 0.0) -> dict[str, Any]:
    return {
        "mass_kg": mass,
        "initial_bond_j_kg": 10e6,
        "initial_energy_j_kg": 30e6,
        "initial_ionized_nuclei_fraction": 0.0,
        "max_energy_error_j_kg": error,
        "parent_returns": [{"time_s": 193e-6, "net_bond_return_j_kg": 2e6}],
        "trajectory": {
            "accepted_steps": 3,
            "rejected_steps": 1,
            "samples": []
            if net is None
            else [
                {
                    "time_s": 193e-6,
                    "net_bond_return_j_kg": net,
                    "work_j_kg": -3e6,
                    "forcing_j_kg": 1e6,
                }
            ],
        },
    }


def test_net_return_is_signed_and_missing_mass_is_not_renormalized_away() -> None:
    result = summarize([record(2, 4e6), record(1, -2e6), record(1, None)], 193e-6)
    assert result["initial_bond_store_j"] == 40e6
    assert result["resolved_net_return_j"] == 6e6
    assert result["resolved_return_over_total_start_store"] == pytest.approx(0.15)
    assert result["resolved_initial_store_fraction"] == pytest.approx(0.75)
    assert result["including_unresolved_return_upper_j"] == 16e6
    assert result["including_unresolved_return_lower_j"] == pytest.approx(
        16e6 - ew.FULL_ATOMIZATION_ENERGY
    )
    assert result["resolved_work_j"] == -9e6
    assert result["resolved_parent_forcing_j"] == 3e6


def test_failed_energy_gate_remains_in_denominator() -> None:
    result = summarize([record(1, 4e6, error=100), record(1, 2e6)], 193e-6)
    assert result["resolved_mass_fraction"] == 0.5
    assert result["resolved_return_over_total_start_store"] == 0.1


@pytest.mark.parametrize("net", [math.nan, math.inf, 11e6, -ew.FULL_ATOMIZATION_ENERGY])
def test_invalid_completed_bond_store_is_rejected(net: float) -> None:
    with pytest.raises(ValueError, match="bond store"):
        summarize([record(1, net)], 193e-6)


@pytest.mark.parametrize("tolerance", [0, -1, math.nan, math.inf])
def test_invalid_solver_settings_fail_before_parent_io(tolerance: float) -> None:
    with pytest.raises(ValueError, match="solver settings"):
        prepare(Path("does_not_exist.jsonl"), 1, tolerance, 1e-7, ["baseline"])
