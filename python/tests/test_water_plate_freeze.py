"""Frozen coefficients must reproduce the existing species thermodynamics exactly."""

from itertools import pairwise

import pytest

pytest.importorskip("CoolProp")

from puffsat import eos_water as ew
from puffsat.water_plate.flow import CsvValue
from puffsat.water_plate.freeze import gas_parameters, sample_times
from puffsat.water_plate.freeze_report import comparison, valid_quad
from puffsat.water_plate.selective_report import ordered_comparison


@pytest.mark.parametrize("rho,reference", [(0.01, 1e5), (1.0, 16000.0), (300.0, 8000.0)])
def test_exported_parcel_thermal_modes_reproduce_frozen_eos(rho: float, reference: float) -> None:
    gas = gas_parameters(rho, reference)
    frozen = ew.frozen_composition(rho, reference)
    for temp in (1000.0, 3000.0, reference, 1e5):
        p, e = ew.pressure_energy_frozen(rho, temp, frozen)
        assert rho * gas.gas_constant * temp == pytest.approx(p, rel=1e-12)
        assert gas.thermal(temp) == pytest.approx(e, rel=1e-12)
    assert gas.chemical_energy == pytest.approx(gas.bond_store + gas.ion_store)


def controls() -> list[dict[str, CsvValue]]:
    rows: list[dict[str, CsvValue]] = []
    for case, branch, j in (
        ("combined", "equilibrium", 300.0),
        ("combined", "parcel_frozen", 280.0),
        ("unsprayed", "equilibrium", 200.0),
        ("unsprayed", "parcel_frozen", 190.0),
    ):
        rows.append(
            {
                "name": case + "_h1_n40_freeze_153.5us",
                "branch": branch,
                "status": "completed_time_window",
                "start_time_s": 0.0001535,
                "time_s": 0.0005,
                "wall_impulse_ns": j,
                "initial_wall_impulse_ns": 100.0,
                "max_energy_error_j": 0.001,
                "restart_physical_energy_j": 1e6,
                "momentum_error_ns": 0.0,
                "max_matching_pressure_correction_fraction": 0.001,
                "abs_matching_energy_offsets_j": 1000.0,
                "max_start_pressure_jump_fraction": 0.0,
                "max_start_temperature_jump_fraction": 0.0,
            }
        )
    return rows


def test_reaction_effect_subtracts_the_unsprayed_response_and_cold_control_once() -> None:
    rows = controls()
    assert valid_quad(rows, 0.1)
    result = comparison(rows, 1000.0, 0.1)
    assert result["reaction_effect_combined_ns"] == 20.0
    assert result["reaction_effect_unsprayed_ns"] == 10.0
    assert result["reaction_effect_on_gain_ns"] == 10.0
    assert result["equilibrium_gain_ns"] == 99.5
    assert result["frozen_gain_ns"] == 89.5


@pytest.mark.parametrize(
    "field,value",
    [
        ("status", "invalid_domain"),
        ("time_s", 0.0004),
        ("max_energy_error_j", 1e4),
        ("max_energy_error_j", float("nan")),
        ("max_start_pressure_jump_fraction", float("nan")),
        ("max_matching_pressure_correction_fraction", 0.02),
        ("abs_matching_energy_offsets_j", 1e4),
        ("restart_physical_energy_j", 1.00001e6),
        ("name", "combined_h1_n80_freeze_153.5us"),
    ],
)
def test_invalid_or_unmatched_switches_cannot_report_a_gain(field: str, value: CsvValue) -> None:
    rows = controls()
    rows[1][field] = value
    assert not valid_quad(rows, 0.1)
    assert "reaction_effect_on_gain_ns" not in comparison(rows, 1000.0, 0.1)


@pytest.mark.parametrize("start,end", [(153.5e-6, 0.0005), (0.0003, 0.001), (0.0001, 0.0002)])
def test_restart_samples_cover_only_the_requested_window(start: float, end: float) -> None:
    samples = sample_times(start, end)
    assert samples[0] == start
    assert samples[-1] == end
    assert all(a < b for a, b in pairwise(samples))


@pytest.mark.parametrize("start,end", [(0.001, 0.001), (0.002, 0.001), (0, 1), (1, float("nan"))])
def test_invalid_restart_window_is_rejected(start: float, end: float) -> None:
    with pytest.raises(ValueError, match="switch"):
        sample_times(start, end)


def molecular_controls() -> list[dict[str, CsvValue]]:
    rows = controls()
    for row in rows:
        if row["branch"] == "parcel_frozen":
            row["branch"] = "molecular_frozen"
            row["wall_impulse_ns"] = float(row["wall_impulse_ns"]) + (
                8.0 if str(row["name"]).startswith("combined") else 3.0
            )
    return rows


def test_ordered_differences_telescope_without_confusing_frozen_closures() -> None:
    molecular = molecular_controls()
    assert not valid_quad(molecular, 0.1)
    result = ordered_comparison(controls(), molecular, 1000.0, 0.1)
    assert result["all_chemistry_effect_ns"] == 10.0
    assert result["molecular_effect_with_ionization_active_ns"] == 5.0
    assert result["ionization_effect_with_molecules_frozen_ns"] == 5.0


def test_ordered_difference_preserves_negative_relative_ionization_response() -> None:
    molecular = molecular_controls()
    molecular[1]["wall_impulse_ns"] = 278.0
    result = ordered_comparison(controls(), molecular, 1000.0, 0.1)
    assert result["molecular_effect_with_ionization_active_ns"] == 15.0
    assert result["ionization_effect_with_molecules_frozen_ns"] == -5.0
    assert result["all_chemistry_effect_ns"] == 10.0


@pytest.mark.parametrize(
    "field,value",
    [
        ("time_s", 0.0006),
        ("wall_impulse_ns", 300.01),
        ("name", "combined_h1_n80_freeze_153.5us"),
        ("status", "invalid_domain"),
    ],
)
def test_ordered_difference_rejects_mismatched_controls(field: str, value: CsvValue) -> None:
    rows = molecular_controls()
    rows[0][field] = value
    assert "all_chemistry_effect_ns" not in ordered_comparison(controls(), rows, 1000.0, 0.1)
