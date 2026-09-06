"""Closed-parcel energy tests for the local equilibrium mixing screen."""

import pytest

pytest.importorskip("CoolProp")

import CoolProp.CoolProp as CP

from puffsat import eos_water
from puffsat.water_plate.thermodynamics import energy_reference, mixed_state


def test_liquid_vapor_alignment_preserves_measured_phase_energy_difference() -> None:
    reference = energy_reference()
    real_difference = float(CP.PropsSI("U", "D", 0.01, "T", 400.0, "Water")) - float(
        CP.PropsSI("U", "P", 101325.0, "T", 293.15, "Water")
    )
    gas = eos_water.pressure_energy(0.01, 400.0)[1]
    assert gas - reference.liquid_energy_j_kg == pytest.approx(real_difference)
    assert 2.4e6 < real_difference < 2.6e6
    assert reference.liquid_energy_j_kg < 0  # Relative to bound gas, not an invalid state.


@pytest.mark.parametrize("k", [5.0, 8.5, 10.0])
@pytest.mark.parametrize("injection_speed", [0.0, 100.0])
def test_local_merge_conserves_total_energy_and_signed_momentum(
    k: float, injection_speed: float
) -> None:
    row = mixed_state(45_580.0, k, 1.0, injection_speed_m_s=injection_speed)
    mass = row.projectile_kg * (1 + k)
    assert mass * row.merge_velocity_m_s == pytest.approx(
        row.projectile_kg * (45_580.0 - k * injection_speed)
    )
    final_energy = mass * (
        eos_water.pressure_energy(row.mixed_density_kg_m3, row.temperature_k)[1]
        + 0.5 * row.merge_velocity_m_s**2
    )
    assert final_energy == pytest.approx(row.total_input_energy_j, rel=1e-8)
    assert row.thermal_energy_j_kg > 0
    assert row.ionization_energy_j_kg >= 0
    assert 0 <= row.bond_fraction <= 1


def test_no_merge_heating_when_no_water_is_added() -> None:
    row = mixed_state(45_580.0, 0.0, 0.1241518740905039)
    assert row.merge_heat_j == 0
    assert row.temperature_k == pytest.approx(400.0, rel=1e-8)


def test_matching_point_uncertainty_is_small_but_not_hidden() -> None:
    cold = energy_reference(match_temperature_k=400.0)
    warm = energy_reference(match_temperature_k=600.0)
    assert abs(cold.liquid_energy_j_kg - warm.liquid_energy_j_kg) < 20_000.0
