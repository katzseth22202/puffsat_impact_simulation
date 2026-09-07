"""Acceptance checks for the cold-liquid/hot-gas flow EOS bridge."""

import numpy as np
import pytest

pytest.importorskip("CoolProp")

import CoolProp.CoolProp as CP

from puffsat import eos_water
from puffsat.water_plate import flow_eos
from puffsat.water_plate.thermodynamics import energy_reference


def test_cold_and_hot_limits_preserve_the_same_phase_energy_bill() -> None:
    reference = energy_reference()
    for rho, temp in ((4.81, 293.15), (0.01, 400.0)):
        p, e = flow_eos.pressure_energy(rho, temp)
        assert p == pytest.approx(float(CP.PropsSI("P", "D", rho, "T", temp, "Water")))
        assert e - flow_eos.ENERGY_SHIFT == pytest.approx(
            float(CP.PropsSI("U", "D", rho, "T", temp, "Water")) + reference.offset_j_kg
        )
    p, e = flow_eos.pressure_energy(1.0, 16000.0)
    expected_p, expected_e = eos_water.pressure_energy(1.0, 16000.0)
    residual_p, residual_e, _ = flow_eos.dense_residual(1.0, 16000.0)
    assert p == pytest.approx(expected_p + residual_p)
    assert e - flow_eos.ENERGY_SHIFT == pytest.approx(expected_e + residual_e)
    assert abs(residual_p / expected_p) < 1e-3
    assert abs(residual_e / expected_e) < 1e-3


@pytest.mark.parametrize("rho", [0.001, 0.1, 1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0])
def test_bridge_has_positive_heat_capacity_and_acoustic_speed(rho: float) -> None:
    for temp in (300.0, 600.0, 1000.0, 1200.0, 1500.0, 1800.0, 10000.0):
        assert (
            flow_eos.pressure_energy(rho, temp + 0.01)[1] > flow_eos.pressure_energy(rho, temp)[1]
        )
        assert flow_eos.sound_speed(rho, temp) > 0


@pytest.mark.parametrize("rho", [1.0, 303.822, 1000.0, 2000.0])
def test_join_preserves_first_derivatives_and_sound_speed(rho: float) -> None:
    temp, delta = flow_eos.COLD_MAX_T, 0.01
    left = flow_eos.pressure_energy(rho, temp - delta)
    middle = flow_eos.pressure_energy(rho, temp)
    right = flow_eos.pressure_energy(rho, temp + delta)
    for i in (0, 1):
        assert (middle[i] - left[i]) / delta == pytest.approx(
            (right[i] - middle[i]) / delta, rel=1e-3
        )
    assert flow_eos.sound_speed(rho, temp - delta) == pytest.approx(
        flow_eos.sound_speed(rho, temp + delta), rel=1e-3
    )


def test_liquid_source_uses_internal_energy_not_a_cold_vapor_state() -> None:
    rho = 4.81
    temp = flow_eos.liquid_source_temperature(rho)
    energy = flow_eos.pressure_energy(rho, temp)[1] - flow_eos.ENERGY_SHIFT
    assert energy == pytest.approx(energy_reference().liquid_energy_j_kg, abs=0.01)
    assert 273.16 < temp < 293.15  # Small equilibrium flash paid by sensible heat.


@pytest.mark.parametrize("rho,temp", [(0.01, 400.0), (1.0, 6000.0), (10.0, 16000.0)])
def test_gas_potential_reproduces_existing_pressure_and_energy(rho: float, temp: float) -> None:
    dr, dt = rho * 1e-4, temp * 1e-4
    potential = flow_eos.gas_free_energy
    expected_p, expected_e = eos_water.pressure_energy(rho, temp)
    assert rho**2 * (potential(rho + dr, temp) - potential(rho - dr, temp)) / (
        2 * dr
    ) == pytest.approx(expected_p, rel=1e-5)
    energy = potential(rho, temp) - temp * (
        potential(rho, temp + dt) - potential(rho, temp - dt)
    ) / (2 * dt)
    assert energy == pytest.approx(expected_e, rel=1e-5)


@pytest.mark.parametrize("rho", [0.01, 1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0])
def test_seam_obeys_maxwell_energy_pressure_identity(rho: float) -> None:
    for temp in (1100.0, 1624.95, 3000.0, 16000.0):
        pressure = flow_eos.pressure_energy(rho, temp)[0]
        _, er, pt, _ = flow_eos.derivatives(rho, temp)
        assert rho**2 * er == pytest.approx(pressure - temp * pt, abs=1e-4 * pressure)


def test_dense_transition_grid_rejects_the_previous_negative_cv_regression() -> None:
    # Sparse temperature endpoints missed the old blended-potential instability.
    # Cover the failed point and the full density/transition box, including midpoints.
    for rho in np.unique(np.append(np.geomspace(1e-5, 2000, 72), 303.822)):
        for temp in np.unique(np.append(np.linspace(1000, 1800, 33), 1624.95)):
            assert flow_eos.sound_speed(float(rho), float(temp)) > 0


@pytest.mark.parametrize("rho", [0.01, 303.822, 1000.0, 2000.0])
def test_extended_heat_capacity_is_the_integrated_model(rho: float) -> None:
    for temp in (1001.0, 1624.95, 10000.0):
        dt = temp * 1e-4
        gas_cv = (
            eos_water.pressure_energy(rho, temp + dt)[1]
            - eos_water.pressure_energy(rho, temp - dt)[1]
        ) / (2 * dt)
        expected = (
            gas_cv
            + flow_eos.residual_reference(rho).cv
            * (flow_eos.COLD_MAX_T / temp) ** flow_eos.RESIDUAL_POWER
        )
        assert flow_eos.derivatives(rho, temp)[3] == pytest.approx(expected, rel=1e-6)
