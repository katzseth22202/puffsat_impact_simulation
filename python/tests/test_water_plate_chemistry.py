"""Acceptance checks for the conditional chemistry/expansion screen."""

import math

import pytest

pytest.importorskip("CoolProp")

from puffsat import eos_water, expansion, recombination
from puffsat.water_plate import chemistry
from puffsat.water_plate.thermodynamics import mixed_state


def test_hemisphere_clock_conserves_mass_and_matches_density_derivative() -> None:
    mass, rho0, speed = 237.5, 1.0, 10_000.0
    clock = chemistry.hemisphere_clock(mass, rho0, rho0 / 8, speed)
    r0 = (3 * mass / (2 * math.pi * rho0)) ** (1 / 3)
    assert clock.radius_m == pytest.approx(2 * r0)
    assert clock.time_s == pytest.approx(r0 / speed)
    assert (2 * math.pi / 3) * clock.radius_m**3 * rho0 / 8 == pytest.approx(mass)
    dt = 1e-6 * clock.time_s
    density_before = mass / ((2 * math.pi / 3) * (clock.radius_m - speed * dt) ** 3)
    density_after = mass / ((2 * math.pi / 3) * (clock.radius_m + speed * dt) ** 3)
    assert clock.expansion_time_s == pytest.approx(
        2 * dt / math.log(density_before / density_after), rel=1e-6
    )
    slower = chemistry.hemisphere_clock(mass, rho0, rho0 / 8, speed / 2)
    assert slower.time_s == pytest.approx(2 * clock.time_s)
    assert slower.expansion_time_s == pytest.approx(2 * clock.expansion_time_s)


def test_release_crossing_is_log_interpolated_and_missing_is_explicit() -> None:
    states = [
        expansion.ExpansionState(1, 10000, 1e7, 1e8),
        expansion.ExpansionState(0.01, 2500, 1e4, 2e7),
    ]
    crossing = chemistry.release_crossing(states, [100.0, 0.0], 0.5)
    assert crossing is not None
    assert crossing.rho == pytest.approx(0.1)
    assert crossing.temp == pytest.approx(5000)
    assert chemistry.release_crossing(states, [100.0, 80.0], 0.5) is None


def test_store_accounting_and_reaction_clocks_do_not_credit_photons_as_heat() -> None:
    mixed = mixed_state(45_580.0, 8.5, 1.0)
    initial = expansion.ExpansionState(
        mixed.mixed_density_kg_m3,
        mixed.temperature_k,
        mixed.pressure_pa,
        mixed.mixed_internal_energy_j_kg,
    )
    row = chemistry.diagnose(initial, mixed)
    assert row.bond_returned_fraction == pytest.approx(0, abs=1e-8)
    assert row.ion_returned_fraction == pytest.approx(0, abs=1e-8)
    assert row.internal_work_j_kg == 0
    assert row.bond_energy_j_kg == pytest.approx(mixed.bond_energy_j_kg)
    assert row.ion_energy_j_kg == pytest.approx(mixed.ionization_energy_j_kg)
    comp = eos_water.composition(initial.rho, initial.temp)
    assert 1 / row.tau_ion_three_body_s == pytest.approx(
        recombination.three_body_coefficient(initial.temp) * comp.n_e**2
    )
    assert row.tau_ion_radiative_s > row.tau_ion_three_body_s
    assert row.tau_bond_equilibrium_oh_s > row.tau_bond_h_proxy_s


@pytest.mark.parametrize("args", [(0, 1, 0.1, 1), (1, 1, 2, 1), (1, 1, 0.1, 0)])
def test_clock_rejects_invalid_inputs(args: tuple[float, float, float, float]) -> None:
    with pytest.raises(ValueError):
        chemistry.hemisphere_clock(*args)
