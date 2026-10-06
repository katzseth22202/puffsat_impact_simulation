"""Limits the arriving-gas collisionality check must reproduce (R3)."""

from __future__ import annotations

import math

import pytest

from puffsat.water_plate import collisionality as co


@pytest.mark.parametrize("b", [1e-12, 1e-11, 1e-10])
def test_deflection_reproduces_rutherford(b: float) -> None:
    """Bare Coulomb has `tan(chi/2) = A / (2 E b)`; the integrator must match it."""
    pot = co.Potential(1, 1, screened=False)
    energy = 100.0 * co.EV
    exact = 2.0 * math.atan(co.COULOMB_J_M / (2.0 * energy * b))
    assert co.deflection(b, energy, pot) == pytest.approx(exact, rel=1e-8)


def test_head_on_collision_backscatters() -> None:
    pot = co.Potential(co.Z_O, co.Z_AR)
    assert co.deflection(0.0, 100.0 * co.EV, pot) == pytest.approx(math.pi, rel=1e-6)


def test_cross_section_shrinks_with_energy_and_sits_near_atomic_size() -> None:
    pot = co.Potential(co.Z_O, co.Z_AR)
    low = co.momentum_transfer_cross_section(100.0 * co.EV, pot)
    high = co.momentum_transfer_cross_section(300.0 * co.EV, pot)
    assert high < low
    assert 5e-21 < high < low < 5e-20


def test_hard_sphere_at_speed_is_the_smaller_value() -> None:
    """`pi r0^2` must undercut `sigma_mt`, or it is not the generous case."""
    pot = co.Potential(co.Z_O, co.Z_AR)
    energy = 200.0 * co.EV
    assert co.hard_sphere_at_speed(energy, pot) < co.momentum_transfer_cross_section(energy, pot)


def test_gaussian_fractions_limits() -> None:
    assert co.fraction_below_gaussian_1d(1.0) == 1.0
    assert co.fraction_below_gaussian_3d(1.0) == 1.0
    # At rho0/e^(1/2) a 1-D Gaussian is at one sigma: erfc(1/sqrt 2) of the mass lies outside.
    assert co.fraction_below_gaussian_1d(math.exp(-0.5)) == pytest.approx(math.erfc(math.sqrt(0.5)))
    assert co.fraction_below_exponential(1e-5) == 1e-5


def test_design_point_is_collisional_by_five_decades() -> None:
    for speed in co.SPEEDS_M_S:
        row = co.evaluate(speed)
        assert row.kn_neutral_generous < 1e-5
        assert row.relaxation_lengths_in_column > 1e5
