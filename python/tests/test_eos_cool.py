"""Analytic-limit acceptance tests for the cool-gas two-phase water EOS (Rung C, C1).

These pin the low-v physics that governs `e_eff` at 3.2 km/s: the superheated-vapor ideal-gas limit,
the saturation-pressure collapse and latent-heat content in the two-phase dome (the bulk
condensation channel, ADR-0004), the two-phase sound speed (where CoolProp's analytic value is
undefined), the condensed-fraction `liquid_frac` the wall-sticking sink will read, and the loader
invariants (positivity, monotone `e(T)`). CoolProp (the `sci` extra) is required, so the module is
skipped when it is absent.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("CoolProp")

import CoolProp.CoolProp as CP

from puffsat import eos_cool as ec

R_WATER = 461.5  # specific gas constant of water vapor [J/kg/K]


def test_superheated_vapor_is_ideal_gas() -> None:
    """Far below saturation the vapor obeys `p ≈ rho R_water T` and carries no condensate."""
    rho, temp = 0.05, 600.0  # rho << sat-vapor density at 600 K, so single-phase superheated vapor
    p, e = ec.pressure_energy(rho, temp)
    np.testing.assert_allclose(p, rho * R_WATER * temp, rtol=0.02)
    assert e > 0.0
    assert ec.liquid_fraction(rho, temp) == 0.0


def test_two_phase_pressure_is_saturation_pressure() -> None:
    """Inside the dome the equilibrium pressure collapses to `p_sat(T)` (bulk collapse)."""
    temp = 450.0
    rho_g = CP.PropsSI("D", "T", temp, "Q", 1, "Water")
    rho_f = CP.PropsSI("D", "T", temp, "Q", 0, "Water")
    rho_mid = 0.5 * (rho_g + rho_f)  # deep in the two-phase region
    p_sat = CP.PropsSI("P", "T", temp, "Q", 0, "Water")
    p, _ = ec.pressure_energy(rho_mid, temp)
    np.testing.assert_allclose(p, p_sat, rtol=1e-6)
    assert 0.0 < ec.liquid_fraction(rho_mid, temp) < 1.0


def test_latent_heat_lives_in_energy() -> None:
    """Crossing the dome at fixed T, `e` rises by water's (large) internal-energy latent heat."""
    temp = 400.0
    rho_g = CP.PropsSI("D", "T", temp, "Q", 1, "Water")
    rho_f = CP.PropsSI("D", "T", temp, "Q", 0, "Water")
    _, e_vap = ec.pressure_energy(rho_g, temp)  # saturated vapor
    _, e_liq = ec.pressure_energy(rho_f, temp)  # saturated liquid
    # Internal-energy latent heat of water near 400 K is ~2 MJ/kg; assert it is unmistakably there.
    assert e_vap - e_liq > 1.5e6


def test_sound_speed_matches_coolprop_single_phase_and_is_low_in_dome() -> None:
    """`c_s` matches CoolProp's `A` in single phase; small but positive (Wood dip) in the dome."""
    rho, temp = 2.727, 800.0  # superheated vapor (single phase)
    np.testing.assert_allclose(
        ec.sound_speed(rho, temp), CP.PropsSI("A", "D", rho, "T", temp, "Water"), rtol=1e-3
    )

    temp = 500.0
    rho_g = CP.PropsSI("D", "T", temp, "Q", 1, "Water")
    rho_f = CP.PropsSI("D", "T", temp, "Q", 0, "Water")
    cs_two_phase = ec.sound_speed(0.5 * (rho_g + rho_f), temp)
    assert 0.0 < cs_two_phase < CP.PropsSI("A", "D", rho_g, "T", temp, "Water")  # Wood dip


def test_liquid_fraction_spans_the_phases() -> None:
    """`liquid_frac` is 0 in vapor/supercritical, 1 in compressed liquid, and (0,1) in the dome."""
    # supercritical (T > Tcrit = 647 K): single phase, no distinct condensate.
    assert ec.liquid_fraction(50.0, 800.0) == 0.0
    # compressed liquid: rho above the saturated-liquid density at a sub-critical T.
    temp = 400.0
    rho_f = CP.PropsSI("D", "T", temp, "Q", 0, "Water")
    assert ec.liquid_fraction(rho_f * 1.01, temp) == 1.0
    # two-phase: lever-rule fraction.
    rho_g = CP.PropsSI("D", "T", temp, "Q", 1, "Water")
    q = CP.PropsSI("Q", "D", 0.5 * (rho_g + rho_f), "T", temp, "Water")
    np.testing.assert_allclose(ec.liquid_fraction(0.5 * (rho_g + rho_f), temp), 1.0 - q, rtol=1e-9)


def test_grid_fields_positive_and_energy_monotone() -> None:
    """`p, e, c_s > 0` and `e(T)` strictly increasing at fixed rho (the loader's requirements)."""
    rho_grid = np.array([0.05, 0.5, 5.0])
    t_grid = np.geomspace(300.0, 1700.0, 24)
    p, e, cs, lf = ec.eos_grid_lowv(rho_grid, t_grid)

    assert np.all(p > 0.0)
    assert np.all(e > 0.0)
    assert np.all(cs > 0.0)
    assert np.all((lf >= 0.0) & (lf <= 1.0))
    assert np.all(np.diff(e, axis=1) > 0.0)  # e increases with T along every density row
