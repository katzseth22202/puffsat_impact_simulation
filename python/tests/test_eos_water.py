"""Analytic-limit acceptance tests for the equilibrium water EOS (B5c-1).

These pin the physics that governs `e_eff`: the cold-vapor ideal-gas limit, the Saha and
dissociation building blocks against their closed forms, and the integrated invariants the Rust
table loader requires (positivity, monotone `e(T)`) plus a positive sound speed.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np

from puffsat import eos_water as ew

R_WATER = ew.K_B / ew.M_H2O  # specific gas constant of water vapor [J/kg/K] = 461.5


def test_saha_constant_matches_closed_form() -> None:
    """`ln_k_saha` reproduces `K = (g_i g_e/g_n)(2 pi m_e k T/h^2)^{3/2} exp(-IP/kT)`."""
    temp = 12_000.0
    expected = (
        (ew.G_HP * ew.G_E / ew.G_H)
        * (2.0 * np.pi * ew.M_E * ew.K_B * temp / ew.H_PLANCK**2) ** 1.5
        * np.exp(-ew.IP_H / (ew.K_B * temp))
    )
    got = np.exp(ew.ln_k_saha(ew.IP_H, ew.G_HP, ew.G_H, temp))
    np.testing.assert_allclose(got, expected, rtol=1e-12)


def test_low_temperature_is_ideal_water_vapor() -> None:
    """At 400 K the gas is essentially undissociated H2O and obeys `p = rho R_water T`."""
    rho, temp = 0.32, 400.0
    comp = ew.composition(rho, temp)
    n_f = rho / ew.M_H2O
    # Dissociation and ionization are exponentially suppressed: vapor is ~pure H2O.
    assert comp.n_h2o / n_f > 0.999
    assert comp.n_h / n_f < 1e-6
    assert comp.n_e / n_f < 1e-12

    p, e = ew.pressure_energy(rho, temp)
    np.testing.assert_allclose(p, rho * R_WATER * temp, rtol=1e-3)
    assert e > 0.0
    # Translational + rotational gives ~3 R_water T; vibration is nearly frozen at 400 K.
    np.testing.assert_allclose(e, 3.0 * R_WATER * temp, rtol=0.05)


def test_dissociation_grows_with_temperature() -> None:
    """The dissociated fraction 1 - n_H2O/n_f rises monotonically from ~0 (cold) toward ~1 (hot)."""
    rho = 0.32
    n_f = rho / ew.M_H2O
    temps = [800.0, 2000.0, 4000.0, 8000.0]
    frac = [1.0 - ew.composition(rho, t).n_h2o / n_f for t in temps]
    assert frac[0] < 1e-3  # cold: bound
    assert frac[-1] > 0.9  # hot: dissociated
    assert all(b > a for a, b in pairwise(frac))


def test_ionization_and_charge_neutrality_at_high_temperature() -> None:
    """At 30 kK the plasma is meaningfully ionized and electrically neutral."""
    rho, temp = 0.32, 30_000.0
    comp = ew.composition(rho, temp)
    n_f = rho / ew.M_H2O
    assert comp.n_e / n_f > 0.1  # appreciable free-electron density
    np.testing.assert_allclose(comp.n_hp + comp.n_op, comp.n_e, rtol=1e-6)


def test_energy_monotone_and_fields_positive() -> None:
    """`e(T)` is strictly increasing at fixed rho; p, e, c_s are all positive across the grid.

    Positivity + monotone `e(T)` are exactly what the Rust loader's log-interpolation requires
    (ADR-0007); strictly increasing `e(T)` is also the physical statement `c_v > 0`.
    """
    rho_grid = np.array([0.16, 0.32, 0.64])
    t_grid = np.geomspace(300.0, 50_000.0, 40)
    p, e, cs = ew.eos_grid(rho_grid, t_grid)

    assert np.all(p > 0.0)
    assert np.all(e > 0.0)
    assert np.all(cs > 0.0)
    # e increases with T along every density row.
    assert np.all(np.diff(e, axis=1) > 0.0)


def test_sound_speed_cold_vapor_is_physical() -> None:
    """Cold-vapor sound speed sits near the ideal `sqrt(gamma R_water T)` (gamma ~ 1.3 for H2O)."""
    cs = ew.sound_speed(0.32, 400.0)
    ideal = np.sqrt(1.3 * R_WATER * 400.0)
    assert 0.7 * ideal < cs < 1.3 * ideal
