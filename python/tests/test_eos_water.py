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
    """At 30 kK the plasma is meaningfully ionized and electrically neutral (charge balance now
    sums the full O ladder, charge-weighted)."""
    rho, temp = 0.32, 30_000.0
    comp = ew.composition(rho, temp)
    n_f = rho / ew.M_H2O
    assert comp.n_e / n_f > 0.1  # appreciable free-electron density
    charge = comp.n_hp + sum((k + 1) * n for k, n in enumerate(comp.n_o_ions))
    np.testing.assert_allclose(charge, comp.n_e, rtol=1e-6)


def test_oxygen_ladder_climbs_and_banks_energy_at_very_high_temperature() -> None:
    """Jupiter-retrograde regime (69 km/s stagnation): at 2e5 K and dilute density the dominant O
    charge state is well past O+, charge neutrality holds across the ladder, and the specific
    energy carries a multi-ionization chemical sink well above the single-stage model's ceiling."""
    rho, temp = 0.1, 2.0e5
    comp = ew.composition(rho, temp)
    n_f = rho / ew.M_H2O
    # H fully stripped; O climbed past the first stage.
    assert comp.n_hp / (2.0 * n_f) > 0.99
    stages = np.array(comp.n_o_ions)
    assert stages.argmax() >= 2, f"expected dominant O charge >= 3, got {stages.argmax() + 1}"
    charge = comp.n_hp + sum((k + 1) * n for k, n in enumerate(comp.n_o_ions))
    np.testing.assert_allclose(charge, comp.n_e, rtol=1e-6)
    # The energy at 2e5 K must exceed the single-stage ceiling: thermal (all 13 particles were it
    # fully stripped would be ~1.6 GJ/kg) plus > 100 eV/molecule of ladder energy.
    _, e = ew.pressure_energy(rho, temp)
    e_single_stage_ceiling = (
        1.5 * ew.K_B * temp * 6.0 + 2.0 * ew.IP_H + ew.IP_O + ew.D_AT
    ) / ew.M_H2O
    assert e > e_single_stage_ceiling


def test_energy_monotone_to_megakelvin() -> None:
    """`e(T)` stays strictly increasing (c_v > 0) through the extended Jupiter-table T range —
    what the Rust loader's monotone inversion requires up to the new grid top."""
    rho_grid = np.array([1e-3, 0.07, 2.0])
    t_grid = np.geomspace(300.0, 1.2e6, 50)
    p, e, cs = ew.eos_grid(rho_grid, t_grid)
    assert np.all(p > 0.0)
    assert np.all(e > 0.0)
    assert np.all(cs > 0.0)
    assert np.all(np.diff(e, axis=1) > 0.0)


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


# ---- Frozen-composition EOS (sudden-freeze bounding runs, frozen-recombination check) ----------


def test_frozen_fractions_conserve_elements_and_charge() -> None:
    """The frozen fractions inherit the equilibrium invariants: O nuclei sum to 1 per formula
    unit, H:O stays 2:1, and the gas is electrically neutral."""
    y = ew.frozen_composition(1.0, 12_000.0)
    np.testing.assert_allclose(y.y_h2o + y.y_o + y.y_op, 1.0, rtol=1e-8)
    np.testing.assert_allclose(2.0 * y.y_h2o + y.y_h + y.y_hp, 2.0, rtol=1e-8)
    np.testing.assert_allclose(y.y_e, y.y_hp + y.y_op, rtol=1e-6)
    # The reference state is chosen hot enough to be meaningfully dissociated.
    assert y.y_h2o < 0.5


def test_frozen_at_reference_state_matches_equilibrium() -> None:
    """Freezing the composition at `(rho*, T*)` reproduces the equilibrium `p` and `e` *at that
    state* exactly — the splice's continuity requirement at the freeze instant."""
    rho, temp = 1.0, 12_000.0
    y = ew.frozen_composition(rho, temp)
    p_eq, e_eq = ew.pressure_energy(rho, temp)
    p_fr, e_fr = ew.pressure_energy_frozen(rho, temp, y)
    np.testing.assert_allclose(p_fr, p_eq, rtol=1e-8)
    np.testing.assert_allclose(e_fr, e_eq, rtol=1e-8)


def test_frozen_pure_h2o_has_no_chemical_sink() -> None:
    """`PURE_H2O_FROZEN` (freeze *before* the plate) is chemistry-free water vapor: it matches
    the equilibrium EOS cold, and stores far less energy hot (no dissociation/ionization sink)."""
    p_fr, e_fr = ew.pressure_energy_frozen(0.32, 400.0, ew.PURE_H2O_FROZEN)
    p_eq, e_eq = ew.pressure_energy(0.32, 400.0)
    np.testing.assert_allclose(p_fr, p_eq, rtol=1e-3)
    np.testing.assert_allclose(e_fr, e_eq, rtol=1e-3)

    _, e_hot_fr = ew.pressure_energy_frozen(0.32, 30_000.0, ew.PURE_H2O_FROZEN)
    _, e_hot_eq = ew.pressure_energy(0.32, 30_000.0)
    assert e_hot_fr < 0.5 * e_hot_eq


def test_frozen_composition_locks_chemical_energy_on_cooling() -> None:
    """Cooling a frozen dissociated gas returns only its thermal energy; the equilibrium path
    returns the chemical energy too — the whole point of the pessimistic bound."""
    rho, t_hot, t_cold = 1.0, 12_000.0, 600.0
    y = ew.frozen_composition(rho, t_hot)

    _, e_hot_eq = ew.pressure_energy(rho, t_hot)
    _, e_cold_eq = ew.pressure_energy(rho, t_cold)
    _, e_hot_fr = ew.pressure_energy_frozen(rho, t_hot, y)
    _, e_cold_fr = ew.pressure_energy_frozen(rho, t_cold, y)

    released_eq = e_hot_eq - e_cold_eq
    released_fr = e_hot_fr - e_cold_fr
    assert released_fr > 0.0
    # Equilibrium recovers the (large) chemical store on top of the thermal energy.
    assert released_eq > 1.5 * released_fr
    # The locked chemical energy survives in the cold frozen state.
    e_chem = (y.y_hp * ew.IP_H + y.y_op * ew.IP_O + (1.0 - y.y_h2o) * ew.D_AT) / ew.M_H2O
    assert e_chem > 0.0
    assert e_cold_fr > e_chem


def test_frozen_grid_positive_and_monotone() -> None:
    """The frozen EOS obeys the same loader invariants as the equilibrium one (ADR-0007):
    positive `p`, `e`, `c_s` and strictly increasing `e(T)` at fixed `rho`."""
    y = ew.frozen_composition(1.0, 12_000.0)
    rho_grid = np.array([0.16, 0.32, 0.64])
    t_grid = np.geomspace(300.0, 50_000.0, 40)
    p, e, cs = ew.eos_grid_frozen(rho_grid, t_grid, y)

    assert np.all(p > 0.0)
    assert np.all(e > 0.0)
    assert np.all(cs > 0.0)
    assert np.all(np.diff(e, axis=1) > 0.0)


def test_frozen_sound_speed_is_ideal_mixture_like() -> None:
    """At fixed composition the gas is an ideal mixture, so `c_s^2 = gamma p/rho` with
    `1 < gamma <= 5/3` (monatomic ceiling)."""
    y = ew.frozen_composition(1.0, 12_000.0)
    for temp in (2_000.0, 12_000.0, 40_000.0):
        p, _ = ew.pressure_energy_frozen(1.0, temp, y)
        cs = ew.sound_speed_frozen(1.0, temp, y)
        gamma = cs * cs * 1.0 / p
        assert 1.0 < gamma <= 5.0 / 3.0 + 1e-6
