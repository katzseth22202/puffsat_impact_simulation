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
    n_o = y.y_h2o + y.y_o + sum(y.y_o_ions) + y.y_oh + 2.0 * y.y_o2
    n_h = 2.0 * y.y_h2o + y.y_h + y.y_hp + y.y_oh + 2.0 * y.y_h2
    np.testing.assert_allclose(n_o, 1.0, rtol=1e-8)
    np.testing.assert_allclose(n_h, 2.0, rtol=1e-8)
    charge = y.y_hp + sum((k + 1) * n for k, n in enumerate(y.y_o_ions))
    np.testing.assert_allclose(y.y_e, charge, rtol=1e-6)
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


def test_frozen_fractions_carry_the_full_o_ladder_when_hot() -> None:
    """Freezing at a 69 km/s-class turnaround state (multi-charge O) must keep every ladder
    stage: O nuclei sum to 1 per formula unit *including O2+ and up*, and neutrality counts
    each ion's charge. The single-stage bookkeeping loses both."""
    y = ew.frozen_composition(0.5, 60_000.0)
    o_nuclei = y.y_h2o + y.y_o + sum(y.y_o_ions)
    np.testing.assert_allclose(o_nuclei, 1.0, rtol=1e-6)
    charge = y.y_hp + sum((k + 1) * frac for k, frac in enumerate(y.y_o_ions))
    np.testing.assert_allclose(y.y_e, charge, rtol=1e-5)
    # The state is chosen hot enough that O2+ (and beyond) actually holds population.
    assert sum(y.y_o_ions[1:]) > 0.1


def test_frozen_at_hot_reference_state_matches_equilibrium() -> None:
    """Splice continuity must hold at a multi-charge freeze state (the 69 km/s bracket), not
    just the transitional-grid ~12 kK states — every ladder stage's cumulative ionization
    energy has to be carried frozen."""
    rho, temp = 0.5, 60_000.0
    y = ew.frozen_composition(rho, temp)
    p_eq, e_eq = ew.pressure_energy(rho, temp)
    p_fr, e_fr = ew.pressure_energy_frozen(rho, temp, y)
    np.testing.assert_allclose(p_fr, p_eq, rtol=1e-8)
    np.testing.assert_allclose(e_fr, e_eq, rtol=1e-8)


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


# --- Molecular intermediates: OH, H2, O2 (2026-08-25) -------------------------------------------
#
# Q-P showed the 2000-6000 K transition region is load-bearing: the cold leg's exhaust freezes at
# 3908 K, and the water-reformation rate is limited by OH, which the atoms-only species set could
# only proxy with `n_H`. These pin the three added equilibria against their closed forms, check the
# element bookkeeping the extra species change, and -- the point of the exercise -- check that both
# endpoints the original model got right are still right.


def test_diatomic_constant_matches_closed_form() -> None:
    """`ln_k_diatomic` reproduces `K = z_A z_B / z_AB * exp(-D0/kT)` for OH <=> H + O."""
    temp = 3000.0
    d = ew.OH
    z_a = (2.0 * np.pi * ew.M_H * ew.K_B * temp / ew.H_PLANCK**2) ** 1.5 * ew.G_H
    z_b = (2.0 * np.pi * ew.M_O * ew.K_B * temp / ew.H_PLANCK**2) ** 1.5 * ew.G_O
    z_rot = temp / (d.symmetry * d.theta_rot)
    z_vib = 1.0 / (1.0 - np.exp(-d.theta_vib / temp))
    z_ab = (
        (2.0 * np.pi * d.mass * ew.K_B * temp / ew.H_PLANCK**2) ** 1.5
        * d.degeneracy
        * z_rot
        * z_vib
    )
    expected = z_a * z_b / z_ab * np.exp(-d.d0 / (ew.K_B * temp))
    np.testing.assert_allclose(np.exp(ew.ln_k_diatomic(d, temp)), expected, rtol=1e-12)


def test_first_oh_bond_of_water_is_the_difference_of_two_independent_sources() -> None:
    """`D_AT - D0(OH)` must be water's first O-H bond, ~492 kJ/mol.

    `D_AT` and `D0_OH` come from the same JANAF 0 K heats of formation but enter the code by
    different routes, so this is the check that a sign or a wrong dissociation limit cannot hide.
    """
    first_bond_kj = (ew.D_AT - ew.OH.d0) * ew.N_A / 1e3
    assert 488.0 < first_bond_kj < 496.0, first_bond_kj
    # And the two O-H bonds must sum back to the atomization energy.
    np.testing.assert_allclose(first_bond_kj + ew.OH.d0 * ew.N_A / 1e3, ew.D_AT * ew.N_A / 1e3)


def test_intermediates_peak_in_the_transition_band_and_vanish_at_both_ends() -> None:
    """OH is a real population between ~2000 and ~6000 K, and negligible outside it.

    This is the whole reason the species set was extended: at 400 K everything is bound H2O and at
    15 kK everything is atoms, so the original two-endpoint model was right at both ends and wrong
    only in between -- which is exactly where the cold leg's exhaust freezes.
    """
    rho = 0.32
    n_f = rho / ew.M_H2O
    y_oh = {t: ew.composition(rho, t).n_oh / n_f for t in (400.0, 2000.0, 3000.0, 4000.0, 15_000.0)}
    assert y_oh[400.0] < 1e-8
    assert y_oh[15_000.0] < 1e-3
    assert max(y_oh[2000.0], y_oh[3000.0], y_oh[4000.0]) > 0.02
    # The band is a peak, not a plateau: OH is 2+ decades down at both ends of the sweep.
    assert y_oh[15_000.0] < 0.01 * max(y_oh[2000.0], y_oh[3000.0], y_oh[4000.0])


def test_element_conservation_holds_with_the_intermediates_carrying_nuclei() -> None:
    """H:O = 2:1 and the O-nucleus count, summed over *all* species including OH, H2, O2."""
    rho = 0.32
    n_f = rho / ew.M_H2O
    for temp in (400.0, 1500.0, 2500.0, 3500.0, 5000.0, 9000.0, 25_000.0):
        c = ew.composition(rho, temp)
        n_h_nuclei = 2.0 * c.n_h2o + c.n_oh + 2.0 * c.n_h2 + c.n_h + c.n_hp
        n_o_nuclei = c.n_h2o + c.n_oh + 2.0 * c.n_o2 + c.n_o + sum(c.n_o_ions)
        np.testing.assert_allclose(n_o_nuclei, n_f, rtol=1e-8, err_msg=f"O nuclei at {temp} K")
        np.testing.assert_allclose(n_h_nuclei, 2.0 * n_f, rtol=1e-8, err_msg=f"H:O at {temp} K")


def test_hot_endpoint_is_unmoved_by_the_intermediates() -> None:
    """Across the whole `f(v)` stagnation regime the intermediates hold a negligible *energy* share.

    They do not vanish there -- at 0.32 kg/m^3 and 15 kK about 1.4e-3 of the oxygen is still in a
    molecule, because `n_H n_O` is enormous even when the Boltzmann factor is not -- but what the
    restitution depends on is the specific energy, and the bond energy they hold back is under
    0.1% of it everywhere the tables are used. This is the check that the extension refines the
    transition band without disturbing the shipped regime.
    """
    n_f_per_rho = 1.0 / ew.M_H2O
    for rho in (0.032, 0.32, 3.2):
        for temp in (15_000.0, 20_000.0, 30_000.0, 40_000.0):
            c = ew.composition(rho, temp)
            _, e = ew.pressure_energy(rho, temp)
            e_bond_held = (c.n_oh * ew.OH.d0 + c.n_h2 * ew.H2.d0 + c.n_o2 * ew.O2.d0) / rho
            assert e_bond_held / e < 5e-3, f"rho={rho}, T={temp}: {e_bond_held / e:.2e}"
            assert (c.n_oh + c.n_h2 + c.n_o2) / (rho * n_f_per_rho) < 2e-2


def test_chemical_energy_reduces_to_the_atomization_form_when_intermediates_vanish() -> None:
    """Where OH, H2, O2 are absent the generalized chemical energy must equal `(n_f - n_H2O) D_AT`.

    The generalization is `n_f D_AT - n_H2O D_AT - sum n_i D0_i`; this pins the reduction, so a
    misplaced term shows up as a hot-state energy error rather than silently rescaling `e_eff`.
    """
    n_f = 1.0e25
    atoms_only = ew.Composition(
        n_h2o=0.1 * n_f, n_h=1.8 * n_f, n_o=0.9 * n_f, n_hp=0.0, n_o_ions=(0.0,) * 8, n_e=0.0
    )
    np.testing.assert_allclose(
        ew._bond_energy_held(atoms_only, n_f), (n_f - atoms_only.n_h2o) * ew.D_AT, rtol=1e-12
    )
    # And the one state that can be checked against a hand number: all the water as OH + H is
    # water's *first* O-H bond, 492 kJ/mol, not the full 918 the atoms-only form would charge.
    all_oh = ew.Composition(
        n_h2o=0.0, n_h=n_f, n_o=0.0, n_hp=0.0, n_o_ions=(0.0,) * 8, n_e=0.0, n_oh=n_f
    )
    np.testing.assert_allclose(
        ew._bond_energy_held(all_oh, n_f) * ew.N_A / (n_f * 1e3), 492.1, rtol=1e-3
    )


def test_intermediates_split_water_earlier_but_hold_the_atoms_longer() -> None:
    """The physical signature of the added species, in the direction that matters for the freeze.

    Breaking H2O to OH + H costs 492 kJ/mol, not the full 918, so water comes apart at a *lower*
    temperature than the atoms-only model said -- while the O-H bond that survives keeps the gas
    from being fully atomic until hotter. At 2500 K the free-atom fraction must therefore sit
    below the bound-H2O deficit: some of the broken water is held as OH, not as atoms.
    """
    rho, temp = 0.32, 2500.0
    n_f = rho / ew.M_H2O
    c = ew.composition(rho, temp)
    broken = 1.0 - c.n_h2o / n_f
    free_o = (c.n_o + sum(c.n_o_ions)) / n_f
    assert broken > free_o, "intermediates must hold part of the broken water"
    assert c.n_oh > c.n_o, "OH is the reservoir doing the holding at 2500 K"
