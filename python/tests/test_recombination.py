"""Tests for the freeze-out criterion: does recombination keep up with the expansion?

`expansion.py` produces the cooling history on **both** branches of ADR-0026 and cannot choose
between them -- equilibrium recombination holds the plume at 16 224 K at the nozzle exit, frozen
recombination drops it to 5 297 K. That factor of 3 spans the whole conductivity cliff, and it is
the largest single uncertainty left on the nozzle side.

Choosing is a race between two clocks. `tau_rec` is how long the plasma takes to put itself back
together; `tau_exp` is how long the gas takes to thin out. Their ratio is a Damkohler number, and
freeze-out is where it passes 1 -- the classical Bray criterion.

The rate coefficients are literature values, and the tests below pin them to published literals
rather than to the implementation.
"""

from __future__ import annotations

import math

import pytest

from puffsat import eos_water, recombination


def test_radiative_recombination_matches_the_published_case_b_coefficient() -> None:
    """Slice 1: the anchor literal every nebular and plasma text carries.

    Hydrogenic case-B radiative recombination at 10^4 K is `2.59e-13 cm^3/s`
    (Osterbrock & Ferland, *AGN2*, table 2.1) -- that is `2.59e-19 m^3/s` in SI. Case B rather
    than case A because recombinations straight to the ground state emit a photon that is
    immediately reabsorbed by a neighbour, so they return no energy and must not be counted.
    """
    assert recombination.radiative_recombination(1.0e4) == pytest.approx(2.59e-19, rel=0.02)


def test_radiative_recombination_falls_with_temperature_at_the_published_power() -> None:
    """Slice 1b: the scaling, checked as a ratio so the anchor value cancels out.

    Case B goes as `T^-0.75` over the range that matters here, so quadrupling the temperature
    must multiply the coefficient by `4^-0.75 = 0.3536`. Hot plasma recombines *reluctantly*:
    a fast electron flies past the ion instead of being captured.
    """
    ratio = recombination.radiative_recombination(4.0e4) / recombination.radiative_recombination(
        1.0e4
    )
    assert ratio == pytest.approx(4.0**-0.75, rel=0.02)


def test_three_body_recombination_matches_the_published_coefficient_in_ev() -> None:
    """Slice 2: `X+ + e + e -> X + e`, where the second electron carries off the binding energy.

    The published coefficient is `8.75e-27 T_e^-4.5 cm^6/s` with **`T_e` in electronvolts**
    (Zel'dovich & Raizer, *Physics of Shock Waves and High-Temperature Hydrodynamic Phenomena*).
    Reading that `T_e` as kelvin instead of eV changes the rate by ~1e11 and silently inverts the
    verdict, so the unit is pinned here by an independent order-of-magnitude argument as well as
    by the literal.

    **The independent argument (Thomson).** Capture happens when an electron comes within the
    Coulomb radius `b = e^2/(4 pi eps0 k T)`, which at 1 eV is 1.44 nm. A three-body rate
    coefficient has units of `m^6/s`, i.e. `(cross-section)(volume)(speed) ~ b^5 v_e`. With
    `v_e ~ 5.9e5 m/s` at 1 eV that is `(1.44e-9)^5 * 5.9e5 = 3.6e-39 m^6/s` -- the same few times
    `1e-39 m^6/s` the eV reading gives, and nowhere near the kelvin reading.
    """
    k_at_1ev = recombination.three_body_coefficient(1.0 * recombination.EV_IN_KELVIN)
    assert k_at_1ev == pytest.approx(8.75e-39, rel=0.02)
    assert 1.0e-39 < k_at_1ev < 1.0e-38  # the Thomson estimate, order of magnitude


def test_three_body_recombination_is_linear_in_the_third_body() -> None:
    """Slice 2b: it needs a spectator electron to carry the energy away, so the effective
    two-body rate is proportional to how many spectators there are. That linearity is what makes
    it dominate in a dense plasma and vanish in a thin one."""
    at_1e24 = recombination.three_body_recombination(8000.0, 1.0e24)
    at_2e24 = recombination.three_body_recombination(8000.0, 2.0e24)
    assert at_2e24 == pytest.approx(2.0 * at_1e24, rel=1e-9)


def test_three_body_dominates_radiative_at_plume_conditions() -> None:
    """Slice 2c: which channel actually returns the ionisation store here.

    Textbook expectation is that three-body wins in dense, cool plasma and radiative wins in thin,
    hot plasma. The plume sits at ~1e25 electrons/m^3, which is dense -- about a fifth of
    atmospheric number density -- so three-body should dominate by orders of magnitude. If it did
    not, the whole freeze question would be decided by the (much slower) radiative channel.
    """
    temp, n_e = 20000.0, 1.0e25
    assert recombination.three_body_recombination(
        temp, n_e
    ) > 100.0 * recombination.radiative_recombination(temp)


def test_expansion_time_is_the_e_folding_time_of_density() -> None:
    """Slice 3: the clock recombination is racing against.

    The expansion timescale is `rho / |d rho/dt|`, i.e. the time for density to fall by a factor
    `e`. Hand cases: density falling by exactly `e` over 1 ms gives 1 ms, and falling by 4 over
    2 ms gives `2/ln 4 = 1.4427 ms` -- the log matters, and using the naive `2 ms` would overstate
    the time available for chemistry by 39%.
    """
    assert recombination.expansion_time(0.0, 1.0, 1.0e-3, 1.0 / math.e) == pytest.approx(1.0e-3)
    assert recombination.expansion_time(0.0, 0.32, 2.0e-3, 0.08) == pytest.approx(
        2.0e-3 / math.log(4.0)
    )


def test_damkohler_verdict_splits_at_one() -> None:
    """Slice 3b: the Bray criterion. `Da = tau_exp / tau_rec`, and the verdict is which branch of
    ADR-0026's bracket the gas actually follows -- not a matter of taste."""
    assert recombination.verdict(1.0e3) == "equilibrium"
    assert recombination.verdict(1.0e-3) == "frozen"
    assert recombination.verdict(1.0) == "freezing"


def test_atom_recombination_uses_the_three_body_water_reformation_rate() -> None:
    """Slice 4: the dissociation store, which is larger than the ionisation store and returns by a
    different and much slower route.

    Re-forming water needs a three-body collision, `H + OH + M -> H2O + M`. Baulch et al.'s
    evaluation gives `k0 = 6.1e-26 T^-2 cm^6 molecule^-2 s^-1`; `1 cm^6 = 1e-12 m^6`, so at
    3000 K that is `6.1e-38 / 3000^2 = 6.78e-45 m^6/s`.
    """
    assert recombination.atom_three_body_coefficient(3000.0) == pytest.approx(6.78e-45, rel=0.02)
    # T^-2, so doubling the temperature quarters it.
    ratio = recombination.atom_three_body_coefficient(
        6000.0
    ) / recombination.atom_three_body_coefficient(3000.0)
    assert ratio == pytest.approx(0.25, rel=1e-9)


def test_atom_recombination_is_far_slower_than_ionic_recombination() -> None:
    """Slice 4b: the ordering that decides which store freezes first.

    The ionic channel has a spectator *electron* as third body and a huge Coulomb cross-section;
    the atomic channel needs a neutral third body and a much smaller one. So the dissociation
    store -- the larger of the two -- is the one at risk, and a study that only checked the
    ionisation freeze would be checking the easy half.
    """
    temp, n_e, n_third, n_atom = 5000.0, 5.0e21, 1.0e25, 5.0e24
    ionic = 1.0 / (recombination.three_body_recombination(temp, n_e) * n_e)
    atomic = recombination.atom_recombination_time(temp, n_third, n_atom)
    assert atomic > 10.0 * ionic  # 79x at this state; an order of magnitude is the claim


def test_freeze_station_flips_to_frozen_when_the_expansion_outruns_the_chemistry() -> None:
    """Slice 5: the verdict at one station, and the control that shows it can say "frozen".

    A criterion that only ever returns "equilibrium" would be worthless, so the same state is run
    at two expansion rates. The plume's own rate keeps up; the same gas blown out a million times
    faster does not. The margin is reported in decades of rate coefficient, because that is the
    uncertainty that actually threatens the verdict -- evaluated three-body rates carry a factor
    of a few, so a verdict with less than ~0.5 decades of margin is not a verdict.
    """
    slow = recombination.freeze_station(
        time=1.0e-3,
        temp=6000.0,
        rho=0.1,
        n_e=5.0e21,
        n_h=5.0e24,
        n_oh=5.0e24,
        n_third=1.0e25,
        tau_expansion=1.4e-3,
    )
    assert slow.verdict_dissociation_h_limited == "equilibrium"
    assert slow.margin_decades > 0.5

    fast = recombination.freeze_station(
        time=1.0e-9,
        temp=6000.0,
        rho=0.1,
        n_e=5.0e21,
        n_h=5.0e24,
        n_oh=5.0e24,
        n_third=1.0e25,
        tau_expansion=1.4e-9,
    )
    assert fast.verdict_dissociation_h_limited == "frozen"
    assert fast.margin_decades < 0.0

    # Da is a pure ratio, so a millionfold faster expansion is a millionfold smaller Da.
    assert fast.da_dissociation_h_limited == pytest.approx(
        slow.da_dissociation_h_limited * 1e-6, rel=1e-9
    )


def test_binding_damkohler_ignores_stations_whose_store_is_already_spent() -> None:
    """Slice 6: `min Da` is the wrong statistic, and using it would have inverted this answer.

    A store that has already been returned cannot freeze. On the coldest leg the ionisation
    Damkohler number falls to 0.067 -- but at a station holding 0.01% of the reservoir store,
    where freezing it releases nothing. The number that governs the energy release is the
    *smallest* `Da` among stations that still hold a meaningful fraction of the store.

    Here: three stations, the last one starved. `min Da` says 0.05 (frozen); the binding value
    ignores it and says 500 (equilibrium), which is what the physics supports.
    """
    stations = [(5.0e3, 0.80), (5.0e2, 0.20), (5.0e-2, 0.0001)]
    assert recombination.binding_damkohler(stations, threshold=0.01) == pytest.approx(500.0)
    assert min(da for da, _ in stations) == pytest.approx(0.05)


def test_binding_damkohler_needs_at_least_one_station_holding_the_store() -> None:
    """Slice 6b: if nothing holds the store the question is empty, and saying "equilibrium"
    would be a false pass. It has to refuse rather than return a default."""
    with pytest.raises(ValueError):
        recombination.binding_damkohler([(1.0e5, 1e-9)], threshold=0.01)


# --- The OH bracket (2026-08-25) ----------------------------------------------------------------
#
# `H + OH + M -> H2O + M` was rated with `n_atom = n_H`, because the species set had no OH and the
# code said so: an *over*estimate of the rate wherever OH is scarce. It is scarce -- equilibrium OH
# is 2e-5 of n_H on the hot legs and 4e-2 on the cold one -- so the shipped rate was 30x to 50000x
# too fast. Now that both densities exist the race is run at both edges and the answer is a bracket.


def test_the_two_edges_of_the_bracket_are_the_two_partner_densities() -> None:
    """`Da` is linear in the partner density, so the bracket's width *is* the `n_OH/n_H` ratio.

    The H-limited edge assumes an OH is always waiting (OH formation never limits); the OH-limited
    edge assumes OH never rises above its equilibrium value. Neither is the truth, which needs a
    reaction network -- what the bracket buys is that the truth is between two computed numbers
    instead of past one.
    """
    station = recombination.freeze_station(
        time=1.0e-3,
        temp=6000.0,
        rho=0.1,
        n_e=5.0e21,
        n_h=5.0e24,
        n_oh=5.0e20,
        n_third=1.0e25,
        tau_expansion=1.4e-3,
    )
    ratio = 5.0e20 / 5.0e24
    assert station.da_dissociation_oh_limited == pytest.approx(
        station.da_dissociation_h_limited * ratio, rel=1e-9
    )
    assert station.tau_dissociation_oh_limited == pytest.approx(
        station.tau_dissociation_h_limited / ratio, rel=1e-9
    )


def test_the_conservative_edge_can_freeze_where_the_optimistic_one_does_not() -> None:
    """The bracket has to be able to straddle the verdict, or it is not doing any work.

    Same state, same expansion: with an OH always available the chemistry keeps up; with only the
    equilibrium OH it does not. That is the entire content of the finding, at one station.
    """
    station = recombination.freeze_station(
        time=1.0e-3,
        temp=6000.0,
        rho=0.1,
        n_e=5.0e21,
        n_h=5.0e24,
        n_oh=5.0e19,
        n_third=1.0e25,
        tau_expansion=1.4e-3,
    )
    assert station.verdict_dissociation_h_limited == "equilibrium"
    assert station.verdict_dissociation_oh_limited == "frozen"


def test_the_held_bond_energy_is_not_the_dissociated_fraction() -> None:
    """The store is what no bond is paying for, and OH pays for half of one.

    Charging `(1 - n_H2O/n_f) * D_AT` was right when the only alternatives to H2O were free atoms.
    With the intermediates it overcharges: at 0.32 kg/m^3 and 3000 K the old measure calls the gas
    14.4% dissociated and would strand 7.3 MJ/kg, while the bonds actually missing are worth
    2.5 MJ/kg -- a factor of 3.
    """
    rho, temp = 0.32, 3000.0
    held = eos_water.bond_energy_held(rho, temp)
    n_f = rho / eos_water.M_H2O
    old_measure = (
        (1.0 - eos_water.composition(rho, temp).n_h2o / n_f) * eos_water.D_AT / (eos_water.M_H2O)
    )
    assert held < 0.4 * old_measure
    assert held / eos_water.FULL_ATOMIZATION_ENERGY == pytest.approx(0.049, abs=0.005)


def test_a_fully_atomised_gas_still_strands_the_whole_ceiling() -> None:
    """The generalisation must not move the case it generalises: where nothing is molecular, the
    held bond energy is the full atomization energy and the old form was right."""
    rho, temp = 2.354e-2, 16_062.0
    assert eos_water.bond_energy_held(rho, temp) == pytest.approx(
        eos_water.FULL_ATOMIZATION_ENERGY, rel=1e-3
    )
