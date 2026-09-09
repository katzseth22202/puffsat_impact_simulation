"""N10 item 3: the propellant ladder.

The strongest tests here are the two *external* ones -- the ask's own stated slug ratios for water
and hydrogen, reproduced from a fixed point that was never told about them -- and the scaling law,
which says the whole ladder is one variable in disguise.
"""

from __future__ import annotations

import math

import pytest

from puffsat.walled_nozzle import chamber, propellants


@pytest.fixture(scope="module")
def rungs() -> list[propellants.Rung]:
    return propellants.ladder()


def test_the_asks_own_slug_ratios_are_reproduced(rungs: list[propellants.Rung]) -> None:
    """**The cross-check that makes the rest of the ladder credible.** N10 item 3 states
    `k = 37.70` for water and `k = 7.99` for pure hydrogen. Nothing in the fixed point was told
    either number -- it solves `(1+k) e(rho, T) = w^2/2` on each fluid's own EOS -- and it lands
    within 5% of both."""
    water = next(r for r in rungs if r.name == "water")
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    assert water.slug_ratio == pytest.approx(37.70, rel=0.06)
    assert hydrogen.slug_ratio == pytest.approx(7.99, rel=0.06)


def test_the_ladder_sorts_by_mean_atomised_particle_mass(
    rungs: list[propellants.Rung],
) -> None:
    """Lighter stores more per kilogram, so it needs less mass to absorb the pulse, so `k` falls
    and the exhaust speed rises. One variable explains the whole ordering."""
    ordered = sorted(rungs, key=lambda r: r.mean_atomised_mass)
    assert [r.name for r in ordered] == ["hydrogen", "methane", "water"]
    assert all(a.slug_ratio < b.slug_ratio for a, b in zip(ordered, ordered[1:], strict=False))
    assert all(
        a.isp_effective > b.isp_effective for a, b in zip(ordered, ordered[1:], strict=False)
    )


def test_effective_isp_tracks_one_over_root_mean_mass(rungs: list[propellants.Rung]) -> None:
    """`Isp ~ 1/sqrt(m_bar)` to within 10% across a factor of six in particle mass. This is the
    rocket equation's molecular-weight law reappearing in a chamber that is heated by an impact
    rather than by combustion, which is why the ordering is not negotiable by nozzle design."""
    methane = next(r for r in rungs if r.name == "methane")
    for r in rungs:
        predicted = math.sqrt(methane.mean_atomised_mass / r.mean_atomised_mass)
        assert r.isp_effective / methane.isp_effective == pytest.approx(predicted, rel=0.10)


def test_water_is_below_methane_and_hydrogen_far_above_it(
    rungs: list[propellants.Rung],
) -> None:
    """The solved part of the ladder. Water is *below* methane and hydrogen is nearly double it,
    both on their own EOS -- so the carbon exposure of N10 item 5 cannot be escaped by switching
    to water without paying for it in specific impulse.

    **Ammonia is deliberately not asserted against methane here**: its conversion fraction is
    assumed rather than solved, and `test_ammonia_ordering_is_not_established` records why that
    makes the comparison unsafe in the direction that matters.
    """
    methane = next(r for r in rungs if r.name == "methane")
    water = next(r for r in rungs if r.name == "water")
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    assert water.isp_effective < methane.isp_effective
    assert hydrogen.isp_effective > 1.7 * methane.isp_effective
    # Ammonia clears water at every plausible conversion; that much is safe.
    for conv in (0.35, 0.406, 0.45):
        assert propellants.ammonia_rung(conv).isp_effective > water.isp_effective


def test_ammonia_ordering_is_not_established(rungs: list[propellants.Rung]) -> None:
    """**Ammonia only sits below methane because it was handed methane's conversion fraction.**

    Methane converts just 0.406 for a reason W9 explains: its carbon is parked in C3 at the exit
    with only 30% of atomisation returned. Ammonia has no C3 trap -- its store comes back as N2
    (941 kJ/mol) and H2 -- so a higher conversion is physically plausible. The crossover is 0.601,
    a factor of 1.48 above methane, and hydrogen already reaches 0.532 at this same area ratio,
    so that is a stretch rather than an impossibility.

    This test pins the crossover so nobody reads the W8 table as settling ammonia.
    """
    methane = next(r for r in rungs if r.name == "methane")
    u = propellants.ammonia_energy()
    k = chamber.SPECIFIC_PULSE_ENERGY / u - 1.0
    needed = (methane.isp_effective * propellants.G0 * k / (1.0 + k)) ** 2 / (2.0 * u)
    assert needed == pytest.approx(0.601, abs=0.01)
    assert needed > methane.conversion
    # Not out of reach: a solved rung in this very table beats it on conversion... nearly.
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    assert hydrogen.conversion > 0.5


def test_effective_exceeds_total_and_rewards_a_small_slug_ratio(
    rungs: list[propellants.Rung],
) -> None:
    """`(1+k)/k` is where `k` enters a second time. It is 1.13 for hydrogen and 1.03 for water, so
    the effective convention *widens* the gap the fixed point already opened."""
    for r in rungs:
        assert r.isp_effective > r.isp_total
    ordered = sorted(rungs, key=lambda r: r.slug_ratio)
    boosts = [r.isp_effective / r.isp_total for r in ordered]
    assert all(a > b for a, b in zip(boosts, boosts[1:], strict=False))


def test_conversion_is_a_fraction_and_below_the_ideal_speed(
    rungs: list[propellants.Rung],
) -> None:
    """Sanity: no nozzle beats full conversion, and `w/sqrt(1+k)` is that ceiling."""
    for r in rungs:
        assert 0.0 < r.conversion < 1.0
        assert r.exhaust_speed < r.ideal_exhaust_speed
        assert r.exhaust_speed == pytest.approx(
            r.ideal_exhaust_speed * math.sqrt(r.conversion), rel=1e-6
        )


def test_the_ammonia_rung_is_flagged_and_the_solved_ones_are_not() -> None:
    """Ammonia is assembled from the ask's atomisation energy plus an exact translational term,
    not solved on an EOS. The flag is what keeps it out of a table of solved results."""
    assert propellants.ammonia_rung(0.4).estimated
    assert not any(r.estimated for r in propellants.ladder())


def test_the_ionisation_allowance_is_small_enough_not_to_matter() -> None:
    """A previous estimate assumed 8 MJ/kg of ionisation for ammonia. The solved methane chamber
    puts it at 0.45% of `u`, i.e. ~0.6 MJ/kg, and the assumption moved `k` by ~8%. This pins the
    allowance small so that error cannot return."""
    assert propellants.IONISATION_SHARE < 0.02
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    from puffsat import eos_methane

    held = state.store_charged * eos_methane.FULL_ATOMIZATION_ENERGY
    sensible = state.sensible_fraction * state.energy
    solved_share = (state.energy - held - sensible) / state.energy
    assert solved_share < propellants.IONISATION_SHARE


# ---- W9: the acetylene route -------------------------------------------------------------------


def test_acetylene_returns_most_of_the_carbon_store_without_condensing() -> None:
    """**The finding that reframes N10 item 5.** `2 CH4 -> C2H2 + 3 H2` is stoichiometrically
    exact and needs no rate constant: acetylene is nearly as strongly bound per atom as methane
    was, so reaching it returns 89% of full atomisation against 52.6% for the H2 channel alone.

    That is 84% of the 43.3% carbon store recovered **in the gas phase**, leaving only ~7 points
    of atomisation -- about 5% of the chamber energy budget -- genuinely hostage to soot
    nucleation. It also removes most of the two-phase lag exposure, because carbon leaving as
    C2H2 is gas and stays momentum-coupled.
    """
    from puffsat import eos_methane as m

    acetylene_route = (m.C2H2.d_at + 3.0 * m.H2.d_at) / (2.0 * m.CH4.d_at)
    h2_only = 4.0 * m.H2.d_at / (2.0 * m.CH4.d_at)
    assert acetylene_route == pytest.approx(0.890, abs=0.005)
    assert h2_only == pytest.approx(0.526, abs=0.005)
    # The carbon store is 43.3% of atomisation; acetylene recovers most of it.
    recovered = (acetylene_route - h2_only) / (
        m.CARBON_CONDENSATION_ENERGY / m.FULL_ATOMIZATION_ENERGY
    )
    assert 0.80 < recovered < 0.88


def test_the_equilibrium_carbon_path_runs_through_c3_to_acetylene() -> None:
    """C -> C3 -> C2H2 as the gas cools, and `held` falls to the 89% return once acetylene wins.

    Worth pinning because **C3 is this EOS's most exposed species** (its quasilinear bend is
    treated harmonically), and it is the one carrying the carbon through 4000-6000 K.
    """
    from puffsat import eos_methane as m

    rho = 1.0
    frac = {}
    for temp in (6000.0, 4000.0, 3000.0):
        c = m.composition(rho, temp)
        frac[temp] = (
            c.density_of(m.C3) * 3.0 / c.n_c_nuclei,
            c.density_of(m.C2H2) * 2.0 / c.n_c_nuclei,
        )
    assert frac[6000.0][0] > 0.5  # C3 dominates hot
    assert frac[3000.0][1] > 0.95  # acetylene dominates cold
    assert frac[4000.0][1] > frac[6000.0][1]
    held = m.bond_energy_held(rho, 3000.0) / m.FULL_ATOMIZATION_ENERGY
    assert held == pytest.approx(1.0 - 0.890, abs=0.02)


def test_the_flown_nozzle_leaves_the_acetylene_energy_on_the_table() -> None:
    """**W9 arriving at W6's conclusion from the carbon side.** At the ask's 7 m^2 throat the
    exit is 5584 K and the carbon is still ~99% free atoms, so essentially none of the acetylene
    energy is collected. The limit is expansion ratio, not nucleation and not rate constants."""
    from puffsat import eos_methane as m
    from puffsat.walled_nozzle import freeze

    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    exit_station = freeze.history(state, throat_area=7.0, steps=64)[-1]
    c = m.composition(exit_station.rho, exit_station.temp)
    assert exit_station.temp > 5000.0
    # Almost no acetylene: the carbon has not reached the species that banks the store.
    assert c.density_of(m.C2H2) * 2.0 / c.n_c_nuclei < 0.05
    # It is parked in C3 instead -- *not* in free atoms, which is a better place to be but still
    # 59 points of atomisation short of the acetylene equilibrium.
    assert c.density_of(m.C3) * 3.0 / c.n_c_nuclei > 0.5
    assert c.carbon_free < 0.35
    assert exit_station.store_held > 0.6
