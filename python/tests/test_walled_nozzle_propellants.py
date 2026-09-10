"""N10 item 3: the propellant ladder.

The strongest tests here are the two *external* ones -- the ask's own stated slug ratios for water
and hydrogen, reproduced from a fixed point that was never told about them -- and the scaling law,
which says the whole ladder is one variable in disguise.
"""

from __future__ import annotations

import math
from itertools import pairwise

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
    assert all(a.slug_ratio < b.slug_ratio for a, b in pairwise(ordered))
    assert all(a.isp_effective > b.isp_effective for a, b in pairwise(ordered))


def test_the_scaling_law_is_a_property_of_the_gross_convention(
    rungs: list[propellants.Rung],
) -> None:
    """`Isp ~ 1/sqrt(m_bar)` to within 10% across a factor of six in particle mass -- **on the
    credit-only column**. This is the rocket equation's molecular-weight law reappearing in a
    chamber heated by an impact rather than by combustion, which is why the ordering is not
    negotiable by nozzle design.

    **The corrected column does not obey it, and that is a result rather than a defect.** The
    head-on debit `w/(k g0)` also runs as `1/k`, so it falls hardest on exactly the light fluids
    the law rewards, and the ladder comes out compressed. Asserting the law on `isp_effective`
    would be asserting that the debit does not exist.
    """
    methane = next(r for r in rungs if r.name == "methane")
    for r in rungs:
        predicted = math.sqrt(methane.mean_atomised_mass / r.mean_atomised_mass)
        gross = r.isp_carried_gross / methane.isp_carried_gross
        assert gross == pytest.approx(predicted, rel=0.10)
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    corrected = hydrogen.isp_effective / methane.isp_effective
    assert corrected < hydrogen.isp_carried_gross / methane.isp_carried_gross
    assert corrected > 1.0


def test_water_is_below_methane_and_hydrogen_far_above_it(
    rungs: list[propellants.Rung],
) -> None:
    """The solved part of the ladder. Water is *below* methane and hydrogen is half again above
    it, both on their own EOS -- so the carbon exposure of N10 item 5 cannot be escaped by
    switching to water without paying for it in specific impulse.

    **The margins are the corrected ones and they are narrower than the credit-only column
    showed**: hydrogen leads methane by 1.55x rather than 1.88x, and methane leads water by 1.23x
    rather than 1.42x. The debit compresses the ladder because it punishes a small `k`. The
    *ordering* is untouched, which is the part the verdict rests on.

    **Ammonia is deliberately not asserted against methane here**: its conversion fraction is
    assumed rather than solved, and `test_ammonia_ordering_is_not_established` records why that
    makes the comparison unsafe in the direction that matters.
    """
    methane = next(r for r in rungs if r.name == "methane")
    water = next(r for r in rungs if r.name == "water")
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    assert water.isp_effective < methane.isp_effective
    assert hydrogen.isp_effective > 1.5 * methane.isp_effective
    assert methane.isp_effective / water.isp_effective == pytest.approx(1.23, abs=0.02)
    # **Ammonia no longer clears water at the bottom of its assumed conversion range.** It did on
    # the credit-only column. Ammonia's `k` is 28.5 against water's 39.5, so it pays a bigger
    # debit, and at 0.35 conversion that is enough to put it just under -- 594 s against 600 s.
    assert propellants.ammonia_rung(0.35).isp_effective < water.isp_effective
    assert propellants.ammonia_rung(0.35).isp_carried_gross > water.isp_carried_gross
    for conv in (0.406, 0.45):
        assert propellants.ammonia_rung(conv).isp_effective > water.isp_effective


def test_ammonia_ordering_is_not_established(rungs: list[propellants.Rung]) -> None:
    """**Ammonia only sits below methane because it was handed methane's conversion fraction.**

    Methane converts just 0.406 for a reason W9 explains: its carbon is parked in C3 at the exit
    with only 30% of atomisation returned. Ammonia has no C3 trap -- its store comes back as N2
    (941 kJ/mol) and H2 -- so a higher conversion is physically plausible.

    **The correction moves the crossover a long way, and towards ammonia.** On the credit-only
    column ammonia had to reach 0.601 to match methane -- a factor 1.48 above methane's own, a
    real stretch. With the head-on debit charged it needs only **0.477**, a factor 1.18, which
    hydrogen already beats at this same area ratio. The reason is that inverting the corrected form
    adds back the debit before dividing: `u_e = (Isp_eff k g0 + w)/(1+k)`, and ammonia's leaner
    charge pays more of that debit than methane does, so less exhaust speed is asked of it.

    This test pins the crossover so nobody reads the W8 table as settling ammonia -- and after the
    correction it settles it even less than before.
    """
    methane = next(r for r in rungs if r.name == "methane")
    u = propellants.ammonia_energy()
    k = chamber.SPECIFIC_PULSE_ENERGY / u - 1.0
    speed = (methane.isp_effective * k * propellants.G0 + chamber.CLOSING_SPEED) / (1.0 + k)
    needed = speed**2 / (2.0 * u)
    assert needed == pytest.approx(0.477, abs=0.01)
    assert needed > methane.conversion
    # The credit-only crossover, kept so the size of the shift is on the record.
    gross_speed = methane.isp_carried_gross * propellants.G0 * k / (1.0 + k)
    assert gross_speed**2 / (2.0 * u) == pytest.approx(0.601, abs=0.01)
    # Not out of reach: a solved rung in this very table beats it on conversion... nearly.
    hydrogen = next(r for r in rungs if r.name == "hydrogen")
    assert hydrogen.conversion > 0.5


def test_the_credit_exceeds_total_but_the_corrected_figure_does_not(
    rungs: list[propellants.Rung],
) -> None:
    """**The correction that reverses the sign of a whole column.**

    `(1+k)/k` is where `k` enters a second time, and on its own it *widens* the gap the fixed
    point opened: 1.13x for hydrogen against 1.03x for water. But `k` enters a third time through
    the head-on debit `w/(k g0)`, and since the credit is only `u_e/(k g0)`, the debit wins
    wherever `u_e < w` -- which is every rung in this study, by a factor of four or more.

    So the corrected figure of merit sits **below** the real Isp, not above it. A reader who takes
    the credit and forgets the debit reads a number too high by `(w - u_e)/(k g0)`.
    """
    for r in rungs:
        assert r.isp_carried_gross > r.isp_total
        assert r.isp_effective < r.isp_total
        assert r.exhaust_speed < chamber.CLOSING_SPEED
        assert r.isp_carried_gross - r.isp_effective == pytest.approx(r.momentum_debit)
    ordered = sorted(rungs, key=lambda r: r.slug_ratio)
    boosts = [r.isp_carried_gross / r.isp_total for r in ordered]
    assert all(a > b for a, b in pairwise(boosts))
    # The debit runs the same way as the credit, which is why it cannot be waved off as small.
    debits = [r.momentum_debit for r in ordered]
    assert all(a > b for a, b in pairwise(debits))


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


def test_the_free_impactor_bonus_rewards_the_smallest_slug(
    rungs: list[propellants.Rung],
) -> None:
    """**Why the effective convention favours hydrogen, stated as its own result.**

    The vehicle carries the slug; the impactor arrives from outside at 75 km/s and was never
    lifted. So the impulse per kilogram *carried* exceeds the impulse per kilogram *expelled* by
    `(1+k)/k = 1 + 1/k` -- a bonus that grows as the required slug shrinks. Hydrogen soaks the
    pulse up in 194 kg and collects +12.9%; water needs 989 kg and collects +2.5%.

    This is a second advantage for light propellants, independent of the chemistry that set `k`
    in the first place, and it is why W8 prints both columns rather than picking one.

    **It is only half the mass ledger.** The head-on momentum debit is also `1/k`-shaped and is
    larger, so this bonus never survives into the figure of merit intact.
    """
    by_k = sorted(rungs, key=lambda r: r.slug_ratio)
    assert [r.name for r in by_k] == ["hydrogen", "methane", "water"]
    bonuses = [r.free_impactor_bonus for r in by_k]
    assert all(a > b for a, b in pairwise(bonuses))
    assert by_k[0].free_impactor_bonus == pytest.approx(0.129, abs=0.003)
    assert by_k[-1].free_impactor_bonus == pytest.approx(0.025, abs=0.003)
    for r in rungs:
        assert r.isp_carried_gross / r.isp_true == pytest.approx(1.0 + r.free_impactor_bonus)
        assert r.carried_slug_mass == pytest.approx(r.slug_ratio * chamber.IMPACTOR_MASS)
        # And the bonus is not the whole mass ledger: the debit runs as 1/k too, and is bigger.
        assert r.momentum_debit > r.free_impactor_bonus * r.isp_true
    assert by_k[0].carried_slug_mass == pytest.approx(194.2, abs=1.0)


# ---- The head-on momentum debit ---------------------------------------------------------------


def test_the_ladder_uses_the_companion_head_on_form(rungs: list[propellants.Rung]) -> None:
    """`Isp_eff = w (eta_jet sqrt(1+k) - 1) / (k g0)` -- the companion's own head-on expression,
    and the same `sqrt(1+k) - 1` the tamper study's `beta_ideal` carries at `eta_jet = 1`."""
    for r in rungs:
        paper = (
            chamber.CLOSING_SPEED
            * (r.eta_jet * math.sqrt(1.0 + r.slug_ratio) - 1.0)
            / (r.slug_ratio * propellants.G0)
        )
        assert r.isp_effective == pytest.approx(paper, rel=1e-9)


def test_eta_jet_is_the_square_root_of_the_conversion(rungs: list[propellants.Rung]) -> None:
    """The bridge to the companion's convention: `(1+k) u = w^2/2` makes `eta_jet^2` the
    conversion fraction, so a conversion quoted *as* an `eta_jet` is wrong by a square root."""
    for r in rungs:
        assert r.eta_jet == pytest.approx(math.sqrt(r.conversion), rel=1e-6)


def test_the_debit_is_one_projectiles_momentum(rungs: list[propellants.Rung]) -> None:
    """`w/(k g0)`: the same momentum for every rung, spread over whatever slug was carried."""
    for r in rungs:
        assert r.momentum_debit * r.slug_ratio * propellants.G0 == pytest.approx(
            chamber.CLOSING_SPEED, rel=1e-9
        )


def test_every_rung_clears_the_thrust_floor(rungs: list[propellants.Rung]) -> None:
    """Below `eta_jet = 1/sqrt(1+k)` the burn is a brake. The ladder clears it, but on this
    convention that has to be checked rather than assumed -- it is what `isp_effective > 0` means.
    """
    for r in rungs:
        assert r.eta_jet > r.isp.thrust_floor
        assert not r.isp.pushes_backwards
        assert r.isp_effective > 0.0


def test_the_ordering_survives_the_debit(rungs: list[propellants.Rung]) -> None:
    """**The verdict that had to be re-checked, not assumed.** The debit punishes small `k`, which
    is hydrogen's whole advantage, so the correction could in principle have reordered the ladder.
    It does not: it compresses hydrogen's lead over methane from 1.88x to 1.55x."""
    ordered = sorted(rungs, key=lambda r: r.slug_ratio)
    assert [r.name for r in ordered] == ["hydrogen", "methane", "water"]
    assert all(a.isp_effective > b.isp_effective for a, b in pairwise(ordered))
    hydrogen, methane = ordered[0], ordered[1]
    assert hydrogen.isp_carried_gross / methane.isp_carried_gross == pytest.approx(1.88, abs=0.02)
    assert hydrogen.isp_effective / methane.isp_effective == pytest.approx(1.55, abs=0.02)
