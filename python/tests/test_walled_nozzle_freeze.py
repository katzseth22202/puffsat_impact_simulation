"""N10.1-3: the rate literature and the freeze race.

The rate constants are literature values, so the tests that matter are not "does the arithmetic
work" but **"is the number we typed the number the source published, and does it still say what
we claim it says"**. Hence the cross-evaluation agreement tests: two independent evaluations of
the same channel landing within a few percent is the evidence that the channel is really known,
and if someone edits a constant that agreement is what breaks first.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import eos_methane
from puffsat.walled_nozzle import chamber, freeze, rates

# ---- The literature constants ------------------------------------------------------------------


def test_published_values_are_carried_in_si_without_a_lost_factor() -> None:
    """`1 cm^6 = 1e-12 m^6`. Baulch publishes 8.85e-33 cm^6 molecule^-2 s^-1."""
    assert rates.H_H2(rates.T_REFERENCE) == pytest.approx(8.85e-45)
    assert rates.N_N2(rates.T_REFERENCE) == pytest.approx(4.16e-45)


def test_two_independent_evaluations_of_the_same_channel_agree() -> None:
    """**The reason `H + H + M` can be trusted where the carbon channels cannot.**

    Baulch et al. 1992 and Cohen & Westberg 1983 evaluated `M = H2` separately and landed within
    2% with the same exponent. That is a far stronger statement than either uncertainty bar.
    """
    assert rates.H_H2.n == rates.H_H2_COHEN.n
    ratio = rates.H_H2_COHEN(3000.0) / rates.H_H2(3000.0)
    assert 0.95 < ratio < 1.05


def test_stated_uncertainty_converts_to_decades() -> None:
    """A factor of 3.16 is exactly half a decade -- which is how it enters the margin."""
    assert rates.H_H2.decades() == pytest.approx(0.5, abs=0.005)
    assert rates.H_H2_COHEN.decades() < rates.H_H2.decades()
    # Silence must not buy headroom.
    assert rates.H_H2_HOT.decades() == 0.0


def test_atomic_hydrogen_is_the_best_third_body_in_every_study() -> None:
    """The load-bearing find: `k_H > k_H2 > k_Ar`, in all three shock tubes, by an order of
    magnitude at the top. At 10 kK the charge is 93-98% dissociated, so this is the case that
    actually applies -- and a margin computed on H2 or an inert understates the rate."""
    for name, h_over_ar, h2_over_ar in rates.THIRD_BODY_EFFICIENCY:
        assert h_over_ar > h2_over_ar > 1.0, name
        assert h_over_ar > 6.0, name


def test_the_hot_measurement_is_faster_than_the_extrapolated_cold_fit() -> None:
    """Hurle et al. measured 2500-7000 K directly and found no temperature dependence.

    Extrapolating Baulch's `T^-0.6` up to 5000 K therefore *understates* the rate by about 3x.
    This is why using the cold evaluations outside their window is conservative here rather than
    optimistic -- an important direction to have right when the result is a margin.
    """
    assert rates.H_H2_HOT(5000.0) > 2.0 * rates.H_H2(5000.0)
    assert not rates.H_H2_HOT.extrapolated(5000.0)
    assert rates.H_H2.extrapolated(6000.0)


def test_byron_needs_no_extrapolation_at_the_chamber_temperature() -> None:
    """Nitrogen's one high-temperature measurement straddles the flown 10 kK chamber."""
    assert rates.N_ATOMIC_N.t_min <= 8000.0 <= rates.N_ATOMIC_N.t_max


# ---- The mixture rule --------------------------------------------------------------------------


def test_mixture_coefficient_is_bracketed_by_its_own_third_bodies() -> None:
    """A weighted mean cannot escape the range of the things it averages."""
    comp = eos_methane.composition(1.0, 6000.0)
    k = freeze.mixture_coefficient(6000.0, comp)
    assert rates.H_AR_JACOBS(6000.0) <= k <= rates.H_ATOMIC_H(6000.0)


def test_a_fully_dissociated_gas_recombines_at_the_atomic_rate() -> None:
    """The limit that matters: at chamber conditions almost every third body *is* an H atom, so
    the mixture rule must collapse onto `H_ATOMIC_H` rather than onto an inert average."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    comp = eos_methane.composition(state.rho, state.temp)
    assert comp.hydrogen_dissociated > 0.9
    k = freeze.mixture_coefficient(state.temp, comp)
    assert k / rates.H_ATOMIC_H(state.temp) > 0.5


def test_recombination_time_is_quadratic_in_density() -> None:
    """`tau = 1/(k n_partner n_M)`: thin the gas 10x and the channel slows 100x. This quadratic
    sensitivity is the entire physical case for a walled throat over a magnetic one."""
    fast = freeze.recombination_time(1e-44, 1e24, 1e24)
    thin = freeze.recombination_time(1e-44, 1e23, 1e23)
    assert thin / fast == pytest.approx(100.0)
    assert freeze.recombination_time(1e-44, 0.0, 1e24) == math.inf


def test_verdict_bands() -> None:
    assert freeze.verdict(1e3) == "equilibrium"
    assert freeze.verdict(1.0) == "freezing"
    assert freeze.verdict(1e-3) == "frozen"


# ---- The race ----------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def flown() -> list[freeze.FreezeStation]:
    """The flown 200 m^3 chamber through the 7 m^2 throat, at a coarse but sufficient grid."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    return freeze.history(state, throat_area=7.0, steps=64)


def test_the_chemistry_never_freezes_inside_the_bore(flown: list[freeze.FreezeStation]) -> None:
    """**The headline answer to N10 items 1-3.** Every station in ADR-0016's nozzle keeps
    `Da` above the equilibrium threshold, so the equilibrium expansion is self-consistent and
    the frozen branch does not apply. Rubbia's 2700 s, which W2 showed *requires* recombination,
    is therefore not in tension with this geometry."""
    assert freeze.freeze_point(flown) is None
    assert min(s.damkohler for s in flown) > freeze.DA_EQUILIBRIUM
    assert all(s.verdict == "equilibrium" for s in flown)


def test_the_verdict_survives_the_rate_being_wrong_by_a_decade(
    flown: list[freeze.FreezeStation],
) -> None:
    """The margin, which is the form the answer had to take. The stated uncertainty on the
    evaluated coefficient is half a decade; the margin here is larger, so the verdict does not
    rest on the rate constant being right."""
    margin = min(s.margin_decades for s in flown)
    assert margin > rates.H_H2.decades()
    assert margin > 1.0


def test_the_store_is_returned_monotonically(flown: list[freeze.FreezeStation]) -> None:
    """Bond energy held can only fall along an equilibrium expansion: the gas is cooling, so
    atoms are pairing up, never the reverse."""
    held = [s.store_held for s in flown]
    assert all(b <= a + 1e-9 for a, b in pairwise(held))
    assert 0.0 < freeze.returned_fraction(flown) < 1.0


def test_refining_the_grid_does_not_invent_a_freeze_at_the_throat() -> None:
    """**Regression.** `cooling_history` interpolates exact stations onto the throat and the exit
    ratio, and at fine step counts one can land on top of its neighbour. That zero-width interval
    has `dt = 0`, so `tau_exp = 0` and `Da = 0` -- which reads as "frozen at the throat", the
    densest and fastest-recombining point in the nozzle. The artefact appears only when the grid
    is refined, so a coarse test would never see it."""
    state = chamber.solve_chamber(400.0, 14.0, chamber.FLOWN_TEMPERATURE)
    for steps in (64, 192, 384):
        st = freeze.history(state, throat_area=7.0, steps=steps)
        assert all(s.tau_expansion > 0.0 for s in st)
        assert freeze.freeze_point(st) is None, steps


def test_a_narrower_throat_expands_further_and_returns_more() -> None:
    """**Item 4's answer.** The throat is the knob, but not in the direction the ask expected: a
    narrower throat gives a larger area ratio, a colder exit, and more of the store handed back.
    The chemistry has margin to spare at both settings, so the binding constraint is the
    expansion ratio rather than the rate."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    narrow = freeze.history(state, throat_area=2.0, steps=64)
    wide = freeze.history(state, throat_area=7.0, steps=64)
    assert narrow[-1].temp < wide[-1].temp
    assert freeze.returned_fraction(narrow) > freeze.returned_fraction(wide)
    assert freeze.freeze_point(narrow) is None


# ---- The ammonia rung and the carbon non-answer -------------------------------------------------


def test_the_nitrogen_channel_is_uncertain_by_more_than_a_decade() -> None:
    """Byron's shock tube and the 2025 ab initio master equation disagree by ~1.3 decades. That
    spread, not either value, is the honest state of `N + N + M` -- so an ammonia claim that
    needs Byron's end is not supported."""
    cmp = freeze.nitrogen_comparison(6000.0)
    assert cmp["nitrogen_spread_decades"] > 1.0
    assert cmp["ratio_byron"] < 1.0 < cmp["ratio_abinitio"]


def test_the_asks_eight_times_slower_sits_inside_the_nitrogen_band() -> None:
    """The ask states "nitrogen is about eight times slower". It is inside the band, but only
    near the ab initio end; Byron's fit makes nitrogen *faster* than hydrogen. So the claim is
    defensible but not robust, and should not be quoted without the spread."""
    cmp = freeze.nitrogen_comparison(6000.0)
    assert cmp["ratio_byron"] < 8.0 < cmp["ratio_abinitio"] * 1.5


def test_carbon_has_no_three_body_answer_and_says_so() -> None:
    """`C + H + M` has no record in the evaluated literature at all, and `C + C + M` has one
    measurement. The module states that rather than inventing a coefficient."""
    assert "nucleation" in rates.CARBON_CHANNEL_NOTE
    assert rates.C_AR.category == "Experiment"
    assert rates.C_AR.t_max - rates.C_AR.t_min <= 1000.0


def test_total_and_effective_isp_are_not_the_same_number(
    flown: list[freeze.FreezeStation],
) -> None:
    """**Regression on a mistake this study actually made.** The exit station's flow speed over
    `g0` is a *total* specific impulse: it counts every kilogram expelled and carries no
    divergence, drift or vessel-mass normalisation. The companion's effective-Isp column is per
    kilogram of *launched slug*, a factor `(1+k)/k` larger for the same exhaust speed, and then
    additionally normalised by conventions this repository does not own.

    An earlier draft compared the two directly and reported a spurious 0.4% agreement with
    ADR-0016's 1080 s. The two differ by 5% from the mass convention alone, so any claim that
    they agree closely is two errors cancelling. This pins the gap so the confusion cannot
    silently return.
    """
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    total = flown[-1].speed / 9.80665
    effective = total * (1.0 + state.slug_ratio) / state.slug_ratio
    assert effective / total == pytest.approx(1.051, abs=0.002)
    assert effective - total > 50.0


def test_the_energy_conversion_fraction_is_convention_free(
    flown: list[freeze.FreezeStation],
) -> None:
    """What this study actually owns: the share of the wall-cap energy density that leaves as
    directed kinetic energy. It needs no mass convention and no launch ledger, which is why W5
    now leads with it rather than with an Isp."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    ideal = math.sqrt(2.0 * state.energy)
    converted = (flown[-1].speed / ideal) ** 2
    assert 0.35 < converted < 0.45
    # A narrower throat expands further and converts more -- the W6 lever, stated without Isp.
    narrow = freeze.history(state, throat_area=2.0, steps=64)
    assert (narrow[-1].speed / ideal) ** 2 > converted + 0.08
