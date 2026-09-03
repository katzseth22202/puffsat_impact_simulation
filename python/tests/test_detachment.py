"""Acceptance tests for detachment along the column (Rung 2, N6 and Rung 5's `beta`).

The rung's exit criterion is not an analytic solution -- it is a *reconciliation*. The paper reports
`M_A` = 1.63/2.06 at one station and the solved history gives 0.35/0.58, so the tests have to show
that the gap is understood rather than merely observed: the paper's own route must reproduce its
own number, and the correction factor must close the difference exactly.

These read `data/results/cooling_history.csv`, so they are skipped rather than failed when it has
not been generated (`make analysis-expansion`).
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import detachment

pytestmark = pytest.mark.skipif(
    not detachment.DEFAULT_HISTORY.exists(),
    reason="needs data/results/cooling_history.csv (make analysis-expansion)",
)


@pytest.fixture(scope="module")
def stations() -> list[detachment.Station]:
    return detachment.load_stations()


@pytest.fixture(scope="module")
def cases(
    stations: list[detachment.Station],
) -> dict[tuple[float, str], list[detachment.Station]]:
    return detachment.by_case(stations)


# ---- The reconciliation ---------------------------------------------------------------------


def test_the_papers_own_route_reproduces_the_papers_own_number() -> None:
    """`M_A` = 1.63 must fall out of the bag field and bag density it is quoted from.

    Without this the disagreement below would be an assertion. 1.67 against a printed 1.63 is the
    rounding of the paper's own two-significant-figure inputs.
    """
    assert detachment.paper_alfven_mach() == pytest.approx(1.63, rel=0.04)


def test_the_state_mixing_factor_closes_the_gap_exactly(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """**The test that carries the finding.** The whole discrepancy is one density ratio.

    `v_A ~ rho^-1/2`, so dividing the post-expansion speed by a pre-expansion Alfven speed inflates
    `M_A` by `sqrt(rho_bag/rho_exit)`. If that factor did not reproduce the paper's 2.06 from the
    solved 0.58, some *other* difference would be hiding in the gap and the correction would not
    be safe to report.
    """
    hot_exit = cases[(75.0, "equilibrium")][-1]
    corrected = hot_exit.alfven_mach_design * detachment.state_mixing_factor(hot_exit.rho)
    assert corrected == pytest.approx(2.06, rel=0.02)


def test_the_flow_is_sub_alfvenic_at_every_station_of_every_case(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """N6's answer: the crossing the paper places at the exit does not happen inside the column."""
    for curve in cases.values():
        assert max(s.alfven_mach_design for s in curve) < 1.0


def test_the_alfven_mach_rises_monotonically_down_the_column(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """It is heading the right way, which is why the downstream crossing argument survives.

    A non-monotonic `M_A` would mean the flow re-magnetises somewhere, which is the failure mode
    N6 asked to be flagged.
    """
    for curve in cases.values():
        for a, b in pairwise(curve):
            assert b.alfven_mach_design >= a.alfven_mach_design - 1e-9


# ---- Why: the field is over-strength for the flow it steers ------------------------------------


def test_beta_is_far_below_one_everywhere_so_standoff_does_not_describe_the_expansion(
    stations: list[detachment.Station],
) -> None:
    """`eq:alfven_from_standoff` assumes `beta = 1`; the solved expansion runs at 0.01-0.08.

    This is the mechanism behind the `M_A` gap, not a separate finding: an over-strength field
    raises `v_A` and holds the flow sub-Alfvenic. The field is graded against the collision's
    snowplow pressure and then asked to steer the much thinner expansion that follows.
    """
    assert max(s.beta for s in stations) < 0.1
    assert min(s.beta for s in stations) > 0.0


def test_a_weaker_realizable_field_raises_the_alfven_mach(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """Rung 1's shortfall helps here, and the sign matters.

    Less field means a lower `v_A` and so a higher `M_A`. It is the one place the winding's
    inability to make the chamber gradient works in the design's favour -- and it is nowhere near
    enough to reach 1.
    """
    curve = cases[(75.0, "equilibrium")]
    throat = curve[0]
    assert throat.b_built_t < throat.b_design_t
    assert throat.alfven_mach_built > throat.alfven_mach_design
    assert throat.alfven_mach_built < 1.0


# ---- The downstream argument, quantified ---------------------------------------------------------


def test_downstream_crossing_is_the_identity_when_already_super_alfvenic() -> None:
    """A flow that has already crossed does not have to travel to cross."""
    assert detachment.downstream_crossing(1.0) == 1.0
    assert detachment.downstream_crossing(2.5) == 1.0


def test_downstream_crossing_follows_the_three_halves_power() -> None:
    """`M_A ~ R^{3/2}` on the `R^-5` against `R^-6` scalings, so crossing is at `M_A^{-2/3}`."""
    assert detachment.downstream_crossing(0.58) == pytest.approx((1.0 / 0.58) ** (2 / 3))
    assert detachment.downstream_crossing(0.125) == pytest.approx(4.0)


# ---- The physical-wall question ------------------------------------------------------------


def test_the_frozen_branch_ablates_sub_micron_per_pulse(
    stations: list[detachment.Station],
) -> None:
    """**Seth's correction, pinned.** A steady-state flux ratio is the wrong test for a 2 ms pulse.

    The first screen compared incident flux against graphite's re-radiation ceiling and called the
    hot frozen leg an ablator at 3.4x. Asked as depth per pulse instead, 3.4x over the ceiling for
    1.7 ms is 0.42 microns -- a centimetre of liner lasts ~24 000 pulses. The verdict flips.
    """
    frozen = [c for c in detachment.ablation_cases(stations) if c.branch == "frozen"]
    assert len(frozen) == 4
    assert max(c.depth_per_pulse_m for c in frozen) < 1e-6
    assert min(c.pulses_per_liner for c in frozen) > 10_000


def test_the_hot_equilibrium_legs_do_not_survive_even_as_a_transient(
    stations: list[detachment.Station],
) -> None:
    """The transient framing rescues the frozen branch and does not rescue this one.

    50 microns a pulse eats a centimetre of liner in about 200 pulses, against a vehicle life the
    paper puts in the thousands.
    """
    hot = next(
        c
        for c in detachment.ablation_cases(stations)
        if c.branch == "equilibrium" and c.closing_speed_km_s == pytest.approx(75.0)
    )
    assert hot.depth_per_pulse_m > 1e-5
    assert hot.pulses_per_liner < 1000.0


def test_every_exit_is_supersonic_so_a_diverging_wall_would_work_on_it(
    stations: list[detachment.Station],
) -> None:
    """The gas dynamics were never the blocker: a de Laval's diverging section is supersonic."""
    assert all(c.sonic_mach > 1.0 for c in detachment.ablation_cases(stations))


def test_ablation_is_zero_below_the_re_radiation_ceiling() -> None:
    """A wall that can shed what arrives does not ablate at all, however long the pulse."""
    ceiling = detachment.SIGMA_SB * detachment.GRAPHITE_SUBLIMATION_K**4
    assert detachment.ablation_depth(ceiling * 0.9, 1.0) == 0.0
    assert detachment.ablation_depth(ceiling * 2.0, 1e-3) > 0.0


def test_ablation_depth_is_linear_in_both_flux_and_time() -> None:
    """The transient scaling is the whole point: halving the pulse halves the damage."""
    ceiling = detachment.SIGMA_SB * detachment.GRAPHITE_SUBLIMATION_K**4
    base = detachment.ablation_depth(ceiling + 100e6, 2e-3)
    assert detachment.ablation_depth(ceiling + 100e6, 1e-3) == pytest.approx(base / 2)
    assert detachment.ablation_depth(ceiling + 200e6, 2e-3) == pytest.approx(base * 2)


def test_the_standoff_flux_is_far_below_what_would_ablate_anything() -> None:
    """On the paper's own booking a liner behind a field never reaches sublimation from radiation.

    Which makes the 4.9 kg/pulse the paper books for liner loss a charge against something other
    than the radiative flash -- worth asking the paper which.
    """
    ceiling = detachment.SIGMA_SB * detachment.GRAPHITE_SUBLIMATION_K**4
    assert detachment.standoff_flux_w_m2() < ceiling / 100.0
    assert detachment.ablation_depth(detachment.standoff_flux_w_m2(), 2e-3) == 0.0


# ---- The field window ------------------------------------------------------------------------


def test_the_closed_form_reproduces_the_solved_alfven_mach_on_the_frozen_branch(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """**The test the window rests on.** `M_A = M_sonic sqrt(gamma/2) sqrt(beta)`.

    Field and density cancel out of the ratio, leaving only the sonic Mach number and `beta`. If
    this did not hold, "lower the field to detach" would not follow, because `M_A` would not be a
    function of `beta` alone at fixed `M_sonic`.

    Frozen only: recombination moves the equilibrium branch's effective exponent off 5/3, which
    the next test bounds rather than ignores.
    """
    for (_, branch), curve in cases.items():
        if branch != "frozen":
            continue
        s = curve[-1]
        predicted = detachment.alfven_mach_closed_form(s.mach, s.beta)
        assert predicted == pytest.approx(s.alfven_mach_design, rel=1e-3)


def test_the_closed_form_over_predicts_equilibrium_by_a_bounded_amount(
    cases: dict[tuple[float, str], list[detachment.Station]],
) -> None:
    """Recombination moves the effective gamma, so `5/3` is the wrong exponent on that branch.

    Bounded rather than hidden: if the deviation grew past ~20% the closed form would stop being
    a safe way to reason about the equilibrium branch at all.
    """
    for (_, branch), curve in cases.items():
        if branch != "equilibrium":
            continue
        s = curve[-1]
        ratio = detachment.alfven_mach_closed_form(s.mach, s.beta) / s.alfven_mach_design
        assert 1.0 <= ratio < 1.2


def test_at_standoff_the_closed_form_is_the_papers_own_relation() -> None:
    """`beta = 1` gives `M_A = M_sonic/1.095` -- "the Alfven surface a tenth past the sonic throat".

    So the paper's physics is right and only its `beta` is wrong. That distinction is the whole
    reason this is a fixable problem rather than a broken argument.
    """
    for mach in (1.5, 2.7, 3.4):
        assert detachment.alfven_mach_closed_form(mach, 1.0) == pytest.approx(
            mach / 1.0954, rel=1e-3
        )


def test_containment_and_release_fields_invert_their_own_definitions() -> None:
    """`B_contain` must give `beta = 1`, and `B_release` must give `M_A = 1`. Otherwise nothing."""
    p, rho, v = 7.25e5, 2.512e-2, 16.19e3
    b_c = detachment.containment_field(p)
    assert p / (b_c**2 / (2 * detachment.MU0)) == pytest.approx(1.0)
    b_r = detachment.release_field(v, rho)
    assert v / (b_r / math.sqrt(detachment.MU0 * rho)) == pytest.approx(1.0)


def test_a_window_exists_on_every_leg_and_the_design_sits_above_it(
    stations: list[detachment.Station],
) -> None:
    """**The finding.** There is room between containing the plume and releasing it -- and the
    flown design is not in it, but above it.

    A window exists whenever the flow is supersonic past `M_sonic` = 1.095, which every exit is by
    a wide margin. So the design's failure to detach is a positioning error, not a conflict.
    """
    windows = detachment.field_windows(stations)
    assert len(windows) == 8
    assert all(w.window_exists for w in windows)
    assert all(w.design_above_window for w in windows)
    assert all(w.alfven_mach_at_standoff > 2.0 for w in windows)


def test_reaching_the_window_needs_a_harder_flare_than_the_design_has(
    stations: list[detachment.Station],
) -> None:
    """Flux conservation ties the exit field to the flare, so the fix is geometric."""
    assert all(w.required_area_ratio > 4.0 for w in detachment.field_windows(stations))


def test_one_flare_serves_the_whole_fleet_and_the_flown_value_is_below_it(
    stations: list[detachment.Station],
) -> None:
    """The actionable form: a single `A/A*` interval that every leg both contains and releases in.

    Narrower than any single leg's window, because it is their intersection -- and the flown 4
    is well below it, which is why no leg detaches today.
    """
    cw = detachment.common_window(stations)
    assert cw.exists
    assert cw.area_ratio_min > cw.flown_area_ratio
    assert 10.0 < cw.area_ratio_min < cw.area_ratio_max < 20.0


def test_buying_detachment_lets_the_snowplow_past(
    stations: list[detachment.Station],
) -> None:
    """The cost, and it is not optional: one static field cannot do both jobs at one station.

    A field weakened to the release ceiling is 3-8x under the collision's snowplow pressure there.
    That is what a physical wall in the diverging section would exist to take.
    """
    windows = detachment.field_windows(stations)
    assert all(w.snowplow_overshoot > 1.0 for w in windows)
    assert max(w.snowplow_overshoot for w in windows) > 3.0
