"""Acceptance tests for the snowplow containment budget -- and for withdrawing P9's fix.

The question: P9 proposed flaring to `A/A*` = 11.3-14.8 to drop the exit field into a window that
both contains the expansion and lets it detach. The cost was left uncomputed. It is computed here,
and **it sinks the proposal**: weakening the field for the expansion un-contains the *collision*,
and contact with a 94 600 K front is five orders of magnitude past what graphite can shed.

The tests are structured around that reversal. The budget is the paper's own booking turned into a
share of pulse energy; the scaling is checked in closed form so the conclusion cannot rest on one
arithmetic slip; and the last test records the reason the proposal is not merely unaffordable but
unnecessary.
"""

from __future__ import annotations

import pytest

from puffsat import detachment, snowplow

# ---- The budget, from the paper's own booking ---------------------------------------------------


def test_the_liner_budget_is_the_papers_own_number_in_energy_terms() -> None:
    """4.9 kg/pulse x 59.7 MJ/kg = 292 MJ, which is 0.46% of a 62.9 GJ pulse.

    Everything downstream is measured against this, so it has to be the paper's booking and not an
    assumption introduced here.
    """
    assert snowplow.liner_energy_budget() == pytest.approx(2.925e8, rel=0.01)
    assert snowplow.liner_budget_share() == pytest.approx(0.00465, rel=0.02)


def test_the_budget_inverts_back_to_the_booked_ablation() -> None:
    """Spending the allowed share must reproduce 4.9 kg, or the conversion is wrong."""
    kg = snowplow.ablation_if_contacted(snowplow.liner_budget_share())
    assert kg == pytest.approx(snowplow.BOOKED_ABLATION_KG, rel=1e-9)


# ---- What weakening the field does to the collision ---------------------------------------------


def test_the_flown_flare_is_standoff_by_construction() -> None:
    """`A/A*` = 4 is where the graded profile was derived, so `beta` there is 1 by definition."""
    assert snowplow.snowplow_beta(snowplow.FLOWN_AREA_RATIO) == pytest.approx(1.0)
    assert snowplow.containment_radius_ratio(1.0) == 1.0


def test_the_front_grows_as_the_closed_form_says() -> None:
    """**The scaling the conclusion rests on**, checked against its closed form rather than trusted.

    `beta = (A/A* / 4)^2` and the front expands as `beta^{1/(2 gamma)}`, so the growth is
    `(A/A* / 4)^{0.6}` at `gamma = 5/3`. If this did not hold in closed form the whole argument
    would be one arithmetic slip away from nothing.
    """
    for ar in (4.0, 6.0, 8.0, 11.3, 14.8, 20.0):
        growth = snowplow.containment_radius_ratio(snowplow.snowplow_beta(ar))
        assert growth == pytest.approx((ar / 4.0) ** 0.6, rel=1e-9)


def test_flaring_past_the_flown_value_puts_the_front_on_the_wall_in_both_bore_readings() -> None:
    """**The finding, and it is not what P9 assumed.**

    The front already fills the bore for the last three quarters of the crossing, so once the field
    is weakened it must exceed whatever bore contains it. A flared bore does not rescue it: the
    wall grows as `sqrt(A/A*)` while the front grows as `sqrt(A/A*) x (A/A*/4)^{0.6}`.

    So P7's bore ambiguity, which decides several other things, does *not* decide this one.
    """
    cases = {c.area_ratio: c for c in snowplow.containment_cases()}
    assert not cases[4.0].touches_cylinder and not cases[4.0].touches_flared
    for ar in (6.0, 8.0, 11.3, 14.8, 20.0):
        assert cases[ar].touches_cylinder, f"cylinder reading at {ar}"
        assert cases[ar].touches_flared, f"flared reading at {ar}"


def test_the_window_flare_puts_the_front_well_past_the_wall() -> None:
    """At P9's own recommended 11.3 the front is 86% past the bore, not marginally over."""
    case = next(c for c in snowplow.containment_cases() if c.area_ratio == 11.3)
    assert case.snowplow_beta == pytest.approx(7.98, abs=0.02)
    assert case.radius_growth > 1.8


# ---- Why contact is disqualifying rather than merely costly ------------------------------------


def test_the_front_radiates_five_orders_past_what_graphite_can_shed() -> None:
    """94 600 K goes as `T^4`: 4.5 TW/m^2 against graphite's 13.1 MW/m^2 at sublimation.

    This is why the question is only *whether* contact happens. There is no regime in which a
    liner takes this flux and survives, so no dwell-time refinement changes the verdict.
    """
    ratio = snowplow.front_radiative_flux() / (snowplow.SIGMA_SB * 3900.0**4)
    assert ratio > 1e5


def test_even_one_percent_of_the_pulse_blows_the_ablation_budget() -> None:
    """The budget is 0.46%, so a contacting front has almost no room before it is over."""
    assert snowplow.ablation_if_contacted(0.01) > 2.0 * snowplow.BOOKED_ABLATION_KG
    assert snowplow.ablation_if_contacted(0.10) > 20.0 * snowplow.BOOKED_ABLATION_KG


# ---- The reason the proposal is unnecessary, not merely unaffordable ----------------------------


@pytest.mark.skipif(
    not detachment.DEFAULT_HISTORY.exists(),
    reason="needs data/results/cooling_history.csv (make analysis-expansion)",
)
def test_detachment_happens_downstream_without_any_flare_change() -> None:
    """**Why P9 should be withdrawn rather than merely costed.**

    P3 established that the paper's downstream argument survives: past the last coil the plume's
    pressure falls as `R^-5` against a vacuum field's `R^-6`, so `M_A` crosses 1 a short way out.
    At the flown flare that lands 1.4-2.0 exit radii downstream. The nozzle does not have to
    achieve detachment *inside* itself, so it does not have to buy the window that costs the liner.
    """
    cases = detachment.by_case(detachment.load_stations())
    crossings = [
        detachment.downstream_crossing(curve[-1].alfven_mach_design) for curve in cases.values()
    ]
    assert all(1.0 < c < 2.5 for c in crossings), crossings
