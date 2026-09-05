"""Acceptance tests for the snowplow front's spreading (R9).

The exit criterion is that the module **reproduces the two numbers the paper already carries and
cites to this repository** -- a 94 600 K shocked layer and a 21.1 km/s spreading speed -- from the
strong-shock piston relation and the water EOS, without either being fitted. Those numbers were
quoted to us; until this module they had never been derived here.

`test_the_piston_relation_reproduces_the_papers_front_temperature` is the one that carries the
argument. If `e = v^2/2` on `eos_water` did not land on 94 600 K, then either the paper's figure
came from somewhere else or this closure is the wrong one, and every contact station below would
be built on it.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import coupling, front

# ---- The closure ---------------------------------------------------------------------------------


def test_the_piston_relation_reproduces_the_papers_front_temperature() -> None:
    """`e = v^2/2` at the cold leg's 45.58 km/s must give the paper's 94 600 K."""
    temp, _c_s = front.shock_state(45.58e3)
    assert temp == pytest.approx(94_600.0, rel=0.005)


def test_the_same_state_reproduces_the_papers_spreading_speed() -> None:
    """And its sound speed must be the 21.1 km/s `sec:needle_through_fog` spreads at."""
    _temp, c_s = front.shock_state(45.58e3)
    assert c_s == pytest.approx(21.1e3, rel=0.01)


def test_the_sound_speed_barely_depends_on_the_compression_ratio() -> None:
    """The compression is bracketed rather than known, so the answer must not turn on it."""
    speeds = [p.sound_speed_m_s for p in front.shock_sound_speeds()]
    assert (max(speeds) - min(speeds)) / min(speeds) < 0.03


def test_the_hot_leg_shocks_harder() -> None:
    """Monotonicity, and the reason the hot leg's front opens at a wider angle."""
    cold, _ = front.shock_state(45.58e3)
    hot, _ = front.shock_state(75.0e3)
    assert hot > cold


def test_the_spread_to_speed_ratio_falls_as_the_projectile_slows() -> None:
    """**The whole correction.** An ideal gas gives a constant ratio and an exact cone.

    Water's dissociation and ionisation thresholds make the post-shock temperature rise more
    slowly than `v^2`, so `c_exp/v` falls -- the front opens more slowly than a fixed cone and
    touches the wall *later*. If this ratio were constant the correction would not exist.
    """
    ratios = [front.shock_state(v)[1] / v for v in (10e3, 20e3, 30e3, 45.58e3)]
    assert all(b > a for a, b in pairwise(ratios))
    assert ratios[0] < 0.35
    assert ratios[-1] > 0.44


def test_the_tabulated_closure_matches_the_direct_solve() -> None:
    """The integrator interpolates for speed; interpolation must not cost accuracy."""
    import numpy as np

    speeds = [45.58e3 * f for f in np.geomspace(0.02, 1.0, 40)]
    table = front.spread_speed_table(speeds)
    for v in (5e3, 12e3, 25e3, 40e3):
        assert table(v) == pytest.approx(front.shock_state(v)[1], rel=0.05)


# ---- The trajectory ------------------------------------------------------------------------------


def test_the_front_grows_monotonically_and_the_projectile_slows() -> None:
    """Both halves of the system, asserted separately so a sign error cannot hide."""
    run = front.integrate(45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    radii = [s.radius_m for s in run.stations]
    speeds = [s.speed_m_s for s in run.stations]
    assert all(b >= a for a, b in pairwise(radii))
    assert all(b <= a for a, b in pairwise(speeds))


def test_the_closed_front_touches_later_than_a_frozen_one() -> None:
    """The artifact being avoided: freezing `c_exp` while `v` falls opens the cone too fast."""
    closed = front.integrate(45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    frozen = front.integrate(
        45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100, fixed_cone=True
    )
    assert closed.contact_x_m is not None
    assert frozen.contact_x_m is not None
    assert closed.contact_x_m > frozen.contact_x_m


def test_the_straight_cone_is_nearly_right_after_all() -> None:
    """ADR-0012's construction survives: the closed integration moves contact by well under a metre.

    This is the result, not a formality. The two effects that could have moved it -- deceleration
    and the weakening shock -- very nearly cancel, which is why a fixed-angle cone was a good
    approximation for a reason ADR-0012 did not give.
    """
    run = front.integrate(45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=8000, sample_every=200)
    assert run.contact_x_m is not None
    assert abs(run.contact_x_m - run.cone_contact_x_m) < 0.5


def test_a_wider_wall_is_touched_later() -> None:
    """R15's whole lever: measuring to the liner rather than the bag bore buys column length."""
    bore = front.integrate(45.58e3, 1.9, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    liner = front.integrate(
        45.58e3, 1.9, wall_m=front.LINER_CHAMBER_M, steps=4000, sample_every=100
    )
    assert bore.contact_x_m is not None
    assert liner.contact_x_m is not None
    assert liner.contact_x_m > bore.contact_x_m


def test_a_faster_spread_touches_sooner_and_demands_more_field() -> None:
    """The asymmetry R9 flags: the conservative coupling assumption is the risky field one."""
    slow = front.integrate(45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    fast = front.integrate(45.58e3, 1.9, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    assert slow.contact_x_m is not None
    assert fast.contact_x_m is not None
    assert fast.contact_x_m < slow.contact_x_m
    assert front.field_demanded_at(fast.contact_x_m) > front.field_demanded_at(slow.contact_x_m)


def test_swept_mass_is_capped_at_the_bag_and_not_at_the_wall() -> None:
    """The clearance gap is vacuum, which is exactly why R15 answers 'the liner'.

    Integrating to the liner must not sweep more mass than integrating to the bag bore.
    """
    bore = front.integrate(45.58e3, 1.0, wall_m=front.BAG_BORE_M, steps=4000, sample_every=100)
    liner = front.integrate(
        45.58e3, 1.0, wall_m=front.LINER_CHAMBER_M, steps=4000, sample_every=100
    )
    assert liner.slug_ratio == pytest.approx(bore.slug_ratio, rel=1e-9)


def test_the_field_demand_follows_the_flown_profile() -> None:
    """A shelf starting further out sits lower, which is the entire saving R9 is chasing."""
    assert front.field_demanded_at(3.0) > front.field_demanded_at(6.0)
    assert front.field_demanded_at(6.0) == pytest.approx(9.0, rel=0.05)
    assert not math.isnan(front.field_demanded_at(0.1))


def test_the_projectile_radius_is_the_paper_s_compact_slug() -> None:
    """A 25 kg ice sphere, which is what makes the front have to grow at all."""
    assert pytest.approx(0.187, abs=0.005) == front.PROJECTILE_RADIUS_M
    assert coupling.BORE_RADIUS > front.PROJECTILE_RADIUS_M
