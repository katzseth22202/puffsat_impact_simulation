"""Limits the chamber-wall fatigue and spike analyses must reproduce."""

from __future__ import annotations

import math

import pytest

from puffsat.walled_nozzle import chamber_fatigue as cf
from puffsat.walled_nozzle import wall_waves as ww


@pytest.mark.parametrize("ratio", [0.0, 0.4, 0.75, 1.0])
def test_breathing_overshoot_matches_the_parents_ramp_formula(ratio: float) -> None:
    """Undamped ramp-and-hold: `D = 1 + |sin(pi r)|/(pi r)`, and 2 for an instant step."""
    period = 1e-3
    _, x = cf.breathing_response(period, 0.0, ratio * period, 1e3, duration=20 * period)
    expected = 2.0 if ratio == 0 else 1 + abs(math.sin(math.pi * ratio)) / (math.pi * ratio)
    assert x.max() == pytest.approx(expected, abs=2e-3)


def test_hydrogen_grows_cracks_faster_than_air() -> None:
    for dk in (8.0, 15.0, 30.0):
        assert cf.rate_hydrogen(dk, 0.0, 818.0) > cf.rate_air(dk, 0.0, 0.0)


def test_rainflow_counts_every_ring_of_a_ringdown() -> None:
    """A decaying ring leaves half-cycles (ASTM residue); ten periods must count as ~ten cycles."""
    import numpy as np

    t = np.linspace(0, 10, 2001)
    series = np.exp(-0.1 * t) * np.sin(2 * math.pi * t)
    assert sum(c[2] for c in cf.rainflow(series)) == pytest.approx(10.0, abs=0.6)


def test_slow_rise_removes_the_fatigue_problem() -> None:
    """A rise of one breathing period gives `D = 1` and leaves the steel near 100 MPa."""
    chamber = cf.HYDROGEN
    fast = cf.evaluate(chamber, "hydrogen", 0.4e-3, chamber.efold_s[0], 0.02)
    slow = cf.evaluate(chamber, "hydrogen", chamber.breathing_period(), chamber.efold_s[0], 0.02)
    assert fast.pulses_a1 < cf.PULSES_PER_DEPARTURE
    assert slow.peak_mpa < 150.0
    assert slow.pulses_a1 > cf.PULSES_FIVE_FLIGHTS


def test_steel_to_carbon_step_reflects_81_percent_as_tension() -> None:
    """The parent's figure: (Z_s - Z_c)/(Z_s + Z_c) ~ 0.81 for a spike short against the steel."""
    stack = ww.Stack("t", (ww.Layer("A", 0.05, 5900.0, 7830.0, math.inf), ww.carbon(0.5)))
    (_, steel), _ = ww.response(stack, 1e9, 2e-6)
    assert steel.max() / 1e9 == pytest.approx(0.81, abs=0.01)


def test_thin_aluminium_is_invisible_to_a_long_spike() -> None:
    """A 5 mm layer crossed in 0.8 us cannot split a 30 us step: carbon tension barely moves."""
    with_al, without = (s for s in ww.stacks() if s.name.startswith("methane"))
    carbon_with = next(r for r in ww.evaluate(with_al, 30e-6) if r.layer == "IM7/8552")
    carbon_without = next(r for r in ww.evaluate(without, 30e-6) if r.layer == "IM7/8552")
    assert carbon_with.tension_per_gpa_mpa == pytest.approx(
        carbon_without.tension_per_gpa_mpa, rel=0.1
    )


def test_autofrettage_radial_bias() -> None:
    assert ww.interface_bias() == pytest.approx(-16.7e6, rel=0.01)
