"""Closed forms the zone-by-zone wall study must reproduce."""

from __future__ import annotations

import math

import pytest

from puffsat.walled_nozzle import chamber_fatigue as cf
from puffsat.walled_nozzle import wall_zones as wz


def test_penny_crack_life_matches_the_paris_closed_form() -> None:
    """On the air curve `C dK^m` with `dK = k sqrt(a)`, `N = (a0^(1-m/2) - ac^(1-m/2)) /
    ((m/2 - 1) C k^m)`."""
    sigma = 400.0  # MPa
    c, m = 1.36e-10, 2.25
    k = 2.0 / math.pi * sigma * math.sqrt(math.pi)
    a_c = min((wz.K_FAIL * math.pi / (2.0 * sigma)) ** 2 / math.pi, 0.5 * 0.05)
    expected = (wz.CRACK_START ** (1 - m / 2) - a_c ** (1 - m / 2)) / ((m / 2 - 1) * c * k**m)
    assert wz.crack_pulses(sigma * 1e6, 0.05, cf.rate_air, 150.0) == pytest.approx(
        expected, rel=0.01
    )


def test_a_crack_already_past_critical_fails_at_once() -> None:
    assert wz.crack_pulses(3e9, 0.1, cf.rate_air, 150.0) == 0.0


def test_multilayer_keeps_hoop_and_mass_and_loses_bending() -> None:
    mono = wz._section("steel", 0.08)
    multi = wz._section("multilayer", 0.08)
    assert multi[0] == mono[0]
    assert multi[2] == mono[2]
    assert multi[1] == pytest.approx(mono[1] / wz.SHELLS**2)


def test_maraging_is_scored_only_behind_a_barrier() -> None:
    """Maraging is 'extreme' in hydrogen: no hydrogen life is credited, and it needs the barrier."""
    assert not wz.steel_of("maraging").hydrogen_qualified
    assert wz.steel_of("steel").hydrogen_qualified
    assert wz.steel_of("maraging").spall > wz.steel_of("steel").spall
