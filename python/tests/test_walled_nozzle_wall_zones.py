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


def test_doubling_the_pulse_scales_lengths_and_times_and_keeps_pressure() -> None:
    """Cube-root scaling: lengths and times x2^(1/3), volume x2, the same quasi-static pressure."""
    import numpy as np

    h = wz.History(
        np.array([0.5, 1.0]), np.array([1.4, 1.4]), np.array([0.0, 1e-6]), np.ones((2, 2)), 40.0
    )
    big = h.scaled(2.0)
    k = 2.0 ** (1.0 / 3.0)
    assert big.qsp == pytest.approx(h.qsp)
    assert big.volume == pytest.approx(80.0)
    assert big.r[0] == pytest.approx(1.4 * k)
    assert big.time[1] == pytest.approx(1e-6 * k)
    assert big.length_scale == pytest.approx(k)


def test_hoop_crack_life_of_one_cycle_matches_the_paris_closed_form() -> None:
    """One rest -> peak -> rest cycle per pulse on the air curve `C dK^m`, `dK = 0.713 sigma
    sqrt(pi a)`: `N = (a0^(1-m/2) - ac^(1-m/2)) / ((m/2 - 1) C k^m)`, failing at 80% of the
    layer when that comes before `K_fail`."""
    import numpy as np

    sigma, layer = 300.0, 0.010  # MPa, m
    c, m = 1.36e-10, 2.25
    k = 0.713 * sigma * math.sqrt(math.pi)
    a_c = min((wz.K_FAIL / (0.713 * sigma)) ** 2 / math.pi, 0.8 * layer)
    expected = (wz.HOOP_CRACK_START ** (1 - m / 2) - a_c ** (1 - m / 2)) / ((m / 2 - 1) * c * k**m)
    stress = np.array([-562.5e6, sigma * 1e6, -562.5e6])
    got = wz.hoop_crack_pulses(stress, layer, cf.rate_air, 142.0)
    assert got == pytest.approx(expected, rel=0.01)


def test_a_wall_that_never_goes_into_tension_never_cracks() -> None:
    import numpy as np

    stress = np.array([-562.5e6, -100e6, -300e6, -50e6, -562.5e6])
    assert wz.hoop_crack_pulses(stress, 0.01, cf.rate_hydrogen, 142.0) == math.inf


def test_ringing_shortens_the_hoop_crack_life() -> None:
    """Ten extra tensile cycles about the peak can only grow the crack faster."""
    import numpy as np

    one = np.array([-562.5e6, 400e6, -562.5e6])
    ring = np.r_[-562.5e6, np.tile([400e6, 100e6], 10), -562.5e6]
    assert wz.hoop_crack_pulses(ring, 0.01, cf.rate_hydrogen, 142.0) < wz.hoop_crack_pulses(
        one, 0.01, cf.rate_hydrogen, 142.0
    )


def test_autofrettage_rest_is_minus_three_quarters_yield() -> None:
    assert wz.rest_stress("wrap") == pytest.approx(-562.5e6)
    assert wz.rest_stress("steel") == pytest.approx(-562.5e6)
    assert wz.rest_stress("maraging") == pytest.approx(-0.6e9)
