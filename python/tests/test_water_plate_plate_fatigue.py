"""Limits the floor's stress-wave, spall and fatigue analysis must reproduce."""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.water_plate import plate_fatigue as pf

L = 0.030
P = 1.0e9


def _load(time: list[float], pressure: list[float]) -> pf.Load:
    return pf.Load("test", np.array(time), np.array(pressure))


def test_rainflow_reproduces_astm_e1049_example() -> None:
    """ASTM E1049 Fig. 6: ranges 3 (x0.5), 4 (x1.5), 6 (x0.5), 8 (x1.0), 9 (x0.5)."""
    series = np.array([-2, 1, -3, 5, -1, 3, -4, 4, -2], dtype=float)
    counts: dict[float, float] = {}
    for rng, _, count in pf.rainflow(series, 0.0):
        counts[rng] = counts.get(rng, 0.0) + count
    assert counts == {3.0: 0.5, 4.0: 1.5, 6.0: 0.5, 8.0: 1.0, 9.0: 0.5}


def test_held_load_relaxes_to_rigid_body_acceleration() -> None:
    """A held face pressure on a free plate settles to `sigma = -P (1 - x/L)`."""
    time, sigma = pf.stress_history(_load([0.0, 3e-3], [P, P]), L, 3.0)
    late = sigma[int(2.5e-3 / time[1])]
    expected = -P * (1.0 - np.linspace(0.0, 1.0, late.size))
    assert np.max(np.abs(late - expected)) < 0.03 * P


def test_short_pulse_reflects_as_full_tension() -> None:
    """A square pulse shorter than the round trip reflects off the free back as tension P."""
    tau = 0.2 * 2.0 * L / pf.longitudinal_speed()
    load = _load([0.0, tau, tau * 1.0001, 1e-4], [P, P, 0.0, 0.0])
    _, sigma = pf.stress_history(load, L, 1000.0)
    assert sigma.max() == pytest.approx(P, rel=0.02)


def test_rise_of_one_round_trip_suppresses_ringing() -> None:
    """A rise lasting one thickness-mode period leaves a slowly released load without tension."""
    period = 2.0 * L / pf.longitudinal_speed()
    times = [0.0, 1e-9, 0.5e-3, 1.5e-3]
    sharp = pf.stress_history(_load(times, [0.0, P, P, 0.0]), L, 1000.0)[1]
    times[1] = period
    ramped = pf.stress_history(_load(times, [0.0, P, P, 0.0]), L, 1000.0)[1]
    assert sharp.max() > 0.4 * P
    assert ramped.max() < 0.05 * P


def test_life_curve_hits_its_anchor_points() -> None:
    assert pf.cycles_to_failure(550.0, 0.0, pf.SMOOTH) == pytest.approx(1e7, rel=1e-9)
    assert pf.cycles_to_failure(810.0, 0.0, pf.SMOOTH) == pytest.approx(1e5, rel=1e-9)
    # A compressive mean is read as zero mean.
    assert pf.cycles_to_failure(245.0, -300.0, pf.NOTCHED) == pytest.approx(1e7, rel=1e-9)


def test_critical_radius_inverts_penny_k() -> None:
    a = pf.critical_radius(pf.MARAGING_300, P)
    assert pf.penny_k(P, a) == pytest.approx(pf.MARAGING_300.k_ic, rel=1e-12)
    assert pf.critical_radius(pf.MARAGING_350, P) < a


def test_spall_strength_from_pullback() -> None:
    """`0.5 rho c_b W` with c_b ~ 4.5 km/s: 227.5 m/s is ~4.1 GPa."""
    assert pf.spall_strength(227.5) == pytest.approx(4.12e9, rel=0.01)


def test_compressive_only_cycling_grows_no_crack() -> None:
    cycles = [(1e9, -1e9, 1.0)]
    assert pf.pulses_to_critical(cycles, 0.5e-3, 5e-3) == 100_000
    assert not math.isinf(pf.critical_radius(pf.MARAGING_300, 1e9))
