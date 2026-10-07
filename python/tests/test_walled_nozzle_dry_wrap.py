"""Limits the dry-wrap layer chain must reproduce."""

from __future__ import annotations

import numpy as np
import pytest

from puffsat.walled_nozzle import dry_wrap as dw


def _wall(bonded: bool) -> dw.Wall:
    return dw.Wall(radius=2.0, wrap=dw.VECTRAN, wrap_t=0.04, e_z=3e9, n_s=4, n_w=10, bonded=bonded)


def test_slow_pressure_gives_the_static_hoop_strain() -> None:
    """A ramp far slower than the ring's period: the added wrap strain is `p r / sum(E t) / ...`,
    the thin-ring static value `p r / sum(E_h t)` (no overshoot)."""
    wall = dw.Wall(radius=2.0, wrap=dw.VECTRAN, wrap_t=0.04, e_z=3e9, n_s=2, n_w=6)
    _, _, eh, _ = wall.layers()
    ts = wall.layers()[0]
    p = 10e6
    load = np.minimum(np.arange(30_000) / 25_000.0, 1.0) * p  # 25 ms ramp, period ~2.5 ms
    r = dw.simulate(wall, load, 1e-6, 0.03)
    static = p * wall.radius / float(np.sum(eh * ts))
    assert r.peak_wrap_strain - r.rest_wrap_strain == pytest.approx(static, rel=0.05)


def test_short_pulse_lifts_a_dry_wrap_but_not_a_bonded_one() -> None:
    spike = np.r_[np.full(20, 2e9), np.zeros(5000)]  # 2 GPa for 20 us
    dry = dw.simulate(_wall(bonded=False), spike, 1e-6, 2e-3)
    bonded = dw.simulate(_wall(bonded=True), spike, 1e-6, 2e-3)
    assert dry.max_gap_mm > 0.01
    assert bonded.max_gap_mm == 0.0
    assert dry.outer_speed > 0.0


def test_sizing_holds_the_design_peak_to_the_liner_swing() -> None:
    qsp, r = 12.4e6, 2.2
    t_w = dw.sized_wrap(qsp, r, dw.KEVLAR)
    stiffness = dw.STEEL_E * dw.STEEL_T + dw.KEVLAR.hoop_modulus * t_w
    assert dw.DESIGN_D * qsp * r / stiffness == pytest.approx(dw.SWING, rel=1e-12)


def _strip(n: int = 60) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = np.arange(n) * 0.02
    return s, np.full(n, 2.0), np.full(n, 0.04)


def test_coupled_chain_under_uniform_load_matches_the_single_ring() -> None:
    """A load uniform along the wall bends nothing: away from the clamps every station is the
    single ring."""
    s, radius, wrap_t = _strip()
    spike = np.r_[np.full(20, 2e9), np.zeros(3000)]
    field = np.repeat(spike[:, None], s.size, axis=1)
    wall = dw.Wall(radius=2.0, wrap=dw.VECTRAN, wrap_t=0.04, e_z=3e9, n_s=2, n_w=8)
    ring = dw.simulate(wall, spike, 1e-6, 1.5e-3)
    coupled = dw.simulate_coupled(dw.VECTRAN, 3e9, s, radius, wrap_t, field, 1e-6, 1.5e-3, 2, 8)
    mid = coupled[s.size // 2 - 5 : s.size // 2 + 5]
    assert mid.max() == pytest.approx(ring.peak_wrap_strain, rel=0.02)


def test_liner_bending_barely_changes_a_narrow_bands_fling() -> None:
    """The spike hands its momentum to the wrap in microseconds; the liner's bending acts over
    ~0.4 ms, and lifted dry layers are free rings. So coupling leaves a narrow band's fling
    within ~10% of its free ring (the liner's rebound adds a few percent)."""
    s, radius, wrap_t = _strip(100)
    spike = np.r_[np.full(20, 2e9), np.zeros(3000)]
    field = np.zeros((spike.size, s.size))
    band = (s > 0.95) & (s < 1.05)
    field[:, band] = spike[:, None]
    wall = dw.Wall(radius=2.0, wrap=dw.VECTRAN, wrap_t=0.04, e_z=3e9, n_s=2, n_w=8)
    ring = dw.simulate(wall, spike, 1e-6, 1.5e-3)
    coupled = dw.simulate_coupled(dw.VECTRAN, 3e9, s, radius, wrap_t, field, 1e-6, 1.5e-3, 2, 8)
    assert coupled.max() == pytest.approx(ring.peak_wrap_strain, rel=0.1)
