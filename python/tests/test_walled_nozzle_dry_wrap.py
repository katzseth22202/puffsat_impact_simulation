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
