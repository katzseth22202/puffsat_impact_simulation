"""Limits the shell (beam-on-elastic-foundation) model must reproduce."""

from __future__ import annotations

import numpy as np
import pytest

from puffsat.walled_nozzle import shell_response as sr

STEEL = sr.WallSpec("steel 20 mm", 20.0, 200e9 * 0.02, sr.plate_bending(200e9, 0.02), 7830 * 0.02)
S = np.arange(0.0, 4.0, 0.01)
R = np.full(S.size, 1.5)


def _uniform(history: np.ndarray) -> np.ndarray:
    return np.repeat(history[:, None], S.size, axis=1)


def test_slow_uniform_pressure_gives_the_static_hoop_strain() -> None:
    """A free ring settles at the membrane value `p R / sum(E t)`. The clamped shell overshoots
    it just inside each clamp by the classic edge factor `1 + exp(-pi)` ~ 1.043 (Timoshenko)."""
    p = 5e6
    hist = np.minimum(np.arange(40) / 30.0, 1.0) * p  # 30 ms ramp, ring period ~1.9 ms
    membrane = p * 1.5 / STEEL.hoop_stiffness
    ring = sr.respond(STEEL, S, R, _uniform(hist), 1e-3, 0.04, bending=False)
    shell = sr.respond(STEEL, S, R, _uniform(hist), 1e-3, 0.04)
    assert ring == pytest.approx(membrane, rel=0.02)
    assert shell == pytest.approx((1.0 + np.exp(-np.pi)) * membrane, rel=0.02)


def test_sudden_uniform_step_doubles_the_strain_in_a_ring() -> None:
    p = 5e6
    hist = np.full(20, p)
    saved = sr.ZETA
    sr.ZETA = 0.0
    try:
        peak = sr.respond(STEEL, S, R, _uniform(hist), 1e-4, 2e-3, bending=False)
    finally:
        sr.ZETA = saved
    assert peak == pytest.approx(2.0 * p * 1.5 / STEEL.hoop_stiffness, rel=0.02)


def test_neighbours_hold_back_a_narrow_struck_band() -> None:
    field = np.zeros((40, S.size))
    band = (S > 1.95) & (S < 2.05)  # 10 cm band, against a ~0.1 m characteristic length
    field[:5, band] = 200e6  # a 50 us hammer tap
    ring = sr.respond(STEEL, S, R, field, 1e-5, 2e-3, bending=False)
    shell = sr.respond(STEEL, S, R, field, 1e-5, 2e-3, bending=True)
    assert shell < 0.7 * ring
