"""Unbonded (contact) interfaces in the wall-wave model: limits they must reproduce."""

from __future__ import annotations

import math

import pytest

from puffsat.walled_nozzle import wall_waves as ww

STEEL = (5900.0, 7830.0)
L1, L2 = 0.04, 0.02


def _two(bonded: bool) -> ww.Stack:
    inner = ww.Layer("inner", L1, *STEEL, math.inf, bonded=bonded)
    outer = ww.Layer("outer", L2, *STEEL, math.inf)
    return ww.Stack("two", (inner, outer))


def _square(amplitude: float, duration: float) -> ww.Load:
    return lambda t: amplitude if t < duration else 0.0


def test_compression_held_on_the_face_does_not_open_the_contact() -> None:
    """A slow ramp-and-hold is compressive throughout, so bonded and unbonded agree exactly."""
    ramp = 200e-6

    def load(t: float) -> float:
        return 1e8 * min(t / ramp, 1.0)

    bonded = ww.simulate(_two(True), load, 400e-6)
    contact = ww.simulate(_two(False), load, 400e-6)
    for (_, a), (_, b) in zip(bonded.peak, contact.peak, strict=True):
        assert b == pytest.approx(a, abs=1e3)
    assert contact.gaps[0] == 0.0


def test_a_short_pulse_is_trapped_in_the_outer_shell() -> None:
    """The momentum trap. A pulse shorter than the outer shell's round trip reflects off its free
    back as tension, which an unbonded contact cannot pass. The outer shell leaves with all the
    impulse, the inner shell is left at rest, and no tension returns into the inner shell."""
    p, tau = 1e9, 0.5 * 2 * L2 / STEEL[0]
    impulse = p * tau
    run = ww.simulate(_two(False), _square(p, tau), 30e-6, cells_thin=32)
    inner, outer = run.momentum
    assert outer == pytest.approx(impulse, rel=0.02)
    assert abs(inner) < 0.02 * impulse
    assert run.gaps[0] > 0.0
    assert run.peak[0][1].max() < 0.05 * p


def test_bonded_layers_ring_and_pass_tension_back() -> None:
    """The same pulse in bonded layers reflects tension through the whole stack."""
    p, tau = 1e9, 0.5 * 2 * L2 / STEEL[0]
    run = ww.simulate(_two(True), _square(p, tau), 30e-6, cells_thin=32)
    assert run.peak[0][1].max() > 0.9 * p
    assert sum(run.momentum) == pytest.approx(p * tau, rel=0.02)


def test_a_preload_holds_the_contact_shut_against_small_tension() -> None:
    """Shrink-fit contact pressure above the arriving tension keeps the shells together."""
    p, tau = 1e7, 0.5 * 2 * L2 / STEEL[0]
    held = ww.simulate(_two(False), _square(p, tau), 30e-6, cells_thin=32, preload=2e7)
    bonded = ww.simulate(_two(True), _square(p, tau), 30e-6, cells_thin=32)
    assert held.gaps[0] == 0.0
    assert held.momentum[0] == pytest.approx(bonded.momentum[0], abs=1e-6 * p * tau)
