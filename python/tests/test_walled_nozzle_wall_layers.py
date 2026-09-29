"""Tests for the coated-wall conduction solve (`wall_layers`).

The known answer is the semi-infinite solid under a convective boundary: with gas at `T_g`
and coefficient `h`, the surface follows

    T_s = T_i + (T_g - T_i) [1 - exp(beta^2) erfc(beta)],   beta = h sqrt(alpha t) / k

which the solver must reproduce on bare steel, where the 30 mm slab is effectively semi-infinite
over a 250 ms pulse.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.walled_nozzle import wall_layers as wl


def _constant_gas(temp: float, h: float, emissivity: float = 0.0) -> wl.GasHistory:
    return wl.GasHistory(
        time=np.array([0.0, 1.0]),
        temp=np.array([temp, temp]),
        pressure=np.array([1.0e7, 1.0e7]),
        emissivity=np.array([emissivity, emissivity]),
        h0=h,
    )


def test_bare_steel_matches_the_convective_semi_infinite_solution() -> None:
    """Surface temperature after one 250 ms period against the erfc solution, to 1%."""
    t_gas, h = 1500.0, 5.0e4  # stays below the 1750 K bare-steel cap
    run = wl.run(wl.BARE, _constant_gas(t_gas, h), pulses=1)
    alpha = wl.STEEL_K / (wl.STEEL_RHO * wl.STEEL_C)
    beta = h * math.sqrt(alpha * wl.PULSE_PERIOD) / wl.STEEL_K
    exact = wl.INITIAL_TEMP + (t_gas - wl.INITIAL_TEMP) * (
        1.0 - math.exp(beta**2) * math.erfc(beta)
    )
    assert run.peak_surface == pytest.approx(exact, rel=0.01)


def test_melt_cap_holds_and_books_removal() -> None:
    """A load far beyond what 1 mm of graphite can conduct pins the surface at the melt point
    and removes material, rather than letting the surface run away."""
    run = wl.run(wl.graphite(1e-3, 3.0), _constant_gas(7000.0, 2.0e5, 1.0), pulses=1)
    assert run.peak_surface == pytest.approx(wl.CARBON_MELT, abs=1.0)
    assert run.removed_per_pulse > 0.0


def test_a_low_conductivity_coat_protects_the_steel_better() -> None:
    """Same load, same thickness: the 3 W/m/K coat must leave the steel cooler than 40 W/m/K."""
    gas = _constant_gas(5000.0, 2.0e4, 0.5)
    low = wl.run(wl.graphite(1e-3, 3.0), gas, pulses=1)
    high = wl.run(wl.graphite(1e-3, 40.0), gas, pulses=1)
    assert low.peak_steel < high.peak_steel
