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


# --- N18 item 1: chemical ablation and blowing under hydrogen ------------------------------------


def _blown_gas(g0: float) -> wl.GasHistory:
    gas = _constant_gas(5500.0, 3.0e4)
    return wl.GasHistory(gas.time, gas.temp, gas.pressure, gas.emissivity, gas.h0, g0)


def test_chemistry_off_is_the_old_solver() -> None:
    """`g0 = 0` must leave every number the methane study printed untouched."""
    gas = _constant_gas(5500.0, 3.0e4, 0.5)
    r = wl.run(wl.pitch(), gas, pulses=1)
    assert r.chem_removed_per_pulse == 0.0
    assert r.min_blowing_correction == 1.0


def test_bare_steel_does_not_react() -> None:
    r = wl.run(wl.BARE, _blown_gas(0.4), pulses=1)
    assert r.chem_removed_per_pulse == 0.0


def test_pinned_surface_loses_carbon_at_the_closed_form_rate() -> None:
    """At the 3900 K cap and constant p, mdot = g ln(1 + B'(3900, p)) exactly, so the carbon
    removed is that times the time spent at the cap, over the density -- to the step at which
    the surface first reaches the cap."""
    from puffsat.walled_nozzle import blowing

    g0 = 0.4
    r = wl.run(wl.pitch(), _blown_gas(g0), pulses=1)
    p = 1.0e7
    rate = g0 * blowing.mass_flux_factor(blowing.blowing_parameter(3900.0, p))
    expected = rate * r.time_at_cap / wl.pitch().rho
    assert r.time_at_cap > 0.05
    # The surface is below the cap (and B' lower) during the climb, so expect a bit more than
    # the pinned-only integral, but within the few ms of that climb.
    assert r.chem_removed_per_pulse >= expected * 0.999
    assert r.chem_removed_per_pulse < expected * 1.15


def test_blowing_shields_the_surface_and_chemistry_is_stable_under_step_refinement() -> None:
    gas = _blown_gas(0.4)
    coarse = wl.run(wl.pitch(), gas, pulses=1)
    fine = wl.run(wl.pitch(), gas, pulses=1, dt_min=2.5e-7, dt_max=5e-5)
    assert coarse.min_blowing_correction < 0.5
    assert coarse.chem_removed_per_pulse == pytest.approx(fine.chem_removed_per_pulse, rel=0.02)


def test_chemical_heat_sink_cuts_the_thermal_term() -> None:
    gas = _blown_gas(0.4)
    base = wl.run(wl.pitch(), gas, pulses=1)
    sink = wl.run(wl.pitch(), gas, pulses=1, chem_heat=wl.CHEM_HEAT_C2H2)
    assert sink.removed_per_pulse < base.removed_per_pulse
