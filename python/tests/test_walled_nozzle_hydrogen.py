"""Acceptance tests for pure-hydrogen equilibrium and the Project 242 verdict (ask N10).

The exit criterion is the one the ask states: reproduce Project 242's thin-end case and say
whether its 2700 s baseline is reachable without recombination. The answer is known to be
*decidable* before the solver exists -- a frozen expansion has a hard ceiling -- which is what
makes this a test rather than a plot.
"""

from __future__ import annotations

import math

import pytest

from puffsat import eos_water
from puffsat.walled_nozzle import hydrogen


def test_full_dissociation_energy_is_the_textbook_value() -> None:
    """`D0(H2)/M(H2)` = 214.3 MJ/kg, from `eos_water`'s own JANAF-derived constant."""
    assert pytest.approx(214.3e6, rel=1e-3) == hydrogen.FULL_DISSOCIATION_ENERGY


@pytest.mark.parametrize("temp", [2000.0, 5000.0, 10000.0, 20000.0])
@pytest.mark.parametrize("rho", [1e-3, 1e-1, 1.0])
def test_nuclei_and_charge_are_conserved(rho: float, temp: float) -> None:
    """Checked on the reconstructed densities, so a bug in the state assembly cannot hide behind
    a converged residual."""
    s = hydrogen.state(rho, temp)
    nuclei = s.n_h + s.n_hp + 2.0 * s.n_h2
    assert nuclei == pytest.approx(2.0 * rho / hydrogen.M_H2, rel=1e-8)
    assert s.n_e == pytest.approx(s.n_hp, rel=1e-8)


def test_cold_endpoint_is_molecular_and_hot_endpoint_is_atomic() -> None:
    """Both ends are known without a solver."""
    assert hydrogen.state(1e-2, 1000.0).dissociated < 1e-6
    assert hydrogen.state(1e-2, 20000.0).dissociated > 0.999


def test_dissociation_matches_a_hand_solved_binary_equilibrium() -> None:
    """At 10 kK ionisation is a part per thousand, so the full solve must reduce to
    `H2 <=> 2H` solved by hand off `eos_water`'s own constant."""
    rho, temp = 1e-2, 10000.0
    k = math.exp(eos_water.ln_k_diatomic(eos_water.H2, temp))
    n = 2.0 * rho / hydrogen.M_H2
    hand = (-k + math.sqrt(k**2 + 8.0 * k * n)) / 4.0 / n
    assert hydrogen.state(rho, temp).dissociated == pytest.approx(hand, rel=1e-3)


def test_at_pressure_inverts_the_state() -> None:
    """The Rubbia case is quoted at `(p, T)`, so the inversion is load-bearing."""
    s = hydrogen.at_pressure(3.0e5, 10000.0)
    assert s.pressure == pytest.approx(3.0e5, rel=1e-6)


# --- The verdict the ask asks for ----------------------------------------------------------------


@pytest.mark.parametrize("pressure", [1.0e5, 3.0e5, 1.0e6])
def test_project_242s_baseline_is_unreachable_frozen(pressure: float) -> None:
    """**This is the answer to the Project 242 tension.**

    Their prose says molecular recombination "is not" fast; their carried baseline is 2700 s.
    A perfectly frozen expansion from equilibrium hydrogen at 10 kK and a few bar cannot exceed
    ~2100-2200 s however good the nozzle, because the dissociation store is simply not available
    to it. 2700 s is above that ceiling, so the arithmetic is right and the prose is wrong.
    """
    b = hydrogen.isp_bracket(pressure, hydrogen.RUBBIA_TEMPERATURE)
    assert b.isp_frozen < hydrogen.RUBBIA_ISP
    assert b.isp_equilibrium > hydrogen.RUBBIA_ISP  # and it is reachable, so it is not absurd
    assert 0.4 < b.implied_return(hydrogen.RUBBIA_ISP) < 0.8


def test_the_asks_own_enthalpy_figures_are_close_but_low() -> None:
    """The ask carries 351 MJ/kg as what 2700 s needs, 206 as sensible-only hydrogen and 421 as
    full recovery. The equilibrium solve gives 217-239 frozen and 430-453 total, so the ask's
    bracket is right in shape and about 5-10% low at both ends."""
    b = hydrogen.isp_bracket(3.0e5, hydrogen.RUBBIA_TEMPERATURE)
    assert b.frozen_enthalpy / 1e6 == pytest.approx(206.0, rel=0.15)
    assert b.enthalpy / 1e6 == pytest.approx(421.0, rel=0.10)
    assert pytest.approx(351.0, rel=0.02) == 0.5 * (hydrogen.RUBBIA_ISP * hydrogen.G0) ** 2 / 1e6


def test_hydrogen_is_essentially_fully_dissociated_at_a_few_bar() -> None:
    """The premise underneath everything above: at 10 kK and a few bar there is nothing left to
    dissociate, so the store is fully charged and the only question is whether it returns."""
    assert hydrogen.isp_bracket(3.0e5, 10000.0).dissociated > 0.999


def test_the_walled_chambers_density_suppresses_dissociation_but_not_much() -> None:
    """Ask N10 exercises the mechanism across four orders of magnitude in three-body rate by
    running Rubbia's few bar against our few hundred. The charge barely moves over that span --
    86-100% -- which is the same finding the methane chamber gives."""
    thin = hydrogen.isp_bracket(3.0e5, 10000.0).dissociated
    dense = hydrogen.state(2.0, 10000.0).dissociated
    assert dense < thin
    assert dense > 0.85
