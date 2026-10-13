"""Tests for the finite-rate nozzle (S26 of the paper repo).

The known answers are the two limits of the march. A rate a million times the measured one must
reproduce the branch's own shifting-equilibrium isentrope; a zero rate must freeze the particle
count per kilogram at the throat. Between them the exit speed must rise monotonically with the
rate. The constrained equilibrium the exchange model rests on must be ordinary equilibrium at
`phi = 1`, conserve nuclei at any `phi`, and give `n_H^2 / n_H2 = K phi` exactly.
"""

from __future__ import annotations

import math
from functools import cache

import pytest

from puffsat import eos_methane
from puffsat.walled_nozzle import finite_rate as fr
from puffsat.walled_nozzle import near_term

METHANE_7000 = near_term.PAIRINGS[1]
THROAT = 0.198
AREA_RATIO = 30.0


@cache
def _solution() -> near_term.ChamberSolution:
    return near_term.solve_charge(METHANE_7000, "janaf")


@cache
def _setup(model: str) -> tuple[fr.Model, near_term.Isentrope]:
    sol = _solution()
    mod = fr.model(model, sol)
    tab = near_term.isentrope(mod.branch, sol.rho, sol.pairing.temp, rho_ratio_end=1e-5, steps=300)
    return mod, tab


def _march(model: str, scale: float, step: float = fr.DEFAULT_STEP) -> fr.Expansion:
    mod, tab = _setup(model)
    return fr.march(mod, tab, 0, THROAT, AREA_RATIO, rate_scale=scale, step=step)


def test_constrained_composition_at_phi_one_is_equilibrium() -> None:
    """`phi = 1` removes the constraint, so the composition is the ordinary equilibrium one."""
    feed = _solution().feed
    rho, temp = 0.05, 4500.0
    eq = eos_methane.composition(rho, temp, feed)
    con, _ = eos_methane.constrained_composition(rho, temp, feed, 0.0)
    assert con.n_h == pytest.approx(eq.n_h, rel=1e-8)
    for a, b in zip(con.n_mol, eq.n_mol, strict=True):
        assert a == pytest.approx(b, rel=1e-6, abs=1e-30 * eq.n_neutral_heavy)


def test_constrained_composition_conserves_nuclei_and_counts_particles_up_with_phi() -> None:
    """Any `phi` keeps the feed's C and H nuclei. `phi > 1` penalises every molecule by
    `phi^(atoms - 1)`, so the particle count per kilogram rises with it."""
    feed = _solution().feed
    rho, temp = 0.05, 4500.0
    n_f = rho / feed.formula_mass
    totals = []
    for ln_phi in (-0.5, 0.0, 0.5, 1.5):
        comp, _ = eos_methane.constrained_composition(rho, temp, feed, ln_phi)
        assert comp.n_c_nuclei == pytest.approx(n_f, rel=1e-8)
        assert comp.n_h_nuclei == pytest.approx(feed.hc_ratio * n_f, rel=1e-8)
        ln_k2 = eos_methane.ln_k_molecule_thermo(eos_methane.H2, temp, feed.thermo)
        ratio = comp.n_h**2 / comp.density_of(eos_methane.H2)
        assert math.log(ratio) == pytest.approx(ln_k2 + ln_phi, abs=1e-6)
        totals.append(comp.n_total)
    assert totals == sorted(totals)


def test_frozen_carbon_pool_scales_only_the_hydrogen_constant() -> None:
    """In the frozen-carbon pool `phi` moves H against H2 and leaves the carbon per kg alone."""
    sol = _solution()
    br = near_term.FrozenCarbon(sol.feed, sol.rho, sol.pairing.temp)
    rho, temp = 0.1 * sol.rho, 4000.0
    base, _ = br.state(rho, temp)
    shifted, _ = br.state(rho, temp, ln_phi=1.0)
    ln_k2 = eos_methane.ln_k_molecule_thermo(eos_methane.H2, temp, sol.feed.thermo)
    ratio = shifted.n_h**2 / shifted.density_of(eos_methane.H2)
    assert math.log(ratio) == pytest.approx(ln_k2 + 1.0, abs=1e-6)
    assert shifted.n_h > base.n_h
    assert shifted.n_c == base.n_c
    assert shifted.n_h_nuclei == pytest.approx(base.n_h_nuclei, rel=1e-10)


@pytest.mark.parametrize("model", fr.MODELS)
def test_fast_rate_reproduces_the_shifting_equilibrium_exit(model: str) -> None:
    """A rate 1e6 times the measured one keeps `phi = 1` all the way, which is the branch's own
    isentrope: the exit speed must match `near_term.nozzle_exit` from the same chamber state."""
    _, tab = _setup(model)
    ex = near_term.nozzle_exit(tab, 0, AREA_RATIO)
    assert ex is not None
    fast = _march(model, 1e6)
    assert fast.speed == pytest.approx(ex.speed, rel=3e-3)
    assert abs(fast.ln_phi_exit) < 1e-2


@pytest.mark.parametrize("model", fr.MODELS)
def test_zero_rate_freezes_the_particle_count_from_the_throat(model: str) -> None:
    """No sink, no change: particles per kilogram leave the exit as they left the throat."""
    frozen = _march(model, 0.0)
    assert frozen.particles_exit == pytest.approx(frozen.particles_throat, rel=1e-8)
    assert frozen.speed < _march(model, 1e6).speed


@pytest.mark.parametrize("model", fr.MODELS)
def test_recombination_and_vacuum_thrust_rise_with_the_rate(model: str) -> None:
    """More rate returns more of the bond store and more vacuum thrust per unit mass flow.

    The bare exit *speed* need not rise: at a fixed area ratio a gas that freezes early keeps a
    higher gamma and converts its thermal enthalpy faster, so it can leave slightly faster but
    colder and at lower pressure. `u + p / (rho u)`, which counts the exit pressure, is ordered.
    """
    runs = [_march(model, s) for s in (0.0, 0.01, 0.1, 1.0, 1e6)]
    bonds = [r.bond_returned for r in runs]
    assert bonds == sorted(bonds)
    thrusts = [r.effective_speed for r in runs]
    assert thrusts == sorted(thrusts)
    assert runs[0].speed < min(r.speed for r in runs[1:])


def test_march_converges_under_step_halving() -> None:
    """Backward Euler is first order; the default step must already sit inside 0.2%."""
    coarse = _march("frozen_carbon", 1.0, step=2.0 * fr.DEFAULT_STEP)
    fine = _march("frozen_carbon", 1.0, step=0.5 * fr.DEFAULT_STEP)
    mid = _march("frozen_carbon", 1.0)
    assert mid.speed == pytest.approx(fine.speed, rel=2e-3)
    assert abs(mid.speed - fine.speed) < abs(coarse.speed - fine.speed) + 1.0


def test_wall_angle_sets_only_the_clock() -> None:
    """A slower-opening cone is longer by `tan(15)/tan(10)` and gives the gas more time."""
    mod, tab = _setup("frozen_carbon")
    ten = fr.march(mod, tab, 0, THROAT, AREA_RATIO, half_angle=math.radians(10.0))
    fifteen = fr.march(mod, tab, 0, THROAT, AREA_RATIO, half_angle=math.radians(15.0))
    ratio = math.tan(math.radians(15.0)) / math.tan(math.radians(10.0))
    assert ten.length == pytest.approx(fifteen.length * ratio, rel=1e-9)
    assert ten.transit_time > fifteen.transit_time
    assert ten.bond_returned >= fifteen.bond_returned
    assert ten.effective_speed >= fifteen.effective_speed


def test_required_carbon_rate_is_the_bray_threshold() -> None:
    """`k_req tau_exp n_partner n_M = 1` at the half-release station, by definition."""
    _, tab = _setup("exchange_carbon")
    req = fr.required_carbon_rate(_solution(), tab, THROAT, AREA_RATIO)
    assert req.coefficient > 0.0
    assert req.coefficient * req.tau_expansion * req.n_partner * req.n_third == pytest.approx(1.0)
    assert req.carbon_free_throat > req.carbon_free_half > req.carbon_free_exit
