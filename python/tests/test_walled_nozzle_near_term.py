"""Tests for the near-term chamber (asks A1-A4 of the methane-7000 K list).

The known answers: the fixed point closes on the rod's energy; the isentrope conserves energy
exactly (`integral h drho = rho0 e0` for an emptying vessel); each constrained branch conserves
its nuclei; the graphite saturation reproduces JANAF; and the Isp ledger with a plug reduces to
the paper's eq:eta_isp.
"""

from __future__ import annotations

import math
from functools import cache

import numpy as np
import pytest

from puffsat import eos_methane
from puffsat.walled_nozzle import chamber, near_term

METHANE_7000 = near_term.PAIRINGS[1]
METHANE_10000 = near_term.PAIRINGS[2]


@cache
def _solution(thermo: str = "janaf") -> near_term.ChamberSolution:
    return near_term.solve_charge(METHANE_7000, thermo)


@cache
def _isentrope(branch_name: str) -> tuple[near_term.Branch, near_term.Isentrope]:
    sol = _solution()
    br = near_term.branch(branch_name, sol)
    return br, near_term.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=1e-4, steps=160)


def test_fixed_point_shares_the_rod_energy_over_everything_expelled() -> None:
    """`(charge + rod + plug) u = m_rod w^2 / 2`, at the chamber temperature asked for."""
    sol = _solution()
    assert sol.total_mass * sol.energy == pytest.approx(near_term.PULSE_ENERGY, rel=1e-9)
    assert sol.rho == pytest.approx(sol.total_mass / near_term.VOLUME, rel=1e-12)


def test_chamber_is_a_sphere_with_the_papers_wall() -> None:
    """20 m^3 and 35.6 m^2 are the same sphere; the radius is what the Bartz station uses."""
    sol = _solution()
    assert sol.wall_area == pytest.approx(near_term.WALL_AREA, rel=2e-3)
    assert sol.radius == pytest.approx(1.684, abs=1e-3)


def test_isentrope_hands_out_exactly_the_stored_energy() -> None:
    """Along an isentrope `d(rho e)/drho = h`, so an emptying vessel delivers `rho0 e0` in total.
    This is the identity `eta_blowdown` rests on; it checks the integrator and the EOS together.
    """
    _, tab = _isentrope("equilibrium")
    rho, h = tab.rho[::-1], tab.enthalpy[::-1]
    integral = float(np.sum(0.5 * (h[1:] + h[:-1]) * np.diff(rho)))
    expected = tab.rho[0] * tab.energy[0] - tab.rho[-1] * tab.energy[-1]
    assert integral == pytest.approx(expected, rel=2e-4)


def test_blowdown_efficiency_sits_below_the_peak_convention() -> None:
    """Late parcels leave an emptier, cooler chamber: the blowdown can never beat the peak state,
    and the impulse-averaged speed can never beat the peak speed."""
    br, tab = _isentrope("equilibrium")
    for eff in near_term.efficiencies(_solution(), br, tab, (4.0, 14.0)):
        assert 0.0 < eff.eta_blowdown < eff.eta_peak < 1.0
        assert eff.speed_impulse < eff.speed_peak


def test_frozen_carbon_is_the_floor_and_keeps_its_carbon() -> None:
    """Freezing the carbon strands its bonds, so it must not beat equilibrium -- and the frozen
    species must scale with density alone."""
    br_eq, tab_eq = _isentrope("equilibrium")
    br_fr, tab_fr = _isentrope("frozen_carbon")
    sol = _solution()
    eq = near_term.efficiencies(sol, br_eq, tab_eq, (14.0,))[0]
    fr = near_term.efficiencies(sol, br_fr, tab_fr, (14.0,))[0]
    assert fr.eta_peak < eq.eta_peak
    rho = sol.rho / 10.0
    comp, _ = br_fr.state(rho, 3000.0)
    comp0 = eos_methane.composition(sol.rho, sol.pairing.temp, sol.feed)
    i_c3 = eos_methane.MOLECULES.index(eos_methane.C3)
    assert comp.n_mol[i_c3] == pytest.approx(comp0.n_mol[i_c3] / 10.0, rel=1e-12)
    n_f = rho / sol.feed.formula_mass
    assert comp.n_h_nuclei == pytest.approx(sol.feed.hc_ratio * n_f, rel=1e-9)
    assert comp.n_c_nuclei == pytest.approx(n_f, rel=1e-9)


def test_frozen_chamber_state_is_the_equilibrium_one() -> None:
    """At the state it was frozen at, the frozen branch *is* equilibrium."""
    sol = _solution()
    br = near_term.branch("frozen_carbon", sol)
    p, e = br(sol.rho, sol.pairing.temp)
    p_eq, e_eq = eos_methane.pressure_energy(sol.rho, sol.pairing.temp, sol.feed)
    assert p == pytest.approx(p_eq, rel=1e-8)
    assert e == pytest.approx(e_eq, rel=1e-8)


@pytest.mark.parametrize(
    ("temp", "janaf_log_kf"), [(3000.0, -4.296), (4000.0, -1.207), (5000.0, 0.632)]
)
def test_graphite_saturation_reproduces_janaf(temp: float, janaf_log_kf: float) -> None:
    """Monatomic C over graphite: JANAF C-003 `log Kf` is `log10(p_C / 1 bar)` directly."""
    p = near_term.graphite_saturation_density(temp) * eos_methane.K_B * temp
    assert math.log10(p / 1e5) == pytest.approx(janaf_log_kf, abs=0.01)


def test_graphite_ceiling_condenses_to_saturation_and_conserves_nuclei() -> None:
    """Below ~4000 K the chamber's carbon is far supersaturated: gas-phase C sits at the
    saturation density and the rest is condensate, with every nucleus accounted for."""
    sol = _solution()
    br = near_term.branch("graphite_ceiling", sol)
    rho, temp = sol.rho / 50.0, 3000.0
    gas, n_s = br.state(rho, temp)
    n_f = rho / sol.feed.formula_mass
    assert n_s > 0.0
    assert gas.n_c == pytest.approx(near_term.graphite_saturation_density(temp), rel=1e-8)
    assert gas.n_c_nuclei + n_s == pytest.approx(n_f, rel=1e-9)
    assert gas.n_h_nuclei == pytest.approx(sol.feed.hc_ratio * n_f, rel=1e-8)


def test_condensation_releases_energy_at_fixed_state() -> None:
    """At the same `(rho, T)`, condensing carbon can only lower the stored energy."""
    sol = _solution()
    rho, temp = sol.rho / 50.0, 3000.0
    e_ceiling = near_term.branch("graphite_ceiling", sol)(rho, temp)[1]
    e_gas = eos_methane.pressure_energy(rho, temp, sol.feed)[1]
    assert e_ceiling < e_gas


def test_plug_ledger_is_the_papers_eta_isp() -> None:
    """`I = [(k+1+P) v - w] / ((k+P) g0)` (paper eq:eta_isp), and `P = 0` is the old ledger."""
    k, p, v = 27.5, near_term.PLUG_RATIO, 9000.0
    ledger = chamber.isp_ledger(v, k, p)
    w = chamber.CLOSING_SPEED
    assert ledger.isp_effective == pytest.approx(((k + 1 + p) * v - w) / ((k + p) * chamber.G0))
    old = chamber.isp_ledger(v, k)
    assert old.isp_effective == pytest.approx(((1 + k) * v - w) / (k * chamber.G0), rel=1e-15)


def test_h_minus_saha_matches_a_hand_value() -> None:
    """At 7000 K the H- Saha factor is lambda_e^3 / 4 * exp(0.754 eV / kT) = 9.6e-27 m^3."""
    comp = eos_methane.Composition(
        n_c=0.0, n_h=1.0, n_hp=0.0, n_c_ions=(0.0,) * 6, n_e=1.0, n_mol=(0.0,) * 9
    )
    assert near_term.h_minus_density(comp, 7000.0) == pytest.approx(9.6e-27, rel=0.02)


def test_band_fraction_is_the_planck_integral() -> None:
    """0.25-1.55 um holds 90% of a 7000 K blackbody -- the share the H- floor can speak for."""
    assert near_term.band_fraction(7000.0) == pytest.approx(0.896, abs=0.005)
    assert near_term._planck_fraction(1e9) == pytest.approx(1.0, abs=1e-5)  # series tail
