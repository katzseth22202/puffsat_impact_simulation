"""Fast checks of the hot-hydrogen asks' pieces (`walled_nozzle.hot_hydrogen`)."""

from __future__ import annotations

import math
from itertools import pairwise

import numpy as np
import pytest

from puffsat import eos_methane, tops
from puffsat.walled_nozzle import hot_hydrogen as hh
from puffsat.walled_nozzle import near_term as nt


def test_charge_holds_the_rods_energy_at_any_closing_speed() -> None:
    """`solve_charge(energy=...)`: the whole chamber shares the pulse at the held temperature."""
    for w in hh.SPEEDS:
        sol = hh.chamber(7000.0, w)
        assert sol.pairing.temp == 7000.0
        assert sol.total_mass * sol.energy == pytest.approx(0.5 * nt.ROD_MASS * w**2, rel=1e-9)
    # k follows w: a faster rod holds the temperature with more charge.
    assert hh.chamber(7000.0, 74e3).charge_mass < hh.chamber(7000.0, 81e3).charge_mass


def test_today_is_the_fixed_60_kg_charge_at_4884_k() -> None:
    sol = hh.chamber(None)
    assert sol.charge_mass == pytest.approx(hh.TODAY_CHARGE, rel=1e-9)
    assert sol.pairing.temp == pytest.approx(4884.0, abs=5.0)


def test_hotter_needs_less_charge_but_the_polyethylene_floors_it() -> None:
    """The rod and plug's 6.9 kg do not store hydrogen's energy, so k falls slower than u rises."""
    k = [hh.chamber(t).slug_ratio for t in hh.TEMPERATURES]
    assert all(a > b for a, b in pairwise(k))
    assert 7.0 < hh.chamber(10000.0).slug_ratio < 8.0  # not the parent's 5.5


def test_hc_ratio_counts_the_polyethylene() -> None:
    sol = hh.chamber(None)
    n_c = (nt.ROD_MASS + nt.PLUG_MASS) / (eos_methane.M_C + 2 * eos_methane.M_H)
    n_h = sol.charge_mass / eos_methane.M_H + 2 * n_c
    assert hh.hc_ratio(sol) == pytest.approx(n_h / n_c)


def test_tops_interpolation_reproduces_its_nodes() -> None:
    hc, path = hh.TOPS_PULLS[0]
    g = tops.parse_tops_gray(path.read_text())
    i, j = 20, 1
    got = hh.tops_rosseland(hc, float(g.rho_grid[i]), float(g.t_grid[j]))
    assert got == pytest.approx(float(g.kappa_rosseland[i, j]), rel=1e-9)
    assert math.isnan(hh.tops_rosseland(hc, 0.3, 5000.0))  # below the OPLIB floor


def test_copper_rise_matches_the_constant_flux_closed_form() -> None:
    """`dT = 2 q sqrt(t / pi) / e` for a step in flux on a semi-infinite face."""
    q, t_end = 1.0e8, 0.05
    times = np.linspace(0.0, t_end, 400)
    got = hh.copper_rise(list(times), [q] * len(times))
    eff = math.sqrt(hh.COPPER[0] * hh.COPPER[1] * hh.COPPER[2])
    assert got == pytest.approx(2.0 * q * math.sqrt(t_end / math.pi) / eff, rel=1e-9)


def test_coolant_capacity_is_hydrogen_cp_times_the_rise() -> None:
    """77 K -> 700 K of H2 gas: ~14.4 kJ/kg/K on average, so ~9 MJ/kg."""
    assert 8.0e6 < hh.coolant_capacity() < 10.0e6
