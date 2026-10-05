"""Limits the argon arm's atomic spray EOS must reproduce (ADR-0055)."""

from __future__ import annotations

import math
from itertools import pairwise

import numpy as np
import pytest

from puffsat import eos_water as ew
from puffsat.water_plate import spray_eos as se

MIX = se.SprayMixture(k=10.0)


def test_cold_mixture_is_an_ideal_monatomic_gas() -> None:
    rho, temp = 1.0, 1000.0
    p, e = se.pressure_energy(MIX, rho, temp)
    n_atoms = sum(MIX.per_kg.values())
    assert p == pytest.approx(n_atoms * rho * ew.K_B * temp, rel=1e-9)
    assert e == pytest.approx(1.5 * n_atoms * ew.K_B * temp, rel=1e-9)


def test_weak_first_stage_matches_hand_saha_for_pure_argon() -> None:
    # Pure argon, first stage only in play: x^2/(1-x) = K1/n_Ar.
    mix = se.SprayMixture(k=1e12)
    rho, temp = 0.1, 9000.0
    n_ar = rho / se.M_AR
    k1 = math.exp(ew.ln_k_saha(se.IP_AR_EV[0] * ew.EV, 6.0, 1.0, temp))
    x = (-k1 + math.sqrt(k1 * k1 + 4.0 * k1 * n_ar)) / (2.0 * n_ar)
    s = se.state(mix, rho, temp)
    assert s.fractions["Ar"][1] == pytest.approx(x, rel=1e-6)


def test_charge_neutrality_holds() -> None:
    for temp in (5e3, 3e4, 2e5):
        s = se.state(MIX, 0.5, temp)
        supply = sum(s.n_element[n] * s.mean_charge(n) for n in s.n_element)
        assert supply == pytest.approx(s.n_e, rel=1e-9)


def test_very_hot_gas_is_fully_stripped() -> None:
    s = se.state(MIX, 1e-3, 1e8)
    assert s.mean_charge("Ar") == pytest.approx(18.0, abs=1e-3)
    assert s.mean_charge("O") == pytest.approx(8.0, abs=1e-3)
    assert s.mean_charge("H") == pytest.approx(1.0, abs=1e-6)


def test_energy_rises_with_temperature() -> None:
    temps = np.geomspace(300.0, 1e6, 80)
    e = [se.pressure_energy(MIX, 1.0, float(t))[1] for t in temps]
    assert all(b > a for a, b in pairwise(e))


def test_no_argon_reproduces_hot_water_less_its_bond_energy() -> None:
    # At 30 kK water holds no molecules, so the two EOSs differ only by the atomization energy.
    rho, temp = 1.0, 3.0e4
    p_w, e_w = ew.pressure_energy(rho, temp)
    p_s, e_s = se.pressure_energy(se.SprayMixture(k=0.0), rho, temp)
    assert p_s == pytest.approx(p_w, rel=1e-3)
    assert e_s + ew.D_AT / ew.M_H2O == pytest.approx(e_w, rel=1e-3)
