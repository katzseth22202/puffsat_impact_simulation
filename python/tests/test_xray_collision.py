"""X-ray collision study: closed-form limits of the kinematics, EOS and one-zone fireball."""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.xray_collision import eos
from puffsat.xray_collision.fireball import Kinematics, collide, planck_above, shocked_state
from puffsat.xray_collision.materials import ARGON, IRON


@pytest.mark.parametrize("k", [1.0, 3.0, 10.0])
def test_heat_is_reduced_mass_energy(k: float) -> None:
    """A resting target of k m heats by k/(1+k) of the impactor KE; the rest is CM motion."""
    kin = Kinematics(2.0, 2.0 * k, 67.0e3)
    v_cm = kin.v_rel / (1.0 + k)
    assert kin.heat_fraction == pytest.approx(k / (1.0 + k))
    assert kin.heat + 0.5 * kin.m_total * v_cm**2 == pytest.approx(kin.ke_ship_frame)


def test_saha_limits() -> None:
    """Neutral when cold and dense; fully stripped when very hot and dilute; p = (1+Z) n k T."""
    temps = np.array([3000.0, 1.0e9])
    p, _, z = eos.saha(ARGON, 1.0e3, temps[:1])
    assert z[0] < 1e-6
    p, _, z = eos.saha(ARGON, 1.0e-6, temps[1:])
    assert z[0] == pytest.approx(18.0, rel=1e-3)
    n = 1.0e-6 / (ARGON.mass_amu * eos.AMU)
    assert p[0] == pytest.approx(19.0 * n * eos.K_B * 1.0e9, rel=1e-3)


def test_lowering_only_raises_ionization() -> None:
    temps = np.array([10.0, 30.0]) * eos.EV / eos.K_B
    _, _, z_ideal = eos.saha(IRON, 3.0e4, temps)
    _, _, z_low = eos.saha(IRON, 3.0e4, temps, lowering=True)
    assert np.all(z_low > z_ideal)


def test_table_inverts_energy() -> None:
    table = eos.build_table(ARGON)
    for rho, t_ev in [(1.0e-3, 2.0), (7.0e3, 30.0), (1.0, 8.0)]:
        temp = t_ev * eos.EV / eos.K_B
        _, e, _ = eos.saha(ARGON, rho, np.array([temp]))
        assert table.temp(rho, float(e[0])) == pytest.approx(temp, rel=0.02)


def test_planck_tail() -> None:
    assert planck_above(0.0) == 1.0
    assert planck_above(4.965114) == pytest.approx(0.25, abs=1e-3)  # shortward of Wien lambda_max
    assert planck_above(10.0) == pytest.approx(15 / math.pi**4 * math.exp(-10) * 1366.0, rel=0.01)


def test_hugoniot_holds_the_jump_energy() -> None:
    """The shocked state carries e = u_p^2/2 and compresses ~4-6x (ionizing ideal plasma)."""
    table = eos.build_table(ARGON)
    s = shocked_state(table, 33.5e3)
    assert 4.0 < s.rho / ARGON.rho0 < 6.5
    p = table.state(s.rho, s.temp)[0]
    u_s = 33.5e3 * s.rho / (s.rho - ARGON.rho0)
    assert p == pytest.approx(ARGON.rho0 * u_s * 33.5e3, rel=1e-3)


def test_energy_closes() -> None:
    r = collide(ARGON, 1.0e-2, 3.0)
    assert r.closure == pytest.approx(1.0, abs=1e-4)
    assert r.heat_fraction == pytest.approx(0.75)


def test_opaque_limit_is_adiabatic() -> None:
    """With the opacity scaled to near zero emission vanishes and the heat ends as expansion."""
    r = collide(IRON, 1.0e-2, 3.0, kappa_scale=1e-12)
    assert r.f_rad < 1e-6
    assert 0.3 * (4.0 * r.m_impactor) * r.v_exp**2 / r.heat == pytest.approx(1.0, abs=0.02)
