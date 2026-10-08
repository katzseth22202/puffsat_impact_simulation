"""X-ray collision study: closed-form limits of the kinematics, EOS and one-zone fireball."""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.xray_collision import eos, multigroup
from puffsat.xray_collision.fireball import (
    Kinematics,
    collide,
    piston,
    planck_above,
    shocked_state,
)
from puffsat.xray_collision.materials import ARGON, CARBON, IRON


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


def test_radiation_adds_a_t4() -> None:
    """With trapped radiation the table carries e += a T^4/rho and p += a T^4/3."""
    plain = eos.build_table(CARBON)
    rad = eos.build_table(CARBON, radiation=True)
    rho, temp = 40.0, 470.0 * eos.EV / eos.K_B
    u = eos.A_RAD * temp**4
    assert rad.state(rho, temp)[0] - plain.state(rho, temp)[0] == pytest.approx(u / 3.0, rel=0.02)


def test_near_sun_reproduces_the_parents_fireball() -> None:
    """1 kg + 1 kg carbon foam at 10 kg/m^3 and 618 km/s: kT ~ 470 eV, mostly radiated."""
    r = collide(CARBON, (3.0 / (4.0 * math.pi * 10.0)) ** (1 / 3), 1.0, 618.0e3, rho0=10.0)
    assert r.kt_start_ev == pytest.approx(470.0, rel=0.05)
    assert r.zbar_shock == pytest.approx(6.0, abs=0.05)
    assert r.f_rad > 0.5
    assert r.closure == pytest.approx(1.0, abs=1e-4)


def test_piston_puts_all_heat_in_the_target() -> None:
    """A snowplow dissipates u^2/2 per swept kg: the target alone carries Q, so it runs hotter."""
    foam = collide(CARBON, (3.0 / (4.0 * math.pi * 10.0)) ** (1 / 3), 1.0, 618.0e3, rho0=10.0)
    plow = piston(CARBON, 1.0, 10.0, 1.0)
    assert plow.heat == pytest.approx(foam.heat)
    assert plow.kt_start_ev > 1.8 * foam.kt_start_ev
    assert plow.closure == pytest.approx(1.0, abs=1e-4)


def test_multigroup_pull_sees_argons_windows() -> None:
    """Cold argon is clear below its 11.6 eV resonance line and opaque above its 15.76 eV edge."""
    grp = multigroup.load_tops_groups()
    assert len(grp.edges_ev) == 121
    assert grp.edges_ev[0] == pytest.approx(1.0) and grp.edges_ev[-1] == pytest.approx(300.0)
    _, kp = grp.kappa(multigroup.COLD_KT_EV, 1.0e-3)
    lower = grp.edges_ev[:-1]
    assert np.max(kp[(lower >= 2.0) & (lower < 11.0)]) < 0.1  # [m^2/kg]
    assert np.min(kp[(lower >= 16.0) & (lower < 30.0)]) > 1e4  # ~40 Mb per atom


@pytest.mark.parametrize("kt_ev", [0.3, 4.0, 80.0])
def test_group_fractions_hold_all_of_sigma_t4(kt_ev: float) -> None:
    edges = multigroup.load_tops_groups().edges_ev
    assert multigroup.group_fractions(edges, kt_ev).sum() == pytest.approx(1.0, abs=1e-9)


def test_multigroup_closes_and_tracks_gray() -> None:
    """Per-group escape conserves energy, brackets itself, and stays near the gray answer."""
    rho0, m = 1.0e-3, 1.0
    planck = multigroup.collide_groups(rho0, m, emission="planck")
    ross = multigroup.collide_groups(rho0, m, emission="rosseland")
    gray = collide(ARGON, (3.0 * m / (4.0 * math.pi * rho0)) ** (1 / 3), 2.0, 66.0e3, rho0=rho0)
    assert planck.closure == pytest.approx(1.0, abs=1e-4)
    assert ross.f_rad <= planck.f_rad
    assert planck.f_rad == pytest.approx(gray.f_rad, abs=0.1)
    bands = (
        planck.f_below_resonance
        + planck.f_resonance_to_edge
        + planck.f_edge_to_euv
        + planck.f_above_euv
    )
    assert bands == pytest.approx(planck.f_rad, rel=1e-9)
    assert planck.f_below_resonance < planck.f_through_cold_target < planck.f_rad
