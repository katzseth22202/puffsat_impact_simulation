"""X-ray collision study: closed-form limits of the kinematics, EOS and one-zone fireball."""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.xray_collision import anvil, eos, multigroup
from puffsat.xray_collision.fireball import (
    Kinematics,
    collide,
    piston,
    planck_above,
    shocked_state,
)
from puffsat.xray_collision.materials import ARGON, CARBON, IRON, IRON_HOT


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


def test_steel_target_runs_cooler_than_carbon_at_400_kms() -> None:
    """More energy per iron atom, but split over ~23 free electrons and paid into ionization:
    at k = 1 the iron fireball starts cooler than carbon's, and radiates about as much."""
    fe = piston(IRON_HOT, 1.0, 1.0, 1.0, 400.0e3)
    c = piston(CARBON, 1.0, 1.0, 1.0, 400.0e3)
    assert fe.heat == pytest.approx(c.heat)
    assert fe.kt_start_ev < 0.7 * c.kt_start_ev
    assert fe.zbar_shock > 20.0
    assert fe.light_of_ship_ke == pytest.approx(c.light_of_ship_ke, abs=0.05)
    assert fe.closure == pytest.approx(1.0, abs=1e-4)


def test_hot_iron_pull_reaches_kev() -> None:
    """The 3 keV pull is not held at the 60 eV edge: hot dilute iron is far more transparent."""
    hot = eos.build_table(IRON_HOT)
    kr_hot = hot.state(1.0, 300.0 * eos.EV / eos.K_B)[2]
    kr_held = eos.build_table(IRON).state(1.0, 300.0 * eos.EV / eos.K_B)[2]
    assert kr_hot < 0.5 * kr_held


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


def test_anvil_split_reproduces_the_parents_pass_through_table() -> None:
    """templateArxiv.tex tab:pass_through_isp: 71% of 127.3 GJ radiated, 80% delivered, the plate
    catching 80% of the forward momentum, a 2 kg target plus the case as fixed spent mass."""
    energy = 0.8 * 0.71 * 127.3e9
    for case, cap, isp, spent in (
        (0.0, 150e3, 5040, 5.6),
        (5.0, 150e3, 2930, 14.0),
        (10.0, 150e3, 2240, 24.0),
        (0.0, 100e3, 4330, 7.4),
    ):
        s = anvil.best_split(energy, 2.0 + case, cap, 0.8)
        assert s.isp == pytest.approx(isp, rel=0.003)
        assert s.spent == pytest.approx(spent, abs=0.06)
        assert s.u_forward <= cap * (1 + 1e-9)


def test_pulse_split_is_symmetric_and_best_split_equally() -> None:
    """eq. pulse_split: swapping the masses changes nothing; at a fixed total, equal is best."""
    assert anvil.pulse_split(1e9, 1.0, 3.0) == pytest.approx(anvil.pulse_split(1e9, 3.0, 1.0))
    assert anvil.pulse_split(1e9, 2.0, 2.0) > anvil.pulse_split(1e9, 1.0, 3.0)


def test_foam_on_anvil_is_the_flipped_snowplow() -> None:
    """Heat is the reduced-mass share k/(1+k) of the projectile KE with k = anvil/foam, and the
    light share counts dead sheet mass against the projectile."""
    bare = anvil.foam_on_anvil(50.0)
    assert bare.heat_share == pytest.approx(50.0 / 60.0, rel=1e-9)
    sheeted = anvil.foam_on_anvil(50.0, sheet_mass=1.8)
    assert sheeted.light_share < bare.light_share
    assert bare.rt_graphite > bare.pulse_time > bare.rt_iron


def test_case_delivery_limits() -> None:
    assert anvil.case_delivery(1.0, 1.0, 5.0, 0.0) == pytest.approx(1.0)
    assert anvil.case_delivery(0.0, 1.0, 3.0, 0.0) == pytest.approx(0.25)


def test_entry_hole_plug_needs_a_foam_longer_than_the_case() -> None:
    """A foam shorter than the case has passed the hole before any parcel stagnates, so all of
    the light meets an open hole. Narrowing the foam lengthens it, plugs the hole, and helps."""
    kernel = anvil.emission_kernel()
    wide = anvil.hole_delivery(1.0, 0.95, kernel, case_length=1.5)
    assert wide.foam_length < 1.5
    assert wide.open_share == pytest.approx(1.0, abs=1e-3)
    assert wide.delivery_plugged == pytest.approx(wide.delivery_static, abs=1e-3)
    narrow = anvil.hole_delivery(0.5, 0.95, kernel, case_length=0.5)
    assert narrow.open_share < 0.25
    assert narrow.delivery_plugged > wide.delivery_plugged + 0.2


def test_emission_kernel_is_fast_against_the_pulse() -> None:
    """A stagnated parcel sends half its light out within ns, so emission tracks arrival."""
    t, cum = anvil.emission_kernel()
    assert cum[-1] == pytest.approx(1.0)
    assert float(np.interp(0.5, cum, t)) < 0.1e-6


def test_narrower_foam_carries_heavier_sheets() -> None:
    """Same volume, side area 4 V / d: halving the diameter roughly doubles the sunward C-C."""
    wide, narrow = anvil.sheet_mass(1.0), anvil.sheet_mass(0.5)
    assert 1.4 < wide < 1.7
    assert narrow > 1.6 * wide
    assert anvil.sheet_mass(0.5) < anvil.PROJECTILE_MASS


def test_jupiter_fall_boost_and_opaque_charge() -> None:
    """A fall from Jupiter to 9.86 R_sun needs ~20.8 km/s to cross Earth at 100 km/s, and a
    near-solid polyethylene charge fireball is so opaque it radiates essentially nothing."""
    from puffsat.xray_collision import charge

    assert charge.boost_dv() == pytest.approx(20.76e3, rel=0.01)
    assert charge.boost_dv(aphelion=charge.AU) > charge.boost_dv()
    fire = charge.charge_fireball(30.0)
    assert fire.f_rad < 1e-3
    assert fire.tau0 > 1e5
    assert fire.drift == pytest.approx(400e3 / 31.0)
    assert charge.plate_area(0.0, 0.0) == pytest.approx(math.pi * 100.0)
    assert charge.plate_area(0.30, 4.0) > charge.plate_area(0.30, 0.0)


def test_charge_fireball_recombines_but_keeps_its_bonds() -> None:
    """The dense charge fireball keeps its ionization in equilibrium all the way to the plate
    (Damkohler numbers far above 1), but arrives too hot for H2: the bond toll is paid in full."""
    from puffsat.xray_collision import charge, freeze

    tab = charge.ch2_table()
    f30 = freeze.freeze(30.0, tab)
    assert 0.15 < f30.ion_share_start < 0.25
    assert f30.ion_share_frozen < 0.01
    assert f30.da_ion_min > 1e3
    assert f30.kt_plate_ev > 0.5  # hotter than any H2 window
    f100 = freeze.freeze(100.0, tab)
    assert f100.velocity_loss_bonds_all_frozen > 3.0 * f30.velocity_loss_bonds_all_frozen
    assert f100.velocity_loss_bonds_all_frozen < 0.1


def test_bore_collar_ledger() -> None:
    """The committed collar ledger: the jet's own structural increment is a fraction of a percent
    to ~1% of a 100 t ship, the eroded collar leaves faster than the charge's exhaust when only
    the arriving column heats it, and even the hot-jet erosion ceiling keeps the mass ratio
    under 3."""
    import csv

    from puffsat.xray_collision import collar

    with collar.OUT.open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["grid"] == "x1"]
    assert len(rows) == 4
    for r in rows:
        assert 0.0 < float(r["jet_extra"]) < float(r["collar_extra"])
        assert float(r["collar_extra"]) < 0.05 * collar.FINAL_MASS
        assert float(r["mass_ratio_with_collar_best"]) <= float(r["mass_ratio"])
        assert float(r["mass_ratio_with_collar_worst"]) < 3.0
        assert float(r["eroded_per_final_best"]) < float(r["eroded_per_final_worst"])
        # The rings tile the plate, so the innermost one carries several times the mean.
        assert float(r["ring0_impulse"]) > 5.0 * float(r["plate_mean_impulse"])


def test_overtake_pusher_ledger() -> None:
    """Result 11: with the merged mass drifting onto the plate, k ~ 10 meets a 1/3 propellant
    fraction for the 11 km/s departure on every plate, k = 40 does not, and the finer grid only
    raises the impulse."""
    import csv

    from puffsat.xray_collision import overtake

    with overtake.OUT.open() as fh:
        rows = list(csv.DictReader(fh))
    x1 = {(r["case"], r["kind"], r["sheet_radius"], r["k"]): r for r in rows if r["grid"] == "x1"}
    coarse = {
        (r["case"], r["kind"], r["sheet_radius"], r["k"]): r for r in rows if r["grid"] == "x0.5"
    }
    for key, r in x1.items():
        if r["kind"] == "point" and float(r["k"]) == 10.0:
            assert float(r["mass_ratio_carbon_frozen"]) < 1.5
        if float(r["k"]) == 40.0:
            assert float(r["mass_ratio"]) > 2.0
        assert float(coarse[key]["impulse_per_mv"]) <= float(r["impulse_per_mv"])
        # The impactor's own momentum is delivered, plus the expansion's.
        assert float(r["impulse_per_mv"]) > 2.0


def test_mission_basis_mass_ratio_two() -> None:
    """The recorded mission basis (ADR-0060 amendment): mass ratio 2 at 100 km/s is met at
    k ~ 24 / 40 / 52 on the bowl / skirted bowl / deep paraboloid, inside the simulated range, and
    a deeper periapsis lowers the boost while multiplying the sunlight as 1/r^2."""
    from puffsat.xray_collision import mission

    ks = [mission.k_for_mass_ratio(p, 2.0).k for p in mission.PLATES]
    assert ks == pytest.approx([24.3, 40.4, 51.7], abs=1.0)
    assert all(mission.k_for_mass_ratio(p, 2.0).in_simulated_range for p in mission.PLATES)
    parker, deep = mission.periapsis(9.86), mission.periapsis(4.0)
    assert parker.dv_100 == pytest.approx(20.76e3, rel=0.01)
    assert deep.dv_100 < 0.7 * parker.dv_100
    assert deep.suns / parker.suns == pytest.approx((9.86 / 4.0) ** 2)


def test_impactor_sheet_budget_closed_forms() -> None:
    """The sheet arithmetic quoted in Result 11 follows from its closed forms."""
    from puffsat.xray_collision import overtake

    b = overtake.sheet_budget(0.05)
    assert b.thickness["polyimide (Kapton)"] == pytest.approx(0.05 / 1420.0)
    # Two-face radiative balance: T scales as r^(-1/2).
    assert b.temp_half_au / b.temp_one_au == pytest.approx(2.0**0.5)
    assert b.sail_accel_half_au[1] == pytest.approx(2.0 * b.sail_accel_half_au[0])
    assert b.spread_tilt < 0.05 * overtake.PULSE_TIME
    assert b.coast_drift[0] > 1.0e7  # tens of thousands of km: tracking is required
