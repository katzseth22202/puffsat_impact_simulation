"""Limits the blowdown, cycle, zone-load and cushion models must reproduce (hydrocode asks
2026-10-12)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import cushion as cu
from puffsat.walled_nozzle import cycle
from puffsat.walled_nozzle import near_term as nt
from puffsat.walled_nozzle import wall_zones as wz
from puffsat.walled_nozzle import zone_loads as zl


def test_contours_hold_the_sweeps_volumes() -> None:
    """The ported contour reproduces `euler2d::vessel`: 40 m^3 solved, 106 with the r 3.0 m
    bulge, and x2 at 5 kg by cube-root scaling."""
    a = bd.CH4_THROAT_5KG / bd.K5**2
    assert bd.long_cone(40.0, a).volume == pytest.approx(40.0, rel=1e-3)
    bulged = bd.long_cone(40.0, a, bulge=True)
    assert bulged.volume == pytest.approx(106.0, rel=0.01)
    assert bulged.scaled(bd.K5).volume == pytest.approx(2.0 * bulged.volume, rel=1e-3)
    assert max(bulged.r_wall(z) for z in np.linspace(0, bulged.length, 400)) == pytest.approx(3.0)
    # The sphere's head opens to the 7.5 cm port and its nose closes on the 25 cm throat, not on
    # the axis: 20.7 m^3, as the sweep's own sphere reports.
    assert bd.sphere(20.0, 0.198).volume == pytest.approx(20.72, rel=2e-3)


def test_closed_form_is_the_isentropic_gamma_law() -> None:
    """`closed_form_time` inverts `m/m0 = (1 + q lambda t)^(-1/q)`; at small t it is the
    linear drain `t = (1 - m/m0) / lambda`."""
    g, a, v, at = 1.2, 4000.0, 100.0, 1.0
    lam = at / v * a * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    assert bd.closed_form_time(g, a, v, at, 0.999) == pytest.approx(0.001 / lam, rel=2e-3)
    t = bd.closed_form_time(g, a, v, at, 0.05)
    assert (1.0 + 0.1 * lam * t) ** (-10.0) == pytest.approx(0.05)


def test_choked_flux_matches_the_exact_real_gas_sonic_flux() -> None:
    """The local-gamma_s formula against `near_term.nozzle_exit` on the real isentrope, within
    the 1% the module reports."""
    sol = bd.solve(bd.methane_specs()[0])
    br = nt.branch("equilibrium", sol)
    tab = nt.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=0.05, steps=120)
    temps, ratios = bd.calibration(br, tab, every=20)
    assert temps.size >= 3
    assert np.all(np.abs(ratios - 1.0) < 0.01)


def test_fixed_charge_solve_returns_the_charge() -> None:
    """The fixed-charge bisection lands on the charge it was asked for, and the 212 m^3 chamber
    is the 7000 K one."""
    for spec in bd.methane_specs():
        sol = bd.solve(spec)
        assert sol.charge_mass == pytest.approx(bd.CH4_CHARGE_25, rel=1e-4)
    assert bd.solve(bd.methane_specs()[0]).pairing.temp == pytest.approx(7000.0, rel=2e-3)


def test_time_to_interpolates() -> None:
    rows = [bd.Row(t, 1.0 - t, 1.0, 1.0 - t, 0.0, 0.0, 0.0, 0.0) for t in (0.0, 0.5, 1.0)]
    assert bd.time_to(rows, lambda r: r.mass_fraction, 0.25) == pytest.approx(0.75)
    assert math.isnan(bd.time_to(rows, lambda r: r.mass_fraction, -1.0))


def test_fill_energy_inverts() -> None:
    for gas, temp in (("methane_7000", 420.0), ("methane_7000", 1100.0), ("hydrogen_5500", 77.0)):
        e = cycle.fill_energy(gas, 0.5, temp)
        assert cycle.temperature_for(gas, 0.5, e) == pytest.approx(temp, rel=1e-6)


def test_methane_fill_heat_capacity_is_methanes() -> None:
    """Between 400 and 450 K pure CH4 has c_v ~ 1.9-2.0 kJ/kg/K (c_p 2.4-2.5 less R = 0.52)."""
    cv = (
        cycle.fill_energy("methane_7000", 1.0, 450.0)
        - cycle.fill_energy("methane_7000", 1.0, 400.0)
    ) / 50.0
    assert 1.8e3 < cv < 2.1e3


def _history(pressure: np.ndarray, dt: float, qsp: float) -> wz.History:
    n = pressure.shape[1]
    z = np.linspace(0.1, 2.0, n)
    return wz.History(
        z,
        np.full(n, 1.4),
        np.arange(pressure.shape[0], dtype=np.float64) * dt,
        pressure,
        100.0,
        1.0,
        qsp,
    )


def test_spike_impulse_and_width() -> None:
    """A 10 µs top-hat of 1 GPa over the QSP: 10 kPa s above the QSP, 10 µs wide."""
    dt, qsp = 1e-6, 10e6
    p = np.full((400, 1), qsp)
    p[200:210, 0] += 1e9
    impulse, width = zl._spike(_history(p, dt, qsp), 0)
    assert impulse == pytest.approx(1e4, rel=1e-9)
    assert width == pytest.approx(10e-6)


def _drive(p_t: np.ndarray, dt: float, qsp: float, n: int = 3) -> cu.Drive:
    pressure = np.repeat(p_t[:, None], n, axis=1)
    h = _history(pressure, dt, qsp)
    return cu.Drive(
        "test", h.z, h.r, ["head"] * n, h.time, pressure, qsp, float(h.time[-1]), h, 1e9
    )


def test_tile_stays_on_its_stops_below_the_precharge() -> None:
    qsp = 10e6
    d = _drive(np.full(2000, 0.8 * qsp), 1e-5, qsp)
    m = cu.move(d, cu.Design(10.0, 6.0, 0.85, 0.0))
    assert np.all(m.x_max == 0.0)
    assert np.all(m.back_max == pytest.approx(0.85 * qsp))


def test_small_impulse_rings_at_the_cushion_frequency() -> None:
    """A tile at its precharge, kicked by a small impulse I: undamped, it travels I / (mu w)
    with `w^2 = gamma p0 / (d mu)`, the linear gas spring."""
    qsp, dt = 10e6, 1e-6
    design = cu.Design(10.0, 8.0, 1.0, 0.0)
    p = np.full(20_000, qsp)
    p[10] += 50e6  # 50 MPa for 1 µs: 50 Pa s
    m = cu.move(_drive(p, dt, qsp), design)
    omega = math.sqrt(cu.GAMMA_K * qsp / (0.08 * design.mu))
    assert m.x_max[0] == pytest.approx(50.0 / (design.mu * omega), rel=0.03)


def test_backwall_overshoot_of_a_slow_load_is_one() -> None:
    """A backwall load that rises far slower than the ring's period is static: D ~ 1."""
    qsp = 10e6
    t = np.arange(0.0, 0.2, 1e-4)
    back = np.minimum(t / 0.1, 1.0) * qsp
    m = cu.Motion(
        np.zeros(1),
        np.zeros(1),
        np.zeros(1),
        back[None, :].max(axis=1),
        back[:, None],
        t,
        np.zeros(1, dtype=bool),
        np.zeros(1),
        0.0,
        0.0,
        np.zeros(1),
    )
    d = _drive(np.full(10, qsp), 1e-3, qsp, n=1)
    assert cu.backwall_overshoot(d, m, 0) == pytest.approx(1.0, abs=0.02)


def test_uniform_spike_does_not_bend_a_tile() -> None:
    """A spike uniform along the tile drives only its rigid translation: no bending."""
    dt, qsp, n = 1e-6, 10e6, 60
    p = np.full((300, n), qsp)
    p[100:110, :] += 1e9
    h = wz.History(
        np.linspace(0.0, 0.59, n),
        np.full(n, 1.4),
        np.arange(300, dtype=np.float64) * dt,
        p,
        100.0,
        1.0,
        qsp,
    )
    d = cu.Drive("uniform", h.z, h.r, ["head"] * n, h.time, p, qsp, float(h.time[-1]), h, 1e9)
    assert cu.bending_stress(d, cu.Design(10.0, 6.0, 0.85, 0.0), "head") < 1e3
