"""A light foam projectile stopped on a dense anvil inside a radiation case, and its thrust ledger.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.anvil

The collision (Result 8 of docs/xray_collision_study.md). A carbon-foam projectile of mass m at
`v_rel` stops on an anvil of mass M at rest. By Galilean invariance this is `fireball.piston` with
the roles swapped: the dense anvil "sweeps" the light foam, so the whole heat Q = 1/2 mu v^2 lands
in the foam, which stagnates on the anvil's face. Shares are quoted against the projectile's own
ship-frame kinetic energy 1/2 m v^2, not the flipped call's.

The plate stage (the parent's §pass_through_pusher). The case delivers a share `eta_case` of the
light into motion of the channel filler and the propellant plate. The filler blows off astern with
mass m_b; the propellant plate goes forward with m_f. Momentum balance gives the parent's
eq. pulse_split,

    p = sqrt(2 E m_f m_b / (m_f + m_b)),

and the propellant must strike the pusher plate no faster than the impingement cap `u_cap`. The
ship receives J = c p, where c is the parent's 0.8 ("the plate catches 80% of the forward
momentum", no bounce credit) or up to 2 f with a rebound (`f` as in this repo's f(v) study).
Specific impulse is J over everything the ship spends per pulse, g0 included.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.xray_collision.fireball import CollisionResult, piston
from puffsat.xray_collision.materials import CARBON

G0 = 9.80665
OUT_DIR = Path("data/results/xray_collision")
V_REL = 400.0e3  # [m/s], a 200 km/s periapsis met head-on
PROJECTILE_MASS = 10.0  # [kg], foam plus its C-C sheets
FOAM_RHO0 = 10.0  # [kg/m^3]
FOAM_DIAMETER = 1.0  # [m], threads the bore
ANVIL_MASS = 50.0  # [kg], the settled design (k = 5)
SHEET_MASS = 1.8  # [kg], 0.5 mm C-C at 1.8 g/cc over the sunward half of the side
GRAPHITE_RHO = 2200.0  # [kg/m^3]
IRON_RHO = 7874.0  # [kg/m^3]
STAGNATION = 1.5  # face pressure / rho0 v^2 (the f(v) study's reflected stagnation is ~1.2)
RT_BUBBLE = 0.05  # classical bubble-front constant, h = RT_BUBBLE A g t^2 with A = 1


@dataclass(frozen=True)
class FoamOnAnvil:
    foam_mass: float  # [kg] that radiates (projectile less its dense sheets)
    anvil_mass: float  # [kg]
    foam_rho0: float  # [kg/m^3]
    kappa_scale: float
    heat_share: float  # Q / projectile KE, over the radiating foam and the anvil
    kt_start_ev: float
    tau0: float
    f_rad: float  # of Q
    light_share: float  # light / projectile KE (the sheets' share counted as lost)
    hard_share: float  # light above 1 keV / projectile KE
    kt_emit_ev: float
    debris_speed: float  # [m/s], centre of mass of projectile + anvil, astern
    pulse_time: float  # [s], the foam cylinder's arrival time at v_rel
    rt_graphite: float  # [s], bubble burn-through of a graphite anvil body
    rt_iron: float  # [s], same areal mass in iron


def rt_burnthrough(anvil_mass: float, area: float, rho_anvil: float, pressure: float) -> float:
    """Time for Rayleigh-Taylor bubbles of the light foam plasma to eat through an anvil of
    areal mass sigma = M/A pushed at `pressure`: thickness d = sigma/rho, acceleration
    g = P/sigma, and h = RT_BUBBLE g t^2 reaches d at t = sqrt(d / (RT_BUBBLE g))."""
    sigma = anvil_mass / area
    return math.sqrt(sigma / rho_anvil / (RT_BUBBLE * pressure / sigma))


def foam_on_anvil(
    anvil_mass: float = ANVIL_MASS,
    foam_rho0: float = FOAM_RHO0,
    *,
    projectile_mass: float = PROJECTILE_MASS,
    sheet_mass: float = 0.0,
    v_rel: float = V_REL,
    kappa_scale: float = 1.0,
) -> FoamOnAnvil:
    """Foam (projectile less `sheet_mass` of dense C-C) stopped on the anvil. The sheets are dead
    mass: they collide dense-on-dense, stay opaque, and radiate ~0.4% (Result 4's graphite row)."""
    foam = projectile_mass - sheet_mass
    r: CollisionResult = piston(
        CARBON, anvil_mass, foam_rho0, foam / anvil_mass, v_rel, kappa_scale=kappa_scale
    )
    ke = 0.5 * projectile_mass * v_rel**2
    area = math.pi * FOAM_DIAMETER**2 / 4.0
    length = foam / (foam_rho0 * area)
    p_face = STAGNATION * foam_rho0 * v_rel**2
    return FoamOnAnvil(
        foam_mass=foam,
        anvil_mass=anvil_mass,
        foam_rho0=foam_rho0,
        kappa_scale=kappa_scale,
        heat_share=r.heat / ke,
        kt_start_ev=r.kt_start_ev,
        tau0=r.tau0,
        f_rad=r.f_rad,
        light_share=r.f_rad * r.heat / ke,
        hard_share=r.f_hard * r.heat / ke,
        kt_emit_ev=r.kt_emit_ev,
        debris_speed=projectile_mass * v_rel / (projectile_mass + anvil_mass),
        pulse_time=length / v_rel,
        rt_graphite=rt_burnthrough(anvil_mass, area, GRAPHITE_RHO, p_face),
        rt_iron=rt_burnthrough(anvil_mass, area, IRON_RHO, p_face),
    )


def case_delivery(albedo: float, a_filler: float, a_wall: float, a_hole: float) -> float:
    """Steady cavity balance: the share of the light absorbed by the filler, against the walls'
    net loss (1 - albedo) and an open hole, A_f / (A_f + (1 - alpha) A_w + A_h)."""
    return a_filler / (a_filler + (1.0 - albedo) * a_wall + a_hole)


def pulse_split(energy: float, m_f: float, m_b: float) -> float:
    """Momentum each way when `energy` is shared by a forward m_f and a backward m_b (eq.
    pulse_split): p = sqrt(2 E m_f m_b / (m_f + m_b))."""
    return math.sqrt(2.0 * energy * m_f * m_b / (m_f + m_b))


@dataclass(frozen=True)
class Split:
    m_forward: float  # [kg] propellant plate
    m_backward: float  # [kg] filler blow-off
    u_forward: float  # [m/s] propellant speed at the plate
    impulse: float  # [N s] on the ship
    spent: float  # [kg] fixed spent mass + filler + propellant
    isp: float  # [s]


def best_split(energy: float, fixed_spent: float, u_cap: float, catch: float) -> Split:
    """Filler and propellant masses that maximize specific impulse J / (spent g0), with
    J = catch * p and the propellant at u <= u_cap. Parametrized by (m_f, u): p = m_f u, the
    forward energy is m_f u^2 / 2, and the backward mass takes the rest at p^2 / (2 E_b)."""
    m_f = np.geomspace(1e-3, 1e3, 3000)[:, None] * energy / u_cap**2
    u = np.linspace(1e-3, 1.0, 1000)[None, :] * u_cap
    e_b = energy - 0.5 * m_f * u**2
    ok = e_b > 0.0
    p = m_f * u
    m_b = np.where(ok, p**2 / (2.0 * np.where(ok, e_b, 1.0)), np.inf)
    spent = fixed_spent + m_f + m_b
    isp = np.where(ok, catch * p / (spent * G0), 0.0)
    i, j = np.unravel_index(int(np.argmax(isp)), isp.shape)
    mf, uu = float(m_f[i, 0]), float(u[0, j])
    pp = mf * uu
    return Split(mf, float(m_b[i, j]), uu, catch * pp, float(spent[i, j]), float(isp[i, j]))


@dataclass(frozen=True)
class Ledger:
    """One pulse, per kg of projectile, ship frame."""

    walls: float  # [kg/kg] case walls and liner beyond the anvil
    eta_case: float
    u_cap: float  # [m/s]
    catch: float  # J / p_forward: 0.8 = parent, 1.6 = f(v) bounce at f = 0.8
    light: float  # [J/kg]
    delivered: float  # [J/kg] into filler + propellant motion
    forward_ke: float  # [J/kg] propellant plate KE toward the plate
    m_forward: float
    m_backward: float
    u_forward: float
    impulse: float  # [N s/kg]
    spent: float  # [kg/kg]
    isp: float  # [s]
    jet_efficiency: float  # (1/2 spent v_e^2) / projectile KE, v_e = J / spent


def ledger(
    collision: FoamOnAnvil,
    walls: float,
    eta_case: float,
    u_cap: float,
    catch: float,
    projectile_mass: float = PROJECTILE_MASS,
    v_rel: float = V_REL,
) -> Ledger:
    ke = 0.5 * v_rel**2  # per kg of projectile
    light = collision.light_share * ke
    delivered = eta_case * light
    fixed = collision.anvil_mass / projectile_mass + walls
    s = best_split(delivered, fixed, u_cap, catch)
    v_e = s.impulse / s.spent
    return Ledger(
        walls=walls,
        eta_case=eta_case,
        u_cap=u_cap,
        catch=catch,
        light=light,
        delivered=delivered,
        forward_ke=0.5 * s.m_forward * s.u_forward**2,
        m_forward=s.m_forward,
        m_backward=s.m_backward,
        u_forward=s.u_forward,
        impulse=s.impulse,
        spent=s.spent,
        isp=s.isp,
        jet_efficiency=0.5 * s.spent * v_e**2 / ke,
    )


ANVIL_MASSES = (20.0, 30.0, 40.0, 50.0, 60.0, 100.0)
FOAM_RHO0S = (1.0, 3.0, 10.0, 30.0, 100.0)
WALLS = (0.0, 1.5, 5.0, 10.0)  # [kg per kg of projectile]; 1.5 fills the parent's 6.5 cap
ETA_CASE = (0.5, 0.65, 0.8)
U_CAPS = (150.0e3, 100.0e3)
CATCHES = (0.8, 1.2, 1.6)


def _write(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"python: wrote {len(rows)} rows -> {path}")


def main() -> None:
    print("== foam on anvil, 400 km/s, 10 kg projectile (no sheets), opacity 0.1 / 1 / 10x")
    rows = []
    for m in ANVIL_MASSES:
        for rho in FOAM_RHO0S:
            for ks in (0.1, 1.0, 10.0):
                c = foam_on_anvil(m, rho, kappa_scale=ks)
                rows.append(asdict(c))
                if ks == 1.0:
                    print(
                        f"  anvil {m:5.0f} rho0 {rho:5g} heat {c.heat_share:.2f} "
                        f"kT0 {c.kt_start_ev:4.0f} tau0 {c.tau0:7.2g} f_rad {c.f_rad:.2f} "
                        f"light {c.light_share:.2f} >1keV {c.hard_share:.2f} "
                        f"pulse {1e6 * c.pulse_time:4.1f} us RT C {1e6 * c.rt_graphite:4.1f} us "
                        f"Fe {1e6 * c.rt_iron:4.1f} us"
                    )
    _write(OUT_DIR / "anvil_sweep.csv", rows)

    print("\n== thrust ledger: 50 kg anvil, sheets 1.8 kg, per kg of projectile")
    design = foam_on_anvil(sheet_mass=SHEET_MASS)
    print(f"  light share {design.light_share:.3f} of the projectile KE")
    led = []
    for walls in WALLS:
        for eta in ETA_CASE:
            for cap in U_CAPS:
                for catch in CATCHES:
                    g = ledger(design, walls, eta, cap, catch)
                    led.append(asdict(g))
                    if eta == 0.8 and cap == 150.0e3:
                        print(
                            f"  walls {walls:4.1f} catch {catch:.1f} Isp {g.isp:5.0f} s "
                            f"spent {g.spent:5.2f} kg J {g.impulse / 1e3:5.0f} kN s "
                            f"m_f {g.m_forward:.2f} m_b {g.m_backward:.2f} "
                            f"jet eff {g.jet_efficiency:.3f}"
                        )
    _write(OUT_DIR / "anvil_ledger.csv", led)


if __name__ == "__main__":
    main()
