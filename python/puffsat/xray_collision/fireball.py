"""One-zone collision fireball: shocked state, expansion, and the share of heat that radiates.

Kinematics (ship frame). An impactor of mass m at v_rel hits a target of k m at rest. The merged
mass keeps centre-of-mass speed v_rel/(1+k), which carries 1/(1+k) of the impactor's kinetic
energy away as bulk motion. Only

    Q = 1/2 mu v_rel^2,   mu = k m / (1+k)       (= k/(1+k) of the impactor's KE)

can become heat.

Shocked state. For one material on both sides the contact moves at v_rel/2, so the first-shocked
matter has particle speed u_p = v_rel/2 relative to its upstream state. Its density comes from the
Rankine-Hugoniot jumps, p = rho0 U_s u_p with U_s = u_p rho/(rho - rho0), with e - e0 = u_p^2/2
from the EOS. The one-zone fireball starts at that density, at rest in the CM frame, holding the
mass-averaged heat Q / ((1+k) m). For k = 1 that average equals u_p^2/2. For k > 1 the shock
decays into the target, and the average stands in for a graded heating.

Expansion. A homologous uniform sphere of the merged mass:

    E_kin = 3/10 M Rdot^2,  dE_kin/dt = p dV/dt  ->  Rddot = (20 pi / 3) R^2 p / M
    dE_int/dt = -p dV/dt - L

Luminosity. L = 4 pi R^2 sigma T^4 g, with

    g = 1 / (1 + 1/(4/3 kP rho R) + 3/4 kR rho D),   D = 2R

Thin, this is volume emission 4 kP rho sigma T^4 V. Thick, it is the parent's
sigma T^4 / (1 + 3/4 tau) (templateArxiv.tex eq. cooling_race). The spectrum is a blackbody at
the zone temperature, which overstates its hardness while the zone is opaque (the true
photosphere is cooler). That bias favours X-rays, so the X-ray shares here are upper bounds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat.xray_collision.eos import EV, K_B, EosTable, build_table
from puffsat.xray_collision.materials import Material

SIGMA_SB = 5.670374419e-8
V_REL_ROCKET = 67.0e3  # 56 km/s impactor + 11 km/s ship
T_STOP_EV = 0.5  # the TOPS/OPLIB floor; below it argon and iron are neutral and transparent
STEP_FRACTION = 0.02  # RK4 step as a fraction of min(hydro, radiative) time
EUV_EV = 30.0
XRAY_EV = 100.0


@dataclass(frozen=True)
class Kinematics:
    """Ship-frame energy ledger for an impactor `m_impactor` at `v_rel` on a resting target."""

    m_impactor: float
    m_target: float
    v_rel: float

    @property
    def m_total(self) -> float:
        return self.m_impactor + self.m_target

    @property
    def ke_ship_frame(self) -> float:
        return 0.5 * self.m_impactor * self.v_rel**2

    @property
    def heat(self) -> float:
        return 0.5 * self.m_impactor * self.m_target / self.m_total * self.v_rel**2

    @property
    def heat_fraction(self) -> float:
        """Share of the ship-frame kinetic energy available as heat, k/(1+k)."""
        return self.m_target / self.m_total


def planck_above(x: float) -> float:
    """Fraction of blackbody flux carried by photons above x = E/kT (series of 15/pi^4 ...)."""
    if x <= 0.0:
        return 1.0
    if x > 200.0:
        return 0.0
    total = 0.0
    for n in range(1, 80):
        total += math.exp(-n * x) * (x**3 / n + 3 * x**2 / n**2 + 6 * x / n**3 + 6 / n**4)
    return min(15.0 / math.pi**4 * total, 1.0)


@dataclass(frozen=True)
class ShockState:
    rho: float  # [kg/m^3]
    temp: float  # [K]
    p: float  # [Pa]
    zbar: float


def shocked_state(table: EosTable, u_p: float, rho0: float | None = None) -> ShockState:
    """Hugoniot point from rest at density `rho0` (default: the condensed material) and particle
    speed `u_p`, by bisection in rho. A porous start (foam, spray) is treated as fully collapsing
    onto the plasma EOS, with no strength or pore physics."""
    rho0 = table.material.rho0 if rho0 is None else rho0
    e = 0.5 * u_p**2

    def excess(rho: float) -> float:
        p_eos = table.state(rho, table.temp(rho, e))[0]
        return p_eos - rho0 * rho * u_p**2 / (rho - rho0)

    lo, hi = 1.001 * rho0, 40.0 * rho0
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if excess(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    rho = math.sqrt(lo * hi)
    temp = table.temp(rho, e)
    p, zbar, _, _ = table.state(rho, temp)
    return ShockState(rho, temp, p, zbar)


@dataclass(frozen=True)
class CollisionResult:
    """One collision. Energy shares `f_*` are fractions of the heat Q, not of the ship-frame KE."""

    material: str
    impactor_radius: float  # [m]
    mass_ratio: float  # target / impactor
    v_rel: float  # [m/s]
    lowering: bool
    kappa_scale: float
    rho0: float  # starting density of both bodies [kg/m^3]
    m_impactor: float  # [kg]
    heat: float  # Q [J]
    heat_fraction: float  # Q / ship-frame KE
    rho_shock: float
    kt_shock_ev: float  # first-shocked matter, e = u_p^2/2
    kt_start_ev: float  # one-zone start, e = Q/M (equal to kt_shock_ev for k = 1)
    zbar_shock: float
    r0: float  # merged fireball radius at the start [m]
    tau0: float  # Rosseland optical depth across it at the start
    race0: float  # parent's t_cool/t_exp at the start
    f_rad: float
    f_euv: float  # photons above EUV_EV
    f_xray: float  # photons above XRAY_EV
    kt_emit_ev: float  # luminosity-weighted zone kT
    v_exp: float  # final edge speed [m/s]
    closure: float  # (E_rad + E_kin + E_int) / Q, should be 1
    t_end: float  # [s]

    @property
    def light_of_ship_ke(self) -> float:
        """Radiated energy as a share of the impactor's ship-frame kinetic energy."""
        return self.f_rad * self.heat_fraction


def _derivs(
    y: tuple[float, ...], table: EosTable, mass: float, kappa_scale: float
) -> tuple[tuple[float, ...], float, float, float]:
    """d/dt of (R, Rdot, E_int, E_rad, E_euv, E_xray, int L kT); also (T, p, hydro time)."""
    r, rdot, e_int = y[0], y[1], y[2]
    rho = mass / (4.0 / 3.0 * math.pi * r**3)
    temp = table.temp(rho, max(e_int, 0.0) / mass)
    p, _, kr, kp = table.state(rho, temp)
    kr *= kappa_scale
    kp *= kappa_scale
    g = 1.0 / (1.0 + 1.0 / (4.0 / 3.0 * kp * rho * r) + 0.75 * kr * rho * 2.0 * r)
    lum = 4.0 * math.pi * r**2 * SIGMA_SB * temp**4 * g
    kt_ev = K_B * temp / EV
    d = (
        rdot,
        20.0 * math.pi / 3.0 * r**2 * p / mass,
        -p * 4.0 * math.pi * r**2 * rdot - lum,
        lum,
        lum * planck_above(EUV_EV / kt_ev),
        lum * planck_above(XRAY_EV / kt_ev),
        lum * kt_ev,
    )
    t_hydro = r / max(abs(rdot), math.sqrt(p / rho))
    return d, temp, lum, t_hydro


def collide(
    material: Material,
    impactor_radius: float,
    mass_ratio: float = 3.0,
    v_rel: float = V_REL_ROCKET,
    *,
    lowering: bool = False,
    kappa_scale: float = 1.0,
    rho0: float | None = None,
) -> CollisionResult:
    """Ball of `material` with radius `impactor_radius` at `v_rel` into `mass_ratio` times its mass
    at rest. Integrates the one-zone fireball (RK4) until it cools to T_STOP_EV.

    `rho0` overrides the starting density of both bodies (porous foam or spray); the default is
    the condensed material."""
    table = build_table(material, lowering)
    rho_start = material.rho0 if rho0 is None else rho0
    m_imp = 4.0 / 3.0 * math.pi * impactor_radius**3 * rho_start
    kin = Kinematics(m_imp, mass_ratio * m_imp, v_rel)
    mass = kin.m_total
    shock = shocked_state(table, 0.5 * v_rel, rho_start)
    r = (3.0 * mass / (4.0 * math.pi * shock.rho)) ** (1.0 / 3.0)
    q = kin.heat
    y: tuple[float, ...] = (r, 0.0, q, 0.0, 0.0, 0.0, 0.0)

    temp0 = table.temp(shock.rho, q / mass)
    p0, z0, kr0, _ = table.state(shock.rho, temp0)
    tau0 = kappa_scale * kr0 * shock.rho * 2.0 * r
    n_tot = (1.0 + z0) * shock.rho / (material.mass_amu * 1.66053906660e-27)
    flux0 = SIGMA_SB * temp0**4 / (1.0 + 0.75 * tau0)
    race0 = n_tot * K_B * temp0 * math.sqrt(5.0 / 3.0 * p0 / shock.rho) / (2.0 * flux0)

    t = 0.0
    while True:
        k1, temp, lum, t_hydro = _derivs(y, table, mass, kappa_scale)
        if temp < T_STOP_EV * EV / K_B or y[2] < 1e-6 * q:
            break
        dt = STEP_FRACTION * min(t_hydro, y[2] / max(lum, 1e-300))
        k2 = _derivs(
            tuple(a + 0.5 * dt * b for a, b in zip(y, k1, strict=True)), table, mass, kappa_scale
        )[0]
        k3 = _derivs(
            tuple(a + 0.5 * dt * b for a, b in zip(y, k2, strict=True)), table, mass, kappa_scale
        )[0]
        k4 = _derivs(
            tuple(a + dt * b for a, b in zip(y, k3, strict=True)), table, mass, kappa_scale
        )[0]
        y = tuple(
            a + dt / 6.0 * (b1 + 2.0 * b2 + 2.0 * b3 + b4)
            for a, b1, b2, b3, b4 in zip(y, k1, k2, k3, k4, strict=True)
        )
        t += dt

    e_kin = 0.3 * mass * y[1] ** 2
    return CollisionResult(
        material=material.name,
        impactor_radius=impactor_radius,
        mass_ratio=mass_ratio,
        v_rel=v_rel,
        lowering=lowering,
        kappa_scale=kappa_scale,
        rho0=rho_start,
        m_impactor=m_imp,
        heat=q,
        heat_fraction=kin.heat_fraction,
        rho_shock=shock.rho,
        kt_shock_ev=K_B * shock.temp / EV,
        kt_start_ev=K_B * temp0 / EV,
        zbar_shock=z0,
        r0=r,
        tau0=tau0,
        race0=race0,
        f_rad=y[3] / q,
        f_euv=y[4] / q,
        f_xray=y[5] / q,
        kt_emit_ev=y[6] / max(y[3], 1e-300),
        v_exp=y[1],
        closure=(y[3] + e_kin + y[2]) / q,
        t_end=t,
    )
