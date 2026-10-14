"""The chamber's blowdown, real gas and wall heat: hydrocode asks B1 and B6 (2026-10-12).

`aim_is_all_you_need` (ADR 0053) times the cycle with a closed-form quasi-steady adiabatic
blowdown: a choked throat, discharge coefficient one, an effective gamma. This module solves the
same 0-D blowdown on the real equilibrium gas and adds the wall's heat loss:

    dM/dt = -C A* F*(rho, e),        dU/dt = -C A* F*(rho, e) h(rho, e) - A_w q(rho, T, p),

with `U = M e` and `h = e + p/rho` the stagnation enthalpy of gas at rest in the chamber.

- **Choked flux.** `F* = sqrt(gamma_s rho p) (2/(gamma_s+1))^((gamma_s+1)/(2(gamma_s-1)))`
  with the equilibrium isentropic exponent `gamma_s = rho a^2 / p`. Checked against the exact
  real-gas sonic flux on the adiabatic isentrope (`near_term.nozzle_exit`); `calibration`
  reports the ratio, and the integration multiplies by it.
- **Wall heat.** `near_term.wall_heat`'s bracket, applied along the trajectory. Low edge: the H-
  in-band radiation floor plus frozen-`c_p` Bartz convection. High edge: an opaque blackbody
  plus equilibrium-`c_p` Bartz. Convection goes as `p^0.8 (T - T_w)` from its peak value.
- **Discharge coefficient `C`.** 1 here. The whole blowdown was run in 2-D
  (`sweep --vessel-blowdown`) to measure it, but the ghost-cell wall does not conserve mass over
  tens of milliseconds (`early_outflow`), so the 2-D run gives only the blast's head start in the
  first milliseconds.
- **Recombination.** The gas is in equilibrium. `equilibrium_margin` checks that `H + H + M`
  keeps up as the chamber empties, in decades of margin against the rate's stated uncertainty.

The 5 kg chambers are solved on the 2.5 kg rod's intensive state. A 5 kg chamber at `V` holds
the same gas as a 2.5 kg chamber at `V/2` (same charge per rod kilogram, same energy per
volume); the clock uses the real `V`, `A*` and wall area.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat import eos_methane
from puffsat.walled_nozzle import chamber, rates
from puffsat.walled_nozzle import near_term as nt

OUTPUT_DIR = nt.OUTPUT_DIR / "blowdown"
BLOWDOWN_2D = nt.OUTPUT_DIR / "vessel_blowdown.jsonl"

R_CH4 = eos_methane.K_B / eos_methane.M_CH4
R_H2 = eos_methane.K_B / (2.0 * eos_methane.M_H)
FILL_TEMP_CH4 = 420.0
"""The warm methane charge, through the regen passages (parent ADR 0053) [K]."""
FILL_TEMP_H2 = 77.0
"""The cold hydrogen gas charge (parent ADR 0053) [K]."""
GATE_THRUST_COST = 0.01
"""Share of exhaust momentum lost at the port (parent `chamber_isp.GATE_THRUST_COST`)."""
SIGMA_SB = 5.670374419e-8


# ---- Chamber geometry: `euler2d::vessel` ported, for wall area and zones ---------------------


@dataclass(frozen=True)
class Contour:
    """The `euler2d::vessel::Vessel` contour (domed head, cylinder, straight cone, bulge)."""

    r_c: float
    z_c: float
    z_hi: float
    r_t: float
    h_head: float = 0.7
    r_port: float = 0.075
    #: `(dr, z0, z1, ramp, ramp_out)`; `dr = 0` for none.
    bulge: tuple[float, float, float, float, float] = (0.0, 0.0, 0.0, 1.0, 1.0)
    curve: float = 0.0
    scale: float = 1.0
    """Every length times this: `2^(1/3)` for the 5 kg chamber."""

    def r_wall(self, z: float) -> float:
        z = z / self.scale
        if self.h_head > 0.0 and z < self.h_head:
            u = 1.0 - max(z, 0.0) / self.h_head
            r = self.r_port + (self.r_c - self.r_port) * math.sqrt(1.0 - u * u)
        elif z <= self.z_c:
            r = self.r_c
        else:
            s = min(max((z - self.z_c) / (self.z_hi - self.z_c), 0.0), 1.0)
            g = (1.0 - self.curve) * (1.0 - s) + self.curve * math.sqrt(1.0 - s * s)
            r = self.r_t + (self.r_c - self.r_t) * g
        dr, z0, z1, ramp, ramp_out = self.bulge
        if dr > 0.0:
            if z < z0:
                x = 1.0 - (z0 - z) / ramp
            elif z > z1:
                x = 1.0 - (z - z1) / ramp_out
            else:
                x = 1.0
            x = min(max(x, 0.0), 1.0)
            r += dr * 0.5 * (1.0 - math.cos(math.pi * x))
        return r * self.scale

    @property
    def length(self) -> float:
        return self.z_hi * self.scale

    def _samples(self, n: int = 4000) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        z = (np.arange(n) + 0.5) * self.length / n
        return z, np.array([self.r_wall(float(x)) for x in z])

    @property
    def volume(self) -> float:
        z, r = self._samples()
        return float(np.sum(math.pi * r * r) * (z[1] - z[0]))

    def area(self, z_lo: float = 0.0, z_hi: float = math.inf) -> float:
        """Wall area of revolution between two stations [m^2] (port face excluded)."""
        z, r = self._samples()
        dz = float(z[1] - z[0])
        slope = np.gradient(r, dz)
        keep = (z >= z_lo) & (z < z_hi)
        return float(np.sum(2.0 * math.pi * r[keep] * np.sqrt(1.0 + slope[keep] ** 2)) * dz)

    def scaled(self, factor: float) -> Contour:
        return dataclasses.replace(self, scale=self.scale * factor)


def long_cone(volume: float, throat_area: float, bulge: bool = False) -> Contour:
    """The sweep's `vessel_sized(1.4, cone, 0, 0.7, volume, r_t)` [with the r 3.0 m bulge].

    The cylinder is solved for `volume` before the bulge is added, as the sweep does (its r 3.0 m
    chamber is the 40 m^3 long cone plus the bulge).
    """
    r_t = math.sqrt(throat_area / math.pi)
    cone = (1.4 - r_t) / math.tan(math.radians(12.0))

    def make(length: float) -> Contour:
        return Contour(r_c=1.4, z_c=0.7 + length, z_hi=0.7 + length + cone, r_t=r_t)

    lo, hi = 0.0, 40.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if make(mid).volume < volume:
            lo = mid
        else:
            hi = mid
    c = make(0.5 * (lo + hi))
    if bulge:
        c = dataclasses.replace(c, bulge=(1.6, 1.4, 2.8, 0.5, 3.0))
    return c


def sphere(volume: float, throat_area: float) -> Contour:
    """The parent's sphere: a head and an elliptical nose each one radius deep."""
    r_s = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    return Contour(
        r_c=r_s,
        z_c=r_s,
        z_hi=2.0 * r_s,
        r_t=math.sqrt(throat_area / math.pi),
        h_head=r_s,
        curve=1.0,
    )


# ---- The chambers the asks name -------------------------------------------------------------


@dataclass(frozen=True)
class Spec:
    """One chamber: its gas, rod, real volume and throat, and how its charge is fixed."""

    name: str
    pairing: str
    #: Rod per pulse [kg]; the gas state is solved on the 2.5 kg rod at `volume * 2.5 / rod`.
    rod: float
    volume: float
    throat_area: float
    contour: Contour
    #: Charge per 2.5 kg of rod [kg] (fixed charge); `None` to fix the temperature instead.
    charge_25: float | None = None
    fill_temp: float = FILL_TEMP_CH4

    @property
    def energy_ratio(self) -> float:
        return self.rod / nt.ROD_MASS

    @property
    def proxy_volume(self) -> float:
        return self.volume / self.energy_ratio


CH4_THROAT_5KG = 1.33
"""The survivable chamber's throat (parent `tab:one_chamber`) [m^2]."""
CH4_CHARGE_25 = 54.22
"""Methane per 2.5 kg of rod in the 212 m^3 / 5 kg chamber at 7000 K (`near_term` at 106 m^3)."""
K5 = 2.0 ** (1.0 / 3.0)


def methane_specs() -> list[Spec]:
    """B1/B3/B6: the bulged 212 m^3 and unbulged 106 and 53 m^3, one throat, one charge."""
    a25 = CH4_THROAT_5KG / K5**2
    return [
        Spec(
            "CH4 212 m^3 bulged",
            "methane_7000",
            5.0,
            212.0,
            CH4_THROAT_5KG,
            long_cone(40.0, a25, bulge=True).scaled(K5),
            CH4_CHARGE_25,
        ),
        Spec(
            "CH4 106 m^3 unbulged",
            "methane_7000",
            5.0,
            106.0,
            CH4_THROAT_5KG,
            long_cone(53.0, a25).scaled(K5),
            CH4_CHARGE_25,
        ),
        Spec(
            "CH4 53 m^3 unbulged",
            "methane_7000",
            5.0,
            53.0,
            CH4_THROAT_5KG,
            long_cone(26.5, a25).scaled(K5),
            CH4_CHARGE_25,
        ),
    ]


H2_CHARGE_25 = 59.913
"""The parent's hydrogen charge: k = 24 on the 2.5 kg rod, 5500 K at 20 m^3."""


COMPACT = (
    # (label as in `sweep --vessel-blast --compact`, cylinder change [m], cone [deg], taper [m],
    #  throat over the 0.838 m^2 one), all at the 2.5 kg rod.
    ("compact: 2 m taper", 0.0, 12.0, 2.0, 1.0),
    ("compact: cylinder 1.5 m shorter", -1.5, 12.0, 3.0, 1.0),
    ("compact: 30 deg cone", 0.0, 30.0, 3.0, 1.0),
    ("compact: throat x2", 0.0, 12.0, 3.0, 2.0),
    ("compact: all (2 m taper, cylinder -1.5 m, 20 deg cone, throat x2)", -1.5, 20.0, 2.0, 2.0),
)


def compact_contour(dz_c: float, cone_deg: float, taper: float, throat_area: float) -> Contour:
    """The sweep's `--compact` chamber at the 2.5 kg rod: the r 3.0 m bulge on a trimmed base."""
    a25 = CH4_THROAT_5KG / K5**2
    base = long_cone(40.0, a25)
    r_t = math.sqrt(throat_area / math.pi)
    z_c = base.z_c + dz_c
    return Contour(
        r_c=1.4,
        z_c=z_c,
        z_hi=z_c + (1.4 - r_t) / math.tan(math.radians(cone_deg)),
        r_t=r_t,
        bulge=(1.6, 1.4, 2.8, 0.5, taper),
    )


def compact_specs(pairing: str = "methane_7000") -> list[Spec]:
    """The compact bulged chambers at 5 kg (methane) or 2.5 kg (hydrogen), fixed charge."""
    a25 = CH4_THROAT_5KG / K5**2
    out = []
    for label, dz, cone, taper, k in COMPACT:
        c = compact_contour(dz, cone, taper, k * a25)
        if pairing.startswith("hydrogen"):
            out.append(
                Spec(f"H2 {label}", pairing, 2.5, c.volume, k * a25, c, H2_CHARGE_25, FILL_TEMP_H2)
            )
        else:
            out.append(
                Spec(
                    f"CH4 {label}",
                    pairing,
                    5.0,
                    2.0 * c.volume,
                    k * CH4_THROAT_5KG,
                    c.scaled(K5),
                    CH4_CHARGE_25,
                )
            )
    return out


def hydrogen_specs(throat_bulged: float = 0.198) -> list[Spec]:
    """B5: the 20 m^3 sphere, and the r 3.0 m bulge (106 m^3 at 2.5 kg) at the same charge."""
    return [
        Spec(
            "H2 20 m^3 sphere",
            "hydrogen_5500",
            2.5,
            20.0,
            0.198,
            sphere(20.0, 0.198),
            H2_CHARGE_25,
            FILL_TEMP_H2,
        ),
        Spec(
            f"H2 106 m^3 bulged, A* {throat_bulged:.3g}",
            "hydrogen_5500",
            2.5,
            106.0,
            throat_bulged,
            long_cone(40.0, throat_bulged, bulge=True),
            H2_CHARGE_25,
            FILL_TEMP_H2,
        ),
    ]


# ---- The gas --------------------------------------------------------------------------------


def pairing_named(name: str) -> nt.Pairing:
    return next(p for p in nt.PAIRINGS if p.name == name)


def solve(spec: Spec) -> nt.ChamberSolution:
    """The peak chamber state on the 2.5 kg rod at the proxy volume.

    Fixed charge: the temperature is found by bisection so that the rod's energy heats the
    charge (with rod and plug) to it, then `near_term.solve_charge` is run at that temperature
    and returns the same charge.
    """
    pairing = pairing_named(spec.pairing)
    if spec.charge_25 is None:
        return nt.solve_charge(pairing, "janaf", spec.proxy_volume)
    feed = nt.feed_for(pairing, spec.charge_25)
    total = spec.charge_25 + nt.ROD_MASS + nt.PLUG_MASS
    rho, target = total / spec.proxy_volume, nt.PULSE_ENERGY / total
    lo, hi = 2000.0, 20000.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if eos_methane.pressure_energy(rho, mid, feed)[1] < target:
            lo = mid
        else:
            hi = mid
    fixed = dataclasses.replace(pairing, temp=0.5 * (lo + hi), paper_charge=spec.charge_25)
    return nt.solve_charge(fixed, "janaf", spec.proxy_volume)


def gamma_s(br: nt.Branch, rho: float, temp: float) -> tuple[float, float, float]:
    """`(gamma_s, p, a)`: the equilibrium isentropic exponent `rho a^2 / p`."""
    p = br(rho, temp)[0]
    a = eos_methane.sound_speed(rho, temp, br.feed)
    return rho * a * a / p, p, a


def choked_flux(br: nt.Branch, rho: float, temp: float) -> float:
    """Sonic mass flux [kg/m^2/s] of gas at rest at `(rho, T)`, on its local `gamma_s`."""
    g, p, _ = gamma_s(br, rho, temp)
    return float(math.sqrt(g * rho * p) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0))))


def calibration(br: nt.Branch, tab: nt.Isentrope, every: int = 25) -> tuple[Vec, Vec]:
    """`(T, exact / formula)` along the adiabatic isentrope, down to 0.1% of the density."""
    temps, ratios = [], []
    for j in range(0, len(tab.rho), every):
        if tab.rho[j] < 1e-3 * tab.rho[0]:
            break
        ex = nt.nozzle_exit(tab, j, 2.0)
        if ex is None:
            break
        temps.append(float(tab.temp[j]))
        ratios.append(ex.flux_star / choked_flux(br, float(tab.rho[j]), float(tab.temp[j])))
    order = np.argsort(temps)
    return np.asarray(temps)[order], np.asarray(ratios)[order]


Vec = NDArray[np.float64]


# ---- Wall heat ------------------------------------------------------------------------------


@dataclass(frozen=True)
class HeatModel:
    """`q(rho, T, p)` [W/m^2] on one edge of `near_term`'s bracket."""

    edge: str  # "none", "low", "high"
    conv_peak: float
    p0: float
    t0: float
    beam_length: float
    wall_temp: float = nt.PITCH_SURFACE_TEMPERATURE

    def __call__(self, br: nt.Branch, rho: float, temp: float, p: float) -> float:
        if self.edge == "none":
            return 0.0
        conv = (
            self.conv_peak
            * (p / self.p0) ** 0.8
            * max(0.0, temp - self.wall_temp)
            / (self.t0 - self.wall_temp)
        )
        if self.edge == "high":
            return float(conv + SIGMA_SB * temp**4)
        comp, _ = br.state(rho, temp)
        tau = nt.optical_depth_floor(comp, temp, self.beam_length)
        return float(conv + SIGMA_SB * temp**4 * nt.band_fraction(temp) * (1.0 - math.exp(-tau)))


def heat_model(edge: str, spec: Spec, sol: nt.ChamberSolution, tab: nt.Isentrope) -> HeatModel:
    # Bartz on the proxy's throat: its `d*^-0.2` makes the 5 kg chamber's convection ~5% lower.
    lo, hi = nt.bartz_liner(
        sol, tab, spec.throat_area / spec.energy_ratio ** (2 / 3), nt.PITCH_SURFACE_TEMPERATURE
    )
    return HeatModel(
        edge,
        hi if edge == "high" else lo,
        sol.pressure,
        sol.pairing.temp,
        4.0 * spec.volume / spec.contour.area(),
    )


# ---- The blowdown ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Row:
    time: float
    mass_fraction: float
    temp: float
    pressure: float
    #: Specific energy and stagnation enthalpy [J/kg].
    energy: float
    enthalpy: float
    #: Wall heat flux [W/m^2] and the wall heat so far, as a share of the pulse.
    q_wall: float
    wall_share: float


#: A wall heat flux `q(branch, rho, T, p)` [W/m^2]: `HeatModel`, or any law of the same shape.
HeatFlux = Callable[[nt.Branch, float, float, float], float]


def integrate(
    spec: Spec,
    sol: nt.ChamberSolution,
    br: nt.Branch,
    heat: HeatFlux,
    cal: tuple[Vec, Vec],
    discharge: float = 1.0,
    floor: float = 1.0e-3,
    step: float = 0.01,
) -> list[Row]:
    """March the 0-D blowdown until `floor` of the mass is left. Heun's method in time, the step
    a fixed share `step` of the mass e-fold time."""
    volume, area_w = spec.volume, spec.contour.area()
    m0 = sol.rho * volume
    pulse = nt.PULSE_ENERGY * spec.energy_ratio
    m, u = m0, sol.rho * volume * sol.energy
    temp = sol.pairing.temp

    def rates_at(m: float, u: float, guess: float) -> tuple[float, float, float, float, float]:
        rho, e = m / volume, u / m
        t = nt._temperature_near(rho, e, br, guess)
        flux = choked_flux(br, rho, t) * float(np.interp(t, cal[0], cal[1]))
        p = br(rho, t)[0]
        mdot = discharge * spec.throat_area * flux
        q = heat(br, rho, t, p)
        return -mdot, -mdot * (e + p / rho) - area_w * q, t, p, q

    time, lost = 0.0, 0.0
    dm, du, temp, p, q = rates_at(m, u, temp)
    rows = [Row(0.0, 1.0, temp, p, u / m, u / m + p * volume / m, q, 0.0)]
    while m > floor * m0:
        dt = step * m / -dm
        m1, u1 = m + dt * dm, u + dt * du
        dm1, du1, t1, _, q1 = rates_at(m1, u1, temp)
        m_new = m + 0.5 * dt * (dm + dm1)
        u_new = u + 0.5 * dt * (du + du1)
        lost += 0.5 * dt * area_w * (q + q1)
        m, u, time = m_new, u_new, time + dt
        dm, du, temp, p, q = rates_at(m, u, t1)
        rows.append(Row(time, m / m0, temp, p, u / m, u / m + p * volume / m, q, lost / pulse))
    return rows


def time_to(rows: Sequence[Row], key: Callable[[Row], float], level: float) -> float:
    """First time `key(row)` falls to `level`, interpolated [s]; NaN if it never does."""
    for a, b in pairwise(rows):
        if key(b) <= level < key(a) or key(b) <= level == key(a):
            f = (key(a) - level) / (key(a) - key(b))
            return a.time + f * (b.time - a.time)
    return math.nan


def state_at(rows: Sequence[Row], t: float) -> Row:
    """The row nearest `t` from below."""
    best = rows[0]
    for r in rows:
        if r.time > t:
            break
        best = r
    return best


def closed_form_time(
    gamma: float, speed: float, volume: float, throat: float, left: float
) -> float:
    """The parent's `chamber_isp.blowdown_time` [s]: `(left^-q - 1)/(q lambda)`."""
    q = 0.5 * (gamma - 1.0)
    choke = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    lam = throat / volume * speed * choke
    return float((left ** (-q) - 1.0) / (q * lam))


def fill_pressure(spec: Spec, sol: nt.ChamberSolution) -> float:
    """The next charge, cold, in the chamber [Pa]: methane at 420 K or hydrogen at 77 K."""
    gas_r = R_H2 if spec.pairing.startswith("hydrogen") else R_CH4
    return sol.charge_mass * spec.energy_ratio * gas_r * spec.fill_temp / spec.volume


# ---- Residue ---------------------------------------------------------------------------------


def species(br: nt.Branch, rho: float, temp: float) -> dict[str, float]:
    """Mole fractions of the leading species at `(rho, T)`."""
    comp, _ = br.state(rho, temp)
    n = {"H": comp.n_h, "C": comp.n_c, "H+": comp.n_hp, "C+": sum(comp.n_c_ions), "e-": comp.n_e}
    for mol, x in zip(eos_methane.MOLECULES, comp.n_mol, strict=True):
        n[mol.name] = x
    total = sum(n.values())
    return {k: v / total for k, v in sorted(n.items(), key=lambda kv: -kv[1])}


def equilibrium_margin(br: nt.Branch, rho: float, temp: float, residence: float) -> float:
    """Decades by which `H + H + M` recombination outruns the chamber's emptying.

    The recombination time of the free H is `1 / (2 k n_H n_M)`; it is compared with the mass
    e-fold time. Baulch's `H + H + H2` (`rates.H_H2`, the conservative fit, half a decade stated)
    stands for every third body. A margin above `rates.H_H2.decades()` is a verdict.
    """
    comp, _ = br.state(rho, temp)
    k = rates.H_H2(temp)
    tau = 1.0 / (2.0 * k * max(comp.n_h, 1.0) * comp.n_total)
    return math.log10(residence / tau)


# ---- Impulse per pulse ----------------------------------------------------------------------


@dataclass(frozen=True)
class Impulse:
    area_ratio: float
    eta_blowdown: float
    eta_net: float
    #: Share of the exhaust momentum the wall heat removes (low, high edge).
    heat_debit_low: float
    heat_debit_high: float
    #: Net impulse per pulse [N s], head-on debit taken, gate charged: adiabatic and with the
    #: wall heat at each edge.
    impulse: float
    impulse_low: float
    impulse_high: float


def exhaust_speeds(tab: nt.Isentrope, ar: float, every: int = 4) -> tuple[Vec, Vec, Vec]:
    """`(mass fraction, h0, exhaust speed)` of the parcel leaving at each isentrope station."""
    frac, h0, v = [], [], []
    for j in range(0, len(tab.rho), every):
        ex = nt.nozzle_exit(tab, j, ar)
        if ex is None:
            break
        frac.append(float(tab.rho[j] / tab.rho[0]))
        h0.append(float(tab.enthalpy[j]))
        v.append(ex.speed)
    return np.asarray(frac), np.asarray(h0), np.asarray(v)


def heat_debit(rows: Sequence[Row], speeds: tuple[Vec, Vec, Vec]) -> float:
    """Share of the adiabatic exhaust momentum lost to the wall.

    A parcel leaving at mass fraction `x` with stagnation enthalpy `h0` is given the adiabatic
    parcel's speed at the same `x`, scaled by `sqrt(h0 / h0_adiabatic)`. Parcels the adiabatic
    nozzle could not resolve are counted at zero, as in `near_term.efficiencies`.
    """
    frac, h_ad, v_ad = (a[::-1] for a in speeds)  # ascending in mass fraction for np.interp
    lossy = adiabatic = 0.0
    for a, b in pairwise(rows):
        dm = a.mass_fraction - b.mass_fraction
        x = 0.5 * (a.mass_fraction + b.mass_fraction)
        if x < frac[0]:
            continue
        h = 0.5 * (a.enthalpy + b.enthalpy)
        v = float(np.interp(x, frac, v_ad))
        ref = float(np.interp(x, frac, h_ad))
        adiabatic += dm * v
        lossy += dm * v * math.sqrt(max(h, 0.0) / ref)
    return 1.0 - lossy / adiabatic


def net_impulse(spec: Spec, sol: nt.ChamberSolution, speed: float) -> float:
    """`m_rod (1 - gate)(1 + k + P) u_e - m_rod w` [N s] for a mean exhaust speed `u_e`."""
    k = sol.charge_mass / nt.ROD_MASS
    expelled = 1.0 + k + nt.PLUG_RATIO
    return spec.rod * ((1.0 - GATE_THRUST_COST) * expelled * speed - chamber.CLOSING_SPEED)


# ---- One chamber, everything -----------------------------------------------------------------


@dataclass(frozen=True)
class Result:
    spec: Spec
    sol: nt.ChamberSolution
    wall_area: float
    p_fill: float
    gamma_peak: float
    calibration_span: tuple[float, float]
    #: `{edge: rows}` for "none", "low", "high".
    runs: dict[str, list[Row]]
    impulses: list[Impulse]
    margin_at_fill: float

    def t95(self, edge: str) -> float:
        return time_to(self.runs[edge], lambda r: r.mass_fraction, 0.05)

    def t_fill(self, edge: str) -> float:
        return time_to(self.runs[edge], lambda r: r.pressure, self.p_fill)

    def t_efold(self, edge: str) -> float:
        return time_to(self.runs[edge], lambda r: r.mass_fraction, math.exp(-1.0))


def run_spec(
    spec: Spec, area_ratios: Sequence[float] = (100.0, 300.0), report: nt.Report = print
) -> Result:
    sol = solve(spec)
    br = nt.branch("equilibrium", sol)
    tab = nt.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=1.0e-8, steps=900)
    cal = calibration(br, tab)
    runs = {
        edge: integrate(spec, sol, br, heat_model(edge, spec, sol, tab), cal)
        for edge in ("none", "low", "high")
    }
    p_fill = fill_pressure(spec, sol)
    effs = {e.area_ratio: e for e in nt.efficiencies(sol, br, tab, area_ratios)}
    impulses = []
    for ar in area_ratios:
        eff = effs[ar]
        # The freeze clock needs only the throat's size (the nozzle's length); the state is
        # intensive, so the real throat serves the proxy's gas.
        fc = dataclasses.replace(
            nt.freeze_check(sol, tab, spec.throat_area, ar), eta_blowdown=eff.eta_blowdown
        )
        eta_net = fc.eta_net
        speed = eff.speed_impulse * math.sqrt(eta_net / eff.eta_blowdown)
        speeds = exhaust_speeds(tab, ar)
        lo, hi = heat_debit(runs["low"], speeds), heat_debit(runs["high"], speeds)
        impulses.append(
            Impulse(
                ar,
                eff.eta_blowdown,
                eta_net,
                lo,
                hi,
                net_impulse(spec, sol, speed),
                net_impulse(spec, sol, speed * (1.0 - lo)),
                net_impulse(spec, sol, speed * (1.0 - hi)),
            )
        )
    fill_row = state_at(runs["none"], time_to(runs["none"], lambda r: r.pressure, p_fill))
    res = Result(
        spec,
        sol,
        spec.contour.area(),
        p_fill,
        gamma_s(br, sol.rho, sol.pairing.temp)[0],
        (float(cal[1].min()), float(cal[1].max())),
        runs,
        impulses,
        equilibrium_margin(
            br,
            fill_row.mass_fraction * sol.rho,
            fill_row.temp,
            max(time_to(runs["none"], lambda r: r.mass_fraction, math.exp(-1.0)), 1e-3),
        ),
    )
    report(summary_line(res))
    return res


def summary_line(r: Result) -> str:
    s = r.spec
    t95 = "/".join(f"{1e3 * r.t95(e):4.0f}" for e in ("none", "low", "high"))
    fill = "/".join(f"{1e3 * r.t_fill(e):4.0f}" for e in ("none", "low", "high"))
    return (
        f"{s.name:28s} V={s.volume:5.0f} A*={s.throat_area:5.3f} T0={r.sol.pairing.temp:5.0f} K "
        f"p0={r.sol.pressure / 1e5:5.0f} bar g_s={r.gamma_peak:.3f} "
        f"cal {r.calibration_span[0]:.3f}-{r.calibration_span[1]:.3f} | 95% out {t95} ms | "
        f"to {r.p_fill / 1e5:.2f} bar {fill} ms (adiabatic/low/high)"
    )


# ---- The 2-D check: the early outflow, and why only the early part ------------------------


@dataclass(frozen=True)
class EarlyOutflow:
    """The 2-D blowdown against the 0-D clock on the same gamma-law gas and pulse."""

    case: str
    volume: float
    throat_area: float
    time: float
    #: Share of the mass out through the throat in 2-D, and on the uniform 0-D clock.
    out_2d: float
    out_0d: float
    #: Mass the immersed wall has lost (negative) or made, as a share of the start.
    wall_leak: float


def early_outflow(
    path: Path = BLOWDOWN_2D, times: Sequence[float] = (2e-3, 5e-3, 10e-3, 20e-3)
) -> list[EarlyOutflow]:
    """Read `vessel-blowdown` runs at a few times.

    The 0-D reference starts uniform with the run's mass and total energy, so the difference is
    the blast's head start through the throat. `wall_leak = (inside + out - start) / start` is the
    ghost-cell wall's conservation error, which grows with time on curved and sloped walls: a few
    percent by 20 ms in the rocket chambers and a third of the charge by 30 ms in the sphere. The
    throat's outflow also stalls past ~10-20 ms. So the 2-D solver gives the first milliseconds
    only, and the clock is the 0-D one.
    """
    out = []
    for line in path.open():
        rec = json.loads(line)
        if not rec.get("throat_out_kg"):
            continue
        g, v, a = rec["gamma"], rec["volume_m3"], rec["throat_area_m2"]
        m0, e0 = rec["mass0_kg"], rec["energy0_j"]
        p0 = (g - 1.0) * e0 / v
        speed = math.sqrt(g * p0 / (m0 / v))
        q = 0.5 * (g - 1.0)
        lam = a / v * speed * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        t = np.asarray(rec["time_s"])
        inside = np.asarray(rec["mass_kg"])
        thrown = np.asarray(rec["throat_out_kg"])
        for when in times:
            if when > t[-1]:
                continue
            i = int(np.searchsorted(t, when))
            out.append(
                EarlyOutflow(
                    rec["case"],
                    v,
                    a,
                    float(t[i]),
                    float(thrown[i]) / m0,
                    1.0 - (1.0 + q * lam * float(t[i])) ** (-1.0 / q),
                    (float(inside[i]) + float(thrown[i]) - m0) / m0,
                )
            )
    return out


# ---- Output ---------------------------------------------------------------------------------


def write(results: Sequence[Result], out: Path = OUTPUT_DIR) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "summary.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "chamber",
                "volume_m3",
                "throat_m2",
                "wall_area_m2",
                "charge_kg",
                "temp_k",
                "pressure_bar",
                "gamma_s",
                "fill_bar",
                "edge",
                "t_efold_ms",
                "t95_ms",
                "t_fill_ms",
                "wall_heat_share",
                "residue_temp_k",
                "residue_mass_share",
                "residue_kg",
            ]
        )
        for r in results:
            for edge, rows in r.runs.items():
                tf = time_to(rows, lambda x: x.pressure, r.p_fill)
                res = state_at(rows, tf) if math.isfinite(tf) else rows[-1]
                w.writerow(
                    [
                        r.spec.name,
                        r.spec.volume,
                        r.spec.throat_area,
                        f"{r.wall_area:.1f}",
                        f"{r.sol.charge_mass * r.spec.energy_ratio:.1f}",
                        f"{r.sol.pairing.temp:.0f}",
                        f"{r.sol.pressure / 1e5:.1f}",
                        f"{r.gamma_peak:.3f}",
                        f"{r.p_fill / 1e5:.2f}",
                        edge,
                        f"{1e3 * time_to(rows, lambda x: x.mass_fraction, math.exp(-1)):.1f}",
                        f"{1e3 * time_to(rows, lambda x: x.mass_fraction, 0.05):.1f}",
                        f"{1e3 * tf:.1f}",
                        f"{rows[-1].wall_share:.4f}",
                        f"{res.temp:.0f}",
                        f"{res.mass_fraction:.4f}",
                        f"{res.mass_fraction * r.sol.total_mass * r.spec.energy_ratio:.2f}",
                    ]
                )
    with (out / "impulse.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "chamber",
                "volume_m3",
                "area_ratio",
                "eta_blowdown",
                "eta_net",
                "heat_debit_low",
                "heat_debit_high",
                "impulse_kns",
                "impulse_low_kns",
                "impulse_high_kns",
            ]
        )
        for r in results:
            for i in r.impulses:
                w.writerow(
                    [
                        r.spec.name,
                        r.spec.volume,
                        i.area_ratio,
                        f"{i.eta_blowdown:.4f}",
                        f"{i.eta_net:.4f}",
                        f"{i.heat_debit_low:.4f}",
                        f"{i.heat_debit_high:.4f}",
                        f"{i.impulse / 1e3:.0f}",
                        f"{i.impulse_low / 1e3:.0f}",
                        f"{i.impulse_high / 1e3:.0f}",
                    ]
                )
    for r in results:
        slug = r.spec.name.replace(" ", "_").replace("^", "").replace(",", "").replace("*", "")
        with (out / f"history_{slug}.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "edge",
                    "time_ms",
                    "mass_fraction",
                    "temp_k",
                    "pressure_bar",
                    "q_wall_mw_m2",
                    "wall_heat_share",
                ]
            )
            for edge, rows in r.runs.items():
                for row in rows[::4]:
                    w.writerow(
                        [
                            edge,
                            f"{1e3 * row.time:.2f}",
                            f"{row.mass_fraction:.5f}",
                            f"{row.temp:.0f}",
                            f"{row.pressure / 1e5:.3f}",
                            f"{row.q_wall / 1e6:.2f}",
                            f"{row.wall_share:.4f}",
                        ]
                    )


def main() -> None:
    results = [run_spec(s) for s in methane_specs() + hydrogen_specs()]
    write(results)
    for r in results:
        fill = state_at(r.runs["none"], r.t_fill("none"))
        mix = species(nt.branch("equilibrium", r.sol), fill.mass_fraction * r.sol.rho, fill.temp)
        top = ", ".join(f"{k} {v:.2f}" for k, v in list(mix.items())[:5])
        print(
            f"{r.spec.name:28s} residue at fill: {fill.temp:.0f} K, "
            f"{100 * fill.mass_fraction:.1f}% of the mass; {top}; "
            f"H+H+M margin {r.margin_at_fill:.1f} decades"
        )
        for i in r.impulses:
            print(
                f"    AR {i.area_ratio:3.0f}: eta_bd {i.eta_blowdown:.3f} eta_net {i.eta_net:.3f} "
                f"I {i.impulse / 1e3:.0f} kN s (wall heat {i.impulse_low / 1e3:.0f}-"
                f"{i.impulse_high / 1e3:.0f})"
            )
    quick = BLOWDOWN_2D.with_name("vessel_blowdown_quick.jsonl")
    if quick.exists():
        for o in early_outflow(quick):
            print(
                f"2-D {o.case:38s} {1e3 * o.time:4.0f} ms: out {o.out_2d:.3f} (0-D "
                f"{o.out_0d:.3f}), wall leak {o.wall_leak:+.3f}"
            )


if __name__ == "__main__":
    main()
