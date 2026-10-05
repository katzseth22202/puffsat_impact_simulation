"""A coated steel wall under the near-term methane pulse: 1-D transient conduction (A5 follow-up).

The paper puts 0.2 mm of pitch over Cr-Mo steel and rejects graphite because "it conducts too
well": a thin graphite layer passes the heat to the steel, which fails first. The design question
is whether ~1 mm of *sprayed* graphite over bare steel survives methane at 7000 K. The chemistry
is already answered: methane's boundary layer is carbon-saturated, so it plates carbon rather
than eating it (`near_term`, section 6.5 of the asks). This module answers the thermal half.

**Gas side, from the solved pulse.** The chamber state against time is `near_term.blowdown` on
the real-gas choked clock. Convection is Bartz, scaled over the pulse as `p^0.8` and driven by the
*actual* surface temperature, so a hot coating is shielded by its own temperature:

    q_conv = h0 (p/p0)^0.8 (T_gas - T_s)

`h0` is bracketed exactly as `near_term.bartz_liner` brackets it. Radiation is the exchange
between a grey wall and a gas of emissivity `eps_g`:

    q_rad = eps_w eps_g sigma (T_gas^4 - T_s^4)

`eps_g = 1` at the high edge (opaque). At the low edge `eps_g` is the H- in-band floor
`band_fraction * (1 - exp(-tau))`.

**Solid side.** A node-based finite volume, backward Euler, with a node on the coating/steel
interface and perfect thermal contact. That is the conservative choice for the steel: any
contact resistance would insulate it. The back face is adiabatic. Graphite's `c_p(T)` is the
NASA-Glenn fit. Its conductivity is **bracketed, because it is the whole answer**:

- ~3 W/m/K for a porous sprayed coat, or pyrolytic graphite across its planes;
- ~10 W/m/K for an intermediate case;
- ~40 W/m/K for dense isotropic graphite at temperature, which is the paper's "conducts too well".

**The surface is capped at carbon's melting point.** At 496 bar the chamber is above carbon's
triple point (107 atm, Haaland 1976), so a carbon surface melts rather than sublimes, and the
boundary gas is near saturation, which suppresses sublimation. Heat arriving beyond the cap melts
graphite at its heat of fusion. That is the *low* removal heat (~9 MJ/kg, against sublimation's
59), so the melted thickness is an upper edge. The melt point (~4800 K) and the heat of fusion
(~105 kJ/mol) are literature values that were not checked here.
"""

from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import blowing, near_term

Vec = NDArray[np.float64]

SIGMA_SB = 5.670374419e-8

# --- Materials ---------------------------------------------------------------------------------

#: Cr-Mo steel, as the paper takes it: k = 35 W/m/K, rho c chosen so sqrt(k rho c) = 11400.
STEEL_K = 35.0
STEEL_RHO = 7850.0
STEEL_C = 11400.0**2 / (STEEL_K * STEEL_RHO)
#: Steel temperatures [K] the answer is read against: the paper's pitch design holds 450-490 K,
#: Cr-Mo tempering and strength loss begin near ~800 K, and it melts near ~1750 K.
STEEL_TARGET, STEEL_STRENGTH, STEEL_MELT = 490.0, 800.0, 1750.0

#: Sprayed-graphite density [kg/m^3], between a porous spray (~1500) and dense graphite (~1800).
GRAPHITE_RHO = 1700.0
#: Through-thickness conductivity bracket [W/m/K]; see the module docstring.
GRAPHITE_K = (3.0, 10.0, 40.0)
#: Grey-body emissivity of the carbon surface.
WALL_EMISSIVITY = 0.9
#: Carbon melt point [K] and heat of fusion [J/kg] (literature; not checked here).
CARBON_MELT = 4800.0
CARBON_FUSION = 105.0e3 / 12.011e-3

#: The paper's pitch model, for the cross-check: k 0.2, rho 1300, c 1600, surface held at 3900 K,
#: and heat beyond it removes pitch at graphite's sublimation heat, 59.3 MJ/kg.
PITCH = (0.2, 1300.0, 1600.0, 3900.0, 59.3e6)

INITIAL_TEMP = 300.0
STEEL_DEPTH = 0.030
PULSE_PERIOD = 0.25  # 4 Hz

OUTPUT = Path("data/results/walled_nozzle/near_term/wall_layers.csv")
OUTPUT_HYDROGEN = Path("data/results/walled_nozzle/near_term/wall_layers_hydrogen.csv")

#: Hot band is 3x the throat area (the paper's "0.594 m^2" at 0.198 m^2); the whole wall is 35.6.
HOT_BAND_FACTOR = 3.0
#: Heat the surface chemistry absorbs [J/kg C] in the sensitivity run: graphite -> C2H2 at 298 K,
#: `dHf(C2H2) = 226.7 kJ/mol` per 2 C (JANAF). A bracket for the dominant 3900 K product, not a
#: solved enthalpy of the equilibrium mixture at the wall.
CHEM_HEAT_C2H2 = 226.7e3 / (2.0 * 12.011e-3)


@dataclass(frozen=True)
class Coating:
    """A surface layer: thickness, properties and the surface temperature cap."""

    name: str
    thickness: float
    k: float
    rho: float
    #: Constant heat capacity [J/kg/K], or `None` for graphite's c_p(T).
    c: float | None
    cap: float
    #: Heat [J/kg] that removes material once the surface sits at `cap`.
    removal_heat: float
    #: Whether the surface is carbon that hot hydrogen can eat (pitch, graphite). Bare steel is not.
    reacts: bool = False


def graphite(thickness: float, k: float) -> Coating:
    """Sprayed graphite of a given thickness and through-thickness conductivity."""
    return Coating(
        f"graphite {thickness * 1e3:.1f} mm k={k:g}",
        thickness,
        k,
        GRAPHITE_RHO,
        None,
        CARBON_MELT,
        CARBON_FUSION,
        reacts=True,
    )


def pitch(thickness: float = 0.2e-3) -> Coating:
    """The paper's pitch model."""
    k, rho, c, cap, heat = PITCH
    return Coating(
        f"pitch {thickness * 1e3:.1f} mm (paper model)", thickness, k, rho, c, cap, heat, True
    )


BARE = Coating("bare steel", 0.0, STEEL_K, STEEL_RHO, STEEL_C, STEEL_MELT, 0.0)


# --- Gas-side boundary --------------------------------------------------------------------------


@dataclass(frozen=True)
class GasHistory:
    """The chamber against time, and the heat-transfer coefficient that goes with it."""

    time: Vec
    temp: Vec
    pressure: Vec
    #: Gas emissivity seen by the wall.
    emissivity: Vec
    #: Bartz coefficient at the peak state [W/m^2/K].
    h0: float
    #: Mass-transfer coefficient `g = h_g / c_p` at the peak state [kg/m^2/s]; zero switches the
    #: chemical ablation off (methane, whose boundary layer plates carbon rather than eating it).
    g0: float = 0.0

    def at(self, t: float) -> tuple[float, float, float]:
        """`(T_gas, h, eps_g)` at time `t` within a pulse; zero load once the chamber is empty."""
        if t > self.time[-1]:
            return INITIAL_TEMP, 0.0, 0.0
        temp = float(np.interp(t, self.time, self.temp))
        p = float(np.interp(t, self.time, self.pressure))
        eps = float(np.interp(t, self.time, self.emissivity))
        return temp, self.h0 * (p / float(self.pressure[0])) ** 0.8, eps

    def ablation_at(self, t: float) -> tuple[float, float]:
        """`(g, p)` at time `t`: the mass-transfer coefficient, scaled as `p^0.8` like `h`, and
        the chamber pressure [Pa] that sets B'. Zero `g` once the chamber is empty."""
        if t > self.time[-1] or self.g0 == 0.0:
            return 0.0, float(self.pressure[-1])
        p = float(np.interp(t, self.time, self.pressure))
        return self.g0 * (p / float(self.pressure[0])) ** 0.8, p


def gas_history(
    throat_area: float,
    edge: str,
    reference_wall: float = 3000.0,
    pairing: near_term.Pairing = near_term.PAIRINGS[1],
    ablation: bool = False,
) -> GasHistory:
    """A near-term pulse (JANAF, equilibrium blowdown) as a wall boundary condition; methane
    7000 K by default.

    `ablation=True` adds `g0 = h0 / c_p` with the `c_p` of the same Bartz edge (frozen for low,
    equilibrium for high), so the `c_p` bracket cancels as the companion's H1 requires. Pair it
    with `reference_wall = 3900` for hydrogen, the temperature the companion took `g` against.

    `edge="low"` takes the low Bartz corner and the H- radiation floor; `"high"` takes the high
    corner and an opaque gas. `h0` is backed out of `bartz_liner` at `reference_wall`, since the
    correlation's own wall-temperature factor is weak and the driving difference is applied live.
    """
    sol = near_term.solve_charge(pairing)
    br = near_term.Branch(sol.feed)
    tab = near_term.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=1e-4, steps=400)
    rows = near_term.blowdown(sol, br, tab, throat_area)
    q_lo, q_hi = near_term.bartz_liner(sol, tab, throat_area, reference_wall)
    h0 = (q_lo if edge == "low" else q_hi) / (sol.pairing.temp - reference_wall)
    cp_fr, cp_eq = near_term.liner_heat_capacities(sol)
    g0 = h0 / (cp_fr if edge == "low" else cp_eq) if ablation else 0.0
    temp = np.array([r.temp for r in rows])
    if edge == "low":
        eps = np.array(
            [near_term.band_fraction(r.temp) * (1.0 - math.exp(-r.tau_floor)) for r in rows]
        )
    else:
        eps = np.ones(len(rows))
    return GasHistory(
        time=np.array([r.time for r in rows]),
        temp=temp,
        pressure=np.array([r.pressure for r in rows]),
        emissivity=eps,
        h0=h0,
        g0=g0,
    )


# --- The solid ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class Grid:
    """Nodes from the surface (0) to the back face, with the interface node index."""

    x: Vec
    interface: int
    coating: Coating

    def is_coating(self) -> NDArray[np.bool_]:
        """Per *segment* (between node i and i+1): whether it lies in the coating."""
        mid = 0.5 * (self.x[1:] + self.x[:-1])
        return np.asarray(mid < self.coating.thickness)


def build_grid(coating: Coating, n_coat: int = 40, first: float = 5e-6) -> Grid:
    """Uniform nodes through the coating, then geometric growth into the steel."""
    xs = [0.0]
    if coating.thickness > 0.0:
        xs += list(np.linspace(0.0, coating.thickness, n_coat + 1)[1:])
        dx = coating.thickness / n_coat
    else:
        dx = first
    interface = len(xs) - 1
    while xs[-1] < coating.thickness + STEEL_DEPTH:
        dx = min(dx * 1.12, 1e-3)
        xs.append(xs[-1] + dx)
    return Grid(np.array(xs), interface, coating)


def _segment_props(grid: Grid, temps: Vec) -> tuple[Vec, Vec]:
    """`(k, rho c)` per segment, graphite c_p evaluated at the segment's mean temperature."""
    in_coat = grid.is_coating()
    c = grid.coating
    tmid = 0.5 * (temps[1:] + temps[:-1])
    k = np.where(in_coat, c.k, STEEL_K)
    if c.c is None:
        cp = np.array([near_term.graphite_heat_capacity(float(t)) for t in tmid])
        rc_coat = c.rho * cp
    else:
        rc_coat = np.full(len(tmid), c.rho * c.c)
    rc = np.where(in_coat, rc_coat, STEEL_RHO * STEEL_C)
    return k, rc


def _blowing(g_mass: float, surface_temp: float, pressure: float) -> tuple[float, float]:
    """`(phi, mdot)`: the blowing correction on the convective heat and the carbon mass flux
    [kg/m^2/s] at a surface temperature and chamber pressure [Pa]."""
    b = blowing.blowing_parameter(surface_temp, pressure)
    return blowing.blowing_correction(b), g_mass * blowing.mass_flux_factor(b)


def _thomas(a: Vec, b: Vec, c: Vec, d: Vec) -> Vec:
    """Tridiagonal solve: `a` sub-, `b` main, `c` super-diagonal."""
    n = len(b)
    cp, dp = np.empty(n), np.empty(n)
    cp[0], dp[0] = c[0] / b[0], d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


@dataclass(frozen=True)
class WallRun:
    """One coating under one pulse edge, over one or more pulses."""

    coating: str
    throat_area: float
    edge: str
    pulses: int
    peak_surface: float
    peak_steel: float
    #: Steel temperature at the interface at the end of the last cycle -- the ratchet.
    end_steel: float
    #: Coating removed per pulse [m] (melt for graphite, ablation for the paper's pitch).
    removed_per_pulse: float
    #: Heat taken into the wall per pulse [J/m^2], and the part still in it at cycle end.
    absorbed_per_pulse: float
    #: Coating eaten by hydrogen chemistry per pulse [m]; `removed_per_pulse` is the thermal part.
    chem_removed_per_pulse: float = 0.0
    #: Time the surface spends within 1 K of its cap in the last pulse [s].
    time_at_cap: float = 0.0
    #: Fraction of the convective flux that survives the blown film, at the peak (1 = no blowing).
    min_blowing_correction: float = 1.0

    @property
    def verdict(self) -> str:
        """Against the paper's own steel targets."""
        if self.peak_steel > STEEL_MELT:
            return "steel melts"
        if self.peak_steel > STEEL_STRENGTH:
            return "steel loses strength"
        if self.peak_steel > STEEL_TARGET:
            return "above paper's steel target"
        return "within paper's steel target"


def run(
    coating: Coating,
    gas: GasHistory,
    pulses: int = 1,
    throat_area: float = 0.0,
    edge: str = "",
    dt_min: float = 1e-6,
    dt_max: float = 2e-4,
    chem_heat: float = 0.0,
) -> WallRun:
    """March the wall through `pulses` cycles at 4 Hz, adiabatic between pulses.

    Adiabatic between pulses is the no-quench case, so the steel's end-of-cycle temperature shows
    the ratchet the paper's between-pulse methane spray has to remove. Each step: build the
    backward-Euler system with properties lagged one Picard pass, linearise the surface flux about
    the current iterate, solve, and if the surface would exceed the cap, pin it there and book the
    excess arriving flux as removed material.

    **Chemical ablation (N18 item 1).** If the gas carries `g0 > 0` and the coating `reacts`, the
    surface also loses carbon continuously, whatever its temperature:

        mdot = g ln(1 + B'(T_s, p))
        q_conv = [ln(1 + B') / B'] h (T_gas - T_s) - mdot * chem_heat

    The blowing factor cuts the convective heat as it cuts the mass flux. B' depends on the very
    surface temperature being solved, so it is lagged one Picard pass like the properties (five
    passes, not three, because B' is steep: x8 to x30 between 2500 and 3900 K). `chem_heat` [J/kg C]
    is the heat the surface chemistry absorbs; default zero, so that is an *omission*, not a
    result -- see the sensitivity in `main`. The grid does not move as material is removed, as
    for the thermal removal already here.
    """
    grid = build_grid(coating)
    n = len(grid.x)
    temps = np.full(n, INITIAL_TEMP)
    dxs = np.diff(grid.x)
    peak_s = peak_st = INITIAL_TEMP
    removed = absorbed = chem_removed = time_at_cap = 0.0
    min_corr = 1.0
    chem = coating.reacts and gas.g0 > 0.0
    for pulse in range(pulses):
        t, dt = 0.0, dt_min
        while t < PULSE_PERIOD - 1e-12:
            dt = min(dt, PULSE_PERIOD - t)
            tg, h, eps = gas.at(t + dt)
            g_mass, p_gas = gas.ablation_at(t + dt) if chem else (0.0, 0.0)
            new = temps.copy()
            pinned = False
            phi, mdot = 1.0, 0.0
            for _ in range(5 if chem else 3):
                k, rc = _segment_props(grid, 0.5 * (temps + new))
                g = k / dxs
                cap = np.zeros(n)
                cap[:-1] += 0.5 * rc * dxs
                cap[1:] += 0.5 * rc * dxs
                ts = new[0]
                rad = WALL_EMISSIVITY * eps * SIGMA_SB
                if chem and g_mass > 0.0:
                    phi, mdot = _blowing(g_mass, min(ts, coating.cap), p_gas)
                q = phi * h * (tg - ts) + rad * (tg**4 - ts**4) - mdot * chem_heat
                dq = -phi * h - 4.0 * rad * ts**3
                a = np.zeros(n)
                b = cap / dt
                c = np.zeros(n)
                d = cap / dt * temps
                a[1:] = -g
                c[:-1] = -g
                b[:-1] += g
                b[1:] += g
                if pinned:
                    b[0], c[0], d[0] = 1.0, 0.0, coating.cap
                else:
                    b[0] -= dq
                    d[0] += q - dq * ts
                new = _thomas(a, b, c, d)
                if not pinned and new[0] > coating.cap:
                    pinned = True
                    new[0] = coating.cap
            if pinned:
                ts = coating.cap
                if chem and g_mass > 0.0:
                    phi, mdot = _blowing(g_mass, ts, p_gas)
                q_in = (
                    phi * h * (tg - ts)
                    + WALL_EMISSIVITY * eps * SIGMA_SB * (tg**4 - ts**4)
                    - mdot * chem_heat
                )
                k, rc = _segment_props(grid, new)
                conducted = k[0] / dxs[0] * (new[0] - new[1])
                stored = 0.5 * rc[0] * dxs[0] * (new[0] - temps[0]) / dt
                excess = max(0.0, q_in - conducted - stored)
                if coating.removal_heat > 0.0 and pulse == pulses - 1:
                    removed += excess * dt / (coating.removal_heat * coating.rho)
            else:
                if chem and g_mass > 0.0:
                    phi, mdot = _blowing(g_mass, min(float(new[0]), coating.cap), p_gas)
                q_in = (
                    phi * h * (tg - new[0])
                    + WALL_EMISSIVITY * eps * SIGMA_SB * (tg**4 - new[0] ** 4)
                    - mdot * chem_heat
                )
            if pulse == pulses - 1:
                absorbed += max(0.0, q_in) * dt
                chem_removed += mdot * dt / coating.rho
                if mdot > 0.0:
                    min_corr = min(min_corr, phi)
                if new[0] >= coating.cap - 1.0:
                    time_at_cap += dt
            temps = new
            peak_s = max(peak_s, float(temps[0]))
            peak_st = max(peak_st, float(temps[grid.interface]))
            t += dt
            quiet = t > gas.time[-1]
            dt = min(dt * 1.05, 5e-3 if quiet else dt_max)
    return WallRun(
        coating=coating.name,
        throat_area=throat_area,
        edge=edge,
        pulses=pulses,
        peak_surface=peak_s,
        peak_steel=peak_st,
        end_steel=float(temps[grid.interface]),
        removed_per_pulse=removed,
        absorbed_per_pulse=absorbed,
        chem_removed_per_pulse=chem_removed,
        time_at_cap=time_at_cap,
        min_blowing_correction=min_corr,
    )


def _scaled(gas: GasHistory, factor: float) -> GasHistory:
    """The same pulse with the convective load (and the mass transfer that rides on it) scaled,
    for a hot spot. Radiation is not scaled."""
    return GasHistory(
        gas.time, gas.temp, gas.pressure, gas.emissivity, gas.h0 * factor, gas.g0 * factor
    )


def hydrogen_main() -> None:
    """N18 item 1: the paper's 0.2 mm pitch over bare Cr-Mo steel under the hydrogen 5500 K pulse.

    Four cases per throat and flux edge: chemistry off (the thermal term alone), chemistry on,
    chemistry on with the C2H2 heat sink, and chemistry on under a 2x hot spot.
    """
    coat = pitch()
    pairing = near_term.PAIRINGS[0]
    cases = (
        ("thermal only", 1.0, False, 0.0),
        ("chem + thermal", 1.0, True, 0.0),
        ("chem + thermal, C2H2 heat sink", 1.0, True, CHEM_HEAT_C2H2),
        ("chem + thermal, 2x hot spot", 2.0, True, 0.0),
    )
    print(f"{coat.name} on Cr-Mo steel, {pairing.name} pulse (single pulse from 300 K)\n")
    rows: list[tuple[str, WallRun, float]] = []
    for throat in near_term.THROAT_AREAS:
        for edge in ("low", "high"):
            gas = gas_history(throat, edge, 3900.0, pairing, ablation=True)
            print(f"A* = {throat} m^2, {edge} flux edge (h0 = {gas.h0:.0f}, g0 = {gas.g0:.3f})")
            for label, factor, chem, heat in cases:
                g = _scaled(gas, factor)
                if not chem:
                    g = GasHistory(g.time, g.temp, g.pressure, g.emissivity, g.h0)
                r = run(coat, g, 1, throat, edge, chem_heat=heat)
                rows.append((label, r, factor))
                total = r.removed_per_pulse + r.chem_removed_per_pulse
                band = total * coat.rho * HOT_BAND_FACTOR * throat
                whole = total * coat.rho * near_term.WALL_AREA
                print(
                    f"  {label:32} chem {r.chem_removed_per_pulse * 1e6:5.1f} + thermal "
                    f"{r.removed_per_pulse * 1e6:5.1f} = {total * 1e6:5.1f} um  band {band:.3f} kg"
                    f"  wall {whole:.2f} kg  surf {r.peak_surface:.0f} K "
                    f"({r.time_at_cap * 1e3:.0f} ms at cap)  steel {r.peak_steel:.0f} K"
                )
    OUTPUT_HYDROGEN.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_HYDROGEN.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "case", "throat_m2", "edge", "chem_um", "thermal_um", "total_um", "hot_band_kg",
                "whole_wall_kg", "peak_surface_k", "time_at_cap_ms", "peak_steel_k",
                "min_blowing_correction",
            ]
        )  # fmt: skip
        for label, r, _ in rows:
            total = r.removed_per_pulse + r.chem_removed_per_pulse
            w.writerow(
                [
                    label, r.throat_area, r.edge, f"{r.chem_removed_per_pulse * 1e6:.2f}",
                    f"{r.removed_per_pulse * 1e6:.2f}", f"{total * 1e6:.2f}",
                    f"{total * coat.rho * HOT_BAND_FACTOR * r.throat_area:.4f}",
                    f"{total * coat.rho * near_term.WALL_AREA:.3f}", f"{r.peak_surface:.0f}",
                    f"{r.time_at_cap * 1e3:.1f}", f"{r.peak_steel:.0f}",
                    f"{r.min_blowing_correction:.3f}",
                ]
            )  # fmt: skip
    print(f"\nwrote {OUTPUT_HYDROGEN}")


def main() -> None:
    """Graphite 1 mm (three conductivities), the paper's pitch, and bare steel, both flux edges."""
    coatings = [*(graphite(1e-3, k) for k in GRAPHITE_K), pitch(), BARE]
    rows: list[WallRun] = []
    print("Methane 7000 K pulse on a coated Cr-Mo wall (single pulse from 300 K)\n")
    for throat in near_term.THROAT_AREAS:
        for edge in ("low", "high"):
            gas = gas_history(throat, edge)
            print(f"A* = {throat} m^2, {edge} flux edge (h0 = {gas.h0:.0f} W/m^2/K)")
            for coat in coatings:
                r = run(coat, gas, 1, throat, edge)
                rows.append(r)
                print(
                    f"  {r.coating:30} surface {r.peak_surface:6.0f} K  steel {r.peak_steel:6.0f} K"
                    f"  removed {r.removed_per_pulse * 1e6:6.1f} um  absorbed "
                    f"{r.absorbed_per_pulse / 1e6:5.2f} MJ/m2  -> {r.verdict}"
                )
    print("\nRatchet: 10 pulses at 4 Hz, no between-pulse cooling, 0.149 m^2 throat, high edge")
    gas = gas_history(near_term.THROAT_AREAS[0], "high")
    for coat in (graphite(1e-3, 3.0), graphite(1e-3, 40.0), pitch()):
        r = run(coat, gas, 10, near_term.THROAT_AREAS[0], "high-x10")
        rows.append(r)
        print(
            f"  {r.coating:30} peak steel {r.peak_steel:6.0f} K, steel at cycle end "
            f"{r.end_steel:6.0f} K -> {r.verdict}"
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "coating", "throat_m2", "edge", "pulses", "peak_surface_k", "peak_steel_k",
                "end_steel_k", "removed_um_per_pulse", "absorbed_mj_m2_per_pulse", "verdict",
            ]
        )  # fmt: skip
        for r in rows:
            w.writerow(
                [
                    r.coating, r.throat_area, r.edge, r.pulses, f"{r.peak_surface:.0f}",
                    f"{r.peak_steel:.0f}", f"{r.end_steel:.0f}", f"{r.removed_per_pulse * 1e6:.2f}",
                    f"{r.absorbed_per_pulse / 1e6:.3f}", r.verdict,
                ]
            )  # fmt: skip
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    if "--hydrogen" in sys.argv[1:]:
        hydrogen_main()
    else:
        main()
