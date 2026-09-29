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
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import near_term

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
    )


def pitch(thickness: float = 0.2e-3) -> Coating:
    """The paper's pitch model."""
    k, rho, c, cap, heat = PITCH
    return Coating(f"pitch {thickness * 1e3:.1f} mm (paper model)", thickness, k, rho, c, cap, heat)


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

    def at(self, t: float) -> tuple[float, float, float]:
        """`(T_gas, h, eps_g)` at time `t` within a pulse; zero load once the chamber is empty."""
        if t > self.time[-1]:
            return INITIAL_TEMP, 0.0, 0.0
        temp = float(np.interp(t, self.time, self.temp))
        p = float(np.interp(t, self.time, self.pressure))
        eps = float(np.interp(t, self.time, self.emissivity))
        return temp, self.h0 * (p / float(self.pressure[0])) ** 0.8, eps


def gas_history(throat_area: float, edge: str, reference_wall: float = 3000.0) -> GasHistory:
    """The methane 7000 K pulse (JANAF, equilibrium blowdown) as a wall boundary condition.

    `edge="low"` takes the low Bartz corner and the H- radiation floor; `"high"` takes the high
    corner and an opaque gas. `h0` is backed out of `bartz_liner` at `reference_wall`, since the
    correlation's own wall-temperature factor is weak and the driving difference is applied live.
    """
    sol = near_term.solve_charge(near_term.PAIRINGS[1])
    br = near_term.Branch(sol.feed)
    tab = near_term.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=1e-4, steps=400)
    rows = near_term.blowdown(sol, br, tab, throat_area)
    q_lo, q_hi = near_term.bartz_liner(sol, tab, throat_area, reference_wall)
    h0 = (q_lo if edge == "low" else q_hi) / (sol.pairing.temp - reference_wall)
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
) -> WallRun:
    """March the wall through `pulses` cycles at 4 Hz, adiabatic between pulses.

    Adiabatic between pulses is the no-quench case, so the steel's end-of-cycle temperature shows
    the ratchet the paper's between-pulse methane spray has to remove. Each step: build the
    backward-Euler system with properties lagged one Picard pass, linearise the surface flux about
    the current iterate, solve, and if the surface would exceed the cap, pin it there and book the
    excess arriving flux as removed material.
    """
    grid = build_grid(coating)
    n = len(grid.x)
    temps = np.full(n, INITIAL_TEMP)
    dxs = np.diff(grid.x)
    peak_s = peak_st = INITIAL_TEMP
    removed = absorbed = 0.0
    for pulse in range(pulses):
        t, dt = 0.0, dt_min
        while t < PULSE_PERIOD - 1e-12:
            dt = min(dt, PULSE_PERIOD - t)
            tg, h, eps = gas.at(t + dt)
            new = temps.copy()
            pinned = False
            for _ in range(3):
                k, rc = _segment_props(grid, 0.5 * (temps + new))
                g = k / dxs
                cap = np.zeros(n)
                cap[:-1] += 0.5 * rc * dxs
                cap[1:] += 0.5 * rc * dxs
                ts = new[0]
                rad = WALL_EMISSIVITY * eps * SIGMA_SB
                q = h * (tg - ts) + rad * (tg**4 - ts**4)
                dq = -h - 4.0 * rad * ts**3
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
                q_in = h * (tg - ts) + WALL_EMISSIVITY * eps * SIGMA_SB * (tg**4 - ts**4)
                k, rc = _segment_props(grid, new)
                conducted = k[0] / dxs[0] * (new[0] - new[1])
                stored = 0.5 * rc[0] * dxs[0] * (new[0] - temps[0]) / dt
                excess = max(0.0, q_in - conducted - stored)
                if coating.removal_heat > 0.0 and pulse == pulses - 1:
                    removed += excess * dt / (coating.removal_heat * coating.rho)
            else:
                q_in = h * (tg - new[0]) + WALL_EMISSIVITY * eps * SIGMA_SB * (tg**4 - new[0] ** 4)
            if pulse == pulses - 1:
                absorbed += max(0.0, q_in) * dt
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
    )


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
    main()
