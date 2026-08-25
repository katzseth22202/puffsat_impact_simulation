"""Cooling history `T(t)` of the shocked plume expanding through the magnetic nozzle.

The companion repo defers the field-leak fraction to "a cooling history, then a quadrature", and
owns the quadrature. This module owns the history. It is a **quasi-1D steady isentropic
expansion** of the post-impact plume, run on the real water EOS rather than a constant-`gamma`
adiabat, because the whole question turns on a chemistry effect a constant `gamma` cannot carry:

**Equilibrium recombination buffers the cooling; frozen recombination does not.** On the way down
the equilibrium gas hands its dissociation and ionisation store back to the thermal pool, holding
`gamma_eff` near 1.15 and keeping the plume hot. If the composition freezes, that store leaves as
inert enthalpy and the temperature falls at the frozen mixture's much stiffer rate. So the plume's
temperature at the end of the field region is bracketed, not single-valued -- and it is the *same*
bracket ADR-0026 already runs on the plate side, reached from the other direction.

**Parametrised by density, so nothing is root-found on the nozzle.** `eos_water` exposes `(p, e)`
from `(rho, T)` and has no entropy function, so the isentrope is obtained by integrating the
adiabat `de = -p d(1/rho)` directly and inverting for `T` at each step. Then `u` follows from
energy conservation, the area ratio from mass conservation, and the time from `dx/u`.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from functools import cache
from itertools import pairwise
from pathlib import Path

import numpy as np

from puffsat import conductivity, eos_water, tops

#: `(rho [kg/m^3], T [K]) -> (p [Pa], e [J/kg])`. Both `eos_water.pressure_energy` and a
#: `pressure_energy_frozen` closure satisfy it, which is what makes the bracket one code path.
Eos = Callable[[float, float], tuple[float, float]]

#: Bracket for the temperature inversion [K]. The floor is set by `eos_water`, whose partition
#: functions underflow below ~50 K, not by the physics -- the plume never approaches it. The
#: ceiling is far above anything reached (item 1's hottest is 26 200 K at 75 km/s).
T_FLOOR = 50.0
T_CEILING = 1.0e6


def temperature_at(rho: float, e_target: float, eos: Eos, tol: float = 1.0e-10) -> float:
    """Invert the caloric EOS: the `T` at which `eos(rho, T)` has specific energy `e_target`.

    Bisection rather than Newton: `e(T)` is monotone in `T` at fixed density for both EOS, but
    through the dissociation and ionisation knees its derivative swings by orders of magnitude,
    which is where a secant method wanders. Monotonicity is all bisection needs.
    """
    lo, hi = T_FLOOR, T_CEILING
    if eos(rho, hi)[1] < e_target:
        raise ValueError(f"e = {e_target:.4e} J/kg is above the EOS ceiling at rho = {rho:.4e}")
    if eos(rho, lo)[1] > e_target:
        raise ValueError(f"e = {e_target:.4e} J/kg is below the EOS floor at rho = {rho:.4e}")
    while hi - lo > tol * hi:
        mid = 0.5 * (lo + hi)
        if eos(rho, mid)[1] < e_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class ExpansionState:
    """One point on the isentrope."""

    rho: float
    temp: float
    pressure: float
    energy: float

    @property
    def enthalpy(self) -> float:
        """Specific enthalpy `h = e + p/rho` [J/kg] -- the quantity the nozzle conserves."""
        return self.energy + self.pressure / self.rho


def expand(
    rho_0: float, temp_0: float, rho_end: float, eos: Eos, steps: int = 256
) -> list[ExpansionState]:
    """Integrate the adiabat `de = -p d(1/rho)` from `(rho_0, temp_0)` down to `rho_end`.

    This is the isentrope: adiabatic *and* reversible, so no entropy function is needed. The
    grid is uniform in specific volume `v = 1/rho` because that is the integration variable.

    Heun's predictor-corrector (second order): the predictor takes an Euler step on `de/dv = -p`,
    the corrector re-evaluates `p` at the predicted state and averages. Each evaluation costs a
    temperature inversion, which for the real EOS is the dominant cost -- hence second order
    rather than a higher-order scheme with more stages.
    """
    if not 0.0 < rho_end < rho_0:
        raise ValueError(f"expansion needs 0 < rho_end ({rho_end}) < rho_0 ({rho_0})")
    v_0, v_end = 1.0 / rho_0, 1.0 / rho_end
    ratio = (v_end / v_0) ** (1.0 / steps)

    p, e = eos(rho_0, temp_0)
    states = [ExpansionState(rho=rho_0, temp=temp_0, pressure=p, energy=e)]
    v = v_0
    for _ in range(steps):
        v_next = v * ratio
        dv = v_next - v
        rho_next = 1.0 / v_next
        e_predict = e - p * dv
        p_predict = eos(rho_next, temperature_at(rho_next, e_predict, eos))[0]
        e = e - 0.5 * (p + p_predict) * dv
        temp = temperature_at(rho_next, e, eos)
        p, _ = eos(rho_next, temp)
        states.append(ExpansionState(rho=rho_next, temp=temp, pressure=p, energy=e))
        v = v_next
    return states


#: `(rho, T) -> c_s [m/s]`. `eos_water.sound_speed` and `sound_speed_frozen` both satisfy it.
SoundSpeed = Callable[[float, float], float]


@dataclass(frozen=True)
class NozzlePoint:
    """One station in the nozzle: the isentrope point plus the flow quantities it implies."""

    rho: float
    temp: float
    pressure: float
    energy: float
    #: Flow speed [m/s] from `h + u^2/2 = h0`, measured from the plume at rest in the bag.
    speed: float
    sound_speed: float
    #: `A/A*`, the area referenced to the sonic throat -- the ratio the paper's field profile
    #: gives directly, since flux conservation makes `A/A* = B*/B`.
    area_ratio: float

    @property
    def enthalpy(self) -> float:
        """Specific enthalpy `h = e + p/rho` [J/kg]."""
        return self.energy + self.pressure / self.rho

    @property
    def mach(self) -> float:
        """`M = u/c_s`. Below 1 upstream of the throat, above it in the expansion."""
        return self.speed / self.sound_speed


def nozzle_history(
    rho_0: float,
    temp_0: float,
    rho_end: float,
    eos: Eos,
    sound_speed: SoundSpeed,
    steps: int = 256,
) -> list[NozzlePoint]:
    """Quasi-1D steady isentropic expansion from a plume at rest at `(rho_0, temp_0)`.

    Parametrised by density, so no equation is root-found on the nozzle:

    1. the adiabat gives `e`, `T` and `p` at each density (`expand`);
    2. `u = sqrt(2 (h0 - h))` -- total enthalpy is conserved, the plume starting at rest;
    3. `rho u A = const` fixes the area, referenced to the sonic point `M = 1` so the result is
       `A/A*`, which is what the paper's `B*/B` field profile reports.

    The area reference is the throat rather than the reservoir because the reservoir's area is
    infinite (`u -> 0`), and because the paper states the nozzle as a field ratio -- `20 T -> 5 T`
    is `A/A* = 4` by flux conservation, with no further assumption.
    """
    states = expand(rho_0, temp_0, rho_end, eos, steps)
    h_0 = states[0].enthalpy

    raw: list[tuple[ExpansionState, float, float, float]] = []
    for st in states:
        # Clamp at zero: the reservoir point is h == h0 exactly and can land a hair negative.
        speed = math.sqrt(max(0.0, 2.0 * (h_0 - st.enthalpy)))
        c_s = sound_speed(st.rho, st.temp)
        raw.append((st, speed, c_s, st.rho * speed))

    mass_flux_star = _sonic_mass_flux(raw)
    return [
        NozzlePoint(
            rho=st.rho,
            temp=st.temp,
            pressure=st.pressure,
            energy=st.energy,
            speed=speed,
            sound_speed=c_s,
            area_ratio=math.inf if flux == 0.0 else mass_flux_star / flux,
        )
        for st, speed, c_s, flux in raw
    ]


def _sonic_mass_flux(raw: list[tuple[ExpansionState, float, float, float]]) -> float:
    """Mass flux `rho u` at the throat, from `rho` and `u` interpolated on the `M = 1` crossing.

    `rho` and `u` are interpolated **separately** and then multiplied, rather than interpolating
    the product: `rho u` is at a maximum here, so a chord across the product is biased low and
    would put the neighbouring samples at `A/A* < 1` -- i.e. *behind* the throat, and at negative
    distance along the bore. Interpolating the state and forming the flux from it is also exactly
    what `_at_throat` does, so the throat station lands at `A/A* = 1` by construction instead of
    by luck.
    """
    for (state_a, u_a, c_a, _), (state_b, u_b, c_b, _) in pairwise(raw):
        if u_a - c_a <= 0.0 <= u_b - c_b:
            span = (u_b - c_b) - (u_a - c_a)
            w = 0.0 if span == 0.0 else -(u_a - c_a) / span
            rho = state_a.rho + w * (state_b.rho - state_a.rho)
            return rho * (u_a + w * (u_b - u_a))
    raise ValueError("the expansion never reaches M = 1; extend rho_end further down")


@dataclass(frozen=True)
class CoolingRow:
    """One row of the cooling history `T(t)` -- the deliverable the companion repo consumes."""

    #: Time since the plume passed the sonic throat [s].
    time: float
    #: Distance along the field region [m].
    x: float
    area_ratio: float
    rho: float
    temp: float
    pressure: float
    #: Specific internal energy [J/kg] -- the reservoir the radiative cooling time draws down.
    energy: float
    speed: float
    mach: float


def cooling_history(
    rho_0: float,
    temp_0: float,
    eos: Eos,
    sound_speed: SoundSpeed,
    area_ratio_end: float,
    length: float,
    steps: int = 256,
    expansion_ratio: float = 512.0,
) -> list[CoolingRow]:
    """`T(t)` from the sonic throat to `area_ratio_end`, over a field region `length` long.

    **`T` against area ratio is the robust half of this result; the clock is the assumed half.**
    The temperature at a station depends only on how far the flux tube has opened, which the
    paper states directly as a field ratio (`A/A* = B*/B`, so `20 T -> 5 T` is `A/A* = 4`). The
    *time* additionally needs the area profile along the bore, which the paper does not give;
    a linear opening over `length` is assumed here, and it is the only place a shape enters.
    Consumers that can work in area ratio should, and are then free of that assumption.

    Positions are interpolated onto the requested exit ratio so the endpoint is exact rather than
    whichever density sample happened to land nearest.
    """
    if area_ratio_end <= 1.0:
        raise ValueError(f"area_ratio_end must exceed 1 (the throat), got {area_ratio_end}")
    points = nozzle_history(rho_0, temp_0, rho_0 / expansion_ratio, eos, sound_speed, steps)

    # Keep the supersonic branch: from the throat (M = 1) out to the requested area ratio.
    throat = _at_throat(points)
    supersonic = [throat, *(pt for pt in points if pt.mach > throat.mach)]
    if supersonic[-1].area_ratio < area_ratio_end:
        raise ValueError(
            f"expansion_ratio = {expansion_ratio} does not reach A/A* = {area_ratio_end}"
        )
    kept = [pt for pt in supersonic if pt.area_ratio < area_ratio_end]
    kept.append(_at_area_ratio(supersonic, area_ratio_end))

    rows: list[CoolingRow] = []
    time = 0.0
    for i, pt in enumerate(kept):
        x = length * (pt.area_ratio - 1.0) / (area_ratio_end - 1.0)
        if i > 0:
            # Trapezoid on 1/u, which is smooth and monotone across a supersonic expansion.
            time += (x - rows[-1].x) * 0.5 * (1.0 / pt.speed + 1.0 / kept[i - 1].speed)
        rows.append(
            CoolingRow(
                time=time,
                x=x,
                area_ratio=pt.area_ratio,
                rho=pt.rho,
                temp=pt.temp,
                pressure=pt.pressure,
                energy=pt.energy,
                speed=pt.speed,
                mach=pt.mach,
            )
        )
    return rows


def _blend(lo: NozzlePoint, hi: NozzlePoint, w: float) -> NozzlePoint:
    """Linear blend of two stations, weight `w` from `lo` to `hi`."""

    def mix(a: float, b: float) -> float:
        return a + w * (b - a)

    return NozzlePoint(
        rho=mix(lo.rho, hi.rho),
        temp=mix(lo.temp, hi.temp),
        pressure=mix(lo.pressure, hi.pressure),
        energy=mix(lo.energy, hi.energy),
        speed=mix(lo.speed, hi.speed),
        sound_speed=mix(lo.sound_speed, hi.sound_speed),
        area_ratio=mix(lo.area_ratio, hi.area_ratio),
    )


def _at_throat(points: list[NozzlePoint]) -> NozzlePoint:
    """The `M = 1` station, interpolated on the Mach crossing.

    It cannot be interpolated on area ratio: `A/A*` is 1 *at* the throat and rises on **both**
    sides, so it does not bracket the throat. `M` is monotone through it and does.
    """
    for lo, hi in pairwise(points):
        if lo.mach <= 1.0 <= hi.mach:
            span = hi.mach - lo.mach
            throat = _blend(lo, hi, 0.0 if span == 0.0 else (1.0 - lo.mach) / span)
            # `A/A* = 1` at `M = 1` by the definition of `A*`. Interpolation lands a hair above
            # it because `A/A*` has a *minimum* there, so the exact value is restored rather than
            # carried as a spurious few-micron offset in `x`.
            return replace(throat, area_ratio=1.0)
    raise ValueError("the expansion never reaches M = 1; extend the expansion ratio")


def _at_area_ratio(points: list[NozzlePoint], target: float) -> NozzlePoint:
    """Interpolate the supersonic branch onto an exact area ratio, where `A/A*` is monotone."""
    for lo, hi in pairwise(points):
        if lo.area_ratio <= target <= hi.area_ratio:
            span = hi.area_ratio - lo.area_ratio
            return _blend(lo, hi, 0.0 if span == 0.0 else (target - lo.area_ratio) / span)
    raise ValueError(f"A/A* = {target} is outside the supersonic branch")


#: Stefan-Boltzmann constant [W m^-2 K^-4].
SIGMA_SB = 5.670374419e-8


@dataclass(frozen=True)
class RadiationCheck:
    """Whether the isentrope is allowed to ignore radiation at one station.

    The expansion above is adiabatic. That is a *claim*, and this is the test of it: if the
    radiative cooling time is short against the transit, the plume does not follow the isentrope
    and the history is wrong in the direction of extra cooling.
    """

    #: Rosseland (transport) optical depth across the flux tube.
    optical_depth: float
    #: Planck (emission) optical depth. TOPS puts this ~100x the Rosseland value for water here,
    #: which is why the regime cannot be read off the Rosseland depth alone.
    planck_depth: float
    #: Which of the three limits binds: "emission", "free-streaming" or "diffusion".
    regime: str
    cooling_time_thin: float
    cooling_time_thick: float
    cooling_time_free_streaming: float
    #: The binding branch -- the *longest* of the three, since the smallest loss rate wins.
    cooling_time: float


def radiation_check(
    temp: float,
    rho: float,
    energy: float,
    kappa_planck: float,
    kappa_rosseland: float,
    radius: float,
) -> RadiationCheck:
    """Radiative cooling time of the plume at one station: flux-limited, two-mean, gray.

    Three loss rates compete and the **smallest binds**, which is ADR-0006's gray two-mean
    flux-limited form applied to a whole flux tube instead of a cell:

    - **Emission-limited** `4 kappa_P rho sigma T^4` per unit volume -- every photon the gas
      makes escapes. Valid only where the *Planck* depth is below 1.
    - **Free-streaming** `(2/R) sigma T^4` -- the flux limiter. Nothing radiates faster than a
      blackbody at its own temperature through its own surface, whatever the opacity says.
    - **Diffusion** `(2/R) 4 sigma T^4 / (3 tau_R)` -- photons random-walk out, so *more*
      opacity now means *slower* cooling. The opposite sign from the emission limit, which is
      why the two cannot be interpolated with one formula.

    Taking the minimum rather than branching on `tau_Rosseland` is not a detail: near the
    crossover the Rosseland depth reads "thin" while the Planck depth is ~100x larger and the gas
    is deeply thick to its own lines. Branching on the wrong mean there overstates the loss by
    two to three orders of magnitude and reports the plume radiating away several times its own
    internal energy.
    """
    tau_r = kappa_rosseland * rho * radius
    tau_p = kappa_planck * rho * radius
    flux_bb = SIGMA_SB * temp**4
    reservoir = rho * energy  # [J/m^3]

    loss_emission = 4.0 * kappa_planck * rho * flux_bb
    loss_free = (2.0 / radius) * flux_bb
    loss_diffusion = (2.0 / radius) * 4.0 * flux_bb / (3.0 * tau_r) if tau_r > 0.0 else math.inf

    times = {
        "emission": reservoir / loss_emission,
        "free-streaming": reservoir / loss_free,
        "diffusion": reservoir / loss_diffusion if loss_diffusion > 0.0 else math.inf,
    }
    regime = max(times, key=lambda k: times[k])
    return RadiationCheck(
        optical_depth=tau_r,
        planck_depth=tau_p,
        regime=regime,
        cooling_time_thin=times["emission"],
        cooling_time_thick=times["diffusion"],
        cooling_time_free_streaming=times["free-streaming"],
        cooling_time=times[regime],
    )


#: Post-impact plume states from the companion repo's item 1 (Saha + energy conservation at the
#: flown bag density), cross-checked against `eos_water` to 1-3% in `T`. Keyed by closing speed.
PLUME_STATES: tuple[tuple[float, float], ...] = (
    (75.0, 26200.0),
    (65.0, 22400.0),
    (56.53, 19400.0),
    (45.58, 14700.0),
)
#: Flown bag: 213 kg over 660 m^3.
BAG_RHO = 0.323
#: Throat radius [m] -- the paper's 3.0 m bore, where the field is 20 T.
THROAT_RADIUS = 3.0
#: `A/A* = B*/B`: the paper's `20 T -> 5 T` nozzle. Flux conservation, no further assumption.
AREA_RATIO_EXIT = 4.0
#: Field-region length [m]: the 23 m capsule (`8 pi r^3 = 660 m^3` at aspect 4 gives 23.8 m).
FIELD_LENGTH = 23.8
#: Seed mass fraction the paper flies.
SEED_FRACTION = 0.01

DEFAULT_COOLING_PATH = Path("data/results/cooling_history.csv")

CSV_HEADER = (
    "closing_speed_km_s,branch,time_ms,x_m,area_ratio,rho,temp_k,pressure_pa,speed_m_s,mach,"
    "radius_m,optical_depth,radiation_regime,t_rad_over_t_transit,v_l,sigma_s_m,rm,"
    "leak_fraction,gamma_t"
)


def plume_radius(area_ratio: float) -> float:
    """Flux-tube radius [m]: `A/A*` is an area ratio, so the radius goes as its square root."""
    return THROAT_RADIUS * math.sqrt(area_ratio)


@dataclass(frozen=True)
class IsentropicExponents:
    """The two exponents of an isentrope, measured between adjacent stations rather than assumed.

    `gamma_pressure` (`Gamma_1`) governs the sound speed; **`gamma_temperature` (`Gamma_3`) is the
    one the field-retention question turns on** (Q-J), because what decides whether `Rm` rises or
    collapses through the push is how fast `T` falls as the plume thins, not how fast `p` does.
    Quoting the pressure exponent for that argument overstates it by ~0.06.
    """

    #: `dln p / dln rho`.
    gamma_pressure: float
    #: `1 + dln T / dln rho`, so an ideal constant-`gamma` gas returns its own `gamma`.
    gamma_temperature: float


def isentropic_exponents(lo: CoolingRow, hi: CoolingRow) -> IsentropicExponents:
    """Local isentropic exponents across one station pair.

    A **diagnostic of the real EOS**, not an input to it: nothing in this module assumes a
    `gamma`. Reading it back out is how the chemistry's buffering becomes a number -- 5/3 means
    no store is being returned, and the further below it the exponent sits the harder the
    chemistry is working to hold the temperature up.
    """
    if lo.rho == hi.rho:
        raise ValueError("stations must differ in density")
    span = math.log(hi.rho / lo.rho)
    return IsentropicExponents(
        gamma_pressure=math.log(hi.pressure / lo.pressure) / span,
        gamma_temperature=1.0 + math.log(hi.temp / lo.temp) / span,
    )


@dataclass(frozen=True)
class HistoryRow:
    """A cooling-history station with everything the companion repo's quadrature consumes."""

    closing_speed: float
    branch: str
    row: CoolingRow
    radiation: RadiationCheck
    #: Transit time from the throat to this station [s] -- the clock radiation is judged against.
    transit_time: float
    #: `v L` from the solution itself: local flow speed times local flux-tube radius. This is the
    #: product `tab:seed_window` needs and does not state, and it is now an output.
    v_l: float
    sigma: float
    rm: float
    leak_fraction: float


def history(
    closing_speed: float, temp_0: float, frozen: bool, steps: int = 320, stride: int = 4
) -> list[HistoryRow]:
    """One cooling history, on the equilibrium or frozen branch of the ADR-0026 bracket.

    The default resolution is set by the **radiated-energy integral**, not by the temperature
    history: `T(t)` is converged at a quarter of it, but the radiative loss is concentrated in a
    narrow window around the opacity crossover, and at 8 stations that window is under-sampled by
    a factor ~3. At these settings the radiated fraction is converged to three digits.
    """
    if frozen:
        y = eos_water.frozen_composition(BAG_RHO, temp_0)
        eos: Eos = lambda r, t: eos_water.pressure_energy_frozen(r, t, y)  # noqa: E731
        c_s: SoundSpeed = lambda r, t: eos_water.sound_speed_frozen(r, t, y)  # noqa: E731
    else:
        eos, c_s = eos_water.pressure_energy, eos_water.sound_speed

    rows = cooling_history(BAG_RHO, temp_0, eos, c_s, AREA_RATIO_EXIT, FIELD_LENGTH, steps=steps)
    sampled = [*rows[::stride], rows[-1]] if rows[-1] not in rows[::stride] else rows[::stride]

    out: list[HistoryRow] = []
    for row in sampled:
        radius = plume_radius(row.area_ratio)
        kap_p, kap_r = _tops_kappa(row.rho, row.temp)
        rad = radiation_check(row.temp, row.rho, row.energy, kap_p, kap_r, radius)
        v_l = row.speed * radius
        sig = conductivity.sigma(row.temp, row.rho, SEED_FRACTION)
        rm = conductivity.magnetic_reynolds(sig, v_l)
        out.append(
            HistoryRow(
                closing_speed=closing_speed,
                branch="frozen" if frozen else "equilibrium",
                row=row,
                radiation=rad,
                transit_time=row.time,
                v_l=v_l,
                sigma=sig,
                rm=rm,
                leak_fraction=min(1.0, 1.0 / rm) if rm > 0.0 else 1.0,
            )
        )
    return out


DEFAULT_TOPS_PATH = Path("data/tables/tops/tops_water_gray.html")


@cache
def _tops_grid(path: Path = DEFAULT_TOPS_PATH) -> tops.TopsGray:
    """The TOPS gray pull, parsed once. Cached because every station queries it."""
    return tops.load_tops_gray(path)


def _tops_kappa(rho: float, temp: float) -> tuple[float, float]:
    """`(kappa_Planck, kappa_Rosseland)` [m^2/kg] at `(rho, T)`, log-log bilinear on the TOPS grid.

    Clamped at the grid edges. The TOPS floor is 5802 K, below which this returns the floor
    value -- the cool end of a frozen expansion runs off the bottom of the table, and the check
    there is a lower bound on the opacity rather than the opacity.
    """
    grid = _tops_grid()
    out: list[float] = []
    for field in (grid.kappa_planck, grid.kappa_rosseland):
        i = float(
            np.clip(
                np.interp(math.log(rho), np.log(grid.rho_grid), np.arange(len(grid.rho_grid))),
                0,
                len(grid.rho_grid) - 1,
            )
        )
        j = float(
            np.clip(
                np.interp(math.log(temp), np.log(grid.t_grid), np.arange(len(grid.t_grid))),
                0,
                len(grid.t_grid) - 1,
            )
        )
        i0, j0 = min(int(i), len(grid.rho_grid) - 2), min(int(j), len(grid.t_grid) - 2)
        fi, fj = i - i0, j - j0
        ln = np.log(field[i0 : i0 + 2, j0 : j0 + 2])
        val = (
            ln[0, 0] * (1 - fi) * (1 - fj)
            + ln[1, 0] * fi * (1 - fj)
            + ln[0, 1] * (1 - fi) * fj
            + ln[1, 1] * fi * fj
        )
        out.append(float(np.exp(val)))
    return out[0], out[1]


def residence_weighted_leak(samples: Sequence[tuple[float, float]]) -> float:
    """Time-average of the leak fraction over the field residence -- item 10's quadrature.

    `samples` are `(time [s], leak fraction)` in time order. The companion repo's framing is
    exact: "weight `1/Rm(T)` by how long the plume spends at each `T` while the field is doing
    work, and integrate". A station average would not do, because the plume crosses the cool,
    leaky end of the history quickly and the hot end slowly -- so the two averages differ by
    much more than their uncertainty.

    Trapezoid in time. Returns the leak that a uniform history of the same duration would give.
    """
    if len(samples) < 2:
        raise ValueError("the quadrature needs at least two stations")
    total = samples[-1][0] - samples[0][0]
    if total <= 0.0:
        raise ValueError("the history has no duration")
    integral = sum(
        (t_b - t_a) * 0.5 * (leak_a + leak_b) for (t_a, leak_a), (t_b, leak_b) in pairwise(samples)
    )
    return integral / total


def radiated_fraction(samples: Sequence[tuple[float, float]]) -> float:
    """Fraction of internal energy lost to radiation over the transit -- `integral dt / t_rad`.

    `samples` are `(time [s], radiative cooling time [s])` in time order. This is the gate on the
    whole model: the expansion above is adiabatic, and that is only allowed while this number is
    small. It is *not* well estimated by a station average, because the radiative loss is
    concentrated in a narrow window around the opacity crossover -- under-sampling that window
    understates the loss by a factor of ~3.
    """
    if len(samples) < 2:
        raise ValueError("the quadrature needs at least two stations")
    return sum(
        (t_b - t_a) * 0.5 * (1.0 / rad_a + 1.0 / rad_b)
        for (t_a, rad_a), (t_b, rad_b) in pairwise(samples)
    )


def write_cooling_history(rows: Sequence[HistoryRow], path: Path = DEFAULT_COOLING_PATH) -> None:
    """Write the cooling history the companion repo's quadrature consumes.

    `gamma_t` is carried per station because Q-J's field-retention question reduces to it: an
    exponent near 1.15 means the chemistry is buffering and `Rm` rises through the push, one near
    5/3 means it is spent and `Rm` collapses. It is blank on the first station of each branch,
    which has no predecessor to difference against.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [CSV_HEADER]
    previous: HistoryRow | None = None
    for h in rows:
        r, rad = h.row, h.radiation
        same_curve = (
            previous is not None
            and previous.branch == h.branch
            and previous.closing_speed == h.closing_speed
            and previous.row.rho != r.rho
        )
        gamma_t = (
            f"{isentropic_exponents(previous.row, r).gamma_temperature:.5f}"
            if same_curve and previous is not None
            else ""
        )
        previous = h
        lines.append(
            f"{h.closing_speed:g},{h.branch},{r.time * 1e3:.6f},{r.x:.4f},{r.area_ratio:.5f},"
            f"{r.rho:.6e},{r.temp:.2f},{r.pressure:.6e},{r.speed:.2f},{r.mach:.4f},"
            f"{plume_radius(r.area_ratio):.4f},{rad.optical_depth:.6e},{rad.regime},"
            f"{rad.cooling_time / max(h.transit_time, 1e-12):.4e},{h.v_l:.6e},"
            f"{h.sigma:.6e},{h.rm:.6e},{h.leak_fraction:.6f},{gamma_t}"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    """Run both branches of the bracket at every plume state and report what the threads need."""
    all_rows: list[HistoryRow] = []
    print(
        f"{'w':>7} {'branch':>11} {'T_exit':>8} {'t':>7} {'min t_rad/t':>12} "
        f"{'v L':>9} {'leak_exit':>10} {'leak_avg':>9}"
    )
    for speed, temp_0 in PLUME_STATES:
        for frozen in (False, True):
            rows = history(speed, temp_0, frozen)
            all_rows.extend(rows)
            total = rows[-1].row.time
            ratio = min(h.radiation.cooling_time / total for h in rows)
            radiated = radiated_fraction([(h.row.time, h.radiation.cooling_time) for h in rows])
            avg = residence_weighted_leak([(h.row.time, h.leak_fraction) for h in rows])
            print(
                f"{speed:7.2f} {rows[-1].branch:>11} {rows[-1].row.temp:8.0f} "
                f"{total * 1e3:6.2f}ms {ratio:12.1f} {radiated:9.4f} {rows[-1].v_l:9.3g} "
                f"{rows[-1].leak_fraction:10.4f} {avg:9.4f}"
            )
    write_cooling_history(all_rows)
    print(f"\nwrote {len(all_rows)} stations -> {DEFAULT_COOLING_PATH}")


if __name__ == "__main__":
    main()
