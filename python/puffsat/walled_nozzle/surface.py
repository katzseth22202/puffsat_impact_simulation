"""N16: one conversion-fraction surface over chamber `T`, chamber volume and throat area.

The ask asks for three dials at once and says why: "everything that matters reduces to the exit
temperature and the exit density. The throat sets the first. The chamber volume sets the second.
The chamber temperature sets both." Sweeping them separately answers none of them.

**The trick that makes a 128-point grid cheap is that the throat is not a third dimension.**
A chamber at `(rho_c, T_c)` has exactly one isentrope, and the throat area only decides *where
along it* the gas stops -- the exit area ratio is `A_bore / A*` with `A_bore` fixed at the 3 m
bore. So one deep expansion per `(fluid, T_c, V)` is run out to the widest area ratio in the
grid, and every throat is read off it. That is 32 expansions rather than 256, and it also
guarantees the throat comparison is exactly like-for-like: the same gas, the same adiabat, one
cut point moved.

**The clock is re-derived on a cone rather than inherited.** `expansion.cooling_history` opens
the area *linearly in `A`* over a length the caller supplies, which was right when the geometry
was a field profile of fixed length (`freeze.py`). It is wrong for this ask: at fixed length a
larger area ratio means a faster opening, so `Da` would fall with area ratio for a reason that is
pure bookkeeping. Here the nozzle is a physical cone from the throat radius out to the 3 m bore
at `CONE_HALF_ANGLE_DEG`, so

    x(A) = (sqrt(A / pi) - r*) / tan(theta),      L = (R_bore - r*) / tan(theta)

and the length is set by where the cone starts. Across the throats the ask now recommends
(0.5 m^2 and below) `r*` is small against the 3 m bore and the length is nearly constant --
9.7 to 10.7 m, within 11% -- so those columns differ only in where they start. The flown 7 m^2
throat is genuinely a *shorter* nozzle (5.6 m), because a 1.5 m throat radius is already half the
bore. That is a real property of a cone and not a bookkeeping artefact, and it is why the clock is
rebuilt on the geometry rather than held at one length: a wide throat has less nozzle to run the
chemistry in, and the grid should say so rather than hide it.

**What this module reports and what it does not.** It reports the equilibrium branch: exit `T`,
exit density, conversion fraction, store returned, and the freeze race that says whether the
equilibrium branch is legitimate. It does *not* re-solve a frozen expansion where the race is
lost -- W5 found the flow does not freeze inside this geometry, and the grid's job is to check
that at the new densities, not to price a branch that does not occur. Where a point does freeze,
the freeze station and the store still held there are reported and the conversion is flagged.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from puffsat import eos_methane, eos_water, expansion, recombination
from puffsat.walled_nozzle import chamber, freeze, hydrogen, propellants, rates

#: Bore area [m^2] -- the exit plane of every nozzle in this grid. ADR-0016's 3 m bore radius.
BORE_AREA = math.pi * chamber.BORE_RADIUS**2

#: Half-angle [deg] of the conical nozzle the clock is integrated along. 15 degrees is the
#: ordinary bell-nozzle value; `Da` scales linearly in the resulting length, so the sensitivity is
#: transparent and `main` prints it at 10 and 20 degrees as well.
CONE_HALF_ANGLE_DEG = 15.0

#: The grid the ask names.
TEMPERATURES = (6000.0, 8000.0, 10000.0, 12000.0)
VOLUMES = (50.0, 100.0, 200.0, 400.0)
#: Throat areas [m^2], the flown 7 down to 0.05, log-spaced. Area ratios 4.04 to 565.
THROATS = (7.0, 4.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05)

#: Pulse period [s] the blowdown has to fit inside (ADR-0016's 400 ms).
PULSE_PERIOD = 0.400

#: Below this the equation of state has no condensed carbon and `eos_methane` is quoting a
#: gas-phase equilibrium that reality would partly replace with soot (ADR-0050 weakness 4).
CARBON_FLOOR_K = 4000.0

DEFAULT_OUTPUT = Path("data/results/walled_nozzle/surface.csv")
HELD_OUTPUT = Path("data/results/walled_nozzle/held_surface.csv")
TEMPERATURE_OUTPUT = Path("data/results/walled_nozzle/exit_temperature.csv")

Eos = Callable[[float, float], tuple[float, float]]
SoundSpeed = Callable[[float, float], float]


# ---- Fluids ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Fluid:
    """One working fluid, with everything the grid needs to run it end to end."""

    name: str
    pressure_energy: Eos
    sound_speed: SoundSpeed
    #: `(rho, T) -> ` bond energy still held [J/kg].
    held: Callable[[float, float], float]
    #: Bond energy of a fully atomised kilogram [J/kg]: the denominator of `store_held`.
    full_atomisation: float
    #: `(rho, T, tau_exp) -> (Da optimistic, Da conservative)`. Equal where the channel has no
    #: scarcity bracket (methane's `H + H + M`, where both partners are the same species).
    damkohler: Callable[[float, float, float], tuple[float, float]]
    #: `True` where the fluid is assembled rather than solved on its own equilibrium EOS.
    estimated: bool = False


def _methane_damkohler(rho: float, temp: float, tau_exp: float) -> tuple[float, float]:
    """`H + H + M -> H2 + M` on the solved methane composition (the channel W5 established)."""
    comp = eos_methane.composition(rho, temp)
    k = freeze.mixture_coefficient(temp, comp)
    tau = freeze.recombination_time(k, comp.n_h, comp.n_neutral_heavy)
    da = tau_exp / tau if tau > 0.0 else math.inf
    return da, da


def _water_damkohler(rho: float, temp: float, tau_exp: float) -> tuple[float, float]:
    """`H + OH + M -> H2O + M`, carried as the scarcity bracket `recombination.py` established.

    The optimistic edge assumes an OH is always available so the H-atom density paces the return;
    the conservative edge holds OH at its equilibrium value, which is far scarcer. The truth is
    between and locating it needs an evolved network. **The conservative edge is the one quoted as
    the verdict here**, because a walled nozzle's whole argument is that it holds density up long
    enough, and that claim should be made on the slow edge.
    """
    comp = eos_water.composition(rho, temp)
    n_third = comp.n_neutral_heavy
    tau_h = recombination.atom_recombination_time(temp, n_third, comp.n_h)
    tau_oh = recombination.atom_recombination_time(temp, n_third, comp.n_oh)
    return (
        tau_exp / tau_h if tau_h > 0 else math.inf,
        tau_exp / tau_oh if tau_oh > 0 else math.inf,
    )


def _hydrogen_held(rho: float, temp: float) -> float:
    """Dissociation energy still held by pure hydrogen [J/kg]."""
    return hydrogen.state(rho, temp).chemical_energy


def _hydrogen_damkohler(rho: float, temp: float, tau_exp: float) -> tuple[float, float]:
    """`H + H + M -> H2 + M` on pure hydrogen, on the same Jacobs coefficients methane uses."""
    s = hydrogen.state(rho, temp)
    n_third = s.n_h + s.n_h2
    if n_third <= 0.0:
        return math.inf, math.inf
    k = (rates.H_ATOMIC_H(temp) * s.n_h + rates.H_H2_JACOBS(temp) * s.n_h2) / n_third
    tau = freeze.recombination_time(k, s.n_h, n_third)
    da = tau_exp / tau if tau > 0.0 else math.inf
    return da, da


METHANE = Fluid(
    "methane",
    eos_methane.pressure_energy,
    eos_methane.sound_speed,
    eos_methane.bond_energy_held,
    eos_methane.FULL_ATOMIZATION_ENERGY,
    _methane_damkohler,
)
WATER = Fluid(
    "water",
    eos_water.pressure_energy,
    eos_water.sound_speed,
    eos_water.bond_energy_held,
    eos_water.FULL_ATOMIZATION_ENERGY,
    _water_damkohler,
)
HYDROGEN = Fluid(
    "hydrogen",
    hydrogen.pressure_energy,
    hydrogen.sound_speed,
    _hydrogen_held,
    hydrogen.FULL_DISSOCIATION_ENERGY,
    _hydrogen_damkohler,
)
GRID_FLUIDS = (METHANE, WATER)


# ---- The mixture the ask names, assembled rather than solved ------------------------------------

#: Mass fraction of liquid hydrogen in the configuration N16 asks to price first.
HYDROGEN_FRACTION = 0.10


def mixture(water_share: float = 1.0 - HYDROGEN_FRACTION) -> Fluid:
    """Water with a hydrogen dilution, as an **ideal mixture of two independent subsystems**.

    Each component is solved on its own equilibrium EOS at its own partial density and the common
    temperature; pressures add and specific energies are mass-weighted. That is exact for the
    thermal terms and for each component's internal chemistry, and it is **wrong about exactly one
    thing**: the free hydrogen the H2 contributes does not shift water's own dissociation
    equilibrium, because the two solves never see each other. Le Chatelier says the real mixture
    holds slightly *more* of water's bond energy than this gives, so the mixture's conversion here
    is an optimistic edge and is labelled an estimate everywhere it appears.

    An honest fix is a joint H/O equilibrium solve over one nuclei budget, which is a new EOS
    rather than a composition of two.
    """
    h_share = 1.0 - water_share

    def pe(rho: float, temp: float) -> tuple[float, float]:
        p_w, e_w = eos_water.pressure_energy(water_share * rho, temp)
        if h_share <= 0.0:
            return p_w, e_w
        p_h, e_h = hydrogen.pressure_energy(h_share * rho, temp)
        return p_w + p_h, water_share * e_w + h_share * e_h

    def held(rho: float, temp: float) -> float:
        water_part = water_share * eos_water.bond_energy_held(water_share * rho, temp)
        if h_share <= 0.0:
            return water_part
        return water_part + h_share * _hydrogen_held(h_share * rho, temp)

    def cs(rho: float, temp: float) -> float:
        return eos_water.sound_speed_fd(pe, rho, temp)

    def da(rho: float, temp: float, tau_exp: float) -> tuple[float, float]:
        return _water_damkohler(water_share * rho, temp, tau_exp)

    full = (
        water_share * eos_water.FULL_ATOMIZATION_ENERGY
        + h_share * hydrogen.FULL_DISSOCIATION_ENERGY
    )
    return Fluid(f"water{water_share:.0%}+H2{h_share:.0%}", pe, cs, held, full, da, estimated=True)


# ---- Geometry: a cone rather than a length ------------------------------------------------------


def throat_radius(throat_area: float) -> float:
    """Radius [m] of a circular throat of this area."""
    return math.sqrt(throat_area / math.pi)


def station_x(
    area: float, throat_area: float, half_angle_deg: float = CONE_HALF_ANGLE_DEG
) -> float:
    """Axial distance [m] from the throat to the station of physical area `area`, on a cone."""
    return (math.sqrt(area / math.pi) - throat_radius(throat_area)) / math.tan(
        math.radians(half_angle_deg)
    )


def nozzle_length(throat_area: float, half_angle_deg: float = CONE_HALF_ANGLE_DEG) -> float:
    """Throat-to-bore length [m] of the conical nozzle this grid assumes."""
    return station_x(BORE_AREA, throat_area, half_angle_deg)


# ---- One point of the grid ----------------------------------------------------------------------


@dataclass(frozen=True)
class GridPoint:
    """One `(fluid, T_c, V, A*)` cell, on the equilibrium branch."""

    fluid: str
    estimated: bool
    temp_c: float
    volume: float
    throat_area: float
    area_ratio: float
    slug_ratio: float
    rho_c: float
    pressure_c: float
    energy_c: float
    #: Fraction of full atomisation charged in the chamber -- the ceiling on what can return.
    store_charged: float
    exit_temp: float
    exit_rho: float
    exit_speed: float
    #: `u_e^2 / (2u)`: the share of the chamber's energy density that left as directed motion.
    #: **Legitimate only where the freeze race is won.** Past the freeze station the equilibrium
    #: branch keeps releasing chemistry the real flow no longer has time to release, so read
    #: `conversion_capped` wherever `verdict` is not "equilibrium".
    conversion: float
    #: Lower bound on the conversion once the freeze is charged for: the equilibrium branch minus
    #: the chemical energy it releases *downstream* of the freeze station. It is a bound rather
    #: than the answer because the frozen gas is also colder, so part of that energy would not
    #: have reached the jet as thermal enthalpy either -- the true value sits between this and
    #: `conversion`, and the two coincide wherever nothing freezes.
    conversion_capped: float
    #: Fraction of the *charged* store handed back before the exit plane.
    store_returned: float
    store_held_exit: float
    #: Minimum Damkoehler anywhere in the nozzle, conservative edge, and its margin in decades
    #: against the equilibrium threshold.
    min_damkohler: float
    margin_decades: float
    verdict: str
    #: Area ratio where `Da` first falls through the threshold, or `None`.
    freeze_area_ratio: float | None
    freeze_temp: float | None
    nozzle_length: float
    #: Chamber emptying time [s] against the 400 ms pulse period.
    blowdown_s: float
    #: `True` where the exit is below the temperature at which `eos_methane` stops having a
    #: condensed carbon phase to offer. Water carries `False` by construction -- it has none.
    below_carbon_floor: bool

    @property
    def exit_speed_capped(self) -> float:
        """Exhaust speed [m/s] on `conversion_capped` -- the bound, not the equilibrium branch."""
        return math.sqrt(2.0 * self.conversion_capped * self.energy_c)

    @property
    def isp(self) -> chamber.IspLedger:
        """Every specific impulse for this cell, on the **capped** conversion.

        The algebra -- and in particular the sign of the head-on momentum term -- lives in
        `chamber.IspLedger`, shared with the propellant ladder so the two cannot diverge.
        """
        return chamber.isp_ledger(self.exit_speed_capped, self.slug_ratio)

    @property
    def exhaust_speed_ideal(self) -> float:
        """`w / sqrt(1+k)` [m/s] -- ADR-0016's full-conversion identity."""
        return self.isp.ideal_exhaust_speed

    @property
    def isp_true(self) -> float:
        """**Real Isp** [s]: per kilogram of everything expelled, on the capped conversion."""
        return self.isp.isp_true

    @property
    def isp_carried_gross(self) -> float:
        """Per kilogram of carried slug [s], **credit only, before the momentum debit**.

        **Not the figure of merit on a head-on burn** and must not be quoted as one. Kept because
        it is the quantity W5 and W8 first published, and because `isp_effective` is built on it.
        """
        return self.isp.isp_carried_gross

    @property
    def isp_effective(self) -> float:
        """**THE figure of merit** [s]: `[ (1+k) u_e - w ] / (k g0)`.

        The free impactor credited *and* the head-on momentum debited -- the companion's own
        head-on form. See `chamber.IspLedger` for the derivation and the sign warning.
        """
        return self.isp.isp_effective

    @property
    def eta_jet(self) -> float:
        """The companion's jet efficiency: exhaust speed over the loss-free one-axis speed.

        `u_e / (w / sqrt(1+k))`, which is identically **`sqrt(conversion)`** -- because
        `(1+k) u = w^2/2` makes `w^2/(1+k) = 2u`, so the ratio of the squares is `u_e^2/(2u)`.
        That identity is the bridge between this repository's convention-free conversion fraction
        and the companion's `eta_jet`, and it is pinned in the tests.
        """
        return self.isp.eta_jet

    @property
    def thrust_floor(self) -> float:
        """`1/sqrt(1+k)`: the `eta_jet` below which this burn pushes the ship *backwards*."""
        return self.isp.thrust_floor

    @property
    def momentum_debit(self) -> float:
        """The incoming momentum charged against thrust, in seconds [s].

        `w / (k g0)`. It is the *same* momentum for every configuration, so it falls hardest on a
        lean charge: methane at 12 kK pays 422 s of it, water at 6 kK only 85 s.
        """
        return self.isp.momentum_debit

    @property
    def blowdown_fits(self) -> bool:
        """Whether the chamber empties inside the pulse period."""
        return self.blowdown_s <= PULSE_PERIOD


@dataclass(frozen=True)
class Station:
    """One re-timed station of the isentrope, with the freeze race called."""

    area_ratio: float
    x: float
    time: float
    rho: float
    temp: float
    pressure: float
    speed: float
    tau_expansion: float
    damkohler: float
    damkohler_optimistic: float
    store_held: float


def isentrope(
    fluid: Fluid,
    rho_c: float,
    temp_c: float,
    area_ratio_max: float,
    steps: int = 220,
    expansion_ratio: float = 4.0e5,
) -> list[expansion.CoolingRow]:
    """The supersonic branch out to `area_ratio_max`, geometry-free.

    `length` is passed as 1.0 and the returned clock is discarded -- `retime` rebuilds it on the
    cone. What survives from `cooling_history` is the part that depends on no geometry at all:
    the state at each area ratio.
    """
    return expansion.cooling_history(
        rho_c,
        temp_c,
        fluid.pressure_energy,
        fluid.sound_speed,
        area_ratio_max,
        1.0,
        steps=steps,
        expansion_ratio=expansion_ratio,
    )


def _mix(a: float, b: float, w: float) -> float:
    """Linear blend from `a` to `b` at weight `w`."""
    return a + w * (b - a)


def _interpolate(rows: list[expansion.CoolingRow], ratio: float) -> expansion.CoolingRow:
    """The row at `ratio`, linearly interpolated between the two samples bracketing it."""
    for lo, hi in pairwise(rows):
        if lo.area_ratio <= ratio <= hi.area_ratio:
            span = hi.area_ratio - lo.area_ratio
            w = 0.0 if span <= 0.0 else (ratio - lo.area_ratio) / span
            return expansion.CoolingRow(
                time=0.0,
                x=0.0,
                area_ratio=ratio,
                rho=_mix(lo.rho, hi.rho, w),
                temp=_mix(lo.temp, hi.temp, w),
                pressure=_mix(lo.pressure, hi.pressure, w),
                energy=_mix(lo.energy, hi.energy, w),
                speed=_mix(lo.speed, hi.speed, w),
                mach=_mix(lo.mach, hi.mach, w),
            )
    raise ValueError(f"A/A* = {ratio} is outside the isentrope [1, {rows[-1].area_ratio}]")


def retime(
    fluid: Fluid,
    rows: list[expansion.CoolingRow],
    throat_area: float,
    exit_ratio: float,
    half_angle_deg: float = CONE_HALF_ANGLE_DEG,
) -> list[Station]:
    """Cut the isentrope at `exit_ratio` and put a conical clock on it.

    `tau_exp = dt / d ln(rho)` across each interval -- the e-folding time of the density, which is
    what a three-body rate races. Degenerate intervals are dropped for the reason `freeze._stations`
    records: a zero-width step reads as `Da = 0` at the densest station in the nozzle.
    """
    kept = [r for r in rows if r.area_ratio < exit_ratio]
    kept.append(_interpolate(rows, exit_ratio))
    tan = math.tan(math.radians(half_angle_deg))
    r_star = throat_radius(throat_area)

    xs = [(math.sqrt(throat_area * r.area_ratio / math.pi) - r_star) / tan for r in kept]
    times = [0.0]
    for i in range(1, len(kept)):
        times.append(
            times[-1] + (xs[i] - xs[i - 1]) * 0.5 * (1.0 / kept[i].speed + 1.0 / kept[i - 1].speed)
        )

    out: list[Station] = []
    for i in range(1, len(kept)):
        prev, row = kept[i - 1], kept[i]
        if not (prev.rho > row.rho > 0.0) or times[i] <= times[i - 1]:
            continue
        tau = (times[i] - times[i - 1]) / math.log(prev.rho / row.rho)
        da_opt, da_con = fluid.damkohler(row.rho, row.temp, tau)
        out.append(
            Station(
                area_ratio=row.area_ratio,
                x=xs[i],
                time=times[i],
                rho=row.rho,
                temp=row.temp,
                pressure=row.pressure,
                speed=row.speed,
                tau_expansion=tau,
                damkohler=da_con,
                damkohler_optimistic=da_opt,
                store_held=fluid.held(row.rho, row.temp) / fluid.full_atomisation,
            )
        )
    return out


def blowdown_time(rho_c: float, sound_c: float, volume: float, throat_area: float) -> float:
    """Chamber emptying time [s] through a choked throat -- `chamber.acoustic_check`'s clock.

    `t ~ M / (C rho* a* A*)` with `C` the standard choked-flow coefficient at `gamma = 1.2`. The
    same 0.65 `chamber.py` uses, restated here rather than imported so the two studies' geometries
    cannot silently merge (ADR-0050).
    """
    gamma = 1.2
    coeff = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    flux = coeff * rho_c * sound_c * throat_area
    return math.inf if flux <= 0.0 else rho_c * volume / flux


def grid_point(
    fluid: Fluid,
    temp_c: float,
    volume: float,
    throat_area: float,
    rows: list[expansion.CoolingRow],
    slug_ratio: float,
    rho_c: float,
    energy_c: float,
    sound_c: float,
    charged: float,
    half_angle_deg: float = CONE_HALF_ANGLE_DEG,
) -> GridPoint:
    """Assemble one cell from a pre-computed isentrope."""
    exit_ratio = BORE_AREA / throat_area
    stations = retime(fluid, rows, throat_area, exit_ratio, half_angle_deg)
    worst = min(stations, key=lambda s: s.damkohler)
    frozen = next((s for s in stations if s.damkohler < freeze.DA_EQUILIBRIUM), None)
    cutoff = frozen or stations[-1]
    exit_row = stations[-1]
    return GridPoint(
        fluid=fluid.name,
        estimated=fluid.estimated,
        temp_c=temp_c,
        volume=volume,
        throat_area=throat_area,
        area_ratio=exit_ratio,
        slug_ratio=slug_ratio,
        rho_c=rho_c,
        pressure_c=fluid.pressure_energy(rho_c, temp_c)[0],
        energy_c=energy_c,
        store_charged=charged,
        exit_temp=exit_row.temp,
        exit_rho=exit_row.rho,
        exit_speed=exit_row.speed,
        conversion=exit_row.speed**2 / (2.0 * energy_c),
        conversion_capped=max(
            0.0,
            exit_row.speed**2 / 2.0
            - max(0.0, cutoff.store_held - exit_row.store_held) * fluid.full_atomisation,
        )
        / energy_c,
        store_returned=0.0 if charged <= 0.0 else max(0.0, (charged - cutoff.store_held) / charged),
        store_held_exit=exit_row.store_held,
        min_damkohler=worst.damkohler,
        margin_decades=math.log10(worst.damkohler / freeze.DA_EQUILIBRIUM)
        if worst.damkohler > 0.0
        else -math.inf,
        verdict=freeze.verdict(worst.damkohler),
        freeze_area_ratio=None if frozen is None else frozen.area_ratio,
        freeze_temp=None if frozen is None else frozen.temp,
        nozzle_length=nozzle_length(throat_area, half_angle_deg),
        blowdown_s=blowdown_time(rho_c, sound_c, volume, throat_area),
        below_carbon_floor=fluid.name == "methane" and exit_row.temp < CARBON_FLOOR_K,
    )


def column(
    fluid: Fluid,
    temp_c: float,
    volume: float,
    throats: tuple[float, ...] = THROATS,
    half_angle_deg: float = CONE_HALF_ANGLE_DEG,
) -> list[GridPoint]:
    """Every throat at one `(fluid, T_c, V)`, off a single isentrope."""
    k, rho_c, energy_c = propellants.solve_slug_ratio(fluid.pressure_energy, temp_c, volume)
    sound_c = fluid.sound_speed(rho_c, temp_c)
    charged = fluid.held(rho_c, temp_c) / fluid.full_atomisation
    rows = isentrope(fluid, rho_c, temp_c, BORE_AREA / min(throats))
    return [
        grid_point(
            fluid, temp_c, volume, a, rows, k, rho_c, energy_c, sound_c, charged, half_angle_deg
        )
        for a in throats
    ]


def surface(
    fluids: tuple[Fluid, ...] = GRID_FLUIDS,
    temperatures: tuple[float, ...] = TEMPERATURES,
    volumes: tuple[float, ...] = VOLUMES,
    throats: tuple[float, ...] = THROATS,
) -> list[GridPoint]:
    """The whole grid."""
    return [
        p
        for fluid in fluids
        for temp in temperatures
        for volume in volumes
        for p in column(fluid, temp, volume, throats)
    ]


# ---- The curve the ask calls the most valuable thing here ----------------------------------------


@dataclass(frozen=True)
class HeldPoint:
    """`held` at one `(T, rho)`, plus the local density exponent."""

    fluid: str
    temp: float
    rho: float
    held: float
    #: `d ln(held) / d ln(rho)` evaluated locally by central difference. **The paper side is
    #: carrying a global -0.21 from a three-point fit; this is the real thing, and it is not
    #: constant.**
    exponent: float


#: Density decade the exit states of this grid actually occupy [kg/m^3].
HELD_DENSITIES = (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
HELD_TEMPERATURES = (2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 6000.0, 8000.0)


def held_exponent(fluid: Fluid, rho: float, temp: float, step: float = 1.05) -> float:
    """`d ln(held) / d ln(rho)` by central difference in log density."""
    lo = fluid.held(rho / step, temp) / fluid.full_atomisation
    hi = fluid.held(rho * step, temp) / fluid.full_atomisation
    if lo <= 0.0 or hi <= 0.0:
        return math.nan
    return math.log(hi / lo) / (2.0 * math.log(step))


def held_surface(
    fluids: tuple[Fluid, ...] = GRID_FLUIDS,
    temperatures: tuple[float, ...] = HELD_TEMPERATURES,
    densities: tuple[float, ...] = HELD_DENSITIES,
) -> list[HeldPoint]:
    """`held(T, rho)` computed directly, retiring the paper side's `held_W9(T) * rho^-0.21`."""
    return [
        HeldPoint(
            fluid=f.name,
            temp=t,
            rho=r,
            held=f.held(r, t) / f.full_atomisation,
            exponent=held_exponent(f, r, t),
        )
        for f in fluids
        for t in temperatures
        for r in densities
    ]


@dataclass(frozen=True)
class TemperaturePoint:
    """`T_e/T_c` at one area ratio, with the local slope the ask asks for."""

    fluid: str
    temp_c: float
    area_ratio: float
    exit_temp: float
    ratio: float
    exit_rho: float
    #: `d ln(T) / d ln(A/A*)`. Ideal monatomic frozen flow is -0.364 (well, `-(gamma-1)*2/(...)`;
    #: the paper side quotes -0.364 and that is the number being compared against).
    slope: float


#: The exponent the paper side is extrapolating, for the comparison.
PAPER_TEMPERATURE_EXPONENT = -0.138
IDEAL_TEMPERATURE_EXPONENT = -0.364


def temperature_curve(
    fluid: Fluid,
    temp_c: float,
    volume: float = 100.0,
    area_ratio_max: float = 1000.0,
) -> list[TemperaturePoint]:
    """`T_e(A/A*)` out to 1000, with the local slope -- "where does the curve steepen back?"."""
    _k, rho_c, _u = propellants.solve_slug_ratio(fluid.pressure_energy, temp_c, volume)
    rows = isentrope(fluid, rho_c, temp_c, area_ratio_max, expansion_ratio=2.0e6)
    out: list[TemperaturePoint] = []
    for prev, row, nxt in zip(rows, rows[1:], rows[2:], strict=False):
        if prev.area_ratio <= 1.0 or prev.temp <= 0.0:
            continue
        slope = math.log(nxt.temp / prev.temp) / math.log(nxt.area_ratio / prev.area_ratio)
        out.append(
            TemperaturePoint(
                fluid=fluid.name,
                temp_c=temp_c,
                area_ratio=row.area_ratio,
                exit_temp=row.temp,
                ratio=row.temp / temp_c,
                exit_rho=row.rho,
                slope=slope,
            )
        )
    return out


# ---- Output --------------------------------------------------------------------------------------

CSV_HEADER = (
    "fluid,estimated,temp_c_k,volume_m3,throat_area_m2,area_ratio,slug_ratio,rho_c,pressure_c_bar,"
    "u_j_kg,store_charged,exit_temp_k,exit_rho,exit_speed_m_s,conversion,conversion_capped,"
    "store_returned,"
    "store_held_exit,min_damkohler,margin_decades,verdict,freeze_area_ratio,freeze_temp_k,"
    "nozzle_length_m,blowdown_ms,blowdown_fits,eta_jet,thrust_floor,isp_true_s,"
    "isp_effective_s,isp_carried_gross_s,momentum_debit_s,below_carbon_floor\n"
)


def write(points: list[GridPoint], path: Path = DEFAULT_OUTPUT) -> None:
    """Write the grid. Committed: the asks document is read in the paper repository."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(CSV_HEADER)
        for p in points:
            fh.write(
                f"{p.fluid},{p.estimated},{p.temp_c:g},{p.volume:g},{p.throat_area:g},"
                f"{p.area_ratio:.3f},{p.slug_ratio:.4f},{p.rho_c:.6e},{p.pressure_c / 1e5:.4f},"
                f"{p.energy_c:.6e},{p.store_charged:.5f},{p.exit_temp:.2f},{p.exit_rho:.6e},"
                f"{p.exit_speed:.2f},{p.conversion:.5f},{p.conversion_capped:.5f},"
                f"{p.store_returned:.5f},"
                f"{p.store_held_exit:.5f},{p.min_damkohler:.6e},{p.margin_decades:.4f},"
                f"{p.verdict},"
                f"{'' if p.freeze_area_ratio is None else f'{p.freeze_area_ratio:.3f}'},"
                f"{'' if p.freeze_temp is None else f'{p.freeze_temp:.1f}'},"
                f"{p.nozzle_length:.3f},{p.blowdown_s * 1e3:.2f},{p.blowdown_fits},"
                f"{p.eta_jet:.5f},{p.thrust_floor:.5f},{p.isp_true:.1f},"
                f"{p.isp_effective:.1f},{p.isp_carried_gross:.1f},"
                f"{p.momentum_debit:.1f},{p.below_carbon_floor}\n"
            )


def write_held(points: list[HeldPoint], path: Path = HELD_OUTPUT) -> None:
    """Write the `held(T, rho)` surface -- the curve N16 calls the most valuable thing here."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("fluid,temp_k,rho,held,d_ln_held_d_ln_rho\n")
        for p in points:
            fh.write(f"{p.fluid},{p.temp:g},{p.rho:g},{p.held:.6f},{p.exponent:.4f}\n")


def write_temperature(points: list[TemperaturePoint], path: Path = TEMPERATURE_OUTPUT) -> None:
    """Write `T_e(A/A*)` and its local slope."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("fluid,temp_c_k,area_ratio,exit_temp_k,ratio,exit_rho,slope\n")
        for p in points:
            fh.write(
                f"{p.fluid},{p.temp_c:g},{p.area_ratio:.4f},{p.exit_temp:.2f},{p.ratio:.5f},"
                f"{p.exit_rho:.6e},{p.slope:.4f}\n"
            )


def _print_grid(points: list[GridPoint]) -> None:
    """The grid, one block per `(fluid, T_c)`, throats across and volumes down."""
    keys = sorted({(p.fluid, p.temp_c) for p in points}, key=lambda kv: (kv[0], kv[1]))
    for fluid, temp in keys:
        block = [p for p in points if p.fluid == fluid and p.temp_c == temp]
        throats = sorted({p.throat_area for p in block}, reverse=True)
        print(f"\n--- {fluid}, chamber {temp:.0f} K: conversion fraction ---")
        print("  V [m3] " + "".join(f"{a:>9.2f}" for a in throats) + "   <- throat A* [m^2]")
        for volume in sorted({p.volume for p in block}):
            row = {p.throat_area: p for p in block if p.volume == volume}
            print(
                f"{volume:8.0f} "
                + "".join(f"{row[a].conversion:9.3f}" for a in throats if a in row)
            )
        print("  exit T [K]")
        for volume in sorted({p.volume for p in block}):
            row = {p.throat_area: p for p in block if p.volume == volume}
            print(
                f"{volume:8.0f} " + "".join(f"{row[a].exit_temp:9.0f}" for a in throats if a in row)
            )
        print("  conversion, capped at the freeze station")
        for volume in sorted({p.volume for p in block}):
            row = {p.throat_area: p for p in block if p.volume == volume}
            print(
                f"{volume:8.0f} "
                + "".join(f"{row[a].conversion_capped:9.3f}" for a in throats if a in row)
            )
        print("  Isp EFFECTIVE [s] -- free impactor credited, head-on momentum debited.")
        for volume in sorted({p.volume for p in block}):
            row = {p.throat_area: p for p in block if p.volume == volume}
            print(
                f"{volume:8.0f} "
                + "".join(f"{row[a].isp_effective:9.0f}" for a in throats if a in row)
            )
        print("  freeze margin [decades, conservative edge]")
        for volume in sorted({p.volume for p in block}):
            row = {p.throat_area: p for p in block if p.volume == volume}
            print(
                f"{volume:8.0f} "
                + "".join(f"{row[a].margin_decades:9.2f}" for a in throats if a in row)
            )


def _print_held(points: list[HeldPoint]) -> None:
    """The held surface and the exponent the paper side assumed."""
    for fluid in sorted({p.fluid for p in points}):
        block = [p for p in points if p.fluid == fluid]
        densities = sorted({p.rho for p in block})
        print(f"\n--- {fluid}: held(T, rho), fraction of full atomisation still locked ---")
        print("   T [K] " + "".join(f"{r:>8g}" for r in densities) + "   <- rho [kg/m^3]")
        for temp in sorted({p.temp for p in block}):
            row = {p.rho: p for p in block if p.temp == temp}
            print(f"{temp:8.0f} " + "".join(f"{row[r].held:8.3f}" for r in densities if r in row))
        print("   d ln(held) / d ln(rho)  -- the paper side assumes a flat -0.21")
        for temp in sorted({p.temp for p in block}):
            row = {p.rho: p for p in block if p.temp == temp}
            print(
                f"{temp:8.0f} " + "".join(f"{row[r].exponent:8.3f}" for r in densities if r in row)
            )


def main() -> None:
    """Run N16's grid, the held surface, and the exit-temperature curve."""
    print("N16: conversion fraction over chamber T, chamber volume and throat area")
    print(
        f"3 m bore ({BORE_AREA:.2f} m^2 exit), conical nozzle at {CONE_HALF_ANGLE_DEG:.0f} deg, "
        f"{chamber.IMPACTOR_MASS:.0f} kg at {chamber.CLOSING_SPEED / 1e3:.0f} km/s\n"
    )
    print(f"{'A* [m^2]':>9} {'A/A*':>8} {'nozzle L [m]':>13}")
    for a in THROATS:
        print(f"{a:9.2f} {BORE_AREA / a:8.1f} {nozzle_length(a):13.2f}")

    points = surface()
    _print_grid(points)

    print("\n\n=== The curve N16 calls the most valuable thing here: held(T, rho) ===")
    held = held_surface()
    _print_held(held)

    print("\n\n=== Exit temperature against area ratio, and where it steepens ===")
    curves: list[TemperaturePoint] = []
    for fluid in GRID_FLUIDS:
        for temp_c in (8000.0, 10000.0):
            curve = temperature_curve(fluid, temp_c)
            curves.extend(curve)
            print(f"\n--- {fluid.name}, chamber {temp_c:.0f} K, 100 m^3 ---")
            print(f"{'A/A*':>8} {'T_e [K]':>9} {'T_e/T_c':>9} {'rho_e':>10} {'slope':>8}")
            for target in (4.0, 10.0, 30.0, 100.0, 140.0, 300.0, 565.0, 1000.0):
                near = min(curve, key=lambda p: abs(math.log(p.area_ratio / target)))
                print(
                    f"{near.area_ratio:8.1f} {near.exit_temp:9.0f} {near.ratio:9.4f} "
                    f"{near.exit_rho:10.3e} {near.slope:8.3f}"
                )

    print(
        f"\npaper-side fit slope {PAPER_TEMPERATURE_EXPONENT:+.3f}, "
        f"ideal frozen {IDEAL_TEMPERATURE_EXPONENT:+.3f}"
    )

    print("\n\n=== The configuration N16 asks to price first ===")
    named = mixture()
    first = column(named, 8000.0, 100.0, throats=(0.1,))[0]
    beside = column(METHANE, 8000.0, 100.0, throats=(0.1,))[0]
    plain_water = column(WATER, 8000.0, 100.0, throats=(0.1,))[0]
    print(
        f"{'fluid':22} {'k':>7} {'p_c bar':>9} {'exit T':>7} {'exit rho':>10} {'conv':>6} "
        f"{'Isp eff':>8} {'margin':>7}"
    )
    for p in (first, plain_water, beside):
        flag = " *ESTIMATE" if p.estimated else ""
        print(
            f"{p.fluid:22} {p.slug_ratio:7.2f} {p.pressure_c / 1e5:9.0f} {p.exit_temp:7.0f} "
            f"{p.exit_rho:10.3e} {p.conversion:6.3f} {p.isp_effective:8.0f} "
            f"{p.margin_decades:7.2f}{flag}"
        )

    write(points)
    write_held(held)
    write_temperature(curves)
    print(f"\nwrote {DEFAULT_OUTPUT}, {HELD_OUTPUT}, {TEMPERATURE_OUTPUT}")


if __name__ == "__main__":
    main()
