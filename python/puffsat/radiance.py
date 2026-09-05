"""`phi`: what share of the pulse's radiation is emitted while the plume is still inside? (R12)

The paper needs one number per leg and branch, and it says why it cannot get it from its side:
every corrected liner figure it has carries `phi` as a multiplier. Its own summary of the problem
is exact --

> only radiation emitted **while the plume is still inside** lands on the wall. Transit is a few
> milliseconds and the plume radiates for far longer downstream, into nothing.

-- and at `phi` = 1 on the equilibrium branch the liner sits past graphite's 3900 K, which is the
only place in this whole exercise where the passive-structure claim actually fails. So the number
decides a hardware question (how many refractory shields, if any) and not just a margin.

# What is integrated, and against what

Radiated power per unit mass is `e/t_rad`, with `e` the internal energy and `t_rad` the radiative
cooling time `radiation_check` already computes on the right side of the optically-thick /
free-streaming crossover. Then

    phi = integral_0^{t_exit} (e/t_rad) dt  /  integral_0^{infinity} (e/t_rad) dt.

**The clock is a parcel's own transit, not the leading edge's.** R12 asks for "before the plume's
leading edge passes the exit plane", but in a steady quasi-1D expansion every parcel has the same
history, so the meaningful residence is how long *a parcel* spends inside -- which is exactly what
`expansion.cooling_history`'s clock measures, 1.7 to 2.8 ms. Reading it as the leading edge would
be the same integral with a different name, since the leading edge is just the first parcel.

**The denominator is continued down the free jet**, on the same isentrope, with
`fireball.py`'s jet clock: past the lip the plume opens at `DIVERGENCE_HALF_ANGLE_DEG` and the
time to reach a radius is `dr/(u tan theta)`. The tail's convergence is reported rather than
assumed, because a `phi` computed against a truncated denominator is biased *high*, which is the
direction that would falsely condemn the liner.

# The one assumption that matters

**Radiation is treated as a perturbation on an adiabatic expansion.** The plume's state comes from
the isentrope, and the radiated energy is integrated against it without letting the loss feed back
on the temperature. That is the same assumption `expansion.py` makes and `radiation_check` exists
to police; it holds while the radiated fraction is small. Where it is not small the number is
reported anyway, flagged, because a `phi` from an invalid history is still the right *ratio* as
long as both halves of it are wrong the same way -- and the numerator and denominator share an
isentrope.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from puffsat import eos_water, expansion

DEFAULT_OUTPUT = Path("data/results/nozzle_phi.csv")

JET_HALF_ANGLE_DEG = 45.0
"""Free-jet divergence past the lip. `fireball.DIVERGENCE_HALF_ANGLE_DEG`, kept in step with it."""

TAIL_EXPANSION_RATIO = 8192.0
"""How far down the isentrope the denominator is carried. `fireball.EXPANSION_RATIO`."""


@dataclass(frozen=True)
class RadiatingStation:
    """One station with the two quantities the quadrature needs."""

    time_s: float
    area_ratio: float
    radius_m: float
    rho: float
    temp_k: float
    energy_j_kg: float
    cooling_time_s: float
    regime: str
    inside: bool

    @property
    def power_w_kg(self) -> float:
        """`e/t_rad` [W/kg] -- radiated power per unit mass at this station."""
        if self.cooling_time_s <= 0.0 or not math.isfinite(self.cooling_time_s):
            return 0.0
        return self.energy_j_kg / self.cooling_time_s


def _eos_for(branch: str, temp_0: float) -> tuple[expansion.Eos, expansion.SoundSpeed]:
    """The `(p, e)` and sound-speed closures for a branch of the ADR-0026 bracket."""
    if branch == "frozen":
        y = eos_water.frozen_composition(expansion.BAG_RHO, temp_0)
        return (
            lambda r, t: eos_water.pressure_energy_frozen(r, t, y),
            lambda r, t: eos_water.sound_speed_frozen(r, t, y),
        )
    return eos_water.pressure_energy, eos_water.sound_speed


def history(
    temp_0: float,
    *,
    branch: str = "equilibrium",
    area_ratio_exit: float = expansion.AREA_RATIO_EXIT,
    length_m: float = expansion.FIELD_LENGTH,
    half_angle_deg: float = JET_HALF_ANGLE_DEG,
    expansion_ratio: float = TAIL_EXPANSION_RATIO,
    steps: int = 320,
    stride: int = 2,
) -> list[RadiatingStation]:
    """Nozzle plus free jet on one isentrope, with the radiative cooling time at every station.

    The in-nozzle leg reuses `expansion.cooling_history` unchanged so the two clocks cannot drift;
    the free leg continues from the lip on the same isentrope, exactly as `fireball._history`
    does.
    """
    eos, c_s = _eos_for(branch, temp_0)
    kappa = _kappa_for(branch, temp_0)

    inside = expansion.cooling_history(
        expansion.BAG_RHO, temp_0, eos, c_s, area_ratio_exit, length_m, steps=steps
    )
    points = _longest_tail(temp_0, eos, c_s, expansion_ratio, steps)
    throat = expansion._at_throat(points)
    supersonic = [throat, *(p for p in points if p.mach > throat.mach)]

    out: list[RadiatingStation] = []
    for row in [*inside[::stride], inside[-1]]:
        out.append(
            _station(
                row.time,
                row.area_ratio,
                row.rho,
                row.temp,
                row.energy,
                kappa,
                inside=True,
            )
        )

    tan_theta = math.tan(math.radians(half_angle_deg))
    prev_radius = expansion.plume_radius(inside[-1].area_ratio)
    prev_speed = inside[-1].speed
    time = inside[-1].time
    for pt in (p for p in supersonic if p.area_ratio > area_ratio_exit):
        radius = expansion.plume_radius(pt.area_ratio)
        time += (radius - prev_radius) * 0.5 * (1.0 / pt.speed + 1.0 / prev_speed) / tan_theta
        prev_radius, prev_speed = radius, pt.speed
        out.append(_station(time, pt.area_ratio, pt.rho, pt.temp, pt.energy, kappa, inside=False))
    return out


def _longest_tail(
    temp_0: float,
    eos: expansion.Eos,
    c_s: expansion.SoundSpeed,
    expansion_ratio: float,
    steps: int,
) -> list[expansion.NozzlePoint]:
    """Continue the isentrope as far down as the EOS can follow it, halving on failure.

    **The frozen branch runs out before the equilibrium one does**, and that is physics rather
    than a numerical problem: with no recombination store to hand back, a frozen plume is at
    3000 K by the nozzle exit and falls below `eos_water`'s ~50 K partition-function floor well
    before the equilibrium branch does. A plume that cold has no internal energy left to radiate,
    so truncating there costs the denominator almost nothing -- and `PhiResult.last_decade_share`
    is reported precisely so the reader can check that rather than take it on trust.
    """
    ratio = expansion_ratio
    while ratio >= 128.0:
        try:
            return expansion.nozzle_history(
                expansion.BAG_RHO, temp_0, expansion.BAG_RHO / ratio, eos, c_s, steps
            )
        except ValueError:
            ratio *= 0.5
    raise ValueError(f"the isentrope from T0 = {temp_0} K will not solve past A/A* = {ratio}")


#: `(rho, T) -> (kappa_Planck, kappa_Rosseland)` [m^2/kg].
Kappa = Callable[[float, float], tuple[float, float]]


def _kappa_for(branch: str, temp_0: float) -> Kappa:
    """The opacity closure. Both branches read the same TOPS gray table `expansion.py` uses."""
    _ = branch, temp_0
    return expansion._tops_kappa


def _station(
    time: float,
    area_ratio: float,
    rho: float,
    temp: float,
    energy: float,
    kappa: Kappa,
    *,
    inside: bool,
) -> RadiatingStation:
    """Attach the radiative cooling time to one isentrope point."""
    radius = expansion.plume_radius(area_ratio)
    kap_p, kap_r = kappa(rho, temp)
    check = expansion.radiation_check(temp, rho, energy, kap_p, kap_r, radius)
    return RadiatingStation(
        time_s=time,
        area_ratio=area_ratio,
        radius_m=radius,
        rho=rho,
        temp_k=temp,
        energy_j_kg=energy,
        cooling_time_s=check.cooling_time,
        regime=check.regime,
        inside=inside,
    )


# ---- The quadrature ------------------------------------------------------------------------------


def radiated_energy(stations: list[RadiatingStation]) -> float:
    """`integral (e/t_rad) dt` [J/kg] over a run of stations, by trapezoid in time."""
    return sum(
        (b.time_s - a.time_s) * 0.5 * (a.power_w_kg + b.power_w_kg)
        for a, b in pairwise(stations)
        if b.time_s > a.time_s
    )


@dataclass(frozen=True)
class PhiResult:
    """`phi` for one leg and branch, with the convergence evidence beside it."""

    closing_speed_km_s: float
    branch: str
    phi: float
    inside_j_kg: float
    total_j_kg: float
    transit_ms: float
    tail_ms: float
    last_decade_share: float
    """Share of the *total* contributed by the final decade of area ratio -- the tail's weight.

    A large value means the denominator has not converged and `phi` is an upper bound.
    """
    inside_fraction_of_internal: float
    """`E_rad,inside / e_throat` -- whether the perturbative treatment is defensible."""


def phi(
    closing_speed: float, temp_0: float, *, branch: str = "equilibrium", **kwargs: float
) -> PhiResult:
    """The in-nozzle share of one leg's radiated energy."""
    stations = history(temp_0, branch=branch, **kwargs)  # type: ignore[arg-type]
    inside = [s for s in stations if s.inside]
    e_inside = radiated_energy(inside)
    e_total = radiated_energy(stations)

    decade = [s for s in stations if s.area_ratio >= 0.1 * stations[-1].area_ratio]
    tail = radiated_energy(decade) if len(decade) > 1 else 0.0

    return PhiResult(
        closing_speed_km_s=closing_speed,
        branch=branch,
        phi=e_inside / e_total if e_total > 0.0 else float("nan"),
        inside_j_kg=e_inside,
        total_j_kg=e_total,
        transit_ms=inside[-1].time_s * 1e3,
        tail_ms=(stations[-1].time_s - inside[-1].time_s) * 1e3,
        last_decade_share=tail / e_total if e_total > 0.0 else float("nan"),
        inside_fraction_of_internal=e_inside / stations[0].energy_j_kg,
    )


CSV_HEADER = (
    "closing_speed_km_s",
    "branch",
    "phi",
    "radiated_inside_J_kg",
    "radiated_total_J_kg",
    "transit_ms",
    "tail_ms",
    "last_decade_share",
    "inside_fraction_of_internal",
)


def write_phi(rows: list[PhiResult], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for r in rows:
            writer.writerow(
                [
                    f"{r.closing_speed_km_s:g}",
                    r.branch,
                    f"{r.phi:.5f}",
                    f"{r.inside_j_kg:.6e}",
                    f"{r.total_j_kg:.6e}",
                    f"{r.transit_ms:.4f}",
                    f"{r.tail_ms:.4f}",
                    f"{r.last_decade_share:.5f}",
                    f"{r.inside_fraction_of_internal:.6f}",
                ]
            )


LINER_LOAD_GW_AT_UNITY = 3.32
"""Hot-branch liner load [GW] at `phi` = 1, from R12's own table. Scales linearly in `phi`."""

GRAPHITE_CEILING_GW = 8.28
"""What the graphite liner sheds, R12's figure."""


def main() -> None:
    """R12: the in-nozzle share of the radiated energy, per leg and branch."""
    parser = argparse.ArgumentParser(description="phi, the in-nozzle radiated share (R12)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rows: list[PhiResult] = []
    for speed, temp_0 in expansion.PLUME_STATES:
        for branch in ("equilibrium", "frozen"):
            rows.append(phi(speed, temp_0, branch=branch))

    print("== R12: what share of the radiation lands while the plume is still inside ==")
    print("Numerator: a parcel's own nozzle transit. Denominator: continued down the free jet")
    print(f"to A/A* = {TAIL_EXPANSION_RATIO:g} at a {JET_HALF_ANGLE_DEG:g} deg half-angle.\n")
    print(
        f"{'leg':>8} {'branch':>12} {'phi':>8} {'transit':>9} {'jet tail':>10} "
        f"{'last decade':>12} {'E_rad/e_0 inside':>17}"
    )
    for r in rows:
        print(
            f"{r.closing_speed_km_s:8.2f} {r.branch:>12} {r.phi:8.4f} {r.transit_ms:8.3f}ms "
            f"{r.tail_ms:9.2f}ms {100 * r.last_decade_share:11.2f}% "
            f"{100 * r.inside_fraction_of_internal:16.3f}%"
        )

    print("\n== What it does to the liner (R12's own scaling, hot 3.6% branch) ==")
    print(f"{'leg':>8} {'branch':>12} {'phi':>8} {'liner load':>12} {'of ceiling':>11}")
    for r in rows:
        load = LINER_LOAD_GW_AT_UNITY * r.phi
        print(
            f"{r.closing_speed_km_s:8.2f} {r.branch:>12} {r.phi:8.4f} {load:11.3f}GW "
            f"{load / GRAPHITE_CEILING_GW:10.3f}x"
        )

    write_phi(rows, args.output)
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
