"""How fast does the snowplow front spread, and where does it first touch the wall? (R9)

R9 is worth 3 T of peak field. ADR-0012 caps the nozzle's peak by observing that the field's job
at the chamber is wall standoff and the front does not reach the wall for the first several
metres, so **the station where it first touches sets the cap**. The paper brackets the spreading
speed and the bracket moves that station from 6.0 m (sound speed, asking 9.0 T) to 3.26 m
(1.9x sound speed, asking 11.8 T). It adopted the safe end and would like the cheap one.

# Both of the paper's front numbers are the strong-shock piston solution, and they check out

`sec:needle_through_fog` carries two figures this repository is cited for: a 94 600 K shocked
layer and a 21.1 km/s spreading speed. Neither was ever reproduced here. They are both **the
strong-shock piston solution on `eos_water`**, and they are exact.

A piston driven into cold gas at speed `v` leaves the shocked gas moving *at the piston speed*,
so energy conservation across the shock gives its specific internal energy directly:

    e_shocked = v^2 / 2

At the cold leg's 45.58 km/s and four-fold compression of the flown bag that is 1038.8 MJ/kg,
which `eos_water` inverts to **94 632 K**, whose sound speed is **21.12 km/s**. Against the
paper's 94 600 K and 21.1 km/s: 0.03% and 0.1%. The compression ratio barely enters -- `c_s` runs
20.8 to 21.3 km/s across a 2x to 16x bracket -- so the agreement is not a tuned one.

# The correction that follows, and it goes R9's way

`coupling.snowplow` integrates the growing front the paper's own `m(x) = m_0 + rho A x` does not:

    m v = m_0 w        momentum
    dm/dx = rho pi min(r, R_bag)^2
    dr/dx = c_exp / v

ADR-0012 draws a **straight cone** off the third line, `tan(theta) = c_exp/w` fixed at the entry
speed. Two things move once the front is integrated instead, and **they are not independent**:

1. The projectile decelerates. It ends carrying `k + 1` times its own mass, so `v` falls by that
   factor and a fixed `c_exp` would open the cone faster and faster.
2. **But `c_exp` is the shocked layer's sound speed, and the shock is driven by the same `v`.** A
   weaker shock is a cooler layer with a lower sound speed. Both halves of `c_exp/v` fall together,
   which is why the cone is nearly right at all -- and why holding `c_exp` fixed while letting `v`
   fall is an artifact rather than a correction.

**What survives is the part that does not scale.** For an ideal gas `e ~ v^2` would give
`c_exp/v` exactly constant and the cone would be exact. Water is not ideal: dissociation and
ionisation absorb energy at thresholds, so the post-shock temperature rises *more slowly* than
`v^2`, and `c_exp/v` **falls** as the projectile slows -- 0.463 at 45.58 km/s to 0.29 near
10 km/s. The front therefore opens more slowly than the cone says and **touches the wall later**,
which is the direction R9 wants.

So `c_exp` is not a parameter here. It is closed at every step from the current speed, on the same
EOS that produced the 21.1 km/s the paper quotes.

# The asymmetry R9 flags, restated

The paper takes the sound speed and calls it conservative. It is -- for the coupling argument, since
a slower front sweeps less and understates `k`. For the field cap it runs the other way, because a
faster front needs field sooner. **The same assumption is conservative for one conclusion and
anti-conservative for the other**, and that is worth carrying whatever the number turns out to be.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat import coupling, eos_water, expansion, field, snowplow

DEFAULT_OUTPUT = Path("data/results/nozzle_front.csv")

BAG_BORE_M = 3.02
"""`eq:bore_from_length`'s bag bore [m] -- the surface ADR-0012's cap table measures to."""

LINER_CHAMBER_M = 3.50
"""ADR-0011's liner at the chamber [m] -- the surface ADR-0013 proposes measuring to."""

PROJECTILE_MASS_KG = coupling.PROJECTILE_MASS
PROJECTILE_RADIUS_M = coupling.projectile_radius(PROJECTILE_MASS_KG)

SHOCK_TEMPERATURE_K = snowplow.FRONT_TEMPERATURE_K
"""Shocked-layer temperature [K] -- this repository's own 94 600 K, which the paper cites."""

PAPER_SPREAD_M_S = 21.1e3
"""The spreading speed `sec:needle_through_fog` adopts, called the shocked layer's sound speed."""

SPREAD_BRACKET: tuple[float, ...] = (1.0, 1.6, 1.9)
"""The paper's own multipliers on the sound speed, from its lateral-venting argument."""

SHOCK_COMPRESSION = 4.0
"""Density jump across the front. The strong-shock ideal value; `c_s` is weakly sensitive to it
(20.8-21.3 km/s over 2x to 16x), which `main` shows before using it.
"""


# ---- Is 21.1 km/s the sound speed? ----------------------------------------------------------


@dataclass(frozen=True)
class SoundSpeedProbe:
    """The shocked layer's sound speed at one assumed compression."""

    compression: float
    rho: float
    temp_k: float
    sound_speed_m_s: float


def shock_sound_speeds(
    *,
    ambient_rho: float = coupling.BAG_RHO,
    temp_k: float = SHOCK_TEMPERATURE_K,
    compressions: tuple[float, ...] = (2.0, 4.0, 6.0, 10.0, 16.0),
) -> list[SoundSpeedProbe]:
    """`c_s` of the shocked layer across a compression bracket, on the real water EOS.

    The compression ratio is bracketed rather than assumed: a strong shock in an ideal
    `gamma = 5/3` gas gives 4, but dissociation and ionisation are energy sinks that raise it, and
    radiative losses raise it further. The sound speed is only weakly sensitive to it, which is
    the point of showing the bracket.
    """
    return [
        SoundSpeedProbe(
            compression=c,
            rho=ambient_rho * c,
            temp_k=temp_k,
            sound_speed_m_s=eos_water.sound_speed(ambient_rho * c, temp_k),
        )
        for c in compressions
    ]


def shock_state(
    speed: float, *, ambient_rho: float = coupling.BAG_RHO, compression: float = SHOCK_COMPRESSION
) -> tuple[float, float]:
    """`(T, c_s)` of the layer shocked by a piston advancing at `speed`, on the real EOS.

    Strong-shock piston relation `e = v^2/2`, inverted for temperature at the compressed density.
    This is where the paper's 94 600 K and 21.1 km/s come from; see the module docstring.
    """
    rho = ambient_rho * compression
    energy = 0.5 * speed * speed
    temp = expansion.temperature_at(rho, energy, eos_water.pressure_energy)
    return temp, eos_water.sound_speed(rho, temp)


def spread_speed_table(
    speeds: Sequence[float],
    *,
    ambient_rho: float = coupling.BAG_RHO,
    compression: float = SHOCK_COMPRESSION,
) -> Callable[[float], float]:
    """Tabulate `c_s(v)` once and return a fast interpolant.

    The integrator takes 20 000 steps and the EOS inversion is this package's hot spot, so the
    closure is tabulated on a log grid and interpolated rather than solved per step. The relation
    is smooth and monotone in `v`, so interpolation costs nothing in accuracy -- pinned by a test
    against the direct solve.
    """
    grid = sorted(speeds)
    values = [shock_state(v, ambient_rho=ambient_rho, compression=compression)[1] for v in grid]

    def interpolate(v: float) -> float:
        return float(np.interp(v, grid, values))

    return interpolate


# ---- The front's trajectory --------------------------------------------------------------------


@dataclass(frozen=True)
class FrontStation:
    """One step of the integrated snowplow."""

    x_m: float
    radius_m: float
    speed_m_s: float
    swept_kg: float
    slug_ratio: float


@dataclass(frozen=True)
class FrontRun:
    """One integration, with the two contact stations R9 and R15 disagree about."""

    closing_speed_m_s: float
    spread_m_s: float
    wall_m: float
    stations: tuple[FrontStation, ...]
    contact_x_m: float | None
    """Station where `r_front` first reaches `wall_m`, or `None` if it never does."""
    cone_contact_x_m: float
    """Where ADR-0012's fixed-angle cone puts the same contact -- the number being corrected."""
    slug_ratio: float
    exit_speed_m_s: float


def integrate(
    closing_speed: float,
    spread_multiple: float,
    *,
    wall_m: float,
    spread_of_speed: Callable[[float], float] | None = None,
    fixed_cone: bool = False,
    front_radius_0: float = PROJECTILE_RADIUS_M,
    bulk_density: float = coupling.BAG_RHO,
    bag_radius: float = coupling.BORE_RADIUS,
    length: float = coupling.BORE_LENGTH,
    projectile_mass: float = PROJECTILE_MASS_KG,
    steps: int = 20000,
    sample_every: int = 200,
) -> FrontRun:
    """Integrate `coupling.snowplow`'s system, keeping the trajectory rather than the endpoint.

    **`c_exp` is closed on the current speed**, not held fixed: `spread_multiple` times the sound
    speed of the layer the *current* speed shocks. `fixed_cone=True` freezes it at the entry value,
    which reproduces ADR-0012's construction and is what the comparison in `main` is against.

    **Swept area is capped at the bag radius, not at `wall_m`.** The mist is inside the bag; the
    clearance gap between bag and liner is vacuum. So a front wider than the bag sweeps no more
    mass, which is exactly the reason R15 answers "the liner" -- see that item.
    """
    closure = spread_of_speed or spread_speed_table(
        [closing_speed * f for f in np.geomspace(0.02, 1.0, 40)]
    )
    entry_spread = spread_multiple * closure(closing_speed)

    mass, speed, radius, x = projectile_mass, closing_speed, front_radius_0, 0.0
    dx = length / steps
    contact: float | None = None
    stations: list[FrontStation] = []

    for i in range(steps + 1):
        if i % sample_every == 0 or (contact is None and radius >= wall_m):
            stations.append(
                FrontStation(
                    x_m=x,
                    radius_m=radius,
                    speed_m_s=speed,
                    swept_kg=mass - projectile_mass,
                    slug_ratio=(mass - projectile_mass) / projectile_mass,
                )
            )
        if contact is None and radius >= wall_m:
            contact = x
        if i == steps:
            break
        mass += bulk_density * math.pi * min(radius, bag_radius) ** 2 * dx
        speed = projectile_mass * closing_speed / mass
        spread = entry_spread if fixed_cone else spread_multiple * closure(speed)
        radius += spread / speed * dx
        x += dx

    cone_slope = entry_spread / closing_speed
    cone_contact = (wall_m - front_radius_0) / cone_slope if cone_slope > 0.0 else math.inf

    return FrontRun(
        closing_speed_m_s=closing_speed,
        spread_m_s=entry_spread,
        wall_m=wall_m,
        stations=tuple(stations),
        contact_x_m=contact,
        cone_contact_x_m=cone_contact,
        slug_ratio=(mass - projectile_mass) / projectile_mass,
        exit_speed_m_s=speed,
    )


def field_demanded_at(x_m: float, cap_t: float | None = None) -> float:
    """The flown profile's field [T] at a station -- what a shelf starting there has to hold.

    This is ADR-0012's construction: hold the flown profile from the contact station outward and
    flatten everything upstream of it into a shelf at the profile's value there.
    """
    profile = field.fit_profile() if cap_t is None else field.capped_profile(cap_t)
    return profile.field(max(x_m, profile.z_ref_m))


CSV_HEADER = (
    "closing_speed_km_s",
    "spread_multiple",
    "spread_km_s",
    "wall",
    "wall_m",
    "contact_x_m",
    "cone_contact_x_m",
    "field_demanded_T",
    "slug_ratio",
    "exit_speed_km_s",
)


def write_runs(rows: list[tuple[str, float, FrontRun]], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for wall_name, multiple, run in rows:
            contact = run.contact_x_m
            writer.writerow(
                [
                    f"{run.closing_speed_m_s / 1e3:g}",
                    f"{multiple:g}",
                    f"{run.spread_m_s / 1e3:.3f}",
                    wall_name,
                    f"{run.wall_m:.3f}",
                    "" if contact is None else f"{contact:.4f}",
                    f"{run.cone_contact_x_m:.4f}",
                    "" if contact is None else f"{field_demanded_at(contact):.4f}",
                    f"{run.slug_ratio:.4f}",
                    f"{run.exit_speed_m_s / 1e3:.4f}",
                ]
            )


def main() -> None:
    """R9: the front's spread, the contact station, and the field cap it sets."""
    parser = argparse.ArgumentParser(description="Snowplow front spreading and wall contact (R9)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    print("== Is 21.1 km/s the shocked layer's sound speed? (on the real water EOS) ==")
    print(f"shocked layer at {SHOCK_TEMPERATURE_K:.0f} K, ambient {coupling.BAG_RHO} kg/m^3\n")
    print(f"{'compression':>12} {'rho':>10} {'c_s':>10} {'vs paper 21.1':>14}")
    for probe in shock_sound_speeds():
        print(
            f"{probe.compression:12.1f} {probe.rho:9.3f} {probe.sound_speed_m_s / 1e3:9.2f}k "
            f"{probe.sound_speed_m_s / PAPER_SPREAD_M_S:13.3f}x"
        )

    walls = (("bag bore", BAG_BORE_M), ("liner", LINER_CHAMBER_M))
    legs = (45.58e3, 75.0e3)

    print("\n== The closure: c_exp/v is not constant, and that is the whole correction ==")
    print("An ideal gas would give e ~ v^2 and a constant ratio, i.e. an exact cone. Water's")
    print("dissociation and ionisation thresholds make T rise more slowly than v^2.\n")
    print(f"{'v [km/s]':>9} {'T shocked':>11} {'c_s':>9} {'c_s/v':>8}")
    for v in (45.58e3, 30e3, 20e3, 10e3, 5e3, 75e3):
        temp, c_s = shock_state(v)
        print(f"{v / 1e3:9.2f} {temp:10.0f}K {c_s / 1e3:8.2f}k {c_s / v:8.4f}")

    rows: list[tuple[str, float, FrontRun]] = []
    print("\n== Where the front first touches, and the shelf it sets ==")
    print("Three constructions, and the middle one is the artifact to avoid:")
    print("  'cone'   ADR-0012's straight line at the entry angle")
    print("  'frozen' integrated, but holding c_exp at its entry value while v falls")
    print("  'closed' integrated with c_exp closed on the current speed -- the answer\n")
    print(
        f"{'leg':>7} {'spread':>8} {'wall':>10} {'cone':>8} {'frozen':>9} {'closed':>9} "
        f"{'vs cone':>9} {'cone T':>8} {'closed T':>9}"
    )
    for closing in legs:
        closure = spread_speed_table([closing * f for f in np.geomspace(0.02, 1.0, 40)])
        for multiple in SPREAD_BRACKET:
            for wall_name, wall in walls:
                run = integrate(closing, multiple, wall_m=wall, spread_of_speed=closure)
                cone = integrate(
                    closing, multiple, wall_m=wall, spread_of_speed=closure, fixed_cone=True
                )
                rows.append((wall_name, multiple, run))
                contact = run.contact_x_m
                frozen_c = cone.contact_x_m
                straight = cone.cone_contact_x_m
                shown = "never" if contact is None else f"{contact:8.2f}m"
                frozen_shown = "never" if frozen_c is None else f"{frozen_c:8.2f}m"
                moves = "-" if contact is None else f"{contact - straight:+8.2f}m"
                cone_field = f"{field_demanded_at(straight):7.2f}T"
                closed_field = "-" if contact is None else f"{field_demanded_at(contact):8.2f}T"
                print(
                    f"{closing / 1e3:7.2f} {multiple:7.1f}x {wall_name:>10} "
                    f"{straight:7.2f}m {frozen_shown:>9} {shown:>9} {moves:>9} "
                    f"{cone_field:>8} {closed_field:>9}"
                )

    print("\n== The cap, read off the binding case at each bracket end ==")
    for multiple in SPREAD_BRACKET:
        for wall_name, wall in walls:
            group = [
                r
                for r in rows
                if r[0] == wall_name and r[1] == multiple and r[2].contact_x_m is not None
            ]
            if not group:
                continue
            worst = max(group, key=lambda r: field_demanded_at(r[2].contact_x_m or 0.0))
            contact = worst[2].contact_x_m or 0.0
            print(
                f"  spread {multiple:.1f}x, wall = {wall_name:>8} at {wall:.2f} m "
                f"-> contact {contact:5.2f} m, cap {field_demanded_at(contact):5.2f} T"
            )

    write_runs(rows, args.output)
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
