"""Velikhov stability along the nozzle, station by station (ADR-0038).

`conductivity.electrothermal_loop` decides one state. This module walks it along the cooling
history `expansion.history()` produces, which is what turns three exit verdicts into the
leg-level statement ADR-0038 reports: *where* a leg crosses into instability and how much of its
transit it spends there.

**The field is not an assumption.** Flux conservation gives `B*/B = A/A*`, the same relation that
makes the paper's `20 T -> 5 T` nozzle exactly area ratio 4, so the local field follows from the
area ratio the expansion already solves for.

**The current-layer thickness is the open question, and it is exposed rather than buried.**
`use_skin_depth` selects between the conservative skin estimate (shipped) and the full flux-tube
radius that this plume's low plasma beta argues for. They disagree on the verdict. See ADR-0038
Addendum 3 -- settling it is an MHD problem this repository does not solve.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache
from itertools import pairwise
from pathlib import Path

from puffsat import conductivity, expansion

#: Field at the sonic throat [T] -- the paper's ~20 T at the 3 m bore.
B_THROAT = 20.0

DEFAULT_SCAN_PATH = Path("data/results/electrothermal_scan.csv")


def local_field(area_ratio: float) -> float:
    """Local field [T] from the area ratio: `B*/B = A/A*` by flux conservation."""
    if area_ratio <= 0.0:
        raise ValueError("area_ratio must be positive")
    return B_THROAT / area_ratio


@dataclass(frozen=True)
class StabilityRow:
    """One station of one leg, with the Velikhov verdict at it."""

    closing_speed: float
    #: Time since the sonic throat [s].
    time: float
    area_ratio: float
    temp: float
    rho: float
    b_field: float
    #: Flux-tube radius [m], the outer bound on the current-layer thickness.
    length_scale: float
    loop: conductivity.ElectrothermalLoop

    @property
    def t_e(self) -> float:
        return self.loop.balance.t_e

    @property
    def elevation(self) -> float:
        """`T_e - T_gas` [K] -- what actually sets `beta_cr`."""
        return self.loop.balance.elevation

    @property
    def beta(self) -> float:
        return self.loop.screen.hall_parameter

    @property
    def critical_hall_parameter(self) -> float:
        return self.loop.critical_hall_parameter

    @property
    def e_folding_time(self) -> float:
        return self.loop.e_folding_time

    @property
    def unstable(self) -> bool:
        return self.loop.unstable


@cache
def _history(
    closing_speed: float, temp_0: float, frozen: bool, steps: int, stride: int
) -> tuple[expansion.HistoryRow, ...]:
    """`expansion.history` memoised -- both readings of the current layer walk the same history.

    The cooling history is the expensive half of this calculation and it does not depend on the
    stability model at all, so recomputing it per reading doubles the cost of every comparison
    ADR-0038 Addendum 3 asks for.
    """
    return tuple(expansion.history(closing_speed, temp_0, frozen, steps=steps, stride=stride))


def scan(
    closing_speed: float,
    temp_0: float,
    *,
    frozen: bool = False,
    use_skin_depth: bool = True,
    steps: int = 160,
    stride: int = 2,
) -> list[StabilityRow]:
    """Velikhov stability at every station of one leg's cooling history.

    The diffusion time fed to the skin estimate is the **whole leg transit**, not the time elapsed
    since the throat: the plume has been inside the field since the bag, so the field has had at
    least that long to diffuse across it. Using the elapsed-since-throat time would send the skin
    to zero at the throat and the current density to infinity with it, which is an artifact of
    where the integration starts rather than physics.

    **The resolution trade runs the opposite way to `expansion.history`'s.** There the integration
    step is set by the radiative-loss integral and the station spacing is free; here the loss
    integral is not consumed, so `steps` can halve (the crossing is unchanged to four digits from
    160 to 320), while `stride` must be *finer* -- the deliverable is a threshold crossing, and a
    coarse station grid brackets it no better than the spacing. At these settings the crossing and
    the dwell are converged; at `expansion.history`'s own default they are not.
    """
    rows = _history(closing_speed, temp_0, frozen, steps, stride)
    transit = rows[-1].row.time

    out: list[StabilityRow] = []
    for entry in rows:
        row = entry.row
        radius = expansion.plume_radius(row.area_ratio)
        b_field = local_field(row.area_ratio)
        loop = conductivity.electrothermal_loop(
            row.temp,
            row.rho,
            expansion.SEED_FRACTION,
            b_field,
            radius,
            transit,
            use_skin_depth=use_skin_depth,
        )
        out.append(
            StabilityRow(
                closing_speed=closing_speed,
                time=row.time,
                area_ratio=row.area_ratio,
                temp=row.temp,
                rho=row.rho,
                b_field=b_field,
                length_scale=radius,
                loop=loop,
            )
        )
    return out


def crossing(rows: list[StabilityRow]) -> StabilityRow | None:
    """The first unstable station, or `None` if the leg is stable throughout."""
    return next((r for r in rows if r.unstable), None)


def unstable_dwell(rows: list[StabilityRow]) -> float:
    """Time [s] the plume spends unstable, on the station clock.

    An interval counts when its **earlier** endpoint is unstable, so the instability is taken to
    begin at the first station that reports it. That is the conservative half of the bracket the
    station spacing leaves: the true crossing lies between the last stable station and the first
    unstable one, so this form approaches it from below and the `hi.unstable` form from above.
    They close on each other under refinement -- 1.698 vs 1.752 ms at the shipped resolution.

    Reported against the transit rather than against the e-folding time, because the two differ by
    ~3 orders: a station that is unstable at all is unstable hundreds of times over, so the
    question is how long the plume is exposed, not whether it has time to grow once.
    """
    return sum(hi.time - lo.time for lo, hi in pairwise(rows) if lo.unstable)


def write_scan(rows: list[StabilityRow], path: Path = DEFAULT_SCAN_PATH) -> None:
    """One CSV over all legs and both readings of the current-layer thickness."""
    lines = [
        "closing_speed_km_s,current_layer,time_s,area_ratio,temp_K,rho_kg_m3,b_field_T,"
        "length_scale_m,current_length_scale_m,t_e_K,elevation_K,beta,beta_crit,"
        "e_folding_s,unstable"
    ]
    for r in rows:
        layer = "skin" if r.loop.balance.current_length_scale < r.length_scale else "flux-tube"
        e_fold = "" if math.isinf(r.e_folding_time) else f"{r.e_folding_time:.6e}"
        lines.append(
            f"{r.closing_speed:g},{layer},{r.time:.6e},{r.area_ratio:.6f},{r.temp:.2f},"
            f"{r.rho:.6e},{r.b_field:.4f},{r.length_scale:.4f},"
            f"{r.loop.balance.current_length_scale:.4f},{r.t_e:.2f},{r.elevation:.2f},"
            f"{r.beta:.4f},{r.critical_hall_parameter:.6g},{e_fold},{int(r.unstable)}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    """ADR-0038's leg-level verdict, under both readings of the current-layer thickness."""
    everything: list[StabilityRow] = []
    for use_skin in (True, False):
        label = "skin (shipped, conservative)" if use_skin else "flux-tube radius (low-beta)"
        print(f"python: current layer = {label}")
        header = (
            f"  {'leg':>8} {'T exit':>8} {'dT_e':>8} {'beta':>7} "
            f"{'beta_cr':>10} {'unstable':>9} {'transit':>9} {'e-fold':>9}"
        )
        print(header)
        for speed, temp_0 in expansion.PLUME_STATES:
            rows = scan(speed, temp_0, use_skin_depth=use_skin)
            everything.extend(rows)
            last = rows[-1]
            dwell = unstable_dwell(rows)
            worst = min((r.e_folding_time for r in rows if r.unstable), default=math.inf)
            e_fold = "--" if math.isinf(worst) else f"{worst * 1e6:.1f}us"
            print(
                f"  {speed:6g}   {last.temp:8.0f} {last.elevation:8.1f} {last.beta:7.2f} "
                f"{last.critical_hall_parameter:10.4g} {dwell * 1e3:8.2f}ms "
                f"{last.time * 1e3:8.2f}ms {e_fold:>9}"
            )
    write_scan(everything)
    print(f"python: wrote {DEFAULT_SCAN_PATH}")
    print(
        "  the two readings disagree on the cold leg and agree on the hot legs; ADR-0038 "
        "Addendum 3 is why the disagreement is reported rather than resolved."
    )


if __name__ == "__main__":
    main()
