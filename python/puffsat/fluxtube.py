"""Flux-tube accounting against the real Biot-Savart solve (R2, R5, and R1's divergence term).

The paper asks for two things from the same object, and they turn out to be the same run:

- **R2** -- its flux-tube accounting redone against a real field solve rather than the paraxial
  `B_r = -(r/2) dB_z/dz` it used out to 6 m against a 3.5 m coil, which is past where paraxial is
  defensible.
- **R5** -- what share of the plume is born on a tube that will not fit inside the winding.

And R1 needs a third number from it, `<cos theta>`, which is the factor its corrected `eta_geom`
turns on. All three are properties of the same traced tubes, so they are computed together.

# Why flux tubes are the right object even though `mu` is not

`continuum.py` retires the *single-particle* invariant: the plume collides millions of times
crossing the bore, so no particle remembers its own `mu`. That says nothing against **flux
tubes**, which are a statement about the *fluid*. The relevant number is the magnetic Reynolds
number, and `data/results/cooling_history.csv` carries it: `Rm` runs into the hundreds through
the column, so the field is frozen into the plasma and a fluid element stays on its tube.

So the geometry survives the regime correction even though the energetics do not. A parcel does
**not** convert its gyration to axial motion by falling down `B` (that was `jet.py`'s error), but
it **does** travel along the tube it was born on, and the tube's inclination to the thrust axis is
a real loss of axial momentum. That is R1's divergence term, and it is the half of R1 that stands
after its own amendment.

# Where the tubes are launched, and why the answer is not the chamber plane

This looks like a detail and is not: it is a factor of five on R5's number, and it is where the
two repositories' models differ without either saying so.

**A steady nozzle has one inlet**, so every gram enters at the chamber and rides a tube that falls
the full `B_chamber/B_exit`. Launching there gives 43.8% of the plume missing a straight 3.5 m
winding. **But the plume is not fed through an inlet.** It is a 23.8 m column that already fills
the bore when it starts to leave, so a parcel at `z` = 20 m never sees the chamber field at all --
it rides the tube through its own starting point, which has only `B(20)/B_exit` left to fall. That
is the paper's own station-weighted picture, and for *geometry* it is right.

**The same weighting is wrong for energetics, and R4 rests on it.** Which tube a parcel rides is
fixed by where it starts. How much of its thermal energy becomes directed motion is fixed by the
total pressure drop it falls through, and in a collisional fluid that is not tied to its birth
station's field at all. `mu` was the thing that tied them, and `continuum.py` retires it. So
station weighting survives here and does not survive in R4.

Launches therefore cover the **whole bore volume**, `(r0, z0)` on a grid, weighted by `r0 dr0 dz0`
for uniform initial density. The uniform-density assumption is the weak step, and it is the same
one on both sides: the real mass-versus-station profile at the end of the snowplow is what R4 asks
for, and nobody has it.

# What else is assumed, and stated because it is a choice

- **The trace is exact within its own model** -- RK4 on `dr/dz = B_r/B_z` through the superposed
  loop fields, no paraxial expansion. It is a *vacuum* field: the plume's own diamagnetism is
  absent, and at low `beta` that is a small correction (`field.py`'s module docstring says why it
  is the optimistic direction).
- **`theta` is the field line's inclination**, `atan(B_r/B_z)`, which is what the paraxial formula
  approximates. Taking it from the trace rather than from `d ln B/dz` is the whole point of R2.
- **The trace stops at the exit plane**, so the `<cos theta>` reported here is the divergence
  *inside* the winding, which is near zero because a solenoid bore is very nearly axial. **The
  loss R1 is about happens downstream**, past the last coil, where the field fans; that is
  `extension.py`'s job and R11's question. What this module establishes is that the inside
  contributes nothing, so the whole term is a detachment-surface question.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat import field

DEFAULT_OUTPUT = Path("data/results/nozzle_fluxtube.csv")

BAG_RADIUS_M = 3.0
"""Bag bore radius [m] -- where the tubes are launched. `expansion.CHAMBER_RADIUS`."""

COLUMN_LENGTH_M = 23.8
"""Column length [m], the sim's value, which R14 adopts paper-side."""

LINER_GAP_M = field.LINER_GAP_M
"""Bore-to-winding standoff [m]: graphite liner, vacuum gap, aluminium shell."""

STRAIGHT_COIL_RADIUS_M = 3.5
"""The winding P8 assumed: constant radius, just outside the bag."""

N_COILS = 48
"""Coil count for the traces. Well above where ripple matters (see `field.ripple_sweep`), so the
traced geometry is the smooth grading rather than an artifact of discretisation.
"""


# ---- The two windings ------------------------------------------------------------------------


def flared_radius(
    z: float,
    profile: field.Profile,
    *,
    bag_radius: float = BAG_RADIUS_M,
    gap: float = LINER_GAP_M,
    z_chamber: float = 0.0,
) -> float:
    """ADR-0011's contour: the winding follows the plume's bounding flux tube, plus clearance.

    `r(z) = r_bag sqrt(B_chamber/B(z)) + gap`. At ADR-0012's 12 T cap this runs 3.50 m at the
    chamber to 5.17 m at the exit, which is the paper's stated flare and the check that this
    reproduces its geometry rather than inventing one.
    """
    b_chamber = profile.field(max(z_chamber, profile.z_ref_m))
    b_here = profile.field(max(z, profile.z_ref_m))
    return bag_radius * math.sqrt(b_chamber / b_here) + gap


def build_stacks(
    *,
    cap_t: float = field.ADR0012_CAP_T,
    length: float = COLUMN_LENGTH_M,
    n_coils: int = N_COILS,
) -> dict[str, tuple[field.CoilStack, field.Profile]]:
    """The windings the asks are about: P8's straight one, and the current flared/capped one."""
    flown = field.fit_profile()
    capped = field.capped_profile(cap_t)

    def flare(z: float) -> float:
        return flared_radius(min(max(z, 0.0), length), capped)

    return {
        "straight-flown": (
            field.build_winding(
                n_coils, coil_radius=STRAIGHT_COIL_RADIUS_M, length=length, profile=flown
            ),
            flown,
        ),
        "straight-capped": (
            field.build_winding(
                n_coils, coil_radius=STRAIGHT_COIL_RADIUS_M, length=length, profile=capped
            ),
            capped,
        ),
        "flared-capped": (
            field.build_winding(n_coils, length=length, profile=capped, radius_at=flare),
            capped,
        ),
    }


# ---- One traced tube -------------------------------------------------------------------------


@dataclass(frozen=True)
class Tube:
    """One flux tube, launched somewhere in the bore and followed to the exit.

    `escaped` marks a tube that reached the winding contour before the exit. Its trace past that
    point is not physical -- outside the winding the field reverses and the tracer happily follows
    it out to tens of metres -- so an escaped tube contributes to the missing fraction and is
    excluded from the `<cos theta>` average rather than being allowed to poison it.
    """

    r_launch_m: float
    z_launch_m: float
    z: tuple[float, ...]
    r: tuple[float, ...]
    b_t: tuple[float, ...]
    theta_rad: tuple[float, ...]
    escaped: bool = False
    weight: float = 1.0

    @property
    def r_exit_m(self) -> float:
        return self.r[-1]

    @property
    def theta_exit_rad(self) -> float:
        return self.theta_rad[-1]

    def radius_at(self, z: float) -> float:
        """Tube radius [m] at a station, by linear interpolation on the trace."""
        return float(np.interp(z, np.array(self.z), np.array(self.r)))

    def theta_at(self, z: float) -> float:
        """Field-line inclination [rad] at a station."""
        return float(np.interp(z, np.array(self.z), np.array(self.theta_rad)))

    @property
    def max_radius_m(self) -> float:
        return max(self.r)


def trace_tube(
    stack: field.CoilStack,
    r_launch: float,
    *,
    z_start: float = 0.0,
    z_end: float = COLUMN_LENGTH_M,
    steps: int = 400,
    weight: float = 1.0,
    clearance: float = 0.0,
) -> Tube:
    """Follow one tube from its launch point, recording radius, `|B|` and inclination.

    Uses `field.trace_line`'s RK4, then recomputes `(B_r, B_z)` on the traced path for `theta`.
    The trace is truncated at the winding contour: past it the field reverses, and following a
    tube through that would report a divergence angle that belongs to the return flux rather than
    to any plume.
    """
    if z_start >= z_end:
        b_r, b_z = stack.field(max(r_launch, 1e-9), z_start)
        return Tube(
            r_launch_m=r_launch,
            z_launch_m=z_start,
            z=(z_start,),
            r=(r_launch,),
            b_t=(math.hypot(b_r, b_z),),
            theta_rad=(math.atan2(b_r, b_z),),
            weight=weight,
        )

    try:
        path = field.trace_line(stack, r_launch, z_start, z_end=z_end, steps=steps)
    except ValueError:
        path = [(z_start, r_launch, stack.magnitude(max(r_launch, 1e-9), z_start))]

    zs: list[float] = []
    rs: list[float] = []
    bs: list[float] = []
    thetas: list[float] = []
    escaped = False
    for z, r, b in path:
        if r >= _winding_radius(stack, z) - clearance:
            escaped = True
            break
        b_r, b_z = stack.field(max(r, 1e-9), z)
        zs.append(float(z))
        rs.append(float(r))
        bs.append(float(b))
        thetas.append(math.atan2(b_r, b_z))
    if not zs:
        escaped = True
        zs, rs, bs, thetas = [z_start], [r_launch], [path[0][2]], [0.0]

    return Tube(
        r_launch_m=r_launch,
        z_launch_m=z_start,
        z=tuple(zs),
        r=tuple(rs),
        b_t=tuple(bs),
        theta_rad=tuple(thetas),
        escaped=escaped,
        weight=weight,
    )


def trace_bundle(
    stack: field.CoilStack,
    *,
    bag_radius: float = BAG_RADIUS_M,
    n_radial: int = 12,
    n_axial: int = 12,
    length: float = COLUMN_LENGTH_M,
    z_end: float = COLUMN_LENGTH_M,
    steps: int = 400,
) -> list[Tube]:
    """A bundle spanning the **bore volume**, weighted for uniform initial density.

    Radial launch points are the midpoints of `n_radial` equal-area annuli,
    `r_i = r_bag sqrt((i + 1/2)/n)`, so they carry equal radial weight. Axial launch points are
    evenly spaced midpoints, equal weight again. Every tube therefore has weight 1 and the
    averages are plain means -- the quadrature is in the placement, not in the weights.

    See the module docstring for why the launch is over the volume rather than at the chamber:
    the plume fills the column before it leaves, so a parcel at `z` never falls the chamber's
    field.
    """
    radii = [bag_radius * math.sqrt((i + 0.5) / n_radial) for i in range(n_radial)]
    stations = [length * (j + 0.5) / n_axial for j in range(n_axial)]
    return [
        trace_tube(stack, r, z_start=z0, z_end=z_end, steps=steps) for z0 in stations for r in radii
    ]


def chamber_bundle(
    stack: field.CoilStack,
    *,
    bag_radius: float = BAG_RADIUS_M,
    n_radial: int = 12,
    z_end: float = COLUMN_LENGTH_M,
    steps: int = 400,
) -> list[Tube]:
    """The steady-nozzle idealization: every gram enters at the chamber plane.

    Kept for contrast, because it is what a nozzle calculation would ordinarily assume and it is
    the pessimistic bound on the geometric load.
    """
    radii = [bag_radius * math.sqrt((i + 0.5) / n_radial) for i in range(n_radial)]
    return [trace_tube(stack, r, z_start=0.0, z_end=z_end, steps=steps) for r in radii]


# ---- The three numbers the asks want ----------------------------------------------------------


def mean_cos_theta(tubes: list[Tube], z: float) -> float:
    """Mass-weighted `<cos theta>` at a station -- R1's divergence factor.

    Averaged over the tubes that **reach** the station and have not escaped, since a tube that
    left the bore carries no exhaust. Equal weights, because the launch grid is equal-weight.
    """
    live = [t for t in tubes if not t.escaped and t.z[-1] >= z - 1e-9 and t.z[0] <= z + 1e-9]
    if not live:
        return float("nan")
    return sum(math.cos(t.theta_at(z)) for t in live) / len(live)


def mean_cos_theta_exit(tubes: list[Tube]) -> float:
    """`<cos theta>` at the exit, over the tubes that get there."""
    live = [t for t in tubes if not t.escaped]
    if not live:
        return float("nan")
    return sum(math.cos(t.theta_exit_rad) for t in live) / len(live)


def clearing_fraction(tubes: list[Tube]) -> float:
    """Share of the plume whose whole tube fits inside the winding -- R5's geometric load."""
    return sum(1 for t in tubes if not t.escaped) / len(tubes)


def _winding_radius(stack: field.CoilStack, z: float) -> float:
    """The winding's radius at a station, interpolated between coils."""
    zs = np.array([c.z_m for c in stack.coils])
    rs = np.array([c.radius_m for c in stack.coils])
    order = np.argsort(zs)
    return float(np.interp(z, zs[order], rs[order]))


@dataclass(frozen=True)
class BundleResult:
    """What one winding does to the plume, under one launch model."""

    winding: str
    launch: str
    n_tubes: int
    mean_cos_theta_exit: float
    max_theta_deg: float
    clearing_fraction: float
    missing_fraction: float
    r_exit_bounding_m: float
    winding_r_exit_m: float


def evaluate(
    name: str,
    stack: field.CoilStack,
    *,
    launch: str = "column",
    n_radial: int = 12,
    n_axial: int = 12,
    z_end: float = COLUMN_LENGTH_M,
) -> tuple[BundleResult, list[Tube]]:
    """Trace a bundle through one winding and reduce it to the numbers R1, R2 and R5 want."""
    if launch == "column":
        tubes = trace_bundle(stack, n_radial=n_radial, n_axial=n_axial, z_end=z_end)
    elif launch == "chamber":
        tubes = chamber_bundle(stack, n_radial=n_radial, z_end=z_end)
    else:
        raise ValueError(f"unknown launch model {launch!r}")
    cleared = clearing_fraction(tubes)
    live = [t for t in tubes if not t.escaped]
    return (
        BundleResult(
            winding=name,
            launch=launch,
            n_tubes=len(tubes),
            mean_cos_theta_exit=mean_cos_theta_exit(tubes),
            max_theta_deg=math.degrees(max((abs(t.theta_exit_rad) for t in live), default=0.0)),
            clearing_fraction=cleared,
            missing_fraction=1.0 - cleared,
            r_exit_bounding_m=max((t.r_exit_m for t in live), default=0.0),
            winding_r_exit_m=_winding_radius(stack, z_end),
        ),
        tubes,
    )


CSV_HEADER = (
    "winding",
    "launch",
    "n_tubes",
    "mean_cos_theta_exit",
    "max_theta_deg",
    "clearing_fraction",
    "missing_fraction",
    "r_exit_bounding_m",
    "winding_r_exit_m",
)


def write_results(rows: list[BundleResult], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for r in rows:
            writer.writerow(
                [
                    r.winding,
                    r.launch,
                    r.n_tubes,
                    f"{r.mean_cos_theta_exit:.6f}",
                    f"{r.max_theta_deg:.3f}",
                    f"{r.clearing_fraction:.4f}",
                    f"{r.missing_fraction:.4f}",
                    f"{r.r_exit_bounding_m:.4f}",
                    f"{r.winding_r_exit_m:.4f}",
                ]
            )


def main() -> None:
    """R2 and R5, on the real solve; and R1's `<cos theta>`."""
    parser = argparse.ArgumentParser(description="Flux-tube accounting on the Biot-Savart field")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--radial", type=int, default=12)
    parser.add_argument("--axial", type=int, default=12)
    args = parser.parse_args()

    stacks = build_stacks()
    rows: list[BundleResult] = []
    traced: dict[str, list[Tube]] = {}
    for name, (stack, _profile) in stacks.items():
        for launch in ("column", "chamber"):
            result, tubes = evaluate(
                name, stack, launch=launch, n_radial=args.radial, n_axial=args.axial
            )
            rows.append(result)
            traced[f"{name}/{launch}"] = tubes

    print("== R2: flux tubes on the real Biot-Savart solve, not a paraxial expansion ==")
    print(f"{N_COILS} coils, column {COLUMN_LENGTH_M} m, tubes truncated at the winding contour")
    print("'column' launches over the whole bore (the plume fills it); 'chamber' is the")
    print("steady-nozzle idealization, every gram entering at z = 0.\n")
    print(
        f"{'winding':>16} {'launch':>8} {'<cos theta>':>12} {'max theta':>10} "
        f"{'clears':>8} {'MISSING':>9} {'bounding r':>11} {'winding r':>10}"
    )
    for r in rows:
        print(
            f"{r.winding:>16} {r.launch:>8} {r.mean_cos_theta_exit:12.4f} "
            f"{r.max_theta_deg:9.2f}d {100 * r.clearing_fraction:7.1f}% "
            f"{100 * r.missing_fraction:8.1f}% {r.r_exit_bounding_m:10.2f}m "
            f"{r.winding_r_exit_m:9.2f}m"
        )

    print("\n== The divergence profile down the column (flared-capped, column launch) ==")
    tubes = traced["flared-capped/column"]
    print(f"{'z [m]':>7} {'<cos theta>':>12} {'theta of the widest live tube':>32}")
    for z in (1.0, 3.0, 6.0, 12.0, 18.0, COLUMN_LENGTH_M):
        live = [t for t in tubes if not t.escaped and t.z[0] <= z <= t.z[-1]]
        widest = max(live, key=lambda t: t.radius_at(z)) if live else None
        tail = (
            f"{math.degrees(widest.theta_at(z)):9.2f}d at r = {widest.radius_at(z):.2f} m"
            if widest
            else "-"
        )
        print(f"{z:7.1f} {mean_cos_theta(tubes, z):12.4f} {tail:>32}")

    print("\n== R5: what misses the winding, per launch model ==")
    for r in rows:
        kg = 213.0 * r.missing_fraction
        print(
            f"  {r.winding:>16} / {r.launch:<8}: {100 * r.missing_fraction:5.1f}%, "
            f"{kg:6.1f} kg/pulse against the 4.9 kg booked ({kg / 4.9:5.1f}x)"
        )

    write_results(rows, args.output)
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()


# ---- Downstream of the last coil: where R1's divergence loss actually happens ------------------


@dataclass(frozen=True)
class FanStation:
    """The plume's divergence at one station past the exit plane."""

    z_m: float
    radii_out: float
    """Distance past the exit plane, in exit radii -- the units P3 states detachment in."""
    mean_cos_theta: float
    max_theta_deg: float
    live_fraction: float
    """Share of tubes still running forward here. The rest have turned; see `downstream_fan`."""


def downstream_fan(
    stack: field.CoilStack,
    *,
    z_exit: float = COLUMN_LENGTH_M,
    r_exit: float,
    n_radial: int = 10,
    out_to_radii: float = 2.5,
    steps: int = 1200,
) -> list[FanStation]:
    """Trace the plume's tubes past the last coil and report how far they have fanned.

    **This is where R1's divergence term lives.** Inside the winding the bore is nearly axial and
    `<cos theta>` is within half a percent of 1; past the last coil the vacuum field opens toward
    its own return path and the tubes fan hard. P3 puts detachment at 1.44-2.00 exit radii out,
    so that window is where `<cos theta>` has to be read.

    A tube is dropped once its field line **turns** -- `B_z <= 0`, or the traced radius stops
    increasing monotonically, which is the tracer following the line back toward the return flux.
    A plume cannot follow a line back, so a turned tube is not a divergence measurement; it is a
    statement that detachment must have happened upstream of there. The surviving share is
    reported at every station so the reader can see the estimate thinning.
    """
    launch = [r_exit * math.sqrt((i + 0.5) / n_radial) for i in range(n_radial)]
    z_end = z_exit + out_to_radii * r_exit

    traces: list[list[tuple[float, float, float]]] = []
    for r0 in launch:
        try:
            path = field.trace_line(stack, r0, z_exit, z_end=z_end, steps=steps)
        except ValueError:
            traces.append([])
            continue
        kept: list[tuple[float, float, float]] = []
        for z, r, _b in path:
            b_r, b_z = stack.field(max(r, 1e-9), z)
            if b_z <= 0.0:
                break
            if kept and r < kept[-1][1]:
                break
            kept.append((z, r, math.atan2(b_r, b_z)))
        traces.append(kept)

    out: list[FanStation] = []
    n_stations = 26
    for j in range(n_stations):
        radii_out = out_to_radii * j / (n_stations - 1)
        z = z_exit + radii_out * r_exit
        thetas = [
            float(np.interp(z, [p[0] for p in t], [p[2] for p in t]))
            for t in traces
            if t and t[0][0] <= z <= t[-1][0]
        ]
        if not thetas:
            out.append(FanStation(z, radii_out, float("nan"), float("nan"), 0.0))
            continue
        out.append(
            FanStation(
                z_m=z,
                radii_out=radii_out,
                mean_cos_theta=sum(math.cos(t) for t in thetas) / len(thetas),
                max_theta_deg=math.degrees(max(abs(t) for t in thetas)),
                live_fraction=len(thetas) / len(traces),
            )
        )
    return out
