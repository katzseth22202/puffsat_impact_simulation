"""Can a field ripple hold plasma in the bore? (R10, and the correct form of N5)

R10 accepts P8's ripple sweep and disputes its **pass criterion**: P8 reported whether an off-axis
`|B|` minimum *exists* and concluded ">= 36 coils", while a magnetic trap only holds what its
mirror ratio can hold, so the test should be `R > 1/sin^2(theta)`. R10 is right that existence is
the wrong test. **It is right for the wrong reason, and the correct reason retires the criterion
rather than sharpening it.**

# The loss cone is a collisionless object

Mirror trapping is `mu` conservation seen from the other side: a particle turns around where
`v_perp^2 = v^2` because its own `mu` is fixed, and the loss cone `sin^2 theta < 1/R` is the set of
pitch angles for which that never happens. `continuum.py` shows the plume collides millions of
times crossing the bore, so no particle keeps its `mu` and there is no pitch angle to test.

The fluid statement is cleaner still. `(J x B) . B = 0` identically, so **the Lorentz force has no
component along a field line**. In a collisional plasma the pressure is a scalar, the parallel
momentum equation is `rho Du/Dt = -dp/ds` with no magnetic term, and a `|B|` variation exerts no
force along the flow at all. The mirror force that appears in kinetic theory,
`-(p_perp - p_par) d ln B/ds`, is proportional to the anisotropy, and collisions destroy the
anisotropy on the mean-free-path timescale -- about 1e-10 s here, against a 2 ms transit.

**So there is no trap, at any coil count, and P8's constraint is withdrawn rather than relaxed.**
R10's own threshold (`1/(1 - alpha)` = 1.096, from P1's `alpha` = 0.088) inherits the same problem
twice over: it is a loss-cone threshold, and it is built from the anisotropy of a *free* expansion
with no nozzle rather than from anything in the bore.

# What can still stop the flow, and it is a real constraint

A flux tube's cross-section is `A ~ 1/|B|`, so a local `|B|` **maximum** is a local area
**minimum** -- a throat. Supersonic flow through a contraction decelerates, and if the contraction
is deep enough the flow **chokes**: it reaches `M` = 1 and a shock stands there. That is a genuine
fluid mechanism, it is what a ripple can actually do, and it has a sharp criterion. Isentropic flow
at Mach `M` survives any contraction up to

    A/A*(M) = (1/M) [ (2/(gamma+1)) (1 + (gamma-1)/2 M^2) ]^{(gamma+1)/(2(gamma-1))}

and the ripple's contraction ratio is exactly the mirror ratio `B_max/B_min` that `field.py`
already reports. So the test becomes **`R < A/A*(M_local)`**, and the same sweep answers it with a
different comparison.

**The binding station is the chamber, not the wall.** `A/A*` goes to 1 at the sonic point, so the
margin *vanishes* at the throat and grows fast downstream -- 1.03 at `M` = 1.2, but 4.0 at
`M` = 2.7 and 9.7 at `M` = 3.4. That inverts P8's picture, where the danger was near the winding
at the cold end. And it lands exactly on R10's own caveat: ADR-0012's cap puts a **flat 12 T
shelf** at the chamber, with no background gradient to swamp ripple, in the one place where the
flow has no margin to spare.

**Only ripples downstream of the sonic point are tested.** Upstream the flow is subsonic, where a
contraction *accelerates* it; the worst a subsonic contraction can do is move the sonic point,
not block the column. It is supersonic flow that a contraction decelerates and can shock, so each
minimum is judged at the Mach number of **its own station**, read off the solved cooling history,
rather than at one nominal value for the whole winding.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from puffsat import field

DEFAULT_OUTPUT = Path("data/results/nozzle_residence.csv")

GAMMA_BRACKET: tuple[float, ...] = (1.15, 1.25, 1.667)
"""Effective `gamma` bracket. 1.15 is the equilibrium branch's buffered value, 5/3 the frozen
monatomic limit; the choking margin is weakly sensitive across it, which `main` shows.
"""

ALPHA_FREE_EXPANSION = 0.088
"""P1's axial share, which R10 uses to set its threshold. Recorded to say why it does not apply."""


def area_mach(mach: float, gamma: float) -> float:
    """`A/A*` at a Mach number -- the deepest contraction isentropic flow survives without choking.

    The standard isentropic area-Mach relation. At `M` = 1 it is 1 (no margin at all); it rises
    steeply on both sides.
    """
    if mach <= 0.0:
        raise ValueError("Mach number must be positive")
    term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach)
    return float(term ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))) / mach)


def loss_cone_threshold(alpha: float) -> float:
    """`1/(1 - alpha)` -- R10's threshold, computed so the answer can say what it would have been.

    Not used as a criterion here. See the module docstring.
    """
    if alpha >= 1.0:
        return math.inf
    return 1.0 / (1.0 - alpha)


@dataclass(frozen=True)
class Minimum:
    """One local `|B|` minimum: where it is, and how deep."""

    z_m: float
    mirror_ratio: float


def local_minima(
    stack: field.CoilStack, radius: float, *, length: float, n: int = 600
) -> list[Minimum]:
    """Every interior local minimum of `|B|` along a line of constant radius, with its depth.

    `field.scan_axial` reports only the deepest, which was enough for an existence test and is
    not enough for a criterion that varies along the column.
    """
    import numpy as np

    zs = [float(z) for z in np.linspace(0.0, length, n)]
    mags = [stack.magnitude(radius, z) for z in zs]
    out: list[Minimum] = []
    for i in range(1, len(mags) - 1):
        if mags[i] < mags[i - 1] and mags[i] < mags[i + 1]:
            left = max(mags[:i])
            right = max(mags[i + 1 :])
            out.append(Minimum(z_m=zs[i], mirror_ratio=min(left, right) / mags[i]))
    return out


def mach_profile(path: Path = Path("data/results/cooling_history.csv")) -> Callable[[float], float]:
    """`M(z)` from the solved cooling history -- the margin's own argument.

    Takes the **lowest** Mach number across the leg/branch cases at each station, so the margin
    reported is the tightest of the flown cases rather than an average of them.
    """
    import numpy as np

    rows = list(csv.DictReader(path.open()))
    by_x: dict[float, float] = {}
    for r in rows:
        x = round(float(r["x_m"]), 3)
        m = float(r["mach"])
        by_x[x] = min(by_x.get(x, math.inf), m)
    xs = sorted(by_x)
    ms = [by_x[x] for x in xs]

    def interpolate(z: float) -> float:
        return float(np.interp(z, xs, ms))

    return interpolate


@dataclass(frozen=True)
class RippleVerdict:
    """One ripple measurement, judged by choking rather than by trapping."""

    winding: str
    n_coils: int
    radius_m: float
    z_m: float
    """Station of the binding minimum -- where the margin has to be read."""
    mach_here: float
    mirror_ratio: float
    n_local_minima: int
    choking_margin: float
    """`A/A*(M)` at the station the ripple sits at -- the contraction the flow can survive."""
    chokes: bool
    would_trap_under_r10: bool
    """Whether R10's loss-cone criterion would have flagged it. Recorded, not used."""


def sweep(
    *,
    coil_counts: tuple[int, ...] = (12, 18, 24, 36, 48, 72),
    radii: tuple[float, ...] = (0.0, 1.0, 2.0, 2.5),
    gamma: float = 1.25,
    cap_t: float = field.ADR0012_CAP_T,
    length: float = 23.8,
    windings: tuple[str, ...] = ("straight-flown", "flared-capped"),
    mach_at: Callable[[float], float] | None = None,
) -> list[RippleVerdict]:
    """Re-run P8's sweep on the current winding, judged by choking at each ripple's own station.

    The binding minimum for a winding is the one with the **smallest margin**, `A/A*(M(z))/R`,
    not the deepest one: a deep ripple far downstream is harmless and a shallow one near the
    sonic point is not.
    """
    # Imported here rather than at module scope: `fluxtube` imports `field`, and this module is
    # about the criterion rather than about the geometry, so the dependency is one-way on demand.
    from puffsat import fluxtube

    mach = mach_at or mach_profile()
    stacks = fluxtube.build_stacks(cap_t=cap_t, length=length)
    threshold = loss_cone_threshold(ALPHA_FREE_EXPANSION)

    out: list[RippleVerdict] = []
    for name in windings:
        for n in coil_counts:
            _stack, profile = stacks[name]
            if name == "flared-capped":

                def flare(z: float, profile: field.Profile = profile) -> float:
                    return fluxtube.flared_radius(min(max(z, 0.0), length), profile)

                stack = field.build_winding(n, length=length, profile=profile, radius_at=flare)
            else:
                stack = field.build_winding(
                    n,
                    coil_radius=fluxtube.STRAIGHT_COIL_RADIUS_M,
                    length=length,
                    profile=profile,
                )
            for r in radii:
                minima = local_minima(stack, r, length=length)
                supersonic = [m for m in minima if mach(m.z_m) > 1.0]
                if not supersonic:
                    out.append(
                        RippleVerdict(
                            winding=name,
                            n_coils=n,
                            radius_m=r,
                            z_m=float("nan"),
                            mach_here=float("nan"),
                            mirror_ratio=1.0,
                            n_local_minima=len(minima),
                            choking_margin=float("inf"),
                            chokes=False,
                            would_trap_under_r10=False,
                        )
                    )
                    continue

                def margin_of(m: Minimum) -> float:
                    return area_mach(mach(m.z_m), gamma) / m.mirror_ratio

                binding = min(supersonic, key=margin_of)
                local_mach = mach(binding.z_m)
                out.append(
                    RippleVerdict(
                        winding=name,
                        n_coils=n,
                        radius_m=r,
                        z_m=binding.z_m,
                        mach_here=local_mach,
                        mirror_ratio=binding.mirror_ratio,
                        n_local_minima=len(minima),
                        choking_margin=area_mach(local_mach, gamma),
                        chokes=binding.mirror_ratio >= area_mach(local_mach, gamma),
                        would_trap_under_r10=max(m.mirror_ratio for m in minima) >= threshold,
                    )
                )
    return out


CSV_HEADER = (
    "winding",
    "n_coils",
    "radius_m",
    "z_m",
    "mach_here",
    "mirror_ratio",
    "n_local_minima",
    "choking_margin",
    "chokes",
    "would_trap_under_r10",
)


def write_verdicts(rows: list[RippleVerdict], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for r in rows:
            writer.writerow(
                [
                    r.winding,
                    r.n_coils,
                    f"{r.radius_m:g}",
                    f"{r.z_m:.4f}",
                    f"{r.mach_here:.4f}",
                    f"{r.mirror_ratio:.6f}",
                    r.n_local_minima,
                    f"{r.choking_margin:.4f}",
                    int(r.chokes),
                    int(r.would_trap_under_r10),
                ]
            )


def main() -> None:
    """R10: the depth criterion, corrected for the regime that actually applies."""
    parser = argparse.ArgumentParser(description="Ripple, residence and choking (R10)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    print("== The margin a supersonic flow has against a contraction ==")
    print("A ripple's contraction ratio is its mirror ratio. The flow chokes if the ratio")
    print("exceeds A/A*(M). The margin vanishes at the sonic point and grows fast after it.\n")
    print(f"{'M':>6} " + " ".join(f"{'g=' + format(g, '.3g'):>10}" for g in GAMMA_BRACKET))
    for m in (1.01, 1.05, 1.1, 1.2, 1.5, 2.0, 2.7, 3.4):
        print(f"{m:6.2f} " + " ".join(f"{area_mach(m, g):10.3f}" for g in GAMMA_BRACKET))

    rows = sweep()
    print("\n== P8's sweep, judged by choking at each ripple's OWN station ==")
    print("`M(z)` is the tightest of the flown legs at that station, from the cooling history.")
    print("A blank row is a winding whose only minima sit in subsonic flow, where a")
    print("contraction cannot choke anything.\n")
    print(
        f"{'winding':>16} {'coils':>6} {'r [m]':>6} {'minima':>7} {'binding z':>10} "
        f"{'M there':>8} {'mirror R':>9} {'margin':>8} {'chokes?':>8} {'R10?':>6}"
    )
    for r in rows:
        if r.n_local_minima == 0:
            continue
        z = "-" if math.isnan(r.z_m) else f"{r.z_m:9.2f}m"
        m_here = "-" if math.isnan(r.mach_here) else f"{r.mach_here:8.3f}"
        margin = "inf" if math.isinf(r.choking_margin) else f"{r.choking_margin:8.3f}"
        print(
            f"{r.winding:>16} {r.n_coils:6d} {r.radius_m:6.1f} {r.n_local_minima:7d} "
            f"{z:>10} {m_here:>8} {r.mirror_ratio:9.4f} {margin:>8} "
            f"{('YES' if r.chokes else 'no'):>8} "
            f"{('yes' if r.would_trap_under_r10 else 'no'):>6}"
        )

    print("\n== The coil count each criterion demands ==")
    for winding in ("straight-flown", "flared-capped"):
        group = [r for r in rows if r.winding == winding]
        counts = sorted({r.n_coils for r in group})
        choke_ok = next(
            (n for n in counts if not any(r.chokes for r in group if r.n_coils == n)),
            None,
        )
        p8_ok = next(
            (n for n in counts if not any(r.n_local_minima > 0 for r in group if r.n_coils == n)),
            None,
        )
        r10_ok = next(
            (n for n in counts if not any(r.would_trap_under_r10 for r in group if r.n_coils == n)),
            None,
        )
        print(
            f"  {winding:>16}: choking (correct) >= {choke_ok}, "
            f"R10's loss cone >= {r10_ok}, P8's existence >= {p8_ok}"
        )

    print(
        "\nR10's loss-cone threshold would have been "
        f"{loss_cone_threshold(ALPHA_FREE_EXPANSION):.3f}"
        " -- it does not apply; see the module docstring."
    )

    write_verdicts(rows, args.output)
    print(f"wrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
