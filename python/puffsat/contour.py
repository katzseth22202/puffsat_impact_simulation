"""The continuous survivable-density contour the cloud schedule flies (Q15/Q13, design SS7/SS12.1).

The discrete construction evaluates 27 fixed cloud shapes per velocity and keeps the best survivor.
That works, but it lets the *shape grid* decide where the answer steps rather than the physics: at
45 -> 46 km/s the optimum jumps `r_foot/R` 0.5 -> 0.7 and `f` drops 0.817 -> 0.791, purely because
no intermediate footprint exists to be chosen.

Q15 replaces it with the contour itself:

```
rho(v) = min(rho_ceiling(v), rho_max),   rho_ceiling = P_limit / (c_stag v^2)
```

solved back through the Sigma contract (ADR-0003) for the cloud shape that lands on it. Which of
the two limits binds is reported, not smoothed over -- "as dense as the box can build" and "as
dense as the plate can take" are different engineering situations with different levers.

Fixing `rho` is one equation in two shape parameters, so a *one-parameter family* of shapes hits any
given contour density. The point is therefore **chosen**: `e_eff` is already pinned by `rho`, so the
choice maximizes `eta_capture` along that family. In practice this always lands at the shortest
feasible `L/D` -- `eta` falls monotonically with cloud length -- which is why the discrete
construction kept picking `L/D = 0.3`.

**The `eta_capture` interpolation is only claimed along this contour.** Across the full shape box
`eta` spans 0.79-0.99 (see `ETA_BOX_RANGE`); a bilinear fit over that whole range would be
advertising a precision the 3x3x3 sweep does not support. Along the schedule the contour flies it
stays inside a ~0.05 band, which is the only regime the fit is used in and the only one it is
offered for.

Pure post-processing: no new physics, no new sweep. It reads the geometry sweep's `eta_capture` and
the heavy-plate sweep's `e_eff(rho)` and reduces them onto one curve.
"""

from __future__ import annotations

import itertools
import json
import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from puffsat.analysis import impact_density

# The heavy-plate scenario's pulse mass and plate radius (design SS12.1). Defined here rather than
# imported from `heavyplate` to keep the dependency one-way -- `heavyplate` consumes this module --
# and deliberately NOT taken from `analysis`, whose 25 kg / 5 m belong to the core envelope study.
# Reading the wrong pair silently moves every contour density by a factor of ~7, which is why a test
# pins these against `heavyplate`'s rather than trusting the comment.
PULSE_MASS_KG = 100.0
PLATE_RADIUS_M = 15.0

# The assumed cloud shape box (design SS13). The contour is only solved inside it: a shape outside
# is not a shape this study claims can be delivered.
L_OVER_D_BOX = (0.3, 1.0)
R_FOOT_BOX = (0.3, 0.7)

# Plate curvature the headline contour flies. Shallow concave is the recovery lever for the rebound
# axiality the survivability-driven cloud stretch forces you to give up (design SS7); the deep dish
# is foreclosed by ADR-0021, and flat is reported alongside as the conservative floor.
D_OVER_D_HEADLINE = 0.10
D_OVER_D_FLAT = 0.0

DEFAULT_GEOMETRY_PATH = Path("data/results/sweep_geometry_m40.jsonl")

# Full-shape-box eta range, for the caveat above. Filled from the geometry sweep on first use.
ETA_BOX_RANGE = (0.7924, 0.9920)


def rho_ceiling(v: float, c_stag: float, p_limit: float) -> float:
    """Densest cloud the facesheet survives at closing speed `v`: `rho = P_limit/(c_stag v^2)`.

    The stagnation law of design SS7 / ADR-0010, inverted. Inverse-square in `v`, which is why the
    contour falls steeply across the extended range and the cloud must be stretched to follow it.
    """
    return p_limit / (c_stag * v * v)


def l_over_d_for(
    rho: float,
    r_foot_over_r: float,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
) -> float:
    """Cloud aspect ratio that puts the Sigma contract at `rho` for a given footprint.

    The exact inverse of `analysis.impact_density`: from `rho = m/(2 pi (L/D) r_foot^3)`,
    `L/D = m/(2 pi rho r_foot^3)`.
    """
    r_foot = r_foot_over_r * plate_radius
    return mass / (2.0 * math.pi * rho * r_foot * r_foot * r_foot)


def _geometry_table(path: Path = DEFAULT_GEOMETRY_PATH) -> list[dict[str, float]]:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    return [{k: float(v) for k, v in r.items() if isinstance(v, (int, float))} for r in rows]


_GEO_CACHE: list[dict[str, float]] | None = None


def _geo(path: Path = DEFAULT_GEOMETRY_PATH) -> list[dict[str, float]]:
    global _GEO_CACHE
    if _GEO_CACHE is None:
        _GEO_CACHE = _geometry_table(path)
    return _GEO_CACHE


def eta_at(
    l_over_d: float,
    r_foot_over_r: float,
    d_over_d: float,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> float:
    """`eta_capture` at an arbitrary shape, bilinear over the geometry sweep's `(L/D, r_foot/R)`
    nodes at fixed curvature.

    Bilinear rather than anything cleverer because the sweep is 3x3 per curvature: a higher-order
    fit would be inventing structure between four points. Clamped at the box edges -- the contour is
    never solved outside the box, so a clamp can only be reached by round-off.

    **Only claimed along the contour** (see the module docstring): across the full box `eta` spans
    `ETA_BOX_RANGE`, and this fit is not offered as a general surrogate for the 2D track.
    """
    rows = [r for r in _geo(path) if abs(r["d_over_d"] - d_over_d) < 1e-9]
    if not rows:
        raise KeyError(f"no geometry rows at d/D = {d_over_d}")
    xs = sorted({r["l_over_d"] for r in rows})
    ys = sorted({r["r_foot_over_r"] for r in rows})
    node = {(r["l_over_d"], r["r_foot_over_r"]): r["eta_capture"] for r in rows}

    def bracket(vals: list[float], t: float) -> tuple[float, float, float]:
        t = min(max(t, vals[0]), vals[-1])
        for lo, hi in itertools.pairwise(vals):
            if lo <= t <= hi:
                return lo, hi, 0.0 if hi == lo else (t - lo) / (hi - lo)
        return vals[-2], vals[-1], 1.0

    x0, x1, tx = bracket(xs, l_over_d)
    y0, y1, ty = bracket(ys, r_foot_over_r)
    a = node[(x0, y0)] * (1 - ty) + node[(x0, y1)] * ty
    b = node[(x1, y0)] * (1 - ty) + node[(x1, y1)] * ty
    return a * (1 - tx) + b * tx


@dataclass(frozen=True)
class ContourPoint:
    """The cloud the schedule flies at one closing speed."""

    v: float
    rho_ceiling: float
    rho_contour: float
    #: True when survivability binds; False when the shape box does. Different levers.
    ceiling_limited: bool
    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    eta_capture: float
    #: `e_eff` at the contour density, or None when no restitution source was supplied.
    e_eff: float | None = None
    #: `f = eta*(1 + e_eff)/2` (ADR-0001), or None likewise.
    f: float | None = None


def rho_max_achievable(mass: float = PULSE_MASS_KG, plate_radius: float = PLATE_RADIUS_M) -> float:
    """Densest cloud the shape box can deliver: shortest `L/D` at the tightest footprint."""
    return impact_density(L_OVER_D_BOX[0], R_FOOT_BOX[0], mass, plate_radius)


def contour_point(
    v: float,
    c_stag: float,
    p_limit: float,
    d_over_d: float = D_OVER_D_HEADLINE,
    samples: int = 512,
    path: Path = DEFAULT_GEOMETRY_PATH,
    e_eff_at: Callable[[float], float] | None = None,
) -> ContourPoint:
    """The contour point at `v`: the highest-`eta` shape on the survivable iso-density curve.

    Sampled rather than solved analytically. The objective is `eta` bilinear over a 3x3 grid
    restricted to a curve, so it is piecewise-smooth with kinks at the grid lines -- a gradient
    method would stop at one of those kinks, and 512 samples across a one-parameter family costs
    nothing in post-processing.

    `e_eff_at` is injected rather than imported so this module stays free of the heavy-plate
    analysis that consumes it; pass `heavyplate.e_eff_interpolator_at_v(rows, v)` to get `f`.
    """
    ceiling = rho_ceiling(v, c_stag, p_limit)
    box = rho_max_achievable()
    rho = min(ceiling, box)

    best: tuple[float, float, float] | None = None  # (eta, l_over_d, r_foot_over_r)
    for i in range(samples):
        rf = R_FOOT_BOX[0] + i * (R_FOOT_BOX[1] - R_FOOT_BOX[0]) / (samples - 1)
        lod = l_over_d_for(rho, rf)
        if not (L_OVER_D_BOX[0] <= lod <= L_OVER_D_BOX[1]):
            continue
        eta = eta_at(lod, rf, d_over_d, path)
        if best is None or eta > best[0]:
            best = (eta, lod, rf)
    if best is None:
        raise ValueError(
            f"the iso-density curve at rho={rho:.4g} (v={v:.0f} m/s) misses the shape box"
        )
    e_eff = e_eff_at(rho) if e_eff_at is not None else None
    return ContourPoint(
        v=v,
        rho_ceiling=ceiling,
        rho_contour=rho,
        ceiling_limited=ceiling <= box,
        d_over_d=d_over_d,
        l_over_d=best[1],
        r_foot_over_r=best[2],
        eta_capture=best[0],
        e_eff=e_eff,
        f=None if e_eff is None else best[0] * (1.0 + e_eff) / 2.0,
    )
