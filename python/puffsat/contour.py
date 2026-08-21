"""The continuous survivable-density contour the cloud schedule flies (Q15/Q13, design SS7/SS12.1).

The discrete construction evaluates 27 fixed cloud shapes per velocity and keeps the best survivor.
That works, but it lets the *shape grid* decide where the answer steps rather than the physics: at
45 -> 46 km/s the optimum jumps `r_foot/R` 0.5 -> 0.7 and `f` drops 0.817 -> 0.791, purely because
no intermediate footprint exists to be chosen.

Q15 replaces it with a **constrained maximum over the shape box**: the Sigma contract (ADR-0003)
already fixes `rho` from the shape, so sweep the box, drop the shapes the facesheet cannot survive,
and keep the best of what is left.

```
survives(shape)  <=>  c_stag * rho(shape) * v^2 * focusing(shape)  <=  P_limit
```

Which limit binds is reported, not smoothed over -- "as dense as the box can build" and "as dense
as the plate can take" are different engineering situations with different levers.

**Corrected 2026-08-21.** This module originally set `rho = min(rho_ceiling, rho_max)` first and
then maximized `eta_capture` along that iso-density curve. That is valid only if the survivable
density is a property of `v` alone -- and it is not, because the facesheet sees the plane-wave
stagnation load *concentrated by the dish* (`focusing_at`, ADR-0010), which runs 1.15-2.20 across
the box. The old construction therefore flew over its own pressure limit wherever survivability
bound: 470-565 MPa against a 400 MPa baseline across 28-60 km/s.

`heavyplate.heavyplate_frontier` and `analysis.survivability_frontier` had always applied the
factor to the discrete grid, so the two constructions disagreed; the discrete one was right. The
corrected contour reproduces its answer at 45 km/s (`rho = 0.1251` against the grid's 0.1258, both
at 399 MPa), which is what identified the contour as the faulty side.

Stating the search over the 2D box rather than a 1D iso-density family also removes a circularity
(`rho` -> shape -> `focusing` -> `rho`) that would otherwise need a fixed point, and lets the
optimizer trade `eta_capture` against `focusing` -- a flatter-`eta` shape that concentrates less can
fly denser and win on `e_eff`.

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
from functools import cache
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


@cache
def _nodes(
    field: str, d_over_d: float, path: Path
) -> tuple[tuple[float, ...], tuple[float, ...], dict[tuple[float, float], float]]:
    """The `(L/D, r_foot/R)` node table for one field at one curvature.

    Cached: the contour search is 2D over the shape box, so this is hit ~10^6 times per curve and
    rebuilding the dict each time dominated the runtime.
    """
    rows = [r for r in _geo(path) if abs(r["d_over_d"] - d_over_d) < 1e-9]
    if not rows:
        raise KeyError(f"no geometry rows at d/D = {d_over_d}")
    return (
        tuple(sorted({r["l_over_d"] for r in rows})),
        tuple(sorted({r["r_foot_over_r"] for r in rows})),
        {(r["l_over_d"], r["r_foot_over_r"]): r[field] for r in rows},
    )


def _bilinear(
    field: str,
    l_over_d: float,
    r_foot_over_r: float,
    d_over_d: float,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> float:
    """One geometry-sweep `field` at an arbitrary shape, bilinear over the `(L/D, r_foot/R)` nodes
    at fixed curvature.

    Bilinear rather than anything cleverer because the sweep is 3x3 per curvature: a higher-order
    fit would be inventing structure between four points. Clamped at the box edges -- the contour is
    never solved outside the box, so a clamp can only be reached by round-off.
    """
    xs, ys, node = _nodes(field, d_over_d, path)

    def bracket(vals: tuple[float, ...], t: float) -> tuple[float, float, float]:
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


def eta_at(
    l_over_d: float,
    r_foot_over_r: float,
    d_over_d: float,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> float:
    """`eta_capture` at an arbitrary shape (see `_bilinear` for the scheme).

    **Only claimed along the contour** (see the module docstring): across the full box `eta` spans
    `ETA_BOX_RANGE`, and this fit is not offered as a general surrogate for the 2D track.
    """
    return _bilinear("eta_capture", l_over_d, r_foot_over_r, d_over_d, path)


def focusing_at(
    l_over_d: float,
    r_foot_over_r: float,
    d_over_d: float,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> float:
    """Concave local-peak concentration over the flat reference at the same shape (ADR-0010).

    `analysis.survivability_frontier` and `heavyplate.heavyplate_frontier` both classify the
    facesheet against `c_stag rho v^2 * focusing`, where `focusing = peak_local_pressure(d/D) /
    peak_local_pressure(flat)`. A dish concentrates the stagnation load into a smaller patch, so the
    *local* peak the facesheet must survive is above the plane-wave value even though the total
    impulse is unchanged.

    It runs **1.15-2.20** across the shape box and climbs steeply with footprint, so leaving it out
    is not a small approximation: the pre-2026-08-21 contour did, and flew 470-565 MPa against the
    400 MPa limit it claimed to respect. Flat is 1.0 by construction.
    """
    if abs(d_over_d) < 1e-12:
        return 1.0
    concave = _bilinear("peak_local_pressure", l_over_d, r_foot_over_r, d_over_d, path)
    flat = _bilinear("peak_local_pressure", l_over_d, r_foot_over_r, 0.0, path)
    return concave / flat


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
    #: Concave local-peak concentration at this shape (ADR-0010); 1.0 for a flat plate.
    focusing: float = 1.0
    #: Peak facesheet pressure actually flown [Pa], `c_stag rho v^2 * focusing`, <= `p_limit`.
    peak_pressure: float = 0.0
    #: `e_eff` at the contour density, or None when no restitution source was supplied.
    e_eff: float | None = None
    #: `f = eta*(1 + e_eff)/2` (ADR-0001), or None likewise.
    f: float | None = None


def rho_max_achievable(mass: float = PULSE_MASS_KG, plate_radius: float = PLATE_RADIUS_M) -> float:
    """Densest cloud the shape box can deliver: shortest `L/D` at the tightest footprint."""
    return impact_density(L_OVER_D_BOX[0], R_FOOT_BOX[0], mass, plate_radius)


def survivable(
    l_over_d: float,
    r_foot_over_r: float,
    v: float,
    c_stag: float,
    p_limit: float,
    d_over_d: float,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> tuple[bool, float, float, float]:
    """Does the facesheet survive this shape? Returns `(ok, rho, focusing, peak)`.

    The Sigma contract fixes `rho` from the shape, the stagnation law turns that into a plane-wave
    peak, and `focusing_at` concentrates it (ADR-0010). Identical to the test
    `heavyplate.heavyplate_frontier` applies to the discrete grid -- that agreement is the point.
    """
    rho = impact_density(l_over_d, r_foot_over_r, mass, plate_radius)
    focusing = focusing_at(l_over_d, r_foot_over_r, d_over_d, path)
    peak = c_stag * rho * v * v * focusing
    return peak <= p_limit, rho, focusing, peak


def contour_point(
    v: float,
    c_stag: float,
    p_limit: float,
    d_over_d: float = D_OVER_D_HEADLINE,
    samples: int = 192,
    path: Path = DEFAULT_GEOMETRY_PATH,
    e_eff_at: Callable[[float], float] | None = None,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
) -> ContourPoint:
    """The best cloud the plate survives at `v`: a constrained maximum over the shape box.

    **Reformulated 2026-08-21 (the focusing fix).** The original construction picked
    `rho = min(rho_ceiling, rho_max)` first and then maximized `eta` along that iso-density curve.
    That is only valid when the survivable density is a property of `v` alone. It is not: the
    facesheet sees `c_stag rho v^2 * focusing(shape)`, and `focusing` runs 1.15-2.20 across the box
    (`focusing_at`). Fixing `rho` up front therefore chose a density the plate could not take.

    Making `rho` shape-dependent makes the old formulation circular (`rho` -> shape -> `focusing`
    -> `rho`). Rather than iterate that to a fixed point, the search is stated the way it should
    always have been: **the Sigma contract already determines `rho` from the shape**, so sweep the
    2D shape box directly, drop the shapes the facesheet cannot survive, and keep the best of what
    is left. No circularity, and `eta` versus `focusing` is now a trade the optimizer can actually
    make -- a slightly flatter-`eta` shape that concentrates less can fly denser and win on `e_eff`.

    **Objective.** With `e_eff_at` supplied it maximizes `f` directly, which is the deliverable.
    Without one it maximizes `rho` -- the monotone surrogate, since a denser cloud is optically
    thicker and bounces better -- and reports `f = None` rather than inventing a restitution.

    `e_eff_at` is injected rather than imported so this module stays free of the heavy-plate
    analysis that consumes it; pass `heavyplate.e_eff_interpolator_at_v(rows, v)` to get `f`.
    """

    def score(rho: float, eta: float) -> float:
        if e_eff_at is None:
            return rho
        return eta * (1.0 + e_eff_at(rho)) / 2.0

    def smallest_surviving_l_over_d(rf: float) -> float | None:
        """The best `L/D` at this footprint, or None if nothing at it survives.

        Both factors of the objective fall as the cloud lengthens -- `eta_capture` drops (a longer
        column splats with more radial relief) and `rho` drops as `1/(L/D)` through the Sigma
        contract, taking `e_eff` with it. So the optimum at fixed `rf` is always the *shortest*
        cloud the facesheet survives, and `peak` falls monotonically with `L/D`, which makes that a
        bisection rather than a scan. `test_peak_and_score_fall_with_cloud_length` pins both
        monotonicities, since the bisection is only valid while they hold.
        """
        lo, hi = L_OVER_D_BOX
        if survivable(lo, rf, v, c_stag, p_limit, d_over_d, mass, plate_radius, path)[0]:
            return lo  # even the shortest cloud survives: the box binds here, not the plate
        if not survivable(hi, rf, v, c_stag, p_limit, d_over_d, mass, plate_radius, path)[0]:
            return None  # nothing at this footprint survives
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if survivable(mid, rf, v, c_stag, p_limit, d_over_d, mass, plate_radius, path)[0]:
                hi = mid
            else:
                lo = mid
        return hi

    # The unconstrained optimum over the box, for reporting which limit binds. Because the score
    # falls monotonically with `L/D`, the best shape ignoring pressure always sits on the box's own
    # `L/D` floor -- so this is one scan over `rf`, not a second search.
    best_free = max(
        score(
            impact_density(L_OVER_D_BOX[0], rf, mass, plate_radius),
            eta_at(L_OVER_D_BOX[0], rf, d_over_d, path),
        )
        for rf in (
            R_FOOT_BOX[0] + i * (R_FOOT_BOX[1] - R_FOOT_BOX[0]) / (samples - 1)
            for i in range(samples)
        )
    )

    best: tuple[float, float, float, float, float, float] | None = None
    for i in range(samples):
        rf = R_FOOT_BOX[0] + i * (R_FOOT_BOX[1] - R_FOOT_BOX[0]) / (samples - 1)
        lod = smallest_surviving_l_over_d(rf)
        if lod is None:
            continue
        _ok, rho, _focusing, peak = survivable(
            lod, rf, v, c_stag, p_limit, d_over_d, mass, plate_radius, path
        )
        eta = eta_at(lod, rf, d_over_d, path)
        sc = score(rho, eta)
        if best is None or sc > best[0]:
            best = (sc, lod, rf, eta, rho, peak)
    if best is None:
        raise ValueError(
            f"no shape in the box survives {p_limit / 1e6:.0f} MPa at v={v:.0f} m/s "
            f"(d/D={d_over_d}) -- the plate cannot fly this velocity"
        )
    _, l_over_d, r_foot_over_r, eta, rho, peak = best
    e_eff = e_eff_at(rho) if e_eff_at is not None else None
    return ContourPoint(
        v=v,
        rho_ceiling=rho_ceiling(v, c_stag, p_limit),
        rho_contour=rho,
        # Survivability binds when the pressure constraint is *active* at the optimum -- the cloud
        # is as dense as the plate can take. Otherwise the shape box is what limits it: the box
        # cannot build anything denser, and the facesheet has margin to spare.
        #
        # Survivability binds when the pressure limit actually cost something -- the best shape the
        # box can build is not one the facesheet survives. Comparing the two optima catches the
        # case where the constraint pushed the search off the footprint floor as well as the
        # cloud-length floor; testing only the latter mislabels everything from 28 to 43 km/s.
        ceiling_limited=best[0] < best_free - 1e-12,
        d_over_d=d_over_d,
        l_over_d=l_over_d,
        r_foot_over_r=r_foot_over_r,
        eta_capture=eta,
        focusing=focusing_at(l_over_d, r_foot_over_r, d_over_d, path),
        peak_pressure=peak,
        e_eff=e_eff,
        f=None if e_eff is None else eta * (1.0 + e_eff) / 2.0,
    )


# --- Q7: per-pulse radiative fluence -> sacrificial-layer recession -----------------------------
#
# A *diagnostic*, not a gate (Q7). It converts the wall's radiative fluence into how deep the
# renewable sacrificial layer recedes per shot, so the between-pulse MEMS renewal rate (ADR-0014)
# has a number to be sized against.
#
# **It is an upper bound, and a loose one.** The quasi-steady balance sends *all* intercepted flux
# into ablation enthalpy, crediting nothing to two effects that are both large here:
#
#   1. **Re-radiation.** An ablating surface at these fluences runs hot and radiates a large
#      fraction straight back out. Not modelled.
#   2. **Vapor shielding.** ADR-0014's own curtain (`kappa_vapor`) attenuates the incoming flux
#      before it reaches the surface; the ablating-bounce runs show that is not a small correction.
#      This diagnostic reads the *rigid-wall* fluence, which is the unshielded case.
#
# Read the numbers as "no more than this", never as an estimate.

# Silicone-class ablator density [kg/m^3]. ADR-0014's own sanity figure -- "a few um of ablator
# ~ 0.4 kg" over the R = 5 m plate -- implies ~1200, which is where this comes from.
RHO_ABLATOR = 1200.0
# Effective heat of ablation [J/kg]. ADR-0014 parameterises it over the silicone literature range
# and requires the sensitivity be reported rather than a single value quoted.
Q_STAR_BRACKET = (2.0e6, 10.0e6)


def ablation_depth(fluence: float, q_star: float, rho_ablator: float = RHO_ABLATOR) -> float:
    """Recession depth [m] per pulse from a radiative `fluence` [J/m^2].

    ADR-0014's quasi-steady surface energy balance: all intercepted flux goes into the ablation
    enthalpy, so `depth = fluence / (rho_ablator * Q*)`. Conservative by construction -- it credits
    nothing to re-radiation or to conduction into the substrate, both of which reduce the recession.
    """
    return fluence / (rho_ablator * q_star)


@dataclass(frozen=True)
class AblationBand:
    """Recession bracket over the literature `Q*` range."""

    fluence: float
    q_star_lo: float
    q_star_hi: float
    depth_max: float
    depth_min: float


def ablation_bracket(fluence: float, rho_ablator: float = RHO_ABLATOR) -> AblationBand:
    """Recession depth over `Q_STAR_BRACKET`, deep end first.

    The *softest* ablator recedes most, so `depth_max` pairs with `q_star_lo`. Reporting the pair
    rather than a midpoint is ADR-0014's explicit requirement: `Q*` is a parameter of the model, not
    a measurement of this system.
    """
    return AblationBand(
        fluence=fluence,
        q_star_lo=Q_STAR_BRACKET[0],
        q_star_hi=Q_STAR_BRACKET[1],
        depth_max=ablation_depth(fluence, Q_STAR_BRACKET[0], rho_ablator),
        depth_min=ablation_depth(fluence, Q_STAR_BRACKET[1], rho_ablator),
    )


# --- The reported band on f --------------------------------------------------------------------


@dataclass(frozen=True)
class FBand:
    """`f` with its combined uncertainty band."""

    f: float
    freeze_delta: float
    opacity_delta: float
    half_width: float
    lo: float
    hi: float
    #: False when the freeze bracket at this velocity was carried forward rather than measured.
    freeze_measured: bool


def f_band(
    f: float,
    freeze_delta: float,
    opacity_delta: float,
    freeze_measured: bool = True,
) -> FBand:
    """Combine the freeze-timing (Q4) and opacity (Q10/Q18) brackets into one band on `f`.

    **In quadrature, not linearly.** The two are independent: freeze timing is *when* the
    composition stops equilibrating during re-expansion, opacity accuracy is *how well TOPS knows
    kappa*. Adding them linearly would claim a correlation that does not exist and overstate the
    band; in practice the freeze bracket dominates by an order of magnitude anyway.

    `freeze_measured` is carried through rather than folded into the number. Q4 measures the freeze
    bracket at three anchors only -- it tracks an ionization staircase, not a ramp, so it cannot be
    interpolated -- and a reader has to be able to tell a measured bracket from a carried-forward
    one.
    """
    half = math.hypot(freeze_delta, opacity_delta)
    return FBand(
        f=f,
        freeze_delta=freeze_delta,
        opacity_delta=opacity_delta,
        half_width=half,
        lo=f - half,
        hi=f + half,
        freeze_measured=freeze_measured,
    )


# --- The assembled deliverable ------------------------------------------------------------------

DEFAULT_CONTOUR_PATH = Path("data/results/frontier_contour_heavyplate.csv")

CSV_HEADER = (
    "v,rho_ceiling,rho_contour,binds,d_over_d,l_over_d,r_foot_over_r,eta_capture,"
    "focusing,peak_mpa,"
    "e_eff,f,freeze_delta,freeze_measured,opacity_delta,f_lo,f_hi,"
    "wall_fluence,ablation_depth_max_um,ablation_depth_min_um"
)


@dataclass(frozen=True)
class ContourRow:
    """One fully-resolved point on the reported curve."""

    point: ContourPoint
    band: FBand
    ablation: AblationBand


def build_curve(
    velocities: list[float],
    c_stag_at: Callable[[float], float],
    e_eff_at: Callable[[float], Callable[[float], float]],
    fluence_at: Callable[[float], Callable[[float], float]],
    freeze_delta_at: Callable[[float], tuple[float, bool]],
    opacity_delta_at: Callable[[float], float],
    p_limit: float,
    d_over_d: float = D_OVER_D_HEADLINE,
    path: Path = DEFAULT_GEOMETRY_PATH,
) -> list[ContourRow]:
    """Assemble the contour curve with its bands and the Q7 recession diagnostic.

    Every input is a callable rather than a dataframe so this stays pure post-processing over
    whatever the caller has already loaded, and so the module keeps no dependency on the analysis
    that consumes it.
    """
    out: list[ContourRow] = []
    for v in velocities:
        e_of_rho = e_eff_at(v)
        pt = contour_point(
            v,
            c_stag=c_stag_at(v),
            p_limit=p_limit,
            d_over_d=d_over_d,
            path=path,
            e_eff_at=e_of_rho,
        )
        assert pt.f is not None  # e_eff_at was supplied, so f is populated
        freeze_delta, measured = freeze_delta_at(v)
        band = f_band(pt.f, freeze_delta, opacity_delta_at(v), freeze_measured=measured)
        ablation = ablation_bracket(fluence_at(v)(pt.rho_contour))
        out.append(ContourRow(point=pt, band=band, ablation=ablation))
    return out


def write_curve(rows: list[ContourRow], path: Path = DEFAULT_CONTOUR_PATH) -> None:
    """Write the contour deliverable CSV."""
    lines = [CSV_HEADER]
    for r in rows:
        p, b, a = r.point, r.band, r.ablation
        lines.append(
            f"{p.v},{p.rho_ceiling:.6e},{p.rho_contour:.6e},"
            f"{'survivability' if p.ceiling_limited else 'shape_box'},"
            f"{p.d_over_d},{p.l_over_d:.6f},{p.r_foot_over_r:.6f},{p.eta_capture:.6f},"
            f"{p.focusing:.6f},{p.peak_pressure / 1e6:.3f},"
            f"{p.e_eff:.6f},{p.f:.6f},{b.freeze_delta:.6f},{b.freeze_measured},"
            f"{b.opacity_delta:.6e},{b.lo:.6f},{b.hi:.6f},"
            f"{a.fluence:.6e},{a.depth_max * 1e6:.3f},{a.depth_min * 1e6:.3f}"
        )
    path.write_text("\n".join(lines) + "\n")


def pinned_shape(
    pin_v: float,
    c_stag: float,
    p_limit: float,
    d_over_d: float = D_OVER_D_HEADLINE,
    path: Path = DEFAULT_GEOMETRY_PATH,
    e_eff_at: Callable[[float], float] | None = None,
) -> ContourPoint:
    """The single cloud shape flown at every velocity when the schedule is *not* available (Q13).

    Design SS7 treats the cloud as a per-shot schedule `shape(v)` -- the plate is one built object,
    the cloud is not. The pinned curve is the companion that asks what happens without that freedom.

    The shape is set by the **most demanding velocity in the range**, since one shape must survive
    everywhere; pinning it anywhere else yields a curve that fails at the top end, which is not a
    curve anyone can fly. Below `pin_v` it is then more dilute than survivability requires, and a
    dilute cloud radiates away more of the bounce -- so the pinned curve sits at or below the
    floating one, and the gap is the value of being able to schedule.
    """
    return contour_point(
        pin_v,
        c_stag=c_stag,
        p_limit=p_limit,
        d_over_d=d_over_d,
        path=path,
        e_eff_at=e_eff_at,
    )


def pinned_point(
    v: float,
    shape: ContourPoint,
    c_stag: float,
    p_limit: float,
    e_eff_at: Callable[[float], float] | None = None,
) -> ContourPoint:
    """Evaluate an already-pinned `shape` at velocity `v`.

    The shape's density and `eta_capture` are carried over unchanged -- that is what "pinned" means
    -- while `rho_ceiling` and `e_eff` are re-read at `v`. Comparing `rho_contour` against
    `rho_ceiling` on the result tells you whether the pinned shape still survives there; sized at
    the range's worst velocity it always does, with margin to spare lower down, and that unused
    margin is exactly the restitution the schedule would have recovered.
    """
    e_eff = e_eff_at(shape.rho_contour) if e_eff_at is not None else None
    return ContourPoint(
        v=v,
        rho_ceiling=rho_ceiling(v, c_stag, p_limit),
        rho_contour=shape.rho_contour,
        ceiling_limited=False,
        d_over_d=shape.d_over_d,
        l_over_d=shape.l_over_d,
        r_foot_over_r=shape.r_foot_over_r,
        eta_capture=shape.eta_capture,
        e_eff=e_eff,
        f=None if e_eff is None else shape.eta_capture * (1.0 + e_eff) / 2.0,
    )
