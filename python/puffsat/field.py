"""Rung 1 -- the magnetostatic field of the graded column (N5, and `p_design(z)` for N3).

**Why this module has to exist.** Neither the paper nor this repo has a field model that can answer
N5's question. The paper's entire magnet specification is four numbers -- 20 T at 1 m, 12 T at 3 m,
9 T at 6 m, 5 T at the exit -- with no coil count, no radii, no currents and no turn count; its
`eps_b` integral models the magnet as a single current sheet. This repo's only field model is
`electrothermal.local_field`, a two-point `20 T -> 5 T` flux-conservation relation. **Flux
conservation gives `|B|` along a flux tube and says nothing about topology**, so it structurally
cannot answer "is there a local minimum off-axis" -- it assumes the tube whose existence is the
question.

**The question N5 actually asks.** `sec:jet_efficiency` defends residence with "the graded profile
of 20 T at the chamber falling to 5 T at the exit has no local minimum anywhere along it, so
nothing can sit in it. Every gram has a downhill path out." That claim is about the *on-axis*
profile. A local `|B|` minimum is a magnetic mirror: plasma in one is trapped unless its pitch
angle falls inside the loss cone. The radial direction is where a real winding is least monotonic,
because **discrete coils ripple** -- `|B|` peaks at each coil and dips between them, and the dip
deepens fast as you approach the winding.

**So the answer is contingent on a design decision the paper has not made**, and this module is
built to say so quantitatively rather than to guess a winding and report one number. The winding
is a swept input: many coils approaches a continuous current sheet with no ripple, few coils is
chunky and traps. What N5 returns is a *constraint on admissible windings*.

## The tension inside the paper

`sec:jet_efficiency` argues residence from "a solenoid's winding is **continuous**, so the bore is
walled by field rather than fenced by it" -- contrasting Schilling's 32-strut cage. But a *uniform*
solenoid cannot produce a graded profile at all; grading it means varying turn density along the
column, which is exactly what makes the winding discrete-able. The paper's own two statements about
the magnet pull in opposite directions, and nothing in it picks a resolution.

## What is assumed here, and stated because it is a choice

- **The profile between the four stations.** They fit `B = 19.80 z^-0.4405` to within 1.7%, so the
  four numbers are a smooth power law rather than four independent requirements. Interpolating on
  that fit is safe; `design_field` does it and `PROFILE_FIT_MAX_ERROR` pins the quality.
- **Inboard of the first station.** The power law diverges at `z -> 0`, so `B` is held at its
  1 m value for `z < 1 m`. The paper says "20 T at the chamber" elsewhere, which is consistent.
- **Coil radius.** 3.5 m, just outside the 3.0 m bore -- the paper says only that "a 3 m bore
  means coils about 7 m across". Swept, because ripple depth depends on standoff from the winding.

## What this module does *not* do

No plasma. This is vacuum magnetostatics: `curl B = 0` outside the conductors, so the plume's own
diamagnetism is absent. At the standoff condition the plume is a `beta = 1` object and *will*
push the field out, which reduces `|B|` in the bore and deepens any ripple rather than smoothing
it. So a trap found here is real and a trap not found here is not proof of absence -- the vacuum
field is the optimistic case.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MU0 = 4.0e-7 * math.pi
"""Vacuum permeability [H/m]."""

GRADED_STATIONS: tuple[tuple[float, float], ...] = (
    (1.0, 20.0),
    (3.0, 12.0),
    (6.0, 9.0),
    (23.0, 5.0),
)
"""The paper's `(z [m], B [T])` stations (tex:1129) -- the whole magnet specification."""

COLUMN_LENGTH_M = 23.0
"""Column length [m], the paper's value. `expansion.FIELD_LENGTH` is 23.8; see the ledger."""

BORE_RADIUS_M = 3.0
"""Bore *radius* [m]. Confirmed a radius by the paper's own 28 m^2 cross-section."""

COIL_RADIUS_M = 3.5
"""Default winding radius [m] -- just outside the bore. Swept; the paper pins only ~7 m across."""

PROFILE_FIT_MAX_ERROR = 0.02
"""The power-law fit reproduces all four stations to better than this. Pinned by test."""

DEFAULT_OUTPUT_DIR = Path("data/results")


# ---- Complete elliptic integrals, by AGM (no SciPy dependency) -----------------------------------


def elliptic_k_e(m: float) -> tuple[float, float]:
    """Complete elliptic integrals `K(m)`, `E(m)` of the first and second kind, parameter `m = k^2`.

    Arithmetic-geometric mean, which converges quadratically -- six iterations is machine precision
    over the range a coil field needs. With `a0 = 1`, `b0 = sqrt(1-m)`, `c0 = sqrt(m)` and
    `a_{n+1} = (a_n+b_n)/2`, `b_{n+1} = sqrt(a_n b_n)`, `c_{n+1} = (a_n-b_n)/2`,

        K = pi/(2 a_N),    E = K (1 - sum_n 2^{n-1} c_n^2).

    Implemented rather than imported because the project's core dependency set is numpy alone
    (`pyproject.toml`), and the tamper ledger sets the precedent of writing the one numerical
    routine a module needs instead of pulling SciPy in for it.
    """
    if not 0.0 <= m < 1.0:
        raise ValueError("m must lie in [0, 1); m = 1 is the singular limit of a filamentary loop")
    a, b, c = 1.0, math.sqrt(1.0 - m), math.sqrt(m)
    total = 0.5 * c * c
    power = 1.0
    for _ in range(64):
        if abs(c) < 1e-16:
            break
        a, b, c = 0.5 * (a + b), math.sqrt(a * b), 0.5 * (a - b)
        total += power * c * c
        power *= 2.0
    k = math.pi / (2.0 * a)
    return k, k * (1.0 - total)


# ---- The field of one circular loop --------------------------------------------------------------


def loop_field_on_axis(radius: float, current: float, z: float) -> float:
    """On-axis axial field of a current loop [T]: `mu0 I a^2 / (2 (a^2+z^2)^{3/2})`.

    Closed form, and the acceptance test for `loop_field`: the off-axis expression must reduce to
    this as `r -> 0`, which it does because `K(0) = E(0) = pi/2`.
    """
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    return float(MU0 * current * radius * radius / (2.0 * (radius * radius + z * z) ** 1.5))


def loop_field(radius: float, current: float, r: float, z: float) -> tuple[float, float]:
    """Off-axis field `(B_r, B_z)` [T] of a circular loop, via complete elliptic integrals.

    With `Q = (a+r)^2 + z^2` and `m = 4 a r / Q`,

        B_z = mu0 I /(2 pi sqrt(Q)) [ K(m) + E(m) (a^2 - r^2 - z^2)/((a-r)^2 + z^2) ]
        B_r = mu0 I z/(2 pi r sqrt(Q)) [ -K(m) + E(m) (a^2 + r^2 + z^2)/((a-r)^2 + z^2) ]

    `B_r` vanishes on the axis by symmetry, where the bracket is an exact `0/0`; it is returned as
    zero below a radial tolerance rather than expanded, because every field line this module traces
    starts at finite radius.

    The filament is singular *on* the conductor (`r = a`, `z = 0`), which is physical for a
    zero-thickness loop and not a regime any plume occupies -- the bore stops well inside it.
    """
    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if r < 0.0:
        raise ValueError("r must be non-negative")
    q = (radius + r) ** 2 + z * z
    denom = (radius - r) ** 2 + z * z
    if denom <= 1e-12 * radius * radius:
        raise ValueError("evaluated on the conductor, where a filamentary loop's field diverges")
    m = 4.0 * radius * r / q
    k_int, e_int = elliptic_k_e(m)
    pref = MU0 * current / (2.0 * math.pi * math.sqrt(q))
    b_z = pref * (k_int + e_int * (radius * radius - r * r - z * z) / denom)
    if r < 1e-9:
        return 0.0, b_z
    b_r = pref * (z / r) * (-k_int + e_int * (radius * radius + r * r + z * z) / denom)
    return b_r, b_z


# ---- A stack of coils ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Coil:
    """One circular loop of the winding: axial position, radius, and total ampere-turns."""

    z_m: float
    radius_m: float
    current_a: float


@dataclass(frozen=True)
class CoilStack:
    """A discretised winding. `n_coils` is the sweep axis N5's answer depends on."""

    coils: tuple[Coil, ...]

    @property
    def n_coils(self) -> int:
        return len(self.coils)

    @property
    def total_ampere_turns(self) -> float:
        return sum(c.current_a for c in self.coils)

    def field(self, r: float, z: float) -> tuple[float, float]:
        """Superposed `(B_r, B_z)` [T] at a point. Linear superposition -- vacuum, no materials."""
        b_r = 0.0
        b_z = 0.0
        for c in self.coils:
            dr, dz = loop_field(c.radius_m, c.current_a, r, z - c.z_m)
            b_r += dr
            b_z += dz
        return b_r, b_z

    def magnitude(self, r: float, z: float) -> float:
        """`|B|` [T] -- the quantity whose local minima are magnetic traps."""
        b_r, b_z = self.field(r, z)
        return math.hypot(b_r, b_z)

    def on_axis(self, z: float) -> float:
        """On-axis `B_z` [T] by the closed form, so the fit never depends on the elliptic path."""
        return sum(loop_field_on_axis(c.radius_m, c.current_a, z - c.z_m) for c in self.coils)


# ---- The design profile the winding has to reproduce ---------------------------------------------


@dataclass(frozen=True)
class ProfileFit:
    """The power law through the paper's four stations, and how well it holds."""

    b0_t: float
    exponent: float
    max_rel_error: float
    z_ref_m: float
    """First station of the stations this was fitted to -- where the flat-inboard hold begins.

    Carried on the fit rather than read from `GRADED_STATIONS`, because a fit to a *scaled* column
    has its own first station. Hardcoding the module constant silently clamped a 16 m column's
    profile against a 23 m column's first station and understated its demanded field by 15%.
    """

    def field(self, z: float) -> float:
        """`B(z)` [T], held flat inboard of the first station where the power law diverges."""
        return float(self.b0_t * max(z, self.z_ref_m) ** self.exponent)


def fit_profile(stations: tuple[tuple[float, float], ...] = GRADED_STATIONS) -> ProfileFit:
    """Least-squares power law `B = B0 z^n` through the paper's stations.

    **A finding in its own right.** The four values are not four independent requirements: they lie
    on `B = 19.80 z^-0.4405` to within 1.7%. Since standoff sets `B^2 = 2 mu0 p`, that is a
    snowplow pressure falling as `p ~ z^-0.88`, i.e. very nearly `1/z` -- which is what a front
    spreading in one dimension down a fixed bore should do, and a check the paper never states.
    """
    xs = [math.log(z) for z, _ in stations]
    ys = [math.log(b) for _, b in stations]
    n = len(xs)
    x_bar, y_bar = sum(xs) / n, sum(ys) / n
    s_xy = sum((x - x_bar) * (y - y_bar) for x, y in zip(xs, ys, strict=True))
    s_xx = sum((x - x_bar) ** 2 for x in xs)
    slope = s_xy / s_xx
    b0 = math.exp(y_bar - slope * x_bar)
    worst = max(abs(b0 * z**slope - b) / b for z, b in stations)
    return ProfileFit(
        b0_t=b0, exponent=slope, max_rel_error=worst, z_ref_m=min(z for z, _ in stations)
    )


def design_field(z: float, fit: ProfileFit | None = None) -> float:
    """The design `B(z)` [T] the winding is fitted to reproduce."""
    return (fit or fit_profile()).field(z)


def design_pressure(z: float, fit: ProfileFit | None = None) -> float:
    """`p_design(z) = B(z)^2/(2 mu0)` [Pa] -- what Rung 5 divides the solved expansion by.

    Beta against this profile is 1 at every station by construction, because the profile *is* the
    standoff condition. Only `p_actual/p_design` carries information.
    """
    b = design_field(z, fit)
    return b * b / (2.0 * MU0)


# ---- Fitting a winding to the design profile -----------------------------------------------------


def _extended_design_field(z: float, fit: ProfileFit, length: float) -> float:
    """The design profile continued past both ends of the column, so the winding can overhang.

    Held flat at the first station's value inboard and at the exit value outboard. The overhang
    exists only to move end sag off the region the plume occupies; what it carries there does not
    have to mean anything.
    """
    return fit.field(min(max(z, 0.0), length))


def build_winding(
    n_coils: int,
    *,
    coil_radius: float = COIL_RADIUS_M,
    length: float = COLUMN_LENGTH_M,
    overhang: float = 0.25,
    stations: tuple[tuple[float, float], ...] = GRADED_STATIONS,
) -> CoilStack:
    """Discretise the thin-solenoid current density into `n_coils` rings, then fit one scalar gain.

    **Why not a free least-squares fit for every current.** That was the first attempt and it is
    wrong: with many coils the normal equations are near-singular, the solution oscillates in sign
    (the end rings came out at -29 MA), and the resulting `|B|` ripple is an artifact of the fit
    rather than a property of the winding. The ripple sweep would then have been measuring its own
    regularisation. `test_ripple_deepens_toward_the_winding_and_with_fewer_coils` caught it.

    What replaces it is physical and cannot oscillate. An infinite solenoid gives `B = mu0 K`, so
    the surface current density that produces the design profile is `K(z) = B(z)/mu0`; lumping it
    into `n_coils` rings gives `I_i = K(z_i) dz`. Every current is positive, the grading is
    monotone, and **`n_coils` becomes a clean sweep axis**: large `n` approaches the continuous
    sheet the paper argues from, small `n` is the chunky winding that grading actually implies.

    The single scalar gain absorbs the finite-length end correction, which is a uniform factor
    deep inside a long solenoid. It is one parameter fitted to many targets, so it cannot
    manufacture structure.
    """
    if n_coils < 2:
        raise ValueError("need at least two coils")
    pad = overhang * length
    z_lo, z_hi = -pad, length + pad
    z_coils = np.linspace(z_lo, z_hi, n_coils)
    dz = (z_hi - z_lo) / (n_coils - 1)
    fit = fit_profile(stations)
    currents = [_extended_design_field(float(z), fit, length) / MU0 * dz for z in z_coils]

    trial = CoilStack(
        tuple(
            Coil(z_m=float(z), radius_m=coil_radius, current_a=i)
            for z, i in zip(z_coils, currents, strict=True)
        )
    )
    z_targets = np.linspace(0.0, length, 240)
    produced = np.array([trial.on_axis(float(z)) for z in z_targets])
    wanted = np.array([fit.field(float(z)) for z in z_targets])
    gain = float(produced @ wanted / (produced @ produced))
    return CoilStack(
        tuple(
            Coil(z_m=c.z_m, radius_m=c.radius_m, current_a=c.current_a * gain) for c in trial.coils
        )
    )


def profile_error(
    stack: CoilStack,
    length: float = COLUMN_LENGTH_M,
    n: int = 200,
    z_min: float = 0.0,
    stations: tuple[tuple[float, float], ...] = GRADED_STATIONS,
) -> float:
    """Worst relative error of a built winding against the design profile, on axis, past `z_min`.

    `z_min` exists because the error is not uniform: it is concentrated in the first two metres,
    where the design profile demands a gradient shorter than the coil radius. See
    `chamber_realizability`.
    """
    fit = fit_profile(stations)
    worst = 0.0
    for z in np.linspace(z_min, length, n):
        target = fit.field(float(z))
        worst = max(worst, abs(stack.on_axis(float(z)) - target) / target)
    return worst


@dataclass(frozen=True)
class Realizability:
    """How well a physical winding can reproduce the paper's profile at the chamber end."""

    coil_radius_m: float
    delivered_t: float
    """`|B|` the winding actually produces at the first station."""
    demanded_t: float
    shortfall: float
    """Fractional undershoot at the first station."""
    beta_at_chamber: float
    """`p_design/(B_delivered^2/2mu0)` -- how far past standoff the chamber sits from this alone."""


def chamber_realizability(
    n_coils: int = 72, *, coil_radius: float = COIL_RADIUS_M, station: float = 1.0
) -> Realizability:
    """**A finding, not a diagnostic.** Can a solenoid make the gradient the paper's profile needs?

    A solenoid of radius `a` smooths field structure over a length of order `a`: the field at a
    point is set by conductors within roughly one radius of it. The paper's profile demands 20 T
    at 1 m falling to 12 T at 3 m -- a 40% drop over 2 m -- from coils about 3.5 m in radius. The
    demanded gradient scale is *shorter than the coil radius*, so an all-positive winding cannot
    deliver it and undershoots at the chamber.

    That undershoot is not free. The field is there to stand off the plume, `B^2/2mu0 = p`, so
    delivering less field means the standoff condition fails where the pressure is highest.
    `beta_at_chamber` is how far past 1 the chamber station sits from **realizability alone**,
    before the solved expansion says anything.

    Matching the profile exactly requires counter-wound coils -- which is what an unconstrained
    least-squares fit reaches for, and what no one builds. The lever that genuinely helps is a
    *smaller* coil radius, which trades against bore.
    """
    stack = build_winding(n_coils, coil_radius=coil_radius)
    demanded = design_field(station)
    delivered = stack.on_axis(station)
    return Realizability(
        coil_radius_m=coil_radius,
        delivered_t=delivered,
        demanded_t=demanded,
        shortfall=(demanded - delivered) / demanded,
        beta_at_chamber=(demanded / delivered) ** 2 if delivered > 0.0 else math.inf,
    )


# ---- Ripple, field lines, and the trap question (N5) ---------------------------------------------


@dataclass(frozen=True)
class RippleResult:
    """`|B|` structure along one path -- the object N5's claim is actually about."""

    radius_m: float
    n_coils: int
    b_max_t: float
    b_min_t: float
    n_local_minima: int
    """Interior local minima of `|B|`. Any nonzero count is a magnetic trap."""
    mirror_ratio: float
    """Deepest `B_max/B_min` across a local minimum. 1.0 means no trap."""
    trapped_fraction: float
    """Loss-cone estimate `1 - sqrt(1 - 1/R)` -- share of an isotropic population held."""


def _local_minima(values: list[float]) -> list[int]:
    """Indices of strict interior local minima."""
    return [
        i
        for i in range(1, len(values) - 1)
        if values[i] < values[i - 1] and values[i] < values[i + 1]
    ]


def trapped_fraction(mirror_ratio: float) -> float:
    """Share of an isotropic velocity distribution confined by a mirror of ratio `R`.

    Particles escape only if their pitch angle lies inside the loss cone, `sin^2 theta < 1/R`, so
    the trapped share is `1 - sqrt(1 - 1/R)`. This is the quantity that makes a ripple matter or
    not: `R = 1.01` traps 10% and `R = 2` traps 29%.
    """
    if mirror_ratio <= 1.0:
        return 0.0
    return 1.0 - math.sqrt(1.0 - 1.0 / mirror_ratio)


def scan_axial(
    stack: CoilStack, radius: float, *, length: float = COLUMN_LENGTH_M, n: int = 400
) -> RippleResult:
    """Sample `|B|` along a line of constant radius and report its minima.

    Constant radius rather than a traced field line: inside the bore the field is nearly axial, so
    the two nearly coincide, and a straight sample cannot smuggle in a tracing error. `trace_line`
    does the exact version for the cases this flags.
    """
    zs = [float(z) for z in np.linspace(0.0, length, n)]
    mags = [stack.magnitude(radius, z) for z in zs]
    minima = _local_minima(mags)
    ratio = 1.0
    for i in minima:
        left = max(mags[:i]) if i > 0 else mags[i]
        right = max(mags[i + 1 :]) if i + 1 < len(mags) else mags[i]
        ratio = max(ratio, min(left, right) / mags[i])
    return RippleResult(
        radius_m=radius,
        n_coils=stack.n_coils,
        b_max_t=max(mags),
        b_min_t=min(mags),
        n_local_minima=len(minima),
        mirror_ratio=ratio,
        trapped_fraction=trapped_fraction(ratio),
    )


def trace_line(
    stack: CoilStack, r0: float, z0: float, *, z_end: float, steps: int = 2000
) -> list[tuple[float, float, float]]:
    """Trace a field line from `(r0, z0)` to `z_end`, returning `(z, r, |B|)`.

    RK4 on `dr/dz = B_r/B_z`. Valid while the line stays field-aligned in `z`, which holds inside
    a solenoid bore; it is not a general tracer and will fail near a field null, which is exactly
    where it should stop rather than produce a plausible wrong answer.
    """
    h = (z_end - z0) / steps

    def slope(r: float, z: float) -> float:
        b_r, b_z = stack.field(max(r, 1e-9), z)
        if abs(b_z) < 1e-12:
            raise ValueError(
                f"field line turned at z = {z:.3f}; the tracer is not valid through a null"
            )
        return b_r / b_z

    out: list[tuple[float, float, float]] = []
    r, z = r0, z0
    for _ in range(steps + 1):
        out.append((z, r, stack.magnitude(max(r, 1e-9), z)))
        k1 = slope(r, z)
        k2 = slope(r + 0.5 * h * k1, z + 0.5 * h)
        k3 = slope(r + 0.5 * h * k2, z + 0.5 * h)
        k4 = slope(r + h * k3, z + h)
        r += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        z += h
        if r < 0.0:
            r = -r
    return out


def ripple_sweep(
    coil_counts: tuple[int, ...] = (12, 18, 24, 36, 48, 72, 120, 200),
    radii: tuple[float, ...] = (0.0, 1.0, 2.0, 2.5),
    *,
    coil_radius: float = COIL_RADIUS_M,
) -> list[RippleResult]:
    """Sweep the winding family and report the trap structure at each radius.

    This is N5's deliverable. Many coils approaches a continuous sheet and should show no interior
    minimum at any radius; few coils ripples, and the ripple deepens toward the winding. The output
    is the **constraint on admissible windings** -- the coil count above which the claim holds.
    """
    results: list[RippleResult] = []
    for n in coil_counts:
        stack = build_winding(n, coil_radius=coil_radius)
        for r in radii:
            results.append(scan_axial(stack, r))
    return results


def write_ripple(rows: list[RippleResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "n_coils",
                "radius_m",
                "b_max_T",
                "b_min_T",
                "n_local_minima",
                "mirror_ratio",
                "trapped_fraction",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.n_coils,
                    f"{r.radius_m:g}",
                    f"{r.b_max_t:.4f}",
                    f"{r.b_min_t:.4f}",
                    r.n_local_minima,
                    f"{r.mirror_ratio:.6f}",
                    f"{r.trapped_fraction:.6f}",
                ]
            )


# ---- Can the chamber shortfall be fixed? (N3) ----------------------------------------------------
#
# The shortfall is a SHAPE mismatch, not a strength one: the winding overshoots mid-column while
# undershooting at the chamber, so the demanded B(1)/B(3) of 1.62 arrives as 1.18. Two levers that
# look obvious are both *scales*, and a scale cannot change a ratio:
#
#   - more current      -- fixes the chamber only by over-provisioning everywhere, and pays for it
#                          in stored energy, hence in virial structure mass;
#   - a thinner plume   -- lowers demanded and delivered field together, leaving the ratio
#                          untouched, while widening the bore and so making the mismatch worse.
#
# What does move the ratio is the geometry: the coil radius over the length the field must change
# in. `column_length_scan` is the paper's own dial on that, because eq:bore_from_length ties a
# longer column to a *narrower* bore.

V_STANDOFF_M3 = 660.0
"""Standoff volume [m^3]. Set by `PV = nR_gT` and shape-independent (paper tex:1095)."""

LINER_GAP_M = 0.5
"""Bore-to-winding standoff [m]: graphite liner, vacuum gap, aluminium shell."""


@dataclass(frozen=True)
class ScaleUpCost:
    """What buying the chamber field with more current costs, and what it buys you elsewhere."""

    current_scale: float
    energy_scale: float
    """Stored field energy goes as `B^2`, and the virial structure floor is linear in it."""
    overshoot: tuple[tuple[float, float], ...]
    """`(z, delivered/demanded)` once scaled. Every station past the chamber over-provisions."""


def scale_up_cost(
    n_coils: int = 72, *, coil_radius: float = COIL_RADIUS_M, station: float = 1.0
) -> ScaleUpCost:
    """Uniformly scale the winding until the chamber reaches standoff, and price it.

    Physically safe -- `beta < 1` elsewhere just means the plume is held more firmly than needed --
    but not free. Field energy goes as `B^2`, and `sec:minimum_nozzle`'s virial floor is linear in
    contained energy, so the structure mass moves with `energy_scale`.

    **Better superconductors do not buy this back.** The magnet has two mass terms with different
    physics: conductor mass, which better tape reduces, and structure mass, which comes from having
    to mechanically react the magnetic pressure whatever carries the current. The term that grows
    here is the second one.
    """
    rz = chamber_realizability(n_coils, coil_radius=coil_radius, station=station)
    scale = rz.demanded_t / rz.delivered_t
    stack = build_winding(n_coils, coil_radius=coil_radius)
    return ScaleUpCost(
        current_scale=scale,
        energy_scale=scale * scale,
        overshoot=tuple(
            (z, scale * stack.on_axis(z) / design_field(z)) for z, _ in GRADED_STATIONS
        ),
    )


@dataclass(frozen=True)
class ColumnOption:
    """One column length, the bore `eq:bore_from_length` gives it, and the shortfall."""

    length_m: float
    bore_m: float
    coil_m: float
    gradient_m: float
    """Distance between the first two field stations -- the length the field must change in."""
    brush_ratio: float
    """`coil radius / gradient length` -- the controlling parameter.

    A wide brush cannot draw a narrow line, and everything that helps reduces this number.
    """
    shortfall: float
    beta_at_chamber: float


def column_length_scan(
    lengths: tuple[float, ...] = (16.0, 23.0, 32.0, 50.0, 75.0), *, n_coils: int = 72
) -> list[ColumnOption]:
    """Does a longer column fix the shape? It attacks the brush ratio from both ends at once.

    `eq:bore_from_length` gives `r = sqrt(V/(pi L))`, so a longer column has a **narrower** bore --
    hence narrower coils -- while the pressure profile it has to reproduce stretches over a longer
    distance. Both move `brush_ratio` the right way.

    **Assumption, stated because it is one.** The four field stations are taken to scale with column
    length, since they come from a snowplow front spreading down the column and so are naturally a
    function of `z/L`. The paper derives them only at 23 m.

    It helps substantially and does not fully close: `beta` at the chamber runs 1.65 at the flown
    23 m to 1.26 at 50 m. And it is not free -- `tab:axial_bag`'s conductor column runs 1.19 at
    23 m against 1.76 at 50 m, so the length that fixes the field shape is the one the paper
    rejected on launch-envelope grounds.
    """
    out: list[ColumnOption] = []
    for length in lengths:
        bore = math.sqrt(V_STANDOFF_M3 / (math.pi * length))
        coil = bore + LINER_GAP_M
        scale = length / COLUMN_LENGTH_M
        stations = tuple((z * scale, b) for z, b in GRADED_STATIONS)
        fit = fit_profile(stations)
        stack = build_winding(n_coils, coil_radius=coil, length=length, stations=stations)
        first_z = stations[0][0]
        demanded = fit.field(first_z)
        delivered = stack.on_axis(first_z)
        gradient = stations[1][0] - stations[0][0]
        out.append(
            ColumnOption(
                length_m=length,
                bore_m=bore,
                coil_m=coil,
                gradient_m=gradient,
                brush_ratio=coil / gradient,
                shortfall=(demanded - delivered) / demanded,
                beta_at_chamber=(demanded / delivered) ** 2,
            )
        )
    return out


def main() -> None:
    """Fit the winding, check it against the paper's stations, and sweep for traps."""
    parser = argparse.ArgumentParser(description="Rung 1: the graded column's magnetostatic field")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    fit = fit_profile()
    print("== The paper's four stations are one power law ==")
    print(f"B(z) = {fit.b0_t:.3f} * z^({fit.exponent:.4f}) T,  worst error {fit.max_rel_error:.1%}")
    print(f"{'z [m]':>7} {'paper [T]':>10} {'fit [T]':>9} {'p_design [MPa]':>16}")
    for z, b in GRADED_STATIONS:
        print(f"{z:7.1f} {b:10.1f} {fit.field(z):9.2f} {design_pressure(z, fit) / 1e6:16.2f}")
    print(f"standoff implies snowplow pressure ~ z^{2 * fit.exponent:.2f}, i.e. nearly 1/z")

    print("\n== A built winding reproduces it, except at the chamber ==")
    for n in (24, 72, 200):
        stack = build_winding(n)
        print(
            f"{n:4d} coils: worst error {profile_error(stack):6.2%} over the whole column, "
            f"{profile_error(stack, z_min=2.0):6.2%} past 2 m, "
            f"{stack.total_ampere_turns / 1e6:7.3f} MA-turns"
        )

    print("\n== The chamber gradient is sharper than a solenoid can make ==")
    print(
        f"{'coil radius':>12} {'demanded':>9} {'delivered':>10} {'shortfall':>10} "
        f"{'beta at chamber':>16}"
    )
    for a in (1.5, 2.5, 3.5, 5.0):
        rz = chamber_realizability(coil_radius=a)
        print(
            f"{rz.coil_radius_m:11.1f}m {rz.demanded_t:8.2f}T {rz.delivered_t:9.2f}T "
            f"{rz.shortfall:10.1%} {rz.beta_at_chamber:16.2f}"
        )
    print("A solenoid smooths field structure over about its own radius; the profile demands a")
    print("2 m gradient from 3.5 m coils. Undershooting B is overshooting beta, at the station")
    print("where the design pressure is highest (159 MPa).")

    print("\n== N5: is there an off-axis local minimum? ==")
    print(
        f"{'n_coils':>8} {'r [m]':>7} {'B_max':>8} {'B_min':>8} {'minima':>7} "
        f"{'mirror R':>10} {'trapped':>9}"
    )
    rows = ripple_sweep()
    for r in rows:
        print(
            f"{r.n_coils:8d} {r.radius_m:7.1f} {r.b_max_t:8.2f} {r.b_min_t:8.2f} "
            f"{r.n_local_minima:7d} {r.mirror_ratio:10.4f} {r.trapped_fraction:8.2%}"
        )

    trapping = [r for r in rows if r.n_local_minima > 0]
    if trapping:
        worst = max(trapping, key=lambda r: r.mirror_ratio)
        counts = sorted({r.n_coils for r in trapping})
        print(f"\nInterior minima appear at coil counts {counts}.")
        print(
            f"Worst: {worst.n_coils} coils at r = {worst.radius_m:g} m, mirror ratio "
            f"{worst.mirror_ratio:.3f}, holding {worst.trapped_fraction:.1%} of an isotropic "
            f"population."
        )
    else:
        print("\nNo interior local minimum at any swept coil count or radius.")

    print("\n== Can the chamber shortfall be fixed? ==")
    cost = scale_up_cost()
    print(
        f"more current: {cost.current_scale:.3f}x current reaches standoff at the chamber, "
        f"but costs {cost.energy_scale:.3f}x field energy"
    )
    print(
        "   and over-provisions everywhere else: "
        + ", ".join(
            f"{z:g} m {ratio:+.0%}" for z, ratio in ((z, r - 1.0) for z, r in cost.overshoot)
        )
    )
    print("   virial structure is linear in energy, and better superconductors reduce the")
    print("   conductor term, not that one.")
    print("a thinner plume: lowers demanded and delivered together, so the ratio does not move,")
    print("   and the wider bore it needs makes the brush wider still. Counterproductive.")
    print()
    print(
        f"{'L [m]':>7} {'bore':>7} {'coil':>7} {'gradient':>9} {'brush ratio':>12} "
        f"{'shortfall':>10} {'beta':>6}"
    )
    for opt in column_length_scan():
        print(
            f"{opt.length_m:7.1f} {opt.bore_m:6.2f}m {opt.coil_m:6.2f}m {opt.gradient_m:8.2f}m "
            f"{opt.brush_ratio:12.2f} {opt.shortfall:10.1%} {opt.beta_at_chamber:6.2f}"
        )
    print("The controlling parameter is the brush ratio. Lengthening helps and does not close it;")
    print("the real fix is a bore that converges, which eq:bore_from_length's cylinder forbids.")

    out = args.output_dir / "nozzle_field_ripple.csv"
    write_ripple(rows, out)
    print(f"\nwrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
