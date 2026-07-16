"""Pulse-shape sensitivity assembly (design §13, ADR-0028): do slight pulse-shape changes move
`f` and the per-pulse impulse only slightly?

The mission argument this feeds: delivery shape dispersion is easy for the pushed vehicle's
guidance to absorb **iff** the impulse response is smooth — errors that propagate linearly can be
corrected; a cliff cannot. The measured object is **raw `f(shape)` at the fixed best-survivable
baseline design** (`d/D = 0.1` headline + flat cross-check, `L/D = 0.3`, `r_foot/R = 0.5`,
M = 20), *not* the survivability-constrained frontier (whose argmax is intrinsically
discontinuous, ADR-0028 decision 1). Survivability enters only as a separate **margin check**.

Inputs (the `--shape` sweep, ADR-0028):

- **2D** (`sweep_shape_geometry.jsonl`): `eta_capture` per shape sample on one frozen plate/grid,
  plus refined-resolution repeats of a noise-floor subset.
- **1D** (`sweep_shape_sigma.jsonl`): fresh equilibrium `e_eff` at each sample's physical
  `(rho, L)` via the Σ contract at fixed pulse mass — the taper axis at its mass-weighted mean Σ
  with Σ-hi/lo90 bound rows (the named halt condition); divergence reuses the nominal `e_eff`
  (Σ untouched).

Deliverables: `data/results/shape_sensitivity.csv` + `shape_sensitivity.png`, the normalized
shape sensitivity `S_x = (Δf/f)/(Δx/x)` per axis (CONTEXT; one-sided axes are normalized by the
box extent), a cliff detector (local second-difference outliers above the measured noise floor),
the Σ-profile bound, and the survivability margin over the whole box. The **shape box is an
assumption, not derived dispersion** — the deferred cloud-schedule study owns the real delivery
numbers, and every quoted result carries that caveat.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, fields, replace
from itertools import pairwise
from pathlib import Path
from statistics import median, pstdev

from puffsat.analysis import P_LIMIT_BASELINE, P_LIMIT_HIGHV, _write_csv, reconcile_f

DEFAULT_SHAPE_2D_PATH = Path("data/results/sweep_shape_geometry.jsonl")
DEFAULT_SHAPE_1D_PATH = Path("data/results/sweep_shape_sigma.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/shape_sensitivity.csv")
DEFAULT_PLOT_DIR = Path("data/results")
# The §13 three-point frozen-chemistry spot-check at the dip anchor (`--frozen-shape` sweep).
DEFAULT_FROZEN_SWEEP_PATH = Path("data/results/sweep_frozen_shape.jsonl")
DEFAULT_FROZEN_SUMMARY_PATH = Path("data/results/shape_frozen_spotcheck.csv")

# The two quoted velocity anchors [m/s] — must match the Rust `SHAPE_V`.
V_ANCHORS = (11_000.0, 16_000.0)
# The nominal shape coordinates — must match the Rust `shape_nominal()`.
NOM_RFOOT_OVER_R = 0.5
NOM_L_OVER_D = 0.3
# One-sided axes have no nominal to normalize `Δx/x` by; `S` uses the box extent instead
# (CONTEXT: the quoted `S` then reads "fraction of the assumed box").
BOX_MAX = {"taper_frac": 0.30, "alpha_div": 0.10}
AXES = ("r_foot_over_r", "l_over_d", "taper_frac", "alpha_div")

# Exit criteria (design §13): bounded S (≲ a few), no cliffs above noise, positive survivability
# margin, Σ-profile bound small, frozen spot-check slope comparable to equilibrium.
S_BOUND = 3.0
# A cliff = a segment whose Δf departs from the axis's median-slope prediction by more than ALL
# of: the measured 2*sigma_noise floor (below-noise structure is a pass — indistinguishable from
# flat, design §13), a CLIFF_FACTOR outlier margin over the axis's own median departure, and an
# absolute floor of CLIFF_REL_FLOOR·f (below the Δf ~ 0.005-0.02 scale the study resolves).
CLIFF_FACTOR = 5.0
CLIFF_REL_FLOOR = 0.005
# The Σ-profile bound tripwire in f units (ADR-0028: non-small ⇒ the deferred Σ-resolved
# e_eff(rho) work has become load-bearing and the study halts there as a finding).
SIGMA_BOUND_DELTA_F_MAX = 0.02
# Frozen spot-check pass band: |Δe_eff-across-the-Σ-box| difference vs equilibrium (ADR-0026).
FROZEN_SLOPE_TOL = 0.05
# Two-resolution validity protocol for the coupled 1D rows (2026-07-16 finding: the coupled
# radiation operator has a resolution-onset radiative-collapse instability whose onset dx
# coarsens with rho·v² — see the Rust `Shape1DRecord` docs). The 300-cell headline is accepted
# only if it agrees with the 150-cell stable-window run; otherwise the coarse value stands in,
# provided it is consistent with the (verified-smooth) EOS-only reference.
FINE_COARSE_TOL = 0.01
COARSE_EOS_TOL = 0.02


@dataclass(frozen=True)
class Shape2DRow:
    """One 2D shape-sweep row (`crates/sweep --shape`): the sample, the plate/resolution, and the
    `eta_capture` pieces with the *measured* initialized `p_in` (the §13 normalization)."""

    axis: str
    d_over_d: float
    r_foot_over_r: float
    l_over_d: float
    taper_frac: float
    alpha_div: float
    mach: float
    resolution_scale: float
    eta_capture: float
    restitution_free: float
    restitution_confined: float
    incident_momentum: float
    peak_local_pressure: float


@dataclass(frozen=True)
class Shape1DRow:
    """One 1D Σ-contract row (`crates/sweep --shape`): the sample it serves, its Σ role, and the
    equilibrium coupled-bounce result at the physical `(rho, L)`."""

    v: float
    axis: str
    sigma_role: str
    r_foot_over_r: float
    l_over_d: float
    taper_frac: float
    sigma: float
    rho_impact: float
    length: float
    e_eff: float  # coupled, 300 cells (the headline when the validity protocol accepts it)
    e_eff_coarse: float  # coupled, 150 cells (stable-window cross-check/fallback)
    e_eff_eos: float  # EOS-only reference (verified smooth over the whole box)
    peak_wall_pressure: float
    peak_wall_pressure_coarse: float
    incident_momentum: float
    wall_impulse: float


def resolve_1d(row: Shape1DRow) -> tuple[float, float, bool]:
    """Apply the two-resolution validity protocol to one coupled 1D row: return
    `(e_eff, peak_wall_pressure, solver_valid)`. The 300-cell headline is used when physical and
    in agreement with the 150-cell run; a fine-run radiative collapse falls back to the coarse
    stable-window value (flagged `solver_valid = False`) if that is physical and consistent with
    the EOS-only reference; anything else is a hard error — the row cannot be quoted."""
    fine_ok = 0.0 < row.e_eff < 1.0 and abs(row.e_eff - row.e_eff_coarse) <= FINE_COARSE_TOL
    if fine_ok:
        return row.e_eff, row.peak_wall_pressure, True
    coarse_ok = (
        0.0 < row.e_eff_coarse < 1.0 and abs(row.e_eff_coarse - row.e_eff_eos) <= COARSE_EOS_TOL
    )
    if coarse_ok:
        return row.e_eff_coarse, row.peak_wall_pressure_coarse, False
    raise ValueError(
        f"1D row v={row.v} axis={row.axis} rho={row.rho_impact}: both coupled resolutions "
        f"outside validity (fine {row.e_eff}, coarse {row.e_eff_coarse}, eos {row.e_eff_eos})"
    )


@dataclass(frozen=True)
class ShapePoint:
    """One assembled operating point of the shape box: the sample's coordinates, the two ADR-0003
    factors, `f`, its normalized sensitivity `S` (NaN at the nominal), and the survivability
    margin inputs."""

    v: float
    d_over_d: float
    axis: str
    x: float  # the perturbed coordinate's value (the nominal's own coordinate for "nominal")
    rel_delta: float  # Δx/x_nom (two-sided) or x/x_box (one-sided); 0 at the nominal
    eta_capture: float
    e_eff: float
    f: float
    s: float  # S_x = (Δf/f)/rel_delta; NaN at the nominal
    peak_compressive: float  # 1D physical peak x the sample's concave focusing factor
    survives_baseline: bool
    cliff: bool
    sigma_noise_f: float  # the (v, d/D) assembly's measured noise floor (constant per assembly)
    # False when the 300-cell coupled run radiatively collapsed and the 150-cell stable-window
    # value stands in (the two-resolution validity protocol; see `resolve_1d`).
    solver_valid: bool


def read_shape_2d(path: Path = DEFAULT_SHAPE_2D_PATH) -> list[Shape2DRow]:
    """Parse the 2D shape sweep JSONL (one JSON object per line; blank lines tolerated)."""
    rows: list[Shape2DRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            Shape2DRow(
                axis=str(d["axis"]),
                **{f.name: float(d[f.name]) for f in fields(Shape2DRow) if f.name != "axis"},
            )
        )
    return rows


def read_shape_1d(path: Path = DEFAULT_SHAPE_1D_PATH) -> list[Shape1DRow]:
    """Parse the 1D Σ-contract shape sweep JSONL (one JSON object per line; blanks tolerated)."""
    rows: list[Shape1DRow] = []
    str_fields = ("axis", "sigma_role")
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            Shape1DRow(
                axis=str(d["axis"]),
                sigma_role=str(d["sigma_role"]),
                **{
                    f.name: float(d[f.name]) for f in fields(Shape1DRow) if f.name not in str_fields
                },
            )
        )
    return rows


def _sample_key(
    axis: str, rff: float, lod: float, taper: float, alpha: float
) -> tuple[str, float, float, float, float]:
    """The sample identity both JSONL files share (verbatim round-trip of the same Rust floats)."""
    return (axis, rff, lod, taper, alpha)


def _coordinate(axis: str, rff: float, lod: float, taper: float, alpha: float) -> float:
    """The perturbed coordinate's value for a sample of `axis` (the nominal reports its own)."""
    return {
        "nominal": rff,
        "r_foot_over_r": rff,
        "l_over_d": lod,
        "taper_frac": taper,
        "alpha_div": alpha,
    }[axis]


def _rel_delta(axis: str, x: float) -> float:
    """`Δx/x` for the two-sided axes, `x/x_box` for the one-sided ones, 0 at the nominal."""
    if axis == "nominal":
        return 0.0
    if axis in BOX_MAX:
        return x / BOX_MAX[axis]
    nominal = {"r_foot_over_r": NOM_RFOOT_OVER_R, "l_over_d": NOM_L_OVER_D}[axis]
    return (x - nominal) / nominal


def e_eff_lookup(
    rows1d: list[Shape1DRow], v: float
) -> dict[tuple[str, float, float, float], Shape1DRow]:
    """The per-sample 1D result at anchor `v`, keyed by `(axis, r_foot/R, L/D, taper)` — the
    headline roles only (`sample`/`taper_mean`; the Σ bound rows are read separately)."""
    out: dict[tuple[str, float, float, float], Shape1DRow] = {}
    for r in rows1d:
        if r.v == v and r.sigma_role in ("sample", "taper_mean"):
            out[(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac)] = r
    return out


def assemble(
    rows2d: list[Shape2DRow],
    rows1d: list[Shape1DRow],
    v: float,
    d_over_d: float,
    *,
    resolution_scale: float = 1.0,
) -> list[ShapePoint]:
    """Assemble `f = eta_capture·(1 + e_eff)/2` per shape sample at one `(v, d/D)`: `eta` from the
    fixed-grid 2D run, `e_eff` from the fresh Σ-contract 1D run (the divergence axis reuses the
    nominal's — Σ untouched), and the survivability peak = the 1D physical `peak_wall_pressure`
    scaled by the sample's concave focusing factor (its dished local peak over its flat
    counterpart's, both measured on the same frozen grid)."""
    e1d = e_eff_lookup(rows1d, v)
    nominal_1d = e1d[("nominal", NOM_RFOOT_OVER_R, NOM_L_OVER_D, 0.0)]
    flat_local = {
        _sample_key(
            r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, r.alpha_div
        ): r.peak_local_pressure
        for r in rows2d
        if r.d_over_d == 0.0 and r.resolution_scale == resolution_scale
    }
    sigma_noise = noise_floor(rows2d, rows1d, v, d_over_d)
    points: list[ShapePoint] = []
    for r in rows2d:
        if r.d_over_d != d_over_d or r.resolution_scale != resolution_scale:
            continue
        one_d = (
            nominal_1d
            if r.axis == "alpha_div"
            else e1d[(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac)]
        )
        e_eff, peak_1d, solver_valid = resolve_1d(one_d)
        key = _sample_key(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, r.alpha_div)
        focusing = r.peak_local_pressure / flat_local[key] if d_over_d > 0.0 else 1.0
        peak = peak_1d * focusing
        x = _coordinate(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, r.alpha_div)
        points.append(
            ShapePoint(
                v=v,
                d_over_d=d_over_d,
                axis=r.axis,
                x=x,
                rel_delta=_rel_delta(r.axis, x),
                eta_capture=r.eta_capture,
                e_eff=e_eff,
                f=reconcile_f(r.eta_capture, e_eff),
                s=math.nan,  # filled by `with_sensitivities`
                peak_compressive=peak,
                survives_baseline=peak < P_LIMIT_BASELINE,
                cliff=False,  # filled by `with_sensitivities`
                sigma_noise_f=sigma_noise,
                solver_valid=solver_valid,
            )
        )
    return with_sensitivities(points)


def noise_floor(
    rows2d: list[Shape2DRow], rows1d: list[Shape1DRow], v: float, d_over_d: float
) -> float:
    """The measured noise floor `sigma_noise` in `f` units at one `(v, d/D)`: the max population std
    of `f` over the ≥ 3 resolution repeats of each noise-subset sample (`e_eff` is a 1D quantity,
    so only `eta_capture` varies with the 2D grid). Structure below `2*sigma_noise` is reported as
    below-noise — itself a pass (design §13)."""
    e1d = e_eff_lookup(rows1d, v)
    nominal_1d = e1d[("nominal", NOM_RFOOT_OVER_R, NOM_L_OVER_D, 0.0)]
    by_sample: dict[tuple[str, float, float, float, float], list[float]] = {}
    for r in rows2d:
        if r.d_over_d != d_over_d:
            continue
        key = _sample_key(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, r.alpha_div)
        by_sample.setdefault(key, []).append(r.eta_capture)
    sigma = 0.0
    for key, etas in by_sample.items():
        if len(etas) < 3:
            continue
        axis, rff, lod, taper, _ = key
        one_d = nominal_1d if axis == "alpha_div" else e1d[(axis, rff, lod, taper)]
        e_eff, _, _ = resolve_1d(one_d)
        fs = [reconcile_f(eta, e_eff) for eta in etas]
        sigma = max(sigma, pstdev(fs))
    return sigma


def _axis_curve(points: list[ShapePoint], axis: str) -> list[ShapePoint]:
    """The axis's curve including the nominal, sorted by coordinate (the nominal sits at its own
    coordinate for the two-sided axes and at 0 for the one-sided ones)."""
    nominal = next(p for p in points if p.axis == "nominal")
    if axis in BOX_MAX:
        nominal_x = 0.0
    else:
        nominal_x = NOM_RFOOT_OVER_R if axis == "r_foot_over_r" else NOM_L_OVER_D
    anchor = replace(nominal, axis=axis, x=nominal_x)
    curve = [p for p in points if p.axis == axis] + [anchor]
    return sorted(curve, key=lambda p: p.x)


def cliff_flags(curve: list[ShapePoint], sigma_noise: float, f_ref: float) -> list[bool]:
    """Flag segments (attributed to each segment's right point, so one flag per `curve[1:]`)
    whose Δf departs from the axis's median-slope prediction beyond the threshold — the
    uneven-grid form of a second-difference outlier vs the local slope (design §13). A flagged
    cliff must additionally survive grid refinement before it is called physical."""
    slopes = [(b.f - a.f) / (b.x - a.x) for a, b in pairwise(curve)]
    if not slopes:
        return []
    med_slope = median(slopes)
    excess = [
        abs((s - med_slope) * (b.x - a.x))
        for s, (a, b) in zip(slopes, pairwise(curve), strict=True)
    ]
    threshold = max(
        2.0 * sigma_noise,
        CLIFF_FACTOR * median(excess),
        CLIFF_REL_FLOOR * abs(f_ref),
    )
    return [e > threshold for e in excess]


def with_sensitivities(points: list[ShapePoint]) -> list[ShapePoint]:
    """Fill each non-nominal point's `S = (Δf/f_nom)/rel_delta` and its cliff flag (computed on
    the full axis curve including the nominal)."""
    nominal = next(p for p in points if p.axis == "nominal")
    flagged: set[tuple[str, float]] = set()
    for axis in AXES:
        curve = _axis_curve(points, axis)
        flags = cliff_flags(curve, nominal.sigma_noise_f, nominal.f)
        for right, flag in zip(curve[1:], flags, strict=True):
            if flag:
                flagged.add((axis, right.x))
    out: list[ShapePoint] = []
    for p in points:
        s = math.nan
        if p.axis != "nominal" and p.rel_delta != 0.0:
            s = ((p.f - nominal.f) / nominal.f) / p.rel_delta
        out.append(replace(p, s=s, cliff=(p.axis, p.x) in flagged))
    return out


def refinement_verdict(
    rows2d: list[Shape2DRow],
    rows1d: list[Shape1DRow],
    points: list[ShapePoint],
    v: float,
    d_over_d: float,
) -> list[tuple[str, float, bool | None]]:
    """The §13 rule: a flagged cliff must survive grid refinement before it is called physical.
    For each flagged point that has refined-resolution repeats, substitute its refined `f` into
    the base curve and re-run the detector — the flag *survives* only if it re-fires at every
    refined resolution. Returns `(axis, x, verdict)` per flag; `verdict = None` means the sample
    has no refined repeats yet (add it to the Rust noise subset and re-run)."""
    e1d = e_eff_lookup(rows1d, v)
    nominal_1d = e1d[("nominal", NOM_RFOOT_OVER_R, NOM_L_OVER_D, 0.0)]
    nominal = next(p for p in points if p.axis == "nominal")
    verdicts: list[tuple[str, float, bool | None]] = []
    for p in (q for q in points if q.cliff):
        one_d = (
            nominal_1d
            if p.axis == "alpha_div"
            else next(
                r
                for r in rows1d
                if r.v == v
                and r.axis == p.axis
                and r.sigma_role in ("sample", "taper_mean")
                and math.isclose(
                    _coordinate(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, 0.0), p.x
                )
            )
        )
        e_eff, _, _ = resolve_1d(one_d)
        refined = [
            r
            for r in rows2d
            if r.d_over_d == d_over_d
            and r.resolution_scale > 1.0
            and r.axis == p.axis
            and math.isclose(
                _coordinate(r.axis, r.r_foot_over_r, r.l_over_d, r.taper_frac, r.alpha_div), p.x
            )
        ]
        if not refined:
            verdicts.append((p.axis, p.x, None))
            continue
        survives = True
        curve = _axis_curve(points, p.axis)
        idx = next(i for i, q in enumerate(curve) if math.isclose(q.x, p.x))
        for r in refined:
            substituted = list(curve)
            substituted[idx] = replace(curve[idx], f=reconcile_f(r.eta_capture, e_eff))
            flags = cliff_flags(substituted, nominal.sigma_noise_f, nominal.f)
            if not flags[idx - 1]:
                survives = False
        verdicts.append((p.axis, p.x, survives))
    return verdicts


def max_abs_s_per_axis(points: list[ShapePoint]) -> dict[str, float]:
    """The quotable `max |S_x|` per axis over one assembly's box samples."""
    out: dict[str, float] = {}
    for axis in AXES:
        ss = [abs(p.s) for p in points if p.axis == axis and not math.isnan(p.s)]
        if ss:
            out[axis] = max(ss)
    return out


def sigma_profile_bound(
    rows1d: list[Shape1DRow], points: list[ShapePoint], v: float
) -> tuple[float, float]:
    """The taper Σ-profile bound (ADR-0028 halt condition): the spread of `e_eff` between the
    widest taper's Σ-hi/lo90 bound rows and its mass-weighted-mean row, and its `f` impact at
    that sample's `eta_capture`. Returns `(Δe_eff_max, Δf)`."""
    widest = BOX_MAX["taper_frac"]
    roles = {
        r.sigma_role: resolve_1d(r)[0]
        for r in rows1d
        if r.v == v and r.sigma_role.startswith("sigma_")
    }
    mean = next(
        resolve_1d(r)[0]
        for r in rows1d
        if r.v == v and r.sigma_role == "taper_mean" and abs(r.taper_frac - widest) < 1e-12
    )
    delta_e = max(abs(roles["sigma_hi"] - mean), abs(roles["sigma_lo90"] - mean))
    eta = next(
        (p.eta_capture for p in points if p.axis == "taper_frac" and abs(p.x - widest) < 1e-12),
        1.0,
    )
    return delta_e, eta * delta_e / 2.0


@dataclass(frozen=True)
class FrozenShapeRow:
    """One `--frozen-shape` spot-check row (the ADR-0026 three-curve schema at one Σ point)."""

    v: float
    rho_impact: float
    e_eff_eq: float
    e_eff_frozen_rebound: float
    e_eff_frozen_all: float
    swap_energy_jump_frac: float


def read_frozen_shape(path: Path = DEFAULT_FROZEN_SWEEP_PATH) -> list[FrozenShapeRow]:
    """Parse the `--frozen-shape` sweep JSONL (one JSON object per line; blanks tolerated)."""
    rows: list[FrozenShapeRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(FrozenShapeRow(**{f.name: float(d[f.name]) for f in fields(FrozenShapeRow)}))
    return rows


def frozen_slope_check(rows: list[FrozenShapeRow]) -> tuple[float, float, bool]:
    """Compare the frozen-rebound `Δe_eff` across the Σ box against the equilibrium one (the §13
    spot-check: the smoothness claim must not silently inherit an equilibrium-only slope).
    Returns `(Δe_eq, Δe_frozen, pass)` where Δ is endpoint-to-endpoint over the box."""
    ordered = sorted(rows, key=lambda r: r.rho_impact)
    d_eq = ordered[-1].e_eff_eq - ordered[0].e_eff_eq
    d_frozen = ordered[-1].e_eff_frozen_rebound - ordered[0].e_eff_frozen_rebound
    return d_eq, d_frozen, abs(d_frozen - d_eq) <= FROZEN_SLOPE_TOL


def write_summary(points: list[ShapePoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write every assembled point to the committed CSV (header = `ShapePoint` field names)."""
    header = [f.name for f in fields(ShapePoint)]
    ordered = sorted(points, key=lambda p: (p.v, -p.d_over_d, p.axis, p.x))
    _write_csv(header, ([getattr(p, name) for name in header] for p in ordered), path)


def plot_shape(points: list[ShapePoint], out_dir: Path = DEFAULT_PLOT_DIR) -> list[Path]:
    """Render `f` vs each shape axis (2x2 panels), one curve per `(anchor, plate)`, with the
    `+/-2*sigma_noise` band on the nominal. matplotlib is imported lazily (`sci` extra)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    labels = {
        "r_foot_over_r": r"footprint $r_\mathrm{foot}/R$",
        "l_over_d": r"aspect $L/D$",
        "taper_frac": r"edge taper (fraction of $r_\mathrm{foot}$)",
        "alpha_div": r"radial divergence $\alpha$",
    }
    fig, axes_grid = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
    series = sorted({(p.v, p.d_over_d) for p in points}, key=lambda t: (t[0], -t[1]))
    for ax, axis in zip(axes_grid.flat, AXES, strict=False):
        for v, dd in series:
            sel = [p for p in points if p.v == v and p.d_over_d == dd]
            curve = _axis_curve(sel, axis)
            plate = "concave" if dd > 0.0 else "flat"
            (line,) = ax.plot(
                [p.x for p in curve],
                [p.f for p in curve],
                "o-",
                label=f"{v / 1000:.0f} km/s, {plate}",
            )
            nominal = next(p for p in sel if p.axis == "nominal")
            if axis in BOX_MAX:
                x_nom = 0.0
            else:
                x_nom = NOM_RFOOT_OVER_R if axis == "r_foot_over_r" else NOM_L_OVER_D
            ax.errorbar(
                [x_nom],
                [nominal.f],
                yerr=[2.0 * nominal.sigma_noise_f],
                fmt="s",
                color=line.get_color(),
                capsize=4,
            )
        ax.set_xlabel(labels[axis])
        ax.grid(True, alpha=0.3)
    for ax in axes_grid[:, 0]:
        ax.set_ylabel(r"$f = \eta_\mathrm{capture}\,(1+e_\mathrm{eff})/2$")
    axes_grid.flat[0].legend(loc="best", fontsize="small")
    fig.suptitle(
        "Pulse-shape sensitivity at the fixed baseline design (raw $f$, one frozen grid;\n"
        "the shape box is an assumption — square markers: nominal $\\pm 2\\sigma_\\mathrm{noise}$)"
    )
    fig.tight_layout()
    path = out_dir / "shape_sensitivity.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _run_frozen(frozen_path: Path, summary_path: Path) -> None:
    """The `--frozen` path: the three-point dip-anchor spot-check — the frozen-rebound slope
    across the Σ box must be comparably gentle to the equilibrium one (design §13)."""
    rows = read_frozen_shape(frozen_path)
    d_eq, d_frozen, ok = frozen_slope_check(rows)
    header = [f.name for f in fields(FrozenShapeRow)]
    ordered = sorted(rows, key=lambda r: r.rho_impact)
    _write_csv(header, ([getattr(r, name) for name in header] for r in ordered), summary_path)
    jump_max = max(abs(r.swap_energy_jump_frac) for r in rows)
    print("python: shape frozen spot-check at the dip anchor (3 Σ points, ADR-0026):")
    for r in ordered:
        print(
            f"    rho={r.rho_impact:.3f}: e_eff eq={r.e_eff_eq:.4f} "
            f"frozen-rebound={r.e_eff_frozen_rebound:.4f} frozen-all={r.e_eff_frozen_all:.4f}"
        )
    print(
        f"python: Δe_eff across the Σ box — equilibrium {d_eq:+.4f}, frozen-rebound "
        f"{d_frozen:+.4f} (|diff| ≤ {FROZEN_SLOPE_TOL}: {'PASS' if ok else 'FAIL'}); "
        f"splice jump max {jump_max:.2e}"
    )
    print(f"python: wrote {summary_path}")


def main() -> None:
    """Assemble the shape-sensitivity study; write the CSV + figure; print the §13 verdicts."""
    parser = argparse.ArgumentParser(
        description="Pulse-shape sensitivity: raw f over the shape box at the fixed design."
    )
    parser.add_argument("--shape-2d", type=Path, default=DEFAULT_SHAPE_2D_PATH)
    parser.add_argument("--shape-1d", type=Path, default=DEFAULT_SHAPE_1D_PATH)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument(
        "--frozen",
        action="store_true",
        help="run the three-point dip-anchor frozen spot-check instead of the headline assembly",
    )
    parser.add_argument("--frozen-sweep", type=Path, default=DEFAULT_FROZEN_SWEEP_PATH)
    args = parser.parse_args()

    if args.frozen:
        summary = (
            args.summary if args.summary != DEFAULT_SUMMARY_PATH else DEFAULT_FROZEN_SUMMARY_PATH
        )
        _run_frozen(args.frozen_sweep, summary)
        return

    rows2d = read_shape_2d(args.shape_2d)
    rows1d = read_shape_1d(args.shape_1d)
    plates = sorted({r.d_over_d for r in rows2d}, reverse=True)

    all_points: list[ShapePoint] = []
    print("python: pulse-shape sensitivity at the fixed baseline design (shape box = assumption):")
    for v in V_ANCHORS:
        for dd in plates:
            pts = assemble(rows2d, rows1d, v, dd)
            all_points.extend(pts)
            nominal = next(p for p in pts if p.axis == "nominal")
            s_axis = max_abs_s_per_axis(pts)
            cliffs = [p for p in pts if p.cliff]
            margin = min((P_LIMIT_BASELINE - p.peak_compressive) / P_LIMIT_BASELINE for p in pts)
            relaxed_limit = max(P_LIMIT_HIGHV)
            margin_relaxed = min((relaxed_limit - p.peak_compressive) / relaxed_limit for p in pts)
            plate = "concave d/D=0.1" if dd > 0.0 else "flat"
            print(
                f"  v={v / 1000:4.0f} km/s, {plate}: f_nom={nominal.f:.3f} "
                f"(eta={nominal.eta_capture:.3f}, e_eff={nominal.e_eff:.3f}), "
                f"sigma_noise={nominal.sigma_noise_f:.4f}"
            )
            print(
                "      max|S|: "
                + ", ".join(f"{a}={s:.2f}" for a, s in s_axis.items())
                + f" (bound {S_BOUND}); cliffs {len(cliffs)}; min survivability margin "
                f"{margin:+.1%} of {P_LIMIT_BASELINE / 1e6:.0f} MPa "
                f"({margin_relaxed:+.1%} of the relaxed {relaxed_limit / 1e6:.0f} MPa)"
            )
            verdicts = refinement_verdict(rows2d, rows1d, pts, v, dd)
            for (axis, x, survives), p in zip(verdicts, cliffs, strict=True):
                state = {
                    True: "SURVIVES refinement — physical",
                    False: "does not survive refinement — grid noise, not physical",
                    None: "UNREFINED — add the sample to the noise subset and re-run",
                }[survives]
                print(f"      CLIFF FLAG: axis={axis} x={x:.3f} f={p.f:.4f} -> {state}")
            fallbacks = sorted(
                {(p.axis, p.x) for p in pts if not p.solver_valid}, key=lambda t: (t[0], t[1])
            )
            if fallbacks:
                print(
                    "      SOLVER WINDOW: 300-cell coupled run radiatively collapsed at "
                    + ", ".join(f"{a}={x:.3f}" for a, x in fallbacks)
                    + " — 150-cell stable-window values stand in (see Shape1DRecord docs)"
                )

    for v in V_ANCHORS:
        pts_v = [p for p in all_points if p.v == v and p.d_over_d > 0.0]
        delta_e, delta_f = sigma_profile_bound(rows1d, pts_v, v)
        verdict = "ok" if delta_f <= SIGMA_BOUND_DELTA_F_MAX else "NON-SMALL — halt (ADR-0028)"
        print(
            f"python: taper Σ-profile bound at v={v / 1000:.0f} km/s: Δe_eff={delta_e:.4f} "
            f"⇒ Δf={delta_f:.4f} ({verdict})"
        )

    write_summary(all_points, args.summary)
    figs = plot_shape(all_points, args.plot_dir)
    print(f"python: wrote {args.summary} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


if __name__ == "__main__":
    main()
