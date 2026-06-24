"""Frontier extraction + plots for the e_eff sweeps (B5d-2, extended for Rung C and the
transitional anchor).

Reads a sweep's JSONL (`crates/sweep` output, ADR-0019), builds the `e_eff(rho_impact)` frontier and
the per-density loss decomposition (each ADR-0016 channel as a fraction of the total loss), writes a
summary CSV, and renders the figures. The extraction is **stdlib-only and unit-tested** (no science
deps); plotting imports matplotlib lazily (the `sci` extra), so this module imports and tests fine
without it.

Channels: 1a radiative-to-wall, 1b escape-to-space, 2 conductive (deferred to B-flux), and
3 condensation (Rung C low-v). The `--sweep`/`--summary`/`--plot-dir`/`--tag` CLI options point it
at either the high-v (`sweep.jsonl`) or low-v (`sweep_lowv.jsonl`) results.

**`--axis v` (transitional anchor, ADR-0012):** instead of the `e_eff(rho)` frontier, build the
`e_eff(v)` curves from the two transitional sweeps (`sweep_transitional_eos.jsonl` and
`..._rad.jsonl`): the EOS-only (opacity-independent) curve and the radiation-on (interim-opacity)
curve, the gap between them being the radiative-uncertainty band. The dip locator reports whether
the EOS-only `e_eff(v)` has an interior minimum below the swept endpoints — the dissociation/
ionization specific-heat floor the rung exists to find.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass, fields
from itertools import pairwise
from pathlib import Path

DEFAULT_SWEEP_PATH = Path("data/results/sweep.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/frontier.csv")
DEFAULT_PLOT_DIR = Path("data/results")

# Transitional anchor (ADR-0012): the two `crates/sweep --transitional` outputs and their frontier.
DEFAULT_TRANS_EOS_PATH = Path("data/results/sweep_transitional_eos.jsonl")
DEFAULT_TRANS_RAD_PATH = Path("data/results/sweep_transitional_rad.jsonl")
DEFAULT_TRANS_SUMMARY_PATH = Path("data/results/frontier_transitional.csv")

# Design reference anchors from the prior rungs, marked on the `e_eff(v)` figure for context: the
# 3.2 km/s adiabatic low-v point (Rung C) and the 16 km/s radiation-on point (Rung B, B5d-1). Each
# is `(v [m/s], e_eff)`.
ANCHOR_LOWV = (3200.0, 0.74)
ANCHOR_HIGHV = (16000.0, 0.63)

# Geometry sweep (Rung D follow-on): the `crates/sweep --geometry` output and its f-reconciliation.
DEFAULT_GEOMETRY_PATH = Path("data/results/sweep_geometry.jsonl")
DEFAULT_GEOMETRY_SUMMARY_PATH = Path("data/results/frontier_geometry.csv")
DEFAULT_SURVIVABILITY_SUMMARY_PATH = Path("data/results/frontier_survivability.csv")
DEFAULT_MARGIN_SUMMARY_PATH = Path("data/results/frontier_margin.csv")

# 1D `e_eff` anchors for the `f = eta_capture·(1+e_eff)/2` reconciliation (ADR-0003). `EEFF_DIP` is
# the transitional EOS worst case (the conservative floor, ADR-0012); `EEFF_HIGHV` the 16 km/s point
# (Rung B). `eta_capture` is geometry-dominated and ~velocity-independent, so the first bracket
# pairs it with these two `e_eff` scenarios; the full `Sigma`-resolved `e_eff(rho(r_foot))` lookup
# is the deferred refinement (the dual-curve `f(v)` deliverable, ADR-0013).
EEFF_DIP = 0.57
EEFF_HIGHV = 0.63
USEFUL_F_GATE = 0.8  # the useful-`f` gate (ADR-0009), marked on the figure

# --- Rung S: survivability frontier (design §7, ADR-0010/0011) ---
# Peak facesheet pressure is the stagnation pressure `c_stag·rho·v²`; `c_stag` (~2.0 at 16 km/s)
# is backed out of the 1D `peak_wall_force` so the frontier inherits the kernel's measured number,
# not an assumed ideal. The structural limits are the SiC+Ti facesheet's: a compressive `P_limit`
# (the §5 band) and a reflected-tensile spall limit at the SiC-Ti interface (ADR-0011), where
# `|R| ~ 0.15` of the incident compressive returns as tension.
PULSE_MASS_KG = 25.0  # gas delivered per PuffSat (design §2)
PLATE_RADIUS_M = 5.0  # plate radius R, fixed (design §2/§7)
V_DIP = 11_000.0  # transitional worst-case velocity (~11 km/s, ADR-0012), paired with EEFF_DIP
P_LIMIT_BASELINE = 400.0e6  # conservative floor of the §5 SiC+Ti band (Pa)
P_LIMIT_HIGHV = (700.0e6, 900.0e6)  # relaxed limits swept at the 16 km/s anchor (design §7)
REFLECT_FRAC = 0.15  # |R| at the SiC-Ti impedance step -> reflected tensile fraction (ADR-0011)
SIC_SPALL_LO = 0.3e9  # SiC dynamic spall strength, conservative end (ADR-0011)
SIC_SPALL_HI = 1.0e9  # SiC dynamic spall strength, upper end

# --- closed-form `f`-margin exploration (design §7, ADR-0010 amendment) ---
# Peak facesheet pressure `c_stag·rho·v²` is intensive (set by the gas at the wall); via the Σ
# contract `rho = m/(2π·(L/D)·(r_foot/R)³·R³)` it scales analytically as `rho ∝ m/R³`, while
# `eta_capture(r_foot/R)` is scale-invariant. So a wider plate (R↑) or a smaller pulse (m↓) only
# relaxes the pressure ceiling, admitting denser/higher-`eta` shapes and buying `f` back above the
# 0.8 gate (ADR-0009). R/m are otherwise PINNED — R by the vehicle dry-mass budget, m by the
# per-pulse thrust x pulse-rate — so this maps the f-side of a system trade, not an optimum.
MARGIN_PLATE_RADII_M = (5.0, 5.5, 6.0, 6.5, 7.0)  # current R=5 m up to +40% (design §2/§7)
MARGIN_MASSES_KG = (25.0, 20.0, 15.0)  # current m=25 kg down to 15 kg (more, smaller pulses)


@dataclass(frozen=True)
class SweepRow:
    """One sweep result row (the fields the analysis needs from the JSONL schema)."""

    rho_impact: float
    v: float
    e_eff: float
    peak_wall_force: float
    loss_radiative_wall: float
    loss_escape_space: float
    loss_conductive: float
    loss_condensation: float


@dataclass(frozen=True)
class FrontierPoint:
    """One point on the `e_eff(rho)` frontier with its loss decomposition (fractions sum to 1 when
    there is any loss, else all 0)."""

    rho_impact: float
    e_eff: float
    total_loss: float
    frac_radiative_wall: float
    frac_escape_space: float
    frac_conductive: float
    frac_condensation: float


def read_sweep(path: Path = DEFAULT_SWEEP_PATH) -> list[SweepRow]:
    """Parse the JSONL sweep results (one JSON object per line; blank lines tolerated)."""
    rows: list[SweepRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            SweepRow(
                rho_impact=float(d["rho_impact"]),
                v=float(d["v"]),
                e_eff=float(d["e_eff"]),
                peak_wall_force=float(d.get("peak_wall_force", 0.0)),
                loss_radiative_wall=float(d["loss_radiative_wall"]),
                loss_escape_space=float(d["loss_escape_space"]),
                loss_conductive=float(d["loss_conductive"]),
                loss_condensation=float(d.get("loss_condensation", 0.0)),
            )
        )
    return rows


def frontier(rows: list[SweepRow]) -> list[FrontierPoint]:
    """Build the `e_eff(rho)` frontier (ascending in rho) with its loss decomposition."""
    points: list[FrontierPoint] = []
    for r in sorted(rows, key=lambda row: row.rho_impact):
        losses = (
            r.loss_radiative_wall,
            r.loss_escape_space,
            r.loss_conductive,
            r.loss_condensation,
        )
        total = sum(losses)
        fa, fb, fc, fd = tuple(x / total for x in losses) if total > 0.0 else (0.0, 0.0, 0.0, 0.0)
        points.append(
            FrontierPoint(
                rho_impact=r.rho_impact,
                e_eff=r.e_eff,
                total_loss=total,
                frac_radiative_wall=fa,
                frac_escape_space=fb,
                frac_conductive=fc,
                frac_condensation=fd,
            )
        )
    return points


def _write_csv(header: list[str], rows: Iterable[Iterable[object]], path: Path) -> None:
    """Write a CSV with `header` then `rows`, creating the parent directory. Shared by the rho and
    velocity frontier writers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def write_summary(points: list[FrontierPoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write the frontier to a CSV (header = the `FrontierPoint` field names), ascending in rho."""
    header = [f.name for f in fields(FrontierPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def plot_frontier(
    points: list[FrontierPoint], out_dir: Path = DEFAULT_PLOT_DIR, tag: str = ""
) -> list[Path]:
    """Render `e_eff` vs rho and the stacked loss decomposition. Returns the saved figure paths.

    `tag` prefixes the filenames so the low-v figures don't overwrite the high-v ones. matplotlib is
    imported lazily (the `sci` extra) on the headless `Agg` backend so `make analysis` needs no
    display. The numeric inputs are plain Python floats built above, so no `Any` from the untyped
    matplotlib API escapes this function.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    rho = [p.rho_impact for p in points]
    e_eff = [p.e_eff for p in points]
    saved: list[Path] = []

    # --- e_eff(rho) frontier ---
    fig, ax = plt.subplots()
    ax.plot(rho, e_eff, "o-")
    ax.set_xlabel(r"$\rho_\mathrm{impact}$ [kg/m$^3$]")
    ax.set_ylabel(r"$e_\mathrm{eff}$")
    ax.set_title(r"Restitution frontier $e_\mathrm{eff}(\rho)$")
    ax.grid(True, alpha=0.3)
    frontier_path = out_dir / f"{tag}e_eff_frontier.png"
    fig.savefig(frontier_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(frontier_path)

    # --- stacked loss decomposition (fractions of total loss) ---
    fa = [p.frac_radiative_wall for p in points]
    fb = [p.frac_escape_space for p in points]
    fc = [p.frac_conductive for p in points]
    fd = [p.frac_condensation for p in points]
    base_b = list(fa)
    base_c = [a + b for a, b in zip(fa, fb, strict=True)]
    base_d = [a + b + c for a, b, c in zip(fa, fb, fc, strict=True)]
    fig, ax = plt.subplots()
    width = 0.6 * min((b - a for a, b in pairwise(rho)), default=1.0)
    ax.bar(rho, fa, width=width, label="1a radiative→wall")
    ax.bar(rho, fb, width=width, bottom=base_b, label="1b escape→space")
    ax.bar(rho, fc, width=width, bottom=base_c, label="2 conductive")
    ax.bar(rho, fd, width=width, bottom=base_d, label="3 condensation")
    ax.set_xlabel(r"$\rho_\mathrm{impact}$ [kg/m$^3$]")
    ax.set_ylabel("fraction of total loss")
    ax.set_title("Loss-channel decomposition (ADR-0016)")
    ax.legend(loc="upper right", fontsize="small")
    losses_path = out_dir / f"{tag}loss_decomposition.png"
    fig.savefig(losses_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(losses_path)

    return saved


@dataclass(frozen=True)
class TransitionalPoint:
    """One velocity on the transitional `e_eff(v)` frontier (ADR-0012). `e_eff_eos`/`e_eff_rad` are
    the rho-mean restitutions of the EOS-only and radiation-on sweeps; `rad_band = e_eff_eos -
    e_eff_rad` is the radiative-uncertainty band; `e_eff_eos_min`/`_max` bracket the spread over the
    swept densities (the shaded band on the figure)."""

    v: float
    e_eff_eos: float
    e_eff_rad: float
    rad_band: float
    e_eff_eos_min: float
    e_eff_eos_max: float


def _mean_e_eff_by_v(rows: list[SweepRow]) -> dict[float, list[float]]:
    """Group `e_eff` by impact speed `v` (the transitional sweep keys; each `v` has one row per
    swept density)."""
    groups: dict[float, list[float]] = {}
    for r in rows:
        groups.setdefault(r.v, []).append(r.e_eff)
    return groups


def transitional_frontier(
    eos_rows: list[SweepRow], rad_rows: list[SweepRow]
) -> list[TransitionalPoint]:
    """Build the rho-mean `e_eff_eos(v)` and `e_eff_rad(v)` curves (ascending in v) from the two
    transitional sweeps over the same `v x rho` grid. A velocity present in the EOS sweep but absent
    from the radiation sweep gets `nan` for the radiation-on value (the band is then undefined)."""
    eos_by_v = _mean_e_eff_by_v(eos_rows)
    rad_by_v = _mean_e_eff_by_v(rad_rows)
    points: list[TransitionalPoint] = []
    for v in sorted(eos_by_v):
        es = eos_by_v[v]
        rs = rad_by_v.get(v, [])
        e_eos = sum(es) / len(es)
        e_rad = sum(rs) / len(rs) if rs else float("nan")
        points.append(
            TransitionalPoint(
                v=v,
                e_eff_eos=e_eos,
                e_eff_rad=e_rad,
                rad_band=e_eos - e_rad,
                e_eff_eos_min=min(es),
                e_eff_eos_max=max(es),
            )
        )
    return points


def locate_dip(points: list[TransitionalPoint]) -> TransitionalPoint | None:
    """Return the interior velocity of minimum EOS-only `e_eff` if it is a genuine dip — strictly
    below *both* swept endpoints — else `None` (the floor sits at an endpoint, so any transitional
    dip is purely radiative and needs the deferred real opacity table). `points` must be ascending
    in v, as `transitional_frontier` returns."""
    if len(points) < 3:
        return None
    dip = min(points[1:-1], key=lambda p: p.e_eff_eos)
    if dip.e_eff_eos < points[0].e_eff_eos and dip.e_eff_eos < points[-1].e_eff_eos:
        return dip
    return None


def write_transitional_summary(
    points: list[TransitionalPoint], path: Path = DEFAULT_TRANS_SUMMARY_PATH
) -> None:
    """Write the `e_eff(v)` frontier to a CSV (header = the `TransitionalPoint` field names),
    ascending in v."""
    header = [f.name for f in fields(TransitionalPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def plot_transitional(
    points: list[TransitionalPoint],
    dip: TransitionalPoint | None,
    out_dir: Path = DEFAULT_PLOT_DIR,
    tag: str = "",
) -> list[Path]:
    """Render the `e_eff(v)` overlay: the EOS-only and radiation-on curves with the EOS spread over
    rho shaded, the prior-rung anchors marked, and the located dip (if any) annotated. Returns the
    saved figure path. matplotlib is imported lazily (the `sci` extra) on the headless `Agg`
    backend; all numeric inputs are plain Python floats, so no `Any` escapes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    v = [p.v / 1000.0 for p in points]  # km/s
    e_eos = [p.e_eff_eos for p in points]
    e_rad = [p.e_eff_rad for p in points]
    lo = [p.e_eff_eos_min for p in points]
    hi = [p.e_eff_eos_max for p in points]

    fig, ax = plt.subplots()
    ax.fill_between(v, lo, hi, alpha=0.15, color="C0", label=r"EOS-only spread over $\rho$")
    ax.plot(v, e_eos, "o-", color="C0", label=r"$e_\mathrm{eff}$ EOS-only ($\rho$-mean)")
    ax.plot(v, e_rad, "s--", color="C1", label=r"$e_\mathrm{eff}$ radiation-on (interim opacity)")
    ax.scatter(
        [ANCHOR_LOWV[0] / 1000.0],
        [ANCHOR_LOWV[1]],
        marker="*",
        s=140,
        color="k",
        zorder=5,
        label="3.2 km/s anchor (Rung C)",
    )
    ax.scatter(
        [ANCHOR_HIGHV[0] / 1000.0],
        [ANCHOR_HIGHV[1]],
        marker="D",
        s=55,
        color="k",
        zorder=5,
        label="16 km/s anchor (Rung B)",
    )
    if dip is not None:
        ax.axvline(dip.v / 1000.0, color="C3", ls=":", alpha=0.7)
        ax.annotate(
            f"EOS dip {dip.e_eff_eos:.3f}\n@ {dip.v / 1000.0:.0f} km/s",
            xy=(dip.v / 1000.0, dip.e_eff_eos),
            xytext=(6, 12),
            textcoords="offset points",
            fontsize="small",
            color="C3",
        )
    ax.set_xlabel(r"impact speed $v$ [km/s]")
    ax.set_ylabel(r"$e_\mathrm{eff}$")
    ax.set_title(r"Transitional restitution $e_\mathrm{eff}(v)$ (ADR-0012)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize="small")
    path = out_dir / f"{tag}transitional_e_eff_v.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _run_transitional(
    eos_path: Path, rad_path: Path, summary_path: Path, plot_dir: Path, tag: str
) -> None:
    """The `--axis v` path: read the two transitional sweeps, build the `e_eff(v)` frontier, locate
    the EOS dip, write the CSV, and render the overlay figure."""
    points = transitional_frontier(read_sweep(eos_path), read_sweep(rad_path))
    dip = locate_dip(points)
    write_transitional_summary(points, summary_path)
    figs = plot_transitional(points, dip, plot_dir, tag)
    if dip is not None:
        print(
            f"python: EOS-only dip e_eff={dip.e_eff_eos:.4f} at v={dip.v / 1000.0:.1f} km/s "
            f"(radiative band {dip.rad_band:+.4f}); below both swept endpoints."
        )
    else:
        print("python: no interior EOS dip — the e_eff(v) floor is at a swept endpoint.")
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


@dataclass(frozen=True)
class GeoRow:
    """One geometry-sweep result row (`crates/sweep --geometry`, ADR-0023): the case and its
    `eta_capture` with the two restitution ratios it was formed from."""

    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    mach: float
    eta_capture: float
    restitution_free: float
    restitution_confined: float
    peak_force: float
    peak_local_pressure: float  # peak local plate pressure of the free run (Rung S focusing factor)


@dataclass(frozen=True)
class GeometryPoint:
    """One reconciled operating point: the geometry case, its `eta_capture`, the axial column
    density it implies (the `Sigma` contract), and `f = eta_capture·(1+e_eff)/2` at the two
    `e_eff` anchors."""

    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    mach: float
    eta_capture: float
    # Axial column density per unit cloud density, `Sigma/rho = L = 2·(L/D)·r_foot` (`r_foot = 1`).
    # The `Sigma = m/(pi r_foot^2)` contract (ADR-0003) is set by `L/D` alone — the footprint
    # cancels for a uniform cylinder, so `r_foot/R` is purely the `eta_capture` lever, not a
    # `Sigma` knob.
    sigma_over_rho: float
    f_dip: float  # f at the transitional worst-case e_eff (the conservative floor)
    f_highv: float  # f at the 16 km/s e_eff


def read_geometry(path: Path = DEFAULT_GEOMETRY_PATH) -> list[GeoRow]:
    """Parse the geometry sweep JSONL (one JSON object per line; blank lines tolerated)."""
    rows: list[GeoRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            GeoRow(
                d_over_d=float(d["d_over_d"]),
                l_over_d=float(d["l_over_d"]),
                r_foot_over_r=float(d["r_foot_over_r"]),
                mach=float(d["mach"]),
                eta_capture=float(d["eta_capture"]),
                restitution_free=float(d["restitution_free"]),
                restitution_confined=float(d["restitution_confined"]),
                peak_force=float(d["peak_force"]),
                peak_local_pressure=float(d.get("peak_local_pressure", 0.0)),
            )
        )
    return rows


def reconcile_f(eta_capture: float, e_eff: float) -> float:
    """The paper's fudge factor `f = eta_capture·(1 + e_eff)/2` (ADR-0003)."""
    return eta_capture * (1.0 + e_eff) / 2.0


def geometry_frontier(rows: list[GeoRow]) -> list[GeometryPoint]:
    """Reconcile each geometry case into `f` at the two `e_eff` anchors, sorted by
    `(mach, L/D, r_foot/R, d/D)` for a stable CSV/plot order."""
    points: list[GeometryPoint] = []
    for r in sorted(rows, key=lambda x: (x.mach, x.l_over_d, x.r_foot_over_r, x.d_over_d)):
        points.append(
            GeometryPoint(
                d_over_d=r.d_over_d,
                l_over_d=r.l_over_d,
                r_foot_over_r=r.r_foot_over_r,
                mach=r.mach,
                eta_capture=r.eta_capture,
                sigma_over_rho=2.0 * r.l_over_d,
                f_dip=reconcile_f(r.eta_capture, EEFF_DIP),
                f_highv=reconcile_f(r.eta_capture, EEFF_HIGHV),
            )
        )
    return points


def write_geometry_summary(
    points: list[GeometryPoint], path: Path = DEFAULT_GEOMETRY_SUMMARY_PATH
) -> None:
    """Write the reconciled geometry frontier to a CSV (header = `GeometryPoint` field names)."""
    header = [f.name for f in fields(GeometryPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def plot_geometry(
    points: list[GeometryPoint], out_dir: Path = DEFAULT_PLOT_DIR, tag: str = ""
) -> list[Path]:
    """Render `eta_capture` and `f` vs footprint `r_foot/R`, one curve per curvature `d/D`, for a
    representative `(Mach, L/D)` slice — the flat plate is the floor, the shallow-concave plates the
    recovery lever. `f` uses the transitional worst-case `e_eff` (the conservative floor), with the
    useful-`f` gate marked. matplotlib is imported lazily (`sci` extra); all inputs are plain
    floats, so no `Any` escapes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    # A representative slice: the higher Mach anchor and the median L/D.
    machs = sorted({p.mach for p in points})
    lds = sorted({p.l_over_d for p in points})
    sel_mach = machs[-1]
    sel_ld = lds[len(lds) // 2]
    sliced = [p for p in points if p.mach == sel_mach and p.l_over_d == sel_ld]
    curvatures = sorted({p.d_over_d for p in sliced})

    fig, (ax_eta, ax_f) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
    for dd in curvatures:
        row = sorted((p for p in sliced if p.d_over_d == dd), key=lambda p: p.r_foot_over_r)
        rf = [p.r_foot_over_r for p in row]
        label = "flat" if dd == 0.0 else f"d/D = {dd:.2f}"
        ax_eta.plot(rf, [p.eta_capture for p in row], "o-", label=label)
        ax_f.plot(rf, [p.f_dip for p in row], "o-", label=label)
    ax_eta.axhline(1.0, color="grey", ls=":", alpha=0.6, label="1D limit (flat ceiling)")
    ax_eta.set_ylabel(r"$\eta_\mathrm{capture}$")
    ax_eta.set_title(
        rf"Geometry sweep: $M={sel_mach:.0f}$, $L/D={sel_ld:.2f}$ "
        r"(concave lifts $\eta$ over the flat floor)"
    )
    ax_eta.grid(True, alpha=0.3)
    ax_eta.legend(loc="best", fontsize="small")
    ax_f.axhline(
        USEFUL_F_GATE, color="C3", ls="--", alpha=0.7, label=rf"useful-$f$ gate {USEFUL_F_GATE}"
    )
    ax_f.set_xlabel(r"footprint $r_\mathrm{foot}/R$")
    ax_f.set_ylabel(rf"$f = \eta\,(1+e_\mathrm{{eff}})/2$, $e_\mathrm{{eff}}={EEFF_DIP}$")
    ax_f.grid(True, alpha=0.3)
    ax_f.legend(loc="best", fontsize="small")
    path = out_dir / f"{tag}geometry_eta_f.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _run_geometry(geometry_path: Path, summary_path: Path, plot_dir: Path, tag: str) -> None:
    """The `--axis geometry` path: read the geometry sweep, reconcile `f`, write the CSV, render the
    figure, and report the flat-floor vs best-concave `f` bracket."""
    points = geometry_frontier(read_geometry(geometry_path))
    write_geometry_summary(points, summary_path)
    figs = plot_geometry(points, plot_dir, tag)
    flat = [p for p in points if p.d_over_d == 0.0]
    concave = [p for p in points if p.d_over_d > 0.0]
    if flat and concave:
        best_flat = max(p.f_dip for p in flat)
        best_concave = max(p.f_dip for p in concave)
        best_concave_hi = max(p.f_highv for p in concave)
        print(
            f"python: f (dip e_eff={EEFF_DIP}) — flat floor up to {best_flat:.3f}, "
            f"concave up to {best_concave:.3f}; at 16 km/s (e_eff={EEFF_HIGHV}) concave up to "
            f"{best_concave_hi:.3f} (gate {USEFUL_F_GATE})."
        )
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


# --- Rung S: survivability frontier primitives (design §7, ADR-0010/0011) ---


@dataclass(frozen=True)
class SurvivabilityVerdict:
    """A peak facesheet load classified against the SiC+Ti structural limits (ADR-0010/0011)."""

    peak_compressive: float  # peak facesheet pressure (Pa)
    reflected_tensile: float  # reflected tension at the SiC-Ti interface (Pa)
    survives_compressive: bool
    survives_spall: bool


def peak_facesheet_pressure(rho: float, v: float, c_stag: float) -> float:
    """Peak facesheet pressure = the stagnation pressure `c_stag·rho·v²` (design §7): a cold
    coasting cloud's ram pressure recompressed at the wall."""
    return c_stag * rho * v * v


def stagnation_coefficient(rows: list[SweepRow], v: float, *, tol: float = 1.0) -> float:
    """Back the stagnation coefficient `c_stag = peak_wall_force/(rho·v²)` out of the 1D sweep,
    averaged over the densities at velocity `v` — so survivability uses the kernel's measured number
    rather than an assumed ideal (`peak ≈ 2.0·rho·v²` at 16 km/s)."""
    coeffs = [
        r.peak_wall_force / (r.rho_impact * r.v * r.v)
        for r in rows
        if abs(r.v - v) <= tol and r.rho_impact > 0.0
    ]
    if not coeffs:
        raise ValueError(f"no sweep rows at v={v}")
    return sum(coeffs) / len(coeffs)


def impact_density(
    l_over_d: float,
    r_foot_over_r: float,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
) -> float:
    """The Σ contract (ADR-0003) made physical: `Σ = m/(π r_foot²) = rho·L`, with
    `L = 2·(L/D)·r_foot` and `r_foot = (r_foot/R)·R`, so `rho = m / (2π·(L/D)·r_foot³)`. A disk
    (small L/D) is dense (high peak pressure); a cylinder dilute; a tighter footprint is denser
    (∝ 1/(r_foot/R)³)."""
    r_foot = r_foot_over_r * plate_radius
    return mass / (2.0 * math.pi * l_over_d * r_foot * r_foot * r_foot)


def density_ceiling(v: float, c_stag: float, p_limit: float) -> float:
    """The densest cloud whose peak pressure stays under `p_limit`: inverts the pressure law to
    `rho = p_limit/(c_stag·v²)`."""
    return p_limit / (c_stag * v * v)


def reflected_tensile(peak_compressive: float) -> float:
    """The tensile wave reflected into the brittle SiC at the SiC-Ti impedance step: `|R|·peak`,
    `R ~ -0.15` (ADR-0011)."""
    return REFLECT_FRAC * peak_compressive


def classify_survivability(
    peak_compressive: float, p_limit: float, spall_strength: float = SIC_SPALL_LO
) -> SurvivabilityVerdict:
    """Classify a peak facesheet load against the two structural limits (ADR-0010/0011): the
    compressive `p_limit` and the reflected-tensile SiC spall strength."""
    tensile = reflected_tensile(peak_compressive)
    return SurvivabilityVerdict(
        peak_compressive=peak_compressive,
        reflected_tensile=tensile,
        survives_compressive=peak_compressive < p_limit,
        survives_spall=tensile < spall_strength,
    )


@dataclass(frozen=True)
class SurvivabilityPoint:
    """One geometry case resolved to physical survivability at a velocity anchor: the Σ-bridge
    density, the peak facesheet pressure it implies, its `f`, and whether it clears the baseline /
    relaxed structural limits (ADR-0010/0011)."""

    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    mach: float
    v: float
    e_eff: float
    eta_capture: float
    rho_impact: float
    focusing_factor: float  # concave local-peak concentration over the flat reference (Rung S)
    peak_compressive: float  # the plane-wave stagnation peak scaled by focusing_factor
    reflected_tensile: float
    f: float
    survives_baseline: bool  # peak < P_LIMIT_BASELINE (400 MPa) and spall OK
    survives_relaxed: bool  # peak < max(P_LIMIT_HIGHV) (900 MPa) and spall OK


def survivability_frontier(
    rows: list[GeoRow],
    anchors: Iterable[tuple[float, float, float]],
    *,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
    spall_strength: float = SIC_SPALL_LO,
) -> list[SurvivabilityPoint]:
    """Resolve each geometry case to physical survivability at each `(v, e_eff, c_stag)` anchor: the
    Σ contract gives `rho`, the stagnation law the plane-wave peak, the concave focusing factor (the
    case's local peak over its flat `d/D=0` counterpart at the same `L/D, r_foot/R, mach`) the local
    concentration, and `classify_survivability` the verdict against the 400 MPa baseline and the
    relaxed 900 MPa high-v limit (design §7)."""
    relaxed_limit = max(P_LIMIT_HIGHV)
    flat_local = {
        (r.l_over_d, r.r_foot_over_r, r.mach): r.peak_local_pressure
        for r in rows
        if r.d_over_d == 0.0 and r.peak_local_pressure > 0.0
    }
    points: list[SurvivabilityPoint] = []
    for v, e_eff, c_stag in anchors:
        for r in sorted(rows, key=lambda x: (x.mach, x.l_over_d, x.r_foot_over_r, x.d_over_d)):
            rho = impact_density(r.l_over_d, r.r_foot_over_r, mass, plate_radius)
            ref = flat_local.get((r.l_over_d, r.r_foot_over_r, r.mach))
            focusing = r.peak_local_pressure / ref if ref else 1.0
            peak = peak_facesheet_pressure(rho, v, c_stag) * focusing
            base = classify_survivability(peak, P_LIMIT_BASELINE, spall_strength)
            relaxed = classify_survivability(peak, relaxed_limit, spall_strength)
            points.append(
                SurvivabilityPoint(
                    d_over_d=r.d_over_d,
                    l_over_d=r.l_over_d,
                    r_foot_over_r=r.r_foot_over_r,
                    mach=r.mach,
                    v=v,
                    e_eff=e_eff,
                    eta_capture=r.eta_capture,
                    rho_impact=rho,
                    focusing_factor=focusing,
                    peak_compressive=peak,
                    reflected_tensile=base.reflected_tensile,
                    f=reconcile_f(r.eta_capture, e_eff),
                    survives_baseline=base.survives_compressive and base.survives_spall,
                    survives_relaxed=relaxed.survives_compressive and relaxed.survives_spall,
                )
            )
    return points


def best_survivable_f(points: list[SurvivabilityPoint], *, relaxed: bool = False) -> float | None:
    """The highest `f` among the cases that clear the limit (baseline by default, relaxed 900 MPa if
    `relaxed`). `None` if nothing survives — the answer to 'what is the best survivable f?'."""
    survivors = [p for p in points if (p.survives_relaxed if relaxed else p.survives_baseline)]
    return max((p.f for p in survivors), default=None)


def write_survivability_summary(
    points: list[SurvivabilityPoint], path: Path = DEFAULT_SURVIVABILITY_SUMMARY_PATH
) -> None:
    """Write the survivability frontier to a CSV (header = `SurvivabilityPoint` field names)."""
    header = [f.name for f in fields(SurvivabilityPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


@dataclass(frozen=True)
class MarginPoint:
    """The best survivable `f` at a scaled plate radius `R` and pulse mass `m`, per velocity anchor:
    one cell of the closed-form margin map. `headroom` is the pressure-ceiling relief vs the fixed
    baseline `(R=5 m, m=25 kg)`: `(R/R0)**3 * (m0/m)`, the exact `1/rho` factor by which every
    shape's peak `2*rho*v**2` drops, since `rho` scales as `m/R**3` (the Sigma contract). A higher
    headroom admits denser shapes and lifts the best survivable `f`."""

    v: float
    plate_radius: float
    mass: float
    headroom: float
    best_f_baseline: float | None  # best survivable f at the 400 MPa baseline
    best_f_relaxed: float | None  # best survivable f at the relaxed 900 MPa high-v limit


def margin_map(
    rows: list[GeoRow],
    anchors: Iterable[tuple[float, float, float]],
    *,
    plate_radii: Iterable[float] = MARGIN_PLATE_RADII_M,
    masses: Iterable[float] = MARGIN_MASSES_KG,
    base_radius: float = PLATE_RADIUS_M,
    base_mass: float = PULSE_MASS_KG,
) -> list[MarginPoint]:
    """The closed-form `f`-margin map: best survivable `f` over the `(plate radius R, pulse mass m)`
    grid, per velocity anchor. A pure rescaling of `survivability_frontier`: only `impact_density`
    (which scales as `m/R**3`) changes between cells, so the geometry/anchor data is reused with no
    kernel reruns. This is the f-SIDE of a system trade; the cost-side (plate-mass and pulse-count)
    is out of scope. It quantifies margin above the already-passing `f ~ 0.8` (de-risk, not a
    necessity). See design §7 / ADR-0010 (amendment)."""
    anchors = list(anchors)
    points: list[MarginPoint] = []
    for radius in plate_radii:
        for mass in masses:
            headroom = (radius / base_radius) ** 3 * (base_mass / mass)
            resolved = survivability_frontier(rows, anchors, mass=mass, plate_radius=radius)
            for v in sorted({p.v for p in resolved}):
                at_v = [p for p in resolved if p.v == v]
                points.append(
                    MarginPoint(
                        v=v,
                        plate_radius=radius,
                        mass=mass,
                        headroom=headroom,
                        best_f_baseline=best_survivable_f(at_v, relaxed=False),
                        best_f_relaxed=best_survivable_f(at_v, relaxed=True),
                    )
                )
    return points


def write_margin_summary(
    points: list[MarginPoint], path: Path = DEFAULT_MARGIN_SUMMARY_PATH
) -> None:
    """Write the closed-form margin map to a CSV (header = `MarginPoint` field names)."""
    header = [f.name for f in fields(MarginPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def plot_survivability(
    points: list[SurvivabilityPoint], out_dir: Path = DEFAULT_PLOT_DIR, tag: str = ""
) -> list[Path]:
    """Render the survivability story at the highest-v anchor: peak facesheet pressure vs `L/D` (one
    curve per footprint, with the 400/700/900 MPa `P_limit` lines) over `f` vs `L/D` with survivable
    points filled and non-survivable open. The high-`f` short-disk / tight-footprint corner sits
    above the pressure limit — the elongated, wider-footprint shapes are the survivable ones. The
    deepest curvature is shown for `f`; peak pressure is curvature-independent (the Σ bridge). Lazy
    matplotlib (`sci` extra); all inputs are plain floats/bools, so no `Any` escapes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    sel_v = max(p.v for p in points)
    sel_mach = max(p.mach for p in points if p.v == sel_v)
    sel_dd = max(p.d_over_d for p in points)
    sliced = [p for p in points if p.v == sel_v and p.mach == sel_mach and p.d_over_d == sel_dd]
    footprints = sorted({p.r_foot_over_r for p in sliced})

    fig, (ax_p, ax_f) = plt.subplots(2, 1, sharex=True, figsize=(6, 7))
    for rf in footprints:
        row = sorted((p for p in sliced if p.r_foot_over_r == rf), key=lambda p: p.l_over_d)
        ld = [p.l_over_d for p in row]
        ax_p.plot(ld, [p.peak_compressive / 1e6 for p in row], "o-", label=rf"$r_f/R={rf:.1f}$")
        surv = [p.survives_baseline for p in row]
        fvals = [p.f for p in row]
        ax_f.plot(ld, fvals, "-", color="grey", alpha=0.4)
        ax_f.scatter(
            ld,
            fvals,
            marker="o",
            facecolors=["C0" if s else "none" for s in surv],
            edgecolors="C0",
        )
    for plim, style in ((400.0, "--"), (700.0, ":"), (900.0, "-.")):
        ax_p.axhline(plim, color="C3", ls=style, alpha=0.7, label=f"{plim:.0f} MPa")
    ax_p.set_yscale("log")
    ax_p.set_ylabel("peak facesheet pressure (MPa)")
    ax_p.set_title(
        rf"Survivability at $v={sel_v / 1000:.0f}$ km/s "
        r"(high-$f$ disk/tight-footprint fails on pressure)"
    )
    ax_p.grid(True, alpha=0.3, which="both")
    ax_p.legend(loc="best", fontsize="small", ncol=2)
    ax_f.axhline(USEFUL_F_GATE, color="C3", ls="--", alpha=0.7, label=r"useful-$f$ gate")
    ax_f.set_xlabel(r"$L/D$")
    ax_f.set_ylabel(rf"$f$ ($e_\mathrm{{eff}}={sliced[0].e_eff}$); filled = survives 400 MPa")
    ax_f.grid(True, alpha=0.3)
    ax_f.legend(loc="best", fontsize="small")
    path = out_dir / f"{tag}survivability.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _run_survivability(
    geometry_path: Path,
    highv_sweep_path: Path,
    dip_sweep_path: Path,
    summary_path: Path,
    plot_dir: Path,
    tag: str,
) -> None:
    """The `--axis survivability` path: back `c_stag` out of the 1D sweeps at each anchor, resolve
    every geometry case to a peak facesheet pressure via the Σ contract, classify against the SiC+Ti
    limits, write the CSV + figure, and report the best survivable `f` per anchor."""
    c_stag_hi = stagnation_coefficient(read_sweep(highv_sweep_path), ANCHOR_HIGHV[0])
    c_stag_dip = stagnation_coefficient(read_sweep(dip_sweep_path), V_DIP)
    anchors = [(V_DIP, EEFF_DIP, c_stag_dip), (ANCHOR_HIGHV[0], EEFF_HIGHV, c_stag_hi)]
    points = survivability_frontier(read_geometry(geometry_path), anchors)
    write_survivability_summary(points, summary_path)
    figs = plot_survivability(points, plot_dir, tag)

    hi = [p for p in points if p.v == ANCHOR_HIGHV[0]]
    dip = [p for p in points if p.v == V_DIP]

    def _fmt(x: float | None) -> str:
        return f"{x:.3f}" if x is not None else "none survive"

    best_dip = best_survivable_f(dip)
    best_hi_base = best_survivable_f(hi)
    best_hi_relax = best_survivable_f(hi, relaxed=True)
    worst = max(points, key=lambda p: p.peak_compressive)
    print(
        f"python: best survivable f — dip ({V_DIP / 1000:.0f} km/s, 400 MPa) {_fmt(best_dip)}; "
        f"16 km/s baseline 400 MPa {_fmt(best_hi_base)}, relaxed 900 MPa {_fmt(best_hi_relax)}."
    )
    print(
        f"python: the f-max corner (L/D={worst.l_over_d}, r_foot/R={worst.r_foot_over_r}) peaks at "
        f"{worst.peak_compressive / 1e6:.0f} MPa at {worst.v / 1000:.0f} km/s — "
        f"{'survives' if worst.survives_relaxed else 'fails'} the relaxed limit."
    )
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


def _run_margin(
    geometry_path: Path,
    highv_sweep_path: Path,
    dip_sweep_path: Path,
    summary_path: Path,
) -> None:
    """The `--axis margin` path: the closed-form `f`-margin exploration (design §7, ADR-0010
    amendment). Reuses the survivability anchors (c_stag from the 1D sweeps) and rescales the
    frontier over the `(R, m)` grid — no kernel reruns — to report how much survivable `f` a wider
    plate or a smaller pulse buys above the already-passing baseline."""
    c_stag_hi = stagnation_coefficient(read_sweep(highv_sweep_path), ANCHOR_HIGHV[0])
    c_stag_dip = stagnation_coefficient(read_sweep(dip_sweep_path), V_DIP)
    anchors = [(V_DIP, EEFF_DIP, c_stag_dip), (ANCHOR_HIGHV[0], EEFF_HIGHV, c_stag_hi)]
    points = margin_map(read_geometry(geometry_path), anchors)
    write_margin_summary(points, summary_path)

    def _fmt(x: float | None) -> str:
        return f"{x:.3f}" if x is not None else "  none "

    for v in sorted({p.v for p in points}):
        print(f"python: f-margin map at {v / 1000:.0f} km/s (best survivable f, 400 MPa baseline):")
        print("python:   R\\m (kg) |" + "".join(f"  {m:>5.0f}" for m in MARGIN_MASSES_KG))
        for radius in MARGIN_PLATE_RADII_M:
            cells = "".join(
                _fmt(
                    next(
                        p.best_f_baseline
                        for p in points
                        if p.v == v and p.plate_radius == radius and p.mass == m
                    )
                )
                for m in MARGIN_MASSES_KG
            )
            print(f"python:   R={radius:>4.1f} m | {cells}")
    base = next(
        p
        for p in points
        if p.plate_radius == PLATE_RADIUS_M and p.mass == PULSE_MASS_KG and p.v == ANCHOR_HIGHV[0]
    )
    best = max(points, key=lambda p: (p.v == ANCHOR_HIGHV[0], p.best_f_baseline or 0.0))
    print(
        f"python: 16 km/s baseline rises from f={_fmt(base.best_f_baseline)} (R=5 m, m=25 kg) to "
        f"f={_fmt(best.best_f_baseline)} (R={best.plate_radius:.1f} m, m={best.mass:.0f} kg, "
        f"headroom {best.headroom:.1f}x). R and m are pinned by external budgets — this is the "
        f"f-side of a system trade. Wrote {summary_path}."
    )


def main() -> None:
    """Read a sweep, extract the frontier, write the CSV summary, and render the figures. `--axis
    rho` (default) builds the `e_eff(rho)` frontier + loss decomposition for the high-v/low-v
    sweeps; `--axis v` builds the transitional `e_eff(v)` overlay (ADR-0012) from the two
    transitional sweeps; `--axis geometry` reconciles `f = eta_capture·(1+e_eff)/2` from the
    geometry sweep (ADR-0003, Rung D follow-on); `--axis survivability` resolves each geometry case
    to a peak facesheet pressure via the Σ contract and reports the best survivable `f` per anchor
    against the SiC+Ti limits (Rung S, ADR-0010/0011)."""
    parser = argparse.ArgumentParser(description="e_eff frontier extraction + loss decomposition.")
    parser.add_argument(
        "--axis",
        choices=["rho", "v", "geometry", "survivability", "margin"],
        default="rho",
        help="frontier axis",
    )
    parser.add_argument(
        "--sweep", type=Path, default=None, help="input sweep JSONL (EOS-only file when --axis v)"
    )
    parser.add_argument(
        "--sweep-rad",
        type=Path,
        default=DEFAULT_TRANS_RAD_PATH,
        help="radiation-on transitional sweep JSONL (--axis v only)",
    )
    parser.add_argument("--summary", type=Path, default=None, help="output CSV")
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR, help="figure directory")
    parser.add_argument("--tag", default="", help="figure filename prefix (e.g. 'lowv_')")
    args = parser.parse_args()
    plot_dir: Path = args.plot_dir
    tag: str = args.tag

    if args.axis == "v":
        eos_path: Path = args.sweep or DEFAULT_TRANS_EOS_PATH
        rad_path: Path = args.sweep_rad
        summary_path: Path = args.summary or DEFAULT_TRANS_SUMMARY_PATH
        _run_transitional(eos_path, rad_path, summary_path, plot_dir, tag)
        return

    if args.axis == "geometry":
        geometry_path: Path = args.sweep or DEFAULT_GEOMETRY_PATH
        summary_path = args.summary or DEFAULT_GEOMETRY_SUMMARY_PATH
        _run_geometry(geometry_path, summary_path, plot_dir, tag)
        return

    if args.axis == "survivability":
        geometry_path = args.sweep or DEFAULT_GEOMETRY_PATH
        summary_path = args.summary or DEFAULT_SURVIVABILITY_SUMMARY_PATH
        _run_survivability(
            geometry_path,
            DEFAULT_SWEEP_PATH,
            DEFAULT_TRANS_EOS_PATH,
            summary_path,
            plot_dir,
            tag,
        )
        return

    if args.axis == "margin":
        geometry_path = args.sweep or DEFAULT_GEOMETRY_PATH
        summary_path = args.summary or DEFAULT_MARGIN_SUMMARY_PATH
        _run_margin(geometry_path, DEFAULT_SWEEP_PATH, DEFAULT_TRANS_EOS_PATH, summary_path)
        return

    sweep_path: Path = args.sweep or DEFAULT_SWEEP_PATH
    summary_path = args.summary or DEFAULT_SUMMARY_PATH
    rows = read_sweep(sweep_path)
    points = frontier(rows)
    write_summary(points, summary_path)
    figs = plot_frontier(points, plot_dir, tag)
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


if __name__ == "__main__":
    main()
