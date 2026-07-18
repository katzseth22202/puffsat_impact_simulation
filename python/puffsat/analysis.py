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
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from _typeshed import DataclassInstance

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

# Ablating-wall recovery sweep (Rung E, ADR-0014): the `crates/sweep --ablating` output and its
# tau-bracketed e_eff recovery + the 16 km/s f-gate call.
DEFAULT_ABLATING_PATH = Path("data/results/sweep_ablating.jsonl")
DEFAULT_ABLATING_SUMMARY_PATH = Path("data/results/frontier_ablating.csv")

# Frozen-recombination bounding sweep (audit finding 3): the `crates/sweep --frozen` output and its
# freeze-timing bracket around the equilibrium e_eff(v) curve.
DEFAULT_FROZEN_PATH = Path("data/results/sweep_frozen.jsonl")
DEFAULT_FROZEN_SUMMARY_PATH = Path("data/results/frontier_frozen.csv")

# 1D `e_eff` anchors for the `f = eta_capture·(1+e_eff)/2` reconciliation (ADR-0003). `EEFF_DIP` is
# the transitional EOS worst case (the conservative floor, ADR-0012); `EEFF_HIGHV` the 16 km/s point
# (Rung B). `eta_capture` is geometry-dominated and ~velocity-independent, so the first bracket
# pairs it with these two `e_eff` scenarios; the full `Sigma`-resolved `e_eff(rho(r_foot))` lookup
# is the deferred refinement (the dual-curve `f(v)` deliverable, ADR-0013).
EEFF_DIP = 0.57
EEFF_HIGHV = 0.63
USEFUL_F_GATE = 0.8  # the useful-`f` gate (ADR-0009), marked on the figure

# --- Rung S: survivability frontier (design §7, ADR-0010/0011) ---
# Peak facesheet pressure is the stagnation pressure `c_stag·rho·v²`; `c_stag` (≈1.1, the
# reflected-shock (gamma_eff+1)/2) is backed out of the 1D `peak_wall_pressure` (physical EOS p,
# AV excluded — ADR-0010 correction: the earlier peak_wall_force-based ≈2.0 was the artificial-
# viscosity spike) so the frontier inherits the kernel's measured number, not an assumed ideal.
# The structural limits are the SiC+Ti facesheet's: a compressive `P_limit`
# (the §5 band), a reflected-tensile spall limit at the SiC-Ti interface (ADR-0011), where
# `|R| ~ 0.15` of the incident compressive returns as tension, and a Ti back-face free-surface
# spall limit (ADR-0011 amendment), where the compression transmitted through the SiC-Ti step
# reflects off the solid Ti layer's back surface as tension.
PULSE_MASS_KG = 25.0  # gas delivered per PuffSat (design §2)
PLATE_RADIUS_M = 5.0  # plate radius R, fixed (design §2/§7)
V_DIP = 11_000.0  # transitional worst-case velocity (~11 km/s, ADR-0012), paired with EEFF_DIP
P_LIMIT_BASELINE = 400.0e6  # conservative floor of the §5 SiC+Ti band (Pa)
P_LIMIT_HIGHV = (700.0e6, 900.0e6)  # relaxed limits swept at the 16 km/s anchor (design §7)
REFLECT_FRAC = 0.15  # |R| at the SiC-Ti impedance step -> reflected tensile fraction (ADR-0011)
SIC_SPALL_LO = 0.3e9  # SiC dynamic spall strength, conservative end (ADR-0011)
SIC_SPALL_HI = 1.0e9  # SiC dynamic spall strength, upper end
# Ti back-face free-surface spall (ADR-0011 amendment). The compression transmitted through the
# SiC-Ti step (stress convention `T = 1 + R = 1 - |R| ~ 0.85` of the incident) reaches the solid Ti
# layer's back face and reflects there as tension (`R ~ -1`: a free — or, per ADR-0011's "no voids
# behind the SiC" corollary, low-impedance-terminated — back surface), so the peak Ti back-face
# tension is `~0.85*peak`. Checked against Ti dynamic spall strength (~2.5-4.5 GPa; ductile, ~8x the
# brittle SiC), this gate is looser than the SiC-interface one and never controls the frontier — the
# SiC spalls first despite seeing far less tension (0.15 vs 0.85 of the peak).
TI_TRANSMIT_FRAC = 1.0 - REFLECT_FRAC  # compressive fraction transmitted into Ti (T = 1 + R)
TI_SPALL_LO = 2.5e9  # Ti-6Al-4V / CP-Ti dynamic spall strength, conservative end
TI_SPALL_HI = 4.5e9  # Ti dynamic spall strength, upper end

# Physical ceiling for eta_capture: shallow-concave over-collimation reaches ~1.01 (Rung D-cc), so
# anything past this is a solver blow-up (one M=40 case once returned 7.6), not physics. Shared by
# the jupiter/heavyplate special-scenario geometry gates.
ETA_PHYSICAL_MAX = 1.2

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
    #: Physical peak wall pressure (EOS p only, artificial viscosity excluded) — the facesheet
    #: survivability load (ADR-0010 correction). 0.0 when reading pre-fix JSONL (stale data).
    peak_wall_pressure: float
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


def read_jsonl_rows[T: DataclassInstance](
    cls: type[T],
    path: Path,
    *,
    str_fields: tuple[str, ...] = (),
    defaults: dict[str, float] | None = None,
) -> list[T]:
    """Parse a sweep-result JSONL into instances of the frozen dataclass `cls` (one JSON object
    per line; blank lines tolerated) — the read + skip-blank + `json.loads` + per-field cast shape
    every `read_*_sweep` in this codebase shares.

    Every field is read by name and cast to `float`, except those named in `str_fields` (cast to
    `str` — the shape sweeps' `axis`/`sigma_role` columns). `defaults` supplies a fallback for
    fields that may be absent from older JSONL (back-compat, e.g. `peak_wall_pressure` before
    ADR-0010); a field with neither a value in the row nor a listed default raises `KeyError`, so
    genuine schema drift still fails loudly rather than silently defaulting.
    """
    defaults = defaults or {}
    rows: list[T] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        kwargs: dict[str, object] = {}
        for f in fields(cls):
            raw = d.get(f.name, defaults[f.name]) if f.name in defaults else d[f.name]
            kwargs[f.name] = str(raw) if f.name in str_fields else float(raw)
        rows.append(cls(**kwargs))
    return rows


def read_sweep(path: Path = DEFAULT_SWEEP_PATH) -> list[SweepRow]:
    """Parse the JSONL sweep results (one JSON object per line; blank lines tolerated)."""
    return read_jsonl_rows(
        SweepRow,
        path,
        defaults={
            "peak_wall_force": 0.0,
            "peak_wall_pressure": 0.0,
            "loss_condensation": 0.0,
        },
    )


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
    return read_jsonl_rows(GeoRow, path, defaults={"peak_local_pressure": 0.0})


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
    back_face_tensile: float  # tension at the solid Ti layer's back face (Pa)
    survives_compressive: bool
    survives_spall: bool  # SiC-interface reflected-tensile spall (brittle, binds first)
    survives_back_spall: bool  # Ti back-face free-surface spall (ductile, looser gate)


def peak_facesheet_pressure(rho: float, v: float, c_stag: float) -> float:
    """Peak facesheet pressure = the stagnation pressure `c_stag·rho·v²` (design §7): a cold
    coasting cloud's ram pressure recompressed at the wall."""
    return c_stag * rho * v * v


def stagnation_coefficient(rows: list[SweepRow], v: float, *, tol: float = 1.0) -> float:
    """Back the stagnation coefficient `c_stag = peak_wall_pressure/(rho·v²)` out of the 1D sweep,
    averaged over the densities at velocity `v` — the kernel's measured **physical** peak (EOS `p`
    of the wall cell, artificial viscosity excluded; ≈(gamma_eff+1)/2 ≈ 1.1). The earlier
    `peak_wall_force`-based value ≈ 2.0 was the AV first-impact spike, a numerical artifact
    (ADR-0010 correction, 2026-07)."""
    coeffs = [
        r.peak_wall_pressure / (r.rho_impact * r.v * r.v)
        for r in rows
        if abs(r.v - v) <= tol and r.rho_impact > 0.0
    ]
    if not coeffs:
        raise ValueError(f"no sweep rows at v={v}")
    c_stag = sum(coeffs) / len(coeffs)
    if c_stag <= 0.0:
        raise ValueError(
            f"peak_wall_pressure is 0 at v={v}: stale pre-fix JSONL — rerun `make sweep` "
            "(and `make sweep-transitional`) to regenerate with the physical peak recorded"
        )
    return c_stag


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


def plate_mass(radius: float, d_over_d: float, areal_density: float = 45.0) -> float:
    """Plate mass [kg]: areal density x disk area, with the shallow dish's extra-area factor
    `1 + (2 d/D)²` (spherical-cap area `π(a² + d²)`, `d = (d/D)·2R`). `areal_density` defaults to
    45 kg/m² — the baseline stack's central estimate (3-4 t at R = 5 m, design §2), shared by the
    jupiter/heavyplate special-scenario plate sizing."""
    return areal_density * math.pi * radius * radius * (1.0 + (2.0 * d_over_d) ** 2)


class _LogInterp:
    """Piecewise-linear interpolation in `ln x`, clamped at the ends (small, typed, stdlib) — the
    `e_eff(rho)` / `e_eff(v)` slice interpolator shared by the jupiter/heavyplate frontiers."""

    def __init__(self, xs: list[float], ys: list[float]) -> None:
        if len(xs) != len(ys) or len(xs) < 2:
            raise ValueError("need >= 2 matching points")
        self._lx = [math.log(x) for x in xs]
        self._ys = ys

    def __call__(self, x: float) -> float:
        lx = math.log(x)
        if lx <= self._lx[0]:
            return self._ys[0]
        if lx >= self._lx[-1]:
            return self._ys[-1]
        for i in range(1, len(self._lx)):
            if lx <= self._lx[i]:
                t = (lx - self._lx[i - 1]) / (self._lx[i] - self._lx[i - 1])
                return self._ys[i - 1] * (1.0 - t) + self._ys[i] * t
        return self._ys[-1]


def density_ceiling(v: float, c_stag: float, p_limit: float) -> float:
    """The densest cloud whose peak pressure stays under `p_limit`: inverts the pressure law to
    `rho = p_limit/(c_stag·v²)`."""
    return p_limit / (c_stag * v * v)


def reflected_tensile(peak_compressive: float) -> float:
    """The tensile wave reflected into the brittle SiC at the SiC-Ti impedance step: `|R|·peak`,
    `R ~ -0.15` (ADR-0011)."""
    return REFLECT_FRAC * peak_compressive


def back_face_tensile(peak_compressive: float) -> float:
    """The tension at the solid Ti layer's back face (ADR-0011 amendment): the compression
    transmitted through the SiC-Ti step (`(1-|R|)·peak`, stress convention `T = 1 + R`) reflects at
    the free / low-impedance-terminated back surface as tension of the same amplitude (`R ~ -1`)."""
    return TI_TRANSMIT_FRAC * peak_compressive


def classify_survivability(
    peak_compressive: float,
    p_limit: float,
    spall_strength: float = SIC_SPALL_LO,
    ti_spall_strength: float = TI_SPALL_LO,
) -> SurvivabilityVerdict:
    """Classify a peak facesheet load against the structural limits (ADR-0010/0011): the compressive
    `p_limit`, the reflected-tensile SiC spall strength at the SiC-Ti interface, and the Ti
    back-face free-surface spall strength. The brittle SiC-interface spall binds first (low
    strength); the ductile Ti back-face check is the looser confirmatory gate."""
    tensile = reflected_tensile(peak_compressive)
    back_tensile = back_face_tensile(peak_compressive)
    return SurvivabilityVerdict(
        peak_compressive=peak_compressive,
        reflected_tensile=tensile,
        back_face_tensile=back_tensile,
        survives_compressive=peak_compressive < p_limit,
        survives_spall=tensile < spall_strength,
        survives_back_spall=back_tensile < ti_spall_strength,
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
    reflected_tensile: float  # SiC-Ti interface reflected tension
    back_face_tensile: float  # Ti back-face free-surface tension
    f: float
    survives_baseline: bool  # peak < P_LIMIT_BASELINE (400 MPa) and both spall checks OK
    survives_relaxed: bool  # peak < max(P_LIMIT_HIGHV) (900 MPa) and both spall checks OK


def survivability_frontier(
    rows: list[GeoRow],
    anchors: Iterable[tuple[float, float, float]],
    *,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
    spall_strength: float = SIC_SPALL_LO,
    ti_spall_strength: float = TI_SPALL_LO,
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
            base = classify_survivability(peak, P_LIMIT_BASELINE, spall_strength, ti_spall_strength)
            relaxed = classify_survivability(peak, relaxed_limit, spall_strength, ti_spall_strength)
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
                    back_face_tensile=base.back_face_tensile,
                    f=reconcile_f(r.eta_capture, e_eff),
                    survives_baseline=base.survives_compressive
                    and base.survives_spall
                    and base.survives_back_spall,
                    survives_relaxed=relaxed.survives_compressive
                    and relaxed.survives_spall
                    and relaxed.survives_back_spall,
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


# --- Rung E: ablating-wall recovery (the tau-bracket + the 16 km/s f-gate call, ADR-0014) ---


@dataclass(frozen=True)
class AblatingRow:
    """One ablating-sweep result row (`crates/sweep --ablating`, ADR-0014): the case axes plus the
    rigid floor, the ablating restitution, and the ablation bookkeeping."""

    v: float
    rho_impact: float
    opacity_scale: float
    q_star: float
    kappa_vapor: float
    e_eff_rigid: float
    e_eff_ablating: float
    recovery: float
    ablated_mass: float
    ablated_fraction: float
    loss_radiative_wall: float
    loss_escape_space: float
    loss_ablation: float
    peak_wall_force: float


def read_ablating(path: Path = DEFAULT_ABLATING_PATH) -> list[AblatingRow]:
    """Parse the ablating sweep JSONL (one JSON object per line; blank lines tolerated)."""
    return read_jsonl_rows(AblatingRow, path)


@dataclass(frozen=True)
class AblatingPoint:
    """The rho-mean ablating recovery at one `(v, opacity_scale, Q*)` case: the rigid floor, the
    ablating best estimate, and the recovery, averaged over the impact-density grid (the
    single-anchor convention the geometry/survivability reconciliation already uses, ADR-0013)."""

    v: float
    opacity_scale: float
    q_star: float
    e_eff_rigid: float
    e_eff_ablating: float
    recovery: float
    ablated_fraction: float


def _mean(values: list[float]) -> float:
    """Arithmetic mean of a non-empty list."""
    return sum(values) / len(values)


def ablating_points(rows: list[AblatingRow]) -> list[AblatingPoint]:
    """Collapse the impact-density axis: rho-mean the rigid floor, the ablating `e_eff`, and the
    recovery at each `(v, opacity_scale, Q*)`, sorted for a stable CSV/plot order."""
    keys = sorted({(r.v, r.opacity_scale, r.q_star) for r in rows})
    points: list[AblatingPoint] = []
    for v, scale, q_star in keys:
        grp = [r for r in rows if r.v == v and r.opacity_scale == scale and r.q_star == q_star]
        points.append(
            AblatingPoint(
                v=v,
                opacity_scale=scale,
                q_star=q_star,
                e_eff_rigid=_mean([r.e_eff_rigid for r in grp]),
                e_eff_ablating=_mean([r.e_eff_ablating for r in grp]),
                recovery=_mean([r.recovery for r in grp]),
                ablated_fraction=_mean([r.ablated_fraction for r in grp]),
            )
        )
    return points


def write_ablating_summary(
    points: list[AblatingPoint], path: Path = DEFAULT_ABLATING_SUMMARY_PATH
) -> None:
    """Write the rho-mean ablating recovery to a CSV (header = `AblatingPoint` field names)."""
    header = [f.name for f in fields(AblatingPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def best_f_at(
    geo: list[GeoRow], v: float, e_eff: float, c_stag: float, *, relaxed: bool = False
) -> float | None:
    """Best survivable `f` at velocity `v` with restitution `e_eff` — resolve the geometry cases to
    survivability at the `(v, e_eff, c_stag)` anchor and take the survivable maximum (Rung S)."""
    points = survivability_frontier(geo, [(v, e_eff, c_stag)])
    return best_survivable_f(points, relaxed=relaxed)


def plot_ablating(
    points: list[AblatingPoint], out_dir: Path = DEFAULT_PLOT_DIR, tag: str = ""
) -> list[Path]:
    """Render the ablating recovery vs the opacity scale (the tau-bracket), one curve per `Q*`, one
    panel per velocity anchor: recovery grows as the gas un-traps (scale drops) and radiation
    reaches the vapor curtain (ADR-0012). matplotlib is lazy (`sci` extra); inputs are floats."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    velocities = sorted({p.v for p in points})
    q_stars = sorted({p.q_star for p in points})
    fig, axes = plt.subplots(1, len(velocities), figsize=(5 * len(velocities), 4.5), squeeze=False)
    for ax, v in zip(axes[0], velocities, strict=True):
        for q_star in q_stars:
            row = sorted(
                (p for p in points if p.v == v and p.q_star == q_star),
                key=lambda p: p.opacity_scale,
            )
            ax.plot(
                [p.opacity_scale for p in row],
                [p.recovery for p in row],
                "o-",
                label=rf"$Q^*={q_star / 1e6:.0f}$ MJ/kg",
            )
        ax.set_xscale("log")
        ax.axhline(0.0, color="grey", ls=":", alpha=0.6)
        ax.set_xlabel("opacity scale (1 = interim Kramers)")
        ax.set_ylabel(r"$\Delta e_\mathrm{eff}$ (ablating $-$ rigid floor)")
        ax.set_title(rf"$v={v / 1000:.0f}$ km/s")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize="small")
    fig.suptitle("Ablating-wall recovery vs optical depth (the tau-bracket, ADR-0014)")
    fig.tight_layout()
    path = out_dir / f"{tag}ablating_recovery.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return [path]


def _run_ablating(
    ablating_path: Path,
    geometry_path: Path,
    highv_sweep_path: Path,
    dip_sweep_path: Path,
    summary_path: Path,
    plot_dir: Path,
    tag: str,
) -> None:
    """The `--axis ablating` path: read the ablating sweep, rho-mean the recovery, write the CSV +
    figure, and report (1) the dip-fill as a tau-bracket against the EOS floor and (2) the 16 km/s
    `f >= 0.8`-at-a-survivable-shape call (the user-deferred recovery-lever decision, ADR-0009)."""
    rows = read_ablating(ablating_path)
    points = ablating_points(rows)
    write_ablating_summary(points, summary_path)
    figs = plot_ablating(points, plot_dir, tag)

    geo = read_geometry(geometry_path)
    c_stag_hi = stagnation_coefficient(read_sweep(highv_sweep_path), ANCHOR_HIGHV[0])
    c_stag_dip = stagnation_coefficient(read_sweep(dip_sweep_path), V_DIP)

    def _fmt(x: float | None) -> str:
        return f"{x:.3f}" if x is not None else "none survive"

    # Dip fill (the EOS worst case is not radiatively driven, so the shield has little to recover).
    dip = [p for p in points if p.v == V_DIP]
    dip_rigid = _mean([p.e_eff_rigid for p in dip])
    dip_abl = [p.e_eff_ablating for p in dip]
    print(
        f"python: dip ({V_DIP / 1000:.0f} km/s) — EOS floor e_eff={EEFF_DIP}; radiation-on rigid "
        f"{dip_rigid:.3f}; ablating recovers to [{min(dip_abl):.3f}, {max(dip_abl):.3f}] "
        f"(recovery [{min(p.recovery for p in dip):+.4f}, {max(p.recovery for p in dip):+.4f}]). "
        "The dip is EOS-dominated, not radiatively fillable."
    )

    # 16 km/s f-gate call: the ablating e_eff bracket over (scale, Q*) vs the rigid floor. The
    # bracket's optimistic end is the low-Q*/high-tau corner (most ablation, most shielding); its
    # conservative end is high-Q*/low-tau. Report which end clears the gate, with the plate ablation
    # the clearing corner costs — clearing at the optimistic end alone is a marginal, Q*-dependent
    # clear, not a robust one.
    hi = [p for p in points if p.v == ANCHOR_HIGHV[0]]
    e_rigid = _mean([p.e_eff_rigid for p in hi])
    p_lo = min(hi, key=lambda p: p.e_eff_ablating)
    p_hi = max(hi, key=lambda p: p.e_eff_ablating)
    e_abl_lo, e_abl_hi = p_lo.e_eff_ablating, p_hi.e_eff_ablating
    f_rigid = best_f_at(geo, ANCHOR_HIGHV[0], e_rigid, c_stag_hi)
    f_abl_lo = best_f_at(geo, ANCHOR_HIGHV[0], e_abl_lo, c_stag_hi)
    f_abl_hi = best_f_at(geo, ANCHOR_HIGHV[0], e_abl_hi, c_stag_hi)
    f_abl_hi_relaxed = best_f_at(geo, ANCHOR_HIGHV[0], e_abl_hi, c_stag_hi, relaxed=True)
    clears_lo = f_abl_lo is not None and f_abl_lo >= USEFUL_F_GATE
    clears_hi = f_abl_hi is not None and f_abl_hi >= USEFUL_F_GATE
    if clears_lo:
        verdict = "CLEARS across the bracket (robust)"
    elif clears_hi:
        pct = p_hi.ablated_fraction * 100
        q_mj = p_hi.q_star / 1e6
        verdict = (
            f"STRADDLES the gate — clears only at the optimistic end (f {_fmt(f_abl_hi)}, "
            f"ablating {pct:.1f}% of the plate at Q*={q_mj:.0f} MJ/kg); "
            f"the conservative end lands {_fmt(f_abl_lo)}, just under"
        )
    else:
        verdict = "does NOT clear at either end of the bracket"
    print(
        f"python: 16 km/s — rigid floor e_eff={e_rigid:.3f} (best survivable f {_fmt(f_rigid)}) -> "
        f"ablating e_eff [{e_abl_lo:.3f}, {e_abl_hi:.3f}] (best survivable f "
        f"[{_fmt(f_abl_lo)}, {_fmt(f_abl_hi)}]). Gate {USEFUL_F_GATE} at the 400 MPa baseline: "
        f"{verdict}. Relaxed 900 MPa: {_fmt(f_abl_hi_relaxed)}."
    )
    # The dip f for completeness (ablating barely moves it, so survivability still binds).
    f_dip_abl = best_f_at(geo, V_DIP, max(dip_abl), c_stag_dip)
    print(
        f"python: dip best survivable f at the ablating e_eff {_fmt(f_dip_abl)} "
        f"(gate {USEFUL_F_GATE})."
    )
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


# ---- Frozen-recombination bounding sweep (audit finding 3): freeze-timing brackets --------------


@dataclass(frozen=True)
class FrozenRow:
    """One frozen-sweep result row (`crates/sweep --frozen`): the three e_eff curves at one case."""

    v: float
    rho_impact: float
    e_eff_eq: float
    e_eff_frozen_rebound: float
    e_eff_frozen_all: float
    rho_star: float
    t_star: float
    swap_energy_jump_frac: float


@dataclass(frozen=True)
class FrozenPoint:
    """One `e_eff(v)` point with its freeze-timing bracket (rho-means): the equilibrium curve, the
    sudden-freeze-at-turnaround pessimistic bound, and the pure-H2O no-chemistry optimistic bound.
    `delta_frozen = e_eff_eq - e_eff_frozen_rebound >= 0` is the chemistry-return content of the
    equilibrium rebound — the amount at risk if recombination freezes."""

    v: float
    e_eff_eq: float
    e_eff_frozen_rebound: float
    e_eff_frozen_all: float
    delta_frozen: float
    e_eff_frozen_min: float
    e_eff_frozen_max: float


@dataclass(frozen=True)
class FrozenDipImpact:
    """The worst-case (dip) points of the equilibrium and frozen-rebound curves and the implied
    upper bound on the `f` shift, `delta_f_max = eta*(dip_eq - dip_frozen)/2` at `eta = 1`."""

    dip_eq: FrozenPoint
    dip_frozen: FrozenPoint
    delta_f_max: float


def read_frozen(path: Path = DEFAULT_FROZEN_PATH) -> list[FrozenRow]:
    """Parse the frozen-sweep JSONL (one JSON object per line; blank lines tolerated)."""
    return read_jsonl_rows(FrozenRow, path)


def frozen_frontier(rows: list[FrozenRow]) -> list[FrozenPoint]:
    """Rho-mean the three curves per velocity (ascending in v), carrying the frozen-rebound spread
    over rho as `[min, max]`."""
    groups: dict[float, list[FrozenRow]] = {}
    for r in rows:
        groups.setdefault(r.v, []).append(r)
    points: list[FrozenPoint] = []
    for v in sorted(groups):
        rs = groups[v]
        e_eq = _mean([r.e_eff_eq for r in rs])
        e_frozen = _mean([r.e_eff_frozen_rebound for r in rs])
        points.append(
            FrozenPoint(
                v=v,
                e_eff_eq=e_eq,
                e_eff_frozen_rebound=e_frozen,
                e_eff_frozen_all=_mean([r.e_eff_frozen_all for r in rs]),
                delta_frozen=e_eq - e_frozen,
                e_eff_frozen_min=min(r.e_eff_frozen_rebound for r in rs),
                e_eff_frozen_max=max(r.e_eff_frozen_rebound for r in rs),
            )
        )
    return points


def frozen_dip_impact(points: list[FrozenPoint]) -> FrozenDipImpact:
    """Locate each curve's minimum over the swept velocities and bound the `f` impact of freezing:
    `f = eta·(1 + e_eff)/2`, so the dip-to-dip `e_eff` drop costs at most `eta·delta_e/2 ≤
    delta_e/2` of `f` (`eta ≈ 1` at the concave dip operating point)."""
    dip_eq = min(points, key=lambda p: p.e_eff_eq)
    dip_frozen = min(points, key=lambda p: p.e_eff_frozen_rebound)
    delta = dip_eq.e_eff_eq - dip_frozen.e_eff_frozen_rebound
    return FrozenDipImpact(dip_eq=dip_eq, dip_frozen=dip_frozen, delta_f_max=delta / 2.0)


def write_frozen_summary(
    points: list[FrozenPoint], path: Path = DEFAULT_FROZEN_SUMMARY_PATH
) -> None:
    """Write the freeze-timing bracket to a CSV (header = the `FrozenPoint` field names)."""
    header = [f.name for f in fields(FrozenPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def plot_frozen(
    points: list[FrozenPoint],
    impact: FrozenDipImpact,
    out_dir: Path = DEFAULT_PLOT_DIR,
    tag: str = "",
) -> list[Path]:
    """Render the freeze-timing bracket around the equilibrium `e_eff(v)` curve: pure-H2O
    no-chemistry above, sudden-freeze-at-turnaround below (with its rho-spread shaded), the two
    dips annotated. matplotlib is imported lazily (the `sci` extra) on the headless Agg backend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    v = [p.v / 1000.0 for p in points]
    fig, ax = plt.subplots()
    ax.fill_between(
        v,
        [p.e_eff_frozen_min for p in points],
        [p.e_eff_frozen_max for p in points],
        alpha=0.15,
        color="C3",
        label=r"frozen-rebound spread over $\rho$",
    )
    ax.plot(
        v,
        [p.e_eff_frozen_all for p in points],
        "^--",
        color="C2",
        label="frozen throughout (no chemistry — freeze before plate)",
    )
    ax.plot(
        v,
        [p.e_eff_eq for p in points],
        "o-",
        color="C0",
        label="equilibrium (study curve)",
    )
    ax.plot(
        v,
        [p.e_eff_frozen_rebound for p in points],
        "v-",
        color="C3",
        label="sudden freeze at turnaround (freeze after plate)",
    )
    dip = impact.dip_frozen
    ax.annotate(
        f"frozen dip {dip.e_eff_frozen_rebound:.3f}\n@ {dip.v / 1000.0:.0f} km/s",
        xy=(dip.v / 1000.0, dip.e_eff_frozen_rebound),
        xytext=(6, -24),
        textcoords="offset points",
        fontsize="small",
        color="C3",
    )
    ax.set_xlabel(r"impact speed $v$ [km/s]")
    ax.set_ylabel(r"$e_\mathrm{eff}$")
    ax.set_title(r"Freeze-timing bracket on $e_\mathrm{eff}(v)$ (EOS-only)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize="small")
    path = out_dir / f"{tag}frozen_e_eff_v.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _run_frozen(frozen_path: Path, summary_path: Path, plot_dir: Path, tag: str) -> None:
    """The `--axis frozen` path: read the frozen sweep, build the freeze-timing bracket, write the
    CSV, render the overlay, and report the worst-case dip shift and its `f` impact bound."""
    rows = read_frozen(frozen_path)
    points = frozen_frontier(rows)
    impact = frozen_dip_impact(points)
    write_frozen_summary(points, summary_path)
    figs = plot_frozen(points, impact, plot_dir, tag)
    jump_max = max(abs(r.swap_energy_jump_frac) for r in rows)
    print(
        f"python: equilibrium dip e_eff={impact.dip_eq.e_eff_eq:.4f} "
        f"@ {impact.dip_eq.v / 1000.0:.0f} km/s; sudden-freeze dip "
        f"e_eff={impact.dip_frozen.e_eff_frozen_rebound:.4f} "
        f"@ {impact.dip_frozen.v / 1000.0:.0f} km/s "
        f"(delta_e={impact.dip_eq.e_eff_eq - impact.dip_frozen.e_eff_frozen_rebound:+.4f})."
    )
    print(
        f"python: worst-case f impact of frozen recombination <= {impact.delta_f_max:.4f} "
        f"(delta_f = eta*delta_e/2 at eta=1); splice energy-jump diagnostic max "
        f"{jump_max:.2e} of incident KE."
    )
    print(f"python: wrote {summary_path} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


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
        choices=["rho", "v", "geometry", "survivability", "margin", "ablating", "frozen"],
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

    if args.axis == "ablating":
        ablating_path: Path = args.sweep or DEFAULT_ABLATING_PATH
        summary_path = args.summary or DEFAULT_ABLATING_SUMMARY_PATH
        _run_ablating(
            ablating_path,
            DEFAULT_GEOMETRY_PATH,
            DEFAULT_SWEEP_PATH,
            DEFAULT_TRANS_EOS_PATH,
            summary_path,
            plot_dir,
            tag,
        )
        return

    if args.axis == "frozen":
        frozen_path: Path = args.sweep or DEFAULT_FROZEN_PATH
        summary_path = args.summary or DEFAULT_FROZEN_SUMMARY_PATH
        _run_frozen(frozen_path, summary_path, plot_dir, tag)
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
