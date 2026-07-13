"""Heavy-plate 16-28 km/s special-scenario analysis: `f(v)` + facesheet survivability at a pinned
30 m / ≤ 40 t plate (design §12.1, ADR-0027).

The scenario scales the vehicle up from the core envelope: a **100 kg** PuffSat pulse on a
**tripled** pusher plate (`R = 15 m`, 30 m diameter) of mass **≤ 40 t**, swept **16-28 km/s at
0.5 km/s** (25 anchors; 16 km/s overlaps the core anchor as a consistency check). The question
(design §12.1): does the fudge factor stay near `f ≈ 0.8`, and does the facesheet survive, at this
larger, faster, heavier configuration?

Unlike `jupiter.py` — which pins the velocity and sweeps the plate radius — here **`R` and `m` are
pinned** and the **velocity** is the swept axis. The Sigma contract (ADR-0003)
`rho = m/(2pi·(L/D)·r_foot³)` with `r_foot = (r_foot/R)·R` is therefore velocity-independent: each
geometry shape has one impact density, and `f(v)` moves only through `e_eff(v, rho)` (and the
survivability peak `∝ v²`).

Pipeline (reuses the core machinery, design §12.1):

- **`e_eff(v, rho)`** from the `--heavyplate` coupled sweep on the reused Jupiter extended-grid
  table (real TOPS/OPLIB opacity), equilibrium headline at a fixed representative stretched-cloud
  length, interpolated per velocity over `rho` (log-rho). The `L`-sensitivity and opacity tau-checks
  are reported as diagnostics (`tau >> 1` => both flat).
- **`eta_capture`** from the `M = 40` geometry sweep (scale-invariant, Mach-converged), flat +
  shallow-concave, with the concave local-focusing factor from `peak_local_pressure`.
- **Peak facesheet pressure** `c_stag·rho·v²·focusing` (ADR-0010), classified against the SiC+Ti
  ladder (400 baseline / 900 relaxed MPa) and the SiC-Ti spall limit.
- **Deliverable** `f = eta_capture·(1 + e_eff)/2` as a `f(v)` dual curve (flat + best concave),
  with the `f = 0.8` reference line (ADR-0009, a reference anchor, not a mission gate).

The whole-plate structural go/no-go (rigid-during-pulse, membrane/bending, spall) is the closed-form
companion in `puffsat.structure` (ADR-0027), decoupled from this `f(v)` frontier.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, fields
from pathlib import Path

from puffsat.analysis import (
    P_LIMIT_BASELINE,
    P_LIMIT_HIGHV,
    SIC_SPALL_LO,
    GeoRow,
    _write_csv,
    classify_survivability,
    impact_density,
    peak_facesheet_pressure,
    read_geometry,
    reconcile_f,
)

DEFAULT_HEAVYPLATE_SWEEP_PATH = Path("data/results/sweep_heavyplate.jsonl")
DEFAULT_GEOMETRY_M40_PATH = Path("data/results/sweep_geometry_m40.jsonl")
DEFAULT_GEOMETRY_M20_PATH = Path("data/results/sweep_geometry.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/frontier_heavyplate.csv")
DEFAULT_PLOT_DIR = Path("data/results")
# The 16-28 km/s freeze-timing bracket (ADR-0026 instrument): the `--frozen-heavyplate` EOS-only
# sweep and the frontier that translates its e_eff delta onto the headline survivable f.
DEFAULT_FROZEN_SWEEP_PATH = Path("data/results/sweep_frozen_heavyplate.jsonl")
DEFAULT_FROZEN_SUMMARY_PATH = Path("data/results/frontier_frozen_heavyplate.csv")

PULSE_MASS_KG = 100.0  # gas per shot (design §12.1: the `m` in `rho ∝ m/R³`)
PLATE_RADIUS_M = 15.0  # tripled radius, 30 m diameter (design §12.1)
PLATE_MASS_CEILING_KG = 40_000.0  # ≤ 40 t ceiling (design §12.1)

# Fixed representative stretched-cloud length [m] the headline `e_eff` slice is read at — must match
# the Rust `HEAVY_LENGTH`. The `L`-spot rows (other lengths) are diagnostics, not the headline.
LENGTH_ANCHOR = 10.0
# The three freeze-timing / L-sensitivity bracket anchors [m/s] (design §12.1) — must match the Rust
# `HEAVY_V_ANCHORS`.
V_ANCHORS = (16_000.0, 22_000.0, 28_000.0)
# The opacity τ-check velocity [m/s] (design §12.1) — must match the Rust `HEAVY_TAU_V`.
V_TAU_CHECK = 28_000.0

# Plate areal density [kg/m²]: baseline stack 3-4 t at R = 5 m (design §2) ⇒ 38-51; 45 central. At
# R = 15 m this is ~32 t flat (`45·π·15²`), inside the ≤ 40 t ceiling with margin.
AREAL_DENSITY = 45.0
AREAL_DENSITY_BAND = (38.0, 51.0)

# Physical ceiling for eta_capture: shallow-concave over-collimation reaches ~1.01 (Rung D-cc), so
# anything past this is a solver blow-up, not physics (a stray M = 40 case once returned 7.6).
ETA_PHYSICAL_MAX = 1.2

FLOAT_TOL = 1e-6  # grid-value match tolerance (the sweep writes exact round grid values)


@dataclass(frozen=True)
class HeavyPlateRow:
    """One heavy-plate sweep row (the fields this analysis needs from the JSONL schema)."""

    v: float
    rho_impact: float
    length: float
    opacity_scale: float
    e_eff: float
    peak_wall_pressure: float
    incident_momentum: float
    wall_impulse: float
    loss_radiative_wall: float
    loss_escape_space: float


def read_heavyplate_sweep(path: Path = DEFAULT_HEAVYPLATE_SWEEP_PATH) -> list[HeavyPlateRow]:
    """Parse the `--heavyplate` sweep JSONL (one JSON object per line; blank lines tolerated)."""
    rows: list[HeavyPlateRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(HeavyPlateRow(**{f.name: float(d[f.name]) for f in fields(HeavyPlateRow)}))
    return rows


class _LogInterp:
    """Piecewise-linear interpolation in `ln x`, clamped at the ends (small, typed, stdlib)."""

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


def headline_rows(rows: list[HeavyPlateRow]) -> list[HeavyPlateRow]:
    """The headline slice: the representative length at real opacity (`κ = 1`) — the `e_eff(v, rho)`
    grid the frontier is built from (L-spot and τ-check rows are diagnostics, excluded here)."""
    return [
        r
        for r in rows
        if abs(r.length - LENGTH_ANCHOR) <= FLOAT_TOL and abs(r.opacity_scale - 1.0) <= FLOAT_TOL
    ]


def sweep_velocities(rows: list[HeavyPlateRow]) -> list[float]:
    """The distinct swept velocities present in the headline slice (ascending) — the `f(v)` axis."""
    return sorted({r.v for r in headline_rows(rows)})


def e_eff_interpolator_at_v(rows: list[HeavyPlateRow], v: float) -> _LogInterp:
    """`e_eff(rho)` along the headline slice at velocity `v`, log-rho linear interpolation."""
    slice_rows = sorted(
        (r for r in headline_rows(rows) if abs(r.v - v) <= FLOAT_TOL),
        key=lambda r: r.rho_impact,
    )
    if not slice_rows:
        raise ValueError(f"no headline sweep rows at v={v}")
    return _LogInterp([r.rho_impact for r in slice_rows], [r.e_eff for r in slice_rows])


def stagnation_coefficient_at_v(rows: list[HeavyPlateRow], v: float) -> float:
    """`c_stag = peak_wall_pressure/(rhov²)` averaged over the headline rho slice at velocity `v`
    (ADR-0010: the physical EOS peak, AV excluded). Weakly v-dependent (~1.1-1.2)."""
    coeffs = [
        r.peak_wall_pressure / (r.rho_impact * v * v)
        for r in headline_rows(rows)
        if abs(r.v - v) <= FLOAT_TOL and r.rho_impact > 0.0
    ]
    if not coeffs:
        raise ValueError(f"no headline sweep rows at v={v}")
    c = sum(coeffs) / len(coeffs)
    if c <= 0.0:
        raise ValueError("non-positive c_stag — stale or empty sweep JSONL?")
    return c


def plate_mass(radius: float, d_over_d: float, areal_density: float = AREAL_DENSITY) -> float:
    """Plate mass [kg]: areal density x disk area, with the shallow dish's extra-area factor
    `1 + (2 d/D)²` (spherical-cap area `π(a² + d²)`, `d = (d/D)·2R`)."""
    return areal_density * math.pi * radius * radius * (1.0 + (2.0 * d_over_d) ** 2)


@dataclass(frozen=True)
class HeavyPlatePoint:
    """One (velocity x cloud shape) case resolved to survivability and `f` at the pinned plate."""

    v: float
    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    rho_impact: float
    e_eff: float
    eta_capture: float
    focusing_factor: float
    plate_mass_t: float  # tonnes at AREAL_DENSITY (band scales linearly); should be ≤ 40 t
    peak_compressive: float
    f: float
    survives_baseline: bool
    survives_relaxed: bool
    within_mass_budget: bool


def heavyplate_frontier(
    geo_rows: list[GeoRow],
    sweep_rows: list[HeavyPlateRow],
    velocities: list[float],
    *,
    mass: float = PULSE_MASS_KG,
    plate_radius: float = PLATE_RADIUS_M,
) -> list[HeavyPlatePoint]:
    """Resolve every (velocity x geometry case) to survivability and `f` at the pinned 30 m plate.

    Per case: the Sigma contract gives the (velocity-independent) `rho(shape)`; the per-velocity
    headline slice gives `e_eff(rho)`; the geometry row gives `eta_capture` and the focusing factor;
    the stagnation law + focusing gives the peak, classified against the 400 MPa baseline and the
    relaxed 900 MPa limit; the plate mass is checked against the ≤ 40 t ceiling."""
    relaxed_limit = max(P_LIMIT_HIGHV)
    flat_local = {
        (r.l_over_d, r.r_foot_over_r): r.peak_local_pressure
        for r in geo_rows
        if r.d_over_d == 0.0 and r.peak_local_pressure > 0.0
    }
    points: list[HeavyPlatePoint] = []
    for v in velocities:
        e_of_rho = e_eff_interpolator_at_v(sweep_rows, v)
        c_stag = stagnation_coefficient_at_v(sweep_rows, v)
        for r in sorted(geo_rows, key=lambda x: (x.l_over_d, x.r_foot_over_r, x.d_over_d)):
            rho = impact_density(r.l_over_d, r.r_foot_over_r, mass, plate_radius)
            ref = flat_local.get((r.l_over_d, r.r_foot_over_r))
            focusing = r.peak_local_pressure / ref if ref else 1.0
            peak = peak_facesheet_pressure(rho, v, c_stag) * focusing
            base = classify_survivability(peak, P_LIMIT_BASELINE, SIC_SPALL_LO)
            relaxed = classify_survivability(peak, relaxed_limit, SIC_SPALL_LO)
            mass_kg = plate_mass(plate_radius, r.d_over_d)
            e_eff = e_of_rho(rho)
            points.append(
                HeavyPlatePoint(
                    v=v,
                    d_over_d=r.d_over_d,
                    l_over_d=r.l_over_d,
                    r_foot_over_r=r.r_foot_over_r,
                    rho_impact=rho,
                    e_eff=e_eff,
                    eta_capture=r.eta_capture,
                    focusing_factor=focusing,
                    plate_mass_t=mass_kg / 1000.0,
                    peak_compressive=peak,
                    f=reconcile_f(r.eta_capture, e_eff),
                    survives_baseline=base.survives_compressive and base.survives_spall,
                    survives_relaxed=relaxed.survives_compressive and relaxed.survives_spall,
                    within_mass_budget=mass_kg <= PLATE_MASS_CEILING_KG,
                )
            )
    return points


def best_at_v(
    points: list[HeavyPlatePoint], v: float, *, concave: bool | None = None
) -> HeavyPlatePoint | None:
    """The highest-`f` baseline-surviving, within-mass-budget case at velocity `v`; `concave`
    filters d/D > 0 (True) / d/D = 0 (False) / any (None)."""
    candidates = [
        p
        for p in points
        if abs(p.v - v) <= FLOAT_TOL
        and p.survives_baseline
        and p.within_mass_budget
        and (concave is None or (p.d_over_d > 0.0) == concave)
    ]
    return max(candidates, key=lambda p: p.f, default=None)


def write_summary(points: list[HeavyPlatePoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write the heavy-plate frontier to a CSV (header = `HeavyPlatePoint` field names)."""
    header = [f.name for f in fields(HeavyPlatePoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def _e_eff_at(
    rows: list[HeavyPlateRow], v: float, rho: float, length: float, scale: float
) -> float:
    """The `e_eff` of the exact grid row `(v, rho, L, κ)` (the sweep writes exact grid values)."""
    for r in rows:
        if (
            abs(r.v - v) <= FLOAT_TOL
            and abs(r.rho_impact - rho) <= FLOAT_TOL
            and abs(r.length - length) <= FLOAT_TOL
            and abs(r.opacity_scale - scale) <= FLOAT_TOL
        ):
            return r.e_eff
    raise KeyError(f"no sweep row at v={v}, rho={rho}, L={length}, scale={scale}")


def length_sensitivity(rows: list[HeavyPlateRow], rho: float) -> dict[float, float]:
    """Max `e_eff` spread across the available lengths at each anchor velocity (at fixed `rho`, real
    opacity): the design-§12.1 check that `τ ≫ 1` makes `e_eff` `L`-insensitive."""
    lengths = sorted({r.length for r in rows if abs(r.opacity_scale - 1.0) <= FLOAT_TOL})
    spreads: dict[float, float] = {}
    for v in V_ANCHORS:
        vals = []
        for length in lengths:
            try:
                vals.append(_e_eff_at(rows, v, rho, length, 1.0))
            except KeyError:
                continue
        if len(vals) >= 2:
            spreads[v] = max(vals) - min(vals)
    return spreads


def opacity_sensitivity(rows: list[HeavyPlateRow], rho: float) -> float:
    """Max `e_eff` spread across the opacity scales at the τ-check velocity (fixed `rho`, headline
    length): small ⇒ `τ ≫ 1` at the dilute top and the equilibrium headline is opacity-robust."""
    scales = sorted({r.opacity_scale for r in rows if abs(r.v - V_TAU_CHECK) <= FLOAT_TOL})
    vals = []
    for scale in scales:
        try:
            vals.append(_e_eff_at(rows, V_TAU_CHECK, rho, LENGTH_ANCHOR, scale))
        except KeyError:
            continue
    return max(vals) - min(vals) if len(vals) >= 2 else 0.0


def plot_f_of_v(
    points: list[HeavyPlatePoint], out_dir: Path = DEFAULT_PLOT_DIR, tag: str = "heavyplate_"
) -> list[Path]:
    """Render the `f(v)` deliverable: the best-surviving flat floor and best-surviving concave `f`
    versus closing speed, with the `f = 0.8` reference line and the 16 km/s core-study overlap
    marked. matplotlib is imported lazily (`sci` extra); all inputs are plain floats, so no `Any`
    escapes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def best_f(v: float, *, concave: bool) -> float:
        p = best_at_v(points, v, concave=concave)
        return p.f if p is not None else float("nan")

    out_dir.mkdir(parents=True, exist_ok=True)
    velocities = sorted({p.v for p in points})
    v_km = [v / 1000.0 for v in velocities]
    flat_f = [best_f(v, concave=False) for v in velocities]
    conc_f = [best_f(v, concave=True) for v in velocities]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(v_km, flat_f, "o-", label="flat plate (best survivable)")
    ax.plot(v_km, conc_f, "s-", label="shallow concave (best survivable)")
    ax.axhline(0.8, color="C3", ls="--", alpha=0.7, label=r"$f = 0.8$ reference (ADR-0009)")
    ax.axvline(16.0, color="grey", ls=":", alpha=0.6, label="core-study anchor (16 km/s)")
    ax.set_xlabel(r"closing speed $v$ [km/s]")
    ax.set_ylabel(r"$f = \eta_\mathrm{capture}\,(1 + e_\mathrm{eff})/2$")
    ax.set_title(r"Heavy-plate $f(v)$ — 100 kg pulse, 30 m / $\leq$40 t plate")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize="small")
    path = out_dir / f"{tag}f_v.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


# ---- Freeze-timing bracket at 16-28 km/s (ADR-0026 instrument) -----------------------------------
#
# The high-v turnaround is ionized, so freezing the composition locks an unrecoverable
# ionization-energy sink; the sudden-freeze bound sits below equilibrium and the penalty grows with
# v. The EOS-only `--frozen-heavyplate` sweep (at the 16 / 22 / 28 km/s anchors) measures that
# delta; here it is translated onto the coupled headline f at each anchor's design point.


@dataclass(frozen=True)
class FrozenHeavyRow:
    """One heavy-plate freeze-bracket row (`crates/sweep --frozen-heavyplate`): the three EOS-only
    e_eff curves at one `(v, rho)` — equilibrium, sudden-freeze at the turnaround (pessimistic), and
    pure-H2O no-chemistry — plus the splice-consistency diagnostic."""

    v: float
    rho_impact: float
    e_eff_eq: float
    e_eff_frozen_rebound: float
    e_eff_frozen_all: float
    swap_energy_jump_frac: float


def read_frozen_heavyplate(path: Path = DEFAULT_FROZEN_SWEEP_PATH) -> list[FrozenHeavyRow]:
    """Parse the `--frozen-heavyplate` sweep JSONL (one JSON object per line; blanks tolerated)."""
    rows: list[FrozenHeavyRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(FrozenHeavyRow(**{f.name: float(d[f.name]) for f in fields(FrozenHeavyRow)}))
    return rows


@dataclass(frozen=True)
class FrozenHeavyPoint:
    """The freeze-timing bracket at one `(v, rho)`, translated onto the survivable f at that
    anchor's design collimation. The EOS-only frozen delta is applied to the coupled (radiation-on)
    headline e_eff, per ADR-0026 (the delta is composition kinetics, ~indep. of the coupling)."""

    v: float
    rho_impact: float
    e_eff_eq: float
    e_eff_frozen_rebound: float
    e_eff_frozen_all: float
    delta_frozen: float  # e_eff_eq - e_eff_frozen_rebound (chemistry-return content forfeited)
    e_eff_coupled: float  # radiation-on headline e_eff(v, rho) at the headline length, κ = 1
    eta_design: float  # the anchor's design-point collimation
    f_eq: float  # headline f at eta_design (no freeze penalty)
    f_frozen_rebound: float  # f under the sudden-freeze pessimistic bound
    f_frozen_all: float  # f under the pure-H2O no-chemistry bound


def frozen_heavyplate_bracket(
    frozen_rows: list[FrozenHeavyRow],
    sweep_rows: list[HeavyPlateRow],
    eta_design: dict[float, float],
) -> list[FrozenHeavyPoint]:
    """Per `(v, rho)`, translate the EOS-only freeze bracket onto f at that anchor's collimation.
    The frozen delta `e_eff_eq - e_eff_frozen_*` (EOS-only) is subtracted from the coupled headline
    `e_coupled(v, rho)` before the `f = eta·(1 + e_eff)/2` reconciliation."""
    points: list[FrozenHeavyPoint] = []
    interps = {v: e_eff_interpolator_at_v(sweep_rows, v) for v in eta_design}
    for r in sorted(frozen_rows, key=lambda x: (x.v, x.rho_impact)):
        if r.v not in eta_design:
            continue
        eta = eta_design[r.v]
        ec = interps[r.v](r.rho_impact)
        delta_frozen = r.e_eff_eq - r.e_eff_frozen_rebound
        delta_all = r.e_eff_eq - r.e_eff_frozen_all
        points.append(
            FrozenHeavyPoint(
                v=r.v,
                rho_impact=r.rho_impact,
                e_eff_eq=r.e_eff_eq,
                e_eff_frozen_rebound=r.e_eff_frozen_rebound,
                e_eff_frozen_all=r.e_eff_frozen_all,
                delta_frozen=delta_frozen,
                e_eff_coupled=ec,
                eta_design=eta,
                f_eq=reconcile_f(eta, ec),
                f_frozen_rebound=reconcile_f(eta, ec - delta_frozen),
                f_frozen_all=reconcile_f(eta, ec - delta_all),
            )
        )
    return points


def write_frozen_summary(
    points: list[FrozenHeavyPoint], path: Path = DEFAULT_FROZEN_SUMMARY_PATH
) -> None:
    """Write the freeze-timing bracket to a CSV (header = `FrozenHeavyPoint` field names)."""
    header = [f.name for f in fields(FrozenHeavyPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def _load_geometry(geo_path: Path) -> list[GeoRow]:
    """Load the geometry sweep (M = 40, falling back to the M = 20 production anchor), keeping only
    the physical `eta_capture` rows (dropping strong-shock solver blow-ups)."""
    if not geo_path.exists():
        geo_path = DEFAULT_GEOMETRY_M20_PATH
    geo_rows = [r for r in read_geometry(geo_path) if r.mach in (20.0, 40.0)]
    mach = max(r.mach for r in geo_rows)
    return [r for r in geo_rows if r.mach == mach and 0.0 < r.eta_capture <= ETA_PHYSICAL_MAX]


def _design_collimation(points: list[HeavyPlatePoint]) -> dict[float, float]:
    """The design-point `eta_capture` at each anchor velocity: the best-surviving concave case (the
    shape the frontier selects), falling back to the best-surviving flat case."""
    out: dict[float, float] = {}
    for v in V_ANCHORS:
        design = best_at_v(points, v, concave=True) or best_at_v(points, v)
        if design is not None:
            out[v] = design.eta_capture
    return out


def _run_frozen(frozen_path: Path, sweep_path: Path, geo_path: Path, summary_path: Path) -> None:
    """The `--frozen` path: resolve the headline design collimation per anchor, translate the
    EOS-only freeze bracket onto f there, write the CSV, and print the per-anchor bracket."""
    sweep_rows = read_heavyplate_sweep(sweep_path)
    geo_rows = _load_geometry(geo_path)
    points = heavyplate_frontier(geo_rows, sweep_rows, sweep_velocities(sweep_rows))
    eta_design = _design_collimation(points)

    frozen_rows = read_frozen_heavyplate(frozen_path)
    bracket = frozen_heavyplate_bracket(frozen_rows, sweep_rows, eta_design)
    write_frozen_summary(bracket, summary_path)

    jump_max = max((abs(r.swap_energy_jump_frac) for r in frozen_rows), default=0.0)
    print("python: heavy-plate 16-28 km/s freeze-timing bracket on f (per anchor design point):")
    for v in V_ANCHORS:
        anchor_pts = [p for p in bracket if abs(p.v - v) <= FLOAT_TOL]
        if not anchor_pts:
            continue
        # Report at the densest bracket rho (the design point trends dense/survivable).
        p = max(anchor_pts, key=lambda x: x.rho_impact)
        print(
            f"    v={v / 1000:4.0f} km/s (eta_design={p.eta_design:.3f}): "
            f"f = [{p.f_frozen_rebound:.3f}, {p.f_eq:.3f}] "
            f"(sudden-freeze .. equilibrium; delta_e={p.delta_frozen:+.3f})"
        )
    print(f"python: splice energy-jump max {jump_max:.2e} of incident KE; wrote {summary_path}")


def main() -> None:
    """Run the heavy-plate `f(v)` + facesheet-survivability frontier; print the design answers."""
    parser = argparse.ArgumentParser(
        description="Heavy-plate 16-28 km/s f(v) + facesheet survivability at the pinned plate."
    )
    parser.add_argument("--sweep", type=Path, default=DEFAULT_HEAVYPLATE_SWEEP_PATH)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY_M40_PATH)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument(
        "--frozen",
        action="store_true",
        help="translate the 16-28 km/s EOS-only freeze-timing bracket onto the survivable f "
        "(ADR-0026) -> frontier_frozen_heavyplate.csv, instead of the headline f(v) run",
    )
    parser.add_argument(
        "--frozen-sweep",
        type=Path,
        default=DEFAULT_FROZEN_SWEEP_PATH,
        help="the --frozen-heavyplate sweep JSONL (with --frozen)",
    )
    args = parser.parse_args()

    if args.frozen:
        summary = (
            args.summary if args.summary != DEFAULT_SUMMARY_PATH else DEFAULT_FROZEN_SUMMARY_PATH
        )
        _run_frozen(args.frozen_sweep, args.sweep, args.geometry, summary)
        return

    sweep_rows = read_heavyplate_sweep(args.sweep)
    geo_rows = _load_geometry(args.geometry)
    mach = max(r.mach for r in geo_rows)
    velocities = sweep_velocities(sweep_rows)

    points = heavyplate_frontier(geo_rows, sweep_rows, velocities)
    write_summary(points, args.summary)
    figs = plot_f_of_v(points, args.plot_dir)

    # Design-§12.1 diagnostics: e_eff should be flat in L (τ ≫ 1) and flat in opacity at the top.
    rho_probe = 0.08  # a representative mid-grid density present in every slice
    l_spreads = length_sensitivity(sweep_rows, rho_probe)
    tau_spread = opacity_sensitivity(sweep_rows, rho_probe)

    print(
        f"geometry anchor: M = {mach:.0f} ({len(geo_rows)} cases); plate R = {PLATE_RADIUS_M} m, "
        f"m = {PULSE_MASS_KG} kg, ≤ {PLATE_MASS_CEILING_KG / 1000:.0f} t"
    )
    print(
        f"L-sensitivity at rho={rho_probe} (τ ≫ 1 ⇒ ~0): "
        + ", ".join(f"{v / 1000:.0f} km/s Δe_eff={d:+.4f}" for v, d in sorted(l_spreads.items()))
    )
    print(
        f"opacity τ-check at 28 km/s, rho={rho_probe}: Δe_eff={tau_spread:+.4f} over 0.1x-10x κ "
        f"({'τ ≫ 1 confirmed' if abs(tau_spread) < 0.02 else 'opacity-sensitive — τ ~ 1'})"
    )

    # The 16 km/s overlap with the core envelope study (consistency check).
    overlap = best_at_v(points, 16_000.0, concave=True) or best_at_v(points, 16_000.0)
    if overlap is not None:
        print(
            f"16 km/s overlap check vs core study (f ~ 0.8): best survivable f = {overlap.f:.3f} "
            f"(e_eff={overlap.e_eff:.3f}, eta={overlap.eta_capture:.3f})"
        )

    print(f"\n{'v [km/s]':>9} | {'flat best f (shape)':<34} | concave best f (shape)")
    for v in velocities:
        flat = best_at_v(points, v, concave=False)
        conc = best_at_v(points, v, concave=True)

        def fmt(p: HeavyPlatePoint | None) -> str:
            if p is None:
                return "none survives".ljust(34)
            return (
                f"{p.f:.3f} (L/D={p.l_over_d:.1f}, rf/R={p.r_foot_over_r:.1f}, "
                f"d/D={p.d_over_d:.2f}, {p.plate_mass_t:.0f} t)"
            ).ljust(34)

        print(f"{v / 1000:9.1f} | {fmt(flat)} | {fmt(conc)}")

    print(f"\npython: wrote {args.summary} and {len(figs)} figure(s): " + ", ".join(map(str, figs)))


if __name__ == "__main__":
    main()
