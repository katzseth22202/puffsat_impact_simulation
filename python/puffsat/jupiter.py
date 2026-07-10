"""Jupiter-retrograde 69 km/s special-scenario analysis: plate sizing + survivable `f`.

The scenario ("Sorry No ISRU" launch-capability study): a **100 kg** PuffSat pulse arriving at
**69 km/s**, on a pusher plate whose mass budget is **up to 100 t** — the plate radius `R` is a
free design variable (the baseline study pinned `R = 5 m`), so the question is how *wide* the
plate must be to survive, whether the budget could be spent lighter, and whether the shallow
concave shape (Rung D-cc) is still worth its focusing penalty.

Physics of the sizing (all inherited from Rung S, ADR-0010):

- Peak facesheet pressure is **intensive**: `peak = c_stag·rho·v² · focusing`, blind to plate
  thickness/width as acreage. At 69 km/s, `v²` is 18.6x the 16 km/s value, so the survivable
  impact density collapses to `rho <= P_limit/(c_stag v²) ~ 0.07 kg/m³` (400 MPa baseline).
- The Σ contract (ADR-0003) `rho = m/(2π·(L/D)·(r_foot/R)³·R³)` is the only way to *reach* that
  dilution: `rho ∝ m/R³`, so **width is the whole game** — a wider plate dilutes the cloud
  cubically, while extra thickness only buys the relaxed material allowable (400 → 900 MPa,
  a factor 2.25 in `rho`, equivalent to just 1.31x in `R`).
- `e_eff(rho)` comes from the 69 km/s coupled sweep (`--jupiter`), interpolated per shape at its
  Σ-implied density, and bracketed over the table's opacity scale (the stagnated slab sits near
  `tau ~ 1` here, so unlike the 16 km/s anchor the radiative loss IS opacity-sensitive).
- `eta_capture` comes from the geometry sweep at the M = 40 strong-shock spot check (falling
  back to the production M = 20 anchor), scale-invariant in `R`.

Plate mass model: the baseline stack (SiC facesheet + solid Ti impedance layer + truss) is
~3-4 t at R = 5 m (design §2), i.e. an areal density of ~38-51 kg/m²; we carry 45 kg/m² central
with that band, and a concave dish pays its small extra-area factor `1 + (2·d/D)²`.
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

DEFAULT_JUPITER_SWEEP_PATH = Path("data/results/sweep_jupiter.jsonl")
DEFAULT_GEOMETRY_M40_PATH = Path("data/results/sweep_geometry_m40.jsonl")
DEFAULT_GEOMETRY_M20_PATH = Path("data/results/sweep_geometry.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/frontier_jupiter.csv")

V_JUPITER = 69_000.0  # m/s — Jupiter-retrograde encounter
PULSE_MASS_JUPITER_KG = 100.0
PLATE_BUDGET_KG = 100_000.0

# Realistic slug column length [m] for the e_eff anchor (the survivable cloud is ~10 m long);
# the 1 m production-convention rows are kept as a sensitivity check.
LENGTH_ANCHOR = 12.0

# Plate areal density [kg/m²]: baseline stack 3-4 t at R = 5 m (design §2) => 38-51; 45 central.
AREAL_DENSITY = 45.0
AREAL_DENSITY_BAND = (38.0, 51.0)

# Physical ceiling for eta_capture: shallow-concave over-collimation reaches ~1.01 (Rung D-cc),
# so anything past this is a solver blow-up (one M=40 case returned 7.6), not physics.
ETA_PHYSICAL_MAX = 1.2


@dataclass(frozen=True)
class JupiterRow:
    """One 69 km/s sweep row (the fields this analysis needs from the JSONL schema)."""

    rho_impact: float
    length: float
    opacity_scale: float
    e_eff: float
    peak_wall_pressure: float
    loss_radiative_wall: float
    loss_escape_space: float


def read_jupiter_sweep(path: Path = DEFAULT_JUPITER_SWEEP_PATH) -> list[JupiterRow]:
    """Parse the `--jupiter` sweep JSONL (one JSON object per line; blank lines tolerated)."""
    rows: list[JupiterRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            JupiterRow(
                rho_impact=float(d["rho_impact"]),
                length=float(d["length"]),
                opacity_scale=float(d["opacity_scale"]),
                e_eff=float(d["e_eff"]),
                peak_wall_pressure=float(d["peak_wall_pressure"]),
                loss_radiative_wall=float(d["loss_radiative_wall"]),
                loss_escape_space=float(d["loss_escape_space"]),
            )
        )
    return rows


def e_eff_interpolator(rows: list[JupiterRow], length: float, opacity_scale: float) -> _LogInterp:
    """`e_eff(rho)` along one `(length, opacity_scale)` slice, log-rho linear interpolation."""
    slice_rows = sorted(
        (r for r in rows if r.length == length and r.opacity_scale == opacity_scale),
        key=lambda r: r.rho_impact,
    )
    if not slice_rows:
        raise ValueError(f"no sweep rows at length={length}, opacity_scale={opacity_scale}")
    return _LogInterp(
        [r.rho_impact for r in slice_rows],
        [r.e_eff for r in slice_rows],
    )


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


def stagnation_coefficient_jupiter(rows: list[JupiterRow]) -> float:
    """`c_stag = peak_wall_pressure/(rho v²)` averaged over the sweep (ADR-0010: physical EOS
    peak, AV excluded). Opacity/length move it only weakly; the mean is the sizing anchor."""
    coeffs = [r.peak_wall_pressure / (r.rho_impact * V_JUPITER**2) for r in rows]
    c = sum(coeffs) / len(coeffs)
    if c <= 0.0:
        raise ValueError("non-positive c_stag — stale or empty sweep JSONL?")
    return c


def plate_mass(radius: float, d_over_d: float, areal_density: float = AREAL_DENSITY) -> float:
    """Plate mass [kg]: areal density x disk area, with the shallow dish's extra-area factor
    `1 + (2 d/D)²` (spherical-cap area `π(a² + d²)`, `d = (d/D)·2R`)."""
    return areal_density * math.pi * radius * radius * (1.0 + (2.0 * d_over_d) ** 2)


@dataclass(frozen=True)
class JupiterPoint:
    """One (plate radius x cloud shape) case resolved to survivability and `f` at 69 km/s."""

    plate_radius: float
    plate_mass_t: float  # tonnes, at AREAL_DENSITY (band scales linearly)
    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    rho_impact: float
    e_eff: float  # at the kappa-scale anchor (1.0); bracket carried separately
    e_eff_lo: float  # opacity-scale bracket, low end
    e_eff_hi: float  # opacity-scale bracket, high end
    eta_capture: float
    focusing_factor: float
    peak_compressive: float
    f: float
    f_lo: float
    f_hi: float
    survives_baseline: bool
    survives_relaxed: bool


def jupiter_frontier(
    geo_rows: list[GeoRow],
    sweep_rows: list[JupiterRow],
    radii: list[float],
    *,
    mass: float = PULSE_MASS_JUPITER_KG,
) -> list[JupiterPoint]:
    """Resolve every (plate radius x geometry case) to survivability and `f` at 69 km/s.

    Per case: the Σ contract gives `rho(R, shape)`; the 69 km/s sweep gives `e_eff(rho)` (at the
    realistic slug length, bracketed over opacity scale); the geometry row gives `eta_capture`
    and the concave focusing factor; the stagnation law + focusing gives the peak, classified
    against the 400 MPa baseline and the relaxed 900 MPa limit."""
    c_stag = stagnation_coefficient_jupiter(sweep_rows)
    e_mid = e_eff_interpolator(sweep_rows, LENGTH_ANCHOR, 1.0)
    e_lo = e_eff_interpolator(sweep_rows, LENGTH_ANCHOR, 0.1)
    e_hi = e_eff_interpolator(sweep_rows, LENGTH_ANCHOR, 10.0)
    relaxed_limit = max(P_LIMIT_HIGHV)

    flat_local = {
        (r.l_over_d, r.r_foot_over_r): r.peak_local_pressure
        for r in geo_rows
        if r.d_over_d == 0.0 and r.peak_local_pressure > 0.0
    }
    points: list[JupiterPoint] = []
    for radius in radii:
        for r in sorted(geo_rows, key=lambda x: (x.l_over_d, x.r_foot_over_r, x.d_over_d)):
            rho = impact_density(r.l_over_d, r.r_foot_over_r, mass, radius)
            ref = flat_local.get((r.l_over_d, r.r_foot_over_r))
            focusing = r.peak_local_pressure / ref if ref else 1.0
            peak = peak_facesheet_pressure(rho, V_JUPITER, c_stag) * focusing
            base = classify_survivability(peak, P_LIMIT_BASELINE, SIC_SPALL_LO)
            relaxed = classify_survivability(peak, relaxed_limit, SIC_SPALL_LO)
            # The opacity brackets are physical bounds, so order them.
            brackets = sorted((e_lo(rho), e_mid(rho), e_hi(rho)))
            points.append(
                JupiterPoint(
                    plate_radius=radius,
                    plate_mass_t=plate_mass(radius, r.d_over_d) / 1000.0,
                    d_over_d=r.d_over_d,
                    l_over_d=r.l_over_d,
                    r_foot_over_r=r.r_foot_over_r,
                    rho_impact=rho,
                    e_eff=brackets[1],
                    e_eff_lo=brackets[0],
                    e_eff_hi=brackets[2],
                    eta_capture=r.eta_capture,
                    focusing_factor=focusing,
                    peak_compressive=peak,
                    f=reconcile_f(r.eta_capture, brackets[1]),
                    f_lo=reconcile_f(r.eta_capture, brackets[0]),
                    f_hi=reconcile_f(r.eta_capture, brackets[2]),
                    survives_baseline=base.survives_compressive and base.survives_spall,
                    survives_relaxed=relaxed.survives_compressive and relaxed.survives_spall,
                )
            )
    return points


def best_at_radius(
    points: list[JupiterPoint], radius: float, *, concave: bool | None = None
) -> JupiterPoint | None:
    """The highest-`f` baseline-surviving case at one plate radius; `concave` filters d/D > 0
    (True) / d/D = 0 (False) / any (None)."""
    candidates = [
        p
        for p in points
        if p.plate_radius == radius
        and p.survives_baseline
        and (concave is None or (p.d_over_d > 0.0) == concave)
    ]
    return max(candidates, key=lambda p: p.f, default=None)


def write_summary(points: list[JupiterPoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write the Jupiter frontier to a CSV (header = `JupiterPoint` field names)."""
    header = [f.name for f in fields(JupiterPoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def main() -> None:
    """Run the plate-sizing frontier and print the design answers."""
    parser = argparse.ArgumentParser(description="Jupiter 69 km/s plate sizing + survivable f.")
    parser.add_argument("--sweep", type=Path, default=DEFAULT_JUPITER_SWEEP_PATH)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY_M40_PATH)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_PATH)
    args = parser.parse_args()

    sweep_rows = read_jupiter_sweep(args.sweep)
    geo_path: Path = args.geometry
    if not geo_path.exists():
        geo_path = DEFAULT_GEOMETRY_M20_PATH
    geo_rows = [r for r in read_geometry(geo_path) if r.mach in (20.0, 40.0)]
    mach = max(r.mach for r in geo_rows)
    geo_rows = [r for r in geo_rows if r.mach == mach]
    blown_up = [r for r in geo_rows if not 0.0 < r.eta_capture <= ETA_PHYSICAL_MAX]
    for r in blown_up:
        print(
            f"WARNING: dropping unphysical geometry row eta={r.eta_capture:.3f} "
            f"(d/D={r.d_over_d:.2f}, L/D={r.l_over_d:.1f}, r_foot/R={r.r_foot_over_r:.1f}, "
            f"M={r.mach:.0f}) — strong-shock solver blow-up, not physics"
        )
    geo_rows = [r for r in geo_rows if 0.0 < r.eta_capture <= ETA_PHYSICAL_MAX]

    radii = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0]
    points = jupiter_frontier(geo_rows, sweep_rows, radii)
    write_summary(points, args.summary)

    c_stag = stagnation_coefficient_jupiter(sweep_rows)
    rho_ceiling = P_LIMIT_BASELINE / (c_stag * V_JUPITER**2)
    print(f"c_stag = {c_stag:.3f}; rho ceiling (400 MPa) = {rho_ceiling:.4f} kg/m^3")
    print(f"geometry anchor: M = {mach:.0f} ({len(geo_rows)} cases)")
    print(f"{'R [m]':>6} {'mass [t]':>9} | flat best f (shape)      | concave best f (shape)")
    for radius in radii:
        flat = best_at_radius(points, radius, concave=False)
        conc = best_at_radius(points, radius, concave=True)

        def fmt(p: JupiterPoint | None) -> str:
            if p is None:
                return "  none survives      "
            return (
                f"{p.f:.3f} [{p.f_lo:.3f},{p.f_hi:.3f}] "
                f"(L/D={p.l_over_d:.1f}, rf/R={p.r_foot_over_r:.1f}, d/D={p.d_over_d:.2f})"
            )

        mass_t = plate_mass(radius, 0.0) / 1000.0
        print(f"{radius:6.1f} {mass_t:9.1f} | {fmt(flat)} | {fmt(conc)}")


if __name__ == "__main__":
    main()
