"""The thick Sun disk: a 1 cm, 22 kg/m^2 graphite plate in place of the 3.6 mm, 8 kg/m^2 disk.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.thick_disk

The request is `docs/thick_disk_impact_request.md` (2026-10-27). The Sun boost flies graphite disks
head-on at ~392 km/s through the bore of Result 22's δ/D 0.5 paraboloid into a polyethylene charge
`k = 20` times their mass at 1.5 F. The pieces:

1. **Coupling** (`coupling_rows`, items 1-2). The 2-D runs (`sweep --thick-disk [coarse|fine]`).
   The gamma gas has no material, so the disk's areal density `sigma` enters only as its charge's
   column, `L = k sigma / rho`, the same for every disk of one thickness. A heavier 1 cm disk is
   only wider: its `L / R` falls. The shocked gas leaks sideways while the disk sweeps the column,
   and that leak, not tunnelling, is what a thick disk costs.
2. **The bowl** (`bowl_rows`, item 2). Each pulse's impulse on the paper's 5 t bowl: its free kick,
   the lossless spring's stroke at 5 t (`m s = J tan(theta/4) / (2 theta f)`, Result 17) and the
   bowl that an 8 m stroke needs.
3. **Off centre** (`offset_rows`, item 3). The 2-D code is axisymmetric. A disk offset `d` from a
   charge of radius `R_c` strikes it over the two circles' lens (`shaped_plug.passed_share`); the
   lune that misses flies on, away from the plate, and is lost. The struck share's impulse is
   bracketed: the lone lens (`J ~ phi`), or the lens with the unstruck charge joining it at no
   cost in energy (`J ~ sqrt(E M)`). The fireball's centroid off the axis turns the impulse; the
   ballistic tracer (`iapetus_bag.trace`, drifting off the plate) gives that side impulse.
   Wide charges (`R_c > R`) catch the whole disk out to `R_c - R` at the 2-D runs' price.
4. **A split disk** (`split_rows`, item 4). Two halves pushed apart across the crack by `g`; each
   half's rim overhangs the charge. Along the axis a gap is a delay of `g / v` against the sweep.
5. **An alumina coat** (`alumina_rows`, item 5). Its share of the pulse, and two bounds: its
   atomization and first ionization against the fireball's heat, and its radiation scaled from
   argon's 15% at the Sun bowl (Result 22) by its share of the atoms.

The parent's S62 (`docs/thin_disk_wide_charge_asks_for_impact_sim.md`) asks the same of the thin
1 kg disk:

6. **The thin disk off centre** (`offset_rows(..., 1.0, THIN)`, `side_rows`). Its bore clearance,
   and the side impulse traced at each offset up to it, with the slope through the origin and each
   point's departure from it (linearity).
7. **4 R_sun** (`speed_rows`). The 2-D gas has no speed scale, so every ratio above carries over;
   the pulse, side impulse and bowl scale with the closing speed (`mission.periapsis`), the peak
   with its square. Polyethylene's frozen ion and carbon toll is the one speed-dependent loss
   (`freeze.freeze`); it is the same on every charge width, so it moves no ratio.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

from puffsat.xray_collision import charge, freeze, gas_spring, mission, shaped_plug, spring_trim
from puffsat.xray_collision import iapetus_bag as ib

OUT_DIR = charge.OUT_DIR
JSONL = OUT_DIR / "thick_disk_2d.jsonl"
JSONL_COARSE = OUT_DIR / "thick_disk_2d_coarse.jsonl"
JSONL_FINE = OUT_DIR / "thick_disk_2d_fine.jsonl"
#: The thin disk on the wide charges (`sweep --thick-disk thin-wide [coarse]`), merged by grid.
JSONL_THIN_WIDE = OUT_DIR / "thick_disk_2d_thinwide.jsonl"
JSONL_THIN_WIDE_COARSE = OUT_DIR / "thick_disk_2d_thinwide_coarse.jsonl"
JSONL_THIN_WIDE_FINE = OUT_DIR / "thick_disk_2d_thinwide_fine.jsonl"
GRIDS = (
    ("x1", (JSONL, JSONL_THIN_WIDE)),
    ("x0.5", (JSONL_COARSE, JSONL_THIN_WIDE_COARSE)),
    ("x2", (JSONL_FINE, JSONL_THIN_WIDE_FINE)),
)

# --- The case (the request, and the child paper's Sun boost) ---------------------------------

V = 392.0e3  # [m/s], closing head-on at perihelion
K = 20.0  # the Sun bowl's charge (Result 22)
GRAPHITE_RHO = 2200.0  # [kg/m^3], the companion's molded graphite
THIN, THICK = 8.0, 22.0  # [kg/m^2]: the 3.6 mm disk and the 1 cm plate
BORE_R = 0.25  # [m], the Sun bowl's on-axis bore in the runs (`CHARGE_HOLE_R`)
SUN_PLATE_R = 10.0  # [m], the 20 m bowl (`CHARGE_PLATE_R`)
#: The paper's Sun boost: 7.2 t of disks give 9.8 km/s. Its per-kg kick is the baseline's.
DISKS_TOTAL, DV_TOTAL = 7200.0, 9.8e3
DV_PER_KG = DV_TOTAL / DISKS_TOTAL  # [m/s per kg of disk], at the baseline's J/mv
BOWL = 5000.0  # [kg], the paper's Sun bowl
STROKE = gas_spring.STROKE  # [m], Seth's
#: `(label, disk kg, cadence [Hz], charge radius over the disk's)`: the paper's pulse, the
#: proposal's, and the proposal on the x1.3 charge.
PULSES = (
    ("1 kg at 2 Hz (paper)", 1.0, 2.0, 1.0),
    ("2.5 kg at 1 Hz (proposed)", 2.5, 1.0, 1.0),
    ("2.5 kg at 1 Hz, charge x1.3", 2.5, 1.0, 1.3),
)
#: The parent's 618 km/s, 1 kg X-ray pusher disk: the per-pulse energy the 2.5 kg disk matches.
PARENT_V, PARENT_DISK = 618.0e3, 1.0


def disk_radius(mass: float, areal: float) -> float:
    """[m]"""
    return math.sqrt(mass / (math.pi * areal))


def on_grid(radius: float, dr: float) -> float:
    """[m] The radius a run actually resolves. `DiskCharge` fills (and the bore's bookkeeping
    takes) every cell whose centre lies inside `radius`, rescaling density to keep the mass:
    `n = ceil(radius / dr - 1/2)` cells, so `radius` rounds to the nearest cell edge."""
    return math.ceil(radius / dr - 0.5) * dr


def column_length(areal: float, k: float = K, rho: float = 950.0) -> float:
    """[m] A charge as wide as its disk: `k sigma / rho`, whatever the disk's mass."""
    return k * areal / rho


# --- 1. Coupling -----------------------------------------------------------------------------


def _f(r: dict[str, object], key: str) -> float:
    return float(str(r[key]))


@cache
def _runs(path: Path) -> tuple[dict[str, object], ...]:
    return ib.runs(path)


def _find(
    rows: tuple[dict[str, object], ...], kind: str, mass: float, areal: float, wide: float = 1.0
) -> dict[str, object] | None:
    for r in rows:
        if str(r["kind"]) != kind:
            continue
        if kind == "point":
            return r
        same_disk = abs(_f(r, "disk_kg") - mass) < 1e-9 and abs(_f(r, "areal") - areal) < 1e-6
        ratio = _f(r, "charge_radius") / _f(r, "disk_radius")
        if same_disk and abs(ratio - wide) < 1e-6:
            return r
    return None


@dataclass(frozen=True)
class CouplingRow:
    grid: str
    case: str
    kind: str
    disk_kg: float
    areal: float
    disk_radius: float
    charge_radius: float
    bore: float  # [m]
    length: float  # [m]
    length_over_radius: float  # the column's L / R_c
    cells_across: float  # disk radius over the cell
    dr: float  # [m] the run's cell
    disk_radius_on_grid: float  # [m] what the run resolves (`on_grid`)
    charge_radius_on_grid: float
    wide_on_grid: float  # the charge's resolved radius over the disk's: the width actually run
    impulse_per_mv: float  # J / (m v), the bore's part removed
    over_point: float  # over the point fireball's on the same grid
    over_baseline: float  # over the 1 kg, 8 kg/m^2 disk's on the same grid
    over_own_matched: float  # over the same disk on a charge as wide as itself (the wide gain)
    pulse_ns: float  # [N s]
    dv_per_pulse: float  # [m/s] the paper's per-kg kick, scaled by over_baseline
    peak_gpa: float  # peak local pressure outside the bore (not converged; Result 22)
    peak_radius: float
    peak_over_baseline: float  # over the thin disk's peak on the same grid (grid-to-grid ratio)
    energy_gj: float  # the disk's kinetic energy in the bowl's frame
    heat_gj: float  # its inelastic share, the fireball's heat: k / (1 + k)


def coupling_rows() -> list[CouplingRow]:
    out = []
    for grid, paths in GRIDS:
        rows = tuple(r for path in paths for r in _runs(path))
        if not rows:
            continue
        point = _find(rows, "point", 0.0, 0.0)
        base = _find(rows, "matched", 1.0, THIN)
        for r in rows:
            kind = str(r["kind"])
            m = 1.0 if kind == "point" else _f(r, "disk_kg")
            j = _f(r, "impulse_per_mv")
            ratio = j / _f(base, "impulse_per_mv") if base else math.nan
            peak = _f(r, "peak_per_mv2") * m * V * V
            base_peak = _f(base, "peak_per_mv2") * V * V if base else math.nan
            rc = _f(r, "charge_radius")
            own = (
                _find(rows, "matched", m, _f(r, "areal"))
                if kind == "wide"
                else (r if kind == "matched" else None)
            )
            out.append(
                CouplingRow(
                    grid=grid,
                    case=str(r["case"]),
                    kind=kind,
                    disk_kg=m,
                    areal=_f(r, "areal"),
                    disk_radius=_f(r, "disk_radius"),
                    charge_radius=rc,
                    bore=_f(r, "bore") if "bore" in r else BORE_R,
                    length=_f(r, "length"),
                    length_over_radius=_f(r, "length") / rc,
                    cells_across=_f(r, "disk_radius") / _f(r, "dr"),
                    dr=_f(r, "dr"),
                    disk_radius_on_grid=on_grid(_f(r, "disk_radius"), _f(r, "dr")),
                    charge_radius_on_grid=on_grid(rc, _f(r, "dr")),
                    wide_on_grid=on_grid(rc, _f(r, "dr"))
                    / on_grid(_f(r, "disk_radius"), _f(r, "dr")),
                    impulse_per_mv=j,
                    over_point=j / _f(point, "impulse_per_mv") if point else math.nan,
                    over_baseline=ratio,
                    over_own_matched=j / _f(own, "impulse_per_mv") if own else math.nan,
                    pulse_ns=j * m * V,
                    dv_per_pulse=DV_PER_KG * m * ratio,
                    peak_gpa=peak / 1e9,
                    peak_radius=_f(r, "peak_radius"),
                    peak_over_baseline=peak / base_peak,
                    energy_gj=0.5 * m * V * V / 1e9,
                    heat_gj=0.5 * m * V * V * K / (1.0 + K) / 1e9,
                )
            )
    return out


@dataclass(frozen=True)
class ConvergenceRow:
    case: str
    disk_kg: float
    areal: float
    wide: float
    coarse: float  # over the 1 kg, 8 kg/m^2 disk on each grid (0.05 / 0.025 / 0.0125 m)
    production: float
    fine: float
    change_last: float  # the last refinement's change, relative


def convergence_rows(coupling: list[CouplingRow]) -> list[ConvergenceRow]:
    """Each matched or wide case's ratio to the thin baseline on every grid that ran it."""
    out = []
    keys = sorted(
        {
            (c.case, c.disk_kg, c.areal, round(c.charge_radius / c.disk_radius, 6))
            for c in coupling
            if c.kind != "point"
        },
        key=lambda t: (t[3], t[1], t[2]),
    )
    for case, m, areal, wide in keys:
        by = {
            c.grid: c.over_baseline
            for c in coupling
            if c.case == case and c.disk_kg == m and c.areal == areal
        }
        vals = [by.get(g, math.nan) for g in ("x0.5", "x1", "x2")]
        done = [v for v in vals if not math.isnan(v)]
        change = done[-1] / done[-2] - 1.0 if len(done) >= 2 else math.nan
        out.append(
            ConvergenceRow(
                case=case,
                disk_kg=m,
                areal=areal,
                wide=wide,
                coarse=vals[0],
                production=vals[1],
                fine=vals[2],
                change_last=change,
            )
        )
    return out


def production_grid(rows: list[CouplingRow]) -> str:
    """The finest grid that ran the case set's matched disks."""
    grids = {r.grid for r in rows}
    return next(g for g in ("x2", "x1", "x0.5") if g in grids)


def ratio_of(rows: list[CouplingRow], grid: str, mass: float, areal: float, wide: float) -> float:
    """J / (m v) of a case over the matched 1 kg, 8 kg/m^2 disk's, on one grid."""
    for r in rows:
        if r.grid == grid and r.kind != "point":
            same = abs(r.disk_kg - mass) < 1e-9 and abs(r.areal - areal) < 1e-6
            if same and abs(r.charge_radius / r.disk_radius - wide) < 1e-6:
                return r.over_baseline
    return math.nan


@dataclass(frozen=True)
class WidthRow:
    disk: str
    wide: float  # nominal charge radius over the disk's
    grid: str
    dr: float
    wide_on_grid: float  # the width the run resolved
    length_over_radius_on_grid: float  # the column's L over the resolved charge radius
    over_own_matched: float  # J over the same disk's matched charge on the same grid


def width_rows(coupling: list[CouplingRow]) -> list[WidthRow]:
    """Each disk's gain over its matched charge against the width each grid resolved (S62:
    nominal widths a cell apart run as one width on the coarser grids)."""
    out = []
    for c in coupling:
        if c.kind == "point":
            continue
        out.append(
            WidthRow(
                disk=f"{c.disk_kg:.1f} kg at {c.areal:g} kg/m^2",
                wide=round(c.charge_radius / c.disk_radius, 6),
                grid=c.grid,
                dr=c.dr,
                wide_on_grid=c.wide_on_grid,
                length_over_radius_on_grid=c.length / c.charge_radius_on_grid,
                over_own_matched=c.over_own_matched,
            )
        )
    return sorted(out, key=lambda w: (w.disk, w.dr, w.wide))


# --- 2. The bowl -----------------------------------------------------------------------------


def spring_stroke(pulse: float, cadence: float, mass: float, theta: float) -> float:
    """[m] The lossless spring's stroke (Result 17): `s = J tan(theta/4) / (2 theta f m)`."""
    return pulse * math.tan(theta / 4.0) / (2.0 * theta * cadence * mass)


def spring_mass(pulse: float, cadence: float, stroke: float, theta: float) -> float:
    """[kg] The bowl that holds `stroke`."""
    return pulse * math.tan(theta / 4.0) / (2.0 * theta * cadence * stroke)


@dataclass(frozen=True)
class BowlRow:
    pulse_case: str
    grid: str
    disk_kg: float
    cadence: float
    pulse_ns: float
    mean_force_mn: float
    free_kick: float  # [m/s] the 5 t bowl, unsprung, per pulse
    theta: str
    stroke_5t: float  # [m] at the paper's 5 t
    bowl_t_8m: float  # [t] at Seth's 8 m stroke
    dv_per_pulse: float  # [m/s] to the ship


def bowl_rows(coupling: list[CouplingRow]) -> list[BowlRow]:
    grid = production_grid(coupling)
    out = []
    for label, mass, cadence, wide in PULSES:
        areal = THIN if mass == 1.0 else THICK
        ratio = ratio_of(coupling, grid, mass, areal, wide)
        base = next(
            r.impulse_per_mv
            for r in coupling
            if r.grid == grid and r.kind == "matched" and r.disk_kg == 1.0 and r.areal == THIN
        )
        pulse = base * ratio * mass * V
        for name in ("2pi/3", "5pi/6"):
            theta = spring_trim.THETAS[name]
            out.append(
                BowlRow(
                    pulse_case=label,
                    grid=grid,
                    disk_kg=mass,
                    cadence=cadence,
                    pulse_ns=pulse,
                    mean_force_mn=pulse * cadence / 1e6,
                    free_kick=pulse / BOWL,
                    theta=name,
                    stroke_5t=spring_stroke(pulse, cadence, BOWL, theta),
                    bowl_t_8m=spring_mass(pulse, cadence, STROKE, theta) / 1e3,
                    dv_per_pulse=DV_PER_KG * mass * ratio,
                )
            )
    return out


# --- 3. Off centre ---------------------------------------------------------------------------


def lens_share(offset: float, charge_r: float, disk_r: float) -> float:
    """The share of the disk's face over the charge's when their centres are `offset` apart."""
    return 1.0 - shaped_plug.passed_share(offset, charge_r, disk_r)


def struck_bracket(phi: float, k: float = K) -> tuple[float, float]:
    """Impulse over the centred strike's when `phi` of the disk strikes `phi` of the charge.
    Low: the lens alone, `J ~ phi`. High: the unstruck charge joins it at no energy cost,
    `J ~ sqrt(E M)` with `E ~ phi` and `M ~ phi (1 + k) + (1 - phi) k`."""
    high = math.sqrt(phi * (phi * (1.0 + k) + (1.0 - phi) * k) / (1.0 + k))
    return phi, high


#: The ballistic bowl: Result 22's δ/D 0.5 paraboloid with the burst at 1.5 F.
_SHAPE = ib.Shape(0.5, 0.0)
_HEIGHT = 1.5 * _SHAPE.focus
#: Misses at which the tracer is run [m, on the 20 m bowl], interpolated between.
TRACE_MISSES = (0.0, 0.025, 0.05, 0.1, 0.2, 0.4)


@cache
def traced() -> tuple[np.ndarray, np.ndarray]:
    """`(axial over centred, lateral over axial)` at each of `TRACE_MISSES`. The tracer's plate is
    Result 20's 9.5 m floor; the bowl is the same shape at 10 m, so the miss is scaled."""
    scale = ib.PLATE_R / SUN_PLATE_R
    axial, lateral = [], []
    j0 = 0.0
    for miss in TRACE_MISSES:
        j, _, _ = ib.trace(_SHAPE, _HEIGHT, miss * scale, K, away=True)
        if miss == 0.0:
            j0 = float(-j[2])
        axial.append(float(-j[2]) / j0)
        lateral.append(float(j[0]) / float(-j[2]))
    return np.array(axial), np.array(lateral)


def off_focus(centroid: float) -> tuple[float, float]:
    """`(axial factor, lateral over axial)` for the fireball's centroid `centroid` [m] off the
    bowl's axis."""
    axial, lateral = traced()
    misses = np.array(TRACE_MISSES)
    c = min(abs(centroid), misses[-1])
    return float(np.interp(c, misses, axial)), float(np.interp(c, misses, lateral))


def coupling_at(offset: float, charge_r: float, disk_r: float) -> tuple[float, float, float]:
    """`(low, high, lateral over axial)`: impulse over the centred strike's, for a disk `offset`
    off the charge (and bowl) axis. The burst sits at the struck mass's centroid: `offset / 2` for
    the matched charge, whose struck lens is symmetric; for a wide charge, `offset` itself, exact
    while the disk is covered and an upper bound on the side impulse past that."""
    phi = lens_share(offset, charge_r, disk_r)
    low, high = struck_bracket(phi)
    centroid = 0.5 * offset if charge_r <= disk_r + 1e-12 else offset
    axial, lateral = off_focus(centroid)
    return low * axial, high * axial, lateral


@dataclass(frozen=True)
class OffsetRow:
    grid: str
    charge: str  # `x1.00` matched, or the wide charge's radius over the disk's
    disk_radius: float
    charge_radius: float
    centred_over_matched: float  # the 2-D runs: this charge's centred J over the matched one's
    margin: float  # [m] full overlap out to here, `R_c - R`
    bore: float  # [m] the run's bore: 0.25 m, or the wide charge's radius
    bore_clearance: float  # [m] the disk clips the bore past this offset; it caps the tolerance
    drop: float  # 0.05 or 0.20, below the matched charge's centred impulse
    offset_low: float  # [m] the offset at that drop, lone-lens bound
    offset_high: float  # [m] with the unstruck charge joining
    offset_low_over_rc: float
    offset_high_over_rc: float
    #: [m] the usable tolerance: the lone-lens offset, capped where the disk clips the bore. NaN
    #: if the centred strike is already below the drop.
    tolerance: float
    side_impulse_low: float  # lateral over axial at the usable tolerance
    side_ns_low: float  # [N s] at the disk's pulse at 392 km/s, at the usable tolerance


def _solve(target: float, f: object, hi: float) -> float:
    """The smallest offset in [0, hi] where `f(offset) <= target` (bisection; `f` falls)."""
    assert callable(f)
    if f(hi) > target:
        return math.nan
    lo = 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


def offset_rows(
    coupling: list[CouplingRow], mass: float = 2.5, areal: float = THICK
) -> list[OffsetRow]:
    grid = production_grid(coupling)
    r = disk_radius(mass, areal)
    matched = ratio_of(coupling, grid, mass, areal, 1.0)
    pulse = next(
        c.pulse_ns
        for c in coupling
        if c.grid == grid and c.kind == "matched" and c.disk_kg == mass and c.areal == areal
    )
    out = []
    wides = sorted({c.charge_radius / c.disk_radius for c in coupling if c.grid == grid})
    for wide in [w for w in wides if w >= 1.0 - 1e-9]:
        centred = ratio_of(coupling, grid, mass, areal, round(wide, 6)) / matched
        if math.isnan(centred):
            continue
        rc = r * wide
        bore = next(
            c.bore
            for c in coupling
            if c.grid == grid
            and c.kind != "point"
            and c.disk_kg == mass
            and c.areal == areal
            and abs(c.charge_radius / c.disk_radius - wide) < 1e-6
        )
        for drop in (0.05, 0.20):
            target = 1.0 - drop

            def low(d: float, rc: float = rc, c: float = centred) -> float:
                return c * coupling_at(d, rc, r)[0]

            def high(d: float, rc: float = rc, c: float = centred) -> float:
                return c * coupling_at(d, rc, r)[1]

            if centred <= target:
                d_low = d_high = math.nan
            else:
                d_low = _solve(target, low, rc + r)
                d_high = _solve(target, high, rc + r)
            usable = min(d_low, bore - r) if not math.isnan(d_low) else math.nan
            lat = coupling_at(usable, rc, r)[2] if not math.isnan(usable) else math.nan
            out.append(
                OffsetRow(
                    grid=grid,
                    charge=f"x{wide:.2f}",
                    disk_radius=r,
                    charge_radius=rc,
                    centred_over_matched=centred,
                    margin=rc - r,
                    bore=bore,
                    bore_clearance=bore - r,
                    drop=drop,
                    offset_low=d_low,
                    offset_high=d_high,
                    offset_low_over_rc=d_low / rc,
                    offset_high_over_rc=d_high / rc,
                    tolerance=usable,
                    side_impulse_low=lat,
                    side_ns_low=lat * pulse * low(usable) if not math.isnan(usable) else math.nan,
                )
            )
    return out


# --- 6. Side impulse against offset (S62) ----------------------------------------------------

#: Offsets [m] at which the side impulse is traced directly (not interpolated), to the bore.
SIDE_OFFSETS = (0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.095)


@cache
def side_traced(centroid: float) -> tuple[float, float]:
    """`(axial over centred, lateral over axial)`, traced at the fireball's `centroid` [m]."""
    scale = ib.PLATE_R / SUN_PLATE_R
    j0, _, _ = ib.trace(_SHAPE, _HEIGHT, 0.0, K, away=True)
    j, _, _ = ib.trace(_SHAPE, _HEIGHT, centroid * scale, K, away=True)
    return float(j[2] / j0[2]), float(j[0]) / float(-j[2])


@dataclass(frozen=True)
class SideRow:
    grid: str
    disk: str
    charge: str
    bore_clearance: float  # [m]
    offset: float  # [m], at or inside the bore's clearance
    struck_share: float  # of the disk, over the charge
    low: float  # impulse over this charge's centred strike: lone lens / unstruck charge joining
    high: float
    over_matched_low: float  # the same, over the disk's matched charge centred
    lateral: float  # side impulse over axial, traced at the fireball's centroid
    slope: float  # [1/m] this charge's least-squares `lateral / offset` through the origin
    off_line: float  # this point's departure from the slope, relative
    side_ns: float  # [N s] at 392 km/s
    side_ns_4rsun: float  # [N s] at 4 R_sun's closing speed


def side_rows(coupling: list[CouplingRow], mass: float, areal: float) -> list[SideRow]:
    """Side impulse at each `SIDE_OFFSETS` inside the bore's clearance, on each charge the finest
    grid ran for this disk. The burst sits at the struck mass's centroid, as in `coupling_at`."""
    grid = production_grid(coupling)
    r = disk_radius(mass, areal)
    base_j = next(
        c.impulse_per_mv
        for c in coupling
        if c.grid == grid and c.kind == "matched" and c.disk_kg == mass and c.areal == areal
    )
    v4 = closing(4.0)
    out = []
    for o in offset_rows(coupling, mass, areal):
        if o.drop != 0.05:
            continue
        wide = o.charge_radius / r
        rows = []
        for d in [x for x in SIDE_OFFSETS if x <= o.bore_clearance + 1e-9]:
            phi = lens_share(d, o.charge_radius, r)
            lo, hi = struck_bracket(phi)
            centroid = 0.5 * d if wide <= 1.0 + 1e-9 else d
            axial, lateral = side_traced(centroid)
            j = base_j * o.centred_over_matched * lo * axial * mass
            rows.append((d, phi, lo * axial, hi * axial, lateral, j))
        ds = np.array([x[0] for x in rows])
        lats = np.array([x[4] for x in rows])
        slope = float(ds @ lats / (ds @ ds)) if len(rows) else math.nan
        for d, phi, lo, hi, lateral, j in rows:
            out.append(
                SideRow(
                    grid=grid,
                    disk=f"{mass:.1f} kg at {areal:.0f} kg/m^2",
                    charge=o.charge,
                    bore_clearance=o.bore_clearance,
                    offset=d,
                    struck_share=phi,
                    low=lo,
                    high=hi,
                    over_matched_low=lo * o.centred_over_matched,
                    lateral=lateral,
                    slope=slope,
                    off_line=lateral / (slope * d) - 1.0,
                    side_ns=lateral * j * V,
                    side_ns_4rsun=lateral * j * v4,
                )
            )
    return out


# --- 7. 4 R_sun (S62 item 4) -----------------------------------------------------------------


def closing(r_sun: float) -> float:
    """[m/s] Head-on closing speed at perihelion `r_sun` (`mission.periapsis`); 392 km/s at the
    headline's 9.86 is `V`."""
    return V if r_sun == 9.86 else mission.periapsis(r_sun).closing_speed


@dataclass(frozen=True)
class SpeedRow:
    depth: str
    closing: float  # [m/s]
    case: str
    grid: str
    over_matched: float  # the 2-D ratio, speed-free
    impulse_per_mv: float
    pulse_ns: float  # [N s]
    energy_gj: float
    peak_over_392: float  # the peak scales as v^2 (the ratio to 392 km/s; absolute not converged)
    clearance: float  # [m], geometric
    side_ns_at_clearance: float  # [N s], the lone-lens bound at the bore's edge
    pe_velocity_factor: float  # polyethylene's sqrt(1 - frozen ions - carbon toll) (`freeze`)
    bowl_t_8m_2hz: str  # [t] the lossless spring at 8 m and 2 Hz, theta 2pi/3-5pi/6


@cache
def pe_velocity_factor(w: float, mass: float = 1.0) -> float:
    """Polyethylene's frozen-ion and carbon-bond toll on the rebound speed (Result 22's pricing)."""
    f = freeze.freeze(K, charge.ch2_table(), 3.75, v_rel=w, projectile=mass)
    carbon = (freeze.CH2_ATOMIZATION - freeze.H2_BOND) / freeze.CH2_ATOMIZATION
    return math.sqrt(max(1.0 - f.ion_share_frozen - f.bond_share * carbon, 0.0))


def speed_rows(
    coupling: list[CouplingRow], mass: float = 1.0, areal: float = THIN
) -> list[SpeedRow]:
    """The thin disk's matched and wide charges at 9.86 and 4 R_sun."""
    grid = production_grid(coupling)
    offsets = [o for o in offset_rows(coupling, mass, areal) if o.drop == 0.05]
    out = []
    for depth in (9.86, 4.0):
        w = closing(depth)
        for o in offsets:
            wide = o.charge_radius / o.disk_radius
            ratio = ratio_of(coupling, grid, mass, areal, round(wide, 6))
            row = next(
                c
                for c in coupling
                if c.grid == grid
                and c.kind != "point"
                and c.disk_kg == mass
                and c.areal == areal
                and abs(c.charge_radius / c.disk_radius - wide) < 1e-6
            )
            pulse = row.impulse_per_mv * mass * w
            bowls = [
                spring_mass(pulse, 2.0, STROKE, spring_trim.THETAS[n]) / 1e3
                for n in ("2pi/3", "5pi/6")
            ]
            out.append(
                SpeedRow(
                    depth=f"{depth:g} R_sun",
                    closing=w,
                    case=row.case,
                    grid=grid,
                    over_matched=ratio / ratio_of(coupling, grid, mass, areal, 1.0),
                    impulse_per_mv=row.impulse_per_mv,
                    pulse_ns=pulse,
                    energy_gj=0.5 * mass * w * w / 1e9,
                    peak_over_392=(w / V) ** 2,
                    clearance=o.bore_clearance,
                    side_ns_at_clearance=o.side_ns_low * w / V,
                    pe_velocity_factor=pe_velocity_factor(w, mass),
                    bowl_t_8m_2hz=f"{min(bowls):.1f}-{max(bowls):.1f}",
                )
            )
    return out


# --- 4. A split disk -------------------------------------------------------------------------


def split_missed(gap: float, disk_r: float, charge_r: float, n: int = 801) -> float:
    """The share of a disk that misses the charge when it has cracked along a diameter and the
    halves have moved `gap` apart across the crack, the pair centred on the charge."""
    xs = np.linspace(-disk_r, disk_r, n)
    x, y = np.meshgrid(xs, xs)
    inside = x * x + y * y <= disk_r * disk_r
    shifted = x + np.sign(x) * 0.5 * gap
    missed = inside & (shifted * shifted + y * y > charge_r * charge_r)
    return float(missed.sum()) / float(inside.sum())


@dataclass(frozen=True)
class SplitRow:
    charge: str
    charge_radius: float
    gap_cm: float
    missed_share: float  # of the disk, flying on past the charge
    low: float  # impulse over the whole disk's, the two bounds of item 3
    high: float
    axial_delay_over_sweep: float  # a gap along the axis instead: g / v over the sweep time


def split_rows() -> list[SplitRow]:
    r = disk_radius(2.5, THICK)
    length = column_length(THICK)
    sweep = length * (1.0 + 0.5 * K) / V  # Result 5's `L (1 + k/2) / v`
    out = []
    for wide in (1.0, 1.15, 1.3):
        rc = r * wide
        for gap_cm in (1.0, 2.0, 3.0, 5.0, 10.0):
            g = gap_cm / 100.0
            miss = split_missed(g, r, rc)
            low, high = struck_bracket(1.0 - miss)
            out.append(
                SplitRow(
                    charge=f"x{wide:.2f}",
                    charge_radius=rc,
                    gap_cm=gap_cm,
                    missed_share=miss,
                    low=low,
                    high=high,
                    axial_delay_over_sweep=(g / V) / sweep,
                )
            )
    return out


# --- 5. An alumina coat ----------------------------------------------------------------------

ALUMINA_RHO = 3950.0  # [kg/m^3]
ALUMINA_M = 0.10196  # [kg/mol]
#: Atomization of Al2O3(s) to 2 Al + 3 O, from JANAF heats of formation at 298 K [J/mol]:
#: -1675.7 (Al2O3) against +330.0 (Al) and +249.2 (O) kJ/mol.
ALUMINA_ATOMIZE = (2.0 * 330.0 + 3.0 * 249.2 + 1675.7) * 1e3
EV_PER_MOL = 96485.33
#: First ionization [eV]: Al 5.986, O 13.618 (NIST ASD).
ALUMINA_ION = (2.0 * 5.986 + 3.0 * 13.618) * EV_PER_MOL
#: Argon's charge on the Sun bowl radiates its stagnation heat into the face: 0.85 of
#: polyethylene's impulse (Result 22, `snow_charge_material.csv`).
ARGON_SUN_LOSS = 0.15
#: Allowance for alumina's ions radiating more per ion than argon's at the same state.
RADIATION_ALLOWANCE = 3.0
CH2_M, C_M = 0.014027, 0.012011  # [kg/mol]


@dataclass(frozen=True)
class AluminaRow:
    coat_um: float
    coat_areal: float  # [kg/m^2]
    share_of_disk: float
    share_of_pulse_mass: float
    atom_fraction: float  # of the pulse's atoms
    toll_bound: float  # atomization + first ionization over the fireball's heat, all stranded
    radiation_bound: float  # argon's loss times the atom fraction times the allowance


def alumina_rows() -> list[AluminaRow]:
    out = []
    heat_per_kg = 0.5 * V * V * K / (1.0 + K) / (1.0 + K)  # [J/kg] of the merged mass
    for um in (50.0, 100.0, 200.0):
        areal = ALUMINA_RHO * um * 1e-6
        share_disk = areal / THICK
        pulse_mass = 1.0 + K  # per kg of disk
        coat = share_disk  # per kg of disk
        atoms_coat = 5.0 * coat / ALUMINA_M
        atoms_rest = 1.0 / C_M + 3.0 * K / CH2_M
        frac = atoms_coat / (atoms_coat + atoms_rest)
        toll = coat * (ALUMINA_ATOMIZE + ALUMINA_ION) / ALUMINA_M
        out.append(
            AluminaRow(
                coat_um=um,
                coat_areal=areal,
                share_of_disk=share_disk,
                share_of_pulse_mass=coat / pulse_mass,
                atom_fraction=frac,
                toll_bound=toll / (heat_per_kg * pulse_mass),
                radiation_bound=ARGON_SUN_LOSS * frac * RADIATION_ALLOWANCE,
            )
        )
    return out


# --- Main ------------------------------------------------------------------------------------

OUT_COUPLING = OUT_DIR / "thick_disk_coupling.csv"
OUT_BOWL = OUT_DIR / "thick_disk_bowl.csv"
OUT_OFFSET = OUT_DIR / "thick_disk_offset.csv"
OUT_SPLIT = OUT_DIR / "thick_disk_split.csv"
OUT_ALUMINA = OUT_DIR / "thick_disk_alumina.csv"
OUT_CONVERGENCE = OUT_DIR / "thick_disk_convergence.csv"
#: S62: the thin 1 kg disk's offset table, side impulse against offset (both disks), 4 R_sun.
OUT_THIN_OFFSET = OUT_DIR / "thin_disk_offset.csv"
OUT_SIDE = OUT_DIR / "thin_disk_side.csv"
OUT_WIDTH = OUT_DIR / "thin_disk_width.csv"
OUT_SPEED = OUT_DIR / "thin_disk_speed.csv"


def main() -> None:
    print(
        f"python: columns {column_length(THIN):.3f} m (8 kg/m^2) and {column_length(THICK):.3f} m"
        f" (22 kg/m^2); 2.5 kg at 1 cm is {disk_radius(2.5, THICK) * 200:.1f} cm across;"
        f" energy {0.5 * 2.5 * V * V / 1e9:.0f} GJ against the parent's"
        f" {0.5 * PARENT_DISK * PARENT_V**2 / 1e9:.0f} GJ"
    )
    coupling = coupling_rows()
    for c in coupling:
        print(
            f"  {c.grid:>4} {c.case:<36} L/R {c.length_over_radius:5.2f}"
            f" cells {c.cells_across:4.1f}"
            f" J/mv {c.impulse_per_mv:5.3f} /point {c.over_point:5.3f} /base {c.over_baseline:5.3f}"
            f" dv {c.dv_per_pulse:4.2f} m/s peak {c.peak_gpa:7.2f} GPa at {c.peak_radius:4.2f} m"
            f" ({c.peak_over_baseline:4.2f} of thin)"
        )
    convergence = convergence_rows(coupling)
    for v in convergence:
        print(
            f"  converge {v.case:<36} {v.coarse:6.3f} {v.production:6.3f} {v.fine:6.3f}"
            f" (last {v.change_last:+.1%})"
        )
    for c in coupling:
        if c.kind == "wide":
            print(
                f"  wide gain {c.grid:>4} {c.case:<36} over own matched {c.over_own_matched:5.3f}"
            )
    bowls = bowl_rows(coupling)
    for b in bowls:
        print(
            f"  bowl {b.pulse_case:<26} {b.theta:>5}: J {b.pulse_ns / 1e6:5.2f} MN s,"
            f" mean {b.mean_force_mn:4.2f} MN, free 5 t kick {b.free_kick:4.0f} m/s,"
            f" stroke at 5 t {b.stroke_5t:5.1f} m, bowl at 8 m {b.bowl_t_8m:5.1f} t,"
            f" ship dv {b.dv_per_pulse:4.2f}"
        )
    offsets = offset_rows(coupling)
    for o in offsets:
        print(
            f"  offset {o.charge} (centred {o.centred_over_matched:5.3f}, margin"
            f" {o.margin * 100:4.1f} cm, bore clears {o.bore_clearance * 100:4.1f} cm)"
            f" -{o.drop:.0%}:"
            f" {o.offset_low * 100:5.2f}-{o.offset_high * 100:5.2f} cm"
            f" ({o.offset_low_over_rc:5.3f}-{o.offset_high_over_rc:5.3f} R_c),"
            f" usable {o.tolerance * 100:5.2f} cm,"
            f" side {o.side_impulse_low:6.4f} of axial, {o.side_ns_low / 1e3:5.1f} kN s"
        )
    widths = width_rows(coupling)
    for w in widths:
        print(
            f"  width {w.disk:<18} x{w.wide:4.2f} dr {w.dr:6.4f}: on grid x{w.wide_on_grid:5.3f}"
            f" L/R {w.length_over_radius_on_grid:5.3f} over own matched {w.over_own_matched:5.3f}"
        )
    thin_offsets = offset_rows(coupling, 1.0, THIN)
    for o in thin_offsets:
        print(
            f"  thin offset {o.charge} (centred {o.centred_over_matched:5.3f}, bore {o.bore:5.3f},"
            f" clears {o.bore_clearance * 100:4.2f} cm) -{o.drop:.0%}:"
            f" {o.offset_low * 100:5.2f}-{o.offset_high * 100:5.2f} cm,"
            f" usable {o.tolerance * 100:5.2f} cm,"
            f" side {o.side_impulse_low:6.4f} of axial, {o.side_ns_low / 1e3:5.1f} kN s"
        )
    sides = side_rows(coupling, 1.0, THIN) + side_rows(coupling, 2.5, THICK)
    for d in sides:
        print(
            f"  side {d.disk:<16} {d.charge} d {d.offset * 100:4.1f} cm:"
            f" struck {d.struck_share:5.3f}"
            f" J {d.low:5.3f}-{d.high:5.3f} (over matched {d.over_matched_low:5.3f})"
            f" lateral {d.lateral:6.4f} slope {d.slope:5.3f}/m off-line {d.off_line:+6.1%}"
            f" {d.side_ns / 1e3:5.2f} kN s ({d.side_ns_4rsun / 1e3:5.2f} at 4 R_sun)"
        )
    speeds = speed_rows(coupling)
    for p in speeds:
        print(
            f"  speed {p.depth:<10} {p.closing / 1e3:5.1f} km/s {p.case:<34} x{p.over_matched:5.3f}"
            f" J {p.pulse_ns / 1e6:5.3f} MN s E {p.energy_gj:5.1f} GJ peak x{p.peak_over_392:4.2f}"
            f" clears {p.clearance * 100:4.2f} cm side {p.side_ns_at_clearance / 1e3:5.2f} kN s"
            f" PE {p.pe_velocity_factor:6.4f} bowl {p.bowl_t_8m_2hz} t"
        )
    splits = split_rows()
    for s in splits:
        print(
            f"  split {s.charge} gap {s.gap_cm:4.1f} cm: missed {s.missed_share:5.3f},"
            f" impulse {s.low:5.3f}-{s.high:5.3f}, axial delay/sweep {s.axial_delay_over_sweep:.4f}"
        )
    alumina = alumina_rows()
    for a in alumina:
        print(
            f"  alumina {a.coat_um:4.0f} um: {a.share_of_disk:.2%} of disk,"
            f" {a.share_of_pulse_mass:.3%} of pulse, atoms {a.atom_fraction:.2e},"
            f" toll <= {a.toll_bound:.1e},"
            f" radiation <= {a.radiation_bound:.1e}"
        )
    for path, rows in (
        (OUT_COUPLING, coupling),
        (OUT_BOWL, bowls),
        (OUT_OFFSET, offsets),
        (OUT_SPLIT, splits),
        (OUT_ALUMINA, alumina),
        (OUT_CONVERGENCE, convergence),
        (OUT_THIN_OFFSET, thin_offsets),
        (OUT_WIDTH, widths),
        (OUT_SIDE, sides),
        (OUT_SPEED, speeds),
    ):
        if rows:
            ib._write(path, list(rows))
    print(f"python: wrote {OUT_DIR}/thick_disk_*.csv")


if __name__ == "__main__":
    main()
