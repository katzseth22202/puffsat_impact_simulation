"""Result 21: the shaped plug spreads the bag, the gun is the backup (parent S52).

    PYTHONPATH=python uv run python -m puffsat.xray_collision.shaped_plug

Result 20 (`iapetus_bag.py`) put the ship's ice charge, the parent's **plug**, in the bag's own
footprint, a column 1.84 m long at k = 8.5. The bag ploughs along it and bursts along a line, which
a paraboloid focuses poorly. Seth (2026-10-21) asked for a plug shaped to spread the bag, with the
gun's bullet fired only on a predicted miss. The parent's sizing is a squat lens about the bag's
own depth (Birkhoff's penetration law: water stops in ~1.04 of its length of ice). A lens wider than
the bag also catches a bag that strays sideways.

The pieces, all on Result 20's delta/D 0.5 paraboloid with the plug at 2 F unless stated:

1. **The lens sweep** (`lens_rows`, S52 items 1 and 7). `sweep --shaped-plug` runs ice lenses
   whose rim is 0.25-2x the bag's 0.363 m depth at k = 1-15, Result 20's compact column at the
   same k as the control, the lens at 1.5 F and 2.5 F, a dished lens and the point fireball at
   k = 1 and 3. Each run goes through Result 20's water corrections (`iapetus_bag.sweep_case`:
   gamma, rebound, frozen-rebound edge) and its plate grading (`plate_case`).
2. **Punch-through** (`punch_rows`, item 2). The thinnest lens that still holds the face below
   the stopped lenses' level: past it the bag's core reaches the plate as a jet.
3. **The safe overlap** (`overlap_rows`, item 3). The 2-D code is axisymmetric, so an off-centre
   strike is stood in for by an open core through the lens. The core lets a share of the bag's
   face through as a compact jet on the axis, where the axis focuses it; that makes it an upper
   bound on the face peak of the real crescent beside the plug. The largest share that holds
   maraging's 2.3 GPa, put through the overlap geometry of two circles (`passed_share`), gives
   the **catch radius**, and a two-axis Gaussian track error turns that into the backup's rate.
4. **The backup pulse** (`backup_rows`, items 4 and 5). The bullet's cloud onto a centred lens,
   and its misses met nearer and further up-track than Result 20's 31 m. The plug's share of the
   spread cloud, against ice's vaporization energy and its station-keeping, says whether it
   survives. The gun's lead time follows from the meeting distance and its muzzle speed.
5. **The plate** (item 6). Result 20's graded floor and graphite liner, for the lens's load.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from puffsat.xray_collision import charge
from puffsat.xray_collision import iapetus_bag as ib

# --- The case -------------------------------------------------------------------------------

OUT_DIR = charge.OUT_DIR
LENS = OUT_DIR / "shaped_plug.jsonl"
LENS_COARSE = OUT_DIR / "shaped_plug_coarse.jsonl"
POINT = OUT_DIR / "shaped_plug_point.jsonl"
POINT_COARSE = OUT_DIR / "shaped_plug_point_coarse.jsonl"
OVERLAP = OUT_DIR / "shaped_plug_overlap.jsonl"
OVERLAP_COARSE = OUT_DIR / "shaped_plug_overlap_coarse.jsonl"
OVERLAP_FINE = OUT_DIR / "shaped_plug_overlap_fine.jsonl"
BACKUP = OUT_DIR / "shaped_plug_backup.jsonl"
BACKUP_COARSE = OUT_DIR / "shaped_plug_backup_coarse.jsonl"
BACKUP_FAR = OUT_DIR / "shaped_plug_backup_far.jsonl"
BACKUP_FAR_COARSE = OUT_DIR / "shaped_plug_backup_far_coarse.jsonl"

#: The paraboloid Result 20 recommends, and its plug height.
PLATE = "paraboloid 0.50 at 2 F"
FOCUS = ib.PLATE_R / (8.0 * 0.5)  # [m]
RIM_HEIGHT = 0.5 * 2.0 * ib.PLATE_R  # [m] the 0.5 bowl's rim above its vertex
#: The intact bag: a 25 kg water sphere (the parent's 0.363 m), and its radius in the plane.
BAG_DEPTH = ib.sphere_diameter(ib.BAG, ib.WATER_RHO)
BAG_RADIUS = 0.5 * BAG_DEPTH
#: Ice's energy to go from ~150 K to steam [J/kg]: ~0.28 MJ/kg to the melt, 0.33 latent, 0.42 to
#: 373 K, 2.26 latent.
ICE_VAPORIZE = 3.3e6
#: The plug's last steering before the bag arrives (the parent's ~100 ms).
PLUG_STEER = 0.1  # [s]
#: Muzzle speeds [m/s] for the trigger window: Result 20's ~6 km firing range implies ~230 m/s.
MUZZLE = (230.0, 500.0, 1000.0, 2000.0, 3000.0)
#: One-sigma track errors per axis [m] at the plug, for the backup's rate.
TRACK_SIGMA = (0.01, 0.02, 0.03, 0.05, 0.1)
#: The backup rates the parent's ledger can carry: its 5% placeholder, and 1%.
BACKUP_RATES = (0.05, 0.01)
#: Halving the cell raised the centred share-2 lens's peak 0.69 -> 0.85 GPa (+23%), so a meeting
#: distance sets the trigger only if its production-grid miss holds 2.3 GPa by this margin.
PEAK_MARGIN = 1.2


def _f(r: dict[str, object], key: str) -> float:
    return float(str(r[key]))


@cache
def _runs(path: Path) -> tuple[dict[str, object], ...]:
    return ib.runs(path)


def as_sweep_record(r: dict[str, object]) -> dict[str, object]:
    """A shaped-plug run in Result 20's schema, so `iapetus_bag.sweep_case` reads it: the plug's
    radius stands in the sheet's column and its depth in the charge length."""
    return {
        **r,
        "case": PLATE,
        "charge": "ice",
        "gamma": 1.4,
        "k": _f(r, "k") + _f(r, "bullet_k"),
        "sheet_radius": _f(r, "plug_radius"),
        "charge_length": _f(r, "plug_depth"),
        "height": _f(r, "height"),
    }


# --- 1. The lens sweep (items 1 and 7) ------------------------------------------------------


@dataclass(frozen=True)
class LensRow:
    kind: str  # lens, compact (Result 20's column), dish, point (the ideal)
    k: float
    depth_share: float  # rim depth over the bag's 0.363 m
    column_ratio: float  # the plug's column over the code bag's
    plug_radius: float
    plug_depth: float
    height_focal: float
    impulse_per_mv: float  # B: Result 20's water corrections applied
    eta: float
    eta_frozen: float
    leverage: float
    leverage_frozen: float
    peak_gpa: float
    peak_radius: float
    exit_direction: float
    ringing_uniform: float
    ringing_graded: float
    graded_t: float
    graded_centre_steel_mm: float
    catch_radius_geometric: float  # plug minus bag radius: no part of the bag passes


def _lens_row(r: dict[str, object]) -> LensRow:
    rec = as_sweep_record(r)
    s = ib.sweep_case(rec, "ice")
    p = ib.plate_case(s, rec)
    return LensRow(
        kind=str(r["kind"]),
        k=s.k,
        depth_share=_f(r, "depth_share"),
        column_ratio=_f(r, "column_ratio"),
        plug_radius=_f(r, "plug_radius"),
        plug_depth=_f(r, "plug_depth"),
        height_focal=_f(r, "height_focal"),
        impulse_per_mv=s.impulse_per_mv,
        eta=s.eta,
        eta_frozen=s.eta_frozen,
        leverage=s.leverage,
        leverage_frozen=s.leverage_frozen,
        peak_gpa=s.peak_pa / 1e9,
        peak_radius=_f(r, "peak_radius"),
        exit_direction=s.exit_direction,
        ringing_uniform=p.ringing_uniform,
        ringing_graded=p.ringing_graded,
        graded_t=p.graded_t,
        graded_centre_steel_mm=p.graded_centre_steel_mm,
        catch_radius_geometric=max(_f(r, "plug_radius") - BAG_RADIUS, 0.0),
    )


def _point_row(r: dict[str, object]) -> LensRow:
    s = ib.sweep_case(r, "ice")
    p = ib.plate_case(s, r)
    return LensRow(
        kind="point",
        k=s.k,
        depth_share=0.0,
        column_ratio=0.0,
        plug_radius=0.0,
        plug_depth=0.0,
        height_focal=2.0,
        impulse_per_mv=s.impulse_per_mv,
        eta=s.eta,
        eta_frozen=s.eta_frozen,
        leverage=s.leverage,
        leverage_frozen=s.leverage_frozen,
        peak_gpa=s.peak_pa / 1e9,
        peak_radius=_f(r, "peak_radius"),
        exit_direction=s.exit_direction,
        ringing_uniform=p.ringing_uniform,
        ringing_graded=p.ringing_graded,
        graded_t=p.graded_t,
        graded_centre_steel_mm=p.graded_centre_steel_mm,
        catch_radius_geometric=math.nan,
    )


def lens_rows(fine: bool = True) -> list[LensRow]:
    """Every lens, column and dish run, and the point fireball at k = 1-15 (this result's k = 1
    and 3 with Result 20's 2, 4, 8.5 and 15)."""
    out = [_lens_row(r) for r in _runs(LENS if fine else LENS_COARSE)]
    points = [r for r in _runs(POINT if fine else POINT_COARSE)]
    points += [
        r
        for r in ib.runs(ib.JSONL if fine else ib.JSONL_COARSE)
        if r["kind"] == "point"
        and r["case"] == PLATE
        and _f(r, "gamma") == 1.4
        and _f(r, "k") <= 15.0
    ]
    out += [_point_row(r) for r in points]
    return sorted(out, key=lambda x: (x.kind, x.height_focal, x.k, x.depth_share))


@dataclass(frozen=True)
class BestLensRow:
    k: float
    best_share: float  # the lens depth share with the most leverage at 2 F, inside 2.3 GPa
    leverage: float
    leverage_frozen: float
    eta: float
    peak_gpa: float
    plug_radius: float
    compact_leverage: float
    compact_leverage_frozen: float
    point_leverage: float
    point_leverage_frozen: float
    gain_over_compact: float
    share_of_ideal: float  # (lens - compact) / (point - compact)


def best_lens_rows(rows: list[LensRow], max_share: float = 99.0) -> list[BestLensRow]:
    """Per k, the 2 F lens with the most leverage among those whose face peak holds maraging's
    2.3 GPa (`max_share` excludes nothing by default). A k with no such lens is left out."""
    out = []
    for k in sorted({r.k for r in rows if r.kind == "lens"}):
        lenses = [
            r
            for r in rows
            if r.kind == "lens"
            and r.k == k
            and r.height_focal == 2.0
            and r.depth_share <= max_share
            and r.peak_gpa <= ib.MARAGING / 1e9
        ]
        if not lenses:
            continue
        best = max(lenses, key=lambda r: r.leverage)
        comp = [r for r in rows if r.kind == "compact" and r.k == k]
        pt = [r for r in rows if r.kind == "point" and r.k == k]
        c = comp[0] if comp else None
        p = pt[0] if pt else None
        cl = c.leverage if c else math.nan
        pl = p.leverage if p else math.nan
        out.append(
            BestLensRow(
                k=k,
                best_share=best.depth_share,
                leverage=best.leverage,
                leverage_frozen=best.leverage_frozen,
                eta=best.eta,
                peak_gpa=best.peak_gpa,
                plug_radius=best.plug_radius,
                compact_leverage=cl,
                compact_leverage_frozen=c.leverage_frozen if c else math.nan,
                point_leverage=pl,
                point_leverage_frozen=p.leverage_frozen if p else math.nan,
                gain_over_compact=best.leverage / cl if c else math.nan,
                share_of_ideal=(
                    (best.leverage - cl) / (pl - cl) if c and p and pl != cl else math.nan
                ),
            )
        )
    return out


@dataclass(frozen=True)
class PlateLensRow:
    """The Earth bag plate for one plug (item 6): Result 20's spring law and grading, and the
    graphite liner it erodes, a pulse and over the 2000-pulse push."""

    kind: str
    k: float
    depth_share: float
    pulse_ns: float
    ship_t_1g: float
    spring_t: float
    centre_over_mean: float
    ringing_uniform: float
    ringing_graded: float
    graded_t: float
    graded_centre_steel_mm: float
    inner_energy_mj_m2: float
    graphite_kg_per_pulse: tuple[float, float]
    graphite_t_per_push: tuple[float, float]
    leverage: float
    leverage_with_graphite: float
    heat_valid: bool


def plate_rows(kinds: tuple[str, ...] = ("lens", "compact")) -> list[PlateLensRow]:
    out = []
    for r in _runs(LENS):
        if r["kind"] not in kinds or _f(r, "height_focal") != 2.0:
            continue
        rec = as_sweep_record(r)
        s = ib.sweep_case(rec, "ice")
        p = ib.plate_case(s, rec)
        e = ib.erosion_case(s, rec)
        g = e.graphite_kg_per_pulse
        out.append(
            PlateLensRow(
                kind=str(r["kind"]),
                k=s.k,
                depth_share=_f(r, "depth_share"),
                pulse_ns=p.pulse_ns,
                ship_t_1g=p.ship_t_1g,
                spring_t=p.spring_t,
                centre_over_mean=p.centre_over_mean,
                ringing_uniform=p.ringing_uniform,
                ringing_graded=p.ringing_graded,
                graded_t=p.graded_t,
                graded_centre_steel_mm=p.graded_centre_steel_mm,
                inner_energy_mj_m2=e.inner_energy_mj_m2,
                graphite_kg_per_pulse=g,
                graphite_t_per_push=(g[0] * ib.PULSES / 1e3, g[1] * ib.PULSES / 1e3),
                leverage=s.leverage,
                leverage_with_graphite=e.leverage_with_graphite,
                heat_valid=e.heat_valid,
            )
        )
    return sorted(out, key=lambda x: (x.kind, x.k, x.depth_share))


# --- 2. Punch-through (item 2) --------------------------------------------------------------

#: A lens has stopped the bag if its face peak stays inside maraging's 2.3 GPa. A thinner lens
#: makes a thin hot pancake that blows along its axis, and the jet lands on the vertex: the peak
#: sits on the axis and rises as the lens thins.


@dataclass(frozen=True)
class PunchRow:
    k: float
    thinnest_stopped_share: float  # rim depth over the bag's depth; nan if the deepest fails
    thinnest_column_ratio: float
    plug_radius: float  # the widest plug that stops the bag
    catch_radius_geometric: float
    peak_gpa_stopped: float  # at the thinnest stopped lens
    next_thinner_share: float  # the first lens past the limit (nan if none ran)
    peak_gpa_past: float
    peak_gpa_deep: float  # the deepest lens's, the reference


def punch_rows(rows: list[LensRow]) -> list[PunchRow]:
    out = []
    for k in sorted({r.k for r in rows if r.kind == "lens"}):
        lenses = sorted(
            (r for r in rows if r.kind == "lens" and r.k == k and r.height_focal == 2.0),
            key=lambda r: r.depth_share,
        )
        deep = lenses[-1]
        limit = ib.MARAGING / 1e9
        # Walk down from the deepest: the last lens still under the limit, with every deeper one.
        stopped: LensRow | None = None
        past: LensRow | None = None
        for r in reversed(lenses):
            if r.peak_gpa <= limit:
                stopped = r
            else:
                past = r
                break
        if stopped is None:
            # Even the deepest lens run sends a jet over the limit: no depth reliably stops it.
            out.append(
                PunchRow(
                    k=k,
                    thinnest_stopped_share=math.nan,
                    thinnest_column_ratio=math.nan,
                    plug_radius=math.nan,
                    catch_radius_geometric=math.nan,
                    peak_gpa_stopped=math.nan,
                    next_thinner_share=deep.depth_share,
                    peak_gpa_past=deep.peak_gpa,
                    peak_gpa_deep=deep.peak_gpa,
                )
            )
            continue
        out.append(
            PunchRow(
                k=k,
                thinnest_stopped_share=stopped.depth_share,
                thinnest_column_ratio=stopped.column_ratio,
                plug_radius=stopped.plug_radius,
                catch_radius_geometric=stopped.catch_radius_geometric,
                peak_gpa_stopped=stopped.peak_gpa,
                next_thinner_share=past.depth_share if past else math.nan,
                peak_gpa_past=past.peak_gpa if past else math.nan,
                peak_gpa_deep=deep.peak_gpa,
            )
        )
    return out


# --- 3. The safe overlap (item 3) -----------------------------------------------------------


def passed_share(offset: float, plug_radius: float, bag_radius: float = BAG_RADIUS) -> float:
    """The share of a bag's face (a disk of `bag_radius`) that misses a plug face (a disk of
    `plug_radius`) when their centres are `offset` apart: one minus the two circles' lens-shaped
    intersection over the bag's area."""
    a, b, d = bag_radius, plug_radius, offset
    if d >= a + b:
        return 1.0
    if d <= abs(b - a):
        return 0.0 if b >= a else 1.0 - (b / a) ** 2
    alpha = math.acos((d * d + a * a - b * b) / (2.0 * d * a))
    beta = math.acos((d * d + b * b - a * a) / (2.0 * d * b))
    inter = (
        a * a * alpha
        + b * b * beta
        - 0.5 * math.sqrt((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b))
    )
    return 1.0 - inter / (math.pi * a * a)


def catch_radius(plug_radius: float, safe_share: float, bag_radius: float = BAG_RADIUS) -> float:
    """[m] The largest offset of the bag's centre whose passed share stays at `safe_share`."""
    lo, hi = 0.0, plug_radius + bag_radius
    if passed_share(0.0, plug_radius, bag_radius) > safe_share:
        return 0.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if passed_share(mid, plug_radius, bag_radius) > safe_share:
            hi = mid
        else:
            lo = mid
    return lo


def miss_rate(catch: float, sigma: float) -> float:
    """The chance a bag's offset at the plug exceeds `catch` for an isotropic Gaussian error of
    `sigma` per axis (the Rayleigh tail)."""
    return math.exp(-0.5 * (catch / sigma) ** 2)


def sigma_for_rate(catch: float, rate: float) -> float:
    """[m] The per-axis track error at the plug whose Rayleigh tail past `catch` is `rate`."""
    return catch / math.sqrt(2.0 * math.log(1.0 / rate))


@dataclass(frozen=True)
class OverlapRow:
    grid: str
    kind: str  # hole (an open core passes the bag's middle) or rim (a narrow plug passes its rim)
    k_lens: float  # the solid lens this plug stands in for
    k: float  # the plug run's own
    depth_share: float
    plug_radius: float
    hole: float
    passed_share: float  # of the bag's face, meeting no plug
    impulse_per_mv_lossless: float
    peak_gpa: float
    peak_radius: float
    holds_maraging: bool


def overlap_rows() -> list[OverlapRow]:
    out = []
    sources = (("x1", OVERLAP), ("x0.5", OVERLAP_COARSE), ("x2", OVERLAP_FINE))
    for grid, path in sources:
        for r in _runs(path):
            peak = _f(r, "peak_per_mv2") * ib.BAG * ib.W_REF**2
            out.append(
                OverlapRow(
                    grid=grid,
                    kind=str(r["kind"]),
                    k_lens=_f(r, "k_lens"),
                    k=_f(r, "k"),
                    depth_share=_f(r, "depth_share"),
                    plug_radius=_f(r, "plug_radius"),
                    hole=_f(r, "hole"),
                    passed_share=_f(r, "passed_share"),
                    impulse_per_mv_lossless=_f(r, "impulse_per_mv"),
                    peak_gpa=peak / 1e9,
                    peak_radius=_f(r, "peak_radius"),
                    holds_maraging=peak <= ib.MARAGING,
                )
            )
    return out


def safe_share(rows: list[OverlapRow], solid_peak_gpa: float) -> float:
    """The passed share at which the face peak reaches 2.3 GPa, interpolated linearly between the
    solid lens (share 0) and the runs of one bracket of one lens on one grid (`rows`, already
    filtered). The grid's own solid run is the reference where it ran one, else
    `solid_peak_gpa`. 0 if the solid lens already exceeds the limit; held at the largest run if
    every one holds."""
    own = [r for r in rows if r.passed_share == 0.0]
    pts = sorted(
        [(0.0, own[0].peak_gpa if own else solid_peak_gpa)]
        + [(r.passed_share, r.peak_gpa) for r in rows if r.passed_share > 0.0]
    )
    lim = ib.MARAGING / 1e9
    if pts[0][1] > lim:
        return 0.0
    for (s0, p0), (s1, p1) in itertools.pairwise(pts):
        if p0 <= lim < p1:
            return s0 + (lim - p0) / (p1 - p0) * (s1 - s0)
    return pts[-1][0]


@dataclass(frozen=True)
class CatchRow:
    grid: str
    k: float
    depth_share: float
    plug_radius: float
    safe_share_core: float  # the open-core bracket
    safe_share_rim: float  # the narrow-plug bracket
    safe_share: float  # the lesser, which sets the catch
    catch_radius: float  # [m] offset of the bag's centre the plug takes; 0 if centred fails
    catch_radius_geometric: float  # with no part of the bag allowed past
    miss_rates: tuple[float, ...]  # at TRACK_SIGMA
    sigma_for_rates: tuple[float, ...]  # [m] the track error that holds each of BACKUP_RATES


def catch_rows(lenses: list[LensRow], overlaps: list[OverlapRow]) -> list[CatchRow]:
    out = []
    for grid, k, share in sorted({(o.grid, o.k_lens, o.depth_share) for o in overlaps}):
        lens = [
            r
            for r in lenses
            if r.kind == "lens"
            and r.height_focal == 2.0
            and abs(r.k - k) < 1e-6
            and r.depth_share == share
        ]
        if not lens:
            continue
        mine = [o for o in overlaps if (o.grid, o.k_lens, o.depth_share) == (grid, k, share)]
        solid = [o for o in mine if o.passed_share == 0.0]
        core = safe_share([o for o in mine if o.kind == "hole"], lens[0].peak_gpa)
        rim = safe_share(solid + [o for o in mine if o.kind == "rim"], lens[0].peak_gpa)
        s = min(core, rim)
        # A lens that exceeds the limit even when centred catches nothing.
        solid_ok = (solid[0].peak_gpa if solid else lens[0].peak_gpa) <= ib.MARAGING / 1e9
        c = catch_radius(lens[0].plug_radius, s) if solid_ok else 0.0
        out.append(
            CatchRow(
                grid=grid,
                k=k,
                depth_share=share,
                plug_radius=lens[0].plug_radius,
                safe_share_core=core,
                safe_share_rim=rim,
                safe_share=s,
                catch_radius=c,
                catch_radius_geometric=lens[0].catch_radius_geometric,
                miss_rates=tuple(miss_rate(c, sig) for sig in TRACK_SIGMA),
                sigma_for_rates=tuple(sigma_for_rate(c, q) for q in BACKUP_RATES),
            )
        )
    return out


@dataclass(frozen=True)
class GeometricCatchRow:
    """The catch with no part of the bag allowed past, for the recommended 2 F lens at each k:
    the trigger that needs no overlap credit, since the 2-D overlap bound does not converge."""

    k: float
    depth_share: float
    plug_radius: float
    peak_gpa: float
    catch_radius: float  # [m] plug radius less the bag's
    miss_rates: tuple[float, ...]  # at TRACK_SIGMA
    sigma_for_rates: tuple[float, ...]  # [m] at BACKUP_RATES


def geometric_catch_rows(lenses: list[LensRow], share: float = 2.0) -> list[GeometricCatchRow]:
    out = []
    for r in lenses:
        if r.kind != "lens" or r.height_focal != 2.0 or r.depth_share != share:
            continue
        c = r.catch_radius_geometric if r.peak_gpa <= ib.MARAGING / 1e9 else 0.0
        out.append(
            GeometricCatchRow(
                k=r.k,
                depth_share=share,
                plug_radius=r.plug_radius,
                peak_gpa=r.peak_gpa,
                catch_radius=c,
                miss_rates=tuple(miss_rate(c, sig) for sig in TRACK_SIGMA),
                sigma_for_rates=tuple(sigma_for_rate(c, q) for q in BACKUP_RATES),
            )
        )
    return sorted(out, key=lambda x: x.k)


# --- 4. The backup pulse and the trigger window (items 4 and 5) -----------------------------


@dataclass(frozen=True)
class BackupRow:
    grid: str
    kind: str  # backup (bullet, then a centred lens) or miss (bullet, no plug)
    k: float  # the plug's
    bullet_g: float
    meeting_m: float  # bullet meets the bag this far up-track of the plug
    impulse_mn_s: float  # lossless gamma gas
    impulse_per_mv: float  # with Result 20's water corrections, where a plug is in place
    leverage: float
    peak_gpa: float
    holds_maraging: bool
    inner_energy_mj_m2: float  # arriving inside 0.6 m (an upper estimate: see heat_valid)
    liner_graphite_kg: tuple[float, float]  # the face heat taken by eroding graphite
    liner_t_per_push: tuple[float, float]  # [t] at the first of BACKUP_RATES of the push's pulses
    heat_valid: bool


def _liner(rec: dict[str, object], k: float) -> tuple[float, tuple[float, float], bool]:
    """The inner ring's arriving energy and the graphite a pulse erodes, as Result 20's
    `erosion_case` charges it, with water's wall share at the merged gas's arrival speed."""
    bands = rec["bands"]
    assert isinstance(bands, list)
    w = ib.W_REF
    u = w * math.sqrt(k) / (1.0 + k) if k > 0.0 else 0.0
    v_arr = u + w / (1.0 + k)
    column = (1.0 + k) * ib.BAG / (2.0 * math.pi * ib.PLATE_R**2)
    _, wall = ib.rebound("ice", v_arr, column)
    t = FOCUS * 2.0 / max(u, v_arr)
    threshold = (
        ib.FACE_MELT_RISE * math.sqrt(math.pi * ib.FACE_K * ib.FACE_RHO * ib.FACE_CP * t) / 2.0
    )
    excess, arriving = 0.0, 0.0
    for b in bands:
        assert isinstance(b, dict)
        lo, hi = float(b["r_lo"]), min(float(b["r_hi"]), ib.PLATE_R)
        if hi <= lo:
            continue
        area = math.pi * (hi * hi - lo * lo)
        arriving += area * float(b["energy_in"])
        q = float(b["energy_in"]) * ib.BAG * w * w * wall
        excess += area * max(q - threshold, 0.0)
    inner = float(bands[0]["energy_in"]) * ib.BAG * w * w
    return inner, (excess / ib.Q_GRAPHITE[0], excess / ib.Q_GRAPHITE[1]), arriving / 0.5 <= 1.0


def backup_rows() -> list[BackupRow]:
    """This result's runs (the far backup's where it has run), and Result 20's 100 g miss at
    31 m on the production grid."""
    out = []
    srcs = [
        (grid, r)
        for grid, path in (
            ("x1", BACKUP),
            ("x0.5", BACKUP_COARSE),
            ("x1", BACKUP_FAR),
            ("x0.5", BACKUP_FAR_COARSE),
        )
        if path.exists()
        for r in _runs(path)
    ]
    srcs += [
        ("x1", {**r, "kind": "miss", "bullet_k": r["bumper_k"], "bullet_gap": r["bumper_gap"]})
        for r in _runs(ib.BULLET_KEY)
        if _f(r, "k") == 0.0 and _f(r, "bumper_radius") < _f(r, "bag_radius")
    ]
    for grid, r in srcs:
        k = _f(r, "k")
        kb = _f(r, "bullet_k")
        peak = _f(r, "peak_per_mv2") * ib.BAG * ib.W_REF**2
        if r["kind"] == "backup":
            s = ib.sweep_case(as_sweep_record(r), "ice")
            b, lev = s.impulse_per_mv, s.leverage
        else:
            b, lev = math.nan, math.nan
        inner, graphite, valid = _liner(r, k + kb)
        out.append(
            BackupRow(
                grid=grid,
                kind=str(r["kind"]),
                k=k,
                bullet_g=1e3 * kb * ib.BAG,
                meeting_m=_f(r, "bullet_gap"),
                impulse_mn_s=_f(r, "impulse_per_mv") * ib.BAG * ib.W_REF / 1e6,
                impulse_per_mv=b,
                leverage=lev,
                peak_gpa=peak / 1e9,
                holds_maraging=peak <= ib.MARAGING,
                inner_energy_mj_m2=inner / 1e6,
                liner_graphite_kg=graphite,
                liner_t_per_push=(
                    BACKUP_RATES[0] * ib.PULSES * graphite[0] / 1e3,
                    BACKUP_RATES[0] * ib.PULSES * graphite[1] / 1e3,
                ),
                heat_valid=valid,
            )
        )
    return sorted(out, key=lambda x: (x.kind, x.grid, x.bullet_g, x.meeting_m, x.k))


@dataclass(frozen=True)
class PlugFateRow:
    """The plug on a backup pulse: the bullet's cloud has spread to `cloud_radius` by the plug,
    and a plug inside it takes its share of the cloud's mass and energy (a uniform cloud, a
    centred hit: the upper end)."""

    k: float
    plug_radius: float
    cloud_radius: float
    intercepted_share: float
    energy_mj_per_kg: float  # into the plug, against ice's ~3.3 MJ/kg to vaporize
    kick_m_s: float  # the plug's velocity from the intercepted momentum
    spent: bool


def plug_fate(k: float, plug_radius: float, cloud_radius: float = 2.5) -> PlugFateRow:
    share = min(1.0, (plug_radius / cloud_radius) ** 2)
    mass = k * ib.BAG
    energy = share * 0.5 * ib.BAG * ib.W_REF**2 / mass
    kick = share * ib.BAG * ib.W_REF / mass
    # Spent if it vaporizes, or is knocked further off station than its steering can return in
    # the half second to the next bag (a few m/s of authority is generous).
    return PlugFateRow(
        k, plug_radius, cloud_radius, share, energy / 1e6, kick, energy > ICE_VAPORIZE or kick > 5.0
    )


@dataclass(frozen=True)
class TriggerRow:
    kind: str  # miss (no plug in the cloud) or backup (the centred 2x lens in it)
    k: float  # the plug's (0 for a miss)
    bullet_g: float
    meeting_m: float  # the shortest meeting distance whose pulse holds 2.3 GPa / PEAK_MARGIN
    muzzle_m_s: float
    path_m: float  # from the rim to the meeting point
    lead_ms: float  # the gun fires this long before the bag reaches the plug
    range_km: float  # the bag's distance when the gun fires
    after_plug_steer: bool  # fires inside the plug's last ~100 ms, so it sees the steered plug


def gun_path(meeting: float) -> float:
    """[m] From the paraboloid's rim gun to the meeting point on the axis, `meeting` up-track of
    the plug at 2 F."""
    return math.hypot(ib.PLATE_R, 2.0 * FOCUS + meeting - RIM_HEIGHT)


def trigger_rows(backups: list[BackupRow]) -> list[TriggerRow]:
    """For each bullet, with no plug in its cloud (a miss) and with the centred lens at each k (a
    backup), the shortest meeting distance whose pulse held 2.3 GPa with PEAK_MARGIN to spare
    (production grid), and the gun's lead time at each muzzle speed: `t = path / v_b + D / w`.
    A case none of whose meetings held is left out."""
    out = []
    prod = [b for b in backups if b.grid == "x1"]
    for kind, k, g in sorted({(b.kind, b.k, b.bullet_g) for b in prod}):
        held = [
            b.meeting_m
            for b in prod
            if (b.kind, b.k, b.bullet_g) == (kind, k, g)
            and b.peak_gpa * PEAK_MARGIN <= ib.MARAGING / 1e9
        ]
        if not held:
            continue
        d = min(held)
        for v in MUZZLE:
            lead = gun_path(d) / v + d / ib.W_REF
            out.append(
                TriggerRow(
                    kind,
                    k,
                    g,
                    d,
                    v,
                    gun_path(d),
                    1e3 * lead,
                    ib.W_REF * lead / 1e3,
                    lead < PLUG_STEER,
                )
            )
    return out


# --- Output ---------------------------------------------------------------------------------

OUT_LENS = OUT_DIR / "shaped_plug_lens.csv"
OUT_BEST = OUT_DIR / "shaped_plug_best.csv"
OUT_PUNCH = OUT_DIR / "shaped_plug_punch.csv"
OUT_PLATE = OUT_DIR / "shaped_plug_plate.csv"
OUT_OVERLAP = OUT_DIR / "shaped_plug_overlap.csv"
OUT_CATCH = OUT_DIR / "shaped_plug_catch.csv"
OUT_CATCH_GEOMETRIC = OUT_DIR / "shaped_plug_catch_geometric.csv"
OUT_BACKUP = OUT_DIR / "shaped_plug_backup.csv"
OUT_FATE = OUT_DIR / "shaped_plug_fate.csv"
OUT_TRIGGER = OUT_DIR / "shaped_plug_trigger.csv"
OUT_CONVERGENCE = OUT_DIR / "shaped_plug_convergence.csv"


@dataclass(frozen=True)
class ConvergenceRow:
    kind: str
    k: float
    depth_share: float
    height_focal: float
    impulse_per_mv: float
    impulse_per_mv_coarse: float
    peak_gpa: float
    peak_gpa_coarse: float


def convergence_rows(fine: list[LensRow], coarse: list[LensRow]) -> list[ConvergenceRow]:
    out = []
    for c in coarse:
        m = [
            f
            for f in fine
            if (f.kind, f.k, f.depth_share, f.height_focal)
            == (c.kind, c.k, c.depth_share, c.height_focal)
        ]
        if m:
            f = m[0]
            out.append(
                ConvergenceRow(
                    c.kind,
                    c.k,
                    c.depth_share,
                    c.height_focal,
                    f.impulse_per_mv,
                    c.impulse_per_mv,
                    f.peak_gpa,
                    c.peak_gpa,
                )
            )
    return out


def main() -> None:
    lenses = lens_rows()
    for r in lenses:
        print(
            f"  {r.kind:<7} k {r.k:5.2f} share {r.depth_share:4.2f} (col {r.column_ratio:4.2f})"
            f" R {r.plug_radius:4.2f} m h {r.height_focal:3.1f} F | B {r.impulse_per_mv:5.3f}"
            f" eta {r.eta:.3f} ({r.eta_frozen:.3f}) lev {r.leverage:5.2f}"
            f" ({r.leverage_frozen:5.2f}) | peak {r.peak_gpa:5.2f} GPa at r {r.peak_radius:4.2f}"
            f" | dir {r.exit_direction:.3f} | ring even {r.ringing_uniform:5.1%} graded"
            f" {r.ringing_graded:5.1%} ({r.graded_t:4.1f} t, centre {r.graded_centre_steel_mm:4.0f}"
            " mm)"
        )
    punch = punch_rows(lenses)
    for p in punch:
        print(
            f"  punch k {p.k:5.2f}: thinnest stopped share {p.thinnest_stopped_share:4.2f}"
            f" (column {p.thinnest_column_ratio:4.2f}, R {p.plug_radius:4.2f} m, geometric catch"
            f" {p.catch_radius_geometric:4.2f} m) peak {p.peak_gpa_stopped:5.2f} GPa | next"
            f" {p.next_thinner_share:4.2f}: {p.peak_gpa_past:5.2f} GPa | deep"
            f" {p.peak_gpa_deep:5.2f}"
        )
    best = best_lens_rows(lenses)
    for bl in best:
        print(
            f"  best k {bl.k:5.2f}: share {bl.best_share:4.2f} R {bl.plug_radius:4.2f} m ->"
            f" {bl.leverage:5.2f}x ({bl.leverage_frozen:5.2f}) eta {bl.eta:.3f} peak"
            f" {bl.peak_gpa:5.2f} GPa | compact {bl.compact_leverage:5.2f}x"
            f" ({bl.compact_leverage_frozen:5.2f}) point {bl.point_leverage:5.2f}x"
            f" ({bl.point_leverage_frozen:5.2f}) | x{bl.gain_over_compact:.3f},"
            f" {bl.share_of_ideal:.0%} of the way to the point"
        )
    plates = plate_rows()
    for pr in plates:
        print(
            f"  plate {pr.kind:<7} k {pr.k:5.2f} share {pr.depth_share:4.2f}: J"
            f" {pr.pulse_ns / 1e6:4.2f} MN s, ship {pr.ship_t_1g:4.0f} t | spring"
            f" {pr.spring_t:4.1f} t, centre x{pr.centre_over_mean:4.1f} | ring even"
            f" {pr.ringing_uniform:5.1%} graded {pr.ringing_graded:5.1%} ({pr.graded_t:4.1f} t,"
            f" {pr.graded_centre_steel_mm:4.0f} mm)"
            f" | inner {pr.inner_energy_mj_m2:5.2f} MJ/m^2, graphite"
            f" {pr.graphite_kg_per_pulse[0]:4.1f}-{pr.graphite_kg_per_pulse[1]:4.1f} kg/pulse"
            f" ({pr.graphite_t_per_push[0]:4.1f}-{pr.graphite_t_per_push[1]:4.1f} t/push) | lev"
            f" {pr.leverage:5.2f} -> {pr.leverage_with_graphite:5.2f}"
            f"{'' if pr.heat_valid else ' (heat over-counted)'}"
        )
    overlaps = overlap_rows()
    for o in overlaps:
        print(
            f"  overlap {o.grid:>4} {o.kind:<4} lens k {o.k_lens:4.1f} (k {o.k:5.2f}) share"
            f" {o.depth_share:4.2f} R {o.plug_radius:4.2f}"
            f" hole {o.hole:5.3f} passes {o.passed_share:5.1%} | J/mv"
            f" {o.impulse_per_mv_lossless:5.3f} peak {o.peak_gpa:6.2f} GPa at r"
            f" {o.peak_radius:4.2f} | {'holds' if o.holds_maraging else 'exceeds'}"
        )
    catches = catch_rows(lenses, overlaps)
    for c in catches:
        rates = ", ".join(
            f"{1e2 * s:.0f} cm {m:.1e}" for s, m in zip(TRACK_SIGMA, c.miss_rates, strict=True)
        )
        print(
            f"  catch {c.grid:>4} k {c.k:5.2f} share {c.depth_share:4.2f} R"
            f" {c.plug_radius:4.2f}: safe share core {c.safe_share_core:5.1%} rim"
            f" {c.safe_share_rim:5.1%} -> catch {c.catch_radius:4.2f} m (geometric"
            f" {c.catch_radius_geometric:4.2f}) | backup rate {rates} | sigma for 5% / 1%"
            f" {1e2 * c.sigma_for_rates[0]:.1f} / {1e2 * c.sigma_for_rates[1]:.1f} cm"
        )
    geometric = geometric_catch_rows(lenses)
    for gc in geometric:
        print(
            f"  catch geometric k {gc.k:5.2f} share {gc.depth_share:4.2f} R {gc.plug_radius:4.2f}"
            f" peak {gc.peak_gpa:4.2f} GPa: catch {gc.catch_radius:4.2f} m | sigma for 5% / 1%"
            f" {1e2 * gc.sigma_for_rates[0]:.1f} / {1e2 * gc.sigma_for_rates[1]:.1f} cm"
        )
    backups = backup_rows()
    for b in backups:
        print(
            f"  {b.kind:<6} {b.grid:>4} k {b.k:5.2f} bullet {b.bullet_g:5.0f} g at"
            f" {b.meeting_m:5.1f} m | {b.impulse_mn_s:5.2f} MN s, B {b.impulse_per_mv:5.3f},"
            f" lev {b.leverage:5.2f} | peak {b.peak_gpa:5.2f} GPa"
            f" {'holds' if b.holds_maraging else 'exceeds'} | inner {b.inner_energy_mj_m2:6.1f}"
            f" MJ/m^2, graphite {b.liner_graphite_kg[0]:5.1f}-{b.liner_graphite_kg[1]:5.1f} kg"
            f" ({'valid' if b.heat_valid else 'over-counted'}), at a 5% rate"
            f" {b.liner_t_per_push[0]:4.1f}-{b.liner_t_per_push[1]:4.1f} t a push"
        )
    fates = [
        plug_fate(r.k, r.plug_radius)
        for r in lenses
        if r.kind == "lens" and r.depth_share == 2.0 and r.height_focal == 2.0
    ]
    for f in fates:
        print(
            f"  plug fate k {f.k:5.2f} R {f.plug_radius:4.2f}: takes {f.intercepted_share:5.1%} of"
            f" a {f.cloud_radius:.1f} m cloud, {f.energy_mj_per_kg:5.1f} MJ/kg, kick"
            f" {f.kick_m_s:6.1f} m/s -> {'spent' if f.spent else 'survives'}"
        )
    triggers = trigger_rows(backups)
    for t in triggers:
        print(
            f"  trigger {t.kind:<6} k {t.k:4.1f} {t.bullet_g:4.0f} g meets at {t.meeting_m:5.1f} m,"
            " muzzle"
            f" {t.muzzle_m_s:5.0f} m/s: path {t.path_m:5.1f} m, fires {t.lead_ms:6.1f} ms /"
            f" {t.range_km:5.2f} km out{' (after the plug steers)' if t.after_plug_steer else ''}"
        )
    conv = convergence_rows(lenses, lens_rows(fine=False)) if LENS_COARSE.exists() else []
    ib._write(OUT_LENS, list(lenses))
    ib._write(OUT_PUNCH, list(punch))
    ib._write(OUT_BEST, list(best))
    for path, rows in (
        (OUT_PLATE, plates),
        (OUT_OVERLAP, overlaps),
        (OUT_CATCH, catches),
        (OUT_CATCH_GEOMETRIC, geometric),
        (OUT_BACKUP, backups),
        (OUT_FATE, fates),
        (OUT_TRIGGER, triggers),
        (OUT_CONVERGENCE, conv),
    ):
        if rows:
            ib._write(path, list(rows))
    print(f"python: wrote {OUT_DIR}/shaped_plug_*.csv")


if __name__ == "__main__":
    main()
