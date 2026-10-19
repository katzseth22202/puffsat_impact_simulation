"""Result 15: the growth wave above 100 km/s -- Kapton, graphite and a bare hit on the bowl.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.fast_wave

The parent's S40 (`docs/fast_growth_wave_asks_for_impact_sim.md`, commit d9c7fa1) notes that
Result 13 ran the 2 Hz, 1 g push at 100 km/s, while three of its four closures cross 400 km at
132-208 km/s. Four things are computed here, each on a model the study already has.

1. **The sheets at speed.** `coated_sheet.py`'s plume (the parent's expanding clouds) with each
   closure's `V_CLOSE`, `K_PUSH` and `VE_PUSH`, conducted through bare Kapton. The plate side is
   `c(k, h) sigma w^2` with `c` re-run at each closure's k (`sweep --fast-wave-2d`). A gamma-law
   gas has no speed scale, so one run per k serves every speed.
2. **The Kapton window.** The thickness range that holds the face under 400 C to the last 10 ms
   (thermal floor) and the plate under its allowable (pressure ceiling), at 2, 4 and 8 m standoff.
   Longer standoffs carry `c(6.74, h)` from Result 13, scaled by `c(k, 2) / c(6.74, 2)`.
3. **Graphite impactors.**
   - *Heating:* the same conduction with graphite's specific heat rising with T (NIST-JANAF) and
     its face sublimation (C atoms, Hertz-Knudsen; the real vapour has C2/C3 too, so the loss is
     taken x `VAPOUR_FACTOR`).
   - *The overtake law:* Result 11's flat-plate point rows re-priced at the closing speed. The
     rebound is taken at the fireball's arrival speed and the bond toll at the heat it really
     has, and the impactor's own atomization is charged. The calcs repo's `push_ve` instead
     carries the 100 km/s law linearly in `w`.
   - *Compact disks:* the disk-into-compact-charge and point-fireball runs on the flown 9.5 m push
     plate at 2.7-6 m, against the sheet: the peak, the ring loading and the impulse.
4. **A bare compact impactor on the bowl.** The geometry first: whether the impactor's line,
   with no charge at 1.5 F, meets the δ/D 0.5 paraboloid or leaves through its mouth. Where it
   strikes, the planar-impact peak (linear `U_s-u_p` Hugoniots), whether the struck patch
   vaporizes, the hole's impulse loss (the bowl's ring loads), and the buffer spray that would
   bring the peak down to a normal pulse's.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from puffsat.xray_collision import (
    cadence_plate,
    charge,
    coated_sheet,
    departure_hole,
    freeze,
    overtake,
)
from puffsat.xray_collision.eos import EosTable

G0 = 9.80665
CADENCE = 2.0  # [Hz], the flown push (ADR-0064 amendment)
G_MEAN = 1.0
PULSE = coated_sheet.LAUNCH_UNIT * G_MEAN * G0 / CADENCE  # [N s], 7.35 MN s


@dataclass(frozen=True)
class Closure:
    """One row of the parent's S40 table: the growth stream's speed at 400 km, the push's k and
    the calcs repo's `push_ve` at that k and speed."""

    name: str
    w: float  # [m/s]
    k: float
    v_e: float  # [m/s] per kg of charge


REFERENCE = Closure("Result 13 (100 km/s)", 100.0e3, 6.74, 38.4e3)
CLOSURES = (
    Closure("9.86 R_sun, 60 (headline)", 99.2e3, 7.07, 36.7e3),
    Closure("4 R_sun, 90 (headline)", 137.1e3, 6.67, 53.0e3),
    Closure("9.86 R_sun, 90", 132.1e3, 7.93, 44.7e3),
    Closure("4 R_sun, 150", 207.8e3, 8.21, 68.4e3),
)
SPREADS = coated_sheet.SPREADS
KAPTON_EPS = 0.8
#: Face allowables [Pa] (`docs/spray_plate_steel_face_limits.md`): the two unprestressed solid
#: faces of Result 13, and the prestressed maraging face that note recommends as the ceiling.
FACES = {
    "HY-100 / 4340": 1.5e9,
    "maraging 300": 2.5e9,
    "maraging, prestressed": 3.5e9,
}
STANDOFFS = (2.0, 4.0, 8.0)  # [m]

JSONL = charge.OUT_DIR / "fast_wave_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "fast_wave_2d_coarse.jsonl"
OUT_KAPTON = charge.OUT_DIR / "fast_wave_kapton.csv"
OUT_WINDOW = charge.OUT_DIR / "fast_wave_window.csv"
OUT_GRAPHITE = charge.OUT_DIR / "fast_wave_graphite.csv"
OUT_GRAPHITE_WINDOW = charge.OUT_DIR / "fast_wave_graphite_window.csv"
OUT_LAW = charge.OUT_DIR / "fast_wave_law.csv"
OUT_ALUMINIZED = charge.OUT_DIR / "fast_wave_aluminized.csv"
OUT_COMPACT = charge.OUT_DIR / "fast_wave_compact.csv"
OUT_BOWL = charge.OUT_DIR / "fast_wave_bowl.csv"
OUT_BUFFER = charge.OUT_DIR / "fast_wave_buffer.csv"


# --- The 2-D runs ---------------------------------------------------------------------------


@cache
def _runs(grid: str = "x1") -> tuple[dict[str, object], ...]:
    path = JSONL if grid == "x1" else JSONL_COARSE
    return tuple(charge.load_jsonl(path)) if path.exists() else ()


def _f(row: dict[str, object], key: str) -> float:
    return float(str(row[key]))


def sheet_coeff(k: float, grid: str = "x1") -> float:
    """`c(k, 2 m)`: the planar sheet's `P_peak / (sigma w^2)`, the mean over the two radii run."""
    vals = [
        _f(r, "peak_per_mv2") * math.pi * _f(r, "sheet_radius") ** 2
        for r in _runs(grid)
        if r["kind"] == "sheet" and _f(r, "k") == k
    ]
    return float(np.mean(vals))


def sheet_coeff_at(k: float, height: float) -> float:
    """`c(k, h)`: the 2 m value at k, carried to `height` by Result 13's `c(6.74, h)` shape."""
    if height == 2.0:
        return sheet_coeff(k)
    runs = cadence_plate.load_runs()
    return (
        sheet_coeff(k)
        * cadence_plate.peak_coeff(runs, height)
        / cadence_plate.peak_coeff(runs, 2.0)
    )


def sheet_capture(height: float) -> float:
    """Result 13's capture at `height` against 2 m, for the 5 m sheet on its r + 1.5 m plate with
    a 2 m skirt (k = 6.74): what a longer standoff costs in impulse."""
    runs = [
        r
        for r in cadence_plate.load_runs()
        if r.grid == "x1" and r.sheet_radius == 5.0 and r.skirt == 2.0
    ]
    at = {r.height: r.impulse_per_mv for r in runs}
    return at[height] / at[2.0]


# --- 1. The sheets at speed -----------------------------------------------------------------


def push_at(c: Closure) -> coated_sheet.Push:
    return coated_sheet.push(c.name, CADENCE, G_MEAN, c.w, c.k, c.v_e)


@cache
def flux(c: Closure, u: float) -> NDArray[np.float64]:
    return coated_sheet.flux(push_at(c), u, coated_sheet._times())


def dose(c: Closure, u: float) -> float:
    q = flux(c, u)
    return float(np.sum(0.5 * (q[1:] + q[:-1]) * coated_sheet.DT))


@cache
def kapton(c: Closure, u: float, thickness: float) -> coated_sheet.Heating:
    h, _ = coated_sheet.heat(flux(c, u), thickness, KAPTON_EPS, KAPTON_EPS, 0.0)
    return h


def kapton_holds(c: Closure, u: float, thickness: float) -> bool:
    """Result 13's rule: the face stays under 400 C until the last 10 ms."""
    return kapton(c, u, thickness).rating_before <= coated_sheet.LAST


def kapton_floor(c: Closure, u: float, hi: float = 2.0e-3) -> float:
    """The thinnest Kapton [m] that holds (thicker films run cooler); NaN if even `hi` fails."""
    if not kapton_holds(c, u, hi):
        return math.nan
    lo = 10.0e-6
    if kapton_holds(c, u, lo):
        return lo
    while hi / lo > 1.02:
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if kapton_holds(c, u, mid) else (mid, hi)
    return hi


def sheet_peak(c: Closure, areal: float, height: float = 2.0) -> float:
    """[Pa] on the push plate's face."""
    return sheet_coeff_at(c.k, height) * areal * c.w**2


@dataclass(frozen=True)
class KaptonRow:
    closure: str
    w: float  # [km/s]
    k: float
    u: float  # [km/s]
    dose: float  # [kJ/m^2]
    kapton_um: float
    face_max: float  # [K]
    rating_before_ms: float
    fail_before_ms: float
    skin_um: float
    holds: bool
    peak_2m_gpa: float


KAPTON_THICKNESSES = (100e-6, 150e-6, 250e-6, 400e-6)


def kapton_rows() -> list[KaptonRow]:
    out = []
    for c in (REFERENCE, *CLOSURES):
        for u in SPREADS:
            d = dose(c, u)
            for th in KAPTON_THICKNESSES:
                h = kapton(c, u, th)
                out.append(
                    KaptonRow(
                        c.name,
                        c.w / 1e3,
                        c.k,
                        u / 1e3,
                        d / 1e3,
                        1e6 * th,
                        h.face_max,
                        1e3 * h.rating_before,
                        1e3 * h.fail_before,
                        1e6 * h.fail_depth,
                        h.rating_before <= coated_sheet.LAST,
                        sheet_peak(c, coated_sheet.KAPTON_RHO * th) / 1e9,
                    )
                )
    return out


#: A ~100 nm vapour-deposited aluminium coat (Seth's practical option) sets that face's infrared
#: emissivity to ~0.04 and adds its film's heat capacity (coated_sheet.AL_FILM).
AL_EPS = 0.04
ALUMINIZED = {
    "bare": (KAPTON_EPS, KAPTON_EPS, 0.0),
    "Al on the back": (KAPTON_EPS, AL_EPS, coated_sheet.AL_FILM),
    "Al on the plume face": (AL_EPS, KAPTON_EPS, coated_sheet.AL_FILM),
    "Al on both faces": (AL_EPS, AL_EPS, 2.0 * coated_sheet.AL_FILM),
}


@dataclass(frozen=True)
class AluminizedRow:
    closure: str
    u: float  # [km/s]
    kapton_um: float
    coat: str
    face_max: float  # [K]
    rise_over_bare: float  # [K]
    holds: bool


def aluminized_rows(c: Closure = CLOSURES[0]) -> list[AluminizedRow]:
    """Kapton with an aluminium coat on either face, in the plume of the 9.86 R_sun loop (the
    closure where Kapton is flown). The plume's heat is particles stopping in the top ~0.1 um, so
    the metal reflects none of it; it only lowers that face's infrared emission."""
    out = []
    for u in SPREADS:
        q = flux(c, u)
        for th in (100e-6, 150e-6):
            bare = coated_sheet.heat(q, th, KAPTON_EPS, KAPTON_EPS, 0.0)[0].face_max
            for coat, (ef, eb, films) in ALUMINIZED.items():
                h, _ = coated_sheet.heat(q, th, ef, eb, films)
                out.append(
                    AluminizedRow(
                        c.name,
                        u / 1e3,
                        1e6 * th,
                        coat,
                        h.face_max,
                        h.face_max - bare,
                        h.rating_before <= coated_sheet.LAST,
                    )
                )
    return out


# --- 2. The Kapton window -------------------------------------------------------------------


@dataclass(frozen=True)
class WindowRow:
    closure: str
    u: float  # [km/s]
    face: str
    height: float  # [m] standoff
    capture: float  # impulse against the 2 m standoff (Result 13, k = 6.74)
    floor_um: float  # thermal: thinnest that holds
    ceiling_um: float  # pressure: thickest the face takes
    window: bool


def window_rows() -> list[WindowRow]:
    out = []
    for c in (REFERENCE, *CLOSURES):
        for u in SPREADS:
            floor = kapton_floor(c, u)
            for face, allow in FACES.items():
                for h in STANDOFFS:
                    ceil = allow / (sheet_coeff_at(c.k, h) * c.w**2) / coated_sheet.KAPTON_RHO
                    out.append(
                        WindowRow(
                            c.name,
                            u / 1e3,
                            face,
                            h,
                            sheet_capture(h),
                            1e6 * floor,
                            1e6 * ceil,
                            bool(floor <= ceil),
                        )
                    )
    return out


# --- 3a. Graphite sheets --------------------------------------------------------------------

GRAPHITE_RHO = 2200.0  # [kg/m^3], dense graphite (materials.CARBON)
#: Through-thickness conductivity [W/m/K]: pyrolytic / flexible graphite sheet across its planes
#: is ~2-10, so 5 is taken; a film this thin is near-isothermal either way.
GRAPHITE_K = 5.0
GRAPHITE_EPS = 0.85  # the parent's
#: NIST-JANAF graphite (C, ref), `C_p` [J/mol/K] at T [K] (Chase 1998), to ~2%.
_JANAF_T = (298.15, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0)
_JANAF_CP = (8.517, 11.81, 14.62, 16.84, 19.83, 21.61, 23.67, 24.59, 25.2, 25.6)
C_MOLAR = 12.011e-3  # [kg/mol]


def graphite_cp(t: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.interp(t, _JANAF_T, _JANAF_CP) / C_MOLAR


GRAPHITE = coated_sheet.Solid(
    "graphite", GRAPHITE_RHO, GRAPHITE_K, graphite_cp, t_rating=2500.0, t_fail=3000.0
)
#: Monatomic carbon over graphite: Delta H = 716.68 kJ/mol and Delta S = 158.10 - 5.74 J/mol/K at
#: 298 K (NIST-JANAF), held constant. C2 and C3 dominate the real vapour above ~2500 K, so the
#: loss is multiplied by VAPOUR_FACTOR as a bound.
SUBLIMATION_DH = 716.68e3
SUBLIMATION_DS = 152.36
VAPOUR_FACTOR = 10.0
R_GAS = 8.314462618


def sublimation_flux(t: NDArray[np.float64]) -> NDArray[np.float64]:
    """[kg/m^2/s] off a graphite face at `t`, Hertz-Knudsen into vacuum, times VAPOUR_FACTOR."""
    p = 1.0e5 * np.exp(-SUBLIMATION_DH / (R_GAS * t) + SUBLIMATION_DS / R_GAS)
    return VAPOUR_FACTOR * p * np.sqrt(C_MOLAR / (2.0 * math.pi * R_GAS * t))


@dataclass(frozen=True)
class GraphiteRow:
    closure: str
    w: float
    u: float
    graphite_um: float
    areal: float  # [kg/m^2]
    face_max: float
    sublimated_nm: float
    peak_2m_gpa: float
    lumped_710_k: float  # the parent's check: lumped, c_p held at 710, same radiation


GRAPHITE_THICKNESSES = (10e-6, 25e-6, 50e-6, 100e-6, 250e-6, 1000e-6)


def lumped(q: NDArray[np.float64], areal: float, cp: float = 710.0) -> float:
    """The parent's lumped film: one temperature, `c_p` fixed, both faces at GRAPHITE_EPS."""
    t, top = coated_sheet.T_START, coated_sheet.T_START
    for qi in q[1:]:
        t += coated_sheet.DT * (qi - 2.0 * GRAPHITE_EPS * coated_sheet.SIGMA * t**4) / (areal * cp)
        top = max(top, t)
    return top


def graphite_rows() -> list[GraphiteRow]:
    out = []
    for c in CLOSURES:
        for u in SPREADS:
            q = flux(c, u)
            for th in GRAPHITE_THICKNESSES:
                h, face = coated_sheet.heat(q, th, GRAPHITE_EPS, GRAPHITE_EPS, 0.0, solid=GRAPHITE)
                lost = float(np.sum(sublimation_flux(face)) * coated_sheet.DT) / GRAPHITE_RHO
                areal = GRAPHITE_RHO * th
                out.append(
                    GraphiteRow(
                        c.name,
                        c.w / 1e3,
                        u / 1e3,
                        1e6 * th,
                        areal,
                        h.face_max,
                        1e9 * lost,
                        sheet_peak(c, areal) / 1e9,
                        lumped(q, areal),
                    )
                )
    return out


@dataclass(frozen=True)
class GraphiteWindowRow:
    closure: str
    face: str
    height: float  # [m]
    ceiling_um: float  # the thickest graphite film the face takes
    areal: float  # [kg/m^2]
    sheet_radius: float  # [m] at the pulse's impactor mass (graphite law)
    face_max_u5: float  # [K] at that thickness, the slowest spread
    sublimated_nm_u5: float


def graphite_window_rows() -> list[GraphiteWindowRow]:
    """For each closure and face: the pressure ceiling on a graphite film, its footprint, and
    its heating there at u = 5 km/s. Thinner films only run hotter, so this is the thermal
    worst case inside the window only if the ceiling is taken; the floor is set by sublimation,
    which stays negligible below ~2600 K."""
    out = []
    for c in CLOSURES:
        m = PULSE / (c.k * push_ve(c.k, c.w, IMPACTOR_TOLL["graphite"])[0])
        q = flux(c, SPREADS[0])
        for face, allow in FACES.items():
            for h in (2.0, 4.0):
                areal = allow / (sheet_coeff_at(c.k, h) * c.w**2)
                th = areal / GRAPHITE_RHO
                heat, hist = coated_sheet.heat(
                    q, th, GRAPHITE_EPS, GRAPHITE_EPS, 0.0, solid=GRAPHITE
                )
                lost = float(np.sum(sublimation_flux(hist)) * coated_sheet.DT) / GRAPHITE_RHO
                out.append(
                    GraphiteWindowRow(
                        c.name,
                        face,
                        h,
                        1e6 * th,
                        areal,
                        math.sqrt(m / (math.pi * areal)),
                        heat.face_max,
                        1e9 * lost,
                    )
                )
    return out


# --- 3b. The overtake law at speed ----------------------------------------------------------

#: Result 11's flat plate with a 4 m skirt, the plate the calcs repo's `push_ve` is copied from.
PUSH_PLATE = "flat + skirt 4 m, 2.7 m off"
#: The impactor's own atomization [J/kg]. Graphite: Delta_f H(C, g) = 716.68 kJ/mol (NIST-JANAF).
#: Kapton (C22H10N2O5, 382.3 g/mol): mean bond enthalpies over the repeat unit -- 18 aromatic C-C
#: at 509.5 kJ/mol (benzene's atomization less six C-H), 10 C-H at 413, 4 C=O at 745, 4 C-C at
#: 348, 6 C-N at 305, 2 C-O at 358 -- 20.2 MJ/mol. Results 10-14 never charged it (the toll
#: was the charge's polyethylene alone), so the law is shown both ways.
IMPACTOR_TOLL = {
    "none (calcs law's basis)": 0.0,
    "Kapton": (18 * 509.5 + 10 * 413.0 + 4 * 745.0 + 4 * 348.0 + 6 * 305.0 + 2 * 358.0)
    * 1e3
    / 0.38232,
    "graphite": 716.68e3 / C_MOLAR,
}


@cache
def _table() -> EosTable:
    return charge.ch2_table()


def _point(k: float, key: str) -> float:
    rows = [
        r
        for r in charge.load_jsonl(overtake.JSONL)
        if r["case"] == PUSH_PLATE and r["kind"] == "point"
    ]
    ks = np.log([_f(r, "k") for r in rows])
    return float(np.interp(math.log(k), ks, [_f(r, key) for r in rows]))


def calcs_push_ve(k: float, w: float) -> float:
    """The calcs repo's `push_ve` [m/s]: the 100 km/s carbon-atomic rows, linear in log k (its
    `_log_interp`), scaled by `w / 100 km/s`."""
    rows = [
        r
        for r in csv.DictReader(overtake.OUT.open())
        if r["case"] == PUSH_PLATE and r["kind"] == "point" and r["grid"] == "x1"
    ]
    ks = np.log([float(r["k"]) for r in rows])
    ve = [float(r["v_e_carbon_frozen"]) for r in rows]
    return float(np.interp(math.log(k), ks, ve)) * w / 100.0e3


@cache
def push_ve(k: float, w: float, toll: float) -> tuple[float, float, float]:
    """`(v_e [m/s], e_eff, velocity loss)`: Result 11's flat-plate point fireball priced at `w`,
    carbon atomic (H2 re-formed), plus the impactor's own atomization `toll` [J/kg]."""
    per_mv = _point(k, "impulse_per_mv")
    peak = _point(k, "peak_per_mv2") * w**2
    m_max = overtake.FACESHEET / peak
    column = (1.0 + k) * min(m_max, 20.0) / (2.0 * math.pi * charge.PLATE_R**2)
    fire = charge.charge_fireball(k, table=_table(), v_rel=w)
    e = charge.rebound(fire.u_rms + fire.drift, column)
    confined = _f(charge.load_jsonl(overtake.JSONL)[0], "confined_ratio")
    j = per_mv * w * (1.0 + e) / confined
    fr = freeze.freeze(k, _table(), v_rel=w)
    carbon_part = (freeze.CH2_ATOMIZATION - freeze.H2_BOND) / freeze.CH2_ATOMIZATION
    heat = 0.5 * k / (1.0 + k) * w**2  # [J per kg of impactor]
    share = fr.bond_share * carbon_part + toll / heat
    loss = 1.0 - math.sqrt(max(1.0 - share, 0.0))
    return j / k * (1.0 - loss), e, loss


@dataclass(frozen=True)
class LawRow:
    closure: str
    w: float  # [km/s]
    k: float
    calcs_ve: float  # [km/s], `push_ve` as the parent quotes it
    impactor: str
    ve: float  # [km/s], re-priced at w
    ratio: float  # re-priced over calcs
    speed_factor: float  # re-priced at w over (re-priced at 100 km/s) x w / 100 km/s
    e_eff: float
    loss: float
    impactor_per_pulse: float  # [kg] at 7.35 MN s, re-priced


def law_rows() -> list[LawRow]:
    out = []
    for c in (REFERENCE, *CLOSURES):
        calc = calcs_push_ve(c.k, c.w)
        for name, toll in IMPACTOR_TOLL.items():
            ve, e, loss = push_ve(c.k, c.w, toll)
            linear = push_ve(c.k, 100.0e3, toll)[0] * c.w / 100.0e3
            out.append(
                LawRow(
                    c.name,
                    c.w / 1e3,
                    c.k,
                    calc / 1e3,
                    name,
                    ve / 1e3,
                    ve / calc,
                    ve / linear,
                    e,
                    loss,
                    PULSE / (c.k * ve),
                )
            )
    return out


# --- 3c. Compact graphite on the push plate -------------------------------------------------

#: The calcs repo's survival floor (cadence_plate.calcs_floor), on pulse energy: 30 kg/m^2 over
#: the 20 m skirted plate, rated for a 12.5 kg sheet at 100 km/s, as energy^(2/3).
FLOOR_RATED = 0.5 * cadence_plate.FLAT_PLATE_PULSE * 100.0e3**2
DISK_RADIUS = 0.25  # [m], the runs' compact disk (MIX_RADIUS)


@dataclass(frozen=True)
class CompactRow:
    closure: str
    kind: str  # "sheet (5 m)", "compact" or "point"
    height: float
    impactor: float  # [kg] per pulse
    capture: float  # J/mv against the 5 m sheet at 2 m on its own plate
    peak_gpa: float
    ring_peak_to_mean: float  # hottest ring's impulse per area over the plate's mean
    face: str  # lightest face that takes the peak
    calcs_floor_t: float  # [t], the point-blast law at this pulse's energy
    inertia_t: float  # [t], the absorber's bound (Result 13: independent of the impactor)


def _ring_peak_to_mean(row: dict[str, object]) -> float:
    rp = _f(row, "plate_radius")
    bands = row["bands"]
    assert isinstance(bands, list)
    loads = [float(b["impulse"]) for b in bands if isinstance(b, dict) and float(b["r_lo"]) < rp]
    return max(loads) / (_f(row, "impulse_per_mv") / (math.pi * rp**2))


def compact_rows() -> list[CompactRow]:
    out = []
    inertia = cadence_plate.inertia_plate(
        PULSE, coated_sheet.LAUNCH_UNIT / 1.330, cadence_plate.G_CAP, 2.0
    )
    for c in CLOSURES:
        ve, _, _ = push_ve(c.k, c.w, IMPACTOR_TOLL["graphite"])
        m = PULSE / (c.k * ve)
        floor = cadence_plate.calcs_floor(cadence_plate.FLAT_PLATE_PULSE) * (
            0.5 * m * c.w**2 / FLOOR_RATED
        ) ** (2.0 / 3.0)
        rows = [r for r in _runs() if _f(r, "k") == c.k]
        sheet = next(r for r in rows if r["kind"] == "sheet" and _f(r, "sheet_radius") == 5.0)
        ref = _f(sheet, "impulse_per_mv")
        for r in rows:
            if r["kind"] == "sheet" and _f(r, "sheet_radius") != 5.0:
                continue
            if r["kind"] == "sheet":
                peak = sheet_peak(c, coated_sheet.KAPTON_RHO * 100e-6)
                label = "sheet (5 m, 100 um Kapton)"
            else:
                peak = _f(r, "peak_per_mv2") * m * c.w**2
                label = str(r["kind"])
            out.append(
                CompactRow(
                    c.name,
                    label,
                    _f(r, "height"),
                    m,
                    _f(r, "impulse_per_mv") / ref,
                    peak / 1e9,
                    _ring_peak_to_mean(r),
                    next((n for n, a in FACES.items() if peak <= a), "none"),
                    floor / 1e3,
                    inertia / 1e3,
                )
            )
    return out


# --- 4. A bare compact impactor on the bowl -------------------------------------------------

BOWL = departure_hole.DESIGN_PLATE  # "paraboloid 0.50 at 1.5 F"
BOWL_D_OVER_D = 0.5
BOWL_R = charge.PLATE_R
BOWL_DEPTH = BOWL_D_OVER_D * 2.0 * BOWL_R  # [m], the rim's height above the vertex
BOWL_FOCUS = BOWL_R / (8.0 * BOWL_D_OVER_D)
BOWL_K = 20.0  # the departure's k (Result 14: k ~ 20-30)
BOWL_HOLE = 0.5  # [m], the hole's radius (Result 14's largest)
#: A compact graphite disk that threads the hole with 10 cm to spare.
BOWL_DISK_RADIUS = 0.4
STEEL_RHO = 7850.0
#: Linear Hugoniots `U_s = c0 + s u_p` [m/s]. Iron: Marsh (LASL 1980) as tabulated by Meyers,
#: *Dynamic Behavior of Materials*, Table 4.1. Graphite: a bracket, not a fit -- the lower bound
#: `c0 = 0, s = 1` (`P = rho u_p^2`) and a stiff `c0 = 4 km/s, s = 1.7`.
IRON_HUGONIOT = (3574.0, 1.92)
GRAPHITE_HUGONIOTS = {"lower bound": (0.0, 1.0), "stiff": (4000.0, 1.7)}
STEEL_VAPORIZE = 7.0e6  # [J/kg] to vaporize iron from 300 K (~6.3 heat + ~0.7 warming), rounded up
INJECTION_SPEEDS = (20.0, 50.0, 100.0)  # [m/s], the buffer spray's (not set by any study yet)


def charge_height(times_focus: float = 1.5) -> float:
    return times_focus * BOWL_FOCUS


def exit_point(theta: float, h: float) -> tuple[float, float]:
    """`(r, z)` where the impactor's line, through the charge's place at height `h` with no charge
    there, meets the paraboloid `z = r^2 / 4F` again on the far side (z may exceed the rim)."""
    phi = math.radians(180.0 - theta)
    # Beyond the charge: z = h + t cos phi, r = t sin phi.
    a, b, cc = math.sin(phi) ** 2 / (4.0 * BOWL_FOCUS), -math.cos(phi), -h
    t = (-b + math.sqrt(b * b - 4.0 * a * cc)) / (2.0 * a)
    return t * math.sin(phi), h + t * math.cos(phi)


def critical_theta(h: float) -> float:
    """[deg] above which the bare impactor leaves through the mouth: the line through the rim."""
    return 180.0 - math.degrees(math.atan2(BOWL_R, BOWL_DEPTH - h))


def graze(theta: float, r: float) -> float:
    """[rad] between the impactor's path and the wall at radius `r` (90 deg is normal)."""
    phi = math.radians(180.0 - theta)
    slope = r / (2.0 * BOWL_FOCUS)  # dz/dr of the wall
    cross = (math.cos(phi) - slope * math.sin(phi)) / math.hypot(1.0, slope)
    return math.asin(abs(cross))


def impact_pressure(v: float, projectile: tuple[float, float]) -> float:
    """[Pa] Planar impedance match of graphite at `v` onto iron (both linear Hugoniots)."""
    c0p, sp = projectile
    c0t, st = IRON_HUGONIOT

    def mismatch(u_t: float) -> float:
        u_p = v - u_t
        return STEEL_RHO * (c0t + st * u_t) * u_t - GRAPHITE_RHO * (c0p + sp * u_p) * u_p

    u_t = float(brentq(mismatch, 0.0, v))
    return STEEL_RHO * (c0t + st * u_t) * u_t


def ring_loss(theta: float, k: float, r: float, area: float) -> float:
    """Impulse share landing on `area` [m^2] of the bowl at radius `r` (the 2-D ring loads),
    as `departure_hole.hole_loss` prices the entry hole."""
    row = departure_hole._runs("x1")[(BOWL, theta, k)]
    bands = row["bands"]
    assert isinstance(bands, list)
    band = next(
        b for b in bands if isinstance(b, dict) and float(b["r_lo"]) <= r < float(b["r_hi"])
    )
    return float(band["impulse"]) * area / _f(row, "impulse_per_mv")


@dataclass(frozen=True)
class BowlRow:
    route: str
    theta: float
    w: float  # [km/s]
    charge_at_focus: float  # the charge's height in focal lengths
    critical_theta: float
    strikes: bool
    exit_r: float
    exit_z: float
    rim_clearance: float  # [m] between the line and the rim at the rim's height (if it misses)
    disk_clears: bool  # the clearance exceeds the disk's radius
    graze_deg: float
    v_normal: float  # [km/s]
    impactor: float  # [kg]
    disk_mm: float  # its thickness at BOWL_DISK_RADIUS
    wall_mm: float  # the bowl's steel at its mass (Result 14)
    peak_low_gpa: float
    peak_stiff_gpa: float
    patch_mj_per_kg: float  # dissipated energy per kg of disk plus struck patch
    patch_vaporizes: bool
    hole_m2: float  # the footprint, the least hole
    impulse_loss_footprint: float
    impulse_loss_twice: float  # a hole of twice the footprint's diameter
    debris_speed: float  # [km/s], the merged patch's speed (tangential plus shared normal)


def _stage(route: str, k: float) -> dict[str, str]:
    rows = list(csv.DictReader(departure_hole.OUT_STAGE.open()))
    return next(
        r for r in rows if r["route"] == route and float(r["k"]) == k and r["carbon"] == "atomic"
    )


def bowl_rows() -> list[BowlRow]:
    out = []
    area = charge.plate_area(BOWL_D_OVER_D, 0.0)
    for route in (r for r in departure_hole.ROUTES if r.sign == "visible"):
        stage = _stage(f"{route.name}, {route.sign}", BOWL_K)
        m = float(stage["impactor_per_pulse"])
        wall = float(stage["plate"]) / (area * STEEL_RHO)
        disk = m / (GRAPHITE_RHO * math.pi * BOWL_DISK_RADIUS**2)
        for f in (1.5, 2.0, 2.5):
            h = charge_height(f)
            r, z = exit_point(route.theta, h)
            hits = z <= BOWL_DEPTH
            phi = math.radians(180.0 - route.theta)
            clear = BOWL_R - (BOWL_DEPTH - h) * math.tan(phi)
            row = dict(
                route=f"{route.name}, {route.sign}",
                theta=route.theta,
                w=route.w / 1e3,
                charge_at_focus=f,
                critical_theta=critical_theta(h),
                strikes=hits,
                exit_r=r,
                exit_z=z,
                rim_clearance=0.0 if hits else clear,
                disk_clears=(not hits) and clear > BOWL_DISK_RADIUS,
                impactor=m,
                disk_mm=1e3 * disk,
                wall_mm=1e3 * wall,
            )
            if hits:
                a = graze(route.theta, r)
                vn = route.w * math.sin(a)
                footprint = math.pi * BOWL_DISK_RADIUS**2 / math.sin(a)
                # Disk and struck patch merge inelastically along the normal.
                s_disk = GRAPHITE_RHO * disk * math.sin(a)  # [kg/m^2] of wall footprint
                s_wall = STEEL_RHO * wall
                dissipated = 0.5 * s_disk * s_wall / (s_disk + s_wall) * vn**2
                per_kg = dissipated / (s_disk + s_wall)
                v_merged = math.hypot(route.w * math.cos(a), vn * s_disk / (s_disk + s_wall))
                row |= dict(
                    graze_deg=math.degrees(a),
                    v_normal=vn / 1e3,
                    peak_low_gpa=impact_pressure(vn, GRAPHITE_HUGONIOTS["lower bound"]) / 1e9,
                    peak_stiff_gpa=impact_pressure(vn, GRAPHITE_HUGONIOTS["stiff"]) / 1e9,
                    patch_mj_per_kg=per_kg / 1e6,
                    patch_vaporizes=per_kg > STEEL_VAPORIZE,
                    hole_m2=footprint,
                    impulse_loss_footprint=ring_loss(route.theta, BOWL_K, r, footprint),
                    impulse_loss_twice=ring_loss(route.theta, BOWL_K, r, 4.0 * footprint),
                    debris_speed=v_merged / 1e3,
                )
            else:
                row |= dict(
                    graze_deg=math.nan,
                    v_normal=math.nan,
                    peak_low_gpa=math.nan,
                    peak_stiff_gpa=math.nan,
                    patch_mj_per_kg=math.nan,
                    patch_vaporizes=False,
                    hole_m2=0.0,
                    impulse_loss_footprint=0.0,
                    impulse_loss_twice=0.0,
                    debris_speed=route.w / 1e3,
                )
            out.append(BowlRow(**row))  # type: ignore[arg-type]
    return out


@dataclass(frozen=True)
class BufferRow:
    route: str
    k_buffer: float
    normal_peak_gpa: float  # a normal pulse's peak on the bowl (charge at 1.5 F, k = BOWL_K)
    height_needed: float  # [m] off the struck wall for the buffer's peak to match it
    methane: float  # [kg] of spray
    lead_s: tuple[float, ...]  # [s] at INJECTION_SPEEDS


def _point_peak(k: float, h: float) -> float:
    """Point fireball `peak_per_mv2` on the flat plate at standoff `h`: a power law in h through
    the 2.7-6 m runs (extrapolated outside them)."""
    rows = [r for r in _runs() if r["kind"] == "point" and _f(r, "k") == k]
    hs = np.log([_f(r, "height") for r in rows])
    ps = np.log([_f(r, "peak_per_mv2") for r in rows])
    slope, icept = np.polyfit(hs, ps, 1)
    return float(math.exp(icept + slope * math.log(h)))


def buffer_rows() -> list[BufferRow]:
    out = []
    for route in (r for r in departure_hole.ROUTES if r.sign == "visible"):
        name = f"{route.name}, {route.sign}"
        m = float(_stage(name, BOWL_K)["impactor_per_pulse"])
        normal = (
            _f(departure_hole._runs("x1")[(BOWL, route.theta, BOWL_K)], "peak_per_mv2")
            * m
            * route.w**2
        )
        for kb in (2.0, 5.0):
            # The full closing speed drives the buffer fireball onto the wall: an upper bound,
            # since only the normal share closes on it.
            def excess(
                x: float, kb: float = kb, load: float = m * route.w**2, cap: float = normal
            ) -> float:
                return _point_peak(kb, x) * load - cap

            h = float(brentq(excess, 0.5, 200.0))
            out.append(
                BufferRow(
                    name,
                    kb,
                    normal / 1e9,
                    h,
                    kb * m,
                    tuple(h / v for v in INJECTION_SPEEDS),
                )
            )
    return out


# --- Output ---------------------------------------------------------------------------------


def _write(path: Path, rows: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


def main() -> None:
    print("  c(k, 2 m) = P / (sigma w^2), planar sheet on the flown 9.5 m plate (x1 / x0.5):")
    for c in (REFERENCE, *CLOSURES):
        coarse = sheet_coeff(c.k, "x0.5") if _runs("x0.5") else math.nan
        print(f"    k {c.k:4.2f}: {sheet_coeff(c.k):.3f} / {coarse:.3f}")
    kap = kapton_rows()
    for r in kap:
        print(
            f"  {r.closure:<26} u {r.u:2.0f} dose {r.dose:4.0f} kJ/m^2 | {r.kapton_um:3.0f} um:"
            f" face {r.face_max:5.0f} K, skin {r.skin_um:3.0f} um,"
            f" {'holds' if r.holds else 'fails'} | {r.peak_2m_gpa:4.2f} GPa at 2 m"
        )
    alu = aluminized_rows()
    for ar in alu:
        print(
            f"  aluminized {ar.closure:<26} u {ar.u:2.0f} {ar.kapton_um:3.0f} um {ar.coat:<22}"
            f" face {ar.face_max:5.0f} K (+{ar.rise_over_bare:3.0f} K)"
            f" {'holds' if ar.holds else 'fails'}"
        )
    win = window_rows()
    for wr in win:
        if wr.u == 8.0 or wr.window:
            print(
                f"  window {wr.closure:<26} u {wr.u:2.0f} {wr.face:<22} h {wr.height:.0f} m"
                f" (capture {wr.capture:.2f}) | Kapton {wr.floor_um:5.0f}-{wr.ceiling_um:5.0f} um"
                f" {'OPEN' if wr.window else 'shut'}"
            )
    gra = graphite_rows()
    for g in gra:
        print(
            f"  graphite {g.closure:<26} u {g.u:2.0f} {g.graphite_um:5.0f} um"
            f" ({g.areal:.3f} kg/m^2): face {g.face_max:5.0f} K (lumped, c_p 710:"
            f" {g.lumped_710_k:5.0f} K), sublimed {g.sublimated_nm:8.2e} nm |"
            f" {g.peak_2m_gpa:5.2f} GPa at 2 m"
        )
    gwin = graphite_window_rows()
    for gw in gwin:
        print(
            f"  graphite window {gw.closure:<26} {gw.face:<22} h {gw.height:.0f} m: <="
            f" {gw.ceiling_um:4.0f} um ({gw.areal:.3f} kg/m^2, R {gw.sheet_radius:4.1f} m) |"
            f" u = 5: face {gw.face_max_u5:5.0f} K, sublimed {gw.sublimated_nm_u5:.1e} nm"
        )
    law = law_rows()
    for lr in law:
        print(
            f"  law {lr.closure:<26} w {lr.w:5.1f} k {lr.k:4.2f} | calcs {lr.calcs_ve:5.1f} km/s"
            f" | {lr.impactor:<24} {lr.ve:5.1f} km/s (x{lr.ratio:.3f};"
            f" speed x{lr.speed_factor:.3f}), e {lr.e_eff:.3f},"
            f" toll {lr.loss:.1%} | {lr.impactor_per_pulse:4.1f} kg/pulse"
        )
    com = compact_rows()
    for cr in com:
        print(
            f"  push plate {cr.closure:<26} {cr.kind:<27} h {cr.height:3.1f}"
            f" | {cr.impactor:4.1f} kg | capture {cr.capture:.2f}"
            f" | peak {cr.peak_gpa:6.2f} GPa -> {cr.face:<22}"
            f" | ring peak/mean {cr.ring_peak_to_mean:5.2f} | floor {cr.calcs_floor_t:4.1f} t,"
            f" inertia {cr.inertia_t:4.1f} t"
        )
    bowl = bowl_rows()
    for b in bowl:
        head = (
            f"  bowl {b.route:<20} {b.theta:5.1f} deg, charge at {b.charge_at_focus:.1f} F"
            f" (critical {b.critical_theta:5.1f} deg):"
        )
        if not b.strikes:
            print(
                f"{head} leaves through the mouth, {b.rim_clearance:.1f} m inside the rim"
                f" ({'clears' if b.disk_clears else 'clips the lip with'} the"
                f" {BOWL_DISK_RADIUS:.1f} m disk)"
            )
            continue
        print(
            f"{head} strikes r {b.exit_r:.2f} z {b.exit_z:.2f} at {b.graze_deg:.0f} deg,"
            f" {b.v_normal:.0f} km/s normal | {b.impactor:.1f} kg, {b.disk_mm:.1f} mm disk on"
            f" {b.wall_mm:.1f} mm wall | {b.peak_low_gpa:.0f}-{b.peak_stiff_gpa:.0f} GPa |"
            f" {b.patch_mj_per_kg:.0f} MJ/kg ({'vaporizes' if b.patch_vaporizes else 'holds'})"
            f" | hole >= {b.hole_m2:.2f} m^2: {b.impulse_loss_footprint:.2%}-"
            f"{b.impulse_loss_twice:.2%} of the impulse | debris {b.debris_speed:.0f} km/s"
        )
    buf = buffer_rows()
    for br in buf:
        print(
            f"  buffer {br.route:<20} k {br.k_buffer:.0f}: normal pulse {br.normal_peak_gpa:.3f}"
            f" GPa -> {br.height_needed:4.1f} m off the wall, {br.methane:4.0f} kg CH4, lead "
            + "/".join(f"{t:.2f}" for t in br.lead_s)
            + f" s at {'/'.join(f'{v:.0f}' for v in INJECTION_SPEEDS)} m/s"
        )
    _write(OUT_KAPTON, list(kap))
    _write(OUT_WINDOW, list(win))
    _write(OUT_GRAPHITE, list(gra))
    _write(OUT_GRAPHITE_WINDOW, list(gwin))
    _write(OUT_LAW, list(law))
    _write(OUT_ALUMINIZED, list(alu))
    _write(OUT_COMPACT, list(com))
    _write(OUT_BOWL, list(bowl))
    _write(OUT_BUFFER, list(buf))
    print(f"python: wrote {OUT_KAPTON.parent}/fast_wave_*.csv")


if __name__ == "__main__":
    main()
