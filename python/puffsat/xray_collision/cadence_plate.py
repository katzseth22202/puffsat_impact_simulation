"""Result 13: the plate for the 1 Hz growth push, and how thick its Kapton sheets may be.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.cadence_plate

The parent (`sec:charge_pusher_plume`, asks S36-S37) slows the 400 km overtake push from 4 Hz at
12 MN s to 1 Hz at 14.7 MN s, which holds ~1 g on the 1500 t unit, so that the 100 km/s sheets
survive the earlier pulses' gas. This module asks what that does to the flat plate, and whether a
thicker sheet on a smaller plate pays.

Two bounds set the plate (the calcs repo's `solar_dive_growth.fly_unit`):

- **Inertia.** A plate of mass `m_p` takes impulse `I` and is stopped over stroke `s`. A
  constant-force absorber then pushes the unit (mass `M`) at `I^2 / (2 m_p s M)`, so a peak cap
  `a` needs `m_p >= I^2 / (2 s M a)`. A linear spring doubles it. This depends on the impulse per
  pulse, not on the plate's size or the mean thrust: at a fixed mean, `I` grows as 1/cadence.
- **Face survival.** The planar sheet's peak is `c(h) sigma v^2` (`sweep --sheet-plate-2d`, fit
  here): set by the sheet's mass per area `sigma` and the standoff `h`, independent of its radius.
  The calcs repo's floor (30 kg/m^2 over the area, scaled as pulse energy^(2/3)) is the point
  fireball's law, and it is printed for comparison only.

The sheet's least thickness is the parent's plume dose absorbed in its own heat capacity: Kapton
from 300 K to its 400 C rating. The doses are the parent's estimates (full absorption), not ours.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.xray_collision import charge

G0 = 9.80665
V_CLOSE = 100.0e3  # [m/s]
# --- The push (parent sec:charge_pusher_plume and S33 table; calcs solar_dive_growth).
LAUNCH_UNIT = 1500.0e3  # [kg]
K_PUSH = 6.74  # the flown slug ratio
VE_PUSH = 38.4e3  # [m/s] per kg of charge, calcs `law.push_ve(6.74, 100)` as the parent quotes it
#: (label, cadence [Hz], pulse impulse [N s], push mass ratio).
PUSHES = (
    ("4 Hz, 12 MN s (priced)", 4.0, 12.0e6, 1.325),
    ("2 Hz, 14.7 MN s", 2.0, 14.7e6, 1.330),
    ("2 Hz, 7.35 MN s (1 g)", 2.0, 7.35e6, 1.330),
    ("1 Hz, 14.7 MN s (1 g)", 1.0, 14.7e6, 1.330),
)
G_CAP = 100.0  # peak g of the growth units' cargo (calcs G_CAP, author 2026-10-16)
STROKES = (2.0, 4.0, 10.0)  # [m]; 2 m is calcs ABSORBER_STROKE
# --- The calcs repo's survival floor, for comparison (solar_dive_growth._plate_for).
SURVIVAL_AREAL = 30.0  # [kg/m^2]
FLAT_PLATE_AREA = 565.487  # [m^2], the 20 m plate with its 4 m skirt
FLAT_PLATE_PULSE = 12.5  # [kg] of sheet it was rated for at 400 MPa
# --- Plume dose on the sheet [J/m^2] at one-axis spread u = 5 / 8 / 12 km/s (parent's table and
# S33.4; full absorption, no radiation). 2 Hz is quoted at u = 8 only.
DOSES = {
    "4 Hz, 12 MN s (priced)": (635.0e3, 251.0e3, 114.0e3),
    "2 Hz, 14.7 MN s": (math.nan, 75.0e3, math.nan),
    "1 Hz, 14.7 MN s (1 g)": (47.0e3, 19.0e3, 8.0e3),
    # Not the parent's: coated_sheet.py's re-implementation of its model (which reproduces the
    # parent's other rows), at 2 Hz and 1 g.
    "2 Hz, 7.35 MN s (1 g)": (96.0e3, 38.0e3, 17.0e3),
}
KAPTON_RHO = 1420.0  # [kg/m^3]
KAPTON_CP = 1090.0  # [J/kg/K], DuPont Kapton HN at room temperature
KAPTON_RISE = 673.15 - 300.0  # [K], 300 K to the 400 C rating
# --- Face allowables [Pa]: the sandwich facesheet of Result 11, and a solid steel face without
# prestress (docs/spray_plate_steel_face_limits.md: HY-100 / 4340 Rc 38-45, maraging 300).
FACES = {"sandwich facesheet": 0.4e9, "HY-100 / 4340": 1.5e9, "maraging 300": 2.5e9}
STEEL_RHO = 7850.0
STEEL_C = 5900.0  # [m/s], longitudinal wave speed: 2L/c is the face's ring period

JSONL = charge.OUT_DIR / "sheet_plate_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "sheet_plate_2d_coarse.jsonl"
OUT_PLATE = charge.OUT_DIR / "cadence_plate.csv"
OUT_SHEET = charge.OUT_DIR / "cadence_sheet.csv"
OUT_RUNS = charge.OUT_DIR / "sheet_plate_runs.csv"


def inertia_plate(impulse: float, ship: float, g_cap: float, stroke: float) -> float:
    """Least plate [kg] that holds `ship` under `g_cap` on a constant-force absorber."""
    return float(impulse**2 / (2.0 * stroke * ship * g_cap * G0))


def calcs_floor(sheet_per_pulse: float) -> float:
    """The calcs repo's survival floor [kg]: the point-fireball E^(2/3) law (comparison only)."""
    return float(SURVIVAL_AREAL * FLAT_PLATE_AREA * (sheet_per_pulse / FLAT_PLATE_PULSE) ** (2 / 3))


def least_sheet(dose: float) -> float:
    """Least sheet [kg/m^2] that holds `dose` inside Kapton's rating."""
    return dose / (KAPTON_CP * KAPTON_RISE)


def plate_area(radius: float, skirt: float) -> float:
    """Floor plus skirt wall [m^2]."""
    return math.pi * radius**2 + 2.0 * math.pi * radius * skirt


@dataclass(frozen=True)
class Run:
    """One 2-D sheet-on-plate run, reduced."""

    grid: str
    sheet_radius: float
    height: float
    plate_radius: float
    skirt: float
    impulse_per_mv: float
    capture: float  # of the flown plate's J/mv on the same grid
    peak_coeff: float  # P_peak / (sigma v^2)
    exit_direction: float
    ring_peak_to_mean: float  # impulse per area: the hottest ring over the plate's mean


def load_runs() -> list[Run]:
    out: list[Run] = []
    for path, grid in ((JSONL, "x1"), (JSONL_COARSE, "x0.5")):
        if not path.exists():
            continue
        rows = charge.load_jsonl(path)
        ref = next(r for r in rows if float(str(r["plate_radius"])) == 10.0)
        ref_j = float(str(ref["impulse_per_mv"]))
        for r in rows:
            rs = float(str(r["sheet_radius"]))
            rp = float(str(r["plate_radius"]))
            per_mv = float(str(r["impulse_per_mv"]))
            bands = r["bands"]
            assert isinstance(bands, list)
            loads = []
            for b in bands:
                assert isinstance(b, dict)
                lo, hi = float(b["r_lo"]), min(float(b["r_hi"]), rp)
                if hi > lo:
                    loads.append(float(b["impulse"]))
            mean = per_mv / (math.pi * rp**2)  # floor-area mean, per kg of impactor
            out.append(
                Run(
                    grid=grid,
                    sheet_radius=rs,
                    height=float(str(r["height"])),
                    plate_radius=rp,
                    skirt=float(str(r["skirt"])),
                    impulse_per_mv=per_mv,
                    capture=per_mv / ref_j,
                    peak_coeff=float(str(r["peak_per_mv2"])) * math.pi * rs**2,
                    exit_direction=float(str(r["exit_direction"])),
                    ring_peak_to_mean=max(loads) / mean if loads else math.nan,
                )
            )
    return out


def peak_coeff(runs: list[Run], height: float, grid: str = "x1") -> float:
    """`c(h)`: the mean over sheet radii of `P_peak / (sigma v^2)` at standoff `height`."""
    return float(np.mean([r.peak_coeff for r in runs if r.grid == grid and r.height == height]))


@dataclass(frozen=True)
class PlateRow:
    push: str
    cadence: float
    pulse: float  # [N s]
    sheet_per_pulse: float  # [kg]
    end_mass: float  # [kg], the unit at the end of the push (the inertia bound's worst case)
    calcs_floor: float  # [kg], point-fireball law, comparison only
    inertia_2m: float  # [kg], at G_CAP
    inertia_4m: float
    inertia_10m: float
    inertia_16g_10m: float  # [kg], the parent's 16 g case
    least_sheet_u5: float  # [kg/m^2]
    least_sheet_u8: float
    least_sheet_u12: float


def plate_rows() -> list[PlateRow]:
    out = []
    for label, f, pulse, mr in PUSHES:
        end = LAUNCH_UNIT / mr
        doses = DOSES[label]
        out.append(
            PlateRow(
                push=label,
                cadence=f,
                pulse=pulse,
                sheet_per_pulse=pulse / (K_PUSH * VE_PUSH),
                end_mass=end,
                calcs_floor=calcs_floor(pulse / (K_PUSH * VE_PUSH)),
                inertia_2m=inertia_plate(pulse, end, G_CAP, STROKES[0]),
                inertia_4m=inertia_plate(pulse, end, G_CAP, STROKES[1]),
                inertia_10m=inertia_plate(pulse, end, G_CAP, STROKES[2]),
                inertia_16g_10m=inertia_plate(pulse, end, 16.0, 10.0),
                least_sheet_u5=least_sheet(doses[0]),
                least_sheet_u8=least_sheet(doses[1]),
                least_sheet_u12=least_sheet(doses[2]),
            )
        )
    return out


@dataclass(frozen=True)
class SheetRow:
    """A sheet thickness on a 1 g push: its footprint, its face pressure, and the plate."""

    push: str
    areal: float  # [kg/m^2]
    kapton_um: float
    sheet_radius: float  # [m] at the pulse's sheet mass
    height: float  # [m]
    peak: float  # [Pa], production grid
    face: str  # the lightest face that holds it ("none" if no face does)
    plate_radius: float  # [m], sheet + 1.5 m, with a 2 m skirt
    plate_area: float  # [m^2]
    plate_mass: float  # [kg], max(inertia at 2 m / 100 g, 30 kg/m^2 floor)
    plate_steel_mm: float  # that mass as a solid steel plate over the area
    ring_period_us: float  # 2L/c of that plate: the face's pressure must rise slower


SHEET_AREALS = (0.05, 0.088, 0.12, 0.142, 0.2, 0.35, 0.72, 1.42)
#: The 1 g pushes the sheet table is run for.
SHEET_PUSHES = ("1 Hz, 14.7 MN s (1 g)", "2 Hz, 7.35 MN s (1 g)")


def sheet_rows(runs: list[Run]) -> list[SheetRow]:
    out = []
    for label, _, pulse, mr in (p for p in PUSHES if p[0] in SHEET_PUSHES):
        out += _sheet_rows(runs, label, pulse, LAUNCH_UNIT / mr)
    return out


def _sheet_rows(runs: list[Run], label: str, pulse: float, end: float) -> list[SheetRow]:
    mass = pulse / (K_PUSH * VE_PUSH)
    inertia = inertia_plate(pulse, end, G_CAP, STROKES[0])
    out = []
    for areal in SHEET_AREALS:
        rs = math.sqrt(mass / (math.pi * areal))
        for h in sorted({r.height for r in runs if r.grid == "x1"}):
            p = peak_coeff(runs, h) * areal * V_CLOSE**2
            face = next((n for n, a in FACES.items() if p <= a), "none")
            rp = rs + 1.5
            area = plate_area(rp, 2.0)
            m_p = max(inertia, SURVIVAL_AREAL * area)
            thick = m_p / (STEEL_RHO * area)
            out.append(
                SheetRow(
                    push=label,
                    areal=areal,
                    kapton_um=1e6 * areal / KAPTON_RHO,
                    sheet_radius=rs,
                    height=h,
                    peak=p,
                    face=face,
                    plate_radius=rp,
                    plate_area=area,
                    plate_mass=m_p,
                    plate_steel_mm=1e3 * thick,
                    ring_period_us=1e6 * 2.0 * thick / STEEL_C,
                )
            )
    return out


def _write(path: Path, rows: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


def main() -> None:
    runs = load_runs()
    print("  2-D sheet on plate, k = 6.74 (capture against the flown 20 m skirted plate):")
    for r in runs:
        if r.grid == "x1":
            print(
                f"    sheet R {r.sheet_radius:3.1f} h {r.height:3.1f}"
                f" | plate R {r.plate_radius:4.1f}"
                f" skirt {r.skirt:3.1f} | J/mv {r.impulse_per_mv:5.3f} ({r.capture:4.0%})"
                f" dir {r.exit_direction:.2f} | P/(sigma v^2) {r.peak_coeff:.3f}"
                f" | ring peak/mean {r.ring_peak_to_mean:4.2f}"
            )
    for h in sorted({r.height for r in runs}):
        print(
            f"  c(h = {h:.0f} m) = {peak_coeff(runs, h):.3f} (x1),"
            f" {peak_coeff(runs, h, 'x0.5'):.3f}"
            " (x0.5): the peak grows on refinement, so these are lower values"
        )
    plates = plate_rows()
    for p in plates:
        print(
            f"  {p.push:<24} sheet {p.sheet_per_pulse:5.1f} kg/pulse | calcs floor"
            f" {p.calcs_floor / 1e3:5.1f} t | inertia at {G_CAP:.0f} g:"
            f" 2 m {p.inertia_2m / 1e3:5.1f} t,"
            f" 4 m {p.inertia_4m / 1e3:5.1f} t, 10 m {p.inertia_10m / 1e3:5.1f} t"
            f" (16 g, 10 m: {p.inertia_16g_10m / 1e3:5.1f} t) | least sheet u = 5/8/12:"
            f" {p.least_sheet_u5:.3f}/{p.least_sheet_u8:.3f}/{p.least_sheet_u12:.3f} kg/m^2"
        )
    sheets = sheet_rows(runs)
    for s in sheets:
        print(
            f"  {s.push[:4]} sheet {s.areal:5.3f} kg/m^2 ({s.kapton_um:5.0f} um)"
            f" R {s.sheet_radius:5.2f} m"
            f" h {s.height:.0f} | peak {s.peak / 1e9:5.2f} GPa -> {s.face:<18} | plate R"
            f" {s.plate_radius:5.2f} + 2 m skirt, {s.plate_mass / 1e3:4.1f} t,"
            f" {s.plate_steel_mm:5.1f} mm steel, 2L/c {s.ring_period_us:4.1f} us"
        )
    _write(OUT_RUNS, list(runs))
    _write(OUT_PLATE, list(plates))
    _write(OUT_SHEET, list(sheets))
    print(f"python: wrote {OUT_RUNS}, {OUT_PLATE}, {OUT_SHEET}")


if __name__ == "__main__":
    main()
