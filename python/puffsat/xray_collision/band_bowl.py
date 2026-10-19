"""Result 16: Kapton at 115-122 km/s, the bowl at 2.5 F, and the bowl and push on the Jupiter chain.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.band_bowl

The parent's follow-ups S44 (`docs/fast_growth_wave_asks_for_impact_sim.md`, commit 5d1379c) and
S42 (`docs/plate_jupiter_loop_asks.md`) ask four things, each on a model the study already has.

1. **Kapton in the 100-130 km/s band (S44.1).** Three printed closures sit between Result 13's
   99 km/s and Result 15's 132 km/s. Result 15's window is re-run at each closure's own k
   (`sweep --band-bowl-2d` re-runs the planar sheet's `c(k, 2 m)` at k = 6.20, 6.53 and 7.45):
   the thinnest Kapton that holds its face under 400 C to the last 10 ms (`fast_wave.kapton_floor`)
   against the thickest the face takes at the flown 2 m standoff.
2. **The bowl with its charge at 2.5 F (S44.2).** Result 15 moved the charge out so a bare
   impactor clears the rim at 117 deg. The point fireball is re-run on that bowl at 117, 129.9 and
   180 deg, and Result 14's `hole_ve`, hole loss and stage are re-priced on it.
3. **The bowl nearly head-on (S42.1).** The flown 1.5 F bowl at 160-176 deg, between Result 14's
   129.9 and 180 deg runs, at the calcs repo's dollar-optimal k ~ 12.5. The calcs repo
   interpolated in `cos theta` between 129.9 and 180; the runs check it.
4. **The push at 54-69 km/s (S42.2).** Result 15's re-priced overtake law (`fast_wave.push_ve`)
   at the Jupiter chain's slower growth streams, against the calcs repo's linear-in-w law, and
   the sheet's heating there at 2 Hz / 1 g.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np

from puffsat.xray_collision import charge, coated_sheet, departure_hole, fast_wave

G0 = 9.80665
SHEET_JSONL = charge.OUT_DIR / "band_sheet_2d.jsonl"
SHEET_JSONL_COARSE = charge.OUT_DIR / "band_sheet_2d_coarse.jsonl"

#: S44 item 1: the ledger's own k and `push_ve` on Kapton.
BAND = (
    fast_wave.Closure("9.86 R_sun, X 75", 114.8e3, 7.45, 40.9e3),
    fast_wave.Closure("4 R_sun, X 90, 5-day gap", 119.0e3, 6.20, 48.4e3),
    fast_wave.Closure("4 R_sun, X 75", 121.8e3, 6.53, 47.8e3),
)
#: The faces S44 names: the unprestressed solid faces of Result 13.
FACES = {k: v for k, v in fast_wave.FACES.items() if k in ("HY-100 / 4340", "maraging 300")}
STANDOFF = 2.0  # [m], the flown push

FAR = "paraboloid 0.50 at 2.5 F"
FLOWN = departure_hole.DESIGN_PLATE  # "paraboloid 0.50 at 1.5 F"
FAR_THETAS = (117.0, 129.9)
HEAD_THETAS = (160.0, 164.0, 168.0, 176.0)
HEAD_KS = (7.0, 10.0, 12.5, 14.0)
#: The Jupiter-only chain's departure streams (calcs ADR 0075: 160-180 deg at 66-76 km/s; the
#: parent's S42 ballpark, 70-85 km/s).
HEAD_SPEEDS = (66.0e3, 76.0e3, 85.0e3)
#: S42.2: the chain's growth stream at 400 km (calcs: 54-65 km/s; parent: 57-69) and its k.
PUSH_SPEEDS = (54.0e3, 57.0e3, 61.0e3, 65.0e3, 69.0e3)
PUSH_KS = (6.0, 6.5, 6.8, 8.0)

OUT_WINDOW = charge.OUT_DIR / "band_window.csv"
OUT_FAR = charge.OUT_DIR / "band_far_bowl.csv"
OUT_HEAD = charge.OUT_DIR / "band_head_on.csv"
OUT_PUSH = charge.OUT_DIR / "band_slow_push.csv"


# --- 1. Kapton in the band ------------------------------------------------------------------


@cache
def _sheets(grid: str = "x1") -> tuple[dict[str, object], ...]:
    path = SHEET_JSONL if grid == "x1" else SHEET_JSONL_COARSE
    return tuple(charge.load_jsonl(path))


def sheet_coeff(k: float, grid: str = "x1") -> float:
    """`c(k, 2 m)`: the planar sheet's `P_peak / (sigma w^2)`, the mean over the two radii run."""
    vals = [
        float(str(r["peak_per_mv2"])) * math.pi * float(str(r["sheet_radius"])) ** 2
        for r in _sheets(grid)
        if float(str(r["k"])) == k
    ]
    return float(np.mean(vals))


@dataclass(frozen=True)
class WindowRow:
    closure: str
    w: float  # [km/s]
    k: float
    c_2m: float
    c_2m_coarse: float
    u: float  # [km/s]
    dose: float  # [kJ/m^2]
    floor_um: float  # thinnest Kapton that holds
    face_100um_k: float  # the 100 um sheet's face peak [K]
    face: str
    ceiling_um: float  # thickest the face takes at 2 m
    window: bool


def window_rows() -> list[WindowRow]:
    out = []
    for c in BAND:
        coeff = sheet_coeff(c.k)
        for u in coated_sheet.SPREADS:
            floor = fast_wave.kapton_floor(c, u)
            for face, allow in FACES.items():
                ceil = allow / (coeff * c.w**2) / coated_sheet.KAPTON_RHO
                out.append(
                    WindowRow(
                        c.name,
                        c.w / 1e3,
                        c.k,
                        coeff,
                        sheet_coeff(c.k, "x0.5"),
                        u / 1e3,
                        fast_wave.dose(c, u) / 1e3,
                        1e6 * floor,
                        fast_wave.kapton(c, u, 100e-6).face_max,
                        face,
                        1e6 * ceil,
                        bool(floor <= ceil),
                    )
                )
    return out


# --- 2. The bowl at 2.5 F -------------------------------------------------------------------


def _hole_ring(plate: str, theta: float, k: float) -> tuple[float, float, float]:
    """`(r_hole, peak on the hole's band [Pa per (m w^2)], energy in on it [per m w^2])`."""
    row = departure_hole._runs("x1")[(plate, theta, k)]
    r_hole, _ = departure_hole.hole_geometry(plate, theta)
    bands = row["bands"]
    assert isinstance(bands, list)
    band = next(
        b for b in bands if isinstance(b, dict) and float(b["r_lo"]) <= r_hole < float(b["r_hi"])
    )
    return r_hole, float(band["peak_pressure"]), float(band["energy_in"])


@dataclass(frozen=True)
class FarRow:
    route: str
    theta: float
    k: float
    capture_15: float  # J / (m w) on the flown 1.5 F bowl
    capture_25: float
    hole_ve_15: float  # [km/s], carbon atomic
    hole_ve_25: float
    isp_25: float  # [s]
    hole_r_15: float  # [m], where the line through the charge meets the bowl
    hole_r_25: float
    hole_share_15: float  # impulse on a 0.5 m hole
    hole_share_25: float
    hole_gas_15: float  # gas share on that hole
    hole_gas_25: float
    rim_peak_ratio: float  # peak on the hole's band, 2.5 F over 1.5 F
    rim_energy_ratio: float  # energy arriving per area there, 2.5 F over 1.5 F
    plate_peak_ratio: float  # the bowl's peak anywhere, 2.5 F over 1.5 F
    mass_ratio_15: float
    mass_ratio_25: float


def far_rows() -> list[FarRow]:
    out = []
    routes = {r.theta: r for r in departure_hole.ROUTES if r.sign == "visible"}
    runs = departure_hole._runs("x1")
    for theta in FAR_THETAS:
        route = routes[theta]
        for k in (10.0, 20.0, 30.0):
            v15 = departure_hole.hole_ve(FLOWN, theta, k, route.w)
            v25 = departure_hole.hole_ve(FAR, theta, k, route.w)
            s15, g15 = departure_hole.hole_loss(FLOWN, theta, k, 0.5)
            s25, g25 = departure_hole.hole_loss(FAR, theta, k, 0.5)
            r15, p15, e15 = _hole_ring(FLOWN, theta, k)
            r25, p25, e25 = _hole_ring(FAR, theta, k)
            peak = float(str(runs[(FAR, theta, k)]["peak_per_mv2"])) / float(
                str(runs[(FLOWN, theta, k)]["peak_per_mv2"])
            )
            out.append(
                FarRow(
                    f"{route.name}, {route.sign}",
                    theta,
                    k,
                    departure_hole.capture(FLOWN, theta, k),
                    departure_hole.capture(FAR, theta, k),
                    v15 / 1e3,
                    v25 / 1e3,
                    v25 / G0,
                    r15,
                    r25,
                    s15,
                    s25,
                    g15,
                    g25,
                    p25 / p15,
                    e25 / e15,
                    peak,
                    departure_hole.stage(route, k).mass_ratio,
                    departure_hole.stage(route, k, bowl=FAR).mass_ratio,
                )
            )
    return out


# --- 3. The bowl nearly head-on -------------------------------------------------------------


def capture_interpolated(theta: float, k: float) -> float:
    """What the calcs repo used: linear in `cos theta` between the 129.9 and 180 deg runs only."""
    lo = departure_hole._capture_at(FLOWN, 129.9, k, "x1")
    hi = departure_hole._capture_at(FLOWN, 180.0, k, "x1")
    c, c_lo, c_hi = (math.cos(math.radians(t)) for t in (theta, 129.9, 180.0))
    return lo + (hi - lo) * (c - c_lo) / (c_hi - c_lo)


@dataclass(frozen=True)
class HeadRow:
    theta: float
    k: float
    capture_run: float
    capture_coarse: float
    capture_interpolated: float
    error: float  # interpolated over run, less 1
    bowl_ceiling: float  # cos theta + sqrt(1 + k)
    share_of_ceiling: float
    hole_r: float  # [m]
    hole_share: float  # impulse on a 0.5 m hole
    hole_gas: float
    ve_66: float  # [km/s], hole_ve at HEAD_SPEEDS, carbon atomic
    ve_76: float
    ve_85: float


def head_rows() -> list[HeadRow]:
    out = []
    for theta in HEAD_THETAS:
        for k in HEAD_KS:
            run = departure_hole._capture_at(FLOWN, theta, k, "x1")
            interp = capture_interpolated(theta, k)
            ceiling = departure_hole.bowl_ceiling(k, theta)
            share, gas = departure_hole.hole_loss(FLOWN, theta, k, 0.5)
            r_hole, _ = departure_hole.hole_geometry(FLOWN, theta)
            ves = [departure_hole.hole_ve(FLOWN, theta, k, w) / 1e3 for w in HEAD_SPEEDS]
            out.append(
                HeadRow(
                    theta,
                    k,
                    run,
                    departure_hole._capture_at(FLOWN, theta, k, "x0.5"),
                    interp,
                    interp / run - 1.0,
                    ceiling,
                    run / ceiling,
                    r_hole,
                    share,
                    gas,
                    *ves,
                )
            )
    return out


# --- 4. The push at 54-69 km/s --------------------------------------------------------------


@dataclass(frozen=True)
class PushRow:
    w: float  # [km/s]
    k: float
    calcs_ve: float  # [km/s], the 100 km/s law x w / 100
    ve: float  # [km/s], re-priced at w (Kapton's atomization charged)
    ratio: float
    e_eff: float
    toll: float  # velocity share lost to bonds
    impactor_kg: float  # per 7.35 MN s pulse, re-priced
    dose_u8: float  # [kJ/m^2] at 2 Hz / 1 g, u = 8 km/s
    face_100um_u5: float  # [K]
    face_100um_u8: float
    floor_u8_um: float  # thinnest Kapton that holds at u = 8
    floor_u5_um: float


def push_rows() -> list[PushRow]:
    out = []
    toll = fast_wave.IMPACTOR_TOLL["Kapton"]
    for w in PUSH_SPEEDS:
        for k in PUSH_KS:
            ve, e, loss = fast_wave.push_ve(k, w, toll)
            c = fast_wave.Closure(f"Jupiter chain {w / 1e3:.0f}", w, k, ve)
            out.append(
                PushRow(
                    w / 1e3,
                    k,
                    fast_wave.calcs_push_ve(k, w) / 1e3,
                    ve / 1e3,
                    ve / fast_wave.calcs_push_ve(k, w),
                    e,
                    loss,
                    fast_wave.PULSE / (k * ve),
                    fast_wave.dose(c, 8.0e3) / 1e3,
                    fast_wave.kapton(c, 5.0e3, 100e-6).face_max,
                    fast_wave.kapton(c, 8.0e3, 100e-6).face_max,
                    1e6 * fast_wave.kapton_floor(c, 8.0e3),
                    1e6 * fast_wave.kapton_floor(c, 5.0e3),
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
    win = window_rows()
    for r in win:
        print(
            f"  window {r.closure:<26} {r.w:5.1f} km/s k {r.k:4.2f} (c {r.c_2m:.3f}, coarse"
            f" {r.c_2m_coarse:.3f}) u {r.u:2.0f}: dose {r.dose:4.0f} kJ/m^2, 100 um face"
            f" {r.face_100um_k:4.0f} K | Kapton {r.floor_um:5.0f}-{r.ceiling_um:4.0f} um on"
            f" {r.face:<14} {'OPEN' if r.window else 'shut'}"
        )
    far = far_rows()
    for f in far:
        print(
            f"  2.5 F {f.route:<20} {f.theta:5.1f} deg k {f.k:3.0f} | J/mw {f.capture_15:.3f} ->"
            f" {f.capture_25:.3f} | v_e {f.hole_ve_15:5.2f} -> {f.hole_ve_25:5.2f} km/s"
            f" ({f.isp_25:4.0f} s) | MR {f.mass_ratio_15:.3f} -> {f.mass_ratio_25:.3f} | hole at r"
            f" {f.hole_r_15:.2f} -> {f.hole_r_25:.2f} m: impulse {f.hole_share_15:.2%} ->"
            f" {f.hole_share_25:.2%}, gas {f.hole_gas_15:.2%} -> {f.hole_gas_25:.2%} | rim peak"
            f" x{f.rim_peak_ratio:.2f}, energy x{f.rim_energy_ratio:.2f} | bowl peak"
            f" x{f.plate_peak_ratio:.2f}"
        )
    head = head_rows()
    for h in head:
        print(
            f"  head-on {h.theta:5.1f} deg k {h.k:4.1f} | J/mw run {h.capture_run:.3f} (coarse"
            f" {h.capture_coarse:.3f}), calcs interp {h.capture_interpolated:.3f}"
            f" ({h.error:+.2%}) | {h.share_of_ceiling:.1%} of ceiling | hole r {h.hole_r:.2f} m:"
            f" {h.hole_share:.2%} impulse, {h.hole_gas:.2%} gas | v_e"
            f" {h.ve_66:4.1f}/{h.ve_76:4.1f}/{h.ve_85:4.1f} km/s at 66/76/85"
        )
    push = push_rows()
    for p in push:
        print(
            f"  push {p.w:4.0f} km/s k {p.k:3.1f} | calcs {p.calcs_ve:5.2f}, re-priced"
            f" {p.ve:5.2f} km/s (x{p.ratio:.3f}; e {p.e_eff:.3f}, toll {p.toll:.1%}) |"
            f" {p.impactor_kg:4.1f} kg | dose u8 {p.dose_u8:3.0f} kJ/m^2, 100 um face"
            f" {p.face_100um_u5:4.0f} K (u5) / {p.face_100um_u8:4.0f} K (u8) | floor u8"
            f" {p.floor_u8_um:4.0f} um, u5 {p.floor_u5_um:4.0f} um"
        )
    _write(OUT_WINDOW, list(win))
    _write(OUT_FAR, list(far))
    _write(OUT_HEAD, list(head))
    _write(OUT_PUSH, list(push))
    print(f"python: wrote {OUT_WINDOW.parent}/band_*.csv")


if __name__ == "__main__":
    main()
