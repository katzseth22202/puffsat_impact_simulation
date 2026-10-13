"""A wall of tiles on a sealed gas cushion: hydrocode ask B4 (2026-10-12; parent ADR 0054).

The parent proposes an Orion-style cushioned wall: inner tiles of 10 mm Cr-Mo on a sealed gas
layer (6-8 cm, precharged to 60-80 bar, orifice-damped) against a backwall that carries the
pressure. Its one-degree-of-freedom estimate assumed a uniform 2-5 kPa s spike. Here each tile
is driven by the **solved** wall history at its station (`zone_loads`, rod and plug as material,
5 kg), continued past the 1.2 ms solve by the real-gas blowdown (`blowdown`, adiabatic).

# The tile

Radial, one degree of freedom per station `x` (outward, into the cushion), areal mass `mu`:

    mu x'' = p_c(t) - p_k(x) - c x',      p_k = p_0 (d / (d - x))^gamma_k,    x >= 0.

- The tile sits on **inner stops** at `x = 0` until the chamber beats the precharge. When it
  comes back onto them its speed is recorded as a stop impact (a slam, inelastic).
- The cushion is adiabatic (milliseconds); methane precharge, `gamma_k = 1.3`.
- `c` is the parent's linear orifice damping. It also loads the backwall.
- **Bottoming:** the tile reaching `d - 5 mm` is a failure; the cushion gas there is at
  `p_0 (d/5 mm)^gamma`, which no seal holds.

# What it reports, per zone

- **travel** `x_max` and the stop-impact speed;
- **backwall:** the cushion's peak `p_k + c x'` over the QSP, and the backwall ring's own
  overshoot `D_back` under that load (`chamber_fatigue.breathing_response_to`, the dry wrap's
  breathing period at the backwall radius);
- **tile stresses:** through-thickness tension from the spike (`wall_waves`, a free back face,
  as `wall_zones.spall`), and bending from the spike's unevenness across one tile (free-free
  beam, the first three modes, impulsive: the spike lasts microseconds, the modes milliseconds).
  Each is turned into a crack life on the CC2938 hydrogen curve: the tile faces the gas;
- **seal gaps:** the hoop gap a radial travel opens (`x w / R`), the step between neighbouring
  tiles along the meridian, their sliding speed, and the largest pressure across the seal;
- **energy:** the work the gas does on the tiles over the pulse, as a share of the pulse
  energy: what the cushion takes out of the exhaust.

# Limits

- Stations are independent, coupled only through the bending estimate.
- The tile's own bending under the cushion is impulsive only; the quasi-static load across a
  tile is uniform (same chamber and cushion pressure on both faces).
- Sloshing after the 1.2 ms solve is not carried: the tail is the 0-D real-gas decay.
- Linear damping, as the parent's model; an orifice is quadratic.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import chamber_fatigue as cf
from puffsat.walled_nozzle import dry_wrap as dw
from puffsat.walled_nozzle import near_term as nt
from puffsat.walled_nozzle import wall_zones as wz
from puffsat.walled_nozzle import zone_loads as zl

OUTPUT = nt.OUTPUT_DIR / "cushion.csv"
STEEL_RHO = 7830.0
STEEL_E = 200e9
POISSON = 0.3
GAMMA_K = 1.3
BOTTOM_CLEARANCE = 5e-3
T_END = 60e-3
DT_TAIL = 5e-6
TILE_LENGTH = 0.5
"""Default tile pitch along the meridian and around the hoop [m] (`Design.pitch`)."""
ZETA_BACK = 0.05
CRACK_START = 0.5e-3

_SPALL: dict[tuple[str, float, str], float] = {}
_IMPULSE: dict[str, NDArray[np.float64]] = {}


@dataclass(frozen=True)
class Design:
    tile_mm: float
    gap_cm: float
    #: Precharge over the chamber's QSP (the parent's 60-80 bar on 94 is 0.64-0.85).
    precharge: float
    damping: float  # Pa s/m
    #: Tile pitch along the meridian and around the hoop [m].
    pitch: float = TILE_LENGTH

    @property
    def mu(self) -> float:
        return STEEL_RHO * self.tile_mm * 1e-3

    def label(self) -> str:
        size = "" if self.pitch == TILE_LENGTH else f" {100 * self.pitch:.0f} cm"
        return (
            f"{self.tile_mm:.0f} mm{size} tile, {self.gap_cm:.0f} cm at {self.precharge:.2f} "
            f"QSP, c {self.damping:.0e}"
        )


def designs() -> list[Design]:
    out = []
    for tile in (10.0, 25.0, 50.0):
        for gap in (6.0, 8.0, 15.0, 30.0):
            for pre in (0.64, 0.85):
                for c in (0.0, 1e5, 1e6):
                    out.append(Design(tile, gap, pre, c))
    # Tiles no wider than the spike's footprint see a near-uniform load and barely bend.
    for pitch in (0.1, 0.2):
        for tile in (10.0, 25.0, 50.0):
            out.append(Design(tile, 30.0, 0.85, 1e5, pitch))
    return out


@dataclass(frozen=True)
class Drive:
    """The chamber pressure at every station of one chamber, at 5 kg, on a fixed time grid."""

    case: str
    z: NDArray[np.float64]
    r: NDArray[np.float64]
    zone: list[str]
    time: NDArray[np.float64]
    pressure: NDArray[np.float64]  # [t, station]
    qsp: float
    #: The solve's own grid ends here; the tail is the 0-D decay.
    solve_end: float
    h: wz.History
    pulse_energy: float


def drive(case: zl.Case) -> Drive:
    """The solved history, then the real-gas adiabatic blowdown `p(t)/p0` times the QSP."""
    h = zl.history(case)
    res = bd.run_spec(case.spec, area_ratios=(100.0,), report=lambda _: None)
    rows = res.runs["none"]
    t0d = np.array([r.time for r in rows])
    p0d = np.array([r.pressure for r in rows]) / rows[0].pressure
    t_tail = np.arange(float(h.time[-1]) + DT_TAIL, T_END, DT_TAIL)
    tail = h.qsp * np.interp(t_tail, t0d, p0d)
    # Blend the solve's last 0.2 ms into the tail's level, so the splice is not a step.
    p = h.pressure.copy()
    blend = h.time > h.time[-1] - 0.2e-3 * h.length_scale
    w = (h.time[blend] - h.time[blend][0]) / (h.time[-1] - h.time[blend][0])
    level = h.qsp * float(np.interp(h.time[-1], t0d, p0d))
    p[blend] = p[blend] * (1.0 - w[:, None]) + level * w[:, None]
    zones = [wz.zone_of(float(z) / h.length_scale) for z in h.z]
    return Drive(
        case.name,
        h.z,
        h.r,
        zones,
        np.r_[h.time, t_tail],
        np.vstack([p, np.repeat(tail[:, None], h.z.size, axis=1)]),
        h.qsp,
        float(h.time[-1]),
        h,
        nt.PULSE_ENERGY * case.spec.energy_ratio,
    )


@dataclass(frozen=True)
class Motion:
    x_max: NDArray[np.float64]
    v_max: NDArray[np.float64]
    stop_speed: NDArray[np.float64]
    back_max: NDArray[np.float64]  # peak backwall pressure [Pa]
    back: NDArray[np.float64]  # backwall pressure history [t, station], decimated
    back_time: NDArray[np.float64]
    bottomed: NDArray[np.bool_]
    seal_dp: NDArray[np.float64]  # largest chamber-minus-cushion pressure [Pa]
    step_max: float  # largest radial step between tiles one pitch apart [m]
    slide_max: float  # largest relative radial speed between them [m/s]
    work: NDArray[np.float64]  # work done by the gas on the tile [J/m^2]


def move(d: Drive, design: Design) -> Motion:
    """Integrate every station's tile at once (symplectic Euler; the spike is impulsive)."""
    n = d.z.size
    gap = design.gap_cm * 1e-2
    p0 = design.precharge * d.qsp
    x = np.zeros(n)
    v = np.zeros(n)
    x_max = np.zeros(n)
    v_max = np.zeros(n)
    stop_speed = np.zeros(n)
    back_max = np.zeros(n)
    bottomed = np.zeros(n, dtype=bool)
    seal_dp = np.zeros(n)
    work = np.zeros(n)
    s = np.r_[0.0, np.cumsum(np.hypot(np.diff(d.z), np.diff(d.r)))]
    pitch = max(1, round(design.pitch / max(float(np.median(np.diff(s))), 1e-6)))
    step_max = slide_max = 0.0
    back_hist, back_t = [], []
    limit = gap - BOTTOM_CLEARANCE
    for k in range(1, d.time.size):
        dt = float(d.time[k] - d.time[k - 1])
        pc = d.pressure[k]
        pk = p0 * (gap / np.maximum(gap - x, BOTTOM_CLEARANCE)) ** GAMMA_K
        force = pc - pk - design.damping * v
        v_new = v + dt * force / design.mu
        x_new = x + dt * v_new
        on_stop = x_new <= 0.0
        stop_speed = np.where(on_stop & (v_new < 0.0), np.maximum(stop_speed, -v_new), stop_speed)
        x_new = np.where(on_stop, 0.0, x_new)
        v_new = np.where(on_stop, 0.0, v_new)
        work += pc * (x_new - x)
        x, v = x_new, v_new
        bottomed |= x >= limit
        x = np.minimum(x, limit)
        x_max = np.maximum(x_max, x)
        v_max = np.maximum(v_max, np.abs(v))
        back = p0 * (gap / (gap - x)) ** GAMMA_K + design.damping * np.maximum(v, 0.0)
        back = np.where(x > 0.0, back, p0)
        back_max = np.maximum(back_max, back)
        seal_dp = np.maximum(seal_dp, pc - back)
        if pitch < n:
            step_max = max(step_max, float(np.abs(x[pitch:] - x[:-pitch]).max()))
            slide_max = max(slide_max, float(np.abs(v[pitch:] - v[:-pitch]).max()))
        if k % 4 == 0:
            back_hist.append(back.copy())
            back_t.append(float(d.time[k]))
    return Motion(
        x_max,
        v_max,
        stop_speed,
        back_max,
        np.array(back_hist),
        np.array(back_t),
        bottomed,
        seal_dp,
        step_max,
        slide_max,
        work,
    )


def backwall_overshoot(d: Drive, m: Motion, j: int) -> float:
    """The backwall ring's peak over its static response to the QSP, under station `j`'s load.

    Its period is the dry wrap's breathing period at the backwall radius: `2 pi R / c_hoop`,
    with the wrap's hoop wave speed.
    """
    radius = float(d.r[j])
    period = 2.0 * math.pi * radius / math.sqrt(dw.KEVLAR.hoop_modulus / dw.KEVLAR.density)
    grid = np.arange(0.0, float(m.back_time[-1]), period / 200.0)
    f = np.interp(grid, m.back_time, m.back[:, j]) / d.qsp
    return float(cf.breathing_response_to(period, ZETA_BACK, grid, f).max())


# Free-free beam modes: (beta L, sigma).
_MODES = ((4.7300408, 0.982502215), (7.8532046, 1.000777312), (10.9956078, 0.999966450))


def _phi(
    xi: NDArray[np.float64], bl: float, sg: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Mode shape and its second derivative over `xi = s/L` (the latter per `L^2`)."""
    a = bl * xi
    phi = np.cosh(a) + np.cos(a) - sg * (np.sinh(a) + np.sin(a))
    d2 = bl**2 * (np.cosh(a) - np.cos(a) - sg * (np.sinh(a) - np.sin(a)))
    return phi, d2


def bending_stress(d: Drive, design: Design, zone: str) -> float:
    """Worst bending stress [Pa] in a tile of the zone from the spike's uneven impulse.

    Each tile is a free-free plate strip of length `design.pitch`; the spike (microseconds) gives
    each mode a velocity `int I phi ds / (mu int phi^2 ds)`; the amplitude is that over the
    mode's frequency (plate bending plus the cushion as an elastic foundation), and the surface
    stress `E t / (2 (1 - nu^2)) q phi''`. The first three modes are added in magnitude, after
    the rigid translation and rotation (a least-squares line through the tile's impulse) are taken
    out. Whole tiles are slid along the zone a station at a time.
    """
    h = d.h
    s = np.r_[0.0, np.cumsum(np.hypot(np.diff(h.z), np.diff(h.r)))]
    if d.case not in _IMPULSE:
        _IMPULSE[d.case] = np.array([zl._spike(h, j)[0] for j in range(h.z.size)])
    impulse = _IMPULSE[d.case]
    t = design.tile_mm * 1e-3
    rigidity = STEEL_E * t**3 / (12.0 * (1.0 - POISSON**2))
    foundation = GAMMA_K * design.precharge * d.qsp / (design.gap_cm * 1e-2)
    idx = np.flatnonzero([z == zone for z in d.zone])
    worst = 0.0
    for start in idx:
        length = design.pitch
        if s[start] + length > s[-1]:
            break  # only whole tiles: on a partial one the modes are not orthogonal to rigid motion
        window = (s >= s[start]) & (s <= s[start] + length)
        if window.sum() < 6:
            continue
        xi = (s[window] - s[start]) / length
        # Rigid translation and rotation carry no bending: take them out (a least-squares line).
        flex = impulse[window] - np.polyval(np.polyfit(xi, impulse[window], 1), xi)
        total = 0.0
        for bl, sg in _MODES:
            phi, d2 = _phi(xi, bl, sg)
            q_dot = float(np.sum(flex * phi)) / (design.mu * float(np.sum(phi**2)))
            # The cushion is an elastic foundation under the plate: `k = gamma p0 / d`.
            omega = math.sqrt(
                (bl**2 / length**2) ** 2 * rigidity / design.mu + foundation / design.mu
            )
            curvature = q_dot / omega * float(np.abs(d2).max()) / length**2
            total += abs(curvature) * STEEL_E * t / (2.0 * (1.0 - POISSON**2))
        worst = max(worst, total)
    return worst


def surface_crack_pulses(stress_range: float, thickness: float, p_bar: float) -> float:
    """Pulses for a semicircular surface crack from 0.5 mm to half the tile, CC2938 hydrogen."""
    if stress_range <= 0.0:
        return math.inf
    a, n, end = CRACK_START, 0.0, 0.5 * thickness
    while a < end:
        dk = cf.surface_k(stress_range, a)
        if dk >= wz.K_FAIL:
            break
        da = cf.rate_hydrogen(dk, 0.0, p_bar)
        step = max(1e-6 * thickness, 0.02 * a)
        n += step / max(da, 1e-30)
        a += step
        if n > 1e9:
            return math.inf
    return n


@dataclass(frozen=True)
class Row:
    case: str
    design: str
    zone: str
    qsp_bar: float
    travel_cm: float
    bottomed: bool
    stop_speed: float
    back_over_qsp: float
    backwall_d: float
    spall_over: float
    bending_mpa: float
    bending_pulses: float
    hoop_gap_mm: float
    step_mm: float
    slide_speed: float
    seal_dp_gpa: float
    energy_share: float

    @property
    def survives(self) -> bool:
        return (
            not self.bottomed
            and self.bending_pulses >= wz.REQUIRED_PULSES
            and self.spall_over < 1.0
        )


def evaluate(d: Drive, design: Design) -> list[Row]:
    m = move(d, design)
    t = design.tile_mm * 1e-3
    out = []
    s = np.r_[0.0, np.cumsum(np.hypot(np.diff(d.z), np.diff(d.r)))]
    ds = np.gradient(s)
    work_total = float(np.sum(2.0 * math.pi * d.r * ds * m.work))
    for zone in dict.fromkeys(d.zone):
        idx = np.flatnonzero([z == zone for z in d.zone])
        j = int(idx[int(np.argmax(m.x_max[idx]))])
        key = (d.case, design.tile_mm, zone)
        if key not in _SPALL:
            _SPALL[key] = wz.spall("steel", t, d.h, zone)[0]
        tension = _SPALL[key]
        bend = bending_stress(d, design, zone)
        out.append(
            Row(
                case=d.case,
                design=design.label(),
                zone=zone,
                qsp_bar=d.qsp / 1e5,
                travel_cm=float(m.x_max[idx].max()) * 100,
                bottomed=bool(m.bottomed[idx].any()),
                stop_speed=float(m.stop_speed[idx].max()),
                back_over_qsp=float(m.back_max[idx].max()) / d.qsp,
                backwall_d=backwall_overshoot(d, m, j),
                spall_over=tension / wz.CRMO.spall,
                bending_mpa=bend / 1e6,
                bending_pulses=surface_crack_pulses(bend, t, d.qsp / 1e5),
                hoop_gap_mm=float(m.x_max[idx].max()) * design.pitch / float(d.r[idx].max()) * 1e3,
                step_mm=m.step_max * 1e3,
                slide_speed=m.slide_max,
                seal_dp_gpa=float(m.seal_dp[idx].max()) / 1e9,
                energy_share=work_total / d.pulse_energy,
            )
        )
    return out


def write(rows: Sequence[Row]) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = [*Row.__dataclass_fields__, "survives"]
    with OUTPUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in rows:
            w.writerow([getattr(r, f) for f in fields])


def main() -> None:
    rows: list[Row] = []
    for case in zl.cases():
        if not case.path.exists() or "0.25 cm" in case.name:
            continue
        d = drive(case)
        print(f"\n{case.name}: QSP {d.qsp / 1e5:.0f} bar")
        for design in designs():
            got = evaluate(d, design)
            rows.extend(got)
            ok = all(r.survives for r in got)
            worst = max(got, key=lambda r: r.travel_cm)
            print(
                f"  {design.label():42s} {'SURVIVES' if ok else 'fails   '} travel "
                f"{worst.travel_cm:5.1f} cm{' BOTTOMS' if any(r.bottomed for r in got) else ''} "
                f"back {max(r.back_over_qsp for r in got):5.2f} D_back "
                f"{max(r.backwall_d for r in got):5.2f} bend {max(r.bending_mpa for r in got):6.0f}"
                f" MPa ({min(r.bending_pulses for r in got):9.3g} pulses) stop "
                f"{max(r.stop_speed for r in got):5.1f} m/s gap "
                f"{max(r.hoop_gap_mm for r in got):4.1f} mm step {worst.step_mm:4.1f} mm "
                f"slide {worst.slide_speed:5.1f} m/s energy {100 * got[0].energy_share:.2f}%"
            )
    write(rows)
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
