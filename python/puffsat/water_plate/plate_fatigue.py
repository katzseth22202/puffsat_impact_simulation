"""Spall and fatigue of the spray plate's maraging floor over 3,000 pulses.

The face pressure enters a ~30 mm steel floor as a stress wave. It crosses in ~5 µs, reflects off
the free back face as tension, and rings in the thickness mode (period `2L/c` ~ 10.6 µs) until
damping removes it. This module drives that wave with the recorded face histories. It asks three
questions:

1. **Spall.** Does the peak interior tension reach the measured spall strength?
2. **Initiation (stress-life).** Rainflow-count the normal-stress history at each depth and sum
   Miner damage against NI76's smooth and notched constant-life diagrams for maraging 300.
3. **Damage tolerance.** Grow an embedded crack lying parallel to the face (the spall plane,
   which `sigma_x` opens) by the Paris law. Start it at an inspectable size and check whether it
   reaches the critical size set by `K_Ic` within 3,000 pulses.

# Two loads

- **Merged** (the design intent): the 2-D bowl (d/D 0.30, 2 m skirt) scaled to the 12 MN s pulse,
  as `skirt_hoop.py` scales it. This is effective-gamma gas, sampled every ~3 µs.
- **Unmerged** (the merge fails): the 1-D stratified real-EOS pulse on a bare face at the same
  design point (`cargo run --release -p sweep -- --spray-face-history`). It rises in ~0.5 µs.

# The wave model

The model is linear elastic, in uniaxial strain, with longitudinal speed
`c = sqrt(E(1-nu)/((1+nu)(1-2nu)rho))`. The rightward and leftward stress waves are carried on a
grid with `dx = c dt`, so advection is exact. The front is a traction boundary,
`sigma(0) = -p(t)`. The back is free, `sigma(L) = 0`, which is the worst case for spall. Damping
is a quality factor `Q` on the thickness mode, swept because it is unknown: the film, the joints
and the taper all add damping. Above the HEL (~3.1-3.6 GPa calculated for 300) the face yields and
the elastic model overstates the ringing.

**Not modelled.** In-plane bending of the plate (the taper is meant to remove it), thermal
stress at the face, and the interface to the gas spring behind the floor (treated as free).

# Material data

- `K_Ic`: 66.5-68.7 MPa m^0.5 for 300 (19 mm plate), 33 for 350 (NI76 Table 12).
- Fatigue: smooth and Kt 2.4-2.5 constant-life lower bands at 1e4/1e5/1e7 cycles for 300
  (NI76 Figs. 14, 15). None given for 350.
- Spall strength: 4.1-5.4 GPa for 18Ni-9Co-5Mo, from `0.5 rho c_b W` (Mescheryakov et al.
  2017, Table 5). None found for 350.
- Paris law: `1.36e-10 dK^2.25` m/cycle, the martensitic-steel upper bound (Barsom and Rolfe,
  not re-fetched).

The fatigue curves are read off NI76's figures by eye from the lower edge of each band, to about
+/-20 MPa. The Paris constants are the standard martensitic-steel upper bound. No crack-growth
threshold is used, so every tensile cycle grows the crack. That choice is conservative.
"""

from __future__ import annotations

import csv
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

FACE_HISTORY = Path("data/results/water_plate/spray_face_history.jsonl")
SHAPE_HISTORY = Path("data/results/water_plate/spray_2d_shape.jsonl")
OUTPUT = Path("data/results/water_plate/plate_fatigue.csv")
RISE_OUTPUT = Path("data/results/water_plate/plate_fatigue_rise.csv")
RISE_TIMES = (0.5e-6, 2e-6, 5e-6, 10e-6, 20e-6)
"""Front rise times [s] for the merged pulse's sharpness sweep."""

PULSES = 3000
"""Baseline plate life in pulses (two 1,500-pulse pushes), chosen as conservative."""
PULSE_IMPULSE = 12.0e6
"""The 150 t design pulse [N s] (ADR-0055, Q59-Q62)."""
CLOSING_SPEED = 45_580.0
DESIGN_SHAPE = (0.30, 2.0)
"""(d/D, skirt m) of the design bowl."""

# ---- Steel ----------------------------------------------------------------------------------

E_STEEL = 190.0e9
NU_STEEL = 0.30
RHO_STEEL = 8000.0
K_BULK = 164.0e9
"""Bulk modulus [Pa] for the spall conversion (Thomas et al. 2018 value for steel, as in the
face-limits note)."""
FLOOR_THICKNESSES = (0.030, 0.036)
"""Floor thickness [m] of the d/D 0.30 and 0.25 bowls with a 2 m skirt (handoff, Q59-Q63)."""
Q_FACTORS = (30.0, 100.0, 1000.0)

SPALL_PULLBACK_M_S = (227.5, 239.3, 271.5, 295.6, 281.2)
"""Pull-back velocities [m/s] for 02N18K9M5-VI, which is 18Ni-9Co-5Mo with sigma_0.2 2.0 GPa, the
maraging 300 chemistry (Mescheryakov et al., Mater. Phys. Mech. 32, 152 (2017), Table 5)."""


def longitudinal_speed() -> float:
    """Uniaxial-strain (longitudinal) wave speed [m/s]."""
    modulus = E_STEEL * (1 - NU_STEEL) / ((1 + NU_STEEL) * (1 - 2 * NU_STEEL))
    return math.sqrt(modulus / RHO_STEEL)


def spall_strength(pullback: float) -> float:
    """Acoustic spall strength `0.5 rho c_b W` [Pa]."""
    return 0.5 * RHO_STEEL * math.sqrt(K_BULK / RHO_STEEL) * pullback


@dataclass(frozen=True)
class Grade:
    name: str
    k_ic: float  # [Pa m^0.5]


MARAGING_300 = Grade("maraging 300", 66.5e6)
MARAGING_350 = Grade("maraging 350", 33.0e6)

PARIS_C = 1.36e-10
"""Paris coefficient [m/cycle at dK in MPa m^0.5] (Barsom and Rolfe, martensitic steels)."""
PARIS_M = 2.25

# Constant-life lower bands, NI76 Fig. 14 (smooth) and Fig. 15 (Kt 2.4-2.5), 18Ni1700/1900.
# {cycles: (mean MPa, amplitude MPa)}.
SMOOTH: dict[float, tuple[tuple[float, ...], tuple[float, ...]]] = {
    1e4: ((0, 200, 400, 600, 800, 1000, 1200), (1140, 990, 860, 740, 620, 500, 400)),
    1e5: ((0, 200, 400, 600, 800, 1000), (810, 690, 560, 440, 320, 210)),
    1e7: ((0, 200, 400, 600, 800, 1000), (550, 400, 270, 180, 110, 60)),
}
NOTCHED: dict[float, tuple[tuple[float, ...], tuple[float, ...]]] = {
    1e4: ((0, 100, 200, 300, 400, 500, 600), (560, 515, 470, 435, 400, 365, 330)),
    1e5: ((0, 100, 200, 300, 400, 500), (365, 310, 265, 225, 195, 170)),
    1e7: ((0, 100, 200, 300, 400), (245, 200, 160, 130, 110)),
}


# ---- Loads ----------------------------------------------------------------------------------


@dataclass(frozen=True)
class Load:
    label: str
    time: NDArray[np.float64]
    pressure: NDArray[np.float64]


def merged_load(path: Path = SHAPE_HISTORY) -> Load:
    """The design bowl's hardest-hit floor band, scaled to the 12 MN s pulse."""
    for line in path.open():
        rec = json.loads(line)
        if (rec["d_over_d"], rec["skirt_m"]) != DESIGN_SHAPE:
            continue
        unit = float(rec["unit_m"])
        tau = unit / CLOSING_SPEED
        scale = PULSE_IMPULSE / (2.0 * math.pi * float(rec["wall_impulse"])) / (unit * unit * tau)
        bands = np.asarray(rec["band_pressure"], dtype=np.float64) * scale
        worst = int(np.argmax(bands.max(axis=0)))
        time = np.asarray(rec["time"], dtype=np.float64) * tau
        return Load("merged (2-D bowl)", time, bands[:, worst])
    raise ValueError(f"design shape {DESIGN_SHAPE} not in {path}")


def unmerged_loads(path: Path = FACE_HISTORY) -> list[Load]:
    """The 1-D stratified (failed-merge) face histories at the design point."""
    out = []
    for line in path.open():
        rec = json.loads(line)
        time = np.asarray(rec["time_s"], dtype=np.float64)
        keep = np.concatenate(([True], np.diff(time) > 0.0))
        label = f"unmerged {rec['material']} {rec['w'] / 1e3:.2f} km/s"
        out.append(Load(label, time[keep], np.asarray(rec["pressure_pa"], dtype=np.float64)[keep]))
    return out


# ---- Elastic wave through the thickness -----------------------------------------------------


def stress_history(
    load: Load, thickness: float, q_factor: float, cells: int = 48, ring_floor: float = 0.01
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Normal stress `sigma_x(t)` [Pa, tension positive] at `cells + 1` depths.

    Right- and left-going waves `r`, `s` with `sigma = r + s`, advected one cell per step. After
    the load ends the run continues until the ringing falls to `ring_floor` of the peak load.
    """
    c = longitudinal_speed()
    dt = thickness / (cells * c)
    period = 2.0 * thickness / c
    decay = math.exp(-math.pi / q_factor / (2 * cells))  # per step; e^-pi/Q per period
    ring_periods = math.log(1.0 / ring_floor) * q_factor / math.pi
    t_end = float(load.time[-1]) + ring_periods * period
    steps = int(t_end / dt) + 1
    p = np.interp(np.arange(steps) * dt, load.time, load.pressure, right=0.0)

    r = np.zeros(cells + 1)
    s = np.zeros(cells + 1)
    out = np.empty((steps, cells + 1))
    for n in range(steps):
        r[1:] = r[:-1] * decay
        s[:-1] = s[1:] * decay
        r[0] = -p[n] - s[0]  # traction at the face: sigma(0) = -p
        s[-1] = -r[-1]  # free back: sigma(L) = 0
        out[n] = r + s
    return np.arange(steps, dtype=np.float64) * dt, out


# ---- Rainflow (ASTM E1049, three-point) -----------------------------------------------------


def reversals(series: NDArray[np.float64], tolerance: float) -> list[float]:
    """Turning points, ignoring excursions smaller than `tolerance` (hysteresis filter)."""
    d = np.diff(series)
    turning = np.flatnonzero(np.sign(d[1:]) * np.sign(d[:-1]) < 0) + 1
    series = series[np.concatenate(([0], turning, [len(series) - 1]))]
    out = [float(series[0])]
    extreme = float(series[0])
    direction = 0
    for value in series[1:]:
        v = float(value)
        if direction == 0:
            if abs(v - out[0]) > tolerance:
                direction = 1 if v > out[0] else -1
                extreme = v
        elif (v - extreme) * direction > 0.0:
            extreme = v
        elif abs(v - extreme) > tolerance:
            out.append(extreme)
            direction = -direction
            extreme = v
    if direction != 0:
        out.append(extreme)
    return out


def rainflow(series: NDArray[np.float64], tolerance: float) -> list[tuple[float, float, float]]:
    """`(range, mean, count)` cycles, count 1 or 0.5."""
    stack: list[float] = []
    cycles: list[tuple[float, float, float]] = []
    for point in reversals(series, tolerance):
        stack.append(point)
        while len(stack) >= 3:
            x = abs(stack[-1] - stack[-2])
            y = abs(stack[-2] - stack[-3])
            if x < y:
                break
            if len(stack) == 3:
                cycles.append((y, 0.5 * (stack[-2] + stack[-3]), 0.5))
                stack.pop(0)
            else:
                cycles.append((y, 0.5 * (stack[-2] + stack[-3]), 1.0))
                last = stack.pop()
                stack.pop()
                stack.pop()
                stack.append(last)
    for a, b in pairwise(stack):
        cycles.append((abs(b - a), 0.5 * (a + b), 0.5))
    return cycles


# ---- Stress-life ----------------------------------------------------------------------------


def cycles_to_failure(
    amplitude_mpa: float,
    mean_mpa: float,
    curves: Mapping[float, tuple[tuple[float, ...], tuple[float, ...]]],
) -> float:
    """Life at `(amplitude, mean)` from a constant-life diagram, log-log between curves.

    A compressive mean is taken as zero (conservative). Beyond 1e7 the 1e5-1e7 slope is
    extended with no endurance limit (conservative for gigacycle ringing).
    """
    mean = max(mean_mpa, 0.0)
    lives = sorted(curves)
    amps = [float(np.interp(mean, *curves[n])) for n in lives]
    if amplitude_mpa <= 0.0:
        return math.inf
    log_n = [math.log10(n) for n in lives]
    log_a = [math.log10(max(a, 1.0)) for a in amps]
    la = math.log10(amplitude_mpa)
    if la >= log_a[0]:
        i = 0
    elif la <= log_a[-1]:
        i = len(lives) - 2
    else:
        i = next(k for k in range(len(lives) - 1) if log_a[k] >= la >= log_a[k + 1])
    slope = (log_n[i + 1] - log_n[i]) / (log_a[i + 1] - log_a[i])
    return float(10.0 ** (log_n[i] + slope * (la - log_a[i])))


def miner_damage(
    cycles: list[tuple[float, float, float]],
    curves: Mapping[float, tuple[tuple[float, ...], tuple[float, ...]]],
) -> float:
    return sum(
        count / cycles_to_failure(0.5 * rng / 1e6, mean / 1e6, curves)
        for rng, mean, count in cycles
    )


# ---- Damage tolerance -----------------------------------------------------------------------


def penny_k(stress: float, radius: float) -> float:
    """`K = (2/pi) sigma sqrt(pi a)` for an embedded penny crack [Pa m^0.5]."""
    return 2.0 / math.pi * stress * math.sqrt(math.pi * radius)


def critical_radius(grade: Grade, peak_tension: float) -> float:
    """Penny radius [m] at which the peak tension reaches `K_Ic`."""
    if peak_tension <= 0.0:
        return math.inf
    return math.pi * grade.k_ic**2 / (4.0 * peak_tension**2)


def tensile_cycles(cycles: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    """Cycles whose maximum is tensile; fully compressive cycles barely drive fatigue."""
    return [c for c in cycles if c[1] + 0.5 * c[0] > 0.0]


def pulses_to_critical(
    cycles: list[tuple[float, float, float]], a0: float, a_crit: float, cap: int = 100_000
) -> int:
    """Pulses for an embedded penny crack to grow from `a0` to `a_crit`, tensile part only.

    0 means the starting flaw is already critical: it fails by fast fracture on the first peak.
    """
    if a0 >= a_crit:
        return 0
    steps = []
    for rng, mean, count in cycles:
        hi = mean + 0.5 * rng
        lo = max(mean - 0.5 * rng, 0.0)
        if hi > 0.0:
            steps.append((hi - lo, count))
    if not steps:
        return cap
    a = a0
    for pulse in range(1, cap + 1):
        for d_sigma, count in steps:
            a += count * PARIS_C * (penny_k(d_sigma, a) / 1e6) ** PARIS_M
        if a >= a_crit:
            return pulse
    return cap


def sharpened(load: Load, rise: float = 0.5e-6) -> Load:
    """The same pulse with its rise compressed to `rise` [s]: zero, then a ramp to the peak."""
    i = int(np.argmax(load.pressure))
    t_peak = float(load.time[i])
    time = np.concatenate(([0.0, t_peak - rise], load.time[i:]))
    pressure = np.concatenate(([0.0, 0.0], load.pressure[i:]))
    return Load(load.label + ", 0.5 us front", time, pressure)


# ---- Run ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Row:
    load: str
    thickness_mm: float
    q_factor: float
    peak_face_gpa: float
    peak_tension_gpa: float
    tension_over_spall: float
    cycles_per_pulse: float
    miner_smooth: float
    miner_notched: float
    miner_notched_tensile: float
    pulses_notched: float
    a_crit_300_mm: float
    a_crit_350_mm: float
    pulses_crack_300: int
    pulses_crack_350: int


def evaluate(load: Load, thickness: float, q_factor: float, a0: float = 0.5e-3) -> Row:
    """Every verdict for one load, floor and damping, at the worst depth for each."""
    _, sigma = stress_history(load, thickness, q_factor)
    peak_face = float(load.pressure.max())
    peak_t = float(sigma.max())
    tol = 0.005 * peak_face
    smooth = notched = notched_t = 0.0
    worst: list[tuple[float, float, float]] = []
    for col in range(0, sigma.shape[1], 4):
        cycles = rainflow(sigma[:, col], tol)
        d_n = miner_damage(cycles, NOTCHED)
        if d_n > notched:
            worst = cycles
        smooth = max(smooth, miner_damage(cycles, SMOOTH))
        notched = max(notched, d_n)
        notched_t = max(notched_t, miner_damage(tensile_cycles(cycles), NOTCHED))
    # Crack growth on the depth with the most tensile cycling.
    by_tension = max(
        (rainflow(sigma[:, col], tol) for col in range(0, sigma.shape[1], 4)),
        key=lambda cs: sum(c[0] * c[2] for c in tensile_cycles(cs)),
    )
    spall = min(spall_strength(w) for w in SPALL_PULLBACK_M_S)
    a300 = critical_radius(MARAGING_300, peak_t)
    a350 = critical_radius(MARAGING_350, peak_t)
    return Row(
        load=load.label,
        thickness_mm=thickness * 1e3,
        q_factor=q_factor,
        peak_face_gpa=peak_face / 1e9,
        peak_tension_gpa=peak_t / 1e9,
        tension_over_spall=peak_t / spall,
        cycles_per_pulse=sum(c[2] for c in worst),
        miner_smooth=smooth * PULSES,
        miner_notched=notched * PULSES,
        miner_notched_tensile=notched_t * PULSES,
        pulses_notched=1.0 / notched_t if notched_t > 0 else math.inf,
        a_crit_300_mm=min(a300 * 1e3, 1e3),
        a_crit_350_mm=min(a350 * 1e3, 1e3),
        pulses_crack_300=pulses_to_critical(by_tension, a0, a300),
        pulses_crack_350=pulses_to_critical(by_tension, a0, a350),
    )


def main() -> None:
    merged = merged_load()
    loads = [merged, sharpened(merged), *unmerged_loads()]
    spall = [spall_strength(w) / 1e9 for w in SPALL_PULLBACK_M_S]
    print(
        f"steel c_L = {longitudinal_speed():.0f} m/s; spall strength {min(spall):.1f}-"
        f"{max(spall):.1f} GPa; {PULSES} pulses; starting flaw radius 0.5 mm\n"
    )
    rows = [evaluate(ld, h, q) for ld in loads for h in FLOOR_THICKNESSES for q in Q_FACTORS]
    print(
        f"{'load':40s} {'mm':>3} {'Q':>5} {'face':>5} {'tens':>5} {'/spall':>6} {'cyc/p':>6} "
        f"{'D smth':>8} {'D notch':>8} {'D n,ten':>8} {'N notch':>8} {'ac300':>6} {'ac350':>6} "
        f"{'Ncr300':>7} {'Ncr350':>7}"
    )
    for r in rows:
        print(
            f"{r.load:40s} {r.thickness_mm:3.0f} {r.q_factor:5.0f} {r.peak_face_gpa:5.2f} "
            f"{r.peak_tension_gpa:5.2f} {r.tension_over_spall:6.2f} {r.cycles_per_pulse:6.0f} "
            f"{r.miner_smooth:8.1e} {r.miner_notched:8.1e} {r.miner_notched_tensile:8.1e} "
            f"{r.pulses_notched:8.2g} {r.a_crit_300_mm:6.1f} {r.a_crit_350_mm:6.1f} "
            f"{r.pulses_crack_300:7d} {r.pulses_crack_350:7d}"
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Row.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([getattr(r, f) for f in fields])
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")

    # How sharp may the merged front be? The 2-D grid (12.5 cm cells) smears it to ~36 us.
    print("\nMerged pulse, 30 mm floor, front rise time against damping:")
    rise_rows = []
    for rise in RISE_TIMES:
        for q in (100.0, 1000.0):
            r = evaluate(sharpened(merged, rise), 0.030, q)
            rise_rows.append((rise * 1e6, q, r))
            print(
                f"  rise {rise * 1e6:5.1f} us  Q {q:5.0f}: tension {r.peak_tension_gpa:.2f} GPa"
                f"  D notch,tensile {r.miner_notched_tensile:.2e}"
                f"  pulses to critical crack 300 / 350: {r.pulses_crack_300} / {r.pulses_crack_350}"
            )
    with RISE_OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rise_us", "q_factor", *fields])
        for rise_us, q, r in rise_rows:
            writer.writerow([rise_us, q, *(getattr(r, f) for f in fields)])


if __name__ == "__main__":
    main()
