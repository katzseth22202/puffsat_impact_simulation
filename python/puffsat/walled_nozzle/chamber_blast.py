"""Read the solved blast-in-vessel wall histories against the chamber wall's three failure modes.

`cargo run --release -p sweep -- --chamber-blast` solves a 1-D spherical blast in the 20 m^3
chamber (`euler2d::sphere`). It covers both chambers and γ 1.2 / 1.4, a dry fill against the
parent's near-wall liquid share, and one charge against split deliveries. Each record holds the
wall-pressure history. This module reads it three ways:

1. **Scale.** The γ-law gas sets its own quasi-static pressure (QSP, the wall mean over 3-6 ms).
   Every history is divided by it and multiplied by the solved real-gas QSP: 496 bar methane,
   818 bar hydrogen (`chambers.csv`). Spike-to-QSP is the γ-law result; the level is the real
   chamber's.
2. **Breathing fatigue.** The normalised history drives the wall's breathing mode
   (`chamber_fatigue.breathing_response_to`), continued by the chamber's blowdown after 6 ms.
   The steel is stressed as the parent sizes it. A 1 mm crack grows on the CC2938 hydrogen curve
   for the hydrogen chamber, and on the air curve (with the hydrogen bound) for methane.
3. **Spall and delamination.** The scaled history, first 0.6 ms, drives the through-thickness
   wave model (`wall_waves.response_to`) on the parent's stack with aluminium. Reported: the
   peak tension in the steel and in the carbon, and the through-thickness strength the overwrap
   would need, against the resin's 64 MPa.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat.walled_nozzle import chamber_fatigue as cf
from puffsat.walled_nozzle import wall_waves as ww

HISTORIES = Path("data/results/walled_nozzle/near_term/chamber_blast.jsonl")
OUTPUT = Path("data/results/walled_nozzle/near_term/chamber_blast.csv")

REAL_QSP = {"methane 7000 K": 495.9e5, "hydrogen 5500 K": 818.1e5}
CHAMBERS = {"methane 7000 K": cf.METHANE, "hydrogen 5500 K": cf.HYDROGEN}
QSP_WINDOW = (3.0e-3, 6.0e-3)
WAVE_WINDOW = 0.6e-3
ZETA = 0.02


@dataclass(frozen=True)
class Row:
    chamber: str
    gamma: float
    fill: str
    liquid_fraction: float
    shell_m: float
    graded: bool
    delivery: str
    spike_over_qsp: float
    spike_mpa: float
    arrival_us: float
    rise_10_90_us: float
    overshoot: float
    steel_peak_mpa: float
    pulses_h2: int
    pulses_air: int
    steel_radial_tension_mpa: float
    carbon_tension_mpa: float
    carbon_over_resin: float


def rise_time(t: np.ndarray, p: np.ndarray, qsp: float) -> float:
    """10-90% rise of the 50 µs running mean toward the QSP [s]."""
    width = max(1, round(50e-6 / float(t[1] - t[0])))
    smooth = np.convolve(p, np.ones(width) / width, mode="same")
    above10 = np.flatnonzero(smooth >= 0.1 * qsp)
    above90 = np.flatnonzero(smooth >= 0.9 * qsp)
    if above10.size == 0 or above90.size == 0:
        return math.nan
    return float(t[above90[0]] - t[above10[0]])


def evaluate(rec: dict[str, object]) -> Row:
    chamber_name = str(rec["chamber"])
    chamber = CHAMBERS[chamber_name]
    t = np.asarray(rec["time_s"], dtype=np.float64)
    p = np.asarray(rec["pressure_pa"], dtype=np.float64)
    window = (t >= QSP_WINDOW[0]) & (t <= QSP_WINDOW[1])
    qsp = float(p[window].mean())
    scale = REAL_QSP[chamber_name] / qsp

    # Breathing: the normalised history on the SDOF grid, then the chamber's blowdown.
    period = chamber.breathing_period()
    grid = np.arange(0.0, cf.PULSE_PERIOD, period / 400.0)
    f_hist = np.interp(grid, t, p / qsp)
    after = grid > t[-1]
    tail = float(np.mean(p[-200:]) / qsp)
    f_hist[after] = tail * np.exp(-(grid[after] - t[-1]) / min(chamber.efold_s))
    x = cf.breathing_response_to(period, ZETA, grid, f_hist)
    stress = cf.REST_STRESS + cf.STRESS_PER_STATIC * x
    cycles = cf.rainflow(stress)
    h2 = cf.pulses_to_failure(cycles, 1e-3, cf.rate_hydrogen, REAL_QSP[chamber_name] / 1e5)
    air = cf.pulses_to_failure(cycles, 1e-3, cf.rate_air, 0.0)

    # Spall and delamination: the scaled history through the parent's stack with aluminium.
    stack = next(
        s for s in ww.stacks() if s.name.startswith(chamber_name.split()[0]) and "with" in s.name
    )
    start = float(rec["first_arrival_s"]) - 20e-6

    def load(time: float) -> float:
        return float(np.interp(start + time, t, p)) * scale

    peaks = ww.response_to(stack, load, WAVE_WINDOW)
    steel = max(float(arr.max()) for layer, arr in peaks if layer.name == "Cr-Mo")
    carbon = max(float(arr.max()) for layer, arr in peaks if layer.name == "IM7/8552")

    return Row(
        chamber=chamber_name,
        gamma=float(rec["gamma"]),  # type: ignore[arg-type]
        fill=str(rec["fill"]),
        liquid_fraction=float(rec["liquid_fraction"]),  # type: ignore[arg-type]
        shell_m=float(rec["shell_m"]),  # type: ignore[arg-type]
        graded=bool(rec["graded"]),
        delivery=str(rec["delivery"]),
        spike_over_qsp=float(p.max()) / qsp,
        spike_mpa=float(p.max()) * scale / 1e6,
        arrival_us=float(rec["first_arrival_s"]) * 1e6,  # type: ignore[arg-type]
        rise_10_90_us=rise_time(t, p, qsp) * 1e6,
        overshoot=float(x.max()),
        steel_peak_mpa=float(stress.max()) / 1e6,
        pulses_h2=h2,
        pulses_air=air,
        steel_radial_tension_mpa=steel / 1e6,
        carbon_tension_mpa=carbon / 1e6,
        carbon_over_resin=carbon / ww.RESIN_TENSION,
    )


def main() -> None:
    rows = [evaluate(json.loads(line)) for line in HISTORIES.open()]
    print(
        f"{'chamber':16s} {'g':>3} {'fill':12s} {'liq':>3} {'shell':>5} {'grd':>3} "
        f"{'delivery':8s} {'spk/Q':>6} {'spk MPa':>7} {'rise':>5} {'D':>5} {'N_H2':>8} "
        f"{'N_air':>8} {'steel r':>7} {'carbon':>7} {'C/64':>5}"
    )
    for r in rows:
        print(
            f"{r.chamber:16s} {r.gamma:3.1f} {r.fill:12s} {r.liquid_fraction:3.1f} "
            f"{r.shell_m:5.1f} {int(r.graded):3d} {r.delivery:8s} {r.spike_over_qsp:6.2f} "
            f"{r.spike_mpa:7.0f} {r.rise_10_90_us:5.0f} {r.overshoot:5.2f} {r.pulses_h2:8d} "
            f"{r.pulses_air:8d} {r.steel_radial_tension_mpa:7.0f} {r.carbon_tension_mpa:7.0f} "
            f"{r.carbon_over_resin:5.2f}"
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Row.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([getattr(r, f) for f in fields])
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
