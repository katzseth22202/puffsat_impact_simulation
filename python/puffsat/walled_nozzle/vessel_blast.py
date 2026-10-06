"""The blast in rocket-shaped 20 m^3 chambers, read station by station along the wall.

`cargo run --release -p sweep -- --vessel-blast` solves the rod's energy, released one port-face
radius inside the port, in four axisymmetric contours (`euler2d::vessel`):
- cylinder + 45 deg cone;
- cylinder + elliptical nose;
- curved cone;
- straight cone.

It records the wall pressure along the side wall and the port face. This module reads each
station:

- **Spike.** Peak over the quasi-static pressure (QSP), scaled to the solved real-gas QSP. The
  QSP is the mean wall pressure over the cylinder/cone wall, 2.5-4 ms.
- **Hoop fatigue.** The local wall is taken as a ring of the local radius, as `skirt_hoop.py` does
  for the plate's skirt. Its period is `2 pi r / sqrt(sum(E t)/sum(rho t))`, with the parent's
  layups. It is driven by the station's own `p/QSP`, with blowdown after 4 ms. The steel is sized
  as the parent sizes it, and a 1 mm crack grows as in `chamber_fatigue.py`. Port-face stations are
  a plate in bending, not a ring, so they get the spike and wave checks only.
- **Spall and delamination.** The five stations with the highest spike in each case go through the
  layered wave model, the parent's overwrap wall with aluminium.
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

HISTORIES = Path("data/results/walled_nozzle/near_term/vessel_blast.jsonl")
OUTPUT = Path("data/results/walled_nozzle/near_term/vessel_blast.csv")
REAL_QSP = {"methane 7000 K": 495.9e5, "hydrogen 5500 K": 818.1e5}
CHAMBERS = {"methane 7000 K": cf.METHANE, "hydrogen 5500 K": cf.HYDROGEN}
QSP_WINDOW = (2.5e-3, 4.0e-3)
ZETA = 0.02
WAVE_STATIONS = 5


def ring_period(chamber: cf.Chamber, radius: float) -> float:
    """Hoop ring period [s] of the chamber's layup at a local radius."""
    stiff = sum(layer.modulus * layer.thickness for layer in chamber.layers)
    mass = sum(layer.density * layer.thickness for layer in chamber.layers)
    return 2.0 * math.pi * radius / math.sqrt(stiff / mass)


@dataclass(frozen=True)
class Row:
    shape: str
    chamber: str
    gamma: float
    qsp_model_mpa: float
    worst_side_spike_mpa: float
    worst_side_z: float
    worst_face_spike_mpa: float
    side_spike_over_qsp_median: float
    worst_overshoot: float
    worst_steel_peak_mpa: float
    min_pulses_h2: int
    min_pulses_air: int
    worst_carbon_tension_mpa: float
    carbon_over_resin: float
    worst_steel_radial_mpa: float


def evaluate(rec: dict[str, object]) -> Row:
    name = str(rec["chamber"])
    chamber = CHAMBERS[name]
    t = np.asarray(rec["time_s"], dtype=np.float64)
    p = np.asarray(rec["pressure_pa"], dtype=np.float64)  # [bin, station]
    zs = np.asarray(rec["station_z"], dtype=np.float64)
    rs = np.asarray(rec["station_r"], dtype=np.float64)
    side = zs > 0.0
    window = (t >= QSP_WINDOW[0]) & (t <= QSP_WINDOW[1])
    qsp = float(p[np.ix_(window, side)].mean())
    scale = REAL_QSP[name] / qsp
    peaks = p.max(axis=0)

    # Hoop fatigue along the side wall.
    worst_d, worst_peak, n_h2, n_air = 0.0, -math.inf, 10**6, 10**6
    for j in np.flatnonzero(side):
        period = ring_period(chamber, float(rs[j]))
        grid = np.arange(0.0, cf.PULSE_PERIOD, period / 200.0)
        f = np.interp(grid, t, p[:, j] / qsp)
        after = grid > t[-1]
        tail = float(np.mean(p[-100:, j]) / qsp)
        f[after] = tail * np.exp(-(grid[after] - t[-1]) / min(chamber.efold_s))
        x = cf.breathing_response_to(period, ZETA, grid, f)
        stress = cf.REST_STRESS + cf.STRESS_PER_STATIC * x
        cycles = cf.rainflow(stress)
        worst_d = max(worst_d, float(x.max()))
        worst_peak = max(worst_peak, float(stress.max()))
        n_h2 = min(n_h2, cf.pulses_to_failure(cycles, 1e-3, cf.rate_hydrogen, REAL_QSP[name] / 1e5))
        n_air = min(n_air, cf.pulses_to_failure(cycles, 1e-3, cf.rate_air, 0.0))

    # Spall and delamination at the hardest-hit stations.
    stack = next(s for s in ww.stacks() if s.name.startswith(name.split()[0]) and "with" in s.name)
    carbon, steel = 0.0, 0.0
    for j in np.argsort(peaks)[-WAVE_STATIONS:]:
        series = p[:, j] * scale
        start = float(t[int(np.argmax(series > 2.0 * series[0]))]) - 20e-6

        def load(time: float, series: np.ndarray = series, start: float = start) -> float:
            return float(np.interp(start + time, t, series))

        out = ww.response_to(stack, load, 0.6e-3)
        carbon = max(carbon, max(float(a.max()) for la, a in out if la.name == "IM7/8552"))
        steel = max(steel, max(float(a.max()) for la, a in out if la.name == "Cr-Mo"))

    js = int(np.argmax(np.where(side, peaks, -math.inf)))
    return Row(
        shape=str(rec["shape"]),
        chamber=name,
        gamma=float(rec["gamma"]),  # type: ignore[arg-type]
        qsp_model_mpa=qsp / 1e6,
        worst_side_spike_mpa=float(peaks[js]) * scale / 1e6,
        worst_side_z=float(zs[js]),
        worst_face_spike_mpa=float(peaks[~side].max()) * scale / 1e6 if (~side).any() else math.nan,
        side_spike_over_qsp_median=float(np.median(peaks[side])) / qsp,
        worst_overshoot=worst_d,
        worst_steel_peak_mpa=worst_peak / 1e6,
        min_pulses_h2=n_h2,
        min_pulses_air=n_air,
        worst_carbon_tension_mpa=carbon / 1e6,
        carbon_over_resin=carbon / ww.RESIN_TENSION,
        worst_steel_radial_mpa=steel / 1e6,
    )


def main() -> None:
    rows = [evaluate(json.loads(line)) for line in HISTORIES.open()]
    print(
        f"{'shape':27s} {'chamber':16s} {'g':>3} {'side spike':>10} {'at z':>5} "
        f"{'face spike':>10} {'med/QSP':>7} {'D':>5} {'N_H2':>7} {'N_air':>7} {'carbon':>7} "
        f"{'x64':>5} {'steel r':>7}"
    )
    for r in rows:
        print(
            f"{r.shape:27s} {r.chamber:16s} {r.gamma:3.1f} {r.worst_side_spike_mpa:9.0f}M "
            f"{r.worst_side_z:5.2f} {r.worst_face_spike_mpa:9.0f}M "
            f"{r.side_spike_over_qsp_median:7.1f} {r.worst_overshoot:5.2f} {r.min_pulses_h2:7d} "
            f"{r.min_pulses_air:7d} {r.worst_carbon_tension_mpa:6.0f}M {r.carbon_over_resin:5.1f} "
            f"{r.worst_steel_radial_mpa:6.0f}M"
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
