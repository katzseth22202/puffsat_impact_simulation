"""Does the spray plate's skirt wrap delaminate? (ADR-0055 skirt; the chamber study's question.)

The skirt is a thin maraging liner wrapped in fibre (`skirt_hoop.py` sizes the wrap for hoop
strain). The walled-nozzle chamber found that a bonded wrap fails through its thickness: a sharp
spike crosses the steel and reflects off the wrap's free outer surface as tension, which the
resin between plies cannot hold (64 MPa). This applies the same layered wave model
(`walled_nozzle.wall_waves.response_to`) to the skirt. The load is the solved skirt-wall pressure
of the 150 t design bowl (d/D 0.30, 2 m skirt; `make water-plate-argon-2d-shape`), scaled to the
12 MN s pulse as `plate_fatigue.merged_load` scales the floor.

Each skirt row gets its own wrap: the `skirt_hoop` thickness for 0.7% hoop strain over a 4 mm
liner. The bonded wrap's through-thickness impedance:
- carbon/epoxy: 4.8 MPa s/m, as the parent;
- aramid/epoxy: ~1380 kg/m^3 at ~2.0 km/s, ~2.8 MPa s/m (an assumed value).

The 64 MPa transverse resin limit is applied to both. Aramid bonds worse than carbon, so for
aramid it is generous.

**Merged loads only.** The skirt histories come from the effective-gamma merged solve. A failed
merge would load the skirt harder, and is not checked here.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat.walled_nozzle import wall_waves as ww
from puffsat.water_plate import skirt_hoop as sh

SHAPES = Path("data/results/water_plate/spray_2d_shape.jsonl")
OUTPUT = Path("data/results/water_plate/skirt_delamination.csv")
DESIGN = (0.30, 2.0)
PULSE_IMPULSE = 12.0e6
LINER = 0.004
ALLOWABLE = 0.007
RESIN = 64.0e6

WRAPS = {
    "carbon": (sh.CARBON, 4.8e6 / 1579.0, 1579.0),
    "aramid": (sh.ARAMID, 2000.0, 1380.0),
}


@dataclass(frozen=True)
class Row:
    wrap: str
    row: int
    z_m: float
    peak_pressure_mpa: float
    rise_us: float
    wrap_mm: float
    wrap_tension_mpa: float
    over_resin: float


def skirt_histories() -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """`(time [s], pressure[t, row] [Pa], row z [m], radius [m])` for the design bowl."""
    for line in SHAPES.open():
        rec = json.loads(line)
        if (rec["d_over_d"], rec["skirt_m"]) != DESIGN:
            continue
        unit = float(rec["unit_m"])
        tau = unit / sh.CLOSING_SPEED
        scale = PULSE_IMPULSE / (2.0 * math.pi * float(rec["wall_impulse"])) / (unit * unit * tau)
        time = np.asarray(rec["time"], dtype=np.float64) * tau
        p = np.asarray(rec["skirt_pressure"], dtype=np.float64) * scale
        z = np.asarray(rec["skirt_z"], dtype=np.float64) * unit
        return time, p, z, float(rec["r_plate"]) * unit
    raise ValueError(f"design shape {DESIGN} not in {SHAPES}")


def rise_time(time: np.ndarray, series: np.ndarray) -> float:
    """10-90% rise of the row's first peak [s]."""
    peak = float(series.max())
    i10 = int(np.argmax(series >= 0.1 * peak))
    i90 = int(np.argmax(series >= 0.9 * peak))
    return float(time[i90] - time[i10])


def evaluate() -> list[Row]:
    time, p, z, radius = skirt_histories()
    out: list[Row] = []
    for name, (material, speed, density) in WRAPS.items():
        for j in range(p.shape[1]):
            series = p[:, j]
            liner = sh.Layer(sh.MARAGING, LINER)
            t_wrap = sh.wrap_thickness(time, series, radius, liner, material, ALLOWABLE)
            layers = [ww.Layer("maraging liner", LINER, 5654.0, 8000.0, math.inf)]
            if t_wrap > 0.0:
                layers.append(ww.Layer("wrap", t_wrap, speed, density, RESIN))
            stack = ww.Stack(name, tuple(layers))

            def load(t: float, series: np.ndarray = series) -> float:
                return float(np.interp(t, time, series))

            peaks = ww.response_to(stack, load, float(time[-1]))
            tension = max((float(a.max()) for la, a in peaks if la.name == "wrap"), default=0.0)
            out.append(
                Row(
                    wrap=name,
                    row=j,
                    z_m=float(z[j] - z[0]),
                    peak_pressure_mpa=float(series.max()) / 1e6,
                    rise_us=rise_time(time, series) * 1e6,
                    wrap_mm=t_wrap * 1e3,
                    wrap_tension_mpa=tension / 1e6,
                    over_resin=tension / RESIN,
                )
            )
    return out


def main() -> None:
    rows = evaluate()
    print(
        f"{'wrap':6s} {'row':>3} {'z m':>5} {'peak MPa':>8} {'rise us':>7} {'wrap mm':>7} "
        f"{'tension MPa':>11} {'x resin':>7}"
    )
    for r in rows:
        print(
            f"{r.wrap:6s} {r.row:3d} {r.z_m:5.2f} {r.peak_pressure_mpa:8.1f} {r.rise_us:7.0f} "
            f"{r.wrap_mm:7.1f} {r.wrap_tension_mpa:11.1f} {r.over_resin:7.2f}"
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
