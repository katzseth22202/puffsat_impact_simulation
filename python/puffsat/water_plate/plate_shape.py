"""The 150 t plate's shape: dish depth against skirt height (ADR-0055, Q62-Q63).

For each shape of `make water-plate-argon-2d-shape` (k = 8.52), scale the recorded wall histories
to the real pulse and size the plate inside the parent ledger's 150 t:

* **the skirt** is fixed to the vehicle (the floor slides inside it, Q58): a 2 mm maraging liner
  with a carbon hoop wrap sized row by row (`skirt_hoop`);
* **the bowl's band** is the floor's own steel, carrying the hoop load of each riser ring with the
  steel of that ring as its liner; carbon wrap is added where the steel alone strains too far;
* **the floor** takes what is left of 150 t, spread over the dish's area. It and the band's wrap
  are the sprung mass, so the gas spring's stroke follows from `J = 8 m nu s` at 4 Hz.

η_jet carries the full-physics 1-D share at k = 8.5 (`spray_k.jsonl`) times each shape's 2-D
containment, unmixed and premixed. Effective-gamma 2-D, scaled to the real pulse.
"""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.water_plate import skirt_hoop as sh

SHAPES = Path("data/results/water_plate/spray_2d_shape.jsonl")
K_SWEEP = Path("data/results/water_plate/spray_k.jsonl")
OUTPUT = Path("data/results/water_plate/plate_shape.csv")

#: The pulse (Q59), the plate budget (the parent ledger's `PLATE_MASS`) and the gas spring's rate.
PULSE_IMPULSE = 12.0e6
PLATE_BUDGET = 150.0e3
RATE = 4.0
#: Closing speed setting the time unit: the cold end, where a pulse lasts longest.
CLOSING_SPEED = 45580.0
STRAIN = 0.007
SKIRT_LINER = sh.Layer(sh.MARAGING, 0.002)
WRAP = sh.CARBON
#: Thinnest floor the thermal run covers (15 mm held 518 K, Q57).
MIN_FLOOR = 0.015
MAX_STROKE = 4.0


def paraboloid_area(radius: float, depth: float) -> float:
    """Surface area of `z = depth (r / radius)^2` out to `radius` [m^2]; `pi radius^2` when flat."""
    if depth <= 0.0:
        return math.pi * radius**2
    rim = (radius**2 + 4.0 * depth**2) ** 1.5
    return float(math.pi * radius / (6.0 * depth**2) * (rim - radius**3))


def share_1d(k: float) -> dict[str, float]:
    """The full-physics 1-D share of the ceiling at `k` (argon, cold end), by mixing bound."""
    out: dict[str, float] = {}
    with K_SWEEP.open() as fh:
        for line in fh:
            r = json.loads(line)
            if r["k"] == k and r["material"] == "argon" and r["w"] == CLOSING_SPEED:
                out[r["bound"]] = float(r["share_1d"])
    return out


def eta_jet(share: float, k: float) -> float:
    s = math.sqrt(1.0 + k)
    return (share * (1.0 + s) - 1.0) / s


@dataclass(frozen=True)
class Sized:
    d_over_d: float
    skirt_m: float
    containment: float
    eta_unmixed: float
    eta_premixed: float
    skirt_t: float
    band_wrap_t: float
    floor_t: float
    floor_mm: float
    stroke_m: float
    worst_skirt_kPa_s: float
    worst_band_kPa_s: float
    feasible: bool


def size(rec: dict[str, object], shares: dict[str, float]) -> Sized:
    unit = float(rec["unit_m"])  # type: ignore[arg-type]
    radius = float(rec["r_plate"]) * unit  # type: ignore[arg-type]
    d_over_d = float(rec["d_over_d"])  # type: ignore[arg-type]
    tau = unit / CLOSING_SPEED
    scale = PULSE_IMPULSE / (2.0 * math.pi * float(rec["wall_impulse"]))  # type: ignore[arg-type]
    to_pa = scale / (unit * unit * tau)
    time = np.asarray(rec["time"], dtype=np.float64) * tau

    # The skirt: rows of 2 mm liner + carbon at the real radius.
    skirt_z = np.asarray(rec["skirt_z"], dtype=np.float64) * unit
    skirt_p = np.asarray(rec["skirt_pressure"], dtype=np.float64) * to_pa
    skirt_t, worst_skirt = 0.0, 0.0
    if skirt_z.size:
        dz_s = float(skirt_z[1] - skirt_z[0])
        band = 2.0 * math.pi * radius * dz_s
        for j in range(skirt_z.size):
            wrap = sh.wrap_thickness(time, skirt_p[:, j], radius, SKIRT_LINER, WRAP, STRAIN)
            skirt_t += band * (SKIRT_LINER.thickness * sh.MARAGING.density + wrap * WRAP.density)
        impulse = np.sum(skirt_p[:-1] * np.diff(time)[:, None], axis=0)
        worst_skirt = float(impulse.max())

    # The band: each riser ring's own floor steel is its liner; iterate the floor thickness.
    band_z = np.asarray(rec["band_z"], dtype=np.float64) * unit
    band_r = np.asarray(rec["band_r"], dtype=np.float64) * unit
    band_p = np.asarray(rec["band_pressure"], dtype=np.float64) * to_pa
    depth = d_over_d * 2.0 * radius
    area = paraboloid_area(radius, depth)
    worst_band = 0.0
    if band_z.size:
        impulse = np.sum(band_p[:-1] * np.diff(time)[:, None], axis=0)
        worst_band = float(impulse.max())
    floor_mass = PLATE_BUDGET - skirt_t
    wrap_mass = 0.0
    for _ in range(6):
        t_floor = (floor_mass - wrap_mass) / (sh.MARAGING.density * area)
        wrap_mass = 0.0
        if band_z.size:
            dz_b = float(band_z[1] - band_z[0]) if band_z.size > 1 else 0.0
            for j in range(band_z.size):
                r = float(band_r[j])
                # Meridian length per unit height: the ring's steel per dz is t ds/dz.
                drdz = radius**2 / (2.0 * depth * r) if r > 0.0 else 0.0
                liner = sh.Layer(sh.MARAGING, t_floor * math.sqrt(1.0 + drdz**2))
                wrap = sh.wrap_thickness(time, band_p[:, j], r, liner, WRAP, STRAIN)
                wrap_mass += 2.0 * math.pi * r * dz_b * wrap * WRAP.density
    t_floor = (floor_mass - wrap_mass) / (sh.MARAGING.density * area)
    stroke = PULSE_IMPULSE / (8.0 * floor_mass * RATE)
    containment = float(rec["impulse_vs_1d"])  # type: ignore[arg-type]
    k = float(rec["k"])  # type: ignore[arg-type]
    return Sized(
        d_over_d=d_over_d,
        skirt_m=float(rec["skirt_m"]),  # type: ignore[arg-type]
        containment=containment,
        eta_unmixed=eta_jet(shares["stratified"] * containment, k),
        eta_premixed=eta_jet(shares["premixed"] * containment, k),
        skirt_t=skirt_t / 1e3,
        band_wrap_t=wrap_mass / 1e3,
        floor_t=(floor_mass - wrap_mass) / 1e3,
        floor_mm=t_floor * 1e3,
        stroke_m=stroke,
        worst_skirt_kPa_s=worst_skirt / 1e3,
        worst_band_kPa_s=worst_band / 1e3,
        feasible=t_floor >= MIN_FLOOR and stroke <= MAX_STROKE,
    )


def main() -> None:
    shares = share_1d(8.5)
    with SHAPES.open() as fh:
        records = [json.loads(line) for line in fh]
    with ProcessPoolExecutor() as pool:
        rows = list(pool.map(size, records, [shares] * len(records)))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        for r in rows:
            writer.writerow(asdict(r))
    for r in sorted(rows, key=lambda r: -r.eta_unmixed):
        print(
            f"dish {r.d_over_d:.2f} + skirt {r.skirt_m:2.0f} m: containment {r.containment:.3f}, "
            f"eta_jet {r.eta_unmixed:.3f}-{r.eta_premixed:.3f}; skirt {r.skirt_t:5.1f} t, band "
            f"wrap {r.band_wrap_t:5.1f} t, floor {r.floor_t:5.1f} t ({r.floor_mm:4.1f} mm), "
            f"stroke {r.stroke_m:.2f} m; worst skirt {r.worst_skirt_kPa_s:5.1f}, band "
            f"{r.worst_band_kPa_s:5.1f} kPa s; {'ok' if r.feasible else 'OUT'}"
        )


if __name__ == "__main__":
    main()
