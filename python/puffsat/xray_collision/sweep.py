"""Run the collision cases and write `data/results/xray_collision/*.csv`.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.sweep

`baseline.csv`: two equal 10 kg balls (k = 1) of argon and of iron at 67 km/s, both Saha
variants, opacity scaled 0.01-10x. It is the first rough look: is a dense kinetic collision an
X-ray source at rocket speeds?
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Iterable
from pathlib import Path

from puffsat.xray_collision.fireball import CollisionResult, collide
from puffsat.xray_collision.materials import ARGON, IRON, Material

OUT_DIR = Path("data/results/xray_collision")
BALL_MASS = 10.0  # [kg], baseline
KAPPA_SCALES = (0.01, 0.1, 1.0, 10.0)


def radius_for_mass(mat: Material, mass: float) -> float:
    return float((3.0 * mass / (4.0 * math.pi * mat.rho0)) ** (1.0 / 3.0))


def write_csv(path: Path, rows: Iterable[CollisionResult]) -> list[CollisionResult]:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [f.name for f in dataclasses.fields(CollisionResult)] + ["light_of_ship_ke"]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in rows:
            w.writerow([*(getattr(r, f.name) for f in dataclasses.fields(r)), r.light_of_ship_ke])
    print(f"python: wrote {len(rows)} rows -> {path}")
    return rows


def show(rows: Iterable[CollisionResult]) -> None:
    print(
        f"  {'mat':6} {'EOS':7} {'k':>3} {'rho0':>8} {'r_imp':>8} {'kappa':>6} {'kT0':>5} "
        f"{'tau0':>8} {'race0':>8} {'f_rad':>8} {'f>30eV':>8} {'f>100eV':>8} {'kT_em':>6} "
        f"{'light/KE':>8}"
    )
    for r in rows:
        print(
            f"  {r.material:6} {'lowered' if r.lowering else 'ideal':7} {r.mass_ratio:>3g} "
            f"{r.rho0:>8.2g} {r.impactor_radius:>8.0e} {r.kappa_scale:>6g} {r.kt_start_ev:>5.1f} "
            f"{r.tau0:>8.2e} {r.race0:>8.2e} {r.f_rad:>8.2e} {r.f_euv:>8.2e} {r.f_xray:>8.2e} "
            f"{r.kt_emit_ev:>6.2f} {r.light_of_ship_ke:>8.2e}"
        )


def baseline() -> list[CollisionResult]:
    """Two equal 10 kg balls, k = 1."""
    return [
        collide(m, radius_for_mass(m, BALL_MASS), 1.0, lowering=low, kappa_scale=ks)
        for m in (ARGON, IRON)
        for low in (False, True)
        for ks in KAPPA_SCALES
    ]


def main() -> None:
    print("== baseline: two 10 kg balls, k = 1, 67 km/s")
    show(write_csv(OUT_DIR / "baseline.csv", baseline()))


if __name__ == "__main__":
    main()
