"""Argon gas bags at 66 km/s, 1 moving : 2 at rest: does most of the heat leave as light?

    PYTHONPATH=python uv run python -m puffsat.xray_collision.argon_bags

Writes `data/results/xray_collision/argon_bags_k2.csv`. Two sweeps, both with the impactor and
target at the same density at the moment of impact:

- density, for a 10 kg impactor into 20 kg, from 1e-2 down to 1e-5 kg/m^3, with the gray
  model alongside. Denser starts (an intact 10 atm bag at 16 kg/m^3, and 1 and 0.1 kg/m^3)
  shock past the multigroup pull's 0.1 kg/m^3 top and are run gray only (`emission` "gray");
- bag mass from 0.1 kg to 1000 kg at 1e-2, 1e-3 and 1e-4 kg/m^3.

Each row is run with the Planck and the Rosseland group mean in emission (`multigroup`). It also
prints the density ripple of a cubic lattice of small bags burst before impact, each spread into
a Gaussian cloud of width sigma at spacing d.
"""

from __future__ import annotations

import csv
import dataclasses
import itertools
import math
from pathlib import Path

import numpy as np

from puffsat.xray_collision.fireball import collide
from puffsat.xray_collision.materials import ARGON
from puffsat.xray_collision.multigroup import (
    V_REL_ARGON_BAGS,
    Emission,
    MultigroupResult,
    collide_groups,
)

OUT = Path("data/results/xray_collision/argon_bags_k2.csv")
MASS_RATIO = 2.0  # 1 moving : 2 at rest in the ship frame
BAG_MASS = 10.0  # [kg]
DENSE = (16.0, 1.0, 1.0e-1)  # [kg/m^3] gray only
DENSITIES = (1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 1.0e-5)  # [kg/m^3]
MASSES = (0.1, 1.0, 10.0, 100.0, 1000.0)  # [kg]
MASS_SWEEP_DENSITIES = (1.0e-2, 1.0e-3, 1.0e-4)
TOPS_RHO_MAX = 0.1  # [kg/m^3] the multigroup pull's top density
EMISSIONS: tuple[Emission, ...] = ("planck", "rosseland")


def lattice_ripple(sigma_over_d: float, n: int = 21, reach: int = 3) -> tuple[float, float]:
    """(max, min) of rho/<rho> - 1 for unit-weight Gaussians of width sigma on a cubic lattice."""
    x = np.linspace(0.0, 1.0, n)
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    rho = np.zeros_like(xx)
    for i, j, k in itertools.product(range(-reach, reach + 1), repeat=3):
        rho += np.exp(-((xx - i) ** 2 + (yy - j) ** 2 + (zz - k) ** 2) / (2.0 * sigma_over_d**2))
    mean = float(rho.mean())
    return float(rho.max()) / mean - 1.0, float(rho.min()) / mean - 1.0


def gray_f_rad(rho0: float, m_impactor: float) -> float:
    r = (3.0 * m_impactor / (4.0 * math.pi * rho0)) ** (1.0 / 3.0)
    return collide(ARGON, r, MASS_RATIO, V_REL_ARGON_BAGS, rho0=rho0).f_rad


def run() -> list[dict[str, float | str | bool]]:
    rows: list[dict[str, float | str | bool]] = []

    def add(res: MultigroupResult, sweep: str) -> None:
        row: dict[str, float | str | bool] = {"sweep": sweep, **dataclasses.asdict(res)}
        row["light_of_ship_ke"] = res.light_of_ship_ke
        row["f_rad_gray"] = gray_f_rad(res.rho0, res.m_impactor)
        rows.append(row)

    for rho0 in DENSE:
        r = (3.0 * BAG_MASS / (4.0 * math.pi * rho0)) ** (1.0 / 3.0)
        g = collide(ARGON, r, MASS_RATIO, V_REL_ARGON_BAGS, rho0=rho0)
        rows.append(
            {
                "sweep": "density",
                "rho0": rho0,
                "m_impactor": BAG_MASS,
                "mass_ratio": MASS_RATIO,
                "v_rel": V_REL_ARGON_BAGS,
                "emission": "gray",
                "r_impactor": r,
                "heat": g.heat,
                "heat_fraction": g.heat_fraction,
                "rho_shock": g.rho_shock,
                "kt_start_ev": g.kt_start_ev,
                "light_of_ship_ke": g.light_of_ship_ke,
                "f_rad_gray": g.f_rad,
            }
        )
    for rho0 in DENSITIES:
        for em in EMISSIONS:
            add(collide_groups(rho0, BAG_MASS, MASS_RATIO, emission=em), "density")
    for rho0 in MASS_SWEEP_DENSITIES:
        for m in MASSES:
            for em in EMISSIONS:
                add(collide_groups(rho0, m, MASS_RATIO, emission=em), "mass")
    return rows


def main() -> None:
    rows = run()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[-1]), restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"python: wrote {len(rows)} rows -> {OUT}")
    print(
        f"  {'sweep':7} {'emis':9} {'M':>6} {'rho0':>7} {'R_imp':>6} {'kT0':>5} {'gray':>6} "
        f"{'f_rad':>6} {'cold':>6} {'<11.6':>6} {'-15.8':>6} {'-30':>6} {'>30':>6} "
        f"{'t_half':>8} {'KE frac':>7}"
    )
    for r in rows:
        if r["emission"] == "gray":
            print(
                f"  {r['sweep']!s:7} {'gray':9} {r['m_impactor']:>6g} {r['rho0']:>7.0e} "
                f"{r['r_impactor']:>6.2f} {r['kt_start_ev']:>5.1f} {r['f_rad_gray']:>6.3f}"
            )
            continue
        print(
            f"  {r['sweep']!s:7} {r['emission']!s:9} {r['m_impactor']:>6g} {r['rho0']:>7.0e} "
            f"{r['r_impactor']:>6.2f} {r['kt_start_ev']:>5.1f} {r['f_rad_gray']:>6.3f} "
            f"{r['f_rad']:>6.3f} {r['f_through_cold_target']:>6.3f} "
            f"{r['f_below_resonance']:>6.3f} {r['f_resonance_to_edge']:>6.3f} "
            f"{r['f_edge_to_euv']:>6.3f} {r['f_above_euv']:>6.3f} {r['t_half']:>8.1e} "
            f"{r['light_of_ship_ke']:>7.3f}"
        )
    if any(float(r["rho_shock"]) > TOPS_RHO_MAX for r in rows if r["emission"] != "gray"):
        print(f"  note: shocked states above {TOPS_RHO_MAX} kg/m^3 hold the pull's top density")
    print("  burst-bag lattice, sigma/d -> density ripple (max, min)")
    for s in (0.3, 0.35, 0.4, 0.45, 0.5):
        hi, lo = lattice_ripple(s)
        print(f"    {s:.2f}: {hi:+.2f} {lo:+.2f}")


if __name__ == "__main__":
    main()
