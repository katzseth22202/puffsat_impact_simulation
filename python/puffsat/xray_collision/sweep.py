"""Run the collision cases and write `data/results/xray_collision/*.csv`.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.sweep

`baseline.csv`: two equal 10 kg balls (k = 1) of argon and of iron at 67 km/s, both Saha
variants, opacity scaled 0.01-10x. It is the first rough look: is a dense kinetic collision an
X-ray source at rocket speeds?

`sweep_k3.csv`: a 3:1 target (25% of the ship-frame KE is lost to centre-of-mass motion), impactor
radius 1 mm, 1 cm and 10 cm as asked, extended down to 0.1 um to find where shrinking stops
helping. Both Saha variants, opacity 0.1-10x. Goal: 80% of the heat leaving as light.

`porous_k3.csv`: the same 3:1 collision started from foam or spray density rather than condensed
matter, to find where the 80% line actually is.

`near_sun.csv`: the parent's own pass-through pusher case, a 1 kg carbon projectile on a 1 kg
carbon target at 618 km/s. It runs over foam densities of 10, 30 and 100 kg/m^3 and solid
graphite, opacity 0.1-10x, with and without trapped radiation in the EOS. It checks the parent's
assumed 60% radiated fraction.
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Iterable
from pathlib import Path

from puffsat.xray_collision.fireball import V_REL_NEAR_SUN, CollisionResult, collide
from puffsat.xray_collision.materials import ARGON, CARBON, IRON, Material

OUT_DIR = Path("data/results/xray_collision")
BALL_MASS = 10.0  # [kg], baseline
KAPPA_SCALES = (0.01, 0.1, 1.0, 10.0)
MASS_RATIO = 3.0
ASKED_RADII = (1.0e-3, 1.0e-2, 1.0e-1)  # [m]
EXTENDED_RADII = (1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7)  # [m]
POROUS_RHO0 = (100.0, 10.0, 1.0, 0.1, 0.01, 1.0e-3)  # [kg/m^3]
POROUS_RADII = (1.0e-3, 1.0e-2, 1.0e-1, 1.0)  # [m]
TARGET_SHARE = 0.8
NEAR_SUN_MASS = 1.0  # [kg] each, the parent's pass-through projectile and target
NEAR_SUN_RHO0 = (10.0, 30.0, 100.0, CARBON.rho0)  # [kg/m^3]: the parent's foams and graphite


def radius_for_mass(mat: Material, mass: float, rho0: float | None = None) -> float:
    rho = mat.rho0 if rho0 is None else rho0
    return float((3.0 * mass / (4.0 * math.pi * rho)) ** (1.0 / 3.0))


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


def sweep_k3() -> list[CollisionResult]:
    """Condensed balls on a 3:1 target across impactor radius."""
    return [
        collide(m, r, MASS_RATIO, lowering=low, kappa_scale=ks)
        for m in (ARGON, IRON)
        for low in (False, True)
        for ks in (0.1, 1.0, 10.0)
        for r in ASKED_RADII + EXTENDED_RADII
    ]


def porous_k3() -> list[CollisionResult]:
    """Foam/spray starts on a 3:1 target: start density x radius, ideal Saha, TOPS opacity."""
    return [
        collide(m, r, MASS_RATIO, rho0=rho0)
        for m in (ARGON, IRON)
        for rho0 in POROUS_RHO0
        for r in POROUS_RADII
    ]


def near_sun() -> list[CollisionResult]:
    """The parent's 618 km/s carbon pass-through collision, k = 1."""
    return [
        collide(
            CARBON,
            radius_for_mass(CARBON, NEAR_SUN_MASS, rho0),
            1.0,
            V_REL_NEAR_SUN,
            rho0=rho0,
            kappa_scale=ks,
            radiation=rad,
        )
        for rad in (False, True)
        for rho0 in NEAR_SUN_RHO0
        for ks in (0.1, 1.0, 10.0)
    ]


def show_near_sun(rows: Iterable[CollisionResult]) -> None:
    print(
        f"  {'a T^4':6} {'rho0':>6} {'kappa':>5} {'R0 cm':>6} {'kT0':>5} {'tau0':>8} "
        f"{'race0':>8} {'f_rad':>6} {'>100eV':>13} {'>1keV':>6} {'kT_emit':>11}"
    )
    for r in rows:
        print(
            f"  {'in' if r.radiation else 'out':6} {r.rho0:>6g} {r.kappa_scale:>5g} "
            f"{100 * r.r0:>6.1f} {r.kt_start_ev:>5.0f} {r.tau0:>8.2e} {r.race0:>8.2e} "
            f"{r.f_rad:>6.3f} {r.f_xray_ph:>6.3f}-{r.f_xray:<6.3f} {r.f_hard:>6.3f} "
            f"{r.kt_emit_ph_ev:>5.0f}-{r.kt_emit_ev:<5.0f}"
        )


def reaching(rows: Iterable[CollisionResult]) -> list[CollisionResult]:
    return [r for r in rows if r.f_rad >= TARGET_SHARE]


def main() -> None:
    print("== baseline: two 10 kg balls, k = 1, 67 km/s")
    show(write_csv(OUT_DIR / "baseline.csv", baseline()))

    print(f"\n== k = {MASS_RATIO:g}, condensed balls, radius sweep (opacity x1 shown)")
    rows = write_csv(OUT_DIR / "sweep_k3.csv", sweep_k3())
    show(r for r in rows if r.kappa_scale == 1.0)
    hits = reaching(rows)
    best = max(rows, key=lambda r: r.f_rad)
    print(f"  cases at >= {TARGET_SHARE:.0%} of the heat as light: {len(hits)} of {len(rows)}")
    print(f"  best: {best.material} r = {best.impactor_radius:.0e} m, f_rad = {best.f_rad:.2e}")

    print(f"\n== k = {MASS_RATIO:g}, porous start (foam / spray), ideal Saha, opacity x1")
    rows = write_csv(OUT_DIR / "porous_k3.csv", porous_k3())
    show(rows)
    hits = reaching(rows)
    print(f"  cases at >= {TARGET_SHARE:.0%} of the heat as light: {len(hits)} of {len(rows)}")
    for r in hits:
        print(
            f"    {r.material} rho0 = {r.rho0:g} kg/m^3, r = {r.impactor_radius:g} m, "
            f"m = {r.m_impactor:.1e} kg: f_rad = {r.f_rad:.2f}, above 100 eV {r.f_xray:.1e}"
        )

    print("\n== near Sun: 1 kg carbon on 1 kg carbon at 618 km/s (the parent's pass-through case)")
    print("  (>100 eV and kT_emit: photosphere-zone bounds)")
    show_near_sun(write_csv(OUT_DIR / "near_sun.csv", near_sun()))


if __name__ == "__main__":
    main()
