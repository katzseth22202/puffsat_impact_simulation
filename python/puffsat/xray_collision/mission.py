"""The charge pusher's mission arithmetic: periapsis choice, and the k that meets a mass ratio.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.mission

Two tables back decisions taken after Result 10, so that none rests on arithmetic done by hand:

- **Periapsis** (`mission_periapsis.csv`). For a fall from Jupiter to periapsis `r_p`, it gives the
  periapsis speed, the head-on closing speed with a mirror-image retrograde projectile, the
  impulsive boost to an Earth crossing at 100 or 150 km/s, and the sunlight on a shield there.
  `flat_plate_temp` is a black flat plate facing the Sun and radiating from that face only. It is
  an upper value: Parker Solar Probe's white-coated shield runs ~1650 K at 9.86 R_sun against
  ~1840 K here.
- **k for a target mass ratio** (`mission_k.csv`). For each plate of Result 10's ledger
  (`charge_ledger.csv`, 1 kg projectile, production grid), `v_e` with every polyethylene bond
  charged is fitted as a power of k, and solved for `v_e = dv / ln(MR / arc)`. `arc` is the
  finite-arc penalty, the 10 h arc's mass ratio over the impulsive one at the solved `v_e`
  (`charge.finite_arc_mass_ratio`, ~1%). The projectile bill is `(MR - 1)/k` per kg of final mass,
  and its cost per kg of propellant is `PROJECTILE_COST / k`.

The decisions they back: the 400 km/s stage at mass ratio 2 (Seth, 2026-10-15), at $1000/kg of
projectile, and Parker's 9.86 R_sun rather than 4 R_sun as the baseline periapsis.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.xray_collision import charge

PROJECTILE_COST = 1000.0  # [$/kg], Seth's mass-production basis
SIGMA = 5.670374419e-8
SOLAR_CONSTANT = 1361.0  # [W/m^2] at 1 AU
PERIAPSES_RSUN = (9.86, 8.0, 6.0, 4.0)
CROSSINGS = (100.0e3, 150.0e3)
TARGET_MRS = (1.5, 2.0, 3.0)
PLATES = (
    "paraboloid 0.30 at focus (parent bowl)",
    "bowl 0.30 + skirt 4 m at focus",
    "paraboloid 0.50 at 1.5 focus",
)
OUT_PERIAPSIS = charge.OUT_DIR / "mission_periapsis.csv"
OUT_K = charge.OUT_DIR / "mission_k.csv"


@dataclass(frozen=True)
class Periapsis:
    r_sun: float
    periapsis_speed: float  # [m/s]
    closing_speed: float  # [m/s], head-on with a mirror-image retrograde projectile
    dv_100: float  # [m/s] for a 100 km/s Earth crossing
    dv_150: float
    flux: float  # [W/m^2]
    suns: float
    flat_plate_temp: float  # [K], black, Sun-facing face only


def periapsis(r_sun: float) -> Periapsis:
    rp = r_sun * charge.R_SUN
    a, gm = charge.APHELION, charge.GM_SUN
    v_p = math.sqrt(2.0 * gm * a / (rp * (a + rp)))
    flux = SOLAR_CONSTANT * (charge.AU / rp) ** 2
    return Periapsis(
        r_sun=r_sun,
        periapsis_speed=v_p,
        closing_speed=2.0 * v_p,
        dv_100=charge.boost_dv(periapsis=rp, crossing=CROSSINGS[0]),
        dv_150=charge.boost_dv(periapsis=rp, crossing=CROSSINGS[1]),
        flux=flux,
        suns=flux / SOLAR_CONSTANT,
        flat_plate_temp=(flux / SIGMA) ** 0.25,
    )


@dataclass(frozen=True)
class KForMassRatio:
    plate: str
    crossing: float  # [m/s]
    mass_ratio: float
    v_e_needed: float  # [m/s], bonds charged, arc included
    k: float
    in_simulated_range: bool  # k inside the ledger's 20-100
    fit_exponent: float  # v_e ~ k^exponent
    projectile_per_final: float  # [kg per kg of final mass]
    cost_per_kg_propellant: float  # [$]


def _ledger(plate: str) -> tuple[np.ndarray, np.ndarray]:
    """`(k, v_e with bonds charged)` for one plate: production grid, 1 kg projectile."""
    seen: dict[float, float] = {}
    with (charge.OUT_DIR / "charge_ledger.csv").open() as fh:
        for r in csv.DictReader(fh):
            if (
                r["case"] == plate
                and r["grid"] == "x1"
                and abs(float(r["projectile_mass"]) - 1.0) < 1e-9
            ):
                seen.setdefault(float(r["k"]), float(r["isp_bonds_frozen"]) * charge.G0)
    ks = np.array(sorted(seen))
    return ks, np.array([seen[k] for k in ks])


def k_for_mass_ratio(
    plate: str, mass_ratio: float, crossing: float = CROSSINGS[0]
) -> KForMassRatio:
    ks, ve = _ledger(plate)
    slope, icept = np.polyfit(np.log(ks), np.log(ve), 1)
    dv = charge.boost_dv(crossing=crossing)
    need = dv / math.log(mass_ratio)
    # The finite arc costs ~1% of mass ratio; fold it in once at the impulsive solution.
    k0 = math.exp((math.log(need) - icept) / slope)
    v0 = math.exp(icept + slope * math.log(k0))
    arc = (
        charge.finite_arc_mass_ratio(v0, charge.ARC_HOURS * 3600.0) / math.exp(dv / v0)
        if crossing == CROSSINGS[0]
        else 1.0
    )
    need = dv / math.log(mass_ratio / arc)
    k = math.exp((math.log(need) - icept) / slope)
    return KForMassRatio(
        plate=plate,
        crossing=crossing,
        mass_ratio=mass_ratio,
        v_e_needed=need,
        k=k,
        in_simulated_range=bool(ks[0] <= k <= ks[-1]),
        fit_exponent=float(slope),
        projectile_per_final=(mass_ratio - 1.0) / k,
        cost_per_kg_propellant=PROJECTILE_COST / k,
    )


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    per = [periapsis(r) for r in PERIAPSES_RSUN]
    for p in per:
        print(
            f"  {p.r_sun:5.2f} R_sun: v_p {p.periapsis_speed / 1e3:4.0f} km/s, closing "
            f"{p.closing_speed / 1e3:4.0f} | dv to 100 km/s {p.dv_100 / 1e3:5.1f}, to 150 "
            f"{p.dv_150 / 1e3:5.1f} km/s | {p.suns:5.0f} suns,"
            f" flat plate {p.flat_plate_temp:5.0f} K"
        )
    ks = [
        k_for_mass_ratio(plate, mr, x) for x in CROSSINGS for plate in PLATES for mr in TARGET_MRS
    ]
    for c in ks:
        flag = "" if c.in_simulated_range else "  (extrapolated)"
        print(
            f"  {c.crossing / 1e3:3.0f} km/s {c.plate[:30]:<30} MR {c.mass_ratio:3.1f}: "
            f"v_e {c.v_e_needed / 1e3:5.1f} km/s -> k {c.k:6.1f}, "
            f"{c.projectile_per_final:.3f} kg/kg final, ${c.cost_per_kg_propellant:5.0f}/kg{flag}"
        )
    OUT_PERIAPSIS.parent.mkdir(parents=True, exist_ok=True)
    _write(OUT_PERIAPSIS, [asdict(p) for p in per])
    _write(OUT_K, [asdict(c) for c in ks])
    print(f"python: wrote {len(per)} rows -> {OUT_PERIAPSIS}, {len(ks)} rows -> {OUT_K}")


if __name__ == "__main__":
    main()
