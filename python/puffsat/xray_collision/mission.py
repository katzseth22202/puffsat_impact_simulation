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

- **The breeding cycle** (`mission_breeding.csv`). The Sun wave's expensive impactor only
  redirects a Sun ship, and that ship's *whole* remaining mass becomes the next Earth wave's
  100 km/s impactors. Per kg of Earth impactor, the departure (Result 11: mass ratio 1.5 at
  `k_E`, from `overtake_ledger.csv`) sends off `D = k_E/(MR_E - 1)` kg. A Sun ship of departing mass
  `S` returns `S (1 - f_dead)/MR_S` kg of impactor, so one kg back needs `S = MR_S/(1 - f_dead)`.
  - The growth per generation is `G = D/S`.
  - The share of departures left for payload is `1 - 1/G`.
  - The expensive Sun-wave impactor per kg of payload departure is
    `S (1 - 1/MR_S)/k_S / (D - S)`.
  - The mass that must already be in Earth orbit per kg of payload departure is
    `D MR_E / (D - S)`.

  `f_dead` is the Sun ship's mass that cannot be released as impactor. It is bracketed between 0
  and its plate: the δ/D 0.5 paraboloid's floor at 30 kg/m^2 on a 100 t final mass.

The decisions they back: the 400 km/s stage at mass ratio 2 (Seth, 2026-10-15), at $1000/kg of
projectile; Parker's 9.86 R_sun rather than 4 R_sun as the baseline periapsis; and the Sun wave as
the redirector of mass that becomes the Earth wave's impactors.
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
                k = float(r["k"])
                # The real collision, not the point fireball (charge.mixing_factor).
                seen.setdefault(
                    k, float(r["isp_bonds_frozen"]) * charge.G0 * charge.mixing_factor(plate, k)
                )
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


#: Result 11's recommended departure plate, and the Sun-stage plate it shares.
EARTH_PLATE = "paraboloid 0.50 at 1.5 F"
SUN_PLATE = "paraboloid 0.50 at 1.5 focus"
MR_EARTH = 1.5
MR_SUN = 2.0
FINAL_MASS = 100.0e3  # [kg], the handoffs' nominal final mass, for the plate's dead share
OUT_BREEDING = charge.OUT_DIR / "mission_breeding.csv"


def earth_k(frozen: bool) -> float:
    """Result 11's k for `MR_EARTH` on `EARTH_PLATE` from the committed `overtake_ledger.csv`
    (point fireball, production grid): carbon left atomic if `frozen`, else bonds returned."""
    rows = []
    with (charge.OUT_DIR / "overtake_ledger.csv").open() as fh:
        ledger = list(csv.DictReader(fh))
    for r in ledger:
        if r["case"] == EARTH_PLATE and r["kind"] == "point" and r["grid"] == "x1":
            k = float(r["k"])
            v = float(r["v_e_carbon_frozen" if frozen else "v_e"])
            rows.append((k, v * overtake_mixing(ledger, k)))
    slope, icept = np.polyfit(np.log([k for k, _ in rows]), np.log([v for _, v in rows]), 1)
    need = DV_EARTH / math.log(MR_EARTH)
    return float(math.exp((math.log(need) - icept) / slope))


DV_EARTH = 11.0e3  # [m/s], Result 11's departure


def overtake_mixing(ledger: list[dict[str, str]], k: float) -> float:
    """Result 11's real compact collision over its point fireball on the parent bowl (the only
    plate with both, k = 10 and 20), interpolated in ln k and held outside: the overtake's
    counterpart of `charge.mixing_factor`."""
    plate = "parent bowl at focus"
    x1 = [r for r in ledger if r["case"] == plate and r["grid"] == "x1"]
    pts = []
    for r in x1:
        if r["kind"] == "compact":
            point = next(p for p in x1 if p["kind"] == "point" and p["k"] == r["k"])
            pts.append((float(r["k"]), float(r["impulse_per_mv"]) / float(point["impulse_per_mv"])))
    pts.sort()
    return float(np.interp(math.log(k), [math.log(a) for a, _ in pts], [b for _, b in pts]))


@dataclass(frozen=True)
class Breeding:
    bonds: str  # "carbon atomic" or "bonds returned" (the Earth stage's bracket)
    dead_fraction: float
    k_earth: float
    k_sun: float
    departing_per_impactor: float  # D [kg per kg of Earth impactor]
    sun_ship_per_impactor: float  # S [kg departing per kg of impactor it returns]
    growth: float  # G = D / S
    payload_share: float  # 1 - 1/G of each departure
    sun_impactor_per_payload: float  # [kg per kg of payload departure]
    sun_impactor_cost_per_payload: float  # [$ per kg of payload departure]
    orbit_mass_per_payload: float  # [kg in Earth orbit per kg of payload departure]


def breeding_cycle(frozen: bool, dead_fraction: float) -> Breeding:
    k_e = earth_k(frozen)
    k_s = k_for_mass_ratio(SUN_PLATE, MR_SUN).k
    d = k_e / (MR_EARTH - 1.0)
    s = MR_SUN / (1.0 - dead_fraction)
    sun = s * (1.0 - 1.0 / MR_SUN) / k_s / (d - s)
    return Breeding(
        bonds="carbon atomic" if frozen else "bonds returned",
        dead_fraction=dead_fraction,
        k_earth=k_e,
        k_sun=k_s,
        departing_per_impactor=d,
        sun_ship_per_impactor=s,
        growth=d / s,
        payload_share=1.0 - s / d,
        sun_impactor_per_payload=sun,
        sun_impactor_cost_per_payload=PROJECTILE_COST * sun,
        orbit_mass_per_payload=d * MR_EARTH / (d - s),
    )


def plate_dead_fraction() -> float:
    """The Sun ship's δ/D 0.5 paraboloid floor at `charge.SURVIVAL_AREAL`, as a share of
    `FINAL_MASS`: the mass that cannot be released as impactor if the plate is kept whole."""
    return charge.SURVIVAL_AREAL * charge.plate_area(0.50, 0.0) / FINAL_MASS


OUT_REAL = charge.OUT_DIR / "mission_real_collision.csv"


def real_collision_ledger() -> list[dict[str, object]]:
    """Result 10's ledger on the real-collision basis: per plate and simulated k, the point
    fireball's `v_e` (bonds charged), the mixing factor, and the mass ratio and projectile bill
    for the 100 km/s crossing with the 10 h arc."""
    dv = charge.boost_dv()
    out: list[dict[str, object]] = []
    for plate in PLATES:
        ks, ve = _ledger(plate)
        for k, v in zip(ks, ve, strict=True):
            f = charge.mixing_factor(plate, float(k))
            v_point = float(v) / f
            mr = charge.finite_arc_mass_ratio(float(v), charge.ARC_HOURS * 3600.0)
            out.append(
                {
                    "plate": plate,
                    "k": float(k),
                    "v_e_point": v_point,
                    "mixing_factor": f,
                    "v_e_real": float(v),
                    "mass_ratio_point": math.exp(dv / v_point),
                    "mass_ratio_real": mr,
                    "projectile_per_final_real": (mr - 1.0) / float(k),
                }
            )
    return out


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
    breed = [
        breeding_cycle(frozen, dead)
        for frozen in (True, False)
        for dead in (0.0, plate_dead_fraction())
    ]
    for b in breed:
        print(
            f"  breeding ({b.bonds}, dead {b.dead_fraction:.2f}): k_E {b.k_earth:4.1f}, k_S"
            f" {b.k_sun:4.1f} | {b.departing_per_impactor:4.1f} kg departs per kg impactor,"
            f" G {b.growth:4.1f}, payload {b.payload_share:.0%} | Sun impactor"
            f" {1e3 * b.sun_impactor_per_payload:.2f} g (${b.sun_impactor_cost_per_payload:.2f})"
            f" and {b.orbit_mass_per_payload:.2f} kg in orbit per kg of payload departing"
        )
    real = real_collision_ledger()
    for r in real:
        print(
            f"  real collision {str(r['plate'])[:30]:<30} k {float(str(r['k'])):5.0f}: mixing"
            f" {float(str(r['mixing_factor'])):.2f}, mass ratio"
            f" {float(str(r['mass_ratio_point'])):.2f}"
            f" (point) -> {float(str(r['mass_ratio_real'])):.2f} (real, 10 h arc),"
            f" {float(str(r['projectile_per_final_real'])):.3f} kg/kg final"
        )
    OUT_PERIAPSIS.parent.mkdir(parents=True, exist_ok=True)
    _write(OUT_BREEDING, [asdict(b) for b in breed])
    _write(OUT_REAL, real)
    _write(OUT_PERIAPSIS, [asdict(p) for p in per])
    _write(OUT_K, [asdict(c) for c in ks])
    print(
        f"python: wrote {len(per)} rows -> {OUT_PERIAPSIS}, {len(ks)} rows -> {OUT_K},"
        f" {len(breed)} rows -> {OUT_BREEDING}"
    )


if __name__ == "__main__":
    main()
