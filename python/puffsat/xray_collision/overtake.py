"""Result 11: the Earth-departure overtake pusher -- what charge ratio k gives the ship 11 km/s?

    PYTHONPATH=python uv run python -m puffsat.xray_collision.overtake

An impactor arriving at 100 km/s from behind meets a polyethylene charge held a few metres off the
plate. Unlike Result 10 there is no bore: the merged mass drifts *onto* the plate at `v/(1+k)`, so
the impactor's own momentum is delivered as well as the expansion's. The plate impulse is the
2-D effective-gamma run (`sweep --overtake-2d`), with three geometries:

- the point fireball (drift toward the plate), on the parent bowl, the deep paraboloid at 1.5 F
  and a flat plate with a 4 m skirt;
- the real compact collision on the parent bowl, as its check;
- a planar collision: a wide impactor meeting a wide, thin polyethylene sheet in front of the flat
  skirted plate.

As in Result 10, polyethylene's radiative rebound is swapped in for the lossless gamma one,
`(1 + e_eff)/(1 + e_gamma)`, at the gas's arrival speed `u_rms + w`. The bond toll comes from
`freeze.freeze` at 100 km/s. Unlike the 400 km/s fireball, this one cools below 2500 K before the
plate and re-forms H2 with decades of margin, so only carbon's share is in doubt: it is bracketed
between condensed (no toll) and atomic.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache

import numpy as np

from puffsat.xray_collision import charge, freeze
from puffsat.xray_collision.eos import EosTable

V_CLOSE = 100.0e3  # [m/s], the impactor's speed relative to the departing ship
DV = 11.0e3  # [m/s], the departure boost to a high apoapsis
TARGET_MR = 1.5  # propellant 1/3 of the initial mass
PROJECTILE_COST = 1000.0  # [$/kg]
FACESHEET = 400.0e6  # [Pa]
JSONL = charge.OUT_DIR / "overtake_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "overtake_2d_coarse.jsonl"
#: The confinement follow-up (`sweep --overtake-2d confine`): deeper bowls, skirts, and the planar
#: sheet inside the δ/D 0.5 paraboloid.
JSONL_CONFINE = charge.OUT_DIR / "overtake_2d_confine.jsonl"
JSONL_CONFINE_COARSE = charge.OUT_DIR / "overtake_2d_confine_coarse.jsonl"
OUT = charge.OUT_DIR / "overtake_ledger.csv"
OUT_SHEET = charge.OUT_DIR / "overtake_sheet.csv"


@cache
def _fire(k: float) -> charge.ChargeFireball:
    return charge.charge_fireball(k, table=_table(), v_rel=V_CLOSE)


@cache
def _freeze(k: float) -> freeze.FreezeResult:
    return freeze.freeze(k, _table(), v_rel=V_CLOSE)


@cache
def _table() -> EosTable:
    return charge.ch2_table()


@dataclass(frozen=True)
class OvertakeCase:
    case: str
    kind: str
    grid: str
    k: float
    sheet_radius: float
    capture: float  # J / (M u_rms), lossless gamma gas (the point fireball's yardstick)
    impulse_per_mv: float  # J / (m v), lossless
    ceiling_share: float  # of the perfect plate's m v (1 + sqrt(1+k)), lossless
    eta: float  # (J/mv - 1) / sqrt(1+k): the share of the expansion's ideal momentum kept
    arrival_speed: float  # [m/s], u_rms + w
    e_eff: float
    v_e: float  # [m/s], per kg of charge, bonds returned
    v_e_carbon_frozen: float  # [m/s], H2 re-formed, carbon left atomic
    mass_ratio: float
    mass_ratio_carbon_frozen: float
    impactor_per_final: float  # [kg per kg of final mass], at the carbon-frozen mass ratio
    cost_per_kg_propellant: float  # [$], projectile cost over k
    max_impactor: float  # [kg per pulse] for the facesheet on a 20 m plate
    max_areal: float  # [kg/m^2] of impactor sheet at the facesheet limit (sheets only, else 0)
    exit_direction: float
    exit_uniformity: float


def cases(rows: list[dict[str, object]], grid: str) -> list[OvertakeCase]:
    out = []
    for r in rows:
        k = float(str(r["k"]))
        fire = _fire(k)
        loss = _freeze(k).velocity_loss_bonds_carbon_frozen
        per_mv = float(str(r["impulse_per_mv"]))
        peak = float(str(r["peak_per_mv2"])) * V_CLOSE**2  # [Pa per kg of impactor]
        m_max = FACESHEET / peak
        mass = (1.0 + k) * min(m_max, 20.0)
        column = mass / (2.0 * math.pi * charge.PLATE_R**2)
        v_arr = fire.u_rms + fire.drift
        e = charge.rebound(v_arr, column)
        j = per_mv * V_CLOSE * (1.0 + e) / float(str(r["confined_ratio"]))
        v_e = j / k
        v_e_c = v_e * (1.0 - loss)
        mr = math.exp(DV / v_e)
        mr_c = math.exp(DV / v_e_c)
        out.append(
            OvertakeCase(
                case=str(r["case"]),
                kind=str(r["kind"]),
                grid=grid,
                k=k,
                sheet_radius=float(str(r["sheet_radius"])),
                capture=float(str(r["impulse_per_mu"])),
                ceiling_share=per_mv / (1.0 + math.sqrt(1.0 + k)),
                eta=(per_mv - 1.0) / math.sqrt(1.0 + k),
                impulse_per_mv=per_mv,
                arrival_speed=v_arr,
                e_eff=e,
                v_e=v_e,
                v_e_carbon_frozen=v_e_c,
                mass_ratio=mr,
                mass_ratio_carbon_frozen=mr_c,
                impactor_per_final=(mr_c - 1.0) / k,
                cost_per_kg_propellant=PROJECTILE_COST / k,
                max_impactor=m_max,
                max_areal=(
                    m_max / (math.pi * float(str(r["sheet_radius"])) ** 2)
                    if r["kind"] == "sheet"
                    else 0.0
                ),
                exit_direction=float(str(r["exit_direction"])),
                exit_uniformity=float(str(r["exit_uniformity"])),
            )
        )
    return out


def k_for(rows: list[OvertakeCase], target: float = TARGET_MR, frozen: bool = True) -> float:
    """The k at which a plate's point-fireball rows reach `target`: a power-law fit of `v_e`
    against k (log-linear), solved for `v_e = DV / ln target`."""
    ks = np.log([c.k for c in rows])
    ve = np.log([c.v_e_carbon_frozen if frozen else c.v_e for c in rows])
    slope, icept = np.polyfit(ks, ve, 1)
    need = math.log(DV / math.log(target))
    return float(math.exp((need - icept) / slope))


#: Impactor-sheet materials: density [kg/m^3]. Polyimide (Kapton) is the choice (Seth,
#: 2026-10-15); carbon fibre is a fibre-plus-matrix laminate; Kevlar needs a UV coating.
SHEET_MATERIALS = {"polyimide (Kapton)": 1420.0, "carbon fibre laminate": 1600.0, "Kevlar": 1440.0}
SIGMA = 5.670374419e-8
SOLAR_CONSTANT = 1361.0  # [W/m^2] at 1 AU
C_LIGHT = 2.998e8
#: The sheet coasts from 0.5 AU to 1 AU. Its outward speed sets the coast time; 100 km/s is the
#: order of the Earth-crossing speed. It is a stated assumption, not a trajectory.
COAST = (0.5, 1.0)  # [AU]
COAST_SPEED = 100.0e3  # [m/s]
PULSE_TIME = 100.0e-6  # [s], the plate pulse's order: metres at 20-45 km/s


@dataclass(frozen=True)
class SheetBudget:
    areal: float  # [kg/m^2], the impactor sheet
    thickness: dict[str, float]  # [m] per material
    temp_half_au: float  # [K], black, face-on, radiating from both faces
    temp_one_au: float
    sail_accel_half_au: tuple[float, float]  # [m/s^2], absorbing / perfectly reflecting
    coast_time: float  # [s]
    coast_drift: tuple[float, float]  # [m], 1/2 <a> t^2 over the coast, absorbing / reflecting
    spread_sag: float  # [s], arrival spread from a 10 cm sag
    spread_tilt: float  # [s], from a 1 degree tilt across a 7 m radius


def sheet_budget(
    areal: float, sag: float = 0.1, tilt_deg: float = 1.0, radius: float = 7.0
) -> SheetBudget:
    """What a sheet of `areal` kg/m^2 must survive and how well it can be aimed, on its coast
    from 0.5 to 1 AU. The temperature is a black sheet facing the Sun that radiates from both
    faces; the sail push is `(1 + reflectivity) S / (c areal)`, averaged over `1/r^2` along the
    coast (`<1/r^2> = 1 / (r0 r1)`); the arrival spread is `sag / v` and `2 R tan(tilt) / v`."""
    r0, r1 = COAST
    s_half = SOLAR_CONSTANT / r0**2

    def temp(flux: float) -> float:
        return float((flux / (2.0 * SIGMA)) ** 0.25)

    a_half = s_half / (C_LIGHT * areal)
    mean = SOLAR_CONSTANT / (r0 * r1) / (C_LIGHT * areal)
    t = (r1 - r0) * charge.AU / COAST_SPEED
    return SheetBudget(
        areal=areal,
        thickness={name: areal / rho for name, rho in SHEET_MATERIALS.items()},
        temp_half_au=temp(s_half),
        temp_one_au=temp(SOLAR_CONSTANT),
        sail_accel_half_au=(a_half, 2.0 * a_half),
        coast_time=t,
        coast_drift=(0.5 * mean * t * t, mean * t * t),
        spread_sag=sag / V_CLOSE,
        spread_tilt=2.0 * radius * math.tan(math.radians(tilt_deg)) / V_CLOSE,
    )


def main() -> None:
    out: list[OvertakeCase] = []
    for path, grid in (
        (JSONL, "x1"),
        (JSONL_COARSE, "x0.5"),
        (JSONL_CONFINE, "x1"),
        (JSONL_CONFINE_COARSE, "x0.5"),
    ):
        if path.exists():
            out += cases(charge.load_jsonl(path), grid)
    print(
        f"  departure: {DV / 1e3:.0f} km/s from a {V_CLOSE / 1e3:.0f} km/s overtake; "
        f"mass ratio {TARGET_MR} needs v_e >= {DV / math.log(TARGET_MR) / 1e3:.1f} km/s"
    )
    for c in out:
        print(
            f"  {c.grid:>4} {c.case:<28} {c.kind:<7} R {c.sheet_radius:3.0f} k {c.k:3.0f}"
            f" | J/mv {c.impulse_per_mv:5.2f} ({c.ceiling_share:.0%} of ceiling) e {c.e_eff:.2f}"
            f" v_e {c.v_e / 1e3:5.1f}-{c.v_e_carbon_frozen / 1e3:5.1f} km/s"
            f" MR {c.mass_ratio:.3f}-{c.mass_ratio_carbon_frozen:.3f}"
            f" | {c.impactor_per_final:.3f} kg/kg final, ${c.cost_per_kg_propellant:.0f}/kg"
            f" | <= {c.max_impactor:5.1f} kg/pulse | dir {c.exit_direction:.2f}"
        )
    groups = sorted({(c.case, c.kind, c.sheet_radius) for c in out if c.kind != "compact"})
    for plate, kind, radius in groups:
        rows = [
            c
            for c in out
            if (c.case, c.kind, c.sheet_radius) == (plate, kind, radius) and c.grid == "x1"
        ]
        if len(rows) >= 2:
            label = plate if kind == "point" else f"{plate}, {radius:.1f} m sheet"
            print(
                f"  {label:<40}: mass ratio {TARGET_MR} at k = {k_for(rows, frozen=False):4.1f}"
                f" (bonds back) to {k_for(rows):4.1f} (carbon atomic)"
            )
    sheets = [c for c in out if c.kind == "sheet" and c.grid == "x1"]
    if sheets:
        lo, hi = 0.05, max(c.max_areal for c in sheets)
        print(
            f"  sheet limit at {FACESHEET / 1e6:.0f} MPa: {min(c.max_areal for c in sheets):.3f}-"
            f"{hi:.3f} kg/m^2 (production grid, an upper value; {lo} kg/m^2 the safer figure)"
        )
        budgets = [sheet_budget(a) for a in (lo, hi)]
        for b in budgets:
            thick = ", ".join(f"{n} {1e6 * t:.0f} um" for n, t in b.thickness.items())
            print(
                f"    {b.areal:.3f} kg/m^2: {thick} | {b.temp_half_au:.0f} K at 0.5 AU,"
                f" {b.temp_one_au:.0f} K at 1 AU | sail {1e3 * b.sail_accel_half_au[0]:.2f}-"
                f"{1e3 * b.sail_accel_half_au[1]:.2f} mm/s^2 at 0.5 AU, drift "
                f"{b.coast_drift[0] / 1e3:.0f}-{b.coast_drift[1] / 1e3:.0f} km over "
                f"{b.coast_time / 86400:.1f} d | arrival spread {1e6 * b.spread_sag:.1f} us (sag),"
                f" {1e6 * b.spread_tilt:.1f} us (tilt) vs ~{1e6 * PULSE_TIME:.0f} us"
            )
        with OUT_SHEET.open("w", newline="") as fh:
            sheet_rows: list[dict[str, object]] = []
            for b in budgets:
                row: dict[str, object] = {
                    key: val for key, val in asdict(b).items() if key != "thickness"
                }
                for name, t in b.thickness.items():
                    row[f"thickness_{name}"] = t
                sheet_rows.append(row)
            w = csv.DictWriter(fh, fieldnames=list(sheet_rows[0]))
            w.writeheader()
            w.writerows(sheet_rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(out[0])))
        w.writeheader()
        w.writerows(asdict(c) for c in out)
    print(f"python: wrote {len(out)} rows -> {OUT}")


if __name__ == "__main__":
    main()
