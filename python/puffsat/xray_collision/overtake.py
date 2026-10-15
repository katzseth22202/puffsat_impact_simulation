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
OUT = charge.OUT_DIR / "overtake_ledger.csv"


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
    arrival_speed: float  # [m/s], u_rms + w
    e_eff: float
    v_e: float  # [m/s], per kg of charge, bonds returned
    v_e_carbon_frozen: float  # [m/s], H2 re-formed, carbon left atomic
    mass_ratio: float
    mass_ratio_carbon_frozen: float
    impactor_per_final: float  # [kg per kg of final mass], at the carbon-frozen mass ratio
    cost_per_kg_propellant: float  # [$], projectile cost over k
    max_impactor: float  # [kg per pulse] for the facesheet on a 20 m plate
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


def main() -> None:
    out: list[OvertakeCase] = []
    for path, grid in ((JSONL, "x1"), (JSONL_COARSE, "x0.5")):
        if path.exists():
            out += cases(charge.load_jsonl(path), grid)
    print(
        f"  departure: {DV / 1e3:.0f} km/s from a {V_CLOSE / 1e3:.0f} km/s overtake; "
        f"mass ratio {TARGET_MR} needs v_e >= {DV / math.log(TARGET_MR) / 1e3:.1f} km/s"
    )
    for c in out:
        print(
            f"  {c.grid:>4} {c.case:<28} {c.kind:<7} R {c.sheet_radius:3.0f} k {c.k:3.0f}"
            f" | J/Mu {c.capture:.3f} J/mv {c.impulse_per_mv:5.2f} e {c.e_eff:.2f}"
            f" v_e {c.v_e / 1e3:5.1f}-{c.v_e_carbon_frozen / 1e3:5.1f} km/s"
            f" MR {c.mass_ratio:.3f}-{c.mass_ratio_carbon_frozen:.3f}"
            f" | {c.impactor_per_final:.3f} kg/kg final, ${c.cost_per_kg_propellant:.0f}/kg"
            f" | <= {c.max_impactor:5.1f} kg/pulse | dir {c.exit_direction:.2f}"
        )
    for plate in sorted({c.case for c in out if c.kind == "point"}):
        rows = [c for c in out if c.case == plate and c.kind == "point" and c.grid == "x1"]
        if len(rows) >= 2:
            print(
                f"  {plate:<28}: mass ratio {TARGET_MR} at k = {k_for(rows, frozen=False):4.1f}"
                f" (bonds back) to {k_for(rows):4.1f} (carbon atomic)"
            )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(out[0])))
        w.writeheader()
        w.writerows(asdict(c) for c in out)
    print(f"python: wrote {len(out)} rows -> {OUT}")


if __name__ == "__main__":
    main()
