"""A dense projectile stopped in a cheap polyethylene charge: the point fireball onto the plate.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.charge

Result 10 of docs/xray_collision_study.md. A graphite disk at 400 km/s threads the bore and stops in
a polyethylene charge of k times its mass, astern of the plate. The fireball is dense and opaque,
so it radiates almost nothing (the one-zone model on the CH2 table): its heat becomes expansion.
It drifts away from the plate at w = v/(1+k).

The plate impulse comes from the 2-D runs (`cargo run --release -p sweep -- --charge-2d`), which
are lossless effective-gamma Euler. Following the f(v) study's 2-D/1-D split, their gas-dynamic
rebound `1 + e_gamma` (the confined plane wave at the same gamma) is swapped for the polyethylene
table's radiative `1 + e_eff` at the speed and column the plate actually sees. The open bore's
share is subtracted. The mission is the craft falling from Jupiter to a 9.86 solar radii
periapsis, boosted there to cross Earth's orbit at 100 km/s.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat.xray_collision import fireball as fb
from puffsat.xray_collision.eos import EV, K_B, EosTable
from puffsat.xray_collision.materials import CARBON

G0 = 9.80665
V_REL = 400.0e3  # [m/s]
OUT_DIR = Path("data/results/xray_collision")
CH2_TABLE = Path("data/tables/ch2.json")
CHARGE_2D = OUT_DIR / "charge_2d.jsonl"
REBOUND_CSV = OUT_DIR / "propellant_rebound.csv"
PLATE_R = 10.0  # [m], the parent's 20 m bowl
GM_SUN = 1.32712440018e20
R_SUN = 6.957e8
AU = 1.495978707e11
PERIAPSIS = 9.86 * R_SUN  # Parker Solar Probe's closest approach
APHELION = 5.2 * AU  # falling from Jupiter
EARTH_CROSSING = 100.0e3  # [m/s]
SURVIVAL_AREAL = 30.0  # [kg/m^2], the parent's "tens of kilograms per square metre"


def boost_dv(aphelion: float = APHELION, periapsis: float = PERIAPSIS) -> float:
    """Impulsive periapsis boost turning a fall from `aphelion` into a 100 km/s Earth crossing."""
    v_p = math.sqrt(2.0 * GM_SUN * aphelion / (periapsis * (aphelion + periapsis)))
    need = math.sqrt(EARTH_CROSSING**2 + 2.0 * GM_SUN * (1.0 / periapsis - 1.0 / AU))
    return need - v_p


def ch2_table() -> EosTable:
    """The polyethylene table (`make xray-collision-ch2-table`) as a fireball `EosTable`."""
    d = json.loads(CH2_TABLE.read_text())
    n = d["shape"]
    f = {k: np.array(v).reshape(n) for k, v in d["fields"].items()}
    mat = dataclasses.replace(CARBON, name="ch2", mass_amu=14.027 / 3.0)
    return EosTable(
        mat,
        np.log(d["rho_grid"]),
        np.log(d["T_grid"]),
        np.log(f["p"]),
        f["e"],
        np.zeros(n),
        np.log(f["kappa_rosseland"]),
        np.log(f["kappa_planck"]),
    )


@dataclass(frozen=True)
class ChargeFireball:
    k: float
    heat: float  # [J per kg of projectile]
    kt_start_ev: float
    tau0: float
    f_rad: float
    u_rms: float  # [m/s], from the heat left after radiation
    drift: float  # [m/s], w = v/(1+k), away from the plate


def charge_fireball(
    k: float, charge_rho0: float = 950.0, table: EosTable | None = None
) -> ChargeFireball:
    """One-zone fireball of a 1 kg projectile in k kg of near-solid polyethylene."""
    tab = table or ch2_table()
    m, big = 1.0, k
    q = 0.5 * m * big / (m + big) * V_REL**2
    shock = fb.shocked_state(tab, math.sqrt(2.0 * q / big), charge_rho0)
    fire = fb._expand(tab, big, shock.rho, q, 1.0)
    rad = fire.y[3]
    return ChargeFireball(
        k=k,
        heat=q,
        kt_start_ev=fire.temp0 * K_B / EV,
        tau0=fire.tau0,
        f_rad=rad / q,
        u_rms=math.sqrt(2.0 * (q - rad) / (m + big)),
        drift=V_REL / (1.0 + k),
    )


def rebound(v: float, column: float) -> float:
    """Polyethylene `e_eff` at plate speed `v` and arriving column [kg/m^2] (L = 1 m rows of the
    `--propellant-rebound` sweep, interpolated in v and log column, held at the grid's edges)."""
    rows = [
        r
        for r in csv.DictReader(REBOUND_CSV.open())
        if float(r["opacity_scale"]) == 1.0
        and float(r["length"]) == 1.0
        and float(r["converged"]) == 1.0
    ]
    vs = sorted({float(r["v"]) for r in rows})
    cols = sorted({float(r["column"]) for r in rows})
    grid = np.array(
        [
            [
                next(
                    float(r["e_eff"])
                    for r in rows
                    if float(r["v"]) == vv and float(r["column"]) == cc
                )
                for cc in cols
            ]
            for vv in vs
        ]
    )
    lc = np.log(cols)
    x = float(np.clip(math.log(column), lc[0], lc[-1]))
    per_v = [float(np.interp(x, lc, row)) for row in grid]
    return float(np.interp(v, vs, per_v))


def plate_area(d_over_d: float, skirt: float, r: float = PLATE_R) -> float:
    """Wetted area: the paraboloid z = a r^2 (or a flat disc) plus a straight skirt."""
    if d_over_d > 0.0:
        a = d_over_d * 2.0 * r / r**2
        bowl = math.pi / (6.0 * a**2) * (math.pow(1.0 + 4.0 * a**2 * r**2, 1.5) - 1.0)
    else:
        bowl = math.pi * r**2
    return bowl + 2.0 * math.pi * r * skirt


@dataclass(frozen=True)
class ChargeCase:
    case: str
    grid: str  # "x1" production, "x2" the twice-finer convergence check
    k: float
    projectile_mass: float  # [kg]
    capture: float  # (J - J_hole) / (M u_rms), lossless gamma gas
    hole_share: float  # J_hole / J
    e_eff: float  # the polyethylene rebound swapped in
    v_e: float  # [m/s], impulse per kg of charge
    isp: float  # [s]
    mass_ratio: float
    projectile_per_final: float  # [kg per kg of final mass]
    impulse_per_projectile: float  # [N s per kg]
    peak_pressure: float  # [Pa]
    plate_area: float  # [m^2]
    plate_floor_mass: float  # [kg], SURVIVAL_AREAL x area


def charge_cases(
    rows: Sequence[Mapping[str, object]],
    projectile_mass: float,
    fireballs: Mapping[float, ChargeFireball],
    dv: float,
) -> list[ChargeCase]:
    out = []
    for r in rows:
        k = float(str(r["k"]))
        fire = fireballs[k]
        mass = (1.0 + k) * projectile_mass
        wall, hole = float(str(r["wall_impulse"])), float(str(r["hole_impulse"]))
        column = mass / (2.0 * math.pi * PLATE_R**2)  # about half the fireball over the bowl
        e = rebound(fire.u_rms, column)
        j = (wall - hole) * (1.0 + e) / float(str(r["confined_ratio"])) * mass * fire.u_rms
        v_e = j / (k * projectile_mass)
        mr = math.exp(dv / v_e)
        area = plate_area(float(str(r["d_over_d"])), float(str(r["skirt"])))
        out.append(
            ChargeCase(
                case=str(r["case"]),
                grid=str(r.get("grid", "x1")),
                k=k,
                projectile_mass=projectile_mass,
                capture=wall - hole,
                hole_share=hole / wall,
                e_eff=e,
                v_e=v_e,
                isp=v_e / G0,
                mass_ratio=mr,
                projectile_per_final=(mr - 1.0) / k,
                impulse_per_projectile=j / projectile_mass,
                peak_pressure=float(str(r["peak_local_pressure"])) * mass * fire.u_rms**2,
                plate_area=area,
                plate_floor_mass=SURVIVAL_AREAL * area,
            )
        )
    return out


def _write(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"python: wrote {len(rows)} rows -> {path}")


def load_jsonl(path: Path) -> list[dict[str, object]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def main() -> None:
    dv = boost_dv()
    print(
        f"== boost at 9.86 R_sun from a Jupiter fall: dv {dv / 1e3:.2f} km/s; "
        f"mass ratio 3 needs Isp {dv / math.log(3.0) / G0:.0f} s"
    )
    tab = ch2_table()
    fireballs = {k: charge_fireball(k, table=tab) for k in (20.0, 30.0, 40.0)}
    for f in fireballs.values():
        print(
            f"  k {f.k:3.0f}: kT0 {f.kt_start_ev:5.1f} eV tau0 {f.tau0:.1e} f_rad {f.f_rad:.4f} "
            f"u_rms {f.u_rms / 1e3:5.1f} km/s drift {f.drift / 1e3:4.1f} km/s "
            f"(w/u_rms {f.drift / f.u_rms:.3f})"
        )
    paths = sorted(OUT_DIR.glob("charge_2d*.jsonl"))
    rows = [r | {"grid": "x2" if "_x2" in p.name else "x1"} for p in paths for r in load_jsonl(p)]
    if not rows:
        print("(no charge_2d*.jsonl: run `cargo run --release -p sweep -- --charge-2d`)")
        return
    table_rows = []
    for m_p in (1.0, 10.0):
        print(f"\n== projectile {m_p:g} kg")
        for c in charge_cases(rows, m_p, fireballs, dv):
            table_rows.append(dataclasses.asdict(c))
            print(
                f"  {c.case:<40} {c.grid} k {c.k:2.0f} capture {c.capture:.3f} "
                f"hole {c.hole_share:.3f} "
                f"e {c.e_eff:.2f} Isp {c.isp:5.0f} s MR {c.mass_ratio:4.2f} "
                f"proj/final {c.projectile_per_final:.3f} peak {c.peak_pressure / 1e6:6.1f} MPa "
                f"area {c.plate_area:5.0f} m^2"
            )
    _write(OUT_DIR / "charge_ledger.csv", table_rows)


if __name__ == "__main__":
    main()
