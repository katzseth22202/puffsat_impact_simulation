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
# Polyethylene's atomization energy, 2 C-H (4.3 eV) + 1 C-C (3.6 eV) per CH2: the toll if no bond
# re-forms in the expanding fireball (the CH2 table counts atoms only).
PE_BONDS = (2.0 * 4.3 + 3.6) * 1.602176634e-19 / (14.027 * 1.66053906660e-27)  # ~84 MJ/kg
ARC_HOURS = 10.0  # nominal burn arc for the finite-arc mass ratio


def boost_dv(
    aphelion: float = APHELION, periapsis: float = PERIAPSIS, crossing: float = EARTH_CROSSING
) -> float:
    """Impulsive periapsis boost turning a fall from `aphelion` into an Earth crossing at
    `crossing` (100 km/s by default)."""
    v_p = math.sqrt(2.0 * GM_SUN * aphelion / (periapsis * (aphelion + periapsis)))
    need = math.sqrt(crossing**2 + 2.0 * GM_SUN * (1.0 / periapsis - 1.0 / AU))
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
    drift: float  # [m/s], w = v/(1+k): away from the plate (Result 10), onto it (Result 11)


def charge_fireball(
    k: float,
    charge_rho0: float = 950.0,
    table: EosTable | None = None,
    v_rel: float = V_REL,
) -> ChargeFireball:
    """One-zone fireball of a 1 kg projectile closing at `v_rel` on k kg of near-solid
    polyethylene."""
    tab = table or ch2_table()
    m, big = 1.0, k
    q = 0.5 * m * big / (m + big) * v_rel**2
    # Shocked solid polyethylene (~4x 950 kg/m^3) lies above the CH2 table's 1284 kg/m^3 ceiling,
    # so the fireball starts at the ceiling, already a little expanded. It radiates nothing either
    # way (it is ~1e6 optical depths thick), so u_rms does not depend on this choice.
    ceiling = float(np.exp(tab.ln_rho[-1]))
    shock = fb.shocked_state(tab, math.sqrt(2.0 * q / big), min(charge_rho0, ceiling / 4.5))
    fire = fb._expand(tab, big, min(shock.rho, ceiling), q, 1.0)
    rad = fire.y[3]
    return ChargeFireball(
        k=k,
        heat=q,
        kt_start_ev=fire.temp0 * K_B / EV,
        tau0=fire.tau0,
        f_rad=rad / q,
        u_rms=math.sqrt(2.0 * (q - rad) / (m + big)),
        drift=v_rel / (1.0 + k),
    )


REBOUND_LOW_JSONL = OUT_DIR / "propellant_rebound_low.jsonl"
REBOUND_LOW_CSV = OUT_DIR / "propellant_rebound_low.csv"


def rebound_rows() -> list[dict[str, str]]:
    """The 60-150 km/s rebound rows plus the 20-50 km/s ones (`--propellant-rebound low`), whose
    CSV this writes from the Rust JSONL when it is newer."""
    from puffsat.xray_collision import anvil

    if REBOUND_LOW_JSONL.exists() and (
        not REBOUND_LOW_CSV.exists()
        or REBOUND_LOW_CSV.stat().st_mtime < REBOUND_LOW_JSONL.stat().st_mtime
    ):
        _write(REBOUND_LOW_CSV, anvil.load_rebound(REBOUND_LOW_JSONL))
    rows = list(csv.DictReader(REBOUND_CSV.open()))
    if REBOUND_LOW_CSV.exists():
        rows += list(csv.DictReader(REBOUND_LOW_CSV.open()))
    return rows


def rebound(v: float, column: float) -> float:
    """Polyethylene `e_eff` at plate speed `v` and arriving column [kg/m^2] (L = 1 m rows of the
    `--propellant-rebound` sweeps, interpolated in v and log column, held at the grid's edges)."""
    rows = [
        r
        for r in rebound_rows()
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
    isp_bonds_frozen: float  # [s], if no polyethylene bond re-forms
    mass_ratio_arc: float  # over a ARC_HOURS burn centred on periapsis
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
        heat_per_charge = fire.heat / (k + 1.0)  # per kg of the merged fireball
        v_e_bonds = v_e * math.sqrt(max(1.0 - PE_BONDS * k / (k + 1.0) / heat_per_charge, 0.0))
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
                isp_bonds_frozen=v_e_bonds / G0,
                mass_ratio_arc=finite_arc_mass_ratio(v_e, ARC_HOURS * 3600.0),
                projectile_per_final=(mr - 1.0) / k,
                impulse_per_projectile=j / projectile_mass,
                peak_pressure=float(str(r["peak_local_pressure"])) * mass * fire.u_rms**2,
                plate_area=area,
                plate_floor_mass=SURVIVAL_AREAL * area,
            )
        )
    return out


def _orbit_rhs(y: np.ndarray, accel: float) -> np.ndarray:
    """Heliocentric planar motion with a thrust acceleration `accel` along the velocity."""
    x, yy, vx, vy = y
    r3 = (x * x + yy * yy) ** 1.5
    v = math.hypot(vx, vy)
    return np.array([vx, vy, -GM_SUN * x / r3 + accel * vx / v, -GM_SUN * yy / r3 + accel * vy / v])


def _crossing_speed(v_e: float, burn: float, mass_ratio: float) -> float:
    """Speed at 1 AU after a burn of `burn` seconds centred on periapsis, constant mass flow,
    exhaust `v_e`, burning the ship down by `mass_ratio`, from the Jupiter fall."""
    v_p = math.sqrt(2.0 * GM_SUN * APHELION / (PERIAPSIS * (APHELION + PERIAPSIS)))
    if burn <= 0.0:
        v_after = v_p + v_e * math.log(mass_ratio)
        return math.sqrt(max(v_after**2 - 2.0 * GM_SUN * (1.0 / PERIAPSIS - 1.0 / AU), 0.0))
    n = 600
    dt = burn / 2.0 / n
    y = np.array([PERIAPSIS, 0.0, 0.0, v_p])
    for _ in range(n):  # coast back from periapsis to the burn's start (time-reversed RK4)
        k1 = _orbit_rhs(y, 0.0)
        k2 = _orbit_rhs(y - 0.5 * dt * k1, 0.0)
        k3 = _orbit_rhs(y - 0.5 * dt * k2, 0.0)
        k4 = _orbit_rhs(y - dt * k3, 0.0)
        y = y - dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    mdot = (1.0 - 1.0 / mass_ratio) / burn  # per unit initial mass
    m = 1.0
    for _ in range(2 * n):
        a1 = mdot * v_e / m
        a2 = mdot * v_e / (m - 0.5 * dt * mdot)
        a3 = mdot * v_e / (m - dt * mdot)
        k1 = _orbit_rhs(y, a1)
        k2 = _orbit_rhs(y + 0.5 * dt * k1, a2)
        k3 = _orbit_rhs(y + 0.5 * dt * k2, a2)
        k4 = _orbit_rhs(y + dt * k3, a3)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        m -= dt * mdot
    energy = 0.5 * (y[2] ** 2 + y[3] ** 2) - GM_SUN / math.hypot(y[0], y[1])
    return math.sqrt(max(2.0 * (energy + GM_SUN / AU), 0.0))


def finite_arc_mass_ratio(v_e: float, burn: float) -> float:
    """Mass ratio that crosses 1 AU at 100 km/s with a `burn`-second arc centred on periapsis."""
    impulsive = math.exp(boost_dv() / v_e)
    lo, hi = impulsive, 1.5 * impulsive
    for _ in range(18):
        mid = math.sqrt(lo * hi)
        if _crossing_speed(v_e, burn, mid) < EARTH_CROSSING:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


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
    paths = sorted(OUT_DIR.glob("charge_2d*.jsonl"))
    rows = [r | {"grid": "x2" if "_x2" in p.name else "x1"} for p in paths for r in load_jsonl(p)]
    ks = sorted({float(str(r["k"])) for r in rows} | {20.0, 30.0, 40.0})
    fireballs = {k: charge_fireball(k, table=tab) for k in ks}
    for f in fireballs.values():
        print(
            f"  k {f.k:3.0f}: kT0 {f.kt_start_ev:5.1f} eV tau0 {f.tau0:.1e} f_rad {f.f_rad:.4f} "
            f"u_rms {f.u_rms / 1e3:5.1f} km/s drift {f.drift / 1e3:4.1f} km/s "
            f"(w/u_rms {f.drift / f.u_rms:.3f})"
        )
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
                f"e {c.e_eff:.2f} Isp {c.isp:5.0f} s ({c.isp_bonds_frozen:5.0f}) "
                f"MR {c.mass_ratio:4.2f} (arc {c.mass_ratio_arc:4.2f}) "
                f"proj/final {c.projectile_per_final:.3f} peak {c.peak_pressure / 1e6:6.1f} MPa "
                f"area {c.plate_area:5.0f} m^2"
            )
    _write(OUT_DIR / "charge_ledger.csv", table_rows)


if __name__ == "__main__":
    main()
