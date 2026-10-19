"""Result 14: a second, staged plate for the 3S departure -- the impactor threads a hole.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.departure_hole

The parent's S38 (`docs/departure_occultation_ask_for_aim_repo.md`, commit 8413f3f) finds that the
flown mirror sign of the solar-dive departure puts periapsis behind Earth: the returning sheets
would pass through Earth. The visible sign meets the stream at **129.9 deg** to the thrust axis
(9.86 R_sun; 117.0 at 4 R_sun) instead of 89.2 (76.5). Its stream closes faster, at 84.0 against
75.0 km/s (107.9 against 98.8).

Above 90 deg the impactor comes from the ship's side, so Seth's design copies the Sun plate. The
push plate is dropped, and a second plate takes over. Compact impactors thread a narrow hole in a
bowl and stop in a compact polyethylene charge at its focus. The hole slides as the stream's
direction turns over the burn, and its rim takes ablative spray.

Five things are computed here:

1. **The ceiling.** A flat frictionless plate keeps at most `m w (cos theta + sqrt(1 + k -
   sin^2 theta))` (Result 12). A bowl, which can take tangential force, keeps at most
   `m w (cos theta + sqrt(1 + k))`. Both lose the impactor's own momentum share `cos theta`.
2. **The bowl's capture** (`sweep --hole-plate-2d`). Point fireballs at the focus drift along the
   axis at `-cos(theta)/sqrt(k)` of `u_rms` toward the plate, which is the stream's axial share.
   The sideways share is not run (the runs are axisymmetric). It carries `sin^2 theta / (1 + k)`
   of the energy, 2-6% at k = 10-30.
3. **v_e and Isp**, per kg of charge. Each run's lossless impulse takes polyethylene's rebound
   at the fireball's `u_rms` (as Results 10-12 do) and the real disk-into-charge collision's share
   of the point fireball (`charge.mixing_factor` of the Sun plate). The bond toll is Result 11's
   at the closing speed: H2 re-forms, and carbon is bracketed between atomic and condensed. The
   flown design at 89.2 deg is Result 12's canted flat plate (`growth.earth_ve`).
4. **The hole.** The impactor's line through the charge meets the bowl off-axis. The hole loses the
   impulse and the gas that land on its area there (the ring loads at that radius).
5. **The burn.** A planar two-body burn, centred on the 10-day orbit's 600 km periapsis, with
   constant thrust held in the one fixed direction that loses least. Its loss is the extra ideal
   burn needed to reach the impulsive burn's energy (the calcs repo's `finite_burn_loss` model,
   re-implemented here in SI units).
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from puffsat.xray_collision import charge, freeze, growth, mission

G0 = 9.80665
MU = mission.GM_EARTH
R_P = mission.R_EARTH + mission.PARKING_PERIAPSIS_ALT


@dataclass(frozen=True)
class Route:
    """A departure route from the parent's S38 table (its planar check of the calcs routes)."""

    name: str
    sign: str
    theta: float  # [deg], mean stream angle to the thrust axis
    w: float  # [m/s], mean closing speed
    dv: float  # [m/s], the impulsive burn


ROUTES = (
    Route("9.86 R_sun", "flown (behind Earth)", 89.2, 75.0e3, 4.21e3),
    Route("9.86 R_sun", "visible", 129.9, 84.0e3, 4.21e3),
    Route("4 R_sun", "flown (behind Earth)", 76.5, 98.8e3, 4.44e3),
    Route("4 R_sun", "visible", 117.0, 107.9e3, 4.44e3),
)
KS = (5.0, 10.0, 20.0, 30.0, 40.0, 60.0)
DESIGN_PLATE = growth.SUN_PLATE.replace("focus", "F")  # "paraboloid 0.50 at 1.5 F"
PLATES = (
    "parent bowl at focus",
    "paraboloid 0.50 at focus",
    "paraboloid 0.50 at 1.5 F",
    "paraboloid 0.75 at 2 F",
)
HOLE_RADII = (0.25, 0.5)  # [m]
JSONL = charge.OUT_DIR / "hole_plate_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "hole_plate_2d_coarse.jsonl"
#: Result 16's runs (`sweep --band-bowl-2d`): the bowl with its charge at 2.5 F, and the 1.5 F
#: bowl at 160-176 deg. Read alongside when present; Result 14's own rows use exact angles.
BAND_JSONL = charge.OUT_DIR / "band_bowl_2d.jsonl"
BAND_JSONL_COARSE = charge.OUT_DIR / "band_bowl_2d_coarse.jsonl"

# --- The departure stage (Seth, 2026-10-18): 1 Hz, ~0.25 g mean, the push plate dropped.
CADENCE = 1.0  # [Hz]
G_MEAN = 0.25
#: The stack at the departure: the 1500 t unit after the 2 Hz / 1 g push (mass ratio 1.330),
#: the apoapsis methalox (~0.24 km/s at 3.7 km/s; parent S38 item 3) and its dropped 12.2 t plate
#: (Result 13). An input, not a ledger: the calcs repo owns the unit's books.
PUSHED = 1500.0e3 / 1.330
METHALOX = (0.24e3, 3.7e3)  # [m/s] burn, exhaust speed
PUSH_PLATE = 12.2e3
STACK = PUSHED * math.exp(-METHALOX[0] / METHALOX[1]) - PUSH_PLATE
G_CAP = 100.0
STROKE = 2.0  # [m]
#: The calcs repo's paraboloid survival floor (solar_dive_growth._plate_for): 30 kg/m^2 over
#: 533 m^2, rated for 1 kg at 400 km/s, scaled as pulse energy^(2/3). A point fireball, so the
#: point-blast law does apply here.
SUN_PLATE_AREA = 533.041
SUN_PLATE_RATED = 0.5 * 1.0 * 400.0e3**2  # [J]
BORE_COLLAR = 1.0e3  # [kg], Result 10's ~0.5-1 t (upper)
BURN_TIMES = (300.0, 480.0, 600.0, 1000.0, 1600.0, 2400.0)  # [s]; 480 s is the flown 8 minutes

OUT_LAW = charge.OUT_DIR / "departure_hole_law.csv"
OUT_HOLE = charge.OUT_DIR / "departure_hole_loss.csv"
OUT_BURN = charge.OUT_DIR / "departure_hole_burn.csv"
OUT_STAGE = charge.OUT_DIR / "departure_hole_stage.csv"


def flat_ceiling(k: float, theta: float) -> float:
    t = math.radians(theta)
    return math.cos(t) + math.sqrt(1.0 + k - math.sin(t) ** 2)


def bowl_ceiling(k: float, theta: float) -> float:
    return math.cos(math.radians(theta)) + math.sqrt(1.0 + k)


@cache
def _runs(grid: str) -> dict[tuple[str, float, float], dict[str, object]]:
    paths = (JSONL, BAND_JSONL) if grid == "x1" else (JSONL_COARSE, BAND_JSONL_COARSE)
    return {
        (str(r["case"]), float(str(r["theta_deg"])), float(str(r["k"]))): r
        for path in paths
        if path.exists()
        for r in charge.load_jsonl(path)
    }


def _capture_at(plate: str, theta: float, k: float, grid: str) -> float:
    runs = _runs(grid)
    ks = sorted({kk for (p, t, kk) in runs if p == plate and t == theta})
    vals = [float(str(runs[(plate, theta, kk)]["impulse_per_mv"])) for kk in ks]
    return float(np.interp(math.log(k), np.log(ks), vals))


def capture(plate: str, theta: float, k: float, grid: str = "x1") -> float:
    """Plate impulse per `m w`, log-interpolated in k between the runs and linear in
    `cos theta` between their angles (the drift goes as `cos theta`); linear past the ends."""
    thetas = sorted({t for (p, t, _) in _runs(grid) if p == plate})
    if theta in thetas:
        return _capture_at(plate, theta, k, grid)
    c = math.cos(math.radians(theta))
    cs = [math.cos(math.radians(t)) for t in thetas]  # descending in theta -> ascending in cos
    order = sorted(range(len(thetas)), key=lambda i: cs[i])
    xs = [cs[i] for i in order]
    ys = [_capture_at(plate, thetas[i], k, grid) for i in order]
    i = int(np.clip(np.searchsorted(xs, c) - 1, 0, len(xs) - 2))
    return ys[i] + (ys[i + 1] - ys[i]) * (c - xs[i]) / (xs[i + 1] - xs[i])


def bond_loss(k: float, w: float) -> float:
    """Velocity share lost if carbon stays atomic (H2 re-formed): Result 11's `freeze` toll."""
    return freeze.freeze(k, charge.ch2_table(), v_rel=w).velocity_loss_bonds_carbon_frozen


def hole_ve(
    plate: str, theta: float, k: float, w: float, grid: str = "x1", frozen: bool = True
) -> float:
    """The hole design's `v_e` [m/s] per kg of charge, priced as Results 10-12 price a plate:
    the run's lossless `J/(m w)` with polyethylene's rebound swapped in at the fireball's `u_rms`
    (`(1 + e_eff)/confined_ratio`), reduced to the real disk-into-charge collision
    (`charge.mixing_factor` of the Sun plate, held at its k = 30 value below 30), and charged
    the carbon-atomic toll at `w` if `frozen`."""
    row = next(iter(_runs(grid).values()))
    confined = float(str(row["confined_ratio"]))
    u_rms = w * math.sqrt(k) / (1.0 + k)
    column = (1.0 + k) * growth.REBOUND_IMPACTOR / (2.0 * math.pi * charge.PLATE_R**2)
    e = charge.rebound(u_rms, column)
    j = capture(plate, theta, k, grid) * w * (1.0 + e) / confined
    v = j / k * charge.mixing_factor(growth.SUN_PLATE, k)
    return v * (1.0 - bond_loss(k, w)) if frozen else v


def flat_ve(theta: float, k: float, w: float, frozen: bool = True) -> float:
    """Result 12's canted flat plate (square charge) at `theta`, closing `w`; carbon atomic if
    `frozen`. Its share of the ceiling was fitted at 101 deg and is carried to `theta`."""
    fit = growth.earth_fit("square charge")
    c = growth.Closure("route", 0.0, 0.0, 0.0, 0.0, (theta, theta), (w, w), 0.0, 0.0)
    return growth.earth_ve(fit, c, k, frozen)[0]


def hole_geometry(plate: str, theta: float) -> tuple[float, float]:
    """`(r, z)` [m] where the impactor's line through the charge meets the bowl `z = r^2/4F`."""
    runs = _runs("x1")
    row = next(r for (p, t, _), r in runs.items() if p == plate and t == theta)
    f, h = float(str(row["focus"])), float(str(row["height"]))
    phi = math.radians(180.0 - theta)  # the line's angle from the axis, toward the vertex
    # z = h - t cos phi, r = t sin phi on z = r^2 / 4F.
    a, b, c = math.sin(phi) ** 2 / (4.0 * f), math.cos(phi), -h
    t = (-b + math.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
    return t * math.sin(phi), h - t * math.cos(phi)


def hole_loss(plate: str, theta: float, k: float, radius: float) -> tuple[float, float]:
    """`(impulse share, mass share)` that lands on a hole of `radius` at the hole's ring."""
    row = _runs("x1")[(plate, theta, k)]
    r_hole, _ = hole_geometry(plate, theta)
    bands = row["bands"]
    assert isinstance(bands, list)
    band = next(
        b for b in bands if isinstance(b, dict) and float(b["r_lo"]) <= r_hole < float(b["r_hi"])
    )
    area = math.pi * radius**2
    j = float(str(row["impulse_per_mv"]))
    return float(band["impulse"]) * area / j, float(band["mass_in"]) * area / (1.0 + k)


# --- The finite burn ------------------------------------------------------------------------


_RTOL, _ATOL = 1.0e-10, 1.0e-3


def _gravity(r: NDArray[np.float64]) -> NDArray[np.float64]:
    return -MU * r / float(np.linalg.norm(r)) ** 3


def burn_loss(
    dv: float, exhaust: float, seconds: float, period: float = mission.PARKING_PERIOD
) -> float:
    """Fixed-direction finite-burn loss [m/s]: constant mass flow, centred on periapsis, held in
    the one direction that loses least (the calcs repo's `finite_burn_loss`, SI units)."""
    if seconds <= 0.0:
        return 0.0
    a = (MU * (period / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
    v_p = math.sqrt(MU * (2.0 / R_P - 1.0 / a))
    target = 0.5 * (v_p + dv) ** 2 - MU / R_P

    def coast(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.concatenate((y[2:], _gravity(y[:2])))

    back = solve_ivp(coast, (0.0, -0.5 * seconds), [R_P, 0.0, 0.0, v_p], rtol=_RTOL, atol=_ATOL)
    start = np.asarray(back.y[:, -1], dtype=np.float64)

    def energy(ideal: float, phi: float) -> float:
        spent = 1.0 - math.exp(-ideal / exhaust)
        heading = np.array([math.cos(phi), math.sin(phi)])

        def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
            thrust = exhaust * spent / seconds / (1.0 - spent * t / seconds)
            return np.concatenate((y[2:], _gravity(y[:2]) + thrust * heading))

        end = solve_ivp(rhs, (0.0, seconds), start, rtol=_RTOL, atol=_ATOL).y[:, -1]
        return float(0.5 * end[2:] @ end[2:] - MU / np.linalg.norm(end[:2]))

    def needed(phi: float) -> float:
        return float(brentq(lambda x: energy(x, phi) - target, dv, 2.0 * dv + 1.0, xtol=1e-3))

    along = 0.5 * math.pi
    best = minimize_scalar(needed, bounds=(along - 0.5, along + 0.5), method="bounded")
    return float(best.fun) - dv


# --- The stage ------------------------------------------------------------------------------


@dataclass(frozen=True)
class Stage:
    route: str
    k: float
    v_e: float  # [m/s]
    isp: float  # [s]
    impulsive_dv: float
    loss: float  # [m/s]
    mass_ratio: float
    burn_s: float
    pulses: float
    pulse_impulse: float  # [N s], mean
    impactor_per_pulse: float  # [kg]
    plate_floor: float  # [kg], the point-fireball survival floor
    plate_inertia: float  # [kg], at G_CAP over STROKE
    plate: float  # [kg], the larger, plus the bore collar
    departing: float  # [kg] after the burn, plate dropped
    impactor: float  # [kg] of impactor spent
    departing_per_impactor: float
    carbon: str


def stage(route: Route, k: float, frozen: bool = True, bowl: str = DESIGN_PLATE) -> Stage:
    """Fly the departure on the hole plate: iterate the burn time, its loss and the mass ratio."""
    v_e = hole_ve(bowl, route.theta, k, route.w, frozen=frozen)
    dv = route.dv
    for _ in range(4):
        mr = math.exp(dv / v_e)
        burn = dv / (G_MEAN * G0)  # the mean acceleration is the ideal burn over its length
        dv = route.dv + burn_loss(route.dv, v_e, burn)
    mr = math.exp(dv / v_e)
    charge_mass = STACK * (1.0 - 1.0 / mr)
    impulse = charge_mass * v_e
    pulses = burn * CADENCE
    per_pulse = impulse / pulses
    sheet = per_pulse / (k * v_e)
    floor = (
        growth_survival() * SUN_PLATE_AREA * (0.5 * sheet * route.w**2 / SUN_PLATE_RATED) ** (2 / 3)
    )
    end = STACK / mr
    inertia = per_pulse**2 / (2.0 * STROKE * end * G_CAP * G0)
    plate = max(floor, inertia) + BORE_COLLAR
    return Stage(
        f"{route.name}, {route.sign}",
        k,
        v_e,
        v_e / G0,
        route.dv,
        dv - route.dv,
        mr,
        burn,
        pulses,
        per_pulse,
        sheet,
        floor,
        inertia,
        plate,
        end - plate,
        charge_mass / k,
        (end - plate) / (charge_mass / k),
        "atomic" if frozen else "condensed",
    )


def growth_survival() -> float:
    return float(charge.SURVIVAL_AREAL)


def _write(path: object, rows: list[dict[str, object]]) -> None:
    from pathlib import Path

    assert isinstance(path, Path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    law: list[dict[str, object]] = []
    print("  impulse per m w: flat ceiling / bowl ceiling / hole design (paraboloid 0.5 at 1.5 F)")
    for k in KS:
        row: dict[str, object] = {"k": k}
        for t in (90.0, 89.2, 117.0, 129.9, 180.0):
            row[f"flat_ceiling_{t:g}"] = flat_ceiling(k, t)
            row[f"bowl_ceiling_{t:g}"] = bowl_ceiling(k, t)
        for t in (89.2, 117.0, 129.9, 180.0):
            for plate in PLATES:
                row[f"capture {plate} {t:g}"] = capture(plate, t, k)
        for r in ROUTES:
            tag = f"{r.name} {r.sign}"
            row[f"hole_ve {tag}"] = hole_ve(DESIGN_PLATE, r.theta, k, r.w)
            row[f"hole_ve_condensed {tag}"] = hole_ve(DESIGN_PLATE, r.theta, k, r.w, frozen=False)
            row[f"flat_ve {tag}"] = flat_ve(r.theta, k, r.w)
            row[f"flat_ve_condensed {tag}"] = flat_ve(r.theta, k, r.w, False)
        law.append(row)
        loss_flat = 1.0 - flat_ceiling(k, 129.9) / flat_ceiling(k, 90.0)
        loss_bowl = 1.0 - bowl_ceiling(k, 129.9) / bowl_ceiling(k, 90.0)
        run = 1.0 - capture(DESIGN_PLATE, 129.9, k) / capture(DESIGN_PLATE, 89.2, k)
        print(
            f"  k {k:3.0f} | 129.9 vs 90 deg per m w: flat ceiling -{loss_flat:.1%}, bowl ceiling"
            f" -{loss_bowl:.1%}, bowl run (vs 89.2) -{run:.1%} | per kg at 84 vs 75 km/s:"
            f" bowl run {84.0 / 75.0 * (1 - run) - 1:+.1%}"
        )
        for r in ROUTES:
            hv, fv = hole_ve(DESIGN_PLATE, r.theta, k, r.w), flat_ve(r.theta, k, r.w)
            print(
                f"      {r.name:<10} {r.sign:<21} {r.theta:5.1f} deg {r.w / 1e3:5.1f} km/s |"
                f" hole bowl v_e {hv / 1e3:5.2f} km/s ({hv / G0:4.0f} s)"
                f" | canted flat v_e {fv / 1e3:5.2f} km/s ({fv / G0:4.0f} s)"
            )
    holes: list[dict[str, object]] = []
    for plate in PLATES:
        for t in (180.0, 129.9, 117.0):
            r_h, z_h = hole_geometry(plate, t) if t < 180 else (0.0, 0.0)
            for k in (10.0, 20.0, 30.0, 40.0, 60.0):
                for rad in HOLE_RADII:
                    if t == 180.0:
                        # On axis: the bore's own ring (the first band).
                        row0 = _runs("x1")[(plate, t, k)]
                        b0 = row0["bands"]
                        assert isinstance(b0, list) and isinstance(b0[0], dict)
                        j0 = float(str(row0["impulse_per_mv"]))
                        share = float(b0[0]["impulse"]) * math.pi * rad**2 / j0
                        mass = float(b0[0]["mass_in"]) * math.pi * rad**2 / (1.0 + k)
                    else:
                        share, mass = hole_loss(plate, t, k, rad)
                    holes.append(
                        {
                            "plate": plate,
                            "theta": t,
                            "k": k,
                            "hole_radius": rad,
                            "hole_r": r_h,
                            "hole_z": z_h,
                            "impulse_share": share,
                            "mass_share": mass,
                        }
                    )
    for h in holes:
        if h["k"] in (10.0, 30.0) and h["hole_radius"] == 0.5:
            print(
                f"  hole {h['plate']:<26} {h['theta']:5.1f} deg k {h['k']:3.0f}: at r"
                f" {h['hole_r']:4.2f} m, a 0.5 m hole takes {h['impulse_share']:.2%} of the"
                f" impulse and {h['mass_share']:.2%} of the gas"
            )
    burns: list[dict[str, object]] = []
    for t in BURN_TIMES:
        loss = burn_loss(4.21e3, 12.0e3, t)
        burns.append({"burn_s": t, "dv": 4.21e3, "exhaust": 12.0e3, "loss": loss})
        print(f"  4.21 km/s burn over {t:5.0f} s, fixed direction: loss {loss:5.0f} m/s")
    stages = [
        asdict(stage(r, k, fz))
        for r in ROUTES[1::2]
        for fz in (True, False)
        for k in (5.0, 10.0, 15.0, 20.0, 30.0)
    ]
    for s in stages:
        print(
            f"  {s['route']:<22} k {s['k']:3.0f}: v_e {s['v_e'] / 1e3:5.2f} km/s"
            f" ({s['isp']:4.0f} s)"
            f" | loss {s['loss']:4.0f} m/s, MR {s['mass_ratio']:.3f}, burn {s['burn_s']:5.0f} s"
            f" | {s['pulse_impulse'] / 1e6:4.2f} MN s, {s['impactor_per_pulse']:4.1f} kg per pulse"
            f" | plate floor {s['plate_floor'] / 1e3:4.1f} t, inertia"
            f" {s['plate_inertia'] / 1e3:4.1f} t -> {s['plate'] / 1e3:4.1f} t"
            f" | {s['departing_per_impactor']:4.1f} kg departs per kg impactor (C {s['carbon']})"
        )
    _write(OUT_LAW, law)
    _write(OUT_HOLE, holes)
    _write(OUT_BURN, burns)
    _write(OUT_STAGE, stages)
    print(f"python: wrote {OUT_LAW}, {OUT_HOLE}, {OUT_BURN}, {OUT_STAGE}")


if __name__ == "__main__":
    main()
