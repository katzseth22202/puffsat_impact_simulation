"""The bore collar: how much it weighs, how much the back-splash eats, and what its vapour pushes.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.collar

Result 10's disk-into-charge runs throw a narrow back-splash jet onto the ring just outside the
bore. Its peak pressure (GPa) exceeds any facesheet, and it grows on refinement. The collar is the
sacrificial ring that takes it. It is priced here from three sources:

- **The ring loads** from 2-D (`sweep --charge-mix K collar`). For each ring the run tallies the
  impulse per area, the arriving column, its speed, its energy and the peak pressure, per kg of
  disk. The point fireball's rows give the same rings' load without the jet.
- **The structural mass.** The pulse is microseconds long, far shorter than any structural
  response, so the ring is damaged by its impulse, not its peak (the impulsive branch of a
  pressure-impulse diagram). The collar is sized to be **impulse-matched**: its areal mass is
  the plate's (`charge.SURVIVAL_AREAL`, 30 kg/m^2) scaled by the ring's impulse over a reference,
  so it takes the same velocity kick per pulse as that reference and loads its mount no harder.
  Two references bracket it. Against the plate's *mean* impulse (no structural spreading at all)
  is the upper bound; it charges the smooth fireball's own centre-heavy load as well as the jet's.
  Against the *point fireball's* impulse on the same ring (the load the ledger assumed) charges
  the jet's increment alone, for a plate already built to the smooth load. All the impulses scale
  with the disk's mass, so this mass does not.
- **The erosion and its thrust** from 1-D (`sweep --collar-ablation`). Each ring's flow is re-run
  as the study's radiation-hydrodynamic slug bounce on the CH2 table, against a rigid wall and
  against graphite that ablates at `q_wall / Q*`. The vapour is injected into the gas at the wall.
  The ablated mass per pulse is the erosion allowance, and the extra wall impulse is what the
  vapour adds. Graphite's `Q*` is bracketed between sublimation (59 MJ/kg) and a third of it. The
  vapour's shielding curtain is off: calibrated for ablated mass well below the cloud's, it makes
  energy here (the sweep's `COLLAR_KAPPA_VAPOR`). Each row is checked: its departing momentum must
  not need more energy than the slug brought, `(J - p_in)^2 / 2M <= KE`.
- **The jet's heat.** The 1-D slug arrives cold; the 2-D jet arrives carrying several times its
  kinetic energy as enthalpy. The worst case scales the erosion up by that ratio, a ceiling.

The eroded collar is reaction mass: it leaves with the exhaust. So the ledger treats it as
propellant with its own exhaust speed (`added impulse / eroded mass`), and recomputes the mass
ratio with it charged. Mechanical erosion (spall, crushing) is not modelled; a bonded, backed
collar has no free rear face to spall.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

from puffsat.walled_nozzle import wall
from puffsat.xray_collision import charge

OUT = charge.OUT_DIR / "collar_ledger.csv"
LEDGER = charge.OUT_DIR / "charge_ledger.csv"
ABLATION = charge.OUT_DIR / "collar_ablation.jsonl"
FINAL_MASS = 100.0e3  # [kg], the handoff's nominal final mass
GRAPHITE_DENSITY = wall.GRAPHITE_DENSITY  # [kg/m^3], to turn eroded mass into depth
#: The 2-D case label's plate, to the charge ledger's.
PLATES = {
    "bowl 0.30 at focus": "paraboloid 0.30 at focus (parent bowl)",
    "paraboloid 0.50 at 1.5 F": "paraboloid 0.50 at 1.5 focus",
}
#: `(Q*, kappa_vapor)` read as the best and the worst erosion (the worst also hot-jet scaled).
BEST = (59.2e6, 0.0)
WORST = (20.0e6, 0.0)
#: Rings arriving with less [kg/m^2 per kg of disk] are loaded by gas sliding along the wall and
#: get no 1-D run (the sweep's `COLLAR_MIN_COLUMN`); their erosion is the plate-wide ablation item.
MIN_COLUMN = 1.0e-4


@dataclass(frozen=True)
class Ring:
    r_lo: float
    r_hi: float
    area: float  # [m^2]
    impulse: float  # [Pa s] per kg of disk
    impulse_smooth: float  # the point fireball's, same ring
    peak_pressure: float  # [Pa] per kg of disk
    speed: float  # [m/s] arriving
    column: float  # [kg/m^2] per kg of disk
    areal_mass: float  # [kg/m^2], impulse-matched to the plate's mean
    areal_mass_jet: float  # [kg/m^2], impulse-matched to the point fireball on this ring
    energy_ratio: float  # 2-D arriving energy over the 1-D slug's kinetic energy
    eroded_best: float  # [kg] per kg of disk per pulse
    eroded_worst: float  # Q* 20 MJ/kg, scaled by `energy_ratio`
    added_impulse_best: float  # [N s] per kg of disk
    added_impulse_worst: float  # Q* 20 MJ/kg, unscaled (the conservative pairing)
    energy_ok_best: bool  # the 1-D row's departing momentum is energetically possible
    energy_ok_worst: bool


@dataclass(frozen=True)
class Collar:
    plate: str
    k: float
    grid: str
    rings: tuple[Ring, ...]
    plate_mean_impulse: float  # [Pa s] per kg of disk, the 2-D run's over the projected plate
    collar_mass: float  # [kg], impulse-matched to the mean, over the rings that need it
    collar_extra: float  # [kg] above the plate's own 30 kg/m^2 over the same rings
    collar_radius: float  # [m]
    jet_extra: float  # [kg] above 30 kg/m^2, matched to the point fireball's impulse instead
    jet_radius: float  # [m]
    smooth_extra: float  # [kg] the point fireball's own excess over its own mean (no jet at all)
    eroded_per_final_best: float  # [kg per kg of final mass] over the boost
    eroded_per_final_worst: float
    disks_per_final: float  # [kg of disk per kg of final mass], the ledger's projectile bill
    max_depth_best: float  # [m] of graphite eroded from the worst ring over the boost (FINAL_MASS)
    max_depth_worst: float
    added_share_best: float  # of the pulse's impulse
    added_share_worst: float
    collar_v_e_best: float  # [m/s] the eroded mass's own exhaust speed
    collar_v_e_worst: float
    mass_ratio: float  # the ledger's
    mass_ratio_with_collar_best: float
    mass_ratio_with_collar_worst: float


def _jsonl(path: Path) -> list[dict[str, object]]:
    return charge.load_jsonl(path) if path.exists() else []


def _f(row: Mapping[str, object], key: str) -> float:
    v = row[key]
    assert isinstance(v, (int, float, str))
    return float(v)


def _ledger_row(plate: str, k: float) -> dict[str, str]:
    with LEDGER.open() as fh:
        for row in csv.DictReader(fh):
            if (
                row["case"] == plate
                and row["grid"] == "x1"
                and abs(float(row["k"]) - k) < 1e-9
                and abs(float(row["projectile_mass"]) - 1.0) < 1e-9
            ):
                return row
    raise KeyError((plate, k))


def _ablation(
    rows: list[dict[str, object]], k: float, case: str, r_lo: float
) -> dict[tuple[float, float], dict[str, object]]:
    out = {}
    for r in rows:
        if abs(_f(r, "k") - k) < 1e-9 and r["case"] == case and abs(_f(r, "r_lo") - r_lo) < 1e-9:
            out[(_f(r, "q_star"), _f(r, "kappa_vapor"))] = r
    return out


def collar(
    plate_label: str,
    k: float,
    coarse: bool = False,
    ablation: list[dict[str, object]] | None = None,
) -> Collar:
    """The collar for one plate and k, from the solid-charge 2-D row and its 1-D re-runs."""
    grid = "_coarse" if coarse else ""
    rows = _jsonl(charge.OUT_DIR / f"charge_collar_k{k:.0f}{grid}.jsonl")
    solid = next(r for r in rows if r["case"] == f"{plate_label}: disk into solid charge")
    # The point fireball always runs at the coarse cell, so the coarse pass reuses its row.
    smooth_rows = _jsonl(charge.OUT_DIR / f"charge_collar_k{k:.0f}.jsonl")
    smooth = next(r for r in smooth_rows if r["case"] == f"{plate_label}: point fireball")
    ledger = _ledger_row(PLATES[plate_label], k)
    pulse = float(ledger["impulse_per_projectile"])
    table = ablation if ablation is not None else _jsonl(ABLATION)
    v = charge.V_REL
    rings = []
    solid_bands = solid["bands"]
    smooth_bands = smooth["bands"]
    assert isinstance(solid_bands, list) and isinstance(smooth_bands, list)
    # The mean is the same 2-D run's whole impulse over the same projected plate, so a ring's
    # ratio to it is like for like (the rings tile the plate from the bore's edge to the rim).
    tiled = sum(
        float(b["impulse"]) * v * math.pi * (float(b["r_hi"]) ** 2 - float(b["r_lo"]) ** 2)
        for b in solid_bands
    )
    r_in, r_out = float(solid_bands[0]["r_lo"]), float(solid_bands[-1]["r_hi"])
    mean = tiled / (math.pi * (r_out**2 - r_in**2))
    # The rings tile the plate from the bore's edge to the rim, so their impulses must add up to
    # the run's own plate impulse less the bore's (the record is per M u_rms; u_rms = sqrt(k)/(1+k)
    # of the disk's speed for a merged mass of 1 + k disks).
    plate = (
        (float(str(solid["wall_impulse"])) - float(str(solid["hole_impulse"]))) * math.sqrt(k) * v
    )
    assert abs(tiled / plate - 1.0) < 1e-3, (plate_label, k, tiled, plate)
    mean_smooth = sum(
        float(b["impulse"]) * v * math.pi * (float(b["r_hi"]) ** 2 - float(b["r_lo"]) ** 2)
        for b in smooth_bands
    ) / (math.pi * (r_out**2 - r_in**2))
    for b, s in zip(solid_bands, smooth_bands, strict=True):
        lo, hi = float(b["r_lo"]), float(b["r_hi"])
        area = math.pi * (hi * hi - lo * lo)
        j = float(b["impulse"]) * v
        a = _ablation(table, k, str(solid["case"]), lo)

        def eroded(
            key: tuple[float, float],
            a: dict[tuple[float, float], dict[str, object]] = a,
            area: float = area,
        ) -> float:
            return _f(a[key], "ablated_mass") * area if key in a else 0.0

        def energy_ok(
            key: tuple[float, float],
            a: dict[tuple[float, float], dict[str, object]] = a,
        ) -> bool:
            if key not in a:
                return False
            col, u = _f(a[key], "column"), _f(a[key], "v")
            out = _f(a[key], "impulse_ablating") - col * u
            mass = col + _f(a[key], "ablated_mass")
            return out * out / (2.0 * mass) <= 0.5 * col * u * u

        def added(
            key: tuple[float, float],
            a: dict[tuple[float, float], dict[str, object]] = a,
            area: float = area,
            j: float = j,
        ) -> float:
            # Only an energetically possible row earns thrust; its erosion is charged regardless.
            if key not in a or not energy_ok(key, a):
                return 0.0
            # The vapour pushes on what the arriving column heats: the 1-D slug's own extra
            # impulse per area, not a share of the ring's 2-D impulse, most of which is pressure
            # from gas already stagnated against the wall.
            gain = _f(a[key], "impulse_ablating") - _f(a[key], "impulse_rigid")
            return gain * area

        column = float(b["mass_in"])
        speed = float(b["speed_in"]) * v
        ke = 0.5 * column * speed * speed
        ratio = max(float(b["energy_in"]) * v * v / ke, 1.0) if column > MIN_COLUMN else math.nan

        rings.append(
            Ring(
                r_lo=lo,
                r_hi=hi,
                area=area,
                impulse=j,
                impulse_smooth=float(s["impulse"]) * v,
                peak_pressure=float(b["peak_pressure"]) * v * v,
                speed=speed,
                column=column,
                areal_mass=charge.SURVIVAL_AREAL * j / mean,
                areal_mass_jet=charge.SURVIVAL_AREAL * j / (float(s["impulse"]) * v),
                energy_ratio=ratio,
                eroded_best=eroded(BEST),
                eroded_worst=eroded(WORST) * ratio if column > MIN_COLUMN else 0.0,
                added_impulse_best=added(BEST),
                added_impulse_worst=added(WORST),
                energy_ok_best=energy_ok(BEST),
                energy_ok_worst=energy_ok(WORST),
            )
        )
    # The collar runs out to the last ring whose impulse beats the plate's mean.
    needed = [r for r in rings if r.areal_mass > charge.SURVIVAL_AREAL]
    collar_mass = sum(r.areal_mass * r.area for r in needed)
    extra = sum((r.areal_mass - charge.SURVIVAL_AREAL) * r.area for r in needed)
    jet = [r for r in rings if r.areal_mass_jet > charge.SURVIVAL_AREAL]
    jet_extra = sum((r.areal_mass_jet - charge.SURVIVAL_AREAL) * r.area for r in jet)
    smooth_extra = sum(
        max(charge.SURVIVAL_AREAL * (r.impulse_smooth / mean_smooth - 1.0), 0.0) * r.area
        for r in rings
    )
    per_final = float(ledger["projectile_per_final"])
    v_e = float(ledger["v_e"])
    dv = charge.boost_dv()

    def charged(eroded: float, added: float) -> tuple[float, float, float]:
        """`(added share, collar v_e, mass ratio)` with the eroded collar charged as propellant."""
        share = added / pulse
        v_e_new = v_e * (1.0 + share) / (1.0 + eroded / (1.0 + k))
        return share, added / eroded if eroded > 0 else math.nan, math.exp(dv / v_e_new)

    best = charged(sum(r.eroded_best for r in rings), sum(r.added_impulse_best for r in rings))
    worst = charged(sum(r.eroded_worst for r in rings), sum(r.added_impulse_worst for r in rings))
    return Collar(
        plate=plate_label,
        k=k,
        grid="x0.5" if coarse else "x1",
        rings=tuple(rings),
        plate_mean_impulse=mean,
        collar_mass=collar_mass,
        collar_extra=extra,
        collar_radius=max((r.r_hi for r in needed), default=0.0),
        jet_extra=jet_extra,
        jet_radius=max((r.r_hi for r in jet), default=0.0),
        smooth_extra=smooth_extra,
        eroded_per_final_best=sum(r.eroded_best for r in rings) * per_final,
        eroded_per_final_worst=sum(r.eroded_worst for r in rings) * per_final,
        disks_per_final=per_final,
        max_depth_best=max(r.eroded_best / r.area for r in rings)
        * per_final
        * FINAL_MASS
        / GRAPHITE_DENSITY,
        max_depth_worst=max(r.eroded_worst / r.area for r in rings)
        * per_final
        * FINAL_MASS
        / GRAPHITE_DENSITY,
        added_share_best=best[0],
        added_share_worst=worst[0],
        collar_v_e_best=best[1],
        collar_v_e_worst=worst[1],
        mass_ratio=math.exp(dv / v_e),
        mass_ratio_with_collar_best=best[2],
        mass_ratio_with_collar_worst=worst[2],
    )


def main() -> None:
    ablation = _jsonl(ABLATION)
    out = []
    for k in (30.0, 100.0):
        for plate in PLATES:
            for coarse in (False, True):
                try:
                    c = collar(plate, k, coarse, ablation if not coarse else [])
                except (StopIteration, FileNotFoundError):
                    continue
                print(
                    f"  k {k:3.0f} {plate:<26} {c.grid:>4}: "
                    f"mean {c.plate_mean_impulse:7.0f} Pa s/kg"
                    f" | collar to r {c.collar_radius:.1f} m: {c.collar_mass:5.0f} kg"
                    f" (+{c.collar_extra:4.0f} over the plate; smooth fireball alone"
                    f" +{c.smooth_extra:4.0f}) | jet only: +{c.jet_extra:4.0f} kg to r"
                    f" {c.jet_radius:.1f} m"
                )
                for r in c.rings:
                    print(
                        f"      r {r.r_lo:.2f}-{r.r_hi:.2f}: J {r.impulse:8.0f} Pa s/kg "
                        f"(smooth {r.impulse_smooth:6.0f}) x{r.impulse / c.plate_mean_impulse:5.1f}"
                        f" peak {r.peak_pressure / 1e9:6.2f} GPa  v {r.speed / 1e3:5.0f} km/s"
                        f"  col {r.column:.2e}"
                        f"  -> {r.areal_mass:4.0f} / {r.areal_mass_jet:4.0f} kg/m^2"
                        f" | eroded {1e3 * r.eroded_best:.3g}-{1e3 * r.eroded_worst:.3g} g/kg"
                        f" (heat x{r.energy_ratio:.1f},"
                        f" E-ok {r.energy_ok_best}/{r.energy_ok_worst})"
                        f" over the boost at {FINAL_MASS / 1e3:.0f} t: "
                        f"{r.eroded_best / r.area * c.disks_per_final * FINAL_MASS:.0f}-"
                        f"{r.eroded_worst / r.area * c.disks_per_final * FINAL_MASS:.0f} kg/m^2"
                    )
                if not coarse:
                    print(
                        f"      eroded over the boost {c.eroded_per_final_best:.2e}-"
                        f"{c.eroded_per_final_worst:.2e} kg/kg final "
                        f"({c.eroded_per_final_best * FINAL_MASS:.0f}-"
                        f"{c.eroded_per_final_worst * FINAL_MASS:.0f} kg at "
                        f"{FINAL_MASS / 1e3:.0f} t); vapour adds "
                        f"{c.added_share_best:.2%}-{c.added_share_worst:.2%} of the pulse, "
                        f"v_e {c.collar_v_e_best / 1e3:.0f}-{c.collar_v_e_worst / 1e3:.0f} km/s;"
                        f" MR {c.mass_ratio:.3f} -> {c.mass_ratio_with_collar_best:.3f}-"
                        f"{c.mass_ratio_with_collar_worst:.3f}"
                    )
                row: dict[str, object] = {
                    key: val for key, val in asdict(c).items() if key != "rings"
                }
                for i, r in enumerate(c.rings):
                    for key, val in asdict(r).items():
                        row[f"ring{i}_{key}"] = val
                out.append(row)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"python: wrote {len(out)} rows -> {OUT}")


if __name__ == "__main__":
    main()
