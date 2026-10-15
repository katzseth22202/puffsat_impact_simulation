"""Result 12: optimize k at each wave and price the 3S cycle (Earth -> Jupiter -> Sun -> Earth).

    PYTHONPATH=python uv run python -m puffsat.xray_collision.growth

The **Earth wave** (the Sun ships' released mass) now pushes departing ships onto the 3S Jovian
cycle instead of an 11 km/s escape. They fall from Jupiter to the dive periapsis, where the **Sun
wave** boosts them, and their whole remainder returns as the next Earth wave (ADR-0062).

**Cycle geometry** is the calcs repo's (`aim_is_all_you_need`, commit 868e88a), from
`src.jovian_solar_dive_cycle`:
`solve_synodic_closure(3, dive_solar_radii=R, return_excess=X)` followed by
`departure_nozzle_ledger(c.departure_excess, c.departure_aim_deg, c.return_excess,
c.push_axis_deg, 30.0, 0.6)`. The departure burn starts at the 10-day parking orbit's periapsis.
Its thrust axis meets the stream at an **impact angle** above 90 degrees, so the stream sweeps
*past* the plate (Result 12's canted plate). The values are copied in `CLOSURES`; this repository
cannot import that one.

**Earth stage.** The canted strips (`sweep --canted-2d`, `canted_2d.jsonl`) give `J/(m w)` against
the canted ceiling `cos theta + sqrt(1 + k - sin^2 theta)`. Their share of it, at the best
standoff, is fitted in k and carried to the closure's own angle. Polyethylene's rebound is swapped
in at the sliding fireball's `u_rms`, as in Results 10-11. The bond toll comes from
`freeze.freeze` at the closure's closing speed, bracketed by carbon left atomic or condensed.

**Sun stage.** It uses the δ/D 0.5 paraboloid at 1.5 F from `charge_ledger.csv` (400 km/s,
bonds charged), reduced to the real disk-into-charge collision by `charge.mixing_factor`. Its
lossless capture does not depend on speed (`drift/u_rms = 1/sqrt k` either way), so `v_e` is
carried to the dive's closing speed `v_S` (twice the periapsis speed, `mission.periapsis`) by
`v_S/400 km/s`, the rebound ratio at the new `u_rms`, and the bond factor at the new heat.

**Objective.** Per kg of Earth-wave impactor:
- `D = k_E/(MR_E - 1)` departs;
- a Sun ship returns `(1 - f_dead)/MR_S` of its mass;
- it spends `(1 - 1/MR_S)/k_S` of Sun-wave impactor per kg.

So the impactor returned per kg of *all* impactor spent in a 3S cycle is

`G = D (1 - f_dead)/MR_S / (1 + D (1 - 1/MR_S)/k_S)`.

It is maximized over `k_E` and `k_S`.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache

import numpy as np

from puffsat.xray_collision import charge, freeze, mission
from puffsat.xray_collision.eos import EosTable

CANTED = charge.OUT_DIR / "canted_2d.jsonl"
OUT_STAGE = charge.OUT_DIR / "growth_earth_stage.csv"
OUT = charge.OUT_DIR / "growth_cycle.csv"
SUN_PLATE = "paraboloid 0.50 at 1.5 focus"
SUN_LEDGER_SPEED = charge.V_REL  # [m/s], the closing speed charge_ledger.csv was priced at
#: The Earth stage's rebound is looked up at a column of (1 + k) x this impactor mass over the
#: 20 m plate's area, the overtake ledger's convention (e_eff depends weakly on it).
REBOUND_IMPACTOR = 5.0  # [kg]
#: The searches stay inside the simulated ranges: the canted strips run k = 5-60 and the Sun
#: ledger k = 20-100. An optimum on either edge is flagged rather than extrapolated past it.
K_EARTH = tuple(float(k) for k in np.geomspace(5.0, 60.0, 50))
K_SUN = tuple(float(k) for k in np.geomspace(20.0, 100.0, 40))


@dataclass(frozen=True)
class Closure:
    """One 3S closure from the calcs repo (`aim_is_all_you_need` 868e88a)."""

    name: str
    dive_rsun: float
    return_excess: float  # [km/s], the closure's input
    departure_excess: float  # [km/s]
    departure_aim_deg: float
    departure_dv: float  # [m/s], from the 10-day parking orbit's periapsis
    impact_angle_deg: tuple[float, float]  # stream vs thrust axis, start and end of the burn
    stream_closing: tuple[float, float]  # [m/s], vehicle frame, start and end
    periapsis_boost: float  # [m/s] the Sun stage must give
    cycle_years: float

    @property
    def theta(self) -> float:
        return 0.5 * (self.impact_angle_deg[0] + self.impact_angle_deg[1])

    @property
    def w_earth(self) -> float:
        return 0.5 * (self.stream_closing[0] + self.stream_closing[1])

    @property
    def w_sun(self) -> float:
        return mission.periapsis(self.dive_rsun).closing_speed


CLOSURES = (
    Closure(
        "9.86 R_sun, 100 km/s crossing",
        9.86,
        90.0,
        10.43,
        -13.5,
        4.15e3,
        (100.0, 102.3),
        (102.8e3, 103.6e3),
        20.45e3,
        3.28,
    ),
    Closure(
        "4 R_sun, 150 km/s return",
        4.0,
        150.0,
        10.54,
        -12.1,
        4.23e3,
        (93.9, 95.5),
        (158.2e3, 158.5e3),
        35.04e3,
        3.28,
    ),
)


def canted_ceiling(k: float, theta_deg: float) -> float:
    t = math.radians(theta_deg)
    return math.cos(t) + math.sqrt(1.0 + k - math.sin(t) ** 2)


@cache
def _table() -> EosTable:
    return charge.ch2_table()


@cache
def _bond_loss(k: float, w: float, frozen: bool) -> float:
    f = freeze.freeze(k, _table(), v_rel=w)
    return f.velocity_loss_bonds_carbon_frozen if frozen else 0.0


@dataclass(frozen=True)
class EarthStageFit:
    """The canted plate's share of its ceiling, fitted as `a + b ln k` over the design rows at
    the standoff that gives the most impulse at each k."""

    family: str
    standoff: tuple[float, ...]
    k: tuple[float, ...]
    share: tuple[float, ...]
    a: float
    b: float
    confined_ratio: float

    def at(self, k: float) -> float:
        return self.a + self.b * math.log(k)


#: The canted strip families: the slab and square charges are real sheet-on-charge collisions,
#: the point fireball the merged state they make (Results 10-11's production model).
FAMILIES = ("design", "square charge", "point fireball")
FAMILY_LABEL = {
    "design": "slab charge",
    "square charge": "square charge",
    "point fireball": "point",
}


def earth_fit(family: str = "design") -> EarthStageFit:
    rows = [r for r in charge.load_jsonl(CANTED) if r["family"] == family]
    best: dict[float, dict[str, object]] = {}
    for r in rows:
        k = float(str(r["k"]))
        if k not in best or float(str(r["impulse"])) > float(str(best[k]["impulse"])):
            best[k] = r
    ks = sorted(best)
    share = [float(str(best[k]["ceiling_share"])) for k in ks]
    b, a = np.polyfit(np.log(ks), share, 1)
    return EarthStageFit(
        family=family,
        standoff=tuple(float(str(best[k]["standoff"])) for k in ks),
        k=tuple(ks),
        share=tuple(share),
        a=float(a),
        b=float(b),
        confined_ratio=float(str(rows[0]["confined_ratio"])),
    )


def earth_ve(fit: EarthStageFit, c: Closure, k: float, frozen: bool) -> tuple[float, float]:
    """`(v_e, e_eff)` of the canted Earth stage at `k`."""
    w = c.w_earth
    j = fit.at(k) * canted_ceiling(k, c.theta)
    u_rms = w * math.sqrt(k) / (1.0 + k)
    column = (1.0 + k) * REBOUND_IMPACTOR / (2.0 * math.pi * charge.PLATE_R**2)
    e = charge.rebound(u_rms, column)
    v_e = j * w * (1.0 + e) / fit.confined_ratio / k * (1.0 - _bond_loss(k, w, frozen))
    return v_e, e


@cache
def _sun_rows() -> tuple[tuple[float, float, float], ...]:
    """`(k, v_e with bonds charged, e_eff)` of the Sun plate at 400 km/s, 1 kg projectile."""
    seen: dict[float, tuple[float, float]] = {}
    with (charge.OUT_DIR / "charge_ledger.csv").open() as fh:
        for r in csv.DictReader(fh):
            if (
                r["case"] == SUN_PLATE
                and r["grid"] == "x1"
                and abs(float(r["projectile_mass"]) - 1.0) < 1e-9
            ):
                k = float(r["k"])
                # The real collision, not the point fireball (charge.mixing_factor).
                v = float(r["isp_bonds_frozen"]) * charge.G0 * charge.mixing_factor(SUN_PLATE, k)
                seen.setdefault(k, (v, float(r["e_eff"])))
    return tuple((k, v, e) for k, (v, e) in sorted(seen.items()))


def _bond_factor(k: float, v: float) -> float:
    heat_per_merged = 0.5 * v * v * k / (1.0 + k) / (1.0 + k)
    return math.sqrt(max(1.0 - charge.PE_BONDS * k / (k + 1.0) / heat_per_merged, 0.0))


def sun_ve(k: float, v: float) -> float:
    """The Sun plate's `v_e` (bonds charged) at `k`, carried from 400 km/s to closing speed `v`."""
    rows = _sun_rows()
    ks = np.log([r[0] for r in rows])
    lk = math.log(k)
    v400 = math.exp(float(np.interp(lk, ks, np.log([r[1] for r in rows]))))
    e400 = float(np.interp(lk, ks, [r[2] for r in rows]))
    u_rms = v * math.sqrt(k) / (1.0 + k)
    column = (1.0 + k) / (2.0 * math.pi * charge.PLATE_R**2)  # the ledger's 1 kg convention
    e = charge.rebound(u_rms, column)
    scale = v / SUN_LEDGER_SPEED * (1.0 + e) / (1.0 + e400)
    return v400 * scale * _bond_factor(k, v) / _bond_factor(k, SUN_LEDGER_SPEED)


@dataclass(frozen=True)
class CycleOptimum:
    closure: str
    geometry: str
    bonds: str
    dead_fraction: float
    k_earth: float
    k_sun: float
    v_e_earth: float  # [m/s]
    v_e_sun: float
    mass_ratio_earth: float
    mass_ratio_sun: float
    departing_per_earth_impactor: float  # D
    sun_impactor_per_earth_impactor: float
    growth_per_total_impactor: float  # G, the objective, per 3S
    growth_per_earth_impactor: float  # without charging the Sun wave
    sun_wave_share: float  # of all impactor spent
    doubling_years: float
    earth_orbit_per_departing: float  # MR_E: kg already in Earth orbit per kg departing
    at_edge: str  # which k sits on its search edge ("" if neither): the optimum lies beyond it


def optimize(c: Closure, fit: EarthStageFit, frozen: bool, dead: float) -> CycleOptimum:
    sun = {ks: sun_ve(ks, c.w_sun) for ks in K_SUN}
    best: CycleOptimum | None = None
    for ke in K_EARTH:
        ve, _ = earth_ve(fit, c, ke, frozen)
        if ve <= 0.0:
            continue
        mre = math.exp(c.departure_dv / ve)
        d = ke / (mre - 1.0)
        for ks, vs in sun.items():
            mrs = math.exp(c.periapsis_boost / vs)
            sun_imp = d * (1.0 - 1.0 / mrs) / ks
            ret = d * (1.0 - dead) / mrs
            g = ret / (1.0 + sun_imp)
            if best is None or g > best.growth_per_total_impactor:
                best = CycleOptimum(
                    closure=c.name,
                    geometry=FAMILY_LABEL[fit.family],
                    bonds="carbon atomic" if frozen else "bonds returned",
                    dead_fraction=dead,
                    k_earth=ke,
                    k_sun=ks,
                    v_e_earth=ve,
                    v_e_sun=vs,
                    mass_ratio_earth=mre,
                    mass_ratio_sun=mrs,
                    departing_per_earth_impactor=d,
                    sun_impactor_per_earth_impactor=sun_imp,
                    growth_per_total_impactor=g,
                    growth_per_earth_impactor=ret,
                    sun_wave_share=sun_imp / (1.0 + sun_imp),
                    doubling_years=c.cycle_years * math.log(2.0) / math.log(g)
                    if g > 1.0
                    else math.inf,
                    earth_orbit_per_departing=mre,
                    at_edge=" ".join(
                        name
                        for name, v, grid in (("k_E", ke, K_EARTH), ("k_S", ks, K_SUN))
                        if v in (grid[0], grid[-1])
                    ),
                )
    assert best is not None
    return best


def main() -> None:
    present = {str(r["family"]) for r in charge.load_jsonl(CANTED)}
    fits = [earth_fit(f) for f in FAMILIES if f in present]
    for fit in fits:
        print(
            f"  {FAMILY_LABEL[fit.family]:<14} share of the canted ceiling at the best standoff: "
            + ", ".join(
                f"k {k:.0f} {s:.2f} (h {h:.0f})"
                for k, s, h in zip(fit.k, fit.share, fit.standoff, strict=True)
            )
        )
    stage: list[dict[str, object]] = []
    for fit in fits:
        for c in CLOSURES:
            for k in (5.0, 10.0, 20.0, 30.0, 40.0, 60.0):
                ve_c, e = earth_ve(fit, c, k, True)
                ve_b, _ = earth_ve(fit, c, k, False)
                stage.append(
                    {
                        "geometry": FAMILY_LABEL[fit.family],
                        "closure": c.name,
                        "k": k,
                        "theta_deg": c.theta,
                        "w_earth": c.w_earth,
                        "share": fit.at(k),
                        "e_eff": e,
                        "v_e_carbon_atomic": ve_c,
                        "v_e_bonds_returned": ve_b,
                        "mass_ratio_carbon_atomic": math.exp(c.departure_dv / ve_c),
                    }
                )
    out = [
        optimize(c, fit, frozen, dead)
        for fit in fits
        for c in CLOSURES
        for frozen in (True, False)
        for dead in (0.0, mission.plate_dead_fraction())
    ]
    for o in out:
        print(
            f"  {o.geometry:<13} {o.closure:<30} {o.bonds:<14} dead {o.dead_fraction:.2f} |"
            f" k_E {o.k_earth:5.1f} (v_e {o.v_e_earth / 1e3:5.1f}, MR {o.mass_ratio_earth:.2f},"
            f" D {o.departing_per_earth_impactor:5.1f}) k_S {o.k_sun:5.1f} (MR"
            f" {o.mass_ratio_sun:.2f}) | G {o.growth_per_total_impactor:5.2f} per 3S"
            f" (Sun wave {o.sun_wave_share:.0%}), doubling {o.doubling_years:.2f} yr"
            + (f"  [{o.at_edge} at its edge]" if o.at_edge else "")
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT_STAGE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stage[0]))
        w.writeheader()
        w.writerows(stage)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(out[0])))
        w.writeheader()
        w.writerows(asdict(o) for o in out)
    print(f"python: wrote {len(stage)} rows -> {OUT_STAGE}, {len(out)} rows -> {OUT}")


if __name__ == "__main__":
    main()
