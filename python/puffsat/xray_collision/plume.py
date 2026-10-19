"""Result 18: the push plate's backward plume, and what it does to the next sheet (parent S36).

    PYTHONPATH=python uv run python -m puffsat.xray_collision.plume

On the 2 Hz, 1 g overtake push each ~100 km/s sheet catches up from behind and flies face-on
through the gas of the pulses before it, which the plate threw back up the stream. The parent
models each earlier pulse as a Gaussian cloud receding at U and spreading at u along any axis
(`coated_sheet.flux`); u is a free input and sets most of the dose. This module replaces the
Gaussian with the 2-D exhaust (`sweep --plume-2d`).

1. **The exhaust.** Each run tallies the gas leaving the grid (and what is left on it at the end)
   by angle from the stream axis and terminal speed `sqrt(u^2 + 2h)`, the enthalpy spent along the
   parcel's present direction. The grid reaches 50 m above the plate so the rebound has finished
   expanding (doses averaged over >= 3 deg move < 3% from 50 to 100 m on the coarse grid).
2. **The dose.** Beyond ~1 km the gas is ballistic (its mean free path passes metres there), so
   each earlier pulse is a point source: the density on the axis at distance `z` and age `t` is
   `f(v) / (v^2 t^3)`, with `v = z / t` and `f = dm / (dOmega dv)` on the axis. The next sheet
   subtends only ~0.02 deg at the ~25 km where the dose arrives, so the axis value is what counts.
   An axisymmetric grid resolves it worst there, so `f` is averaged over cones of 1, 3, 5 and
   10 deg and all four are reported. The flux `1/2 rho (v + w)^3` in the sheet's frame is summed
   over earlier pulses exactly as `coated_sheet.flux` sums the Gaussians.
3. **The equivalent u.** The u at which the parent's Gaussian (with its momentum-set U) gives the
   same dose, and Kapton's face on the computed flux (Result 13's rule).
4. **Species.** From Result 11's freeze-out (`freeze.py`): what is still ionized or dissociated,
   and the recombination energy a Kapton face could collect, against the kinetic dose.
5. **Cones.** Convex cones tilted 5-20 deg, apex toward the stream: on-axis dose and impulse.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from puffsat.xray_collision import band_bowl, charge, coated_sheet, fast_wave, freeze

JSONL = charge.OUT_DIR / "plume_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "plume_2d_coarse.jsonl"
CLOSURE = fast_wave.CLOSURES[0]  # 9.86 R_sun headline: 99.2 km/s, k = 7.07 (the runs' k)
#: Every printed closure, each on its own speed, k and v_e. The runs are at k = 7.07; a gamma-law
#: gas has no speed scale, so the exhaust's shape in units of `w` carries to any speed, and to
#: k = 6.2-7.45 to within the weak k-dependence of the sheet runs (Results 15-16).
CLOSURES = (
    CLOSURE,
    band_bowl.BAND[0],
    band_bowl.BAND[1],
    band_bowl.BAND[2],
    fast_wave.CLOSURES[1],  # 4 R_sun headline, 137.1 km/s (flies dry carbon fabric)
    fast_wave.CLOSURES[2],  # 9.86 R_sun, 90: 132.1 km/s
    fast_wave.CLOSURES[3],  # 4 R_sun, 150: 207.8 km/s
)
#: The cases carried to every closure: the flown sheet and the floor-wide sheet of ADR-0068.
CLOSURE_CASES = ("flat + skirt, 8 m sheet", "flat + skirt, 9.5 m sheet")
CAPS = (1.0, 3.0, 5.0, 10.0)  # [deg], the cones the on-axis density is averaged over
DESIGN_CAP = 3.0  # the narrowest cone that has converged in the grid's height
KAPTON = 100e-6  # [m], the flown sheet
#: Speed and angle bins for item 1's table (in the closing speed, and degrees from the axis).
TABLE_SPEEDS = (0.0, 0.2, 0.3, 0.4, 0.5, 2.0)
TABLE_ANGLES = (0.0, 5.0, 10.0, 20.0, 40.0, 90.0, 180.0)

OUT_DOSE = charge.OUT_DIR / "plume_dose.csv"
OUT_TABLE = charge.OUT_DIR / "plume_table.csv"
OUT_SPECIES = charge.OUT_DIR / "plume_species.csv"
OUT_SHEET = charge.OUT_DIR / "plume_sheet.csv"
OUT_PUSH = charge.OUT_DIR / "plume_push.csv"


@cache
def runs(grid: str = "x1") -> tuple[dict[str, object], ...]:
    path = JSONL if grid == "x1" else JSONL_COARSE
    return tuple(charge.load_jsonl(path)) if path.exists() else ()


def _hist(row: dict[str, object]) -> tuple[list[float], NDArray[np.float64], float]:
    angles = row["angles"]
    assert isinstance(angles, list)
    h = np.array(row["exhaust"], dtype=np.float64)
    return [float(a) for a in angles], h, float(str(row["speed_bin"]))


def on_axis(row: dict[str, object], cap: float, w: float) -> tuple[NDArray[np.float64], ...]:
    """`(v [m/s], f)`: the exhaust's `dm / (dOmega dv)` per kg of merged gas, averaged over the
    cone `cap` [deg] about the axis, by terminal speed."""
    angles, h, dv = _hist(row)
    i = angles.index(cap)
    omega = 2.0 * math.pi * (1.0 - math.cos(math.radians(cap)))
    v = np.asarray((np.arange(h.shape[1]) + 0.5) * dv * w, dtype=np.float64)
    f = np.asarray(h[:i].sum(axis=0) / (omega * dv * w), dtype=np.float64)
    return v, f


def point_flux(
    v: NDArray[np.float64],
    f: NDArray[np.float64],
    c: fast_wave.Closure = CLOSURE,
    power: int = 3,
) -> NDArray[np.float64]:
    """Incident kinetic flux [W/m^2] on the sheet, on `coated_sheet`'s time grid, when every
    earlier pulse is a ballistic point source whose on-axis `dm / (dOmega dv)` per kg of cloud is
    `f(v)` (plate frame, `v` [m/s] away from the plate). `power = 2` gives the momentum flux
    `rho (v + w)^2` [Pa] instead: the push that slows the sheet."""
    p = fast_wave.push_at(c)
    before = coated_sheet._times()
    q = np.zeros_like(before)
    tau = 1.0 / p.cadence
    for k in range(1, int(coated_sheet.WINDOW * p.cadence) + 40):
        t = k * tau - before
        ok = (t > 0.0) & (before > 0.0)
        tt = t[ok]
        vz = c.w * before[ok] / tt  # the gas speed that is at the sheet now
        rho = p.cloud * np.interp(vz, v, f, right=0.0) / (vz**2 * tt**3)
        q[ok] += (0.5 if power == 3 else 1.0) * rho * (vz + c.w) ** power
    return q


def flux(row: dict[str, object], cap: float, c: fast_wave.Closure = CLOSURE) -> NDArray[np.float64]:
    """The 2-D exhaust's flux on the sheet, its on-axis density averaged over the cone `cap`."""
    v, f = on_axis(row, cap, c.w)
    return point_flux(v, f, c)


def dose(q: NDArray[np.float64]) -> float:
    return float(np.sum(0.5 * (q[1:] + q[:-1]) * coated_sheet.DT))


def equivalent_u(d: float, c: fast_wave.Closure = CLOSURE) -> float:
    """The Gaussian spread [m/s] whose dose is `d`, with the parent's momentum-set U."""
    p = fast_wave.push_at(c)

    def gap(u: float) -> float:
        return dose(coated_sheet.flux(p, u, coated_sheet._times())) - d

    return float(brentq(gap, 1.0e3, 60.0e3, xtol=10.0))


def kapton_floor(q: NDArray[np.float64]) -> float:
    """The thinnest Kapton [m] whose face holds Result 13's rule on `q`; NaN if 2 mm fails."""

    def holds(th: float) -> bool:
        h, _ = coated_sheet.heat(q, th, fast_wave.KAPTON_EPS, fast_wave.KAPTON_EPS, 0.0)
        return h.rating_before <= coated_sheet.LAST

    lo, hi = 10.0e-6, 2.0e-3
    if not holds(hi):
        return math.nan
    if holds(lo):
        return lo
    while hi / lo > 1.02:
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if holds(mid) else (mid, hi)
    return hi


@dataclass(frozen=True)
class DoseRow:
    closure: str
    case: str
    grid: str
    tilt_deg: float
    impulse_per_mv: float
    impulse_vs_flat: float  # against the flat skirted plate under the same 9.5 m sheet
    within_3deg: float  # mass share of the merged gas
    within_5deg: float
    cap_deg: float
    dose_kj: float
    u_equivalent: float  # [km/s]
    face_100um_k: float
    holds_100um: bool
    kapton_floor_um: float


def dose_rows() -> list[DoseRow]:
    out = []
    for grid in ("x1", "x0.5"):
        rows = runs(grid)
        if not rows:
            continue
        ref = next(
            float(str(r["impulse_per_mv"]))
            for r in rows
            if r["case"] == "flat + skirt, 9.5 m sheet"
        )
        cases = [(CLOSURE, r, cap) for r in rows for cap in CAPS]
        cases += [
            (c, r, DESIGN_CAP)
            for c in CLOSURES[1:]
            for r in rows
            if grid == "x1" and r["case"] in CLOSURE_CASES
        ]
        for c, r, cap in cases:
            angles, h, _ = _hist(r)
            share = {a: float(h[: angles.index(a)].sum()) for a in (3.0, 5.0)}
            q = flux(r, cap, c)
            d = dose(q)
            heat, _ = coated_sheet.heat(q, KAPTON, fast_wave.KAPTON_EPS, fast_wave.KAPTON_EPS, 0.0)
            out.append(
                DoseRow(
                    c.name,
                    str(r["case"]),
                    grid,
                    float(str(r["tilt_deg"])),
                    float(str(r["impulse_per_mv"])),
                    float(str(r["impulse_per_mv"])) / ref,
                    share[3.0],
                    share[5.0],
                    cap,
                    d / 1e3,
                    equivalent_u(d, c) / 1e3,
                    heat.face_max,
                    heat.rating_before <= coated_sheet.LAST,
                    1e6 * kapton_floor(q),
                )
            )
    return out


@dataclass(frozen=True)
class TableRow:
    case: str
    angle_lo: float
    angle_hi: float
    shares: tuple[float, ...]  # by TABLE_SPEEDS bin, of the merged gas


def table_rows() -> list[TableRow]:
    """Item 1: the exhaust's mass by angle from the axis and terminal speed (production grid)."""
    out = []
    for r in runs():
        angles, h, dv = _hist(r)
        speeds = (np.arange(h.shape[1]) + 0.5) * dv
        for a_lo, a_hi in pairwise(TABLE_ANGLES):
            rows = [i for i in range(len(angles) - 1) if a_lo <= angles[i] < a_hi]
            sub = h[rows].sum(axis=0)
            shares = tuple(
                float(sub[(speeds >= s_lo) & (speeds < s_hi)].sum())
                for s_lo, s_hi in pairwise(TABLE_SPEEDS)
            )
            out.append(TableRow(str(r["case"]), a_lo, a_hi, shares))
    return out


@dataclass(frozen=True)
class SpeciesRow:
    ion_share_frozen: float  # of the heat, still ionized at the plate
    h2_margin_decades: float  # H + H + H -> H2 against the expansion
    kt_at_plate_ev: float
    carbon_atomic_mj_per_kg: float  # the bond energy carbon keeps if it stays atomic
    kinetic_mj_per_kg: float  # what the gas brings the sheet, at its mean speed (sheet frame)
    recombination_share: float  # of the kinetic dose, if every carbon atom recombines on the face


def species_row(c: fast_wave.Closure = CLOSURE) -> SpeciesRow:
    """Item 3 on Result 11's freeze-out, and the exhaust's dose-weighted speed."""
    fr = freeze.freeze(c.k, charge.ch2_table(), v_rel=c.w)
    per_ch2 = (freeze.CH2_ATOMIZATION - freeze.H2_BOND) / (14.027 * freeze.AMU)
    # The flown plate's exhaust within the design cone, weighted as the dose weights it.
    row = next(r for r in runs() if r["case"] == "flat + skirt, 8 m sheet")
    v, f = on_axis(row, DESIGN_CAP, c.w)
    weight = f * (v + c.w) ** 3
    rel = float(np.sum(weight * (v + c.w)) / np.sum(weight))
    kinetic = 0.5 * rel**2
    return SpeciesRow(
        fr.ion_share_frozen,
        fr.h2_margin_decades,
        fr.kt_plate_ev,
        per_ch2 / 1e6,
        kinetic / 1e6,
        per_ch2 / kinetic,
    )


def _write(path: Path, rows: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


#: A point mass on the film (the guidance node): its mass, the radius over which it is attached,
#: and the dimple Result 11's flatness tolerance allows (10 cm spreads the arrival ~1 us).
NODE = 0.1  # [kg], a placeholder, not a design value
NODE_RADIUS = 0.05  # [m]
FLATNESS = 0.10  # [m]


@dataclass(frozen=True)
class PushRow:
    """How far the plume slows the floor-wide sheet: it pushes straight back along the flight line,
    and evenly across the sheet (which subtends ~0.02 deg of it), so it delays rather than steers.
    A point mass `m` on the film (a guidance node) is not slowed with it, and the film's tension
    `T` must drag it along: the film dimples into a cone of depth `(m a / 2 pi T) ln(R / r0)`
    under the peak deceleration `a`. The tension that holds that to Result 11's 10 cm is given
    per NODE kg of attached mass, with the film's spin stress and the spin that sets it."""

    closure: str
    momentum_ns_m2: float  # [N s/m^2] taken over the approach
    peak_pa: float
    dv: float  # [m/s] lost by the film
    lag_m: float  # [m] the film arrives short of where it would have been
    delay_us: float  # its arrival delay
    peak_decel: float  # [m/s^2] of the film
    tension_n_m: float  # [N/m] holding a NODE-kg point mass to a 10 cm dimple
    film_stress_mpa: float  # that tension in the film
    spin_rpm: float  # the spin whose hoop stress at the centre gives it


def push_rows() -> list[PushRow]:
    row = next(r for r in runs() if r["case"] == CLOSURE_CASES[1])
    rs = float(str(row["sheet_radius"]))
    out = []
    for c in CLOSURES:
        sigma = fast_wave.PULSE / (c.k * c.v_e) / (math.pi * rs**2)
        v, f = on_axis(row, DESIGN_CAP, c.w)
        p = point_flux(v, f, c, power=2)
        # The time grid runs from WINDOW s before impact to 0; integrate forward in time.
        pt = p[::-1]
        dt = coated_sheet.DT
        impulse = np.concatenate(([0.0], np.cumsum(0.5 * (pt[1:] + pt[:-1]) * dt)))
        dv = impulse / sigma
        lag = float(np.sum(0.5 * (dv[1:] + dv[:-1]) * dt))
        force = NODE * float(p.max()) / sigma
        tension = force * math.log(rs / NODE_RADIUS) / (2.0 * math.pi * FLATNESS)
        th = sigma / coated_sheet.KAPTON_RHO
        # A spinning uniform disk's hoop stress at its centre: (3 + nu) / 8 rho omega^2 R^2.
        omega = math.sqrt(tension / th / ((3.0 + 0.34) / 8.0 * coated_sheet.KAPTON_RHO * rs**2))
        out.append(
            PushRow(
                c.name,
                float(impulse[-1]),
                float(p.max()),
                float(dv[-1]),
                lag,
                1e6 * lag / c.w,
                float(p.max()) / sigma,
                tension,
                tension / th / 1e6,
                60.0 / (2.0 * math.pi) * omega,
            )
        )
    return out


@dataclass(frozen=True)
class SheetRow:
    """Kapton spread over each sheet at each closure: its own heat in this plume, its face's
    pressure ceiling at the flown 2 m standoff, and the impulse it gives."""

    closure: str
    case: str
    sheet_radius: float  # [m]
    kapton_um: float  # the pulse's impactor mass spread over the sheet
    floor_um: float  # thinnest Kapton that holds in this plume (3 deg cone)
    peak_coeff: float  # P / (sigma w^2) on the face, this run
    peak_gpa: float  # at this thickness
    ceiling_hy100_um: float  # thickest the HY-100 / 4340 face takes (1.5 GPa)
    ceiling_maraging_um: float  # maraging 300 (2.5 GPa)
    holds_hy100: bool  # floor <= thickness <= ceiling
    holds_maraging: bool
    impulse_per_mv: float  # lossless gamma-gas J/(m w)
    rebound_e: float  # polyethylene's e_eff at the arriving column (charge.rebound)
    impulse_vs_flown: float  # J/(m w) (1 + e) against the flown 8 m sheet's


def sheet_rows() -> list[SheetRow]:
    """Every closure on the flown and the floor-wide sheet (production grid)."""
    rows = {str(r["case"]): r for r in runs()}
    out = []
    for c in CLOSURES:
        mass = fast_wave.PULSE / (c.k * c.v_e)
        fire = charge.charge_fireball(c.k, table=charge.ch2_table(), v_rel=c.w)
        speed = fire.u_rms + fire.drift

        def impulse(
            r: dict[str, object],
            rs: float,
            mass: float = mass,
            k: float = c.k,
            speed: float = speed,
        ) -> tuple[float, float]:
            column = (1.0 + k) * mass / (math.pi * rs**2)
            e = charge.rebound(speed, column)
            return float(str(r["impulse_per_mv"])) * (1.0 + e), e

        flown = rows[CLOSURE_CASES[0]]
        ref, _ = impulse(flown, float(str(flown["sheet_radius"])))
        for case in CLOSURE_CASES:
            r = rows[case]
            rs = float(str(r["sheet_radius"]))
            sigma = mass / (math.pi * rs**2)
            coeff = float(str(r["peak_per_mv2"])) * math.pi * rs**2
            floor = kapton_floor(flux(r, DESIGN_CAP, c))
            th = sigma / coated_sheet.KAPTON_RHO
            ceil = {
                face: allow / (coeff * c.w**2) / coated_sheet.KAPTON_RHO
                for face, allow in fast_wave.FACES.items()
            }
            j, e = impulse(r, rs)
            out.append(
                SheetRow(
                    c.name,
                    case,
                    rs,
                    1e6 * th,
                    1e6 * floor,
                    coeff,
                    coeff * sigma * c.w**2 / 1e9,
                    1e6 * ceil["HY-100 / 4340"],
                    1e6 * ceil["maraging 300"],
                    bool(floor <= th <= ceil["HY-100 / 4340"]),
                    bool(floor <= th <= ceil["maraging 300"]),
                    float(str(r["impulse_per_mv"])),
                    e,
                    j / ref,
                )
            )
    return out


def main() -> None:
    p = fast_wave.push_at(CLOSURE)
    gauss = [
        dose(coated_sheet.flux(p, u, coated_sheet._times())) / 1e3 for u in coated_sheet.SPREADS
    ]
    print(
        f"  Gaussian (U {p.recession / 1e3:.1f} km/s): dose "
        + " / ".join(f"{d:.0f}" for d in gauss)
        + " kJ/m^2 at u = 5 / 8 / 12"
    )
    rows = dose_rows()
    for d in rows:
        print(
            f"  {d.grid:<4} {d.closure:<26} {d.case:<26} J/mv {d.impulse_per_mv:.3f}"
            f" (x{d.impulse_vs_flat:.3f})"
            f" | <3 deg {d.within_3deg:.3f} | cap {d.cap_deg:4.1f}: dose {d.dose_kj:5.0f} kJ/m^2,"
            f" u_eq {d.u_equivalent:5.1f} km/s, 100 um face {d.face_100um_k:5.0f} K"
            f" ({'holds' if d.holds_100um else 'fails'}), floor {d.kapton_floor_um:5.0f} um"
        )
    table = table_rows()
    for t in table:
        print(
            f"  table {t.case:<26} {t.angle_lo:5.1f}-{t.angle_hi:5.1f} deg: "
            + " ".join(f"{s:.3f}" for s in t.shares)
        )
    sp = species_row()
    print(
        f"  species: ions frozen {sp.ion_share_frozen:.1e} of the heat, H2 margin"
        f" {sp.h2_margin_decades:.1f} decades, kT {sp.kt_at_plate_ev:.2f} eV at the plate; atomic"
        f" carbon keeps {sp.carbon_atomic_mj_per_kg:.0f} MJ/kg against {sp.kinetic_mj_per_kg:.0f}"
        f" MJ/kg kinetic -> {sp.recombination_share:.2%} of the dose"
    )
    _write(OUT_DOSE, list(rows))
    _write(OUT_TABLE, list(table))
    _write(OUT_SPECIES, [sp])
    sheets = sheet_rows()
    for sr in sheets:
        print(
            f"  sheet {sr.closure:<26} {sr.case:<26} {sr.kapton_um:4.0f} um (floor"
            f" {sr.floor_um:5.0f}) | {sr.peak_gpa:4.2f} GPa, ceilings HY-100"
            f" {sr.ceiling_hy100_um:4.0f} / maraging {sr.ceiling_maraging_um:4.0f} um ->"
            f" {'holds' if sr.holds_hy100 else 'fails'} /"
            f" {'holds' if sr.holds_maraging else 'fails'}"
            f" | J/mw {sr.impulse_per_mv:.3f}, e {sr.rebound_e:.3f}, x{sr.impulse_vs_flown:.4f}"
        )
    _write(OUT_SHEET, list(sheets))
    pushes = push_rows()
    for pr in pushes:
        print(
            f"  push {pr.closure:<26} {pr.momentum_ns_m2:.3f} N s/m^2 (peak {pr.peak_pa:5.1f} Pa)"
            f" -> film slows {pr.dv:5.2f} m/s, arrives {pr.lag_m:.3f} m short"
            f" ({pr.delay_us:.1f} us late) | peak {pr.peak_decel:4.0f} m/s^2: a {NODE * 1e3:.0f} g"
            f" node needs {pr.tension_n_m:4.0f} N/m ({pr.film_stress_mpa:.2f} MPa,"
            f" {pr.spin_rpm:4.1f} rpm) for a {100 * FLATNESS:.0f} cm dimple"
        )
    _write(OUT_PUSH, list(pushes))
    print(f"python: wrote {OUT_DOSE.parent}/plume_*.csv")


if __name__ == "__main__":
    main()
