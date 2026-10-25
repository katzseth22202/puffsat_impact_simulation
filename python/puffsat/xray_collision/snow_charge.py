"""Result 22: water or argon snow as the charge, on all three plates (parent S60).

    PYTHONPATH=python uv run python -m puffsat.xray_collision.snow_charge

Seth (2026-10-24) asks what happens if the polyethylene charge becomes water snow, a water-and-
bubbles froth, liquid water, ice, or argon snow, in a thin bag, on the growth push (a Kapton sheet
at 99 km/s, k = 6.4), the departure bowl (a 15 kg impactor at 84 km/s, k = 12.4) and the Sun bowl
(a 1 kg disk at 392 km/s, k = 20). The parent expects two effects of density to pull opposite
ways: pores heat the charge more, and a lighter charge is a longer column.

The pieces:

1. **The density law** (`density_rows`, S60 item 1). `sweep --snow-charge` runs the real disk-
   into-charge collision on the δ/D 0.5 paraboloid with the charge at 1.5 F, its column
   `k m / (rho pi R^2)` long behind the impactor's footprint, at 30-1620 kg/m^3. The gamma gas has
   no material, so bulk density enters only as the column's length. The push sheet's thickness
   is swept on the grid instead (0.1-0.8 m against the real 0.5-2.5 mm).
2. **The material** (`material_rows`, items 2 and 3). Polyethylene, water and argon raced on their
   own equilibrium tables at each plate's heat per kilogram (`iapetus_bag.race`, `freeze.freeze`):
   the effective gamma, the ionization and bonds that freeze, and the plate rebound. The gamma
   enters through each plate's own gamma bracket of 2-D runs.
3. **Porosity** (`hugoniot_rows`, item 2). A 1-D bound. In the snowplow the disk strikes each
   kilogram at its own current speed `V = w / (1 + m)`, so the slowest charge is the last, struck
   at `w / (1 + k)`. Its shock heat `V^2 / 2` is exact (conservation, any EOS). Pore collapse can
   add at most the solid's cold-compression share of it, bounded here with linear Hugoniots and a
   Mie-Gruneisen release at that slowest speed only (the hotter states are past the fits).
4. **The graded charge** (`graded_rows`, item 7). Ice at 917 kg/m^3 under snow at 1/3 or 1/5 of
   it, the dense part 25-75% of k, at either end.
5. **The bag** (`bag_rows`, item 4): the film's share of the charge and the gas it holds.
6. **Argon's flash and the debris** (`argon_flash`, `debris_rows`, items 3 and 6): what a sealed
   argon bag holds as vapour at 60 K, and how long lost ice, snow and argon pieces last.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np

from puffsat import eos_water
from puffsat.xray_collision import charge, freeze, iapetus_bag

OUT_DIR = charge.OUT_DIR
JSONL = OUT_DIR / "snow_charge_2d.jsonl"
JSONL_COARSE = OUT_DIR / "snow_charge_2d_coarse.jsonl"
REBOUND_FAST = OUT_DIR / "iapetus_rebound_fast.jsonl"
REBOUND_SLOW = OUT_DIR / "iapetus_rebound_slow.jsonl"

SIGMA = 5.670374419e-8
R_GAS = 8.314462618
AMU = 1.66053906660e-27
K_BOLTZ = 1.380649e-23
SOLAR_CONSTANT = 1361.0  # [W/m^2] at 1 AU
AU_IN_RSUN = 215.032

PE_SOLID = 950.0
#: The bowls' production cell [m] (`sweep --snow-charge`: half the Sun plate's 0.05 m).
BOWL_DR = 0.025
PE_GAMMA = 1.4  # the 2-D runs' gas, as Results 10-15 price polyethylene


@dataclass(frozen=True)
class Plate:
    name: str
    impactor: float  # [kg]
    w: float  # [m/s] closing speed
    k: float
    disk_radius: float  # [m]
    plate_distance: float  # [m], charge to plate (the race's stopping edge)
    plate_radius: float  # [m], for the arriving column


#: aim `7ce7628`'s floored headline at 9.86 R_sun, X 60 (the S60 table). The push's impactor is
#: the parent's ~200 kg charge over k; its sheet is 0.07 kg/m^2, so its charge 0.45 kg/m^2.
PLATES = (
    Plate("push", 200.0 / 6.4, 99.0e3, 6.4, 5.0, 2.0, 9.5),
    Plate("departure", 15.0, 84.0e3, 12.4, 0.4, 3.75, 10.0),
    Plate("sun", 1.0, 392.0e3, 20.0, 0.25, 3.75, 10.0),
)
PUSH_SHEET = 0.07  # [kg/m^2], the Kapton impactor


def plate(name: str) -> Plate:
    return next(p for p in PLATES if p.name == name)


@dataclass(frozen=True)
class Charge:
    name: str
    material: str  # "pe", "ice" (water in any phase) or "argon"
    rho: float  # bulk density [kg/m^3]
    porous: bool


#: S60 item 1's densities, polyethylene foam as the material control. Argon snow's bulk density is
#: not given by the parent; it is run at water snow's two densities.
CHARGES = (
    Charge("solid polyethylene", "pe", PE_SOLID, False),
    Charge("polyethylene foam 600", "pe", 600.0, True),
    Charge("polyethylene foam 500", "pe", 500.0, True),
    Charge("polyethylene foam 300", "pe", 300.0, True),
    Charge("liquid water", "ice", 1000.0, False),
    Charge("water ice", "ice", 917.0, False),
    Charge("water snow 600", "ice", 600.0, True),
    Charge("water froth 500", "ice", 500.0, True),
    Charge("water snow 300", "ice", 300.0, True),
    Charge("solid argon", "argon", 1620.0, False),
    Charge("argon snow 600", "argon", 600.0, True),
    Charge("argon snow 300", "argon", 300.0, True),
)


# --- 1. The density law ------------------------------------------------------------------------


def _f(r: dict[str, object], key: str) -> float:
    return float(str(r[key]))


@cache
def runs(path: Path = JSONL) -> tuple[dict[str, object], ...]:
    """The `--snow-charge` rows: the whole file, with any plate run apart (`--snow-charge push`,
    its own file) replacing that plate's rows."""
    coarse = "_coarse" if path.name.endswith("_coarse.jsonl") else ""
    stem = path.name.removesuffix(".jsonl").removesuffix("_coarse")
    rows = list(charge.load_jsonl(path)) if path.exists() else []
    for name in ("sun", "departure", "push"):
        p = path.with_name(f"{stem}_{name}{coarse}.jsonl")
        if p.exists():
            rows = [r for r in rows if r["plate"] != name] + list(charge.load_jsonl(p))
    return tuple(rows)


def _uniform(plate_name: str, path: Path = JSONL, gamma: float = 1.4) -> list[tuple[float, float]]:
    """`(rho, J/mv)` of the plate's uniform columns at `gamma`, by density."""
    return sorted(
        (_f(r, "rho"), _f(r, "impulse_per_mv"))
        for r in runs(path)
        if r["plate"] == plate_name
        and r["kind"] == "uniform"
        and abs(_f(r, "gamma") - gamma) < 1e-3
    )


def density_factor(plate_name: str, rho: float, path: Path = JSONL) -> float:
    """`J(rho) / J(950)` on a bowl, interpolated in ln rho (held at the run's edges). The push sheet
    is far thinner than any grid sheet at every density, and its thickness sweep is flat
    (`push_thickness_spread`), so its factor is 1."""
    if plate_name == "push":
        return 1.0
    p = plate(plate_name)
    if rho > PE_SOLID and p.k * p.impactor / (rho * math.pi * p.disk_radius**2) < 4.0 * BOWL_DR:
        # Denser than the reference and under four production cells long (solid argon on the Sun
        # bowl, 2.4 cells): the run reads the grid, not the charge, so it is held at solid's.
        return 1.0
    pts = _uniform(plate_name, path)
    lr = [math.log(a) for a, _ in pts]
    js = [b for _, b in pts]
    ref = float(np.interp(math.log(PE_SOLID), lr, js))
    return float(np.interp(math.log(rho), lr, js)) / ref


def push_sheets(path: Path = JSONL, gamma: float = 1.4) -> list[tuple[float, float]]:
    """`(thickness, J/mv)` of the push's uniform sheets."""
    return sorted(
        (_f(r, "length"), _f(r, "impulse_per_mv"))
        for r in runs(path)
        if r["plate"] == "push"
        and r["kind"] == "sheet"
        and _f(r, "dense_share") == 0.0
        and abs(_f(r, "gamma") - gamma) < 1e-3
    )


def push_thickness_spread(path: Path = JSONL) -> float:
    """Largest departure of a sheet's impulse from the 0.2 m production sheet's."""
    pts = push_sheets(path)
    ref = next(j for t, j in pts if abs(t - 0.2) < 1e-9)
    return max(abs(j / ref - 1.0) for _, j in pts)


def point_reference(plate_name: str, path: Path = JSONL) -> float:
    """The point fireball's `J/mv` at gamma 1.4."""
    return next(
        _f(r, "impulse_per_mv")
        for r in runs(path)
        if r["plate"] == plate_name and r["kind"] == "point" and abs(_f(r, "gamma") - 1.4) < 1e-3
    )


@dataclass(frozen=True)
class DensityRow:
    plate: str
    rho: float
    length: float  # [m]; the sheet's grid thickness on the push
    length_over_radius: float
    impulse_per_mv: float
    factor: float  # over solid polyethylene's 950 (the push: over the 0.2 m sheet)
    coarse_factor: float  # the same on the grid twice as coarse
    real_over_point: float
    peak_per_mv2: float
    peak_over_solid: float  # peak face pressure over solid's (not converged: indicative only)


def density_rows(path: Path = JSONL, coarse: Path = JSONL_COARSE) -> list[DensityRow]:
    out = []
    for p in PLATES:
        if p.name == "push":
            pts = push_sheets(path)
            ref = next(j for t, j in pts if abs(t - 0.2) < 1e-9)
            push_ref_peak = next(
                _f(r, "peak_per_mv2")
                for r in runs(path)
                if r["plate"] == "push"
                and r["kind"] == "sheet"
                and abs(_f(r, "length") - 0.2) < 1e-9
                and _f(r, "dense_share") == 0.0
                and abs(_f(r, "gamma") - 1.4) < 1e-3
            )
            cpts = dict(push_sheets(coarse)) if runs(coarse) else {}
            cref = cpts.get(0.2, math.nan)
            for t, j in pts:
                rec = next(
                    r
                    for r in runs(path)
                    if r["plate"] == "push"
                    and r["kind"] == "sheet"
                    and _f(r, "length") == t
                    and _f(r, "dense_share") == 0.0
                    and abs(_f(r, "gamma") - 1.4) < 1e-3
                )
                out.append(
                    DensityRow(
                        "push",
                        math.nan,
                        t,
                        t / p.disk_radius,
                        j,
                        j / ref,
                        cpts.get(t, math.nan) / cref,
                        math.nan,
                        _f(rec, "peak_per_mv2"),
                        _f(rec, "peak_per_mv2") / push_ref_peak,
                    )
                )
            continue
        pts = _uniform(p.name, path)
        cpts = dict(_uniform(p.name, coarse)) if runs(coarse) else {}
        point = point_reference(p.name, path)
        for rho, j in pts:
            rec = next(
                r
                for r in runs(path)
                if r["plate"] == p.name
                and r["kind"] == "uniform"
                and _f(r, "rho") == rho
                and abs(_f(r, "gamma") - 1.4) < 1e-3
            )
            ln = _f(rec, "length")
            cf = cpts[rho] / cpts[PE_SOLID] if rho in cpts and PE_SOLID in cpts else math.nan
            out.append(
                DensityRow(
                    p.name,
                    rho,
                    ln,
                    ln / p.disk_radius,
                    j,
                    density_factor(p.name, rho, path),
                    cf,
                    j / point,
                    _f(rec, "peak_per_mv2"),
                    _f(rec, "peak_per_mv2") / _solid_peak(p.name, path),
                )
            )
    return out


def breakpoint(plate_name: str, loss: float = 0.02, path: Path = JSONL) -> float:
    """The lowest bulk density [kg/m^3] whose impulse stays within `loss` of solid's (log-linear
    between the runs); the lowest run if every run holds."""
    pts = _uniform(plate_name, path)
    rhos = sorted({a for a, _ in pts}, reverse=True)
    prev = PE_SOLID
    for rho in rhos:
        if rho > PE_SOLID:
            continue
        if density_factor(plate_name, rho, path) < 1.0 - loss:
            # Bisect in ln rho between the last holding run and this one.
            lo, hi = math.log(rho), math.log(prev)
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                if density_factor(plate_name, math.exp(mid), path) < 1.0 - loss:
                    lo = mid
                else:
                    hi = mid
            return math.exp(hi)
        prev = rho
    return rhos[-1]


# --- 2. The material ----------------------------------------------------------------------------


CONFINED = OUT_DIR / "snow_charge_2d_confined.jsonl"


@cache
def confined_ratio(gamma: float) -> float:
    """The confined plane wave's `1 + e` in a gamma gas (`sweep --snow-charge confined`)."""
    rows = charge.load_jsonl(CONFINED)
    return next(_f(r, "confined_ratio") for r in rows if abs(_f(r, "gamma") - gamma) < 1e-3)


def gamma_factor(plate_name: str, gamma: float, path: Path = JSONL) -> float:
    """`(J/cr)(gamma) / (J/cr)(1.4)` of the plate's reference charge (the 950 kg/m^3 column on a
    bowl, the 0.2 m sheet on the push), as Result 20's `gamma_factor`: each run's impulse over its
    own gas's confined bounce `cr`, so the factor carries the approach and the material's
    `1 + e_eff` carries the bounce. Linear in gamma between the bracket's runs and carried on the
    end segment's slope outside them (`gamma_in_bracket` says when)."""
    if plate_name == "push":
        pts = sorted(
            (_f(r, "gamma"), _f(r, "impulse_per_mv"))
            for r in runs(path)
            if r["plate"] == "push"
            and r["kind"] == "sheet"
            and abs(_f(r, "length") - 0.2) < 1e-9
            and _f(r, "dense_share") == 0.0
        )
    else:
        pts = sorted(
            (_f(r, "gamma"), _f(r, "impulse_per_mv"))
            for r in runs(path)
            if r["plate"] == plate_name and r["kind"] == "uniform" and _f(r, "rho") == PE_SOLID
        )
    gs = [a for a, _ in pts]
    # Each run's impulse over its own gas's bounce: the material's `1 + e_eff` replaces the bounce.
    js = [b / confined_ratio(a) for a, b in pts]
    i = 0 if gamma < gs[1] else len(gs) - 2
    if gs[0] <= gamma <= gs[-1]:
        j = float(np.interp(gamma, gs, js))
    else:
        j = js[i] + (js[i + 1] - js[i]) * (gamma - gs[i]) / (gs[i + 1] - gs[i])
    return j / float(np.interp(1.4, gs, js))


def gamma_in_bracket(gamma: float) -> bool:
    return 1.2 <= gamma <= 5.0 / 3.0


@cache
def fast_rebound_rows() -> tuple[dict[str, object], ...]:
    """This result's water and argon rebound runs: 50-100 km/s (`--iapetus-rebound fast`) and
    6-12 km/s (`--iapetus-rebound slow`), beside Result 20's 15-35 km/s."""
    rows: list[dict[str, object]] = []
    for p in (REBOUND_FAST, REBOUND_SLOW):
        if p.exists():
            rows += charge.load_jsonl(p)
    return tuple(r for r in rows if r["converged"] and _f(r, "length") == 1.0)


def rebound(material: str, v: float, column: float) -> float:
    """Plate rebound `e_eff` at arrival speed `v` and column [kg/m^2]: polyethylene from Results
    9-10's sweeps; water and argon from Result 20's 15-35 km/s grid joined to this result's
    6-12 and 50-100 km/s runs (`fast_rebound_rows`), interpolated in v and log column."""
    if material == "pe":
        return charge.rebound(v, column)
    name = "argon" if material == "argon" else "water"
    rows = [r for r in (*iapetus_bag.rebound_rows(), *fast_rebound_rows()) if r["material"] == name]
    vs = sorted({_f(r, "v") for r in rows})
    cols = sorted({_f(r, "rho_impact") for r in rows})
    lc = np.log(cols)
    x = float(np.clip(math.log(column), lc[0], lc[-1]))
    per_v = []
    for v_i in vs:
        pts = sorted(
            (math.log(_f(r, "rho_impact")), _f(r, "e_eff")) for r in rows if _f(r, "v") == v_i
        )
        per_v.append(float(np.interp(x, [a for a, _ in pts], [b for _, b in pts])))
    return float(np.interp(v, vs, per_v))


def rebound_speed_covered(material: str, v: float) -> bool:
    """Whether `v` lies inside the material's rebound grid (otherwise it is held at the edge)."""
    if material == "pe":
        vs = [float(r["v"]) for r in charge.rebound_rows()]
    else:
        name = "argon" if material == "argon" else "water"
        vs = [
            _f(r, "v")
            for r in (*iapetus_bag.rebound_rows(), *fast_rebound_rows())
            if r["material"] == name
        ]
    return min(vs) <= v <= max(vs)


@dataclass(frozen=True)
class MaterialRow:
    plate: str
    material: str
    heat_per_kg: float  # [J/kg] of merged mass
    kt_start_ev: float
    gamma_eff: float
    gamma_factor: float
    ion_share_start: float  # of the heat
    ion_share_frozen: float
    da_ion_min: float
    bond_share_held_at_plate: float  # equilibrium bonds still broken at the plate (water)
    bond_toll_frozen: float  # of the heat: PE carbon atomic / water stranded at its freeze
    bond_toll_returned: float  # PE carbon condensed; water fully re-formed
    da_bond_min: float  # water's least Damkohler number while bonds return; PE's H2 race
    velocity_frozen: float  # sqrt(1 - frozen ions - frozen bonds)
    velocity_returned: float  # sqrt(1 - frozen ions - returned-case toll)
    arrival_speed: float  # [m/s]
    column: float  # [kg/m^2]
    e_eff: float
    e_eff_in_grid: bool
    factor_frozen: float  # gamma x rebound x velocity, over polyethylene's (frozen against frozen)
    factor_returned: float


def arrival(p: Plate) -> tuple[float, float]:
    """`(u_rms, column)`: the merged gas's speed toward the plate and its arriving areal mass."""
    u_rms = p.w * math.sqrt(p.k) / (1.0 + p.k)
    if p.name == "push":
        return u_rms, (1.0 + p.k) * PUSH_SHEET
    return u_rms, (1.0 + p.k) * p.impactor / (2.0 * math.pi * p.plate_radius**2)


@cache
def _material(plate_name: str, material: str) -> MaterialRow:
    p = plate(plate_name)
    r = iapetus_bag.race(material, p.k, p.w, p.plate_distance, bag=p.impactor)
    # Polyethylene keeps the record's pricing (Results 10-15): the 2-D gas's gamma 1.4 with its bond
    # toll charged apart. Its CH2 table counts atoms only, so the gas cooling toward molecules is
    # not in it and a gamma read off it is not meaningful (at the departure it fails outright).
    g = PE_GAMMA if material == "pe" else iapetus_bag.effective_gamma(r)
    gf = gamma_factor(p.name, g)
    if material == "pe":
        f = freeze.freeze(
            p.k, charge.ch2_table(), p.plate_distance, v_rel=p.w, projectile=p.impactor
        )
        carbon_part = (freeze.CH2_ATOMIZATION - freeze.H2_BOND) / freeze.CH2_ATOMIZATION
        ions = f.ion_share_frozen
        frozen_bonds = f.bond_share * carbon_part
        returned_bonds = 0.0
        held = 0.0
        da_b = f.da_h2
    else:
        ions = r.ion_share_frozen
        # Water's bonds are inside its equilibrium table; only a stranded share is a toll.
        frozen_bonds = max(r.stranded_h, r.stranded_oh) if material == "ice" else 0.0
        returned_bonds = 0.0
        held = r.bond_share_end_eq
        da_b = min(r.da_bond_min_h, r.da_bond_min_oh) if material == "ice" else math.nan
    v_frozen = math.sqrt(max(1.0 - ions - frozen_bonds, 0.0))
    v_returned = math.sqrt(max(1.0 - ions - returned_bonds, 0.0))
    u, col = arrival(p)
    e = rebound(material, u, col)
    return MaterialRow(
        plate=p.name,
        material=material,
        heat_per_kg=0.5 * p.k / (1.0 + p.k) ** 2 * p.w**2,
        kt_start_ev=r.kt_start_ev,
        gamma_eff=g,
        gamma_factor=gf,
        ion_share_start=r.ion_share_start,
        ion_share_frozen=ions,
        da_ion_min=r.da_ion_min,
        bond_share_held_at_plate=held,
        bond_toll_frozen=frozen_bonds,
        bond_toll_returned=returned_bonds,
        da_bond_min=da_b,
        velocity_frozen=v_frozen,
        velocity_returned=v_returned,
        arrival_speed=u,
        column=col,
        e_eff=e,
        e_eff_in_grid=rebound_speed_covered(material, u),
        factor_frozen=math.nan,
        factor_returned=math.nan,
    )


def material_rows() -> list[MaterialRow]:
    out = []
    for p in PLATES:
        pe = _material(p.name, "pe")
        for m in ("pe", "ice", "argon"):
            row = _material(p.name, m)

            def lever(x: MaterialRow, v: float) -> float:
                return x.gamma_factor * (1.0 + x.e_eff) * v

            out.append(
                MaterialRow(
                    **(
                        asdict(row)
                        | {
                            "factor_frozen": lever(row, row.velocity_frozen)
                            / lever(pe, pe.velocity_frozen),
                            "factor_returned": lever(row, row.velocity_returned)
                            / lever(pe, pe.velocity_returned),
                        }
                    )
                )
            )
    return out


# --- 3. Porosity: the 1-D bound ----------------------------------------------------------------


@dataclass(frozen=True)
class Hugoniot:
    """A linear `U_s = c0 + s u_p` Hugoniot and a constant `Gamma rho`. These are fits to gas-gun
    data of order 5-10 km/s in `u_p`; they are used only at the snowplow's slowest strike."""

    name: str
    rho0: float
    c0: float
    s: float
    gamma0: float
    source: str


HUGONIOTS = {
    # LASL Shock Hugoniot Data (Marsh 1980), polyethylene at 954 kg/m^3; Gamma_0 is a bracket.
    "pe": Hugoniot("polyethylene", 954.0, 2901.0, 1.481, 1.6, "Marsh 1980 (LASL), linear fit"),
    # Steinberg 1996 (UCRL-MA-106439), water's linear term; ice takes water's above its phases.
    "water": Hugoniot("water", 998.0, 1647.0, 1.921, 0.35, "Steinberg 1996, linear term"),
    "ice": Hugoniot("ice", 917.0, 1647.0, 1.921, 0.35, "water's fit at ice's density"),
    # No fit is taken for solid argon: a generic soft solid, bracketed like the others.
    "argon": Hugoniot("argon", 1620.0, 1300.0, 1.5, 1.0, "generic bracket, not a fit"),
}
#: The bracket every fit is run across: `s` +/- 0.2 and `Gamma_0` 0.5-2.
S_BRACKET = (-0.2, 0.0, 0.2)
GAMMA_BRACKET = (0.5, 2.0)

#: Energy to vaporize from the stored solid [J/kg]. Water: ice at 250 K heated, melted, boiled
#: (steam tables; sublimation in vacuum is less, 2.83 MJ/kg). Polyethylene: the upper end of its
#: heat of gasification (Tewarson, SFPE Handbook, 1.8-3.6 MJ/kg), with full atomization (84 MJ/kg)
#: as a harder bar. Argon: sublimation from 60 K (NIST WebBook, 7.6 kJ/mol near the triple point).
VAPORIZE = {
    "pe": 3.6e6,
    "water": 2100.0 * 23.0 + 334.0e3 + 4186.0 * 100.0 + 2257.0e3,
    "ice": 2100.0 * 23.0 + 334.0e3 + 4186.0 * 100.0 + 2257.0e3,
    "argon": 7.61e3 / 0.039948 + 700.0 * 24.0,
}
ATOMIZE = {
    "pe": charge.PE_BONDS,
    "water": eos_water.FULL_ATOMIZATION_ENERGY,
    "ice": eos_water.FULL_ATOMIZATION_ENERGY,
}


def cold_share(
    h: Hugoniot, up: float, s: float | None = None, gamma0: float | None = None
) -> float:
    """The share of the Hugoniot energy `u_p^2 / 2` a solid holds as cold compression: the energy
    of the isentrope through the unshocked state at the shocked volume, over `E_H`.

    Mie-Gruneisen about the Hugoniot, `P = P_H(V) + (Gamma_0 rho0) (E - E_H(V))`, with
    `P_H = rho0 c0^2 eta / (1 - s eta)^2`, `E_H = P_H (V0 - V) / 2`, `eta = 1 - V/V0`; the
    isentrope `dE/dV = -P` is integrated from `(V0, 0)`. A porous charge struck at the same speed
    gets the same `E_H` almost wholly as heat, so this share bounds what pore collapse can add."""
    s = h.s if s is None else s
    g = h.gamma0 if gamma0 is None else gamma0
    us = h.c0 + s * up
    v0 = 1.0 / h.rho0
    vh = v0 * (1.0 - up / us)
    eh = 0.5 * up * up

    def p_h(v: float) -> float:
        eta = 1.0 - v / v0
        return h.rho0 * h.c0**2 * eta / (1.0 - s * eta) ** 2

    n = 4000
    e, v = 0.0, v0
    dv = (vh - v0) / n
    for _ in range(n):
        # Midpoint step of dE/dV = -P(V, E).
        def rhs(vv: float, ee: float) -> float:
            ph = p_h(vv)
            return -(ph + g * h.rho0 * (ee - 0.5 * ph * (v0 - vv)))

        k1 = rhs(v, e)
        k2 = rhs(v + 0.5 * dv, e + 0.5 * dv * k1)
        e += dv * k2
        v += dv
    return min(max(e / eh, 0.0), 1.0)


@dataclass(frozen=True)
class HugoniotRow:
    plate: str
    material: str
    v_first: float  # [m/s], the first-struck charge's
    v_last: float  # [m/s], the last-struck charge's: w / (1 + k)
    e_last: float  # [J/kg], its shock heat V^2/2 (exact)
    vaporize: float
    margin_vaporize: float  # e_last over the vaporization energy
    margin_atomize: float  # e_last over full atomization (nan for argon)
    mixed_heat: float  # [J/kg], the merged mass's mean heat, once shocks have shared it
    mixed_margin_atomize: float
    cold_share_nominal: float  # at v_last, the nominal fit
    cold_share_max: float  # the bracket's largest
    solid_heat_margin: float  # (1 - cold_share_max) e_last over vaporization
    in_fit_range: bool  # v_last <= 10 km/s


def hugoniot_rows() -> list[HugoniotRow]:
    out = []
    for p in PLATES:
        v_last = p.w / (1.0 + p.k)
        e_last = 0.5 * v_last**2
        mixed = 0.5 * p.k / (1.0 + p.k) ** 2 * p.w**2
        for m in ("pe", "ice", "water", "argon"):
            h = HUGONIOTS[m]
            shares = [cold_share(h, v_last, h.s + ds, g) for ds in S_BRACKET for g in GAMMA_BRACKET]
            nominal = cold_share(h, v_last)
            cmax = max([*shares, nominal])
            atom = ATOMIZE.get(m, math.nan)
            out.append(
                HugoniotRow(
                    plate=p.name,
                    material=m,
                    v_first=p.w,
                    v_last=v_last,
                    e_last=e_last,
                    vaporize=VAPORIZE[m],
                    margin_vaporize=e_last / VAPORIZE[m],
                    margin_atomize=e_last / atom,
                    mixed_heat=mixed,
                    mixed_margin_atomize=mixed / atom,
                    cold_share_nominal=nominal,
                    cold_share_max=cmax,
                    solid_heat_margin=(1.0 - cmax) * e_last / VAPORIZE[m],
                    in_fit_range=v_last <= 10.0e3,
                )
            )
    return out


# --- 4. The graded charge ------------------------------------------------------------------------


@dataclass(frozen=True)
class GradedRow:
    plate: str
    rho_light: float
    rho_dense: float
    dense_share: float
    dense_side: str  # "plate" or "far"
    dense_side_is_entry: bool  # the dense part faces the impactor's entry
    length: float
    impulse_per_mv: float
    over_solid: float  # over the uniform 950 kg/m^3 column
    over_light_uniform: float  # over a uniform column at the light density
    over_mean_uniform: float  # over a uniform column of the same length (mean density)
    peak_over_solid: (
        float  # peak face pressure over the uniform 950 kg/m^3 column's (not converged)
    )


def _solid_peak(plate_name: str, path: Path) -> float:
    """The reference charge's peak: the 950 kg/m^3 column, or the push's 0.8 m uniform sheet."""
    if plate_name == "push":
        return next(
            _f(r, "peak_per_mv2")
            for r in runs(path)
            if r["plate"] == "push"
            and r["kind"] == "sheet"
            and abs(_f(r, "length") - 0.8) < 1e-9
            and _f(r, "dense_share") == 0.0
            and abs(_f(r, "gamma") - 1.4) < 1e-3
        )
    return next(
        _f(r, "peak_per_mv2")
        for r in runs(path)
        if r["plate"] == plate_name
        and r["kind"] == "uniform"
        and _f(r, "rho") == PE_SOLID
        and abs(_f(r, "gamma") - 1.4) < 1e-3
    )


def graded_rows(path: Path = JSONL) -> list[GradedRow]:
    out = []
    for r in runs(path):
        if r["kind"] not in ("graded", "sheet") or _f(r, "dense_share") == 0.0:
            continue
        name = str(r["plate"])
        j = _f(r, "impulse_per_mv")
        side_plate = bool(r["dense_on_plate_side"])
        if name == "push":
            ref = next(t for x, t in push_sheets(path) if abs(x - _f(r, "length")) < 1e-9)
            out.append(
                GradedRow(
                    name,
                    math.nan,
                    math.nan,
                    _f(r, "dense_share"),
                    "plate" if side_plate else "far",
                    # The push impactor arrives from beyond the sheet, heading for the plate.
                    not side_plate,
                    _f(r, "length"),
                    j,
                    j / ref,
                    math.nan,
                    j / ref,
                    _f(r, "peak_per_mv2") / _solid_peak(name, path),
                )
            )
            continue
        p = plate(name)
        pts = _uniform(name, path)
        lr = [math.log(a) for a, _ in pts]
        js = [b for _, b in pts]
        j_solid = float(np.interp(math.log(PE_SOLID), lr, js))
        rho_l = _f(r, "rho")
        j_light = float(np.interp(math.log(rho_l), lr, js))
        mean = p.k * p.impactor / (math.pi * p.disk_radius**2 * _f(r, "length"))
        j_mean = float(np.interp(math.log(mean), lr, js))
        out.append(
            GradedRow(
                name,
                rho_l,
                _f(r, "rho_dense"),
                _f(r, "dense_share"),
                "plate" if side_plate else "far",
                # On both bowls the impactor enters from the plate side.
                side_plate,
                _f(r, "length"),
                j,
                j / j_solid,
                j / j_light,
                j / j_mean,
                _f(r, "peak_per_mv2") / _solid_peak(name, path),
            )
        )
    return out


# --- 5. The bag -----------------------------------------------------------------------------------

FILM = 0.05  # [kg/m^2], the parent's polyethylene film
THIN_FILM = 0.0095  # [kg/m^2], 10 um of polyethylene, as a lighter option
#: Gas held in the bag: saturated water vapour over ice at 273.16 K (611.657 Pa, IAPWS) and over
#: argon frost at 60 K (`argon_sublimation_pressure`).
WATER_TRIPLE_P, WATER_TRIPLE_T = 611.657, 273.16


@dataclass(frozen=True)
class BagRow:
    plate: str
    charge: str
    rho: float
    charge_mass: float  # [kg]
    film_area: float  # [m^2]
    film_share: float  # film mass over the charge mass, at 0.05 kg/m^2
    thin_film_share: float  # the same at 10 um
    gas_share: float  # held vapour over the charge mass


def vapour_density(material: str) -> float:
    """[kg/m^3] of the gas a sealed bag holds over its frozen charge."""
    if material == "argon":
        t = 60.0
        return argon_sublimation_pressure(t) * 0.039948 / (R_GAS * t)
    return WATER_TRIPLE_P * 0.018015 / (R_GAS * WATER_TRIPLE_T)


def bag_rows() -> list[BagRow]:
    out = []
    for p in PLATES:
        for c in CHARGES:
            if c.material == "pe":
                continue  # polyethylene, solid or foam, holds its own shape: no bag
            solid = {"pe": PE_SOLID, "ice": 917.0, "argon": 1620.0}[c.material]
            if p.name == "push":
                mass = p.k * PUSH_SHEET  # per m^2 of sheet, two faces of film
                area = 2.0
            else:
                mass = p.k * p.impactor
                length = mass / (c.rho * math.pi * p.disk_radius**2)
                area = 2.0 * math.pi * p.disk_radius**2 + 2.0 * math.pi * p.disk_radius * length
            pores = max(mass * (1.0 / c.rho - 1.0 / solid), 0.0)
            gas = 0.0 if c.material == "pe" else vapour_density(c.material) * pores
            out.append(
                BagRow(
                    p.name,
                    c.name,
                    c.rho,
                    mass,
                    area,
                    FILM * area / mass,
                    THIN_FILM * area / mass,
                    gas / mass,
                )
            )
    return out


# --- 6. Argon's flash, and the debris -----------------------------------------------------------

#: Argon (NIST WebBook): ITS-90 triple point 83.8058 K, 68.89 kPa; normal boiling point 87.30 K;
#: fusion 1.18 kJ/mol and vaporization 6.43 kJ/mol; liquid c_p ~1.12 kJ/kg/K, solid c_p ~0.7
#: kJ/kg/K near 60-80 K (approximate).
AR_M = 0.039948
AR_T_TRIPLE, AR_P_TRIPLE = 83.8058, 68.89e3
AR_T_BOIL = 87.30
AR_L_FUS = 1.18e3 / AR_M
AR_L_VAP = 6.43e3 / AR_M
AR_L_SUB = AR_L_FUS + AR_L_VAP
AR_CP_LIQ, AR_CP_SOLID = 1120.0, 700.0


def argon_sublimation_pressure(t: float) -> float:
    """[Pa] over argon frost, Clausius-Clapeyron from the triple point at constant `L_sub`."""
    return AR_P_TRIPLE * math.exp(-AR_L_SUB * AR_M / R_GAS * (1.0 / t - 1.0 / AR_T_TRIPLE))


@dataclass(frozen=True)
class ArgonFlash:
    flash_to_triple: (
        float  # vapour share of liquid at its boiling point flashed to the triple point
    )
    flash_to_60k: float  # and on cooling the frost to 60 K by further sublimation
    p_60k: float  # [Pa]
    vapour_density_60k: float  # [kg/m^3]
    vapour_density_triple: float
    volume_per_kg_60k: float  # [m^3] to hold the parent's 22% as vapour, per kg of charge
    volume_per_kg_triple: float
    heat_per_kg: float  # [J/kg] to radiate so the vapour deposits
    radiate_days: dict[str, float]  # shaded sphere at 70 K, emissivity 0.9, per plate's charge


def argon_flash(parent_vapour: float = 0.22, t_bag: float = 70.0) -> ArgonFlash:
    """What a sealed argon bag holds. Liquid at its boiling point sprayed into vacuum flashes to
    solid and vapour at the triple point; cooling the frost to 60 K sublimes more. In a sealed bag
    that vapour stays vapour only until the bag radiates its latent heat, after which it lies on
    the frost: the charge holds only the saturated gas."""
    x1 = (AR_CP_LIQ * (AR_T_BOIL - AR_T_TRIPLE) + AR_L_FUS) / (AR_L_FUS + AR_L_VAP)
    x2 = x1 + (1.0 - x1) * AR_CP_SOLID * (AR_T_TRIPLE - 60.0) / AR_L_SUB
    p60 = argon_sublimation_pressure(60.0)
    rho60 = p60 * AR_M / (R_GAS * 60.0)
    rho_t = AR_P_TRIPLE * AR_M / (R_GAS * AR_T_TRIPLE)
    heat = parent_vapour * AR_L_SUB
    days = {}
    for p in PLATES:
        if p.name == "push":
            continue
        mass = p.k * p.impactor
        r = (3.0 * mass / (4.0 * math.pi * 1620.0)) ** (1.0 / 3.0)
        power = 0.9 * SIGMA * t_bag**4 * 4.0 * math.pi * r * r
        days[p.name] = heat * mass / power / 86400.0
    return ArgonFlash(
        flash_to_triple=x1,
        flash_to_60k=x2,
        p_60k=p60,
        vapour_density_60k=rho60,
        vapour_density_triple=rho_t,
        volume_per_kg_60k=parent_vapour / rho60,
        volume_per_kg_triple=parent_vapour / rho_t,
        heat_per_kg=heat,
        radiate_days=days,
    )


def ice_vapour_pressure(t: float) -> float:
    """[Pa] over hexagonal ice, Murphy & Koop 2005 (QJRMS 131, eq. 7), valid above 110 K."""
    return math.exp(9.550426 - 5723.265 / t + 3.53068 * math.log(t) - 0.00728332 * t)


@dataclass(frozen=True)
class DebrisRow:
    material: str
    albedo: float
    distance_rsun: float
    temperature: float  # [K], a sunlit isothermal sphere cooled by radiation and sublimation
    flux: float  # [kg/m^2/s] sublimed
    lifetimes_days: dict[str, float]  # by radius: the time to sublime a sphere at its bulk density


DEBRIS_RADII = {"0.1 mm": 1e-4, "1 mm": 1e-3, "1 cm": 1e-2, "10 cm": 0.1, "0.4 m": 0.4}


def debris_row(material: str, albedo: float, distance_rsun: float, rho: float) -> DebrisRow:
    """Sublimation of a fast-rotating (isothermal) sunlit sphere: it absorbs `S (1 - A) / 4` per
    m^2 of surface, and loses it as `eps sigma T^4 + L Z(T)`, `Z = P_sat sqrt(M / 2 pi R T)`
    (Hertz-Knudsen, unit sticking). The radius falls at `Z / rho`, so a sphere lasts `r rho / Z`.
    Earth's shadow and infrared are left out."""
    s = SOLAR_CONSTANT * (AU_IN_RSUN / distance_rsun) ** 2
    absorbed = s * (1.0 - albedo) / 4.0
    if material == "argon":
        p_sat, m, latent, eps = argon_sublimation_pressure, AR_M, AR_L_SUB, 0.9
    else:
        p_sat, m, latent, eps = ice_vapour_pressure, 0.018015, 2.83e6, 0.97

    def z(t: float) -> float:
        return p_sat(t) * math.sqrt(m / (2.0 * math.pi * R_GAS * t))

    lo, hi = 20.0 if material == "argon" else 110.0, 400.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if eps * SIGMA * mid**4 + latent * z(mid) > absorbed:
            hi = mid
        else:
            lo = mid
    t = 0.5 * (lo + hi)
    flux = z(t)
    return DebrisRow(
        material,
        albedo,
        distance_rsun,
        t,
        flux,
        {name: r * rho / flux / 86400.0 for name, r in DEBRIS_RADII.items()},
    )


def debris_rows() -> list[DebrisRow]:
    out = []
    for d in (AU_IN_RSUN, 9.86):
        out.append(debris_row("ice", 0.5, d, 917.0))
        out.append(debris_row("ice", 0.05, d, 917.0))  # dark snow
        out.append(debris_row("snow 300", 0.05, d, 300.0))
        out.append(debris_row("argon", 0.5, d, 1620.0))
    return out


#: A lost film fragment's area-to-mass ratio [m^2/kg] at 0.05 kg/m^2: solar pressure and drag
#: act on it ~4000 times as hard as on a 0.2 m solid polyethylene block (190 kg/m^2).
def film_area_to_mass() -> tuple[float, float]:
    block = PE_SOLID * 0.2
    return 1.0 / FILM, 1.0 / block


# --- 7. The low-speed panel --------------------------------------------------------------------

#: Seth (2026-10-24): the 100 km/s Earth runs and the 3S cycle's straight Jupiter return (55 km/s,
#: 44 km/s once up to speed). The 2-D runs are scale-free, so only the chemistry, the rebound and
#: the porosity margin are re-run. The plate is not named, so the push sheet's gamma bracket (the
#: monotonic one) prices gamma, a 15 kg impactor and the bowls' 3.75 m distance set the race.
LOW_SPEEDS = (44.0e3, 55.0e3, 100.0e3)
LOW_K = (5.0, 10.0, 20.0)
LOW_IMPACTOR = 15.0
LOW_DISTANCE = 3.75


@dataclass(frozen=True)
class LowSpeedRow:
    w: float
    k: float
    material: str
    heat_per_kg: float  # [J/kg] of merged mass
    full_bonds_share: float  # full atomization over the heat (the frozen worst case; 0 for argon)
    kt_start_ev: float
    gamma_eff: float  # nan for polyethylene (priced at 1.4) or where the race leaves its table
    gamma_factor: float
    gamma_in_bracket: bool
    ion_share_start: float
    ion_share_frozen: float
    da_ion_min: float
    bond_share_held_at_plate: float
    bond_held_after_push: (
        float  # water: of the stagnation energy, when the rebound layer has thinned
    )
    bond_stranded_in_push: float  # water: left behind by H + OH + M on the rebound
    bond_toll_frozen: float
    velocity_frozen: float
    arrival_speed: float
    e_eff: float
    e_eff_in_grid: bool
    last_struck_heat: float  # [J/kg], (w / (1 + k))^2 / 2
    vaporize_margin: float
    lever_common: float  # (1 + e_eff) x velocity, every material at the 2-D gas's gamma 1.4
    over_argon_common: float
    lever: float  # times the push sheet's gamma factor (polyethylene held at 1.4: favoured)
    over_argon: float
    gamma_factor_bowl: float  # the departure bowl's bracket instead
    over_argon_bowl: float


def _low(w: float, k: float, material: str) -> LowSpeedRow:
    heat = 0.5 * k / (1.0 + k) ** 2 * w * w
    r = iapetus_bag.race(material, k, w, LOW_DISTANCE, bag=LOW_IMPACTOR)
    if material == "pe":
        g = math.nan
        gf = gfb = 1.0
        f = freeze.freeze(k, charge.ch2_table(), LOW_DISTANCE, v_rel=w, projectile=LOW_IMPACTOR)
        carbon_part = (freeze.CH2_ATOMIZATION - freeze.H2_BOND) / freeze.CH2_ATOMIZATION
        ions, bonds, held = f.ion_share_frozen, f.bond_share * carbon_part, 0.0
    else:
        thermal = 1.0 - r.kinetic_share_end
        g = iapetus_bag.effective_gamma(r) if thermal > 0.0 else math.nan
        gf = gamma_factor("push", g) if math.isfinite(g) else math.nan
        gfb = gamma_factor("departure", g) if math.isfinite(g) else math.nan
        ions = r.ion_share_frozen
        bonds = max(r.stranded_h, r.stranded_oh) if material == "ice" else 0.0
        held = r.bond_share_end_eq
    v = math.sqrt(max(1.0 - ions - bonds, 0.0))
    u = w * math.sqrt(k) / (1.0 + k)
    col = (1.0 + k) * LOW_IMPACTOR / (2.0 * math.pi * 10.0**2)
    # Water's bonds through the plate: stagnation re-heats the gas, and the race follows the
    # rebound layer until it has thinned tenfold (`iapetus_bag.race`'s rebound leg).
    after = stranded = math.nan
    if material == "ice":
        rb = iapetus_bag.race(
            material,
            k,
            w,
            LOW_DISTANCE,
            leg="rebound",
            rebound_column=max(col, 0.1),
            bag=LOW_IMPACTOR,
        )
        after, stranded = rb.bond_share_end_eq, max(rb.stranded_h, rb.stranded_oh)
    e = rebound(material, u, col)
    last = 0.5 * (w / (1.0 + k)) ** 2
    vap = VAPORIZE["ice" if material == "ice" else material]
    full = ATOMIZE.get(material, 0.0) * k / (1.0 + k) / heat
    return LowSpeedRow(
        w=w,
        k=k,
        material=material,
        heat_per_kg=heat,
        full_bonds_share=full,
        kt_start_ev=r.kt_start_ev,
        gamma_eff=g,
        gamma_factor=gf,
        gamma_in_bracket=material == "pe" or gamma_in_bracket(g),
        ion_share_start=r.ion_share_start,
        ion_share_frozen=ions,
        da_ion_min=r.da_ion_min,
        bond_share_held_at_plate=held,
        bond_held_after_push=after,
        bond_stranded_in_push=stranded,
        bond_toll_frozen=bonds,
        velocity_frozen=v,
        arrival_speed=u,
        e_eff=e,
        e_eff_in_grid=rebound_speed_covered(material, u),
        last_struck_heat=last,
        vaporize_margin=last / vap,
        lever_common=(1.0 + e) * v,
        over_argon_common=math.nan,
        lever=gf * (1.0 + e) * v,
        over_argon=math.nan,
        gamma_factor_bowl=gfb,
        over_argon_bowl=math.nan,
    )


def low_speed_rows() -> list[LowSpeedRow]:
    out = []
    for w in LOW_SPEEDS:
        for k in LOW_K:
            rows = {m: _low(w, k, m) for m in ("argon", "ice", "pe")}
            ar, ar_c = rows["argon"].lever, rows["argon"].lever_common
            ar_b = rows["argon"].gamma_factor_bowl * ar_c
            out += [
                LowSpeedRow(
                    **(
                        asdict(x)
                        | {
                            "over_argon": x.lever / ar,
                            "over_argon_common": x.lever_common / ar_c,
                            "over_argon_bowl": x.gamma_factor_bowl * x.lever_common / ar_b,
                        }
                    )
                )
                for x in rows.values()
            ]
    return out


# --- Output ---------------------------------------------------------------------------------------


@dataclass(frozen=True)
class LedgerRow:
    plate: str
    charge: str
    rho: float
    density_factor: float
    material_factor_frozen: float
    material_factor_returned: float
    film_share: float
    total_frozen: float  # density x material x the film's dilution
    total_returned: float


def ledger_rows(path: Path = JSONL) -> list[LedgerRow]:
    """Each charge's impulse per kg of impactor (equivalently its v_e at the same k) over solid
    polyethylene's on each plate. The bag's film is charged as dead mass inside the same k (it is
    polyethylene, so it couples like the charge; the dilution only replaces charge material)."""
    mats = {(r.plate, r.material): r for r in material_rows()}
    bags = {(b.plate, b.charge): b for b in bag_rows()}
    out = []
    for p in PLATES:
        for c in CHARGES:
            m = mats[(p.name, c.material)]
            d = density_factor(p.name, c.rho, path)
            film = bags[(p.name, c.name)].film_share if (p.name, c.name) in bags else 0.0
            # The film is polyethylene: its share of the charge couples at factor 1.
            mix_f = (1.0 - film) * m.factor_frozen + film
            mix_r = (1.0 - film) * m.factor_returned + film
            out.append(
                LedgerRow(
                    p.name,
                    c.name,
                    c.rho,
                    d,
                    m.factor_frozen,
                    m.factor_returned,
                    film,
                    d * mix_f,
                    d * mix_r,
                )
            )
    return out


def _write(path: Path, rows: list[object]) -> None:
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    if not dicts:
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0].keys()))
        w.writeheader()
        for d in dicts:
            w.writerow({k: (str(v) if isinstance(v, dict) else v) for k, v in d.items()})


def main() -> None:
    if not runs(JSONL):
        print("(no snow_charge_2d*.jsonl: run `cargo run --release -p sweep -- --snow-charge`)")
        return
    dens = density_rows()
    _write(OUT_DIR / "snow_charge_density.csv", list(dens))
    print("== 1. density law (J over solid 950; push: over the 0.2 m grid sheet)")
    for d in dens:
        print(
            f"  {d.plate:9s} rho {d.rho:6.0f} L {d.length:5.2f} m L/R {d.length_over_radius:5.2f} "
            f"J/mv {d.impulse_per_mv:6.3f} x{d.factor:.3f} (coarse x{d.coarse_factor:.3f}) "
            f"real/point {d.real_over_point:.3f}"
        )
    for name in ("sun", "departure"):
        print(
            f"  {name}: within 2% of solid down to {breakpoint(name):.0f} kg/m^3, "
            f"5% to {breakpoint(name, 0.05):.0f}"
        )
    print(f"  push: thickness 0.1-0.8 m spreads the impulse by {push_thickness_spread():.1%}")

    mats = material_rows()
    _write(OUT_DIR / "snow_charge_material.csv", list(mats))
    print("== 2-3. material (over polyethylene, frozen / returned)")
    for m in mats:
        print(
            f"  {m.plate:9s} {m.material:5s} heat {m.heat_per_kg / 1e6:7.0f} MJ/kg kT0 "
            f"{m.kt_start_ev:5.1f} eV gamma {m.gamma_eff:.3f} (J x{m.gamma_factor:.3f}) ions "
            f"{m.ion_share_start:.2f}->{m.ion_share_frozen:.3f} (Da {m.da_ion_min:.1e}) bonds "
            f"held {m.bond_share_held_at_plate:.3f} toll {m.bond_toll_frozen:.3f} (Da "
            f"{m.da_bond_min:.1e}) e_eff {m.e_eff:.2f}{'' if m.e_eff_in_grid else '*'} at "
            f"{m.arrival_speed / 1e3:.0f} km/s -> x{m.factor_frozen:.3f} / x{m.factor_returned:.3f}"
        )

    hug = hugoniot_rows()
    _write(OUT_DIR / "snow_charge_hugoniot.csv", list(hug))
    print("== porosity bound at the slowest strike")
    for h in hug:
        print(
            f"  {h.plate:9s} {h.material:5s} V_last {h.v_last / 1e3:5.1f} km/s "
            f"E {h.e_last / 1e6:6.1f} MJ/kg = {h.margin_vaporize:6.1f}x vaporize, "
            f"{h.margin_atomize:.2f}x atomize "
            f"(mixed {h.mixed_margin_atomize:.1f}x); cold share {h.cold_share_nominal:.3f} "
            f"(max {h.cold_share_max:.3f}); solid heat {h.solid_heat_margin:.1f}x vaporize"
            f"{'' if h.in_fit_range else ' [fit extrapolated]'}"
        )

    gr = graded_rows()
    _write(OUT_DIR / "snow_charge_graded.csv", list(gr))
    print("== 7. graded")
    for g in sorted(gr, key=lambda x: (x.plate, x.rho_light, x.dense_share, x.dense_side)):
        print(
            f"  {g.plate:9s} {g.rho_light:5.0f}/{g.rho_dense:4.0f} share {g.dense_share:.2f} "
            f"dense on {g.dense_side:5s} (entry {g.dense_side_is_entry!s:5s}) "
            f"L {g.length:4.2f} J/mv "
            f"{g.impulse_per_mv:6.3f} x{g.over_solid:.3f} solid, "
            f"x{g.over_light_uniform:.3f} light, x{g.over_mean_uniform:.3f} same length"
        )

    bags = bag_rows()
    _write(OUT_DIR / "snow_charge_bag.csv", list(bags))
    print("== 4. bag")
    for b in bags:
        print(
            f"  {b.plate:9s} {b.charge:22s} film {b.film_share:.4f} (10 um {b.thin_film_share:.4f})"
            f" gas {b.gas_share:.1e}"
        )

    fl = argon_flash()
    print(
        f"== argon: flash {fl.flash_to_triple:.3f} to the triple point, "
        f"{fl.flash_to_60k:.3f} to 60 K; 60 K frost holds {fl.p_60k:.0f} Pa, "
        f"{fl.vapour_density_60k:.3f} kg/m^3: 22% as vapour "
        f"needs {fl.volume_per_kg_60k:.1f} m^3/kg at 60 K ({fl.volume_per_kg_triple:.3f} at the "
        f"triple point); depositing it radiates {fl.heat_per_kg / 1e3:.0f} kJ/kg: "
        + ", ".join(f"{k} {v:.0f} days" for k, v in fl.radiate_days.items())
    )
    with (OUT_DIR / "snow_charge_argon_flash.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        for k, v in asdict(fl).items():
            w.writerow([k, v])

    deb = debris_rows()
    _write(OUT_DIR / "snow_charge_debris.csv", list(deb))
    print("== 6. debris lifetimes (sunlit)")
    for db in deb:
        print(
            f"  {db.material:8s} A {db.albedo:.2f} at {db.distance_rsun:6.1f} R_sun: "
            f"T {db.temperature:5.1f} K flux {db.flux:.2e} kg/m^2/s; "
            + ", ".join(f"{k} {v:.3g} d" for k, v in db.lifetimes_days.items())
        )
    a_film, a_block = film_area_to_mass()
    print(f"  film {a_film:.0f} m^2/kg against a 0.2 m PE block's {a_block:.4f}")

    low = low_speed_rows()
    _write(OUT_DIR / "snow_charge_low_speed.csv", list(low))
    print("== low-speed panel (15 kg impactor; * outside a grid)")
    for lw in low:
        print(
            f"  {lw.w / 1e3:4.0f} km/s k {lw.k:4.0f} {lw.material:5s} "
            f"heat {lw.heat_per_kg / 1e6:6.1f} MJ/kg full bonds {lw.full_bonds_share:4.0%} "
            f"gamma {lw.gamma_eff:.3f}"
            f"{'' if lw.gamma_in_bracket else '*'} ions {lw.ion_share_start:.2f}->"
            f"{lw.ion_share_frozen:.4f} (Da {lw.da_ion_min:.1e}) "
            f"bonds held {lw.bond_share_held_at_plate:.2f} toll {lw.bond_toll_frozen:.2f} "
            f"v x{lw.velocity_frozen:.3f} e_eff {lw.e_eff:.2f}"
            f"{'' if lw.e_eff_in_grid else '*'} at {lw.arrival_speed / 1e3:.1f} km/s; last-struck "
            f"{lw.vaporize_margin:.1f}x vaporize; "
            f"x{lw.over_argon_common:.3f} argon at common gamma, "
            f"x{lw.over_argon:.3f} with the sheet's gamma, x{lw.over_argon_bowl:.3f} the bowl's"
        )

    led = ledger_rows()
    _write(OUT_DIR / "snow_charge_ledger.csv", list(led))
    print("== ledger: impulse per kg of impactor over solid polyethylene (frozen / returned)")
    for lg in led:
        print(
            f"  {lg.plate:9s} {lg.charge:22s} density x{lg.density_factor:.3f} material "
            f"x{lg.material_factor_frozen:.3f}/x{lg.material_factor_returned:.3f} film "
            f"{lg.film_share:.3f} -> x{lg.total_frozen:.3f} / x{lg.total_returned:.3f}"
        )


if __name__ == "__main__":
    main()
