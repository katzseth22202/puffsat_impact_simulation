"""Result 20: Iapetus water bags on the overtake plate, with carried charge (parent S50).

    PYTHONPATH=python uv run python -m puffsat.xray_collision.iapetus_bag

Water mined on Iapetus comes back as 25 kg bags that overtake an Earth-launched ship at 51.1 km/s
(60 and 68 km/s as sensitivities) and push it to 10.916 km/s. Each bag meets a charge `k` times its
mass held off the flown 9.5 m push plate, as in Result 11's overtake. Unlike the floor-wide Kapton
sheet of Results 13-19, a bag is compact, so the pulse starts from near a point. Seth asked three
things. Should the charge be argon or water, given that water spends energy on bonds that may not
re-form? Is the skirted flat plate enough, or must it be a paraboloid? And what does the plate
weigh at 2 Hz and 1 g? The parent's S50 adds the arrival state, aim and face erosion.

The pieces:

1. **Arrival state** (`arrival_rows`). A bag released ahead breaks into fragments. These couple to a
   charge sheet only if they cover their own area. A fragment of radius `a` sweeps
   `3 k sigma / (4 rho a)` of its mass in charge, which is `k` times the cloud's optical depth.
2. **The chemistry** (`race`). The fireball expands from solid density on the argon or the water
   table. Its ionization store, and water's bonds, are raced against the expansion
   (`recombination.freeze_station`: three-body `H + OH + M`, bracketed between H-limited and
   OH-limited). The same race is run on the gas re-expanding off the plate after it stagnates there.
3. **The k sweep** (`sweep_rows`). The plate impulse is the 2-D effective-gamma run
   (`sweep --iapetus-2d`). Each material's radiative rebound replaces the lossless one,
   `(1 + e_eff)/(1 + e_gamma)`, at the arrival speed `u_rms + w` (`sweep --iapetus-rebound`).
   Any stranded store comes off the expansion's share. The leverage is the parent's overtake law,
   integrated over the push with the measured `B(k) = J/(m w)`.
4. **The plate** (`plate_rows`). Plate mass follows the lossless spring, `m s = J tan(theta/4) /
   (2 theta f)`. The radial impulse profile grades the floor. A free plate given velocity
   `j(r)/sigma(r)` puts `1 - (sum sigma v)^2 / (sum sigma sum sigma v^2)` of the kick's energy into
   bending; for an axisymmetric impulsive kick that share is exact.
5. **Aim** (`aim_rows`). An off-axis burst is traced ballistically: an isotropic point source
   drifting onto the plate, with each ray's normal momentum returned with restitution. The 2-D
   code is axisymmetric and cannot do this. The result is read as Result 19's rocking share.
6. **Erosion** (`erosion_rows`). The inner ring's arriving energy times the wall's radiative share
   from the 1-D runs gives the face heat per pulse, charged at an ablation energy.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np

from puffsat import eos_water
from puffsat import recombination as rc
from puffsat.water_plate import spray_eos
from puffsat.xray_collision import charge
from puffsat.xray_collision.eos import EV, K_B, EosTable
from puffsat.xray_collision.materials import CARBON

# --- The case (parent S50) ------------------------------------------------------------------

BAG = 25.0  # [kg]
W_REF = 51.1e3  # [m/s], closing speed at 200 km on the reference transfer
W_SENS = (60.0e3, 68.0e3)
V_F = 10.916e3  # [m/s], the parent's eccentric-orbit speed at 200 km
K_PARENT, ETA_PARENT = 8.5, 0.775
CADENCE = 2.0  # [Hz]
G0 = 9.80665
STROKE = 8.0  # [m]
THETA = 2.0 * math.pi / 3.0  # the parent's spring phase (Result 17)
PULSES = 2000
PLATE_R = 9.5  # [m], the flown push plate's floor

WATER_RHO, ICE_RHO = 1000.0, 917.0
ARGON_SOLID_RHO, ARGON_LIQUID_RHO = 1616.0, 1395.0
#: The charges the 2-D runs carry (`IAPETUS_CHARGES`), and their tables.
CHARGE_RHO = {"argon": ARGON_SOLID_RHO, "ice": ICE_RHO}
TABLES = {"argon": Path("data/tables/argon.json"), "ice": Path("data/tables/water_jupiter.json")}
MASS_AMU = {"argon": 39.948, "ice": 18.015 / 3.0}
MIXTURES = {
    "argon": spray_eos.AtomMixture((("Ar", 1.0),)),
    "ice": spray_eos.AtomMixture((("H", 2.0), ("O", 1.0))),
}
#: Polyethylene, raced on the CH2 table (atoms only, so its bonds are charged apart, as Result 11)
#: for Result 22's material comparison. Kept out of the maps above, which the 2-D charges iterate.
PE_RHO = 950.0
RACE_RHO = CHARGE_RHO | {"pe": PE_RHO}
RACE_MIXTURES = MIXTURES | {"pe": spray_eos.AtomMixture((("C", 1.0), ("H", 2.0)))}

OUT_DIR = charge.OUT_DIR
JSONL = OUT_DIR / "iapetus_2d.jsonl"
JSONL_COARSE = OUT_DIR / "iapetus_2d_coarse.jsonl"
#: The follow-up: moderate spreads onto a short charge disk (`sweep --iapetus-2d spread`).
SPREAD = OUT_DIR / "iapetus_2d_spread.jsonl"
SPREAD_COARSE = OUT_DIR / "iapetus_2d_spread_coarse.jsonl"
#: The two-piece charge: a thin bumper ahead of a 2.5 m disk (`sweep --iapetus-2d bumper`).
BUMPER = OUT_DIR / "iapetus_2d_bumper.jsonl"
BUMPER_COARSE = OUT_DIR / "iapetus_2d_bumper_coarse.jsonl"
#: The gun's ice bullet fired back along the track (`sweep --iapetus-2d bullet`).
BULLET = OUT_DIR / "iapetus_2d_bullet.jsonl"
BULLET_COARSE = OUT_DIR / "iapetus_2d_bullet_coarse.jsonl"
#: The bullet's four deciding runs on the production grid (`--iapetus-2d bullet key`).
BULLET_KEY = OUT_DIR / "iapetus_2d_bullet_key.jsonl"
REBOUND_JSONL = OUT_DIR / "iapetus_rebound.jsonl"
#: The f(v) study's freeze-timing bracket on the heavy plate at 16-28 km/s (ADR-0026): water's
#: rebound with its composition frozen at the turnaround.
FROZEN_WATER = Path("data/results/frontier_frozen_heavyplate.csv")


# --- 1. Arrival state -----------------------------------------------------------------------


def sphere_diameter(mass: float, rho: float) -> float:
    return float(2.0 * (3.0 * mass / (4.0 * math.pi * rho)) ** (1.0 / 3.0))


def cloud_depth(radius: float, fragment: float, rho: float = ICE_RHO, mass: float = BAG) -> float:
    """Optical depth of a cloud of `mass` in fragments of radius `fragment` spread over a disk of
    `radius`: their summed cross-section over the disk, `3 sigma / (4 rho a)`."""
    sigma = mass / (math.pi * radius**2)
    return 3.0 * sigma / (4.0 * rho * fragment)


def cloud_radius_for_depth(
    fragment: float, depth: float = 1.0, rho: float = ICE_RHO, mass: float = BAG
) -> float:
    """The widest cloud whose fragments still cover its area `depth` times."""
    return math.sqrt(3.0 * mass / (4.0 * math.pi * rho * fragment * depth))


@dataclass(frozen=True)
class ArrivalRow:
    state: str
    fragment_m: float  # 0 for an intact bag
    cloud_radius_m: float
    areal_kg_m2: float
    depth: float  # fragments' cover of the cloud's area
    local_k_at_8_5: float  # charge each fragment sweeps per its own mass, k = 8.5 sheet
    release_lead_s_at_100_m_s: float  # how far ahead to release for this radius at 100 m/s spread
    couples: str


#: The paper's stopping column for an unbroken drop, `2.7 rho_w r` (`sec:icy_puffsat`): the drop
#: stops once it has swept `3 x 2.7 / 4` of its own mass.
STOP_SWEEP = 3.0 * 2.7 / 4.0
#: Fragment radii [m]: the paper's ice-wire grid cap, a millimetre, and two sprays.
FRAGMENTS = (2.5e-3, 1.0e-3, 1.0e-4, 1.0e-5)
CLOUDS = (0.2, 1.6, 4.0, PLATE_R)


def arrival_rows() -> list[ArrivalRow]:
    out = []
    for name, rho in (("intact, liquid", WATER_RHO), ("intact, ice", ICE_RHO)):
        d = sphere_diameter(BAG, rho)
        out.append(
            ArrivalRow(
                name,
                0.0,
                0.5 * d,
                BAG / (math.pi * (0.5 * d) ** 2),
                math.inf,
                K_PARENT,
                0.0,
                "compact: needs a compact charge in its path",
            )
        )
    for a in FRAGMENTS:
        for r in CLOUDS[1:]:
            tau = cloud_depth(r, a)
            out.append(
                ArrivalRow(
                    "released ahead (fragments)",
                    a,
                    r,
                    BAG / (math.pi * r**2),
                    tau,
                    K_PARENT * tau,
                    r / 100.0,
                    "as a sheet"
                    if tau >= 1.0
                    else (
                        "each stops on its own track only"
                        if K_PARENT * tau >= STOP_SWEEP
                        else "fragments punch through the charge"
                    ),
                )
            )
    return out


def sunlit_temperature(au: float, absorb_over_emit: float = 1.0) -> float:
    """A sphere's equilibrium temperature at `au` [K]: `(alpha S / (4 eps sigma))^(1/4)`."""
    s = 1361.0 / au**2
    return float((absorb_over_emit * s / (4.0 * 5.670374419e-8)) ** 0.25)


# --- 2. The chemistry -----------------------------------------------------------------------


@cache
def table(material: str) -> EosTable:
    """The argon or water equilibrium table as a fireball `EosTable` (as `charge.ch2_table`), or
    the CH2 table itself for `pe`."""
    if material == "pe":
        return charge.ch2_table()
    d = json.loads(TABLES[material].read_text())
    n = d["shape"]
    f = {k: np.array(v).reshape(n) for k, v in d["fields"].items()}
    mat = dataclasses.replace(CARBON, name=material, mass_amu=MASS_AMU[material])
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


def ion_energy(material: str, rho: float, temp: float) -> tuple[float, float]:
    """`(e_ion [J/kg], n_e)` in Saha equilibrium (atoms only, as `freeze.ion_energy`)."""
    s = spray_eos.state(RACE_MIXTURES[material], rho, temp)
    e = 0.0
    for name, n in s.n_element.items():
        e += n * float(np.dot(s.fractions[name], spray_eos._cumulative_ip(name)))
    return e / rho, s.n_e


#: Below this the water is all molecular at these densities; the race stops.
T_STOP = 2000.0
#: The rebound layer has delivered its push once it has thinned this far (nine crossing times).
REBOUND_THIN = 0.1


@dataclass(frozen=True)
class Race:
    material: str
    leg: str  # "fireball" (solid density to the plate) or "rebound" (off the stagnated plate)
    k: float
    w: float
    rho_start: float
    kt_start_ev: float
    rho_end: float
    t_end_k: float
    kinetic_share_end: float  # of the heat, as expansion (fireball leg only)
    ion_share_start: float  # of the heat (fireball) or of the stagnation energy (rebound)
    ion_share_frozen: float
    da_ion_min: float
    bond_share_start: float  # water only: of the heat
    bond_share_end_eq: float  # held at the end in equilibrium
    freeze_t_h: float  # [K] where H + OH + M falls behind (Da < 1), H-limited; nan if it never does
    stranded_h: float  # bond energy held there, of the heat
    freeze_t_oh: float  # the OH-limited edge
    stranded_oh: float
    da_bond_min_h: float  # least Da met while the bonds still return (T above T_STOP)
    da_bond_min_oh: float


def race(
    material: str,
    k: float,
    w: float = W_REF,
    plate_distance: float = 4.0,
    leg: str = "fireball",
    rebound_column: float = 0.4,
    bag: float = BAG,
) -> Race:
    """Race the ionization store and (water) the bonds against the expansion.

    `fireball`: the merged mass `(1 + k) BAG` starts at the charge's solid density (held at the
    table's ceiling, as Result 10) with the collision's heat `k/(1+k) w^2/2` per kg of bag, and
    expands as `freeze.freeze`'s homologous sphere until its edge reaches `plate_distance`.

    `rebound`: the gas arrives at `u_rms + w` with the fireball's end density, stagnates at six
    times that density (a strong reflected shock at gamma ~ 1.4) and then re-expands as a planar
    layer of `rebound_column` kg/m^2. That layer's density falls as `1 / (1 + t / tau)`, with
    `tau = column / (rho c_s)` its sound-crossing time.
    """
    tab = table(material)
    mass = (1.0 + k) * bag
    heat_bag = 0.5 * k / (1.0 + k) * w * w  # [J per kg of bag]
    heat = heat_bag / (1.0 + k)  # [J per kg of merged mass]
    ceiling = float(np.exp(tab.ln_rho[-1]))
    if leg == "fireball":
        rho = min(RACE_RHO[material], ceiling)
        e = heat
    else:
        f = race(material, k, w, plate_distance, bag=bag)
        u_rms = math.sqrt(2.0 * heat * f.kinetic_share_end)
        drift = w / (1.0 + k)
        rho_a = f.rho_end
        rho = 6.0 * rho_a
        e = heat * (1.0 - f.kinetic_share_end) + 0.5 * (u_rms + drift) ** 2
    temp = tab.temp(rho, e)
    e_start = e
    e_ion0, _ = ion_energy(material, rho, temp)
    x = e_ion0
    water = material == "ice"
    bond0 = eos_water.bond_energy_held(rho, temp) if water else 0.0
    da_ion_min = math.inf
    da_h_min = da_oh_min = math.inf
    freeze_h = freeze_oh = (math.nan, 0.0)
    rho_start, t_start = rho, temp
    r = (3.0 * mass / (4.0 * math.pi * rho)) ** (1.0 / 3.0)
    rdot, t = 0.0, 0.0
    tau_layer = 0.0
    if leg == "rebound":
        c_s = math.sqrt(1.4 * tab.state(rho, temp)[0] / rho)
        tau_layer = rebound_column / (rho * c_s)
    rho0 = rho
    while True:
        temp = tab.temp(rho, e)
        p = tab.state(rho, temp)[0]
        x_eq, n_e = ion_energy(material, rho, temp)
        if leg == "fireball":
            tau_exp = r / (3.0 * rdot) if rdot > 0.0 else math.inf
        else:
            tau_exp = tau_layer + t
        alpha = rc.three_body_recombination(temp, n_e) + rc.radiative_recombination(temp)
        tau_rec = 1.0 / (alpha * n_e) if n_e > 0.0 else math.inf
        if math.isfinite(tau_exp) and x > 1e-6 * e:
            da_ion_min = min(da_ion_min, tau_exp / tau_rec)
        if water and math.isfinite(tau_exp) and temp > T_STOP:
            comp = eos_water.composition(rho, temp)
            st = rc.freeze_station(
                time=t,
                temp=temp,
                rho=rho,
                n_e=comp.n_e,
                n_h=comp.n_h,
                n_oh=comp.n_oh,
                n_third=comp.n_total,
                tau_expansion=tau_exp,
                bond_energy_fraction=eos_water.bond_energy_held(rho, temp)
                / eos_water.FULL_ATOMIZATION_ENERGY,
            )
            held = eos_water.bond_energy_held(rho, temp) / (heat if leg == "fireball" else e_start)
            # Only a station that still has bonds to return, and is cool enough for water to be
            # re-forming, can strand them.
            if held > 1e-3 and temp < 6000.0:
                da_h_min = min(da_h_min, st.da_dissociation_h_limited)
                da_oh_min = min(da_oh_min, st.da_dissociation_oh_limited)
                if math.isnan(freeze_h[0]) and st.da_dissociation_h_limited < 1.0:
                    freeze_h = (temp, held)
                if math.isnan(freeze_oh[0]) and st.da_dissociation_oh_limited < 1.0:
                    freeze_oh = (temp, held)
        if leg == "fireball":
            edge = r * math.sqrt(5.0 / 3.0)
            if edge >= plate_distance or temp < T_STOP:
                break
            dt = 0.02 * r / max(rdot, math.sqrt(p / rho))
            rddot = 20.0 * math.pi / 3.0 * r**2 * p / mass
            dv = 4.0 * math.pi * r**2 * rdot * dt
            e -= p * dv / mass
            rdot += rddot * dt
            r += rdot * dt
            rho = mass / (4.0 / 3.0 * math.pi * r**3)
        else:
            if rho < REBOUND_THIN * rho0 or temp < T_STOP:
                break
            # Planar adiabatic re-expansion: rho = rho0 / (1 + t / tau), de = p / rho^2 drho.
            rho_new = rho * 0.97
            dt = (rho0 / rho_new - rho0 / rho) * tau_layer
            e += p / rho**2 * (rho_new - rho)
            rho = rho_new
        if math.isfinite(tau_rec):
            x += -(x - x_eq) * (1.0 - math.exp(-dt / tau_rec))
        t += dt
    # Shares are of the collision's heat on the fireball leg, of the stagnation energy after.
    base = heat if leg == "fireball" else e_start
    frozen = max(x - ion_energy(material, rho, temp)[0], 0.0) / base
    kinetic = 0.5 * rdot**2 * 3.0 / 5.0 / heat if leg == "fireball" else math.nan
    return Race(
        material=material,
        leg=leg,
        k=k,
        w=w,
        rho_start=rho_start,
        kt_start_ev=t_start * K_B / EV,
        rho_end=rho,
        t_end_k=temp,
        kinetic_share_end=kinetic,
        ion_share_start=e_ion0 / e_start,
        ion_share_frozen=frozen,
        da_ion_min=da_ion_min,
        bond_share_start=bond0 / base,
        bond_share_end_eq=(eos_water.bond_energy_held(rho, temp) / base if water else 0.0),
        freeze_t_h=freeze_h[0],
        stranded_h=freeze_h[1],
        freeze_t_oh=freeze_oh[0],
        stranded_oh=freeze_oh[1],
        da_bond_min_h=da_h_min,
        da_bond_min_oh=da_oh_min,
    )


def ideal_kinetic_share(gamma: float, rho_start: float, rho_end: float) -> float:
    """The share of a gamma gas's heat that a homologous sphere has turned into expansion by the
    time it thins from `rho_start` to `rho_end` (from rest, adiabatic): `1 - (rho_end /
    rho_start)^(gamma - 1)`."""
    return float(1.0 - (rho_end / rho_start) ** (gamma - 1.0))


def effective_gamma(r: Race) -> float:
    """The gamma whose ideal sphere converts the same share of its heat by the plate as the
    material does on its equilibrium table."""
    thermal = 1.0 - r.kinetic_share_end
    return 1.0 + math.log(thermal) / math.log(r.rho_end / r.rho_start)


# --- 3. The k sweep -------------------------------------------------------------------------


def _f(r: dict[str, object], key: str) -> float:
    return float(str(r[key]))


#: How far over the energy ceiling a lossless run may read before it is called unphysical.
CEILING_SLACK = 1.02


def physical(r: dict[str, object]) -> bool:
    """A run that kept its exhaust moving away from the plate and stayed under the perfect plate's
    `m v (1 + sqrt(1 + k))`: the bag's kinetic energy all leaving as one axial jet, which a lossless
    run cannot beat. Runs that break it created energy numerically -- a thin, light target slab
    on the coarse grid (J/mv 605, exhaust moving backward) and a short ice column on the flat plate
    (J/mv 4.8 against a 3.24 ceiling) -- and are dropped, not used."""
    k = _f(r, "k") + (_f(r, "bumper_k") if "bumper_k" in r else 0.0)
    ceiling = 1.0 + math.sqrt(1.0 + k)
    return _f(r, "exit_direction") > 0.0 and _f(r, "impulse_per_mv") <= CEILING_SLACK * ceiling


@cache
def runs(path: Path = JSONL) -> tuple[dict[str, object], ...]:
    if not path.exists():
        return ()
    rows = charge.load_jsonl(path)
    for r in rows:
        if not physical(r):
            print(f"  dropped an unphysical run in {path.name}: {r['case']} k {r['k']}")
    return tuple(r for r in rows if physical(r))


def bumper_record(r: dict[str, object]) -> dict[str, object]:
    """A bumper run in the main grid's schema. Its `k` is the whole charge the ship carries,
    bumper and disk together, and its case names the bumper."""
    kb, scale = _f(r, "bumper_k"), _f(r, "gap_scale")
    if _f(r, "bumper_k") < 0.05:
        what = "narrow bullet" if _f(r, "bumper_radius") < _f(r, "bag_radius") else "wide bullet"
    else:
        what = "bumper"
    return {
        **r,
        "case": f"{r['case']}, {what} {kb:g} at x{scale:g} gap",
        "kind": "bumper" if _f(r, "k") > 0.0 else "miss",
        "charge": "mixed",
        "gamma": 1.4,
        "k": _f(r, "k") + kb,
        "sheet_radius": _f(r, "charge_radius"),
        "charge_length": 0.0,
    }


def all_runs(fine: bool = True) -> tuple[dict[str, object], ...]:
    """The main grid, the spread follow-up and the bumper runs, on the production (or the
    coarse) grid."""
    bumpers = tuple(
        bumper_record(r)
        for path in ((BUMPER, BULLET, BULLET_KEY) if fine else (BUMPER_COARSE, BULLET_COARSE))
        for r in runs(path)
    )
    return runs(JSONL if fine else JSONL_COARSE) + runs(SPREAD if fine else SPREAD_COARSE) + bumpers


@cache
def rebound_rows() -> tuple[dict[str, object], ...]:
    return tuple(
        r for r in charge.load_jsonl(REBOUND_JSONL) if r["converged"] and _f(r, "length") == 1.0
    )


def rebound(material: str, v: float, column: float) -> tuple[float, float]:
    """`(e_eff, wall share)` of argon or water at arrival speed `v` and arriving column [kg/m^2]
    (1 m slugs, interpolated in v and log column, held at the grid's edges). The wall share is the
    energy radiated into the wall over the slug's incident kinetic energy."""
    name = "argon" if material == "argon" else "water"
    rows = [r for r in rebound_rows() if r["material"] == name]
    vs = sorted({_f(r, "v") for r in rows})
    cols = sorted({_f(r, "rho_impact") for r in rows})

    def grid(key: str) -> np.ndarray:
        out = np.full((len(vs), len(cols)), np.nan)
        for r in rows:
            i, j = vs.index(_f(r, "v")), cols.index(_f(r, "rho_impact"))
            val = _f(r, key)
            if key == "loss_radiative_wall":
                val /= 0.5 * _f(r, "rho_impact") * _f(r, "v") ** 2
            out[i, j] = val
        return out

    lc = np.log(cols)
    x = float(np.clip(math.log(column), lc[0], lc[-1]))
    res = []
    for key in ("e_eff", "loss_radiative_wall"):
        g = grid(key)
        per_v = [float(np.interp(x, lc[~np.isnan(row)], row[~np.isnan(row)])) for row in g]
        res.append(float(np.interp(v, vs, per_v)))
    return res[0], res[1]


@cache
def frozen_water_delta(v: float) -> float:
    """The f(v) study's freeze-timing bracket on the heavy plate: `e_eq - e_frozen_rebound` at its
    densest bag (0.6 kg/m^3), interpolated in v and held at 16-28 km/s."""
    rows = [r for r in csv.DictReader(FROZEN_WATER.open()) if float(r["rho_impact"]) == 0.6]
    vs = [float(r["v"]) for r in rows]
    ds = [float(r["delta_frozen"]) for r in rows]
    return float(np.interp(v, vs, ds))


def gamma_factor(case: str, k: float, gamma: float, path: Path = JSONL) -> float:
    """`(J/cr)(gamma) / (J/cr)(1.4)` of the point fireball: what a softer or stiffer gas does to
    the plate impulse before the rebound swap. Read on the flat 2.7 m plate for flat plates and
    on the 2 F paraboloid for bowls, linear in gamma and in ln k between the bracket's runs."""
    ref = "flat + 2 m skirt, 2.7 m off" if case.startswith("flat") else "paraboloid 0.50 at 2 F"
    pts = [r for r in runs(path) if r["kind"] == "point" and r["case"] == ref]
    ks = sorted({_f(r, "k") for r in pts if _f(r, "gamma") != 1.4})
    gs = sorted({_f(r, "gamma") for r in pts})
    if not ks:
        return 1.0

    def q(kk: float, g: float) -> float:
        r = next(r for r in pts if _f(r, "k") == kk and _f(r, "gamma") == g)
        return _f(r, "impulse_per_mv") / _f(r, "confined_ratio")

    per_k = [float(np.interp(gamma, gs, [q(kk, g) / q(kk, 1.4) for g in gs])) for kk in ks]
    return float(np.interp(math.log(k), np.log(ks), per_k))


def leverage(b: float, k: float, w: float = W_REF) -> float:
    """Kilograms pushed to `V_F` per kilogram of returned bag: the overtake law integrated over
    the push, `M_f/m = k / ((w/(w - v_f))^(k/B) - 1)` (`B/ln(...)` as k -> 0)."""
    x = math.log(w / (w - V_F))
    if k == 0.0:
        return b / x
    return k / math.expm1(k / b * x)


@dataclass(frozen=True)
class SweepRow:
    case: str
    kind: str
    charge: str  # the charge's material: argon or ice (water)
    sheet_radius: float  # [m], released-ahead cloud and charge disk; 0 for compact and point
    k: float
    w: float
    height: float
    charge_length: float
    gamma_eff: float
    gamma_factor: float
    arrival_speed: float
    e_eff: float
    wall_share: float
    impulse_per_mv_lossless: float  # the 2-D run at gamma 1.4
    impulse_per_mv: float  # B: gamma-corrected, rebound swapped, bonds returned
    impulse_per_mv_frozen: float  # water: the f(v) study's frozen-rebound bracket
    eta: float
    eta_frozen: float
    leverage: float
    leverage_frozen: float
    peak_pa: float  # [Pa] for a 25 kg bag
    max_bag_maraging: float  # [kg per pulse] at 2.3 GPa
    exit_direction: float
    exit_uniformity: float


#: Face allowables [Pa] (Result 13, `docs/spray_plate_steel_face_limits.md`).
MARAGING, HY100 = 2.3e9, 1.3e9


def sweep_case(r: dict[str, object], material: str, w: float = W_REF) -> SweepRow:
    k = _f(r, "k")
    fire = race(material, k, w)
    g_eff = effective_gamma(fire)
    gf = gamma_factor(str(r["case"]), k, g_eff)
    u = w * math.sqrt(k) / (1.0 + k)
    v_arr = u + w / (1.0 + k)
    column = (1.0 + k) * BAG / (2.0 * math.pi * PLATE_R**2)
    e, wall = rebound(material, v_arr, column)
    lossless = _f(r, "impulse_per_mv")
    cr = _f(r, "confined_ratio")
    b = lossless * gf * (1.0 + e) / cr
    b_frozen = b
    if material == "ice":
        b_frozen = lossless * gf * (1.0 + e - frozen_water_delta(v_arr)) / cr
    root = math.sqrt(1.0 + k)
    peak = _f(r, "peak_per_mv2") * BAG * w * w
    return SweepRow(
        case=str(r["case"]),
        kind=str(r["kind"]),
        charge=material,
        sheet_radius=_f(r, "sheet_radius"),
        k=k,
        w=w,
        height=_f(r, "height"),
        charge_length=_f(r, "charge_length"),
        gamma_eff=g_eff,
        gamma_factor=gf,
        arrival_speed=v_arr,
        e_eff=e,
        wall_share=wall,
        impulse_per_mv_lossless=lossless,
        impulse_per_mv=b,
        impulse_per_mv_frozen=b_frozen,
        eta=(b - 1.0) / root,
        eta_frozen=(b_frozen - 1.0) / root,
        leverage=leverage(b, k, w),
        leverage_frozen=leverage(b_frozen, k, w),
        peak_pa=peak,
        max_bag_maraging=MARAGING / peak * BAG,
        exit_direction=_f(r, "exit_direction"),
        exit_uniformity=_f(r, "exit_uniformity"),
    )


def sweep_rows(fine: bool = True) -> list[SweepRow]:
    out = []
    for r in all_runs(fine):
        if _f(r, "gamma") != 1.4 or r["kind"] == "miss":
            continue
        mats = [str(r["charge"])] if r["kind"] == "compact" else ["argon", "ice"]
        for mat in mats:
            out.append(sweep_case(r, mat))
            if r["kind"] in ("compact", "sheet") and _f(r, "k") in (8.5, 20.0):
                out += [sweep_case(r, mat, w) for w in W_SENS]
    return out


def best_k(rows: list[SweepRow], frozen: bool = False) -> tuple[float, float]:
    """`(k, leverage)` at the best slug ratio: `B` interpolated in ln k over the rows (one plate,
    one charge, one speed) and the leverage maximized on a fine grid between the first and last
    run."""
    rows = sorted(rows, key=lambda c: c.k)
    ks = np.array([c.k for c in rows])
    bs = np.array([c.impulse_per_mv_frozen if frozen else c.impulse_per_mv for c in rows])
    # Only between the runs: below the first one `np.interp` would hold its B fixed.
    fine = np.linspace(ks[0], ks[-1], 600)
    bf = np.interp(np.log(fine), np.log(ks), bs)
    lev = np.array(
        [leverage(float(b), float(kk), rows[0].w) for b, kk in zip(bf, fine, strict=True)]
    )
    i = int(np.argmax(lev))
    return float(fine[i]), float(lev[i])


# --- 4. The plate ---------------------------------------------------------------------------


def spring_mass(pulse: float, stroke: float = STROKE, theta: float = THETA) -> float:
    """[kg] The paper's lossless-spring plate (Result 17): `m s = J tan(theta/4) / (2 theta f)`."""
    return pulse * math.tan(theta / 4.0) / (2.0 * theta * CADENCE * stroke)


def ring_profile(r: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """`(ring edges [m], impulse per projected area per m v [1/m^2])` inside the floor's radius."""
    prof = r["profile"]
    assert isinstance(prof, list)
    rs = np.array([0.0] + [float(p[0]) for p in prof])
    cum = np.array([0.0] + [float(p[1]) for p in prof])
    keep = rs <= PLATE_R + 1e-9
    rs, cum = rs[keep], cum[keep]
    return rs, np.diff(cum) / (math.pi * np.diff(rs**2))


def ringing_share(areal: np.ndarray, j: np.ndarray, area: np.ndarray) -> float:
    """The share of an impulsive axisymmetric kick's energy that is not rigid translation, for
    ring masses `areal * area` given velocities `j / areal`. A free body's only axisymmetric rigid
    mode is translation, so this is exactly the energy left in bending."""
    m = areal * area
    v = j / areal
    return 1.0 - float(np.sum(m * v)) ** 2 / (float(np.sum(m)) * float(np.sum(m * v * v)))


#: The parent's structural floor where the load is light ("tens of kg/m^2"; charge.py).
AREAL_FLOOR = charge.SURVIVAL_AREAL
STEEL_RHO = 7850.0


@dataclass(frozen=True)
class PlateRow:
    case: str
    kind: str
    sheet_radius: float
    charge: str
    k: float
    impulse_per_mv: float
    pulse_ns: float  # [N s] at the start of the push (closing speed w)
    pulse_end_ns: float  # at its end (w - v_f)
    ship_t_1g: float  # [t] the ship that this pulse holds at 1 g at 2 Hz (start of push)
    spring_t: float  # [t] plate on the lossless spring at an 8 m stroke
    centre_over_mean: float  # impulse per area inside 0.75 m over the floor's mean
    ringing_uniform: float  # bending share of the kick on an even floor
    ringing_graded: float  # on a floor graded to the load, with the 30 kg/m^2 floor
    graded_centre_kg_m2: float
    graded_centre_steel_mm: float
    graded_t: float  # [t] the graded plate (the spring mass, or more if the floor binds)
    rock_per_m_uniform: float  # edge kick per heave kick, per metre of impulse-centroid offset
    rock_per_m_graded: float
    peak_gpa: float
    max_bag_maraging: float
    exit_direction: float


def plate_case(s: SweepRow, rec: dict[str, object]) -> PlateRow:
    pulse = s.impulse_per_mv * BAG * s.w
    m_p = spring_mass(pulse)
    edges, j = ring_profile(rec)
    area = math.pi * np.diff(edges**2)
    mid = 0.5 * (edges[1:] + edges[:-1])
    mean = float(np.sum(j * area) / np.sum(area))
    inner = edges[1:] <= 0.75 + 1e-9
    centre = float(np.sum(j[inner] * area[inner]) / np.sum(area[inner]))
    total_area = float(np.sum(area))
    uniform = np.full_like(j, m_p / total_area)
    # Graded: sigma = max(floor, c j) with c set so the floor's mass is the spring mass.
    floor = AREAL_FLOOR
    if floor * total_area >= m_p:
        graded = np.full_like(j, floor)
    else:
        lo, hi = 0.0, 10.0 * m_p / (mean * total_area)
        for _ in range(200):
            c = 0.5 * (lo + hi)
            if float(np.sum(np.maximum(floor, c * j) * area)) > m_p:
                hi = c
            else:
                lo = c
        graded = np.maximum(floor, 0.5 * (lo + hi) * j)
    m_g = float(np.sum(graded * area))
    # Rocking turns about a diameter: half the polar moment of an axisymmetric floor.
    inertia_u = 0.5 * float(np.sum(uniform * area * mid**2))
    inertia_g = 0.5 * float(np.sum(graded * area * mid**2))
    return PlateRow(
        case=s.case,
        kind=s.kind,
        sheet_radius=s.sheet_radius,
        charge=s.charge,
        k=s.k,
        impulse_per_mv=s.impulse_per_mv,
        pulse_ns=pulse,
        pulse_end_ns=pulse * (s.w - V_F) / s.w,
        ship_t_1g=pulse * CADENCE / G0 / 1e3,
        spring_t=m_p / 1e3,
        centre_over_mean=centre / mean,
        ringing_uniform=ringing_share(uniform, j, area),
        ringing_graded=ringing_share(graded, j, area),
        graded_centre_kg_m2=float(graded[0]),
        graded_centre_steel_mm=1e3 * float(graded[0]) / STEEL_RHO,
        graded_t=m_g / 1e3,
        rock_per_m_uniform=PLATE_R * m_p / inertia_u,
        rock_per_m_graded=PLATE_R * m_g / inertia_g,
        peak_gpa=s.peak_pa / 1e9,
        max_bag_maraging=s.max_bag_maraging,
        exit_direction=s.exit_direction,
    )


@dataclass(frozen=True)
class MissRow:
    """The bag strikes the bumper but misses the disk (or, with no bumper, strikes the plate
    bare, intact or as a cloud it spread itself): what reaches the face, for a 25 kg bag at the
    reference closing speed."""

    case: str
    grid: str  # "x1" production, "x0.5" the coarse grid
    bag_radius: float  # [m], 0.2 intact, or the self-spread cloud's
    bumper_k: float  # 0: the bare bag
    bumper_radius_m: float  # a bumper a cell wider than the bag, or a gun's bullet
    impulse_mn_s: float  # [MN s], lossless gamma gas
    peak_gpa: float
    peak_radius_m: float
    inner_energy_mj_m2: float  # arriving inside 0.6 m
    mid_energy_mj_m2: float  # arriving on the 1-2 m ring
    holds_maraging: bool


def miss_rows() -> list[MissRow]:
    """Both grids, tagged: the bullet's misses have production runs, the bumpers' only coarse."""
    out = []
    tagged = [("x1", p) for p in (BUMPER, BULLET, BULLET_KEY)] + [
        ("x0.5", p) for p in (BUMPER_COARSE, BULLET_COARSE)
    ]
    for grid, r in ((g, row) for g, p in tagged for row in runs(p)):
        if _f(r, "k") != 0.0:
            continue
        bands = r["bands"]
        assert isinstance(bands, list)
        scale = BAG * W_REF * W_REF
        peak = _f(r, "peak_per_mv2") * scale
        out.append(
            MissRow(
                case=str(r["case"]),
                grid=grid,
                bag_radius=_f(r, "bag_radius"),
                bumper_k=_f(r, "bumper_k"),
                bumper_radius_m=_f(r, "bumper_radius") if "bumper_radius" in r else 0.0,
                impulse_mn_s=_f(r, "impulse_per_mv") * BAG * W_REF / 1e6,
                peak_gpa=peak / 1e9,
                peak_radius_m=_f(r, "peak_radius"),
                inner_energy_mj_m2=float(bands[0]["energy_in"]) * scale / 1e6,
                mid_energy_mj_m2=float(bands[2]["energy_in"]) * scale / 1e6,
                holds_maraging=peak <= MARAGING,
            )
        )
    return out


# --- 5. Aim: a ballistic off-axis burst -----------------------------------------------------


@dataclass(frozen=True)
class Shape:
    """A plate as the ray tracer sees it: a paraboloid `z = r^2 / (4 F)` (or a flat floor) out to
    `PLATE_R`, with an optional straight skirt."""

    d_over_d: float
    skirt: float

    @property
    def focus(self) -> float:
        return PLATE_R / (8.0 * self.d_over_d) if self.d_over_d > 0.0 else math.inf

    def floor_z(self, r: float) -> float:
        return r * r / (4.0 * self.focus) if self.d_over_d > 0.0 else 0.0

    @property
    def rim_z(self) -> float:
        return self.floor_z(PLATE_R)


def _hit(shape: Shape, x0: np.ndarray, d: np.ndarray) -> tuple[float, np.ndarray] | None:
    """The first wall the ray `x0 + s d` (s > 0) meets, with the wall's normal into the plate."""
    best: tuple[float, np.ndarray] | None = None
    a = 1.0 / (4.0 * shape.focus) if shape.d_over_d > 0.0 else 0.0
    # Floor: z = a (x^2 + y^2).
    qa = a * (d[0] ** 2 + d[1] ** 2)
    qb = 2.0 * a * (x0[0] * d[0] + x0[1] * d[1]) - d[2]
    qc = a * (x0[0] ** 2 + x0[1] ** 2) - x0[2]
    roots: list[float] = []
    if abs(qa) < 1e-14:
        if abs(qb) > 1e-14:
            roots = [-qc / qb]
    else:
        disc = qb * qb - 4.0 * qa * qc
        if disc >= 0.0:
            sq = math.sqrt(disc)
            roots = [(-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa)]
    for s in sorted(roots):
        if s > 1e-9:
            p = x0 + s * d
            if math.hypot(p[0], p[1]) <= PLATE_R:
                n = np.array([2.0 * a * p[0], 2.0 * a * p[1], -1.0])
                best = (s, n / np.linalg.norm(n))
                break
    # Skirt: x^2 + y^2 = R^2 between the rim and rim + skirt.
    if shape.skirt > 0.0:
        ca = d[0] ** 2 + d[1] ** 2
        cb = 2.0 * (x0[0] * d[0] + x0[1] * d[1])
        cc = x0[0] ** 2 + x0[1] ** 2 - PLATE_R**2
        disc = cb * cb - 4.0 * ca * cc
        if ca > 1e-14 and disc >= 0.0:
            s = (-cb + math.sqrt(disc)) / (2.0 * ca)
            if s > 1e-9 and (best is None or s < best[0]):
                p = x0 + s * d
                if shape.rim_z <= p[2] <= shape.rim_z + shape.skirt:
                    best = (s, np.array([p[0], p[1], 0.0]) / PLATE_R)
    return best


@dataclass(frozen=True)
class AimRow:
    case: str
    k: float
    miss_m: float
    axial: float  # plate impulse per M u_rms, centred burst (the ballistic model's capture)
    tilt_deg_per_m: float  # impulse direction error per metre of lateral miss (toward the miss)
    lateral_per_m: float  # lateral impulse over axial, per metre
    centroid_per_m: float  # the axial impulse's centroid shift per metre of miss
    cm_height: float  # [m] the uniform shell's centre of mass above the vertex
    lever_per_m: float  # angular impulse about the centre of mass over axial impulse, per metre
    rock_share_per_cm: float  # Result 19's edge kick per heave kick, `4 lever / R`, per cm


def trace(
    shape: Shape,
    height: float,
    miss: float,
    k: float,
    e: float = 0.6,
    n: int = 6000,
    away: bool = False,
) -> tuple[np.ndarray, np.ndarray, float]:
    """`(impulse, angular impulse about the vertex, sum of x times axial impulse)` from a point
    burst at `(miss, 0, height)`
    drifting onto the plate at `u_rms/sqrt(k)` (`away`: off it, as on the Sun bowl, whose disk
    enters through the bore), in units of the burst's `M u_rms`. Rays carry
    equal mass in Fibonacci-sphere directions. Each hit returns `(1 + e)` of the normal
    momentum and keeps the tangential; a ray may strike up to four times."""
    imp = np.zeros(3)
    ang = np.zeros(3)
    xmom = 0.0
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        rad = math.sqrt(1.0 - z * z)
        th = golden * i
        drift = (1.0 if away else -1.0) / math.sqrt(k)
        p = np.array([rad * math.cos(th), rad * math.sin(th), z + drift]) / n
        x = np.array([miss, 0.0, height])
        for _ in range(4):
            speed = float(np.linalg.norm(p))
            if speed == 0.0:
                break
            h = _hit(shape, x, p / speed)
            if h is None:
                break
            s, nrm = h
            x = x + s * p / speed
            pn = float(np.dot(p, nrm))
            if pn <= 0.0:
                break
            kick = (1.0 + e) * pn * nrm
            imp += kick
            ang += np.cross(x, kick)
            xmom -= x[0] * kick[2]
            p = p - kick
    return imp, ang, xmom


#: The ballistic plates: `(label, d/D, skirt, height above the vertex)` (the 2-D runs' set).
AIM_SHAPES = (
    ("flat + 2 m skirt, 2.7 m off", 0.0, 2.0, 2.7),
    ("flat + 2 m skirt, 6 m off", 0.0, 2.0, 6.0),
    ("bowl 0.30 at focus", 0.30, 0.0, 1.0),
    ("paraboloid 0.50 at 1.5 F", 0.50, 0.0, 1.5),
    ("paraboloid 0.50 at 2 F", 0.50, 0.0, 2.0),
    ("paraboloid 0.75 at 2 F", 0.75, 0.0, 2.0),
)
AIM_MISS = 0.1  # [m], Result 19's 1-sigma aim per axis


def shell_cm_height(shape: Shape, n: int = 2000) -> float:
    """[m] The centre of mass of the floor (uniform areal density) above the vertex. The skirt
    hangs on the craft's frame (Result 17), so it is left out."""
    r = (np.arange(n) + 0.5) * PLATE_R / n
    if shape.d_over_d == 0.0:
        return 0.0
    slope = r / (2.0 * shape.focus)
    da = 2.0 * math.pi * r * np.sqrt(1.0 + slope**2)
    return float(np.sum(da * r * r / (4.0 * shape.focus)) / np.sum(da))


def aim_rows(k: float = K_PARENT) -> list[AimRow]:
    out = []
    for label, dd, skirt, h in AIM_SHAPES:
        shape = Shape(dd, skirt)
        height = h * shape.focus if dd > 0.0 else h
        j0, _, _ = trace(shape, height, 0.0, k)
        j1, a1, xmom = trace(shape, height, AIM_MISS, k)
        axial = float(-j1[2])
        # About the y axis through the centre of mass: the vertex moment, less the lateral
        # impulse's arm. `a_y = sum (z Jx - x Jz)`, and the plate turns the way the axial
        # centroid's offset turns it.
        z_cm = shell_cm_height(shape)
        a_cm = float(a1[1]) - z_cm * float(j1[0])
        out.append(
            AimRow(
                case=label,
                k=k,
                miss_m=AIM_MISS,
                axial=float(-j0[2]),
                tilt_deg_per_m=math.degrees(math.atan2(float(j1[0]), axial)) / AIM_MISS,
                lateral_per_m=float(j1[0]) / axial / AIM_MISS,
                centroid_per_m=xmom / axial / AIM_MISS,
                cm_height=z_cm,
                lever_per_m=a_cm / axial / AIM_MISS,
                rock_share_per_cm=4.0 * abs(a_cm) / axial / AIM_MISS * 0.01 / PLATE_R,
            )
        )
    return out


# --- 6. Face erosion ------------------------------------------------------------------------

#: Maraging steel's face [SI] and its removal energy: melted off (to the melt, plus latent heat)
#: and boiled off (to the boil, plus both latent heats).
FACE_K, FACE_RHO, FACE_CP = 20.0, 8000.0, 450.0
FACE_MELT_RISE = 1400.0  # [K] from ~300 K to ~1700 K
Q_MELT, Q_BOIL = 1.0e6, 7.5e6  # [J/kg]
#: The sprayed water layer boils off in the steel's place: from ~300 K to steam, 2.6 MJ/kg.
Q_SPRAY = 2.6e6
#: A graphite face eroded instead (Result 10's collar bracket: a third of sublimation to all of it).
Q_GRAPHITE = (20.0e6, 59.0e6)


@dataclass(frozen=True)
class ErosionRow:
    case: str
    kind: str
    sheet_radius: float
    charge: str
    k: float
    pulse_us: float
    inner_energy_mj_m2: float  # arriving on the 0-0.6 m ring, per pulse
    wall_share: float  # of the arriving energy, radiated into the wall (1-D)
    face_heat_mj_m2: float
    melt_threshold_mj_m2: float  # the heat a semi-infinite face takes before its surface melts
    eroded_um_per_pulse: tuple[float, float]  # boiled off / melted off
    eroded_mm_per_push: tuple[float, float]
    spray_kg_per_pulse: float  # water sprayed on the face to boil off the heat over every ring
    hot_radius_m: float  # out to where the face heat exceeds the melt threshold
    leverage: float  # the sweep's, charge only
    leverage_with_spray: float  # the spray carried too, as charge that adds no impulse
    graphite_kg_per_pulse: tuple[float, float]  # the same heat taken by eroding graphite
    leverage_with_graphite: float  # the graphite carried, at the conservative 20 MJ/kg
    arriving_over_ke: float  # the rings' arriving energy over the bag's kinetic energy
    heat_valid: bool  # False where that exceeds 1: the tally has counted recirculating gas


def erosion_case(s: SweepRow, rec: dict[str, object]) -> ErosionRow:
    bands = rec["bands"]
    assert isinstance(bands, list)
    inner = bands[0]
    assert isinstance(inner, dict)
    energy = float(inner["energy_in"]) * BAG * s.w * s.w  # [J/m^2]
    u = s.w * math.sqrt(s.k) / (1.0 + s.k)
    t = max(s.height, 1.0) / u
    heat = energy * s.wall_share
    threshold = FACE_MELT_RISE * math.sqrt(math.pi * FACE_K * FACE_RHO * FACE_CP * t) / 2.0
    excess = max(heat - threshold, 0.0)
    per = (excess / Q_BOIL / FACE_RHO, excess / Q_MELT / FACE_RHO)
    # The sprayed layer: on every ring whose heat beats the threshold, enough water to boil the
    # excess away instead of the steel (`Q_SPRAY`). The rings are projected; the last is clipped
    # to the floor's radius.
    spray, hot = 0.0, 0.0
    # The rings count gas moving toward the wall at every step, so gas that sloshes in the bowl
    # is counted again each time. Where the total passes the bag's own kinetic energy (1/2 in
    # these units) the heat cannot be read from it, and the face-heat columns are left empty.
    arriving = 0.0
    for b in bands:
        assert isinstance(b, dict)
        lo, hi = float(b["r_lo"]), min(float(b["r_hi"]), PLATE_R)
        if hi > lo:
            arriving += math.pi * (hi * hi - lo * lo) * float(b["energy_in"])
    arriving /= 0.5
    valid = arriving <= 1.0
    for b in bands:
        assert isinstance(b, dict)
        lo, hi = float(b["r_lo"]), min(float(b["r_hi"]), PLATE_R)
        if hi <= lo:
            continue
        q = float(b["energy_in"]) * BAG * s.w * s.w * s.wall_share
        if q > threshold:
            spray += math.pi * (hi * hi - lo * lo) * (q - threshold) / Q_SPRAY
            hot = hi
    return ErosionRow(
        case=s.case,
        kind=s.kind,
        sheet_radius=s.sheet_radius,
        charge=s.charge,
        k=s.k,
        pulse_us=1e6 * t,
        inner_energy_mj_m2=energy / 1e6,
        wall_share=s.wall_share,
        face_heat_mj_m2=heat / 1e6,
        melt_threshold_mj_m2=threshold / 1e6,
        eroded_um_per_pulse=(1e6 * per[0], 1e6 * per[1]),
        eroded_mm_per_push=(1e3 * per[0] * PULSES, 1e3 * per[1] * PULSES),
        spray_kg_per_pulse=spray if valid else math.nan,
        hot_radius_m=hot if valid else math.nan,
        leverage=s.leverage,
        leverage_with_spray=(
            leverage(s.impulse_per_mv, s.k + spray / BAG, s.w) if valid else math.nan
        ),
        graphite_kg_per_pulse=(
            (spray * Q_SPRAY / Q_GRAPHITE[0], spray * Q_SPRAY / Q_GRAPHITE[1])
            if valid
            else (math.nan, math.nan)
        ),
        leverage_with_graphite=(
            leverage(s.impulse_per_mv, s.k + spray * Q_SPRAY / Q_GRAPHITE[0] / BAG, s.w)
            if valid
            else math.nan
        ),
        arriving_over_ke=arriving,
        heat_valid=valid,
    )


# --- Output ---------------------------------------------------------------------------------

OUT_ARRIVAL = OUT_DIR / "iapetus_arrival.csv"
OUT_RACE = OUT_DIR / "iapetus_race.csv"
OUT_SWEEP = OUT_DIR / "iapetus_sweep.csv"
OUT_BEST = OUT_DIR / "iapetus_best_k.csv"
OUT_PLATE = OUT_DIR / "iapetus_plate.csv"
OUT_AIM = OUT_DIR / "iapetus_aim.csv"
OUT_EROSION = OUT_DIR / "iapetus_erosion.csv"
OUT_CONVERGENCE = OUT_DIR / "iapetus_convergence.csv"
OUT_MISS = OUT_DIR / "iapetus_miss.csv"


def _write(path: Path, rows: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


@dataclass(frozen=True)
class BestRow:
    case: str
    kind: str
    charge: str
    sheet_radius: float
    w: float
    best_k: float
    leverage: float
    best_k_frozen: float
    leverage_frozen: float
    leverage_at_8_5: float


@dataclass(frozen=True)
class ConvergenceRow:
    case: str
    kind: str
    charge: str
    sheet_radius: float
    k: float
    impulse_per_mv: float
    impulse_per_mv_coarse: float
    peak_per_mv2: float
    peak_per_mv2_coarse: float


def _record(
    case: str, kind: str, mat: str, k: float, sheet: float, fine: bool = True
) -> dict[str, object]:
    return next(
        r
        for r in all_runs(fine)
        if r["case"] == case
        and r["kind"] == kind
        and (kind != "compact" or r["charge"] == mat)
        and _f(r, "k") == k
        and _f(r, "sheet_radius") == sheet
        and _f(r, "gamma") == 1.4
    )


def main() -> None:
    arr = arrival_rows()
    for ar in arr:
        print(
            f"  arrival {ar.state:<28} a {1e3 * ar.fragment_m:6.3f} mm R {ar.cloud_radius_m:4.2f}"
            f" m | {ar.areal_kg_m2:7.2f} kg/m^2, cover {ar.depth:8.3g}, local k"
            f" {ar.local_k_at_8_5:7.3g} | {ar.couples}"
        )
    print(
        f"  sunlit sphere: {sunlit_temperature(1.0):.0f} K at 1 AU, {sunlit_temperature(9.58):.0f}"
        " K at Saturn (black, equal absorptance and emittance)"
    )
    races = [
        race(m, k, leg=leg)
        for m in ("argon", "ice")
        for k in (4.0, 8.5, 15.0, 20.0, 30.0)
        for leg in ("fireball", "rebound")
    ]
    for rr in races:
        print(
            f"  race {rr.material:<5} {rr.leg:<8} k {rr.k:4.1f}: {rr.kt_start_ev:5.2f} eV at"
            f" {rr.rho_start:6.1f} -> {rr.t_end_k:6.0f} K at {rr.rho_end:7.3g} kg/m^3"
            + (
                f", expansion {rr.kinetic_share_end:.3f} (gamma_eff {effective_gamma(rr):.3f})"
                if rr.leg == "fireball"
                else ""
            )
            + f" | ions {rr.ion_share_start:.3f} -> frozen {rr.ion_share_frozen:.1e}"
            f" (Da_min {rr.da_ion_min:.1e})"
            + (
                f" | bonds {rr.bond_share_start:.3f} -> {rr.bond_share_end_eq:.3f} held;"
                f" Da_min {rr.da_bond_min_h:.1e} / {rr.da_bond_min_oh:.1e} (H / OH-limited),"
                f" freezes {rr.freeze_t_h:.0f} / {rr.freeze_t_oh:.0f} K"
                if rr.material == "ice"
                else ""
            )
        )
    sweep = sweep_rows()
    for s in sweep:
        print(
            f"  {s.case:<28} {s.kind:<7} R {s.sheet_radius:3.1f} {s.charge:<5} k {s.k:4.1f}"
            f" w {s.w / 1e3:4.1f}"
            f" | g {s.gamma_eff:.2f} x{s.gamma_factor:.3f} e {s.e_eff:.2f} | J/mv"
            f" {s.impulse_per_mv_lossless:5.2f} -> {s.impulse_per_mv:5.2f}"
            f" ({s.impulse_per_mv_frozen:5.2f}) eta {s.eta:.3f} ({s.eta_frozen:.3f})"
            f" | lev {s.leverage:5.2f} ({s.leverage_frozen:5.2f})"
            f" | peak {s.peak_pa / 1e9:5.2f} GPa,"
            f" <= {s.max_bag_maraging:5.1f} kg | dir {s.exit_direction:.2f}"
        )
    best = []
    groups = sorted({(s.case, s.kind, s.sheet_radius, s.charge, s.w) for s in sweep})
    for case, kind, rs, mat, w in groups:
        rows = [
            s
            for s in sweep
            if (s.case, s.kind, s.sheet_radius, s.charge, s.w) == (case, kind, rs, mat, w)
        ]
        if len(rows) < 3:
            continue
        kb, lb = best_k(rows)
        kf, lf = best_k(rows, frozen=True)
        at = [s.leverage for s in rows if s.k == 8.5]
        best.append(BestRow(case, kind, mat, rs, w, kb, lb, kf, lf, at[0] if at else math.nan))
    for b in best:
        print(
            f"  best {b.case:<28} {b.kind:<7} R {b.sheet_radius:3.1f} {b.charge:<5}"
            f" w {b.w / 1e3:4.1f}: k {b.best_k:4.1f}"
            f" -> {b.leverage:5.2f}x (frozen rebound: k {b.best_k_frozen:4.1f} ->"
            f" {b.leverage_frozen:5.2f}x) | at k 8.5 {b.leverage_at_8_5:5.2f}x"
        )
    print(
        f"  parent: k 8.5, eta 0.775 -> {leverage(1 + ETA_PARENT * math.sqrt(9.5), 8.5):.2f}x;"
        f" bare f 0.8 -> {leverage(1.6, 0.0):.2f}x"
    )
    plates, erosion = [], []
    for s in sweep:
        if s.kind == "point" or s.w != W_REF:
            continue
        rec = _record(s.case, s.kind, s.charge, s.k, s.sheet_radius)
        plates.append(plate_case(s, rec))
        erosion.append(erosion_case(s, rec))
    for p in plates:
        print(
            f"  plate {p.case:<28} {p.kind:<7} R {p.sheet_radius:3.1f} {p.charge:<5}"
            f" k {p.k:4.1f}: J {p.pulse_ns / 1e6:4.2f}-"
            f"{p.pulse_end_ns / 1e6:4.2f} MN s, ship {p.ship_t_1g:4.0f} t at 1 g | spring"
            f" {p.spring_t:4.1f} t | centre x{p.centre_over_mean:4.1f} | ringing even"
            f" {p.ringing_uniform:5.1%}, graded {p.ringing_graded:5.1%} ({p.graded_t:4.1f} t,"
            f" centre {p.graded_centre_steel_mm:5.0f} mm) | rock/m {p.rock_per_m_uniform:.2f} ->"
            f" {p.rock_per_m_graded:.2f} | peak {p.peak_gpa:5.2f} GPa"
        )
    for e in erosion:
        print(
            f"  erosion {e.case:<28} {e.kind:<7} R {e.sheet_radius:3.1f} {e.charge:<5}"
            f" k {e.k:4.1f}: {e.pulse_us:4.0f} us,"
            f" inner {e.inner_energy_mj_m2:6.2f} MJ/m^2 x wall {e.wall_share:.3f} ->"
            f" {e.face_heat_mj_m2:5.3f} vs melt {e.melt_threshold_mj_m2:5.3f} MJ/m^2 |"
            f" {e.eroded_um_per_pulse[0]:.3f}-{e.eroded_um_per_pulse[1]:.3f} um/pulse,"
            f" {e.eroded_mm_per_push[0]:.2f}-{e.eroded_mm_per_push[1]:.2f} mm/push | spray"
            f" {e.spray_kg_per_pulse:5.1f} kg/pulse out to {e.hot_radius_m:.1f} m | lev"
            f" {e.leverage:5.2f} -> {e.leverage_with_spray:5.2f} with the spray carried,"
            f" {e.leverage_with_graphite:5.2f} with graphite ({e.graphite_kg_per_pulse[0]:.1f}-"
            f"{e.graphite_kg_per_pulse[1]:.1f} kg/pulse)"
        )
    misses = miss_rows()
    for mr in misses:
        print(
            f"  miss {mr.grid:>4} {mr.case:<28} bag R {mr.bag_radius:3.1f}"
            f" bumper {mr.bumper_k:5.3f} (r {mr.bumper_radius_m:.2f}):"
            f" {mr.impulse_mn_s:4.2f} MN s, peak"
            f" {mr.peak_gpa:6.2f} GPa at r {mr.peak_radius_m:.2f} m | arriving"
            f" {mr.inner_energy_mj_m2:7.1f} MJ/m^2 inside 0.6 m,"
            f" {mr.mid_energy_mj_m2:6.1f} at 1-2 m"
            f" | {'holds' if mr.holds_maraging else 'exceeds'} maraging"
        )
    aim = aim_rows()
    for am in aim:
        print(
            f"  aim {am.case:<28} k {am.k:4.1f}: capture {am.axial:.3f} | per m miss: tilt"
            f" {am.tilt_deg_per_m:6.2f} deg, lateral {am.lateral_per_m:6.3f}, centroid"
            f" {am.centroid_per_m:6.3f} m, lever about the CM ({am.cm_height:.2f} m up)"
            f" {am.lever_per_m:6.3f} m | rock share {am.rock_share_per_cm:.4f} per cm"
        )
    conv = []
    for rc_ in all_runs(fine=False):
        if _f(rc_, "gamma") != 1.4:
            continue
        try:
            fine = _record(
                str(rc_["case"]),
                str(rc_["kind"]),
                str(rc_["charge"]),
                _f(rc_, "k"),
                _f(rc_, "sheet_radius"),
            )
        except StopIteration:
            continue
        conv.append(
            ConvergenceRow(
                str(rc_["case"]),
                str(rc_["kind"]),
                str(rc_["charge"]),
                _f(rc_, "sheet_radius"),
                _f(rc_, "k"),
                _f(fine, "impulse_per_mv"),
                _f(rc_, "impulse_per_mv"),
                _f(fine, "peak_per_mv2"),
                _f(rc_, "peak_per_mv2"),
            )
        )
    _write(OUT_ARRIVAL, list(arr))
    _write(OUT_RACE, list(races))
    _write(OUT_SWEEP, list(sweep))
    _write(OUT_BEST, list(best))
    _write(OUT_PLATE, list(plates))
    _write(OUT_AIM, list(aim))
    _write(OUT_EROSION, list(erosion))
    if conv:
        _write(OUT_CONVERGENCE, list(conv))
    if misses:
        _write(OUT_MISS, list(misses))
    print(f"python: wrote {OUT_DIR}/iapetus_*.csv")


if __name__ == "__main__":
    main()
