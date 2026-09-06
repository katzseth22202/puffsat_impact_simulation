"""Local equilibrium merge screen with a liquid/gas energy reference bridge.

This closed-parcel calculation does not predict the spatial mixing ratio, density,
peak wall state, cooling history, or impulse. It asks what a parcel would become
if the specified masses mixed inelastically at a specified final gas density.
Liquid energy is aligned to the existing gas EOS at a dilute molecular vapor
state using IAPWS-95/CoolProp. Formation of gas then draws from the conserved
collision energy automatically; dissociation/ionization already live in the EOS.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import CoolProp.CoolProp as CP

from puffsat import eos_water, expansion
from puffsat.water_plate.ledger import ANCHORS, PROJECTILE_KG, WATER_RATIOS, _check

REFERENCE_RHO = 0.01
PROJECTILE_RHO = 0.1241518740905039  # D=10 m, r_foot/R=.75, 100 us, 45.58 km/s input.
STORED_WATER_T = 293.15
STORED_WATER_P = 101_325.0
OUTPUT_DIR = Path("data/results/water_plate")


@dataclass(frozen=True)
class EnergyReference:
    """Add offset_j_kg to CoolProp internal energy to match the bound-gas EOS zero."""

    match_temperature_k: float
    match_density_kg_m3: float
    offset_j_kg: float
    liquid_energy_j_kg: float
    coolprop_version: str


@cache
def energy_reference(match_temperature_k: float = 400.0) -> EnergyReference:
    """Align at dilute superheated vapor; 400/600 K are checked for model mismatch.

    This transfers an internal-energy DIFFERENCE, not a handbook vaporization
    enthalpy. Pressure work for a later injector boundary must be handled there.
    It does not create a globally stitched liquid/plasma EOS.
    """
    if not 350 <= match_temperature_k <= 1000:
        raise ValueError("matching temperature must be in the molecular vapor range 350-1000 K")
    vapor_real = float(CP.PropsSI("U", "D", REFERENCE_RHO, "T", match_temperature_k, "Water"))
    vapor_model = eos_water.pressure_energy(REFERENCE_RHO, match_temperature_k)[1]
    offset = vapor_model - vapor_real
    liquid_real = float(CP.PropsSI("U", "P", STORED_WATER_P, "T", STORED_WATER_T, "Water"))
    return EnergyReference(
        match_temperature_k,
        REFERENCE_RHO,
        offset,
        liquid_real + offset,
        str(CP.get_global_param_string("version")),
    )


@dataclass(frozen=True)
class MixedState:
    """Local equilibrium state, with the input energy and separate chemical stores."""

    case: str
    closing_speed_m_s: float
    k: float
    projectile_kg: float
    projectile_temperature_k: float
    projectile_density_kg_m3: float
    injection_speed_m_s: float
    stored_water_temperature_k: float
    stored_water_pressure_pa: float
    liquid_energy_j_kg: float
    energy_offset_j_kg: float
    reference_match_temperature_k: float
    reference_match_density_kg_m3: float
    coolprop_version: str
    mixed_density_kg_m3: float
    merge_velocity_m_s: float
    merge_heat_j: float
    total_input_energy_j: float
    mixed_internal_energy_j_kg: float
    temperature_k: float
    pressure_pa: float
    bond_energy_j_kg: float
    bond_fraction: float
    ionization_energy_j_kg: float
    thermal_energy_j_kg: float


def mixed_state(
    speed_m_s: float,
    k: float,
    mixed_density_kg_m3: float,
    *,
    case: str = "3-synodic",
    projectile_temperature_k: float = 400.0,
    injection_speed_m_s: float = 0.0,
    match_temperature_k: float = 400.0,
) -> MixedState:
    """Conserve parcel momentum and total energy, then invert the equilibrium EOS.

    Projectile velocity is +w, injected water -u. Relative merge heating is
    (m_i/2)*k/(1+k)*(w+u)^2, distinct from the bulk KE remaining afterward.
    The negative liquid energy on the gas zero is essential: replacing it with
    zero or with cold vapor would grant an uncharged phase-change energy input.
    """
    _check("speed_m_s", speed_m_s, positive=True)
    _check("k", k)
    _check("mixed_density_kg_m3", mixed_density_kg_m3, positive=True)
    _check("projectile_temperature_k", projectile_temperature_k, positive=True)
    _check("injection_speed_m_s", injection_speed_m_s)
    reference = energy_reference(match_temperature_k)
    projectile_e = eos_water.pressure_energy(PROJECTILE_RHO, projectile_temperature_k)[1]
    mass = PROJECTILE_KG * (1 + k)
    velocity = (speed_m_s - k * injection_speed_m_s) / (1 + k)
    heat = 0.5 * PROJECTILE_KG * k / (1 + k) * (speed_m_s + injection_speed_m_s) ** 2
    initial_internal = PROJECTILE_KG * (projectile_e + k * reference.liquid_energy_j_kg)
    energy = (initial_internal + heat) / mass
    total_input = initial_internal + 0.5 * PROJECTILE_KG * (
        speed_m_s**2 + k * injection_speed_m_s**2
    )
    temp = expansion.temperature_at(mixed_density_kg_m3, energy, eos_water.pressure_energy)
    pressure = eos_water.pressure_energy(mixed_density_kg_m3, temp)[0]
    comp = eos_water.composition(mixed_density_kg_m3, temp)
    bonds = eos_water.bond_energy_held(mixed_density_kg_m3, temp)
    ionization = (
        comp.n_hp * eos_water.IP_H
        + sum(n * e for n, e in zip(comp.n_o_ions, eos_water.E_O_CUM, strict=True))
    ) / mixed_density_kg_m3
    return MixedState(
        case,
        speed_m_s,
        k,
        PROJECTILE_KG,
        projectile_temperature_k,
        PROJECTILE_RHO,
        injection_speed_m_s,
        STORED_WATER_T,
        STORED_WATER_P,
        reference.liquid_energy_j_kg,
        reference.offset_j_kg,
        reference.match_temperature_k,
        reference.match_density_kg_m3,
        reference.coolprop_version,
        mixed_density_kg_m3,
        velocity,
        heat,
        total_input,
        energy,
        temp,
        pressure,
        bonds,
        bonds / eos_water.FULL_ATOMIZATION_ENERGY,
        ionization,
        energy - bonds - ionization,
    )


def _write(path: Path, rows: list[MixedState]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    """Write equilibrium screening states and explicit reference/input sensitivity."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = [
        mixed_state(speed, k, rho, case=case)
        for case, speed in ANCHORS
        for k in WATER_RATIOS
        if k > 0
        for rho in (0.1, 1.0, 10.0)
    ]
    _write(OUTPUT_DIR / "thermo_mixing.csv", rows)
    sensitivity = [
        mixed_state(
            45_580.0,
            8.5,
            rho,
            projectile_temperature_k=temp,
            injection_speed_m_s=speed,
            match_temperature_k=match_temp,
        )
        for rho in (0.1, 1.0, 10.0)
        for temp in (400.0, 800.0)
        for speed in (0.0, 100.0)
        for match_temp in (400.0, 600.0)
    ]
    _write(OUTPUT_DIR / "thermo_sensitivity.csv", sensitivity)
    print(f"Wrote {len(rows)} local equilibrium states and {len(sensitivity)} sensitivity rows.")
    print("No cooling history, chemical recovery, or wall impulse is inferred.")


if __name__ == "__main__":
    main()
