"""Step 0: analytic scales, not hydrodynamic performance predictions.

Incoming momentum is a CREDIT in this overtake study. The ceiling assumes cold
stationary carried water, negligible vehicle recoil and injection energy, complete
ejection, and all incoming kinetic energy converted to single-speed axial exhaust.
Gasification and bond energy are separately reported stores/cost scales, not
automatically permanent losses subtracted from that absolute ceiling.

No wall survival, actual ablation, chemistry recovery, or optimum is computed.
See docs/water_injected_overtake_plate_study.md for scope and provenance.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from puffsat.eos_water import FULL_ATOMIZATION_ENERGY

G0 = 9.80665
PROJECTILE_KG = 25.0
# Rounded engineering scale: liquid water near 293 K heated to ordinary boiling
# and vaporized (~0.334 + 2.257 MJ/kg). This is NOT a vacuum saturation model or
# an EOS zero-point offset. Actual source states need an EOS-consistent enthalpy.
GASIFICATION_J_KG = 2.6e6
REFERENCE_K = 8.5
MAX_DIAMETER_M = 30.0
MAX_ASSEMBLY_KG = 50_000.0  # Working interpretation of "50 tons"; no structure solve yet.
ANCHORS = (
    ("3-synodic", 45_580.0),
    ("3-synodic", 56_530.0),
    ("2-synodic", 61_830.0),
    ("2-synodic", 65_130.0),
)
WATER_RATIOS = (0.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 10.0)
COATING_THICKNESSES_M = (0.0, 1e-6, 3e-6, 10e-6, 30e-6, 100e-6)
OUTPUT_DIR = Path("data/results/water_plate")


def _check(name: str, value: float, *, positive: bool = False) -> None:
    if not math.isfinite(value) or value < 0.0 or (positive and value == 0.0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")


@dataclass(frozen=True)
class PulseLedger:
    """SI ledger; water-only ceiling Isp excludes any ablator or other expendables."""

    case: str
    closing_speed_m_s: float
    k: float
    projectile_kg: float
    water_kg: float
    incoming_momentum_ns: float
    incoming_ke_j: float
    merge_speed_m_s: float
    merge_bulk_j: float
    merge_heat_j: float
    gasification_j_kg: float
    gasification_j: float
    gasification_fraction_ke: float
    bond_energy_j_kg: float
    full_dissociation_store_j: float
    full_dissociation_fraction_ke: float
    ceiling_beta: float
    ceiling_impulse_ns: float
    water_only_ceiling_isp_s: float | None


def pulse_ledger(case: str, speed: float, k: float) -> PulseLedger:
    """Complete inelastic merge identity plus a separate axial-exhaust ceiling.

    Merge heat is the free merge BEFORE any further stagnation at the plate:
    E_heat = E_in*k/(1+k). It does not limit the later plate's total energy access.
    The dissociation scale covers all projectile + injected water, assuming
    initially molecular water. No dissociation fraction or ionization is solved.
    """
    _check("speed", speed, positive=True)
    _check("k", k)
    water = PROJECTILE_KG * k
    energy = 0.5 * PROJECTILE_KG * speed**2
    momentum = PROJECTILE_KG * speed
    gasification = water * GASIFICATION_J_KG
    bonds = (PROJECTILE_KG + water) * FULL_ATOMIZATION_ENERGY
    beta = 1.0 + math.sqrt(1.0 + k)
    impulse = beta * momentum
    return PulseLedger(
        case=case,
        closing_speed_m_s=speed,
        k=k,
        projectile_kg=PROJECTILE_KG,
        water_kg=water,
        incoming_momentum_ns=momentum,
        incoming_ke_j=energy,
        merge_speed_m_s=speed / (1.0 + k),
        merge_bulk_j=energy / (1.0 + k),
        merge_heat_j=energy * k / (1.0 + k),
        gasification_j_kg=GASIFICATION_J_KG,
        gasification_j=gasification,
        gasification_fraction_ke=gasification / energy,
        bond_energy_j_kg=FULL_ATOMIZATION_ENERGY,
        full_dissociation_store_j=bonds,
        full_dissociation_fraction_ke=bonds / energy,
        ceiling_beta=beta,
        ceiling_impulse_ns=impulse,
        water_only_ceiling_isp_s=impulse / (G0 * water) if water > 0.0 else None,
    )


@dataclass(frozen=True)
class CoatingLedger:
    """Full-face deposited inventory; full renewal is a scenario, not an ablation result."""

    diameter_m: float
    depth_over_diameter: float
    projected_area_m2: float
    surface_area_m2: float
    thickness_m: float
    layer_density_kg_m3: float
    layer_inventory_kg: float
    reference_water_kg: float
    full_renewal_fraction_reference_water: float


def coating_ledger(
    diameter: float,
    depth_over_diameter: float,
    *,
    thickness_m: float = 100e-6,
    density_kg_m3: float = 1000.0,
) -> CoatingLedger:
    """Inventory on z(r)=delta*(r/R)^2, matching the existing parabolic plate family.

    Integrating 2*pi*r*sqrt(1+(dz/dr)^2) gives A/(pi*R^2) =
    2*((1+s^2)^(3/2)-1)/(3*s^2), s=2*delta/R. The factored expression
    below is continuous and numerically stable at zero depth. Thickness is normal
    to the surface. Density is illustrative; no particular ablator is selected.
    """
    _check("diameter", diameter, positive=True)
    _check("depth_over_diameter", depth_over_diameter)
    _check("thickness_m", thickness_m)
    _check("density_kg_m3", density_kg_m3, positive=True)
    projected_area = math.pi * (diameter / 2.0) ** 2
    q = math.sqrt(1.0 + (4.0 * depth_over_diameter) ** 2)
    surface_area = projected_area * (2.0 / 3.0) * (q * q + q + 1.0) / (q + 1.0)
    inventory = surface_area * thickness_m * density_kg_m3
    reference_water = PROJECTILE_KG * REFERENCE_K
    return CoatingLedger(
        diameter,
        depth_over_diameter,
        projected_area,
        surface_area,
        thickness_m,
        density_kg_m3,
        inventory,
        reference_water,
        inventory / reference_water,
    )


def _write_rows(path: Path, rows: Iterable[PulseLedger | CoatingLedger]) -> None:
    records = [asdict(row) for row in rows]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    """Regenerate the small cited reference tables. No simulation is launched."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _write_rows(
        OUTPUT_DIR / "ledger_pulses.csv",
        (pulse_ledger(case, speed, k) for case, speed in ANCHORS for k in WATER_RATIOS),
    )
    _write_rows(
        OUTPUT_DIR / "ledger_coatings.csv",
        (
            coating_ledger(diameter, depth, thickness_m=thickness)
            for diameter in (10.0, 20.0, MAX_DIAMETER_M)
            for depth in (0.0, 0.10, 0.15)
            for thickness in COATING_THICKNESSES_M
        ),
    )
    print(f"Analytic reference tables written to {OUTPUT_DIR}; no performance prediction.")


if __name__ == "__main__":
    main()
