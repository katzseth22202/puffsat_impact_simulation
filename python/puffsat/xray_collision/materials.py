"""Atomic data for the collision materials: NIST ionization ladders and ground-term weights.

Ionization energies [eV] for stage z -> z+1 and ground-level statistical weights `2J+1` for
stages 0..Z (the last entry is the bare nucleus), from the NIST Atomic Spectra Database
ionization-energy query (physics.nist.gov/cgi-bin/ASD/ie.pl), retrieved 2026-10-08. The argon
ladder matches `puffsat.water_plate.spray_eos.IP_AR_EV` (retrieved 2026-10-05).

Opacities are TOPS (LANL OPLIB) gray Rosseland/Planck means for the pure element, pulled over
rho 1e-9..20-30 g/cc and T 0.5-60 eV (`make xray-collision-tops`).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Material:
    """One element, as a condensed ball before the collision."""

    name: str
    mass_amu: float
    rho0: float  # condensed starting density [kg/m^3]
    ip_ev: tuple[float, ...]  # stage z -> z+1, z = 0..Z-1
    g_ground: tuple[float, ...]  # stages 0..Z
    tops_pull: Path


ARGON = Material(
    name="argon",
    mass_amu=39.948,
    rho0=1395.0,  # liquid at its 87 K normal boiling point
    ip_ev=(
        15.7596119,
        27.62967,
        40.735,
        59.58,
        74.84,
        91.29,
        124.41,
        143.4567,
        422.6,
        479.76,
        540.4,
        619.0,
        685.5,
        755.13,
        855.5,
        918.375,
        4120.66559,
        4426.22407,
    ),
    g_ground=(1, 4, 5, 4, 1, 2, 1, 2, 1, 4, 5, 4, 1, 2, 1, 2, 1, 2, 1),
    tops_pull=Path("data/tables/tops/tops_argon_fireball_gray.html"),
)

IRON = Material(
    name="iron",
    mass_amu=55.845,
    rho0=7874.0,  # solid at room temperature
    ip_ev=(
        7.9024681,
        16.19921,
        30.651,
        54.91,
        75.0,
        98.985,
        124.9671,
        151.06,
        233.6,
        262.1,
        290.9,
        330.8,
        361.0,
        392.2,
        456.2,
        489.312,
        1262.7,
        1357.8,
        1460.0,
        1575.6,
        1687.0,
        1798.4,
        1950.4,
        2045.759,
        8828.1864,
        9277.6886,
    ),
    g_ground=(9, 10, 9, 6, 1, 4, 5, 4, 1, 4, 5, 4, 1, 2, 1, 2, 1, 4, 5, 4, 1, 2, 1, 2, 1, 2, 1),
    tops_pull=Path("data/tables/tops/tops_iron_fireball_gray.html"),
)

MATERIALS = {m.name: m for m in (ARGON, IRON)}
