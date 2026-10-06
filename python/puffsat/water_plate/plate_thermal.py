"""The spray plate's face over a push: ablative film on maraging, propellant-cooled (ADR-0055).

Drives `walled_nozzle.wall_layers` with each pulse's radiant energy on the face, taken from the
solved step 2a' columns (`spray_levers.jsonl`, 4 m cloud, 4 m pulse, bare face). The film is
modelled with the pitch properties the wall solver already carries, as a stand-in for the
carbon-loaded oil. The substrate is maraging 300 (`docs/spray_plate_layer_properties.md`).

Three cycles are compared, each over 1,500 pulses at 4 Hz:
  * none: adiabatic between pulses, so the steel ratchets;
  * respray: the film alone, resprayed cold once a cycle, as the only coolant;
  * water: a water spray from 10 to 100 ms after each pulse, at a floating (Leidenfrost) or wetting
    coefficient, then the cold respray.

The radiant load is a top-hat of duration `tau` from an opaque gas whose temperature carries the
pulse's energy. `tau` is bracketed, since the pulse length is not yet extracted from the runs.
"""

from __future__ import annotations

import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.walled_nozzle import wall_layers as wl

LEVERS = Path("data/results/water_plate/spray_levers.jsonl")
FILMS = Path("data/results/water_plate/spray_film.jsonl")
OUTPUT = Path("data/results/water_plate/plate_thermal.csv")

#: Maraging 300 [NI76]: k 21-28 W/m/K over 20-480 C (mid value), c_p 460 J/kg/K, rho 8000.
#: The melt cap is a placeholder well above anything reached; the solidus is not verified.
MARAGING = wl.Substrate("maraging 300", k=25.0, rho=8000.0, c=460.0, melt=1700.0)
#: Maraging's aging temperature [K] (480 C), the bulk limit.
AGING_LIMIT = 753.0
FOOTPRINT = math.pi * 25.0
#: Water warmed from 300 K and boiled: 0.3 + 2.26 MJ/kg.
WATER_HEAT = 2.56e6
PULSES = 1500
#: Film thickness per arm [m]: the step 1b film masses at ~1,000 kg/m^3 over the footprint.
FILM = {"argon": 150e-6, "water": 30e-6}
#: PuffSat mass per pulse [kg] at each closing speed (the 100 t plate's 8 MN s pulse).
PUFFSAT_MASS = {45580.0: 47.0, 65130.0: 33.0}


#: Carbon-loaded oil film (docs/spray_plate_film_properties.md, central values): k 0.20,
#: rho 960, c 1800, cracking cap ~1300 K (estimate, not verified), removal heat 1.6 MJ/kg.
OIL = (0.20, 960.0, 1800.0, 1300.0, 1.6e6)


@dataclass(frozen=True)
class Case:
    material: str
    w: float
    tau: float
    scheme: str
    h_cool: float
    #: "pitch" (the paper model at the step-1b thickness), or "pitch2x"/"oil2x": the film at twice
    #: its shielded consumption per pulse (spray_film.jsonl, kappa_vapor 5e3).
    film: str = "pitch"
    #: Steel behind the film [m]: the solver's 30 mm, or a thinner floor (70 t averages ~28 mm).
    steel_depth: float = wl.STEEL_DEPTH
    #: PuffSat mass per pulse [kg]; 0 takes the first study's `PUFFSAT_MASS` (8 MN s).
    puffsat_kg: float = 0.0
    #: Pitch thickness [m] for the "pitch" film; 0 takes the step-1b `FILM`.
    film_m: float = 0.0


@dataclass(frozen=True)
class Row:
    material: str
    w: float
    tau: float
    scheme: str
    h_cool: float
    film: str
    steel_depth: float
    puffsat_kg: float
    film_thickness: float
    face_energy: float
    end_steel: float
    peak_steel: float
    film_removed: float
    cooled: float
    respray: float
    water_kg_per_pulse: float
    under_aging_limit: bool


def face_energy(material: str, w: float, puffsat_kg: float = 0.0) -> float:
    """Radiant energy on the face per pulse [J/m^2], from the solved 4 m cloud, 4 m pulse row.

    The face takes the row's share of the pulse's kinetic energy over the footprint, so it scales
    with the PuffSat mass per pulse (`PUFFSAT_MASS` when `puffsat_kg` is 0)."""
    with LEVERS.open() as fh:
        for line in fh:
            r = json.loads(line)
            if (
                r["material"] == material
                and r["w"] == w
                and r["depth"] == 4.0
                and r["pulse_length"] == 4.0
            ):
                sigma_p = (puffsat_kg or PUFFSAT_MASS[w]) / FOOTPRINT
                return float(r["face_radiation_share"]) * 0.5 * sigma_p * w * w
    raise ValueError(f"no levers row for {material} at {w}")


def shielded_film_kg(material: str, w: float, film: str) -> float:
    """Film consumed per pulse [kg] on the shielded face (kappa_vapor 5e3) for `film`."""
    with FILMS.open() as fh:
        for line in fh:
            r = json.loads(line)
            run = r["run"]
            if (
                r["film"] == film
                and run["material"] == material
                and run["w"] == w
                and run["kappa_vapor"] == 5.0e3
            ):
                return float(run["film_per_pulse"])
    raise ValueError(f"no film row for {film} {material} {w}")


def coating(case: Case) -> wl.Coating:
    """The film: pitch at the step-1b thickness, or pitch or oil at twice its shielded use."""
    if case.film == "pitch":
        return wl.pitch(case.film_m or FILM[case.material])
    base = "oil" if case.film == "oil2x" else "pitch"
    k, rho, c, cap, heat = OIL if base == "oil" else wl.PITCH
    thickness = 2.0 * shielded_film_kg(case.material, case.w, base) / (rho * FOOTPRINT)
    return wl.Coating(f"{base} {thickness * 1e6:.0f} um", thickness, k, rho, c, cap, heat, False)


def pulse_gas(energy: float, tau: float) -> wl.GasHistory:
    """An opaque gas whose grey emission delivers `energy` [J/m^2] over `tau` [s]."""
    temp = (energy / (tau * wl.WALL_EMISSIVITY * wl.SIGMA_SB)) ** 0.25
    return wl.GasHistory(
        time=np.array([0.0, tau]),
        temp=np.array([temp, temp]),
        pressure=np.array([1.0e9, 1.0e9]),
        emissivity=np.array([1.0, 1.0]),
        h0=0.0,
    )


def run_case(case: Case) -> Row:
    energy = face_energy(case.material, case.w, case.puffsat_kg)
    respray = None if case.scheme == "none" else wl.Respray(temp=280.0, at=0.20)
    cooling = (
        wl.Cooling(h=case.h_cool, temp=373.0, start=0.010, end=0.100)
        if case.scheme == "water"
        else None
    )
    coat = coating(case)
    run = wl.run(
        coat,
        pulse_gas(energy, case.tau),
        pulses=PULSES,
        substrate=MARAGING,
        cooling=cooling,
        respray=respray,
        steel_depth=case.steel_depth,
    )
    return Row(
        material=case.material,
        w=case.w,
        tau=case.tau,
        scheme=case.scheme,
        h_cool=case.h_cool,
        film=case.film,
        steel_depth=case.steel_depth,
        puffsat_kg=case.puffsat_kg or PUFFSAT_MASS[case.w],
        film_thickness=coat.thickness,
        face_energy=energy,
        end_steel=run.end_steel,
        peak_steel=run.peak_steel,
        film_removed=float(run.removed_per_pulse) * coat.rho * FOOTPRINT,
        cooled=run.cooled_per_pulse,
        respray=run.respray_per_pulse,
        water_kg_per_pulse=max(run.cooled_per_pulse, 0.0) * FOOTPRINT / WATER_HEAT,
        under_aging_limit=run.peak_steel < AGING_LIMIT,
    )


def cases() -> list[Case]:
    out = []
    for material in ("argon", "water"):
        for w in PUFFSAT_MASS:
            out.append(Case(material, w, 1e-3, "none", 0.0))
            out.append(Case(material, w, 1e-3, "respray", 0.0))
            for h in (300.0, 3.0e3, 3.0e4):
                out.append(Case(material, w, 1e-3, "water", h))
    for tau in (3e-4, 3e-3):
        out.append(Case("argon", 45580.0, tau, "respray", 0.0))
    for film in ("pitch2x", "oil2x"):
        for material in ("argon", "water"):
            for w in PUFFSAT_MASS:
                out.append(Case(material, w, 1e-3, "respray", 0.0, film))
                out.append(Case(material, w, 1e-3, "water", 3.0e3, film))
    # The 70 t floor (ADR-0055, Q57): ~28 mm on average, ~15 mm at a thinned rim.
    for depth in (0.028, 0.015):
        for material in ("argon", "water"):
            for w in PUFFSAT_MASS:
                out.append(Case(material, w, 1e-3, "respray", 0.0, "pitch", depth))
    # The 150 t design (Q64): 12 MN s per pulse at k = 8.52 and eta_jet 0.6 (beta = 2.85), on the
    # d/D 0.30 bowl's ~30 mm floor, with the standing film thickened if the steel needs it.
    for w in PUFFSAT_MASS:
        kg = 12.0e6 / (2.85 * w)
        for film_m in (150e-6, 225e-6, 300e-6):
            out.append(Case("argon", w, 1e-3, "respray", 0.0, "pitch", 0.030, kg, film_m))
    return out


def main() -> None:
    with ProcessPoolExecutor() as pool:
        rows = list(pool.map(run_case, cases()))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(rows[0])))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    for r in rows:
        print(
            f"{r.material:5} w={r.w / 1e3:5.2f} {r.film:7} {r.film_thickness * 1e6:4.0f}um "
            f"steel {r.steel_depth * 1e3:2.0f}mm PuffSat {r.puffsat_kg:3.0f}kg "
            f"tau={r.tau * 1e3:3.1f}ms {r.scheme:7} "
            f"h={r.h_cool:7.0f}: face {r.face_energy / 1e6:5.2f} MJ/m2, steel end "
            f"{r.end_steel:5.0f} K peak {r.peak_steel:5.0f} K, film {r.film_removed:5.1f} kg, "
            f"water {r.water_kg_per_pulse:5.1f} kg, under 480C={r.under_aging_limit}"
        )


if __name__ == "__main__":
    main()
