"""Prepare nonlinear neutral parcel histories and independent kernel fixtures.

Cold path only: thermochemistry, EOS/source preparation and result analysis.
Time integration belongs to the Rust water_plate parcels executable.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate import freeze
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate import kinetics_audit as ka
from puffsat.water_plate import thermal_kinetics as tk
from puffsat.water_plate.flow import read_output
from puffsat.water_plate.thermal_kinetics_audit import validate_model_metadata

START = 153.5e-6
ENDS = (193e-6, 250e-6)
VARIANTS = (
    "baseline",
    "li_h2",
    "without_direct_and_steam_h2",
    "association_1e3_slower",
    "association_1e6_slower",
    "frozen",
)


def chemistry(variant: str = "baseline") -> dict[str, object]:
    if variant not in VARIANTS:
        raise ValueError("unknown parcel rate variant")
    reactions = list(kn.REACTIONS)
    if variant == "li_h2":
        reactions = [r for r in reactions if r.source_number not in (12, 13, 14)]
        reactions.append(
            kn.Reaction(
                0,
                (0, 0, 1, 0, 0, 0),
                (2, 0, 0, 0, 0, 0),
                4.577e19,
                -1.4,
                104380,
                (1, 1, 2.5, 1, 1, 12),
            )
        )
    rates = []
    for r in reactions:
        factor = 1.0
        if variant == "frozen" or (
            variant == "without_direct_and_steam_h2" and r.source_number in (14, 15)
        ):
            factor = 0.0
        if variant.startswith("association_") and r.colliders is not None:
            factor = 1e-3 if variant == "association_1e3_slower" else 1e-6
        order = sum(r.reactants) + (r.colliders is not None)
        rates.append(
            {
                "reactants": r.reactants,
                "products": r.products,
                "a": factor * r.a_cgs_mol * (1e-6 / ew.N_A) ** (order - 1),
                "exponent": r.exponent,
                "activation_k": r.activation_cal_mol * 4.184 / (ew.N_A * ew.K_B),
                "colliders": r.colliders,
            }
        )
    return {
        "kb": ew.K_B,
        "molecule_mass": ew.M_H2O,
        "reference_temp": 6000.0,
        "reference_log_n": kn.equilibrium_log_reference(6000).tolist(),
        "ground_k": (kn.GROUND_ENERGIES / ew.K_B).tolist(),
        "linear_cv": [1.5, 1.5, 2.5, 2.5, 2.5, 3.0],
        "vibrations": [
            [],
            [],
            [ew.H2.theta_vib],
            [ew.OH.theta_vib],
            [ew.O2.theta_vib],
            list(ew._THETA_VIB),
        ],
        "reactions": rates,
    }


def fixture() -> dict[str, object]:
    """Independent thermochemistry/source and fixed-energy relaxation references."""
    checks = []
    for temp in (1000.0, 2000.0, 5100.0, 6760.0, 9000.0, 100000.0):
        rho = 71.0
        # Positive off-equilibrium populations, avoiding cancellation-only tests.
        y = np.array([0.4, 0.2, 0.3, 0.1, 0.15, 0.25])
        u, cv = tk.species_thermo(temp)
        n = y * rho / ew.M_H2O
        checks.append(
            {
                "temp": temp,
                "rho": rho,
                "y": y.tolist(),
                "u": (u / ew.K_B).tolist(),
                "cv": (cv / ew.K_B).tolist(),
                "log_reference": kn.equilibrium_log_reference(temp).tolist(),
                "source": (kn.net_production(temp, n) * ew.M_H2O / rho).tolist(),
            }
        )
    temp, rho = 6760.0, 71.0
    comp = ew.composition(rho, temp)
    y_eq = kn.neutral_densities(comp) * ew.M_H2O / rho
    spectator = (comp.n_hp + comp.n_e + sum(comp.n_o_ions)) * ew.M_H2O / rho
    y0 = y_eq.copy()
    # Dissociate 10% of the equilibrium H2O, conserving elements. The initial
    # temperature pays this bond cost at fixed total energy.
    extent = 0.1 * y0[5]
    y0 += extent * np.array([2, 1, 0, 0, 0, -1])

    def energy(y: kn.Vec, t: float) -> float:
        u, _ = tk.species_thermo(t)
        return (float(y @ u) + 1.5 * ew.K_B * spectator * t) / ew.M_H2O

    target = energy(y_eq, temp)
    lo, hi = 1000.0, temp
    if not energy(y0, lo) < target < energy(y0, hi):
        raise ValueError("fixed-energy perturbation is not bracketed")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if energy(y0, mid) < target:
            lo = mid
        else:
            hi = mid
    t0 = 0.5 * (lo + hi)
    tau = kn.relaxation(temp, kn.neutral_densities(comp)).energy_time_s
    return {
        "chemistry": chemistry(),
        "checks": checks,
        "relaxation": {
            "rho": rho,
            "y0": y0.tolist(),
            "temperature0": t0,
            "spectators": spectator,
            "equilibrium_y": y_eq.tolist(),
            "equilibrium_temperature": temp,
            "duration": 1000 * tau,
        },
    }


def prepare(
    parent: Path,
    stride: int,
    tolerance: float,
    maximum_step: float,
    variants: list[str],
    limit: int | None = None,
) -> dict[str, object]:
    if not all(math.isfinite(v) and v > 0 for v in (tolerance, maximum_step)):
        raise ValueError("need finite positive solver settings")
    if not variants or len(set(variants)) != len(variants):
        raise ValueError("need distinct nonempty rate variants")
    provenance = ka.parent_provenance(parent)
    with parent.open() as stream:
        validate_model_metadata(json.loads(next(stream)))
    snapshots, summaries = read_output(parent)
    combined = [s for s in snapshots if s.name.startswith("combined")]
    selected = ka.select_snapshots(combined, stride)
    ka.validate_parent(
        next(r for r in summaries if str(r["name"]).startswith("combined")), selected
    )
    selected = [s for s in selected if s.time_s >= START - 1e-12]
    if not math.isclose(selected[0].time_s, START, abs_tol=1e-12):
        raise ValueError("missing exact common parcel start")
    parcels = []
    for i, cell in enumerate(selected[0].cells):
        if limit is not None and i >= limit:
            break
        rho, temp = cell.rho_kg_m3, cell.temperature_k
        if temp <= 1000:
            raise ValueError("selected initial parcel below neutral thermo domain")
        comp = ew.composition(rho, temp)
        y = kn.neutral_densities(comp) * ew.M_H2O / rho
        spectator = (comp.n_hp + comp.n_e + sum(comp.n_o_ions)) * ew.M_H2O / rho
        initial_bond, _ = ka.cell_state(cell)
        forcing = 0.0
        history = []
        previous = cell
        parent_returns = []
        for snap in selected:
            current = snap.cells[i]
            forcing += (
                current.energy_j_kg
                - previous.energy_j_kg
                + (
                    0.5
                    * (current.pressure_pa + previous.pressure_pa)
                    * (1 / current.rho_kg_m3 - 1 / previous.rho_kg_m3)
                )
            )
            history.append(
                {"time_s": snap.time_s, "rho": current.rho_kg_m3, "forcing_j_kg": forcing}
            )
            if any(math.isclose(snap.time_s, t, abs_tol=1e-12) for t in ENDS):
                parent_returns.append(
                    {
                        "time_s": snap.time_s,
                        "net_bond_return_j_kg": initial_bond - ka.cell_state(current)[0],
                    }
                )
            previous = current
        parcels.append(
            {
                "index": i,
                "mass_kg": cell.mass_kg,
                "initial_y": y.tolist(),
                "initial_temperature": temp,
                "initial_energy": cell.energy_j_kg,
                "initial_pressure": cell.pressure_pa,
                "initial_bond": initial_bond,
                "spectators": spectator,
                "history": history,
                "parent_returns": parent_returns,
                "initial_ionized_nuclei_fraction": (comp.n_hp + sum(comp.n_o_ions))
                * ew.M_H2O
                / (3 * rho),
            }
        )
    return {
        "parent": provenance,
        "stride": stride,
        "kinetics_validated": False,
        "source_url": kn.SOURCE_URL,
        "source_sha256": kn.SOURCE_SHA256,
        "li_source_url": "https://raw.githubusercontent.com/Pele-Suite/PelePhysics/"
        "ab2c54a4183cc97ff6a390b1ddc0da5367535e1f/Mechanisms/LiDryer/mechanism.inp",
        "li_source_sha256": "a883ccfd1e4462bca3e5fb0fe29c126664cbd848e3e0a58859f92d3176f5c485",
        "residual_curve": freeze.residual_curve(),
        "parcels": parcels,
        "settings": {
            "relative_tolerance": tolerance,
            "absolute_population_tolerance": tolerance * 1e-6,
            "maximum_step_s": maximum_step,
            "parent_forcing": True,
        },
        "variants": [{"name": name, "chemistry": chemistry(name)} for name in variants],
        "forcing_modes": [True, False],
        "report_times": list(ENDS),
        "limits": [
            "conditional scenarios, not probabilities or validated rate uncertainty",
            "prescribed parent density; no chemistry feedback on flow or impulse",
            "fixed ions/electrons; omitted ion kinetics and HO2/H2O2",
            "parent residual forcing includes sampling, EOS and discretization defects",
            "no radiation/wall losses; provisional dense residual",
            "net return since 153.5 us, not total release since impact or useful work",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", action="store_true")
    parser.add_argument("--parent", type=Path)
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--tolerance", type=float, default=1e-4)
    parser.add_argument("--maximum-step", type=float, default=2.5e-7)
    parser.add_argument("--variants", nargs="+", choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.fixture:
        payload = fixture()
    else:
        if args.parent is None:
            parser.error("--parent is required unless --fixture is selected")
        if (
            not math.isfinite(args.tolerance)
            or args.tolerance <= 0
            or not math.isfinite(args.maximum_step)
            or args.maximum_step <= 0
            or (args.limit is not None and args.limit <= 0)
        ):
            parser.error("tolerances, maximum step and optional limit must be positive")
        payload = prepare(
            args.parent, args.stride, args.tolerance, args.maximum_step, args.variants, args.limit
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
