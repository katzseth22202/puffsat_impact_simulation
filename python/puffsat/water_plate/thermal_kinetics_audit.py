"""Reaction-control and energy-feedback audit on the saved plate-loading history."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate import kinetics_audit as base
from puffsat.water_plate import thermal_kinetics as tk
from puffsat.water_plate.flow import FlowCell, _write
from puffsat.water_plate.kinetics_audit import Row
from puffsat.water_plate.thermodynamics import OUTPUT_DIR

VARIANTS = ("baseline", "without_direct", "association_1e3_slower", "association_1e6_slower")
MODES = ("isothermal", "isoenergetic")


def scales_for(variant: str) -> kn.Vec:
    if variant not in VARIANTS:
        raise ValueError("unknown rate variant")
    scales = np.ones(10)
    for i, reaction in enumerate(kn.REACTIONS):
        if variant == "without_direct" and reaction.source_number == 15:
            scales[i] = 0
        if reaction.colliders is not None and variant.startswith("association_"):
            scales[i] = 1e-3 if variant == "association_1e3_slower" else 1e-6
    return scales


def validate_model_metadata(metadata: dict[str, object]) -> None:
    table = metadata.get("table_provenance")
    if not isinstance(table, dict) or (
        table.get("dense_hot_cv_residual_power") != 2.0 or table.get("real_fluid_max_k") != 1000.0
    ):
        raise ValueError("thermal audit requires the 1000 K, T^-2 residual parent closure")


class Diagnostic:
    def __init__(self, curve: tk.CvCurve) -> None:
        self.curve = curve

    def __call__(self, cell: FlowCell, comp: ew.Composition | None, weight: float) -> Row:
        row = base.diagnose(cell, comp, weight)
        if comp is None or "energy_relaxation_s" not in row:
            return row
        try:
            u, cv = tk.species_thermo(cell.temperature_k)
            residual = (
                cell.rho_kg_m3
                * self.curve.reference_cv(cell.rho_kg_m3)
                * (1000 / cell.temperature_k) ** 2
            )
            capacity = tk.frozen_cv(comp, cv, residual)
            if not np.isfinite(capacity) or capacity <= 0:
                raise ValueError("nonpositive or nonfinite frozen heat capacity")
            row["frozen_cv_j_kg_k"] = capacity / cell.rho_kg_m3
            row["residual_cv_fraction"] = residual / capacity
            n = kn.neutral_densities(comp)
            for variant in VARIANTS:
                try:
                    result = tk.response(cell.temperature_k, n, u, capacity, scales_for(variant))
                    row[variant + "_thermal_over_isothermal"] = (
                        result.isoenergetic_s / result.isothermal_s
                    )
                    for mode in MODES:
                        prefix = variant + "_" + mode
                        time = (
                            result.isothermal_s if mode == "isothermal" else result.isoenergetic_s
                        )
                        row[prefix + "_time_s"] = time
                        row[prefix + "_time_over_baseline"] = time / float(
                            row["energy_relaxation_s"]
                        )
                        if variant in ("baseline", "without_direct"):
                            sensitivities = (
                                result.isothermal_sensitivities
                                if mode == "isothermal"
                                else result.isoenergetic_sensitivities
                            )
                            for r, s in zip(kn.REACTIONS, sensitivities, strict=True):
                                row[prefix + f"_sensitivity_r{r.source_number}"] = s
                    row["thermal_metric_strength"] = result.thermal_metric_strength
                    row["fixed_energy_susceptibility_ratio"] = result.susceptibility_ratio
                except ValueError as error:
                    row[variant + "_error"] = str(error)
        except ValueError as error:
            row["thermal_input_error"] = str(error)
        return row


def summarize(rows: list[Row], capacity: float = 0) -> Row:
    result = base.summarize(rows, capacity)
    total = sum(float(r["weight_j"]) for r in rows)
    for variant in VARIANTS:
        for mode in MODES:
            prefix = variant + "_" + mode
            time_key = prefix + "_time_over_baseline"
            for row in rows:
                if "da_tracking" in row and time_key in row:
                    row[prefix + "_da_tracking"] = float(row["da_tracking"]) / float(row[time_key])
            resolved = sum(float(r["weight_j"]) for r in rows if time_key in r)
            result[prefix + "_resolved_rate_fraction"] = resolved / total if total else ""
            for threshold in (0, 10):
                covered = sum(
                    float(r["weight_j"])
                    for r in rows
                    if float(r.get(prefix + "_da_tracking", 0)) > threshold
                )
                result[prefix + f"_da_gt{threshold}_fraction"] = covered / total if total else ""
            for key in (prefix + "_da_tracking", time_key):
                for p in (0.1, 0.5, 0.9):
                    result[key + f"_p{int(100 * p)}"] = base.weighted_quantile(rows, key, p)
            if variant in ("baseline", "without_direct"):
                for reaction in kn.REACTIONS:
                    key = prefix + f"_sensitivity_r{reaction.source_number}"
                    numerator = sum(float(r["weight_j"]) * float(r.get(key, 0)) for r in rows)
                    # Conditional mean, with resolved coverage explicitly alongside.
                    result[key + "_mean"] = numerator / resolved if resolved else ""
        for p in (0.1, 0.5, 0.9):
            key = variant + "_thermal_over_isothermal"
            result[key + f"_p{int(100 * p)}"] = base.weighted_quantile(rows, key, p)
    for key in (
        "frozen_cv_j_kg_k",
        "residual_cv_fraction",
        "thermal_metric_strength",
        "fixed_energy_susceptibility_ratio",
    ):
        for p in (0.1, 0.5, 0.9):
            result[key + f"_p{int(100 * p)}"] = base.weighted_quantile(rows, key, p)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--strides", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--residual-knots", type=int, default=512)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "thermal_kinetics_audit.csv")
    args = parser.parse_args()
    for path in args.paths:
        with path.open() as stream:
            validate_model_metadata(json.loads(next(stream)))
    curve = tk.CvCurve.build(args.residual_knots)
    rows = [
        r
        for path in args.paths
        for stride in args.strides
        for r in base.audit(
            path, stride, point_diagnostic=Diagnostic(curve), summary_diagnostic=summarize
        )
    ]
    _write(args.output, rows)
    model = {
        "source_url": kn.SOURCE_URL,
        "source_sha256": kn.SOURCE_SHA256,
        "reactions": [asdict(r) for r in kn.REACTIONS],
        "rate_scales": {name: scales_for(name).tolist() for name in VARIANTS},
        "parents": [base.parent_provenance(p) for p in args.paths],
        "strides": args.strides,
        "window_edges_s": base.EDGES,
        "residual_log_rho": curve.log_rho.tolist(),
        "residual_cv_values": curve.values.tolist(),
        "residual_cv_slopes": curve.slopes.tolist(),
        "kinetics_validated": False,
        "limits": [
            "local linear response only; no nonlinear kinetics or finite-rate hydrodynamics",
            "isoenergetic means fixed volume and total energy, with ions chemically fixed",
            "frozen cv includes neutrals, ion/electron sensible heat and provisional residual",
            "fixed-energy susceptibility/observable transformed consistently, not iso-T weights",
            "sensitivity is -d ln(tau)/d ln(rate), not causal energy or impulse",
            "group slowdowns multiply both directions of ALL six association channels",
            "slowdown factors are parameters, not evaluated kinetic uncertainty bounds",
            "positive interval bond decreases and clocks use the unchanged equilibrium history",
            "missing/unresolved weights stay in denominators; means condition on resolved rates",
            "other source limitations from ADR-0046 still apply; no new evaluated rate set",
        ],
    }
    args.output.with_suffix(".json").write_text(json.dumps(model, indent=2) + "\n")


if __name__ == "__main__":
    main()
