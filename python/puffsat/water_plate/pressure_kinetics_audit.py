"""Water-association pressure-source coverage and channel-removal sensitivity.

Uses the same saved equilibrium histories and weights as kinetics_audit. Other
nine rates, equilibrium composition, clocks and energy susceptibility stay fixed.
No new trajectory, heat release or impulse prediction is made.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate import kinetics_audit as audit_base
from puffsat.water_plate import pressure_kinetics as pk
from puffsat.water_plate.flow import FlowCell, _write
from puffsat.water_plate.kinetics_audit import Row
from puffsat.water_plate.thermodynamics import OUTPUT_DIR

VARIANTS = ("without_direct_water", "source_plog")


def diagnose(cell: FlowCell, comp: ew.Composition | None, weight: float) -> Row:
    row = audit_base.diagnose(cell, comp, weight)
    if comp is None:
        return row
    n = kn.neutral_densities(comp)
    peff = pk.effective_pressure(cell.temperature_k, n)
    row["source_effective_pressure_pa"] = peff
    row["source_water_efficiency"] = pk.water_efficiency(cell.temperature_k)
    row["source_default_collider_fraction"] = float(np.sum(n[:-1]) / np.sum(n))
    row["source_pressure_tabulated"] = (
        pk.PLOG_NODES[0][0] <= peff / pk.ATM_PA <= pk.PLOG_NODES[-1][0]
    )
    if "energy_relaxation_s" not in row:
        return row  # Preserve baseline tiny/unresolved states in every denominator.
    for variant in VARIANTS:
        try:
            scale = 0.0
            if variant == "source_plog":
                # A documented forward ASSOCIATION rate, not a NASA-derived
                # reverse dissociation rate. Match the reverse to our own K.
                source_k = pk.water_coefficient(cell.temperature_k, peff)
                direct = next(r for r in kn.REACTIONS if r.source_number == 15)
                assert direct.colliders is not None
                baseline_k = direct.coefficient(cell.temperature_k) * float(
                    np.dot(direct.colliders, n)
                )
                scale = source_k / baseline_k
                row["source_association_over_baseline"] = scale
            result = kn.relaxation(cell.temperature_k, n, water_association_scale=scale)
            row[variant + "_time_s"] = result.energy_time_s
            row[variant + "_time_over_baseline"] = result.energy_time_s / float(
                row["energy_relaxation_s"]
            )
            row[variant + "_status"] = "resolved"
        except ValueError as error:
            row[variant + "_status"] = str(error)
    return row


def summarize(rows: list[Row], reference_capacity_j: float = 0) -> Row:
    result = audit_base.summarize(rows, reference_capacity_j)
    total = sum(float(r["weight_j"]) for r in rows)
    # The base audit adds its clocks AFTER calling the point diagnostic.
    # Reuse exactly those clocks, never invent one for missing/invalid states.
    for row in rows:
        for variant in VARIANTS:
            ratio_key = variant + "_time_over_baseline"
            if "da_tracking" in row and ratio_key in row:
                row[variant + "_da_tracking"] = float(row["da_tracking"]) / float(row[ratio_key])
    for flag in ("source_pressure_tabulated", "source_above_collider_data_2000k"):
        covered = sum(
            float(r["weight_j"])
            for r in rows
            if (
                bool(r.get(flag, False))
                if flag == "source_pressure_tabulated"
                else float(r["temperature_k"]) > 2000
            )
        )
        result[flag + "_fraction"] = covered / total if total > 0 else ""
    keys = [
        "source_effective_pressure_pa",
        "source_water_efficiency",
        "source_default_collider_fraction",
        "source_association_over_baseline",
    ]
    for variant in VARIANTS:
        keys.extend([variant + "_time_over_baseline", variant + "_da_tracking"])
        for suffix, threshold in (("resolved_clock", 0), ("da_gt10", 10), ("da_gt10000", 10000)):
            value = sum(
                float(r["weight_j"])
                for r in rows
                if float(r.get(variant + "_da_tracking", 0)) > threshold
            )
            result[variant + "_" + suffix + "_fraction"] = value / total if total > 0 else ""
        unresolved = sum(
            float(r["weight_j"])
            for r in rows
            if not math.isfinite(float(r.get(variant + "_time_s", math.nan)))
        )
        result[variant + "_no_rate_fraction"] = unresolved / total if total > 0 else ""
    for key in keys:
        for fraction in (0.1, 0.5, 0.9):
            result[f"{key}_p{int(100 * fraction)}"] = audit_base.weighted_quantile(
                rows, key, fraction
            )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--strides", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "pressure_kinetics_audit.csv")
    args = parser.parse_args()
    rows = [
        r
        for path in args.paths
        for stride in args.strides
        for r in audit_base.audit(
            path, stride, point_diagnostic=diagnose, summary_diagnostic=summarize
        )
    ]
    _write(args.output, rows)
    model = {
        "pressure_source_url": pk.SOURCE_URL,
        "pressure_source_sha256": pk.SOURCE_SHA256,
        "pressure_nodes_atm_a_m3_kmol_s_b_ea_cal_mol": pk.PLOG_NODES,
        "water_efficiency_a_b_ea_cal_mol": pk.WATER_EFFICIENCY,
        "baseline_source_url": kn.SOURCE_URL,
        "baseline_source_sha256": kn.SOURCE_SHA256,
        "baseline_reactions": [asdict(r) for r in kn.REACTIONS],
        "parents": [audit_base.parent_provenance(p) for p in args.paths],
        "strides": args.strides,
        "window_edges_s": audit_base.EDGES,
        "kinetics_validated": False,
        "limits": [
            "saved equilibrium history, isothermal fixed-density fixed-ion local relaxation",
            "positive interval bond-store decreases are neither unique energy nor impulse weights",
            "only direct H+OH association/dissociation changes; other nine rates remain baseline",
            "source_plog uses tabulated effective pressure only; T applicability unestablished",
            "finite PLOG upper node is not a high-pressure limit; no pressure endpoint clamping",
            "neutral ideal pressure times mixture efficiency, not dense EOS pressure",
            "represented neutrals except H2O inherit the default N2 collider efficiency",
            "ions and omitted species do not contribute to collider pressure",
            "reverse rate matched to study equilibrium constants, not imported NASA7",
            "without_direct_water removes both directions of reaction 15 only",
            "removal is a slow endpoint only for this linear network with its other rates fixed",
            "collider data at 300, 1000, 2000 K are not a full reaction validity interval",
            "quantiles use only resolved subsets; coverage fractions retain all positive weights",
            "no kinetic uncertainty envelope or new finite-rate impulse prediction",
        ],
    }
    args.output.with_suffix(".json").write_text(json.dumps(model, indent=2) + "\n")


if __name__ == "__main__":
    main()
