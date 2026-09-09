"""Release-weighted rate/applicability screen on saved planar equilibrium flow.

Weights are positive consecutive parcel bond-store decreases, not chemical
work or unique released energy. Rates/clocks use the right endpoint of each
interval; output thinning checks that sampling approximation. No rate is
declared validated by this screen, including states below NASA's 3500 K limit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections.abc import Callable
from dataclasses import asdict
from itertools import pairwise
from pathlib import Path

from puffsat import eos_water as ew
from puffsat import recombination
from puffsat.water_plate.flow import CsvValue, FlowCell, Snapshot, _write, read_output
from puffsat.water_plate.kinetics import (
    REACTIONS,
    REFERENCE_THERMO_MAX_K,
    SOURCE_SHA256,
    SOURCE_URL,
    neutral_densities,
    relaxation,
)
from puffsat.water_plate.thermodynamics import OUTPUT_DIR

Row = dict[str, CsvValue]
EDGES = (100e-6, 153.5e-6, 193e-6, 250e-6)


def store_change_clock(old: float, new: float, dt: float) -> float:
    """Mean held store / sampled return rate; distinct from an expansion clock."""
    if not (
        math.isfinite(old)
        and math.isfinite(new)
        and math.isfinite(dt)
        and old > new >= 0
        and dt > 0
    ):
        raise ValueError("need finite decreasing nonnegative stores and positive elapsed time")
    return 0.5 * (old + new) * dt / (old - new)


def weighted_quantile(rows: list[Row], key: str, fraction: float) -> CsvValue:
    samples = sorted(
        (float(r[key]), float(r["weight_j"]))
        for r in rows
        if key in r and math.isfinite(float(r[key]))
    )
    total = sum(w for _, w in samples)
    if total <= 0:
        return ""
    cumulative = 0.0
    for value, weight in samples:
        cumulative += weight
        if cumulative >= fraction * total:
            return value
    return samples[-1][0]


def select_snapshots(snapshots: list[Snapshot], stride: int) -> list[Snapshot]:
    if stride < 1:
        raise ValueError("stride must be positive")
    inside = [s for s in snapshots if EDGES[0] - 1e-12 <= s.time_s <= EDGES[-1] + 1e-12]
    for edge in EDGES:
        if not any(math.isclose(s.time_s, edge, abs_tol=1e-12, rel_tol=0) for s in inside):
            raise ValueError("need exact checkpoints at all rate-audit window boundaries")
    return [
        s
        for i, s in enumerate(inside)
        if i % stride == 0
        or any(math.isclose(s.time_s, t, abs_tol=1e-12, rel_tol=0) for t in EDGES)
    ]


def validate_parent(summary: Row, snapshots: list[Snapshot]) -> None:
    """Only accepted fixed-mass flow histories can support this parcel audit."""
    keys = (
        "input_energy_j",
        "max_energy_error_j",
        "momentum_error_ns",
        "wall_impulse_ns",
        "time_s",
    )
    if summary["status"] != "completed_time_window" or not all(
        math.isfinite(float(summary[k])) for k in keys
    ):
        raise ValueError("requires completed parent and finite conservation diagnostics")
    if float(summary["input_energy_j"]) <= 0 or not 0 <= float(
        summary["max_energy_error_j"]
    ) <= 0.005 * float(summary["input_energy_j"]):
        raise ValueError("parent fails energy gate")
    if abs(float(summary["momentum_error_ns"])) > 1e-8 * max(
        1, abs(float(summary["wall_impulse_ns"]))
    ):
        raise ValueError("parent fails momentum gate")
    if float(summary["time_s"]) < EDGES[-1] or not snapshots:
        raise ValueError("parent does not cover audit window")
    initial = snapshots[0].cells
    if not initial or any(not math.isfinite(c.mass_kg) or c.mass_kg <= 0 for c in initial):
        raise ValueError("parent requires positive finite parcel masses")
    if any(a.time_s >= b.time_s for a, b in pairwise(snapshots)):
        raise ValueError("parent times must be strictly increasing")
    for snap in snapshots:
        if len(snap.cells) != len(initial) or any(
            not math.isclose(a.mass_kg, b.mass_kg, rel_tol=1e-12)
            for a, b in zip(initial, snap.cells, strict=True)
        ):
            raise ValueError("audit requires fixed parcel masses and indexing")


def cell_state(cell: FlowCell) -> tuple[float, ew.Composition | None]:
    if cell.temperature_k <= 1000:
        return 0.0, None  # Same cold molecular bookkeeping as flow.cell_stores.
    comp = ew.composition(cell.rho_kg_m3, cell.temperature_k)
    bond = ew._bond_energy_held(comp, cell.rho_kg_m3 / ew.M_H2O) / cell.rho_kg_m3
    return bond, comp


def diagnose(cell: FlowCell, comp: ew.Composition | None, weight: float) -> Row:
    row: Row = {
        "weight_j": weight,
        "temperature_k": cell.temperature_k,
        "pressure_pa": cell.pressure_pa,
        "rho_kg_m3": cell.rho_kg_m3,
        "x_m": cell.x_m,
        "status": "cold_unmodelled",
        "expansion_rate_s": cell.expansion_rate_s,
    }
    if comp is None:
        return row
    p_ideal = comp.n_total * ew.K_B * cell.temperature_k
    row["nonideal_pressure_fraction"] = abs(cell.pressure_pa - p_ideal) / cell.pressure_pa
    row["ionized_nuclei_fraction"] = (comp.n_hp + sum(comp.n_o_ions)) / (
        3 * cell.rho_kg_m3 / ew.M_H2O
    )
    # Retain tiny interval decreases in the denominator but do not solve a
    # spectrum merely to classify floating-point noise in fully atomized gas.
    if weight < 1e-10 * cell.mass_kg * ew.FULL_ATOMIZATION_ENERGY:
        row["status"] = "tiny_interval_change"
        return row
    try:
        result = relaxation(cell.temperature_k, neutral_densities(comp))
    except ValueError as error:
        row["status"] = "unresolved: " + str(error)
        return row
    row.update(
        {
            "status": "resolved_local_network",
            "energy_relaxation_s": result.energy_time_s,
            "water_collider_efficiency": result.water_collider_efficiency,
            "detailed_balance_log_error": result.max_detailed_balance_log_error,
        }
    )
    if not (math.isfinite(cell.expansion_rate_s) and cell.expansion_rate_s > 0):
        row["status"] = "resolved_no_expansion_clock"
        return row
    tau_old = recombination.atom_recombination_time(
        cell.temperature_k, comp.n_neutral_heavy, comp.n_oh
    )
    row["da_energy"] = 1 / (cell.expansion_rate_s * result.energy_time_s)
    row["da_old_oh_proxy"] = 1 / (cell.expansion_rate_s * tau_old)
    row["old_proxy_time_over_network_time"] = tau_old / result.energy_time_s
    return row


def summarize(rows: list[Row], reference_capacity_j: float = 0.0) -> Row:
    total = sum(float(r["weight_j"]) for r in rows)
    result: Row = {
        "positive_interval_bond_decrease_j": total,
        "kinetics_validated": False,
        "weighted_cell_intervals": len(rows),
        "status": "negligible_molecular_return"
        if total < 1e-6 * reference_capacity_j
        else "conditional_rate_screen",
    }
    flags: dict[str, Callable[[Row], bool]] = {
        "above_reference_thermo_3500k": lambda r: (
            float(r["temperature_k"]) > REFERENCE_THERMO_MAX_K
        ),
        "pressure_above_100mpa": lambda r: float(r["pressure_pa"]) > 1e8,
        "pressure_above_1gpa": lambda r: float(r["pressure_pa"]) > 1e9,
        "nonideal_pressure_gt10pct": lambda r: float(r.get("nonideal_pressure_fraction", 0)) > 0.1,
        "ionized_nuclei_gt1pct": lambda r: float(r.get("ionized_nuclei_fraction", 0)) > 0.01,
        "resolved_expansion_clock": lambda r: "da_energy" in r,
        "network_da_gt10": lambda r: "da_energy" in r and float(r["da_energy"]) > 10,
        "network_da_lt0p1": lambda r: "da_energy" in r and float(r["da_energy"]) < 0.1,
        "network_da_gt10000": lambda r: "da_energy" in r and float(r["da_energy"]) > 10000,
        "old_proxy_da_gt10": lambda r: "da_old_oh_proxy" in r and float(r["da_old_oh_proxy"]) > 10,
        "tracking_da_gt10": lambda r: "da_tracking" in r and float(r["da_tracking"]) > 10,
        "tracking_da_lt0p1": lambda r: "da_tracking" in r and float(r["da_tracking"]) < 0.1,
        "tracking_da_gt10000": lambda r: "da_tracking" in r and float(r["da_tracking"]) > 10000,
        "cold_unmodelled": lambda r: r["status"] == "cold_unmodelled",
        "tiny_interval_change": lambda r: r["status"] == "tiny_interval_change",
        "unresolved_network": lambda r: str(r["status"]).startswith("unresolved:"),
        "no_expansion_clock": lambda r: r["status"] == "resolved_no_expansion_clock",
    }
    for name, predicate in flags.items():
        value = sum(float(r["weight_j"]) for r in rows if predicate(r))
        result[name + "_weight_j"] = value
        result[name + "_fraction"] = value / total if total > 0 else ""
    for key in (
        "temperature_k",
        "pressure_pa",
        "rho_kg_m3",
        "x_m",
        "da_energy",
        "da_old_oh_proxy",
        "old_proxy_time_over_network_time",
        "water_collider_efficiency",
        "da_store_change",
        "da_tracking",
    ):
        for quantile in (0.1, 0.5, 0.9):
            result[f"{key}_p{int(100 * quantile)}"] = weighted_quantile(rows, key, quantile)
    return result


def audit(
    path: Path,
    stride: int,
    *,
    point_diagnostic: Callable[[FlowCell, ew.Composition | None, float], Row] = diagnose,
    summary_diagnostic: Callable[[list[Row], float], Row] = summarize,
) -> list[Row]:
    snapshots, summaries = read_output(path)
    by_name = {str(r["name"]): r for r in summaries}
    output = []
    for name in sorted({s.name for s in snapshots if not s.name.startswith("water_only")}):
        selected = select_snapshots([s for s in snapshots if s.name == name], stride)
        validate_parent(by_name[name], selected)
        capacity = sum(c.mass_kg for c in selected[0].cells) * ew.FULL_ATOMIZATION_ENERGY
        previous: list[float] | None = None
        previous_time = selected[0].time_s
        points: list[Row] = []
        for snap in selected:
            states = [cell_state(c) for c in snap.cells]
            if previous is not None:
                for old, (bond, comp), cell in zip(previous, states, snap.cells, strict=True):
                    weight = cell.mass_kg * max(0.0, old - bond)
                    if weight > 0:
                        row = point_diagnostic(cell, comp, weight)
                        row["time_s"] = snap.time_s
                        if old > bond >= 0:
                            row["store_change_time_s"] = store_change_clock(
                                old, bond, snap.time_s - previous_time
                            )
                            if "energy_relaxation_s" in row:
                                row["da_store_change"] = float(row["store_change_time_s"]) / float(
                                    row["energy_relaxation_s"]
                                )
                                if "da_energy" in row:
                                    row["da_tracking"] = min(
                                        float(row["da_energy"]), float(row["da_store_change"])
                                    )
                        points.append(row)
            previous = [s[0] for s in states]
            previous_time = snap.time_s
        for lo, hi in pairwise(EDGES):
            rows = [r for r in points if lo + 1e-12 < float(r["time_s"]) <= hi + 1e-12]
            output.append(
                {
                    "parent": path.stem,
                    "case": name,
                    "stride": stride,
                    "window_start_s": lo,
                    "window_end_s": hi,
                    **summary_diagnostic(rows, capacity),
                }
            )
        output.append(
            {
                "parent": path.stem,
                "case": name,
                "stride": stride,
                "window_start_s": EDGES[0],
                "window_end_s": EDGES[-1],
                **summary_diagnostic(points, capacity),
            }
        )
        print(
            f"{name}, stride {stride}: "
            f"{sum(float(r['weight_j']) for r in points) / 1e9:.6g} GJ positive interval decrease",
            flush=True,
        )
    return output


def parent_provenance(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    with path.open() as stream:
        metadata = json.loads(next(stream))
    if metadata["record"] != "metadata":
        raise ValueError("parent must start with its model metadata")
    return {"path": str(path), "sha256": digest, "metadata": metadata}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--strides", nargs="+", type=int, default=[1, 2])
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR / "kinetics_audit.csv")
    args = parser.parse_args()
    rows = [r for path in args.paths for stride in args.strides for r in audit(path, stride)]
    _write(args.output, rows)
    model = {
        "source_url": SOURCE_URL,
        "source_sha256": SOURCE_SHA256,
        "source_units": "cm, mol, s, cal/mol; converted to molecule, m, s",
        "retained_reactions": [asdict(r) for r in REACTIONS],
        "species_order": ["H", "O", "H2", "OH", "O2", "H2O"],
        "reference_thermo_max_k": REFERENCE_THERMO_MAX_K,
        "parents": [parent_provenance(p) for p in args.paths],
        "strides": args.strides,
        "window_edges_s": EDGES,
        "reverse_rates": "matched study partition functions, not GRI NASA7",
        "kinetics_validated": False,
        "limits": [
            "fixed T and rho; ions held fixed",
            "HO2/H2O2 routes absent",
            "no association falloff",
            "positive interval store decreases are not impulse or unique reaction energy",
            "3500 K is the source thermo-fit ceiling, not a kinetic validation boundary",
            "pressure thresholds are diagnostic bins, not validated rate limits",
            "tracking clock is min(expansion time, sampled mean held store / return rate)",
            "negligible return threshold is 1e-6 of case full-atomization capacity",
        ],
    }
    (OUTPUT_DIR / "kinetics_model.json").write_text(json.dumps(model, indent=2) + "\n")


if __name__ == "__main__":
    main()
