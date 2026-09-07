"""Matched reaction-switch effects; invalid domains never become impulse differences."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import cast

from puffsat.water_plate.flow import CsvValue, _write
from puffsat.water_plate.thermodynamics import OUTPUT_DIR

Row = dict[str, CsvValue]


def read_run(path: Path) -> tuple[list[Row], list[Row], dict[str, object], float, float]:
    history, summaries = [], []
    input_path = Path()
    with path.open() as stream:
        for line in stream:
            raw = json.loads(line)
            if raw["record"] == "metadata":
                input_path = Path(raw["input_path"])
            elif raw["record"] == "summary":
                summaries.append(cast(Row, raw))
            elif raw["record"] == "history":
                history.append(cast(Row, raw))
    inputs = json.loads(input_path.read_text())
    parent = Path(inputs["parent_path"])
    cold_rate = arrival = math.nan
    with parent.open() as stream:
        for line in stream:
            # Only summaries are needed from the large parent trajectory.
            if '"record":"summary"' not in line:
                continue
            raw = json.loads(line)
            if (
                raw["name"].startswith("water_only")
                and raw["status"] == "completed_analytic_wall_window"
            ):
                cold_rate = float(raw["final_wall_pressure_pa"]) * float(raw["area_m2"])
                arrival = float(raw["initial_rarefaction_arrival_s"])
    if not math.isfinite(cold_rate) or not math.isfinite(arrival):
        raise ValueError("missing valid cold-only wall control in parent")
    provenance: dict[str, object] = {
        "input_path": str(input_path),
        "parent_path": str(parent),
        "parent_metadata": inputs["parent_metadata"],
        "residual_knots": len(inputs["residual_curve"]["knots"]),
        "cold_wall_impulse_rate_ns_per_s": cold_rate,
        "cold_head_arrival_s": arrival,
    }
    return history, summaries, provenance, cold_rate, arrival


def valid_quad(rows: list[Row], arrival: float) -> bool:
    """Require matched controls, conservation and small matching corrections."""
    keys = {(str(r["name"]).split("_h")[0], r["branch"]) for r in rows}
    if len(rows) != 4 or keys != {
        (c, b) for c in ("combined", "unsprayed") for b in ("equilibrium", "parcel_frozen")
    }:
        return False
    if len({str(r["name"]).split("_h", 1)[1] for r in rows}) != 1:
        return False
    for row in rows:
        if row["status"] != "completed_time_window":
            return False
        values = [
            float(row[k])
            for k in (
                "time_s",
                "start_time_s",
                "wall_impulse_ns",
                "max_energy_error_j",
                "restart_physical_energy_j",
                "momentum_error_ns",
                "max_matching_pressure_correction_fraction",
                "abs_matching_energy_offsets_j",
                "max_start_pressure_jump_fraction",
                "max_start_temperature_jump_fraction",
            )
        ]
        if not all(math.isfinite(v) for v in values):
            return False
        if not float(row["time_s"]) < arrival:
            return False
        scale = abs(float(row["restart_physical_energy_j"]))
        if float(row["max_energy_error_j"]) > 0.005 * scale:
            return False
        if float(row["abs_matching_energy_offsets_j"]) > 0.005 * scale:
            return False
        if float(row["max_matching_pressure_correction_fraction"]) > 0.01:
            return False
        if abs(float(row["momentum_error_ns"])) > 1e-8 * max(
            1.0, abs(float(row["wall_impulse_ns"]))
        ):
            return False
        if any(
            float(row[k]) > 1e-8
            for k in ("max_start_pressure_jump_fraction", "max_start_temperature_jump_fraction")
        ):
            return False
        if any(
            not math.isclose(float(row[k]), float(rows[0][k]), rel_tol=1e-12)
            for k in ("time_s", "start_time_s")
        ):
            return False
    for case in ("combined", "unsprayed"):
        pair = [r for r in rows if str(r["name"]).startswith(case + "_")]
        if any(
            not math.isclose(float(pair[0][k]), float(pair[1][k]), rel_tol=1e-12)
            for k in ("initial_wall_impulse_ns", "restart_physical_energy_j")
        ):
            return False
    return True


def comparison(rows: list[Row], cold_rate: float, arrival: float) -> Row:
    result: Row = {"switch_time_s": rows[0]["start_time_s"], "status": "invalid_or_unmatched"}
    if not valid_quad(rows, arrival):
        return result
    by = {(str(r["name"]).split("_h")[0], r["branch"]): float(r["wall_impulse_ns"]) for r in rows}
    ce, cf, ue, uf = (
        by[k]
        for k in (
            ("combined", "equilibrium"),
            ("combined", "parcel_frozen"),
            ("unsprayed", "equilibrium"),
            ("unsprayed", "parcel_frozen"),
        )
    )
    end = float(rows[0]["time_s"])
    cold = cold_rate * end
    delta_combined, delta_unsprayed = ce - cf, ue - uf
    result.update(
        {
            "status": "matched",
            "end_time_s": end,
            "combined_equilibrium_ns": ce,
            "combined_frozen_ns": cf,
            "unsprayed_equilibrium_ns": ue,
            "unsprayed_frozen_ns": uf,
            "cold_only_ns": cold,
            "equilibrium_gain_ns": ce - ue - cold,
            "frozen_gain_ns": cf - uf - cold,
            "reaction_effect_combined_ns": delta_combined,
            "reaction_effect_unsprayed_ns": delta_unsprayed,
            "reaction_effect_on_gain_ns": delta_combined - delta_unsprayed,
            "fraction_of_equilibrium_gain": (
                (delta_combined - delta_unsprayed) / (ce - ue - cold) if ce - ue - cold > 0 else ""
            ),
        }
    )
    return result


def report(paths: list[Path], reference: Path) -> None:
    all_runs: list[Row] = []
    comparisons: list[Row] = []
    reference_history: list[Row] = []
    provenance: dict[str, object] = {}
    for path in paths:
        history, summaries, meta, cold_rate, arrival = read_run(path)
        for switch in sorted({float(r["start_time_s"]) for r in summaries}):
            quad = [r for r in summaries if float(r["start_time_s"]) == switch]
            row = comparison(quad, cold_rate, arrival)
            row["run"] = path.stem
            comparisons.append(row)
            if path == reference and row["status"] != "matched":
                raise ValueError("reference must have valid matched branches")
        for row in summaries:
            row["run"] = path.stem
        all_runs.extend(summaries)
        if path == reference:
            reference_history = history
            provenance = meta
    if not reference_history:
        raise ValueError("reference must be included in paths")
    _write(OUTPUT_DIR / "freeze_runs.csv", all_runs)
    _write(OUTPUT_DIR / "freeze_comparison.csv", comparisons)
    (OUTPUT_DIR / "freeze_model.json").write_text(json.dumps(provenance, indent=2) + "\n")
    write_figure(reference_history, reference.stem)
    for row in comparisons:
        print(row)


def write_figure(history: list[Row], name: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    switches = sorted({str(r["name"]).split("_freeze_")[1] for r in history})
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, layout="constrained")
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
    for index, switch in enumerate(switches):
        color = colors[index % len(colors)]

        def selected(case: str, branch: str, current_switch: str = switch) -> list[Row]:
            return [
                r
                for r in history
                if str(r["name"]).startswith(case + "_")
                and str(r["name"]).endswith("_freeze_" + current_switch)
                and r["branch"] == branch
            ]

        ce, cf, ue, uf = (
            selected(c, b)
            for c, b in (
                ("combined", "equilibrium"),
                ("combined", "parcel_frozen"),
                ("unsprayed", "equilibrium"),
                ("unsprayed", "parcel_frozen"),
            )
        )
        times = [1e3 * float(r["time_s"]) for r in ce]
        for rows, style, branch in ((ce, "-", "equilibrium"), (cf, "--", "frozen")):
            axes[0].plot(
                times,
                [float(r["wall_impulse_ns"]) / 1e6 for r in rows],
                style,
                color=color,
                label=f"{branch}, switch {switch}",
            )
        delta = []
        for a, b, c, d in zip(ce, cf, ue, uf, strict=True):
            if any(
                not math.isclose(float(r["time_s"]), float(a["time_s"]), rel_tol=1e-12)
                for r in (b, c, d)
            ):
                raise ValueError("history samples do not match")
            delta.append(
                (
                    float(a["wall_impulse_ns"])
                    - float(b["wall_impulse_ns"])
                    - float(c["wall_impulse_ns"])
                    + float(d["wall_impulse_ns"])
                )
                / 1000
            )
        axes[1].plot(times, delta, color=color, label=f"switch {switch}")
    axes[0].set_ylabel("Combined wall impulse (MN s)")
    axes[1].set_ylabel("Reaction effect on water-layer\nadvantage (kN s)")
    axes[1].set_xlabel("Time (ms)")
    axes[0].set_title(
        "Same-state reaction switch: all later chemistry active versus frozen\n" + name
    )
    for axis in axes:
        axis.legend()
        axis.grid(True, alpha=0.2)
    fig.savefig(OUTPUT_DIR / "freeze_comparison.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    args = parser.parse_args()
    report(args.paths, args.reference)


if __name__ == "__main__":
    main()
