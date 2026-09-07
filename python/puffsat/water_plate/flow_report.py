"""Small, citable evidence from the planar flow histories; no chemistry thrust attribution."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from puffsat.water_plate.flow import CsvValue, _write, matched_gain
from puffsat.water_plate.thermodynamics import OUTPUT_DIR


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def milestones(rows: list[dict[str, str]]) -> list[dict[str, CsvValue]]:
    """Sampled milestones, normalized to this finite window, with explicit time brackets."""
    final = rows[-1]
    final_impulse = float(final["wall_impulse_ns"])
    events = {
        "peak_total_bond_store": max(
            range(len(rows)), key=lambda i: float(rows[i]["bond_store_j"])
        ),
        "peak_total_ion_store": max(range(len(rows)), key=lambda i: float(rows[i]["ion_store_j"])),
        "window_end": len(rows) - 1,
    }
    for fraction in (0.1, 0.5, 0.9):
        events[f"impulse_{fraction:g}_of_window"] = next(
            i
            for i, row in enumerate(rows)
            if float(row["wall_impulse_ns"]) >= fraction * final_impulse
        )
    for store in ("bond", "ion"):
        field = f"{store}_drop_from_sampled_parcel_peaks_j"
        scale = float(final[f"{store}_store_j"]) + float(final[field])
        for fraction in (0.1, 0.5):
            events[f"{store}_parcel_peak_drop_{fraction:g}"] = next(
                (
                    i
                    for i, row in enumerate(rows)
                    if scale > 0 and float(row[field]) >= fraction * scale
                ),
                -1,
            )
    result: list[dict[str, CsvValue]] = []
    for event, index in events.items():
        row: dict[str, CsvValue] = {"event": event, "status": "not_reached"}
        if index >= 0:
            row.update(rows[index])
            row["status"] = "sampled"
            row["window_impulse_fraction"] = float(rows[index]["wall_impulse_ns"]) / final_impulse
            row["previous_sample_time_s"] = rows[max(0, index - 1)]["time_s"]
        result.append(row)
    return result


def report(paths: list[Path], reference: Path) -> None:
    comparisons: list[dict[str, CsvValue]] = []
    reference_valid = False
    for path in paths:
        summaries: list[dict[str, CsvValue]] = [
            dict(row) for row in read_csv(path.with_name(path.stem + "_summary.csv"))
        ]
        gain = matched_gain(summaries)
        sampled = read_csv(path.with_name(path.stem + "_history.csv"))
        if path == reference:
            reference_valid = gain is not None
        for row in summaries:
            case_history = [r for r in sampled if r["name"] == row["name"]]
            if case_history:
                last = case_history[-1]
                row["snapshots"] = len(case_history)
                for store in ("bond", "ion"):
                    row[f"end_{store}_parcel_peak_drop_j"] = last[
                        f"{store}_drop_from_sampled_parcel_peaks_j"
                    ]
                    row[f"peak_total_{store}_store_j"] = max(
                        float(r[f"{store}_store_j"]) for r in case_history
                    )
                for field in (
                    "max_water_density_kg_m3",
                    "max_water_temperature_k",
                    "hot_dense_mass_kg_Tge5000_rhoge10",
                    "dense_closure_mass_kg_Tgt1000_rhoge100",
                    "high_pressure_mass_kg_pgt1GPa",
                ):
                    row[f"sampled_peak_{field}"] = max(float(r[field]) for r in case_history)
            row["run"] = path.stem
            row["matched_gain_ns"] = gain if gain is not None else ""
            if row.get("max_energy_error_j", "") != "":
                row["relative_max_energy_error"] = float(row["max_energy_error_j"]) / abs(
                    float(row["input_energy_j"])
                )
            comparisons.append(row)
    if not reference_valid:
        raise ValueError("reference must be a valid matched comparison included in paths")
    _write(OUTPUT_DIR / "flow_reference_convergence.csv", comparisons)
    history = read_csv(reference.with_name(reference.stem + "_history.csv"))
    combined = [row for row in history if row["name"].startswith("combined_")]
    events = milestones(combined)
    for event in events:
        event["run"] = reference.stem
    _write(OUTPUT_DIR / "flow_reference_milestones.csv", events)
    with reference.open() as stream:
        metadata = json.loads(stream.readline())
    if metadata["record"] != "metadata":
        raise ValueError("reference output must include its table and input provenance")
    metadata["reference_run"] = reference.stem
    (OUTPUT_DIR / "flow_reference_model.json").write_text(json.dumps(metadata, indent=2) + "\n")
    write_figure(history, reference.stem)


def write_figure(history: list[dict[str, str]], label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True, layout="constrained")
    combined = [row for row in history if row["name"].startswith("combined_")]
    for case, color in (("combined", "#0072B2"), ("unsprayed", "#D55E00")):
        selected = [row for row in history if row["name"].startswith(case + "_")]
        axes[0].plot(
            [1000 * float(r["time_s"]) for r in selected],
            [float(r["wall_impulse_ns"]) / 1e6 for r in selected],
            label=case,
            color=color,
        )
    times = [1000 * float(row["time_s"]) for row in combined]
    for store, color in (("bond", "#009E73"), ("ion", "#CC79A7")):
        axes[1].plot(
            times, [float(r[f"{store}_store_j"]) / 1e9 for r in combined], color=color, label=store
        )
        axes[2].plot(
            times,
            [float(r[f"{store}_drop_from_sampled_parcel_peaks_j"]) / 1e9 for r in combined],
            color=color,
            label=store,
        )
    final_j = float(combined[-1]["wall_impulse_ns"])
    t90 = next(
        float(row["time_s"]) for row in combined if float(row["wall_impulse_ns"]) >= 0.9 * final_j
    )
    for axis in axes:
        axis.axvline(1000 * t90, color="#666666", linestyle="--", linewidth=1)
        axis.grid(True, alpha=0.2)
        axis.legend()
    axes[0].set_ylabel("Wall impulse (MN s)")
    axes[1].set_ylabel("Combined: held stores (GJ)")
    axes[2].set_ylabel("Drop from sampled\nparcel peaks (GJ)")
    axes[2].set_xlabel("Time (ms); dashed line: sampled 90% of 1 ms wall impulse")
    axes[0].set_title(
        "Planar stratified-water reference: equilibrium chemistry, fixed lossless wall\n" + label
    )
    fig.savefig(OUTPUT_DIR / "flow_reference.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    args = parser.parse_args()
    report(args.paths, args.reference)


if __name__ == "__main__":
    main()
