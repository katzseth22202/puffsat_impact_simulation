"""Ordered chemistry differences with explicit matched-parent and restart gates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from puffsat.water_plate.flow import _write
from puffsat.water_plate.freeze_report import Row, comparison, read_run, valid_quad
from puffsat.water_plate.thermodynamics import OUTPUT_DIR


def ordered_comparison(
    all_rows: list[Row], molecular_rows: list[Row], cold_rate: float, arrival: float
) -> Row:
    """Eq -> molecular-frozen/ion-active -> all-frozen; order is not unique."""
    result: Row = {"status": "invalid_or_unmatched"}
    if not valid_quad(all_rows, arrival) or not valid_quad(
        molecular_rows, arrival, "molecular_frozen"
    ):
        return result
    for case in ("combined", "unsprayed"):
        a = next(
            r
            for r in all_rows
            if str(r["name"]).startswith(case + "_") and r["branch"] == "equilibrium"
        )
        b = next(
            r
            for r in molecular_rows
            if str(r["name"]).startswith(case + "_") and r["branch"] == "equilibrium"
        )
        if a["name"] != b["name"] or any(
            not math.isclose(float(a[k]), float(b[k]), rel_tol=1e-10)
            for k in (
                "start_time_s",
                "time_s",
                "initial_wall_impulse_ns",
                "restart_physical_energy_j",
                "wall_impulse_ns",
            )
        ):
            return result
    total = comparison(all_rows, cold_rate, arrival)
    molecular = comparison(molecular_rows, cold_rate, arrival, "molecular_frozen")
    delta = float(total["reaction_effect_on_gain_ns"])
    mol = float(molecular["reaction_effect_on_gain_ns"])
    gain = float(total["equilibrium_gain_ns"])
    result.update(
        {
            "status": "matched",
            "switch_time_s": total["switch_time_s"],
            "end_time_s": total["end_time_s"],
            "equilibrium_gain_ns": gain,
            "all_chemistry_effect_ns": delta,
            "molecular_effect_with_ionization_active_ns": mol,
            "ionization_effect_with_molecules_frozen_ns": delta - mol,
            "molecular_fraction_of_equilibrium_gain": mol / gain if gain > 0 else "",
        }
    )
    return result


def matched_inputs(a: dict[str, object], b: dict[str, object]) -> bool:
    """Compare actual conserved arrays and integration grids, not just filenames."""
    if any(
        a[k] != b[k]
        for k in (
            "parent_path",
            "parent_metadata",
            "cold_wall_impulse_rate_ns_per_s",
            "cold_head_arrival_s",
        )
    ):
        return False
    # Untyped JSON is restricted to this comparison at the file boundary.
    left = json.loads(Path(str(a["input_path"])).read_text())
    right = json.loads(Path(str(b["input_path"])).read_text())
    if left["residual_curve"] != right["residual_curve"]:
        return False
    ra = {r["name"]: r for r in left["restarts"]}
    rb = {r["name"]: r for r in right["restarts"]}
    if ra.keys() != rb.keys():
        return False
    keys = (
        "state",
        "start_time_s",
        "initial_wall_impulse_ns",
        "area_m2",
        "storage_shift",
        "timestep_fraction",
        "output_times_s",
    )
    if not all(ra[n][k] == rb[n][k] for n in ra for k in keys):
        return False
    for name in ra:
        pa = [{k: v for k, v in p.items() if k != "ionizing"} for p in ra[name]["parcels"]]
        pb = [{k: v for k, v in p.items() if k != "ionizing"} for p in rb[name]["parcels"]]
        if pa != pb:
            return False
    return True


def report(all_paths: list[Path], molecular_paths: list[Path]) -> None:
    rows: list[Row] = []
    for all_path, molecular_path in zip(all_paths, molecular_paths, strict=True):
        _, a, ma, cold, arrival = read_run(all_path)
        _, b, mb, _, _ = read_run(molecular_path)
        if not matched_inputs(ma, mb):
            raise ValueError(
                "ordered decomposition requires identical parents, states and integration grids"
            )
        switches = sorted({float(r["start_time_s"]) for r in a})
        for switch in switches:
            row = ordered_comparison(
                [r for r in a if float(r["start_time_s"]) == switch],
                [r for r in b if float(r["start_time_s"]) == switch],
                cold,
                arrival,
            )
            row["all_frozen_run"] = all_path.stem
            row["molecular_frozen_run"] = molecular_path.stem
            rows.append(row)
            print(row)
    _write(OUTPUT_DIR / "selective_comparison.csv", rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-frozen", nargs="+", type=Path, required=True)
    parser.add_argument("--molecular-frozen", nargs="+", type=Path, required=True)
    args = parser.parse_args()
    report(args.all_frozen, args.molecular_frozen)


if __name__ == "__main__":
    main()
