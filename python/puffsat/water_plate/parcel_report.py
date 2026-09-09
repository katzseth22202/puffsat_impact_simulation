"""Net bond-store return from nonlinear prescribed-history parcels, not thrust."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from puffsat import eos_water as ew
from puffsat.water_plate.flow import CsvValue, _write
from puffsat.water_plate.parcel import ENDS


def summarize(records: list[dict[str, Any]], end: float) -> dict[str, CsvValue]:
    """Any is confined to the local JSON boundary; missing parcels keep their mass/store."""
    total_mass = sum(float(r["mass_kg"]) for r in records)
    initial = sum(float(r["mass_kg"]) * float(r["initial_bond_j_kg"]) for r in records)
    if not records or total_mass <= 0 or initial <= 0:
        raise ValueError("need a positive ensemble mass and initial store")
    resolved_mass = resolved_store = returned = parent = work = forcing = 0.0
    unresolved_lower = unresolved_upper = 0.0
    max_energy_error = 0.0
    max_energy_relative_error = 0.0
    for record in records:
        mass = float(record["mass_kg"])
        held = float(record["initial_bond_j_kg"])
        samples = record["trajectory"]["samples"]
        found = [s for s in samples if math.isclose(float(s["time_s"]), end, abs_tol=1e-12)]
        parents = [
            s
            for s in record["parent_returns"]
            if math.isclose(float(s["time_s"]), end, abs_tol=1e-12)
        ]
        if len(parents) != 1:
            raise ValueError("missing or duplicate parent endpoint")
        parent += mass * float(parents[0]["net_bond_return_j_kg"])
        error = float(record["max_energy_error_j_kg"])
        relative = error / max(abs(float(record["initial_energy_j_kg"])), 1.0)
        max_energy_error = max(max_energy_error, error)
        max_energy_relative_error = max(max_energy_relative_error, relative)
        if len(found) > 1:
            raise ValueError("duplicate parcel endpoint")
        if found and math.isfinite(relative) and relative <= 1e-6:
            sample = found[0]
            net = float(sample["net_bond_return_j_kg"])
            if (
                not math.isfinite(net)
                or not -1e-3 <= held - net <= ew.FULL_ATOMIZATION_ENERGY + 1e-3
            ):
                raise ValueError("invalid completed parcel bond store")
            resolved_mass += mass
            resolved_store += mass * held
            returned += mass * net
            work += mass * float(sample["work_j_kg"])
            forcing += mass * float(sample["forcing_j_kg"])
        else:
            # Very broad bookkeeping bounds for missing endpoint stores, not
            # a physical rate uncertainty or a confidence interval.
            unresolved_lower += mass * (held - ew.FULL_ATOMIZATION_ENERGY)
            unresolved_upper += mass * held
    return {
        "end_time_s": end,
        "kinetics_validated": False,
        "total_mass_kg": total_mass,
        "initial_bond_store_j": initial,
        "mass_mean_initial_ionized_nuclei_fraction": sum(
            float(r["mass_kg"]) * float(r["initial_ionized_nuclei_fraction"]) for r in records
        )
        / total_mass,
        "resolved_mass_fraction": resolved_mass / total_mass,
        "resolved_initial_store_fraction": resolved_store / initial,
        "resolved_net_return_j": returned,
        "resolved_return_over_total_start_store": returned / initial,
        "including_unresolved_return_lower_j": returned + unresolved_lower,
        "including_unresolved_return_upper_j": returned + unresolved_upper,
        "parent_equilibrium_net_return_j": parent,
        "resolved_work_j": work,
        "resolved_parent_forcing_j": forcing,
        "max_energy_error_j_kg": max_energy_error,
        "max_energy_relative_error": max_energy_relative_error,
        "accepted_steps": sum(int(r["trajectory"]["accepted_steps"]) for r in records),
        "rejected_steps": sum(int(r["trajectory"]["rejected_steps"]) for r in records),
    }


def report(paths: list[Path], output: Path) -> None:
    summaries = []
    provenance = []
    for path in paths:
        with path.open() as stream:
            metadata = json.loads(next(stream))
            records = [json.loads(line) for line in stream]
        input_path = Path(metadata["input_path"])
        with input_path.open() as stream:
            payload = json.load(stream)
        if any(r["record"] != "parcel" for r in records):
            raise ValueError("unexpected output record")
        expected = {int(p["index"]) for p in payload["parcels"]}
        for variant in payload["variants"]:
            for forcing in payload["forcing_modes"]:
                selected = [
                    r
                    for r in records
                    if r["variant"] == variant["name"] and r["parent_forcing"] == forcing
                ]
                if (
                    len(selected) != len(expected)
                    or {int(r["index"]) for r in selected} != expected
                ):
                    raise ValueError("incomplete or duplicated scenario ensemble")
                for end in ENDS:
                    summaries.append(
                        {
                            "run": path.stem,
                            "variant": str(variant["name"]),
                            "parent_forcing": bool(forcing),
                            "history_stride": int(payload["stride"]),
                            "relative_tolerance": float(payload["settings"]["relative_tolerance"]),
                            "maximum_step_s": float(payload["settings"]["maximum_step_s"]),
                            **summarize(selected, end),
                        }
                    )
        hashes = {}
        for p in (path, input_path):
            with p.open("rb") as stream:
                hashes[str(p)] = hashlib.file_digest(stream, "sha256").hexdigest()
        provenance.append(
            {
                "files": hashes,
                "parent": payload["parent"],
                "settings": payload["settings"],
                "stride": payload["stride"],
                "variants": payload["variants"],
                "source_url": payload["source_url"],
                "source_sha256": payload["source_sha256"],
                "li_source_url": payload["li_source_url"],
                "li_source_sha256": payload["li_source_sha256"],
                "limits": payload["limits"],
            }
        )
    _write(output, summaries)
    output.with_suffix(".json").write_text(
        json.dumps(
            {"kinetics_validated": False, "runs": provenance, "net_return_not_impulse": True},
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("data/results/water_plate/parcel_return.csv")
    )
    args = parser.parse_args()
    report(args.paths, args.output)


if __name__ == "__main__":
    main()
