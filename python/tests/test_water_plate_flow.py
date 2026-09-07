"""Acceptance checks for actual-history chemistry bookkeeping, not thrust attribution."""

import json
import math
from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import pytest

pytest.importorskip("CoolProp")

from puffsat.water_plate.flow import (
    CsvValue,
    FlowCell,
    Snapshot,
    cell_stores,
    history_rows,
    matched_gain,
    reaction_damkohler,
)


def gas(temp: float = 10000.0, rate: float = 1000.0) -> FlowCell:
    return FlowCell(1.0, 0.1, 2.0, 10.0, temp, 0.0, 0.0, 0.0, rate)


def test_parcel_peaks_are_separate_by_cell_and_control_and_mass_weighted() -> None:
    hot, cool = gas(16000), gas(5000)
    hb, hi = cell_stores(hot)
    cb, ci = cell_stores(cool)
    snapshots = [
        Snapshot("combined", 0.0, 0.0, 1, [hot, cool]),
        Snapshot("water_only", 0.0, 0.0, 1, [cool]),
        Snapshot("combined", 1e-3, 1e6, 1, [cool, hot]),
    ]
    rows = history_rows(snapshots)
    assert rows[0]["bond_store_j"] == pytest.approx(2 * (hb + cb))
    assert rows[0]["water_bond_store_j"] == pytest.approx(2 * hb)
    assert rows[1]["bond_drop_from_sampled_parcel_peaks_j"] == 0
    assert rows[2]["bond_drop_from_sampled_parcel_peaks_j"] == pytest.approx(2 * (hb - cb))
    assert rows[2]["ion_drop_from_sampled_parcel_peaks_j"] == pytest.approx(2 * (hi - ci))
    assert rows[2]["wall_impulse_ns"] == 1e6  # Measured, not inferred from store decrease.
    assert rows[2]["bond_decrease_midpoint_x_m"] == 1.0


@pytest.mark.parametrize("rate", [math.nan, 0.0, -1000.0])
def test_no_expansion_clock_is_invented_for_missing_stationary_or_compressing_cells(
    rate: float,
) -> None:
    assert reaction_damkohler(gas(rate=rate)) is None


def test_reaction_clocks_scale_with_actual_velocity_divergence() -> None:
    slow = reaction_damkohler(gas(rate=1000))
    fast = reaction_damkohler(gas(rate=10000))
    assert slow is not None and fast is not None
    assert slow[0] == pytest.approx(10 * fast[0])
    assert slow[1] == pytest.approx(10 * fast[1])
    assert reaction_damkohler(gas(temp=300)) is None


def test_freeze_screen_is_weighted_by_held_store_only_in_expanding_gas() -> None:
    expanding = gas(rate=1e50)
    compressed = replace(expanding, expansion_rate_s=-1000)
    bond, ion = cell_stores(expanding)
    row = history_rows([Snapshot("test", 0, 0, 0, [expanding, compressed])])[0]
    assert row["bond_store_j"] == pytest.approx(4 * bond)
    assert row["expanding_gas_bond_store_j"] == pytest.approx(2 * bond)
    assert row["expanding_bond_store_da_oh_lt0p1_j"] == pytest.approx(2 * bond)
    assert row["expanding_ion_store_da_coll_lt0p1_j"] == pytest.approx(2 * ion)
    assert row["expanding_bond_store_da_oh_gt10_j"] == 0


def controls() -> list[dict[str, CsvValue]]:
    return [
        {
            "name": f"{name}_h1_n20",
            "status": "completed_time_window",
            "wall_impulse_ns": impulse,
            "input_energy_j": 1e9,
            "max_energy_error_j": 0.0,
            "momentum_error_ns": 0.0,
            "time_s": 0.001,
            "area_m2": 40.0,
            "layer_length_m": 1.0,
        }
        for name, impulse in (("combined", 3e6), ("unsprayed", 2e6), ("water_only", 100.0))
    ]


def test_gain_requires_all_three_matched_valid_conserved_controls() -> None:
    rows = controls()
    assert matched_gain(rows) == 999900.0
    assert matched_gain(rows[:2]) is None
    for field, value in (
        ("status", "invalid_table_domain"),
        ("time_s", 0.0009),
        ("area_m2", 20.0),
        ("max_energy_error_j", 1e7),
        ("momentum_error_ns", 100.0),
    ):
        changed = [dict(row) for row in rows]
        changed[0][field] = value
        assert matched_gain(changed) is None


def test_analytic_cold_control_is_accepted_only_before_signal_arrival() -> None:
    rows = controls()
    rows[2]["status"] = "completed_analytic_wall_window"
    rows[2]["initial_rarefaction_arrival_s"] = 0.1
    del rows[2]["max_energy_error_j"]
    del rows[2]["momentum_error_ns"]
    assert matched_gain(rows) == 999900.0
    rows[2]["initial_rarefaction_arrival_s"] = 0.001
    assert matched_gain(rows) is None


def test_report_milestones_are_sampled_and_missing_returns_are_explicit() -> None:
    from puffsat.water_plate.flow_report import milestones

    rows = [
        {
            "time_s": str(t),
            "wall_impulse_ns": str(j),
            "bond_store_j": str(b),
            "ion_store_j": str(i),
            "bond_drop_from_sampled_parcel_peaks_j": str(10 - b),
            "ion_drop_from_sampled_parcel_peaks_j": str(20 - i),
        }
        for t, j, b, i in ((0, 0, 10, 20), (1, 5, 10, 10), (2, 10, 8, 5))
    ]
    events = {row["event"]: row for row in milestones(rows)}
    assert events["ion_parcel_peak_drop_0.5"]["window_impulse_fraction"] == 0.5
    assert events["bond_parcel_peak_drop_0.1"]["time_s"] == "2"
    assert events["bond_parcel_peak_drop_0.1"]["previous_sample_time_s"] == "1"
    assert events["bond_parcel_peak_drop_0.5"]["status"] == "not_reached"


def test_control_inputs_share_energy_zero_and_resolve_the_coupling_window(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from puffsat.water_plate import flow, flow_eos

    monkeypatch.setattr(flow, "OUTPUT_DIR", tmp_path)
    configs = json.loads(flow.write_inputs(40, 1.0, coupling_samples=301).read_text())
    assert len(configs) == 3
    for config in configs:
        assert config["energy_shift_j_kg"] == flow_eos.ENERGY_SHIFT
        assert config["timestep_fraction"] == 0.25
        times = config["output_times_s"]
        assert times[0] == 0.0 and times[-1] == 0.001
        assert all(a < b for a, b in pairwise(times))
        assert sum(100e-6 <= t <= 250e-6 for t in times) >= 301
