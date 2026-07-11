"""Tests for the 69 km/s freeze-timing bracket translation (`jupiter.py --frozen`, ADR-0026).

The `--frozen-jupiter` sweep measures an EOS-only e_eff bracket (equilibrium vs sudden-freeze vs
pure-H2O); this module translates that delta onto the coupled headline `f`. These tests pin the
translation math (stdlib + the analysis reconciliation only — no matplotlib/CoolProp)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from puffsat import jupiter


def test_read_frozen_jupiter_round_trips(tmp_path: Path) -> None:
    """The reader parses the `--frozen-jupiter` JSONL and tolerates a trailing blank line."""
    p = tmp_path / "frozen.jsonl"
    p.write_text(
        json.dumps(
            {
                "v": 69_000.0,
                "rho_impact": 0.04,
                "e_eff_eq": 0.70,
                "e_eff_frozen_rebound": 0.52,
                "e_eff_frozen_all": 0.646,
                "rho_star": 0.3,
                "t_star": 150_000.0,
                "swap_energy_jump_frac": 1e-3,
            }
        )
        + "\n\n"
    )
    rows = jupiter.read_frozen_jupiter(p)
    assert len(rows) == 1
    assert rows[0].rho_impact == 0.04
    assert rows[0].e_eff_frozen_rebound == 0.52


def test_frozen_bracket_subtracts_eos_delta_from_coupled_headline() -> None:
    """The freeze bracket applies the EOS-only delta to the *coupled* headline e_eff, not the
    EOS-only `e_eff_eq`: `f_eq` reconciles from `e_coupled`; the sudden-freeze `f` subtracts
    `e_eff_eq - e_eff_frozen_rebound` (positive), so it sits strictly below the headline; the
    pure-H2O bound (smaller delta here) sits between."""
    rows = [
        jupiter.FrozenJupiterRow(0.02, 0.70, 0.52, 0.646, 1e-3),
        jupiter.FrozenJupiterRow(0.08, 0.71, 0.53, 0.646, 1e-3),
    ]
    # Coupled headline e_eff(rho) is below the EOS-only e_eff_eq (radiation removes energy).
    e_coupled = jupiter._LogInterp([0.02, 0.08], [0.60, 0.65])
    eta = 0.9
    pts = jupiter.frozen_jupiter_bracket(rows, e_coupled, eta)

    assert [p.rho_impact for p in pts] == [0.02, 0.08]
    p0 = pts[0]
    # f_eq is reconciled from the coupled e_eff (0.60), not e_eff_eq (0.70).
    assert p0.e_eff_coupled == pytest.approx(0.60)
    assert p0.f_eq == pytest.approx(eta * (1.0 + 0.60) / 2.0)
    # delta_frozen = e_eff_eq - e_eff_frozen_rebound = 0.18, applied to the coupled e_eff.
    assert p0.delta_frozen == pytest.approx(0.18)
    assert p0.f_frozen_rebound == pytest.approx(eta * (1.0 + (0.60 - 0.18)) / 2.0)
    # Freeze penalty strictly lowers f; pure-H2O (delta 0.054) sits between rebound and headline.
    assert p0.f_frozen_rebound < p0.f_frozen_all < p0.f_eq


def test_run_frozen_jupiter_writes_summary_and_translates_design_point(tmp_path: Path) -> None:
    """End-to-end on the committed sweep artifacts: the `--frozen` path writes the CSV (header =
    `FrozenJupiterPoint` fields) and the sudden-freeze design-point f is strictly below the
    equilibrium headline. Skipped if the sweep artifacts are absent (they are gitignored bulk)."""
    frozen = Path("data/results/sweep_frozen_jupiter.jsonl")
    sweep = Path("data/results/sweep_jupiter.jsonl")
    geo = Path("data/results/sweep_geometry_m40.jsonl")
    for pth in (frozen, sweep, geo):
        if not pth.exists():
            pytest.skip(f"missing sweep artifact {pth} (run make sweep-frozen-jupiter)")

    out = tmp_path / "frontier_frozen_jupiter.csv"
    jupiter._run_frozen_jupiter(frozen, sweep, geo, out)
    assert out.exists()

    lines = out.read_text().splitlines()
    header = lines[0].split(",")
    assert header[0] == "rho_impact"
    assert {"f_eq", "f_frozen_rebound", "delta_frozen"} <= set(header)
    # Every swept density loses f to the freeze penalty (delta_frozen > 0 at 69 km/s).
    fe, fr, d = header.index("f_eq"), header.index("f_frozen_rebound"), header.index("delta_frozen")
    for row in lines[1:]:
        cells = row.split(",")
        assert float(cells[d]) > 0.0
        assert float(cells[fr]) < float(cells[fe])
