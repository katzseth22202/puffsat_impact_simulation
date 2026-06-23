"""Tests for the opacity-insensitivity comparison core (B5d-3).

`compare` is the pure, tested heart of the gate: given per-scale sweep rows it measures how far
`e_eff` moves from the 1x reference and decides whether the interim opacity is justified. The scan
itself (`run_scan`, which shells out to the Rust sweep) is the experiment, exercised by
`make sensitivity`, not here.
"""

from __future__ import annotations

import pytest

from puffsat.analysis import SweepRow
from puffsat.sensitivity import PASS_THRESHOLD, compare


def _rows(e_eff_by_rho: dict[float, float], total_loss: float) -> list[SweepRow]:
    """Build sweep rows at `e_eff(rho)`, putting all loss in channel 1a (radiative-wall)."""
    return [
        SweepRow(
            rho_impact=rho,
            v=16_000.0,
            e_eff=e_eff,
            loss_radiative_wall=total_loss,
            loss_escape_space=0.0,
            loss_conductive=0.0,
            loss_condensation=0.0,
        )
        for rho, e_eff in e_eff_by_rho.items()
    ]


def test_compare_measures_small_spread_and_passes() -> None:
    """Near-identical `e_eff` across scales (while losses scale) -> small spread, gate passes."""
    rows_by_scale = {
        0.1: _rows({0.16: 0.628, 0.64: 0.636}, total_loss=10.0),
        1.0: _rows({0.16: 0.629, 0.64: 0.637}, total_loss=100.0),
        10.0: _rows({0.16: 0.630, 0.64: 0.638}, total_loss=1000.0),
    }
    result = compare(rows_by_scale)

    assert result.rho == [0.16, 0.64]  # sorted ascending
    assert result.scales == [0.1, 1.0, 10.0]
    # max excursion from the 1x reference is 0.001 (e.g. 0.630 vs 0.629).
    assert result.max_abs_de_eff == pytest.approx(0.001, abs=1e-9)
    assert result.max_rel_de_eff < PASS_THRESHOLD
    assert result.passes
    # The total loss tracks the opacity scale (1a moved as expected).
    assert result.total_loss_by_scale[10.0][0] == 1000.0
    assert result.total_loss_by_scale[0.1][0] == 10.0


def test_compare_flags_large_spread() -> None:
    """A scale that shifts `e_eff` past the threshold fails the gate."""
    rows_by_scale = {
        1.0: _rows({0.16: 0.60}, total_loss=100.0),
        10.0: _rows({0.16: 0.40}, total_loss=1000.0),  # 33% drop
    }
    result = compare(rows_by_scale)
    assert result.max_abs_de_eff == pytest.approx(0.20)
    assert result.max_rel_de_eff == pytest.approx(0.20 / 0.60)
    assert not result.passes


def test_compare_requires_reference_scale() -> None:
    """Without the 1x reference there is nothing to compare against."""
    with pytest.raises(ValueError, match="reference scale"):
        compare({0.1: _rows({0.16: 0.6}, 1.0), 10.0: _rows({0.16: 0.6}, 1.0)})


def test_compare_requires_matching_rho_grid() -> None:
    """Scales must share the same rho grid, else the per-rho comparison is ill-defined."""
    rows_by_scale = {
        1.0: _rows({0.16: 0.6, 0.64: 0.6}, 1.0),
        10.0: _rows({0.16: 0.6, 0.32: 0.6}, 1.0),  # different second rho
    }
    with pytest.raises(ValueError, match="rho grid"):
        compare(rows_by_scale)
