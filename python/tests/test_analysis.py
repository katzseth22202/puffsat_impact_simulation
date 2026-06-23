"""Tests for the sweep frontier extraction (B5d-2).

The extraction is the load-bearing, science-dep-free core: it parses the JSONL schema the Rust
sweep writes, sorts/normalizes into the `e_eff(rho)` frontier + loss decomposition, and writes the
CSV the figures are built from. Plotting (matplotlib, `sci` extra) gets a smoke test that is skipped
when the extra is absent.
"""

from __future__ import annotations

import csv
import json
from dataclasses import fields
from pathlib import Path

import pytest

from puffsat import analysis as an


def _write_jsonl(path: Path, rows: list[dict[str, float]]) -> None:
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")


def _row(
    rho: float, e_eff: float, la: float, lb: float, lc: float, ld: float = 0.0
) -> dict[str, float]:
    return {
        "rho_impact": rho,
        "v": 16_000.0,
        "e_eff": e_eff,
        "peak_wall_force": 1.0,  # extra fields the reader ignores
        "incident_momentum": 1.0,
        "residual_momentum": 0.5,
        "wall_impulse": 1.5,
        "loss_radiative_wall": la,
        "loss_escape_space": lb,
        "loss_conductive": lc,
        "loss_condensation": ld,
    }


def test_read_sweep_parses_jsonl(tmp_path: Path) -> None:
    """`read_sweep` reads one object per line (tolerating a trailing blank) into `SweepRow`s."""
    path = tmp_path / "sweep.jsonl"
    path.write_text(
        json.dumps(_row(0.32, 0.63, 100.0, 30.0, 0.0)) + "\n\n"  # trailing blank line
    )
    rows = an.read_sweep(path)
    assert len(rows) == 1
    assert rows[0].rho_impact == 0.32
    assert rows[0].e_eff == 0.63
    assert rows[0].loss_radiative_wall == 100.0
    assert rows[0].loss_conductive == 0.0
    assert rows[0].loss_condensation == 0.0


def test_frontier_sorts_and_normalizes(tmp_path: Path) -> None:
    """Frontier is ascending in rho; loss fractions sum to 1 with loss, and are 0 without."""
    path = tmp_path / "sweep.jsonl"
    _write_jsonl(
        path,
        [
            _row(0.64, 0.636, 300.0, 100.0, 0.0),  # out of order on purpose
            _row(0.16, 0.628, 0.0, 0.0, 0.0),  # no loss
            _row(0.32, 0.633, 60.0, 20.0, 20.0),
        ],
    )
    pts = an.frontier(an.read_sweep(path))

    assert [p.rho_impact for p in pts] == [0.16, 0.32, 0.64]  # sorted ascending

    # zero-loss point: total 0, all fractions 0.
    p0 = pts[0]
    assert p0.total_loss == 0.0
    fracs0 = (
        p0.frac_radiative_wall,
        p0.frac_escape_space,
        p0.frac_conductive,
        p0.frac_condensation,
    )
    assert fracs0 == (0.0, 0.0, 0.0, 0.0)

    # loss-bearing points: the four channel fractions sum to 1 and match the channel split.
    for p in pts[1:]:
        assert p.total_loss > 0.0
        s = p.frac_radiative_wall + p.frac_escape_space + p.frac_conductive + p.frac_condensation
        assert s == pytest.approx(1.0)
    mid = pts[1]  # rho=0.32: 60/100, 20/100, 20/100, 0
    assert mid.frac_radiative_wall == pytest.approx(0.6)
    assert mid.frac_escape_space == pytest.approx(0.2)
    assert mid.frac_conductive == pytest.approx(0.2)
    assert mid.frac_condensation == pytest.approx(0.0)


def test_frontier_condensation_channel(tmp_path: Path) -> None:
    """A low-v row whose only loss is condensation (channel 3) reports `frac_condensation == 1`."""
    path = tmp_path / "sweep_lowv.jsonl"
    _write_jsonl(path, [_row(0.32, 0.74, 0.0, 0.0, 0.0, ld=500.0)])
    pts = an.frontier(an.read_sweep(path))
    assert pts[0].frac_condensation == pytest.approx(1.0)
    assert pts[0].frac_radiative_wall == pytest.approx(0.0)


def test_write_summary_has_header_and_rows(tmp_path: Path) -> None:
    """The CSV header is exactly the `FrontierPoint` fields, one row per frontier point."""
    path = tmp_path / "sweep.jsonl"
    _write_jsonl(path, [_row(0.16, 0.628, 1.0, 1.0, 0.0), _row(0.32, 0.633, 2.0, 1.0, 0.0)])
    pts = an.frontier(an.read_sweep(path))

    out = tmp_path / "frontier.csv"
    an.write_summary(pts, out)
    with out.open(newline="") as fh:
        reader = list(csv.reader(fh))
    assert reader[0] == [f.name for f in fields(an.FrontierPoint)]
    assert len(reader) == 1 + len(pts)  # header + points
    assert float(reader[1][0]) == 0.16  # first data row, ascending in rho


def test_plot_frontier_writes_files(tmp_path: Path) -> None:
    """The plotting path renders both figures (skipped when matplotlib is absent)."""
    pytest.importorskip("matplotlib")
    path = tmp_path / "sweep.jsonl"
    _write_jsonl(path, [_row(0.16, 0.628, 1.0, 1.0, 0.0), _row(0.32, 0.633, 2.0, 1.0, 1.0)])
    pts = an.frontier(an.read_sweep(path))

    saved = an.plot_frontier(pts, tmp_path)
    assert len(saved) == 2
    for fig_path in saved:
        assert fig_path.exists()
        assert fig_path.stat().st_size > 0
