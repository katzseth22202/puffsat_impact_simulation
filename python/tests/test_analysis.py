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


def _vrow(v: float, rho: float, e_eff: float) -> dict[str, float]:
    """A transitional-sweep row at a chosen impact speed `v` (lossless; the EOS/rad files differ
    only in `e_eff`)."""
    d = _row(rho, e_eff, 0.0, 0.0, 0.0)
    d["v"] = v
    return d


def _tp(v: float, e_eos: float, e_rad: float = 0.0) -> an.TransitionalPoint:
    """A `TransitionalPoint` with the rho-spread collapsed onto `e_eos` (enough for the dip/plot
    tests, which key on the EOS-only value)."""
    return an.TransitionalPoint(
        v=v,
        e_eff_eos=e_eos,
        e_eff_rad=e_rad,
        rad_band=e_eos - e_rad,
        e_eff_eos_min=e_eos,
        e_eff_eos_max=e_eos,
    )


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


def test_transitional_frontier_sorts_and_means(tmp_path: Path) -> None:
    """The `e_eff(v)` frontier is ascending in v; each point is the rho-mean of its velocity, and
    `rad_band = e_eff_eos - e_eff_rad`. The rho-spread is bracketed by min/max."""
    eos = tmp_path / "eos.jsonl"
    rad = tmp_path / "rad.jsonl"
    _write_jsonl(
        eos,
        [
            _vrow(11_000, 0.16, 0.56),  # out of order on purpose
            _vrow(11_000, 0.64, 0.58),  # 11 km/s rho-mean 0.57
            _vrow(5_000, 0.16, 0.68),
            _vrow(5_000, 0.64, 0.70),  # 5 km/s rho-mean 0.69
        ],
    )
    _write_jsonl(
        rad,
        [
            _vrow(11_000, 0.16, 0.55),
            _vrow(11_000, 0.64, 0.57),  # 11 km/s rho-mean 0.56
            _vrow(5_000, 0.16, 0.67),
            _vrow(5_000, 0.64, 0.69),  # 5 km/s rho-mean 0.68
        ],
    )
    pts = an.transitional_frontier(an.read_sweep(eos), an.read_sweep(rad))

    assert [p.v for p in pts] == [5_000.0, 11_000.0]  # ascending in v
    assert pts[0].e_eff_eos == pytest.approx(0.69)
    assert pts[1].e_eff_eos == pytest.approx(0.57)
    assert pts[1].e_eff_rad == pytest.approx(0.56)
    assert pts[1].rad_band == pytest.approx(0.57 - 0.56)
    assert pts[1].e_eff_eos_min == pytest.approx(0.56)
    assert pts[1].e_eff_eos_max == pytest.approx(0.58)


def test_locate_dip_finds_interior_minimum() -> None:
    """A planted interior minimum below both endpoints is found and returned."""
    pts = [
        _tp(5_000, 0.68),
        _tp(8_000, 0.59),
        _tp(11_000, 0.567),
        _tp(13_000, 0.58),
        _tp(16_000, 0.64),
    ]
    dip = an.locate_dip(pts)
    assert dip is not None
    assert dip.v == 11_000.0
    assert dip.e_eff_eos == pytest.approx(0.567)


def test_locate_dip_none_when_monotonic() -> None:
    """A monotone curve (minimum at an endpoint) reports no interior dip — the floor is at an edge,
    so any transitional dip would be purely radiative."""
    pts = [_tp(5_000, 0.50), _tp(8_000, 0.55), _tp(11_000, 0.60), _tp(16_000, 0.64)]
    assert an.locate_dip(pts) is None


def test_plot_transitional_writes_file(tmp_path: Path) -> None:
    """The `e_eff(v)` overlay renders one figure with the dip annotated (skipped w/o matplotlib)."""
    pytest.importorskip("matplotlib")
    pts = [_tp(5_000, 0.68, 0.67), _tp(11_000, 0.567, 0.56), _tp(16_000, 0.64, 0.63)]
    saved = an.plot_transitional(pts, an.locate_dip(pts), tmp_path)
    assert len(saved) == 1
    assert saved[0].exists()
    assert saved[0].stat().st_size > 0


def _grow(d_over_d: float, l_over_d: float, r_foot_over_r: float, eta: float) -> dict[str, float]:
    """A geometry-sweep row at the given case; the two restitution ratios are stand-ins (the reader
    keeps them but the reconciliation keys on `eta_capture`)."""
    return {
        "d_over_d": d_over_d,
        "l_over_d": l_over_d,
        "r_foot_over_r": r_foot_over_r,
        "mach": 10.0,
        "eta_capture": eta,
        "restitution_free": eta * 1.5,
        "restitution_confined": 1.5,
        "peak_force": 0.6,
    }


def test_reconcile_f_formula() -> None:
    """`f = eta_capture·(1 + e_eff)/2` (ADR-0003)."""
    assert an.reconcile_f(0.9, 0.57) == pytest.approx(0.9 * 1.57 / 2.0)
    assert an.reconcile_f(1.0, 1.0) == pytest.approx(1.0)  # elastic + perfect collimation


def test_geometry_frontier_reconciles_f_and_sigma(tmp_path: Path) -> None:
    """The geometry frontier reconciles `f` at both `e_eff` anchors, reports the `Sigma` contract
    (`= 2·L/D`, footprint-independent), and sorts by `(mach, L/D, r_foot/R, d/D)`."""
    path = tmp_path / "geom.jsonl"
    _write_jsonl(
        path,
        [
            _grow(0.15, 0.6, 0.5, 0.93),  # out of order on purpose
            _grow(0.0, 0.6, 0.5, 0.83),
        ],
    )
    pts = an.geometry_frontier(an.read_geometry(path))
    assert [p.d_over_d for p in pts] == [0.0, 0.15]  # sorted ascending in d/D within the slice
    flat, concave = pts[0], pts[1]
    assert flat.f_dip == pytest.approx(0.83 * (1.0 + an.EEFF_DIP) / 2.0)
    assert concave.f_highv == pytest.approx(0.93 * (1.0 + an.EEFF_HIGHV) / 2.0)
    assert concave.f_dip > flat.f_dip  # the concave plate lifts f over the flat floor
    assert flat.sigma_over_rho == pytest.approx(2.0 * 0.6)  # Sigma set by L/D, not footprint


def test_plot_geometry_writes_file(tmp_path: Path) -> None:
    """The geometry figure renders for a representative slice (skipped without matplotlib)."""
    pytest.importorskip("matplotlib")
    path = tmp_path / "geom.jsonl"
    _write_jsonl(
        path,
        [
            _grow(dd, 0.6, rf, 0.8 + dd + 0.05 * rf)
            for dd in (0.0, 0.10, 0.15)
            for rf in (0.3, 0.5, 0.7)
        ],
    )
    pts = an.geometry_frontier(an.read_geometry(path))
    saved = an.plot_geometry(pts, tmp_path)
    assert len(saved) == 1
    assert saved[0].exists()
    assert saved[0].stat().st_size > 0
