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


def _grow(
    d_over_d: float,
    l_over_d: float,
    r_foot_over_r: float,
    eta: float,
    peak_local: float = 1.0,
) -> dict[str, float]:
    """A geometry-sweep row at the given case; the two restitution ratios are stand-ins (the reader
    keeps them but the reconciliation keys on `eta_capture`). `peak_local` is the free run's local
    peak pressure, the Rung S focusing factor's numerator."""
    return {
        "d_over_d": d_over_d,
        "l_over_d": l_over_d,
        "r_foot_over_r": r_foot_over_r,
        "mach": 10.0,
        "eta_capture": eta,
        "restitution_free": eta * 1.5,
        "restitution_confined": 1.5,
        "peak_force": 0.6,
        "peak_local_pressure": peak_local,
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


# --- Rung S: survivability frontier (peak facesheet pressure vs P_limit) ---


def test_peak_facesheet_pressure_is_c_stag_rho_v2() -> None:
    """Peak facesheet pressure is the stagnation pressure `c_stag·rho·v²` (design §7): a cold
    coasting cloud feels the ram pressure recompressed at the wall."""
    # c_stag ~ 2.0 at 16 km/s (the measured 1D coefficient), rho = 0.64 -> 330 MPa.
    peak = an.peak_facesheet_pressure(0.64, 16_000.0, 2.0)
    assert peak == pytest.approx(2.0 * 0.64 * 16_000.0**2)
    assert peak == pytest.approx(3.2768e8)


def test_stagnation_coefficient_recovered_from_sweep() -> None:
    """`stagnation_coefficient` backs `c_stag = peak_wall_force / (rho·v²)` out of the 1D sweep
    rows, averaging over the densities at one velocity."""
    v = 16_000.0
    rows = [
        an.SweepRow(
            rho_impact=rho,
            v=v,
            e_eff=0.63,
            peak_wall_force=2.0 * rho * v**2,  # exactly c_stag = 2.0
            loss_radiative_wall=0.0,
            loss_escape_space=0.0,
            loss_conductive=0.0,
            loss_condensation=0.0,
        )
        for rho in (0.16, 0.32, 0.64)
    ]
    assert an.stagnation_coefficient(rows, v) == pytest.approx(2.0)


def test_impact_density_sigma_bridge() -> None:
    """The Σ contract maps a geometry case to a physical impact density:
    `rho = m / (2π·(L/D)·(r_foot/R)³·R³)` (m=25 kg, R=5 m). A disk (small L/D) is dense; a
    cylinder dilute; a tighter footprint is denser (∝ 1/(r_foot/R)³)."""
    rho_disk = an.impact_density(0.3, 0.5)  # L/D=0.3, r_foot/R=0.5
    rho_cyl = an.impact_density(1.0, 0.5)
    assert rho_disk == pytest.approx(0.849, abs=2e-3)
    assert rho_cyl == pytest.approx(0.2546, abs=2e-3)
    # density scales as 1/(L/D): the disk is (1.0/0.3)x the cylinder
    assert rho_disk / rho_cyl == pytest.approx(1.0 / 0.3)
    # tighter footprint is denser as 1/(r_foot/R)^3
    assert an.impact_density(0.3, 0.25) / rho_disk == pytest.approx((0.5 / 0.25) ** 3)


def test_density_ceiling_inverts_the_pressure_law() -> None:
    """`density_ceiling` is the densest cloud under `P_limit`: `rho = P_limit/(c_stag·v²)`."""
    rho_ceil = an.density_ceiling(16_000.0, 2.0, an.P_LIMIT_BASELINE)
    assert rho_ceil == pytest.approx(400.0e6 / (2.0 * 16_000.0**2))
    assert rho_ceil == pytest.approx(0.78125, abs=1e-4)


def test_classify_survivability_disk_fails_cylinder_passes() -> None:
    """At 16 km/s baseline 400 MPa: the disk (rho~0.85 -> ~439 MPa) fails the compressive limit
    but its reflected tension (~66 MPa) clears SiC spall; the cylinder (rho~0.25 -> ~130 MPa)
    clears both. This is the punchline: the f>0.8 disk corner is not survivable at baseline."""
    v, c_stag = 16_000.0, 2.0
    peak_disk = an.peak_facesheet_pressure(an.impact_density(0.3, 0.5), v, c_stag)
    peak_cyl = an.peak_facesheet_pressure(an.impact_density(1.0, 0.5), v, c_stag)
    disk = an.classify_survivability(peak_disk, an.P_LIMIT_BASELINE, an.SIC_SPALL_LO)
    cyl = an.classify_survivability(peak_cyl, an.P_LIMIT_BASELINE, an.SIC_SPALL_LO)
    assert disk.peak_compressive == pytest.approx(peak_disk)
    assert disk.reflected_tensile == pytest.approx(an.REFLECT_FRAC * peak_disk)
    assert not disk.survives_compressive  # ~439 MPa > 400 MPa baseline
    assert disk.survives_spall  # reflected ~66 MPa < SiC spall 0.3 GPa
    assert cyl.survives_compressive and cyl.survives_spall
    # the disk survives once the high-v limit is relaxed to 900 MPa
    assert an.classify_survivability(peak_disk, 900.0e6, an.SIC_SPALL_LO).survives_compressive


def test_survivability_frontier_maps_sigma_and_classifies(tmp_path: Path) -> None:
    """Each geometry case is resolved to physical (rho, peak) via the Σ bridge at a velocity anchor
    and classified. At 16 km/s the disk (L/D=0.3) fails the 400 MPa baseline but clears the relaxed
    900; the cylinder (L/D=1.0) clears both."""
    path = tmp_path / "geom.jsonl"
    _write_jsonl(path, [_grow(0.0, 0.3, 0.5, 0.87), _grow(0.0, 1.0, 0.5, 0.73)])
    pts = an.survivability_frontier(an.read_geometry(path), [(16_000.0, 0.63, 2.0)])
    by_ld = {p.l_over_d: p for p in pts}
    disk, cyl = by_ld[0.3], by_ld[1.0]
    assert disk.rho_impact == pytest.approx(an.impact_density(0.3, 0.5))
    assert disk.peak_compressive == pytest.approx(2.0 * disk.rho_impact * 16_000.0**2)
    assert disk.f == pytest.approx(an.reconcile_f(0.87, 0.63))
    assert not disk.survives_baseline  # ~438 MPa > 400
    assert disk.survives_relaxed  # < 900
    assert cyl.survives_baseline and cyl.survives_relaxed


def test_best_survivable_f_excludes_the_unsurvivable_high_f_corner(tmp_path: Path) -> None:
    """The f-maximizing corner (short disk + tight footprint) is the densest case and fails even the
    relaxed limit, so `best_survivable_f` returns a strictly lower f from an intermediate shape."""
    path = tmp_path / "geom.jsonl"
    _write_jsonl(path, [_grow(0.10, 0.3, 0.3, 1.02), _grow(0.10, 0.6, 0.7, 0.87)])
    pts = an.survivability_frontier(an.read_geometry(path), [(16_000.0, 0.63, 2.0)])
    corner = next(p for p in pts if p.l_over_d == 0.3)
    assert not corner.survives_baseline and not corner.survives_relaxed  # ~2 GPa
    best = an.best_survivable_f(pts)
    assert best == pytest.approx(an.reconcile_f(0.87, 0.63))
    assert best is not None and best < corner.f


def test_survivability_applies_concave_focusing_factor(tmp_path: Path) -> None:
    """A concave plate focuses the local peak above the flat reference at the same (L/D, r_foot/R),
    so its survivability peak is scaled up by the focusing factor (concave local / flat local). The
    flat reference itself has focusing 1.0."""
    path = tmp_path / "geom.jsonl"
    _write_jsonl(
        path,
        [
            _grow(0.0, 0.6, 0.5, 0.83, peak_local=1.0),  # flat reference
            _grow(0.15, 0.6, 0.5, 0.93, peak_local=1.5),  # concave focuses 1.5x
        ],
    )
    pts = an.survivability_frontier(an.read_geometry(path), [(16_000.0, 0.63, 2.0)])
    flat = next(p for p in pts if p.d_over_d == 0.0)
    concave = next(p for p in pts if p.d_over_d == 0.15)
    assert flat.focusing_factor == pytest.approx(1.0)
    assert concave.focusing_factor == pytest.approx(1.5)
    assert concave.rho_impact == pytest.approx(flat.rho_impact)  # same L/D, footprint
    assert concave.peak_compressive == pytest.approx(1.5 * flat.peak_compressive)


def test_plot_survivability_writes_file(tmp_path: Path) -> None:
    """The survivability figure renders for the high-v anchor (skipped without matplotlib)."""
    pytest.importorskip("matplotlib")
    path = tmp_path / "geom.jsonl"
    _write_jsonl(
        path,
        [_grow(0.10, ld, rf, 0.9) for ld in (0.3, 0.6, 1.0) for rf in (0.3, 0.5, 0.7)],
    )
    pts = an.survivability_frontier(an.read_geometry(path), [(16_000.0, 0.63, 2.0)])
    saved = an.plot_survivability(pts, tmp_path)
    assert len(saved) == 1
    assert saved[0].exists()
    assert saved[0].stat().st_size > 0


def test_margin_map_widening_or_lightening_buys_f(tmp_path: Path) -> None:
    """The closed-form margin map: a wider plate (R↑) or smaller pulse (m↓) relaxes the pressure
    ceiling (`peak ∝ m/R³`) and admits a denser, higher-`f` shape that fails at the baseline. Three
    shapes at 16 km/s: a dense corner (~2 GPa, never survives), an intermediate disk (~435 MPa at
    R=5/m=25, just over baseline), and a dilute cylinder (~48 MPa, always survives). At the baseline
    only the cylinder survives; widening to R=5.5 m *or* lightening to 20 kg flips the intermediate
    in, lifting the best survivable `f` from the cylinder's to the intermediate's."""
    path = tmp_path / "geom.jsonl"
    _write_jsonl(
        path,
        [
            _grow(0.0, 0.3, 0.3, 1.00),  # dense corner: ~2 GPa, never survives in the grid
            _grow(0.0, 0.3, 0.5, 0.90),  # intermediate: ~435 MPa at R=5/m=25 (just fails 400)
            _grow(0.0, 1.0, 0.7, 0.73),  # dilute cylinder: ~48 MPa, always survives
        ],
    )
    pts = an.margin_map(
        an.read_geometry(path),
        [(16_000.0, 0.63, 2.0)],
        plate_radii=(5.0, 5.5),
        masses=(25.0, 20.0),
    )
    cell = {(p.plate_radius, p.mass): p for p in pts}
    f_cyl = an.reconcile_f(0.73, 0.63)
    f_mid = an.reconcile_f(0.90, 0.63)

    # headroom is exact: (R/R0)^3 * (m0/m)
    assert cell[(5.0, 25.0)].headroom == pytest.approx(1.0)
    assert cell[(5.5, 25.0)].headroom == pytest.approx((5.5 / 5.0) ** 3)
    assert cell[(5.0, 20.0)].headroom == pytest.approx(25.0 / 20.0)

    # at the pinned baseline only the dilute cylinder survives
    assert cell[(5.0, 25.0)].best_f_baseline == pytest.approx(f_cyl)
    # widening the plate OR lightening the pulse flips the intermediate in -> higher f
    assert cell[(5.5, 25.0)].best_f_baseline == pytest.approx(f_mid)
    assert cell[(5.0, 20.0)].best_f_baseline == pytest.approx(f_mid)
    assert f_mid > f_cyl
    # the dense corner never survives even at the most-relaxed grid cell (still ~1.2 GPa)
    assert cell[(5.5, 20.0)].best_f_baseline == pytest.approx(f_mid)


# --- Rung E: ablating-wall recovery (the tau-bracket + the 16 km/s f-gate call) ---


def _arow(
    v: float,
    scale: float,
    q_star: float,
    rho: float,
    e_rigid: float,
    e_abl: float,
    frac: float = 0.02,
) -> dict[str, float]:
    """An ablating-sweep row (`crates/sweep --ablating`): the rigid floor, the ablating restitution,
    and the recovery, at one `(v, opacity_scale, Q*, rho)` case."""
    return {
        "v": v,
        "rho_impact": rho,
        "opacity_scale": scale,
        "q_star": q_star,
        "kappa_vapor": 200.0,
        "e_eff_rigid": e_rigid,
        "e_eff_ablating": e_abl,
        "recovery": e_abl - e_rigid,
        "ablated_mass": frac * rho,
        "ablated_fraction": frac,
        "loss_radiative_wall": 1.0e4,
        "loss_escape_space": 1.0e3,
        "loss_ablation": 2.0e3,
        "peak_wall_force": 1.0,
    }


def test_read_ablating_parses_jsonl(tmp_path: Path) -> None:
    """`read_ablating` reads one object per line (tolerating a trailing blank) into rows."""
    path = tmp_path / "abl.jsonl"
    path.write_text(json.dumps(_arow(16_000.0, 1.0, 5.0e6, 0.32, 0.63, 0.64)) + "\n\n")
    rows = an.read_ablating(path)
    assert len(rows) == 1
    assert rows[0].v == 16_000.0
    assert rows[0].e_eff_rigid == 0.63
    assert rows[0].e_eff_ablating == 0.64
    assert rows[0].recovery == pytest.approx(0.01)


def test_ablating_points_rho_means(tmp_path: Path) -> None:
    """`ablating_points` collapses the impact-density axis: the rho-mean rigid floor, ablating
    `e_eff`, and recovery at each `(v, scale, Q*)`."""
    path = tmp_path / "abl.jsonl"
    _write_jsonl(
        path,
        [
            _arow(16_000.0, 1.0, 5.0e6, 0.16, 0.62, 0.63),
            _arow(16_000.0, 1.0, 5.0e6, 0.48, 0.64, 0.67),  # same case, second density
            _arow(16_000.0, 0.1, 5.0e6, 0.16, 0.62, 0.65),  # different scale
        ],
    )
    pts = an.ablating_points(an.read_ablating(path))
    assert len(pts) == 2  # two (scale) cases, each rho-meaned
    scale1 = next(p for p in pts if p.opacity_scale == 1.0)
    assert scale1.e_eff_rigid == pytest.approx(0.63)  # mean(0.62, 0.64)
    assert scale1.e_eff_ablating == pytest.approx(0.65)  # mean(0.63, 0.67)
    assert scale1.recovery == pytest.approx(0.02)


def test_best_f_at_lifts_with_recovery(tmp_path: Path) -> None:
    """`best_f_at` resolves the geometry cases to survivability at the `(v, e_eff, c_stag)` anchor
    and returns the survivable max `f` — which rises when the ablating recovery lifts `e_eff`, while
    the unsurvivable dense corner stays excluded."""
    path = tmp_path / "geom.jsonl"
    # The dense short-disk corner fails survivability; the elongated, wider-footprint case survives.
    _write_jsonl(path, [_grow(0.10, 0.3, 0.3, 1.02), _grow(0.10, 0.6, 0.7, 0.87)])
    geo = an.read_geometry(path)
    f_floor = an.best_f_at(geo, 16_000.0, 0.63, 2.0)
    f_recovered = an.best_f_at(geo, 16_000.0, 0.70, 2.0)
    assert f_floor == pytest.approx(an.reconcile_f(0.87, 0.63))  # survivable case, not the corner
    assert f_recovered == pytest.approx(an.reconcile_f(0.87, 0.70))
    assert f_recovered is not None and f_floor is not None and f_recovered > f_floor


def test_write_ablating_summary_has_header_and_rows(tmp_path: Path) -> None:
    """The ablating summary CSV has the `AblatingPoint` header and one row per rho-meaned case."""
    path = tmp_path / "abl.jsonl"
    _write_jsonl(
        path,
        [_arow(16_000.0, s, 5.0e6, 0.32, 0.63, 0.64) for s in (0.1, 1.0)],
    )
    pts = an.ablating_points(an.read_ablating(path))
    out = tmp_path / "frontier_ablating.csv"
    an.write_ablating_summary(pts, out)
    with out.open() as fh:
        table = list(csv.reader(fh))
    assert table[0] == [f.name for f in fields(an.AblatingPoint)]
    assert len(table) == 1 + len(pts)


def test_plot_ablating_writes_file(tmp_path: Path) -> None:
    """The ablating recovery figure renders (skipped without matplotlib)."""
    pytest.importorskip("matplotlib")
    path = tmp_path / "abl.jsonl"
    _write_jsonl(
        path,
        [
            _arow(v, s, q, 0.32, 0.62, 0.62 + 0.01 * (1.0 - s))
            for v in (11_000.0, 16_000.0)
            for s in (0.1, 1.0)
            for q in (2.0e6, 5.0e6)
        ],
    )
    pts = an.ablating_points(an.read_ablating(path))
    saved = an.plot_ablating(pts, tmp_path)
    assert len(saved) == 1
    assert saved[0].exists() and saved[0].stat().st_size > 0
