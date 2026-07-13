"""Tests for the ADR-0027 closed-form whole-plate structural bound (`structure.py`, design §12.1).

The three checks (rigid-during-pulse / f-validity, areal-impulse membrane, SiC-Ti spall) are pinned
with synthetic sweep + design-point inputs (stdlib only), so the closed-form arithmetic and verdict
logic are exercised without the expensive real sweep."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from puffsat import structure
from puffsat.heavyplate import (
    PLATE_RADIUS_M,
    HeavyPlatePoint,
    HeavyPlateRow,
    read_heavyplate_sweep,
)


def _row(v: float, rho: float, *, peak: float, impulse: float) -> HeavyPlateRow:
    return HeavyPlateRow(
        v=v,
        rho_impact=rho,
        length=10.0,
        opacity_scale=1.0,
        e_eff=0.7,
        peak_wall_pressure=peak,
        incident_momentum=impulse / 1.7,
        wall_impulse=impulse,
        loss_radiative_wall=1.0e6,
        loss_escape_space=1.0e5,
    )


def _design(peak: float, *, d_over_d: float, r_foot_over_r: float) -> HeavyPlatePoint:
    return HeavyPlatePoint(
        v=20_000.0,
        d_over_d=d_over_d,
        l_over_d=0.6,
        r_foot_over_r=r_foot_over_r,
        rho_impact=0.08,
        e_eff=0.7,
        eta_capture=0.9,
        focusing_factor=1.0,
        plate_mass_t=35.0,
        peak_compressive=peak,
        f=0.765,
        survives_baseline=True,
        survives_relaxed=True,
        within_mass_budget=True,
    )


def test_first_mode_period_is_hundreds_of_ms() -> None:
    """The candidate Ti-sandwich fundamental flexural period is ~0.3 s at R = 15 m, 45 kg/m^2 —
    orders of magnitude above the ~ms bounce, so the rigidity gate has comfortable margin."""
    t1 = structure.first_mode_period()
    assert t1 == pytest.approx(0.295, abs=0.01)
    # Heavier plate -> longer period (softer omega1), so check 1 stays conservative.
    assert structure.first_mode_period(90.0) > t1


def test_pulse_duration_is_impulse_over_peak() -> None:
    """tau_pulse = 2 J / p_peak; a ~2e4 Pa.s impulse at a ~2e7 Pa peak gives a ~ms pulse."""
    assert structure.pulse_duration(2.0e4, 2.0e7) == pytest.approx(2.0e-3)
    with pytest.raises(ValueError, match="non-positive"):
        structure.pulse_duration(2.0e4, 0.0)


def test_structure_point_membrane_and_gates() -> None:
    """The membrane arithmetic and the three gates at a known design point: p_eff = peak*(rf/R)^2,
    N = p_eff*Rc/2 with Rc = R/(4 d/D), implied back-face = N/allowable; the rigidity + spall gates
    pass for a survivable ~ms pulse, and the implied plate mass is reported for the budget check."""
    rows = [
        _row(20_000.0, 0.04, peak=2.0e7, impulse=2.0e4),
        _row(20_000.0, 0.16, peak=2.0e7, impulse=2.0e4),
    ]
    design = _design(3.0e7, d_over_d=0.10, r_foot_over_r=0.5)
    sp = structure.structure_point(design, rows)

    # Check 1 — rigid-during-pulse: tau = 2*2e4/2e7 = 2 ms; T1 ~ 0.295 s; ratio ~ 148 >> 10.
    assert sp.pulse_duration == pytest.approx(2.0e-3)
    assert sp.rigidity_ratio == pytest.approx(structure.first_mode_period() / 2.0e-3)
    assert sp.rigid_ok

    # Check 2 — membrane: p_eff = 3e7*0.25 = 7.5e6; Rc = 15/(4*0.1) = 37.5; N = 7.5e6*37.5/2.
    p_eff = 3.0e7 * 0.25
    r_curv = PLATE_RADIUS_M / (4.0 * 0.10)
    n_expected = p_eff * r_curv / 2.0
    assert sp.membrane_tension == pytest.approx(n_expected)
    assert sp.back_thickness_req == pytest.approx(n_expected / structure.SIGMA_FIBER_WORKING)
    implied_areal = structure.BASE_AREAL_NO_BACK + sp.back_thickness_req * structure.RHO_FIBER
    assert sp.implied_plate_mass_t == pytest.approx(
        implied_areal * math.pi * PLATE_RADIUS_M**2 / 1000.0
    )

    # Check 3 — spall: reflected 0.15*3e7 = 4.5 MPa << 0.3 GPa SiC spall.
    assert sp.reflected_tensile == pytest.approx(0.15 * 3.0e7)
    assert sp.spall_ok
    assert sp.verdict_ok == (sp.rigid_ok and sp.mass_ok and sp.spall_ok)


def test_flat_design_uses_minimum_builtin_dish() -> None:
    """A flat design point (d/D = 0) uses the minimum built-in dish for the membrane radius of
    curvature (no singular Rc), giving a finite membrane tension."""
    rows = [
        _row(20_000.0, 0.04, peak=1.0e7, impulse=1.0e4),
        _row(20_000.0, 0.16, peak=1.0e7, impulse=1.0e4),
    ]
    design = _design(1.0e7, d_over_d=0.0, r_foot_over_r=0.5)
    sp = structure.structure_point(design, rows)
    r_curv = PLATE_RADIUS_M / (4.0 * structure.D_OVER_D_MIN)
    assert sp.membrane_tension == pytest.approx(1.0e7 * 0.25 * r_curv / 2.0)
    assert math.isfinite(sp.membrane_tension)


def test_structure_frontier_end_to_end(tmp_path: Path) -> None:
    """End-to-end on the committed artifacts: the frontier evaluates + writes the CSV (header =
    `StructurePoint` fields). Skipped if the gitignored bulk sweep artifacts are absent."""
    sweep = Path("data/results/sweep_heavyplate.jsonl")
    geo = Path("data/results/sweep_geometry_m40.jsonl")
    for pth in (sweep, geo):
        if not pth.exists():
            pytest.skip(f"missing sweep artifact {pth} (run make sweep-heavyplate)")

    from puffsat.heavyplate import _load_geometry, heavyplate_frontier, sweep_velocities

    rows = read_heavyplate_sweep(sweep)
    frontier = heavyplate_frontier(_load_geometry(geo), rows, sweep_velocities(rows))
    pts = structure.structure_frontier(frontier, rows)
    out = tmp_path / "frontier_structure_heavyplate.csv"
    structure.write_summary(pts, out)
    assert out.exists()
    header = out.read_text().splitlines()[0].split(",")
    assert {"rigidity_ratio", "rigid_ok", "implied_plate_mass_t", "spall_ok", "verdict_ok"} <= set(
        header
    )
    # The rigid-during-pulse / f-validity gate should clear at every surviving anchor.
    assert all(p.rigid_ok for p in pts)
