"""Tests for the heavy-plate 16-28 km/s f(v) + survivability frontier (`heavyplate.py`, §12.1).

The pinned-plate / swept-velocity frontier and the freeze-bracket translation are exercised with
synthetic sweep + geometry rows (stdlib + the analysis reconciliation only — no matplotlib/CoolProp)
so the physics wiring is pinned without the expensive real sweep."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from puffsat import heavyplate
from puffsat.analysis import GeoRow


def _row(
    v: float,
    rho: float,
    e_eff: float,
    *,
    length: float = 10.0,
    scale: float = 1.0,
    peak: float = 2.0e7,
    impulse: float = 2.0e4,
) -> heavyplate.HeavyPlateRow:
    return heavyplate.HeavyPlateRow(
        v=v,
        rho_impact=rho,
        length=length,
        opacity_scale=scale,
        e_eff=e_eff,
        peak_wall_pressure=peak,
        incident_momentum=impulse / (1.0 + e_eff),
        wall_impulse=impulse,
        loss_radiative_wall=1.0e6,
        loss_escape_space=1.0e5,
    )


def _geo(d_over_d: float, l_over_d: float, r_foot_over_r: float, eta: float) -> GeoRow:
    # A flat reference peak of 1.0; concave rows carry a mild focusing concentration.
    local = 1.0 if d_over_d == 0.0 else 1.05
    return GeoRow(
        d_over_d=d_over_d,
        l_over_d=l_over_d,
        r_foot_over_r=r_foot_over_r,
        mach=40.0,
        eta_capture=eta,
        restitution_free=0.9,
        restitution_confined=0.9,
        peak_force=1.0,
        peak_local_pressure=local,
    )


def test_read_heavyplate_sweep_round_trips(tmp_path: Path) -> None:
    """The reader parses the `--heavyplate` JSONL and tolerates a trailing blank line."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(
        json.dumps(
            {
                "v": 22_000.0,
                "rho_impact": 0.08,
                "length": 10.0,
                "opacity_scale": 1.0,
                "e_eff": 0.74,
                "peak_wall_pressure": 4.6e7,
                "incident_momentum": 1.76e4,
                "wall_impulse": 3.06e4,
                "loss_radiative_wall": 1.0e6,
                "loss_escape_space": 1.0e5,
            }
        )
        + "\n\n"
    )
    rows = heavyplate.read_heavyplate_sweep(p)
    assert len(rows) == 1
    assert rows[0].v == 22_000.0
    assert rows[0].wall_impulse == 3.06e4


def test_e_eff_interpolator_is_per_velocity() -> None:
    """`e_eff_interpolator_at_v` selects the requested velocity slice and interpolates log-rho; a
    different velocity's rows do not leak in."""
    rows = [
        _row(16_000.0, 0.04, 0.60),
        _row(16_000.0, 0.16, 0.70),
        _row(28_000.0, 0.04, 0.50),
        _row(28_000.0, 0.16, 0.55),
    ]
    e16 = heavyplate.e_eff_interpolator_at_v(rows, 16_000.0)
    # log(0.08) is the midpoint of log(0.04) and log(0.16) (ratio 2 vs 4), so e_eff is the mean.
    assert e16(0.08) == pytest.approx(0.65)
    # Clamped at the ends.
    assert e16(0.01) == pytest.approx(0.60)
    assert e16(1.0) == pytest.approx(0.70)
    # The 28 km/s slice is distinct.
    assert heavyplate.e_eff_interpolator_at_v(rows, 28_000.0)(0.08) == pytest.approx(0.525)


def test_frontier_reconciles_f_and_flags_mass_budget() -> None:
    """Each (velocity x shape) resolves to `f = eta*(1+e_eff)/2` at the Sigma density, and the
    30 m plate at 45 kg/m^2 (~32 t flat) is within the <= 40 t budget."""
    rows = [_row(16_000.0, 0.04, 0.60), _row(16_000.0, 0.16, 0.70)]
    geo = [_geo(0.0, 0.6, 0.5, 0.80), _geo(0.10, 0.6, 0.5, 0.85)]
    pts = heavyplate.heavyplate_frontier(geo, rows, [16_000.0])
    assert len(pts) == 2
    for p in pts:
        assert p.f == pytest.approx(p.eta_capture * (1.0 + p.e_eff) / 2.0)
        assert p.within_mass_budget  # flat ~32 t, shallow concave ~35 t < 40 t
        assert p.plate_mass_t < 40.0
    # The dilute Sigma density (100 kg over a 15 m footprint) keeps the facesheet peak survivable.
    assert all(p.survives_baseline for p in pts)


def test_best_at_v_picks_highest_surviving_f() -> None:
    """`best_at_v` returns the highest-f surviving case, filterable by curvature."""
    rows = [_row(16_000.0, 0.04, 0.60), _row(16_000.0, 0.16, 0.70)]
    geo = [_geo(0.0, 0.6, 0.5, 0.80), _geo(0.10, 0.6, 0.5, 0.90)]
    pts = heavyplate.heavyplate_frontier(geo, rows, [16_000.0])
    best = heavyplate.best_at_v(pts, 16_000.0)
    assert best is not None and best.d_over_d == 0.10  # higher eta wins
    flat = heavyplate.best_at_v(pts, 16_000.0, concave=False)
    assert flat is not None and flat.d_over_d == 0.0


def test_length_and_opacity_sensitivity_spreads() -> None:
    """The L-sensitivity and opacity-tau diagnostics report the e_eff spread across the spot rows;
    a flat (L,kappa) response gives ~0 (the design-§12.1 tau >> 1 signature)."""
    rows = [
        # headline slice at L = 10, kappa = 1 for the anchor.
        _row(16_000.0, 0.08, 0.66, length=10.0, scale=1.0),
        # L-spot rows: e_eff essentially flat in L.
        _row(16_000.0, 0.08, 0.66, length=6.0, scale=1.0),
        _row(16_000.0, 0.08, 0.665, length=14.0, scale=1.0),
        # tau-check rows at 28 km/s: essentially flat in opacity scale.
        _row(28_000.0, 0.08, 0.60, length=10.0, scale=1.0),
        _row(28_000.0, 0.08, 0.60, length=10.0, scale=0.1),
        _row(28_000.0, 0.08, 0.601, length=10.0, scale=10.0),
    ]
    spreads = heavyplate.length_sensitivity(rows, 0.08)
    assert spreads[16_000.0] == pytest.approx(0.005, abs=1e-9)
    assert heavyplate.opacity_sensitivity(rows, 0.08) == pytest.approx(0.001, abs=1e-9)


def test_frozen_bracket_subtracts_eos_delta_from_coupled_headline() -> None:
    """The freeze bracket applies the EOS-only delta to the coupled headline e_eff(v, rho): f_eq
    reconciles from the coupled e_eff; the sudden-freeze f subtracts the (positive) delta, so it
    sits strictly below the headline; the pure-H2O bound (smaller delta) sits between."""
    sweep = [_row(22_000.0, 0.04, 0.66), _row(22_000.0, 0.16, 0.72)]
    frozen = [heavyplate.FrozenHeavyRow(22_000.0, 0.08, 0.70, 0.52, 0.646, 1e-3)]
    eta_design = {22_000.0: 0.9}
    pts = heavyplate.frozen_heavyplate_bracket(frozen, sweep, eta_design)
    assert len(pts) == 1
    p = pts[0]
    # coupled e_eff(0.08) is the mean of 0.66 and 0.72 = 0.69 (log-midpoint), not the EOS 0.70.
    assert p.e_eff_coupled == pytest.approx(0.69)
    assert p.f_eq == pytest.approx(0.9 * (1.0 + 0.69) / 2.0)
    assert p.delta_frozen == pytest.approx(0.18)
    assert p.f_frozen_rebound == pytest.approx(0.9 * (1.0 + (0.69 - 0.18)) / 2.0)
    assert p.f_frozen_rebound < p.f_frozen_all < p.f_eq


def test_run_headline_end_to_end(tmp_path: Path) -> None:
    """End-to-end on the committed sweep artifacts: the headline run writes the frontier CSV (header
    = `HeavyPlatePoint` fields). Skipped if the gitignored bulk sweep artifacts are absent."""
    sweep = Path("data/results/sweep_heavyplate.jsonl")
    geo = Path("data/results/sweep_geometry_m40.jsonl")
    for pth in (sweep, geo):
        if not pth.exists():
            pytest.skip(f"missing sweep artifact {pth} (run make sweep-heavyplate)")

    rows = heavyplate.read_heavyplate_sweep(sweep)
    geo_rows = heavyplate._load_geometry(geo)
    pts = heavyplate.heavyplate_frontier(geo_rows, rows, heavyplate.sweep_velocities(rows))
    out = tmp_path / "frontier_heavyplate.csv"
    heavyplate.write_summary(pts, out)
    assert out.exists()
    header = out.read_text().splitlines()[0].split(",")
    assert header[0] == "v"
    assert {"f", "peak_compressive", "survives_baseline", "within_mass_budget"} <= set(header)
