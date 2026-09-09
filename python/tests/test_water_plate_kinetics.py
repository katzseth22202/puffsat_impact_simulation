"""Detailed balance, event convention and analytic linear relaxation checks."""

from dataclasses import replace

import numpy as np
import pytest

from puffsat import eos_water as ew
from puffsat.water_plate.flow import FlowCell, Snapshot
from puffsat.water_plate.kinetics import (
    ATOM_COUNTS,
    GROUND_ENERGIES,
    REACTIONS,
    log_equilibrium_constant,
    net_production,
    neutral_densities,
    relaxation,
)
from puffsat.water_plate.kinetics_audit import (
    EDGES,
    Row,
    select_snapshots,
    store_change_clock,
    summarize,
    validate_parent,
    weighted_quantile,
)


def test_molar_cgs_to_molecular_si_and_water_collider_efficiency() -> None:
    water = next(r for r in REACTIONS if r.source_number == 15)
    expected = 2.2e22 / ew.N_A**2 * 1e-12 / 3000**2
    assert water.coefficient(3000) == pytest.approx(expected, rel=1e-14)
    n = np.ones(6) * 1e20
    assert water.forward_flux(3000, n) == pytest.approx(expected * 8.38e60)


def test_every_retained_reaction_conserves_both_elements() -> None:
    for reaction in REACTIONS:
        assert ATOM_COUNTS @ reaction.stoichiometry == pytest.approx(np.zeros(2), abs=0)


@pytest.mark.parametrize(
    "rho,temp",
    [
        (0.1, 3000.0),
        (1.0, 8000.0),
        (100.0, 16000.0),
        (279.0, 9000.0),
        (71.0, 6760.0),
        (11.3, 5100.0),
    ],
)
def test_detailed_balance_and_energy_weighted_relaxation(rho: float, temp: float) -> None:
    n = neutral_densities(ew.composition(rho, temp))
    for reaction in REACTIONS:
        assert reaction.stoichiometry @ np.log(n) == pytest.approx(
            log_equilibrium_constant(reaction, temp), abs=1e-9
        )
    result = relaxation(temp, n)
    assert result.energy_time_s > 0
    assert len(result.mode_rates_s) == 4
    assert sum(result.energy_weights) == pytest.approx(1)
    assert min(result.energy_weights) >= 0
    slow = relaxation(temp, n, rate_scale=0.01)
    assert slow.energy_time_s == pytest.approx(100 * result.energy_time_s, rel=1e-7)
    # Conserved elemental energy zeros must not change the reactive energy projection.
    shifted = GROUND_ENERGIES + ATOM_COUNTS.T @ np.array([ew.D_AT, -2 * ew.D_AT])
    gauge = relaxation(temp, n, shifted)
    assert gauge.energy_time_s == pytest.approx(result.energy_time_s, rel=1e-8)


def test_reversible_jacobian_matches_detailed_balance_outer_products() -> None:
    temp = 5000.0
    n = neutral_densities(ew.composition(1.0, temp))
    analytic = np.zeros((6, 6))
    for reaction in REACTIONS:
        nu = reaction.stoichiometry
        analytic -= reaction.forward_flux(temp, n) * np.outer(nu, nu / n)
    measured = np.zeros((6, 6))
    for j in range(6):
        delta = np.zeros(6)
        delta[j] = 1e-5 * n[j]
        measured[:, j] = (net_production(temp, n + delta) - net_production(temp, n - delta)) / (
            2 * delta[j]
        )
    assert np.linalg.norm(measured - analytic) / np.linalg.norm(analytic) < 1e-7


def test_single_water_association_rate_includes_reverse_and_both_radicals() -> None:
    temp = 5000.0
    n = neutral_densities(ew.composition(1.0, temp))
    water = next(r for r in REACTIONS if r.source_number == 15)
    rate = water.forward_flux(temp, n) * np.sum(water.stoichiometry**2 / n)
    assert water.colliders is not None
    effective_k = water.coefficient(temp) * float(np.dot(water.colliders, n))
    assert rate == pytest.approx(effective_k * (n[0] + n[3] + n[0] * n[3] / n[5]))


def test_nonequilibrium_input_cannot_masquerade_as_a_linear_equilibrium_clock() -> None:
    n = neutral_densities(ew.composition(1.0, 5000.0))
    n[5] *= 2
    with pytest.raises(ValueError, match="not an equilibrium"):
        relaxation(5000, n)


def test_unresolved_release_stays_in_the_rate_screen_denominator() -> None:
    rows: list[Row] = [
        {
            "weight_j": 6.0,
            "temperature_k": 5000.0,
            "pressure_pa": 1e9,
            "status": "resolved_local_network",
            "da_energy": 1e5,
            "da_tracking": 1e4,
        },
        {
            "weight_j": 4.0,
            "temperature_k": 2000.0,
            "pressure_pa": 1e6,
            "status": "unresolved: test",
        },
    ]
    result = summarize(rows)
    assert result["network_da_gt10_fraction"] == 0.6
    assert result["unresolved_network_fraction"] == 0.4
    assert result["above_reference_thermo_3500k_fraction"] == 0.6
    assert result["kinetics_validated"] is False
    assert weighted_quantile(rows, "temperature_k", 0.5) == 5000
    assert summarize([], 100)["network_da_gt10_fraction"] == ""
    assert summarize(rows, 1e8)["status"] == "negligible_molecular_return"


def test_thinning_preserves_every_window_boundary() -> None:
    times = sorted([*EDGES, 150e-6, 180e-6, 230e-6])
    snapshots = [Snapshot("test", t, 0, 0, []) for t in times]
    selected = select_snapshots(snapshots, 3)
    assert set(EDGES) <= {s.time_s for s in selected}
    with pytest.raises(ValueError, match="exact checkpoints"):
        select_snapshots(snapshots[1:], 1)


def test_sampled_store_change_clock_and_invalid_inputs() -> None:
    assert store_change_clock(10, 6, 2) == 4
    assert store_change_clock(10, 0, 2) == 1
    with pytest.raises(ValueError):
        store_change_clock(6, 10, 2)


def test_audit_requires_conserved_fixed_mass_parent() -> None:
    cell = FlowCell(0.5, 1, 1, 1, 5000, 1e7, 1e5, 0, 1)
    snapshots = [Snapshot("test", t, 1, 0, [cell]) for t in (100e-6, 250e-6)]
    summary: Row = {
        "status": "completed_time_window",
        "input_energy_j": 1e8,
        "max_energy_error_j": 100,
        "momentum_error_ns": 0,
        "wall_impulse_ns": 1,
        "time_s": 0.001,
    }
    validate_parent(summary, snapshots)
    for key, value in (
        ("max_energy_error_j", 1e6),
        ("momentum_error_ns", 0.1),
        ("time_s", 0.0002),
        ("input_energy_j", float("nan")),
    ):
        with pytest.raises(ValueError):
            validate_parent({**summary, key: value}, snapshots)
    changed = replace(snapshots[1], cells=[replace(cell, mass_kg=2)])
    with pytest.raises(ValueError, match="fixed parcel"):
        validate_parent(summary, [snapshots[0], changed])
