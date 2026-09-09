"""Pressure-table interpolation and a single-channel local relaxation bound."""

import math
from itertools import pairwise

import numpy as np
import pytest

from puffsat import eos_water as ew
from puffsat.water_plate.kinetics import (
    ATOM_COUNTS,
    GROUND_ENERGIES,
    REACTIONS,
    Vec,
    equilibrium_log_reference,
    neutral_densities,
    relaxation,
)
from puffsat.water_plate.kinetics_audit import Row
from puffsat.water_plate.pressure_kinetics import (
    ATM_PA,
    PLOG_NODES,
    effective_pressure,
    water_coefficient,
    water_efficiency,
)
from puffsat.water_plate.pressure_kinetics_audit import summarize


def test_pressure_nodes_and_logarithmic_interpolation() -> None:
    temp = 6000.0  # Deliberate temperature extrapolation, not validation.
    for pressure_atm, a, b, ea in PLOG_NODES:
        expected = a * temp**b * math.exp(-ea * 4.184 / (ew.N_A * ew.K_B * temp))
        assert water_coefficient(temp, pressure_atm * ATM_PA) == pytest.approx(
            expected / (1000 * ew.N_A), rel=1e-13
        )
    for lo, hi in pairwise(PLOG_NODES):
        p1, p2 = lo[0] * ATM_PA, hi[0] * ATM_PA
        midpoint = water_coefficient(temp, math.sqrt(p1 * p2))
        assert midpoint == pytest.approx(
            math.sqrt(water_coefficient(temp, p1) * water_coefficient(temp, p2)), rel=1e-13
        )


@pytest.mark.parametrize("pressure", [0, -1, 1e-5 * ATM_PA, 1.001e4 * ATM_PA, math.inf])
def test_pressure_table_never_silently_clamps_or_extrapolates(pressure: float) -> None:
    with pytest.raises(ValueError):
        water_coefficient(6000, pressure)


def test_neutral_number_density_sets_effective_pressure_not_dense_flow_pressure() -> None:
    n = np.arange(1, 7, dtype=float) * 1e23
    temp = 6000.0
    expected = ew.K_B * temp * (sum(n[:-1]) + water_efficiency(temp) * n[-1])
    assert effective_pressure(temp, n) == pytest.approx(expected)
    assert water_efficiency(1000) == pytest.approx(8.55 / 1.62, rel=0.002)


@pytest.mark.parametrize("rho,temp", [(1, 5000), (71, 6760), (279, 9000), (11.3, 5100)])
def test_removing_water_association_is_the_slow_endpoint_for_this_network(
    rho: float, temp: float
) -> None:
    n = neutral_densities(ew.composition(rho, temp))
    times = [
        relaxation(temp, n, water_association_scale=s).energy_time_s
        for s in (0, 1e-6, 0.01, 1, 100)
    ]
    assert all(a >= b * (1 - 1e-9) for a, b in pairwise(times))
    assert np.isfinite(times).all()
    assert relaxation(temp, n).energy_time_s == times[3]


@pytest.mark.parametrize("scale", [-1, math.nan, math.inf])
def test_invalid_association_scales_rejected(scale: float) -> None:
    n = neutral_densities(ew.composition(1, 5000))
    with pytest.raises(ValueError, match="water association scale"):
        relaxation(5000, n, water_association_scale=scale)


def test_pressure_variant_missing_rates_remain_in_denominator() -> None:
    rows: list[Row] = [
        {
            "weight_j": 6.0,
            "temperature_k": 6000.0,
            "pressure_pa": 1e8,
            "status": "resolved_local_network",
            "da_tracking": 1e5,
            "source_pressure_tabulated": True,
            "source_plog_time_s": 1e-9,
            "source_plog_time_over_baseline": 10,
        },
        {
            "weight_j": 4.0,
            "temperature_k": 6000.0,
            "pressure_pa": 1e10,
            "status": "resolved_local_network",
            "da_tracking": 1e6,
            "source_pressure_tabulated": False,
        },
    ]
    result = summarize(rows)
    assert result["source_plog_da_tracking_p50"] == 1e4
    assert result["source_plog_da_gt10_fraction"] == 0.6
    assert result["source_plog_no_rate_fraction"] == 0.4
    assert result["source_pressure_tabulated_fraction"] == 0.6
    assert result["source_above_collider_data_2000k_fraction"] == 1
    assert result["kinetics_validated"] is False


@pytest.mark.parametrize(
    "rho,temp,expected",
    [
        (0.01, 2000.0, 2.71122741446056e-20),
        (1.0, 5000.0, 3.1925169668291525e-19),
        (11.3, 5100.0, 1.5000399549235054e-17),
    ],
)
def test_against_cantera_320_linear_burke_forward_rate(
    rho: float, temp: float, expected: float
) -> None:
    # Independent Cantera 3.2.0, pinned SOURCE_SHA256, phase linear-Burke,
    # reaction 275 (zero-based). Set X=n/sum(n), P=kBT sum(n); extract only
    # forward_rate_constants[275]/(1000*N_A), not its NASA equilibrium constant.
    n = neutral_densities(ew.composition(rho, temp))
    assert water_coefficient(temp, effective_pressure(temp, n)) == pytest.approx(
        expected, rel=1e-10
    )


def test_removed_channel_time_matches_conservation_projected_linear_solve() -> None:
    temp = 6760.0
    n = neutral_densities(ew.composition(71, temp))
    root = np.sqrt(n)
    # Independent basis from the nullspace of elemental conservation, not
    # the implementation's chosen reaction columns; solve B rather than eig(B).
    _, _, vt = np.linalg.svd(ATOM_COUNTS * root)
    q = vt[2:].T
    matrix = np.zeros((6, 6))
    for reaction in REACTIONS:
        if reaction.source_number != 15:
            direction = reaction.stoichiometry / root
            matrix += reaction.forward_flux(temp, n) * np.outer(direction, direction)
    energy = q.T @ (GROUND_ENERGIES * root)
    expected = energy @ np.linalg.solve(q.T @ matrix @ q, energy) / (energy @ energy)
    assert relaxation(temp, n, water_association_scale=0).energy_time_s == pytest.approx(
        expected, rel=1e-9
    )


def test_pressure_dependent_matched_reverse_has_the_same_equilibrium_jacobian() -> None:
    temp = 5100.0
    n = neutral_densities(ew.composition(11.3, temp))
    direct = next(r for r in REACTIONS if r.source_number == 15)
    nu = direct.stoichiometry

    def production(populations: Vec) -> Vec:
        forward = (
            water_coefficient(temp, effective_pressure(temp, populations))
            * float(populations[0])
            * float(populations[3])
        )
        deviation = float(nu @ (np.log(populations) - equilibrium_log_reference(temp)))
        return -nu * forward * math.expm1(deviation)

    flux = water_coefficient(temp, effective_pressure(temp, n)) * n[0] * n[3]
    analytic = -flux * np.outer(nu, nu / n)
    measured = np.zeros((6, 6))
    for j in range(6):
        delta = np.zeros(6)
        delta[j] = 1e-5 * n[j]
        measured[:, j] = (production(n + delta) - production(n - delta)) / (2 * delta[j])
    assert np.linalg.norm(measured - analytic) / np.linalg.norm(analytic) < 1e-7
