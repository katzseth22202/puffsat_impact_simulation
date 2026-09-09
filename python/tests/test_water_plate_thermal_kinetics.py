"""Energy-constrained local chemistry and rate-sensitivity acceptance checks."""

import math

import numpy as np
import pytest

from puffsat import eos_water as ew
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate import thermal_kinetics as tk
from puffsat.water_plate.freeze import gas_parameters


@pytest.mark.parametrize("temp", [2000.0, 5100.0, 9000.0])
def test_species_internal_energies_obey_concentration_vant_hoff(temp: float) -> None:
    energy, _ = tk.species_thermo(temp)
    dt = 1e-4 * temp
    for r in kn.REACTIONS:
        measured = (
            kn.log_equilibrium_constant(r, temp + dt) - kn.log_equilibrium_constant(r, temp - dt)
        ) / (2 * dt)
        expected = float(r.stoichiometry @ energy) / (ew.K_B * temp**2)
        assert measured == pytest.approx(expected, rel=1e-6)


def test_fixed_composition_heat_capacity_not_equilibrium_heat_capacity() -> None:
    rho, temp = 71.0, 6760.0
    comp = ew.composition(rho, temp)
    energy, cv = tk.species_thermo(temp)
    gas = gas_parameters(rho, temp)
    dt = 0.1
    expected = rho * (gas.thermal(temp + dt) - gas.thermal(temp - dt)) / (2 * dt)
    assert tk.frozen_cv(comp, cv, 0.0) == pytest.approx(expected, rel=1e-8)
    assert np.all(np.isfinite(energy))


@pytest.mark.parametrize("rho,temp", [(1.0, 5000.0), (71.0, 6760.0), (279.0, 9000.0)])
def test_thermal_response_reduces_to_existing_isothermal_time(rho: float, temp: float) -> None:
    comp = ew.composition(rho, temp)
    n = kn.neutral_densities(comp)
    u, cv = tk.species_thermo(temp)
    result = tk.response(temp, n, u, tk.frozen_cv(comp, cv, 0.0))
    assert result.isothermal_s == pytest.approx(kn.relaxation(temp, n).energy_time_s, rel=1e-8)
    bath = tk.response(temp, n, u, 1e20 * tk.frozen_cv(comp, cv, 0.0))
    assert bath.isoenergetic_s == pytest.approx(result.isothermal_s, rel=1e-8)
    assert sum(result.isothermal_sensitivities) == pytest.approx(1, abs=1e-8)
    assert sum(result.isoenergetic_sensitivities) == pytest.approx(1, abs=1e-8)


def test_rate_sensitivities_are_logarithmic_derivatives() -> None:
    temp = 6760.0
    comp = ew.composition(71, temp)
    n = kn.neutral_densities(comp)
    u, cv = tk.species_thermo(temp)
    capacity = tk.frozen_cv(comp, cv, 0.0)
    base = tk.response(temp, n, u, capacity)
    for i in range(len(kn.REACTIONS)):
        lo, hi = np.ones(10), np.ones(10)
        lo[i], hi[i] = math.exp(-1e-3), math.exp(1e-3)
        left, right = (tk.response(temp, n, u, capacity, s) for s in (lo, hi))
        for name in ("isothermal", "isoenergetic"):
            measured = -math.log(getattr(right, name + "_s") / getattr(left, name + "_s")) / 0.002
            assert measured == pytest.approx(getattr(base, name + "_sensitivities")[i], abs=1e-7)


def test_group_removal_has_a_real_conserved_mode_not_a_fast_clipped_rate() -> None:
    comp = ew.composition(71, 6760)
    n = kn.neutral_densities(comp)
    u, cv = tk.species_thermo(6760)
    scales = np.array([0 if r.colliders is not None else 1 for r in kn.REACTIONS], dtype=float)
    with pytest.raises(ValueError, match=r"rank|spectrum"):
        tk.response(6760, n, u, tk.frozen_cv(comp, cv, 0.0), scales)


@pytest.mark.parametrize("capacity", [0, -1, math.inf, math.nan])
def test_invalid_heat_capacity_rejected(capacity: float) -> None:
    n = kn.neutral_densities(ew.composition(71, 6760))
    u, _ = tk.species_thermo(6760)
    with pytest.raises(ValueError):
        tk.response(6760, n, u, capacity)


def test_residual_curve_matches_source_and_refines() -> None:
    from puffsat.water_plate import flow_eos

    coarse, fine = tk.CvCurve.build(256), tk.CvCurve.build(512)
    errors = [0.0, 0.0]
    for rho in np.geomspace(1, 500, 31):
        exact = flow_eos.residual_reference(float(rho)).cv
        for j, curve in enumerate((coarse, fine)):
            errors[j] = max(errors[j], abs(curve.reference_cv(float(rho)) - exact))
    assert errors[1] < 0.15 * errors[0]
    assert errors[1] < 0.1  # J/kg/K at the 1000 K reference, smaller at loading T.


def test_isoenergetic_jacobian_matches_a_conserved_energy_temperature_solve() -> None:
    rho, temp = 11.3, 5100.0
    comp = ew.composition(rho, temp)
    n = kn.neutral_densities(comp)
    root_n = np.sqrt(n)
    u, cv = tk.species_thermo(temp)
    spectator_cv = 1.5 * ew.K_B * (comp.n_hp + comp.n_e + sum(comp.n_o_ions))
    residual_ref = 1000.0 * rho  # Independent nonzero residual cv at 1000 K.
    capacity = tk.frozen_cv(comp, cv, residual_ref * (1000 / temp) ** 2)
    nu = np.column_stack([r.stoichiometry for r in kn.REACTIONS])
    q, _ = np.linalg.qr(nu[:, [0, 1, 4, 7]] / root_n[:, None], mode="reduced")
    directions = q.T @ (nu / root_n[:, None])
    fluxes = np.array([r.forward_flux(temp, n) for r in kn.REACTIONS])
    matrix = (directions * fluxes) @ directions.T
    a = q.T @ (u * root_n)
    expected = -matrix @ (np.eye(4) + np.outer(a, a) / (ew.K_B * temp**2 * capacity))

    def energy(populations: kn.Vec, temperature: float) -> float:
        internal, _ = tk.species_thermo(temperature)
        return (
            float(populations @ internal)
            + spectator_cv * temperature
            - residual_ref * 1e6 / temperature
        )

    target = energy(n, temp)

    def production(populations: kn.Vec) -> kn.Vec:
        temperature = temp
        for _ in range(8):
            _, local_cv = tk.species_thermo(temperature)
            heat_capacity = (
                float(populations @ local_cv)
                + spectator_cv
                + residual_ref * (1000 / temperature) ** 2
            )
            temperature -= (energy(populations, temperature) - target) / heat_capacity
        assert abs(energy(populations, temperature) - target) < 1e-12 * target
        return q.T @ (kn.net_production(temperature, populations) / root_n)

    measured = np.zeros((4, 4))
    for j in range(4):
        direction = root_n * q[:, j]
        step = 1e-4 * float(np.min(n / np.maximum(np.abs(direction), 1e-100)))
        measured[:, j] = (production(n + step * direction) - production(n - step * direction)) / (
            2 * step
        )
    assert np.linalg.norm(measured - expected) / np.linalg.norm(expected) < 1e-6


def test_missing_thermal_rates_stay_in_the_audit_denominator() -> None:
    from puffsat.water_plate.kinetics_audit import Row
    from puffsat.water_plate.thermal_kinetics_audit import summarize

    common: Row = {"temperature_k": 6000, "pressure_pa": 1e8, "status": "resolved_local_network"}
    rows: list[Row] = [
        {
            **common,
            "weight_j": 6,
            "da_tracking": 100,
            "baseline_isoenergetic_time_over_baseline": 0.5,
            "baseline_isoenergetic_sensitivity_r1": 1.0,
        },
        {**common, "weight_j": 4, "baseline_error": "unresolved"},
    ]
    result = summarize(rows)
    assert result["baseline_isoenergetic_da_gt10_fraction"] == 0.6
    assert result["baseline_isoenergetic_resolved_rate_fraction"] == 0.6
    assert result["baseline_isoenergetic_sensitivity_r1_mean"] == 1.0
    assert result["kinetics_validated"] is False


def test_thermal_audit_rejects_a_different_parent_residual() -> None:
    from puffsat.water_plate.thermal_kinetics_audit import validate_model_metadata

    validate_model_metadata(
        {
            "table_provenance": {
                "dense_hot_cv_residual_power": 2.0,
                "real_fluid_max_k": 1000.0,
            }
        }
    )
    with pytest.raises(ValueError, match="parent closure"):
        validate_model_metadata(
            {
                "table_provenance": {
                    "dense_hot_cv_residual_power": 1.0,
                    "real_fluid_max_k": 1000.0,
                }
            }
        )


def test_conserved_elemental_energy_gauge_does_not_change_thermal_response() -> None:
    comp = ew.composition(71, 6760)
    n = kn.neutral_densities(comp)
    u, cv = tk.species_thermo(6760)
    capacity = tk.frozen_cv(comp, cv, 0.0)
    shifted = u + kn.ATOM_COUNTS.T @ np.array([ew.D_AT, -2 * ew.D_AT])
    result, gauge = (tk.response(6760, n, energy, capacity) for energy in (u, shifted))
    assert gauge.isoenergetic_s == pytest.approx(result.isoenergetic_s, rel=1e-8)
    assert gauge.isoenergetic_sensitivities == pytest.approx(
        result.isoenergetic_sensitivities, abs=1e-8
    )
