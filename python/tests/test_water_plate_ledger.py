"""Conservation and geometric acceptance checks for the overtake reference ledger."""

import math

import pytest

from puffsat.water_plate.ledger import coating_ledger, pulse_ledger


@pytest.mark.parametrize("k", [0.0, 5.0, 8.5, 10.0])
@pytest.mark.parametrize("speed", [45_580.0, 65_130.0])
def test_merge_and_axial_ceiling_conserve_energy_and_momentum(k: float, speed: float) -> None:
    row = pulse_ledger("test", speed, k)
    merged_mass = row.projectile_kg + row.water_kg
    assert merged_mass * row.merge_speed_m_s == pytest.approx(row.incoming_momentum_ns)
    assert row.merge_bulk_j + row.merge_heat_j == pytest.approx(row.incoming_ke_j)
    # Construct the single-speed axial exhaust saturating Cauchy-Schwarz. The
    # overtake vehicle receives an incoming credit plus backward exhaust momentum.
    exhaust_momentum = row.ceiling_impulse_ns - row.incoming_momentum_ns
    assert exhaust_momentum**2 / (2.0 * merged_mass) == pytest.approx(row.incoming_ke_j)
    assert row.ceiling_beta == pytest.approx(1.0 + math.sqrt(1.0 + k))
    if k == 0.0:
        assert row.water_only_ceiling_isp_s is None
        assert row.gasification_j == 0.0
    else:
        assert row.water_only_ceiling_isp_s is not None


def test_reference_water_bill_and_distinct_optimization_metrics() -> None:
    row = pulse_ledger("3-synodic", 45_580.0, 8.5)
    assert row.water_kg == 212.5
    assert row.gasification_j == pytest.approx(552.5e6)
    assert row.gasification_fraction_ke == pytest.approx(0.021275, rel=1e-4)
    low, high = (pulse_ledger("test", 45_580.0, k) for k in (5.0, 10.0))
    assert high.ceiling_impulse_ns > low.ceiling_impulse_ns
    assert high.water_only_ceiling_isp_s is not None
    assert low.water_only_ceiling_isp_s is not None
    assert high.water_only_ceiling_isp_s < low.water_only_ceiling_isp_s


@pytest.mark.parametrize("depth", [0.0, 0.10, 0.15])
def test_coating_area_against_surface_quadrature(depth: float) -> None:
    row = coating_ledger(30.0, depth)
    # Independent midpoint integration of 2*pi*r*sqrt(1+(dz/dr)^2) dr.
    radius, n = 15.0, 10_000
    dr = radius / n
    area = sum(
        2.0
        * math.pi
        * (r := (i + 0.5) * dr)
        * math.sqrt(1.0 + (2.0 * (depth * 30.0) * r / radius**2) ** 2)
        * dr
        for i in range(n)
    )
    assert row.surface_area_m2 == pytest.approx(area, rel=1e-8)
    assert row.layer_inventory_kg == pytest.approx(0.1 * area, rel=1e-8)
    assert row.full_renewal_fraction_reference_water == pytest.approx(
        row.layer_inventory_kg / 212.5
    )


def test_thickness_changes_inventory_not_surface_or_assumed_consumption() -> None:
    bare = coating_ledger(30.0, 0.0, thickness_m=0.0)
    thin = coating_ledger(30.0, 0.0, thickness_m=10e-6)
    thick = coating_ledger(30.0, 0.0)
    assert bare.layer_inventory_kg == 0.0
    assert bare.surface_area_m2 == thin.surface_area_m2 == thick.surface_area_m2
    assert thick.layer_inventory_kg == pytest.approx(70.6858347)
    assert thick.layer_inventory_kg == pytest.approx(10.0 * thin.layer_inventory_kg)


@pytest.mark.parametrize("bad", [-1.0, math.nan, math.inf])
def test_invalid_inputs_are_rejected(bad: float) -> None:
    with pytest.raises(ValueError):
        pulse_ledger("test", 45_580.0, bad)
    with pytest.raises(ValueError):
        pulse_ledger("test", bad, 8.5)
    with pytest.raises(ValueError):
        coating_ledger(30.0, bad)
