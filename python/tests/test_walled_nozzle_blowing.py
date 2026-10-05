"""Tests for the B'(T, p) lookup (`blowing`) copied from the companion's H1 table."""

from __future__ import annotations

import math

import pytest

from puffsat.walled_nozzle import blowing


@pytest.mark.parametrize(
    ("temp", "bar", "expected"),
    [(2500, 818, 0.70), (3100, 62, 1.00), (3900, 62, 5.36), (3900, 818, 5.28), (4100, 10, 13.51)],
)
def test_reproduces_the_table_nodes(temp: float, bar: float, expected: float) -> None:
    assert blowing.blowing_parameter(temp, bar * 1e5) == pytest.approx(expected, rel=1e-6)


def test_rises_with_temperature_and_is_flat_in_pressure_when_hot() -> None:
    p = 300e5
    values = [blowing.blowing_parameter(t, p) for t in (2000, 2500, 3000, 3500, 3900)]
    assert values == sorted(values)
    # 3900 K row spans 5.27-5.55 above 30 bar: the companion's "flat" claim.
    hot = [blowing.blowing_parameter(3900, bar * 1e5) for bar in (30, 100, 300, 818)]
    assert max(hot) / min(hot) < 1.06


def test_cold_tail_is_positive_and_continuous_at_the_table_edge() -> None:
    edge = blowing.blowing_parameter(2500.0, 100e5)
    just_below = blowing.blowing_parameter(2499.0, 100e5)
    assert 0.0 < blowing.blowing_parameter(1000.0, 100e5) < just_below < edge
    assert just_below == pytest.approx(edge, rel=0.01)


def test_blowing_factors_match_the_papers_5p3_to_1p8() -> None:
    assert blowing.mass_flux_factor(5.3) == pytest.approx(1.84, abs=0.005)
    assert blowing.blowing_correction(5.3) == pytest.approx(math.log(6.3) / 5.3)
    # Small B': the correction tends to 1 (no blowing).
    assert blowing.blowing_correction(1e-6) == pytest.approx(1.0, abs=1e-5)
