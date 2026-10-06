"""Closed-form helpers of `plate_shape` (ADR-0055, Q62)."""

from __future__ import annotations

import math

import pytest

from puffsat.water_plate import plate_shape as ps


def test_a_flat_dish_has_the_disk_area() -> None:
    assert ps.paraboloid_area(10.0, 0.0) == pytest.approx(math.pi * 100.0)


def test_a_5_m_deep_20_m_bowl_has_its_hand_worked_area() -> None:
    # pi r / (6 h^2) [(r^2 + 4 h^2)^1.5 - r^3] = pi 10 / 150 (200^1.5 - 1000) = 382.9 m^2.
    assert ps.paraboloid_area(10.0, 5.0) == pytest.approx(382.9, rel=1e-3)


def test_eta_jet_inverts_the_overtake_impulse_model() -> None:
    # beta = eta sqrt(1 + k) + 1 and share = beta / (1 + sqrt(1 + k)): at k = 8 a share of
    # (0.6 * 3 + 1) / 4 = 0.7 gives back 0.6.
    assert ps.eta_jet(0.7, 8.0) == pytest.approx(0.6)
