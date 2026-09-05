"""Acceptance tests for `phi`, the in-nozzle radiated share (R12).

`phi` is a ratio of two quadratures over the same isentrope, so the exit criterion is that the
quadrature is right and that the ratio's own bounds hold: it must lie in `[0, 1]`, it must fall
when the denominator is carried further, and it must be *reported alongside its convergence*,
because a truncated denominator biases `phi` high -- the direction that would falsely condemn the
liner.

`test_phi_falls_as_the_denominator_is_carried_further` is the one that carries the argument. A
`phi` that did not move with the tail length would mean the tail contributes nothing, and the
whole item exists because the paper suspects it contributes a great deal.
"""

from __future__ import annotations

import math

import pytest

from puffsat import radiance

COLD_T0 = 14_700.0
HOT_T0 = 26_200.0


@pytest.fixture(scope="module")
def hot_equilibrium() -> radiance.PhiResult:
    """One solved case, shared: the isentrope is this package's hot spot."""
    return radiance.phi(75.0, HOT_T0, branch="equilibrium", steps=120, stride=4)


# ---- The quadrature ------------------------------------------------------------------------------


def test_the_trapezoid_is_exact_on_a_constant_rate() -> None:
    """A constant `e/t_rad` over a known duration integrates to the obvious answer."""
    stations = [
        radiance.RadiatingStation(
            time_s=t,
            area_ratio=1.0,
            radius_m=3.0,
            rho=0.1,
            temp_k=1e4,
            energy_j_kg=2.0,
            cooling_time_s=1.0,
            regime="diffusion",
            inside=True,
        )
        for t in (0.0, 1.0, 2.0, 3.0)
    ]
    assert radiance.radiated_energy(stations) == pytest.approx(6.0)


def test_an_infinite_cooling_time_radiates_nothing() -> None:
    """A station that cannot cool must contribute zero, not a NaN or an exception."""
    station = radiance.RadiatingStation(
        time_s=0.0,
        area_ratio=1.0,
        radius_m=3.0,
        rho=0.1,
        temp_k=1e4,
        energy_j_kg=1e9,
        cooling_time_s=math.inf,
        regime="free-streaming",
        inside=False,
    )
    assert station.power_w_kg == 0.0


def test_a_single_station_integrates_to_zero() -> None:
    """No interval, no energy. Guards the `pairwise` edge case rather than the physics."""
    stations = [
        radiance.RadiatingStation(0.0, 1.0, 3.0, 0.1, 1e4, 1.0, 1.0, "diffusion", inside=True)
    ]
    assert radiance.radiated_energy(stations) == 0.0


# ---- The ratio's own bounds ----------------------------------------------------------------------


def test_phi_lies_between_zero_and_one(hot_equilibrium: radiance.PhiResult) -> None:
    """It is a share of a total that contains it, so this cannot fail unless the split is wrong."""
    assert 0.0 <= hot_equilibrium.phi <= 1.0


def test_the_inside_leg_is_a_strict_subset_of_the_total(
    hot_equilibrium: radiance.PhiResult,
) -> None:
    """The denominator carries the same stations plus the free jet, so it must be larger."""
    assert hot_equilibrium.total_j_kg > hot_equilibrium.inside_j_kg > 0.0


def test_the_transit_matches_the_solved_cooling_history(
    hot_equilibrium: radiance.PhiResult,
) -> None:
    """The numerator's clock is `expansion.cooling_history`'s, which is 1.7-2.8 ms."""
    assert 1.0 < hot_equilibrium.transit_ms < 4.0


def test_the_free_jet_takes_far_longer_than_the_nozzle(
    hot_equilibrium: radiance.PhiResult,
) -> None:
    """The paper's own framing: 'the plume radiates for far longer downstream, into nothing.'"""
    assert hot_equilibrium.tail_ms > hot_equilibrium.transit_ms


def test_phi_falls_as_the_denominator_is_carried_further() -> None:
    """**The convergence check.** A longer tail can only add to the denominator.

    This is also the reason `last_decade_share` is reported: if `phi` were insensitive to the tail
    length the item would be trivial, and if it never converged the number would be meaningless.
    """
    short = radiance.phi(
        75.0, HOT_T0, branch="equilibrium", steps=120, stride=4, expansion_ratio=512.0
    )
    long = radiance.phi(
        75.0, HOT_T0, branch="equilibrium", steps=120, stride=4, expansion_ratio=4096.0
    )
    assert long.phi <= short.phi + 1e-9
    assert long.total_j_kg >= short.total_j_kg


# ---- The branches --------------------------------------------------------------------------------


def test_both_branches_solve() -> None:
    """The frozen branch cools below the EOS floor sooner; the tail must shorten, not crash."""
    for branch in ("equilibrium", "frozen"):
        result = radiance.phi(45.58, COLD_T0, branch=branch, steps=120, stride=4)
        assert not math.isnan(result.phi)
        assert 0.0 <= result.phi <= 1.0


def test_the_perturbative_assumption_is_reported_not_assumed(
    hot_equilibrium: radiance.PhiResult,
) -> None:
    """`E_rad,inside / e_throat` is carried so the reader can see whether the isentrope is valid."""
    assert hot_equilibrium.inside_fraction_of_internal >= 0.0
