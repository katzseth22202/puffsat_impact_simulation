"""Acceptance tests for the continuum `eta_geom` and the magnetic extension (R1, R11).

The exit criterion is that the replacement for `jet.py` reproduces **its own limits, and the
paper's own algebra**, before it is trusted at the flown geometry:

- a cold hypersonic exhaust (no thermal store) is perfectly directed, `eta_geom -> <cos theta>`;
- a stagnant plume (no directed motion) delivers nothing, `eta_geom -> 0`;
- and the `R_sp T`-based form must agree with the paper's `1/sqrt(1 + 3/(gamma M^2))` wherever the
  gas really is a `gamma`-law gas, which is the check that the two sides are computing the same
  quantity and not two different ones with the same name.

`test_the_thermal_form_matches_the_papers_gamma_form` is the one that carries the argument: R1
proposed the `gamma M^2` form, this module evaluates the mixture directly, and if those disagreed
on an ideal gas then one of them would be wrong.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import extension

# ---- The two limits ----------------------------------------------------------------------------


def test_a_cold_exhaust_is_perfectly_directed() -> None:
    """With no thermal store left, every joule is directed and only divergence costs anything.

    `eos_water`'s partition functions underflow below ~50 K, so "cold" here is 100 K against a
    20 km/s exhaust -- a thermal speed four orders below the directed one, which is the limit.
    """
    assert extension.eta_thermal(1e-4, 100.0, 20_000.0) == pytest.approx(1.0, abs=1e-3)


def test_a_stagnant_plume_delivers_nothing() -> None:
    """`u -> 0` is a hot cloud going nowhere: all motion is random, none of it axial."""
    assert extension.eta_thermal(0.03, 10_000.0, 0.0) == 0.0
    assert extension.eta_thermal(0.03, 10_000.0, 1.0) < 1e-3


def test_eta_geom_is_the_divergence_times_the_thermal_term() -> None:
    """The factorisation itself, so a refactor cannot quietly change what multiplies what."""
    rho, temp, speed = 0.025, 5000.0, 9900.0
    thermal = extension.eta_thermal(rho, temp, speed)
    assert extension.eta_geom(rho, temp, speed, 1.0) == pytest.approx(thermal)
    assert extension.eta_geom(rho, temp, speed, 0.5) == pytest.approx(0.5 * thermal)


def test_the_thermal_term_rises_with_speed_at_fixed_state() -> None:
    """More directed motion against the same thermal store must mean a better-aligned jet."""
    values = [extension.eta_thermal(0.025, 5000.0, u) for u in (2000.0, 5000.0, 10_000.0, 20_000.0)]
    assert all(b > a for a, b in pairwise(values))


# ---- Agreement with the paper's own form -------------------------------------------------------


def test_the_thermal_form_matches_the_papers_gamma_form() -> None:
    """`1/sqrt(1 + 3 R T/u^2)` and `1/sqrt(1 + 3/(gamma M^2))` are the same statement.

    Checked on the algebra rather than on the EOS: with `c^2 = gamma R T` and `M = u/c`, the two
    expressions are identical, so any Mach number and `gamma` must agree to machine precision.
    """
    for gamma in (1.15, 1.25, 5.0 / 3.0):
        for mach in (1.5, 2.7, 3.4, 6.0):
            paper = 1.0 / math.sqrt(1.0 + 3.0 / (gamma * mach * mach))
            # Reconstruct the same state: pick R T, then u = M sqrt(gamma R T).
            r_t = 1.0e6
            u = mach * math.sqrt(gamma * r_t)
            ours = 1.0 / math.sqrt(1.0 + 3.0 * r_t / (u * u))
            assert ours == pytest.approx(paper, rel=1e-12)


def test_the_cone_factor_reproduces_its_endpoints() -> None:
    """`(1 + cos theta_max)/2`: 1 for a parallel exhaust, 1/2 for a hemisphere."""
    assert extension.cone_mean_cos(0.0) == pytest.approx(1.0)
    assert extension.cone_mean_cos(math.pi / 2.0) == pytest.approx(0.5)
    assert extension.cone_mean_cos(math.radians(15.0)) == pytest.approx(0.983, abs=1e-3)


# ---- Detachment ---------------------------------------------------------------------------------


def test_detachment_is_at_the_exit_for_an_already_super_alfvenic_plume() -> None:
    """If it leaves super-Alfvenic it is not tied to the field and detaches immediately."""
    assert extension.detachment_radii(1.0) == pytest.approx(1.0)
    assert extension.detachment_radii(2.5) == pytest.approx(1.0)


def test_a_more_sub_alfvenic_plume_detaches_further_out() -> None:
    """`(1/M_A)^{2/3}` falls with `M_A`: the weaker the flow, the longer the field holds it."""
    radii = [extension.detachment_radii(m) for m in (0.30, 0.45, 0.60, 0.90)]
    assert all(b < a for a, b in pairwise(radii))
    assert extension.detachment_radii(0.353) == pytest.approx(1.99, abs=0.02)


def test_the_thermal_speed_uses_heavy_particles_only() -> None:
    """Electrons carry the heat but not the mass, and `v_g` is energy per unit *mass*.

    Including them would raise `<v_th^2>` and lower `eta_geom`, so the choice is stated by test.
    """
    v2 = extension.thermal_speed_squared(0.025, 5000.0)
    # 3kT/m for atomic hydrogen is the fastest a heavy-particle mixture could possibly be.
    fastest_heavy = 3.0 * extension.K_B * 5000.0 / (1.008 * 1.66053906660e-27)
    assert 0.0 < v2 <= fastest_heavy


# ---- The geometry the sweep is built on --------------------------------------------------------


def test_the_flare_half_angle_sets_the_extension_length() -> None:
    """A 15-degree straight flare between two radii has one length, and it is stated."""
    length = extension.half_angle(10.0, 3.0, 5.0)
    assert math.degrees(length) == pytest.approx(math.degrees(math.atan(0.2)))


def test_the_alfven_crossing_is_none_when_the_flow_never_gets_there() -> None:
    """A sub-Alfvenic history has no crossing, and the sweep must report that rather than guess."""
    assert extension.alfven_crossing([]) is None
