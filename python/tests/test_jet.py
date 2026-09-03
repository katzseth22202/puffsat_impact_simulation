"""Acceptance tests for N2 -- `eta_jet`, and whether the flare turns the pancake.

The rung's exit criterion is that the adiabatic conversion reproduces its own two limits, which
are known before any nozzle exists: **no flare changes nothing**, and **an infinite flare makes
the exhaust perfectly axial**. Everything else is that formula evaluated at the flown geometry.

The test that carries the argument is `test_the_flare_undoes_the_pancake`: Rung 4 measured
`alpha` = 0.088, and if the flare could not lift it back toward 1 then the column geometry would
have cost the design its jet rather than merely its baseline.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import jet

# ---- The adiabatic conversion, against its own limits ------------------------------------------


def test_no_flare_changes_nothing() -> None:
    """`A/A*` = 1 is a straight pipe: `B` never falls, so `mu` conservation moves no energy."""
    for alpha in (0.088, 1.0 / 3.0, 0.9):
        assert jet.alpha_after_expansion(alpha, 1.0) == pytest.approx(alpha)


def test_an_unbounded_flare_makes_the_exhaust_perfectly_axial() -> None:
    """`v_perp^2` goes as `B`, so with `B -> 0` every joule ends up axial whatever it started as."""
    assert jet.alpha_after_expansion(0.088, 1e6) == pytest.approx(1.0, abs=1e-5)
    assert jet.alpha_after_expansion(0.001, 1e6) == pytest.approx(1.0, abs=1e-5)


def test_the_conversion_is_monotone_in_the_flare() -> None:
    """More flare must never give less axial share, or the mechanism is stated backwards."""
    ratios = [1.0, 2.0, 4.0, 8.0, 11.3, 20.0]
    alphas = [jet.alpha_after_expansion(jet.ALPHA_MEASURED, r) for r in ratios]
    assert all(b >= a for a, b in pairwise(alphas))


def test_the_flare_needed_for_a_target_inverts_the_conversion() -> None:
    """`area_ratio_for_alpha` and `alpha_after_expansion` must be one relation."""
    for target in (0.5, 0.772, 0.95):
        ar = jet.area_ratio_for_alpha(jet.ALPHA_MEASURED, target)
        assert jet.alpha_after_expansion(jet.ALPHA_MEASURED, ar) == pytest.approx(target)


# ---- The finding ---------------------------------------------------------------------------------


def test_the_flare_undoes_the_pancake() -> None:
    """**The N2 answer.** A measured `alpha` of 0.088 becomes a directed jet at the flown flare.

    `alpha_0` enters the conversion only through `(1 - alpha_0)`, which is at most 1, so the
    starting anisotropy is bounded in how much it can cost. The flare is not a mitigation for the
    pancake -- it is the mechanism that undoes it.
    """
    flown = jet.alpha_after_expansion(jet.ALPHA_MEASURED, jet.FLOWN_AREA_RATIO)
    window = jet.alpha_after_expansion(jet.ALPHA_MEASURED, jet.WINDOW_AREA_RATIO)
    assert flown == pytest.approx(0.772, abs=1e-3)
    assert window == pytest.approx(0.919, abs=1e-3)
    # Both are far past the isotropic 1/3 the paper assumed at the *throat*.
    assert flown > 2.0 * jet.ALPHA_ISOTROPIC


def test_the_spread_formula_reproduces_the_papers_own_baseline() -> None:
    """At `alpha` = 1/3 the Gaussian family must return the paper's 0.461, or it is a different
    formula from the one `eq:reflection_baseline` uses."""
    assert jet.eta_geom_spread(jet.ALPHA_ISOTROPIC) == pytest.approx(0.4607, abs=1e-4)


def test_the_two_eta_geom_readings_bracket_and_do_not_cross() -> None:
    """Directed exhaust is the optimistic bound, a Gaussian spread the pessimistic one.

    Their ratio is the `sqrt(2/pi)` = 0.798 of a Maxwellian, at every `alpha` -- so the bracket is
    a fixed width in ratio and cannot invert.
    """
    for alpha in (0.088, 0.5, 0.772, 0.95):
        lo, hi = jet.eta_geom_spread(alpha), jet.eta_geom_directed(alpha)
        assert lo < hi
        assert lo / hi == pytest.approx(math.sqrt(2.0 / math.pi), rel=1e-9)


def test_eta_geom_clears_the_asks_calibration_by_a_wide_margin() -> None:
    """The ask: ">= 130% of the baseline supports the paper, < 100% is a serious problem."

    Against the *measured* baseline the flown flare returns ~370%, and against the isotropic one
    ~190%. The geometry is not what threatens the chain.
    """
    hi = jet.eta_geom_directed(jet.alpha_after_expansion(jet.ALPHA_MEASURED, jet.FLOWN_AREA_RATIO))
    measured_baseline = jet.eta_geom_spread(jet.ALPHA_MEASURED)
    isotropic_baseline = jet.eta_geom_spread(jet.ALPHA_ISOTROPIC)
    assert hi / measured_baseline > 3.0
    assert hi / isotropic_baseline > 1.3


def test_the_cold_leg_falls_short_on_chemistry_not_on_geometry() -> None:
    """**The reframing.** Where `eta_jet` misses the target it is `eta_chem` doing it.

    At 75 km/s the chain clears 0.775; at 45.58 it does not, and the only thing that changed is
    `eta_chem` (0.910 -> 0.731). The geometric term is identical at both speeds because the flare
    does not care how hot the pulse was.
    """
    hot = jet.jet_efficiency(closing_speed=75.0, area_ratio=jet.FLOWN_AREA_RATIO)
    cold = jet.jet_efficiency(closing_speed=45.58, area_ratio=jet.FLOWN_AREA_RATIO)
    assert hot.eta_geom_hi == pytest.approx(cold.eta_geom_hi)
    assert hot.clears_target
    assert not cold.clears_target
    assert cold.eta_chem < hot.eta_chem


# ---- The overtake leg's mirror -------------------------------------------------------------------


def test_the_drift_pushes_the_pitch_angle_toward_the_loss_cone() -> None:
    """The drift is motion *along* the field -- the one direction a mirror cannot reverse.

    Ignoring it would have reported `sin^2(theta)` = 0.912 instead of 0.69, overstating the
    mirror's margin by a quarter.
    """
    with_drift = jet.mean_sin2_pitch(jet.ALPHA_MEASURED, 0.105)
    thermal_only = 1.0 - jet.ALPHA_MEASURED
    assert with_drift < thermal_only
    assert with_drift == pytest.approx(0.69, abs=0.01)


def test_the_pancake_reflects_in_essentially_any_mirror() -> None:
    """Leg 1's mirror is not the binding constraint, and the pancake is why.

    Reflection acts on the perpendicular component, so a plume that is nearly all perpendicular
    sits far outside the loss cone. The minimum mirror ratio that still turns the flown plume is
    about 1.45 — below anything the graded column comes near.
    """
    assert jet.reflected_fraction(jet.ALPHA_MEASURED, 4.0) == 1.0
    assert jet.reflected_fraction(jet.ALPHA_MEASURED, 11.3) == 1.0
    minimum = 1.0 / jet.mean_sin2_pitch(jet.ALPHA_MEASURED, 0.105)
    assert minimum == pytest.approx(1.45, abs=0.05)


def test_a_mirror_that_does_not_mirror_reflects_nothing() -> None:
    """`R` = 1 is a straight field: no gradient, so nothing turns around."""
    assert jet.reflected_fraction(jet.ALPHA_MEASURED, 1.0) == 0.0
    assert jet.loss_cone_fraction(1.0) == 1.0


def test_both_legs_share_eta_geom_and_the_asymmetry_lives_elsewhere() -> None:
    """`mu` conversion does not care which way the drift points, and the mirror reflects it all.

    The leg asymmetry is the impulse ledger's `+1` (overtake) against `-1` (head-on), which is
    `aim_is_all_you_need`'s bookkeeping and deliberately not re-derived here.
    """
    head_on = jet.jet_efficiency(leg="head-on", area_ratio=jet.FLOWN_AREA_RATIO)
    overtake = jet.jet_efficiency(leg="overtake", area_ratio=jet.FLOWN_AREA_RATIO)
    assert head_on.eta_geom_hi == pytest.approx(overtake.eta_geom_hi)
    assert overtake.reflected == 1.0


# ---- The model's own assumption ------------------------------------------------------------------


def test_the_guiding_centre_picture_holds_by_orders_of_magnitude() -> None:
    """`mu` is invariant only while the field varies slowly over a gyroradius.

    If this failed, every number in the module would be void — so it is checked rather than
    asserted in prose.
    """
    for b_field, length in ((20.0, 3.0), (9.0, 6.0), (5.0, 23.0)):
        assert jet.adiabaticity_parameter(b_field, length) < 1e-4


def test_a_weaker_field_is_less_adiabatic() -> None:
    """The gyroradius grows as `B` falls, so the invariant is weakest at the exit — the scaling
    has to come out that way or the parameter is not measuring what it claims."""
    assert jet.adiabaticity_parameter(5.0, 3.0) > jet.adiabaticity_parameter(20.0, 3.0)
