"""Acceptance tests for the graded column's magnetostatic field (Rung 1, N5).

The rung's exit criterion is that the off-axis machinery reproduces two things whose answers are
known before the solver exists: **a single loop's closed-form on-axis field**, and **the infinite
solenoid's `mu0 n I`**. Everything else here is downstream of those two.

The last test is the one that carries the physics rather than the algebra: a winding that is dense
enough has no interior `|B|` minimum, and one that is not does -- which is why N5's answer is a
constraint on windings rather than a yes or no.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import field

# ---- The elliptic integrals ----------------------------------------------------------------------


def test_elliptic_integrals_at_zero_are_both_a_quarter_turn() -> None:
    """`K(0) = E(0) = pi/2`; the degenerate case every on-axis reduction passes through."""
    k, e = field.elliptic_k_e(0.0)
    assert k == pytest.approx(math.pi / 2)
    assert e == pytest.approx(math.pi / 2)


def test_elliptic_integrals_match_tabulated_values_at_one_half() -> None:
    """AGM against the standard values, to well past the precision any field needs."""
    k, e = field.elliptic_k_e(0.5)
    assert k == pytest.approx(1.8540746773, rel=1e-10)
    assert e == pytest.approx(1.3506438810, rel=1e-10)


def test_the_singular_parameter_is_refused_rather_than_returned() -> None:
    """`m = 1` is a filament's own surface, where `K` diverges. Refusing beats returning inf."""
    with pytest.raises(ValueError):
        field.elliptic_k_e(1.0)


# ---- The loop, against its closed form -----------------------------------------------------------


def test_loop_field_reduces_to_the_closed_form_on_the_axis() -> None:
    """**The acceptance test.** The elliptic expression must reduce to `mu0 I a^2/(2(a^2+z^2)^1.5)`.

    Checked at a sequence of radii approaching the axis so this pins the *limit*, not one point.
    """
    a, current = 3.5, 1.0e6
    for z in (0.0, 1.0, 5.0, -4.0):
        exact = field.loop_field_on_axis(a, current, z)
        for r in (1e-4, 1e-5, 1e-6):
            _, b_z = field.loop_field(a, current, r, z)
            assert b_z == pytest.approx(exact, rel=1e-6)


def test_the_loop_field_is_symmetric_about_its_own_plane() -> None:
    """`B_z` is even in `z` and `B_r` is odd: the loop has no preferred axial direction."""
    a, current, r = 3.5, 1.0e6, 2.0
    b_r_up, b_z_up = field.loop_field(a, current, r, 1.7)
    b_r_dn, b_z_dn = field.loop_field(a, current, r, -1.7)
    assert b_z_up == pytest.approx(b_z_dn, rel=1e-12)
    assert b_r_up == pytest.approx(-b_r_dn, rel=1e-12)


def test_radial_field_vanishes_on_the_axis() -> None:
    """By symmetry there is nowhere for `B_r` to point on the axis."""
    b_r, _ = field.loop_field(3.5, 1.0e6, 0.0, 2.0)
    assert b_r == 0.0


def test_a_long_dense_stack_reproduces_the_finite_solenoid() -> None:
    """**The second acceptance test.** The centre field of a solenoid of finite length.

        B = mu0 n I L/sqrt(L^2 + 4a^2),      n = (N-1)/L for endpoint-inclusive spacing,

    which tends to `mu0 n I` as `L/a -> infinity`. Both corrections are checked rather than
    tolerated: taking the infinite-solenoid value and `N/L` instead was a 0.25% error, which is
    larger than the numerical error being tested for and would have made this test meaningless.
    """
    a, length, n_coils, current = 1.0, 40.0, 800, 1.0e4
    turns_per_m = (n_coils - 1) / length
    stack = field.CoilStack(
        tuple(
            field.Coil(z_m=-length / 2 + length * i / (n_coils - 1), radius_m=a, current_a=current)
            for i in range(n_coils)
        )
    )
    infinite = field.MU0 * turns_per_m * current
    expected = infinite * length / math.sqrt(length * length + 4.0 * a * a)
    assert stack.on_axis(0.0) == pytest.approx(expected, rel=1e-4)
    assert stack.magnitude(0.5 * a, 0.0) == pytest.approx(expected, rel=2e-3)
    assert expected == pytest.approx(infinite, rel=2e-3)  # the infinite limit is nearby


# ---- The design profile --------------------------------------------------------------------------


def test_the_papers_four_stations_lie_on_one_power_law() -> None:
    """20/12/9/5 T are not four requirements but one curve, `B ~ z^-0.44`, to under 2%.

    This is what makes interpolating between them safe, and it is a finding the paper does not
    state: through standoff it says the snowplow pressure falls as nearly `1/z`.
    """
    fit = field.fit_profile()
    assert fit.max_rel_error < field.PROFILE_FIT_MAX_ERROR
    assert fit.exponent == pytest.approx(-0.4405, abs=1e-3)
    assert 2.0 * fit.exponent == pytest.approx(-0.881, abs=2e-3)


def test_design_pressure_is_the_standoff_inverse_of_the_design_field() -> None:
    """`p = B^2/2mu0` and its inverse must be one relation, or beta is not 1 by construction."""
    for z in (1.0, 3.0, 6.0, 23.0):
        b = field.design_field(z)
        assert field.design_pressure(z) == pytest.approx(b * b / (2.0 * field.MU0))


def test_the_profile_is_held_flat_inboard_of_the_first_station() -> None:
    """The power law diverges at the origin; the chamber does not. Stated assumption, pinned."""
    assert field.design_field(0.1) == pytest.approx(field.design_field(1.0))


# ---- The fitted winding --------------------------------------------------------------------------


def test_a_built_winding_tracks_the_design_profile_past_the_chamber_knee() -> None:
    """Beyond the first two metres the winding follows the profile; that part is well posed."""
    stack = field.build_winding(72)
    assert field.profile_error(stack, z_min=2.0) < 0.15


def test_the_chamber_gradient_is_sharper_than_a_solenoid_can_make() -> None:
    """**The finding.** The paper's profile is not realizable by an all-positive winding.

    A solenoid smooths field structure over about its own radius. The profile demands 20 T at 1 m
    falling to 12 T at 3 m -- a gradient scale near 2 m -- from coils of 3.5 m radius, so the
    winding undershoots at the chamber by about a fifth. Since the field is there to stand off
    the plume, `B^2/2mu0 = p`, an undershoot in `B` is an overshoot in `beta` at the highest
    pressure station in the column.

    The mechanism is pinned by its scaling: the shortfall must grow with coil radius, because the
    smoothing length is the coil radius. If it did not, this would be an artifact.
    """
    at_bore = field.chamber_realizability(coil_radius=3.5)
    assert at_bore.shortfall > 0.15
    assert at_bore.beta_at_chamber > 1.4

    tight = field.chamber_realizability(coil_radius=1.5)
    wide = field.chamber_realizability(coil_radius=5.0)
    assert tight.shortfall < at_bore.shortfall < wide.shortfall


def test_every_current_in_a_built_winding_points_the_same_way() -> None:
    """No reverse windings. This is the property the discarded least-squares fit did not have.

    Its end rings came out at -29 MA, which manufactured `|B|` ripple that was an artifact of the
    regularisation rather than of the magnet. Nobody builds a nozzle with counter-wound end coils,
    so a solution containing them is answering a different question.
    """
    for n in (18, 72, 200):
        currents = [c.current_a for c in field.build_winding(n).coils]
        assert min(currents) > 0.0


def test_more_coils_does_not_degrade_the_on_axis_fit() -> None:
    """Refinement must converge, or the ripple sweep is measuring the fit instead of the winding."""
    coarse = field.profile_error(field.build_winding(36))
    fine = field.profile_error(field.build_winding(120))
    assert fine <= coarse * 1.5


# ---- The trap question ---------------------------------------------------------------------------


def test_the_loss_cone_reproduces_its_two_known_ends() -> None:
    """`R = 1` traps nothing; `R = 2` traps `1 - sqrt(1/2)`."""
    assert field.trapped_fraction(1.0) == 0.0
    assert field.trapped_fraction(2.0) == pytest.approx(1.0 - math.sqrt(0.5))


def test_a_dense_winding_has_no_interior_minimum_on_axis() -> None:
    """The paper's claim, at the continuous limit it is argued from -- and here it holds."""
    stack = field.build_winding(200)
    assert field.scan_axial(stack, 0.0).n_local_minima == 0


def test_ripple_deepens_toward_the_winding_and_with_fewer_coils() -> None:
    """**The physics test.** Whether N5's claim holds is a property of the winding, not the paper.

    Two monotonicities have to hold for the sweep to mean anything: at fixed coil count the ripple
    is worse nearer the conductor, and at fixed radius it is worse with fewer coils. If either
    failed, the mirror ratios reported for N5 would be noise.
    """
    sparse = field.build_winding(18)
    dense = field.build_winding(120)
    # 2.5 m, not the 3.0 m bore edge: a filamentary loop is singular on its own conductor at
    # 3.5 m, so ripple measured hard against the winding is the filament model failing rather
    # than the magnet trapping. A real coil pack has finite cross-section there.
    near = field.scan_axial(sparse, 2.5).mirror_ratio
    far = field.scan_axial(sparse, 0.0).mirror_ratio
    assert near >= far
    assert field.scan_axial(dense, 2.5).mirror_ratio <= near


# ---- Can the shortfall be fixed? -------------------------------------------------------------


def test_scaling_the_current_reaches_standoff_but_over_provisions_everywhere_else() -> None:
    """More current is a scale, and the shortfall is a ratio -- so it buys the chamber by excess.

    The chamber lands exactly on standoff (that is what the scale is solved for) and every other
    station overshoots. Overshoot is physically safe and costs stored energy, which is what the
    virial structure floor is charged on.
    """
    cost = field.scale_up_cost()
    chamber = dict(cost.overshoot)[1.0]
    assert chamber == pytest.approx(1.0, rel=1e-6)
    assert all(ratio > 1.2 for z, ratio in cost.overshoot if z > 1.0)
    assert cost.energy_scale == pytest.approx(cost.current_scale**2)


def test_a_longer_column_lowers_the_brush_ratio_and_the_shortfall_together() -> None:
    """**The physics test for the fix.** Shortfall must track `coil radius / gradient length`.

    If it did not, the brush-width explanation would be a story rather than the mechanism, and
    the column-length lever would not be the right one to reach for.
    """
    rows = field.column_length_scan()
    for shorter, longer in pairwise(rows):
        assert longer.brush_ratio < shorter.brush_ratio
        assert longer.shortfall < shorter.shortfall
        assert longer.beta_at_chamber < shorter.beta_at_chamber


def test_the_column_scan_reproduces_the_flown_geometry_it_shares() -> None:
    """The 23 m row and the standalone chamber check are the same calculation, so they must agree.

    They come down different paths -- one scales the stations and rebuilds, the other uses the
    module defaults -- so this is a real cross-check and not a tautology. It is loose because the
    scan derives its bore from `eq:bore_from_length` (3.02 m) while the default coil radius is a
    round 3.5 m.
    """
    flown = next(r for r in field.column_length_scan() if r.length_m == 23.0)
    direct = field.chamber_realizability(coil_radius=flown.coil_m)
    assert flown.shortfall == pytest.approx(direct.shortfall, abs=0.02)
