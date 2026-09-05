"""Acceptance tests for the residence criterion (R10, and the correct form of N5).

The exit criterion is the **isentropic area-Mach relation**, whose values are known before this
module exists: `A/A*` is exactly 1 at `M` = 1, has its minimum there, and rises on both sides.
Those are the properties the whole argument rests on -- specifically that the margin *vanishes*
at the sonic point, which is what moves the binding station from the winding to the chamber.

`test_the_margin_vanishes_at_the_sonic_point` is the one that carries the argument. If the margin
did not collapse there, ripple would be a downstream problem as P8 assumed rather than a chamber
one, and R10's flat-shelf caveat would not bite.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import residence

# ---- The area-Mach relation --------------------------------------------------------------------


def test_the_margin_vanishes_at_the_sonic_point() -> None:
    """`A/A*` = 1 at `M` = 1 for every `gamma`: sonic flow survives no contraction at all."""
    for gamma in residence.GAMMA_BRACKET:
        assert residence.area_mach(1.0, gamma) == pytest.approx(1.0, rel=1e-12)


def test_the_area_ratio_has_its_minimum_at_the_throat() -> None:
    """It rises on *both* sides of `M` = 1 -- which is why `A/A*` cannot bracket the throat."""
    for gamma in residence.GAMMA_BRACKET:
        assert residence.area_mach(0.5, gamma) > 1.0
        assert residence.area_mach(2.0, gamma) > 1.0


def test_the_margin_grows_steeply_with_mach_number() -> None:
    """The ordering that makes a downstream ripple harmless and an upstream one dangerous."""
    values = [residence.area_mach(m, 1.25) for m in (1.1, 1.5, 2.0, 2.7, 3.4)]
    assert all(b > a for a, b in pairwise(values))
    assert values[0] < 1.02
    assert values[-1] > 5.0


def test_the_area_mach_relation_matches_published_values() -> None:
    """Spot-check against the standard `gamma` = 1.4 table, which is tabulated everywhere.

    This is the external anchor: the rest of the module is arithmetic on top of this relation, so
    if it agrees with the textbook the criterion is being evaluated correctly.
    """
    assert residence.area_mach(2.0, 1.4) == pytest.approx(1.6875, rel=1e-4)
    assert residence.area_mach(3.0, 1.4) == pytest.approx(4.2346, rel=1e-4)
    assert residence.area_mach(0.5, 1.4) == pytest.approx(1.3398, rel=1e-4)


def test_a_non_positive_mach_number_is_refused() -> None:
    """Stagnant flow has no area ratio; returning one would be a silently wrong margin."""
    with pytest.raises(ValueError):
        residence.area_mach(0.0, 1.25)


# ---- The criterion that is being retired -------------------------------------------------------


def test_the_loss_cone_threshold_is_computed_but_flagged() -> None:
    """R10's own number, reproduced so the answer can say what it would have concluded."""
    assert residence.loss_cone_threshold(0.088) == pytest.approx(1.096, abs=1e-3)
    assert residence.loss_cone_threshold(1.0 / 3.0) == pytest.approx(1.5, abs=1e-3)


def test_an_isotropic_plume_needs_a_much_deeper_mirror_than_a_pancake() -> None:
    """The interaction R10 found between P1 and P8, kept because the *direction* is right."""
    assert residence.loss_cone_threshold(1.0 / 3.0) > residence.loss_cone_threshold(0.088)


# ---- Finding the minima ------------------------------------------------------------------------


def test_a_smooth_field_has_no_interior_minima() -> None:
    """A dense winding approaches a current sheet, and a monotone profile has nothing to trap."""
    from puffsat import field

    stack = field.build_winding(120, length=23.8, profile=field.capped_profile())
    assert residence.local_minima(stack, 0.0, length=23.8, n=200) == []


def test_a_sparse_winding_ripples_and_the_ripple_deepens_toward_it() -> None:
    """The physical ordering: `|B|` peaks at each coil and dips between, worse near the wall."""
    from puffsat import field

    stack = field.build_winding(12, coil_radius=3.5, length=23.8)
    near = residence.local_minima(stack, 2.5, length=23.8, n=400)
    far = residence.local_minima(stack, 0.0, length=23.8, n=400)
    assert near
    deepest_near = max(m.mirror_ratio for m in near)
    deepest_far = max((m.mirror_ratio for m in far), default=1.0)
    assert deepest_near > deepest_far


def test_every_reported_minimum_is_deeper_than_flat() -> None:
    """A mirror ratio below 1 would mean the 'minimum' is not one."""
    from puffsat import field

    stack = field.build_winding(12, coil_radius=3.5, length=23.8)
    for m in residence.local_minima(stack, 2.5, length=23.8, n=400):
        assert m.mirror_ratio >= 1.0
        assert 0.0 <= m.z_m <= 23.8


def test_the_mach_profile_starts_sonic_and_rises() -> None:
    """Read off the shipped cooling history: `M` = 1 at the chamber, supersonic downstream."""
    mach = residence.mach_profile()
    assert mach(0.0) == pytest.approx(1.0, abs=0.01)
    assert mach(23.8) > 2.0
    assert not math.isnan(mach(6.0))
