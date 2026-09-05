"""Acceptance tests for the flux-tube accounting (R2, R5).

The exit criterion is an **analytic invariant the tracer cannot fake**: a flux tube conserves
`B r^2`, so a correctly traced line must hold that product constant along its length. That is the
strongest single check available here, because it couples the tracer, the off-axis loop field and
the elliptic integrals in one number and would fail on an error in any of them.

The geometric checks are secondary but load-bearing for the answer: `flared_radius` has to
reproduce ADR-0011's stated 3.50 m to 5.17 m, or this module is answering a question about a
different magnet than the paper flies.
"""

from __future__ import annotations

import math

import pytest

from puffsat import field, fluxtube


@pytest.fixture(scope="module")
def capped_stack() -> field.CoilStack:
    """The winding the paper currently flies: flared contour, 12 T shelf."""
    stack, _profile = fluxtube.build_stacks(n_coils=32)["flared-capped"]
    return stack


# ---- The invariant -----------------------------------------------------------------------------


def test_a_traced_tube_conserves_b_r_squared(capped_stack: field.CoilStack) -> None:
    """`B r^2` is constant along a flux tube. This is the tracer's acceptance test.

    The tolerance is set by the RK4 step and by the fact that `|B|` on the traced path is the
    *magnitude*, which picks up the radial component; deep inside a solenoid bore that is a
    sub-percent effect and the invariant holds tightly.
    """
    tube = fluxtube.trace_tube(capped_stack, 1.5, z_start=2.0, z_end=20.0, steps=400)
    products = [b * r * r for b, r in zip(tube.b_t, tube.r, strict=True)]
    spread = (max(products) - min(products)) / (sum(products) / len(products))
    assert spread < 0.05, f"B r^2 drifts by {100 * spread:.1f}% along the tube"


def test_a_tube_launched_on_axis_stays_on_axis(capped_stack: field.CoilStack) -> None:
    """`B_r` vanishes on the axis by symmetry, so the axis is a field line."""
    tube = fluxtube.trace_tube(capped_stack, 1e-6, z_start=1.0, z_end=20.0, steps=200)
    assert max(tube.r) < 1e-3
    assert max(abs(t) for t in tube.theta_rad) < 1e-3


def test_tubes_widen_where_the_field_falls(capped_stack: field.CoilStack) -> None:
    """The whole geometric argument in one assertion: a weakening field is a widening tube."""
    tube = fluxtube.trace_tube(capped_stack, 2.0, z_start=1.0, z_end=20.0, steps=400)
    assert tube.r[-1] > tube.r[0]
    assert tube.b_t[-1] < tube.b_t[0]


# ---- The winding contour -----------------------------------------------------------------------


def test_the_flare_reproduces_the_paper_s_stated_contour() -> None:
    """ADR-0011 as amended: 3.50 m at the chamber to 5.17 m at the exit, at the 12 T cap."""
    capped = field.capped_profile(field.ADR0012_CAP_T)
    assert fluxtube.flared_radius(0.0, capped) == pytest.approx(3.50, abs=0.01)
    assert fluxtube.flared_radius(23.8, capped) == pytest.approx(5.17, abs=0.05)


def test_the_uncapped_flare_is_the_wider_one() -> None:
    """ADR-0011 before its amendment: a 19.8 T chamber expands its tubes further, to 6.50 m."""
    flown = field.fit_profile()
    assert fluxtube.flared_radius(23.0, flown) == pytest.approx(6.50, abs=0.1)


def test_the_cap_is_a_flat_shelf_that_meets_the_profile() -> None:
    """`min(B_profile, cap)`: flat inboard, on the profile outboard, continuous where they meet."""
    capped = field.capped_profile(12.0)
    assert capped.field(0.5) == pytest.approx(12.0)
    assert capped.field(1.0) == pytest.approx(12.0)
    assert capped.shelf_end_m == pytest.approx(3.12, abs=0.02)
    assert capped.field(capped.shelf_end_m) == pytest.approx(12.0, rel=1e-3)
    assert capped.field(23.8) < 5.0


# ---- What the answer rests on ------------------------------------------------------------------


def test_the_flared_winding_accommodates_every_tube(capped_stack: field.CoilStack) -> None:
    """ADR-0011's whole point: sized to the bounding tube, the geometric load is zero."""
    result, _tubes = fluxtube.evaluate(
        "flared-capped", capped_stack, launch="column", n_radial=6, n_axial=6
    )
    assert result.missing_fraction == pytest.approx(0.0)


def test_a_straight_winding_loses_the_tubes_born_at_the_strong_end() -> None:
    """R5's finding: a cylinder cannot hold a fourfold expansion, and the loss is not 2.3%."""
    stack, _profile = fluxtube.build_stacks(n_coils=32)["straight-flown"]
    result, _tubes = fluxtube.evaluate(
        "straight-flown", stack, launch="column", n_radial=6, n_axial=6
    )
    assert result.missing_fraction > 0.05


def test_the_chamber_launch_is_the_pessimistic_bound() -> None:
    """The two launch models bracket the answer, and the ordering is the argument in R2.

    A parcel born downstream has less field left to fall, so it fans less and is more likely to
    fit. Launching everything at the chamber must therefore lose at least as much.
    """
    stack, _profile = fluxtube.build_stacks(n_coils=32)["straight-flown"]
    column, _ = fluxtube.evaluate("s", stack, launch="column", n_radial=6, n_axial=6)
    chamber, _ = fluxtube.evaluate("s", stack, launch="chamber", n_radial=6)
    assert chamber.missing_fraction >= column.missing_fraction


def test_the_bore_is_axial_so_no_divergence_happens_inside_the_magnet(
    capped_stack: field.CoilStack,
) -> None:
    """R1's divergence term is a *downstream* quantity. If it were not, the answer would differ."""
    result, _tubes = fluxtube.evaluate(
        "flared-capped", capped_stack, launch="column", n_radial=6, n_axial=6
    )
    assert result.mean_cos_theta_exit > 0.99


def test_the_fan_opens_monotonically_past_the_last_coil(capped_stack: field.CoilStack) -> None:
    """Downstream of the winding the lines open, so `<cos theta>` must fall, never rise."""
    rows = [
        s
        for s in fluxtube.downstream_fan(capped_stack, r_exit=4.65, n_radial=6, out_to_radii=1.5)
        if not math.isnan(s.mean_cos_theta)
    ]
    assert rows[0].mean_cos_theta > rows[-1].mean_cos_theta
    assert rows[0].mean_cos_theta > 0.99
