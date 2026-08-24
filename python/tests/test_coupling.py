"""Study 2: does a projectile couple to a droplet cloud the way it couples to a vapour?

The routing document asks this because if the field leak turns out small -- and it did, 0.11 to
2.5% -- the bag need not be pressurised, and 213 kg over 660 m^3 could be a cloud of droplets
rather than a vapour. What is at stake is `k = 8.5`.

Answers computed by hand before the module existed: a 10 um droplet couples in ~1.7 us against a
~2.3 ms transit, and the limiting size is ~1.4 cm.
"""

from __future__ import annotations

import math

import pytest

from puffsat import coupling

TRANSIT = 2.3e-3  # the paper's stated snowplow transit [s]
COLD = 45.58e3


def test_drag_time_is_linear_in_radius_and_inverse_in_relative_speed() -> None:
    """`tau_drag = (8/3) a rho_w / (C_d rho_g dv)` -- the particle relaxation time at high Re.

    Both scalings are the whole reason the answer is comfortable: droplets are small and the
    relative speed is enormous, and the time goes down with both.
    """
    base = coupling.drag_time(1.0e-5, coupling.BAG_RHO, COLD)
    assert base == pytest.approx(1.7e-6, rel=0.1), "10 um droplet, hand value"

    assert coupling.drag_time(1.0e-4, coupling.BAG_RHO, COLD) == pytest.approx(10.0 * base)
    assert coupling.drag_time(1.0e-5, coupling.BAG_RHO, 2.0 * COLD) == pytest.approx(base / 2.0)
    assert coupling.drag_time(1.0e-5, 2.0 * coupling.BAG_RHO, COLD) == pytest.approx(base / 2.0)


def test_vaporisation_is_far_faster_than_drag_so_drag_is_what_binds() -> None:
    """A droplet at 45 km/s is not smoothly dragged -- it is flash-vaporised. Both routes end in
    coupling, so the *slower* one is the criterion, and it is drag by ~3 orders.

    This is worth asserting rather than assuming: if vaporisation were the slower channel the
    limiting size would be set by latent heat instead of by inertia, and the scaling with relative
    speed would be cubic rather than linear.
    """
    a = 1.0e-5
    drag = coupling.drag_time(a, coupling.BAG_RHO, COLD)
    vapour = coupling.vaporisation_time(a, coupling.BAG_RHO, COLD)
    assert vapour < drag / 100.0
    assert coupling.coupling_time(a, coupling.BAG_RHO, COLD) == pytest.approx(drag)


def test_the_limiting_droplet_is_centimetres_which_no_bag_would_produce() -> None:
    """The answer to Study 2 as asked. Coupling fails when a droplet cannot be brought up to
    speed within the transit, and that takes a **1.4 cm** droplet.

    The margin depends on how the cloud is made, and it is worth stating per case rather than as
    one number: condensation from vapour makes sub-micron to micron drops (**4 orders** of margin),
    an injected spray makes tens to hundreds of microns (**2 orders**). Both are comfortable, and
    2 is the honest figure for the case a designer would actually build.
    """
    limit = coupling.limiting_radius(TRANSIT, coupling.BAG_RHO, COLD)
    assert limit == pytest.approx(1.4e-2, rel=0.15)

    condensation = coupling.coupling_number(1.0e-6, TRANSIT, coupling.BAG_RHO, COLD)
    spray = coupling.coupling_number(1.0e-4, TRANSIT, coupling.BAG_RHO, COLD)
    assert condensation > 1.0e4
    assert 1.0e2 < spray < 1.0e3, "a spray clears it by two orders, not four"


def test_discreteness_never_binds_at_bag_scale() -> None:
    """The other way a droplet cloud could differ from a vapour: the front passing *between* the
    drops rather than through them. It cannot -- even centimetre drops are intercepted in millions.
    """
    for radius in (1.0e-5, 1.0e-3, 1.0e-2):
        count = coupling.interception_count(radius, coupling.BAG_RHO, coupling.BORE_RADIUS, 23.8)
        assert count > 1.0e3, f"{radius} m droplets must still be a continuum to the front"


def test_the_full_bore_snowplow_reproduces_the_papers_slug_ratio() -> None:
    """`k = rho A L / m` is the paper's own arithmetic, and the module must reproduce it exactly
    when the front is given the bore radius from the start. That is the calibration.
    """
    assert coupling.full_bore_slug_ratio() == pytest.approx(8.69, rel=0.01)

    at_bore = coupling.snowplow(COLD, front_radius=coupling.BORE_RADIUS, expansion_speed=0.0)
    assert at_bore.slug_ratio == pytest.approx(coupling.full_bore_slug_ratio(), rel=0.01)
    assert at_bore.transit_time == pytest.approx(2.79e-3, rel=0.05), "L(1 + k/2)/v_0"


def test_a_compact_projectile_cannot_sweep_the_bore_and_that_is_an_input_not_an_output() -> None:
    """**The finding that fell out of Study 2, and it is about the vapour case too.**

    `lambda = rho A` assumes the sweeping front spans the full bore from x = 0. A compact 25 kg
    ice projectile is 0.187 m against a 3 m bore -- 16x in radius, 260x in area -- and a front
    that has to grow at any plausible lateral speed does not get there in time. So `k = 8.5` is
    a statement about the projectile's *arrival radius*, which this repository does not own.
    """
    compact = coupling.projectile_radius(coupling.PROJECTILE_MASS)
    assert compact == pytest.approx(0.187, rel=0.02)

    grown = coupling.snowplow(COLD, front_radius=compact, expansion_speed=5.0e3)
    assert grown.slug_ratio < 0.6 * coupling.full_bore_slug_ratio(), "falls short by well over 2x"

    # Monotone in the lateral speed, so it is a bracket rather than a single verdict.
    faster = coupling.snowplow(COLD, front_radius=compact, expansion_speed=8.0e3)
    assert faster.slug_ratio > grown.slug_ratio
    assert faster.slug_ratio < coupling.full_bore_slug_ratio()


def test_the_shortfall_is_worse_at_higher_closing_speed() -> None:
    """The direction matters for the mission: the hot legs are where the front has least time.

    A faster projectile crosses the bore in less time, so its front grows less before it exits --
    `dr/dx = c/v`. So if this bites at all it bites hardest exactly where the paper needs `k` most.
    """
    compact = coupling.projectile_radius(coupling.PROJECTILE_MASS)
    cold = coupling.snowplow(COLD, front_radius=compact, expansion_speed=5.0e3)
    hot = coupling.snowplow(75.0e3, front_radius=compact, expansion_speed=5.0e3)
    assert hot.slug_ratio < cold.slug_ratio
    assert math.isfinite(hot.slug_ratio)


def test_the_requirement_is_stated_as_an_arrival_radius_rather_than_a_shortfall() -> None:
    """The constructive form: `k = 8.5` demands the projectile arrive spanning most of the bore.

    Between 74% and 97% of the bore radius depending on how fast the front spreads -- 2.2 to 2.9 m
    of a 3 m bore. That is a checkable requirement on `aim`'s projectile design, where "the
    snowplow overstates k" is only a complaint.

    Design SS7's own footprint box for the plate-side pulse is `r_foot/R` = 0.3-1.0 at `R` = 5 m,
    i.e. 1.5-5 m, so the upper half of that box clears this. The requirement is plausible; it is
    simply unstated.
    """
    for c_exp in (0.0, 3.0e3, 5.0e3):
        need = coupling.required_arrival_radius(COLD, c_exp)
        assert 0.70 <= need / coupling.BORE_RADIUS <= 1.0, f"c_exp = {c_exp}"

    # A faster-spreading front relaxes the requirement, and a still front makes it hardest.
    assert coupling.required_arrival_radius(COLD, 5.0e3) < coupling.required_arrival_radius(
        COLD, 0.0
    )

    # Verifying it: arriving at the required radius really does deliver the target.
    need = coupling.required_arrival_radius(COLD, 5.0e3)
    got = coupling.snowplow(COLD, need, 5.0e3).slug_ratio
    assert got == pytest.approx(0.95 * coupling.full_bore_slug_ratio(), rel=0.01)
