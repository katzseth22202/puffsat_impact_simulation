"""Tests for the continuous rho-ceiling contour (`contour.py`, Q15/Q13 of the 16-63 km/s extension).

The discrete construction evaluates 27 fixed cloud shapes at each velocity and keeps the best
survivor. That is fine as far as it goes, but the shape grid -- not the physics -- then sets where
the curve steps: at 45 -> 46 km/s the optimum jumps r_foot/R 0.5 -> 0.7 and `f` drops 0.817 ->
0.791, purely because no intermediate footprint exists to be chosen.

Q15 replaces that with the contour the schedule actually flies: `rho(v) = min(rho_ceiling(v),
rho_max)`, solved through the Sigma contract for the shape that hits it. The 27 discrete shapes
become validation points *on* the contour rather than the only places it may be sampled."""

from __future__ import annotations

import pytest

from puffsat import contour
from puffsat.analysis import impact_density


def test_sigma_contract_inverts_exactly() -> None:
    """The contour is built by solving the Sigma contract backwards, so the inversion has to be an
    exact inverse of `analysis.impact_density`, not merely close. Round-tripping is the check that
    matters: any drift here silently moves every contour point off the density it claims."""
    for r_foot_over_r in (0.3, 0.42, 0.5, 0.618, 0.7):
        for rho in (0.02, 0.08, 0.3, 0.58):
            l_over_d = contour.l_over_d_for(rho, r_foot_over_r)
            assert impact_density(
                l_over_d, r_foot_over_r, contour.PULSE_MASS_KG, contour.PLATE_RADIUS_M
            ) == pytest.approx(rho, rel=1e-12)


def test_rho_ceiling_follows_the_stagnation_law() -> None:
    """`rho_ceiling = P_limit/(c_stag v^2)` (design SS7 / ADR-0010). Worked by hand at 400 MPa,
    c_stag = 1.26, v = 45 km/s: 4.0e8/(1.26*2.025e9) = 0.15681 kg/m^3."""
    assert contour.rho_ceiling(45_000.0, c_stag=1.26, p_limit=4.0e8) == pytest.approx(
        0.156_770, rel=1e-5
    )
    # Inverse-square in v: doubling the speed quarters the survivable density.
    lo = contour.rho_ceiling(30_000.0, c_stag=1.26, p_limit=4.0e8)
    hi = contour.rho_ceiling(60_000.0, c_stag=1.26, p_limit=4.0e8)
    assert lo / hi == pytest.approx(4.0, rel=1e-12)


def test_contour_takes_the_shape_box_limit_when_survivability_is_slack() -> None:
    """At low speed the plate survives far denser clouds than the shape box can deliver, so the
    contour is limited by geometry, not survivability -- and the reported point must say which,
    because "as dense as we can build" and "as dense as it can take" carry different engineering
    consequences.

    At 16 km/s the ceiling is ~1.24 kg/m^3 while the densest shape in the box (L/D = 0.3,
    r_foot/R = 0.3) gives ~0.58, so the box binds."""
    pt = contour.contour_point(16_000.0, c_stag=1.26, p_limit=4.0e8)

    assert pt.rho_ceiling > pt.rho_contour
    assert not pt.ceiling_limited
    assert pt.rho_contour == pytest.approx(
        impact_density(0.3, 0.3, contour.PULSE_MASS_KG, contour.PLATE_RADIUS_M), rel=1e-9
    )


def test_contour_takes_the_survivability_limit_when_the_box_is_slack() -> None:
    """At high speed survivability binds: 63 km/s admits only ~0.080 kg/m^3, well inside what the
    box can deliver, so the cloud must be stretched to meet it."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    assert pt.ceiling_limited
    assert pt.rho_contour == pytest.approx(pt.rho_ceiling, rel=1e-12)
    assert pt.rho_contour == pytest.approx(0.080_04, rel=1e-3)
    # The solved shape must actually sit on the contour and inside the box.
    assert impact_density(
        pt.l_over_d, pt.r_foot_over_r, contour.PULSE_MASS_KG, contour.PLATE_RADIUS_M
    ) == pytest.approx(pt.rho_contour, rel=1e-9)
    assert contour.L_OVER_D_BOX[0] <= pt.l_over_d <= contour.L_OVER_D_BOX[1]
    assert contour.R_FOOT_BOX[0] <= pt.r_foot_over_r <= contour.R_FOOT_BOX[1]


def test_contour_maximizes_eta_along_the_iso_density_curve() -> None:
    """One equation (fixed rho) in two shape parameters leaves a one-parameter family, so the
    contour point is *chosen*, not determined -- it is the member of that family with the highest
    `eta_capture`, since `e_eff` is already pinned by rho.

    The chosen point must beat every other feasible point on the same curve."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    feasible = []
    for i in range(200):
        rf = contour.R_FOOT_BOX[0] + i * (contour.R_FOOT_BOX[1] - contour.R_FOOT_BOX[0]) / 199
        lod = contour.l_over_d_for(pt.rho_contour, rf)
        if contour.L_OVER_D_BOX[0] <= lod <= contour.L_OVER_D_BOX[1]:
            feasible.append(contour.eta_at(lod, rf, pt.d_over_d))
    assert feasible, "the iso-density curve must intersect the shape box"
    assert pt.eta_capture >= max(feasible) - 1e-9


def test_eta_interpolation_reproduces_the_geometry_grid_nodes() -> None:
    """The interpolation must return the swept value at a swept point, or the contour is reporting
    an `eta_capture` the 2D track never computed."""
    assert contour.eta_at(0.3, 0.5, 0.10) == pytest.approx(0.9783, abs=5e-4)
    assert contour.eta_at(1.0, 0.7, 0.0) == pytest.approx(0.7924, abs=5e-4)
    assert contour.eta_at(0.3, 0.3, 0.15) == pytest.approx(0.9920, abs=5e-4)


def test_contour_is_continuous_where_the_discrete_grid_stepped() -> None:
    """The point of Q15. Across 45 -> 46 km/s the discrete construction dropped `f` 0.817 -> 0.791
    (a 0.026 step) because the shape optimum jumped r_foot/R 0.5 -> 0.7 with nothing in between.
    The contour has intermediate footprints available, so the same crossing must be smooth."""
    pts = [contour.contour_point(v, c_stag=1.26, p_limit=4.0e8) for v in (45_000.0, 46_000.0)]

    d_rf = abs(pts[1].r_foot_over_r - pts[0].r_foot_over_r)
    d_eta = abs(pts[1].eta_capture - pts[0].eta_capture)
    assert d_rf < 0.05, f"footprint should drift, not jump: {d_rf:.4f}"
    assert d_eta < 0.01, f"eta should drift, not jump: {d_eta:.4f}"


def test_eta_validity_is_asserted_only_along_the_contour() -> None:
    """A caveat with teeth. `eta_capture` spans 0.79-0.99 across the full shape box, so an
    interpolation over the whole box would be advertising a precision it does not have. Along the
    schedule the contour actually flies, it stays in a narrow band -- and that is the only claim
    made."""
    etas = [
        contour.contour_point(v, c_stag=1.26, p_limit=4.0e8).eta_capture
        for v in range(16_000, 64_000, 1_000)
    ]
    assert min(etas) > 0.93, f"contour eta dropped to {min(etas):.4f}"
    assert max(etas) - min(etas) < 0.06, "contour eta band should stay narrow"
    # ... while the full box genuinely is that wide, which is why the caveat exists.
    assert contour.ETA_BOX_RANGE[0] < 0.80 and contour.ETA_BOX_RANGE[1] > 0.98


def test_scenario_constants_match_the_heavy_plate_module() -> None:
    """The contour must fly the heavy plate, not the core envelope. `analysis` defines 25 kg / 5 m
    for the core study and `heavyplate` overrides them to 100 kg / 15 m; importing the wrong pair
    moves every contour density by ~7x while still producing a plausible-looking curve. This is
    what caught exactly that mistake."""
    from puffsat import heavyplate

    assert contour.PULSE_MASS_KG == heavyplate.PULSE_MASS_KG
    assert contour.PLATE_RADIUS_M == heavyplate.PLATE_RADIUS_M


def test_f_on_the_contour_reconciles_eta_and_e_eff() -> None:
    """`f = eta*(1 + e_eff)/2` (ADR-0001), with `e_eff` read at the contour density. The restitution
    source is injected rather than imported so this module does not depend on the analysis that
    consumes it."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8, e_eff_at=lambda rho: 0.60)

    assert pt.e_eff == pytest.approx(0.60)
    assert pt.f == pytest.approx(pt.eta_capture * 1.60 / 2.0, rel=1e-12)


def test_f_is_absent_rather_than_guessed_without_a_restitution_source() -> None:
    """No `e_eff` source means no `f`. Reporting a default would put a fabricated number in the
    deliverable's own column."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    assert pt.e_eff is None
    assert pt.f is None
