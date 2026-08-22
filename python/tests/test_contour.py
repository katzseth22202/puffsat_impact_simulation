"""Tests for the continuous rho-ceiling contour (`contour.py`, Q15/Q13 of the 16-63 km/s extension).

The discrete construction evaluates 27 fixed cloud shapes at each velocity and keeps the best
survivor. That is fine as far as it goes, but the shape grid -- not the physics -- then sets where
the curve steps: at 45 -> 46 km/s the optimum jumps r_foot/R 0.5 -> 0.7 and `f` drops 0.817 ->
0.791, purely because no intermediate footprint exists to be chosen.

Q15 replaces that with the contour the schedule actually flies: `rho(v) = min(rho_ceiling(v),
rho_max)`, solved through the Sigma contract for the shape that hits it. The 27 discrete shapes
become validation points *on* the contour rather than the only places it may be sampled."""

from __future__ import annotations

import math
from itertools import pairwise

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
    """At high speed survivability binds, and the density it admits is **below** the plane-wave
    ceiling by exactly the focusing factor.

    `rho_ceiling` is the plane-wave limit `P_limit/(c_stag v^2)`. The facesheet actually sees that
    load concentrated by the dish (ADR-0010), so the flyable density is `rho_ceiling / focusing`.
    Before 2026-08-21 the contour used the plane-wave value directly and flew ~27% over its own
    limit."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    assert pt.ceiling_limited
    assert pt.focusing > 1.0, "the headline contour is concave, so it must concentrate"
    assert pt.rho_contour == pytest.approx(pt.rho_ceiling / pt.focusing, rel=1e-3)
    assert pt.rho_contour < pt.rho_ceiling
    # The solved shape must sit inside the box and match its own Sigma-contract density.
    assert impact_density(
        pt.l_over_d, pt.r_foot_over_r, contour.PULSE_MASS_KG, contour.PLATE_RADIUS_M
    ) == pytest.approx(pt.rho_contour, rel=1e-9)
    assert contour.L_OVER_D_BOX[0] <= pt.l_over_d <= contour.L_OVER_D_BOX[1]
    assert contour.R_FOOT_BOX[0] <= pt.r_foot_over_r <= contour.R_FOOT_BOX[1]


def test_contour_never_flies_over_the_pressure_limit_it_claims_to_respect() -> None:
    """The regression test for the 2026-08-21 focusing bug.

    The contour is *defined* as the densest cloud the facesheet survives, so the peak it actually
    flies must never exceed `p_limit`. The pre-fix construction compared the **plane-wave** peak
    against the limit while the facesheet sees the **focused** peak, so it sat over the limit at
    every velocity where survivability bound -- 509 MPa against a 400 MPa limit at 45 km/s. Nothing
    in the old test suite looked at the flown pressure, which is why it survived review."""
    for v in range(28_000, 64_000, 4_000):
        pt = contour.contour_point(float(v), c_stag=1.26, p_limit=4.0e8)
        assert pt.peak_pressure <= 4.0e8 * (1.0 + 1e-9), (
            f"contour at {v / 1000:.0f} km/s flies {pt.peak_pressure / 1e6:.0f} MPa "
            f"against a 400 MPa limit"
        )
        # Where survivability binds, the box could have built something denser and was forbidden.
        # (Checked against the box rather than the pressure, so it is not just restating the
        # active-constraint test that sets `ceiling_limited`.)
        if pt.ceiling_limited:
            assert pt.rho_contour < contour.rho_max_achievable()


def test_contour_agrees_with_the_discrete_frontier_it_replaced() -> None:
    """The cross-check that decided the focusing bug was in the contour, not in the fix.

    `heavyplate_frontier` evaluates the same physics on the 27-shape grid and has always applied
    the focusing factor. At 45 km/s its best surviving concave shape is `L/D = 0.3, r_foot/R = 0.5`
    -> `rho = 0.1258`, focusing 1.27, peak 399 MPa against the 400 MPa baseline. The continuous
    contour must land on that point, since the grid node is available to it."""
    pt = contour.contour_point(45_000.0, c_stag=1.2315, p_limit=4.0e8)

    assert pt.rho_contour == pytest.approx(0.1258, rel=0.02)
    assert pt.r_foot_over_r == pytest.approx(0.5, abs=0.02)
    assert pt.l_over_d == pytest.approx(0.3, abs=0.02)
    assert pt.focusing == pytest.approx(1.27, abs=0.02)


def test_contour_keeps_the_best_surviving_shape_in_the_box() -> None:
    """The objective, restated for the 2D search. Fixing `rho` first is no longer possible -- the
    survivable density depends on the shape through `focusing` -- so the contour sweeps the box and
    keeps the best *surviving* shape. Nothing that survives may beat it."""
    pt = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    for i in range(60):
        rf = contour.R_FOOT_BOX[0] + i * (contour.R_FOOT_BOX[1] - contour.R_FOOT_BOX[0]) / 59
        for j in range(60):
            lod = (
                contour.L_OVER_D_BOX[0]
                + j * (contour.L_OVER_D_BOX[1] - contour.L_OVER_D_BOX[0]) / 59
            )
            ok, rho, _, _ = contour.survivable(lod, rf, 63_000.0, 1.26, 4.0e8, pt.d_over_d)
            if ok:
                # No restitution source was supplied, so the objective is density.
                assert rho <= pt.rho_contour + 1e-6


def test_peak_and_score_fall_with_cloud_length() -> None:
    """The two monotonicities `contour_point`'s bisection rests on.

    At fixed footprint, lengthening the cloud (a) drops `rho` as `1/(L/D)` through the Sigma
    contract, which drops the facesheet peak, and (b) drops `eta_capture`, because a longer column
    splats with more radial relief. So the best surviving shape at each footprint is the *shortest*
    cloud that survives, and the search can bisect for it instead of scanning.

    If a future geometry sweep breaks either monotonicity the bisection silently returns the wrong
    shape, so both are pinned here rather than left as a comment."""
    lengths = (0.3, 0.45, 0.6, 0.8, 1.0)
    for rf in (0.35, 0.5, 0.65):
        peaks = [contour.survivable(ld, rf, 45_000.0, 1.2315, 4.0e8, 0.10)[3] for ld in lengths]
        etas = [contour.eta_at(ld, rf, 0.10) for ld in lengths]
        rhos = [contour.survivable(ld, rf, 45_000.0, 1.2315, 4.0e8, 0.10)[1] for ld in lengths]
        assert all(a > b for a, b in pairwise(peaks)), f"peak not monotone at rf={rf}"
        assert all(a > b for a, b in pairwise(etas)), f"eta not monotone at rf={rf}"
        assert all(a > b for a, b in pairwise(rhos)), f"rho not monotone at rf={rf}"


def test_contour_leaves_no_survivable_density_unused() -> None:
    """The bisection is what makes this exact rather than grid-resolution-limited.

    Where survivability binds, the contour must not sit *inside* the limit with room to spare --
    unused survivable density is `e_eff` thrown away. There are exactly two ways to be done: the
    cloud sits on the pressure limit, or it has already been shortened to the box's own `L/D` floor
    and cannot be made denser at that footprint. Anything else means the search stopped early."""
    for v in (28_000.0, 34_000.0, 45_000.0, 55_000.0, 63_000.0):
        pt = contour.contour_point(v, c_stag=1.2315, p_limit=4.0e8)
        assert pt.ceiling_limited, f"survivability should bind at {v / 1000:.0f} km/s"
        on_the_limit = pt.peak_pressure == pytest.approx(4.0e8, rel=1e-6)
        at_the_box_floor = pt.l_over_d == pytest.approx(contour.L_OVER_D_BOX[0], abs=1e-9)
        assert on_the_limit or at_the_box_floor, (
            f"{v / 1000:.0f} km/s: {pt.peak_pressure / 1e6:.1f} MPa at L/D={pt.l_over_d:.4f} "
            f"is neither on the limit nor at the box floor"
        )


def test_focusing_model_assumes_a_flat_plate_sees_the_plane_wave_load() -> None:
    """The hidden premise under `focusing_at`, and therefore under every survivability verdict.

    Survivability is classified as `peak = c_stag*rho*v^2 * focusing`, where the first factor is the
    **1D plane-wave** stagnation pressure (`c_stag` is measured from the 1D kernel) and `focusing`
    is the measured **2D** ratio `P_local(concave)/P_local(flat)`. Multiplying them only yields the
    true concave peak if the flat 2D plate sees the plane-wave load in the first place --
    `P_local(flat) == P_local(plane wave)`. Nothing had ever checked that.

    Seam: the free flat run against the confined (plane-wave) run of the **same kernel**, so scheme
    error is common-mode and cancels -- the same construction ADR-0003 uses for `eta_capture`.
    Tolerance: 10% on peak pressure, which stays well inside the SiC+Ti margin (the 400 MPa
    baseline is the conservative floor of a band running to 700/900 MPa).

    **Measured 2026-08-22: the premise holds, to 1.4%.** The flat plate runs `1.0011-1.0142` x the
    plane-wave peak across the box, rising with both `L/D` and `r_foot/R`. The direction is worth
    naming because it is the opposite of the obvious guess: a finite footprint does *not* relieve
    the local peak below plane-wave -- it sits slightly **above** it, so the survivability model is
    marginally **optimistic** rather than conservative. At <= 1.4% that is negligible against a
    400 MPa baseline drawn from a 400/700/900 band, but it is not zero and it is not the sign one
    would assume."""
    rows = [r for r in contour._geo() if abs(r["d_over_d"]) < 1e-9]
    assert rows, "no flat geometry rows to check the focusing premise against"

    worst = 0.0
    for r in rows:
        ratio = r["peak_local_pressure"] / r["peak_local_pressure_confined"]
        worst = max(worst, abs(ratio - 1.0))
        assert ratio == pytest.approx(1.0, rel=0.10), (
            f"flat plate at L/D={r['l_over_d']}, r_foot/R={r['r_foot_over_r']} sees "
            f"{ratio:.3f}x the plane-wave peak -- the focusing model's premise fails there"
        )
    # Drift guard, distinct from the 10% soundness gate above: the measured worst case is 1.4%, so
    # anything past 3% means the kernel or the sweep has moved and the premise wants re-reading,
    # even though the model would still be sound.
    assert worst < 0.03, f"flat-vs-plane-wave departure drifted to {worst * 100:.1f}% (was 1.4%)"


def test_the_plane_wave_peak_does_not_depend_on_the_cloud_shape() -> None:
    """A physics check that fell out of measuring the premise above, and is worth keeping.

    The confined run is the plane-wave limit, so its peak stagnation pressure is a function of the
    incident Mach number and `gamma` alone -- not of how long or how wide the slug is. The nine
    confined runs behind the flat geometry rows differ in both slab length and domain radius, so if
    the denominator of `focusing` carried any shape dependence it would show up here as a spread.

    It does not: all nine agree to a part in 10^4. That is what licenses treating `focusing` as pure
    geometry, and it is an independent check on the confined boundary condition."""
    peaks = [r["peak_local_pressure_confined"] for r in contour._geo() if abs(r["d_over_d"]) < 1e-9]
    assert len(peaks) >= 9
    assert max(peaks) - min(peaks) < 1e-4 * max(peaks)


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


def test_ablation_depth_is_fluence_over_ablator_enthalpy() -> None:
    """Q7's diagnostic: per-pulse radiative fluence at the wall converted to a sacrificial-layer
    recession depth, `depth = fluence / (rho_ablator * Q*)` (ADR-0014's quasi-steady Q* model).

    Worked by hand at 1.0e7 J/m^2, Q* = 5 MJ/kg, rho = 1200 kg/m^3:
    1.0e7/(1200*5.0e6) = 1.667e-3 m = 1.667 mm."""
    d = contour.ablation_depth(1.0e7, q_star=5.0e6, rho_ablator=1200.0)
    assert d == pytest.approx(1.6667e-3, rel=1e-4)
    # Inverse in Q*: the tougher ablator recedes less for the same fluence.
    lo = contour.ablation_depth(1.0e7, q_star=2.0e6, rho_ablator=1200.0)
    hi = contour.ablation_depth(1.0e7, q_star=10.0e6, rho_ablator=1200.0)
    assert lo / hi == pytest.approx(5.0, rel=1e-12)


def test_ablation_bracket_spans_the_literature_q_star_range() -> None:
    """ADR-0014 parameterises `Q*` over the silicone literature range 2-10 MJ/kg and requires the
    sensitivity be reported rather than a single value quoted. The bracket is therefore the pair,
    ordered deep-first: the *softest* ablator recedes most and is the conservative end."""
    band = contour.ablation_bracket(1.0e7)

    assert band.q_star_lo == 2.0e6
    assert band.q_star_hi == 10.0e6
    assert band.depth_max == pytest.approx(contour.ablation_depth(1.0e7, 2.0e6), rel=1e-12)
    assert band.depth_min == pytest.approx(contour.ablation_depth(1.0e7, 10.0e6), rel=1e-12)
    assert band.depth_max > band.depth_min


def test_f_band_combines_the_freeze_and_opacity_brackets() -> None:
    """`f` is reported as a band over both uncertainties (Q10/Q4), not one of them. They are
    independent -- freeze timing is when the composition stops equilibrating, opacity accuracy is
    how well TOPS knows kappa -- so they add in quadrature rather than linearly."""
    band = contour.f_band(f=0.80, freeze_delta=0.06, opacity_delta=0.005)

    assert band.f == pytest.approx(0.80)
    assert band.half_width == pytest.approx(math.hypot(0.06, 0.005), rel=1e-12)
    assert band.lo == pytest.approx(0.80 - math.hypot(0.06, 0.005), rel=1e-12)
    assert band.hi == pytest.approx(0.80 + math.hypot(0.06, 0.005), rel=1e-12)


def test_f_band_marks_an_unmeasured_freeze_bracket_rather_than_extrapolating() -> None:
    """The freeze bracket is measured at three anchors only (Q4: 16/22/28 km/s), because it tracks
    an ionization staircase rather than a ramp and cannot be interpolated. Above 28 km/s the band
    must therefore be *labelled* as carried forward, not silently extrapolated -- a reader has to
    be able to tell a measured bracket from an assumed one."""
    measured = contour.f_band(f=0.80, freeze_delta=0.06, opacity_delta=0.005, freeze_measured=True)
    carried = contour.f_band(f=0.78, freeze_delta=0.081, opacity_delta=0.005, freeze_measured=False)

    assert measured.freeze_measured is True
    assert carried.freeze_measured is False
    # The carried-forward width is the widest measured one, so it cannot understate the band.
    assert carried.half_width > measured.half_width


def test_recession_scales_with_the_column_the_pulse_delivers() -> None:
    """A scaling check that makes the magnitudes auditable rather than merely plausible.

    Fluence tracks a roughly fixed fraction of incident kinetic energy per unit area, which is
    `0.5*rho*L*v^2`. So at fixed `rho` and `v`, ten times the cloud length delivers ten times the
    fluence and ten times the recession -- which is exactly why the heavy-plate scenario (L = 10 m)
    reports ~12x the core envelope study's per-pulse recession at the same 16 km/s, and why
    ADR-0014's few-micron figure (derived for neither) does not describe either of them."""
    shallow = contour.ablation_depth(1.0e6, q_star=5.0e6)
    deep = contour.ablation_depth(1.0e7, q_star=5.0e6)

    assert deep / shallow == pytest.approx(10.0, rel=1e-12)
    # And the bracket keeps that proportionality at both ends.
    assert contour.ablation_bracket(1.0e7).depth_max / contour.ablation_bracket(
        1.0e6
    ).depth_max == pytest.approx(10.0, rel=1e-12)


def test_pinned_shape_is_the_one_that_survives_the_worst_case() -> None:
    """Q13's companion curve. The headline lets the cloud shape float with velocity -- design SS7's
    schedule `shape(v)`, one shot at a time. The pinned curve answers the opposite question: what
    if you cannot vary it, and must fly one shape everywhere?

    Then the shape is forced by the *most demanding* velocity in the range, because it has to
    survive there. Pinning anywhere else would produce a curve that fails at the top end, which is
    not a curve anyone can fly."""
    pin = contour.pinned_shape(63_000.0, c_stag=1.26, p_limit=4.0e8)
    free = contour.contour_point(63_000.0, c_stag=1.26, p_limit=4.0e8)

    # At the pinning velocity the two coincide -- that is where the constraint binds.
    assert pin.rho_contour == pytest.approx(free.rho_contour, rel=1e-12)
    assert pin.r_foot_over_r == pytest.approx(free.r_foot_over_r, rel=1e-12)


def test_pinned_curve_is_dilute_and_therefore_worse_at_low_speed() -> None:
    """The cost of not scheduling. A shape sized to survive 63 km/s is far more dilute than 16 km/s
    needs (0.081 against a 0.582 contour), and a dilute cloud radiates away more of the bounce --
    so the pinned curve must sit at or below the floating one everywhere, with the gap widest where
    the two densities differ most."""
    e_of_rho = {0.0806: 0.36, 0.5822: 0.64}

    def e_eff_at(rho: float) -> float:
        return e_of_rho[min(e_of_rho, key=lambda k: abs(k - rho))]

    free_16 = contour.contour_point(16_000.0, c_stag=1.26, p_limit=4.0e8, e_eff_at=e_eff_at)
    pinned_16 = contour.pinned_shape(63_000.0, c_stag=1.26, p_limit=4.0e8, e_eff_at=e_eff_at)

    assert pinned_16.rho_contour < free_16.rho_contour
    assert pinned_16.f is not None and free_16.f is not None
    assert pinned_16.f < free_16.f
