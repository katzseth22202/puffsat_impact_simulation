"""Tests for the pulse-shape sensitivity assembly (design §13, ADR-0028): the S metric, the
cliff detector, the noise floor, the Σ-profile bound, and the f assembly on synthetic rows."""

import math

import pytest

from puffsat.analysis import reconcile_f
from puffsat.shape import (
    BOX_MAX,
    NOM_L_OVER_D,
    NOM_RFOOT_OVER_R,
    FrozenShapeRow,
    Shape1DRow,
    Shape2DRow,
    assemble,
    cliff_flags,
    frozen_slope_check,
    max_abs_s_per_axis,
    noise_floor,
    refinement_verdict,
    resolve_1d,
    sigma_profile_bound,
)

V = 11_000.0


def _row2d(
    axis: str,
    *,
    rff: float = NOM_RFOOT_OVER_R,
    lod: float = NOM_L_OVER_D,
    taper: float = 0.0,
    alpha: float = 0.0,
    d_over_d: float = 0.1,
    eta: float = 0.98,
    res: float = 1.0,
    peak_local: float = 500.0,
) -> Shape2DRow:
    return Shape2DRow(
        axis=axis,
        d_over_d=d_over_d,
        r_foot_over_r=rff,
        l_over_d=lod,
        taper_frac=taper,
        alpha_div=alpha,
        mach=20.0,
        resolution_scale=res,
        eta_capture=eta,
        restitution_free=1.4,
        restitution_confined=1.43,
        incident_momentum=0.6,
        peak_local_pressure=peak_local,
    )


def _row1d(
    axis: str,
    *,
    rff: float = NOM_RFOOT_OVER_R,
    lod: float = NOM_L_OVER_D,
    taper: float = 0.0,
    role: str = "sample",
    e_eff: float = 0.57,
    e_eff_coarse: float | None = None,
    e_eff_eos: float | None = None,
    peak: float = 1.2e8,
) -> Shape1DRow:
    return Shape1DRow(
        v=V,
        axis=axis,
        sigma_role=role,
        r_foot_over_r=rff,
        l_over_d=lod,
        taper_frac=taper,
        sigma=1.27,
        rho_impact=0.85,
        length=1.5,
        e_eff=e_eff,
        e_eff_coarse=e_eff if e_eff_coarse is None else e_eff_coarse,
        e_eff_eos=(e_eff + 0.006) if e_eff_eos is None else e_eff_eos,
        peak_wall_pressure=peak,
        peak_wall_pressure_coarse=peak * 0.99,
        incident_momentum=1.4e4,
        wall_impulse=2.2e4,
    )


def _synthetic_box(
    *, eta_slope: float = 0.0, e_slope: float = 0.0
) -> tuple[list[Shape2DRow], list[Shape1DRow]]:
    """A minimal synthetic box: nominal + a footprint axis (two-sided) + one alpha sample, with
    optional linear responses of eta (2D) and e_eff (1D) in the footprint's relative offset.
    Both plates (d/D = 0.1 and the flat 0) are emitted so focusing is computable."""
    rels = (0.8, 0.9, 1.1, 1.2)
    rows2d: list[Shape2DRow] = []
    rows1d: list[Shape1DRow] = [_row1d("nominal")]
    for dd in (0.1, 0.0):
        rows2d.append(_row2d("nominal", d_over_d=dd))
        rows2d.append(_row2d("alpha_div", alpha=0.1, d_over_d=dd))
        for rel in rels:
            rows2d.append(
                _row2d(
                    "r_foot_over_r",
                    rff=NOM_RFOOT_OVER_R * rel,
                    d_over_d=dd,
                    eta=0.98 * (1.0 + eta_slope * (rel - 1.0)),
                )
            )
    for rel in rels:
        rows1d.append(
            _row1d(
                "r_foot_over_r",
                rff=NOM_RFOOT_OVER_R * rel,
                e_eff=0.57 + e_slope * (rel - 1.0),
            )
        )
    return rows2d, rows1d


class TestAssembly:
    def test_f_is_the_adr0003_reconciliation(self) -> None:
        rows2d, rows1d = _synthetic_box()
        points = assemble(rows2d, rows1d, V, 0.1)
        nominal = next(p for p in points if p.axis == "nominal")
        assert nominal.f == reconcile_f(0.98, 0.57)

    def test_alpha_axis_reuses_the_nominal_e_eff(self) -> None:
        rows2d, rows1d = _synthetic_box()
        points = assemble(rows2d, rows1d, V, 0.1)
        alpha = next(p for p in points if p.axis == "alpha_div")
        nominal = next(p for p in points if p.axis == "nominal")
        assert alpha.e_eff == nominal.e_eff

    def test_flat_plate_has_unit_focusing(self) -> None:
        rows2d, rows1d = _synthetic_box()
        flat = assemble(rows2d, rows1d, V, 0.0)
        nominal = next(p for p in flat if p.axis == "nominal")
        assert nominal.peak_compressive == 1.2e8  # the 1D physical peak, unscaled

    def test_concave_focusing_scales_the_1d_peak(self) -> None:
        rows2d, rows1d = _synthetic_box()
        # Make the dished nominal focus 1.5x its flat counterpart.
        rows2d = [
            (
                r
                if not (r.axis == "nominal" and r.d_over_d > 0.0)
                else Shape2DRow(**{**r.__dict__, "peak_local_pressure": 750.0})
            )
            for r in rows2d
        ]
        points = assemble(rows2d, rows1d, V, 0.1)
        nominal = next(p for p in points if p.axis == "nominal")
        assert math.isclose(nominal.peak_compressive, 1.5 * 1.2e8)


class TestSensitivity:
    def test_linear_response_gives_constant_s(self) -> None:
        # f moves only through eta, linearly with slope 0.5 in Δx/x ⇒ S = 0.5 everywhere.
        rows2d, rows1d = _synthetic_box(eta_slope=0.5)
        points = assemble(rows2d, rows1d, V, 0.1)
        ss = [p.s for p in points if p.axis == "r_foot_over_r"]
        assert ss and all(math.isclose(s, 0.5, rel_tol=1e-9) for s in ss)

    def test_one_sided_axis_normalizes_by_box_extent(self) -> None:
        rows2d, rows1d = _synthetic_box()
        # Raise the alpha sample's eta by 2% at the full box extent ⇒ S = 0.02/1.0 = 0.02.
        rows2d = [
            (
                r
                if not (r.axis == "alpha_div" and r.d_over_d > 0.0)
                else Shape2DRow(**{**r.__dict__, "eta_capture": 0.98 * 1.02})
            )
            for r in rows2d
        ]
        points = assemble(rows2d, rows1d, V, 0.1)
        alpha = next(p for p in points if p.axis == "alpha_div")
        assert alpha.rel_delta == 0.1 / BOX_MAX["alpha_div"]
        assert math.isclose(alpha.s, 0.02, rel_tol=1e-9)

    def test_max_abs_s_per_axis(self) -> None:
        rows2d, rows1d = _synthetic_box(eta_slope=0.5)
        points = assemble(rows2d, rows1d, V, 0.1)
        s_axis = max_abs_s_per_axis(points)
        assert math.isclose(s_axis["r_foot_over_r"], 0.5, rel_tol=1e-9)


class TestCliffDetector:
    def test_linear_curve_has_no_flags(self) -> None:
        rows2d, rows1d = _synthetic_box(eta_slope=0.5)
        points = assemble(rows2d, rows1d, V, 0.1)
        assert not any(p.cliff for p in points)

    def test_step_is_flagged_above_noise(self) -> None:
        rows2d, rows1d = _synthetic_box(eta_slope=0.05)
        # Inject a step: every footprint sample past +10% jumps by 5% in eta on both plates
        # (the step lands on the segment ending at x = 0.55).
        rows2d = [
            (
                r
                if not (r.axis == "r_foot_over_r" and r.r_foot_over_r > 0.54)
                else Shape2DRow(**{**r.__dict__, "eta_capture": r.eta_capture * 1.05})
            )
            for r in rows2d
        ]
        points = assemble(rows2d, rows1d, V, 0.1)
        flagged = [p for p in points if p.cliff]
        assert len(flagged) == 1
        assert flagged[0].axis == "r_foot_over_r"
        assert math.isclose(flagged[0].x, 0.55)

    def test_structure_below_noise_floor_is_a_pass(self) -> None:
        # The same step, but with a noise floor dominating it: no flags (design §13 —
        # below-noise structure is indistinguishable from flat).
        curve_points = assemble(*_synthetic_box(eta_slope=0.05), V, 0.1)
        curve = sorted((p for p in curve_points if p.axis == "r_foot_over_r"), key=lambda p: p.x)
        assert cliff_flags(curve, sigma_noise=1.0, f_ref=0.77) == [False] * (len(curve) - 1)


class TestNoiseFloor:
    def test_noise_floor_is_max_pstdev_over_repeats(self) -> None:
        rows2d, rows1d = _synthetic_box()
        # Nominal repeated at 3 resolutions with eta spread 0.98 ± 0.005 (dished plate only).
        for res, eta in ((1.5, 0.985), (2.0, 0.975)):
            rows2d.append(_row2d("nominal", d_over_d=0.1, res=res, eta=eta))
        sigma = noise_floor(rows2d, rows1d, V, 0.1)
        fs = [reconcile_f(eta, 0.57) for eta in (0.98, 0.985, 0.975)]
        mean = sum(fs) / 3.0
        expected = math.sqrt(sum((f - mean) ** 2 for f in fs) / 3.0)
        assert math.isclose(sigma, expected)

    def test_no_repeats_means_zero_floor(self) -> None:
        rows2d, rows1d = _synthetic_box()
        assert noise_floor(rows2d, rows1d, V, 0.1) == 0.0


class TestSigmaBound:
    def test_bound_is_the_worst_endpoint_spread(self) -> None:
        rows2d, rows1d = _synthetic_box()
        widest = BOX_MAX["taper_frac"]
        rows2d.append(_row2d("taper_frac", taper=widest, d_over_d=0.1, eta=1.0))
        rows2d.append(_row2d("taper_frac", taper=widest, d_over_d=0.0))
        rows1d += [
            _row1d("taper_frac", taper=widest, role="taper_mean", e_eff=0.570),
            _row1d("taper_frac", taper=widest, role="sigma_hi", e_eff=0.560),
            _row1d("taper_frac", taper=widest, role="sigma_lo90", e_eff=0.585),
        ]
        points = assemble(rows2d, rows1d, V, 0.1)
        delta_e, delta_f = sigma_profile_bound(rows1d, points, V)
        assert math.isclose(delta_e, 0.015)
        assert math.isclose(delta_f, 1.0 * 0.015 / 2.0)


def _stepped_box() -> tuple[list[Shape2DRow], list[Shape1DRow]]:
    """The step-injected box of the cliff-detector test: eta jumps 5% for every footprint sample
    past +10%, flagging the segment ending at x = 0.55."""
    rows2d, rows1d = _synthetic_box(eta_slope=0.05)
    rows2d = [
        (
            r
            if not (r.axis == "r_foot_over_r" and r.r_foot_over_r > 0.54)
            else Shape2DRow(**{**r.__dict__, "eta_capture": r.eta_capture * 1.05})
        )
        for r in rows2d
    ]
    return rows2d, rows1d


class TestRefinementVerdict:
    def test_unrefined_flag_is_unresolved(self) -> None:
        rows2d, rows1d = _stepped_box()
        points = assemble(rows2d, rows1d, V, 0.1)
        verdicts = refinement_verdict(rows2d, rows1d, points, V, 0.1)
        assert verdicts == [("r_foot_over_r", 0.55, None)]

    def test_flag_that_vanishes_at_refined_resolution_is_grid_noise(self) -> None:
        rows2d, rows1d = _stepped_box()
        # Refined repeats put the sample back on the eta trend: the flag must not survive.
        on_trend = 0.98 * (1.0 + 0.05 * 0.1)
        for res in (1.5, 2.0):
            rows2d.append(_row2d("r_foot_over_r", rff=0.55, res=res, eta=on_trend))
        points = assemble(rows2d, rows1d, V, 0.1)
        verdicts = refinement_verdict(rows2d, rows1d, points, V, 0.1)
        assert verdicts == [("r_foot_over_r", 0.55, False)]

    def test_flag_that_persists_at_refined_resolution_is_physical(self) -> None:
        rows2d, rows1d = _stepped_box()
        jumped = 0.98 * (1.0 + 0.05 * 0.1) * 1.05
        for res in (1.5, 2.0):
            rows2d.append(_row2d("r_foot_over_r", rff=0.55, res=res, eta=jumped))
        points = assemble(rows2d, rows1d, V, 0.1)
        verdicts = refinement_verdict(rows2d, rows1d, points, V, 0.1)
        assert verdicts == [("r_foot_over_r", 0.55, True)]


class TestSolverValidityProtocol:
    def test_agreeing_fine_run_is_the_headline(self) -> None:
        row = _row1d("nominal", e_eff=0.638, e_eff_coarse=0.639)
        e_eff, peak, valid = resolve_1d(row)
        assert (e_eff, peak, valid) == (0.638, 1.2e8, True)

    def test_collapsed_fine_run_falls_back_to_the_stable_window(self) -> None:
        # The 16 km/s radiative-collapse signature: unphysical fine e_eff, healthy coarse one.
        row = _row1d("r_foot_over_r", rff=0.4, e_eff=-0.317, e_eff_coarse=0.640, e_eff_eos=0.646)
        e_eff, peak, valid = resolve_1d(row)
        assert e_eff == 0.640
        assert math.isclose(peak, 1.2e8 * 0.99)
        assert not valid

    def test_both_resolutions_bad_is_a_hard_error(self) -> None:
        row = _row1d("r_foot_over_r", rff=0.4, e_eff=-0.3, e_eff_coarse=0.1, e_eff_eos=0.646)
        with pytest.raises(ValueError, match="outside validity"):
            resolve_1d(row)

    def test_fallback_is_flagged_in_the_assembly(self) -> None:
        rows2d, rows1d = _synthetic_box()
        rows1d = [
            (
                r
                if r.axis != "r_foot_over_r" or not math.isclose(r.r_foot_over_r, 0.4)
                else _row1d(
                    "r_foot_over_r", rff=0.4, e_eff=-0.317, e_eff_coarse=0.57, e_eff_eos=0.576
                )
            )
            for r in rows1d
        ]
        points = assemble(rows2d, rows1d, V, 0.1)
        flagged = [p for p in points if not p.solver_valid]
        assert [(p.axis, p.x) for p in flagged] == [("r_foot_over_r", 0.4)]


class TestFrozenSpotCheck:
    def test_comparable_slopes_pass(self) -> None:
        rows = [
            FrozenShapeRow(V, 0.59, 0.60, 0.52, 0.63, 1e-4),
            FrozenShapeRow(V, 0.85, 0.57, 0.49, 0.60, 1e-4),
            FrozenShapeRow(V, 1.33, 0.55, 0.47, 0.58, 1e-4),
        ]
        d_eq, d_frozen, ok = frozen_slope_check(rows)
        assert math.isclose(d_eq, -0.05) and math.isclose(d_frozen, -0.05)
        assert ok

    def test_wild_frozen_slope_fails(self) -> None:
        rows = [
            FrozenShapeRow(V, 0.59, 0.60, 0.60, 0.63, 1e-4),
            FrozenShapeRow(V, 1.33, 0.55, 0.35, 0.58, 1e-4),
        ]
        _, _, ok = frozen_slope_check(rows)
        assert not ok
