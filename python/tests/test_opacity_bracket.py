"""Tests for the grounded opacity bracket (`opacity_bracket.py`, Q18/Q10).

Replaces the inherited 0.1x-10x opacity bracket -- which was built for a *placeholder* Kramers
opacity -- with one grounded in the real TOPS/OPLIB data's published accuracy, propagated to `f`
through the sweep's own measured sensitivity."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from puffsat import opacity_bracket as ob


def test_published_uncertainty_matches_the_cited_source() -> None:
    """Farag et al. 2024 (ApJ 968:16, arXiv:2406.02845) section 2.1, citing Huebner & Barfield
    2014, states the OPLIB Rosseland-mean uncertainty by dominant process. These are the numbers
    the bracket rests on, so they are pinned verbatim."""
    assert ob.OPLIB_UNCERTAINTY["electron_scattering"] == pytest.approx(0.05)
    assert ob.OPLIB_UNCERTAINTY["free_free"] == pytest.approx(0.10)
    assert ob.OPLIB_UNCERTAINTY["bound_free"] == pytest.approx(0.20)
    assert ob.OPLIB_UNCERTAINTY["bound_bound"] == pytest.approx(0.30)


def test_cool_partially_ionized_states_take_the_widest_bracket() -> None:
    """The uncertainty grows as T falls and bound-bound transitions start contributing. Our
    coldest turnaround (16.5 kK, Zbar = 0.16 -- barely ionized, a forest of bound-bound lines)
    therefore takes 30%, while the 69 km/s states (Zbar = 4.4, lithium-like and above, where the
    paper notes uncertainties fall toward the hydrogenic limit) take 20%."""
    assert ob.bracket_fraction(16_497.0, 0.157) == pytest.approx(0.30)
    assert ob.bracket_fraction(137_796.0, 4.521) == pytest.approx(0.20)


def test_slope_is_measured_from_the_bracket_rows(tmp_path: Path) -> None:
    """`de_eff/dln(kappa)` is measured from the sweep's own opacity-scale rows, not assumed. A
    log-slope over the 0.3x and 3.0x rows: e_eff 0.50 -> 0.56 across ln(10) is 0.06/2.3026."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(
        "\n".join(
            json.dumps(
                {
                    "v": 28_000.0,
                    "rho_impact": 0.01,
                    "length": 10.0,
                    "opacity_scale": s,
                    "e_eff": e,
                }
            )
            for s, e in ((0.3, 0.50), (1.0, 0.53), (3.0, 0.56))
        )
    )
    slopes = ob.measure_slopes(p)

    assert slopes[(28_000.0, 0.01, 10.0)] == pytest.approx(0.06 / 2.302585, rel=1e-6)


def test_stalled_rows_are_rejected_rather_than_fitted(tmp_path: Path) -> None:
    """A real defect found in `sweep_heavyplate.jsonl` (2026-08-17): the v = 28 km/s, rho = 0.6,
    kappa = 10x row returns e_eff = 0.0562 where its neighbours run 0.6749-0.6783 — the radiative
    collapse of 4ddaed5, still live because `water_jupiter.json` never got that commit's rho-grid
    extension and still ceilings at 30 kg/m^3.

    Fitting a slope through it flips the sign (-0.134 against a true +0.001), so a stalled row must
    be rejected, not averaged in. It survived the Q20 regression gate because the baseline carried
    the identical corruption — old-vs-new comparison cannot see physical implausibility."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(
        "\n".join(
            json.dumps(
                {"v": 28e3, "rho_impact": 0.6, "length": 10.0, "opacity_scale": s, "e_eff": e}
            )
            for s, e in ((0.1, 0.6749), (0.3, 0.6763), (1.0, 0.6775), (3.0, 0.6783), (10.0, 0.0562))
        )
    )
    slopes, stalled = ob.measure_slopes_checked(p)

    assert stalled == [(28e3, 0.6, 10.0, 10.0)]
    # The 0.1/10.0 pair is unusable, so it falls back to the intact 0.3/3.0 pair.
    assert slopes[(28e3, 0.6, 10.0)] == pytest.approx((0.6783 - 0.6763) / 2.302585, rel=1e-6)


def test_bracket_translates_a_kappa_error_into_an_f_band() -> None:
    """The deliverable: a fractional opacity uncertainty becomes a band on `f` via the measured
    sensitivity and `f = eta*(1 + e_eff)/2`.

    Worked by hand at the 69 km/s, rho = 0.02, L = 12 m row (slope 0.0711, 20% bracket,
    eta = 0.98): d(e_eff) = 0.0711 * ln(1.2) = 0.012963, and d(f) = 0.49 * 0.012963 = 0.006352."""
    band = ob.f_band(slope=0.0711, fraction=0.20, eta_capture=0.98)

    assert band.delta_e_eff == pytest.approx(0.012963, rel=1e-4)
    assert band.delta_f == pytest.approx(0.006352, rel=1e-4)


def test_grounded_bracket_is_far_tighter_than_the_inherited_one() -> None:
    """The point of Q18. The inherited 0.1x-10x bracket spans ln(10) = 2.303 in either direction;
    the grounded 20-30% spans ln(1.2)-ln(1.3) = 0.18-0.26. That is a 9-13x narrower band on `f`,
    and the old one was sized for a placeholder Kramers opacity that is no longer in use."""
    grounded = ob.f_band(slope=0.0711, fraction=0.20, eta_capture=0.98)
    inherited = ob.f_band(slope=0.0711, fraction=9.0, eta_capture=0.98)  # 1x -> 10x

    assert inherited.delta_f / grounded.delta_f == pytest.approx(12.6, rel=0.05)
