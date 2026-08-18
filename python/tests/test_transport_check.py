"""Tests for the Sn transport-check verdict (`transport_check.py`, Q9/Q21).

The Rust audit reports a *relative bias on a loss channel*. That is not yet an answer: a 40% bias
on a channel carrying 0.01% of the energy budget cannot move `f`, while a 5% bias on a channel
carrying half of it certainly can. This module is the translation into the currency the study
actually reports -- `delta_f` -- and the gate that reads it."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from puffsat import transport_check as tc


def _row(**kw: float | bool) -> str:
    base: dict[str, float | bool] = {
        "v": 45e3,
        "rho_impact": 0.02,
        "length": 10.0,
        "e_eff": 0.65,
        "loss_escape_space": 1.0e6,
        "transport_escape_rosseland": 1.05e6,
        "transport_escape_planck": 1.10e6,
        "relative_bias": 0.05,
        "relative_bias_planck": 0.10,
        "worst_relative_difference": 0.2,
        "worst_optical_depth": 1.5,
        "flux_weighted_optical_depth": 3.0,
        "escape_share_of_ke": 0.02,
        "converged": True,
    }
    base.update(kw)
    return json.dumps(base)


def test_a_bias_is_translated_into_delta_f_by_the_energy_it_misplaces() -> None:
    """The deliverable. A relative bias `beta` on a channel carrying fraction `phi` of the incident
    kinetic energy misplaces `phi*beta` of the budget. Rebound KE is `e_eff^2` of incident, so
    `d(e_eff)/e_eff = (1/2)*d(KE)/KE` gives `d(e_eff) = phi*|beta| / (2*e_eff)`, and
    `f = eta*(1 + e_eff)/2` carries it to `d(f) = eta*d(e_eff)/2`.

    Worked by hand at phi = 0.02, beta = 0.05, e_eff = 0.65, eta = 0.98:
    d(e_eff) = 0.02*0.05/(2*0.65) = 7.6923e-4, and d(f) = 0.49*7.6923e-4 = 3.7692e-4."""
    impact = tc.delta_f_from_bias(bias=0.05, escape_share=0.02, e_eff=0.65, eta_capture=0.98)

    assert impact.delta_e_eff == pytest.approx(7.6923e-4, rel=1e-4)
    assert impact.delta_f == pytest.approx(3.7692e-4, rel=1e-4)


def test_the_gate_reads_the_bias_not_the_translated_number() -> None:
    """Design SS12.1 step 5 states the escalation rule on the bias itself: >10% disagreement
    escalates to a coupled M1. That threshold is honoured verbatim, so the decision does not
    quietly become "whatever we judge small" once the energy weighting is applied."""
    assert tc.verdict(bias=0.05) == "PASS"
    assert tc.verdict(bias=-0.099) == "PASS"
    assert tc.verdict(bias=0.11) == "ESCALATE"
    assert tc.verdict(bias=-0.30) == "ESCALATE"
    # Exactly at the threshold is not "above" it.
    assert tc.verdict(bias=0.10) == "PASS"


def test_non_converged_rows_are_rejected_rather_than_judged(tmp_path: Path) -> None:
    """Same contract as everywhere else since Q6: `converged = false` means *no result*, not a low
    one. A stalled bounce reports a truncated escape integral, so its bias is an artifact of where
    the run stopped -- reading it as a transport verdict would be reading a solver bug as
    physics."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(
        "\n".join(
            [
                _row(v=28e3, converged=True),
                _row(v=45e3, converged=False, relative_bias=0.9),
                _row(v=63e3, converged=True),
            ]
        )
    )
    rows, rejected = tc.load(p)

    assert [r.v for r in rows] == [28e3, 63e3]
    assert [r.v for r in rejected] == [45e3]


def test_mean_selection_spread_is_reported_beside_the_closure_number(tmp_path: Path) -> None:
    """A gray transport solve has one extinction coefficient; FLD has two (Planck for emission,
    Rosseland for diffusion). So no single-mean Sn run is "FLD minus the closure", and the honest
    output carries both tallies. The spread between them is how much of the disagreement is really
    about mean selection rather than about diffusion-versus-transport."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(_row(relative_bias=0.04, relative_bias_planck=0.31))
    rows, _ = tc.load(p)

    assert rows[0].mean_selection_spread == pytest.approx(0.27, rel=1e-9)


def test_summary_reports_the_worst_case_by_delta_f_not_by_bias(tmp_path: Path) -> None:
    """Ranking by raw bias would put a huge relative error on an energetically irrelevant channel
    at the top of the table and bury the row that actually threatens `f`. The summary ranks by the
    translated `delta_f`.

    Here the 55 km/s row has 8x the bias but 1/200th the energy share, so it must rank below."""
    p = tmp_path / "sweep.jsonl"
    p.write_text(
        "\n".join(
            [
                _row(v=55e3, relative_bias=0.08, escape_share=0.0001, escape_share_of_ke=0.0001),
                _row(v=28e3, relative_bias=0.01, escape_share_of_ke=0.20),
            ]
        )
    )
    out = tmp_path / "summary.csv"
    ranked = tc.write_summary(p, out)

    assert [r.row.v for r in ranked] == [28e3, 55e3]
    assert out.read_text().splitlines()[0] == tc.CSV_HEADER
    assert len(out.read_text().splitlines()) == 3
