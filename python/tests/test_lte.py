"""Tests for the LTE validity check (`lte.py`, Q5 of the 16-63 km/s extension).

The sweep's EOS is an equilibrium Saha ladder and its opacity is LTE TOPS/OPLIB data, so both
assume the level populations are collisionally controlled. These tests pin the McWhirter threshold
that decides whether that assumption holds at the states the sweep actually visits."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from puffsat import lte


def test_mcwhirter_threshold_matches_the_hydrogen_reference() -> None:
    """The textbook anchor: hydrogen at 10 000 K, Lyman-alpha (10.2 eV).

    `n_e >= 1.6e12 * sqrt(T) * dE^3 cm^-3` gives 1.6e12 * 100 * 1061.2 = 1.70e17 cm^-3, the
    ~10^17 cm^-3 figure quoted for hydrogen at 10^4 K. Returned in SI (m^-3), so 1.70e23."""
    assert lte.mcwhirter_threshold(10_000.0, 10.2) == pytest.approx(1.698e23, rel=1e-3)


def test_critical_gap_inverts_the_threshold() -> None:
    """`critical_gap` answers the criterion's useful inverse: given the state a sweep actually
    reaches, what is the largest level gap still collisionally controlled there?

    That inversion is the robust way to report LTE validity, because it needs no atomic data —
    the species' resonance gaps enter only afterwards, as a comparison. Round-trip: a gap fed
    through `mcwhirter_threshold` and back must return itself."""
    for temp, gap in ((16_497.0, 10.199), (137_796.0, 12.01), (48_057.0, 14.88)):
        n_e = lte.mcwhirter_threshold(temp, gap)
        assert lte.critical_gap(temp, n_e) == pytest.approx(gap, rel=1e-12)


def test_critical_gap_at_the_coldest_probe_state() -> None:
    """Worked by hand at the 16 km/s, rho=0.01 turnaround (T* = 16497 K, n_e = 1.694e24 m^-3):

        dE_crit = (n_e / (1.6e18 * sqrt(T)))^(1/3)
                = (1.694e24 / (1.6e18 * 128.44))^(1/3)
                = (8.240e3)^(1/3) = 20.2 eV

    Comfortably above the 10.2 eV Lyman-alpha gap that dominates this state, so the criterion is
    satisfied for hydrogen there with room to spare."""
    assert lte.critical_gap(16_497.0, 1.694e24) == pytest.approx(20.2, rel=5e-3)


def test_bare_protons_carry_no_criterion() -> None:
    """H+ is a bare nucleus: no bound levels, so no transition to hold in LTE. It dominates the
    69 km/s turnaround (67% of heavies) and must not be handed a gap of its own."""
    assert "H+" not in lte.RESONANCE_GAP_EV


def test_governing_species_is_the_most_demanding_abundant_one() -> None:
    """The criterion is set by the largest gap among species actually present, not by the most
    abundant one. At the 16 km/s turnaround the mix is H 55% / O 28% / H+ 12% / O+ 5%: hydrogen
    dominates by number, but O+ (14.86 eV) has a larger gap than H I (10.20 eV) and so governs.
    Trace species below the abundance floor are excluded — a 0.01%-abundance ion should not
    condemn the state."""
    abundances = {"H": 0.55, "O": 0.28, "H+": 0.12, "O1+": 0.05, "O2+": 1e-4}
    v = lte.evaluate_state(16_497.0, 1.694e24, abundances)

    assert v.governing == "O1+"
    assert v.governing_gap_ev == pytest.approx(14.86, rel=1e-2)
    # dE_crit = 20.2 eV > 14.86 eV, so the state clears its own governing gap.
    assert v.lte_valid
    # Worked by hand: dE = 12398.42/834.47 = 14.858 eV, dE^3 = 3279.9, sqrt(T) = 128.44, so the
    # threshold is 1.6e18 * 128.44 * 3279.9 = 6.740e23 m^-3 and n_e clears it by 2.51x.
    assert v.margin == pytest.approx(2.513, rel=1e-3)


def test_helium_like_oxygen_fails_by_orders_of_magnitude() -> None:
    """The physically important negative: at the 69 km/s turnaround, O6+ (helium-like, a 574 eV
    K-shell resonance gap) is present at ~1% and misses the criterion by ~4 orders of magnitude.
    High-charge stages in hot, fast plasmas are simply not in LTE, and the report must say so
    rather than average it away."""
    abundances = {"H+": 0.67, "O5+": 0.19, "O4+": 0.10, "O3+": 0.03, "O6+": 0.01}
    v = lte.evaluate_state(137_796.0, 1.901e25, abundances)

    assert v.governing == "O6+"
    assert not v.lte_valid
    assert v.margin < 1e-3
    # The abundant low-charge stages that carry the mass are fine; only the trace tail fails.
    assert v.failing == ("O6+",)


def test_o_ii_and_o_iii_resonance_lines_nearly_coincide() -> None:
    """Regression guard on the atomic data, which is easy to get wrong by a plausible amount.

    O II (83.4466 nm) and O III (83.5292 nm) both resonate in the solar "834 A multiplet", so
    their gaps are within 0.1% of each other and O III's is very slightly *smaller*. An earlier
    draft carried O III at 702.90 A (17.64 eV) — 16% high — which wrongly made O III rather than
    O II the governing species across the 28 km/s band and understated the margin there by 1.7x.
    Both values are NIST ASD (see `todos/nist_resonance/`)."""
    o2 = lte.RESONANCE_GAP_EV["O1+"].gap_ev
    o3 = lte.RESONANCE_GAP_EV["O2+"].gap_ev

    assert o2 == pytest.approx(14.858, rel=1e-3)
    assert o3 == pytest.approx(14.843, rel=1e-3)
    assert o3 < o2, "O III's resonance gap is marginally below O II's, so O II governs"


def test_probe_states_are_evaluated_from_a_probe_jsonl(tmp_path: Path) -> None:
    """The reader turns `(rho*, T*)` turnaround states into verdicts by solving the same
    equilibrium composition the sweep's EOS uses — so the LTE check is applied to the populations
    actually assumed, not to an independent estimate of them."""
    p = tmp_path / "probe.jsonl"
    p.write_text(
        json.dumps({"v": 16_000.0, "rho_impact": 0.01, "rho_star": 0.100, "t_star": 16_497.0})
        + "\n\n"  # trailing blank line tolerated, as elsewhere in the repo
    )
    verdicts = lte.evaluate_probe(p)

    assert len(verdicts) == 1
    row, v = verdicts[0]
    assert row.v == 16_000.0
    assert v.temp_k == 16_497.0
    # The composition solver reproduces the n_e that made dE_crit = 20.2 eV in the test above.
    assert v.critical_gap_ev == pytest.approx(20.2, rel=0.02)
    assert v.governing == "O1+"
    assert v.lte_valid


def test_threshold_is_cubic_in_the_energy_gap() -> None:
    """The gap enters cubed, which is why the criterion is set by the *largest* gap in the level
    structure and not by a typical one: doubling the gap costs 8x in required density."""
    base = lte.mcwhirter_threshold(50_000.0, 5.0)
    assert lte.mcwhirter_threshold(50_000.0, 10.0) == pytest.approx(8.0 * base, rel=1e-12)


def test_threshold_is_square_root_in_temperature() -> None:
    """Temperature enters only as sqrt, so the criterion is far more forgiving in T than in dE."""
    base = lte.mcwhirter_threshold(10_000.0, 9.0)
    assert lte.mcwhirter_threshold(40_000.0, 9.0) == pytest.approx(2.0 * base, rel=1e-12)
