"""Acceptance tests for the `eta_chem(w, k)` surface owed to `aim_is_all_you_need` (Q-R).

The closed form is exact *given* `phi`, so most of these pin the algebra and its limits. The one
that carries physics is the last: `phi` must be solved, not assumed, because the cold high-`k`
corner is where a charged-in-full toll would invent energy the plume never stored.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import eos_water, plume, toll


def test_eta_is_one_when_no_bond_energy_is_stranded() -> None:
    """`phi = 0` is the paper's own assumption -- the loan repaid -- and must cost nothing."""
    assert toll.eta_chem(56.53e3, 8.5, 0.0) == pytest.approx(1.0)


def test_eta_matches_the_energy_ratio_it_is_defined_as() -> None:
    """`eta_chem` is the ratio of real to ideal gas speed, so its square is the energy ratio.

    Stated independently of the closed form: `eta^2` must equal available over total collision
    energy per merged kg. If the two ever disagree the algebra has drifted from the definition.
    """
    w, k, phi = 56.53e3, 8.5, 0.92
    total = 0.5 * w * w / (1.0 + k)
    assert toll.available_energy(w, k, phi) == pytest.approx(total - phi * toll.E_BOND)
    assert toll.eta_chem(w, k, phi) ** 2 == pytest.approx(toll.available_energy(w, k, phi) / total)


def test_eta_falls_with_slug_ratio_and_rises_with_closing_speed() -> None:
    """The two monotonicities the whole finding rests on.

    The toll is a *fixed* subtraction while the budget scales as `w^2/(1+k)`, so diluting into more
    slug spends a fixed bill over a thinner budget and closing faster pays it more easily. This is
    why the cold end of the overtake is the exposed one.
    """
    etas_in_k = [toll.eta_chem(56.53e3, k) for k in (1.0, 2.0, 4.0, 8.5, 16.0)]
    assert all(b < a for a, b in pairwise(etas_in_k))
    etas_in_w = [toll.eta_chem(w * 1e3, 8.5) for w in (45.58, 56.53, 65.13, 75.0)]
    assert all(b > a for a, b in pairwise(etas_in_w))


def test_the_zero_boundary_is_where_the_collision_cannot_pay_the_bond_bill() -> None:
    """Above `k_max` there is no exhaust at all, which is a different statement from a poor one.

    `aim`'s chain optimiser searches `k` to 80 with no such bound. At 45.58 km/s the boundary is
    19.4, so most of that box describes a plume that never expands.
    """
    for w in (45.58e3, 56.53e3, 75.0e3):
        k_max = toll.zero_slug_ratio(w)
        assert toll.eta_chem(w, k_max * 0.99) > 0.0
        assert toll.eta_chem(w, k_max * 1.01) == 0.0
        # At the boundary the available energy is exactly spent.
        assert toll.available_energy(w, k_max, 1.0) == pytest.approx(0.0, abs=1e-6 * toll.E_BOND)


def test_the_bond_ceiling_is_the_eos_atomization_energy_not_a_transcribed_constant() -> None:
    """`E_BOND` must come from the same EOS that computes `phi`, or the two disagree silently.

    The paper quotes 50.4 MJ/kg from its own source; this repo derives 50.94 from JANAF heats of
    formation. They agree to 1%, and the check is that this module uses *ours*, since `phi` is
    measured against it.
    """
    assert toll.E_BOND is eos_water.FULL_ATOMIZATION_ENERGY
    assert 50.0e6 < toll.E_BOND < 51.5e6


def test_a_node_solves_end_to_end_and_agrees_with_the_closed_form_at_its_own_phi() -> None:
    """The solver and the closed form are the same statement, and the node proves the wiring.

    `toll_point` runs the real chain -- dissipated energy, stagnation state, freeze -- so this also
    pins that the stagnation solve is the one `plume.py` publishes rather than a second copy.
    """
    p = toll.toll_point(56.53e3, 8.5)
    assert p.dissipated == pytest.approx(plume.dissipated_energy(56.53e3, 8.5))
    assert p.temp_0 == pytest.approx(19_710.0, rel=1e-3)
    assert p.eta == pytest.approx(toll.eta_chem(56.53e3, 8.5, p.bond_fraction))
    assert p.lights
    # The growth push's hot end: the store is fully held, so the whole 50.9 MJ/kg strands.
    assert p.bond_fraction > 0.99


def test_phi_is_solved_rather_than_assumed_at_the_corner_where_it_matters() -> None:
    """The cold, dilute corner is the one place charging the full bill would invent a toll.

    At a low closing speed with the slug spread thin, the stagnation state is cold enough that
    equilibrium water is not fully atomized -- there is less store to strand, so the real `eta` is
    *above* what `phi = 1` would give. The surface has to be able to say so, which is the whole
    reason this is a solve and not a formula `aim` could evaluate itself.
    """
    cold = toll.toll_point(40.0e3, 0.5)
    assert cold.bond_fraction <= 1.0
    assert cold.eta >= toll.eta_chem(40.0e3, 0.5, 1.0) - 1e-12


def test_available_energy_is_never_claimed_beyond_what_was_dissipated() -> None:
    """The toll comes out of the thermal pool, so it cannot exceed what the collision dissipated.

    A node where the bond bill is larger than the dissipated energy is one the plume could not have
    reached in the first place, and reporting a positive `eta` there would be fiction.
    """
    for w in (40.0e3, 45.58e3, 56.53e3):
        for k in (0.5, 4.0, 8.5, 20.0):
            if toll.eta_chem(w, k) > 0.0:
                assert plume.dissipated_energy(w, k) > toll.E_BOND * (1.0 - 1e-9)


def test_the_closed_form_is_what_aim_should_multiply_and_it_is_not_the_whole_beta() -> None:
    """Guards the insertion point, which is the part of Q-R most likely to be got wrong.

    `eta` scales the gross jet `sqrt(1+k)`. Applying it to the whole head-on `beta` instead would
    also scale the `-1` momentum debit, which is pure momentum conservation. The two differ, and
    the difference is what restores the forward-thrust floor at `eta = 1/sqrt(1+k)`.
    """
    k, eta = 8.5, 0.731
    correct = (eta * math.sqrt(1.0 + k) - 1.0) / k
    wrong = eta * (math.sqrt(1.0 + k) - 1.0) / k
    assert correct < wrong, "scaling the debit too flatters the head-on leg"
    # The floor exists only in the correct form.
    floor = 1.0 / math.sqrt(1.0 + k)
    assert (floor * math.sqrt(1.0 + k) - 1.0) / k == pytest.approx(0.0, abs=1e-12)
    assert floor == pytest.approx(0.324, abs=0.001)


def test_the_ignition_gate_is_what_binds_cold_and_eta_alone_would_mislead() -> None:
    """The finding of the first grid run, and the reason `ignites` ships alongside `eta`.

    `eta_chem` is **not monotone in `k`** on the solved surface, and the closed form at `phi = 1`
    is badly wrong where it turns: at 40 km/s it says there is no exhaust above `k` = 14.7, while
    the solve reads 0.70 at `k` = 20. Both are right about their own question. The plume there
    never fully dissociates, so it strands far less than the full bill -- the toll is
    self-limiting. But it also never becomes a conducting plasma, and a magnetic nozzle has
    nothing to grip. Every node whose `phi` falls meaningfully below 1 is a node that fails the
    paper's own ignition bill, so a consumer must read `ignites` before `eta`.
    """
    cold = toll.toll_point(40.0e3, 20.0)
    assert not cold.ignites, "the corner where phi < 1 must be flagged as not lighting"
    assert cold.bond_fraction < 0.5, "and it is cold because it never dissociated"
    assert cold.eta > 0.6, "so eta recovers there, which is the trap"
    assert toll.eta_chem(40.0e3, 20.0, 1.0) == 0.0, "while phi = 1 would call it dead"

    hot = toll.toll_point(56.53e3, 8.5)
    assert hot.ignites and hot.bond_fraction > 0.99


def test_the_closed_form_is_conservative_wherever_the_plume_ignites() -> None:
    """What makes the deliverable usable, stated at the strength the grid actually supports.

    `phi <= 1` by construction and `eta` falls with `phi`, so `eta_chem(w, k, 1.0)` can never
    overstate the surface -- it is a floor, and `aim` may use it without this repository's solver
    provided it stays inside the igniting region. Across the 71 igniting nodes of the published
    grid `phi` bottoms out at **0.859** and the floor understates the solve by at most **0.057**.

    It is not exact, and the cold anchor is why: at 45.58 km/s and `k` = 8.5 the plume is cool
    enough at the lip that some water has re-formed, so `phi` = 0.927 and the real ceiling is
    0.754 against the floor's 0.731. Conservative in the direction that matters.
    """
    for w, k in ((45.58e3, 8.5), (56.53e3, 8.5), (65.13e3, 8.5), (75.0e3, 8.5), (75.0e3, 16.0)):
        p = toll.toll_point(w, k)
        assert p.ignites
        assert p.bond_fraction > 0.85, f"w={w / 1e3}, k={k}: phi={p.bond_fraction}"
        floor = toll.eta_chem(w, k, 1.0)
        assert p.eta >= floor - 1e-9, "the phi = 1 form must never overstate the surface"
        assert p.eta - floor < 0.06, "and must not be uselessly slack either"

    cold = toll.toll_point(45.58e3, 8.5)
    assert cold.bond_fraction == pytest.approx(0.927, abs=0.01)
    assert cold.eta == pytest.approx(0.754, abs=0.005)


def test_the_temperature_gate_separates_the_mush_and_the_seed_is_why_it_can_sit_at_10_kK() -> None:
    """Seth's gate, 2026-08-25, and the check that 10 000 K is not set too high.

    The worry was that a 10 kK floor might be stricter than the potassium seed requires. It is not,
    and conductivity is not what fails down there at all: the seed holds `sigma` at 1480 S/m at
    10 000 K and still 403 S/m at 4500 K. What fails is dissociation -- the excluded nodes are a
    barely-broken mush at `phi` = 0.38-0.75, while everything admitted sits at `phi >= 0.777`,
    which is the region where the closed-form floor is tight.
    """
    admitted = toll.toll_point(56.53e3, 8.5)
    assert admitted.conducts and admitted.bond_fraction > 0.99

    mush = toll.toll_point(40.0e3, 20.0)
    assert not mush.conducts
    assert mush.temp_0 < 6000.0
    assert mush.bond_fraction < 0.5


def test_admissible_k_is_a_closed_interval_because_the_energy_peaks_at_k_equals_one() -> None:
    """The bound that replaces `aim`'s flat `_K_SEARCH_MAX = 80`.

    `k/(1+k)^2` peaks at `k = 1`, so too much slug spreads a fixed energy too thin and too little
    dissipates almost nothing -- the admissible set is an interval, not a ceiling. The paper says
    the same in `sec:two_leg_nozzle` and gives [0.098, 10.21] at 45.58 km/s against its own
    85.1 MJ/kg bill; the 10 kK gate is slightly looser at [0.081, 12.29], which is the expected
    direction because that bill targets 15 000 K.
    """
    span = toll.admissible_slug_ratios(45.58e3)
    assert span is not None
    lo, hi = span
    assert lo == pytest.approx(0.081, abs=0.005)
    assert hi == pytest.approx(12.29, abs=0.1)
    assert lo < 1.0 < hi, "the peak must be inside the interval"

    for k in (lo * 1.05, 1.0, hi * 0.95):
        assert toll.stagnation_temperature(45.58e3, k, plume.BAG_RHO) >= toll.MIN_STAGNATION_TEMP
    for k in (lo * 0.9, hi * 1.1):
        assert toll.stagnation_temperature(45.58e3, k, plume.BAG_RHO) < toll.MIN_STAGNATION_TEMP

    # The interval widens with closing speed, so the *coldest* pulse a leg sees is what binds it.
    fast = toll.admissible_slug_ratios(75.0e3)
    assert fast is not None and fast[1] > hi


def test_the_flown_slug_ratio_is_inside_every_admissible_interval_it_has_to_serve() -> None:
    """`k` = 8.5 must clear the coldest pulse of every leg it flies, which is the real constraint.

    Seth, 2026-08-25: `k` may differ between the head-on and overtake legs but cannot be
    reconfigured per pulse, so each leg carries one value that has to serve its whole falling
    speed range. The 3-synodic overtake ends at 45.58 km/s and that is what binds.
    """
    overtake = toll.admissible_slug_ratios(45.58e3)
    headon = toll.admissible_slug_ratios(75.0e3)
    assert overtake is not None and headon is not None
    overtake_bound, headon_bound = overtake[1], headon[1]
    assert overtake_bound > 8.5, "the flown k must clear the coldest pulse of the slowest cycle"
    assert headon_bound > overtake_bound, "and the head-on leg has room the overtake does not"
