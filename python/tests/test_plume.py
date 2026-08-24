"""The plume-state table `aim_is_all_you_need` cites (routing document, "Added 2026-08-21").

That repository supplies `w` and `rho = m_slug / V`; this one supplies the solve. The audit that
routed the item had already hand-computed four rows and reported agreement to 1-3% in `T`, so the
acceptance test is written against those rows -- the answer was known before the module existed.
"""

from __future__ import annotations

import pytest

from puffsat import eos_water, expansion, plume

#: The audit's hand table (routing document item 1): closing speed -> (dissipated MJ/kg, T, f).
HAND_TABLE = {
    75.0: (265.0, 26200.0, 0.573),
    65.0: (199.0, 22400.0, 0.371),
    56.53: (150.0, 19400.0, 0.217),
    45.58: (98.0, 14700.0, 0.053),
}


def test_dissipated_energy_reproduces_the_papers_own_arithmetic() -> None:
    """`(1/2) k w^2 / (1+k)^2` -- the energy a merged slug keeps after momentum sharing.

    Closed form, so it is checked against the paper's stated MJ/kg rather than against itself.
    """
    for speed, (mj, _, _) in HAND_TABLE.items():
        assert plume.dissipated_energy(speed * 1.0e3) / 1.0e6 == pytest.approx(mj, rel=0.005)


def test_the_solve_reproduces_the_hand_table_to_the_advertised_few_percent() -> None:
    """The acceptance test: `eos_water` against the audit's single-species Saha, at the four
    anchors the paper's burn envelope is quoted at.

    The solve runs **warmer** at every anchor, and the direction is not an accident: the audit
    charged 54 MJ/kg for vaporisation plus dissociation, while `eos_water`'s bond energy is
    50.9 MJ/kg, leaving ~3 MJ/kg more in the thermal pool. The gap is largest at the cold anchor
    because there is least energy there for it to hide in.
    """
    for speed, (_, hand_t, hand_f) in HAND_TABLE.items():
        state = plume.plume_state(speed * 1.0e3, plume.BAG_RHO)
        assert state.temp == pytest.approx(hand_t, rel=0.035), f"{speed} km/s temperature"
        assert state.ionised_fraction == pytest.approx(hand_f, abs=0.01), f"{speed} km/s ionisation"
        assert state.temp > hand_t, "the solve must run warmer, for the reason in the docstring"


def test_the_energy_reference_is_the_trap_and_double_charging_wrecks_the_cold_anchor() -> None:
    """`eos_water` references `e` to **bound molecular H2O at T -> 0**, so dissociation is already
    inside `e`. Subtracting the audit's 54 MJ/kg as well -- the natural move when porting its
    formula -- double-charges the bond energy.

    It is silently plausible at the hot end (23 000 K, only 11% low) and catastrophic at the cold
    end, where 54 MJ/kg is over half the budget. This test exists because that failure mode looks
    like a small calibration error until it does not.
    """
    speed = 45.58e3
    correct = plume.plume_state(speed, plume.BAG_RHO)
    double_charged = expansion.temperature_at(
        plume.BAG_RHO, plume.dissipated_energy(speed) - 54.0e6, eos_water.pressure_energy
    )
    assert correct.temp == pytest.approx(14700.0, rel=0.035)
    assert double_charged < 6000.0, "double-charging the bond energy loses two thirds of the answer"


def test_a_denser_bag_runs_hotter_and_less_ionised_at_the_same_closing_speed() -> None:
    """The reason the table has to be two-dimensional rather than a list of four numbers.

    Specific dissipated energy does not depend on `rho` at all -- it is set by `w` and `k`. But
    Saha does: compressing the plume pushes recombination, so less of that fixed budget is spent
    stripping electrons and more stays as heat. `aim` sets `rho = m_slug / V` from the bag volume,
    which is a live design variable, so it cannot be handed a single row.
    """
    speed = 56.53e3
    thin = plume.plume_state(speed, 0.05)
    thick = plume.plume_state(speed, 1.0)

    assert thick.temp > thin.temp, "denser is hotter at fixed specific energy"
    assert thick.ionised_fraction < thin.ionised_fraction, "and less ionised"
    assert thick.energy == pytest.approx(thin.energy), "the budget itself is density-independent"
    # 20x in density is worth ~30% in temperature -- not a detail that can be rounded away.
    assert thick.temp / thin.temp > 1.2


def test_the_published_table_is_what_the_cooling_history_is_anchored_on() -> None:
    """`expansion.PLUME_STATES` carries the audit's hand temperatures, and everything downstream
    of it -- the cooling history, Q-M's freeze verdict, ADR-0038's stability scan -- inherits them.

    So the two must not drift apart silently. This asserts they agree to the same few percent the
    acceptance test allows; if the anchors are ever re-based onto this solve, this is the test that
    has to be updated deliberately rather than the one that quietly starts failing.
    """
    for speed, anchor_temp in expansion.PLUME_STATES:
        solved = plume.plume_state(speed * 1.0e3, expansion.BAG_RHO).temp
        assert solved == pytest.approx(anchor_temp, rel=0.035), f"{speed} km/s drifted"


def test_the_table_spans_the_burn_envelope_and_stays_monotone_in_speed() -> None:
    """The shipped grid, and the one shape claim worth making about it.

    More closing speed is more dissipated energy is more temperature, at every bag density. A
    non-monotone row would mean the root-find had picked up a second branch, which is the failure
    this catches.
    """
    rows = plume.table(plume.BURN_ENVELOPE, (0.1, plume.BAG_RHO))
    assert len(rows) == len(plume.BURN_ENVELOPE) * 2

    for rho in (0.1, plume.BAG_RHO):
        temps = [r.temp for r in rows if r.rho == rho]
        assert temps == sorted(temps), f"temperature must rise with closing speed at rho={rho}"
    assert min(plume.BURN_ENVELOPE) <= 45.58e3 <= max(plume.BURN_ENVELOPE)
    assert min(plume.BURN_ENVELOPE) <= 75.0e3 <= max(plume.BURN_ENVELOPE)
