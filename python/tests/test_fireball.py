"""Does recombination freeze below 0.01 kg/m^3? (routing document, `sec:watering_it_down`)

Q-M answered the *nozzle*: equilibrium, by two to five decades. It also found that the
dissociation store never returns there at all, and said its marginal Damkohler number would bind
downstream instead. This is downstream. The expected answer was obtained from a standalone probe
before the module existed, so these assert numbers rather than shapes.
"""

from __future__ import annotations

import math

import pytest

from puffsat import eos_water, expansion, fireball, plume

COLD_ANCHOR = 14700.0  # the 45.58 km/s bag state, the marginal one in Q-M


def test_the_clock_is_continuous_across_the_nozzle_lip() -> None:
    """The fireball history is the nozzle history plus a tail, not a second calculation.

    Inside the field the clock is `expansion.cooling_history`'s -- a linear area opening over the
    23.8 m bore. Outside it there is no bore, so the clock becomes `dR/(u tan theta)`. The two
    must meet: if they did not, the expansion time on the first free station would carry a step
    that is an artifact of the handover rather than physics, and the Damkohler number there is
    exactly what this module reports.
    """
    rows = fireball.history(COLD_ANCHOR)
    lip = [r for r in rows if r.area_ratio <= expansion.AREA_RATIO_EXIT][-1]

    nozzle = expansion.cooling_history(
        expansion.BAG_RHO,
        COLD_ANCHOR,
        eos_water.pressure_energy,
        eos_water.sound_speed,
        expansion.AREA_RATIO_EXIT,
        expansion.FIELD_LENGTH,
    )
    assert lip.time == pytest.approx(nozzle[-1].time, rel=1e-3)
    assert lip.radius == pytest.approx(6.0, rel=1e-6), "A/A* = 4 is a 6 m flux tube"


def test_the_expansion_clock_steps_down_at_the_lip_and_the_aggregate_hides_it() -> None:
    """The trigger, and the measure that shows it.

    The paper's bore opens `A/A*` 1 -> 4 over 23.8 m: ~7.9 m of travel per unit area ratio. A
    45-degree free jet covers the same unit in 0.75 m. So the local expansion time steps down ~8x
    across the lip -- the field was holding the plume in a slow expansion and has let go.

    **Averaged over the whole tail it does not show up at all**: the fireball is *slower* per
    decade than the nozzle (3.8 ms against 3.0), because far out `tau_exp` has grown large again.
    This test asserts both, because taking the aggregate for the answer gets the mechanism
    backwards -- which is how it was first written here.
    """
    stations = fireball.scan(COLD_ANCHOR)
    # The one interval that spans the lip is measured across both clocks and belongs to neither.
    clean = [s for s in stations if not s.straddles_lip]
    inside = [s for s in clean if s.inside]
    outside = [s for s in clean if not s.inside]

    assert inside[-1].tau_expansion > 5.0 * outside[0].tau_expansion, "the local step at the lip"
    assert inside[-1].da_dissociation_h_limited > 10.0 * outside[0].da_dissociation_h_limited

    rows = fireball.history(COLD_ANCHOR)

    def seconds_per_decade(seq: list[fireball.FireballRow]) -> float:
        decades = math.log10(seq[0].rho / seq[-1].rho)
        return (seq[-1].time - seq[0].time) / decades

    nozzle = [r for r in rows if r.area_ratio <= expansion.AREA_RATIO_EXIT]
    jet = [r for r in rows if r.area_ratio >= expansion.AREA_RATIO_EXIT]
    assert seconds_per_decade(jet) > seconds_per_decade(nozzle), "the aggregate runs the other way"


def test_the_dissociation_store_freezes_just_below_the_papers_own_threshold() -> None:
    """The finding. `Da` for the atomic three-body channel falls through 1 at ~1.1e-2 kg/m^3 --
    within 10% of the 0.01 kg/m^3 the paper names -- and it does so while the store is still full.

    This is the case Q-M's `binding_damkohler` was built to distinguish. There, the station with
    the smallest `Da` held 0.01% of the store and the freeze was an artifact. Here the crossing
    station still holds **~92%** of the dissociation store, so the freeze is real and it strands
    most of the bond energy.
    """
    stations = fireball.scan(COLD_ANCHOR)
    freeze = fireball.freeze_state(stations)

    assert freeze is not None, "the cold anchor must freeze"
    assert freeze.rho == pytest.approx(1.1e-2, rel=0.25), "freezes near the paper's 0.01 kg/m^3"
    assert freeze.bond_energy_fraction > 0.85, "with the store still substantially held"
    assert freeze.da_dissociation_h_limited < 1.0


def test_the_verdict_survives_the_jet_divergence_it_is_conditional_on() -> None:
    """`tan theta` is the one free parameter and it scales the clock linearly, so `Da ~ 1/tan`.

    A wider jet expands faster and freezes denser. Over a 15-60 degree bracket -- far wider than
    a magnetic nozzle's real divergence -- the freeze density moves less than a decade and the
    store is still >80% held at both edges. So the *number* is conditional on the angle and the
    *verdict* is not.
    """
    narrow = fireball.freeze_state(fireball.scan(COLD_ANCHOR, half_angle_deg=15.0))
    wide = fireball.freeze_state(fireball.scan(COLD_ANCHOR, half_angle_deg=60.0))

    assert narrow is not None and wide is not None
    assert wide.rho > narrow.rho, "a wider jet freezes at higher density"
    assert wide.rho / narrow.rho < 10.0, "less than a decade over a 4x range in tan"
    assert min(narrow.bond_energy_fraction, wide.bond_energy_fraction) > 0.80


def test_nothing_freezes_inside_the_nozzle_which_is_what_q_m_already_found() -> None:
    """Consistency with Q-M, on the channel Q-M reported as marginal.

    If this module said the dissociation store froze inside the field it would contradict the
    committed answer, and one of the two would be wrong. It does not: every station at or below
    the nozzle exit clears `Da = 1`, and the crossing is outside.
    """
    stations = fireball.scan(COLD_ANCHOR)
    inside = [s for s in stations if s.inside]
    assert inside, "the scan must cover the nozzle as well as the fireball"
    assert all(s.da_dissociation_h_limited > 1.0 for s in inside), (
        "Q-M's equilibrium verdict must hold"
    )


def test_the_stranded_energy_is_reported_against_the_budget_that_produced_it() -> None:
    """What the freeze costs, which is the only reason the question is worth asking.

    The bond energy still held at the crossing never becomes directed kinetic energy -- it leaves
    as inert chemical enthalpy. Reported as a fraction of the dissipated energy that made the
    plume, because that is the quantity an Isp claim is built on.
    """
    stations = fireball.scan(COLD_ANCHOR)
    freeze = fireball.freeze_state(stations)
    assert freeze is not None

    stranded = fireball.stranded_energy(freeze)
    assert stranded == pytest.approx(0.92 * eos_water.FULL_ATOMIZATION_ENERGY, rel=0.15)
    # Against the 45.58 km/s budget it is a large fraction, not a correction.
    assert stranded / plume.dissipated_energy(45.58e3) > 0.3


def test_the_oh_bracket_is_zero_width_on_the_hot_legs() -> None:
    """The finding of 2026-08-25, and it is that the apologised-for proxy did not matter.

    `recombination` rated water reformation with `n_H` because the species set had no OH, and said
    in a comment that this *over*estimates the rate. It does, by 30x to 50000x. But the stranded
    energy is `bond_energy_fraction * FULL_ATOMIZATION_ENERGY` at the freeze station, and on the
    hot legs that fraction is already 1.0000 at the optimistic edge -- the store is fully held when
    the clock runs out either way. So a five-decade error in the rate moves the answer by nothing.

    Where a bracket exists at all is the cold leg, and it is narrow: see the test below.
    """
    for temp_0 in (26_200.0, 22_400.0):
        stations = fireball.scan(temp_0)
        optimistic = fireball.freeze_state(stations)
        conservative = fireball.freeze_state(stations, oh_limited=True)
        assert optimistic is not None and conservative is not None
        # The conservative edge freezes far earlier -- at the first station, inside the nozzle.
        assert conservative.rho > 5.0 * optimistic.rho
        # And it buys nothing, because there was nothing left to return by then either way.
        assert fireball.stranded_energy(conservative) == pytest.approx(
            fireball.stranded_energy(optimistic), rel=1e-3
        )


def test_the_cold_leg_is_the_only_leg_with_a_bracket_and_it_is_narrow() -> None:
    """On the cold anchor the optimistic edge does return part of the store, so the edges differ.

    They differ by less than the jet-divergence bracket the answer already carried: `theta` over
    15-60 degrees moves the optimistic edge 0.83 -> 0.94 of the store, and the conservative edge
    sits at 0.999. So the OH question is a smaller uncertainty than the geometry one, which is
    what makes Q-P's published range survive.
    """
    stations = fireball.scan(COLD_ANCHOR)
    optimistic = fireball.freeze_state(stations)
    conservative = fireball.freeze_state(stations, oh_limited=True)
    assert optimistic is not None and conservative is not None

    assert optimistic.bond_energy_fraction < conservative.bond_energy_fraction
    assert conservative.bond_energy_fraction > 0.99, "nothing comes back on the slow edge"
    assert 0.85 < optimistic.bond_energy_fraction < 0.95
