"""Acceptance tests for the walled-nozzle chamber ledger (ADR-0050, asks N9 item 0 / N10 item 4b).

The exit criterion is stated before the solver runs, as the repo's TDD discipline requires: ask
N9 gives its own geometry -- 474 kg of methane on a 25 kg impactor, `k = 18.96` -- and the
fixed-point identity has to reproduce it without being told. It does, at the chamber and
temperature the ask names, which is what licenses the rest of the table.
"""

from __future__ import annotations

import math

import pytest

from puffsat import eos_methane
from puffsat.walled_nozzle import chamber


def test_the_asks_own_geometry_falls_out_of_the_identity() -> None:
    """`(1+k) u = w^2/2` at the flown 10 kK reproduces N9's stated `k = 18.96` to under a percent.

    **This is the acceptance test.** The ask supplies `k` as an input to its own wall question;
    here it is an output of an equilibrium solve that was never shown it. Agreement means the
    energy identity and the EOS are the same ones the companion repository is using, so a
    disagreement anywhere else in this table is a real disagreement rather than a mismatch of
    conventions.
    """
    state = chamber.solve_chamber(673.0, 23.8, chamber.FLOWN_TEMPERATURE)
    assert state.slug_ratio == pytest.approx(chamber.ASK_SLUG_RATIO, rel=0.01)
    assert state.slug_mass == pytest.approx(chamber.ASK_SLUG_MASS, rel=0.01)


def test_the_fixed_point_satisfies_its_own_energy_identity() -> None:
    """Whatever `k` comes out, `(1+k) u` must be `w^2/2` at the solved density -- the check that
    the iteration converged to the equation rather than merely to a number."""
    for volume, length in chamber.CHAMBERS:
        for temp in chamber.TEMPERATURES:
            state = chamber.solve_chamber(volume, length, temp)
            assert state.rho == pytest.approx(
                (1.0 + state.slug_ratio) * chamber.IMPACTOR_MASS / volume, rel=1e-9
            )
            _, u = eos_methane.pressure_energy(state.rho, temp)
            assert (1.0 + state.slug_ratio) * u == pytest.approx(
                chamber.SPECIFIC_PULSE_ENERGY, rel=1e-6
            )


def test_pulse_energy_is_the_adr_headline() -> None:
    """25 kg at 75 km/s is ADR-0016's 70.3 GJ pulse."""
    assert pytest.approx(70.3e9, rel=0.01) == chamber.PULSE_ENERGY


def test_pressures_bracket_the_adr_quoted_values() -> None:
    """ADR-0016 quotes 636 bar at 200 m^3 and N10 names 206 bar as the expansion's start. Both
    are reproduced within a few percent, which is what says the density and the particle count
    agree even where the composition does not."""
    small = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    large = chamber.solve_chamber(673.0, 23.8, chamber.FLOWN_TEMPERATURE)
    assert small.pressure / 1e5 == pytest.approx(636.0, rel=0.05)
    assert large.pressure / 1e5 == pytest.approx(206.0, rel=0.10)


def test_wall_cap_energy_density_matches_the_adr_ladder() -> None:
    """ADR-0016 orders its family of working fluids on the wall-cap energy density and puts
    methane at 143 MJ/kg at 10 kK. The solve gives 137-141 across the geometry range."""
    for volume, length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, length, chamber.FLOWN_TEMPERATURE)
        assert state.energy / 1e6 == pytest.approx(143.0, rel=0.06)


def test_sensible_fraction_matches_the_adrs_temperature_argument() -> None:
    """ADR-0016 declines a hotter chamber because methane's sensible fraction is 0.273 and only
    that part responds to `T`. The solve gives 0.277-0.287 at 10 kK, so the argument stands on a
    number this repository now owns rather than on an estimate."""
    for volume, length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, length, chamber.FLOWN_TEMPERATURE)
        assert state.sensible_fraction == pytest.approx(0.273, abs=0.02)


def test_the_store_is_essentially_fully_charged_at_the_flown_point() -> None:
    """**The finding that reverses the ask's premise.** N10 item 4b carries 49-76% dissociated;
    equilibrium at the flown 10 kK gives 93-98%, and the store 95-99% charged."""
    for volume, length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, length, chamber.FLOWN_TEMPERATURE)
        assert state.hydrogen_dissociated > 0.93
        assert state.store_charged > 0.94


def test_a_bigger_chamber_charges_more_store_but_barely() -> None:
    """The ask's *direction* is right and its *magnitude* is not: ADR-0016 books +3.4% on Isp
    from 200 to 673 m^3, where the solved charge moves 1.7% and the Isp identity 1.8%."""
    small = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    large = chamber.solve_chamber(673.0, 23.8, chamber.FLOWN_TEMPERATURE)
    assert large.store_charged > small.store_charged
    gain = chamber.isp_scaling(large, small.slug_ratio) - 1.0
    assert 0.0 < gain < 0.034


def test_the_c3_exposure_is_bounded_and_small_where_it_is_flown() -> None:
    """The model's weakest partition function must not be load-bearing at the flown point."""
    for volume, length in chamber.CHAMBERS:
        assert chamber.solve_chamber(volume, length, 10000.0).c3_store_exposure < 0.02
    # And it is honest about where it *is* load-bearing.
    assert chamber.solve_chamber(200.0, 7.1, 8000.0).c3_store_exposure > 0.05


def test_isp_scaling_is_the_identity_it_claims() -> None:
    """`sqrt(1+k)/k`, and unity against itself."""
    state = chamber.solve_chamber(400.0, 14.0, chamber.FLOWN_TEMPERATURE)
    assert chamber.isp_scaling(state, state.slug_ratio) == pytest.approx(1.0)
    k, ref = state.slug_ratio, chamber.ASK_SLUG_RATIO
    expected = (math.sqrt(1.0 + k) / k) / (math.sqrt(1.0 + ref) / ref)
    assert chamber.isp_scaling(state, ref) == pytest.approx(expected)


# --- N9 item 0 ----------------------------------------------------------------------------------


def test_the_chamber_sound_speed_is_methanes_not_hydrogens() -> None:
    """ADR-0016 runs item 0 at "an 11 km/s sound speed". That is *atomic hydrogen* at 10 kK.
    Dissociated methane averages 3.3 amu per particle, not 1, so it runs near 6 km/s -- and item
    0's equilibration count is proportional to it."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    assert 5.5e3 < state.sound_speed < 6.6e3
    hydrogen = math.sqrt(5.0 / 3.0 * eos_methane.K_B * 10000.0 / eos_methane.M_H)
    assert hydrogen == pytest.approx(11.7e3, rel=0.05)


def test_equilibration_count_is_an_area_ratio_and_the_column_length_cancels() -> None:
    """**The result that settles item 0, and it is not the one the ask expected.** Blowdown time
    goes as `M / (rho a A*)` and `M ~ rho V`, so the count is `V / (L A* f) = A_bore / (A* f)`.
    The column length divides out exactly -- so item 0 cannot discriminate between 200, 400 and
    673 m^3, and the short column is admissible on the same terms as the long one."""
    counts = []
    for volume, length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, length, chamber.FLOWN_TEMPERATURE)
        counts.append(chamber.acoustic_check(state, throat_area=7.0).equilibrations)
    assert max(counts) - min(counts) < 0.05 * min(counts)


def test_the_impactor_crosses_the_column_supersonically() -> None:
    """The separation that keeps N9 items 1-3 alive. 75 km/s against a 6 km/s sound speed means
    the column completes under a tenth of an acoustic crossing while the impactor is inside it:
    equilibration sets `k`, but it does **not** soften the wall's arrival transient."""
    for volume, length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, length, chamber.FLOWN_TEMPERATURE)
        check = chamber.acoustic_check(state)
        assert check.crossings_during_impact < 0.1
        assert check.crossing_time < check.acoustic_time


@pytest.mark.parametrize("throat_area", [2.0, 4.0, 7.0])
def test_the_sealed_vessel_verdict_depends_on_the_throat(throat_area: float) -> None:
    """**Item 0's answer is conditional on the throat, and it fails at the wide end.**

    `n_eq` is ~24 at 2 m^2, ~12 at 4 m^2 and ~6.8 at 7 m^2 -- so against this study's own
    `n_eq >= 10` bar the sealed vessel sets `k` only for throats of about 4 m^2 or narrower. At
    7 m^2 the chamber empties before it has equilibrated and `k` reverts to being set by what the
    cone sweeps, which is the coupling problem the walled chamber was adopted to escape.

    None of these is the ~50 ADR-0016 assumes: that used hydrogen's sound speed (W4) *and* a
    30 ms blowdown that only the small-throat end delivers.

    **4 m^2 is the crossover and it is load-bearing**, because it is where W6's Isp
    recommendation and this criterion agree.
    """
    state = chamber.solve_chamber(400.0, 14.0, chamber.FLOWN_TEMPERATURE)
    check = chamber.acoustic_check(state, throat_area)
    assert 5.0 < check.equilibrations < 30.0
    assert check.sealed_vessel_sets_k == (throat_area <= 4.0)


def test_blowdown_stretches_as_the_throat_narrows() -> None:
    """The cost side of W6's recommendation. Choked mass flow goes as `A*`, so the pulse lengthens
    in proportion as the throat shrinks -- 8 ms to 28 ms at 200 m^3 across the ask's 7 -> 2 m^2
    range. That is 3.5x longer for the wall to absorb the same pulse, and it is *not* costed in
    this study: it belongs to N9 items 1-7."""
    state = chamber.solve_chamber(200.0, 7.1, chamber.FLOWN_TEMPERATURE)
    wide = chamber.acoustic_check(state, throat_area=7.0).blowdown_time
    narrow = chamber.acoustic_check(state, throat_area=2.0).blowdown_time
    assert narrow / wide == pytest.approx(7.0 / 2.0, rel=1e-6)
