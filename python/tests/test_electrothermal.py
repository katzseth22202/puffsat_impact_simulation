"""Velikhov stability walked along the cooling history (ADR-0038).

`test_conductivity.py` pins the criterion at single states. This file pins the **leg-level**
claims ADR-0038 actually reports -- which legs are exposed, where the cold leg crosses, and how
much of its transit it spends unstable -- because those are scans over `expansion.history()` and
no single-state test can reach them.
"""

from __future__ import annotations

import pytest

from puffsat import electrothermal, expansion

HOT, MID, COLD = 75.0, 56.53, 45.58


def _temp_for(closing_speed: float) -> float:
    """The bag temperature `expansion.PLUME_STATES` anchors to this closing speed."""
    for speed, temp in expansion.PLUME_STATES:
        if speed == closing_speed:
            return temp
    raise AssertionError(f"no anchor at {closing_speed}")


def test_the_field_profile_is_the_area_ratio_the_paper_states() -> None:
    """`B*/B = A/A*` by flux conservation, so the paper's `20 T -> 5 T` fixes the whole profile.

    This is the one place a field enters the stability calculation, and it is not an assumption:
    the same flux conservation that makes `A/A* = 4` the exit area ratio makes 5 T the exit field.
    """
    assert electrothermal.local_field(1.0) == pytest.approx(20.0)
    assert electrothermal.local_field(4.0) == pytest.approx(5.0)
    assert electrothermal.local_field(2.0) == pytest.approx(10.0)


def test_the_hot_legs_are_stable_at_every_station_not_just_at_the_exit() -> None:
    """ADR-0038's first unqualified claim, and the exit test alone does not establish it.

    The exit is the coldest station, so it is the natural place to look -- but `beta` *rises*
    through the nozzle while `beta_cr` falls, so a leg could in principle cross and recover. It
    does not: on both hot legs `beta_cr` exceeds `beta` at every station by two decades or more.
    """
    for speed in (HOT, MID):
        rows = electrothermal.scan(speed, _temp_for(speed))
        assert rows, "the scan must produce stations"
        assert not any(r.unstable for r in rows), f"{speed} km/s must be stable throughout"
        worst = max(r.beta / r.critical_hall_parameter for r in rows)
        assert worst < 0.01, f"{speed} km/s clears the criterion by >100x, worst was {worst}"


def test_the_cold_leg_crosses_once_and_stays_unstable_to_the_exit() -> None:
    """The finding, and the shape of it matters as much as the fact.

    A leg that flickered in and out of instability would be a numerical artifact; one that crosses
    once and stays over is a monotone consequence of the plume cooling. The dwell is what makes it
    reportable -- microsecond e-folding times inside a stretch measured in milliseconds.
    """
    rows = electrothermal.scan(COLD, _temp_for(COLD))
    flags = [r.unstable for r in rows]

    assert any(flags), "the cold leg is the exposed one"
    assert flags[-1], "and it is unstable at the exit"
    assert flags.count(True) == len(flags) - flags.index(True), "one crossing, no recovery"

    crossing = electrothermal.crossing(rows)
    assert crossing is not None
    assert crossing.area_ratio == pytest.approx(1.968, rel=1e-3), "crosses mid-nozzle"
    assert crossing.temp == pytest.approx(7569.0, rel=1e-3)

    dwell = electrothermal.unstable_dwell(rows)
    assert dwell == pytest.approx(rows[-1].time - crossing.time), "one crossing: dwell is the tail"
    assert dwell == pytest.approx(1.70e-3, rel=0.02), "1.70 ms of a 2.71 ms transit"
    assert rows[-1].time == pytest.approx(2.71e-3, rel=0.01)


def test_the_crossing_is_resolved_rather_than_bracketed_by_the_station_spacing() -> None:
    """The deliverable is a threshold crossing, so the station grid -- not the integration step --
    is what has to be converged, and the default had to be chosen for that.

    At `expansion.history`'s own stride the crossing lands anywhere in `A/A*` 1.97-2.28 and the
    dwell reads 1.8-2.1 ms. This is the check that the shipped default is past that, and it is
    cheap: refining `steps` costs a fresh cooling history, refining `stride` costs ~15 ms a
    station.
    """
    coarse = electrothermal.scan(COLD, _temp_for(COLD), steps=320)
    fine = electrothermal.scan(COLD, _temp_for(COLD))

    c_coarse, c_fine = electrothermal.crossing(coarse), electrothermal.crossing(fine)
    assert c_coarse is not None and c_fine is not None
    assert c_fine.area_ratio == pytest.approx(c_coarse.area_ratio, rel=1e-3)
    assert electrothermal.unstable_dwell(fine) == pytest.approx(
        electrothermal.unstable_dwell(coarse), rel=1e-2
    )


def test_suppressing_the_skin_makes_every_station_stable() -> None:
    """ADR-0038 Addendum 3: the verdict rests on the current-layer thickness, which is unsolved.

    Putting the current across the full flux-tube radius instead of a resistive skin is a factor
    ~3 in thickness, hence ~9x in Joule heating and in the elevation. That is enough to flip every
    unstable station. The test exists so the sensitivity cannot quietly stop being true: this is
    the reason the cold-leg finding is reported as model-dependent rather than as a result.
    """
    skin = electrothermal.scan(COLD, _temp_for(COLD))
    tube = electrothermal.scan(COLD, _temp_for(COLD), use_skin_depth=False)

    assert any(r.unstable for r in skin)
    assert not any(r.unstable for r in tube), "the alternative reading is stable throughout"

    # Same station, both readings: the elevation must fall and beta_cr must rise with it.
    assert tube[-1].elevation < skin[-1].elevation / 5.0
    assert tube[-1].critical_hall_parameter > 4.0 * skin[-1].critical_hall_parameter
