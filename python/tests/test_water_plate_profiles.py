"""Mass-normalization acceptance tests for prescribed overtake inputs."""

import math

import pytest

from puffsat.water_plate.profiles import IncomingPulse, InjectionProfile, TimeBand


def test_incoming_mass_is_fixed_when_footprint_and_duration_change() -> None:
    for radius in (2.5, 3.75, 15.0):
        for duration in (10e-6, 100e-6, 1e-3):
            pulse = IncomingPulse(radius, duration)
            assert pulse.density_kg_m3 * math.pi * radius**2 * pulse.length_m == pytest.approx(25.0)


@pytest.mark.parametrize("power", [0, 2])
@pytest.mark.parametrize("pre_fraction", [0.0, 0.5, 1.0])
def test_space_time_flux_integrates_to_supplied_water(power: int, pre_fraction: float) -> None:
    profile = InjectionProfile(
        5.0, power, (TimeBand(-0.002, -0.001, pre_fraction), TimeBand(0, 1e-4, 1 - pre_fraction))
    )
    # Independent area quadrature with time integration per rectangular band.
    dr = profile.radius_m / 1000
    mass = sum(
        profile.mass_flux_kg_m2_s((i + 0.5) * dr, (band.start_s + band.end_s) / 2)
        * 2
        * math.pi
        * (i + 0.5)
        * dr
        * dr
        * (band.end_s - band.start_s)
        for band in profile.bands
        for i in range(1000)
    )
    assert mass == pytest.approx(212.5, rel=1e-6)
    assert profile.mass_between(-1, 1, 5.0) == pytest.approx(mass, rel=1e-6)
    assert profile.mass_between(-1, 0, 5.0) == pytest.approx(212.5 * pre_fraction)
    assert profile.mass_between(0, 5e-5, 5.0) == pytest.approx(212.5 * (1 - pre_fraction) / 2)
    # The complement lies outside the incoming footprint, not necessarily lost.
    assert profile.mass_between(-1, 1, 3.75) / 212.5 == pytest.approx(0.75 ** (power + 2))


def test_invalid_or_unnormalized_sources_are_rejected() -> None:
    with pytest.raises(ValueError):
        IncomingPulse(0, 1e-4)
    with pytest.raises(ValueError):
        InjectionProfile(5, 0, (TimeBand(0, 1, 0.5),))
    with pytest.raises(ValueError):
        TimeBand(1, 0, 1)
    with pytest.raises(ValueError):
        InjectionProfile(5, 0, (TimeBand(0, 2, 0.5), TimeBand(1, 3, 0.5)))
