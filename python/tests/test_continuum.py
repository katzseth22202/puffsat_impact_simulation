"""Acceptance tests for R1's premise -- is the nozzle expansion collisionless?

The exit criterion is known before the module exists, and it is a *regime* rather than a number:
the mean free path must come out microscopic against a 3 m bore, on the **longest** of the three
collision channels. If it did not, `jet.py`'s guiding-centre model would still stand and P2 would
not have been withdrawn, so the ordering of the channels is as load-bearing as the number.

`test_the_generous_channel_really_is_the_longest` is the one that carries the argument: the
headline uses neutral hard-sphere scattering *because* it gives the longest path, and that claim
has to be true or the verdict is rhetoric rather than a bound.
"""

from __future__ import annotations

import math

import pytest

from puffsat import continuum, expansion

# A representative station: the flown bag, hot enough to be well ionised.
RHO = 0.2
TEMP = 20_000.0


def test_the_generous_channel_really_is_the_longest() -> None:
    """The headline uses the neutral path because it is the longest. Verify, do not assert."""
    paths = continuum.mean_free_paths(RHO, TEMP)
    assert paths.lambda_neutral > paths.lambda_ion_neutral
    assert paths.lambda_neutral > paths.lambda_coulomb
    assert paths.lambda_generous == paths.lambda_neutral


def test_the_combined_path_is_shorter_than_any_single_channel() -> None:
    """Collision rates add, so opening a channel can only shorten the path."""
    paths = continuum.mean_free_paths(RHO, TEMP)
    assert paths.lambda_combined < paths.lambda_neutral
    assert paths.lambda_combined < paths.lambda_ion_neutral
    assert paths.lambda_combined < paths.lambda_coulomb


def test_the_neutral_path_scales_inversely_with_density() -> None:
    """`lambda = 1/(sqrt2 n sigma)`: at fixed composition, halving `rho` doubles it.

    Not exact here, because the composition itself shifts with density -- a thinner gas is more
    dissociated and more ionised, which *removes* neutrals and lengthens the path further. So the
    test is that the scaling holds to within that shift and in the right direction.
    """
    lo = continuum.mean_free_paths(0.1, TEMP).lambda_neutral
    hi = continuum.mean_free_paths(0.2, TEMP).lambda_neutral
    assert lo > hi
    assert 1.8 < lo / hi < 4.0


def test_the_coulomb_logarithm_is_floored_rather_than_negative() -> None:
    """Cold dense gas drives the expression negative, where the channel is irrelevant anyway."""
    assert continuum.coulomb_logarithm(1e28, 300.0) >= 2.0
    assert continuum.coulomb_logarithm(1e20, 20_000.0) > 5.0


# ---- The verdict itself --------------------------------------------------------------------------


def test_the_bore_is_deeply_continuum_at_every_flown_state() -> None:
    """`Kn` must be microscopic at every plume state, or `mu` survives and P2 stands.

    The threshold is not marginal: a guiding-centre treatment wants `Kn` of order 0.1, and the
    flown states sit six or more orders below it.
    """
    for _speed, temp_0 in expansion.PLUME_STATES:
        kn = continuum.knudsen(expansion.BAG_RHO, temp_0)
        assert kn < 1e-5, f"Kn = {kn:.2e} at T0 = {temp_0} K is not obviously continuum"


def test_a_parcel_collides_astronomically_many_times_in_transit() -> None:
    """One collision destroys `mu`. The point is that there are millions, not that there is one."""
    n = continuum.collisions_per_transit(0.025, 5000.0, 0.0, 2.0e-3)
    assert n > 1e5


def test_the_mean_particle_mass_lies_between_the_atomic_and_molecular_limits() -> None:
    """A sanity floor on the composition: water is 18 amu whole and 6 amu fully dissociated.

    Ionisation drives it lower still by adding electrons, so only the upper bound is hard.
    """
    m_bar = continuum._mean_particle_mass(RHO, TEMP)
    assert m_bar < 18.0 * 1.66053906660e-27
    assert m_bar > 0.0


def test_knudsen_is_the_path_over_the_stated_length() -> None:
    """The definition, so a future refactor cannot silently change what is divided by what."""
    paths = continuum.mean_free_paths(RHO, TEMP)
    assert continuum.knudsen(RHO, TEMP, 1.0) == pytest.approx(paths.lambda_generous)
    assert continuum.knudsen(RHO, TEMP, 2.0) == pytest.approx(0.5 * paths.lambda_generous)


def test_an_infinitely_thin_gas_would_be_collisionless() -> None:
    """The control: the verdict is a statement about *these* densities, not about the formula."""
    assert math.isinf(continuum.mean_free_paths(1e-30, TEMP).lambda_generous) or (
        continuum.knudsen(1e-30, TEMP) > continuum.COLLISIONLESS_KN
    )
