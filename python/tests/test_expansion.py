"""Tests for the plume cooling history `T(t)` through the magnetic nozzle.

The companion repo's item 10 needs the field leak weighted by how long the plume spends at each
temperature, which is a quadrature over a cooling history this repo owes. That history is a
quasi-1D steady isentropic expansion on the real water EOS, and the reason it cannot be done with
a constant `gamma` is the reason ADR-0026 exists: **equilibrium recombination hands its energy
back to the thermal pool and buffers the cooling; frozen recombination does not.** The two run as
a bracket, exactly as they do on the plate side.

Every acceptance test below is against a closed form that is known before the solver exists: the
constant-`gamma` adiabat, the isentropic nozzle relations, and the area-Mach law.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import eos_water, expansion

#: A constant-`gamma` ideal gas standing in for the EOS, so the classical relations apply exactly.
#: `R` is water's specific gas constant; the value only has to be self-consistent.
IDEAL_GAMMA = 1.4
IDEAL_R = 461.5


def ideal_eos(rho: float, temp: float) -> tuple[float, float]:
    """`(p, e)` for a calorically perfect gas: `p = rho R T`, `e = R T / (gamma - 1)`."""
    return rho * IDEAL_R * temp, IDEAL_R * temp / (IDEAL_GAMMA - 1.0)


def test_temperature_inverts_the_caloric_equation_of_state() -> None:
    """Slice 1: recover `T` from `(rho, e)` against the closed form `T = e (gamma - 1) / R`.

    The isentrope is integrated in specific energy, so every step has to come back to a
    temperature. For the perfect gas the inverse is exact and independent of the solver:
    at `e = 1.0e6 J/kg`, `T = 1.0e6 * 0.4 / 461.5 = 866.74 K`.
    """
    temp = expansion.temperature_at(0.32, 1.0e6, ideal_eos)
    assert temp == pytest.approx(1.0e6 * (IDEAL_GAMMA - 1.0) / IDEAL_R, rel=1e-9)


def test_adiabat_reproduces_the_constant_gamma_power_law() -> None:
    """Slice 2: `de = -p d(1/rho)` integrated down must give `T/T0 = (rho/rho0)^(gamma - 1)`.

    The classical isentrope is the one closed form the integrator can be held to without any
    plasma chemistry, and it is the reason the isentrope is integrated rather than looked up:
    `eos_water` publishes no entropy function, so the adiabat *is* the definition of the path.

    Expanding by a factor 8 in volume at `gamma = 1.4` cools by `8^-0.4 = 0.43528` exactly.
    """
    states = expansion.expand(0.32, 20000.0, 0.32 / 8.0, ideal_eos, steps=256)

    assert states[0].rho == pytest.approx(0.32)
    assert states[0].temp == pytest.approx(20000.0)
    assert states[-1].rho == pytest.approx(0.04)
    # 2.1e-5 at 256 steps; the rate, not the level, is checked by the next test.
    assert states[-1].temp == pytest.approx(20000.0 * 8.0 ** (1.0 - IDEAL_GAMMA), rel=1e-4)


def test_adiabat_integration_is_second_order_accurate() -> None:
    """Slice 2b: the strongest single correctness signal for the integrator (CLAUDE.md).

    Heun's predictor-corrector on `de/dv = -p` is formally second order, so halving the step must
    quarter the error against the exact power law. Anything materially shallower means the
    corrector is not doing its job.
    """
    exact = 20000.0 * 8.0 ** (1.0 - IDEAL_GAMMA)
    errors = [
        abs(expansion.expand(0.32, 20000.0, 0.04, ideal_eos, steps=n)[-1].temp - exact) / exact
        for n in (16, 32, 64)
    ]
    for coarse, fine in pairwise(errors):
        assert math.log2(coarse / fine) > 1.8


def ideal_sound_speed(rho: float, temp: float) -> float:
    """`c = sqrt(gamma R T)` -- exact for the perfect gas, independent of the EOS callable."""
    return math.sqrt(IDEAL_GAMMA * IDEAL_R * temp)


def area_mach_law(mach: float) -> float:
    """Closed-form `A/A*` for isentropic quasi-1D flow of a perfect gas (any gas-dynamics text)."""
    g = IDEAL_GAMMA
    stagnation = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * mach * mach)
    return float((1.0 / mach) * stagnation ** (0.5 * (g + 1.0) / (g - 1.0)))


def test_nozzle_obeys_the_isentropic_stagnation_and_area_mach_relations() -> None:
    """Slice 3: `u` from energy conservation and `A/A*` from mass conservation, held to the
    textbook perfect-gas relations at *every* point on the path -- not just the endpoints.

    `T0/T = 1 + (gamma-1)/2 M^2` is the independent check on the speed, and the area-Mach law is
    the independent check on the area ratio. Both are closed forms the solver never sees.
    """
    points = expansion.nozzle_history(
        0.32, 20000.0, 0.32 / 64.0, ideal_eos, ideal_sound_speed, steps=256
    )

    assert points[0].speed == pytest.approx(0.0, abs=1e-6)
    for pt in points[1:]:
        mach = pt.speed / ideal_sound_speed(pt.rho, pt.temp)
        stagnation = 1.0 + 0.5 * (IDEAL_GAMMA - 1.0) * mach * mach
        assert 20000.0 / pt.temp == pytest.approx(stagnation, rel=2e-4)
        assert pt.area_ratio == pytest.approx(area_mach_law(mach), rel=2e-3)


def test_nozzle_conserves_total_enthalpy_along_the_whole_path() -> None:
    """Slice 3b: `h + u^2/2 = h0` is the constraint that defines the speed, so violating it means
    the adiabat and the energy equation have come apart. Checked on the *real* water EOS, where
    the chemical energy moves between the thermal and chemical pools as the plume recombines."""
    points = expansion.nozzle_history(
        0.323, 26200.0, 0.323 / 16.0, eos_water.pressure_energy, eos_water.sound_speed, steps=96
    )
    h0 = points[0].enthalpy
    for pt in points:
        assert pt.enthalpy + 0.5 * pt.speed**2 == pytest.approx(h0, rel=1e-6)


def test_frozen_recombination_cools_the_plume_faster_than_equilibrium() -> None:
    """Slice 4: the bracket, and the reason a constant-`gamma` adiabat cannot answer this.

    ADR-0026 states the mechanism for the plate side: the equilibrium EOS *returns* the
    dissociation and ionisation store as the gas re-expands, while a frozen composition carries
    it away as inert enthalpy. The same mechanism runs in the nozzle, so the ordering is fixed
    before either curve is computed -- **equilibrium stays hotter at every density below the
    start**, because recombination keeps refilling the thermal pool it is draining.

    The two EOS agree exactly at the freeze reference state (splice continuity), so the test is
    the *ordering below it*, which is the physical claim.
    """
    rho_0, temp_0 = 0.323, 26200.0
    frozen = eos_water.frozen_composition(rho_0, temp_0)

    equilibrium = expansion.expand(rho_0, temp_0, rho_0 / 32.0, eos_water.pressure_energy, steps=64)
    locked = expansion.expand(
        rho_0,
        temp_0,
        rho_0 / 32.0,
        lambda r, t: eos_water.pressure_energy_frozen(r, t, frozen),
        steps=64,
    )

    assert equilibrium[0].temp == pytest.approx(locked[0].temp)
    for eq, fr in zip(equilibrium[1:], locked[1:], strict=True):
        assert eq.rho == pytest.approx(fr.rho)
        assert eq.temp > fr.temp
    # And the gap is large enough to matter, not a rounding artifact.
    assert equilibrium[-1].temp > 1.5 * locked[-1].temp


def test_cooling_history_spans_the_field_region_and_brackets_its_own_residence_time() -> None:
    """Slice 5: the deliverable. `T(t)` from the sonic throat to the stated exit area ratio.

    The residence time is a quadrature of `dx/u`, and the mean-value theorem bounds it
    independently of how that quadrature is done: the flow only accelerates through a supersonic
    expansion, so `length/u_exit < t_total < length/u_throat`. That bound needs no integrator.
    """
    rows = expansion.cooling_history(
        0.323, 26200.0, ideal_eos, ideal_sound_speed, area_ratio_end=4.0, length=23.8
    )

    assert rows[0].mach == pytest.approx(1.0, rel=1e-3)
    assert rows[0].area_ratio == pytest.approx(1.0, rel=1e-3)
    assert rows[-1].area_ratio == pytest.approx(4.0, rel=1e-3)
    assert rows[0].x == pytest.approx(0.0, abs=1e-9)
    assert rows[-1].x == pytest.approx(23.8)
    assert rows[0].time == pytest.approx(0.0, abs=1e-9)

    for a, b in pairwise(rows):
        assert b.time > a.time
        assert b.x > a.x
        assert b.temp < a.temp
        assert b.speed > a.speed

    assert 23.8 / rows[-1].speed < rows[-1].time < 23.8 / rows[0].speed


def test_residence_time_converges_under_grid_refinement() -> None:
    """Slice 5b: the residence time is what the consumers actually read, so it is the number that
    has to be resolution-independent -- a coarse answer that moves under refinement would put the
    plume on the wrong side of the instability threshold for a numerical reason."""
    times = [
        expansion.cooling_history(
            0.323, 26200.0, ideal_eos, ideal_sound_speed, area_ratio_end=4.0, length=23.8, steps=n
        )[-1].time
        for n in (64, 128, 256)
    ]
    for coarse, fine in pairwise(times):
        assert fine == pytest.approx(coarse, rel=2e-3)


def test_opacity_traps_radiation_when_thick_and_vents_it_when_thin() -> None:
    """Slice 6: the radiation check, held to the *sign* of its two opacity limits.

    This is the non-obvious physical claim, and it is why the check cannot be a single formula:
    raising the opacity makes an emission-limited plume cool *faster* (there is more emitter),
    and a diffusion-limited one cool *slower* (the photons are trapped). A model that got the
    regime wrong would move the cooling time in the wrong direction.
    """
    # Genuinely emission-limited: tau_Planck = 0.07, so every photon made does escape.
    thin = expansion.radiation_check(
        temp=6000.0, rho=0.008, energy=2.0e7, kappa_planck=2.9, kappa_rosseland=0.036, radius=3.0
    )
    assert thin.planck_depth < 1.0
    assert thin.regime == "emission"

    thicker_emitter = expansion.radiation_check(
        temp=6000.0, rho=0.008, energy=2.0e7, kappa_planck=5.8, kappa_rosseland=0.036, radius=3.0
    )
    assert thicker_emitter.cooling_time == pytest.approx(0.5 * thin.cooling_time, rel=1e-9)

    # Deep interior: tau_Rosseland = 450, so the loss is set by how slowly photons diffuse out.
    thick = expansion.radiation_check(
        temp=20000.0, rho=0.3, energy=2.6e8, kappa_planck=1.5e4, kappa_rosseland=500.0, radius=3.0
    )
    assert thick.optical_depth > 1.0
    assert thick.regime == "diffusion"

    more_trapping = expansion.radiation_check(
        temp=20000.0, rho=0.3, energy=2.6e8, kappa_planck=1.5e4, kappa_rosseland=1000.0, radius=3.0
    )
    assert more_trapping.cooling_time == pytest.approx(2.0 * thick.cooling_time, rel=1e-9)


def test_residence_weighted_leak_is_the_time_average_the_quadrature_asks_for() -> None:
    """Slice 7: item 10 is "weight `1/Rm(T)` by how long the plume spends at each `T`".

    Held to two hand-computable cases. A history that leaks at a constant rate must average to
    that rate whatever its timing. A history that leaks 0.02 for the first 1 ms and 0.10 for the
    next 3 ms must average to `(0.02*1 + 0.10*3)/4 = 0.08` -- weighted by *duration*, which is
    the whole point: a naive average over stations would give 0.06 and understate it by 25%.
    """
    constant = [(0.0, 0.05), (0.7e-3, 0.05), (3.0e-3, 0.05)]
    assert expansion.residence_weighted_leak(constant) == pytest.approx(0.05)

    stepped = [(0.0, 0.02), (1.0e-3, 0.02), (1.0e-3, 0.10), (4.0e-3, 0.10)]
    assert expansion.residence_weighted_leak(stepped) == pytest.approx(0.08)


def test_radiative_loss_never_exceeds_free_streaming() -> None:
    """Slice 6b: the invariant that caught the first version of this check.

    No body can lose energy faster than a blackbody at its own temperature radiating from its own
    surface. For a flux tube of radius `R` that caps the volumetric loss at `(2/R) sigma T^4`,
    i.e. floors the cooling time at `rho e R / (2 sigma T^4)`.

    The first version violated it by two to three orders of magnitude near `tau_Rosseland ~ 1`,
    because it picked the regime on the Rosseland depth but computed the emission from the Planck
    mean -- and TOPS puts `kappa_P/kappa_R ~ 100` for water here, so those states are still deeply
    thick to their *own* emission. The symptom was a radiated-energy fraction that grew without
    bound under grid refinement and passed 1, which is how the defect was found.
    """
    # The station that broke it: 45.58 km/s equilibrium at T = 9056 K, where tau_Rosseland is
    # 0.36 (nominally "thin") but tau_Planck is 84 (deeply thick to its own emission).
    rho, radius, energy = 0.0393, 3.7, 5.0e7
    check = expansion.radiation_check(
        temp=9000.0,
        rho=rho,
        energy=energy,
        kappa_planck=575.8,
        kappa_rosseland=2.489,
        radius=radius,
    )
    assert check.optical_depth < 1.0  # Rosseland: nominally thin
    floor = rho * energy * radius / (2.0 * expansion.SIGMA_SB * 9000.0**4)
    assert check.cooling_time >= floor


def test_radiated_fraction_is_the_dwell_weighted_loss_not_a_station_average() -> None:
    """Slice 8: the validity gate on the adiabatic isentrope.

    The fraction of internal energy radiated over the transit is `integral dt / t_rad`, and it is
    the number that says whether the adiabatic history is allowed. Hand case: a plume whose
    radiative cooling time is a steady 10 ms, carried for 2 ms, loses `2/10 = 0.2` of its energy.
    Halving the cooling time for the second half of that transit gives
    `1/10 + 1/5 = 0.3` per the same integral.
    """
    steady = [(0.0, 10.0e-3), (1.0e-3, 10.0e-3), (2.0e-3, 10.0e-3)]
    assert expansion.radiated_fraction(steady) == pytest.approx(0.2)

    faster_late = [(0.0, 10.0e-3), (1.0e-3, 10.0e-3), (1.0e-3, 5.0e-3), (2.0e-3, 5.0e-3)]
    assert expansion.radiated_fraction(faster_late) == pytest.approx(0.3)
