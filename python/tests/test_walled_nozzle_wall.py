"""N9's load case: the front, the gate, and the three loads the wall can be charged with."""

from __future__ import annotations

import math

import pytest

from puffsat import coupling, eos_methane, expansion, front
from puffsat.walled_nozzle import surface, wall

# ---- The integrator, pinned against the one it mirrors -------------------------------------------


def test_the_front_integrator_reproduces_front_integrate() -> None:
    """**The consistency test the module's first claim rests on.**

    `wall.integrate_front` is `coupling.snowplow`'s system with a clock and a general EOS bolted
    on. Run on the *magnetic* geometry, where `front.integrate` is the established answer, the two
    have to put contact at the same station -- otherwise the walled numbers are coming from a
    different model wearing the same name.
    """
    closure = front.spread_speed_table(
        [front.PAPER_SPREAD_M_S * f for f in (0.02, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0)]
    )
    reference = front.integrate(
        45.58e3,
        1.0,
        wall_m=coupling.BORE_RADIUS,
        spread_of_speed=closure,
        bag_radius=coupling.BORE_RADIUS,
        steps=4000,
    )
    _stations, contact, _exit, _slope = wall.integrate_front(
        lambda v: (0.0, closure(v)),
        coupling.BAG_RHO,
        coupling.BORE_LENGTH,
        closing_speed=45.58e3,
        wall_radius=coupling.BORE_RADIUS,
        steps=4000,
        sample_every=10,
    )
    assert reference.contact_x_m is not None
    assert contact is not None
    assert contact.x == pytest.approx(reference.contact_x_m, abs=0.05)


def test_the_clock_is_bounded_by_the_free_flight_time() -> None:
    """The front only ever slows, so its transit cannot be shorter than the undecelerated one."""
    run = wall.run_front(surface.METHANE, 400.0, 10000.0)
    assert run.exit_state.time > run.length / run.closing_speed


def test_impactor_radius_is_the_sphere_identity() -> None:
    r = wall.impactor_radius(25.0, wall.ICE_DENSITY)
    assert (4.0 / 3.0) * math.pi * r**3 * wall.ICE_DENSITY == pytest.approx(25.0, rel=1e-12)


# ---- The gate ------------------------------------------------------------------------------------


def test_the_short_chamber_lets_the_front_out_and_the_long_one_does_not() -> None:
    """The gate N9 asks about, at its two ends. Below it items 1-4 do not arise."""
    assert not wall.run_front(surface.METHANE, 50.0, 10000.0).reaches_wall
    assert wall.run_front(surface.METHANE, 400.0, 10000.0).reaches_wall


def test_the_gate_volume_is_the_crossing_it_claims_to_be() -> None:
    """Bisection is only meaningful if the predicate really flips there."""
    gate = wall.gate_volume(surface.METHANE, 10000.0, 1.0, lo=100.0, hi=400.0, tol=2.0)
    assert not wall.run_front(surface.METHANE, gate * 0.9, 10000.0, steps=4000).reaches_wall
    assert wall.run_front(surface.METHANE, gate * 1.1, 10000.0, steps=4000).reaches_wall


def test_a_faster_spreading_front_touches_sooner() -> None:
    """Monotone in the bracket, which is what makes reporting only its ends honest."""
    slow = wall.gate_volume(surface.METHANE, 10000.0, 1.0, lo=50.0, hi=400.0, tol=4.0)
    fast = wall.gate_volume(surface.METHANE, 10000.0, 1.9, lo=50.0, hi=400.0, tol=4.0)
    assert fast < slow


# ---- The strike ----------------------------------------------------------------------------------


def test_no_contact_means_no_strike() -> None:
    """`None` is the gate's answer, not a failure to compute one."""
    assert wall.strike(wall.run_front(surface.METHANE, 50.0, 10000.0)) is None


def test_the_radiative_fluence_never_exceeds_the_blackbody_cap() -> None:
    """Nothing radiates faster than a blackbody at its own temperature, whatever the opacity says.
    This cap is what makes the item-3 verdict independent of the opacity model."""
    for scale in wall.OPACITY_BRACKET:
        s = wall.strike(wall.run_front(surface.METHANE, 400.0, 10000.0), scale)
        assert s is not None
        assert s.fluence_radiative * 1e6 <= s.blackbody_flux * s.transient * (1.0 + 1e-9)
        assert s.fluence_radiative <= s.fluence_ceiling


def test_more_opacity_delivers_less_once_the_layer_is_thick() -> None:
    """In the diffusion regime photons random-walk out, so opacity and delivery move opposite
    ways -- which is why the free-free *floor* is the conservative opacity for this question."""
    thin = wall.strike(wall.run_front(surface.METHANE, 673.0, 10000.0), 1.0)
    thick = wall.strike(wall.run_front(surface.METHANE, 673.0, 10000.0), 10.0)
    assert thin is not None and thick is not None
    assert thick.optical_depth > thin.optical_depth
    assert thick.fluence_radiative < thin.fluence_radiative


def test_the_opacity_floor_is_at_least_thomson() -> None:
    """Free-free plus Thomson can fall to Thomson but never below it."""
    kappa, _planck = wall.plasma_opacity("methane", 10.0, 200000.0)
    comp = eos_methane.composition(10.0, 200000.0)
    assert kappa >= wall.SIGMA_THOMSON * comp.n_e / 10.0


# ---- Carbon --------------------------------------------------------------------------------------


def test_graphite_vapour_pressure_returns_its_anchor() -> None:
    """Clausius-Clapeyron through the measured 1 atm sublimation point."""
    assert wall.graphite_vapour_pressure(wall.GRAPHITE_1ATM_K) == pytest.approx(wall.ATM, rel=1e-12)


def test_the_deposition_threshold_is_where_the_two_pressures_cross() -> None:
    """The verdict is a crossing, so it has to hold as one."""
    rho, temp = 3.4, 8700.0
    threshold = wall.deposition_threshold(rho, temp)
    assert wall.graphite_vapour_pressure(threshold) == pytest.approx(
        wall.carbon_partial_pressure(rho, temp), rel=0.02
    )


def test_a_wall_above_the_threshold_plates_nothing() -> None:
    rho, temp = 3.4, 8700.0
    threshold = wall.deposition_threshold(rho, temp)
    assert wall.hertz_knudsen_thickness(rho, temp, threshold + 200.0, 0.4) == 0.0
    assert wall.hertz_knudsen_thickness(rho, temp, threshold - 500.0, 0.4) > 0.0


# ---- Convective flux -----------------------------------------------------------------------------


def test_bartz_scales_with_the_area_ratio_to_the_nine_tenths() -> None:
    """The part of the correlation that does not depend on its calibration, and the part the
    narrow-throat conclusion actually rests on."""
    throat = wall.bartz_flux(surface.METHANE, 200.0, 10000.0, 7.0, "throat")
    bore = wall.bartz_flux(surface.METHANE, 200.0, 10000.0, 7.0, "bore")
    ratio = surface.BORE_AREA / 7.0
    assert throat.flux_equilibrium / bore.flux_equilibrium == pytest.approx(ratio**0.9, rel=1e-9)


def test_equilibrium_heat_capacity_exceeds_the_frozen_one() -> None:
    """A dissociating gas absorbs heat into chemistry, so its equilibrium `c_p` has to be the
    larger of the two -- and the factor between them is the width of the Bartz answer."""
    point = wall.bartz_flux(surface.METHANE, 200.0, 10000.0, 7.0, "throat")
    assert point.cp_equilibrium > point.cp_frozen


def test_radiative_equilibrium_temperature_inverts_stefan_boltzmann() -> None:
    temp = 4000.0
    flux = expansion.SIGMA_SB * temp**4
    assert wall.radiative_equilibrium_temperature(flux) == pytest.approx(temp, rel=1e-12)


def test_ablation_depth_is_the_energy_identity() -> None:
    """Fluence over (heat of ablation x density), and nothing else."""
    depth = wall.ablation_depth(1.0e9)
    assert depth * wall.GRAPHITE_DENSITY * wall.GRAPHITE_ABLATION == pytest.approx(1.0e9, rel=1e-12)


# ---- Does the throat self-heal? -----------------------------------------------------------


def test_the_wall_adjacent_gas_is_a_third_state() -> None:
    """The whole item-5 correction in one assertion: the gas touching the wall is at the free
    stream's *pressure* and the wall's *temperature*, which is neither end's density."""
    pressure, wall_temp = 623e5, 3900.0
    rho_w = wall.wall_gas_density(pressure, wall_temp)
    assert eos_methane.pressure_energy(rho_w, wall_temp)[0] == pytest.approx(pressure, rel=1e-6)
    assert rho_w > 10.0  # the throat free stream is near 4.4 kg/m^3; the wall side is far denser


def test_the_equilibrium_gas_is_undersaturated_where_graphite_can_survive() -> None:
    """The finding. Between the cold edge and graphite's own working ceiling the acetylene sink
    holds the monatomic-carbon activity below graphite's, so nothing plates."""
    for temp in (2500.0, 3000.0, 3900.0, 4500.0):
        assert wall.equilibrium_saturation(623e5, temp) < 1.0


def test_a_cold_enough_wall_does_plate() -> None:
    """Solid carbon is the stable phase down there, so the band has to close at the cold end or
    the calculation is not reproducing ordinary thermodynamics."""
    cold, _hot = wall.deposition_window(623e5)
    assert 1500.0 < cold < 2500.0
    assert wall.equilibrium_saturation(623e5, cold - 200.0) > 1.0
    assert wall.equilibrium_saturation(623e5, cold + 200.0) < 1.0


def test_the_two_saturation_limits_disagree_by_orders_of_magnitude() -> None:
    """Which limit the boundary layer sits in is the whole question, so the two edges must be far
    enough apart that picking between them matters."""
    healing = wall.self_healing(100.0, 8000.0, 0.5)
    assert healing.saturation_frozen / healing.saturation_equilibrium > 1e3


def test_the_boundary_layer_thickness_is_free_of_the_heat_capacity_choice() -> None:
    """`delta = k/h` with `k ~ mu c_p` and Bartz's `h ~ c_p`, so `c_p` cancels -- which is why the
    factor-of-two frozen/equilibrium spread in flux does not reach this answer."""
    point = wall.bartz_flux(surface.METHANE, 100.0, 8000.0, 0.5, "throat")
    h_frozen = point.flux_frozen / (8000.0 - wall.WALL_TEMPERATURE)
    k_frozen = point.viscosity * point.cp_frozen / point.prandtl
    layer = wall.boundary_layer(point, 4.4, 8000.0)
    assert k_frozen / h_frozen == pytest.approx(layer.thickness, rel=1e-9)


def test_the_boundary_layer_has_room_to_spare_for_its_chemistry() -> None:
    """The margin is reported in decades below a gas-kinetic collision frequency, so the verdict
    says what would have to be true to flip it rather than resting on a rate nobody has."""
    layer = wall.self_healing(100.0, 8000.0, 0.5).layer
    assert layer.decades_of_margin > 3.0
    assert layer.equilibrium


def test_the_throat_does_not_self_heal() -> None:
    """Item 5's answer, at every throat the ask recommends."""
    for area in (2.0, 0.5, 0.2, 0.1):
        healing = wall.self_healing(100.0, 8000.0, area)
        assert not healing.plates
        assert not healing.self_heals
        assert healing.ablation > 0.0


# ---- The contraction, which is where W10 moves the problem ----------------------------------


def test_the_bulk_end_wall_load_does_not_depend_on_chamber_volume() -> None:
    """The reason shrinking the chamber does not make the strike go away: the same pulse energy
    turns the same corner through the same bore area however long the column in front of it was."""
    small = wall.end_wall_strike(surface.METHANE, 50.0)
    large = wall.end_wall_strike(surface.METHANE, 150.0)
    assert small.incident_bulk == pytest.approx(large.incident_bulk, rel=1e-12)


def test_the_arrival_load_gets_worse_as_the_chamber_shrinks() -> None:
    """A shorter column decelerates the front less, so it arrives faster and on a smaller patch.
    This is the term that runs the opposite way from the side-wall strike."""
    strikes = [wall.end_wall_strike(surface.METHANE, v) for v in (50.0, 100.0, 150.0)]
    assert strikes[0].incident_arrival > strikes[1].incident_arrival > strikes[2].incident_arrival
    assert strikes[0].exit_speed > strikes[2].exit_speed


def test_the_wetted_patch_is_capped_at_the_bore() -> None:
    """Above the gate the cone has already reached the wall, so it cannot wet more than the bore."""
    above = wall.end_wall_strike(surface.METHANE, 400.0)
    assert above.wetted_area == pytest.approx(surface.BORE_AREA, rel=1e-9)


def test_the_front_never_carries_more_than_the_pulse() -> None:
    """Energy conservation on the snowplow: it can only ever hold a share of what arrived."""
    for volume in (50.0, 100.0, 200.0):
        assert 0.0 < wall.end_wall_strike(surface.METHANE, volume).front_share < 1.0


# ---- W24: the throat as a consumable ----------------------------------------------------------


def test_throat_life_follows_the_analytic_area_exponent() -> None:
    """**The strongest check on W24, because the exponent is derived before the solver runs.**

    Bartz gives the throat flux as `D*^-0.2 ~ A*^-0.1`; the blowdown goes as `1/A*`; so the pulse
    fluence and the recession go as `A*^-1.1`. Reaching a *fractional* area growth needs a
    recession proportional to `r* ~ A*^0.5`, so the pulse count goes as `A*^1.6`. If the joined
    Bartz/blowdown/ablation chain ever stops reproducing that, one of the three has drifted.
    """
    lives = {a: wall.throat_life(surface.METHANE, 200.0, 10000.0, a) for a in (7.0, 2.0, 1.0, 0.5)}
    for wide, narrow in ((7.0, 2.0), (2.0, 1.0), (1.0, 0.5)):
        predicted = (wide / narrow) ** wall.THROAT_LIFE_EXPONENT
        actual = lives[wide].pulses_to_area_growth(0.10) / lives[narrow].pulses_to_area_growth(0.10)
        assert actual == pytest.approx(predicted, rel=0.02)


def test_the_cost_of_a_narrow_throat_is_dwell_and_not_flux() -> None:
    """**W24's mechanism, asserted so the write-up cannot be read the other way round.** Halving
    the throat raises the throat flux by only about 7% (`D*^-0.2`) while doubling the blowdown, so
    the fluence -- and the recession with it -- is carried by the time, not the intensity."""
    wide = wall.throat_life(surface.METHANE, 200.0, 10000.0, 2.0)
    narrow = wall.throat_life(surface.METHANE, 200.0, 10000.0, 1.0)
    assert narrow.flux / wide.flux == pytest.approx(2.0**0.1, rel=0.02)
    assert narrow.blowdown / wide.blowdown == pytest.approx(2.0, rel=0.02)
    assert narrow.fluence / wide.fluence == pytest.approx(2.0**1.1, rel=0.03)


def test_the_throat_opens_as_it_erodes() -> None:
    """Recession is radial, so a throat narrowed for performance drifts back up the trade curve.

    **The doubling count is exact in the radius rather than linearised in `dA/A`, and the two
    disagree by 18%.** Dividing a 100% area growth by the per-pulse `dA/A` would give 10x the
    +10% count; the honest answer is `(sqrt(2) - 1)/(sqrt(1.1) - 1)` = 8.49x, because the radius
    has further to go than a linear reading of the area suggests. Linearising would overstate
    every doubling life in W24's table.
    """
    life = wall.throat_life(surface.METHANE, 200.0, 10000.0, 1.0)
    assert life.area_growth_per_pulse == pytest.approx(2.0 * life.recession / life.throat_radius)
    assert life.throat_radius == pytest.approx(math.sqrt(1.0 / math.pi))
    exact = (math.sqrt(2.0) - 1.0) / (math.sqrt(1.1) - 1.0)
    ratio = life.pulses_to_area_growth(1.00) / life.pulses_to_area_growth(0.10)
    assert ratio == pytest.approx(exact, rel=1e-9)
    assert ratio == pytest.approx(8.49, abs=0.01)
    assert ratio < 1.00 / 0.10  # the linearised reading, which is optimistic


def test_the_deep_throat_stops_paying_before_it_stops_costing() -> None:
    """**The finding W24 turns on.** Raw conversion keeps climbing all the way down the ladder,
    but *capped* conversion saturates -- past the freeze station the equilibrium branch releases
    chemistry the real flow no longer has time to release. So below about 0.5 m^2 the effective
    Isp is flat to within a percent while the recession per pulse is still doubling per halving.
    """
    grid = {
        p.throat_area: p
        for p in surface.column(surface.METHANE, 10000.0, 200.0, throats=(0.5, 0.2, 0.1))
    }
    assert grid[0.1].conversion > grid[0.2].conversion > grid[0.5].conversion
    assert grid[0.1].conversion_capped == pytest.approx(grid[0.5].conversion_capped, rel=0.02)
    assert grid[0.1].isp_effective / grid[0.5].isp_effective == pytest.approx(1.0, abs=0.02)
    coarse = wall.throat_life(surface.METHANE, 200.0, 10000.0, 0.5)
    fine = wall.throat_life(surface.METHANE, 200.0, 10000.0, 0.1)
    assert fine.recession > 4.0 * coarse.recession
