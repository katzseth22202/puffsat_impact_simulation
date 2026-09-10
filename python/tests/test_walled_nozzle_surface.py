"""N16's grid: the geometry it assumes, the shortcut it takes, and the surface it reports."""

from __future__ import annotations

import math

import pytest

from puffsat import eos_methane, eos_water, expansion
from puffsat.walled_nozzle import chamber, propellants, surface

# ---- Geometry ------------------------------------------------------------------------------------


def test_nozzle_length_is_the_cone_identity() -> None:
    """`L = (R_bore - r*) / tan(theta)` -- the whole geometric content of the clock."""
    for area in surface.THROATS:
        r_star = math.sqrt(area / math.pi)
        expected = (chamber.BORE_RADIUS - r_star) / math.tan(
            math.radians(surface.CONE_HALF_ANGLE_DEG)
        )
        assert surface.nozzle_length(area) == pytest.approx(expected, rel=1e-12)


def test_the_recommended_throats_all_get_the_same_length_of_nozzle() -> None:
    """What makes the deep half of the throat sweep comparable: below 0.5 m^2 the throat radius is
    small against the 3 m bore, so those columns differ only in where the cone starts.

    The flown 7 m^2 throat deliberately does *not* satisfy this -- a 1.5 m throat radius is half
    the bore, so it is a 5.6 m nozzle against the deep end's 10.7 m. The grid reports that rather
    than hiding it behind a fixed length.
    """
    recommended = [surface.nozzle_length(a) for a in surface.THROATS if a <= 0.5]
    assert max(recommended) / min(recommended) < 1.11
    assert surface.nozzle_length(7.0) < 0.6 * surface.nozzle_length(0.05)


def test_station_x_is_zero_at_the_throat_and_the_full_length_at_the_bore() -> None:
    for area in (7.0, 0.1):
        assert surface.station_x(area, area) == pytest.approx(0.0, abs=1e-12)
        assert surface.station_x(surface.BORE_AREA, area) == pytest.approx(
            surface.nozzle_length(area), rel=1e-12
        )


# ---- The shortcut the module rests on ------------------------------------------------------------


def test_one_deep_isentrope_reproduces_a_dedicated_run_at_each_throat() -> None:
    """**The module's central claim**: the throat is not a third dimension.

    A point read off one deep expansion must match an expansion run directly to that exit ratio,
    because both are the same adiabat cut in the same place. If this drifts, every throat column
    in the grid is quietly interpolated rather than solved.
    """
    temp_c, volume = 8000.0, 100.0
    _k, rho_c, _u = propellants.solve_slug_ratio(eos_methane.pressure_energy, temp_c, volume)
    deep = surface.isentrope(surface.METHANE, rho_c, temp_c, surface.BORE_AREA / 0.05)

    for area in (2.0, 0.2):
        ratio = surface.BORE_AREA / area
        direct = expansion.cooling_history(
            rho_c,
            temp_c,
            eos_methane.pressure_energy,
            eos_methane.sound_speed,
            ratio,
            1.0,
            steps=220,
            expansion_ratio=4.0e5,
        )[-1]
        read = surface._interpolate(deep, ratio)
        assert read.temp == pytest.approx(direct.temp, rel=2e-3)
        assert read.rho == pytest.approx(direct.rho, rel=5e-3)
        assert read.speed == pytest.approx(direct.speed, rel=2e-3)


def test_conversion_cannot_exceed_the_energy_it_is_measured_against() -> None:
    """`u_e^2/2 <= u` and `u_e <= w/sqrt(1+k)` -- ADR-0016's own identity, both ends."""
    for point in surface.column(surface.METHANE, 10000.0, 100.0, throats=(7.0, 0.5, 0.05)):
        assert 0.0 < point.conversion < 1.0
        assert point.conversion_capped <= point.conversion + 1e-12
        assert point.exit_speed < point.exhaust_speed_ideal


def test_blowdown_scales_inversely_with_throat_area() -> None:
    """`t ~ M / (rho a A*)`: halve the throat and the chamber takes twice as long to empty."""
    a = surface.blowdown_time(5.0, 4000.0, 100.0, 1.0)
    b = surface.blowdown_time(5.0, 4000.0, 100.0, 0.5)
    assert b == pytest.approx(2.0 * a, rel=1e-12)


# ---- The held surface ----------------------------------------------------------------------------


def test_held_reproduces_the_published_w9_column() -> None:
    """W9's speciation table was computed at 1 kg/m^3 and the paper repository carries it. The
    surface has to return through those points or the correction it offers is not a correction."""
    published = {6000.0: 0.599, 5000.0: 0.397, 4000.0: 0.206, 3500.0: 0.144, 3000.0: 0.118}
    for temp, held in published.items():
        value = eos_methane.bond_energy_held(1.0, temp) / eos_methane.FULL_ATOMIZATION_ENERGY
        assert value == pytest.approx(held, abs=0.002)


def test_the_density_exponent_vanishes_at_both_ends() -> None:
    """Mass action, not a fit: a gas that is fully atomised has nothing left to shift, and one
    that is fully recombined has nothing left to break. The paper side's flat -0.21 cannot be
    right at either end, and this is the shape that says so."""
    hot = surface.held_exponent(surface.METHANE, 0.001, 8000.0)
    cold = surface.held_exponent(surface.METHANE, 0.3, 2500.0)
    middle = surface.held_exponent(surface.METHANE, 0.1, 4000.0)
    assert abs(hot) < 0.01
    assert abs(cold) < 0.05
    assert middle < -0.15


def test_held_falls_with_density_everywhere_it_is_not_saturated() -> None:
    """Compressing an equilibrium mixture drives it back toward the bound side, so `held` cannot
    rise with density. A positive exponent anywhere would mean the EOS is not at equilibrium."""
    for point in surface.held_surface(fluids=(surface.METHANE,)):
        assert point.exponent <= 1e-6


# ---- The mixture ---------------------------------------------------------------------------------


def test_the_mixture_reduces_to_water_when_undiluted() -> None:
    """A composition of two subsystems has to be exact at the endpoints or it is not a mixture."""
    pure = surface.mixture(water_share=1.0)
    p_mix, e_mix = pure.pressure_energy(1.0, 6000.0)
    p_ref, e_ref = eos_water.pressure_energy(1.0, 6000.0)
    assert p_mix == pytest.approx(p_ref, rel=1e-9)
    assert e_mix == pytest.approx(e_ref, rel=1e-9)


def test_the_mixture_is_labelled_an_estimate() -> None:
    """It omits the cross-chemistry between water's oxygen and the added hydrogen, so every row it
    produces has to carry the flag."""
    assert surface.mixture().estimated
    assert not surface.METHANE.estimated


# ---- The incoming-momentum debit -------------------------------------------------------------


def test_eta_jet_is_the_square_root_of_the_conversion_fraction() -> None:
    """The bridge between this repository's convention and the companion's.

    `(1+k) u = w^2/2` makes `w^2/(1+k) = 2u`, so `u_e^2 / (w/sqrt(1+k))^2 = u_e^2/(2u)`, which is
    the conversion fraction. It means the companion can take `eta_jet` straight out of the grid
    without re-deriving anything, and it means a conversion fraction quoted as an `eta_jet` would
    be wrong by a square root.
    """
    for point in surface.column(surface.METHANE, 10000.0, 100.0, throats=(7.0, 0.5, 0.05)):
        assert point.eta_jet == pytest.approx(math.sqrt(point.conversion_capped), rel=1e-4)


def test_the_effective_isp_is_the_companion_head_on_form() -> None:
    """`Isp_eff = w (eta_jet sqrt(1+k) - 1) / (k g0)`, the paper's own head-on expression.

    Asserted against the independently-built `isp_carried_gross - momentum_debit` so that the two
    routes to it cannot drift apart.
    """
    for point in surface.column(surface.METHANE, 12000.0, 100.0, throats=(2.0, 0.2)):
        paper = (
            chamber.CLOSING_SPEED
            * (point.eta_jet * math.sqrt(1.0 + point.slug_ratio) - 1.0)
            / (point.slug_ratio * propellants.G0)
        )
        assert point.isp_effective == pytest.approx(paper, rel=1e-4)
        assert point.isp_effective == pytest.approx(
            point.isp_carried_gross - point.momentum_debit, rel=1e-12
        )


def test_the_debit_is_the_same_momentum_for_every_configuration() -> None:
    """`w/(k g0)`: one projectile's momentum, spread over whatever slug was carried. It is why the
    debit falls hardest on a lean charge -- nothing about the nozzle changes it."""
    for point in surface.column(surface.METHANE, 10000.0, 100.0, throats=(7.0, 0.05)):
        assert point.momentum_debit * point.slug_ratio * propellants.G0 == pytest.approx(
            chamber.CLOSING_SPEED, rel=1e-9
        )


def test_the_head_on_burn_is_a_debit_and_never_a_credit() -> None:
    """**The sign that must not be inherited from the overtake study.** A head-on arrival opposes
    the ship's motion, so the corrected figure is below the credit-only one -- and, because the
    debit `w/(k g0)` beats the credit `u_e/(k g0)` whenever `u_e < w`, below the real Isp too."""
    for point in surface.surface(
        fluids=(surface.METHANE, surface.WATER),
        temperatures=(8000.0,),
        volumes=(100.0,),
        throats=(0.5,),
    ):
        assert point.isp_effective < point.isp_carried_gross
        assert point.isp_effective < point.isp_true
        assert point.exit_speed_capped < chamber.CLOSING_SPEED


def test_every_cell_clears_the_thrust_floor() -> None:
    """Below `eta_jet = 1/sqrt(1+k)` the nozzle pushes the ship backwards and the cell returns
    nothing. Nothing in this grid is close, but the floor has to be checked rather than assumed."""
    for point in surface.surface(
        fluids=(surface.METHANE, surface.WATER),
        temperatures=(6000.0, 12000.0),
        volumes=(100.0,),
        throats=(7.0, 0.05),
    ):
        assert point.eta_jet > point.thrust_floor
        assert not point.isp.pushes_backwards
        assert point.isp_effective > 0.0


def test_the_grid_and_the_ladder_share_one_isp_ledger() -> None:
    """**The reason the sign cannot be fixed in one study and missed in the other.** A grid cell
    and a propellant rung at the same exhaust speed and slug ratio must agree exactly, because
    both resolve to `chamber.IspLedger` rather than re-deriving the algebra."""
    point = next(iter(surface.column(surface.METHANE, 10000.0, 200.0, throats=(7.0,))))
    rung = propellants.Rung(
        name="probe",
        mean_atomised_mass=3.2,
        slug_ratio=point.slug_ratio,
        rho=point.rho_c,
        energy=point.energy_c,
        exit_temp=point.exit_temp,
        exhaust_speed=point.exit_speed_capped,
        conversion=point.conversion_capped,
        estimated=True,
    )
    assert rung.isp_effective == pytest.approx(point.isp_effective, rel=1e-12)
    assert rung.isp_true == pytest.approx(point.isp_true, rel=1e-12)
    assert rung.momentum_debit == pytest.approx(point.momentum_debit, rel=1e-12)
