"""N9 items 1-7: does the shocked front reach the wall, and what does it deliver when it does?

The section clears the wall against the *equilibrated* 10 kK chamber. That is not what arrives
first. A 25 kg impactor at 75 km/s drives a snowplow front into the pre-charge, and the layer at
its nose is at a few hundred thousand kelvin, not ten thousand.

**The ask's own instruction is to price the strike before designing any mitigation**, and it names
the gate: "if the front exits first below about 170 m^3, say so and items 1 to 4 close for the
short chamber". So this module answers the gate first, and everything downstream of it is
conditional on the answer.

# The three physical claims this module rests on, stated up front

1. **The front is a snowplow whose cone angle closes on its own current speed.** That is
   `coupling.snowplow`'s system, and `front.py` integrates it for the magnetic nozzle; here the
   same system is integrated on the walled geometry, with the pre-charge density and the chamber
   length of ADR-0016 rather than the bag's. `test_walled_nozzle_wall.py` pins this integrator
   against `front.integrate` on the magnetic geometry so the two cannot drift.

2. **The shocked layer's state is the strong-shock piston solution `e = v^2/2`, inverted on the
   real equation of state.** That is exactly how `sec:needle_through_fog`'s 94 600 K was obtained
   for water at 45.58 km/s (W-series, `front.py`), and it is used here unchanged -- only the fluid
   and the speed differ.

3. **The layer is radiatively thick, so what it delivers is set by diffusion and not by its
   temperature.** This is the finding. A blackbody reading of a 200 kK layer gives 9e13 W/m^2 and
   is meaningless; the layer cannot deliver energy it has not got, and it cannot deliver it faster
   than photons random-walk out of an optical depth of thousands. The opacity used is the
   **free-free floor** -- Kramers plus Thomson, computed on the solved composition -- because at
   200 kK a C/H plasma is dominated by free-free, and because bound-free would only make the layer
   *thicker* and the wall load *smaller*. The bracket is carried explicitly (`OPACITY_BRACKET`),
   and the verdict survives all of it.

# What is not computed here

Ablation chemistry, mechanical response of the liner, and the trajectory of ablated wall material
back into the flow. And the carbon result (item 5) is a **thermodynamic** verdict plus a kinetic
ceiling: it says where deposition is favoured and how fast the arriving flux could plate, not what
a real surface with a boundary layer does with it.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat import eos_methane, eos_water, expansion
from puffsat.walled_nozzle import chamber, propellants, surface

#: Wall radius [m]: the chamber *is* the bore, so there is no bag/liner gap here. In the magnetic
#: nozzle the mist sits inside a bag and the clearance outside it is vacuum (`front.integrate`
#: caps swept area at the bag radius); a walled chamber is filled to the wall, so the two radii
#: coincide and every kilogram between the front and the wall is swept.
WALL_RADIUS = chamber.BORE_RADIUS

#: Density of the ice impactor [kg/m^3], for its radius. The front starts at the projectile's own
#: frontal radius, which is where `front.py` starts it.
ICE_DENSITY = 917.0

#: Spreading-speed multiples the paper brackets the front over (`front.SPREAD_BRACKET`): the
#: shocked layer's own sound speed, and two faster readings the literature supports.
SPREAD_BRACKET = (1.0, 1.6, 1.9)

#: Density jump across the front. Strong-shock ideal; `c_s` is weakly sensitive to it.
SHOCK_COMPRESSION = 4.0

#: Chamber volumes [m^3] the ask now wants run, and the ones it was written for.
VOLUMES = (50.0, 100.0, 200.0, 400.0)
LEGACY_VOLUMES = (200.0, 400.0, 673.0)

#: The paper-side estimate of the volume below which the front exits before touching the wall.
PAPER_GATE_VOLUME = 170.0

DEFAULT_OUTPUT = Path("data/results/walled_nozzle/wall_front.csv")
STRIKE_OUTPUT = Path("data/results/walled_nozzle/wall_strike.csv")
BARTZ_OUTPUT = Path("data/results/walled_nozzle/wall_bartz.csv")
THROAT_LIFE_OUTPUT = Path("data/results/walled_nozzle/wall_throat_life.csv")


def chamber_length(volume: float, radius: float = WALL_RADIUS) -> float:
    """Column length [m] of a cylindrical chamber of this volume at the bore radius."""
    return volume / (math.pi * radius**2)


def impactor_radius(mass: float = chamber.IMPACTOR_MASS, density: float = ICE_DENSITY) -> float:
    """Frontal radius [m] of a spherical ice impactor."""
    return float((3.0 * mass / (4.0 * math.pi * density)) ** (1.0 / 3.0))


# ---- The front (items 1 and 2, and the gate) -----------------------------------------------------


@dataclass(frozen=True)
class FrontState:
    """The snowplow front at one axial station."""

    x: float
    time: float
    radius: float
    speed: float
    swept_mass: float
    #: Local slug ratio the front has accumulated so far -- the number `sec:needle_through_fog`
    #: reasons with when it says "roughly 14 kg against a 25 kg impactor".
    slug_ratio: float
    shocked_temp: float
    shocked_rho: float


def spread_closure(
    fluid: surface.Fluid,
    ambient_rho: float,
    closing_speed: float,
    compression: float = SHOCK_COMPRESSION,
    samples: int = 40,
) -> Callable[[float], tuple[float, float]]:
    """Tabulate `(T_shocked, c_s)` against piston speed and return a fast interpolant.

    The integrator takes tens of thousands of steps and the EOS inversion is the hot spot, so the
    relation -- smooth and monotone in `v` -- is tabulated on a log grid once. This is
    `front.spread_speed_table` generalised off `eos_water`, and it returns the temperature too
    because item 2 asks for it at contact.
    """
    rho = ambient_rho * compression
    grid = [closing_speed * f for f in np.geomspace(0.01, 1.0, samples)]
    temps = [expansion.temperature_at(rho, 0.5 * v * v, fluid.pressure_energy) for v in grid]
    speeds = [fluid.sound_speed(rho, t) for t in temps]

    def probe(v: float) -> tuple[float, float]:
        return float(np.interp(v, grid, temps)), float(np.interp(v, grid, speeds))

    return probe


@dataclass(frozen=True)
class FrontRun:
    """One integrated front, with the gate called."""

    fluid: str
    volume: float
    temp_c: float
    ambient_rho: float
    length: float
    closing_speed: float
    spread_multiple: float
    stations: tuple[FrontState, ...]
    #: The station where the front first reaches the wall, or `None` if it leaves first.
    contact: FrontState | None
    #: The state at the chamber's exit plane, whether or not contact happened.
    exit_state: FrontState
    #: Straight-cone contact station [m] at the entry angle -- the construction the ask's own
    #: "6 m of the column" comes from, kept so the correction is visible.
    cone_contact_x: float
    #: Half-angle [deg] of that straight cone.
    cone_half_angle: float

    @property
    def reaches_wall(self) -> bool:
        """Whether items 1-4 are live for this chamber at all."""
        return self.contact is not None


def integrate_front(
    probe: Callable[[float], tuple[float, float]],
    ambient_rho: float,
    length: float,
    *,
    spread_multiple: float = 1.0,
    wall_radius: float = WALL_RADIUS,
    closing_speed: float = chamber.CLOSING_SPEED,
    projectile_mass: float = chamber.IMPACTOR_MASS,
    front_radius_0: float | None = None,
    compression: float = SHOCK_COMPRESSION,
    steps: int = 20000,
    sample_every: int = 250,
) -> tuple[tuple[FrontState, ...], FrontState | None, FrontState, float]:
    """`coupling.snowplow`'s system, carrying the clock and the shocked state.

        m v = m_0 w            momentum, the swept gas moving with the front
        dm/dx = rho pi min(r, R_wall)^2
        dr/dx = c_exp / v      the cone closes on the *current* speed, not the entry speed
        dt/dx = 1 / v

    **`front.integrate` is the same system** and is pinned against this one in the tests; it is
    not reused directly because it caps swept area at a bag radius inside a larger liner, carries
    no clock, and knows only `eos_water`. All three differ here.
    """
    r0 = impactor_radius(projectile_mass) if front_radius_0 is None else front_radius_0
    mass, speed, radius, x, time = projectile_mass, closing_speed, r0, 0.0, 0.0
    dx = length / steps
    contact: FrontState | None = None
    stations: list[FrontState] = []

    def state() -> FrontState:
        temp, _cs = probe(speed)
        rho_s = ambient_rho * compression
        return FrontState(
            x=x,
            time=time,
            radius=radius,
            speed=speed,
            swept_mass=mass - projectile_mass,
            slug_ratio=(mass - projectile_mass) / projectile_mass,
            shocked_temp=temp,
            shocked_rho=rho_s,
        )

    for i in range(steps + 1):
        if i % sample_every == 0 or (contact is None and radius >= wall_radius):
            stations.append(state())
        if contact is None and radius >= wall_radius:
            contact = stations[-1]
        if i == steps:
            break
        mass += ambient_rho * math.pi * min(radius, wall_radius) ** 2 * dx
        speed = projectile_mass * closing_speed / mass
        _t, c_s = probe(speed)
        radius += spread_multiple * c_s / speed * dx
        time += dx / speed
        x += dx

    exit_state = stations[-1]
    _t0, c_entry = probe(closing_speed)
    cone_slope = spread_multiple * c_entry / closing_speed
    return tuple(stations), contact, exit_state, cone_slope


def run_front(
    fluid: surface.Fluid,
    volume: float,
    temp_c: float = 10000.0,
    spread_multiple: float = 1.0,
    closing_speed: float = chamber.CLOSING_SPEED,
    steps: int = 20000,
) -> FrontRun:
    """Solve the chamber charge, then drive the front into it."""
    k, _rho_total, _u = propellants.solve_slug_ratio(fluid.pressure_energy, temp_c, volume)
    # The front ploughs into the *pre-charge*, which is the slug alone -- the impactor is the
    # piston, not part of the target.
    ambient_rho = k * chamber.IMPACTOR_MASS / volume
    length = chamber_length(volume)
    probe = spread_closure(fluid, ambient_rho, closing_speed)
    stations, contact, exit_state, cone_slope = integrate_front(
        probe,
        ambient_rho,
        length,
        spread_multiple=spread_multiple,
        closing_speed=closing_speed,
        steps=steps,
    )
    return FrontRun(
        fluid=fluid.name,
        volume=volume,
        temp_c=temp_c,
        ambient_rho=ambient_rho,
        length=length,
        closing_speed=closing_speed,
        spread_multiple=spread_multiple,
        stations=stations,
        contact=contact,
        exit_state=exit_state,
        cone_contact_x=(WALL_RADIUS - impactor_radius()) / cone_slope,
        cone_half_angle=math.degrees(math.atan(cone_slope)),
    )


def gate_volume(
    fluid: surface.Fluid,
    temp_c: float = 10000.0,
    spread_multiple: float = 1.0,
    lo: float = 5.0,
    hi: float = 2000.0,
    tol: float = 0.5,
) -> float:
    """Chamber volume [m^3] at which the front just reaches the wall as it leaves.

    **This is the number the ask asks for.** Below it the front exits the column before its cone
    touches the liner and items 1-4 close; above it they are live. Bisection on `reaches_wall`,
    which is monotone in volume: a longer column at lower density both gives the cone more room
    and slows the front's deceleration less.
    """
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if run_front(fluid, mid, temp_c, spread_multiple, steps=4000).reaches_wall:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ---- The strike (item 3) -------------------------------------------------------------------------

#: Thomson cross-section [m^2].
SIGMA_THOMSON = 6.6524587e-29

#: Kramers free-free coefficient in the standard astrophysical form
#: `kappa_ff = C * (Zbar^2/Abar) * (1 + X) * rho * T^-3.5`, with `rho` in kg/m^3 and `kappa` in
#: m^2/kg. The cgs coefficient is 3.68e22 cm^2 g^-1 for `rho` in g/cm^3; converting
#: (1 cm^2/g = 0.1 m^2/kg, 1 g/cm^3 = 1e3 kg/m^3) gives 3.68e22 * 0.1 / 1e3 = 3.68e18.
#: Gaunt factor taken as 1, which is right to tens of percent in this regime.
KRAMERS_FF = 3.68e18

#: Factor the opacity is bracketed over. The verdict has to survive all of it, and it does: an
#: optical depth of thousands stays in the diffusion regime even divided by ten.
OPACITY_BRACKET = (0.1, 1.0, 10.0)


def plasma_opacity(fluid_name: str, rho: float, temp: float) -> tuple[float, float]:
    """`(kappa_Rosseland, kappa_Planck)` [m^2/kg] of the shocked layer, on the free-free floor.

    **A floor, deliberately.** Free-free plus Thomson is what a *fully ionised* plasma cannot go
    below; the real gas at 200 kK still holds bound electrons on carbon (mean charge ~3.2 of 6),
    so bound-free and line opacity add to this. Every omitted term makes the layer thicker, and a
    thicker layer delivers *less* to the wall -- so the floor is the conservative choice for a
    question about how hard the wall is hit.

    The Planck mean is returned equal to the Rosseland mean here. For free-free the two differ by
    a factor of order unity (the Kramers Planck/Rosseland ratio is about 4), and `radiation_check`
    only uses the Planck mean to decide whether the emission limit binds -- which at these depths
    it does not by three orders of magnitude.
    """
    if fluid_name == "methane":
        methane = eos_methane.composition(rho, temp)
        n_e, n_hp = methane.n_e, methane.n_hp
        heavy = sum((k + 1) ** 2 * n for k, n in enumerate(methane.n_c_ions)) / 12.011
        n_formula = rho / eos_methane.AMU
        x_h = 4.0 * 1.008 / 16.043
    else:
        water = eos_water.composition(rho, temp)
        n_e, n_hp = water.n_e, water.n_hp
        heavy = sum((k + 1) ** 2 * n for k, n in enumerate(water.n_o_ions)) / 15.999
        n_formula = rho / eos_water.AMU
        x_h = 2.0 * 1.008 / 18.015
    kappa_es = SIGMA_THOMSON * n_e / rho
    z2_over_a = (n_hp + heavy) / n_formula
    kappa_ff = KRAMERS_FF * z2_over_a * (1.0 + x_h) * rho * temp**-3.5
    return kappa_es + kappa_ff, kappa_es + kappa_ff


#: Stanton-number bracket for the grazing layer's convective delivery. Reynolds analogy on a
#: turbulent boundary layer (`St = C_f/2`, `C_f ~ 0.026 Re^-0.2`) puts a hypersonic wall-bounded
#: flow at Reynolds numbers of 1e6-1e8 in this band. **This is the crudest input in the module and
#: it carries the convective answer**, which is exactly what the ask says of item 7.
STANTON_BRACKET = (0.001, 0.01)

#: `sec:minimum_nozzle`'s equilibrium wall load, for the comparison the ask wants.
EQUILIBRIUM_FLUENCE = 1.6
EQUILIBRIUM_FLUX = 54.8e6


@dataclass(frozen=True)
class Strike:
    """What the front delivers to the wall over the contact transient (item 3)."""

    fluid: str
    volume: float
    temp_c: float
    spread_multiple: float
    contact_x: float
    contact_time: float
    contact_speed: float
    swept_mass: float
    slug_ratio: float
    shocked_temp: float
    shocked_rho: float
    #: Radial speed [m/s] the layer arrives at the wall with -- the front's own spreading speed.
    radial_speed: float
    #: Layer depth [m] the wall looks through: contact station to the chamber exit.
    layer_depth: float
    optical_depth: float
    regime: str
    #: Radiative cooling time [s] of the layer, the binding branch of `radiation_check`.
    cooling_time: float
    #: Blackbody flux `sigma T^4` [W/m^2]. Reported because it is the number that makes the
    #: strike look unsurvivable, and because it is not the one that binds.
    blackbody_flux: float
    #: Duration [s] of the transient: the front's remaining transit to the chamber exit.
    transient: float
    #: Radiative fluence [MJ/m^2] delivered over the transient.
    fluence_radiative: float
    #: Stagnation enthalpy the layer *brings* to the wall [MJ/m^2] -- an incident flux, not a
    #: delivered one. What crosses into the surface is this times a Stanton number.
    incident_convective: float
    #: Convective fluence [MJ/m^2] at the two ends of `STANTON_BRACKET`.
    fluence_convective_low: float
    fluence_convective_high: float
    #: Energy per unit wall area held by the shocked column at contact [MJ/m^2] -- the hard
    #: ceiling on any delivery mechanism, and what a "the front dumps everything" reading gives.
    fluence_ceiling: float

    @property
    def radiative_ratio(self) -> float:
        """Radiative fluence against the equilibrium load the section already clears."""
        return self.fluence_radiative / EQUILIBRIUM_FLUENCE

    @property
    def convective_ratio(self) -> float:
        """Convective fluence at the pessimistic Stanton edge, against the same reference."""
        return self.fluence_convective_high / EQUILIBRIUM_FLUENCE


def strike(run: FrontRun, opacity_scale: float = 1.0) -> Strike | None:
    """Price the contact transient for a front that reaches the wall.

    `None` where it does not -- which is the gate, and is a real answer rather than a failure.

    **Radiation is priced on the layer's own clock, not on its temperature.** The delivered
    fluence is the energy the column holds times the share of it that can leave in the transient,
    `1 - exp(-dt/tau_rad)`, with `tau_rad` the binding branch of `expansion.radiation_check` --
    emission-limited, free-streaming or diffusive, whichever is *slowest*. It is then capped at
    `sigma T^4 dt`, which nothing can exceed whatever the opacity says. **That cap alone settles
    the item at every chamber where contact happens**, which is why the verdict does not depend on
    the opacity model.

    **Convection is priced as an incident flux and a Stanton bracket**, because the gas arriving
    at the wall does not deposit its energy -- it stagnates, and what crosses the surface is set by
    the boundary layer. The incident number is reported alongside so the bracket is visible.
    """
    if run.contact is None:
        return None
    c = run.contact
    depth = max(run.length - c.x, 0.5 * WALL_RADIUS)
    kappa_r, kappa_p = plasma_opacity(run.fluid, c.shocked_rho, c.shocked_temp)
    kappa_r *= opacity_scale
    kappa_p *= opacity_scale
    energy = 0.5 * c.speed**2
    check = expansion.radiation_check(
        c.shocked_temp, c.shocked_rho, energy, kappa_p, kappa_r, depth
    )
    blackbody = expansion.SIGMA_SB * c.shocked_temp**4
    transient = max(run.length - c.x, 0.0) / c.speed
    held = c.shocked_rho * energy * depth  # [J/m^2] in the column the wall looks through
    escaped = held * (1.0 - math.exp(-transient / check.cooling_time))
    radiative = min(escaped, blackbody * transient)

    probe = spread_closure(
        _fluid_by_name(run.fluid), run.ambient_rho, run.closing_speed, SHOCK_COMPRESSION
    )
    radial = run.spread_multiple * probe(c.speed)[1]
    incident = c.shocked_rho * radial * energy * transient
    return Strike(
        fluid=run.fluid,
        volume=run.volume,
        temp_c=run.temp_c,
        spread_multiple=run.spread_multiple,
        contact_x=c.x,
        contact_time=c.time,
        contact_speed=c.speed,
        swept_mass=c.swept_mass,
        slug_ratio=c.slug_ratio,
        shocked_temp=c.shocked_temp,
        shocked_rho=c.shocked_rho,
        radial_speed=radial,
        layer_depth=depth,
        optical_depth=check.optical_depth,
        regime=check.regime,
        cooling_time=check.cooling_time,
        blackbody_flux=blackbody,
        transient=transient,
        fluence_radiative=radiative / 1e6,
        incident_convective=incident / 1e6,
        fluence_convective_low=STANTON_BRACKET[0] * incident / 1e6,
        fluence_convective_high=STANTON_BRACKET[1] * incident / 1e6,
        fluence_ceiling=held / 1e6,
    )


def _fluid_by_name(name: str) -> surface.Fluid:
    """The `surface.Fluid` a run was made with, for re-probing the shock closure."""
    return {"methane": surface.METHANE, "water": surface.WATER}[name]


# ---- The film (items 4 and 6) -------------------------------------------------------------------

#: Areal densities [kg/m^2] the ask asks about for a continuously injected cold film.
FILM_AREAL = (0.02, 0.05, 0.1, 0.2)


def film_capacity(areal_density: float, fluid: surface.Fluid, temp: float = 10000.0) -> float:
    """Energy [J/m^2] a cold film absorbs before it is fully atomised and heated to `temp`.

    The film's job is to be destroyed instead of the liner. Its capacity is the specific energy
    the working fluid holds at the chamber temperature -- atomisation plus sensible heat, which is
    exactly `u` from the chamber solve -- times its areal density. Latent heat of the condensed
    phase is neglected and is worth under a percent against a 100 MJ/kg atomisation.
    """
    return areal_density * fluid.pressure_energy(1.0, temp)[1]


# ---- Throat carbon (item 5) ----------------------------------------------------------------------

#: Sublimation enthalpy of graphite to monatomic carbon vapour [J/mol]. This is the same
#: 711.2 kJ/mol the W9 acetylene accounting uses for condensation, taken from the other direction.
GRAPHITE_SUBLIMATION = 711.2e3
#: Graphite's 1 atm sublimation point [K]. Anchors the Clausius-Clapeyron line so the vapour
#: pressure is a two-parameter fit to a measured point rather than a free extrapolation.
GRAPHITE_1ATM_K = 3915.0
R_GAS = 8.314462618
ATM = 101325.0


def graphite_vapour_pressure(temp: float) -> float:
    """Equilibrium vapour pressure [Pa] of monatomic carbon over graphite.

    Clausius-Clapeyron with a constant sublimation enthalpy, anchored at the 1 atm sublimation
    point. Constant `dH` is the approximation; over 3000-5000 K it costs tens of percent in `p`,
    which is a fraction of a percent in the threshold *temperature* the verdict is stated as --
    the exponential is steep enough that the answer is robust to it.
    """
    return ATM * math.exp(-GRAPHITE_SUBLIMATION / R_GAS * (1.0 / temp - 1.0 / GRAPHITE_1ATM_K))


def carbon_partial_pressure(rho: float, temp: float) -> float:
    """Partial pressure [Pa] of *monatomic* carbon in equilibrium methane at `(rho, temp)`.

    Monatomic, because that is the species graphite's vapour pressure is defined against. C2 and
    C3 carry most of the carbon at these temperatures and would deposit too, so this is the
    conservative (lowest) reading of the driving pressure.
    """
    return eos_methane.composition(rho, temp).n_c * eos_methane.K_B * temp


def deposition_threshold(rho: float, temp: float, lo: float = 1500.0, hi: float = 8000.0) -> float:
    """Wall temperature [K] above which a **frozen** boundary layer stops plating carbon.

    **This is one edge of a fork, not the answer.** It compares the *free stream's* monatomic
    carbon against graphite's vapour pressure at the *wall* temperature, which is the right
    comparison only if the gas crosses the boundary layer without its chemistry adjusting -- i.e.
    if the carbon atoms the throat is carrying arrive at the surface still as atoms.

    `equilibrium_saturation` is the other edge, and the two differ by six orders of magnitude in
    driving potential. `self_healing` picks between them and reports which applies.
    """
    driving = carbon_partial_pressure(rho, temp)
    while hi - lo > 0.5:
        mid = 0.5 * (lo + hi)
        if graphite_vapour_pressure(mid) < driving:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def wall_gas_density(pressure: float, wall_temp: float) -> float:
    """Density [kg/m^3] of the gas touching the wall: methane's elements at `(p, T_wall)`.

    Pressure is constant across a boundary layer and the gas at a no-slip surface is at the wall
    temperature, so the state adjacent to the surface is the free stream's *pressure* at the
    wall's *temperature*. That is a different state from either end, and using the wrong one is
    the mistake `deposition_threshold` alone would make.
    """
    return _density_at_pressure(surface.METHANE, pressure, wall_temp)


def equilibrium_saturation(pressure: float, wall_temp: float) -> float:
    """Carbon saturation ratio `S = p_C / p_vap(graphite)` for a boundary layer in equilibrium.

    `S > 1` means the gas is supersaturated and carbon plates onto the surface; `S < 1` means the
    gas can dissolve *more* carbon than it is carrying, so a graphite surface is chemically eroded
    on top of whatever the heat flux does to it.

    **The result is not intuitive and the reason is acetylene.** Between roughly 1,900 and 5,000 K
    a hydrogen-rich carbon gas parks its carbon in acetylene, not in free atoms -- the same
    stoichiometry W9 found from the recovery side -- and that sink is deep enough to hold the
    monatomic-carbon activity *below* graphite's own. It is the reason hydrogen-rich flames do not
    soot, and it is why a graphite throat inside this window cannot collect a deposit.
    """
    comp = eos_methane.composition(wall_gas_density(pressure, wall_temp), wall_temp)
    return comp.n_c * eos_methane.K_B * wall_temp / graphite_vapour_pressure(wall_temp)


def frozen_saturation(rho: float, temp: float, pressure: float, wall_temp: float) -> float:
    """The same ratio if the boundary layer's composition is frozen at the free stream's.

    Mole fractions are what survive a frozen compression or cooling, so the carbon partial
    pressure at the wall is `x_C` times the local pressure.
    """
    comp = eos_methane.composition(rho, temp)
    x_c = comp.n_c / comp.n_total
    return x_c * pressure / graphite_vapour_pressure(wall_temp)


def deposition_window(
    pressure: float, lo: float = 1200.0, hi: float = 5500.0, tol: float = 5.0
) -> tuple[float, float | None]:
    """`(cold edge, hot edge)` of the wall-temperature band where an equilibrium gas does not plate.

    The cold edge always exists: below it graphite's vapour pressure collapses and everything is
    supersaturated against it, which is just "solid carbon is the stable phase down there". The
    hot edge is `None` if `S` is still below 1 at `hi` -- the acetylene sink has not unwound yet
    within the range searched.
    """

    def cold_crossing() -> float:
        a, b = lo, hi
        while b - a > tol:
            mid = 0.5 * (a + b)
            if equilibrium_saturation(pressure, mid) > 1.0:
                a = mid
            else:
                b = mid
        return 0.5 * (a + b)

    cold = cold_crossing()
    if equilibrium_saturation(pressure, hi) <= 1.0:
        return cold, None
    a, b = cold, hi
    while b - a > tol:
        mid = 0.5 * (a + b)
        if equilibrium_saturation(pressure, mid) > 1.0:
            b = mid
        else:
            a = mid
    return cold, 0.5 * (a + b)


#: Graphite density [kg/m^3] and carbon atomic mass [kg], for turning a flux into a thickness.
GRAPHITE_DENSITY = 2200.0
M_CARBON = 12.011 * eos_methane.AMU


def hertz_knudsen_thickness(rho: float, temp: float, wall_temp: float, duration: float) -> float:
    """Ceiling on the graphite laid down in `duration` [s], as a thickness [m].

    Hertz-Knudsen net flux with unit sticking:

        `Gamma = (p_C - p_vap(T_wall)) / sqrt(2 pi m k T)`

    **This is a ceiling and a large one.** It assumes every carbon atom that would arrive by
    free flight sticks, with no boundary layer between the free stream and the surface to slow the
    delivery, and no re-entrainment. A real throat has a boundary layer, and diffusion across it is
    what actually limits deposition. The number's use is that if even this ceiling is small, the
    question is closed; if it is large, the answer is that the boundary layer has to be computed.
    """
    driving = carbon_partial_pressure(rho, temp) - graphite_vapour_pressure(wall_temp)
    if driving <= 0.0:
        return 0.0
    flux = driving / math.sqrt(2.0 * math.pi * M_CARBON * eos_methane.K_B * temp)
    return flux * M_CARBON * duration / GRAPHITE_DENSITY


# ---- Convective flux, the Bartz correlation (item 7) ---------------------------------------------

#: Prandtl number bracket for a dissociating plasma. Bartz's own recommendation for rocket
#: exhausts is 0.8; a partly ionised gas runs lower because electron conduction is fast.
PRANDTL_BRACKET = (0.5, 0.8)
#: Hard-sphere collision cross-section [m^2] for the neutral-atom viscosity estimate, bracketed.
#: Atomic hydrogen's gas-kinetic cross-section is a few 1e-20 m^2 at these temperatures.
CROSS_SECTION_BRACKET = (2.0e-20, 1.0e-19)
#: Wall temperature [K] the flux is quoted against: pyrolytic graphite's working ceiling.
WALL_TEMPERATURE = 3900.0
#: Throat radius-of-curvature to throat-diameter ratio Bartz's geometry factor wants.
CURVATURE_RATIO = 1.5
#: Choked-flow coefficient at `gamma = 1.2`, the same one `chamber.acoustic_check` uses.
CHOKED_COEFF = (2.0 / 2.2) ** (2.2 / 0.4)


def chapman_enskog_viscosity(mean_mass: float, temp: float, cross_section: float) -> float:
    """Chapman-Enskog first-approximation viscosity [Pa s] for a hard-sphere monatomic gas.

    `mu = (5/16) sqrt(pi m k T) / (pi d^2)`, with `pi d^2` supplied as `cross_section`. It is the
    right *form* for a dissociated gas of atoms and the wrong one for a plasma, where Coulomb
    collisions take over; at the 0.5% ionisation of a 10 kK chamber the neutral channel still
    dominates, which is why it is used, and the cross-section is bracketed rather than picked.
    """
    return (5.0 / 16.0) * math.sqrt(math.pi * mean_mass * eos_methane.K_B * temp) / cross_section


def frozen_heat_capacity(
    fluid: surface.Fluid, rho: float, temp: float, step: float = 1.02
) -> float:
    """Equilibrium `c_p` [J/kg/K] by finite difference of enthalpy at constant pressure.

    **Equilibrium, not frozen**, and the distinction is the whole size of the answer: a
    dissociating gas's equilibrium `c_p` is an order of magnitude above its frozen value because
    the chemistry absorbs heat. Bartz's `c_p` is a boundary-layer transport property, and which of
    the two belongs there depends on whether the chemistry keeps up across the boundary layer --
    which is the same Damkoehler question `surface.py` runs on the core flow, asked of a much
    thinner region. Both are carried in `bartz_flux`, and the answer is a factor of ~7 wide.
    """
    p0, _ = fluid.pressure_energy(rho, temp)

    def enthalpy_at(t: float) -> float:
        r = _density_at_pressure(fluid, p0, t)
        p, e = fluid.pressure_energy(r, t)
        return e + p / r

    return (enthalpy_at(temp * step) - enthalpy_at(temp / step)) / (temp * step - temp / step)


def _density_at_pressure(fluid: surface.Fluid, pressure: float, temp: float) -> float:
    """Density [kg/m^3] at which `fluid` has `pressure` at `temp`.

    Bisection in log density: `p` rises monotonically with `rho` at fixed `T` for every EOS here.
    """
    lo, hi = 1e-6, 1e4
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if fluid.pressure_energy(mid, temp)[0] < pressure:
            lo = mid
        else:
            hi = mid
        if hi / lo < 1.0 + 1e-10:
            break
    return math.sqrt(lo * hi)


@dataclass(frozen=True)
class BartzPoint:
    """Convective flux at one station, with the transport inputs that set it."""

    fluid: str
    volume: float
    temp_c: float
    throat_area: float
    #: Area ratio the flux is quoted at: 1 at the throat, `A_bore/A*` at the liner.
    area_ratio: float
    station: str
    pressure_c: float
    c_star: float
    viscosity: float
    cp_equilibrium: float
    cp_frozen: float
    prandtl: float
    #: Flux [W/m^2] with the equilibrium `c_p` -- the upper edge.
    flux_equilibrium: float
    #: Flux [W/m^2] with the frozen `c_p` -- the lower edge.
    flux_frozen: float
    #: Pulse fluence [MJ/m^2] at the upper edge, over the blowdown.
    fluence_equilibrium: float
    blowdown: float


def bartz_flux(
    fluid: surface.Fluid,
    volume: float,
    temp_c: float,
    throat_area: float,
    station: str = "throat",
    prandtl: float = 0.8,
    cross_section: float = 5.0e-20,
    wall_temp: float = WALL_TEMPERATURE,
) -> BartzPoint:
    """Bartz's 1957 correlation at the throat or at the bore, on solved chamber properties.

        h_g = (0.026 / D*^0.2) (mu^0.2 cp / Pr^0.6) (p_c / c*)^0.8 (D* / R_c)^0.1 sigma (A*/A)^0.9

    **Read the size, not the digits.** Bartz was fitted on chemical rockets at tens of bar and a
    few thousand kelvin. This chamber is at a thousand bar and ten thousand kelvin, the boundary
    layer is partly ionised, and the correlation's transport inputs (`mu`, `cp`, `Pr`) are the
    least certain things in it -- which is why all three are bracketed and both `c_p` conventions
    are carried. What the number is good for is the *scaling*: `q ~ p_c^0.8 (A*/A)^0.9`, so
    shrinking the throat from 7 m^2 to 0.1 m^2 raises the throat flux by the area ratio to the
    0.9, and that conclusion does not depend on the correlation's absolute calibration.
    """
    _k, rho_c, _u = propellants.solve_slug_ratio(fluid.pressure_energy, temp_c, volume)
    p_c = fluid.pressure_energy(rho_c, temp_c)[0]
    sound_c = fluid.sound_speed(rho_c, temp_c)
    mdot = CHOKED_COEFF * rho_c * sound_c * throat_area
    c_star = p_c * throat_area / mdot

    mean_mass = rho_c / _number_density(fluid, rho_c, temp_c)
    mu = chapman_enskog_viscosity(mean_mass, temp_c, cross_section)
    cp_eq = frozen_heat_capacity(fluid, rho_c, temp_c)
    cp_fr = 2.5 * eos_methane.K_B / mean_mass  # a monatomic mixture's frozen c_p

    d_star = 2.0 * math.sqrt(throat_area / math.pi)
    area_ratio = 1.0 if station == "throat" else surface.BORE_AREA / throat_area
    sigma_factor = 1.0 / (0.5 * (wall_temp / temp_c) + 0.5) ** 0.68

    def flux(cp: float) -> float:
        h = (
            (0.026 / d_star**0.2)
            * (mu**0.2 * cp / prandtl**0.6)
            * (p_c / c_star) ** 0.8
            * CURVATURE_RATIO**-0.1
            * sigma_factor
            * area_ratio**-0.9
        )
        return float(h * (temp_c - wall_temp))

    blowdown = surface.blowdown_time(rho_c, sound_c, volume, throat_area)
    return BartzPoint(
        fluid=fluid.name,
        volume=volume,
        temp_c=temp_c,
        throat_area=throat_area,
        area_ratio=area_ratio,
        station=station,
        pressure_c=p_c,
        c_star=c_star,
        viscosity=mu,
        cp_equilibrium=cp_eq,
        cp_frozen=cp_fr,
        prandtl=prandtl,
        flux_equilibrium=flux(cp_eq),
        flux_frozen=flux(cp_fr),
        fluence_equilibrium=flux(cp_eq) * blowdown / 1e6,
        blowdown=blowdown,
    )


def _number_density(fluid: surface.Fluid, rho: float, temp: float) -> float:
    """Total particle density [m^-3] from the ideal-mixture pressure, `n = p / kT`."""
    return fluid.pressure_energy(rho, temp)[0] / (eos_methane.K_B * temp)


def radiative_equilibrium_temperature(flux: float) -> float:
    """Surface temperature [K] a wall reaches when it re-radiates the flux it takes.

    `T = (q / sigma)^(1/4)`. It is the *floor* on the surface temperature under a convective load
    -- a real wall also conducts inward, which lowers it, and ablates, which caps it. It is used
    twice here: to test whether the Bartz throat load is compatible with a passive graphite
    surface at all, and to ask whether a throat under load is hot enough to stop carbon sticking.
    """
    return float((flux / expansion.SIGMA_SB) ** 0.25)


#: Heat of ablation of graphite [J/kg]: sublimation to monatomic vapour, `711.2 kJ/mol / 12 g`.
#: Sensible heating to the sublimation point is a further ~8 MJ/kg and is neglected, which makes
#: the ablation depth below an *over*estimate by under 15%.
GRAPHITE_ABLATION = GRAPHITE_SUBLIMATION / (12.011e-3)


def ablation_depth(fluence: float) -> float:
    """Graphite removed [m] by a fluence [J/m^2], with no transpiration credit.

    A blowing boundary layer typically halves the net flux at these blowing rates, so read this as
    the pessimistic edge of a factor-of-two band.
    """
    return fluence / (GRAPHITE_ABLATION * GRAPHITE_DENSITY)


def carbon_budget_fraction(mass: float, volume: float, temp_c: float = 10000.0) -> float:
    """Share of the pulse's carbon inventory a mass [kg] represents.

    The ask's item 5 asks whether the liner is self-healing and states its own threshold: "0.34 to
    3.7% redeposition suffices". This puts an ablated mass on that same scale.
    """
    k, _rho, _u = propellants.solve_slug_ratio(eos_methane.pressure_energy, temp_c, volume)
    slug = k * chamber.IMPACTOR_MASS
    carbon = slug * 12.011 / 16.043
    return mass / carbon


# ---- The strike W10 moves the problem to: the contraction ------------------------------------


@dataclass(frozen=True)
class EndWallStrike:
    """What the downstream end of a **sub-gate** chamber takes (W10's successor question).

    Below the gate volume the front leaves the column before its cone touches the liner -- which
    is the good news of W10. It does not leave the *hardware*. With a 0.05 to 1 m^2 throat against
    a 3 m bore, the convergent section is 96-99% of the bore area, so it is very nearly a flat end
    wall and the front stagnates on it at **normal incidence** rather than grazing past.

    Two loads land there and they scale differently, which is the point of this dataclass:

    - the **arrival**, concentrated on whatever patch the cone has opened to by the exit plane.
      It gets *worse* as the chamber shrinks, because a shorter column decelerates the front less;
    - the **bulk**, the rest of the charge turning the corner behind it. It is
      **volume-independent**, because the same pulse energy passes through the same bore area
      however long the column in front of it was.
    """

    fluid: str
    volume: float
    temp_c: float
    throat_area: float
    exit_radius: float
    exit_speed: float
    swept_mass: float
    #: Patch of the end plane the front has opened to when it gets there [m^2], capped at the bore.
    wetted_area: float
    #: Area of the contraction the whole charge must turn against [m^2].
    contraction_area: float
    #: Kinetic energy the front carries at the exit plane [J].
    front_energy: float
    #: Incident fluence [MJ/m^2] of the arrival, over the wetted patch.
    incident_arrival: float
    #: Incident fluence [MJ/m^2] of the whole pulse over the contraction.
    incident_bulk: float

    @property
    def front_share(self) -> float:
        """Share of the pulse energy the front is still carrying when it arrives."""
        return self.front_energy / chamber.PULSE_ENERGY

    def delivered_arrival(self, stanton: float) -> float:
        """Arrival fluence actually crossing into the surface [MJ/m^2] at this Stanton number."""
        return stanton * self.incident_arrival

    def delivered_bulk(self, stanton: float) -> float:
        """Bulk fluence actually crossing into the surface [MJ/m^2] at this Stanton number."""
        return stanton * self.incident_bulk


def end_wall_strike(
    fluid: surface.Fluid,
    volume: float,
    temp_c: float = 10000.0,
    throat_area: float = 0.5,
    spread_multiple: float = 1.0,
) -> EndWallStrike:
    """Price the contraction, for the chambers W10 clears of the side-wall strike.

    **The bulk term is a ceiling reached from the energy budget rather than a flux.** Everything
    the impactor brought has to turn the corner and leave through the throat, so the contraction
    sees all 70.3 GJ pass across it whatever the chamber's shape; dividing by the contraction area
    gives what would land if none of it turned. What actually crosses the surface is that times a
    Stanton number, the same bracket and the same caveat as `strike`.
    """
    run = run_front(fluid, volume, temp_c, spread_multiple)
    exit_state = run.exit_state
    wetted = math.pi * min(exit_state.radius, WALL_RADIUS) ** 2
    contraction = surface.BORE_AREA - throat_area
    mass = chamber.IMPACTOR_MASS + exit_state.swept_mass
    energy = 0.5 * mass * exit_state.speed**2
    return EndWallStrike(
        fluid=fluid.name,
        volume=volume,
        temp_c=temp_c,
        throat_area=throat_area,
        exit_radius=exit_state.radius,
        exit_speed=exit_state.speed,
        swept_mass=exit_state.swept_mass,
        wetted_area=wetted,
        contraction_area=contraction,
        front_energy=energy,
        incident_arrival=energy / wetted / 1e6,
        incident_bulk=chamber.PULSE_ENERGY / contraction / 1e6,
    )


# ---- What a narrower throat costs: the throat as a consumable (W24) --------------------------

#: Chamber temperatures the trade study is run at. Two, not one, because the throat and the
#: temperature dials are the study's only two live levers (W17, W24) and the useful question is
#: what they are worth *together* -- which is what W8's recommended-geometry rung quotes.
THROAT_TEMPERATURES = (10000.0, 12000.0)

#: Throat areas [m^2] the trade study is run over: the ask's flown 7 down to 0.1. It stops one
#: rung short of the N16 grid's 0.05 deliberately -- methane's blowdown already misses the 400 ms
#: pulse period at 0.1, and the effective Isp has been flat to within a percent since 0.5, so a
#: deeper rung would price a throat that is off the table on two counts at once.
THROAT_LADDER = (7.0, 4.0, 2.0, 1.0, 0.5, 0.2, 0.1)

#: Analytic exponent of throat life in throat area, `life ~ A*^1.6`. Derived, not fitted:
#: Bartz gives the throat flux as `D*^-0.2 ~ A*^-0.1`; the blowdown goes as `1/A*`; so the pulse
#: fluence -- and with it the recession -- goes as `A*^-1.1`. Reaching a *fractional* area growth
#: needs a recession proportional to `r* ~ A*^0.5`, so the pulse count goes as `A*^1.6`.
#: `test_throat_life_follows_the_analytic_area_exponent` pins it.
THROAT_LIFE_EXPONENT = 1.6


@dataclass(frozen=True)
class ThroatLife:
    """How many pulses a throat lasts, given W15's verdict that it does not self-heal.

    **This is the term that prices a narrow throat, and it is not the peak flux.** Bartz's throat
    flux goes as `D*^-0.2`, so halving the throat area raises it only about 7%. What doubles is
    the **blowdown time**, because the same chamber empties through a smaller hole. Fluence is
    flux times time, so the dwell carries the whole cost: *a narrow throat does not heat the
    throat harder, it heats it for twice as long.*

    W15 established that the throat plates nothing back at any survivable wall temperature -- it
    is chemically eroded on top of being thermally ablated -- so recession per pulse is a
    consumable rate rather than a survival question, and the right output is a pulse count.

    Recession is radial and the throat therefore **opens** as it erodes, which is why the metric
    is fractional area growth: the nozzle drifts back up the trade curve it was narrowed to
    escape.
    """

    fluid: str
    volume: float
    temp_c: float
    throat_area: float
    #: Throat radius [m], from `A* = pi r*^2`.
    throat_radius: float
    #: Throat flux [W/m^2] on the equilibrium `c_p` -- Bartz's upper edge.
    flux: float
    #: Chamber emptying time [s]: the dwell that carries the cost.
    blowdown: float
    #: Pulse fluence [MJ/m^2] at the throat.
    fluence: float
    #: Radial graphite recession per pulse [m], no transpiration credit (pessimistic edge).
    recession: float

    @property
    def area_growth_per_pulse(self) -> float:
        """`dA/A = 2 dr/r` per pulse: how fast the throat undoes its own narrowing."""
        return 2.0 * self.recession / self.throat_radius

    def pulses_to_area_growth(self, fraction: float) -> float:
        """Pulses until the throat area has grown by `fraction` (0.10 for +10%, 1.0 for 2x).

        Exact rather than linearised in `dA/A`, since a doubling is not a small perturbation:
        `A ~ r^2`, so the radius has to reach `r* sqrt(1 + fraction)`.
        """
        if self.recession <= 0.0:
            return math.inf
        return (math.sqrt(1.0 + fraction) - 1.0) * self.throat_radius / self.recession


def throat_life(
    fluid: surface.Fluid,
    volume: float = 200.0,
    temp_c: float = 10000.0,
    throat_area: float = 1.0,
) -> ThroatLife:
    """Price one throat area as a consumable, on the Bartz flux and the solved blowdown.

    **Read the scaling, not the digits.** Every caveat on `bartz_flux` applies -- the correlation
    is far outside its fit at a thousand bar and ten thousand kelvin -- and `ablation_depth` takes
    no transpiration credit, so a blowing boundary layer roughly doubles every pulse count here.
    What survives both is the exponent: `life ~ A*^1.6`, a factor of three per halving.
    """
    b = bartz_flux(fluid, volume, temp_c, throat_area, station="throat")
    return ThroatLife(
        fluid=fluid.name,
        volume=volume,
        temp_c=temp_c,
        throat_area=throat_area,
        throat_radius=math.sqrt(throat_area / math.pi),
        flux=b.flux_equilibrium,
        blowdown=b.blowdown,
        fluence=b.fluence_equilibrium,
        recession=ablation_depth(b.fluence_equilibrium * 1e6),
    )


# ---- Does the throat self-heal? (item 5, properly) ---------------------------------------

#: Gas-kinetic bimolecular rate coefficient [m^3/s], order of magnitude, for the *bound* below.
#: **This is not a rate and is never used as one** (`rates.py`'s provenance rule forbids inventing
#: coefficients). It is the collision frequency a bimolecular reaction cannot exceed, so
#: `decades_below_gas_kinetic` reports how far under it the carbon chemistry could run and still
#: finish inside the boundary layer -- a statement about what would have to be true to flip the
#: verdict, not a claim about any particular reaction.
GAS_KINETIC = 1.0e-16


@dataclass(frozen=True)
class BoundaryLayer:
    """The thin layer that decides which saturation limit applies."""

    #: Convective film coefficient [W/m^2/K] from Bartz.
    film_coefficient: float
    thermal_conductivity: float
    #: Thermal thickness [m], `k / h`.
    thickness: float
    #: Diffusion clock across it [s], `delta^2 / alpha`.
    time: float
    #: Particle density [m^-3] used for the collision-frequency bound.
    density: float
    #: Decades below gas-kinetic the carbon chemistry could run and still finish in `time`.
    decades_of_margin: float

    @property
    def equilibrium(self) -> bool:
        """Whether the boundary layer's chemistry keeps up. Read the margin, not the boolean."""
        return self.decades_of_margin > 0.0


def boundary_layer(
    point: BartzPoint, rho: float, temp: float, wall_temp: float = WALL_TEMPERATURE
) -> BoundaryLayer:
    """Thickness and clock of the throat's thermal boundary layer, off the Bartz solve.

    **The thickness is free of Bartz's largest uncertainty.** `delta = k/h` with `k = mu c_p / Pr`,
    and Bartz's `h` is itself proportional to `c_p` -- so `c_p` cancels and the frozen and
    equilibrium branches, which differ by a factor of two in flux, give *identical* thicknesses.
    What is left depends on `mu` and `Pr`, which are bracketed.

    The clock is the diffusion time across it, `delta^2 / alpha` with `alpha = mu / (rho Pr)`.
    """
    h = point.flux_equilibrium / (temp - wall_temp)
    k_th = point.viscosity * point.cp_equilibrium / point.prandtl
    thickness = k_th / h
    alpha = point.viscosity / (rho * point.prandtl)
    time = thickness**2 / alpha
    n = eos_methane.composition(rho, 0.5 * (temp + wall_temp)).n_total
    return BoundaryLayer(
        film_coefficient=h,
        thermal_conductivity=k_th,
        thickness=thickness,
        time=time,
        density=n,
        decades_of_margin=math.log10(time * GAS_KINETIC * n),
    )


@dataclass(frozen=True)
class SelfHealing:
    """Item 5's actual question: does the throat's own carbon replace what the heat takes off?"""

    volume: float
    temp_c: float
    throat_area: float
    throat_temp: float
    throat_rho: float
    throat_pressure: float
    wall_temp: float
    #: Saturation ratio if the boundary layer's chemistry keeps up.
    saturation_equilibrium: float
    #: Saturation ratio if it is frozen at the free stream's composition.
    saturation_frozen: float
    layer: BoundaryLayer
    #: Wall temperature [K] below which even an equilibrium gas plates.
    cold_edge: float
    #: Graphite removed per pulse [m] by the convective load, no transpiration credit.
    ablation: float

    @property
    def plates(self) -> bool:
        """Whether carbon deposits at all, on the limit the boundary layer actually sits in."""
        ratio = self.saturation_equilibrium if self.layer.equilibrium else self.saturation_frozen
        return ratio > 1.0

    @property
    def self_heals(self) -> bool:
        """A throat self-heals only if it plates at all. It does not."""
        return self.plates


def self_healing(
    volume: float = 100.0,
    temp_c: float = 8000.0,
    throat_area: float = 0.5,
    wall_temp: float = WALL_TEMPERATURE,
) -> SelfHealing:
    """Settle item 5: plate, erode, or balance?

    Three states have to be kept apart and conflating any two of them gets the wrong answer:

    - **the free stream** at the throat, hot and carrying free carbon atoms;
    - **the wall-adjacent gas**, at the free stream's *pressure* and the wall's *temperature*;
    - **the surface**, graphite at the wall temperature.

    Deposition is driven by the second against the third, not the first against the third. The
    first-against-third comparison is `deposition_threshold` and it is the frozen edge.
    """
    _k, rho_c, _u = propellants.solve_slug_ratio(eos_methane.pressure_energy, temp_c, volume)
    rows = expansion.cooling_history(
        rho_c,
        temp_c,
        eos_methane.pressure_energy,
        eos_methane.sound_speed,
        surface.BORE_AREA / throat_area,
        1.0,
        steps=140,
        expansion_ratio=4.0e5,
    )
    throat = rows[0]
    point = bartz_flux(surface.METHANE, volume, temp_c, throat_area, "throat", wall_temp=wall_temp)
    layer = boundary_layer(point, throat.rho, throat.temp, wall_temp)
    cold, _hot = deposition_window(throat.pressure)
    return SelfHealing(
        volume=volume,
        temp_c=temp_c,
        throat_area=throat_area,
        throat_temp=throat.temp,
        throat_rho=throat.rho,
        throat_pressure=throat.pressure,
        wall_temp=wall_temp,
        saturation_equilibrium=equilibrium_saturation(throat.pressure, wall_temp),
        saturation_frozen=frozen_saturation(throat.rho, throat.temp, throat.pressure, wall_temp),
        layer=layer,
        cold_edge=cold,
        ablation=ablation_depth(point.fluence_equilibrium * 1e6),
    )


# ---- Output --------------------------------------------------------------------------------------

FRONT_HEADER = (
    "fluid,volume_m3,temp_c_k,spread_multiple,ambient_rho,length_m,cone_half_angle_deg,"
    "cone_contact_x_m,reaches_wall,contact_x_m,contact_time_us,contact_speed_m_s,swept_kg,"
    "local_slug_ratio,shocked_temp_k,shocked_rho,exit_radius_m,exit_speed_m_s\n"
)


def write_fronts(runs: Sequence[FrontRun], path: Path = DEFAULT_OUTPUT) -> None:
    """Write the front ledger. Committed: the asks document is read in the paper repository."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(FRONT_HEADER)
        for r in runs:
            c = r.contact
            fh.write(
                f"{r.fluid},{r.volume:g},{r.temp_c:g},{r.spread_multiple:g},"
                f"{r.ambient_rho:.5f},{r.length:.4f},{r.cone_half_angle:.2f},"
                f"{r.cone_contact_x:.4f},{r.reaches_wall},"
                f"{'' if c is None else f'{c.x:.4f}'},"
                f"{'' if c is None else f'{c.time * 1e6:.2f}'},"
                f"{'' if c is None else f'{c.speed:.1f}'},"
                f"{'' if c is None else f'{c.swept_mass:.2f}'},"
                f"{'' if c is None else f'{c.slug_ratio:.3f}'},"
                f"{'' if c is None else f'{c.shocked_temp:.0f}'},"
                f"{'' if c is None else f'{c.shocked_rho:.4f}'},"
                f"{r.exit_state.radius:.4f},{r.exit_state.speed:.1f}\n"
            )


STRIKE_HEADER = (
    "fluid,volume_m3,temp_c_k,spread_multiple,contact_x_m,transient_us,shocked_temp_k,"
    "shocked_rho,radial_speed_m_s,optical_depth,regime,cooling_time_s,blackbody_flux_w_m2,"
    "fluence_radiative_mj_m2,incident_convective_mj_m2,fluence_conv_low,fluence_conv_high,"
    "fluence_ceiling_mj_m2\n"
)


def write_strikes(strikes: Sequence[Strike], path: Path = STRIKE_OUTPUT) -> None:
    """Write the contact-transient ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(STRIKE_HEADER)
        for s in strikes:
            fh.write(
                f"{s.fluid},{s.volume:g},{s.temp_c:g},{s.spread_multiple:g},{s.contact_x:.4f},"
                f"{s.transient * 1e6:.2f},{s.shocked_temp:.0f},{s.shocked_rho:.4f},"
                f"{s.radial_speed:.1f},{s.optical_depth:.6e},{s.regime},{s.cooling_time:.6e},"
                f"{s.blackbody_flux:.6e},{s.fluence_radiative:.6f},{s.incident_convective:.3f},"
                f"{s.fluence_convective_low:.4f},{s.fluence_convective_high:.4f},"
                f"{s.fluence_ceiling:.2f}\n"
            )


BARTZ_HEADER = (
    "fluid,volume_m3,temp_c_k,throat_area_m2,station,area_ratio,pressure_c_bar,c_star,"
    "viscosity,cp_equilibrium,cp_frozen,prandtl,flux_equilibrium_mw_m2,flux_frozen_mw_m2,"
    "blowdown_ms,fluence_equilibrium_mj_m2,surface_temp_k,ablation_um\n"
)


def write_bartz(points: Sequence[BartzPoint], path: Path = BARTZ_OUTPUT) -> None:
    """Write the convective-flux ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(BARTZ_HEADER)
        for b in points:
            fh.write(
                f"{b.fluid},{b.volume:g},{b.temp_c:g},{b.throat_area:g},{b.station},"
                f"{b.area_ratio:.3f},{b.pressure_c / 1e5:.1f},{b.c_star:.1f},{b.viscosity:.4e},"
                f"{b.cp_equilibrium:.4e},{b.cp_frozen:.4e},{b.prandtl:g},"
                f"{b.flux_equilibrium / 1e6:.2f},{b.flux_frozen / 1e6:.2f},"
                f"{b.blowdown * 1e3:.2f},{b.fluence_equilibrium:.2f},"
                f"{radiative_equilibrium_temperature(b.flux_equilibrium):.0f},"
                f"{ablation_depth(b.fluence_equilibrium * 1e6) * 1e6:.1f}\n"
            )


THROAT_LIFE_HEADER = (
    "fluid,volume_m3,temp_c_k,throat_area_m2,throat_radius_mm,area_ratio,slug_ratio,exit_temp_k,"
    "conversion,conversion_capped,verdict,isp_true_s,isp_effective_s,isp_ceiling_s,"
    "isp_gain_vs_7m2,"
    "flux_throat_mw_m2,blowdown_ms,blowdown_fits,fluence_mj_m2,recession_mm_per_pulse,"
    "area_growth_per_pulse,pulses_to_plus_10pc,pulses_to_2x,below_carbon_floor\n"
)


def write_throat_life(
    fluids: Sequence[surface.Fluid] = (surface.METHANE, surface.HYDROGEN),
    volume: float = 200.0,
    temps: Sequence[float] = THROAT_TEMPERATURES,
    path: Path = THROAT_LIFE_OUTPUT,
) -> None:
    """Write W24's trade table: what a narrower throat buys against what it costs.

    `isp_ceiling_s` is the **zero-loss** effective Isp at that chamber's own `k`: the same ledger
    evaluated at `eta_jet = 1`, i.e. full conversion through a perfect nozzle. It is carried so a
    reader can see how much of the available impulse a given geometry actually collects, and so
    that a proposed number can be checked against what the architecture allows before it is
    argued about. `slug_ratio` is beside it because the ceiling is a function of `k` alone.

    Joined here rather than in `surface.py` because the cost side is this module's -- the Bartz
    flux, the blowdown dwell and the ablation rate -- while the gain side is one column read off
    the N16 grid. Committed for the cross-repo reason ADR-0025 records.

    **The gain column comes from the grid and not from `propellants.solved_rung`, deliberately.**
    That function holds the nozzle length at the ask's 7.1 m and applies no freeze cap, both of
    which are fine at the ask's own wide throat and wrong at a deep one. `surface.column` rebuilds
    the cone geometry per throat and carries `conversion_capped` with a Damkoehler verdict, so it
    is the only source in the study that can be trusted down the ladder.

    `isp_gain_vs_7m2` is taken against the widest throat **at the same temperature**, so the
    column is a throat effect with the temperature dial held still.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(THROAT_LIFE_HEADER)
        for fluid in fluids:
            for temp_c in temps:
                points = {
                    p.throat_area: p
                    for p in surface.column(fluid, temp_c, volume, throats=THROAT_LADDER)
                }
                base = points[max(THROAT_LADDER)].isp_effective
                for area in THROAT_LADDER:
                    g = points[area]
                    life = throat_life(fluid, volume, temp_c, area)
                    # Zero-loss effective Isp at this chamber's own k: the ledger at eta_jet = 1.
                    ceiling = chamber.isp_ledger(g.exhaust_speed_ideal, g.slug_ratio)
                    fh.write(
                        f"{fluid.name},{volume:g},{temp_c:g},{area:g},"
                        f"{life.throat_radius * 1e3:.0f},{g.area_ratio:.1f},"
                        f"{g.slug_ratio:.4f},{g.exit_temp:.0f},"
                        f"{g.conversion:.4f},{g.conversion_capped:.4f},{g.verdict},"
                        f"{g.isp_true:.1f},{g.isp_effective:.1f},"
                        f"{ceiling.isp_effective:.1f},"
                        f"{g.isp_effective / base - 1.0:.4f},"
                        f"{life.flux / 1e6:.0f},{life.blowdown * 1e3:.1f},{g.blowdown_fits},"
                        f"{life.fluence:.1f},{life.recession * 1e3:.3f},"
                        f"{life.area_growth_per_pulse:.5f},"
                        f"{life.pulses_to_area_growth(0.10):.0f},"
                        f"{life.pulses_to_area_growth(1.00):.0f},{g.below_carbon_floor}\n"
                    )


def _report_throat_trade(volume: float = 200.0, temp_c: float = 10000.0) -> None:
    """Print W24: the throat-reduction trade, measured from the ask's own 7 m^2 baseline."""
    print(
        f"\n=== W24: what reducing the throat buys and costs, from the ask's 7 m^2 "
        f"({volume:.0f} m^3, {temp_c:.0f} K) ==="
    )
    for fluid in (surface.METHANE, surface.HYDROGEN):
        points = {
            p.throat_area: p for p in surface.column(fluid, temp_c, volume, throats=THROAT_LADDER)
        }
        base = points[max(THROAT_LADDER)].isp_effective
        print(f"\n  {fluid.name}: baseline A* = 7 m^2 gives {base:.0f} s effective")
        print(
            f"  {'A*':>5} {'A/A*':>6} {'exit T':>7} {'conv/capped':>12} {'Isp eff':>8} "
            f"{'vs 7':>7} {'/halving':>9} {'blowdown':>9} {'recess':>8} {'to +10%':>8} "
            f"{'flags':>17}"
        )
        previous = None
        for area in THROAT_LADDER:
            g = points[area]
            life = throat_life(fluid, volume, temp_c, area)
            step = 0.0 if previous is None else 100.0 * (g.isp_effective / previous - 1.0)
            marginal = "--" if previous is None else f"{step:+.1f}%"
            flags = []
            if g.verdict != "equilibrium":
                flags.append("FREEZE")
            if g.below_carbon_floor:
                flags.append("C-floor")
            if not g.blowdown_fits:
                flags.append("BLOWDOWN")
            print(
                f"  {area:5.2f} {g.area_ratio:6.1f} {g.exit_temp:7.0f} "
                f"{g.conversion:5.3f}/{g.conversion_capped:6.3f} {g.isp_effective:8.0f} "
                f"{100.0 * (g.isp_effective / base - 1.0):+6.1f}% {marginal:>9} "
                f"{life.blowdown * 1e3:8.1f}ms {life.recession * 1e3:7.2f}mm "
                f"{life.pulses_to_area_growth(0.10):8.0f} {','.join(flags):>17}"
            )
            previous = g.isp_effective
    print(
        f"\n  Throat life goes as A*^{THROAT_LIFE_EXPONENT} -- a factor of three per halving --"
        "\n  because the throat FLUX barely moves (Bartz D*^-0.2, ~7% per halving) while the"
        "\n  blowdown DWELL doubles. Recession is the pessimistic edge: no transpiration credit."
    )


def _report_recommended_rung(volume: float = 200.0) -> None:
    """Print the W8 rung: the fluids re-run at the geometry W24 and W17 actually recommend.

    W8's ladder is pinned to the ask's own 7 m^2 throat and 10 kK, which is the right choice for
    answering the ask on its own terms and the wrong one for stating what the architecture can do.
    This block moves both live dials to where this study recommends them and carries the throat
    life alongside, so the specific impulse cannot be quoted without its consumable cost.

    **Every column is like-for-like**: each row is one fluid at one geometry, and methane sits
    beside hydrogen at the same throat and the same chamber temperature.
    """
    print("\n=== The W8 rung at the recommended geometry (W17's temperature, W24's throat) ===")
    print(
        f"  {'fluid':>9} {'T_c':>7} {'A*':>5} {'exit T':>7} {'conv':>6} {'verdict':>11} "
        f"{'Isp true':>9} {'Isp EFF':>8} {'vs ask':>7} {'to +10%':>8}"
    )
    for fluid in (surface.METHANE, surface.HYDROGEN):
        reference = None
        for temp_c in THROAT_TEMPERATURES:
            points = {
                p.throat_area: p
                for p in surface.column(fluid, temp_c, volume, throats=(7.0, 2.0, 1.0))
            }
            for area in (7.0, 2.0, 1.0):
                g = points[area]
                if reference is None:
                    reference = g.isp_effective  # the ask's own point: 7 m^2 at 10 kK
                life = throat_life(fluid, volume, temp_c, area)
                print(
                    f"  {fluid.name:>9} {temp_c:7.0f} {area:5.2f} {g.exit_temp:7.0f} "
                    f"{g.conversion_capped:6.3f} {g.verdict:>11} {g.isp_true:9.0f} "
                    f"{g.isp_effective:8.0f} "
                    f"{100.0 * (g.isp_effective / reference - 1.0):+6.1f}% "
                    f"{life.pulses_to_area_growth(0.10):8.0f}"
                )
    print(
        "\n  'vs ask' is against that fluid's own 7 m^2 / 10 kK point, so it is the gain the two"
        "\n  dials buy together. 'to +10%' is pulses until the throat has eroded 10% wider: the"
        "\n  Isp and the throat life are one choice, not two."
    )


def main() -> None:
    """Answer N9 items 1-7 in the order the ask asks for them: the gate first."""
    print("N9: does the shocked front reach the wall, and what does it deliver?\n")
    print(
        f"{chamber.IMPACTOR_MASS:.0f} kg at {chamber.CLOSING_SPEED / 1e3:.0f} km/s into methane, "
        f"{WALL_RADIUS:.0f} m bore, chamber charged at 10 kK\n"
    )

    print("=== The gate: below what volume does the front leave before it touches? ===")
    print(f"{'spread':>8} {'T_c [K]':>8} {'gate V [m^3]':>13} {'gate L [m]':>11}")
    gates: dict[float, float] = {}
    for mult in SPREAD_BRACKET:
        g = gate_volume(surface.METHANE, 10000.0, mult)
        gates[mult] = g
        print(f"{mult:8.1f} {10000:8.0f} {g:13.0f} {chamber_length(g):11.2f}")
    print(f"\npaper-side estimate: {PAPER_GATE_VOLUME:.0f} m^3")

    print("\n=== Items 1 and 2: contact station, time, and the state there ===")
    runs: list[FrontRun] = []
    strikes: list[Strike] = []
    print(
        f"{'V [m3]':>7} {'L [m]':>7} {'rho_0':>7} {'cone':>6} {'x_c [m]':>9} {'t [us]':>8} "
        f"{'v_c [km/s]':>11} {'swept [kg]':>11} {'k_local':>8} {'T_shock [K]':>12}"
    )
    for volume in (*VOLUMES, 673.0):
        run = run_front(surface.METHANE, volume, 10000.0)
        runs.append(run)
        c = run.contact
        if c is None:
            print(
                f"{volume:7.0f} {run.length:7.2f} {run.ambient_rho:7.3f} "
                f"{run.cone_half_angle:6.1f} {'NEVER':>9}   front leaves at r = "
                f"{run.exit_state.radius:.2f} m of {WALL_RADIUS:.0f}"
            )
            continue
        print(
            f"{volume:7.0f} {run.length:7.2f} {run.ambient_rho:7.3f} {run.cone_half_angle:6.1f} "
            f"{c.x:9.2f} {c.time * 1e6:8.1f} {c.speed / 1e3:11.1f} {c.swept_mass:11.1f} "
            f"{c.slug_ratio:8.2f} {c.shocked_temp:12.0f}"
        )
        st = strike(run)
        if st is not None:
            strikes.append(st)

    print("\n=== Item 3: what the transient delivers, MJ/m^2 ===")
    print(
        f"{'V [m3]':>7} {'tau':>10} {'regime':>15} {'sigT^4 [W/m2]':>14} {'radiative':>10} "
        f"{'convective (St bracket)':>25} {'ceiling':>9}"
    )
    for hit in strikes:
        print(
            f"{hit.volume:7.0f} {hit.optical_depth:10.3g} {hit.regime:>15} "
            f"{hit.blackbody_flux:14.3e} {hit.fluence_radiative:10.4f} "
            f"{hit.fluence_convective_low:11.3f} - {hit.fluence_convective_high:9.3f} "
            f"{hit.fluence_ceiling:9.1f}"
        )
    print(f"\nthe equilibrium load the section already clears: {EQUILIBRIUM_FLUENCE:.1f} MJ/m^2")

    print("\n  opacity bracket at 673 m^3 (the verdict has to survive all of it):")
    deep = next(r for r in runs if r.volume == 673.0)
    for scale in OPACITY_BRACKET:
        probed = strike(deep, scale)
        if probed is None:  # pragma: no cover - `deep` reaches the wall by construction
            continue
        print(
            f"    kappa x{scale:<5} tau = {probed.optical_depth:10.4g}  {probed.regime:>15}  "
            f"radiative = {probed.fluence_radiative:.4f} MJ/m^2"
        )

    print("\n=== W10's successor: what the contraction takes when the front misses the liner ===")
    print(
        f"{'V [m3]':>7} {'exit r':>7} {'v [km/s]':>9} {'% of pulse':>11} "
        f"{'arrival (St bracket)':>22} {'bulk (St bracket)':>20}"
    )
    for volume in (50.0, 100.0, 150.0, 178.0):
        e = end_wall_strike(surface.METHANE, volume)
        print(
            f"{volume:7.0f} {e.exit_radius:7.2f} {e.exit_speed / 1e3:9.1f} "
            f"{e.front_share * 100:10.1f}% "
            f"{e.delivered_arrival(STANTON_BRACKET[0]):10.2f} -"
            f"{e.delivered_arrival(STANTON_BRACKET[1]):9.2f} "
            f"{e.delivered_bulk(STANTON_BRACKET[0]):9.2f} -"
            f"{e.delivered_bulk(STANTON_BRACKET[1]):8.2f}"
        )
    print(
        "\n  All MJ/m^2 against the same 1.6 reference. The bulk column is identical at every"
        "\n  volume -- the same pulse turns the same corner -- while the arrival column gets"
        "\n  WORSE as the chamber shrinks, because a shorter column decelerates the front less."
    )

    print("\n=== Item 4: does a cold film measurably reduce it? ===")
    print(f"{'film [kg/m2]':>13} {'capacity [MJ/m2]':>17}   against the loads above")
    for areal in FILM_AREAL:
        print(f"{areal:13.2f} {film_capacity(areal, surface.METHANE) / 1e6:17.1f}")

    print("\n=== Item 7: convective flux, Bartz, at the ask's own grid ===")
    bartz: list[BartzPoint] = []
    print(
        f"{'V [m3]':>7} {'T_c [K]':>8} {'A* [m2]':>8} {'station':>8} {'p_c [bar]':>10} "
        f"{'q frozen':>10} {'q equil':>10} {'MW/m2':>6} {'T_surf [K]':>11}"
    )
    for volume in LEGACY_VOLUMES:
        for temp in (10000.0, 12000.0, 15000.0):
            for station in ("bore", "throat"):
                b = bartz_flux(surface.METHANE, volume, temp, 7.0, station)
                bartz.append(b)
                print(
                    f"{volume:7.0f} {temp:8.0f} {7.0:8.2f} {station:>8} {b.pressure_c / 1e5:10.0f} "
                    f"{b.flux_frozen / 1e6:10.0f} {b.flux_equilibrium / 1e6:10.0f} {'':>6} "
                    f"{radiative_equilibrium_temperature(b.flux_equilibrium):11.0f}"
                )

    print("\n  and at the throats the ask now recommends, 100 m^3 at 8 kK:")
    print(
        f"{'A* [m2]':>8} {'q frozen':>10} {'q equil':>10} {'blowdown':>10} {'fluence':>10} "
        f"{'ablation':>10} {'% of C':>8}"
    )
    for area in surface.THROATS:
        b = bartz_flux(surface.METHANE, 100.0, 8000.0, area, "throat")
        bartz.append(b)
        removed = ablation_depth(b.fluence_equilibrium * 1e6)
        mass = removed * GRAPHITE_DENSITY * area
        print(
            f"{area:8.2f} {b.flux_frozen / 1e6:10.0f} {b.flux_equilibrium / 1e6:10.0f} "
            f"{b.blowdown * 1e3:9.0f}ms {b.fluence_equilibrium:10.0f} {removed * 1e3:9.2f}mm "
            f"{carbon_budget_fraction(mass, 100.0) * 100:8.2f}"
        )

    print("\n=== Item 5: does the throat self-heal? ===")
    print(
        "  Three states, and conflating any two gets the wrong answer: the free stream, the\n"
        "  wall-adjacent gas (free-stream pressure, wall temperature), and the graphite surface.\n"
    )
    print(
        f"{'A* [m2]':>8} {'throat K':>9} {'p [bar]':>9} {'S equil':>9} {'S frozen':>9} "
        f"{'layer um':>9} {'clock us':>9} {'margin':>8} {'plates?':>8} {'ablation':>10}"
    )
    for area in (2.0, 0.5, 0.2, 0.1):
        h = self_healing(100.0, 8000.0, area)
        print(
            f"{area:8.2f} {h.throat_temp:9.0f} {h.throat_pressure / 1e5:9.0f} "
            f"{h.saturation_equilibrium:9.4g} {h.saturation_frozen:9.4g} "
            f"{h.layer.thickness * 1e6:9.1f} {h.layer.time * 1e6:9.2f} "
            f"{h.layer.decades_of_margin:8.1f} {h.plates!s:>8} {h.ablation * 1e3:9.2f}mm"
        )
    probe = self_healing(100.0, 8000.0, 0.5)
    print(
        f"\n  The boundary layer is {probe.layer.decades_of_margin:.1f} decades from freezing, so"
        f" the equilibrium column is\n  the one that applies. S = "
        f"{probe.saturation_equilibrium:.3f} < 1: the gas can hold "
        f"{1.0 / probe.saturation_equilibrium:.0f}x more carbon than it\n  carries, so it erodes"
        f" graphite rather than plating it. **The throat does not self-heal.**"
    )
    print(f"\n  Saturation against wall temperature at {probe.throat_pressure / 1e5:.0f} bar:")
    print(f"{'T_wall':>9} {'S':>10}   (S > 1 plates)")
    for temp in (1500.0, 1800.0, 1912.0, 2000.0, 2500.0, 3000.0, 3900.0, 4500.0):
        print(f"{temp:9.0f} {equilibrium_saturation(probe.throat_pressure, temp):10.4g}")
    print(
        f"  cold edge of the no-plating band: {probe.cold_edge:.0f} K -- only a wall colder than"
        f" that collects carbon."
    )

    print("\n  The frozen edge, kept because it is what the earlier reading compared:")
    print(
        f"{'T_c [K]':>8} {'throat T':>9} {'rho':>8} {'p_C [Pa]':>11} {'self-clean above':>17} "
        f"{'graphite ceiling':>17}"
    )
    for temp in (8000.0, 10000.0):
        _k, rho_c, _u = propellants.solve_slug_ratio(eos_methane.pressure_energy, temp, 100.0)
        rows = expansion.cooling_history(
            rho_c,
            temp,
            eos_methane.pressure_energy,
            eos_methane.sound_speed,
            4.04,
            1.0,
            steps=120,
        )
        th = rows[0]
        print(
            f"{temp:8.0f} {th.temp:9.0f} {th.rho:8.3f} "
            f"{carbon_partial_pressure(th.rho, th.temp):11.3e} "
            f"{deposition_threshold(th.rho, th.temp):17.0f} {WALL_TEMPERATURE:17.0f}"
        )

    _report_throat_trade()
    _report_recommended_rung()

    write_fronts(runs)
    write_strikes(strikes)
    write_bartz(bartz)
    write_throat_life()
    print(f"\nwrote {DEFAULT_OUTPUT}, {STRIKE_OUTPUT}, {BARTZ_OUTPUT}, {THROAT_LIFE_OUTPUT}")


if __name__ == "__main__":
    main()
