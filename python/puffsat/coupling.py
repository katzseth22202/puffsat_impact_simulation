"""Study 2: does a projectile couple to a droplet cloud the way it couples to a vapour?

The routing document asks this because the field leak turned out small (0.11-2.5%, Q-M), so the
bag need not be pressurised and hot: 213 kg spread over 660 m^3 could as well be a cloud of
droplets. **The bag does not become unnecessary either way** -- spreading the slug over that
volume is a magnetic-field requirement independent of heat -- but the *phase* of what is in it
becomes a free choice, and `k = 8.5` had better not depend on it.

**Two ways a droplet cloud could differ from a vapour, and neither binds.**

1. **Inertia.** A droplet has to be brought up to the front's speed, and it is not a fluid
   element -- it lags. The relaxation time is `tau ~ (8/3) a rho_water / (C_d rho_gas dv)`, linear
   in droplet radius and inverse in relative speed. At 45 km/s a 10 um drop couples in ~1.7 us
   against a ~2.3 ms transit, and the limiting size is **~1.4 cm**. Margin depends on how the
   cloud is made: condensation from vapour gives sub-micron drops (4 orders), an injected spray
   gives tens to hundreds of microns (**2 orders**). Two is the honest figure, since a spray is
   what a designer would actually build.
2. **Discreteness.** The front could pass *between* the drops. It cannot: even centimetre drops
   are intercepted in the millions across a 3 m front over 23.8 m of bore.

Flash vaporisation is the other route to coupling and it is ~1000x faster than drag, so drag is
the criterion -- the slower of the two is what has to fit inside the transit.

**What Study 2 turned up that it did not ask about.** The paper's snowplow (item 11) is
`m(x) = m_0 + rho A x` with `A` the **full bore area**, i.e. it assumes the sweeping front spans
3 m from the moment the projectile enters. A compact 25 kg ice projectile is 0.187 m. Integrating
with a front that has to *grow* -- `dr/dx = c_exp / v`, which is the honest version -- falls short
of the paper's `k` by more than a factor two at any plausible lateral speed, and falls short
hardest at high closing speed. **This is not a droplet finding; it applies to the vapour case
equally, and it says `k = 8.5` is an assumption about the projectile's arrival radius rather than
an output of the snowplow.** That radius is `aim`'s to state. See Q-Q.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: Bag bulk density [kg/m^3]: 213 kg over the ~660 m^3 bore.
BAG_RHO = 0.323

#: Bore radius [m] and length [m] -- the bag is the bore.
BORE_RADIUS = 3.0
BORE_LENGTH = 23.8

#: Projectile mass [kg] on the paper's 5 m plate: 213 kg of slug at `k = 8.5`.
PROJECTILE_MASS = 25.0

#: Solid/liquid water density [kg/m^3]. Ice at 917; the answer is insensitive at this margin.
RHO_WATER = 917.0

#: Latent heat to vaporise water [J/kg]. Sublimation, not the dissociation energy: the droplet
#: only has to *leave the condensed phase* to stop being a lagging inertial body.
H_VAPORISATION = 2.8e6

#: Newtonian drag coefficient for a sphere in hypersonic free-molecular-to-continuum flow. It is
#: O(1) across every regime here, and the verdict has three orders of margin, so it is not swept.
DRAG_COEFFICIENT = 1.0


def drag_time(
    radius: float,
    gas_density: float,
    relative_speed: float,
    droplet_density: float = RHO_WATER,
    drag_coefficient: float = DRAG_COEFFICIENT,
) -> float:
    """Particle relaxation time [s]: how long the flow takes to bring a droplet up to speed.

        tau = m_d dv / F = (8/3) a rho_droplet / (C_d rho_gas dv)

    An e-folding time rather than a completion time -- as the droplet accelerates `dv` falls and
    the drag with it -- which is the conservative reading for a criterion.
    """
    if min(radius, gas_density, relative_speed, droplet_density, drag_coefficient) <= 0.0:
        raise ValueError("all arguments must be positive")
    return (
        (8.0 / 3.0) * radius * droplet_density / (drag_coefficient * gas_density * relative_speed)
    )


def vaporisation_time(
    radius: float,
    gas_density: float,
    relative_speed: float,
    droplet_density: float = RHO_WATER,
    latent_heat: float = H_VAPORISATION,
) -> float:
    """Time [s] to flash the droplet to vapour under the stagnation enthalpy flux `(1/2) rho v^3`.

        tau = m_d h_vap / (q A) = (8/3) a rho_droplet h_vap / (rho_gas dv^3)

    The cubic in relative speed is why this is never the binding channel here: it beats drag by
    ~3 orders at 45 km/s. Carried so that `coupling_time` can take the slower of the two rather
    than assume which one it is.
    """
    if min(radius, gas_density, relative_speed, droplet_density, latent_heat) <= 0.0:
        raise ValueError("all arguments must be positive")
    return (8.0 / 3.0) * radius * droplet_density * latent_heat / (gas_density * relative_speed**3)


def coupling_time(
    radius: float, gas_density: float, relative_speed: float, droplet_density: float = RHO_WATER
) -> float:
    """The slower of the two routes to coupling -- what actually has to fit inside the transit."""
    return max(
        drag_time(radius, gas_density, relative_speed, droplet_density),
        vaporisation_time(radius, gas_density, relative_speed, droplet_density),
    )


def coupling_number(
    radius: float,
    transit_time: float,
    gas_density: float,
    relative_speed: float,
    droplet_density: float = RHO_WATER,
) -> float:
    """`transit / coupling` -- the Damkohler number of this question. Above 1, the droplet joins."""
    return transit_time / coupling_time(radius, gas_density, relative_speed, droplet_density)


def limiting_radius(
    transit_time: float,
    gas_density: float,
    relative_speed: float,
    droplet_density: float = RHO_WATER,
    drag_coefficient: float = DRAG_COEFFICIENT,
) -> float:
    """The droplet radius [m] whose coupling time equals the transit -- Study 2's answer.

    Closed form because drag binds: invert `tau_drag = transit` for `a`. Anything smaller joins
    the slug and `k` is blind to phase; anything larger passes through as an unswept projectile
    of its own.

    Linear in `transit_time`, `gas_density` and `relative_speed`, so the margin is worst on the
    **cold, thin** corner -- which is the one this repository quotes it at.
    """
    return (
        3.0
        * transit_time
        * drag_coefficient
        * gas_density
        * relative_speed
        / (8.0 * droplet_density)
    )


def number_density(radius: float, bulk_density: float, droplet_density: float = RHO_WATER) -> float:
    """Droplets per m^3 at a given bulk density [m^-3]."""
    mass = (4.0 / 3.0) * math.pi * radius**3 * droplet_density
    return bulk_density / mass


def interception_count(
    radius: float,
    bulk_density: float,
    front_radius: float,
    path_length: float,
    droplet_density: float = RHO_WATER,
) -> float:
    """How many droplets the front sweeps over `path_length` -- the discreteness check.

    A front is a continuum to the cloud when this is large: the swept mass is then the column
    integral rather than a Poisson draw, and `k` is a number rather than a distribution.
    """
    volume = math.pi * front_radius**2 * path_length
    return number_density(radius, bulk_density, droplet_density) * volume


def projectile_radius(mass: float, density: float = RHO_WATER) -> float:
    """Radius [m] of a compact spherical projectile of `mass` -- the lower bound on arrival size."""
    return float((3.0 * mass / (4.0 * math.pi * density)) ** (1.0 / 3.0))


def full_bore_slug_ratio(
    bulk_density: float = BAG_RHO,
    bore_radius: float = BORE_RADIUS,
    length: float = BORE_LENGTH,
    projectile_mass: float = PROJECTILE_MASS,
) -> float:
    """`k = rho A L / m` -- the paper's snowplow, front spanning the bore throughout."""
    return bulk_density * math.pi * bore_radius**2 * length / projectile_mass


@dataclass(frozen=True)
class SnowplowResult:
    """One integration of the snowplow with a front that has to grow to the bore."""

    slug_ratio: float
    front_radius_end: float
    transit_time: float
    exit_speed: float


def snowplow(
    closing_speed: float,
    front_radius: float,
    expansion_speed: float,
    bulk_density: float = BAG_RHO,
    bore_radius: float = BORE_RADIUS,
    length: float = BORE_LENGTH,
    projectile_mass: float = PROJECTILE_MASS,
    steps: int = 20000,
) -> SnowplowResult:
    """Perfectly inelastic snowplow with a laterally expanding front.

        m v = m_0 v_0                       momentum, the projectile shares with what it sweeps
        dm/dx = rho pi min(r, R_bore)^2     it can only sweep what it spans
        dr/dx = c_exp / v                   lateral growth, converted to a per-length rate

    `expansion_speed = 0` with `front_radius = bore_radius` reproduces the paper's `rho A L / m`
    exactly, which is the calibration.

    **The deceleration helps.** `dr/dx = c_exp/v` steepens as the projectile slows, so filling the
    bore is a positive feedback -- which is why the shortfall is a factor of two or three rather
    than the factor of 260 the raw area ratio would suggest.
    """
    mass, speed, radius, time = projectile_mass, closing_speed, front_radius, 0.0
    dx = length / steps
    for _ in range(steps):
        mass += bulk_density * math.pi * min(radius, bore_radius) ** 2 * dx
        speed = projectile_mass * closing_speed / mass
        radius += expansion_speed / speed * dx
        time += dx / speed
    return SnowplowResult(
        slug_ratio=(mass - projectile_mass) / projectile_mass,
        front_radius_end=min(radius, bore_radius),
        transit_time=time,
        exit_speed=speed,
    )


def required_arrival_radius(
    closing_speed: float,
    expansion_speed: float,
    target_slug_ratio: float | None = None,
    fraction: float = 0.95,
    tolerance: float = 1.0e-4,
    bore_radius: float = BORE_RADIUS,
) -> float:
    """Front radius [m] the projectile must arrive with to reach `fraction` of the target `k`.

    **This is the constructive form of the finding, and the useful one.** Rather than reporting
    that a compact projectile falls short, it converts the paper's `k` into a stated requirement
    on the projectile: it must arrive spanning **74-97% of the bore radius** depending on how fast
    the front spreads. Bisection, because `snowplow` is monotone in the arrival radius.

    Raises if even a bore-filling arrival cannot reach the target -- that would mean the target is
    not achievable at this geometry rather than that the projectile is too small.
    """
    target = full_bore_slug_ratio() if target_slug_ratio is None else target_slug_ratio
    wanted = fraction * target
    if snowplow(closing_speed, bore_radius, expansion_speed).slug_ratio < wanted:
        raise ValueError(f"even a bore-filling front cannot reach k = {wanted:g}")

    lo, hi = tolerance, bore_radius
    while hi - lo > tolerance:
        mid = 0.5 * (lo + hi)
        if snowplow(closing_speed, mid, expansion_speed).slug_ratio < wanted:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> None:
    """Study 2's answer, and the question it turned up. Analytic throughout -- no artifact."""
    transit = 2.3e-3
    print(f"python: droplet coupling at rho = {BAG_RHO} kg/m^3, transit {transit * 1e3:.1f} ms")
    print(f"  {'droplet':>10} {'tau_drag':>10} {'tau_vap':>10} {'transit/tau':>12} {'verdict':>10}")
    for radius, label in ((1e-6, "1 um"), (1e-5, "10 um"), (1e-4, "100 um"), (1e-3, "1 mm")):
        drag = drag_time(radius, BAG_RHO, 45.58e3)
        vap = vaporisation_time(radius, BAG_RHO, 45.58e3)
        number = coupling_number(radius, transit, BAG_RHO, 45.58e3)
        print(
            f"  {label:>10} {drag * 1e6:9.3f}u {vap * 1e6:9.4f}u {number:12.4g} "
            f"{'couples' if number > 1.0 else 'LAGS':>10}"
        )
    limit = limiting_radius(transit, BAG_RHO, 45.58e3)
    print(
        f"  limiting droplet radius {limit * 1e3:.1f} mm (drag binds; vaporisation is 1000x faster)"
    )
    print(
        f"  discreteness never binds: a {BORE_RADIUS:g} m front over {BORE_LENGTH:g} m intercepts "
        f"{interception_count(1e-2, BAG_RHO, BORE_RADIUS, BORE_LENGTH):.3g} even at 1 cm"
    )

    print()
    print("  the question Study 2 turned up: does the front span the bore at all?")
    print(
        f"  paper's full-bore k = {full_bore_slug_ratio():.2f}; a compact 25 kg projectile "
        f"is {projectile_radius(PROJECTILE_MASS):.3f} m"
    )
    print(
        f"  {'w [km/s]':>9} {'c_exp':>7} {'k_eff':>7} {'of paper':>9} {'r_end':>7} {'transit':>9}"
    )
    for closing in (45.58e3, 75.0e3):
        for c_exp in (3.0e3, 5.0e3, 8.0e3):
            r = snowplow(closing, projectile_radius(PROJECTILE_MASS), c_exp)
            print(
                f"  {closing / 1e3:9.2f} {c_exp / 1e3:6.0f}k {r.slug_ratio:7.2f} "
                f"{r.slug_ratio / full_bore_slug_ratio():9.2f} {r.front_radius_end:7.2f} "
                f"{r.transit_time * 1e3:8.2f}ms"
            )
    print()
    print("  stated constructively -- the arrival radius k = 8.5 REQUIRES:")
    print(f"  {'c_exp':>7} {'w [km/s]':>9} {'r_arrival':>10} {'r/R_bore':>9}")
    for c_exp in (0.0, 3.0e3, 5.0e3):
        for closing in (45.58e3, 75.0e3):
            need = required_arrival_radius(closing, c_exp)
            print(
                f"  {c_exp / 1e3:6.0f}k {closing / 1e3:9.2f} {need:10.3f} {need / BORE_RADIUS:9.2f}"
            )
    print("  -> k = 8.5 is a requirement on the projectile's ARRIVAL radius, not a snowplow")
    print("     output. That radius is aim's to state. See Q-Q.")


if __name__ == "__main__":
    main()
