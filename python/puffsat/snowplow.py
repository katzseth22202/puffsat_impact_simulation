"""Is the field weakening that buys detachment affordable? (the first discovered question)

`docs/nozzle_asks_answered.md` P9 found a way out of P3's detachment failure: the exit field is
roughly 15x over-strength for the expansion it steers, and flaring to `A/A*` between 11.3 and 14.8
drops it into a window that both contains the plume and lets it go. That fix is not free, and the
cost was left uncomputed:

> One static field cannot be strong for the collision and weak for the expansion at the same
> station, and the paper requires it static in time. A field weakened to the release ceiling sits
> **3-8x under the collision's snowplow pressure** there, so the snowplow gets past. Whether that
> is affordable in ablation is a collision-phase calculation nobody has done.

This is that calculation, and it is deliberately posed as a **budget** rather than a simulation.

# Why a budget rather than a model

Predicting how much of an uncontained snowplow reaches the liner needs the collision resolved
against a wall -- the deferred full-bore run. But the question "is it affordable" does not: the
paper already books the liner's ablation allowance, and that allowance converts directly into a
**share of pulse energy the liner can absorb**. If the uncontained snowplow must deliver more than
that share, the fix fails on the budget alone and no simulation is needed to say so.

    liner budget = 4.9 kg/pulse x 59.7 MJ/kg = 293 MJ, against a 62.9 GJ pulse -> 0.47%.

**Under half a percent of the pulse may reach the liner.** That is the number every result here is
measured against, and it is the paper's own booking, not an assumption made here.

# What weakening the field actually does

Flux conservation ties the exit field to the flare, `B = B_throat/(A/A*)`, so flaring past the
flown 4 lowers the field the *collision* sees at that station too. The snowplow's design pressure
is unchanged -- it is set by the swept mass, not by the magnet -- so

    beta_snowplow = p_snowplow / (B_weak^2/2mu0) = (B_design/B_weak)^2 = (A/A*_new / 4)^2.

Above 1 the field cannot stand the front off, and the front expands until pressure balances. For a
cylindrical adiabatic expansion `p ~ r^{-2 gamma}`, so the radius grows by `beta^{1/(2 gamma)}`.
Whether that reaches the wall is a geometry question, and it runs straight into the bore ambiguity
of P7 -- so both readings are carried.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

from puffsat import field

MU0 = field.MU0
GAMMA = 5.0 / 3.0

DEFAULT_OUTPUT = Path("data/results/nozzle_snowplow.csv")

BOOKED_ABLATION_KG = 4.9
"""Carbon the paper books the liner shedding per pulse (`sec:watering_it_down` appendix)."""

SUBLIMATION_ENTHALPY = 716.7e3 / 0.012011
"""Graphite -> C(g), 59.7 MJ/kg."""

PULSE_ENERGY_J = 62.9e9
"""Flown Jupiter-only pulse at 75 km/s (25 kg into 213 kg, reduced mass 22.37 kg)."""

FLOWN_AREA_RATIO = 4.0
BORE_RADIUS_M = 3.02
"""`eq:bore_from_length`'s cylinder. The flared reading is `bore_radius_flared`."""

FRONT_TEMPERATURE_K = 94_600.0
"""Shocked-layer temperature at the snowplow's nose, 45.58 km/s (`sec:needle_through_fog`).

The paper cites **this repository** for it: `v^2/2` = 1040 MJ/kg inverted through the dissociation
and ionisation ladder. Quoted rather than recomputed, and it matters because radiation goes as
`T^4` -- a front this hot radiates 4.5 TW/m^2.
"""

SIGMA_SB = 5.670374419e-8


def liner_energy_budget() -> float:
    """Energy the liner may absorb per pulse [J] before exceeding its booked ablation.

    `4.9 kg x 59.7 MJ/kg`. The paper's own number, turned into the currency this question needs.
    """
    return BOOKED_ABLATION_KG * SUBLIMATION_ENTHALPY


def liner_budget_share(pulse_energy: float = PULSE_ENERGY_J) -> float:
    """The budget as a share of pulse energy -- **0.47%**, and the gate for everything below."""
    return liner_energy_budget() / pulse_energy


def snowplow_beta(area_ratio: float, flown_area_ratio: float = FLOWN_AREA_RATIO) -> float:
    """How far past containment the collision front sits once the field is flared to `A/A*`.

    `B = B_throat/(A/A*)` by flux conservation and the snowplow pressure does not move, so
    `beta = (A/A*_new / A/A*_flown)^2`. The flown flare is the definition of standoff here, so it
    returns 1 by construction.
    """
    return float((area_ratio / flown_area_ratio) ** 2)


def containment_radius_ratio(beta: float, gamma: float = GAMMA) -> float:
    """Radius growth of an uncontained front expanding to pressure balance: `beta^{1/(2 gamma)}`.

    Cylindrical and adiabatic: area goes as `r^2`, `p ~ V^{-gamma}`, so `p ~ r^{-2 gamma}` and the
    front expands until it has shed the factor `beta`. Modest — `beta` = 8 is only 87% in radius,
    because a `5/3` gas stiffens fast.
    """
    if beta <= 1.0:
        return 1.0
    return float(beta ** (1.0 / (2.0 * gamma)))


def bore_radius_flared(area_ratio: float, throat_radius: float = BORE_RADIUS_M / 2.0) -> float:
    """Bore radius under the *flared* reading of P7: `r = r_throat sqrt(A/A*)`.

    The throat is taken at half the cylinder's radius so the flown `A/A*` = 4 fills the 3.02 m bore
    at the exit — the only reading in which the cylinder and the area ratio describe one object.
    """
    return throat_radius * math.sqrt(area_ratio)


@dataclass(frozen=True)
class ContainmentCase:
    """One flare, and what it does to the collision front at the exit station."""

    area_ratio: float
    exit_field_t: float
    snowplow_beta: float
    radius_growth: float
    front_radius_cylinder_m: float
    """Where the front ends up, against `eq:bore_from_length`'s fixed 3.02 m wall."""
    front_radius_flared_m: float
    """The same, against a bore that flares with the field."""
    touches_cylinder: bool
    touches_flared: bool


def containment_cases(
    area_ratios: tuple[float, ...] = (4.0, 6.0, 8.0, 11.3, 14.8, 20.0),
    *,
    throat_field_t: float = 20.0,
) -> list[ContainmentCase]:
    """Sweep the flare and ask, at each, whether the collision front still clears the wall.

    **Both bore readings give the same answer, and it is not the one P9 assumed.** The front
    already fills the bore for the last three quarters of the crossing (`sec:needle_through_fog`:
    a 0.15 m arrival "reaches the wall of a 3 m bore after about 6 m of a 23 m column"), so once
    the field is weakened it must exceed whatever bore contains it. Formally, with the front
    filling the bore before weakening,

        r_front/r_bore = beta^{1/(2 gamma)} = (A/A*_new / 4)^{0.6} > 1 for any flare past 4,

    and the bore's own growth cancels out. A flared bore does not rescue it: the wall grows as
    `sqrt(A/A*)` and the front grows as `sqrt(A/A*) x (A/A*/4)^{0.6}`, always faster.

    **So flaring past the flown 4 puts the collision front on the liner, in every reading of P7.**
    """
    out: list[ContainmentCase] = []
    for ar in area_ratios:
        beta = snowplow_beta(ar)
        growth = containment_radius_ratio(beta)
        # The front fills whatever bore it is in before the field is weakened.
        r_cyl = BORE_RADIUS_M * growth
        r_flare = bore_radius_flared(ar) * growth
        out.append(
            ContainmentCase(
                area_ratio=ar,
                exit_field_t=throat_field_t / ar,
                snowplow_beta=beta,
                radius_growth=growth,
                front_radius_cylinder_m=r_cyl,
                front_radius_flared_m=r_flare,
                touches_cylinder=r_cyl > BORE_RADIUS_M * 1.001,
                touches_flared=r_flare > bore_radius_flared(ar) * 1.001,
            )
        )
    return out


def front_radiative_flux() -> float:
    """`sigma T^4` at the shocked layer's temperature — 4.5 TW/m^2 at 94 600 K.

    The reason contact is not survivable if it happens: this is five orders of magnitude past what
    graphite re-radiates at sublimation (13.1 MW/m^2).
    """
    return SIGMA_SB * FRONT_TEMPERATURE_K**4


def contact_share_allowed(pulse_energy: float = PULSE_ENERGY_J) -> float:
    """Share of pulse energy the liner may take — the same 0.47%, named for the contact question."""
    return liner_budget_share(pulse_energy)


def ablation_if_contacted(contact_share: float, pulse_energy: float = PULSE_ENERGY_J) -> float:
    """Carbon lost per pulse [kg] if `contact_share` of the pulse energy reaches the liner."""
    return contact_share * pulse_energy / SUBLIMATION_ENTHALPY


def main() -> None:
    """Does weakening the field to buy detachment cost more liner than the budget allows?"""
    parser = argparse.ArgumentParser(description="The snowplow containment budget")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    budget = liner_energy_budget()
    print("== The liner's own budget, in the currency this question needs ==")
    print(f"booked ablation      : {BOOKED_ABLATION_KG} kg/pulse")
    print(f"sublimation enthalpy : {SUBLIMATION_ENTHALPY / 1e6:.1f} MJ/kg")
    print(f"energy budget        : {budget / 1e6:.0f} MJ")
    print(f"pulse energy         : {PULSE_ENERGY_J / 1e9:.1f} GJ")
    print(
        f"-> the liner may absorb {100 * liner_budget_share():.2f}% of the pulse. That is the gate."
    )

    print("\n== What flaring does to the COLLISION front at the exit station ==")
    print("Flux conservation lowers the field the collision sees too; its pressure does not move.")
    print(
        f"\n{'A/A*':>7} {'B exit':>8} {'beta_snowplow':>14} {'radius x':>9} "
        f"{'front r (cyl)':>14} {'touches?':>9} {'front r (flared)':>17} {'touches?':>9}"
    )
    rows = containment_cases()
    for c in rows:
        print(
            f"{c.area_ratio:7.1f} {c.exit_field_t:7.2f}T {c.snowplow_beta:14.2f} "
            f"{c.radius_growth:9.2f} {c.front_radius_cylinder_m:13.2f}m "
            f"{('YES' if c.touches_cylinder else 'no'):>9} "
            f"{c.front_radius_flared_m:16.2f}m "
            f"{('YES' if c.touches_flared else 'no'):>9}"
        )

    print("\n== If it touches ==")
    flux = front_radiative_flux()
    graphite = SIGMA_SB * 3900.0**4
    print(
        f"shocked-layer temperature (the paper's, cited to this repo): {FRONT_TEMPERATURE_K:.0f} K"
    )
    print(f"its radiative flux        : {flux / 1e12:.2f} TW/m^2")
    print(f"graphite's own ceiling    : {graphite / 1e6:.1f} MW/m^2")
    print(f"-> {flux / graphite:.0f}x past what the liner can shed. Contact is not")
    print("   a survivable regime; the question is only whether it happens.")

    print("\n== The budget, restated as a limit on contact ==")
    for share in (0.0047, 0.01, 0.05, 0.10):
        kg = ablation_if_contacted(share)
        verdict = "within budget" if kg <= BOOKED_ABLATION_KG * 1.01 else "over"
        print(f"  {100 * share:5.2f}% of the pulse to the liner -> {kg:7.1f} kg/pulse  ({verdict})")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "area_ratio",
                "exit_field_T",
                "snowplow_beta",
                "radius_growth",
                "front_r_cylinder_m",
                "touches_cylinder",
                "front_r_flared_m",
                "touches_flared",
            ]
        )
        for c in rows:
            writer.writerow(
                [
                    f"{c.area_ratio:g}",
                    f"{c.exit_field_t:.4f}",
                    f"{c.snowplow_beta:.4f}",
                    f"{c.radius_growth:.4f}",
                    f"{c.front_radius_cylinder_m:.4f}",
                    int(c.touches_cylinder),
                    f"{c.front_radius_flared_m:.4f}",
                    int(c.touches_flared),
                ]
            )
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
