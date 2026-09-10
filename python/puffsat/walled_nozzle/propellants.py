"""N10 item 3: the propellant ladder at a fixed chamber temperature, on solved chemistry.

The ask asks for "a water slug at `k = 37.70` and for pure hydrogen at `k = 7.99`, so the ladder
in ADR-0016 rests on solved chemistry rather than on the equilibrium assumption". This module
runs that ladder like for like: the same fixed point, the same chamber, the same nozzle, and
**each fluid's own equilibrium EOS**, so the conversion fraction is solved per fluid rather than
borrowed from methane.

**The identity is EOS-agnostic and that is the whole trick.** ADR-0016's chamber closes on

    rho = (1 + k) m_impactor / V ;   u = e(rho, T) ;   k = (w^2 / 2) / u - 1

which mentions no species at all. Point it at `eos_water` or `hydrogen` instead of `eos_methane`
and the ladder falls out. `expansion.cooling_history` is likewise parametrised by the EOS and the
sound speed, so the same nozzle runs every fluid.

**Why lighter wins twice, and then gives some of it back.** A light fluid stores more energy per
kilogram, because energy at 10 kK is mostly `3/2 kT` per *particle* and a kilogram of a light gas
is more particles. Storing more per kilogram means less mass is needed to absorb the pulse, so `k`
falls -- and the exhaust speed `w / sqrt(1+k)` rises. Then `k` enters a second time through the
carried-mass credit `(1+k)/k`, which also rewards a small `k`; on that convention alone the ladder
follows `Isp ~ 1/sqrt(mean atomised particle mass)` to a few percent. But `k` enters a **third**
time, through the head-on momentum debit `w/(k g0)`, and that one *punishes* a small `k`. So the
corrected ladder is compressed relative to the scaling law: hydrogen still wins, by 1.55x over
methane rather than the 1.88x the credit-only column shows.

**Three Isp quantities, and they must not be mixed** (the mistake W5 records). All three come from
`chamber.IspLedger`, which is the single place the head-on sign is written down:

- `isp_true` (alias `isp_total`) -- **real Isp**: exhaust speed over `g0`, per kilogram of
  *everything* expelled. Convention-free.
- `isp_effective` -- **the figure of merit**: per kilogram of *carried slug*, with the free
  impactor credited **and the arriving momentum debited**. Below `isp_true` for every fluid here.
- `isp_carried_gross` -- the credit-only intermediate, `(1+k)/k` times `isp_true`. **Not a figure
  of merit**; retained because it is what W5 and W8 first published and because the `1/sqrt(mean
  mass)` scaling law is a property of this convention.

None carries the launch-ledger normalisations this repository does not own (`eta_geom`, the
drift term, vessel mass), so **these numbers are comparable across fluids but not against the
paper's absolute figures.**

**Ammonia is an estimate and is labelled one everywhere.** There is no `eos_ammonia`; nitrogen is
not in any EOS here. Its `u` is assembled from the ask's own atomisation energy plus an exactly
computed translational term, and its conversion fraction is swept rather than solved. It is
carried because the ask asks for ammonia, not because this study can price it.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_methane, eos_water, expansion
from puffsat.walled_nozzle import chamber, hydrogen

#: Standard gravity, taken from `chamber` so every Isp in the study divides by the same one.
G0 = chamber.G0

#: Nozzle the ladder is run through: ADR-0016's 3 m bore over its widest throat, the
#: configuration the ask itself states. W6 argues for a narrower one; the ladder is run at the
#: ask's own geometry so the comparison is against its numbers rather than against our advice.
BORE_AREA = math.pi * chamber.BORE_RADIUS**2
THROAT_AREA = 7.0
NOZZLE_LENGTH = chamber.CHAMBERS[0][1]
VOLUME = chamber.CHAMBERS[0][0]

DEFAULT_OUTPUT = Path("data/results/walled_nozzle/propellants.csv")

#: `(rho, T) -> (p, e)`, the same alias `expansion` uses.
Eos = Callable[[float, float], tuple[float, float]]


@dataclass(frozen=True)
class Rung:
    """One propellant at the flown chamber temperature."""

    name: str
    #: Mean mass [amu] of a particle once fully atomised -- the variable the ladder sorts on.
    mean_atomised_mass: float
    slug_ratio: float
    rho: float
    #: Wall-cap energy density [J/kg]: what a kilogram holds at the chamber temperature.
    energy: float
    exit_temp: float
    exhaust_speed: float
    #: Share of `u` that leaves as directed kinetic energy. Convention-free.
    conversion: float
    #: `True` where the rung is assembled rather than solved (ammonia).
    estimated: bool

    @property
    def isp(self) -> chamber.IspLedger:
        """Every specific impulse for this rung. The algebra lives in `chamber.IspLedger`."""
        return chamber.isp_ledger(self.exhaust_speed, self.slug_ratio)

    @property
    def isp_total(self) -> float:
        """Per kilogram of everything expelled [s]."""
        return self.isp.isp_true

    @property
    def isp_true(self) -> float:
        """Alias for `isp_total`. **W8 prints this column as "Isp true"**, so the name the paper
        will carry resolves here rather than looking like a missing quantity."""
        return self.isp_total

    @property
    def carried_slug_mass(self) -> float:
        """Propellant the vehicle actually lifted [kg] -- the slug alone. The impactor arrives
        from outside at 75 km/s and was never carried, which is the whole reason the effective
        convention exists."""
        return self.slug_ratio * chamber.IMPACTOR_MASS

    @property
    def free_impactor_bonus(self) -> float:
        """`isp_carried_gross / isp_true - 1 = 1/k`: the free ride from mass the vehicle did not
        carry. **It grows as the required slug shrinks**, so it rewards hydrogen (+12.9%) far
        more than water (+2.5%) -- a second, independent advantage on top of the chemistry.

        **It is not the whole mass ledger**: the head-on debit runs as `1/k` too, and is larger.
        """
        return self.isp.free_impactor_bonus

    @property
    def momentum_debit(self) -> float:
        """`w/(k g0)` [s]: the momentum the arriving projectile costs, charged against thrust."""
        return self.isp.momentum_debit

    @property
    def isp_carried_gross(self) -> float:
        """Per kilogram of carried slug [s], credit only. **Not a figure of merit** -- see the
        module docstring and `chamber.IspLedger`. This is the column W8 first published."""
        return self.isp.isp_carried_gross

    @property
    def isp_effective(self) -> float:
        """**The figure of merit** [s]: carried-mass credit *and* head-on momentum debit."""
        return self.isp.isp_effective

    @property
    def eta_jet(self) -> float:
        """The companion's jet efficiency, `u_e/(w/sqrt(1+k))` -- identically `sqrt(conversion)`."""
        return self.isp.eta_jet

    @property
    def ideal_exhaust_speed(self) -> float:
        """`w / sqrt(1+k)` [m/s]: the speed at full conversion, ADR-0016's own identity."""
        return self.isp.ideal_exhaust_speed


def solve_slug_ratio(pe: Eos, temp: float, volume: float = VOLUME) -> tuple[float, float, float]:
    """`(k, rho, u)` at one chamber, for any EOS. Same contraction `chamber.solve_chamber` uses."""
    k = chamber.ASK_SLUG_RATIO
    rho = (1.0 + k) * chamber.IMPACTOR_MASS / volume
    u = pe(rho, temp)[1]
    for _ in range(400):
        rho = (1.0 + k) * chamber.IMPACTOR_MASS / volume
        u = pe(rho, temp)[1]
        k_new = chamber.SPECIFIC_PULSE_ENERGY / u - 1.0
        if abs(k_new - k) < 1e-12 * k:
            k = k_new
            break
        k = 0.5 * (k + k_new)
    else:  # pragma: no cover - the map is a contraction for every EOS carried here
        raise RuntimeError("propellant fixed point did not converge")
    rho = (1.0 + k) * chamber.IMPACTOR_MASS / volume
    return k, rho, pe(rho, temp)[1]


def solved_rung(
    name: str,
    pe: Eos,
    sound_speed: Callable[[float, float], float],
    mean_atomised_mass: float,
    temp: float = chamber.FLOWN_TEMPERATURE,
    throat_area: float = THROAT_AREA,
    steps: int = 128,
) -> Rung:
    """Run one fluid end to end: fixed point, then the nozzle, on its own EOS.

    **Valid at the ask's own wide throat and not down the ladder.** Two simplifications are fine
    at `A/A*` ~ 4 and wrong deeper: the nozzle length is held at the ask's `NOZZLE_LENGTH` rather
    than rebuilt on the cone, and the conversion is the raw expansion with **no freeze cap and no
    Damkoehler verdict**. For any throat narrower than the ask's, use `surface.column`, which does
    both -- it is what `wall.write_throat_life` and W8's recommended-geometry rung are built on.
    """
    k, rho, u = solve_slug_ratio(pe, temp)
    rows = expansion.cooling_history(
        rho, temp, pe, sound_speed, BORE_AREA / throat_area, NOZZLE_LENGTH, steps=steps
    )
    speed = rows[-1].speed
    return Rung(
        name=name,
        mean_atomised_mass=mean_atomised_mass,
        slug_ratio=k,
        rho=rho,
        energy=u,
        exit_temp=rows[-1].temp,
        exhaust_speed=speed,
        conversion=speed**2 / (2.0 * u),
        estimated=False,
    )


# ---- Ammonia, assembled rather than solved -------------------------------------------------------

#: NH3 molar mass [kg]. Nitrogen appears in no EOS here, which is the whole reason for this branch.
M_NH3 = (14.007 + 3.0 * 1.008) * eos_methane.AMU
#: Atomisation energy [J/kg]. The ask's own figure, and it matches 3 x the N-H bond energy.
ATOMISATION_NH3 = 68.9e6
#: Ionisation share of `u` at 10 kK, taken from the solved methane chamber, where it is **0.45%**
#: -- small enough that assuming it for ammonia costs less than a percent on `k`. Recorded because
#: an earlier estimate assumed 8 MJ/kg here, an order of magnitude too high, and moved the answer.
IONISATION_SHARE = 0.005
#: Fraction of the atomisation store charged at the chamber state, from methane's solved 0.947.
CHARGED_FRACTION = 0.95


def ammonia_energy(temp: float = chamber.FLOWN_TEMPERATURE) -> float:
    """Estimated wall-cap energy density [J/kg] for ammonia at `temp`.

    Chemical + translational + a small ionisation allowance. The translational term is *exact*
    given full atomisation -- `NH3 -> N + 3H` is four particles, so a kilogram holds `4 N_A /
    M_NH3` of them, each with `3/2 kT`. The chemical term is the ask's own number. Only the
    ionisation allowance is assumed, and methane shows it is worth half a percent.
    """
    particles = 4.0 / M_NH3
    sensible = particles * 1.5 * eos_methane.K_B * temp
    chemical = CHARGED_FRACTION * ATOMISATION_NH3
    return (chemical + sensible) / (1.0 - IONISATION_SHARE)


def ammonia_rung(conversion: float, temp: float = chamber.FLOWN_TEMPERATURE) -> Rung:
    """Ammonia at an *assumed* conversion fraction. Sweep it; do not quote one value."""
    u = ammonia_energy(temp)
    k = chamber.SPECIFIC_PULSE_ENERGY / u - 1.0
    return Rung(
        name="ammonia",
        mean_atomised_mass=M_NH3 / 4.0 / eos_methane.AMU,
        slug_ratio=k,
        rho=(1.0 + k) * chamber.IMPACTOR_MASS / VOLUME,
        energy=u,
        exit_temp=math.nan,
        exhaust_speed=math.sqrt(2.0 * conversion * u),
        conversion=conversion,
        estimated=True,
    )


#: Mean atomised particle masses [amu]: `H2 -> 2H`, `CH4 -> C + 4H`, `H2O -> O + 2H`.
MEAN_MASS = {"hydrogen": 1.008, "methane": 16.043 / 5.0, "water": 18.015 / 3.0}


def ladder(temp: float = chamber.FLOWN_TEMPERATURE) -> list[Rung]:
    """The three solved rungs, lightest first."""
    return [
        solved_rung(
            "hydrogen", hydrogen.pressure_energy, hydrogen.sound_speed, MEAN_MASS["hydrogen"], temp
        ),
        solved_rung(
            "methane",
            eos_methane.pressure_energy,
            eos_methane.sound_speed,
            MEAN_MASS["methane"],
            temp,
        ),
        solved_rung(
            "water", eos_water.pressure_energy, eos_water.sound_speed, MEAN_MASS["water"], temp
        ),
    ]


CSV_HEADER = (
    "fluid,estimated,mean_atomised_amu,slug_ratio,carried_slug_kg,rho,u_j_kg,exit_temp_k,"
    "exhaust_speed_m_s,conversion,eta_jet,isp_true_s,isp_effective_s,isp_carried_gross_s,"
    "momentum_debit_s,free_impactor_bonus\n"
)


def write(rungs: list[Rung], path: Path = DEFAULT_OUTPUT) -> None:
    """Write the ladder. Committed: the asks document is read in the paper repository."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(CSV_HEADER)
        for r in rungs:
            fh.write(
                f"{r.name},{r.estimated},{r.mean_atomised_mass:.4f},{r.slug_ratio:.4f},"
                f"{r.carried_slug_mass:.2f},{r.rho:.6e},{r.energy:.6e},{r.exit_temp:.2f},"
                f"{r.exhaust_speed:.2f},{r.conversion:.5f},{r.eta_jet:.5f},{r.isp_true:.1f},"
                f"{r.isp_effective:.1f},{r.isp_carried_gross:.1f},{r.momentum_debit:.1f},"
                f"{r.free_impactor_bonus:.5f}\n"
            )


def _by_name(rungs: list[Rung], name: str) -> Rung:
    """The rung with this name. A helper so the report lines stay inside the line limit."""
    return next(r for r in rungs if r.name == name)


def main() -> None:
    """Report the propellant ladder at the flown chamber temperature."""
    print(
        f"N10 item 3: the propellant ladder at {chamber.FLOWN_TEMPERATURE:.0f} K, "
        f"{chamber.CLOSING_SPEED / 1e3:.0f} km/s\n"
    )
    print(
        f"{VOLUME:.0f} m^3 chamber, {THROAT_AREA:.0f} m^2 throat "
        f"(A/A* = {BORE_AREA / THROAT_AREA:.2f}), {NOZZLE_LENGTH} m\n"
    )
    print(
        f"{'fluid':10} {'m_bar':>6} {'k':>7} {'u MJ/kg':>8} {'exit T':>7} {'u_e':>7} "
        f"{'conv':>6} {'Isp TRUE':>9} {'Isp EFF':>8} {'(gross':>8} {'debit)':>8}"
    )

    rungs = ladder()
    for r in rungs:
        print(
            f"{r.name:10} {r.mean_atomised_mass:6.2f} {r.slug_ratio:7.2f} "
            f"{r.energy / 1e6:8.1f} {r.exit_temp:7.0f} {r.exhaust_speed:7.0f} "
            f"{r.conversion:6.3f} {r.isp_total:9.0f} {r.isp_effective:8.0f} "
            f"{r.isp_carried_gross:8.0f} {-r.momentum_debit:8.0f}"
        )

    print()
    for conv in (0.35, 0.406, 0.45):
        a = ammonia_rung(conv)
        print(
            f"{'ammonia*':10} {a.mean_atomised_mass:6.2f} {a.slug_ratio:7.2f} "
            f"{a.energy / 1e6:8.1f} {'--':>7} {a.exhaust_speed:7.0f} "
            f"{a.conversion:6.3f} {a.isp_total:9.0f} {a.isp_effective:8.0f} "
            f"{a.isp_carried_gross:8.0f} {-a.momentum_debit:8.0f}"
        )
    print("  * ESTIMATE -- no eos_ammonia; conversion assumed, not solved")
    print(
        "  Isp EFF = gross - debit: the free impactor credited, the head-on momentum charged."
        "\n  It is THE figure of merit and it is below Isp TRUE, because the debit w/(k g0)"
        "\n  exceeds the credit u_e/(k g0) whenever u_e < w -- which is everywhere here.\n"
    )

    methane = _by_name(rungs, "methane")
    print("Against 1/sqrt(mean atomised particle mass), normalised to methane:")
    for r in rungs:
        scaling = math.sqrt(methane.mean_atomised_mass / r.mean_atomised_mass)
        print(
            f"  {r.name:10} predicted {scaling:5.2f}   "
            f"gross {r.isp_carried_gross / methane.isp_carried_gross:5.2f}   "
            f"EFFECTIVE {r.isp_effective / methane.isp_effective:5.2f}"
        )
    print(
        "  The scaling law is a property of the GROSS convention. The debit compresses the"
        "\n  ladder, because it falls hardest on the small-k fluids the law rewards."
    )

    print("\nThe ask's own slug ratios, for comparison:")
    print(f"  water     ask k = 37.70   solved {_by_name(rungs, 'water').slug_ratio:.2f}")
    print(f"  hydrogen  ask k =  7.99   solved {_by_name(rungs, 'hydrogen').slug_ratio:.2f}")

    write([*rungs, *(ammonia_rung(c) for c in (0.35, 0.406, 0.45))])
    print(f"\nwrote {DEFAULT_OUTPUT}")


if __name__ == "__main__":
    main()
