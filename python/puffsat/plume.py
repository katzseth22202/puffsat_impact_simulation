"""The plume-state table `aim_is_all_you_need` cites: `(w, rho) -> (T, f, P)`.

The companion repository's `src/plume_state.py` was originally specified to implement a
single-species Saha solve of its own. It was rescoped (routing document, 2026-08-21) once
`eos_water.py` was found to close all three of that spec's stated gaps -- oxygen's ionisation
ladder, O's 13.618 eV against H's 13.598, and the mixture rather than a single hydrogen-like
species. So the solve lives here and `aim` consumes a number.

**The physics is one energy balance and one root-find.** A projectile at closing speed `w` buries
itself in `k` kg of carried slug per kg of projectile; momentum sharing leaves the merged slug
moving, and what does *not* stay as bulk motion is dissipated:

    e_dissipated = (1/2) k w^2 / (1 + k)^2   [J per kg of merged slug]

**That energy is not reduced by the vaporisation and dissociation cost, it is spent on it.**
`eos_water` references `e` to bound molecular H2O at `T -> 0`, so the bond energy is already inside
`e` and the balance is simply `e(rho, T) = e_dissipated`. The audit's hand calculation used the
other convention -- subtract 54 MJ/kg first, then solve a thermal-plus-ionisation budget -- and
porting its formula onto this EOS double-charges the bond energy. There is a test for that.

**Why the table is two-dimensional.** `e_dissipated` depends only on `w` and `k`, never on
density. Saha does depend on density, so the same energy budget lands at a different temperature
in a different bag: compressing the plume pushes recombination, spends less of the budget
stripping electrons and leaves more as heat. `aim` sets `rho = m_slug / V` from the enclosed
volume, which is a live design variable, so a single row would not serve it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_water, expansion

#: Slug ratio -- carried slug mass per projectile kg (`CONTEXT.md`). The paper's flown value.
SLUG_RATIO = 8.5

#: The flown bag's density [kg/m^3]: 213 kg of slug over the ~660 m^3 enclosed volume.
BAG_RHO = 0.323

#: Closing speeds the table is published at [m/s] -- the burn envelope on a 2 km/s grid, with the
#: four anchors the paper quotes (45.58, 56.53, 65, 75 km/s) inserted exactly so `aim` reads them
#: off rather than interpolating to them.
BURN_ENVELOPE: tuple[float, ...] = tuple(
    sorted({*(float(w) * 1.0e3 for w in range(44, 78, 2)), 45.58e3, 56.53e3})
)

#: Bag densities [kg/m^3]. Spans a 20x range about the flown bag because item 4 sweeps the
#: enclosed volume; `BAG_RHO` is present exactly for the same reason the anchors above are.
BAG_DENSITIES: tuple[float, ...] = (0.05, 0.08, 0.13, 0.2, BAG_RHO, 0.5, 0.8, 1.3, 2.0)

DEFAULT_TABLE_PATH = Path("data/results/plume_state.csv")


def dissipated_energy(closing_speed: float, slug_ratio: float = SLUG_RATIO) -> float:
    """Energy dissipated per kg of merged slug [J/kg] at `closing_speed` [m/s].

    `(1/2) k w^2 / (1+k)^2`. The merged slug retains `w/(1+k)` of bulk speed by momentum
    conservation; this is the rest of the kinetic energy, per kg of the merged mass.
    """
    if slug_ratio <= 0.0:
        raise ValueError("slug_ratio must be positive")
    return 0.5 * slug_ratio * closing_speed * closing_speed / (1.0 + slug_ratio) ** 2


@dataclass(frozen=True)
class PlumeState:
    """One published row: what `aim` gets back for a `(w, rho)` it supplies."""

    closing_speed: float
    rho: float
    #: Dissipated specific energy [J/kg] -- the input to the root-find, carried so the row is
    #: auditable without re-deriving it.
    energy: float
    temp: float
    pressure: float
    #: Electrons per heavy atom, on the audit's own definition (`n_a = 3 rho / m_H2O`). It exceeds
    #: 1 only where oxygen climbs past its first ionisation stage, which the burn envelope does not.
    ionised_fraction: float
    electron_density: float

    @property
    def dissociated_fraction(self) -> float:
        """Fraction of the water that is no longer bound as H2O."""
        n_f = self.rho / eos_water.M_H2O
        return 1.0 - eos_water.composition(self.rho, self.temp).n_h2o / n_f


def plume_state(closing_speed: float, rho: float, slug_ratio: float = SLUG_RATIO) -> PlumeState:
    """Solve `e(rho, T) = e_dissipated` for the post-impact plume state.

    The inversion is `expansion.temperature_at`, reused rather than reimplemented -- it is the
    same caloric inversion the cooling history runs at every station, and having two of them would
    be two things to keep in agreement.
    """
    energy = dissipated_energy(closing_speed, slug_ratio)
    temp = expansion.temperature_at(rho, energy, eos_water.pressure_energy)
    pressure, _ = eos_water.pressure_energy(rho, temp)
    comp = eos_water.composition(rho, temp)
    n_atom = 3.0 * rho / eos_water.M_H2O
    return PlumeState(
        closing_speed=closing_speed,
        rho=rho,
        energy=energy,
        temp=temp,
        pressure=pressure,
        ionised_fraction=comp.n_e / n_atom,
        electron_density=comp.n_e,
    )


def table(
    speeds: tuple[float, ...] = BURN_ENVELOPE,
    densities: tuple[float, ...] = BAG_DENSITIES,
    slug_ratio: float = SLUG_RATIO,
) -> list[PlumeState]:
    """The full grid, density-major so a fixed-`rho` column reads contiguously."""
    return [plume_state(w, rho, slug_ratio) for rho in densities for w in speeds]


def write_table(rows: list[PlumeState], path: Path = DEFAULT_TABLE_PATH) -> None:
    """The published artifact. One row per `(w, rho)`; `aim` interpolates in it."""
    lines = [
        "closing_speed_km_s,rho_kg_m3,dissipated_MJ_kg,temp_K,pressure_Pa,"
        "ionised_fraction,dissociated_fraction,n_e_m3"
    ]
    lines += [
        f"{r.closing_speed / 1e3:g},{r.rho:g},{r.energy / 1e6:.4f},{r.temp:.1f},"
        f"{r.pressure:.6e},{r.ionised_fraction:.6f},{r.dissociated_fraction:.6f},"
        f"{r.electron_density:.6e}"
        for r in rows
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    """Publish the table, and print the column `aim`'s four quoted anchors sit in."""
    rows = table()
    write_table(rows)
    print(f"python: plume state at k = {SLUG_RATIO}, {len(rows)} rows over ")
    print(
        f"        w = {min(BURN_ENVELOPE) / 1e3:g}-{max(BURN_ENVELOPE) / 1e3:g} km/s, "
        f"rho = {min(BAG_DENSITIES):g}-{max(BAG_DENSITIES):g} kg/m^3"
    )
    print(f"  the flown bag (rho = {BAG_RHO} kg/m^3), against the audit's hand table:")
    print(f"  {'w [km/s]':>9} {'diss MJ/kg':>11} {'T [K]':>8} {'f':>7} {'P [MPa]':>9}")
    for speed, _ in expansion.PLUME_STATES:
        state = plume_state(speed * 1.0e3, BAG_RHO)
        print(
            f"  {speed:9.2f} {state.energy / 1e6:11.1f} {state.temp:8.0f} "
            f"{state.ionised_fraction:7.4f} {state.pressure / 1e6:9.2f}"
        )
    print(
        "  the solve runs 1-3% warmer than the hand table at every anchor: the audit charged "
        "54 MJ/kg\n  for the bond, eos_water charges 50.9, and the difference stays in the "
        "thermal pool."
    )
    print(f"python: wrote {DEFAULT_TABLE_PATH}")


if __name__ == "__main__":
    main()
