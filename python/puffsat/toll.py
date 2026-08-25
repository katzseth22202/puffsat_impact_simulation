"""`eta_chem(w, k)`: the ceiling frozen chemistry puts on the magnetic nozzle's jet (Q-R).

`aim_is_all_you_need` prices the growth chain against a nozzle recovery it sweeps freely to 0.90.
The paper's own `sec:jet_efficiency` defines that efficiency to include "frozen ionization or
dissociation energy" among its five contributions -- so this module does not add a term, it
**computes one of the five** and turns a swept parameter into one with a floor under it.

**Why it is owed here rather than there.** Q-P established that the dissociation store freezes at
the nozzle lip with 90-100% of it still held, against the paper's assumption that "the energy is a
loan rather than a cost". How much is stranded depends on the stagnation state, on where the store
freezes, and on how much of the water dissociated in the first place -- all of which need the
equation of state. `aim` needs a number; this is the number.

**The closed form, and what it is a ratio of.** `aim`'s ideal gross exhaust momentum is
`m w sqrt(1+k)`: the whole collision energy `1/2 m w^2` placed on one axis, which is the
energy-conservation bound. Paying the toll first leaves `1/2 w^2/(1+k) - phi E_B` per kg of merged
mass, so the gas speed falls from `w/sqrt(1+k)` to `sqrt(w^2/(1+k) - 2 phi E_B)` and

    eta_chem = sqrt(1 - 2 phi E_B (1+k) / w^2)

`phi` is the bond fraction still held when the chemistry freezes, which is what the solver supplies
and what a closed form cannot.

**Where this multiplies, which is the part that is easy to get wrong.** It scales the *gross* jet,
`sqrt(1+k)`, and **not** the `-1` (head-on) or `+1` (overtake) that `aim`'s `beta` carries with it.
Those are the merged slug's bulk drift -- pure momentum conservation, which no chemistry can touch
(Seth, 2026-08-25). See the routing document's Q-R for the two corrected `beta` expressions.

**Density convention, stated because it is a choice.** The surface is computed at the flown bag
density, on the assumption that the bag is *resized* with `k` so the plume density stays put. That
matches the paper, which sweeps bag radius as a live design variable (`tab:bag_sizing`). Holding
the bag *volume* fixed instead would make the plume density scale with `(1 + k)` and move the
freeze; `density_sensitivity` reports how much that would matter.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_water, expansion, fireball, plume

#: Minimum stagnation temperature [K] for a node to count as a realistic operating point
#: (Seth, 2026-08-25). **This is the gate `aim` should use**, and it is better than the energy bill
#: below because it is stated in the variable that actually controls conductivity.
#:
#: **Why 10 000 K is the right height, and not too high.** Conductivity is not what fails at the
#: cold corner -- the 1% potassium seed keeps supplying electrons long after the water stops, so
#: `sigma` is still 1480 S/m at 10 000 K and 403 S/m even at 4500 K. What fails is *dissociation*:
#: below ~10 000 K the plume never breaks up, and the excluded nodes sit at `phi` = 0.38-0.75, a
#: barely-dissociated mush. The gate separates that corner cleanly, and every node it admits has
#: `phi >= 0.777`, which is the region where the closed form is a tight floor.
#:
#: **What it does NOT do.** It is not a stability gate. ADR-0038's electrothermal verdict flips
#: between `T_0` = 15 170 and 19 710 K -- *above* this line and above the flown cold anchor -- so
#: the design already operates on the unstable side at its cold end. That is Q-O, still open, and
#: folding it in here would hide it.
MIN_STAGNATION_TEMP = 10_000.0

#: The paper's ignition bill [J/kg of merged mass]: vaporizing the slug, dissociating it and
#: heating it to 15 000 K (`sec:watering_it_down`). A collision that does not clear this does not
#: produce a conducting plume, and a magnetic nozzle has nothing to grip. It is **not** the same
#: gate as `eta_chem > 0`, and at the cold, high-`k` corner it is the stricter one -- see
#: `TollPoint.ignites`.
IGNITION_ENERGY = 85.1e6

#: Full atomization energy of water [J/kg] -- the most the dissociation store can strand.
#: The paper's `sec:watering_it_down` uses 50.4 MJ/kg; this is 50.94 and they agree to 1%.
E_BOND = eos_water.FULL_ATOMIZATION_ENERGY

#: Closing speeds the surface is published at [m/s]. Spans both legs and both cadence cases, with
#: the anchors the paper and `aim` quote inserted exactly so neither has to interpolate to them:
#: 45.58 and 56.53 (3-synodic overtake, its ends), 61.83 and 65.13 (2-synodic overtake, its top),
#: 75 (head-on Jupiter departure).
SPEEDS: tuple[float, ...] = (
    40.0e3,
    45.58e3,
    50.0e3,
    56.53e3,
    61.83e3,
    65.13e3,
    70.0e3,
    75.0e3,
    80.0e3,
)

#: Slug ratios [-]. Brackets `aim`'s chain-optimiser search box from below and runs up to where the
#: chemistry forbids exhaust entirely on the cold legs; that repo's `_K_SEARCH_MAX` is 80, which is
#: well past the point where there is nothing left to expand.
SLUG_RATIOS: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0, 6.0, 8.5, 12.0, 16.0, 20.0)

#: Integration settings for the freeze solve. Cheaper than `fireball`'s shipped defaults and
#: verified against them: the freeze is triggered at the lip, so integrating three decades past it
#: costs a grid most of its runtime and moves `phi` by under 1e-4.
SCAN_STEPS = 160
SCAN_EXPANSION_RATIO = 512.0

DEFAULT_TABLE_PATH = Path("data/results/eta_chem.csv")


def available_energy(closing_speed: float, slug_ratio: float, bond_fraction: float) -> float:
    """Specific energy left for directed motion [J/kg of merged mass].

    The full collision energy per merged kg is `1/2 w^2/(1+k)` -- bulk drift plus thermal, which is
    what `aim`'s ideal bound places on one axis -- less the bond energy no bond is paying for.
    """
    return 0.5 * closing_speed**2 / (1.0 + slug_ratio) - bond_fraction * E_BOND


def eta_chem(closing_speed: float, slug_ratio: float, bond_fraction: float = 1.0) -> float:
    """`sqrt(1 - 2 phi E_B (1+k)/w^2)`, floored at zero.

    Zero means the collision does not dissipate enough to pay the bond bill at all: there is no
    exhaust, not merely a poor one. `zero_slug_ratio` gives that boundary in closed form.
    """
    if closing_speed <= 0.0 or slug_ratio <= 0.0:
        raise ValueError("closing_speed and slug_ratio must be positive")
    inner = 1.0 - 2.0 * bond_fraction * E_BOND * (1.0 + slug_ratio) / closing_speed**2
    return math.sqrt(inner) if inner > 0.0 else 0.0


def zero_slug_ratio(closing_speed: float, bond_fraction: float = 1.0) -> float:
    """The `k` above which the plume produces no exhaust at all: `1 + k = w^2/(2 phi E_B)`.

    `aim`'s chain optimiser searches `k` to 80 with no such bound, so at the cold end its box runs
    far into territory where the answer is not a poor jet but no jet.
    """
    return closing_speed**2 / (2.0 * bond_fraction * E_BOND) - 1.0


def stagnation_temperature(closing_speed: float, slug_ratio: float, rho_bag: float) -> float:
    """Stagnation temperature [K] for one `(w, k)`.

    The cheap half of a node: no freeze solve, so `admissible_slug_ratios` can bisect on it.
    """
    return expansion.temperature_at(
        rho_bag, plume.dissipated_energy(closing_speed, slug_ratio), eos_water.pressure_energy
    )


def admissible_slug_ratios(
    closing_speed: float,
    rho_bag: float = plume.BAG_RHO,
    min_temp: float = MIN_STAGNATION_TEMP,
) -> tuple[float, float] | None:
    """The `(k_lo, k_hi)` interval where the plume clears `min_temp`, or `None` if none does.

    **This is what replaces `two_wave_growth._K_SEARCH_MAX = 80`.** The dissipated energy per kg
    goes as `k/(1+k)^2`, which *peaks at k = 1* and falls away on both sides, so the admissible set
    is a closed interval rather than everything below a maximum -- the paper says so too
    (`sec:two_leg_nozzle`), and gives [0.098, 10.21] at 45.58 km/s for its own 85.1 MJ/kg bill.
    Too much slug spreads a fixed energy too thin; too little dissipates almost nothing.

    Bisected on each side of the `k = 1` peak, on the cheap stagnation solve only.
    """
    if stagnation_temperature(closing_speed, 1.0, rho_bag) < min_temp:
        return None

    def bisect(lo: float, hi: float) -> float:
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if stagnation_temperature(closing_speed, mid, rho_bag) >= min_temp:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    return bisect(1.0, 1.0e-6), bisect(1.0, 1.0e4)


@dataclass(frozen=True)
class TollPoint:
    """One node of the surface: what `aim` gets back for a `(w, k)` it supplies."""

    closing_speed: float
    slug_ratio: float
    rho_bag: float
    #: Dissipated specific energy [J/kg of merged mass], `1/2 k w^2/(1+k)^2`. Carried so a row is
    #: auditable without re-deriving it.
    dissipated: float
    #: Stagnation temperature [K] from `e(rho, T) = e_dissipated`.
    temp_0: float
    #: Where the store freezes. `None` when nothing freezes inside the integrated range.
    rho_freeze: float | None
    temp_freeze: float | None
    #: Bond fraction still held at the freeze -- the `phi` of the closed form. **Solved, not
    #: assumed**: at the cold, high-`k` corner the plume never fully atomizes, so charging the
    #: whole `E_BOND` there would invent a toll the plume does not owe.
    bond_fraction: float
    available: float
    eta: float

    @property
    def lights(self) -> bool:
        """Whether the collision dissipates enough to pay the bond bill at all."""
        return self.available > 0.0

    @property
    def conducts(self) -> bool:
        """Whether the stagnation state clears `MIN_STAGNATION_TEMP` -- the recommended gate.

        Slightly *looser* than `ignites`: it admits three nodes at 72-75 MJ/kg that the paper's
        85.1 bill rejects. That bill is defined as heating to 15 000 K, which is a design target
        rather than a floor, and this EOS puts those nodes at 10 000-11 000 K with `sigma` of
        1500-2200 S/m. Admitting them is deliberate; both flags ship so a consumer can disagree.
        """
        return self.temp_0 >= MIN_STAGNATION_TEMP

    @property
    def ignites(self) -> bool:
        """Whether the collision clears the paper's 85.1 MJ/kg bill for a conducting plume.

        **This is the gate that binds at the cold, high-`k` corner, and `eta_chem` is not.** There
        the plume never fully dissociates, so it strands less and `eta_chem` recovers to ~0.7 --
        but a barely-dissociated plume is a poor conductor, and "mass the field fails to grip" is a
        *different* term of `eta_jet` that this module does not compute. A consumer optimising on
        `eta` alone would walk straight into that corner and read a flattering number.
        """
        return self.dissipated >= IGNITION_ENERGY


def toll_point(
    closing_speed: float,
    slug_ratio: float,
    rho_bag: float = plume.BAG_RHO,
    half_angle_deg: float = fireball.DIVERGENCE_HALF_ANGLE_DEG,
) -> TollPoint:
    """Solve one `(w, k)` node: stagnation state, freeze, held bond fraction, `eta_chem`.

    The chain is three existing pieces and no new physics: `plume.dissipated_energy` for the
    budget, `expansion.temperature_at` for the stagnation state, and `fireball` for the freeze.
    That is deliberate -- a second copy of any of them would be a second thing to keep in step.
    """
    dissipated = plume.dissipated_energy(closing_speed, slug_ratio)
    temp_0 = expansion.temperature_at(rho_bag, dissipated, eos_water.pressure_energy)
    stations = fireball.scan(
        temp_0,
        half_angle_deg,
        steps=SCAN_STEPS,
        rho_0=rho_bag,
        expansion_ratio=SCAN_EXPANSION_RATIO,
    )
    freeze = fireball.freeze_state(stations)
    # Nothing freezing inside the integrated range means the chemistry kept up throughout it --
    # the paper's assumption, where it happens to hold. The honest `phi` there is what is still
    # held at the last station looked at, not zero: the store is returned by then only if the
    # solve says so. Computed rather than assumed, for the same reason `phi` itself is.
    phi = freeze.bond_energy_fraction if freeze is not None else stations[-1].bond_energy_fraction
    return TollPoint(
        closing_speed=closing_speed,
        slug_ratio=slug_ratio,
        rho_bag=rho_bag,
        dissipated=dissipated,
        temp_0=temp_0,
        rho_freeze=freeze.rho if freeze is not None else None,
        temp_freeze=freeze.temp if freeze is not None else None,
        bond_fraction=phi,
        available=available_energy(closing_speed, slug_ratio, phi),
        eta=eta_chem(closing_speed, slug_ratio, phi),
    )


def surface(
    speeds: tuple[float, ...] = SPEEDS,
    ratios: tuple[float, ...] = SLUG_RATIOS,
    rho_bag: float = plume.BAG_RHO,
) -> list[TollPoint]:
    """The full grid, speed-major so a fixed-`w` row reads contiguously."""
    return [toll_point(w, k, rho_bag) for w in speeds for k in ratios]


def fit_error(points: list[TollPoint]) -> tuple[float, TollPoint | None]:
    """Largest gap between the solved `eta` and the closed form evaluated at `phi = 1`.

    This is the deliverable's honesty check. The closed form is *exact* given `phi`, so the only
    question is whether `aim` can use it with `phi = 1` and skip this repository's solver. Where
    the plume fully atomizes it can; the answer here says where it cannot.
    """
    worst, at = 0.0, None
    for p in points:
        gap = abs(p.eta - eta_chem(p.closing_speed, p.slug_ratio, 1.0))
        if gap > worst:
            worst, at = gap, p
    return worst, at


def density_sensitivity(
    closing_speed: float, slug_ratio: float, densities: tuple[float, ...]
) -> list[TollPoint]:
    """`eta_chem` against bag density, for the convention noted in the module docstring.

    The surface assumes the bag is resized with `k`. If the *volume* is held fixed instead, the
    plume density scales with `(1+k)` and the freeze moves; this says how much that is worth.
    """
    return [toll_point(closing_speed, slug_ratio, rho) for rho in densities]


def main() -> None:
    """Publish the `eta_chem(w, k)` surface `aim` consumes (Q-R)."""
    points = surface()
    by_speed: dict[float, list[TollPoint]] = {}
    for p in points:
        by_speed.setdefault(p.closing_speed, []).append(p)

    print("python: eta_chem(w, k) -- frozen-chemistry ceiling on the nozzle's GROSS jet")
    print(f"  bag density {plume.BAG_RHO} kg/m^3; E_bond {E_BOND / 1e6:.2f} MJ/kg")
    print("  it multiplies sqrt(1+k), NOT the -1/+1 momentum term (routing doc Q-R)")
    header = "  " + " ".join(f"{f'k={k:g}':>7}" for k in SLUG_RATIOS)
    print(f"\n  {'w [km/s]':>9} |{header}")
    for w, rows in by_speed.items():
        print(f"  {w / 1e3:9.2f} | " + " ".join(f"{p.eta:7.3f}" for p in rows))

    print("\n  bond fraction held at the freeze (phi) -- 1.0 means the whole store strands")
    print(f"  {'w [km/s]':>9} |{header}")
    for w, rows in by_speed.items():
        print(f"  {w / 1e3:9.2f} | " + " ".join(f"{p.bond_fraction:7.4f}" for p in rows))

    print(
        f"\n  admissible k interval (T_0 >= {MIN_STAGNATION_TEMP:.0f} K)"
        f" -- replaces aim's flat k <= 80"
    )
    for w in SPEEDS:
        span = admissible_slug_ratios(w)
        shown = f"[{span[0]:.3f}, {span[1]:6.2f}]" if span else "empty -- no k lights this speed"
        print(f"    w = {w / 1e3:5.2f} km/s -> k in {shown}")

    print("\n  does it ignite? (dissipated >= 85.1 MJ/kg -- the paper's own bill, for comparison)")
    print(f"  {'w [km/s]':>9} |{header}")
    for w, rows in by_speed.items():
        print(f"  {w / 1e3:9.2f} | " + " ".join(f"{('y' if p.ignites else 'NO'):>7}" for p in rows))

    print("\n  no exhaust at all above this k (closed form, phi = 1 -- PESSIMISTIC, see Q-R):")
    for w in SPEEDS:
        print(f"    w = {w / 1e3:5.2f} km/s -> k_max = {zero_slug_ratio(w):6.2f}")

    worst, at = fit_error(points)
    print(f"\n  closed form at phi = 1 vs the solve: worst gap {worst:.4f}", end="")
    if at is not None:
        print(f" at w = {at.closing_speed / 1e3:.2f} km/s, k = {at.slug_ratio:g}")
    else:
        print()

    lines = [
        "closing_speed_km_s,slug_ratio,rho_bag_kg_m3,dissipated_MJ_kg,temp_0_K,"
        "rho_freeze_kg_m3,temp_freeze_K,bond_fraction,available_MJ_kg,eta_chem,conducts,ignites"
    ]
    for p in points:
        rf = f"{p.rho_freeze:.6e}" if p.rho_freeze is not None else ""
        tf = f"{p.temp_freeze:.1f}" if p.temp_freeze is not None else ""
        lines.append(
            f"{p.closing_speed / 1e3:g},{p.slug_ratio:g},{p.rho_bag:g},"
            f"{p.dissipated / 1e6:.4f},{p.temp_0:.1f},{rf},{tf},"
            f"{p.bond_fraction:.6f},{p.available / 1e6:.4f},{p.eta:.6f},"
            f"{int(p.conducts)},{int(p.ignites)}"
        )
    DEFAULT_TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_TABLE_PATH.write_text("\n".join(lines) + "\n")
    print(f"python: wrote {DEFAULT_TABLE_PATH}")


if __name__ == "__main__":
    main()
