"""The chamber state: what the walled nozzle actually charges (N10 item 4b, N9 item 0).

**This is the ask's own hinge.** N10 item 4b says the paper-side arithmetic found the slug's
hydrogen "only 49 to 76% dissociated across the candidate chambers, not fully atomised, so the
chemical store is only partly charged", and that the correction moved methane from 1132 s to
1059-1095 s. It also says the finding "makes N9 and N10 one question", because a bigger chamber
charges more of the store while a smaller one is better thermally.

That is an equilibrium composition at `(rho, T)` and nothing else, so it is solved here rather
than estimated, on `eos_methane`.

**The slug ratio is a fixed point, not an input.** ADR-0016's own identity is that the exhaust
speed is `w / sqrt(1+k)`, which is energy conservation: all of the impactor's specific kinetic
energy `w^2/2` is shared over `(1+k)` units of mass, so

    `(1 + k) u(rho, T) = w^2 / 2`

with `u` the wall-cap energy density -- the internal energy a kilogram of working fluid holds at
the temperature the wall survives. But `k` sets the slug mass, which sets `rho` at a fixed
chamber volume, which sets `u`. So the three are solved together here. The check that this is
the right identity is that ask N9's own stated geometry -- 474 kg of methane against 25 kg of
impactor, `k = 18.96` -- falls out of it at 10 kK to better than a percent.

**What this module does not do.** It reports the equilibrium *charge*, not the *return*. Whether
that store comes back on the expansion is the Damkoehler question of N10 items 1-3, and whether
the carbon half of it condenses is item 5. A charge quoted without those is a ceiling.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_methane

# --- ADR-0016 geometry. None of this is the magnetic nozzle's (ADR-0050). ----------------------

#: Impactor mass [kg] and head-on closing speed [m/s] on the departure burn.
IMPACTOR_MASS = 25.0
CLOSING_SPEED = 75.0e3
#: Specific kinetic energy the impactor brings [J/kg of total charge], `w^2/2`.
SPECIFIC_PULSE_ENERGY = 0.5 * CLOSING_SPEED**2
#: Total pulse energy [J] -- the "70.3 GJ pulse" of ADR-0016.
PULSE_ENERGY = IMPACTOR_MASS * SPECIFIC_PULSE_ENERGY

#: Chamber bore **radius** [m]. ADR-0016 says "3 m bore" and means the radius: `pi r^2 L` closes
#: its own three volumes only at `r = 3`, and its 24.9 degree cone reaches "a 3 m wall". This is
#: the same 3.0 m as `expansion.CHAMBER_RADIUS`, but it is restated rather than imported, because
#: sharing a constant across two studies is how a geometry silently migrates (ADR-0050).
BORE_RADIUS = 3.0

#: The three candidate chambers of N9 item 0, as `(volume [m^3], length [m])`. Each closes on
#: `pi * BORE_RADIUS^2 * L` to better than half a percent.
CHAMBERS = ((200.0, 7.1), (400.0, 14.0), (673.0, 23.8))

#: Chamber temperatures the asks name. **10 kK is the flown point** (ADR-0016: Rubbia's radiative
#: ceiling, and going hotter buys 2.5% per 2000 K against much larger uncertainties).
TEMPERATURES = (8000.0, 10000.0, 12000.0)
FLOWN_TEMPERATURE = 10000.0

#: The charge N9 states, for the cross-check: 474 kg of methane on 25 kg of impactor.
ASK_SLUG_MASS = 474.0
ASK_SLUG_RATIO = ASK_SLUG_MASS / IMPACTOR_MASS  # 18.96

#: Pre-charge state ADR-0016 quotes: methane vapour at 1.4 bar and 111 K.
PRECHARGE_TEMPERATURE = 111.0

DEFAULT_OUTPUT = Path("data/results/walled_nozzle/chamber.csv")


@dataclass(frozen=True)
class ChamberState:
    """The solved chamber at one `(volume, temperature)`."""

    volume: float
    length: float
    temp: float
    #: Slug ratio `k = m_slug / m_impactor`, the fixed point of `(1+k) u = w^2/2`.
    slug_ratio: float
    slug_mass: float
    rho: float
    pressure: float
    #: Wall-cap energy density `u` [J/kg]: the internal energy a kilogram holds at `temp`.
    energy: float
    sound_speed: float
    #: Fraction of hydrogen *nuclei* that are free (H or H+) -- the number item 4b turns on.
    hydrogen_dissociated: float
    #: Fraction of carbon nuclei that are free atoms or atomic ions.
    carbon_free: float
    #: Fraction of the full atomisation store actually charged at this state. **This, not the
    #: dissociated fraction, is what the specific impulse depends on**: a nucleus in `C2H2` has
    #: paid back a different amount than one in `H2`.
    store_charged: float
    #: Sensible (thermal) share of `u`. ADR-0016's temperature-scaling argument rests on it
    #: being small -- only this part responds to `T`.
    sensible_fraction: float
    #: Free-electron fraction of all particles. Reported because N11's opacity and the wall's
    #: radiative load both scale on it, and because it is what makes carbon a high-Z worry.
    ionisation_fraction: float
    #: **The model's own worst exposure, stated rather than buried.** `C3` is treated as a
    #: harmonic rotor-oscillator, and its `nu2` bend at 63 cm^-1 is quasilinear -- so floppy that
    #: the harmonic partition function runs to ~90 per mode at 8 kK and certainly overstates it.
    #: This field is the *exact upper bound* on what that costs: the share of the atomisation
    #: store `C3` is withholding, all of which would be charged if `C3` vanished entirely. At the
    #: flown 10 kK it is about a percent; at 8 kK it is not.
    c3_store_exposure: float

    @property
    def precharge_pressure(self) -> float:
        """Pressure [Pa] of the cold methane charge before the impact, as an ideal gas at 111 K.

        Reported because ADR-0016 quotes 1.4 bar as *the* pre-charge, and that value belongs to
        the 200 m^3 chamber alone -- the same mass in 673 m^3 sits near 0.4 bar.
        """
        n = self.slug_mass / self.volume / eos_methane.M_CH4
        return n * eos_methane.K_B * PRECHARGE_TEMPERATURE


def solve_chamber(volume: float, length: float, temp: float, tol: float = 1e-10) -> ChamberState:
    """Solve `(k, rho, u)` together at one chamber volume and temperature.

    Fixed-point iteration on `k`, damped by averaging: `u` falls as the chamber is loaded more
    heavily (higher `rho` suppresses dissociation), so the map is a contraction and converges in
    a handful of steps from any sane start.
    """
    k = ASK_SLUG_RATIO
    comp = None
    for _ in range(200):
        rho = (1.0 + k) * IMPACTOR_MASS / volume
        _, u = eos_methane.pressure_energy(rho, temp)
        k_new = SPECIFIC_PULSE_ENERGY / u - 1.0
        if abs(k_new - k) < tol * k:
            k = k_new
            break
        k = 0.5 * (k + k_new)
    else:  # pragma: no cover - the map is a contraction; this guards a future EOS change
        raise RuntimeError(f"chamber fixed point did not converge at V={volume} m^3, T={temp} K")

    rho = (1.0 + k) * IMPACTOR_MASS / volume
    p, u = eos_methane.pressure_energy(rho, temp)
    comp = eos_methane.composition(rho, temp)
    held = eos_methane.bond_energy_held(rho, temp)
    # Sensible = everything that is not chemical. The ionisation store is charged too, but it is
    # counted with the chemical share because it returns by the same race (N10 items 1-3).
    e_ion = (
        comp.n_hp * eos_methane.IP_H
        + sum(n * eos_methane.E_C_CUM[j] for j, n in enumerate(comp.n_c_ions))
    ) / rho
    return ChamberState(
        volume=volume,
        length=length,
        temp=temp,
        slug_ratio=k,
        slug_mass=k * IMPACTOR_MASS,
        rho=rho,
        pressure=p,
        energy=u,
        sound_speed=eos_methane.sound_speed(rho, temp),
        hydrogen_dissociated=comp.hydrogen_dissociated,
        carbon_free=comp.carbon_free,
        store_charged=held / eos_methane.FULL_ATOMIZATION_ENERGY,
        sensible_fraction=(u - held - e_ion) / u,
        ionisation_fraction=comp.n_e / comp.n_total,
        c3_store_exposure=comp.density_of(eos_methane.C3)
        * eos_methane.C3.d_at
        / (rho / eos_methane.M_CH4 * eos_methane.CH4.d_at),
    )


def isp_scaling(state: ChamberState, reference: float) -> float:
    """Specific impulse relative to a chamber of slug ratio `reference`, on ADR-0016's identity.

    **This is a ratio and not a specific impulse, deliberately.** The companion repository's
    effective-Isp column carries launch-ledger conventions this repository does not own -- vessel
    mass, the drift term of `eq:reflection_baseline`, `eta_geom` -- and reproducing a number
    without owning its normalisation would be arithmetic dressed as a result.

    What *is* owned here is the identity ADR-0016 states: exhaust speed `w/sqrt(1+k)`, so the
    impulse per pulse is `m_p w sqrt(1+k)` and the impulse per kilogram of launched slug goes as
    `sqrt(1+k)/k`. Every normalisation the paper applies is a factor on that, so the ratio
    survives them and can be applied directly to whatever the paper's own column says.
    """

    def per_kg(k: float) -> float:
        return math.sqrt(1.0 + k) / k

    return per_kg(state.slug_ratio) / per_kg(reference)


@dataclass(frozen=True)
class AcousticCheck:
    """N9 item 0: does a *sealed* vessel escape the `sec:needle_through_fog` coupling problem?

    The worry that section carries is that a magnetic nozzle's bag is a free-standing cloud, so
    propellant the spreading cone misses is simply left behind and `k` is set by what the cone
    sweeps rather than by what is loaded. In a closed chamber the unswept mass is still in the
    chamber. The question is whether it has *time* to be reached, and that is two separate
    comparisons rather than one -- which is the finding.
    """

    length: float
    sound_speed: float
    #: Time for a sound wave to cross the column once [s].
    acoustic_time: float
    #: Time for the impactor to cross the column at closing speed [s].
    crossing_time: float
    #: Blowdown duration [s] -- how long the chamber takes to empty through the throat.
    blowdown_time: float
    #: Acoustic crossings completed during the blowdown. This is the number that decides item 0.
    equilibrations: float
    #: Acoustic crossings completed during the impactor's own transit. **Much less than one**,
    #: which is why the wall's load case is still a cone problem.
    crossings_during_impact: float

    @property
    def sealed_vessel_sets_k(self) -> bool:
        """Whether `k` is set by what is loaded rather than by what the cone sweeps.

        A column that equilibrates many times over before it empties has no memory of which
        parcels the cone touched first, so every kilogram loaded is every kilogram expelled hot.
        """
        return self.equilibrations >= 10.0


#: Throat area [m^2] of ADR-0016's chamber. The ask gives "a 2 to 7 m^2 throat" and N10 item 4
#: asks for the sensitivity, so both ends are carried.
THROAT_AREA_RANGE = (2.0, 7.0)


def acoustic_check(state: ChamberState, throat_area: float = THROAT_AREA_RANGE[1]) -> AcousticCheck:
    """Run item 0's timescale comparison at a solved chamber state.

    Blowdown is timed from the choked throat: a sonic throat passes `rho* a* A*` and the
    stagnation state falls behind it, so `t ~ M / (rho a A*)` with an order-unity coefficient. The
    coefficient used is `(2/(gamma+1))^((gamma+1)/(2(gamma-1)))`, the standard choked-flow factor,
    evaluated at `gamma = 1.2` -- appropriate for a dissociating gas, where the compression work
    partly goes into chemistry rather than into temperature. It is 0.65, so the answer is not
    sensitive to it: item 0 turns on a factor of fifty, not a factor of 1.5.
    """
    gamma = 1.2
    choked = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    mass = (1.0 + state.slug_ratio) * IMPACTOR_MASS
    mdot = choked * state.rho * state.sound_speed * throat_area
    acoustic = state.length / state.sound_speed
    blowdown = mass / mdot
    return AcousticCheck(
        length=state.length,
        sound_speed=state.sound_speed,
        acoustic_time=acoustic,
        crossing_time=state.length / CLOSING_SPEED,
        blowdown_time=blowdown,
        equilibrations=blowdown / acoustic,
        crossings_during_impact=(state.length / CLOSING_SPEED) / acoustic,
    )


CSV_HEADER = (
    "volume_m3",
    "length_m",
    "temp_k",
    "slug_ratio",
    "slug_mass_kg",
    "rho_kg_m3",
    "pressure_bar",
    "energy_mj_kg",
    "sound_speed_m_s",
    "hydrogen_dissociated",
    "carbon_free",
    "store_charged",
    "sensible_fraction",
    "ionisation_fraction",
    "c3_store_exposure",
    "precharge_bar",
    "isp_vs_ask_geometry",
    "acoustic_ms",
    "blowdown_ms",
    "equilibrations",
    "crossings_during_impact",
)


def scan() -> list[tuple[ChamberState, AcousticCheck]]:
    """Every chamber at every named temperature, with item 0's timescales alongside."""
    out = []
    for volume, length in CHAMBERS:
        for temp in TEMPERATURES:
            state = solve_chamber(volume, length, temp)
            out.append((state, acoustic_check(state)))
    return out


def write(rows: list[tuple[ChamberState, AcousticCheck]], path: Path = DEFAULT_OUTPUT) -> None:
    """Write the scan to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for state, check in rows:
            writer.writerow(
                [
                    f"{state.volume:.1f}",
                    f"{state.length:.2f}",
                    f"{state.temp:.0f}",
                    f"{state.slug_ratio:.3f}",
                    f"{state.slug_mass:.1f}",
                    f"{state.rho:.4f}",
                    f"{state.pressure / 1e5:.1f}",
                    f"{state.energy / 1e6:.2f}",
                    f"{state.sound_speed:.0f}",
                    f"{state.hydrogen_dissociated:.4f}",
                    f"{state.carbon_free:.4f}",
                    f"{state.store_charged:.4f}",
                    f"{state.sensible_fraction:.4f}",
                    f"{state.ionisation_fraction:.5f}",
                    f"{state.c3_store_exposure:.4f}",
                    f"{state.precharge_pressure / 1e5:.3f}",
                    f"{isp_scaling(state, ASK_SLUG_RATIO):.4f}",
                    f"{check.acoustic_time * 1e3:.3f}",
                    f"{check.blowdown_time * 1e3:.2f}",
                    f"{check.equilibrations:.1f}",
                    f"{check.crossings_during_impact:.4f}",
                ]
            )


def main() -> None:
    """Report the chamber ledger, and the two verdicts it settles."""
    rows = scan()
    write(rows)

    print("N10 item 4b -- the chamber charge, solved rather than assumed")
    print(
        f"{'V':>6} {'T':>7} {'k':>7} {'slug':>7} {'rho':>7} {'p':>8} {'u':>8} "
        f"{'diss H':>7} {'free C':>7} {'charge':>7} {'sens':>6} {'Isp/ask':>8}"
    )
    for state, _ in rows:
        print(
            f"{state.volume:6.0f} {state.temp:7.0f} {state.slug_ratio:7.2f} "
            f"{state.slug_mass:7.1f} {state.rho:7.3f} {state.pressure / 1e5:8.1f} "
            f"{state.energy / 1e6:8.1f} {state.hydrogen_dissociated:7.3f} "
            f"{state.carbon_free:7.3f} {state.store_charged:7.3f} "
            f"{state.sensible_fraction:6.3f} {state.c3_store_exposure:7.4f} "
            f"{isp_scaling(state, ASK_SLUG_RATIO):8.4f}"
        )

    print()
    print("N9 item 0 -- sealed-vessel equilibration, at the flown 10 kK")
    print(
        f"{'V':>6} {'L':>6} {'c_s':>7} {'t_ac':>8} {'t_blow':>8} {'n_eq':>7} "
        f"{'t_cross':>8} {'n_cross':>8} {'k set by load?':>16}"
    )
    for state, _check in rows:
        if state.temp != FLOWN_TEMPERATURE:
            continue
        for area in THROAT_AREA_RANGE:
            c = acoustic_check(state, area)
            print(
                f"{state.volume:6.0f} {state.length:6.1f} {c.sound_speed:7.0f} "
                f"{c.acoustic_time * 1e3:8.3f} {c.blowdown_time * 1e3:8.2f} "
                f"{c.equilibrations:7.1f} {c.crossing_time * 1e3:8.4f} "
                f"{c.crossings_during_impact:8.4f} "
                f"{('yes' if c.sealed_vessel_sets_k else 'NO') + f' (A*={area:.0f} m2)':>16}"
            )

    print()
    print(f"wrote {DEFAULT_OUTPUT}")


if __name__ == "__main__":
    main()
