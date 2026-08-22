"""Freeze-out: does the plume recombine faster than it expands? (Q-M)

`expansion.py` delivers the cooling history on both branches of ADR-0026's bracket and cannot
choose between them. At the nozzle exit the equilibrium branch reads 16 224 K and the frozen
branch 5 297 K -- a factor of 3, spanning the whole conductivity cliff, and the difference between
"no instability exposure anywhere" and "the cold leg is exposed". Two studies are each carrying
that bracket at full width.

**Choosing is a rate comparison, not a new simulation.** Recombination has a timescale and so does
the expansion, and the ratio decides:

    Da = tau_expansion / tau_recombination

`Da >> 1` means the chemistry keeps up and the equilibrium branch is right. `Da << 1` means the
composition is quenched at whatever it was and the frozen branch is right. `Da ~ 1` is the freeze
point -- the classical Bray criterion, and the same reasoning nozzle designers use for
frozen-flow losses.

**Two stores recombine, and they are not the same problem.** The ionisation store returns when
ions capture electrons; the (larger) dissociation store returns when atoms re-form molecules,
which needs a three-body collision. They have very different rates, so they can freeze at
different points and are carried separately here.

**Provenance.** Every coefficient below is a literature value, named in its docstring, with the
uncertainty stated. This module computes no rate coefficients of its own.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise

from puffsat import eos_water, expansion

#: Case-B hydrogenic radiative recombination at 10^4 K [m^3/s]. Osterbrock & Ferland, *Astrophysics
#: of Gaseous Nebulae and Active Galactic Nuclei* (2nd ed.), table 2.1: `2.59e-13 cm^3/s`.
ALPHA_B_10K = 2.59e-19
#: Reference temperature for `ALPHA_B_10K` [K].
T_REF_ALPHA_B = 1.0e4
#: Case-B temperature exponent. `alpha ~ T^-0.75` over ~5e3-2e4 K; a slow electron lingers near
#: the ion long enough to be caught, a fast one does not.
ALPHA_B_EXPONENT = -0.75


def radiative_recombination(t_e: float) -> float:
    """Case-B radiative recombination coefficient `alpha` [m^3/s] at electron temperature `t_e`.

    `X+ + e -> X + h nu`: the electron is captured and the binding energy leaves as a photon.

    **Case B, not case A.** Recombination straight to the ground state emits a photon energetic
    enough to re-ionise a neighbour immediately, so it returns nothing to the thermal pool and
    must not be counted. Case B excludes it. Using case A would overstate the rate by ~60% and
    push the freeze point in the optimistic direction.

    Valid to ~20% over 5e3-2e4 K, which covers the whole expansion. Above that the power law
    stiffens slightly; the plume is there only at the very start, where recombination is fast by
    orders of magnitude and the error cannot change a verdict.
    """
    return float(ALPHA_B_10K * (t_e / T_REF_ALPHA_B) ** ALPHA_B_EXPONENT)


#: Kelvin per electronvolt: `e/k_B`. The three-body coefficient is published in eV and everything
#: else here is in kelvin, so the conversion is named rather than inlined.
EV_IN_KELVIN = 1.602176634e-19 / 1.380649e-23

#: Three-body recombination coefficient at `T_e = 1 eV` [m^6/s]. Zel'dovich & Raizer give
#: `8.75e-27 T_e[eV]^-4.5 cm^6/s`; `1 cm^6 = 1e-12 m^6`.
K_THREE_BODY_1EV = 8.75e-39
#: Its temperature exponent. Steep, because capture needs the electron to linger inside the
#: Coulomb radius and both that radius and the dwell time shrink fast with temperature.
THREE_BODY_EXPONENT = -4.5


def three_body_coefficient(t_e: float) -> float:
    """Three-body recombination coefficient `K` [m^6/s] at electron temperature `t_e` [K].

    `X+ + e + e -> X + e`: one electron is captured and the other carries off the binding energy
    as kinetic energy. Because the energy stays in the gas rather than leaving as a photon, this
    channel returns the ionisation store *fully* -- which is exactly what the equilibrium branch
    of ADR-0026 assumes happens.
    """
    return float(K_THREE_BODY_1EV * (t_e / EV_IN_KELVIN) ** THREE_BODY_EXPONENT)


def three_body_recombination(t_e: float, n_e: float) -> float:
    """Effective two-body rate `alpha` [m^3/s] for the three-body channel at electron density `n_e`.

    The spectator electron is the third body, so the effective rate is `K n_e`: this channel turns
    off in a thin plasma and dominates in a dense one.
    """
    return three_body_coefficient(t_e) * n_e


def expansion_time(t_a: float, rho_a: float, t_b: float, rho_b: float) -> float:
    """Local expansion timescale [s]: the time for density to fall by a factor `e`.

    `tau_exp = rho / |d rho/dt|`, evaluated across a step as `(t_b - t_a) / ln(rho_a/rho_b)`.
    The logarithm is the point: density falling 4x over 2 ms leaves 1.44 ms of chemistry time,
    not 2 ms, because most of the fall happens in the last part of the interval.
    """
    if not (rho_a > 0.0 and rho_b > 0.0):
        raise ValueError("densities must be positive")
    if rho_b >= rho_a:
        raise ValueError("expansion requires rho_b < rho_a")
    return (t_b - t_a) / math.log(rho_a / rho_b)


#: Damkohler number above which the chemistry is taken to keep up, and below which it is quenched.
#: Order-unity by construction -- it *is* the definition of the freeze point (Bray) -- but the
#: transition is gradual, so a band is reported rather than a step.
DA_EQUILIBRIUM = 10.0
DA_FROZEN = 0.1


def verdict(damkohler: float) -> str:
    """Which branch of ADR-0026's bracket the gas follows at this `Da = tau_exp / tau_rec`.

    A band rather than a step at 1: real freeze-out takes about a decade in `Da` to complete, so
    calling `Da = 3` "equilibrium" would overstate what the criterion supports.
    """
    if damkohler >= DA_EQUILIBRIUM:
        return "equilibrium"
    if damkohler <= DA_FROZEN:
        return "frozen"
    return "freezing"


#: Low-pressure three-body rate for `H + OH + M -> H2O + M` at 1 K, in SI [m^6/s], from the fit
#: `6.1e-26 T^-2 cm^6 molecule^-2 s^-1` (Baulch et al., *Evaluated Kinetic Data for Combustion
#: Modelling*). This is the dominant water-reformation channel; `H + H + M -> H2 + M` is roughly
#: an order of magnitude slower and is not the bottleneck.
#:
#: **Uncertainty is the weak point of this module.** Evaluated three-body rates carry a factor
#: ~2-3 at combustion temperatures and more when extrapolated above them, and the third-body
#: efficiency of steam differs from the argon/nitrogen the fit is anchored on. The verdict below
#: is therefore reported with the margin in decades of rate, so a reader can see how much of that
#: uncertainty it survives.
K_ATOM_THREE_BODY = 6.1e-38
ATOM_THREE_BODY_EXPONENT = -2.0


def atom_three_body_coefficient(temp: float) -> float:
    """Three-body coefficient `K` [m^6/s] for water re-formation at gas temperature `temp` [K]."""
    return float(K_ATOM_THREE_BODY * temp**ATOM_THREE_BODY_EXPONENT)


def atom_recombination_time(temp: float, n_third_body: float, n_atom: float) -> float:
    """Time [s] to re-form molecules from atoms, returning the dissociation store.

    `tau = 1 / (K n_M n_atom)`: a three-body reaction is second order in the colliding partners
    and first order in the third body, so this channel is *quadratically* sensitive to the
    expansion -- it slows by 100x for every 10x the plume thins. That is why the dissociation
    store, not the ionisation store, is the one that can freeze.
    """
    rate = atom_three_body_coefficient(temp) * n_third_body * n_atom
    return math.inf if rate <= 0.0 else 1.0 / rate


@dataclass(frozen=True)
class FreezeStation:
    """The freeze verdict at one station of the cooling history."""

    time: float
    temp: float
    rho: float
    n_e: float
    n_atom: float
    n_third: float
    tau_expansion: float
    tau_ionisation: float
    tau_dissociation: float
    da_ionisation: float
    da_dissociation: float
    verdict_ionisation: str
    verdict_dissociation: str
    #: Ionisation energy still held here, as a fraction of the reservoir's store. A station with
    #: nothing left cannot freeze anything, however small its `Da`.
    ionisation_store_fraction: float
    #: Dissociated fraction. Equilibrium water stays fully dissociated across the whole nozzle at
    #: these densities, so this store has nothing to return here -- it returns in the fireball.
    dissociated_fraction: float
    #: How many decades the dissociation rate coefficient could be wrong before the verdict
    #: changes. The rate is the least certain input, so this is the number that says whether the
    #: answer is robust or merely arithmetic.
    margin_decades: float


def freeze_station(
    time: float,
    temp: float,
    rho: float,
    n_e: float,
    n_atom: float,
    n_third: float,
    tau_expansion: float,
    ionisation_store_fraction: float = 1.0,
    dissociated_fraction: float = 1.0,
) -> FreezeStation:
    """Race both recombination channels against the expansion at one state.

    The ionisation store returns by whichever ionic channel is faster -- three-body dominates in
    this density range by orders of magnitude, but both are summed so the model stays right if a
    thinner case is ever run. The dissociation store returns only through the three-body atomic
    channel.
    """
    alpha_ion = three_body_recombination(temp, n_e) + radiative_recombination(temp)
    tau_ion = math.inf if alpha_ion * n_e <= 0.0 else 1.0 / (alpha_ion * n_e)
    tau_atom = atom_recombination_time(temp, n_third, n_atom)

    da_ion = tau_expansion / tau_ion
    da_atom = tau_expansion / tau_atom
    return FreezeStation(
        time=time,
        temp=temp,
        rho=rho,
        n_e=n_e,
        n_atom=n_atom,
        n_third=n_third,
        tau_expansion=tau_expansion,
        tau_ionisation=tau_ion,
        tau_dissociation=tau_atom,
        da_ionisation=da_ion,
        da_dissociation=da_atom,
        verdict_ionisation=verdict(da_ion),
        verdict_dissociation=verdict(da_atom),
        ionisation_store_fraction=ionisation_store_fraction,
        dissociated_fraction=dissociated_fraction,
        margin_decades=math.log10(da_atom / DA_EQUILIBRIUM) if da_atom > 0.0 else -math.inf,
    )


def scan(
    closing_speed: float, temp_0: float, steps: int = 320, stride: int = 4
) -> list[FreezeStation]:
    """Run the equilibrium cooling history and test whether it satisfies its own assumption.

    This is the standard way to settle a freeze question and it is why no finite-rate solver is
    needed: take the equilibrium solution, ask at every station whether the chemistry was fast
    enough to have produced it, and see whether the answer is consistent. A history that says
    "the chemistry kept up" everywhere *is* the answer; one that says otherwise localises where
    a finite-rate calculation would have to start.
    """
    rows = expansion.cooling_history(
        expansion.BAG_RHO,
        temp_0,
        eos_water.pressure_energy,
        eos_water.sound_speed,
        expansion.AREA_RATIO_EXIT,
        expansion.FIELD_LENGTH,
        steps=steps,
    )[::stride]

    reservoir_store = _ionisation_store(expansion.BAG_RHO, temp_0)

    out: list[FreezeStation] = []
    for prev, row in pairwise(rows):
        comp = eos_water.composition(row.rho, row.temp)
        n_third = comp.n_h2o + comp.n_h + comp.n_o + comp.n_hp + sum(comp.n_o_ions) + comp.n_e
        out.append(
            freeze_station(
                time=row.time,
                temp=row.temp,
                rho=row.rho,
                n_e=comp.n_e,
                # H atoms are the partner that has to find an OH; the water-reformation channel is
                # limited by them, and `eos_water` carries no OH so this is the closest proxy the
                # species set allows. It is an *over*estimate of the rate where OH is scarce.
                n_atom=comp.n_h,
                n_third=n_third,
                tau_expansion=expansion_time(prev.time, prev.rho, row.time, row.rho),
                ionisation_store_fraction=_ionisation_store(row.rho, row.temp) / reservoir_store,
                dissociated_fraction=1.0 - comp.n_h2o / (row.rho / eos_water.M_H2O),
            )
        )
    return out


def main() -> None:
    """Report the freeze verdict on every leg of the burn envelope."""
    print(
        f"{'w':>7} {'T_exit':>8} {'binding Da_ion':>15} {'margin':>8} {'verdict':>13} "
        f"{'f_diss exit':>12} {'min Da_diss':>12}"
    )
    for speed, temp_0 in expansion.PLUME_STATES:
        st = scan(speed, temp_0)
        binding = binding_damkohler([(s.da_ionisation, s.ionisation_store_fraction) for s in st])
        print(
            f"{speed:7.2f} {st[-1].temp:8.0f} {binding:15.3g} "
            f"{math.log10(binding / DA_EQUILIBRIUM):8.2f} {verdict(binding):>13} "
            f"{st[-1].dissociated_fraction:12.4f} "
            f"{min(s.da_dissociation for s in st):12.2f}"
        )


def binding_damkohler(stations: Sequence[tuple[float, float]], threshold: float = 0.01) -> float:
    """The `Da` that actually governs the energy release, given `(Da, store fraction)` per station.

    **`min Da` over the history is the wrong statistic and it inverts this answer.** The
    Damkohler number keeps falling as the plume cools, because the recombination rate is
    proportional to the density of things left to recombine -- so the coldest, thinnest stations
    always look "frozen". But they hold nothing: freezing a store that has already been returned
    releases no energy and changes no temperature.

    So the binding value is the smallest `Da` among stations still holding more than `threshold`
    of the reservoir store. Raises if no station qualifies, rather than defaulting to a pass.
    """
    holding = [da for da, fraction in stations if fraction > threshold]
    if not holding:
        raise ValueError(f"no station holds more than {threshold:g} of the store")
    return min(holding)


def _ionisation_store(rho: float, temp: float) -> float:
    """Ionisation energy held per kg [J/kg] -- what recombination has left to give back."""
    comp = eos_water.composition(rho, temp)
    ladder = sum(n * eos_water.E_O_CUM[k] for k, n in enumerate(comp.n_o_ions))
    return (comp.n_hp * eos_water.IP_H + ladder) / rho


if __name__ == "__main__":
    main()
