"""N10 items 1-3: does the walled nozzle's expansion outrun its own recombination?

The chamber solve (`chamber.py`) settled how much of the atomisation store is *charged*: at the
flown 10 kK, 95-99% of it (W1). This module answers the other half -- how much of that charge is
handed back as directed kinetic energy before the gas leaves, and how much is carried out of the
nozzle still locked in broken bonds.

**It is a rate comparison on the equilibrium solution, not a new simulation.** Exactly as
`recombination.py` does for the magnetic nozzle: take the equilibrium expansion, and ask at every
station whether the chemistry was fast enough to have produced it.

    Da = tau_expansion / tau_recombination

`Da >> 1` and the equilibrium history is self-consistent; `Da << 1` and the composition is
quenched at whatever it was when the ratio fell through one. The station where it crosses is
Bray's sudden-freezing point, and the bond energy still held there is the store that is lost.

**Why this study needs its own module rather than reusing `recombination.py`.** That module's
dissociation channel is `H + OH + M -> H2O + M` on the water species set, with an OH-scarcity
bracket that has no analogue here. Methane's recoverable store returns through `H + H + M -> H2
+ M`, where both partners are the same species, so the scarcity bracket collapses and the answer
is sharper. The rate constants are in `rates.py`; none is computed here.

**The third body is the point.** A three-body rate is only as good as the efficiency assumed for
whatever is doing the stabilising, and at 10 kK the methane charge is 93-98% dissociated, so the
third body is *atomic hydrogen* -- an order of magnitude more efficient than argon and several
times more efficient than H2 (`rates.THIRD_BODY_EFFICIENCY`). `mixture_coefficient` weights the
species actually present rather than assuming one.

**What this module does not do.** It does not touch carbon: `rates.CARBON_CHANNEL_NOTE` records
that the gas-phase carbon association channels are essentially unmeasured, so the 43% of
methane's atomisation locked in carbon stays with N10 item 5 and its nucleation machinery. Nor
does it produce an ammonia specific impulse -- that needs an `eos_ammonia` this repository does
not have. What it does give for ammonia is the `N + N + M` rate comparison the ask's "nitrogen is
about eight times slower" rests on.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from puffsat import eos_methane, expansion
from puffsat.walled_nozzle import chamber, rates

#: Damkohler thresholds, the same band `recombination.py` uses. Order-unity by construction --
#: it *is* the definition of the freeze point -- but real freeze-out takes about a decade in `Da`
#: to complete, so a band is reported rather than a step at 1.
DA_EQUILIBRIUM = 10.0
DA_FROZEN = 0.1

#: Nozzle exit area ratio `A/A*` for the flown geometry: ADR-0016's 3 m bore against its 7 m^2
#: throat. `chamber.BORE_RADIUS` is the radius, so the bore is `pi r^2 = 28.3 m^2`.
BORE_AREA = math.pi * chamber.BORE_RADIUS**2

#: Nozzle length [m] used to turn the area profile into a clock. **This is the one assumed
#: shape in the calculation** and it is inherited from `expansion.cooling_history`, which opens
#: the area linearly over the length. Set to the flown chamber's own length so the nozzle is not
#: quietly given more room than the vehicle has; `history` takes it as an argument so the
#: sensitivity can be run, and `Da` scales linearly in it.
NOZZLE_LENGTH = chamber.CHAMBERS[0][1]

DEFAULT_OUTPUT = Path("data/results/walled_nozzle/freeze.csv")


def verdict(damkohler: float) -> str:
    """Which branch the gas follows at this `Da`."""
    if damkohler >= DA_EQUILIBRIUM:
        return "equilibrium"
    if damkohler <= DA_FROZEN:
        return "frozen"
    return "freezing"


def mixture_coefficient(temp: float, comp: eos_methane.Composition) -> float:
    """Efficiency-weighted three-body coefficient [m^6/s] for `H + H + M -> H2 + M`.

    A published three-body rate belongs to the third body it was measured with, and those differ
    by an order of magnitude (`rates.THIRD_BODY_EFFICIENCY`). Rather than pick one and hope, the
    stabilising flux is summed over the species actually present,

        k_eff = ( k_H n_H + k_H2 n_H2 + k_inert n_other ) / n_M

    so that `k_eff n_M` is the true stabilisation rate and the caller can keep using a single
    third-body density. Everything that is neither atomic H nor H2 -- the carbon bearers, the
    hydrocarbons -- is charged at the argon rate, which is the *slowest* of the three measured
    efficiencies and therefore the conservative assignment for an unmeasured collider.

    All three coefficients are taken from Jacobs, Giedt & Cohen 1967, i.e. from a **single**
    experiment. Ratios within one study cancel its systematics; ratios across studies do not.
    """
    n_h = comp.n_h
    n_h2 = comp.density_of(eos_methane.H2)
    n_third = comp.n_neutral_heavy
    n_other = max(0.0, n_third - n_h - n_h2)
    if n_third <= 0.0:
        return 0.0
    flux = (
        rates.H_ATOMIC_H(temp) * n_h
        + rates.H_H2_JACOBS(temp) * n_h2
        + rates.H_AR_JACOBS(temp) * n_other
    )
    return flux / n_third


def recombination_time(coefficient: float, n_partner: float, n_third: float) -> float:
    """Time [s] for the dissociation store to return: `tau = 1 / (k n_partner n_M)`.

    Second order in the recombining atoms and first order in the third body, so this channel is
    *quadratically* sensitive to the expansion -- it slows by 100x for every 10x the gas thins.
    That is exactly why a walled throat, which holds density up, is the design being argued for.
    """
    rate = coefficient * n_partner * n_third
    return math.inf if rate <= 0.0 else 1.0 / rate


@dataclass(frozen=True)
class FreezeStation:
    """One station along the nozzle, with the race called."""

    time: float
    x: float
    area_ratio: float
    rho: float
    temp: float
    pressure: float
    speed: float
    #: Atomic hydrogen density [m^-3] -- both the recombining partner and the best third body.
    n_h: float
    #: All neutral heavy particles [m^-3]: the third-body density.
    n_third: float
    #: Efficiency-weighted `k` [m^6/s] actually used here.
    coefficient: float
    tau_expansion: float
    tau_recombination: float
    damkohler: float
    verdict: str
    #: Bond energy still held at this station as a fraction of full atomisation. **This is the
    #: store that leaves as inert enthalpy if the flow freezes here.**
    store_held: float
    #: How many decades the rate coefficient could be wrong before this station's verdict flips
    #: from "equilibrium" to something else. The number the whole ask turns on.
    margin_decades: float

    @property
    def extrapolated(self) -> bool:
        """Whether this station sits outside the window the hot fits were established over."""
        return rates.H_ATOMIC_H.extrapolated(self.temp)


def station(
    time: float,
    x: float,
    area_ratio: float,
    rho: float,
    temp: float,
    pressure: float,
    speed: float,
    tau_expansion: float,
) -> FreezeStation:
    """Race recombination against the expansion at one nozzle state."""
    comp = eos_methane.composition(rho, temp)
    k = mixture_coefficient(temp, comp)
    tau_rec = recombination_time(k, comp.n_h, comp.n_neutral_heavy)
    da = tau_expansion / tau_rec if tau_rec > 0.0 else math.inf
    return FreezeStation(
        time=time,
        x=x,
        area_ratio=area_ratio,
        rho=rho,
        temp=temp,
        pressure=pressure,
        speed=speed,
        n_h=comp.n_h,
        n_third=comp.n_neutral_heavy,
        coefficient=k,
        tau_expansion=tau_expansion,
        tau_recombination=tau_rec,
        damkohler=da,
        verdict=verdict(da),
        store_held=eos_methane.bond_energy_held(rho, temp) / eos_methane.FULL_ATOMIZATION_ENERGY,
        margin_decades=math.log10(da / DA_EQUILIBRIUM) if da > 0.0 else -math.inf,
    )


def history(
    state: chamber.ChamberState,
    throat_area: float = chamber.THROAT_AREA_RANGE[1],
    length: float = NOZZLE_LENGTH,
    steps: int = 192,
) -> list[FreezeStation]:
    """Freeze stations from the sonic throat to the bore exit, for a solved chamber.

    The expansion itself is `expansion.cooling_history`, which is already field-free: it takes
    the EOS and the sound speed as arguments and knows nothing about a magnetic nozzle. Running
    it on `eos_methane` *is* the "no-field variant" -- no new integrator was required, only a
    different equation of state and a geometry that comes from ADR-0016 rather than the paper's
    field profile.
    """
    rows = expansion.cooling_history(
        state.rho,
        state.temp,
        eos_methane.pressure_energy,
        eos_methane.sound_speed,
        BORE_AREA / throat_area,
        length,
        steps=steps,
    )
    return _stations(rows)


def _stations(rows: list[expansion.CoolingRow]) -> list[FreezeStation]:
    """Race every interval of a cooling history, skipping degenerate ones.

    **Degenerate intervals must be dropped rather than scored.** `cooling_history` interpolates
    an exact station onto the throat and another onto the exit ratio, and either can land on top
    of a neighbouring sample. A zero-width interval has `dt = 0`, hence `tau_exp = 0` and a
    Damkohler of zero -- which would be read as "frozen at the throat", the densest and most
    strongly recombining point in the whole nozzle. It is an artefact of the grid, not physics,
    and it appears only at fine step counts, which is exactly when it is least expected.
    """
    out: list[FreezeStation] = []
    for prev, row in pairwise(rows):
        tau = _tau(prev, row)
        if not (tau > 0.0 and math.isfinite(tau)):
            continue
        out.append(
            station(
                time=row.time,
                x=row.x,
                area_ratio=row.area_ratio,
                rho=row.rho,
                temp=row.temp,
                pressure=row.pressure,
                speed=row.speed,
                tau_expansion=tau,
            )
        )
    return out


def _tau(prev: expansion.CoolingRow, row: expansion.CoolingRow) -> float:
    """Local expansion timescale [s]: the time for density to fall by a factor `e`.

    `tau_exp = rho / |d rho/dt|`, evaluated across a step as `(t_b - t_a) / ln(rho_a/rho_b)`.
    The logarithm is the point: density falling 4x over 2 ms leaves 1.44 ms of chemistry time,
    not 2 ms, because most of the fall happens late in the interval.
    """
    if not (prev.rho > 0.0 and row.rho > 0.0) or row.rho >= prev.rho:
        return math.inf
    if row.time <= prev.time:
        return math.inf
    return (row.time - prev.time) / math.log(prev.rho / row.rho)


def freeze_point(stations: list[FreezeStation]) -> FreezeStation | None:
    """First station whose `Da` has fallen through the equilibrium threshold, or `None`.

    `None` is a real answer and the one the design wants: the chemistry kept up all the way to
    the bore exit, so the store came back.
    """
    return next((s for s in stations if s.damkohler < DA_EQUILIBRIUM), None)


def returned_fraction(stations: list[FreezeStation], charged: float | None = None) -> float:
    """Fraction of the *charged* store handed back before the flow froze or reached the exit.

    **The reference is the chamber, not the throat.** `stations[0]` sits at the sonic throat,
    by which point the gas has already cooled from the reservoir and given some of the store
    back; measuring from there would silently discard the subsonic leg and understate the
    return. Pass `charged` -- the chamber's `store_charged` -- to reference the number the ask
    actually wants, which is the fraction of the *chamber's* charge that reappears as directed
    kinetic energy. Omitting it falls back to the throat and is only right for a nozzle-only
    question.

    Whatever is still held at the freeze station never comes back; if the flow never freezes,
    the exit plane is the cutoff and the remainder leaves as inert enthalpy anyway.
    """
    if not stations:
        return 0.0
    start = stations[0].store_held if charged is None else charged
    end = (freeze_point(stations) or stations[-1]).store_held
    return 0.0 if start <= 0.0 else max(0.0, (start - end) / start)


def freeze_search(
    state: chamber.ChamberState,
    throat_area: float = chamber.THROAT_AREA_RANGE[1],
    length: float = NOZZLE_LENGTH,
    area_ratio_end: float = 400.0,
    steps: int = 160,
    expansion_ratio: float = 4.0e5,
) -> list[FreezeStation]:
    """Keep expanding past the bore until the chemistry *does* freeze, to locate the station.

    **This is a diagnostic, not the flown geometry.** Inside ADR-0016's bore the gas never
    freezes (`history` shows `Da` in the hundreds at the exit plane), which answers "does it
    freeze" but not item 1's literal question, "where is the freezing station". Running the same
    isentrope out to a much larger area ratio finds it, and the bond energy still held there is
    the true ceiling on what any longer nozzle could ever recover.

    The clock is scaled with the area ratio so the assumed opening rate is unchanged: a nozzle
    that opens to 400x over a proportionally longer bore, not one that opens 100x faster.
    """
    scaled_length = length * (area_ratio_end - 1.0) / (BORE_AREA / throat_area - 1.0)
    rows = expansion.cooling_history(
        state.rho,
        state.temp,
        eos_methane.pressure_energy,
        eos_methane.sound_speed,
        area_ratio_end,
        scaled_length,
        steps=steps,
        expansion_ratio=expansion_ratio,
    )
    return _stations(rows)


def nitrogen_comparison(temp: float) -> dict[str, float]:
    """`N + N + M` against `H + H + M` at the same temperature, both on atomic third bodies.

    Tests the ask's stated "nitrogen is about eight times slower". Reported three ways because
    the nitrogen channel is genuinely uncertain: Byron's shock-tube fit is the only direct
    high-temperature measurement, and the modern ab initio master-equation result sits more than
    a decade below it. A claim that survives both ends is safe; one that needs Byron is not.
    """
    h_atomic = rates.H_ATOMIC_H(temp)
    n_atomic = rates.N_ATOMIC_N(temp)
    #: Log-interpolate the two published ab initio points onto `temp`.
    (t_a, k_a), (t_b, k_b) = rates.N_ABINITIO
    w = (math.log(temp) - math.log(t_a)) / (math.log(t_b) - math.log(t_a))
    n_abinitio = math.exp(math.log(k_a) + w * (math.log(k_b) - math.log(k_a)))
    return {
        "k_h_atomic": h_atomic,
        "k_n_byron": n_atomic,
        "k_n_abinitio": n_abinitio,
        "ratio_byron": h_atomic / n_atomic,
        "ratio_abinitio": h_atomic / n_abinitio,
        "nitrogen_spread_decades": math.log10(n_atomic / n_abinitio),
    }


CSV_HEADER = (
    "volume_m3,temp_k,throat_area_m2,time_ms,x_m,area_ratio,rho,temp_station_k,pressure_pa,"
    "speed_m_s,n_h,n_third,k_eff_m6_s,tau_exp_s,tau_rec_s,damkohler,verdict,store_held,"
    "margin_decades\n"
)


def write(
    rows: list[tuple[chamber.ChamberState, float, FreezeStation]],
    path: Path = DEFAULT_OUTPUT,
) -> None:
    """Write the freeze ledger. Committed, because the asks document is read in the paper repo."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(CSV_HEADER)
        for state, throat, s in rows:
            fh.write(
                f"{state.volume:g},{state.temp:g},{throat:g},{s.time * 1e3:.6f},{s.x:.4f},"
                f"{s.area_ratio:.5f},{s.rho:.6e},{s.temp:.2f},{s.pressure:.6e},"
                f"{s.speed:.2f},{s.n_h:.6e},{s.n_third:.6e},{s.coefficient:.6e},"
                f"{s.tau_expansion:.6e},{s.tau_recombination:.6e},{s.damkohler:.6e},"
                f"{s.verdict},{s.store_held:.6f},{s.margin_decades:.4f}\n"
            )


def main() -> None:
    """Report the freeze verdict for the candidate chambers and throat areas."""
    print("N10.1-3: recombination against expansion in the walled nozzle\n")
    print(f"rate: {rates.H_ATOMIC_H.source} (M = H), {rates.H_H2.source} (M = H2)")
    print(f"nozzle: bore {BORE_AREA:.1f} m^2, length {NOZZLE_LENGTH:.1f} m\n")

    ledger: list[tuple[chamber.ChamberState, float, FreezeStation]] = []
    print(
        f"{'V [m3]':>8} {'throat':>7} {'exit T':>8} {'min Da':>10} {'verdict':>12} "
        f"{'returned':>9} {'margin':>8}"
    )
    for volume, _length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, _length, chamber.FLOWN_TEMPERATURE)
        for throat in (2.0, 4.0, 7.0):
            st = history(state, throat_area=throat)
            ledger.extend((state, throat, s) for s in st)
            worst = min(st, key=lambda s: s.damkohler)
            print(
                f"{volume:8.0f} {throat:7.1f} {st[-1].temp:8.0f} {worst.damkohler:10.3g} "
                f"{worst.verdict:>12} {returned_fraction(st, state.store_charged):9.3f} "
                f"{worst.margin_decades:8.2f}"
            )

    print("\nWhere it would freeze if the nozzle kept opening (diagnostic, not the flown bore):")
    print(f"{'V [m3]':>8} {'A/A*':>9} {'T [K]':>8} {'held':>7} {'returned':>9} {'p [bar]':>9}")
    for volume, _length in chamber.CHAMBERS:
        state = chamber.solve_chamber(volume, _length, chamber.FLOWN_TEMPERATURE)
        deep = freeze_search(state)
        pt = freeze_point(deep)
        if pt is None:
            print(f"{volume:8.0f}   no freeze out to A/A* = {deep[-1].area_ratio:.0f}")
            continue
        print(
            f"{volume:8.0f} {pt.area_ratio:9.1f} {pt.temp:8.0f} {pt.store_held:7.3f} "
            f"{returned_fraction(deep, state.store_charged):9.3f} {pt.pressure / 1e5:9.3g}"
        )

    print("\nAmmonia rung -- N + N + M against H + H + M at 6000 K:")
    cmp = nitrogen_comparison(6000.0)
    print(f"  k(H+H+H)      = {cmp['k_h_atomic']:.3e} m^6/s")
    print(f"  k(N+N+N) Byron= {cmp['k_n_byron']:.3e} m^6/s   -> H is {cmp['ratio_byron']:.1f}x")
    print(f"  k(N+N+N) ab-in= {cmp['k_n_abinitio']:.3e} m^6/s -> H is {cmp['ratio_abinitio']:.1f}x")
    print(f"  nitrogen spread: {cmp['nitrogen_spread_decades']:.2f} decades")
    print(f"\nCarbon: {rates.CARBON_CHANNEL_NOTE}")

    write(ledger)
    print(f"\nwrote {DEFAULT_OUTPUT}")


if __name__ == "__main__":
    main()
