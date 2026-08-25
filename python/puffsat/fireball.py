"""Does recombination freeze below 0.01 kg/m^3? The fireball past the nozzle lip.

The routing document carries this as an open item against `sec:watering_it_down`, and Q-M
sharpened it: the *dissociation* store never returns inside the nozzle -- equilibrium water is
fully dissociated at every station there -- so the marginal Damkohler number that store carries
has to bind somewhere downstream or not at all. Downstream is here.

**Two things happen at the lip, and the first is the trigger.** Inside the field the paper's bore
stretches the area ratio 1 -> 4 over 23.8 m, which is ~7.9 m of travel per unit `A/A*`. A 45-degree
free jet covers the same unit in 0.75 m of radius growth. So the *local* expansion time steps down
by ~8x at the lip -- 2.66 ms just inside to 0.33 ms just outside -- and `Da` falls from ~35 to
~1.8 in one step. **The magnetic nozzle is holding the plume in a slow expansion; the freeze
begins where the field lets go.**

Then the rates finish it. Three-body atomic recombination is a density-*cubed* process, so every
decade of density costs two decades of `tau_rec` while `tau_exp` grows only slowly. `Da` collapses
monotonically from there and never recovers.

**The local rate is the right measure and the aggregate one hides this.** Averaged over the whole
tail the fireball is *not* faster per decade than the nozzle (3.8 ms/decade against 3.0), because
far out `tau_exp` has grown large again. The step is local to the lip, which is exactly where it
matters.

    inside  (A/A* <= 4):  the nozzle clock, a linear area opening over `FIELD_LENGTH`
    outside (A/A* >  4):  dt = dR / (u tan theta),  R = plume_radius(A/A*)

**`theta` is the only free parameter and it is stated rather than fitted.** It scales the clock
linearly, so `Da ~ 1/tan theta`: a wider jet expands faster and freezes at higher density. The
verdict is reported across a 15-60 degree bracket for that reason -- the freeze *density* is
conditional on the angle, the freeze itself is not.

**The reformation partner is bracketed, and the bracket turns out not to matter (2026-08-25).**
`H + OH + M -> H2O + M` was rated with `n_H`, because the species set had no OH and the code said
plainly that this *over*estimates the rate. `eos_water` now carries OH, and the proxy was wrong by
30x on the cold leg and 50000x on the hot ones. It changes almost nothing: the stranded energy is
whatever the store still holds when the clock runs out, and on the hot legs it holds **100% at both
edges**. Only the cold anchor has a real bracket -- 0.83-0.94 of the store on the optimistic edge
against 0.999 on the conservative one -- and that is narrower than the jet-divergence bracket the
answer already carried. What settles the middle is still a finite-rate network with OH, H2 and O2
as *evolved* species; what this shows is that the answer is not waiting on it.

**And the store is now measured as held bond energy, not as `1 - n_H2O/n_f`.** With OH, H2 and O2
in the set the old measure overcharges -- it calls a gas fully dissociated when half of it is OH,
which is paying for half a bond. See `eos_water.bond_energy_held`.

**This is a self-consistency test, like Q-M's.** The equilibrium history is run and asked at every
station whether the chemistry could have produced it. Where the answer is no, equilibrium is not
the real history below that point -- the composition holds near where it froze and the temperature
stops being buffered. So the stations past the crossing say *where* a finite-rate calculation would
have to start, not what it would find.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache
from itertools import pairwise
from pathlib import Path

from puffsat import eos_water, expansion, recombination

#: Jet divergence half-angle [deg] past the nozzle lip. 45 degrees makes `tan theta = 1`, i.e. the
#: plume's axial advance matches its radial growth -- the neutral choice, bracketed in `main`.
DIVERGENCE_HALF_ANGLE_DEG = 45.0

#: The density the paper names the question at [kg/m^3].
FREEZE_DENSITY = 0.01

#: How far past the lip to integrate. 8192x the bag density reaches ~4e-5 kg/m^3, three decades
#: below the crossing, which is far enough that the frozen tail is unambiguous.
EXPANSION_RATIO = 8192.0

DEFAULT_SCAN_PATH = Path("data/results/fireball_freeze.csv")

#: Slug ratio, for charging the stranded energy against the budget that produced the plume.
SLUG_RATIO = 8.5


@dataclass(frozen=True)
class FireballRow:
    """One station of the combined nozzle-plus-fireball history."""

    time: float
    radius: float
    area_ratio: float
    rho: float
    temp: float
    speed: float
    #: Bond energy still held here as a fraction of full atomization -- the share of the store
    #: that would be stranded if the chemistry froze at this station.
    bond_energy_fraction: float


@cache
def _cached_history(
    temp_0: float, half_angle_deg: float, rho_0: float, steps: int, expansion_ratio: float
) -> tuple[FireballRow, ...]:
    """`history` memoised. Every angle in the bracket re-walks the same isentrope, and the
    isentrope is the expensive half -- the freeze race itself is arithmetic."""
    return tuple(_history(temp_0, half_angle_deg, rho_0, steps, expansion_ratio))


def history(
    temp_0: float,
    half_angle_deg: float = DIVERGENCE_HALF_ANGLE_DEG,
    rho_0: float = expansion.BAG_RHO,
    steps: int = 320,
    expansion_ratio: float = EXPANSION_RATIO,
) -> list[FireballRow]:
    """The equilibrium history from the sonic throat out into the free jet."""
    return list(_cached_history(temp_0, half_angle_deg, rho_0, steps, expansion_ratio))


@cache
def _isentrope(
    temp_0: float, rho_0: float, steps: int, expansion_ratio: float
) -> tuple[tuple[FireballRow, ...], tuple[expansion.NozzlePoint, ...]]:
    """The angle-independent half of the history: the isentrope and the in-nozzle stations.

    Split out and memoized because `theta` only enters the *free-jet clock*. `main` sweeps three
    divergence angles per leg, and without this split each one re-solves the same nozzle from
    scratch -- the EOS Newton solve is the hot spot in this package, so that was three quarters of
    the run time doing nothing.
    """
    points = expansion.nozzle_history(
        rho_0,
        temp_0,
        rho_0 / expansion_ratio,
        eos_water.pressure_energy,
        eos_water.sound_speed,
        steps,
    )
    throat = expansion._at_throat(points)
    supersonic = (throat, *(pt for pt in points if pt.mach > throat.mach))

    inside = expansion.cooling_history(
        rho_0,
        temp_0,
        eos_water.pressure_energy,
        eos_water.sound_speed,
        expansion.AREA_RATIO_EXIT,
        expansion.FIELD_LENGTH,
        steps=steps,
    )
    rows = tuple(
        FireballRow(
            time=r.time,
            radius=expansion.plume_radius(r.area_ratio),
            area_ratio=r.area_ratio,
            rho=r.rho,
            temp=r.temp,
            speed=r.speed,
            bond_energy_fraction=_bond_fraction(r.rho, r.temp),
        )
        for r in inside
    )
    return rows, supersonic


def _history(
    temp_0: float,
    half_angle_deg: float,
    rho_0: float,
    steps: int,
    expansion_ratio: float,
) -> list[FireballRow]:
    """The equilibrium history from the sonic throat out into the free jet.

    The nozzle leg reuses `expansion.cooling_history` unchanged rather than reproducing its clock,
    so the two cannot drift; the free leg continues from its last station on the same isentrope.
    """
    if not 0.0 < half_angle_deg < 90.0:
        raise ValueError(f"half_angle_deg must be in (0, 90), got {half_angle_deg}")
    tan_theta = math.tan(math.radians(half_angle_deg))

    nozzle_rows, supersonic = _isentrope(temp_0, rho_0, steps, expansion_ratio)
    rows = list(nozzle_rows)

    # The free jet: continue from the lip on the same isentrope, with the jet clock.
    for pt in (p for p in supersonic if p.area_ratio > expansion.AREA_RATIO_EXIT):
        radius = expansion.plume_radius(pt.area_ratio)
        prev = rows[-1]
        step = (radius - prev.radius) * 0.5 * (1.0 / pt.speed + 1.0 / prev.speed) / tan_theta
        rows.append(
            FireballRow(
                time=prev.time + step,
                radius=radius,
                area_ratio=pt.area_ratio,
                rho=pt.rho,
                temp=pt.temp,
                speed=pt.speed,
                bond_energy_fraction=_bond_fraction(pt.rho, pt.temp),
            )
        )
    return rows


def _bond_fraction(rho: float, temp: float) -> float:
    """Share of the full atomization energy still held as broken bonds.

    **Replaces `1 - n_H2O/n_f` (2026-08-25).** Once OH, H2 and O2 are in the species set that
    measure overcharges: it calls a gas "dissociated" when it is half OH, which is paying for half
    a bond. At 0.32 kg/m^3 and 3000 K the two differ by a factor of 3.
    """
    return eos_water.bond_energy_held(rho, temp) / eos_water.FULL_ATOMIZATION_ENERGY


@dataclass(frozen=True)
class FreezeRow:
    """A station of the history together with the freeze race at it.

    The two are carried side by side because "is this inside the field?" is a question about the
    *geometry* -- the area ratio -- and `recombination.FreezeStation` has no business knowing about
    nozzles. Splitting the legs on a density threshold instead looks equivalent and is not: the lip
    density differs from leg to leg, so a fixed threshold lands on the wrong side on some of them.
    """

    #: The station the race is evaluated at.
    row: FireballRow
    #: The station before it. `tau_expansion` is a property of the *interval* between the two,
    #: which is why it is carried: an interval spanning the lip is on neither clock.
    prev: FireballRow
    station: recombination.FreezeStation

    @property
    def inside(self) -> bool:
        """True while the magnetic nozzle is still holding the plume."""
        return self.row.area_ratio <= expansion.AREA_RATIO_EXIT

    @property
    def straddles_lip(self) -> bool:
        """True when this station's expansion time is measured across the field boundary.

        Such an interval carries part of the nozzle's slow clock and part of the jet's fast one,
        so its `tau_expansion` is a blend of two regimes and belongs to neither. There is exactly
        one of these per history, and it is excluded when the two legs are compared.
        """
        exit_ratio = expansion.AREA_RATIO_EXIT
        return (self.prev.area_ratio <= exit_ratio) != (self.row.area_ratio <= exit_ratio)

    @property
    def rho(self) -> float:
        return self.station.rho

    @property
    def temp(self) -> float:
        return self.station.temp

    @property
    def tau_expansion(self) -> float:
        return self.station.tau_expansion

    @property
    def da_dissociation_h_limited(self) -> float:
        """The optimistic edge: an OH is always available. What shipped before 2026-08-25."""
        return self.station.da_dissociation_h_limited

    @property
    def da_dissociation_oh_limited(self) -> float:
        """The conservative edge: OH never exceeds its equilibrium density."""
        return self.station.da_dissociation_oh_limited

    @property
    def bond_energy_fraction(self) -> float:
        return self.station.bond_energy_fraction


def scan(
    temp_0: float,
    half_angle_deg: float = DIVERGENCE_HALF_ANGLE_DEG,
    stride: int = 4,
    steps: int = 320,
    rho_0: float = expansion.BAG_RHO,
    expansion_ratio: float = EXPANSION_RATIO,
) -> list[FreezeRow]:
    """Race atomic recombination against the free expansion at every station.

    `rho_0` and `expansion_ratio` are parameters rather than module constants because `toll.py`
    sweeps `(w, k)` and needs both: the reservoir density is a design variable, and integrating
    three decades past a freeze that happens at the lip is most of the cost of a grid.
    """
    rows = list(_cached_history(temp_0, half_angle_deg, rho_0, steps, expansion_ratio))[::stride]
    out: list[FreezeRow] = []
    for prev, row in pairwise(rows):
        comp = eos_water.composition(row.rho, row.temp)
        station = recombination.freeze_station(
            time=row.time,
            temp=row.temp,
            rho=row.rho,
            n_e=comp.n_e,
            n_h=comp.n_h,
            n_oh=comp.n_oh,
            n_third=comp.n_total,
            tau_expansion=recombination.expansion_time(prev.time, prev.rho, row.time, row.rho),
            bond_energy_fraction=row.bond_energy_fraction,
        )
        out.append(FreezeRow(row=row, prev=prev, station=station))
    return out


def freeze_state(stations: list[FreezeRow], oh_limited: bool = False) -> FreezeRow | None:
    """The first station where the expansion outruns water reformation (`Da < 1`).

    `Da = 1` rather than `recombination.DA_FROZEN`: the deliverable is *where the equilibrium
    assumption stops holding*, which is where the two clocks are equal. `DA_FROZEN = 0.1` is a
    decade further out and describes a store that is already stranded rather than one that is
    beginning to strand.

    `oh_limited` selects the edge of the bracket (see `recombination.freeze_station`). The default
    is the optimistic one, so the shipped headline stays the one that is hardest to argue with.
    """

    def damkohler(s: FreezeRow) -> float:
        return s.da_dissociation_oh_limited if oh_limited else s.da_dissociation_h_limited

    return next((s for s in stations if damkohler(s) < 1.0), None)


def stranded_energy(station: FreezeRow) -> float:
    """Bond energy [J/kg] still held when the store freezes.

    This never becomes directed kinetic energy -- it leaves as inert chemical enthalpy. It is the
    reason the question is worth asking: it is charged against the same dissipated budget an Isp
    claim is built on.
    """
    return station.bond_energy_fraction * eos_water.FULL_ATOMIZATION_ENERGY


def main() -> None:
    """The verdict, across the burn envelope and across the divergence bracket."""
    angles = (15.0, DIVERGENCE_HALF_ANGLE_DEG, 60.0)
    lines = [
        "closing_speed_km_s,half_angle_deg,partner,freeze_rho_kg_m3,freeze_temp_K,freeze_time_s,"
        "bond_energy_fraction,tau_expansion_s,tau_dissociation_s,stranded_MJ_kg,"
        "stranded_fraction_of_budget"
    ]
    print(f"python: fireball freeze-out, jet divergence {angles[0]:g}-{angles[-1]:g} deg")
    print(
        "  partner H = an OH is always available (optimistic); OH = equilibrium OH (conservative)"
    )
    print(
        f"  {'w [km/s]':>9} {'theta':>6} {'partner':>8} {'rho_freeze':>11} {'T [K]':>7} "
        f"{'t [ms]':>8} {'bond held':>10} {'stranded':>10} {'of budget':>10}"
    )
    for speed, temp_0 in expansion.PLUME_STATES:
        budget = 0.5 * SLUG_RATIO * (speed * 1e3) ** 2 / (1.0 + SLUG_RATIO) ** 2
        for angle in angles:
            stations = scan(temp_0, angle)
            for label, oh_limited in (("H", False), ("OH", True)):
                station = freeze_state(stations, oh_limited=oh_limited)
                if station is None:
                    print(
                        f"  {speed:9.2f} {angle:6.0f} {label:>8}   equilibrium holds to the end "
                        f"of the range"
                    )
                    continue
                stranded = stranded_energy(station)
                tau = (
                    station.station.tau_dissociation_oh_limited
                    if oh_limited
                    else station.station.tau_dissociation_h_limited
                )
                print(
                    f"  {speed:9.2f} {angle:6.0f} {label:>8} {station.rho:11.4e} "
                    f"{station.temp:7.0f} {station.station.time * 1e3:8.3f} "
                    f"{station.bond_energy_fraction:10.4f} "
                    f"{stranded / 1e6:9.1f}M {stranded / budget:10.1%}"
                )
                lines.append(
                    f"{speed:g},{angle:g},{label},{station.rho:.6e},{station.temp:.1f},"
                    f"{station.station.time:.6e},"
                    f"{station.bond_energy_fraction:.6f},{station.tau_expansion:.6e},"
                    f"{tau:.6e},{stranded / 1e6:.3f},"
                    f"{stranded / budget:.6f}"
                )
    DEFAULT_SCAN_PATH.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_SCAN_PATH.write_text("\n".join(lines) + "\n")
    print(f"python: wrote {DEFAULT_SCAN_PATH}")
    print(
        "  the freeze density is conditional on theta and the freeze is not: over 15-60 deg the\n"
        "  crossing moves less than a decade and the store is >80% held at both edges."
    )


if __name__ == "__main__":
    main()
