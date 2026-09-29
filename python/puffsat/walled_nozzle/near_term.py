"""The near-term walled chamber: the 2.5 kg rod into a 20 m^3 chamber, solved (asks A1-A4).

The paper (`sec:methane_7000_near_term`) proposes methane at 7000 K on a pitch-lined steel wall
as the near-term walled chamber, and prices it with efficiencies **borrowed** from this study's
ADR-0016 cells: 10 kK, 200 m^3, a 25 kg impactor. This module replaces the borrowing with a
solve on the paper's own geometry:

- a 2.5 kg polyethylene rod at 75 km/s into a 20 m^3 chamber (a sphere: the paper's 35.6 m^2
  of wall is `4 pi r^2` at `r = 1.684 m`);
- a 4.4 kg polyethylene plug, which leaves with the exhaust and is carried propellant;
- throats of 0.149-0.198 m^2 (paper `tab:pulse_size_entrance`), opened as a 15 degree cone.

**The charge is a fixed point, as in `chamber.py`**, with the rod and plug in it:

    (m_charge + m_rod + m_plug) u(rho, T; feed) = m_rod w^2 / 2

where the feed -- and so the H:C ratio -- itself depends on `m_charge` (`eos_methane.Feed`).

**Three carbon branches (A2).** One isentrope each, because in each the whole pulse, chamber
blowdown and nozzle alike, follows a single constrained equilibrium:

- `equilibrium` -- every species in gas-phase equilibrium. The companion's current treatment.
- `frozen_carbon` -- every carbon-bearing species frozen at the *peak chamber* composition per
  kilogram; hydrogen alone re-equilibrates (`H + H + M` is the one channel with measured rates,
  and `freeze.py` shows it keeps up). **The floor that matters for the reactor comparison.**
- `graphite_ceiling` -- gas-phase equilibrium *plus* condensed carbon wherever monatomic carbon
  exceeds its saturation density over graphite. **A thermodynamic ceiling, not a nucleation
  result**: it assumes soot forms and stays in velocity and thermal equilibrium with the gas at no
  kinetic cost, which nothing on a millisecond clock does (ADR-0050, N10 item 5).

**Two efficiencies, and they answer different questions.**

- `eta_peak` -- the companion's convention: a steady expansion from the *peak* chamber state,
  `v^2 / 2u`. Its `v` comes from the stagnation enthalpy `h0 = u + p/rho`, so it credits the flow
  work `p/rho` -- about 14% of `u` at 493 bar -- to every kilogram.
- `eta_blowdown` -- the pulse as it happens. The chamber empties along the same isentrope, each
  parcel leaving with the stagnation enthalpy of the chamber *at the time it leaves*. Energy
  conservation is exact: `integral h drho = rho0 u0` along an isentrope, so the whole pulse hands
  out `U0` and not `M h0`. Late parcels leave a cooler, emptier chamber and expand less well.
  **This is the number the paper's `eta` means** ("the share of the rod's kinetic energy that
  leaves as directed exhaust").

The specific impulse uses the **impulse-averaged** exhaust speed, `integral v dm / M`, which is
what momentum sees. It sits below `w sqrt(eta / (k+1+P))` (the paper's eq:eta_isp) by Jensen's
inequality, because a blowdown's exhaust speed is not uniform.

What this module does **not** do: the wall is adiabatic in the expansion (its heat, from A4, is
reported beside `eta` as a debit rather than folded in); there is no boundary layer, no
two-dimensional loss, no divergence loss; and carbon association has no rate -- the branches
bracket it because the literature has no coefficient to use (`rates.CARBON_CHANNEL_NOTE`).
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat import eos_methane, expansion
from puffsat.eos_water import G_H, IP_H, ln_k_saha
from puffsat.walled_nozzle import chamber, freeze, rates, surface, wall

Vec = NDArray[np.float64]

# ---- The paper's near-term geometry (sec:rod_2p5kg, sec:methane_7000_near_term) --------------

#: Rod (projectile) mass [kg]. Polyethylene; it arrives at the closing speed and is not carried.
ROD_MASS = 2.5
#: Plug mass [kg]. Polyethylene foam; carried, and leaves with the exhaust.
PLUG_MASS = 4.4
#: `P` of paper eq:eta_isp: plug per kilogram of rod.
PLUG_RATIO = PLUG_MASS / ROD_MASS
#: Chamber volume [m^3].
VOLUME = 20.0
#: Chamber wall area [m^2]. `4 pi r^2` for a 20 m^3 sphere, to three figures.
WALL_AREA = 35.6
CLOSING_SPEED = chamber.CLOSING_SPEED
#: Pulse energy [J]: the rod's kinetic energy, 7.03 GJ.
PULSE_ENERGY = 0.5 * ROD_MASS * CLOSING_SPEED**2
#: Throat areas [m^2], the fixed-shape throats of the paper's tab:pulse_size_entrance.
THROAT_AREAS = (0.149, 0.198)
#: Exit area ratios: the ask's 4-100, plus 300-1000, which is in family for vacuum upper stages
#: (RL10B-2 flies 285:1 on a carbon-carbon extension; high-expansion studies run 500-1000:1). On
#: these throats 1000:1 is a 14-16 m exit -- the ratio is conventional, the absolute size is not.
AREA_RATIOS = (4.0, 14.0, 30.0, 100.0, 300.0, 500.0, 1000.0)
#: Nozzle cone half-angle [rad]. **The one assumed shape**: it sets only the freeze clock, never
#: `eta`, which depends on the area ratio alone.
CONE_HALF_ANGLE = math.radians(15.0)
#: Wall surface temperature [K] under the pulse: the paper's charring-pitch surface, held at
#: graphite's 1 atm sublimation point (A5 questions this; see `SURFACE_NOTE`).
PITCH_SURFACE_TEMPERATURE = 3900.0

BRANCHES = ("equilibrium", "frozen_carbon", "graphite_ceiling")

OUTPUT_DIR = Path("data/results/walled_nozzle/near_term")


@dataclass(frozen=True)
class Pairing:
    """One gas at one chamber temperature, as the paper's tab:wall_pairings states it."""

    name: str
    temp: float
    charge: eos_methane.FeedComponent
    #: The paper's own charge [kg] and peak pressure [bar], for the comparison.
    paper_charge: float
    paper_pressure_bar: float


PAIRINGS = (
    Pairing("hydrogen_5500", 5500.0, eos_methane.HYDROGEN_GAS, 59.0, 822.0),
    Pairing("methane_7000", 7000.0, eos_methane.METHANE_GAS, 62.0, 493.0),
    Pairing("methane_10000", 10000.0, eos_methane.METHANE_GAS, 45.0, 629.0),
)


#: Partition-function sets the chamber is solved on. **`janaf` is primary**: the harmonic C3 is
#: over-bound ~30x at 6000 K and C2 under-bound ~8x, which moves the 7000 K charge by 5% and the
#: C3 exposure from 16% to 2% (`eos_methane.janaf_correction`). `harmonic` is carried as the
#: sensitivity, and because every earlier result in this study is on it.
THERMOS = ("janaf", "harmonic")


def feed_for(pairing: Pairing, charge_mass: float, thermo: str = "janaf") -> eos_methane.Feed:
    """The chamber's feed: the charge plus the polyethylene rod and plug."""
    return eos_methane.mixture_feed(
        f"{pairing.name}+PE",
        [(charge_mass, pairing.charge), (ROD_MASS + PLUG_MASS, eos_methane.POLYETHYLENE)],
        thermo,
    )


# ---- A1: the chamber ------------------------------------------------------------------------


@dataclass(frozen=True)
class ChamberSolution:
    """The solved chamber for one pairing."""

    pairing: Pairing
    charge_mass: float
    feed: eos_methane.Feed
    #: Chamber volume [m^3]. A sphere: the wall, radius and beam length all follow from it.
    volume: float
    rho: float
    pressure: float
    #: Specific internal energy [J/kg] -- `PULSE_ENERGY / total mass` at the fixed point.
    energy: float
    sound_speed: float
    hydrogen_dissociated: float
    carbon_free: float
    #: Share of `u` held in broken bonds (ionisation excluded) -- the paper's "bond share".
    bond_share: float
    ionisation_fraction: float
    #: Upper bound on the atomisation store C3 withholds, as a fraction of the feed's full store.
    #: The model's worst exposure: C3's harmonic bend (see `chamber.ChamberState`).
    c3_store_exposure: float

    @property
    def radius(self) -> float:
        """Sphere radius [m]."""
        return float((3.0 * self.volume / (4.0 * math.pi)) ** (1.0 / 3.0))

    @property
    def wall_area(self) -> float:
        """Wall area [m^2], `4 pi r^2` -- 35.6 m^2 at the paper's 20 m^3."""
        return 4.0 * math.pi * self.radius**2

    @property
    def beam_length(self) -> float:
        """Mean beam length [m] for gas radiation to the whole wall, `4 V / A` (Hottel)."""
        return 4.0 * self.volume / self.wall_area

    @property
    def total_mass(self) -> float:
        """Everything that leaves [kg]: charge, rod and plug."""
        return self.charge_mass + ROD_MASS + PLUG_MASS

    @property
    def slug_ratio(self) -> float:
        """`k`: charge per kilogram of rod."""
        return self.charge_mass / ROD_MASS


def solve_charge(
    pairing: Pairing, thermo: str = "janaf", volume: float = VOLUME, tol: float = 1e-11
) -> ChamberSolution:
    """The charge that puts the chamber at `pairing.temp` when the rod's energy is shared.

    Damped fixed point on the total mass, as `chamber.solve_chamber`; the feed is rebuilt each
    step because the H:C ratio moves with the charge.
    """
    polyethylene = ROD_MASS + PLUG_MASS
    total = pairing.paper_charge + polyethylene
    for _ in range(300):
        feed = feed_for(pairing, total - polyethylene, thermo)
        u = eos_methane.pressure_energy(total / volume, pairing.temp, feed)[1]
        new = PULSE_ENERGY / u
        if abs(new - total) < tol * total:
            total = new
            break
        total = 0.5 * (total + new)
    else:  # pragma: no cover - the map is a contraction, as in chamber.py
        raise RuntimeError(f"near-term fixed point did not converge for {pairing.name}")

    charge = total - polyethylene
    feed = feed_for(pairing, charge, thermo)
    rho = total / volume
    p, u = eos_methane.pressure_energy(rho, pairing.temp, feed)
    comp = eos_methane.composition(rho, pairing.temp, feed)
    n_f = rho / feed.formula_mass
    return ChamberSolution(
        pairing=pairing,
        charge_mass=charge,
        feed=feed,
        volume=volume,
        rho=rho,
        pressure=p,
        energy=u,
        sound_speed=eos_methane.sound_speed(rho, pairing.temp, feed),
        hydrogen_dissociated=comp.hydrogen_dissociated,
        carbon_free=comp.carbon_free,
        bond_share=eos_methane.bond_energy_held(rho, pairing.temp, feed) / u,
        ionisation_fraction=comp.n_e / comp.n_total,
        c3_store_exposure=comp.density_of(eos_methane.C3)
        * eos_methane.C3.d_at
        / (n_f * feed.reference_atomization),
    )


# ---- A2: the three carbon branches, each a (rho, T) -> (p, e) law ---------------------------


class Branch:
    """A constrained equilibrium: `(rho, T) -> composition`, priced by one accounting."""

    name = "equilibrium"

    def __init__(self, feed: eos_methane.Feed) -> None:
        self.feed = feed

    def state(self, rho: float, temp: float) -> tuple[eos_methane.Composition, float]:
        """`(gas composition, condensed carbon atoms per m^3)` at `(rho, temp)`."""
        return eos_methane.composition(rho, temp, self.feed), 0.0

    def condensed_energy(self, temp: float) -> float:
        """Energy [J] of one condensed carbon atom, from free ground-state C. None here."""
        return 0.0

    def condensed_thermal(self, temp: float) -> float:
        """Thermal part [J] of one condensed atom's energy, above the solid at 0 K. None here."""
        return 0.0

    def __call__(self, rho: float, temp: float) -> tuple[float, float]:
        """`(p [Pa], e [J/kg])` -- an `expansion.Eos`."""
        comp, n_s = self.state(rho, temp)
        p, e = eos_methane.state_from_composition(comp, rho, temp, self.feed)
        return p, e + n_s * self.condensed_energy(temp) / rho

    def thermal_share(self, rho: float, temp: float) -> float:
        """Thermal energy [J/kg] -- the part a frozen flow could still turn into speed.

        The condensate's own vibrational heat counts as thermal: it is in thermal equilibrium
        with the gas in this branch by assumption.
        """
        comp, n_s = self.state(rho, temp)
        solid = n_s * self.condensed_thermal(temp) / rho
        return eos_methane.thermal_energy(comp, temp, self.feed.thermo) / rho + solid


class FrozenCarbon(Branch):
    """Carbon species frozen per kilogram at a reference state; hydrogen re-equilibrates.

    The hydrogen not locked in a carbon-bearing molecule is split among `H`, `H2` and `H+` by
    mass action and Saha, with the frozen carbon ions' electrons in the charge balance:

        n_H + 2 n_H^2 / K_H2 + K_H n_H / n_e = N_free ;   n_e = Q + K_H n_H / n_e

    which is a quadratic in `n_H` at fixed `n_e` and a quadratic in `n_e` at fixed `n_H`, solved
    alternately. Both maps are contractions; the loop converges in a few passes.
    """

    name = "frozen_carbon"

    def __init__(self, feed: eos_methane.Feed, rho0: float, temp0: float) -> None:
        super().__init__(feed)
        comp0 = eos_methane.composition(rho0, temp0, feed)
        mols = eos_methane.MOLECULES
        self._mol_per_kg = tuple(
            0.0 if m is eos_methane.H2 else n / rho0 for m, n in zip(mols, comp0.n_mol, strict=True)
        )
        self._c_per_kg = comp0.n_c / rho0
        self._ions_per_kg = tuple(n / rho0 for n in comp0.n_c_ions)
        self._charge_per_kg = sum((j + 1) * n for j, n in enumerate(self._ions_per_kg))
        bound_h = sum(m.n_h * n for m, n in zip(mols, self._mol_per_kg, strict=True))
        self._h_free_per_kg = feed.hc_ratio / feed.formula_mass - bound_h
        self._i_h2 = mols.index(eos_methane.H2)

    def state(self, rho: float, temp: float) -> tuple[eos_methane.Composition, float]:
        """Frozen carbon at this density, equilibrium hydrogen at this temperature."""
        n_free = self._h_free_per_kg * rho
        q = self._charge_per_kg * rho
        ln_k2 = eos_methane.ln_k_molecule_thermo(eos_methane.H2, temp, self.feed.thermo)
        k_h = math.exp(min(ln_k_saha(IP_H, 1.0, G_H, temp), 700.0))  # n_H+ n_e / n_H

        n_e = max(q, math.sqrt(k_h * n_free), 1e-300)
        n_h = 0.0
        for _ in range(200):
            b = 1.0 + k_h / n_e
            x = math.log(8.0 * n_free) - ln_k2
            if x > 60.0:  # K_H2 negligible: essentially all H2
                n_h = math.exp(0.5 * (ln_k2 + math.log(n_free / 2.0)))
            else:
                n_h = 2.0 * n_free / (b + math.sqrt(b * b + math.exp(x)))
            n_e_new = max(0.5 * (q + math.sqrt(q * q + 4.0 * k_h * n_h)), 1e-300)
            if abs(n_e_new - n_e) <= 1e-13 * n_e:
                n_e = n_e_new
                break
            n_e = n_e_new
        n_hp = k_h * n_h / n_e
        n_h2 = math.exp(2.0 * math.log(n_h) - ln_k2) if n_h > 0.0 else 0.5 * n_free
        n_mol = tuple(n_h2 if i == self._i_h2 else n * rho for i, n in enumerate(self._mol_per_kg))
        comp = eos_methane.Composition(
            n_c=self._c_per_kg * rho,
            n_h=n_h,
            n_hp=n_hp,
            n_c_ions=tuple(n * rho for n in self._ions_per_kg),
            n_e=n_e,
            n_mol=n_mol,
        )
        return comp, 0.0


# NASA-Glenn 7-coefficient fit for C(gr), McBride/Gordon via Cantera `nasa_condensed.yaml`
# (note "X 4/83"), 200-1000 K and 1000-5000 K. Enthalpy referenced to graphite at 298.15 K.
_GRAPHITE_LOW = (
    -0.310872072, 4.40353686e-03, 1.90394118e-06, -6.38546966e-09, 2.98964248e-12,
    -108.650794, 1.11382953,
)  # fmt: skip
_GRAPHITE_HIGH = (
    1.45571829, 1.71702216e-03, -6.97562786e-07, 1.35277032e-10, -9.67590652e-15,
    -695.138814, -8.52583033,
)  # fmt: skip
_R_GAS = 8.314462618
#: `dHf(C, g)` at 0 K [J/mol] (JANAF C-003) and graphite's `H(298) - H(0)` [J/mol] (JANAF
#: C-002). Together they place graphite at 298 K `711185 - 1051` J/mol below free C at 0 K.
_D0_CARBON = 711.185e3
_GRAPHITE_H298_H0 = 1.051e3


def _graphite_coeffs(temp: float) -> tuple[float, ...]:
    return _GRAPHITE_LOW if temp < 1000.0 else _GRAPHITE_HIGH


def graphite_heat_capacity(temp: float) -> float:
    """`c_p` of graphite [J/kg/K], NASA-Glenn fit, floored at its 200 K edge."""
    a = _graphite_coeffs(max(temp, 200.0))
    t = max(temp, 200.0)
    cp_molar = _R_GAS * (a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] * t**4)
    return cp_molar / 12.011e-3


def graphite_enthalpy(temp: float) -> float:
    """`H(T) - H(298.15)` of graphite [J/mol], NASA-Glenn fit."""
    a = _graphite_coeffs(temp)
    t = temp
    return (
        _R_GAS
        * t
        * (a[0] + a[1] * t / 2 + a[2] * t**2 / 3 + a[3] * t**3 / 4 + a[4] * t**4 / 5 + a[5] / t)
    )


def graphite_entropy(temp: float) -> float:
    """`S(T)` of graphite [J/mol/K], NASA-Glenn fit."""
    a = _graphite_coeffs(temp)
    t = temp
    return _R_GAS * (
        a[0] * math.log(t) + a[1] * t + a[2] * t**2 / 2 + a[3] * t**3 / 3 + a[4] * t**4 / 4 + a[6]
    )


def graphite_chemical_potential(temp: float) -> float:
    """Chemical potential [J per atom] of graphite, measured from free ground-state C at rest."""
    h = graphite_enthalpy(temp) - (_D0_CARBON - _GRAPHITE_H298_H0)
    return (h - temp * graphite_entropy(temp)) / eos_methane.N_A


def graphite_saturation_density(temp: float) -> float:
    """Monatomic carbon density [m^-3] in equilibrium with graphite, on this EOS's own C(g).

    Checked against JANAF's `log Kf` for C(g) over C(ref) in the tests: within 2% to 5000 K and
    4% at 6000 K, where the NASA fit is extrapolated past its 5000 K range.
    """
    ln_n = eos_methane.ln_free_carbon_density_scale(temp)
    return math.exp(ln_n + graphite_chemical_potential(temp) / (eos_methane.K_B * temp))


class GraphiteCeiling(Branch):
    """Gas-phase equilibrium plus condensed carbon wherever the gas is supersaturated.

    Where monatomic carbon would exceed `graphite_saturation_density`, it is pinned there and the
    excess carbon condenses; hydrogen and electrons re-solve on the remaining gas. Melting is not
    carried: JANAF's reference state is graphite to 6000 K, and liquid carbon's 105-120 kJ/mol of
    fusion would only return if the condensate first formed as liquid.
    """

    name = "graphite_ceiling"

    def condensed_energy(self, temp: float) -> float:
        """Energy [J] of one graphite atom relative to free ground-state C at rest."""
        return self.condensed_thermal(temp) - _D0_CARBON / eos_methane.N_A

    def condensed_thermal(self, temp: float) -> float:
        """`H(T) - H(0)` of graphite per atom [J]; the NASA fit is floored at its 200 K edge."""
        t = max(temp, 200.0)
        return (graphite_enthalpy(t) + _GRAPHITE_H298_H0) / eos_methane.N_A

    def state(self, rho: float, temp: float) -> tuple[eos_methane.Composition, float]:
        """Gas equilibrium, then condense any supersaturated carbon."""
        comp = eos_methane.composition(rho, temp, self.feed)
        n_sat = graphite_saturation_density(temp)
        if comp.n_c <= n_sat:
            return comp, 0.0
        n_f = rho / self.feed.formula_mass
        ln_nc = math.log(n_sat)
        ln_h_target = math.log(self.feed.hc_ratio * n_f)
        thermo = self.feed.thermo
        y = np.array([math.log(comp.n_h), math.log(comp.n_e)])

        def resid(v: Vec) -> Vec:
            c = eos_methane.species_from_log(ln_nc, float(v[0]), float(v[1]), temp, thermo)
            charge = c.n_hp + sum((j + 1) * n for j, n in enumerate(c.n_c_ions))
            return np.array([math.log(c.n_h_nuclei) - ln_h_target, float(v[1]) - math.log(charge)])

        r = resid(y)
        for _ in range(200):
            if float(np.max(np.abs(r))) < 1e-11:
                break
            jac = np.empty((2, 2))
            for k in range(2):
                yp = y.copy()
                yp[k] += 1e-6
                jac[:, k] = (resid(yp) - r) / 1e-6
            step = np.linalg.solve(jac, -r)
            big = float(np.max(np.abs(step)))
            if big > 2.0:
                step *= 2.0 / big
            y = y + step
            r = resid(y)
        else:  # pragma: no cover - guards a future EOS change
            raise RuntimeError(f"graphite branch did not converge at rho={rho}, T={temp}")
        gas = eos_methane.species_from_log(ln_nc, float(y[0]), float(y[1]), temp, thermo)
        return gas, max(0.0, n_f - gas.n_c_nuclei)


def branch(name: str, sol: ChamberSolution) -> Branch:
    """The named branch for a solved chamber."""
    if name == "equilibrium":
        return Branch(sol.feed)
    if name == "frozen_carbon":
        return FrozenCarbon(sol.feed, sol.rho, sol.pairing.temp)
    if name == "graphite_ceiling":
        return GraphiteCeiling(sol.feed)
    raise ValueError(f"unknown branch {name!r}")


# ---- The isentrope, once per branch ---------------------------------------------------------


@dataclass(frozen=True)
class Isentrope:
    """One isentrope from the peak chamber state, as arrays (densest first)."""

    rho: Vec
    temp: Vec
    pressure: Vec
    energy: Vec

    @property
    def enthalpy(self) -> Vec:
        """`h = e + p / rho` [J/kg]."""
        return self.energy + self.pressure / self.rho


def _temperature_near(rho: float, e_target: float, eos: expansion.Eos, guess: float) -> float:
    """`expansion.temperature_at`, bracketed around the previous station's temperature.

    Same bisection, same monotonicity argument; the bracket is widened geometrically until it
    holds the target, so the answer is identical and the cost is a third.
    """
    lo, hi = 0.7 * guess, 1.05 * guess
    while eos(rho, lo)[1] > e_target:
        lo *= 0.7
        if lo < expansion.T_FLOOR:
            return expansion.temperature_at(rho, e_target, eos)
    while eos(rho, hi)[1] < e_target:
        hi *= 1.5
    while hi - lo > 1e-10 * hi:
        mid = 0.5 * (lo + hi)
        if eos(rho, mid)[1] < e_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def isentrope(
    eos: expansion.Eos,
    rho0: float,
    temp0: float,
    rho_ratio_end: float = 1.0e-9,
    steps: int = 1100,
    temp_floor: float = 150.0,
) -> Isentrope:
    """Integrate `de = -p d(1/rho)` from the chamber, exactly as `expansion.expand` does.

    Heun's predictor-corrector on a grid uniform in `ln v`. The differences from `expand` are a
    warm-started temperature bracket and a stop at `temp_floor`: a frozen-carbon gas expanded a
    millionfold cools towards the EOS's partition-function floor, and nothing below 150 K matters
    to an exit that the ask stops at `A/A* = 100`.
    """
    ratio = (1.0 / rho_ratio_end) ** (1.0 / steps)
    p, e = eos(rho0, temp0)
    rho, temp = rho0, temp0
    out = [(rho, temp, p, e)]
    for _ in range(steps):
        rho_next = rho / ratio
        dv = 1.0 / rho_next - 1.0 / rho
        e_pred = e - p * dv
        t_pred = _temperature_near(rho_next, e_pred, eos, temp)
        p_pred = eos(rho_next, t_pred)[0]
        e = e - 0.5 * (p + p_pred) * dv
        temp = _temperature_near(rho_next, e, eos, t_pred)
        p = eos(rho_next, temp)[0]
        rho = rho_next
        out.append((rho, temp, p, e))
        if temp < temp_floor:
            break
    arr = np.array(out)
    return Isentrope(rho=arr[:, 0], temp=arr[:, 1], pressure=arr[:, 2], energy=arr[:, 3])


@dataclass(frozen=True)
class NozzleExit:
    """A steady expansion from chamber station `j` to one area ratio."""

    area_ratio: float
    #: Sonic mass flux [kg/m^2/s] -- sets the throat's mass flow at this chamber state.
    flux_star: float
    speed: float
    rho: float
    temp: float
    pressure: float
    enthalpy: float
    #: Index of the sonic station, and of the first station past the exit, in the isentrope.
    i_star: int
    i_exit: int


def nozzle_exit(tab: Isentrope, j: int, area_ratio: float) -> NozzleExit | None:
    """Expand from chamber state `j` (at rest) to `area_ratio` on the same isentrope.

    Every station downstream of `j` on the isentrope is a candidate nozzle state; its speed is
    `sqrt(2 (h_j - h))` and its area `flux* / (rho u)`. The sonic flux is the maximum of `rho u`,
    refined by a parabola through the three samples around it (the maximum is flat, so this is
    good to the square of the step). `None` if the isentrope ends before the exit is reached.
    """
    h = tab.enthalpy
    h0 = float(h[j])
    speed = np.sqrt(np.maximum(0.0, 2.0 * (h0 - h[j:])))
    flux = tab.rho[j:] * speed
    i = int(np.argmax(flux))
    if i == 0 or i >= len(flux) - 1:
        return None
    f0, f1, f2 = float(flux[i - 1]), float(flux[i]), float(flux[i + 1])
    denom = f0 - 2.0 * f1 + f2
    flux_star = f1 - 0.125 * (f2 - f0) ** 2 / denom if denom < 0.0 else f1
    target = flux_star / area_ratio
    beyond = np.nonzero(flux[i:] < target)[0]
    if len(beyond) == 0:
        return None
    b = i + int(beyond[0])
    w = (float(flux[b - 1]) - target) / (float(flux[b - 1]) - float(flux[b]))

    def mix(arr: Vec) -> float:
        return float(arr[j + b - 1] + w * (arr[j + b] - arr[j + b - 1]))

    h_exit = mix(h)
    return NozzleExit(
        area_ratio=area_ratio,
        flux_star=flux_star,
        speed=math.sqrt(max(0.0, 2.0 * (h0 - h_exit))),
        rho=mix(tab.rho),
        temp=mix(tab.temp),
        pressure=mix(tab.pressure),
        enthalpy=h_exit,
        i_star=j + i,
        i_exit=j + b,
    )


@dataclass(frozen=True)
class Efficiency:
    """A1/A2's deliverable at one area ratio on one branch."""

    pairing: str
    #: Partition-function set (`THERMOS`).
    thermo: str
    branch: str
    area_ratio: float
    #: Companion convention: steady expansion from the peak state, `v^2 / 2 u0`.
    eta_peak: float
    #: The pulse as it happens: exhaust kinetic energy over the rod's, whole blowdown.
    eta_blowdown: float
    #: Exhaust speed from the peak state [m/s].
    speed_peak: float
    #: Impulse-averaged exhaust speed over the blowdown [m/s] -- what momentum sees.
    speed_impulse: float
    #: Mass fraction whose exit the isentrope could not reach (counted as zero speed).
    unresolved_mass: float
    exit_temp_peak: float
    exit_pressure_peak: float
    #: At the peak-state exit: thermal share and chemical share of `u0` that leave unconverted.
    exit_thermal_share: float
    exit_chemical_share: float
    #: Condensed carbon mass fraction at the peak-state exit (graphite ceiling only).
    exit_condensed_fraction: float
    slug_ratio: float

    @property
    def isp(self) -> chamber.IspLedger:
        """The ledger on the impulse-averaged speed, plug carried (paper eq:eta_isp)."""
        return chamber.isp_ledger(self.speed_impulse, self.slug_ratio, PLUG_RATIO)

    @property
    def isp_peak(self) -> chamber.IspLedger:
        """The ledger on the peak-state speed -- the companion's convention."""
        return chamber.isp_ledger(self.speed_peak, self.slug_ratio, PLUG_RATIO)


def efficiencies(
    sol: ChamberSolution, br: Branch, tab: Isentrope, area_ratios: Sequence[float] = AREA_RATIOS
) -> list[Efficiency]:
    """`eta_peak`, `eta_blowdown` and both exhaust speeds at every area ratio.

    The blowdown integral runs over chamber states `j` down the same isentrope: the parcel that
    leaves while the chamber is at `rho_j` carries `h_j - h_exit(j)` of kinetic energy, and the
    mass leaving is `V d rho_j`. Parcels leaving after the isentrope can no longer resolve their
    exit are counted at zero speed and reported as `unresolved_mass`.
    """
    rho0, u0 = float(tab.rho[0]), float(tab.energy[0])
    out = []
    for ar in area_ratios:
        exits: list[NozzleExit] = []
        for j in range(len(tab.rho)):
            ex = nozzle_exit(tab, j, ar)
            if ex is None:
                break
            exits.append(ex)
        n = len(exits)
        rho_j = tab.rho[:n]
        ke = np.array([0.5 * ex.speed**2 for ex in exits])
        v = np.array([ex.speed for ex in exits])
        drho = -np.diff(rho_j)
        eta_bd = float(np.sum(0.5 * (ke[1:] + ke[:-1]) * drho)) / (rho0 * u0)
        v_imp = float(np.sum(0.5 * (v[1:] + v[:-1]) * drho)) / rho0
        peak = exits[0]
        comp, n_s = br.state(peak.rho, peak.temp)
        thermal = br.thermal_share(peak.rho, peak.temp)
        e_exit = br(peak.rho, peak.temp)[1]
        out.append(
            Efficiency(
                pairing=sol.pairing.name,
                thermo=sol.feed.thermo,
                branch=br.name,
                area_ratio=ar,
                eta_peak=0.5 * peak.speed**2 / u0,
                eta_blowdown=eta_bd,
                speed_peak=peak.speed,
                speed_impulse=v_imp,
                unresolved_mass=float(rho_j[-1]) / rho0,
                exit_temp_peak=peak.temp,
                exit_pressure_peak=peak.pressure,
                exit_thermal_share=(thermal + peak.pressure / peak.rho) / u0,
                exit_chemical_share=(e_exit - thermal) / u0,
                exit_condensed_fraction=n_s * eos_methane.M_C / peak.rho,
                slug_ratio=sol.slug_ratio,
            )
        )
        del comp
    return out


# ---- The freeze clock on H + H + M, for the real (small) nozzle -----------------------------


@dataclass(frozen=True)
class FreezeCheck:
    """`H + H + M` against the expansion from the peak state, to one area ratio."""

    pairing: str
    #: Partition-function set (`THERMOS`).
    thermo: str
    throat_area: float
    area_ratio: float
    nozzle_length: float
    transit_time: float
    min_damkohler: float
    #: `log10(min Da / DA_EQUILIBRIUM)`: decades the rate could be wrong before it matters.
    margin_decades: float
    #: Whether any station sat outside the fitted window of the dominant rate.
    extrapolated: bool
    #: Area ratio of the Bray sudden-freeze station, the first with `Da < 1` (`inf` if none).
    freeze_area_ratio: float
    #: Atomic-H recombination energy held at the freeze station but returned by the exit in the
    #: equilibrium solution, as a share of `u0`: **the eta this channel can strand**. A negative
    #: margin at a station with no atomic H left costs nothing, which is why this, not the
    #: margin, is the verdict.
    store_at_risk: float
    #: The dominant rate's stated uncertainty [decades] (`rates.RateFit.decades`).
    rate_uncertainty: float
    #: `eta_blowdown` of the equilibrium branch at this area ratio, for the net figure below.
    eta_blowdown: float = math.nan

    @property
    def eta_net(self) -> float:
        """Equilibrium `eta_blowdown` less the stranded H store: **the eta to quote** once the
        H + H + M race is lost. An upper-bound debit (Bray sudden freeze), so a little low."""
        return self.eta_blowdown - self.store_at_risk


def _h_store(comp: eos_methane.Composition, rho: float) -> float:
    """Energy [J/kg] free H atoms would release recombining to H2: `n_H D(H2) / 2 / rho`."""
    return comp.n_h * 0.5 * eos_methane.H2.d_at / rho


def freeze_check(
    sol: ChamberSolution, tab: Isentrope, throat_area: float, area_ratio: float
) -> FreezeCheck:
    """Race `H + H + M -> H2 + M` against a 15 degree cone from the peak state.

    Same Damkoehler as `freeze.station`, on the same literature coefficients, with the clock from
    the cone: `r = r* sqrt(A/A*)`, `x = (r - r*) / tan(theta)`, `t = integral dx / u`.
    """
    ex = nozzle_exit(tab, 0, area_ratio)
    if ex is None:
        raise ValueError("isentrope does not reach the exit")
    h0 = float(tab.enthalpy[0])
    r_star = math.sqrt(throat_area / math.pi)
    idx = range(ex.i_star, ex.i_exit)
    xs, ts, das, extrap = [], [], [], False
    freeze_ar, store_freeze, store_exit = math.inf, math.nan, 0.0
    time = 0.0
    prev: tuple[float, float, float] | None = None  # x, speed, rho
    for i in idx:
        rho, temp = float(tab.rho[i]), float(tab.temp[i])
        speed = math.sqrt(max(0.0, 2.0 * (h0 - float(tab.enthalpy[i]))))
        ar = max(1.0, ex.flux_star / (rho * speed))
        x = (r_star * math.sqrt(ar) - r_star) / math.tan(CONE_HALF_ANGLE)
        if prev is not None:
            dt = (x - prev[0]) * 0.5 * (1.0 / speed + 1.0 / prev[1])
            time += dt
            if dt > 0.0 and prev[2] > rho:
                tau_exp = dt / math.log(prev[2] / rho)
                comp = eos_methane.composition(rho, temp, sol.feed)
                k = freeze.mixture_coefficient(temp, comp)
                tau_rec = freeze.recombination_time(k, comp.n_h, comp.n_neutral_heavy)
                da = tau_exp / tau_rec if tau_rec > 0.0 else math.inf
                das.append(da)
                extrap = extrap or rates.H_ATOMIC_H.extrapolated(temp)
                store_exit = _h_store(comp, rho)
                if da < 1.0 and math.isinf(freeze_ar):
                    freeze_ar, store_freeze = ar, store_exit
        xs.append(x)
        ts.append(time)
        prev = (x, speed, rho)
    length = (r_star * math.sqrt(area_ratio) - r_star) / math.tan(CONE_HALF_ANGLE)
    min_da = min(das) if das else math.inf
    return FreezeCheck(
        pairing=sol.pairing.name,
        thermo=sol.feed.thermo,
        throat_area=throat_area,
        area_ratio=area_ratio,
        nozzle_length=length,
        transit_time=time,
        min_damkohler=min_da,
        margin_decades=math.log10(min_da / freeze.DA_EQUILIBRIUM),
        extrapolated=extrap,
        freeze_area_ratio=freeze_ar,
        store_at_risk=0.0
        if math.isinf(freeze_ar)
        else (store_freeze - store_exit) / float(tab.energy[0]),
        rate_uncertainty=rates.H_ATOMIC_H.decades(),
    )


# ---- A4: the blowdown clock and the wall ----------------------------------------------------

#: H- electron affinity [J] (Lykke, Murray & Lineberger 1991: 0.754195 eV).
H_MINUS_AFFINITY = 0.754195 * eos_methane.EV
#: H- bound-free cross-section floor [m^2] across 0.25-1.55 um. The cross-section peaks at
#: 4e-21 m^2 near 0.85 um and stays above 1e-21 m^2 over this band (Wishart 1979; John 1988),
#: so this is a floor. The band holds 90% of a 7000 K blackbody's flux; nothing is claimed
#: beyond it, where H- free-free, C- (affinity 1.26 eV) and C2/C3/CH bands add opacity.
H_MINUS_BF_FLOOR = 1.0e-21
#: Fraction of a blackbody's flux inside that band, at the chamber temperatures here.
_H_MINUS_BAND = (0.25e-6, 1.55e-6)
_ELECTRON_MASS = 9.1093837e-31
_PLANCK = 6.62607015e-34
_C2 = 1.438776877e-2  # second radiation constant [m K]


def _planck_fraction(lam_t: float) -> float:
    """Fraction of blackbody flux below `lambda T` [m K], the standard series."""
    x = _C2 / lam_t
    total = 0.0
    for n in range(1, 60):
        nx = n * x
        total += math.exp(-nx) / n * (nx**3 + 3 * nx**2 + 6 * nx + 6) / n**3
    return 15.0 / math.pi**4 * total


def band_fraction(temp: float) -> float:
    """Share of a `temp` blackbody's flux inside the H- bound-free floor's band."""
    lo, hi = _H_MINUS_BAND
    return _planck_fraction(hi * temp) - _planck_fraction(lo * temp)


def h_minus_density(comp: eos_methane.Composition, temp: float) -> float:
    """H- density [m^-3] by Saha against H + e- (g_H- = 1, g_H = 2, g_e = 2)."""
    lam3 = (_PLANCK**2 / (2.0 * math.pi * _ELECTRON_MASS * eos_methane.K_B * temp)) ** 1.5
    boltzmann = math.exp(H_MINUS_AFFINITY / (eos_methane.K_B * temp))
    return float(comp.n_h * comp.n_e * lam3 / 4.0 * boltzmann)


def optical_depth_floor(comp: eos_methane.Composition, temp: float, length: float) -> float:
    """In-band optical depth across `length` [m], on H- bound-free alone."""
    return h_minus_density(comp, temp) * H_MINUS_BF_FLOOR * length


@dataclass(frozen=True)
class BlowdownRow:
    """One chamber state during the blowdown."""

    time: float
    mass_fraction: float
    temp: float
    pressure: float
    tau_floor: float


def blowdown(
    sol: ChamberSolution, br: Branch, tab: Isentrope, throat_area: float
) -> list[BlowdownRow]:
    """Chamber state against time: `V d rho / dt = -flux*(rho) A*`, on the isentrope.

    Each chamber state's sonic flux comes from the same steady expansion `nozzle_exit` solves, so
    the clock is the choked mass flow of the real gas rather than a `gamma = 1.2` coefficient.
    """
    rows: list[BlowdownRow] = []
    time = 0.0
    prev: tuple[float, float] | None = None  # rho, flux*
    rho0 = float(tab.rho[0])
    for j in range(len(tab.rho)):
        ex = nozzle_exit(tab, j, 2.0)
        if ex is None:
            break
        rho, temp = float(tab.rho[j]), float(tab.temp[j])
        if prev is not None:
            time += (
                sol.volume
                * (prev[0] - rho)
                * 0.5
                * (1.0 / ex.flux_star + 1.0 / prev[1])
                / throat_area
            )
        comp, _ = br.state(rho, temp)
        rows.append(
            BlowdownRow(
                time=time,
                mass_fraction=rho / rho0,
                temp=temp,
                pressure=float(tab.pressure[j]),
                tau_floor=optical_depth_floor(comp, temp, sol.beam_length),
            )
        )
        prev = (rho, ex.flux_star)
        if rho / rho0 < 1.0e-3:
            break
    return rows


def paper_emptying_time(sol: ChamberSolution, throat_area: float) -> float:
    """The paper's emptying time [s]: `surface.blowdown_time`'s `M / (0.65 rho a A*)`."""
    return surface.blowdown_time(sol.rho, sol.sound_speed, sol.volume, throat_area)


@dataclass(frozen=True)
class WallHeat:
    """A4's deliverable for one pairing and throat."""

    pairing: str
    #: Partition-function set (`THERMOS`).
    thermo: str
    throat_area: float
    volume: float
    wall_area: float
    #: Time for the chamber mass to fall by `e` [s] -- the quantity the paper's emptying time,
    #: `M / (0.65 rho a A*)`, estimates.
    t_efold: float
    #: Time to empty 90% and 99% of the mass [s], on the real-gas clock.
    t90: float
    t99: float
    #: The paper's emptying time [s], for comparison.
    t_paper: float
    #: Effective full-load heating times [s]: `integral (q/q0) dt` for opaque radiation (T^4)
    #: and for the Bartz convective scaling (`p^0.8 (T - T_w)`).
    t_eff_radiation: float
    t_eff_convection: float
    tau_floor_peak: float
    #: Time [s] at which the H- in-band optical depth floor falls below 1.
    t_thin: float
    radiation_peak: float
    #: Bartz convective flux at the liner [W/m^2]: frozen-c_p low edge, equilibrium-c_p high edge.
    convection_low: float
    convection_high: float
    #: Fluence per pulse [J/m^2]: radiation floor (in-band, H- only) + convection low edge,
    #: and opaque blackbody + convection high edge.
    fluence_low: float
    fluence_high: float

    @property
    def wall_energy_share(self) -> tuple[float, float]:
        """Share of the pulse energy the whole wall takes, `(low, high)` -- a debit on `eta`."""
        return (
            self.fluence_low * self.wall_area / PULSE_ENERGY,
            self.fluence_high * self.wall_area / PULSE_ENERGY,
        )


def _fluid(sol: ChamberSolution) -> surface.Fluid:
    """A `surface.Fluid` on this feed, for `wall.frozen_heat_capacity`."""
    feed = sol.feed
    return surface.Fluid(
        name=sol.pairing.name,
        pressure_energy=feed.pressure_energy,
        sound_speed=feed.sound_speed,
        held=lambda r, t: eos_methane.bond_energy_held(r, t, feed),
        full_atomisation=feed.full_atomization_energy,
        damkohler=lambda r, t, tau: (math.nan, math.nan),
    )


def bartz_liner(
    sol: ChamberSolution, tab: Isentrope, throat_area: float, wall_temp: float
) -> tuple[float, float]:
    """`(low, high)` Bartz flux [W/m^2] on the chamber liner at the peak state.

    Same correlation and station convention as `wall.bartz_flux`'s liner (`A/A*` of the bore),
    with the bracket spread to its corners: frozen `c_p` at `Pr = 0.8` and the large cross-section
    for the low edge, equilibrium `c_p` at `Pr = 0.5` and the small cross-section for the high
    edge. `c*` is the real-gas `p0 / flux*` rather than a `gamma = 1.2` coefficient.
    """
    temp = sol.pairing.temp
    ex = nozzle_exit(tab, 0, 2.0)
    assert ex is not None
    c_star = sol.pressure / ex.flux_star
    comp = eos_methane.composition(sol.rho, temp, sol.feed)
    mean_mass = sol.rho / comp.n_total
    cp_eq = wall.frozen_heat_capacity(_fluid(sol), sol.rho, temp)
    cp_fr = 2.5 * eos_methane.K_B / mean_mass
    d_star = 2.0 * math.sqrt(throat_area / math.pi)
    area_ratio = math.pi * sol.radius**2 / throat_area
    sigma = 1.0 / (0.5 * (wall_temp / temp) + 0.5) ** 0.68

    def flux(cp: float, prandtl: float, cross: float) -> float:
        mu = wall.chapman_enskog_viscosity(mean_mass, temp, cross)
        h = wall.bartz_coefficient(d_star, mu, cp, prandtl, sol.pressure, c_star, sigma, area_ratio)
        return h * (temp - wall_temp)

    lo_x, hi_x = wall.CROSS_SECTION_BRACKET
    lo_pr, hi_pr = wall.PRANDTL_BRACKET
    return flux(cp_fr, hi_pr, hi_x), flux(cp_eq, lo_pr, lo_x)


def wall_heat(
    sol: ChamberSolution,
    br: Branch,
    tab: Isentrope,
    throat_area: float,
    wall_temp: float = PITCH_SURFACE_TEMPERATURE,
) -> WallHeat:
    """Integrate radiation and convection over the real blowdown history."""
    rows = blowdown(sol, br, tab, throat_area)
    t0, p0 = rows[0].temp, rows[0].pressure
    sigma_sb = expansion.SIGMA_SB
    q_lo, q_hi = bartz_liner(sol, tab, throat_area, wall_temp)

    def conv_shape(r: BlowdownRow) -> float:
        return float((r.pressure / p0) ** 0.8 * max(0.0, r.temp - wall_temp) / (t0 - wall_temp))

    def rad_floor(r: BlowdownRow) -> float:
        return sigma_sb * r.temp**4 * band_fraction(r.temp) * (1.0 - math.exp(-r.tau_floor))

    t_rad = t_conv = f_lo = f_hi = 0.0
    for a, b in pairwise(rows):
        dt = b.time - a.time
        t_rad += dt * 0.5 * ((a.temp / t0) ** 4 + (b.temp / t0) ** 4)
        t_conv += dt * 0.5 * (conv_shape(a) + conv_shape(b))
        f_lo += dt * 0.5 * (rad_floor(a) + rad_floor(b) + q_lo * (conv_shape(a) + conv_shape(b)))
        f_hi += (
            dt * 0.5 * (sigma_sb * (a.temp**4 + b.temp**4) + q_hi * (conv_shape(a) + conv_shape(b)))
        )

    def time_at(frac: float) -> float:
        return next((r.time for r in rows if r.mass_fraction <= frac), rows[-1].time)

    return WallHeat(
        pairing=sol.pairing.name,
        thermo=sol.feed.thermo,
        throat_area=throat_area,
        volume=sol.volume,
        wall_area=sol.wall_area,
        t_efold=time_at(math.exp(-1.0)),
        t90=time_at(0.1),
        t99=time_at(0.01),
        t_paper=paper_emptying_time(sol, throat_area),
        t_eff_radiation=t_rad,
        t_eff_convection=t_conv,
        tau_floor_peak=rows[0].tau_floor,
        t_thin=next((r.time for r in rows if r.tau_floor < 1.0), math.inf),
        radiation_peak=sigma_sb * t0**4,
        convection_low=q_lo,
        convection_high=q_hi,
        fluence_low=f_lo,
        fluence_high=f_hi,
    )


# ---- A6: the entrance leak against wall radiation, re-solved at 7000 K ----------------------

#: The rod port [m] (paper sec:leak_radiation_trade): set by the rod and its clearance, so it does
#: not grow with the chamber.
ENTRANCE_DIAMETER = 0.109
#: Volume over throat area [m], held as the chamber grows so the cycle timing does not change:
#: 20 m^3 over the 0.149 and 0.198 m^2 throats.
VOLUME_PER_THROAT = tuple(VOLUME / a for a in THROAT_AREAS)
#: Chamber volumes the trade is re-solved over [m^3].
TRADE_VOLUMES = (20.0, 30.0, 45.0, 70.0, 100.0, 150.0)


@dataclass(frozen=True)
class LeakRadiation:
    """One chamber size on the paper's leak-radiation trade, with the gas state solved."""

    pairing: str
    thermo: str
    volume: float
    volume_per_throat: float
    pressure: float
    charge_mass: float
    radius: float
    throat_area: float
    #: Entrance over throat area; leaked gas goes the wrong way, so its thrust cost is `2 f`.
    leak_share: float
    #: Radiated + convected share of the pulse, `(low, high)`, from `wall_heat` on the real
    #: blowdown. Its thrust cost is `r / 2` (the paper's model).
    wall_share: tuple[float, float]
    #: In-band H- optical depth floor across a quarter radius at peak: the paper's opacity edge
    #: is where this falls to one.
    tau_quarter_radius: float

    @property
    def leak_cost(self) -> float:
        """`X = 2 f`."""
        return 2.0 * self.leak_share

    @property
    def wall_cost(self) -> tuple[float, float]:
        """`Y = r / 2`, low and high edges."""
        return 0.5 * self.wall_share[0], 0.5 * self.wall_share[1]

    @property
    def total_cost(self) -> tuple[float, float]:
        """`X + Y`, low and high edges."""
        return self.leak_cost + self.wall_cost[0], self.leak_cost + self.wall_cost[1]


def leak_radiation_trade(
    pairing: Pairing = PAIRINGS[1],
    volumes: Sequence[float] = TRADE_VOLUMES,
    thermo: str = "janaf",
) -> list[LeakRadiation]:
    """Re-solve `sec:leak_radiation_trade` for a pairing: the charge, the blowdown and the wall
    load at every volume, rather than scaling a fixed gas state as `p` and `p^-2/3`.

    The chamber stays at the pairing's temperature as it grows, so its density falls, the gas
    dissociates further and the charge re-solves; the paper's `pV = const` is the ideal limit of
    this. The wall load is the full `wall_heat` over the real blowdown (radiation *and* convection,
    pitch surface at 3900 K), where the paper carried radiation alone.
    """
    entrance = math.pi * (0.5 * ENTRANCE_DIAMETER) ** 2
    out = []
    for volume in volumes:
        sol = solve_charge(pairing, thermo, volume)
        br = Branch(sol.feed)
        tab = isentrope(br, sol.rho, pairing.temp, rho_ratio_end=1e-4, steps=400)
        comp = eos_methane.composition(sol.rho, pairing.temp, sol.feed)
        tau_q = optical_depth_floor(comp, pairing.temp, 0.25 * sol.radius)
        for ratio in VOLUME_PER_THROAT:
            throat = volume / ratio
            heat = wall_heat(sol, br, tab, throat)
            out.append(
                LeakRadiation(
                    pairing=pairing.name,
                    thermo=thermo,
                    volume=volume,
                    volume_per_throat=ratio,
                    pressure=sol.pressure,
                    charge_mass=sol.charge_mass,
                    radius=sol.radius,
                    throat_area=throat,
                    leak_share=entrance / throat,
                    wall_share=heat.wall_energy_share,
                    tau_quarter_radius=tau_q,
                )
            )
    return out


# ---- Driver ---------------------------------------------------------------------------------

Report = Callable[[str], None]


def run(
    pairings: Sequence[Pairing] = PAIRINGS,
    branches: Sequence[str] = BRANCHES,
    thermos: Sequence[str] = THERMOS,
    report: Report = print,
) -> tuple[
    list[ChamberSolution], list[Efficiency], list[FreezeCheck], list[WallHeat], dict[str, Isentrope]
]:
    """Everything, for every partition-function set, pairing and branch."""
    sols, effs, freezes, walls = [], [], [], []
    tables: dict[str, Isentrope] = {}
    for thermo, pairing in ((t, p) for t in thermos for p in pairings):
        sol = solve_charge(pairing, thermo)
        sols.append(sol)
        report(
            f"[{thermo}] {pairing.name}: charge {sol.charge_mass:.2f} kg "
            f"(paper {pairing.paper_charge:.0f}), "
            f"k = {sol.slug_ratio:.2f}, p = {sol.pressure / 1e5:.0f} bar, "
            f"u = {sol.energy / 1e6:.1f} MJ/kg, bond share {sol.bond_share:.3f}"
        )
        for name in branches:
            br = branch(name, sol)
            tab = isentrope(br, sol.rho, pairing.temp)
            tables[f"{thermo}/{pairing.name}/{name}"] = tab
            branch_effs = efficiencies(sol, br, tab)
            for eff in branch_effs:
                effs.append(eff)
                report(
                    f"  {name:17} A/A*={eff.area_ratio:5.0f}  eta_peak {eff.eta_peak:.3f}  "
                    f"eta_bd {eff.eta_blowdown:.3f}  Isp_eff {eff.isp.isp_effective:6.0f} s"
                )
            if name == "equilibrium":
                for throat in THROAT_AREAS:
                    for ar in AREA_RATIOS:
                        eta = next(x.eta_blowdown for x in branch_effs if x.area_ratio == ar)
                        freezes.append(
                            dataclasses.replace(
                                freeze_check(sol, tab, throat, ar), eta_blowdown=eta
                            )
                        )
                    walls.append(wall_heat(sol, br, tab, throat))
    return sols, effs, freezes, walls, tables


def _write(path: Path, header: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


def _cols(*parts: str) -> list[str]:
    """A CSV header from comma-separated chunks, so long headers stay readable in source."""
    return "".join(parts).split(",")


def _g(x: float, fmt: str) -> str:
    return format(x, fmt)


def write_all(
    sols: Sequence[ChamberSolution],
    effs: Sequence[Efficiency],
    freezes: Sequence[FreezeCheck],
    walls: Sequence[WallHeat],
    out: Path = OUTPUT_DIR,
) -> None:
    """Write the four ledgers."""
    _write(
        out / "chambers.csv",
        _cols(
            "pairing,thermo,temp_k,charge_kg,paper_charge_kg,slug_ratio,hc_ratio,rho_kg_m3,"
            "pressure_bar,paper_pressure_bar,energy_mj_kg,sound_speed_m_s,"
            "hydrogen_dissociated,carbon_free,bond_share,ionisation_fraction,"
            "c3_store_exposure"
        ),
        [
            [
                s.pairing.name, s.feed.thermo, _g(s.pairing.temp, ".0f"), _g(s.charge_mass, ".3f"),
                _g(s.pairing.paper_charge, ".1f"), _g(s.slug_ratio, ".4f"),
                _g(s.feed.hc_ratio, ".4f"), _g(s.rho, ".5f"), _g(s.pressure / 1e5, ".1f"),
                _g(s.pairing.paper_pressure_bar, ".0f"), _g(s.energy / 1e6, ".3f"),
                _g(s.sound_speed, ".0f"), _g(s.hydrogen_dissociated, ".4f"),
                _g(s.carbon_free, ".4f"), _g(s.bond_share, ".4f"),
                _g(s.ionisation_fraction, ".3e"), _g(s.c3_store_exposure, ".4f"),
            ]
            for s in sols
        ],
    )  # fmt: skip
    _write(
        out / "efficiency.csv",
        _cols(
            "pairing,thermo,branch,area_ratio,eta_peak,eta_blowdown,speed_peak_m_s,"
            "speed_impulse_m_s,unresolved_mass,isp_true_s,isp_effective_s,"
            "isp_effective_peak_s,exit_temp_k,exit_pressure_bar,exit_thermal_share,"
            "exit_chemical_share,exit_condensed_fraction"
        ),
        [
            [
                e.pairing, e.thermo, e.branch, _g(e.area_ratio, ".0f"), _g(e.eta_peak, ".4f"),
                _g(e.eta_blowdown, ".4f"), _g(e.speed_peak, ".0f"), _g(e.speed_impulse, ".0f"),
                _g(e.unresolved_mass, ".2e"), _g(e.isp.isp_true, ".1f"),
                _g(e.isp.isp_effective, ".1f"), _g(e.isp_peak.isp_effective, ".1f"),
                _g(e.exit_temp_peak, ".0f"), _g(e.exit_pressure_peak / 1e5, ".3f"),
                _g(e.exit_thermal_share, ".4f"), _g(e.exit_chemical_share, ".4f"),
                _g(e.exit_condensed_fraction, ".4f"),
            ]
            for e in effs
        ],
    )  # fmt: skip
    _write(
        out / "freeze.csv",
        _cols(
            "pairing,thermo,throat_m2,area_ratio,nozzle_length_m,transit_ms,min_damkohler,"
            "margin_decades,extrapolated,freeze_area_ratio,store_at_risk,",
            "rate_uncertainty_decades,eta_blowdown,eta_net"
        ),
        [
            [
                f.pairing, f.thermo, _g(f.throat_area, ".3f"), _g(f.area_ratio, ".0f"),
                _g(f.nozzle_length, ".2f"), _g(f.transit_time * 1e3, ".4f"),
                _g(f.min_damkohler, ".3e"), _g(f.margin_decades, ".2f"), f.extrapolated,
                _g(f.freeze_area_ratio, ".2f"), _g(f.store_at_risk, ".4f"),
                _g(f.rate_uncertainty, ".2f"), _g(f.eta_blowdown, ".4f"), _g(f.eta_net, ".4f"),
            ]
            for f in freezes
        ],
    )  # fmt: skip
    _write(
        out / "wall_heat.csv",
        _cols(
            "pairing,thermo,throat_m2,t_efold_ms,t90_ms,t99_ms,t_paper_ms,t_eff_radiation_ms,"
            "t_eff_convection_ms,tau_floor_peak,t_thin_ms,radiation_peak_mw_m2,"
            "convection_low_mw_m2,convection_high_mw_m2,fluence_low_mj_m2,fluence_high_mj_m2,"
            "wall_share_low,wall_share_high"
        ),
        [
            [
                w.pairing, w.thermo, _g(w.throat_area, ".3f"), _g(w.t_efold * 1e3, ".2f"),
                _g(w.t90 * 1e3, ".2f"),
                _g(w.t99 * 1e3, ".2f"), _g(w.t_paper * 1e3, ".2f"),
                _g(w.t_eff_radiation * 1e3, ".2f"), _g(w.t_eff_convection * 1e3, ".2f"),
                _g(w.tau_floor_peak, ".3e"), _g(w.t_thin * 1e3, ".2f"),
                _g(w.radiation_peak / 1e6, ".1f"), _g(w.convection_low / 1e6, ".1f"),
                _g(w.convection_high / 1e6, ".1f"), _g(w.fluence_low / 1e6, ".2f"),
                _g(w.fluence_high / 1e6, ".2f"), _g(w.wall_energy_share[0], ".4f"),
                _g(w.wall_energy_share[1], ".4f"),
            ]
            for w in walls
        ],
    )  # fmt: skip


def write_trade(rows: Sequence[LeakRadiation], out: Path = OUTPUT_DIR) -> None:
    """Write the A6 leak-radiation ledger."""
    _write(
        out / "leak_radiation.csv",
        _cols(
            "pairing,thermo,volume_m3,volume_per_throat_m,pressure_bar,charge_kg,radius_m,",
            "throat_m2,leak_share,leak_cost,wall_share_low,wall_share_high,total_cost_low,",
            "total_cost_high,tau_quarter_radius",
        ),
        [
            [
                r.pairing, r.thermo, _g(r.volume, ".0f"), _g(r.volume_per_throat, ".1f"),
                _g(r.pressure / 1e5, ".1f"), _g(r.charge_mass, ".2f"), _g(r.radius, ".3f"),
                _g(r.throat_area, ".4f"), _g(r.leak_share, ".4f"), _g(r.leak_cost, ".4f"),
                _g(r.wall_share[0], ".4f"), _g(r.wall_share[1], ".4f"),
                _g(r.total_cost[0], ".4f"), _g(r.total_cost[1], ".4f"),
                _g(r.tau_quarter_radius, ".3e"),
            ]
            for r in rows
        ],
    )  # fmt: skip


def main() -> None:
    """Run the near-term chamber asks and write the ledgers."""
    print("Near-term chamber (asks A1-A4): 2.5 kg rod, 20 m^3, 4.4 kg plug, 75 km/s\n")
    sols, effs, freezes, walls, _ = run()
    write_all(sols, effs, freezes, walls)
    print("\nH + H + M freeze margin (15 deg cone), decades above Da = 10:")
    for f in freezes:
        print(
            f"  {f.pairing:14} A*={f.throat_area:.3f} A/A*={f.area_ratio:4.0f} "
            f"L={f.nozzle_length:5.2f} m  margin {f.margin_decades:5.2f}"
            f"  freeze A/A*={f.freeze_area_ratio:6.1f}  at risk {f.store_at_risk:.4f} of u0"
            f"{'  (extrapolated)' if f.extrapolated else ''}"
        )
    print("\nWall heat (A4), pitch surface at 3900 K:")
    for w in walls:
        lo, hi = w.wall_energy_share
        print(
            f"  {w.pairing:14} A*={w.throat_area:.3f}  t90 {w.t90 * 1e3:5.1f} ms  "
            f"(paper {w.t_paper * 1e3:5.1f})  rad {w.radiation_peak / 1e6:5.0f} MW/m2 "
            f"tau>={w.tau_floor_peak:8.2e}  conv {w.convection_low / 1e6:5.0f}-"
            f"{w.convection_high / 1e6:5.0f} MW/m2  fluence {w.fluence_low / 1e6:5.2f}-"
            f"{w.fluence_high / 1e6:5.2f} MJ/m2  wall {lo:.3f}-{hi:.3f} of pulse"
        )
    print("\nA6 leak vs wall (methane 7000 K, JANAF): thrust cost X = 2f, Y = r/2")
    trade = leak_radiation_trade()
    write_trade(trade)
    for r in trade:
        print(
            f"  V={r.volume:5.0f} m3  V/A*={r.volume_per_throat:5.1f}"
            f"  p={r.pressure / 1e5:5.0f} bar"
            f"  X={r.leak_cost:.3f}  Y={r.wall_cost[0]:.3f}-{r.wall_cost[1]:.3f}"
            f"  total={r.total_cost[0]:.3f}-{r.total_cost[1]:.3f}"
            f"  tau(R/4)>={r.tau_quarter_radius:.2e}"
        )
    print(f"\nwrote {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
