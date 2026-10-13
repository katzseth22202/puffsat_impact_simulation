"""Finite-rate recombination in the walled nozzle: the paper repo's S26 (ADR-0057).

The impulses so far (`blowdown`, B6) take the nozzle in shifting equilibrium and then debit a
**Bray sudden freeze** of `H + H + M`: the atomic H held at the first station where the
Damkoehler number falls below one is written off. That is an upper bound on the debit. This
module replaces it with a quasi-1D steady expansion that carries the recombination at its
measured rate, from the throat to the exit, for every parcel of the blowdown.

**The one kinetic variable is the particle count per kilogram, `nu = n_total / rho`.**
Bimolecular exchanges (`CH + H -> C + H2`, `C2 + H -> CH + C`) conserve it, and only three-body
association and its reverse, dissociation, change it. The only three-body channel with evaluated
rates is `H + H + M` (`rates.py`), so it is the only sink, at its measured rate:

    d nu / dt = -k_eff n_M (n_H^2 - K_H2 n_H2) / rho,     k_eff = freeze.mixture_coefficient.

`K_H2` comes from the EOS's own thermochemistry, so the reverse rate is detailed balance and the
fast-rate limit is exactly the branch's isentrope (tested).

**Carbon has no rate (ADR-0050), so it is bracketed by two models**, each a composition law
`(rho, T, phi) -> composition` in which `phi` holds the particle count off equilibrium:

- `frozen_carbon`: every carbon bearer frozen per kilogram at the peak chamber state
  (`near_term.FrozenCarbon`); the free hydrogen alone moves, `n_H^2 / n_H2 = K phi`.
  **The floor**: carbon returns nothing, and cannot take up hydrogen.
- `exchange_carbon`: rate-controlled constrained equilibrium on the particle count
  (`eos_methane.constrained_composition`). Every bimolecular exchange is instantaneous, so carbon
  re-forms molecules freely, but only at the cost of particles that `H + H + M` has removed.
  **The upper realistic case**: it is equilibrium carbon in every respect except an unmeasured
  carbon three-body sink. Full equilibrium (`phi = 1` throughout) is the ceiling above it.

**The march.** Steady, adiabatic, inviscid, along the nozzle axis. Density is the independent
variable, stepped uniformly in `ln rho` from the throat. Per step, two unknowns `(T, ln phi)` and
two equations:

    e = e_i - (p_i + p) / 2 (1/rho - 1/rho_i)            (first law; trapezoid, as `isentrope`)
    nu = nu_i + dt R(rho, T, phi)                         (backward Euler: L-stable, so the
                                                           near-equilibrium throat is not stiff)

with `u = sqrt(2 (h0 - h))`, the area from continuity `A = mdot / (rho u)`, the axial station
from the cone `x = r* (sqrt(A/A*) - 1) / tan(theta)` (as `near_term.freeze_check`), and
`dt = dx (1/u + 1/u_i) / 2`. Solved by Newton with a finite-difference Jacobian.

**What it does not do.** The throat state is taken from the branch's shifting isentrope: the
subsonic gas is in equilibrium, and the march reports the throat's Damkoehler number to show it.
No boundary layer, wall heat or divergence loss; the cone has no throat rounding. The carbon
models are limits, and `required_carbon_rate` reports how fast an unmeasured carbon three-body
channel would have to run to approach the ceiling, beside the one measurement there is.
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat import eos_methane
from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import freeze, rates
from puffsat.walled_nozzle import near_term as nt

Vec = NDArray[np.float64]

MODELS = ("frozen_carbon", "exchange_carbon")
#: Area ratio the march starts from, on the branch's shifting isentrope. At the sonic point the
#: area is flat in density, so the cone's clock `dx / u` is ill-posed there (and the unrounded
#: cone is not the real throat); the gas is in equilibrium to here (`Expansion.throat_damkohler`).
START_AREA_RATIO = 1.05
#: March step in `ln rho`. Step halving moves the exit speed by under 0.1% (tested at 0.2%).
DEFAULT_STEP = 0.01
OUTPUT_DIR = nt.OUTPUT_DIR / "finite_rate"


# ---- The two composition laws ---------------------------------------------------------------


class Model:
    """A composition law `(rho, T, ln phi) -> composition`, and the branch it shifts along."""

    name = ""

    def __init__(self, branch: nt.Branch) -> None:
        self.branch = branch
        self.feed = branch.feed

    def composition(self, rho: float, temp: float, ln_phi: float) -> eos_methane.Composition:
        raise NotImplementedError


class FrozenCarbonModel(Model):
    """Carbon frozen at the peak chamber state; the free H pool held off equilibrium by `phi`."""

    name = "frozen_carbon"

    def __init__(self, branch: nt.FrozenCarbon) -> None:
        super().__init__(branch)
        self._frozen = branch

    def composition(self, rho: float, temp: float, ln_phi: float) -> eos_methane.Composition:
        return self._frozen.state(rho, temp, ln_phi)[0]


class ExchangeCarbonModel(Model):
    """Constrained equilibrium on the particle count; warm-started from the last solve."""

    name = "exchange_carbon"

    def __init__(self, branch: nt.Branch) -> None:
        super().__init__(branch)
        self._guess: Vec | None = None

    def composition(self, rho: float, temp: float, ln_phi: float) -> eos_methane.Composition:
        comp, self._guess = eos_methane.constrained_composition(
            rho, temp, self.feed, ln_phi, self._guess
        )
        return comp


def model(name: str, sol: nt.ChamberSolution) -> Model:
    """The named model for a solved chamber."""
    if name == "frozen_carbon":
        return FrozenCarbonModel(nt.FrozenCarbon(sol.feed, sol.rho, sol.pairing.temp))
    if name == "exchange_carbon":
        return ExchangeCarbonModel(nt.Branch(sol.feed))
    raise ValueError(f"unknown model {name!r}")


def sink_rate(
    comp: eos_methane.Composition, rho: float, temp: float, ln_phi: float, thermo: str
) -> float:
    """`d nu / dt` [1/(kg s)] from `H + H + M`, written as `-k n_M K n_H2 (phi - 1) / rho`.

    Both models satisfy `n_H^2 / n_H2 = K phi` exactly, so this is the same number as
    `-k n_M (n_H^2 - K n_H2) / rho` without its cancellation near equilibrium.
    """
    k = freeze.mixture_coefficient(temp, comp)
    ln_k2 = eos_methane.ln_k_molecule_thermo(eos_methane.H2, temp, thermo)
    n_h2 = comp.density_of(eos_methane.H2)
    return -k * comp.n_neutral_heavy * math.exp(ln_k2) * n_h2 * math.expm1(ln_phi) / rho


# ---- One parcel: the march from the throat ---------------------------------------------------


@dataclass(frozen=True)
class Expansion:
    """One steady expansion from one chamber state to one exit, at one rate."""

    model: str
    throat_area: float
    area_ratio: float
    half_angle: float
    rate_scale: float
    speed: float
    temp: float
    rho: float
    pressure: float
    #: Cone length [m] from the throat to this exit, and the gas's time along it [s].
    length: float
    transit_time: float
    ln_phi_exit: float
    #: Particles per kilogram at the throat and at the exit.
    particles_throat: float
    particles_exit: float
    #: Bond energy held [J/kg] in the chamber and at the exit.
    bond_chamber: float
    bond_exit: float
    #: Local Damkoehler number of `H + H + M` at the first step past the throat.
    throat_damkohler: float
    #: Area ratio of the first station with `Da < 1` (`inf` if none): where the race is lost.
    freeze_area_ratio: float
    hydrogen_free_exit: float

    @property
    def effective_speed(self) -> float:
        """Vacuum thrust per unit mass flow [m/s], `u + p / (rho u)`.

        `eta` counts only `u` (the study's convention). At a fixed area ratio the two can order
        frozen and shifting flow differently, so both are reported.
        """
        return self.speed + self.pressure / (self.rho * self.speed)

    @property
    def bond_returned(self) -> float:
        """Share of the chamber's bond store returned by the exit."""
        return (self.bond_chamber - self.bond_exit) / self.bond_chamber


@dataclass(frozen=True)
class _Station:
    rho: float
    temp: float
    ln_phi: float
    p: float
    e: float
    u: float
    x: float
    ar: float
    nu: float
    time: float
    comp: eos_methane.Composition


def _station(
    mod: Model,
    rho: float,
    temp: float,
    ln_phi: float,
    h0: float,
    flux_star: float,
    r_star: float,
    tan_theta: float,
) -> _Station | None:
    comp = mod.composition(rho, temp, ln_phi)
    p, e = eos_methane.state_from_composition(comp, rho, temp, mod.feed)
    u2 = 2.0 * (h0 - e - p / rho)
    if u2 <= 0.0:
        return None
    u = math.sqrt(u2)
    ar = flux_star / (rho * u)
    x = r_star * (math.sqrt(max(ar, 1.0)) - 1.0) / tan_theta
    return _Station(rho, temp, ln_phi, p, e, u, x, ar, comp.n_total / rho, 0.0, comp)


def _step(
    mod: Model,
    prev: _Station,
    rho: float,
    guess: tuple[float, float],
    h0: float,
    flux_star: float,
    r_star: float,
    tan_theta: float,
    rate_scale: float,
) -> tuple[_Station, float] | None:
    """Advance one step to `rho`. Returns the station and its `dt`.

    Newton on `(T, zeta)`, `zeta = ln K_H2(T) + ln phi = ln(n_H^2 / n_H2)`, not on `(T, ln phi)`.
    Once the pool is far from equilibrium the composition depends on `T` and `ln phi` almost only
    through `zeta`, and in those variables the Jacobian is singular to four figures; in `zeta` the
    composition is nearly fixed and `T` carries the thermal energy alone.
    """
    thermo = mod.feed.thermo
    dv = 1.0 / rho - 1.0 / prev.rho

    def ln_k2(temp: float) -> float:
        return eos_methane.ln_k_molecule_thermo(eos_methane.H2, temp, thermo)

    def resid(temp: float, zeta: float) -> tuple[Vec, _Station, float] | None:
        ln_phi = zeta - ln_k2(temp)
        st = _station(mod, rho, temp, ln_phi, h0, flux_star, r_star, tan_theta)
        if st is None:
            return None
        dt = max(0.0, st.x - prev.x) * 0.5 * (1.0 / st.u + 1.0 / prev.u)
        rate = 0.0
        if rate_scale > 0.0:
            rate = rate_scale * sink_rate(st.comp, rho, temp, ln_phi, thermo)
        r1 = (st.e - (prev.e - 0.5 * (prev.p + st.p) * dv)) / h0
        r2 = (st.nu - prev.nu - dt * rate) / prev.nu
        return np.array([r1, r2]), st, dt

    temp = guess[0]
    zeta = ln_k2(temp) + guess[1]
    for _ in range(60):
        out = resid(temp, zeta)
        if out is None:
            temp *= 1.02  # too cold to flow: back towards the previous station
            continue
        r, st, dt = out
        if float(np.max(np.abs(r))) < 1e-11:
            return dataclasses.replace(st, time=prev.time + dt), dt
        dtemp, dzeta = 1e-6 * temp, 1e-6
        a = resid(temp + dtemp, zeta)
        b = resid(temp, zeta + dzeta)
        if a is None or b is None:
            temp *= 1.001
            continue
        jac = np.column_stack(((a[0] - r) / dtemp, (b[0] - r) / dzeta))
        delta = np.linalg.solve(jac, -r)
        delta[0] = max(-0.1 * temp, min(0.1 * temp, float(delta[0])))
        delta[1] = max(-0.5, min(0.5, float(delta[1])))
        temp, zeta = temp + float(delta[0]), zeta + float(delta[1])
        # The exchange model's inner solve is good to ~1e-11, which a fast rate amplifies in
        # r2; a Newton step that no longer moves the state is converged.
        if abs(float(delta[0])) < 1e-10 * temp and abs(float(delta[1])) < 1e-10:
            break
    out = resid(temp, zeta)
    if out is not None and abs(float(out[0][0])) < 1e-8 and abs(float(out[0][1])) < 1e-6:
        return dataclasses.replace(out[1], time=prev.time + out[2]), out[2]
    return _bracket(resid, prev, rho)


def _bracket(
    resid: Callable[[float, float], tuple[Vec, _Station, float] | None], prev: _Station, rho: float
) -> tuple[_Station, float]:
    """Fallback when Newton cycles: bisect `T` on the energy residual, `zeta` solved inside.

    Newton cycles where `e(T)` is not smooth. `eos_methane.janaf_correction` holds its `delta`
    below 1000 K but drops its slope, so `e` steps by ~2e-4 of `h0` there and the energy equation
    has no root. Bisection lands on the step; the next station starts from the `e` actually
    reached, so the error is spent once and does not accumulate.
    """

    def at(temp: float) -> tuple[Vec, _Station, float]:
        zeta = math.log(
            max(prev.comp.n_h, 1e-300) ** 2 / max(prev.comp.density_of(eos_methane.H2), 1e-300)
        )
        for _ in range(60):
            out = resid(temp, zeta)
            if out is None:
                raise RuntimeError("no flow in the bracket")
            b = resid(temp, zeta + 1e-6)
            if b is None:
                raise RuntimeError("no flow in the bracket")
            slope = (b[0][1] - out[0][1]) / 1e-6
            d = max(-0.5, min(0.5, -float(out[0][1]) / slope))
            zeta += d
            if abs(d) < 1e-11:
                break
        final = resid(temp, zeta)
        assert final is not None
        return final

    lo, hi = 0.8 * prev.temp, prev.temp
    if float(at(lo)[0][0]) > 0.0 or float(at(hi)[0][0]) < 0.0:
        raise RuntimeError(f"finite-rate step did not converge at rho={rho:.3e}")
    while hi - lo > 1e-9 * hi:
        mid = 0.5 * (lo + hi)
        if float(at(mid)[0][0]) < 0.0:
            lo = mid
        else:
            hi = mid
    r, st, dt = at(hi)
    if abs(float(r[0])) > 1e-3:
        raise RuntimeError(f"energy residual {float(r[0]):.2e} at rho={rho:.3e}")
    return dataclasses.replace(st, time=prev.time + dt), dt


def march(
    mod: Model,
    tab: nt.Isentrope,
    j: int,
    throat_area: float,
    area_ratio: float,
    half_angle: float = nt.CONE_HALF_ANGLE,
    rate_scale: float = 1.0,
    step: float = DEFAULT_STEP,
) -> Expansion:
    """Expand chamber state `j` of the branch isentrope `tab` through a cone to `area_ratio`.

    `rate_scale` multiplies the measured coefficient: 0 freezes the particle count at the
    throat, 1e6 recovers shifting equilibrium.
    """
    ex = nt.nozzle_exit(tab, j, area_ratio)
    if ex is None:
        raise ValueError("the branch isentrope does not reach this exit from this state")
    h0 = float(tab.enthalpy[j])
    r_star = math.sqrt(throat_area / math.pi)
    tan_theta = math.tan(half_angle)
    rho_j, temp_j = float(tab.rho[j]), float(tab.temp[j])
    bond_chamber = eos_methane.bond_energy_in(mod.composition(rho_j, temp_j, 0.0), rho_j, mod.feed)

    first: _Station | None = None
    for i in range(ex.i_star, ex.i_exit + 1):
        first = _station(
            mod, float(tab.rho[i]), float(tab.temp[i]), 0.0, h0, ex.flux_star, r_star, tan_theta
        )
        if first is not None and first.ar >= START_AREA_RATIO:
            break
    if first is None:
        raise ValueError("no flow past the throat")
    start = first
    prev, older = start, None
    throat_da, freeze_ar = math.nan, math.inf
    for _ in range(100_000):
        rho = prev.rho * math.exp(-step)
        if older is None:
            guess = (prev.temp, prev.ln_phi)
        else:  # extrapolate T and ln phi linearly in ln rho
            guess = (prev.temp**2 / older.temp, 2.0 * prev.ln_phi - older.ln_phi)
        res = _step(mod, prev, rho, guess, h0, ex.flux_star, r_star, tan_theta, rate_scale)
        if res is None:  # pragma: no cover - _step raises rather than returning None
            break
        st, dt = res
        if dt > 0.0 and rate_scale > 0.0:
            k = rate_scale * freeze.mixture_coefficient(st.temp, st.comp)
            tau_rec = freeze.recombination_time(k, st.comp.n_h, st.comp.n_neutral_heavy)
            da = (dt / step) / tau_rec
            if math.isnan(throat_da):
                throat_da = da
            if da < 1.0 and math.isinf(freeze_ar):
                freeze_ar = st.ar
        if st.ar >= area_ratio:
            w = (math.log(area_ratio) - math.log(prev.ar)) / (math.log(st.ar) - math.log(prev.ar))

            def lerp(a: float, b: float, w: float = w) -> float:
                return a + w * (b - a)

            bond_prev = eos_methane.bond_energy_in(prev.comp, prev.rho, mod.feed)
            bond_st = eos_methane.bond_energy_in(st.comp, st.rho, mod.feed)
            return Expansion(
                model=mod.name,
                throat_area=throat_area,
                area_ratio=area_ratio,
                half_angle=half_angle,
                rate_scale=rate_scale,
                speed=lerp(prev.u, st.u),
                temp=lerp(prev.temp, st.temp),
                rho=math.exp(lerp(math.log(prev.rho), math.log(st.rho))),
                pressure=lerp(prev.p, st.p),
                length=r_star * (math.sqrt(area_ratio) - 1.0) / tan_theta,
                transit_time=lerp(prev.time, st.time),
                ln_phi_exit=lerp(prev.ln_phi, st.ln_phi),
                particles_throat=start.nu,
                particles_exit=lerp(prev.nu, st.nu),
                bond_chamber=bond_chamber,
                bond_exit=lerp(bond_prev, bond_st),
                throat_damkohler=throat_da,
                freeze_area_ratio=freeze_ar,
                hydrogen_free_exit=lerp(
                    prev.comp.hydrogen_dissociated, st.comp.hydrogen_dissociated
                ),
            )
        older, prev = prev, st
    raise RuntimeError("march did not reach the exit")  # pragma: no cover


# ---- The pulse: every parcel of the blowdown -------------------------------------------------


@dataclass(frozen=True)
class Pulse:
    """One nozzle at one rate, on one carbon model, over the whole blowdown."""

    model: str
    throat_area: float
    area_ratio: float
    half_angle: float
    rate_scale: float
    #: Shifting equilibrium on this model's branch (the fine `near_term.efficiencies` value).
    eta_shifting: float
    eta_finite: float
    #: Particle count frozen at the throat (rate zero).
    eta_frozen: float
    speed_shifting: float
    speed_finite: float
    speed_frozen: float
    #: Impulse-averaged vacuum thrust per unit mass flow, `u + p / (rho u)`, finite over shifting.
    thrust_ratio_finite: float
    thrust_ratio_frozen: float
    #: Mass-weighted share of each parcel's chamber bond store its exit returns.
    bond_shifting: float
    bond_finite: float
    bond_frozen: float
    peak: Expansion
    parcels: int


def parcel_indices(tab: nt.Isentrope, count: int = 14, floor: float = 1e-3) -> list[int]:
    """Chamber stations spaced geometrically in density from the peak down to `floor` of it."""
    targets = np.geomspace(1.0, floor, count) * float(tab.rho[0])
    idx = sorted({int(np.argmin(np.abs(tab.rho - t))) for t in targets})
    return idx


def pulse(
    mod: Model,
    tab: nt.Isentrope,
    throat_area: float,
    area_ratio: float,
    half_angle: float = nt.CONE_HALF_ANGLE,
    rate_scale: float = 1.0,
    count: int = 14,
    frozen_cache: dict[tuple[str, int, int, float], Expansion] | None = None,
) -> Pulse:
    """The blowdown-weighted `eta`, impulse-averaged speed and bond share, finite-rate.

    The coarse parcel quadrature is used only as a ratio against the same quadrature on the
    shifting exits, which is then applied to the fine `near_term.efficiencies` value, so the
    quadrature's own error cancels.
    """
    rows: list[tuple[float, ...]] = []
    peak: Expansion | None = None
    for j in parcel_indices(tab, count):
        ex = nt.nozzle_exit(tab, j, area_ratio)
        if ex is None:
            break
        fin = march(mod, tab, j, throat_area, area_ratio, half_angle, rate_scale)
        # A frozen march has no clock, so neither the throat nor the cone angle enters it.
        key = (mod.name, id(tab), j, area_ratio)
        if frozen_cache is not None and key in frozen_cache:
            frz = frozen_cache[key]
        else:
            frz = march(mod, tab, j, throat_area, area_ratio, half_angle, 0.0)
            if frozen_cache is not None:
                frozen_cache[key] = frz
        peak = fin if peak is None else peak
        shift_comp = mod.composition(ex.rho, ex.temp, 0.0)
        bond_shift = eos_methane.bond_energy_in(shift_comp, ex.rho, mod.feed)
        b0 = fin.bond_chamber
        rows.append(
            (
                float(tab.rho[j]),
                ex.speed,
                fin.speed,
                frz.speed,
                b0,
                b0 - bond_shift,
                b0 - fin.bond_exit,
                b0 - frz.bond_exit,
                ex.speed + ex.pressure / (ex.rho * ex.speed),
                fin.effective_speed,
                frz.effective_speed,
            )
        )
    if peak is None or len(rows) < 2:
        raise ValueError("no parcel reaches this exit")
    arr = np.array(rows)
    drho = -np.diff(arr[:, 0])

    def trap(col: Vec) -> float:
        return float(np.sum(0.5 * (col[1:] + col[:-1]) * drho))

    ke = [trap(0.5 * arr[:, c] ** 2) for c in (1, 2, 3)]
    vel = [trap(arr[:, c]) for c in (1, 2, 3)]
    bond = [trap(arr[:, c]) / trap(arr[:, 4]) for c in (5, 6, 7)]
    thrust = [trap(arr[:, c]) for c in (8, 9, 10)]
    eff = _shifting_efficiency(mod, tab, area_ratio)
    eta_fine, v_fine = eff.eta_blowdown, eff.speed_impulse
    return Pulse(
        model=mod.name,
        throat_area=throat_area,
        area_ratio=area_ratio,
        half_angle=half_angle,
        rate_scale=rate_scale,
        eta_shifting=eta_fine,
        eta_finite=eta_fine * ke[1] / ke[0],
        eta_frozen=eta_fine * ke[2] / ke[0],
        speed_shifting=v_fine,
        speed_finite=v_fine * vel[1] / vel[0],
        speed_frozen=v_fine * vel[2] / vel[0],
        thrust_ratio_finite=thrust[1] / thrust[0],
        thrust_ratio_frozen=thrust[2] / thrust[0],
        bond_shifting=bond[0],
        bond_finite=bond[1],
        bond_frozen=bond[2],
        peak=peak,
        parcels=len(rows),
    )


def _shifting_efficiency(mod: Model, tab: nt.Isentrope, area_ratio: float) -> nt.Efficiency:
    """`near_term.efficiencies` needs a solution only for its labels; this is its integral."""
    rho0, u0 = float(tab.rho[0]), float(tab.energy[0])
    exits: list[nt.NozzleExit] = []
    for j in range(len(tab.rho)):
        ex = nt.nozzle_exit(tab, j, area_ratio)
        if ex is None:
            break
        exits.append(ex)
    rho_j = tab.rho[: len(exits)]
    ke = np.array([0.5 * ex.speed**2 for ex in exits])
    v = np.array([ex.speed for ex in exits])
    drho = -np.diff(rho_j)
    return nt.Efficiency(
        pairing="",
        thermo=mod.feed.thermo,
        branch=mod.branch.name,
        area_ratio=area_ratio,
        eta_peak=0.5 * exits[0].speed ** 2 / u0,
        eta_blowdown=float(np.sum(0.5 * (ke[1:] + ke[:-1]) * drho)) / (rho0 * u0),
        speed_peak=exits[0].speed,
        speed_impulse=float(np.sum(0.5 * (v[1:] + v[:-1]) * drho)) / rho0,
        unresolved_mass=float(rho_j[-1]) / rho0,
        exit_temp_peak=exits[0].temp,
        exit_pressure_peak=exits[0].pressure,
        exit_thermal_share=math.nan,
        exit_chemical_share=math.nan,
        exit_condensed_fraction=0.0,
        slug_ratio=math.nan,
    )


# ---- How fast carbon would have to associate -------------------------------------------------


@dataclass(frozen=True)
class RequiredCarbonRate:
    """The three-body coefficient a carbon channel needs to keep up to its half-release station.

    Bray's criterion, `Da = 1`, applied to carbon on the equilibrium expansion: free carbon falls
    from `carbon_free_throat` to `carbon_free_exit` in equilibrium, and the coefficient that keeps
    up at the halfway station returns about half the carbon store. It is a threshold, not a rate
    used for any verdict: the provenance rule (ADR-0050) forbids that.
    """

    area_ratio: float
    temp: float
    carbon_free_throat: float
    carbon_free_half: float
    carbon_free_exit: float
    tau_expansion: float
    #: Free atoms a carbon atom could associate with (C + H) [m^-3], and the third bodies.
    n_partner: float
    n_third: float
    coefficient: float
    #: Slack 1976's `C + C + Ar` at this temperature, and whether that is outside its window.
    slack: float
    slack_extrapolated: bool
    #: `H + H + M`'s effective coefficient at the same station, for scale.
    hydrogen: float

    @property
    def decades_short(self) -> float:
        """`log10(required / Slack)`: positive means the one measurement is too slow."""
        return math.log10(self.coefficient / self.slack)


def required_carbon_rate(
    sol: nt.ChamberSolution,
    tab: nt.Isentrope,
    throat_area: float,
    area_ratio: float,
    half_angle: float = nt.CONE_HALF_ANGLE,
) -> RequiredCarbonRate:
    """The Bray threshold for carbon, from the peak state on the equilibrium isentrope `tab`."""
    ex = nt.nozzle_exit(tab, 0, area_ratio)
    if ex is None:
        raise ValueError("isentrope does not reach the exit")
    h0 = float(tab.enthalpy[0])
    r_star = math.sqrt(throat_area / math.pi)
    tan_theta = math.tan(half_angle)
    stations = []
    for i in range(ex.i_star, ex.i_exit + 1):
        rho, temp = float(tab.rho[i]), float(tab.temp[i])
        u = math.sqrt(max(0.0, 2.0 * (h0 - float(tab.enthalpy[i]))))
        ar = max(1.0, ex.flux_star / (rho * u))
        comp = eos_methane.composition(rho, temp, sol.feed)
        stations.append((rho, temp, u, r_star * (math.sqrt(ar) - 1.0) / tan_theta, ar, comp))
    free = [s[5].carbon_free for s in stations]
    half = 0.5 * (free[0] + free[-1])
    k = next(n for n in range(1, len(free)) if free[n] <= half)
    (rho_a, _, u_a, x_a, _, _), (rho_b, temp, u_b, x_b, ar, comp) = stations[k - 1], stations[k]
    dt = (x_b - x_a) * 0.5 * (1.0 / u_a + 1.0 / u_b)
    tau = dt / math.log(rho_a / rho_b)
    n_partner = comp.n_c + comp.n_h
    n_third = comp.n_neutral_heavy
    return RequiredCarbonRate(
        area_ratio=ar,
        temp=temp,
        carbon_free_throat=free[0],
        carbon_free_half=free[k],
        carbon_free_exit=free[-1],
        tau_expansion=tau,
        n_partner=n_partner,
        n_third=n_third,
        coefficient=1.0 / (tau * n_partner * n_third),
        slack=rates.C_AR(temp),
        slack_extrapolated=rates.C_AR.extrapolated(temp),
        hydrogen=freeze.mixture_coefficient(temp, comp),
    )


# ---- S26: the cases, the sweep and the tables -------------------------------------------------

#: The ask's nozzles on the bulged chamber, plus the B6 row at A/A* 100: `(label, A*, A/A*)`.
BULGED_NOZZLES = (
    ("1.33 m^2, AR 100", 1.33, 100.0),
    ("1.33 m^2, AR 300 (22.5 m exit)", 1.33, 300.0),
    ("2.66 m^2, AR 150 (22.5 m exit)", 2.66, 150.0),
    ("2.66 m^2, AR 300 (31.9 m exit)", 2.66, 300.0),
)
#: The sphere for calibration: the paper's two near-term throats, at A/A* 300.
SPHERE_NOZZLES = (
    ("sphere 0.198 m^2, AR 300", 0.198, 300.0),
    ("sphere 0.149 m^2, AR 300", 0.149, 300.0),
)
HALF_ANGLES_DEG = (10.0, 15.0, 20.0)
#: Half a decade either way: Baulch's stated uncertainty on `H + H + H2`, and wider than the
#: factor ~2.5 spread of the three shock tubes' atomic-H efficiencies.
RATE_SCALES = (10.0**-0.5, 1.0, 10.0**0.5)


@dataclass(frozen=True)
class Job:
    case: str  # "bulged" or "sphere"
    label: str
    model: str
    throat_area: float
    area_ratio: float
    half_angle_deg: float
    rate_scale: float


def jobs() -> list[Job]:
    """Every pulse S26 needs: the main rows, the angle rows and the rate rows."""
    out = []
    asked = BULGED_NOZZLES[1:]
    for m in MODELS:
        for label, a, ar in BULGED_NOZZLES:
            out.append(Job("bulged", label, m, a, ar, 15.0, 1.0))
        for label, a, ar in asked:
            out += [Job("bulged", label, m, a, ar, d, 1.0) for d in (10.0, 20.0)]
            out += [
                Job("bulged", label, m, a, ar, 15.0, s) for s in (RATE_SCALES[0], RATE_SCALES[2])
            ]
        for label, a, ar in SPHERE_NOZZLES:
            out.append(Job("sphere", label, m, a, ar, 15.0, 1.0))
    return out


def bulged_solution() -> nt.ChamberSolution:
    """The 212 m^3 bulged chamber's peak state (B6's: 108.4 kg of methane per 5 kg rod)."""
    return bd.solve(bd.methane_specs()[0])


def sphere_solution() -> nt.ChamberSolution:
    """The 20 m^3 near-term sphere: 68.8 kg at 7000 K, 2.5 kg rod."""
    return nt.solve_charge(nt.PAIRINGS[1], "janaf")


_WORKER: dict[str, tuple[nt.ChamberSolution, dict[str, tuple[Model, nt.Isentrope]]]] = {}


def _setup(case: str) -> tuple[nt.ChamberSolution, dict[str, tuple[Model, nt.Isentrope]]]:
    if case not in _WORKER:
        sol = bulged_solution() if case == "bulged" else sphere_solution()
        tabs = {}
        for m in MODELS:
            mod = model(m, sol)
            tabs[m] = (
                mod,
                nt.isentrope(mod.branch, sol.rho, sol.pairing.temp, rho_ratio_end=1e-8, steps=900),
            )
        _WORKER[case] = (sol, tabs)
    return _WORKER[case]


def run_job(job: Job) -> tuple[Job, Pulse]:
    _, tabs = _setup(job.case)
    mod, tab = tabs[job.model]
    return job, pulse(
        mod, tab, job.throat_area, job.area_ratio, math.radians(job.half_angle_deg), job.rate_scale
    )


@dataclass(frozen=True)
class Bray:
    """B6's sudden-freeze figure for the same nozzle, for the comparison."""

    label: str
    eta_blowdown: float
    eta_net: float
    freeze_area_ratio: float
    speed_impulse: float


def bray_rows() -> list[Bray]:
    out = []
    for case, nozzles in (("bulged", BULGED_NOZZLES), ("sphere", SPHERE_NOZZLES)):
        sol, tabs = _setup(case)
        _, tab = tabs["exchange_carbon"]
        for label, a, ar in nozzles:
            eff = _shifting_efficiency(tabs["exchange_carbon"][0], tab, ar)
            fc = dataclasses.replace(
                nt.freeze_check(sol, tab, a, ar), eta_blowdown=eff.eta_blowdown
            )
            out.append(
                Bray(label, eff.eta_blowdown, fc.eta_net, fc.freeze_area_ratio, eff.speed_impulse)
            )
    return out


def carbon_rows() -> list[tuple[str, RequiredCarbonRate]]:
    out = []
    for case, nozzles in (("bulged", BULGED_NOZZLES), ("sphere", SPHERE_NOZZLES)):
        sol, tabs = _setup(case)
        _, tab = tabs["exchange_carbon"]
        for label, a, ar in nozzles:
            out.append((label, required_carbon_rate(sol, tab, a, ar)))
    return out


def impulse(speed: float, throat_area: float) -> float:
    """Net impulse per 5 kg pulse [N s] at a mean exhaust speed (B6's formula, gate included)."""
    spec = dataclasses.replace(bd.methane_specs()[0], throat_area=throat_area)
    return bd.net_impulse(spec, bulged_solution(), speed)


def write(
    results: Sequence[tuple[Job, Pulse]],
    bray: Sequence[Bray],
    carbon: Sequence[tuple[str, RequiredCarbonRate]],
    out: Path = OUTPUT_DIR,
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "pulses.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["case", "nozzle", "model", "throat_area", "area_ratio", "half_angle_deg", "rate_scale",
             "eta_shifting", "eta_finite", "eta_frozen", "speed_shifting", "speed_finite",
             "speed_frozen", "thrust_ratio_finite", "thrust_ratio_frozen", "bond_shifting",
             "bond_finite", "bond_frozen", "peak_bond_finite", "peak_freeze_area_ratio",
             "peak_throat_damkohler", "peak_length", "peak_transit_time", "peak_exit_temp",
             "peak_hydrogen_free_exit", "parcels"]
        )  # fmt: skip
        for j, p in results:
            w.writerow(
                [j.case, j.label, j.model, j.throat_area, j.area_ratio, j.half_angle_deg,
                 f"{j.rate_scale:.4g}", f"{p.eta_shifting:.4f}", f"{p.eta_finite:.4f}",
                 f"{p.eta_frozen:.4f}", f"{p.speed_shifting:.1f}", f"{p.speed_finite:.1f}",
                 f"{p.speed_frozen:.1f}", f"{p.thrust_ratio_finite:.4f}",
                 f"{p.thrust_ratio_frozen:.4f}", f"{p.bond_shifting:.4f}", f"{p.bond_finite:.4f}",
                 f"{p.bond_frozen:.4f}", f"{p.peak.bond_returned:.4f}",
                 f"{p.peak.freeze_area_ratio:.1f}", f"{p.peak.throat_damkohler:.3g}",
                 f"{p.peak.length:.2f}", f"{p.peak.transit_time:.3e}", f"{p.peak.temp:.0f}",
                 f"{p.peak.hydrogen_free_exit:.3f}", p.parcels]
            )  # fmt: skip
    with (out / "bray.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["nozzle", "eta_blowdown", "eta_net", "freeze_area_ratio", "speed_impulse"])
        for b in bray:
            w.writerow([b.label, f"{b.eta_blowdown:.4f}", f"{b.eta_net:.4f}",
                        f"{b.freeze_area_ratio:.1f}", f"{b.speed_impulse:.1f}"])  # fmt: skip
    with (out / "carbon_threshold.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["nozzle", "area_ratio_half", "temp", "carbon_free_throat",
                    "carbon_free_half", "carbon_free_exit", "tau_expansion", "k_required",
                    "k_slack", "slack_extrapolated", "decades_short", "k_hydrogen"])  # fmt: skip
        for label, c in carbon:
            w.writerow([label, f"{c.area_ratio:.2f}", f"{c.temp:.0f}",
                        f"{c.carbon_free_throat:.3f}", f"{c.carbon_free_half:.3f}",
                        f"{c.carbon_free_exit:.3f}", f"{c.tau_expansion:.3e}",
                        f"{c.coefficient:.3e}", f"{c.slack:.3e}", c.slack_extrapolated,
                        f"{c.decades_short:.2f}", f"{c.hydrogen:.3e}"])  # fmt: skip


def report(
    results: Sequence[tuple[Job, Pulse]],
    bray: Sequence[Bray],
    carbon: Sequence[tuple[str, RequiredCarbonRate]],
) -> None:
    by = {
        (j.case, j.label, j.model, j.half_angle_deg, round(j.rate_scale, 4)): p for j, p in results
    }
    bray_by = {b.label: b for b in bray}
    print("\n## Main rows (15 deg cone, measured rate)")
    print("| nozzle | eta eq | Bray | exch. finite | FC shift | FC finite | frozen | "
          "bond eq / exch / FC / frozen | I eq / Bray / exch / FC / frozen [kN s] |")  # fmt: skip
    for case, nozzles in (("bulged", BULGED_NOZZLES), ("sphere", SPHERE_NOZZLES)):
        for label, a, _ in nozzles:
            ex = by[(case, label, "exchange_carbon", 15.0, 1.0)]
            fc = by[(case, label, "frozen_carbon", 15.0, 1.0)]
            b = bray_by[label]
            imp = ""
            if case == "bulged":
                v_bray = b.speed_impulse * math.sqrt(b.eta_net / b.eta_blowdown)
                vals = (
                    ex.speed_shifting,
                    v_bray,
                    ex.speed_shifting * ex.thrust_ratio_finite,
                    fc.speed_shifting * fc.thrust_ratio_finite,
                    fc.speed_shifting * fc.thrust_ratio_frozen,
                )
                imp = " / ".join(f"{impulse(v, a) / 1e3:.0f}" for v in vals)
            print(
                f"| {label} | {ex.eta_shifting:.3f} | {b.eta_net:.3f} | {ex.eta_finite:.3f} | "
                f"{fc.eta_shifting:.3f} | {fc.eta_finite:.3f} | {fc.eta_frozen:.3f} | "
                f"{ex.bond_shifting:.2f} / {ex.bond_finite:.2f} / {fc.bond_finite:.2f} / "
                f"{fc.bond_frozen:.2f} | {imp} |"
            )
    print("\n## Freeze points and clocks (peak parcel, 15 deg, measured rate)")
    print(
        "| nozzle | model | Da at A/A* 1.05 | Bray AR | finite: first Da<1 AR | length | transit "
        "| exit T | H free at exit |"
    )
    for (_, label, m, ang, sc), p in by.items():
        if ang == 15.0 and sc == 1.0:
            br = bray_by[label].freeze_area_ratio if m == "exchange_carbon" else math.nan
            print(f"| {label} | {m} | {p.peak.throat_damkohler:.3g} | {br:.1f} | "
                  f"{p.peak.freeze_area_ratio:.1f} | {p.peak.length:.1f} m | "
                  f"{1e3 * p.peak.transit_time:.2f} ms | {p.peak.temp:.0f} K | "
                  f"{p.peak.hydrogen_free_exit:.3f} |")  # fmt: skip
    print("\n## Sensitivity: wall angle and rate (thrust ratio to shifting; bond returned)")
    print("| nozzle | model | 10 deg | 15 deg | 20 deg | k/3.16 | k x3.16 |")
    for label, _, _ in BULGED_NOZZLES[1:]:
        for m in MODELS:
            cells = []
            for ang, sc in ((10.0, 1.0), (15.0, 1.0), (20.0, 1.0), (15.0, 0.3162), (15.0, 3.1623)):
                p = by[("bulged", label, m, ang, sc)]
                cells.append(f"{p.thrust_ratio_finite:.4f} ({p.bond_finite:.3f})")
            print(f"| {label} | {m} | " + " | ".join(cells) + " |")
    print("\n## Carbon: the three-body coefficient needed to return half the carbon store")
    print(
        "| nozzle | AR at half | T | free C throat/half/exit | k required | Slack C+C+Ar "
        "| decades short | H+H+M k |"
    )
    for label, c in carbon:
        print(f"| {label} | {c.area_ratio:.1f} | {c.temp:.0f} K | {c.carbon_free_throat:.2f} / "
              f"{c.carbon_free_half:.2f} / {c.carbon_free_exit:.2f} | {c.coefficient:.2e} | "
              f"{c.slack:.2e}{' (extrap.)' if c.slack_extrapolated else ''} | "
              f"{c.decades_short:+.1f} | {c.hydrogen:.2e} |")  # fmt: skip


def main() -> None:
    from concurrent.futures import ProcessPoolExecutor

    todo = jobs()
    # Slowest first: the frozen-carbon pulses take the bisection fallback below 1000 K.
    todo.sort(key=lambda j: (j.model != "frozen_carbon", j.case))
    with ProcessPoolExecutor() as pool:
        results = list(pool.map(run_job, todo))
    bray, carbon = bray_rows(), carbon_rows()
    write(results, bray, carbon)
    report(results, bray, carbon)


if __name__ == "__main__":
    main()
