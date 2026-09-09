"""Pure-hydrogen equilibrium, and the Project 242 tension N10 asks the solver to settle.

Ask N10 sets out a contradiction and asks which side is right:

> Project 242 writes that molecular recombination is *not* fast, but its carried baseline of
> 2700 s needs 351 MJ/kg of stagnation enthalpy where sensible-only hydrogen at 10,000 K gives
> 206. That implies about 67% of the dissociation energy returning at a few bar... Those two
> statements are in tension and the solver should say which is right.

**That part of N10 needs no rate coefficient.** It is an energy budget: what stagnation enthalpy
does equilibrium hydrogen actually hold at Rubbia's stated chamber, and what exhaust speed does
each of the two limiting expansions -- fully frozen and fully equilibrium -- deliver from it?
Those two bracket every finite-rate answer, so if 2700 s sits outside the frozen limit then
recombination *must* be contributing, whatever the prose says.

The species set is `H2`, `H`, `H+`, `e-`: two unknowns (`ln n_H`, `ln n_e`), `H2` by mass action
off `eos_water.H2`, `H+` by Saha, closed by nuclei conservation and charge neutrality. It reuses
`eos_water`'s diatomic constants rather than restating them, so hydrogen cannot drift between
the two modules.

This also supplies ask N10 item 3's pure-hydrogen rung of ADR-0016's ladder.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from puffsat import eos_water
from puffsat.eos_water import H2, K_B, ln_k_diatomic, ln_k_saha, sound_speed_fd

M_H2 = H2.mass
#: Dissociation energy of one kilogram of H2 [J/kg] -- 214.3 MJ/kg. The store in question.
FULL_DISSOCIATION_ENERGY = H2.d0 / M_H2

#: Project 242's operating point, as ADR-0016 reads it: hydrogen at 10 kK, "a few bar", and a
#: 9500 K radiative ceiling. Its carried baseline is 2700 s at 3200 N.
RUBBIA_TEMPERATURE = 10000.0
RUBBIA_PRESSURE_RANGE = (1.0e5, 1.0e6)  # 1 to 10 bar: "a few bar", bracketed rather than picked
RUBBIA_ISP = 2700.0
G0 = 9.80665


@dataclass(frozen=True)
class HydrogenState:
    """Equilibrium pure hydrogen at one `(rho, T)`."""

    rho: float
    temp: float
    n_h: float
    n_h2: float
    n_hp: float
    n_e: float
    pressure: float
    energy: float
    #: Fraction of hydrogen nuclei that are free -- atomic H plus H+.
    dissociated: float

    @property
    def enthalpy(self) -> float:
        """Specific stagnation enthalpy `h = e + p/rho` [J/kg].

        **This is the quantity the ask compares against 351 MJ/kg**, and it is enthalpy rather
        than internal energy because what a nozzle converts to jet kinetic energy is `h`, not `e`
        -- the flow does `p/rho` of pushing work on itself on the way out.
        """
        return self.energy + self.pressure / self.rho

    @property
    def chemical_energy(self) -> float:
        """Dissociation energy currently held [J/kg] -- what recombination has left to return."""
        return (self.n_h + self.n_hp) * 0.5 * H2.d0 / self.rho


def _residual(x: np.ndarray, ln_kd: float, ln_kh: float, n_nuclei: float) -> np.ndarray:
    """Nuclei conservation and charge neutrality, both in log form."""
    ln_nh, ln_ne = float(x[0]), float(x[1])
    ln_nh2 = 2.0 * ln_nh - ln_kd
    ln_hp = ln_kh + ln_nh - ln_ne
    ln_total = float(np.logaddexp(np.logaddexp(ln_nh, ln_hp), math.log(2.0) + ln_nh2))
    return np.array([ln_total - math.log(n_nuclei), ln_ne - ln_hp])


def state(rho: float, temp: float) -> HydrogenState:
    """Equilibrium hydrogen at `(rho [kg/m^3], temp [K])`."""
    ln_kd = ln_k_diatomic(H2, temp)  # n_H^2 / n_H2
    ln_kh = ln_k_saha(eos_water.IP_H, 1.0, eos_water.G_H, temp)
    n_nuclei = 2.0 * rho / M_H2

    # Start from the binary H2 <=> 2H solution, which is exact when ionisation is negligible and
    # a good start when it is not.
    k_d = math.exp(min(ln_kd, 700.0))
    n_h0 = (-k_d + math.sqrt(k_d**2 + 8.0 * k_d * n_nuclei)) / 4.0
    x = np.array([math.log(max(n_h0, 1e-30)), 0.5 * (ln_kh + math.log(max(n_h0, 1e-30)))])

    r = _residual(x, ln_kd, ln_kh, n_nuclei)
    for _ in range(200):
        if float(np.max(np.abs(r))) < 1e-12:
            break
        jac = np.empty((2, 2))
        for k in range(2):
            xp = x.copy()
            xp[k] += 1e-6
            jac[:, k] = (_residual(xp, ln_kd, ln_kh, n_nuclei) - r) / 1e-6
        delta = np.linalg.solve(jac, -r)
        step = float(np.max(np.abs(delta)))
        if step > 2.0:
            delta *= 2.0 / step
        x = x + delta
        r = _residual(x, ln_kd, ln_kh, n_nuclei)
    else:  # pragma: no cover - two-unknown Newton from a near-exact start
        raise RuntimeError(f"hydrogen solve did not converge at rho={rho}, T={temp}")

    n_h = math.exp(float(x[0]))
    n_e = math.exp(float(x[1]))
    n_h2 = math.exp(2.0 * float(x[0]) - ln_kd)
    n_hp = n_e
    n_total = n_h + n_h2 + n_hp + n_e
    p = n_total * K_B * temp
    e = (
        1.5 * K_B * temp * (n_h + n_hp + n_e)
        + n_h2 * (2.5 * K_B * temp + eos_water._e_vib_mode(H2.theta_vib, temp))
        + (n_h + n_hp) * 0.5 * H2.d0
        + n_hp * eos_water.IP_H
    ) / rho
    return HydrogenState(
        rho=rho,
        temp=temp,
        n_h=n_h,
        n_h2=n_h2,
        n_hp=n_hp,
        n_e=n_e,
        pressure=p,
        energy=e,
        dissociated=(n_h + n_hp) / n_nuclei,
    )


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """`(p, e)` at `(rho, temp)` -- the `Eos` signature `expansion` and `sound_speed_fd` want."""
    s = state(rho, temp)
    return s.pressure, s.energy


def sound_speed(rho: float, temp: float) -> float:
    """Equilibrium adiabatic sound speed [m/s]."""
    return sound_speed_fd(pressure_energy, rho, temp)


def at_pressure(pressure: float, temp: float) -> HydrogenState:
    """Equilibrium hydrogen at a stated `(p, T)` -- Rubbia's chamber is quoted that way."""
    lo, hi = 1e-6, 1e3
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if state(mid, temp).pressure < pressure:
            lo = mid
        else:
            hi = mid
        if hi / lo < 1.0 + 1e-12:
            break
    return state(math.sqrt(lo * hi), temp)


@dataclass(frozen=True)
class IspBracket:
    """The two limiting expansions from one chamber state, and where a claimed Isp falls.

    Both limits assume a perfect nozzle -- everything converted, nothing left in the exhaust --
    so they are ceilings on their own branch, not predictions. That is what makes the comparison
    decisive in one direction only: an Isp **above the frozen ceiling** cannot be reached without
    recombination, whatever the prose says. An Isp below it proves nothing either way.
    """

    temp: float
    pressure: float
    dissociated: float
    enthalpy: float
    #: Enthalpy less the dissociation store still held: what a fully frozen expansion can convert.
    frozen_enthalpy: float
    isp_equilibrium: float
    isp_frozen: float

    def implied_return(self, isp: float) -> float:
        """Fraction of the held dissociation store a claimed `isp` requires to come back.

        Solves `isp` against `sqrt(2 (h_frozen + f * h_chem))/g0`. Negative means the claim sits
        below the frozen ceiling and needs nothing; above 1 means it does not close even with
        full recombination.
        """
        needed = 0.5 * (isp * G0) ** 2
        chem = self.enthalpy - self.frozen_enthalpy
        return (needed - self.frozen_enthalpy) / chem if chem > 0.0 else math.inf


def isp_bracket(pressure: float, temp: float) -> IspBracket:
    """Frozen and equilibrium Isp ceilings for hydrogen at a stated chamber `(p, T)`."""
    s = at_pressure(pressure, temp)
    frozen = s.enthalpy - s.chemical_energy
    return IspBracket(
        temp=temp,
        pressure=s.pressure,
        dissociated=s.dissociated,
        enthalpy=s.enthalpy,
        frozen_enthalpy=frozen,
        isp_equilibrium=math.sqrt(2.0 * s.enthalpy) / G0,
        isp_frozen=math.sqrt(2.0 * frozen) / G0,
    )


def main() -> None:
    """Settle the Project 242 tension, and report the pure-hydrogen rung of ADR-0016's ladder."""
    print("N10 -- Project 242: does 2700 s need recombination?")
    print(
        f"{'p [bar]':>8} {'T':>7} {'diss':>7} {'h [MJ/kg]':>10} {'h_frozen':>9} "
        f"{'Isp_eq':>8} {'Isp_fr':>8} {'return needed':>14}"
    )
    for pressure in (1.0e5, 3.0e5, 1.0e6, 1.0e7):
        b = isp_bracket(pressure, RUBBIA_TEMPERATURE)
        print(
            f"{b.pressure / 1e5:8.2f} {b.temp:7.0f} {b.dissociated:7.4f} "
            f"{b.enthalpy / 1e6:10.1f} {b.frozen_enthalpy / 1e6:9.1f} "
            f"{b.isp_equilibrium:8.0f} {b.isp_frozen:8.0f} "
            f"{b.implied_return(RUBBIA_ISP):14.3f}"
        )

    print()
    print("The same bracket at the walled nozzle's own density (ADR-0016's hydrogen rung)")
    for rho in (0.5, 1.0, 2.0):
        s = state(rho, RUBBIA_TEMPERATURE)
        frozen = s.enthalpy - s.chemical_energy
        print(
            f"rho {rho:5.2f} kg/m3  p {s.pressure / 1e5:8.1f} bar  diss {s.dissociated:.4f}  "
            f"h {s.enthalpy / 1e6:6.1f} MJ/kg  Isp_eq {math.sqrt(2 * s.enthalpy) / G0:5.0f} s  "
            f"Isp_fr {math.sqrt(2 * frozen) / G0:5.0f} s"
        )


if __name__ == "__main__":
    main()
