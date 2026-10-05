"""Atomic Saha EOS for the argon arm's premixed spray mixture (ADR-0055).

The mixture is the PuffSat's water and the sprayed argon at injection ratio `k`, carried as atoms:
H, O and Ar, every ionization stage of each, and electrons. A shared electron density closes the
Saha ladders through charge neutrality. Ground-term degeneracies, NIST ionization energies and no
continuum lowering follow the convention of `puffsat.eos_water`.

The water fraction carries no molecules, so its bond energy (50.9 MJ/kg of water) never returns.
The caller charges that energy before the gas reaches this table. That is the floor case the
child repo uses when it charges water as never re-forming.

Energy reference: neutral ground-state atoms at T -> 0 have e = 0, so every grid value is positive,
as the Rust loader requires.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from puffsat import eos_water as ew
from puffsat.eos_water import EV, K_B, M_H2O, Vec

M_AR = 39.948 * ew.AMU

# Argon ionization energies Ar I->II .. Ar XVIII->XIX [eV], NIST Atomic Spectra Database
# ionization-energy query for Ar (physics.nist.gov/cgi-bin/ASD/ie.pl), retrieved 2026-10-05.
IP_AR_EV = (
    15.7596119,
    27.62967,
    40.735,
    59.58,
    74.84,
    91.290,
    124.41,
    143.4567,
    422.60,
    479.76,
    540.4,
    619.0,
    685.5,
    755.13,
    855.5,
    918.375,
    4120.66559,
    4426.22407,
)
# Ground-term degeneracies (2S+1)(2L+1) of Ar I .. Ar XIX, from the ground levels in the same pull:
# 1S, 2P, 3P, 4S, 3P, 2P, 1S, 2S, 1S, 2P, 3P, 4S, 3P, 2P, 1S, 2S, 1S, 2S, bare.
G_AR_LADDER = (
    1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0
)  # fmt: skip

# Element ladders: (atomic mass, ionization energies [J], degeneracies of stages 0..Z).
_LADDERS = {
    "H": (ew.M_H, (ew.IP_H,), (ew.G_H, ew.G_HP)),
    "O": (ew.M_O, ew.IP_O_LADDER, ew.G_O_LADDER),
    "Ar": (M_AR, tuple(ip * EV for ip in IP_AR_EV), G_AR_LADDER),
}


@dataclass(frozen=True)
class SprayMixture:
    """Element numbers per kilogram of mixture at injection ratio `k` (argon kg per water kg)."""

    k: float

    @property
    def per_kg(self) -> dict[str, float]:
        water = 1.0 / (1.0 + self.k) / M_H2O
        return {"H": 2.0 * water, "O": water, "Ar": self.k / (1.0 + self.k) / M_AR}


def _stage_log_fractions(name: str, temp: float, ln_ne: float) -> Vec:
    """ln of each stage's share of its element at fixed `n_e` (normalized in log space)."""
    _, ips, gs = _LADDERS[name]
    ln_ratio = np.array(
        [ew.ln_k_saha(ips[j], gs[j + 1], gs[j], temp) - ln_ne for j in range(len(ips))]
    )
    ln_rel = np.concatenate([[0.0], np.cumsum(ln_ratio)])
    top = float(np.max(ln_rel))
    out: Vec = ln_rel - (top + np.log(np.sum(np.exp(ln_rel - top))))
    return out


def _cumulative_ip(name: str) -> Vec:
    _, ips, _ = _LADDERS[name]
    return np.concatenate([[0.0], np.cumsum(np.array(ips))])


@dataclass(frozen=True)
class SprayState:
    """Equilibrium state at one `(rho, T)`: electron density and per-element stage fractions."""

    n_e: float
    fractions: dict[str, Vec]
    n_element: dict[str, float]

    def mean_charge(self, name: str) -> float:
        f = self.fractions[name]
        return float(np.dot(np.arange(len(f)), f))


def state(mix: SprayMixture, rho: float, temp: float) -> SprayState:
    """Solve charge neutrality `n_e = sum n_el Zbar_el(n_e)` by bisection in `ln n_e`.

    The right side falls as `n_e` rises, so the root is unique and bisection cannot miss it.
    """
    n_el = {name: rho * n for name, n in mix.per_kg.items()}
    n_max = sum(n_el[name] * (len(_LADDERS[name][1])) for name in n_el)
    lo, hi = np.log(n_max) - 700.0, np.log(n_max)

    def excess(ln_ne: float) -> float:
        supply = 0.0
        for name, n in n_el.items():
            f = np.exp(_stage_log_fractions(name, temp, ln_ne))
            supply += n * float(np.dot(np.arange(len(f)), f))
        return float(np.log(max(supply, 1e-300)) - ln_ne)

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if excess(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    ln_ne = 0.5 * (lo + hi)
    fractions = {name: np.exp(_stage_log_fractions(name, temp, ln_ne)) for name in n_el}
    return SprayState(n_e=float(np.exp(ln_ne)), fractions=fractions, n_element=n_el)


def pressure_energy(mix: SprayMixture, rho: float, temp: float) -> tuple[float, float]:
    """`(p [Pa], e [J/kg])`: ideal translational pressure and energy plus ionization energy."""
    s = state(mix, rho, temp)
    n_heavy = sum(s.n_element.values())
    p = (n_heavy + s.n_e) * K_B * temp
    e_ion = sum(
        s.n_element[name] * float(np.dot(s.fractions[name], _cumulative_ip(name)))
        for name in s.n_element
    )
    e = (1.5 * p + e_ion) / rho
    return float(p), float(e)


def sound_speed(mix: SprayMixture, rho: float, temp: float) -> float:
    """Equilibrium adiabatic sound speed [m/s], by the water module's finite-difference rule."""
    return ew.sound_speed_fd(lambda r, t: pressure_energy(mix, r, t), rho, temp)


def eos_grid(mix: SprayMixture, rho_grid: Vec, t_grid: Vec) -> tuple[Vec, Vec, Vec]:
    """`(p, e, c_s)` on the `(rho, T)` grid, shaped `(n_rho, n_T)`, row-major over `(rho, T)`."""
    shape = (len(rho_grid), len(t_grid))
    p, e, cs = np.empty(shape), np.empty(shape), np.empty(shape)
    for i, rho in enumerate(rho_grid):
        for j, temp in enumerate(t_grid):
            p[i, j], e[i, j] = pressure_energy(mix, float(rho), float(temp))
            cs[i, j] = sound_speed(mix, float(rho), float(temp))
    return p, e, cs
