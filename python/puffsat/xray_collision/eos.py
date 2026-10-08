"""Single-element Saha EOS plus TOPS gray opacity, tabulated once on a (rho, T) grid.

Ideal-plasma Saha over every ionization stage with ground-term partition functions:

    n_{z+1} n_e / n_z = 2 (g_{z+1}/g_z) (2 pi m_e k T / h^2)^{3/2} exp(-I_z / kT)
    p = (1 + zbar) n k T,     e = [3/2 (1 + zbar) k T + sum_z f_z sum_{z'<z} I_z'] / m_a

Energy is measured from neutral ground-state atoms. Liquid/solid binding (< 0.1 eV/atom) is
dropped. With `lowering=True` every I_z is reduced by the ion-sphere amount
(z+1) e^2 / (4 pi eps0 a_i), with a_i the Wigner-Seitz radius. That brackets the
under-ionization ideal Saha gives at condensed density. The lowered energy is charged the same
way, which is not thermodynamically consistent and is good only as a bracket.

The table is held in log rho and log T. `temp` inverts e(T) on an interpolated row, and `state`
interpolates bilinearly (log p, zbar, log kappa).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from puffsat import tops
from puffsat.eos_water import Vec
from puffsat.xray_collision.materials import Material

K_B = 1.380649e-23
EV = 1.602176634e-19
AMU = 1.66053906660e-27
M_E = 9.1093837015e-31
H_PLANCK = 6.62607015e-34
E2_4PIEPS0_EV_M = 1.43996448e-9  # e^2 / (4 pi eps0) [eV m]

RHO_RANGE = (1.0e-12, 1.0e5)  # [kg/m^3]
T_RANGE_EV = (0.25, 400.0)
N_RHO = 141
N_T = 221
N_BISECT = 60


def ionization_energies(mat: Material, n_nuc: float, lowering: bool) -> Vec:
    """I_z [eV], optionally ion-sphere lowered and floored at zero."""
    ip = np.asarray(mat.ip_ev, dtype=float)
    if not lowering:
        return ip
    a_i = (3.0 / (4.0 * math.pi * n_nuc)) ** (1.0 / 3.0)
    shift = (np.arange(ip.size) + 1.0) * E2_4PIEPS0_EV_M / a_i
    return np.asarray(np.maximum(ip - shift, 0.0))


def saha(mat: Material, rho: float, temps: Vec, lowering: bool = False) -> tuple[Vec, Vec, Vec]:
    """`(p [Pa], e [J/kg], zbar)` at density `rho` for each temperature in `temps` [K].

    Charge neutrality n_e = zbar n is solved by bisection on ln n_e, vectorized over T."""
    m_a = mat.mass_amu * AMU
    n = rho / m_a
    ip = ionization_energies(mat, n, lowering)
    g = np.asarray(mat.g_ground, dtype=float)
    zs = np.arange(g.size, dtype=float)
    kt = K_B * temps
    ln_lam = 1.5 * np.log(2.0 * math.pi * M_E * kt / H_PLANCK**2)
    # ln[(n_{z+1}/n_z) n_e], shape (n_T, Z)
    ln_s = (
        np.log(2.0 * g[1:] / g[:-1])[None, :] + ln_lam[:, None] - (ip * EV)[None, :] / kt[:, None]
    )

    def fractions(ln_ne: Vec) -> Vec:
        ln_f = np.concatenate(
            [np.zeros((temps.size, 1)), np.cumsum(ln_s - ln_ne[:, None], axis=1)], axis=1
        )
        f = np.exp(ln_f - ln_f.max(axis=1, keepdims=True))
        return np.asarray(f / f.sum(axis=1, keepdims=True))

    lo = np.full(temps.size, math.log(n * 1e-30))
    hi = np.full(temps.size, math.log(n * ip.size))
    for _ in range(N_BISECT):
        mid = 0.5 * (lo + hi)
        too_few = fractions(mid) @ zs * n > np.exp(mid)
        lo = np.where(too_few, mid, lo)
        hi = np.where(too_few, hi, mid)
    f = fractions(0.5 * (lo + hi))
    zbar = f @ zs
    cum_ip = np.concatenate([[0.0], np.cumsum(ip)])
    p = (1.0 + zbar) * n * kt
    e = (1.5 * (1.0 + zbar) * kt + (f @ cum_ip) * EV) / m_a
    return p, e, zbar


@dataclass(frozen=True)
class EosTable:
    """Tabulated EOS and opacity of one material; arrays are `[i_rho, j_T]`."""

    material: Material
    ln_rho: Vec
    ln_t: Vec
    ln_p: Vec
    e: Vec
    zbar: Vec
    ln_kr: Vec
    ln_kp: Vec

    def _rho_weight(self, rho: float) -> tuple[int, float]:
        x = min(max(math.log(rho), float(self.ln_rho[0])), float(self.ln_rho[-1]))
        i = int(np.clip(np.searchsorted(self.ln_rho, x) - 1, 0, self.ln_rho.size - 2))
        return i, (x - float(self.ln_rho[i])) / float(self.ln_rho[i + 1] - self.ln_rho[i])

    def temp(self, rho: float, e: float) -> float:
        """T [K] such that e(rho, T) = e, clamped to the table's T range."""
        i, a = self._rho_weight(rho)
        row = (1.0 - a) * self.e[i] + a * self.e[i + 1]
        return math.exp(float(np.interp(e, row, self.ln_t)))

    def state(self, rho: float, temp: float) -> tuple[float, float, float, float]:
        """`(p [Pa], zbar, kappa_R, kappa_P [m^2/kg])` at (rho, T)."""
        i, a = self._rho_weight(rho)
        y = min(max(math.log(temp), float(self.ln_t[0])), float(self.ln_t[-1]))
        j = int(np.clip(np.searchsorted(self.ln_t, y) - 1, 0, self.ln_t.size - 2))
        b = (y - float(self.ln_t[j])) / float(self.ln_t[j + 1] - self.ln_t[j])

        def bil(t: Vec) -> float:
            return float(
                (1 - a) * (1 - b) * t[i, j]
                + a * (1 - b) * t[i + 1, j]
                + (1 - a) * b * t[i, j + 1]
                + a * b * t[i + 1, j + 1]
            )

        return (
            math.exp(bil(self.ln_p)),
            bil(self.zbar),
            math.exp(bil(self.ln_kr)),
            math.exp(bil(self.ln_kp)),
        )


@lru_cache(maxsize=8)
def build_table(mat: Material, lowering: bool = False) -> EosTable:
    """Tabulate `mat` on the module grid. TOPS means are held at the pull's edges."""
    rho_grid = np.geomspace(*RHO_RANGE, N_RHO)
    t_grid = np.geomspace(T_RANGE_EV[0], T_RANGE_EV[1], N_T) * EV / K_B
    p = np.empty((N_RHO, N_T))
    e = np.empty((N_RHO, N_T))
    z = np.empty((N_RHO, N_T))
    for i, rho in enumerate(rho_grid):
        p[i], e[i], z[i] = saha(mat, float(rho), t_grid, lowering)
    kr, kp = tops.resample(tops.load_tops_gray(mat.tops_pull), rho_grid, t_grid)
    return EosTable(mat, np.log(rho_grid), np.log(t_grid), np.log(p), e, z, np.log(kr), np.log(kp))
