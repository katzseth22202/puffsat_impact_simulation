"""Chemical-equilibrium equation of state for water, cold vapor to plasma (B5c-1).

The high-v PuffSat package needs an *equilibrium* water EOS (design §7): at 16 km/s the stagnated
gas is a 10-40 kK plasma, and the dominant control on the restitution `e_eff` is the
dissociation + ionization energy buffering (design §4) -- the specific heat sink that decides how
much impact kinetic energy disappears into internal modes instead of staying available to rebound.
CoolProp only reaches ~1300 K, so this is built from first principles (statistical mechanics), with
no external data. See `puffsat_impact_sim_design.md` §3-4, §7, and Zel'dovich & Raizer, *Physics of
Shock Waves and High-Temperature Hydrodynamic Phenomena*, Ch. III.

Model (reduced species set): H2O, H, O, H+, O+, e- in chemical + ionization equilibrium at each
`(rho, T)`. Dissociation H2O <=> 2H + O by law of mass action; ionization H <=> H+ + e- and
O <=> O+ + e- by the Saha equation; closed by H:O = 2:1 element conservation and charge neutrality.

Simplification (first-pass, documented): the molecular intermediates OH, H2, O2 are omitted, so the
single effective reaction H2O <=> 2H + O carries all the dissociation. This reshapes the
~2000-6000 K *transition* but preserves both endpoints (cold molecular vapor; hot ionized
plasma) and the total energy invested -- what sets `e_eff`. Adding the intermediates is the natural
refinement if the transition region is ever load-bearing.

Energy reference: bound molecular H2O at T -> 0 has e = 0, so dissociation and ionization only *add*
energy. Every `e` on the grid is therefore strictly positive -- required by the Rust table loader,
which interpolates `ln e` (ADR-0007).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Vec = NDArray[np.float64]

# --- Physical constants (SI; CODATA) ---
K_B = 1.380649e-23  # Boltzmann constant [J/K]
H_PLANCK = 6.62607015e-34  # Planck constant [J s]
M_E = 9.1093837015e-31  # electron mass [kg]
AMU = 1.66053906660e-27  # atomic mass unit [kg]
EV = 1.602176634e-19  # electronvolt [J]
N_A = 6.02214076e23  # Avogadro [1/mol]
C_LIGHT = 2.99792458e8  # speed of light [m/s]

# --- Species masses [kg] ---
M_H = 1.008 * AMU
M_O = 15.999 * AMU
M_H2O = 18.015 * AMU

# --- Ground-state electronic degeneracies (statistical weights) ---
G_H = 2.0  # H 2S_1/2
G_O = 9.0  # O 3P (sum 2J+1 over the ground term)
G_HP = 1.0  # H+ (bare proton)
G_OP = 4.0  # O+ 4S_3/2
G_E = 2.0  # electron spin
G_H2O = 1.0  # H2O ground electronic state

# --- Energetics ---
IP_H = 13.598 * EV  # H first ionization potential [J]
IP_O = 13.618 * EV  # O first ionization potential [J]
# Atomization energy of H2O -> 2H + O at 0 K: 2*dHf(H) + dHf(O) - dHf(H2O,g)
#   = 2(216.0) + 246.8 - (-238.9) kJ/mol = 917.7 kJ/mol  (NIST-JANAF 0 K heats of formation)
D_AT = 917.7e3 / N_A  # [J per molecule]

# --- H2O internal structure ---
# Rotational constants [cm^-1] (nonlinear asymmetric top); symmetry number sigma = 2.
_ROT_CONST_CM = (27.877, 14.512, 9.285)
_SIGMA = 2.0
# Vibrational fundamentals [cm^-1]: bend nu2, symmetric stretch nu1, antisymmetric stretch nu3.
_VIB_CM = (1595.0, 3657.0, 3756.0)
# cm^-1 -> characteristic temperature [K]: Theta = h c (100 * nu~) / k_B.
_CM_TO_K = H_PLANCK * C_LIGHT * 100.0 / K_B
_THETA_ROT = tuple(_CM_TO_K * v for v in _ROT_CONST_CM)
_THETA_VIB = tuple(_CM_TO_K * v for v in _VIB_CM)

_EXP_CAP = 700.0  # clamp exp() arguments to avoid float overflow during the Newton iteration


@dataclass(frozen=True)
class Composition:
    """Equilibrium number densities [m^-3] at a single `(rho, T)`."""

    n_h2o: float
    n_h: float
    n_o: float
    n_hp: float
    n_op: float
    n_e: float


def _ln_trans(mass: float, temp: float) -> float:
    """ln of the translational density of states `(2 pi m k T / h^2)^{3/2}` [ln m^-3]."""
    return float(1.5 * np.log(2.0 * np.pi * mass * K_B * temp / H_PLANCK**2))


def _z_rot(temp: float) -> float:
    """Classical rigid-rotor partition function of H2O (nonlinear), valid for T >> Theta_rot."""
    prod = _THETA_ROT[0] * _THETA_ROT[1] * _THETA_ROT[2]
    return float(np.sqrt(np.pi) / _SIGMA * np.sqrt(temp**3 / prod))


def _z_vib(temp: float) -> float:
    """Harmonic-oscillator vibrational partition function of H2O (zero of energy at the well)."""
    z = 1.0
    for theta in _THETA_VIB:
        z /= 1.0 - np.exp(-theta / temp)
    return float(z)


def _e_vib(temp: float) -> float:
    """Mean vibrational energy of one H2O molecule [J] (harmonic, excludes zero-point)."""
    e = 0.0
    for theta in _THETA_VIB:
        e += theta * K_B / (np.exp(theta / temp) - 1.0)
    return float(e)


def ln_k_saha(ip: float, g_ion: float, g_neutral: float, temp: float) -> float:
    """ln of the Saha constant `K = n_ion n_e / n_neutral` [ln m^-3] for one ionization stage.

    `K = (g_ion g_e / g_neutral) (2 pi m_e k T / h^2)^{3/2} exp(-IP / k T)` -- the ion/neutral
    translational factors cancel (equal masses), leaving the electron de Broglie term.
    """
    return float(np.log(g_ion * G_E / g_neutral) + _ln_trans(M_E, temp) - ip / (K_B * temp))


def ln_k_dissoc(temp: float) -> float:
    """ln of the dissociation constant `K_d = n_H^2 n_O / n_H2O` [ln m^-6] for H2O <=> 2H + O."""
    ln_h = _ln_trans(M_H, temp) + np.log(G_H)
    ln_o = _ln_trans(M_O, temp) + np.log(G_O)
    ln_h2o = _ln_trans(M_H2O, temp) + np.log(G_H2O * _z_rot(temp) * _z_vib(temp))
    return float(2.0 * ln_h + ln_o - ln_h2o - D_AT / (K_B * temp))


def _densities(x: Vec, ln_kd: float, ln_kh: float, ln_ko: float) -> tuple[float, ...]:
    """Reconstruct `(n_H2O, n_H, n_O, n_Hp, n_Op, n_e)` from log-densities `x = [ln n_H, ln n_O,
    ln n_e]` via mass action (n_H2O), Saha (n_Hp, n_Op). Kept in log form so the cold limit, where
    the equilibrium constants are ~1e-300, never underflows."""
    ln_nh, ln_no, ln_ne = float(x[0]), float(x[1]), float(x[2])
    n_h = float(np.exp(min(ln_nh, _EXP_CAP)))
    n_o = float(np.exp(min(ln_no, _EXP_CAP)))
    n_e = float(np.exp(min(ln_ne, _EXP_CAP)))
    n_h2o = float(np.exp(min(2.0 * ln_nh + ln_no - ln_kd, _EXP_CAP)))
    n_hp = float(np.exp(min(ln_kh + ln_nh - ln_ne, _EXP_CAP)))
    n_op = float(np.exp(min(ln_ko + ln_no - ln_ne, _EXP_CAP)))
    return n_h2o, n_h, n_o, n_hp, n_op, n_e


def _residual(x: Vec, ln_kd: float, ln_kh: float, ln_ko: float, n_f: float) -> Vec:
    """Equilibrium residual at log-densities `x = [ln n_H, ln n_O, ln n_e]`.

    Three independent, well-conditioned constraints (the naive H+O number pair degenerates when H2O
    dominates -- both then merely say `n_H2O = n_f`):
    - **O nuclei** conservation (pins `n_H2O` ~ `n_f` in the cold limit),
    - **H:O = 2:1** stoichiometry as a *ratio* (pins the atomic split independent of `n_H2O`),
    - **charge neutrality** `n_e^2 = K_H n_H + K_O n_O`, written in log form (well-defined and with
      slope 2 in `ln n_e` even when ionization is negligible).
    """
    ln_nh, ln_no, ln_ne = float(x[0]), float(x[1]), float(x[2])
    n_h2o, n_h, n_o, n_hp, n_op, _ = _densities(x, ln_kd, ln_kh, ln_ko)
    return np.array(
        [
            (n_h2o + n_o + n_op) / n_f - 1.0,  # O nuclei
            (n_h + n_hp) / (n_o + n_op) - 2.0,  # H:O stoichiometry
            2.0 * ln_ne - np.logaddexp(ln_kh + ln_nh, ln_ko + ln_no),  # charge, log-stable
        ]
    )


def _solve_log_densities(temp: float, n_f: float) -> Vec:
    """Levenberg-Marquardt-damped Newton (log-space) for `[ln n_H, ln n_O, ln n_e]`."""
    ln_kd = ln_k_dissoc(temp)
    ln_kh = ln_k_saha(IP_H, G_HP, G_H, temp)
    ln_ko = ln_k_saha(IP_O, G_OP, G_O, temp)
    ln_nf = float(np.log(n_f))

    # Init: H2O <=> 2H + O with n_H = 2 n_O (so 4 n_O^3 = K_d n_f) capped at full dissociation, then
    # n_e from charge neutrality. This satisfies the ratio + charge residuals exactly at every T, so
    # Newton only has to polish the O-conservation residual.
    ln_no = min(ln_nf, (ln_kd + ln_nf - np.log(4.0)) / 3.0)
    ln_nh = min(np.log(2.0) + ln_nf, np.log(2.0) + ln_no)
    ln_ne = 0.5 * float(np.logaddexp(ln_kh + ln_nh, ln_ko + ln_no))
    x = np.array([ln_nh, ln_no, ln_ne])

    eps = 1e-6
    lam = 1e-10
    for _ in range(200):
        r = _residual(x, ln_kd, ln_kh, ln_ko, n_f)
        if float(np.max(np.abs(r))) < 1e-11:
            break
        jac = np.empty((3, 3))
        for k in range(3):
            xp = x.copy()
            xp[k] += eps
            jac[:, k] = (_residual(xp, ln_kd, ln_kh, ln_ko, n_f) - r) / eps
        # LM-regularized step (J^T J + lam I) delta = -J^T r: robust if J is near-singular.
        jtj = jac.T @ jac
        delta = np.linalg.solve(jtj + lam * np.eye(3), -jac.T @ r)
        step = float(np.max(np.abs(delta)))
        if step > 2.0:  # cap the log-density step so a far init cannot overshoot into overflow
            delta *= 2.0 / step
        x = x + delta
    return x


def composition(rho: float, temp: float) -> Composition:
    """Equilibrium composition at `(rho [kg/m^3], temp [K])`."""
    n_f = rho / M_H2O  # H2O formula units per m^3 (conserves H:O = 2:1)
    ln_kd = ln_k_dissoc(temp)
    ln_kh = ln_k_saha(IP_H, G_HP, G_H, temp)
    ln_ko = ln_k_saha(IP_O, G_OP, G_O, temp)
    x = _solve_log_densities(temp, n_f)
    n_h2o, n_h, n_o, n_hp, n_op, n_e = _densities(x, ln_kd, ln_kh, ln_ko)
    return Composition(n_h2o=n_h2o, n_h=n_h, n_o=n_o, n_hp=n_hp, n_op=n_op, n_e=n_e)


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """Equilibrium `(p [Pa], e [J/kg])` at `(rho, temp)`.

    `p = (sum_i n_i) k T` (ideal mixture, including electrons; plasma coupling Gamma << 1 here).
    `e` = translational + H2O rotational/vibrational thermal energy + chemical (dissociation +
    ionization) energy, referenced to bound H2O = 0 so `e > 0` everywhere.
    """
    comp = composition(rho, temp)
    n_f = rho / M_H2O

    n_total = comp.n_h2o + comp.n_h + comp.n_o + comp.n_hp + comp.n_op + comp.n_e
    p = n_total * K_B * temp

    n_monatomic = comp.n_h + comp.n_o + comp.n_hp + comp.n_op + comp.n_e
    e_thermal = 1.5 * K_B * temp * n_monatomic  # translational, all single particles
    # H2O carries translational (3/2 kT) + rotational (3/2 kT, nonlinear) + vibrational.
    e_thermal += comp.n_h2o * (3.0 * K_B * temp + _e_vib(temp))
    # Chemical energy vs bound H2O: ionization potentials + the dissociated fraction's bond energy.
    e_chem = comp.n_hp * IP_H + comp.n_op * IP_O + (n_f - comp.n_h2o) * D_AT

    e = (e_thermal + e_chem) / rho
    return p, e


def _sound_speed_fd(
    pe: Callable[[float, float], tuple[float, float]], rho: float, temp: float
) -> float:
    """Adiabatic sound speed `c_s = sqrt((dp/drho)_s)` [m/s] for any `(p, e)(rho, T)` EOS, by
    central differences: `c_s^2 = (dp/drho)_T + T (dp/dT)_rho^2 / (rho^2 c_v)`, with
    `c_v = (de/dT)_rho` -- the standard decomposition that reuses the EOS at neighbouring states.
    """
    d = 1e-4
    p_rp, _ = pe(rho * (1.0 + d), temp)
    p_rm, _ = pe(rho * (1.0 - d), temp)
    p_tp, e_tp = pe(rho, temp * (1.0 + d))
    p_tm, e_tm = pe(rho, temp * (1.0 - d))

    dp_drho_t = (p_rp - p_rm) / (2.0 * d * rho)
    dp_dt_rho = (p_tp - p_tm) / (2.0 * d * temp)
    c_v = (e_tp - e_tm) / (2.0 * d * temp)
    cs2 = dp_drho_t + temp * dp_dt_rho**2 / (rho**2 * c_v)
    return float(np.sqrt(max(cs2, 0.0)))


def sound_speed(rho: float, temp: float) -> float:
    """Equilibrium adiabatic sound speed `c_s` [m/s] (see `_sound_speed_fd`)."""
    return _sound_speed_fd(pressure_energy, rho, temp)


# ---- Frozen-composition EOS (sudden-freeze bounding runs) --------------------------------------
#
# The equilibrium EOS above lets recombination return the dissociation/ionization energy during
# re-expansion. If the chemistry *freezes* (the nozzle-flow effect: three-body recombination rates
# collapse as the rebounding gas rarefies), the composition stops tracking equilibrium and that
# chemical energy stays locked. This block is the EOS for that regime: the same species set held at
# **fixed** per-formula-unit fractions `y_i = n_i / n_f`, so the mixture is a plain ideal gas with a
# constant mean molecular weight and a constant (inert) chemical energy offset. Used by the
# frozen-recombination bounding runs: frozen at the turnaround state (pessimistic, freeze *after*
# the plate) and frozen pure H2O (optimistic, freeze *before* the plate — no sink at all).


@dataclass(frozen=True)
class FrozenComposition:
    """Species fractions per H2O formula unit, `y_i = n_i / n_f`, held fixed (no chemistry)."""

    y_h2o: float
    y_h: float
    y_o: float
    y_hp: float
    y_op: float
    y_e: float


PURE_H2O_FROZEN = FrozenComposition(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
"""Undissociated molecular water with the chemistry switched off (freeze-before-the-plate)."""


def frozen_composition(rho_ref: float, t_ref: float) -> FrozenComposition:
    """Freeze the equilibrium composition at the reference state `(rho_ref, t_ref)`."""
    comp = composition(rho_ref, t_ref)
    n_f = rho_ref / M_H2O
    return FrozenComposition(
        y_h2o=comp.n_h2o / n_f,
        y_h=comp.n_h / n_f,
        y_o=comp.n_o / n_f,
        y_hp=comp.n_hp / n_f,
        y_op=comp.n_op / n_f,
        y_e=comp.n_e / n_f,
    )


def pressure_energy_frozen(rho: float, temp: float, y: FrozenComposition) -> tuple[float, float]:
    """Frozen-composition `(p [Pa], e [J/kg])` at `(rho, temp)`.

    Identical energetics to `pressure_energy` except the composition is `y` instead of the
    equilibrium one — so the chemical term is a *constant* specific-energy offset (locked, never
    exchanged with the thermal pool), and `p` is ideal with a fixed mean molecular weight.
    At the freeze reference state the two EOS agree exactly (splice continuity).
    """
    n_f = rho / M_H2O
    y_mono = y.y_h + y.y_o + y.y_hp + y.y_op + y.y_e
    p = n_f * (y.y_h2o + y_mono) * K_B * temp

    e_thermal = n_f * (1.5 * K_B * temp * y_mono + y.y_h2o * (3.0 * K_B * temp + _e_vib(temp)))
    e_chem = n_f * (y.y_hp * IP_H + y.y_op * IP_O + (1.0 - y.y_h2o) * D_AT)
    return p, (e_thermal + e_chem) / rho


def sound_speed_frozen(rho: float, temp: float, y: FrozenComposition) -> float:
    """Frozen adiabatic sound speed `c_s` [m/s] (see `_sound_speed_fd`)."""
    return _sound_speed_fd(lambda r, t: pressure_energy_frozen(r, t, y), rho, temp)


def eos_grid_frozen(rho_grid: Vec, t_grid: Vec, y: FrozenComposition) -> tuple[Vec, Vec, Vec]:
    """Evaluate the frozen `(p, e, c_s)` on the full `(rho, T)` grid, row-major over `(rho, T)`.

    Closed-form per node (no Newton solve), so this is cheap even on the production grid.
    """
    n_rho, n_t = len(rho_grid), len(t_grid)
    p = np.empty((n_rho, n_t))
    e = np.empty((n_rho, n_t))
    cs = np.empty((n_rho, n_t))
    for i, rho in enumerate(rho_grid):
        for j, temp in enumerate(t_grid):
            p[i, j], e[i, j] = pressure_energy_frozen(float(rho), float(temp), y)
            cs[i, j] = sound_speed_frozen(float(rho), float(temp), y)
    return p, e, cs


def eos_grid(rho_grid: Vec, t_grid: Vec) -> tuple[Vec, Vec, Vec]:
    """Evaluate `(p, e, c_s)` on the full `(rho, T)` grid, row-major over `(rho, T)`.

    Returns three `(n_rho, n_T)` arrays. Each `(rho, T)` is solved independently; the
    physics-based initial guess makes the per-node Newton robust across the molecular/plasma span.
    """
    n_rho, n_t = len(rho_grid), len(t_grid)
    p = np.empty((n_rho, n_t))
    e = np.empty((n_rho, n_t))
    cs = np.empty((n_rho, n_t))
    for i, rho in enumerate(rho_grid):
        for j, temp in enumerate(t_grid):
            p[i, j], e[i, j] = pressure_energy(float(rho), float(temp))
            cs[i, j] = sound_speed(float(rho), float(temp))
    return p, e, cs
