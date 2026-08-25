"""Chemical-equilibrium equation of state for water, cold vapor to plasma (B5c-1).

The high-v PuffSat package needs an *equilibrium* water EOS (design §7): at 16 km/s the stagnated
gas is a 10-40 kK plasma, and the dominant control on the restitution `e_eff` is the
dissociation + ionization energy buffering (design §4) -- the specific heat sink that decides how
much impact kinetic energy disappears into internal modes instead of staying available to rebound.
CoolProp only reaches ~1300 K, so this is built from first principles (statistical mechanics), with
no external data. See `puffsat_impact_sim_design.md` §3-4, §7, and Zel'dovich & Raizer, *Physics of
Shock Waves and High-Temperature Hydrodynamic Phenomena*, Ch. III.

Model (species set): H2O, OH, H2, O2, H, O, H+, O ions O+ .. O8+ (the full Saha ladder), e- in
chemical + ionization equilibrium at each `(rho, T)`. Dissociation H2O <=> 2H + O and each diatomic
AB <=> A + B by law of mass action; each ionization stage by the Saha equation; closed by
H:O = 2:1 element conservation and charge neutrality. The multi-stage oxygen ladder
(2026-07, Jupiter-retrograde 69 km/s scenario)
extends the original single-stage model: below ~30 kK only the first stages are populated and the
two models agree, but at ~2.4 GJ/kg stagnation (69 km/s) oxygen climbs to O4+..O6+ and the
~0.4-2.3 GJ/kg of multi-ionization energy is the dominant specific-heat sink — without the ladder
the table would overshoot the stagnation temperature severely.

Molecular intermediates OH, H2, O2 (added 2026-08-25): the transition region became load-bearing.
The original model carried the single effective reaction H2O <=> 2H + O, which preserves both
endpoints (cold molecular vapor; hot ionized plasma) and the total energy invested -- what sets
`e_eff` -- but reshapes the ~2000-6000 K transition. Q-P put a freeze point at 3908 K, inside that
band, and made the water-reformation partner (`H + OH + M -> H2O + M`) a rate-limiting density that
had to be proxied by `n_H`. So the three intermediates are now carried explicitly. They add **no
unknowns**: each is fixed by mass action off the same `(n_H, n_O)` the solve already carries, and
they enter only through the element-conservation residuals and the energy sum. Both endpoints are
unmoved -- see `test_hot_endpoint_is_unmoved_by_the_intermediates`.

Energy reference: bound molecular H2O at T -> 0 has e = 0, so dissociation and ionization only *add*
energy. Every `e` on the grid is therefore strictly positive -- required by the Rust table loader,
which interpolates `ln e` (ADR-0007).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
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

# Ground-term degeneracies of the full oxygen ionization ladder O I .. O IX (bare nucleus):
# 3P, 4S, 3P, 2P, 1S, 2S, 1S, 2S, bare. Ground terms only (no excited-state partition sums,
# no continuum lowering) — adequate for the energy bookkeeping that sets `e_eff`; the ladder's
# role is the ionization-potential sink, not spectroscopic fidelity.
G_O_LADDER = (9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0)

# --- Energetics ---
IP_H = 13.598 * EV  # H first ionization potential [J]
IP_O = 13.618 * EV  # O first ionization potential [J]

# Successive oxygen ionization potentials O I->II .. O VIII->IX [J] (NIST ASD). The K-shell pair
# (739, 871 eV) only opens above ~10^6 K; carried for completeness so the table's top edge stays
# physical rather than an artificial full-stripping plateau.
IP_O_LADDER = tuple(
    ip_ev * EV for ip_ev in (13.618, 35.121, 54.936, 77.414, 113.899, 138.120, 739.327, 871.410)
)
# Cumulative ionization energy of the O ion at charge k (index k-1): sum of the first k potentials.
E_O_CUM = tuple(float(np.cumsum(np.array(IP_O_LADDER))[k]) for k in range(len(IP_O_LADDER)))
N_O_STAGES = len(IP_O_LADDER)  # 8 ion charge states O+ .. O8+
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


# --- Molecular intermediates: OH, H2, O2 --------------------------------------------------------
#
# Added 2026-08-25, because the 2000-6000 K transition became load-bearing (Q-P). They introduce no
# new unknowns: given `(n_H, n_O)` each is fixed by mass action, `n_AB = n_A n_B / K_AB`. What they
# change is the element bookkeeping and the energy sum.
#
# Dissociation energies at 0 K come from the *same* NIST-JANAF 0 K heats of formation that set
# `D_AT` (H 216.035, O 246.844, OH 37.278, H2 0, O2 0 kJ/mol), so the set is internally consistent:
#   D0(OH) = dHf(H) + dHf(O) - dHf(OH) = 425.60 kJ/mol, and `D_AT - D0(OH)` = 492.1 kJ/mol is
#           water's *first* O-H bond -- the cross-check that no sign or dissociation limit is wrong.
#   D0(H2) = 2 dHf(H) = 432.07 kJ/mol;  D0(O2) = 2 dHf(O) = 493.69 kJ/mol.
# Spectroscopic constants: Herzberg, *Molecular Spectra and Molecular Structure* I. Ground
# electronic terms only, matching the H2O and O-ladder treatment -- O2's low-lying a(1)Delta_g and
# b(1)Sigma_g+ are not summed, which slightly understates O2 in a band where it never exceeds a few
# percent of the oxygen.


@dataclass(frozen=True)
class Diatomic:
    """A diatomic intermediate AB and the equilibrium `AB <=> A + B` that sets its density."""

    name: str
    mass: float  # molecular mass [kg]
    degeneracy: float  # ground-term electronic degeneracy of AB
    d0: float  # dissociation energy from v=0 [J per molecule]
    theta_rot: float  # rotational temperature hcB/k [K]
    symmetry: float  # rotational symmetry number (2 for homonuclear)
    theta_vib: float  # vibrational temperature hc*omega_e/k [K]
    mass_a: float  # atom A mass [kg]
    degeneracy_a: float  # atom A ground-term degeneracy
    mass_b: float  # atom B mass [kg]
    degeneracy_b: float  # atom B ground-term degeneracy


OH = Diatomic(
    name="OH",
    mass=M_O + M_H,
    degeneracy=4.0,  # X 2Pi: Lambda-doubling (2) x spin (2)
    d0=425.60e3 / N_A,
    theta_rot=_CM_TO_K * 18.911,
    symmetry=1.0,
    theta_vib=_CM_TO_K * 3737.76,
    mass_a=M_H,
    degeneracy_a=G_H,
    mass_b=M_O,
    degeneracy_b=G_O,
)
H2 = Diatomic(
    name="H2",
    mass=2.0 * M_H,
    degeneracy=1.0,  # X 1Sigma_g+
    d0=432.07e3 / N_A,
    theta_rot=_CM_TO_K * 60.853,
    symmetry=2.0,
    theta_vib=_CM_TO_K * 4401.21,
    mass_a=M_H,
    degeneracy_a=G_H,
    mass_b=M_H,
    degeneracy_b=G_H,
)
O2 = Diatomic(
    name="O2",
    mass=2.0 * M_O,
    degeneracy=3.0,  # X 3Sigma_g-
    d0=493.69e3 / N_A,
    theta_rot=_CM_TO_K * 1.4377,
    symmetry=2.0,
    theta_vib=_CM_TO_K * 1580.19,
    mass_a=M_O,
    degeneracy_a=G_O,
    mass_b=M_O,
    degeneracy_b=G_O,
)
INTERMEDIATES = (OH, H2, O2)


@dataclass(frozen=True)
class Composition:
    """Equilibrium number densities [m^-3] at a single `(rho, T)`.

    `n_o_ions[k]` is the O ion of charge `k+1` (O+ .. O8+); `n_op` remains the singly-charged
    O+ density (`n_o_ions[0]`) for the pre-ladder call sites (frozen-composition bracket)."""

    n_h2o: float
    n_h: float
    n_o: float
    n_hp: float
    n_o_ions: tuple[float, ...]
    n_e: float
    #: Molecular intermediates. Defaulted to zero so a hand-built atoms-only `Composition` stays
    #: constructible, but `composition()` always fills them.
    n_oh: float = 0.0
    n_h2: float = 0.0
    n_o2: float = 0.0

    @property
    def n_op(self) -> float:
        """Singly-ionized oxygen density [m^-3] (the first rung of the ladder)."""
        return self.n_o_ions[0]

    @property
    def n_neutral_heavy(self) -> float:
        """All neutral heavy particles [m^-3] -- every collision partner that is not an ion or an
        electron. This is the third-body density for three-body recombination and the neutral
        target density for electron-neutral collisions."""
        return self.n_h2o + self.n_oh + self.n_h2 + self.n_o2 + self.n_h + self.n_o

    @property
    def n_total(self) -> float:
        """Every particle [m^-3], electrons included -- the ideal-mixture pressure count."""
        return self.n_neutral_heavy + self.n_hp + sum(self.n_o_ions) + self.n_e


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


def _ln_ko_ladder(temp: float) -> tuple[float, ...]:
    """ln Saha constants for every O stage `O^{k} <=> O^{k+1} + e-`, k = 0 .. N_O_STAGES-1."""
    return tuple(
        ln_k_saha(IP_O_LADDER[k], G_O_LADDER[k + 1], G_O_LADDER[k], temp) for k in range(N_O_STAGES)
    )


def _z_rot_linear(theta_rot: float, symmetry: float, temp: float) -> float:
    """Classical rigid-rotor partition function of a linear molecule, `T / (sigma Theta_rot)`.

    Valid for `T >> Theta_rot`; the stiffest case here is H2 at 87.6 K, so the whole regime of
    interest (>= ~400 K) is classical to better than a percent.
    """
    return float(temp / (symmetry * theta_rot))


def _z_vib_mode(theta_vib: float, temp: float) -> float:
    """Harmonic-oscillator partition function of one mode, referenced to `v = 0`."""
    return float(1.0 / (1.0 - np.exp(-theta_vib / temp)))


def _e_vib_mode(theta_vib: float, temp: float) -> float:
    """Mean vibrational energy of one harmonic mode [J], excluding zero-point."""
    return float(theta_vib * K_B / (np.exp(theta_vib / temp) - 1.0))


def ln_k_diatomic(d: Diatomic, temp: float) -> float:
    """ln of `K = n_A n_B / n_AB` [ln m^-3] for the dissociation `AB <=> A + B`.

    Same law of mass action as `ln_k_dissoc`, one bond instead of three: the ratio of the atoms'
    translational-electronic partition functions to the molecule's, times `exp(-D0/kT)`. `D0` is
    measured from `v = 0` and `_z_vib_mode` is referenced there too, so the two conventions agree.
    """
    ln_a = _ln_trans(d.mass_a, temp) + float(np.log(d.degeneracy_a))
    ln_b = _ln_trans(d.mass_b, temp) + float(np.log(d.degeneracy_b))
    ln_ab = _ln_trans(d.mass, temp) + float(
        np.log(
            d.degeneracy
            * _z_rot_linear(d.theta_rot, d.symmetry, temp)
            * _z_vib_mode(d.theta_vib, temp)
        )
    )
    return float(ln_a + ln_b - ln_ab - d.d0 / (K_B * temp))


@dataclass(frozen=True)
class _EqConstants:
    """Every equilibrium constant at one temperature, in log form.

    Bundled because the intermediates took the threaded-constant count from four to seven, and the
    same set is needed by `_densities`, `_residual` and both initial guesses.
    """

    ln_kd: float  # H2O <=> 2H + O
    ln_kh: float  # H <=> H+ + e-
    ln_ko: tuple[float, ...]  # the O ionization ladder
    ln_k_mol: tuple[float, ...]  # one per entry of INTERMEDIATES, in that order


def _constants(temp: float) -> _EqConstants:
    """All equilibrium constants at `temp`."""
    return _EqConstants(
        ln_kd=ln_k_dissoc(temp),
        ln_kh=ln_k_saha(IP_H, G_HP, G_H, temp),
        ln_ko=_ln_ko_ladder(temp),
        ln_k_mol=tuple(ln_k_diatomic(d, temp) for d in INTERMEDIATES),
    )


def _densities(x: Vec, c: _EqConstants) -> Composition:
    """Reconstruct the full composition from log-densities `x = [ln n_H, ln n_O, ln n_e]`.

    Everything else follows by mass action from the three unknowns: `n_H2O` from the atomization
    constant, each intermediate `n_AB = n_A n_B / K_AB`, and each O ion by one more link of the
    chained Saha ladder (`ln K_k - ln n_e` per stage). Kept in log form so the cold limit, where
    the equilibrium constants are ~1e-300, never underflows.
    """
    ln_nh, ln_no, ln_ne = float(x[0]), float(x[1]), float(x[2])
    n_h = float(np.exp(min(ln_nh, _EXP_CAP)))
    n_o = float(np.exp(min(ln_no, _EXP_CAP)))
    n_e = float(np.exp(min(ln_ne, _EXP_CAP)))
    n_h2o = float(np.exp(min(2.0 * ln_nh + ln_no - c.ln_kd, _EXP_CAP)))
    n_hp = float(np.exp(min(c.ln_kh + ln_nh - ln_ne, _EXP_CAP)))
    # Intermediates, in INTERMEDIATES order: OH from (H, O), H2 from (H, H), O2 from (O, O).
    ln_ab = (ln_nh + ln_no, 2.0 * ln_nh, 2.0 * ln_no)
    n_oh, n_h2, n_o2 = (
        float(np.exp(min(ln_ab[i] - c.ln_k_mol[i], _EXP_CAP))) for i in range(len(INTERMEDIATES))
    )
    n_o_ions = []
    ln_stage = ln_no
    for k in range(N_O_STAGES):
        ln_stage = ln_stage + c.ln_ko[k] - ln_ne
        n_o_ions.append(float(np.exp(min(ln_stage, _EXP_CAP))))
    return Composition(
        n_h2o=n_h2o,
        n_h=n_h,
        n_o=n_o,
        n_hp=n_hp,
        n_o_ions=tuple(n_o_ions),
        n_e=n_e,
        n_oh=n_oh,
        n_h2=n_h2,
        n_o2=n_o2,
    )


#: Precomputed logs used inside the residual. It is evaluated four times per Newton iteration (the
#: finite-difference Jacobian), so the scalar work there is worth keeping out of numpy.
_LN2 = math.log(2.0)
_LN_CHARGE = tuple(math.log(k + 1.0) for k in range(N_O_STAGES))


def _ln_sum(terms: Sequence[float]) -> float:
    """`ln(sum_i exp(terms_i))`, evaluated without leaving log space.

    Plain-`math` shifted sum rather than `np.logaddexp.reduce`: the arrays here are 3-10 elements
    and numpy's per-call overhead dominated the equilibrium solve once the intermediates doubled
    the number of reductions.
    """
    top = max(terms)
    if top == -math.inf:
        return -math.inf
    return top + math.log(sum(math.exp(t - top) for t in terms))


def _residual(x: Vec, c: _EqConstants, n_f: float) -> Vec:
    """Equilibrium residual at log-densities `x = [ln n_H, ln n_O, ln n_e]`.

    Three independent, well-conditioned constraints (the naive H+O number pair degenerates when
    H2O dominates -- both then merely say `n_H2O = n_f`):
    - **O nuclei** conservation, now counting the O carried by OH and O2 as well (pins `n_H2O` ~
      `n_f` in the cold limit),
    - **H:O = 2:1** stoichiometry as a ratio over the non-H2O species only, **in log form**. H2O
      is 2:1 by construction, so `N_H = 2 N_O` holds if and only if it holds with H2O's
      contribution removed from both sides -- an exact restatement, and the one that stays
      conditioned when H2O carries essentially all the nuclei. The intermediates are what make it
      a real constraint rather than an identity: OH is 1:1, H2 is 2:0, O2 is 0:2. Log form is not
      cosmetic here: at 50 K the trace pool is ~1e-100 in O2 against ~1e-166 in H2, and the linear
      ratio is a division of two underflowing numbers.
    - **charge neutrality** `n_e = n_H+ + sum_k (k+1) n_O(k+1)+`, written in log form as
      `ln n_e - ln(sum of charge-weighted ion terms)`. Each ion term carries its own `-j ln n_e`
      Saha chain, so in the weak-ionization limit this reduces to the classic
      `2 ln n_e = ln(K_H n_H + K_O n_O)` with slope 2 in `ln n_e` -- well-conditioned everywhere.
    """
    ln_nh, ln_no, ln_ne = float(x[0]), float(x[1]), float(x[2])
    ln_h2o = 2.0 * ln_nh + ln_no - c.ln_kd
    ln_hp = c.ln_kh + ln_nh - ln_ne
    ln_oh, ln_h2, ln_o2 = (
        ln_nh + ln_no - c.ln_k_mol[0],
        2.0 * ln_nh - c.ln_k_mol[1],
        2.0 * ln_no - c.ln_k_mol[2],
    )

    # The O ladder, one Saha link per stage; kept as logs for both the charge sum and the O count.
    ln_stages: list[float] = []
    ln_stage = ln_no
    for k in range(N_O_STAGES):
        ln_stage = ln_stage + c.ln_ko[k] - ln_ne
        ln_stages.append(ln_stage)
    ln_charge = _ln_sum([ln_hp, *(_LN_CHARGE[k] + t for k, t in enumerate(ln_stages))])

    # Nuclei carried by everything except H2O, in log form.
    ln_h_free = _ln_sum([ln_nh, ln_hp, ln_oh, _LN2 + ln_h2])
    ln_o_free = _ln_sum([ln_no, *ln_stages, ln_oh, _LN2 + ln_o2])

    n_o_total = math.exp(min(ln_h2o, _EXP_CAP)) + math.exp(min(ln_o_free, _EXP_CAP))
    return np.array(
        [
            n_o_total / n_f - 1.0,  # O nuclei
            ln_h_free - ln_o_free - _LN2,  # H:O, H2O's exact 2:1 divided out
            ln_ne - ln_charge,  # charge neutrality, log-stable
        ]
    )


def _cold_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Molecular/weakly-ionized init: H2O <=> 2H + O with n_H = 2 n_O (so 4 n_O^3 = K_d n_f)
    capped at full dissociation, then n_e from single-stage charge neutrality capped at the
    full-stripping ceiling (10 e- per formula unit).

    Deliberately still the atoms-only guess: the intermediates shift where the atomic densities
    land but not by orders of magnitude in log space, which is the accuracy a Newton start needs.
    """
    ln_no = min(ln_nf, (c.ln_kd + ln_nf - np.log(4.0)) / 3.0)
    ln_nh = min(np.log(2.0) + ln_nf, np.log(2.0) + ln_no)
    ln_ne = min(
        0.5 * float(np.logaddexp(c.ln_kh + ln_nh, c.ln_ko[0] + ln_no)),
        float(np.log(10.0)) + ln_nf,
    )
    return np.array([ln_nh, ln_no, ln_ne])


def _molecular_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Init for the transition band, where the *intermediates* hold the broken water.

    The cold init assumes water goes straight to atoms, so in the 1500-4000 K band -- where the
    first O-H bond (492 kJ/mol) is breaking but the second (426 kJ/mol) is not -- it can start
    several decades off in `ln n_O`. Here the products are taken to be OH + H instead:
    `H2O <=> OH + H` has `K = K_d / K_OH`, and with `n_OH = n_H` that gives `n_H^2 = K n_f`.
    """
    ln_k_first = c.ln_kd - c.ln_k_mol[0]  # H2O <=> OH + H
    ln_nh = min(np.log(2.0) + ln_nf, 0.5 * (ln_k_first + ln_nf))
    # n_O then follows from OH equilibrium with n_OH ~ n_H: n_O = n_OH K_OH / n_H = K_OH.
    ln_no = min(ln_nf, c.ln_k_mol[0])
    ln_ne = min(
        0.5 * float(np.logaddexp(c.ln_kh + ln_nh, c.ln_ko[0] + ln_no)),
        float(np.log(10.0)) + ln_nf,
    )
    return np.array([ln_nh, ln_no, ln_ne])


def _hot_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Hot-plasma init (Jupiter-retrograde regime): a mean-charge fixed point on the Saha ladder.

    Guess `n_e = q n_f`; a stage is 'climbed' when its Saha constant exceeds `n_e`; the H+
    fraction is `K_H/(K_H + n_e)`; iterate the implied mean charge `q = 2 f_H+ + k_dom`. Then
    back the *neutral* O density down the chain from the dominant stage (population ~ n_f)."""
    q = 1.0
    k_dom = 0
    f_hp = 0.5
    for _ in range(30):
        ln_ne = float(np.log(q)) + ln_nf
        k_dom = sum(1 for lk in c.ln_ko if lk > ln_ne)
        f_hp = float(np.exp(c.ln_kh - np.logaddexp(c.ln_kh, ln_ne)))
        q_new = max(2.0 * f_hp + float(k_dom), 1e-3)
        if abs(q_new - q) < 1e-3 * q:
            q = q_new
            break
        q = 0.5 * (q + q_new)
    ln_ne = float(np.log(q)) + ln_nf
    # Neutral H from the ionized fraction; neutral O backed down the k_dom Saha links.
    ln_nh = float(np.log(2.0)) + ln_nf + float(np.log(max(1.0 - f_hp, 1e-12)))
    ln_no = ln_nf - sum(c.ln_ko[j] - ln_ne for j in range(k_dom))
    return np.array([ln_nh, ln_no, ln_ne])


def _newton_polish(x0: Vec, c: _EqConstants, n_f: float) -> tuple[Vec, float]:
    """LM-damped Newton (log-space) from `x0`; returns `(x, max |residual|)`."""
    x = x0.copy()
    eps = 1e-6
    lam = 1e-10
    r = _residual(x, c, n_f)
    for _ in range(400):
        if float(np.max(np.abs(r))) < 1e-11:
            break
        jac = np.empty((3, 3))
        for k in range(3):
            xp = x.copy()
            xp[k] += eps
            jac[:, k] = (_residual(xp, c, n_f) - r) / eps
        # LM-regularized step (J^T J + lam I) delta = -J^T r: robust if J is near-singular.
        jtj = jac.T @ jac
        delta = np.linalg.solve(jtj + lam * np.eye(3), -jac.T @ r)
        step = float(np.max(np.abs(delta)))
        if step > 2.0:  # cap the log-density step so a far init cannot overshoot into overflow
            delta *= 2.0 / step
        x = x + delta
        r = _residual(x, c, n_f)
    return x, float(np.max(np.abs(r)))


def _solve_log_densities(temp: float, n_f: float) -> Vec:
    """Solve `[ln n_H, ln n_O, ln n_e]` by Newton from each of three physics-based inits.

    The three bracket the three regimes the EOS spans -- bound molecular water, the intermediate
    band where OH holds the broken water, and the stripped plasma -- and the first that converges
    wins. The molecular init is the one the added species made necessary.
    """
    c = _constants(temp)
    ln_nf = float(np.log(n_f))

    best: tuple[Vec, float] | None = None
    for init in (_cold_init(ln_nf, c), _molecular_init(ln_nf, c), _hot_init(ln_nf, c)):
        x, res = _newton_polish(init, c, n_f)
        if res < 1e-8:
            return x
        if best is None or res < best[1]:
            best = (x, res)
    assert best is not None
    raise RuntimeError(
        f"equilibrium solve did not converge at T={temp} K, n_f={n_f} m^-3 "
        f"(best max residual {best[1]:.2e})"
    )


def composition(rho: float, temp: float) -> Composition:
    """Equilibrium composition at `(rho [kg/m^3], temp [K])`."""
    n_f = rho / M_H2O  # H2O formula units per m^3 (conserves H:O = 2:1)
    return _densities(_solve_log_densities(temp, n_f), _constants(temp))


def _intermediate_densities(comp: Composition) -> tuple[float, ...]:
    """The intermediates' densities in `INTERMEDIATES` order, so the two sums below stay aligned
    with the species table rather than repeating `n_oh, n_h2, n_o2` by hand."""
    return (comp.n_oh, comp.n_h2, comp.n_o2)


def _intermediate_thermal_energy(comp: Composition, temp: float) -> float:
    """Thermal energy of the diatomic intermediates [J/m^3].

    Each carries translational `3/2 kT` + rotational `kT` (linear, classical) + one vibrational
    mode. That is `5/2 kT` of sensible heat plus vibration, against the `3/2 kT` of an atom -- so
    the intermediates are also a (small) specific-heat sink in their own right, not only a store
    of bond energy.
    """
    total = 0.0
    for d, n in zip(INTERMEDIATES, _intermediate_densities(comp), strict=True):
        total += n * (2.5 * K_B * temp + _e_vib_mode(d.theta_vib, temp))
    return total


def _bond_energy_held(comp: Composition, n_f: float) -> float:
    """Chemical energy stored as broken H2O bonds [J/m^3], net of the bonds the intermediates hold.

    Fully atomizing the water costs `n_f D_AT`; intact H2O owes none of it, and every intermediate
    hands part of it back. Reduces to the old `(n_f - n_H2O) D_AT` when the intermediates vanish,
    and gives the right answer at the one state that can be checked by hand: all-OH-plus-H is
    `n_f (D_AT - D0_OH)` = water's first O-H bond, 492 kJ/mol.
    """
    held = (n_f - comp.n_h2o) * D_AT
    for d, n in zip(INTERMEDIATES, _intermediate_densities(comp), strict=True):
        held -= n * d.d0
    return held


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """Equilibrium `(p [Pa], e [J/kg])` at `(rho, temp)`.

    `p = (sum_i n_i) k T` (ideal mixture, including electrons; plasma coupling Gamma << 1 here).
    `e` = translational + H2O rotational/vibrational thermal energy + chemical (dissociation +
    ionization) energy, referenced to bound H2O = 0 so `e > 0` everywhere.
    """
    comp = composition(rho, temp)
    n_f = rho / M_H2O

    p = comp.n_total * K_B * temp

    n_monatomic = comp.n_h + comp.n_o + comp.n_hp + sum(comp.n_o_ions) + comp.n_e
    e_thermal = 1.5 * K_B * temp * n_monatomic  # translational, all single particles
    # H2O carries translational (3/2 kT) + rotational (3/2 kT, nonlinear) + vibrational.
    e_thermal += comp.n_h2o * (3.0 * K_B * temp + _e_vib(temp))
    e_thermal += _intermediate_thermal_energy(comp, temp)
    # Chemical energy vs bound H2O: ionization potentials (each O stage carries its *cumulative*
    # ladder energy) + the bond energy not currently paid back by a bond.
    e_ionization_o = sum(n * E_O_CUM[k] for k, n in enumerate(comp.n_o_ions))
    e_chem = comp.n_hp * IP_H + e_ionization_o + _bond_energy_held(comp, n_f)

    e = (e_thermal + e_chem) / rho
    return p, e


#: Bond energy of one kilogram of fully atomized water [J/kg] -- `D_AT / M_H2O` = 50.9 MJ/kg.
#: The ceiling on what the dissociation store can strand, and the denominator `fireball` charges
#: the freeze against.
FULL_ATOMIZATION_ENERGY = D_AT / M_H2O


def bond_energy_held(rho: float, temp: float) -> float:
    """Chemical energy still locked in broken H2O bonds [J/kg] at equilibrium `(rho, temp)`.

    **Not `(1 - n_H2O/n_f) * D_AT/M_H2O`, which is what this replaced.** That form charges the
    full atomization energy for every formula unit that is not intact water, and with OH, H2 and
    O2 in the species set that is wrong by a lot: a gas half in OH has handed back roughly half
    its bond energy while still reading as "fully dissociated". What is stranded is the energy no
    bond is currently paying for, which is exactly `_bond_energy_held`.
    """
    return _bond_energy_held(composition(rho, temp), rho / M_H2O) / rho


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
    """Species fractions per H2O formula unit, `y_i = n_i / n_f`, held fixed (no chemistry).

    `y_o_ions[k]` is the O ion of charge `k+1` (O+ .. O8+), mirroring `Composition.n_o_ions` —
    a 69 km/s-class turnaround state is multi-charged, so the whole ladder freezes, each stage
    locking its *cumulative* ionization energy."""

    y_h2o: float
    y_h: float
    y_o: float
    y_hp: float
    y_o_ions: tuple[float, ...]
    y_e: float
    y_oh: float = 0.0
    y_h2: float = 0.0
    y_o2: float = 0.0

    @property
    def y_op(self) -> float:
        """Singly-ionized oxygen fraction (the first rung of the ladder)."""
        return self.y_o_ions[0]

    @property
    def y_intermediates(self) -> tuple[float, ...]:
        """Intermediate fractions in `INTERMEDIATES` order."""
        return (self.y_oh, self.y_h2, self.y_o2)


PURE_H2O_FROZEN = FrozenComposition(
    y_h2o=1.0, y_h=0.0, y_o=0.0, y_hp=0.0, y_o_ions=(0.0,) * N_O_STAGES, y_e=0.0
)
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
        y_o_ions=tuple(n / n_f for n in comp.n_o_ions),
        y_e=comp.n_e / n_f,
        y_oh=comp.n_oh / n_f,
        y_h2=comp.n_h2 / n_f,
        y_o2=comp.n_o2 / n_f,
    )


def pressure_energy_frozen(rho: float, temp: float, y: FrozenComposition) -> tuple[float, float]:
    """Frozen-composition `(p [Pa], e [J/kg])` at `(rho, temp)`.

    Identical energetics to `pressure_energy` except the composition is `y` instead of the
    equilibrium one — so the chemical term is a *constant* specific-energy offset (locked, never
    exchanged with the thermal pool), and `p` is ideal with a fixed mean molecular weight.
    At the freeze reference state the two EOS agree exactly (splice continuity).
    """
    n_f = rho / M_H2O
    y_o_ions_total = sum(y.y_o_ions)
    y_mono = y.y_h + y.y_o + y.y_hp + y_o_ions_total + y.y_e
    y_inter = sum(y.y_intermediates)
    p = n_f * (y.y_h2o + y_mono + y_inter) * K_B * temp

    e_thermal = n_f * (1.5 * K_B * temp * y_mono + y.y_h2o * (3.0 * K_B * temp + _e_vib(temp)))
    # The frozen intermediates keep their sensible heat (5/2 kT + vibration) even though their
    # bond energy is locked -- the composition is frozen, not the thermal degrees of freedom.
    for d, frac in zip(INTERMEDIATES, y.y_intermediates, strict=True):
        e_thermal += n_f * frac * (2.5 * K_B * temp + _e_vib_mode(d.theta_vib, temp))
    e_ionization_o = sum(frac * E_O_CUM[k] for k, frac in enumerate(y.y_o_ions))
    e_bond = (1.0 - y.y_h2o) * D_AT - sum(
        frac * d.d0 for d, frac in zip(INTERMEDIATES, y.y_intermediates, strict=True)
    )
    e_chem = n_f * (y.y_hp * IP_H + e_ionization_o + e_bond)
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
