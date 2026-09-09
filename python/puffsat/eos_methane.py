"""Chemical-equilibrium equation of state for methane, cold vapour to plasma (N9/N10).

ADR-0016 of the companion repository admits a **walled** de Laval chamber on the head-on
departure burn and names methane as the fluid it would fly. Ask N10 item 4b then found that the
headline turns on a composition rather than on a nozzle: the slug's hydrogen is only *partly*
dissociated at chamber conditions, so the atomisation store is only partly charged, and it
responds to pressure as well as to temperature. That is an equilibrium composition at `(rho, T)`
and nothing else -- so it is solved here, once, and both N9 and N10 read it.

**This is `eos_water.py`'s frame with a different element pair.** Three log-unknowns
(`ln n_C`, `ln n_H`, `ln n_e`); every other species follows by law of mass action; closed by
carbon-nuclei conservation, the H:C = 4:1 stoichiometry and charge neutrality. What is new is
that methane's molecules are not all diatomic, so `eos_water.Diatomic` generalises here to
`Molecule`, which carries a stoichiometry `(n_c, n_h)`, a linear *or* nonlinear rigid rotor and
an arbitrary list of harmonic modes.

Species set: `CH4`, `CH3`, `CH2`, `CH`, `C2H2`, `C2`, `C3`, `H2`, `C`, `H`, `H+`,
`C+ .. C6+` (the full Saha ladder), `e-`.

**Condensed carbon is deliberately absent.** Graphite is not a gas-phase mass-action species and
soot does not form by a three-body collision, so putting it in this solve would silently answer
N10 item 5 by assumption. It gets its own nucleation treatment; see ADR-0050.

Energy reference: bound `CH4` at `T -> 0` has `e = 0`, so dissociation and ionisation only *add*
energy and every `e` is strictly positive -- the same convention `eos_water` uses.

Thermochemistry: NIST-JANAF **0 K** heats of formation throughout, matching `eos_water`.
ADR-0016's paper-side arithmetic used 298 K values, which is why its atomisation energy reads
103.7 MJ/kg against this module's 102.35. That is a reference-state convention, not a
disagreement (ADR-0050).

Spectroscopic constants: Herzberg, *Molecular Spectra and Molecular Structure*, and the NIST
Chemistry WebBook. Ground electronic terms only, as in `eos_water`.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat import eos_water
from puffsat.eos_water import (
    _CM_TO_K,
    _EXP_CAP,
    G_H,
    _ln_trans,
    ln_k_saha,
    sound_speed_fd,
)

# Physical constants, **explicitly re-bound rather than implicitly re-exported**. Callers of this
# EOS (`walled_nozzle.*`) need `K_B` and `IP_H` beside the composition they just solved, and a
# downstream module should not have to know that methane's Boltzmann constant lives in the water
# EOS. Binding them here also makes them part of this module's stated interface under
# `mypy --strict`, which a pass-through import is not.
K_B = eos_water.K_B
N_A = eos_water.N_A
AMU = eos_water.AMU
EV = eos_water.EV
M_H = eos_water.M_H
IP_H = eos_water.IP_H

Vec = NDArray[np.float64]

# --- Species masses [kg] ---
M_C = 12.011 * AMU
M_CH4 = M_C + 4.0 * M_H

# --- Ground-term electronic degeneracies ---
G_C = 9.0  # C 3P (sum 2J+1 over the ground term)
#: Ground-term degeneracies of the carbon ionisation ladder C I .. C VII (bare nucleus):
#: 3P, 2P*, 1S, 2S, 1S, 2S, bare.
G_C_LADDER = (9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0)

#: Successive carbon ionisation potentials C I->II .. C VI->VII [J] (NIST ASD). The K-shell pair
#: (392, 490 eV) opens only far above anything the chamber reaches; carried so the table's top
#: edge stays physical rather than plateauing at an artificial full-stripping ceiling.
IP_C_LADDER = tuple(ip_ev * EV for ip_ev in (11.2603, 24.3833, 47.8878, 64.4939, 392.087, 489.993))
#: Cumulative ionisation energy of the C ion at charge `k+1`.
E_C_CUM = tuple(float(np.cumsum(np.array(IP_C_LADDER))[k]) for k in range(len(IP_C_LADDER)))
N_C_STAGES = len(IP_C_LADDER)  # 6 ion charge states C+ .. C6+

# --- 0 K heats of formation [kJ/mol] (NIST-JANAF) -----------------------------------------------
#
# Every dissociation energy below is derived from *these* numbers rather than quoted separately,
# so the species set is internally consistent by construction: a molecule's atomisation energy is
# `sum(atoms) - itself`, and no sign or dissociation limit can disagree with another entry.
_HF0_H = 216.035
_HF0_C = 711.19
_HF0 = {
    "CH4": -66.63,
    "CH3": 149.8,
    "CH2": 390.7,  # X 3B1, the ground (triplet) state
    "CH": 592.5,
    "C2H2": 228.8,
    "C2": 820.0,
    "C3": 814.0,
    "H2": 0.0,
    "C2H": 566.1,
}


def _atomization(name: str, n_c: int, n_h: int) -> float:
    """Atomisation energy of one molecule [J], from the 0 K heats of formation above.

    `D_at = n_C dHf(C) + n_H dHf(H) - dHf(molecule)`, all at 0 K, converted per molecule.
    """
    return (n_c * _HF0_C + n_h * _HF0_H - _HF0[name]) * 1.0e3 / N_A


@dataclass(frozen=True)
class Molecule:
    """A gas-phase `C_a H_b` species and the equilibrium `M <=> a C + b H` that sets its density.

    Generalises `eos_water.Diatomic` in the two ways methane needs: an arbitrary stoichiometry
    (so `CH4` and `C2H2` are the same code path as `H2`), and a rotor that may be linear (one
    `theta_rot`) or nonlinear (three). Degenerate vibrational modes are listed once per degree of
    freedom rather than carrying a multiplicity field, so the mode count is visibly `3N-6`/`3N-5`.
    """

    name: str
    n_c: int
    n_h: int
    mass: float  # molecular mass [kg]
    degeneracy: float  # ground-term electronic degeneracy
    d_at: float  # atomisation energy from v = 0 [J per molecule]
    theta_rot: tuple[float, ...]  # 1 entry if linear, 3 if nonlinear [K]
    symmetry: float  # rotational symmetry number
    theta_vib: tuple[float, ...]  # one entry per vibrational degree of freedom [K]

    @property
    def linear(self) -> bool:
        """Whether the rotor is linear -- decided by the constant count, not by a second flag
        that could disagree with it."""
        return len(self.theta_rot) == 1


def _molecule(
    name: str,
    n_c: int,
    n_h: int,
    degeneracy: float,
    rot_cm: Sequence[float],
    symmetry: float,
    vib_cm: Sequence[float],
) -> Molecule:
    """Build a `Molecule` from spectroscopic constants in cm^-1, deriving mass and `d_at`."""
    return Molecule(
        name=name,
        n_c=n_c,
        n_h=n_h,
        mass=n_c * M_C + n_h * M_H,
        degeneracy=degeneracy,
        d_at=_atomization(name, n_c, n_h),
        theta_rot=tuple(_CM_TO_K * v for v in rot_cm),
        symmetry=symmetry,
        theta_vib=tuple(_CM_TO_K * v for v in vib_cm),
    )


# H2: X 1Sigma_g+. Same constants as `eos_water.H2`, restated here so this species set is one
# table rather than half an import -- the values are checked equal in the tests.
H2 = _molecule("H2", 0, 2, 1.0, (60.853,), 2.0, (4401.21,))
# CH: X 2Pi -- Lambda-doubling (2) x spin (2).
CH = _molecule("CH", 1, 1, 4.0, (14.457,), 1.0, (2858.5,))
# C2: X 1Sigma_g+. The low-lying a 3Pi_u at 716 cm^-1 is *not* summed (ground terms only), which
# understates C2 in a band where it never exceeds a few percent of the carbon.
C2 = _molecule("C2", 2, 0, 1.0, (1.8198,), 2.0, (1854.71,))
# CH2: X 3B1, bent asymmetric top.
CH2 = _molecule("CH2", 1, 2, 3.0, (73.81, 8.45, 7.18), 2.0, (2806.0, 963.0, 3190.0))
# CH3: X 2A2'', planar oblate symmetric top (D3h, sigma = 6). Modes: nu1, nu2 (umbrella),
# nu3 (x2), nu4 (x2).
CH3 = _molecule(
    "CH3", 1, 3, 2.0, (9.578, 9.578, 4.789), 6.0, (3004.0, 606.0, 3161.0, 3161.0, 1396.0, 1396.0)
)
# CH4: spherical top (Td, sigma = 12). Modes: nu1, nu2 (x2), nu3 (x3), nu4 (x3) -- nine in all.
CH4 = _molecule(
    "CH4",
    1,
    4,
    1.0,
    (5.2412, 5.2412, 5.2412),
    12.0,
    (2917.0, 1534.0, 1534.0, 3019.0, 3019.0, 3019.0, 1306.0, 1306.0, 1306.0),
)
# C2H2: linear (D_inf_h, sigma = 2). Modes: nu1, nu2, nu3, nu4 (x2, bend), nu5 (x2, bend).
C2H2 = _molecule(
    "C2H2", 2, 2, 1.0, (1.1766,), 2.0, (3374.0, 1974.0, 3289.0, 612.0, 612.0, 730.0, 730.0)
)
# C2H: X 2Sigma+, linear. The soot pathway's first radical; carried because it is the one C-H
# species that survives where C2H2 has gone.
C2H = _molecule("C2H", 2, 1, 2.0, (1.4568,), 1.0, (3329.0, 1841.0, 372.0, 372.0))
# C3: linear (D_inf_h, sigma = 2). Its nu2 bend at 63 cm^-1 is very floppy, so the harmonic
# treatment overstates C3 somewhat at the cold end; it is a minor species there.
C3 = _molecule("C3", 3, 0, 1.0, (0.4305,), 2.0, (1224.0, 63.0, 63.0, 2040.0))

#: Every molecule, in the order the mass-action vector and the energy sum both use.
MOLECULES = (CH4, CH3, CH2, CH, C2H2, C2H, C2, C3, H2)

#: Atomisation energy of one kilogram of methane [J/kg] -- 102.35 MJ/kg at the 0 K reference.
#: The ceiling on what the dissociation store can hold or strand.
FULL_ATOMIZATION_ENERGY = CH4.d_at / M_CH4

#: What full recombination returns, split by channel, per kilogram of methane [J/kg]. `4 H -> 2 H2`
#: is a three-body gas-phase reaction and gets the `n^2` help a walled chamber provides;
#: `C(g) -> C(graphite)` is a *nucleation* problem and gets none of it (N10 item 5). Their sum is
#: below `FULL_ATOMIZATION_ENERGY` by methane's own formation enthalpy, which never returns.
H2_RETURN_ENERGY = 2.0 * H2.d_at / M_CH4
CARBON_CONDENSATION_ENERGY = _HF0_C * 1.0e3 / N_A / M_CH4
FORMATION_ENTHALPY = FULL_ATOMIZATION_ENERGY - H2_RETURN_ENERGY - CARBON_CONDENSATION_ENERGY


@dataclass(frozen=True)
class Composition:
    """Equilibrium number densities [m^-3] at a single `(rho, T)`.

    `n_mol[i]` is the density of `MOLECULES[i]`; `n_c_ions[k]` is the C ion of charge `k+1`.
    """

    n_c: float
    n_h: float
    n_hp: float
    n_c_ions: tuple[float, ...]
    n_e: float
    n_mol: tuple[float, ...]

    def density_of(self, mol: Molecule) -> float:
        """Density [m^-3] of a named molecule, by identity rather than by remembered index."""
        return self.n_mol[MOLECULES.index(mol)]

    @property
    def n_neutral_heavy(self) -> float:
        """Every neutral heavy particle [m^-3] -- the third-body density for three-body
        recombination, and the neutral target density for electron-neutral collisions."""
        return self.n_c + self.n_h + sum(self.n_mol)

    @property
    def n_total(self) -> float:
        """Every particle [m^-3], electrons included -- the ideal-mixture pressure count."""
        return self.n_neutral_heavy + self.n_hp + sum(self.n_c_ions) + self.n_e

    @property
    def n_h_nuclei(self) -> float:
        """Hydrogen nuclei per m^3, however bound."""
        return (
            self.n_h
            + self.n_hp
            + sum(m.n_h * n for m, n in zip(MOLECULES, self.n_mol, strict=True))
        )

    @property
    def n_c_nuclei(self) -> float:
        """Carbon nuclei per m^3, however bound."""
        return (
            self.n_c
            + sum(self.n_c_ions)
            + sum(m.n_c * n for m, n in zip(MOLECULES, self.n_mol, strict=True))
        )

    @property
    def hydrogen_dissociated(self) -> float:
        """Fraction of hydrogen nuclei that are *free* -- atomic H plus H+.

        **This is the number ask N10 item 4b turns on.** It is not `1 - n_H2/n_f`: a nucleus
        sitting in CH4 or C2H2 is bound too, and at chamber conditions those are gone, but at the
        cold end they are not.
        """
        total = self.n_h_nuclei
        return 0.0 if total <= 0.0 else (self.n_h + self.n_hp) / total

    @property
    def carbon_free(self) -> float:
        """Fraction of carbon nuclei that are free atoms or atomic ions -- what carbon
        condensation would have to collect, and what N10 item 5 asks about."""
        total = self.n_c_nuclei
        return 0.0 if total <= 0.0 else (self.n_c + sum(self.n_c_ions)) / total


def _z_rot(mol: Molecule, temp: float) -> float:
    """Classical rigid-rotor partition function, linear or nonlinear.

    Valid for `T >> theta_rot`; the stiffest rotor here is H2 at 87.6 K, so everything above a
    few hundred kelvin is classical to better than a percent.
    """
    if mol.linear:
        return float(temp / (mol.symmetry * mol.theta_rot[0]))
    prod = mol.theta_rot[0] * mol.theta_rot[1] * mol.theta_rot[2]
    return float(math.sqrt(math.pi) / mol.symmetry * math.sqrt(temp**3 / prod))


def _e_rot(mol: Molecule, temp: float) -> float:
    """Mean rotational energy of one molecule [J]: `kT` linear, `3/2 kT` nonlinear."""
    return (1.0 if mol.linear else 1.5) * K_B * temp


def _z_vib(mol: Molecule, temp: float) -> float:
    """Harmonic-oscillator vibrational partition function, referenced to `v = 0`."""
    z = 1.0
    for theta in mol.theta_vib:
        z /= 1.0 - math.exp(-min(theta / temp, _EXP_CAP))
    return z


def _e_vib(mol: Molecule, temp: float) -> float:
    """Mean vibrational energy of one molecule [J], excluding zero-point."""
    e = 0.0
    for theta in mol.theta_vib:
        x = theta / temp
        e += theta * K_B / math.expm1(min(x, _EXP_CAP))
    return e


def ln_k_molecule(mol: Molecule, temp: float) -> float:
    """ln of `K = n_C^a n_H^b / n_M` for the atomisation `M <=> a C + b H`.

    Units are `ln m^{-3(a+b-1)}`, which is why the residual works in log space throughout: at
    500 K methane's constant is around `exp(-380)` and nothing linear survives it.
    """
    ln_c = _ln_trans(M_C, temp) + math.log(G_C)
    ln_h = _ln_trans(M_H, temp) + math.log(G_H)
    ln_m = _ln_trans(mol.mass, temp) + math.log(
        mol.degeneracy * _z_rot(mol, temp) * _z_vib(mol, temp)
    )
    return mol.n_c * ln_c + mol.n_h * ln_h - ln_m - mol.d_at / (K_B * temp)


def _ln_kc_ladder(temp: float) -> tuple[float, ...]:
    """ln Saha constants for every C stage `C^k <=> C^{k+1} + e-`."""
    return tuple(
        ln_k_saha(IP_C_LADDER[k], G_C_LADDER[k + 1], G_C_LADDER[k], temp) for k in range(N_C_STAGES)
    )


@dataclass(frozen=True)
class _EqConstants:
    """Every equilibrium constant at one temperature, in log form."""

    ln_kh: float  # H <=> H+ + e-
    ln_kc: tuple[float, ...]  # the C ionisation ladder
    ln_k_mol: tuple[float, ...]  # one per entry of MOLECULES, in that order


def _constants(temp: float) -> _EqConstants:
    """All equilibrium constants at `temp`."""
    return _EqConstants(
        ln_kh=ln_k_saha(IP_H, 1.0, G_H, temp),  # H+ is a bare proton: g = 1
        ln_kc=_ln_kc_ladder(temp),
        ln_k_mol=tuple(ln_k_molecule(m, temp) for m in MOLECULES),
    )


def _ln_species(x: Vec, c: _EqConstants) -> tuple[list[float], list[float], float, float]:
    """Log-densities of every derived species from `x = [ln n_C, ln n_H, ln n_e]`.

    Returns `(ln molecules, ln C-ion stages, ln H+, ln n_e)`. Kept in log form because the cold
    limit runs to `exp(-700)` and the hot limit to `exp(+60)` in the same solve.
    """
    ln_nc, ln_nh, ln_ne = float(x[0]), float(x[1]), float(x[2])
    ln_mol = [
        m.n_c * ln_nc + m.n_h * ln_nh - lk for m, lk in zip(MOLECULES, c.ln_k_mol, strict=True)
    ]
    ln_stages: list[float] = []
    ln_stage = ln_nc
    for k in range(N_C_STAGES):
        ln_stage = ln_stage + c.ln_kc[k] - ln_ne
        ln_stages.append(ln_stage)
    return ln_mol, ln_stages, c.ln_kh + ln_nh - ln_ne, ln_ne


def _exp(ln_n: float) -> float:
    """`exp`, clamped so a far Newton iterate cannot overflow before the step cap catches it."""
    return float(math.exp(min(ln_n, _EXP_CAP)))


def _densities(x: Vec, c: _EqConstants) -> Composition:
    """Reconstruct the full composition from `x = [ln n_C, ln n_H, ln n_e]`."""
    ln_mol, ln_stages, ln_hp, ln_ne = _ln_species(x, c)
    return Composition(
        n_c=_exp(float(x[0])),
        n_h=_exp(float(x[1])),
        n_hp=_exp(ln_hp),
        n_c_ions=tuple(_exp(t) for t in ln_stages),
        n_e=_exp(ln_ne),
        n_mol=tuple(_exp(t) for t in ln_mol),
    )


_LN_CHARGE = tuple(math.log(k + 1.0) for k in range(N_C_STAGES))
#: `ln n_c` and `ln n_h` coefficients of each molecule, precomputed: the residual is evaluated
#: four times per Newton iteration for the finite-difference Jacobian.
_MOL_NC = tuple(float(m.n_c) for m in MOLECULES)
_MOL_NH = tuple(float(m.n_h) for m in MOLECULES)
_LN_MOL_NC = tuple(math.log(m.n_c) if m.n_c else -math.inf for m in MOLECULES)
_LN_MOL_NH = tuple(math.log(m.n_h) if m.n_h else -math.inf for m in MOLECULES)
_I_CH4 = MOLECULES.index(CH4)
_LN4 = math.log(4.0)


def _ln_sum(terms: Sequence[float]) -> float:
    """`ln(sum_i exp(terms_i))` without leaving log space."""
    top = max(terms)
    if top == -math.inf:
        return -math.inf
    return top + math.log(sum(math.exp(t - top) for t in terms))


def _residual(x: Vec, c: _EqConstants, n_f: float) -> Vec:
    """Equilibrium residual at `x = [ln n_C, ln n_H, ln n_e]`.

    Three independent, well-conditioned constraints, chosen the same way `eos_water`'s are:

    - **C nuclei** conservation, counting the carbon in every molecule and every ion stage. This
      is the one that pins `n_CH4 ~ n_f` in the cold limit. **In log form**, unlike
      `eos_water`'s: a far Newton iterate here can put `ln n_C` a hundred decades out (the
      methane constants run to `exp(-180)` at 500 K, against water's `exp(-80)`), and the linear
      ratio overflows to `inf` before the step cap can pull it back. The log form is bounded and
      homogeneous with the other two, which is also what makes the LM Jacobian well-scaled.
    - **H:C = 4:1** as a ratio over the non-`CH4` species only, **in log form**. `CH4` is 4:1 by
      construction, so the stoichiometry holds if and only if it holds with `CH4` removed from
      both sides -- an exact restatement, and the one that stays conditioned when `CH4` carries
      essentially all the nuclei and the naive pair degenerates into "n_CH4 = n_f" twice.
    - **charge neutrality**, log form, each ion term carrying its own `-j ln n_e` Saha chain.
    """
    ln_nc, _, ln_ne = float(x[0]), float(x[1]), float(x[2])
    ln_mol, ln_stages, ln_hp, _ = _ln_species(x, c)

    ln_charge = _ln_sum([ln_hp, *(_LN_CHARGE[k] + t for k, t in enumerate(ln_stages))])

    # Nuclei carried by everything except CH4, in log form.
    ln_c_free = _ln_sum(
        [
            ln_nc,
            *ln_stages,
            *(
                _LN_MOL_NC[i] + ln_mol[i]
                for i in range(len(MOLECULES))
                if i != _I_CH4 and _MOL_NC[i]
            ),
        ]
    )
    ln_h_free = _ln_sum(
        [
            float(x[1]),
            ln_hp,
            *(
                _LN_MOL_NH[i] + ln_mol[i]
                for i in range(len(MOLECULES))
                if i != _I_CH4 and _MOL_NH[i]
            ),
        ]
    )

    ln_c_total = _ln_sum([ln_mol[_I_CH4], ln_c_free])
    return np.array(
        [
            ln_c_total - math.log(n_f),  # C nuclei, in log form
            ln_h_free - ln_c_free - _LN4,  # H:C, CH4's exact 4:1 divided out
            ln_ne - ln_charge,  # charge neutrality
        ]
    )


def _cold_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Molecular init: `CH4 <=> C + 4H` with `n_H = 4 n_C`, so `256 n_C^5 = K n_f`."""
    ln_k = c.ln_k_mol[_I_CH4]
    ln_nc = min(ln_nf, (ln_k + ln_nf - math.log(256.0)) / 5.0)
    ln_nh = min(_LN4 + ln_nf, _LN4 + ln_nc)
    ln_ne = min(
        0.5 * float(np.logaddexp(c.ln_kh + ln_nh, c.ln_kc[0] + ln_nc)), math.log(10.0) + ln_nf
    )
    return np.array([ln_nc, ln_nh, ln_ne])


def _methane_dominant_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Init for the cold limit, where `CH4` holds essentially every nucleus.

    `_cold_init` fails here and the failure is instructive. It assumes the trace atom pool sits
    at the same 4:1 as the parent, but it does not: the free hydrogen goes to `H2` (2 nuclei per
    molecule) and the free carbon to `C2H2` (2 per molecule, carrying 2 hydrogen with it), and
    the stoichiometry constraint is on *that* pool. Imposing it there gives a carbon density that
    depends on the equilibrium constants alone,

        `n_H2 = 4 n_C2H2`  =>  `n_C^2 = K_C2H2 / (4 K_H2)`,

    and `n_H` then follows from `n_CH4 = n_f`. At 500 K the two differ by thirty decades in
    `ln n_C`, which no amount of Newton damping recovers from.
    """
    i_h2, i_c2h2 = MOLECULES.index(H2), MOLECULES.index(C2H2)
    ln_nc = 0.5 * (c.ln_k_mol[i_c2h2] - c.ln_k_mol[i_h2] - _LN4)
    ln_nh = (ln_nf + c.ln_k_mol[_I_CH4] - ln_nc) / 4.0
    ln_ne = min(
        0.5 * float(np.logaddexp(c.ln_kh + ln_nh, c.ln_kc[0] + ln_nc)), math.log(10.0) + ln_nf
    )
    return np.array([ln_nc, ln_nh, ln_ne])


def _pyrolysis_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Init for the band where methane has cracked but hydrogen has not: `H2` + free carbon.

    The cold init sends CH4 straight to atoms, so between roughly 1500 and 5000 K -- where the
    C-H bonds are gone and the H-H bond is not -- it can start decades off in `ln n_H`. Here the
    hydrogen is taken to be all `H2` at `2 n_f`, giving `n_H^2 = K_H2 * 2 n_f`, and the carbon to
    be free at `n_f`.
    """
    i_h2 = MOLECULES.index(H2)
    ln_nh = min(_LN4 + ln_nf, 0.5 * (c.ln_k_mol[i_h2] + math.log(2.0) + ln_nf))
    ln_nc = ln_nf
    ln_ne = min(
        0.5 * float(np.logaddexp(c.ln_kh + ln_nh, c.ln_kc[0] + ln_nc)), math.log(10.0) + ln_nf
    )
    return np.array([ln_nc, ln_nh, ln_ne])


def _hot_init(ln_nf: float, c: _EqConstants) -> Vec:
    """Hot-plasma init: a mean-charge fixed point on the Saha ladder, as `eos_water._hot_init`."""
    q = 1.0
    k_dom = 0
    f_hp = 0.5
    for _ in range(30):
        ln_ne = math.log(q) + ln_nf
        k_dom = sum(1 for lk in c.ln_kc if lk > ln_ne)
        f_hp = float(np.exp(c.ln_kh - np.logaddexp(c.ln_kh, ln_ne)))
        q_new = max(4.0 * f_hp + float(k_dom), 1e-3)
        if abs(q_new - q) < 1e-3 * q:
            q = q_new
            break
        q = 0.5 * (q + q_new)
    ln_ne = math.log(q) + ln_nf
    ln_nh = _LN4 + ln_nf + math.log(max(1.0 - f_hp, 1e-12))
    ln_nc = ln_nf - sum(c.ln_kc[j] - ln_ne for j in range(k_dom))
    return np.array([ln_nc, ln_nh, ln_ne])


def _newton_polish(x0: Vec, c: _EqConstants, n_f: float) -> tuple[Vec, float]:
    """LM-damped Newton in log space from `x0`; returns `(x, max |residual|)`."""
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
        jtj = jac.T @ jac
        delta = np.linalg.solve(jtj + lam * np.eye(3), -jac.T @ r)
        step = float(np.max(np.abs(delta)))
        if step > 2.0:  # cap the log-density step so a far init cannot overshoot into overflow
            delta *= 2.0 / step
        x = x + delta
        r = _residual(x, c, n_f)
    return x, float(np.max(np.abs(r)))


def _solve_log_densities(temp: float, n_f: float) -> Vec:
    """Solve `[ln n_C, ln n_H, ln n_e]` by Newton from each of three physics-based inits.

    The four bracket the regimes this EOS spans -- bound methane at the density where its own
    trace pool sets the atom split, bound methane at the density where it does not, the
    cracked-but-molecular hydrogen band, and the stripped plasma -- and the first that converges
    wins.
    """
    c = _constants(temp)
    ln_nf = math.log(n_f)

    best: tuple[Vec, float] | None = None
    inits = (
        _cold_init(ln_nf, c),
        _methane_dominant_init(ln_nf, c),
        _pyrolysis_init(ln_nf, c),
        _hot_init(ln_nf, c),
    )
    for init in inits:
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
    n_f = rho / M_CH4  # CH4 formula units per m^3 (conserves C:H = 1:4)
    return _densities(_solve_log_densities(temp, n_f), _constants(temp))


def _bond_energy_held(comp: Composition, n_f: float) -> float:
    """Chemical energy stored as broken bonds [J/m^3], net of what the molecules pay back.

    Fully atomising the charge costs `n_f D_at(CH4)`; every molecule present, `CH4` included,
    hands back its own atomisation energy. Intact methane therefore owes nothing, and a gas that
    has cracked to `C + 2 H2` still holds the C-H bonds' worth that the H-H bonds do not cover --
    which is the store the walled nozzle is trying to charge and then get back.
    """
    held = n_f * CH4.d_at
    for mol, n in zip(MOLECULES, comp.n_mol, strict=True):
        held -= n * mol.d_at
    return held


def bond_energy_held(rho: float, temp: float) -> float:
    """Chemical energy still locked in broken methane bonds [J/kg] at equilibrium `(rho, temp)`."""
    return _bond_energy_held(composition(rho, temp), rho / M_CH4) / rho


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """Equilibrium `(p [Pa], e [J/kg])` at `(rho, temp)`.

    `p = (sum_i n_i) k T` (ideal mixture, electrons included). `e` is translational + rotational
    + vibrational thermal energy + chemical (dissociation + ionisation) energy, referenced to
    bound `CH4` = 0 so `e > 0` everywhere.
    """
    comp = composition(rho, temp)
    n_f = rho / M_CH4

    p = comp.n_total * K_B * temp

    n_monatomic = comp.n_c + comp.n_h + comp.n_hp + sum(comp.n_c_ions) + comp.n_e
    e_thermal = 1.5 * K_B * temp * n_monatomic
    for mol, n in zip(MOLECULES, comp.n_mol, strict=True):
        e_thermal += n * (1.5 * K_B * temp + _e_rot(mol, temp) + _e_vib(mol, temp))

    e_chem = comp.n_hp * IP_H
    e_chem += sum(n * E_C_CUM[k] for k, n in enumerate(comp.n_c_ions))
    e_chem += _bond_energy_held(comp, n_f)

    return p, (e_thermal + e_chem) / rho


def sound_speed(rho: float, temp: float) -> float:
    """Equilibrium adiabatic sound speed `c_s` [m/s] (see `eos_water.sound_speed_fd`)."""
    return sound_speed_fd(pressure_energy, rho, temp)
