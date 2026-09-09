"""Conditional neutral-network relaxation, not validated dense-water kinetics.

Forward rates: Cantera v3.2.0 h2o2.yaml (GRI-Mech 3.0 submechanism).
Reverse rates are matched to this study's partition functions, NOT its NASA7
fits. Ions are fixed spectators; HO2/H2O2 routes are absent. No pressure
falloff is supplied for the retained association reactions. This is a local,
isothermal, fixed-density linear response, not a finite-rate hydro solve.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat import eos_water as ew

Vec = NDArray[np.float64]
SPECIES = ("H", "O", "H2", "OH", "O2", "H2O")
SOURCE_URL = "https://raw.githubusercontent.com/Cantera/cantera/v3.2.0/data/h2o2.yaml"
SOURCE_SHA256 = "0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4"
REFERENCE_THERMO_MAX_K = 3500.0  # NOT a validated range of every forward rate.
ATOM_COUNTS = np.array([[1, 0, 2, 1, 0, 2], [0, 1, 0, 1, 2, 1]], dtype=float)
GROUND_ENERGIES = np.array(
    [0.0, ew.D_AT, -ew.H2.d0, ew.D_AT - ew.OH.d0, 2 * ew.D_AT - ew.O2.d0, 0.0]
)


@dataclass(frozen=True)
class Reaction:
    source_number: int
    reactants: tuple[int, ...]
    products: tuple[int, ...]
    a_cgs_mol: float
    exponent: float
    activation_cal_mol: float
    # None is bimolecular; otherwise effective collider densities in SPECIES order.
    colliders: tuple[float, ...] | None = None

    def coefficient(self, temp: float) -> float:
        order = sum(self.reactants) + (self.colliders is not None)
        return float(
            self.a_cgs_mol
            * (1e-6 / ew.N_A) ** (order - 1)
            * temp**self.exponent
            * math.exp(-self.activation_cal_mol * 4.184 / (ew.N_A * ew.K_B * temp))
        )

    @property
    def stoichiometry(self) -> Vec:
        return np.array(
            [p - r for p, r in zip(self.products, self.reactants, strict=True)], dtype=np.float64
        )

    def forward_flux(self, temp: float, densities: Vec) -> float:
        value = self.coefficient(temp)
        for n, order in zip(densities, self.reactants, strict=True):
            value *= float(n) ** order
        if self.colliders is not None:
            value *= float(np.dot(self.colliders, densities))
        return value


# Integer stoichiometry is net of the spectator, including explicit H2/H2O
# collider channels in reactions 13/14. No 1/2 factor: GRI's rate convention
# for 2 H is q=k[H]^2[M], with two H atoms consumed per event.
REACTIONS = (
    Reaction(1, (0, 2, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0), 1.2e17, -1, 0, (1, 1, 2.4, 1, 1, 15.4)),
    Reaction(2, (1, 1, 0, 0, 0, 0), (0, 0, 0, 1, 0, 0), 5e17, -1, 0, (1, 1, 2, 1, 1, 6)),
    Reaction(3, (0, 1, 1, 0, 0, 0), (1, 0, 0, 1, 0, 0), 3.87e4, 2.7, 6260),
    Reaction(11, (1, 0, 0, 0, 1, 0), (0, 1, 0, 1, 0, 0), 2.65e16, -0.6707, 17041),
    Reaction(12, (2, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0), 1e18, -1, 0, (1, 1, 0, 1, 1, 0)),
    Reaction(13, (2, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0), 9e16, -0.6, 0, (0, 0, 1, 0, 0, 0)),
    Reaction(14, (2, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0), 6e19, -1.25, 0, (0, 0, 0, 0, 0, 1)),
    Reaction(15, (1, 0, 0, 1, 0, 0), (0, 0, 0, 0, 0, 1), 2.2e22, -2, 0, (1, 1, 0.73, 1, 1, 3.65)),
    Reaction(21, (0, 0, 1, 1, 0, 0), (1, 0, 0, 0, 0, 1), 2.16e8, 1.51, 3430),
    Reaction(23, (0, 0, 0, 2, 0, 0), (0, 1, 0, 0, 0, 1), 3.57e4, 2.4, -2110),
)


def neutral_densities(comp: ew.Composition) -> Vec:
    return np.array([comp.n_h, comp.n_o, comp.n_h2, comp.n_oh, comp.n_o2, comp.n_h2o])


def log_equilibrium_constant(reaction: Reaction, temp: float) -> float:
    """Net concentration K in molecule/m³ units, from the SAME gas partition functions."""
    # At arbitrary positive atomic reference densities n_H=n_O=1, molecular
    # log densities follow mass action. Element conservation makes their
    # stoichiometric dot product independent of these two arbitrary references.
    return float(np.dot(reaction.stoichiometry, equilibrium_log_reference(temp)))


def equilibrium_log_reference(temp: float) -> Vec:
    return np.array(
        [
            0.0,
            0.0,
            -ew.ln_k_diatomic(ew.H2, temp),
            -ew.ln_k_diatomic(ew.OH, temp),
            -ew.ln_k_diatomic(ew.O2, temp),
            -ew.ln_k_dissoc(temp),
        ]
    )


def net_production(temp: float, densities: Vec) -> Vec:
    """Matched reversible rates for Jacobian tests; not a time integrator."""
    deviation = np.log(densities) - equilibrium_log_reference(temp)
    result = np.zeros(6)
    for reaction in REACTIONS:
        nu = reaction.stoichiometry
        result -= nu * reaction.forward_flux(temp, densities) * math.expm1(float(nu @ deviation))
    return result


@dataclass(frozen=True)
class Relaxation:
    energy_time_s: float
    mode_rates_s: tuple[float, ...]
    energy_weights: tuple[float, ...]
    max_detailed_balance_log_error: float
    chemical_energy_susceptibility_j2_m3: float
    water_collider_efficiency: float


def relaxation(
    temp: float,
    densities: Vec,
    ground_energies: Vec = GROUND_ENERGIES,
    rate_scale: float = 1.0,
    *,
    water_association_scale: float = 1.0,
) -> Relaxation:
    """Bond-energy-weighted integrated relaxation time over four reactive modes.

    At detailed balance, the symmetrized negative Jacobian is
    B = sum q_eq * (nu/sqrt(n)) (nu/sqrt(n))^T. Remove the two elemental null
    modes using an independent stoichiometric basis. Project chemical energy
    onto its eigenmodes, then tau_E = sum(w_i/lambda_i)/sum(w_i).
    No reaction progress, heat release or impulse is inferred from equilibrium
    one-way fluxes (whose forward and reverse values cancel).
    water_association_scale multiplies BOTH directions of reaction 15 only.
    Zero removes that channel, retaining the other nine reactions. It is a
    sensitivity endpoint, not a physical high-pressure rate prescription.
    """
    if not (1000 <= temp <= 2e6 and np.all(np.isfinite(densities)) and np.all(densities > 0)):
        raise ValueError("requires positive finite populations and 1000 K <= T <= 2 MK")
    if not math.isfinite(rate_scale) or rate_scale <= 0:
        raise ValueError("rate scale must be finite and positive")
    if not math.isfinite(water_association_scale) or water_association_scale < 0:
        raise ValueError("water association scale must be finite and nonnegative")
    root_n = np.sqrt(densities)
    deviation = np.log(densities) - equilibrium_log_reference(temp)
    matrix = np.zeros((6, 6))
    error = 0.0
    for reaction in REACTIONS:
        nu = reaction.stoichiometry
        error = max(error, abs(float(nu @ deviation)))
        direction = nu / root_n
        channel_scale = water_association_scale if reaction.source_number == 15 else 1.0
        matrix += (
            rate_scale
            * channel_scale
            * reaction.forward_flux(temp, densities)
            * np.outer(direction, direction)
        )
    if error > 1e-7:
        raise ValueError(f"not an equilibrium state of the matched neutral network: {error}")
    basis = np.column_stack([REACTIONS[i].stoichiometry / root_n for i in (0, 1, 4, 7)])
    q, _ = np.linalg.qr(basis, mode="reduced")
    reduced = q.T @ matrix @ q
    scale = float(np.max(np.abs(reduced)))
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("nonfinite or zero relaxation matrix")
    rates, vectors = np.linalg.eigh(0.5 * (reduced + reduced.T) / scale)
    if rates[0] <= 0 or rates[-1] / rates[0] > 1e12:
        raise ValueError("relaxation spectrum not numerically resolved")
    rates *= scale
    weights = (vectors.T @ q.T @ (ground_energies * root_n)) ** 2
    total = float(np.sum(weights))
    if not math.isfinite(total) or total <= 0:
        raise ValueError("no resolved chemical-energy susceptibility")
    weights /= total
    water = next(r for r in REACTIONS if r.source_number == 15)
    assert water.colliders is not None
    return Relaxation(
        float(np.sum(weights / rates)),
        tuple(float(x) for x in rates),
        tuple(float(x) for x in weights),
        error,
        total,
        float(np.dot(water.colliders, densities)) / float(np.sum(densities)),
    )
