"""Rate control and energy-constrained LOCAL linear chemistry, not finite-rate flow.

Same ten conditional neutral rates/partition functions as kinetics.py. Ions
stay chemically fixed but carry sensible heat. Fixed volume and total energy
permit a temperature perturbation; no expansion or nonlinear ODE is integrated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate import flow_eos
from puffsat.water_plate import kinetics as kn
from puffsat.water_plate.kinetics import Vec


def species_thermo(temp: float) -> tuple[Vec, Vec]:
    """Per-neutral internal energy and frozen cv in the existing elemental gauge."""
    if not math.isfinite(temp) or not 1000 <= temp <= 2e6:
        raise ValueError("thermal screen requires 1000 K <= T <= 2 MK")
    cv = ew.K_B * np.array([1.5, 1.5, 2.5, 2.5, 2.5, 3.0])
    energy = kn.GROUND_ENERGIES + cv * temp
    for i, modes in (
        (2, (ew.H2.theta_vib,)),
        (3, (ew.OH.theta_vib,)),
        (4, (ew.O2.theta_vib,)),
        (5, ew._THETA_VIB),
    ):
        for theta in modes:
            x = theta / temp
            denominator = math.expm1(x)
            energy[i] += ew.K_B * theta / denominator
            cv[i] += ew.K_B * x**2 * math.exp(x) / denominator**2
    return energy, cv


def frozen_cv(comp: ew.Composition, neutral_cv: Vec, residual_cv_j_m3_k: float) -> float:
    """Volumetric heat capacity with ALL populations fixed, including ion spectators."""
    spectators = comp.n_hp + sum(comp.n_o_ions) + comp.n_e
    return (
        float(kn.neutral_densities(comp) @ neutral_cv)
        + 1.5 * ew.K_B * spectators
        + residual_cv_j_m3_k
    )


@dataclass(frozen=True)
class CvCurve:
    """Same log-density cubic-Hermite residual-cv continuation as Step 6."""

    log_rho: Vec
    values: Vec
    slopes: Vec

    @classmethod
    def build(cls, points: int = 512) -> CvCurve:
        if points < 2:
            raise ValueError("need at least two residual knots")
        rho = np.geomspace(1e-5, 2000, points)
        refs = [flow_eos.residual_reference(float(r)) for r in rho]
        return cls(
            np.log(rho),
            np.array([ref.cv for ref in refs]),
            np.array(
                [-1000 * ref.pressure_tt / float(r) for ref, r in zip(refs, rho, strict=True)]
            ),
        )

    def reference_cv(self, rho: float) -> float:
        if not math.isfinite(rho) or rho <= 0:
            raise ValueError("density must be finite and positive")
        x = math.log(rho)
        if not self.log_rho[0] <= x <= self.log_rho[-1]:
            raise ValueError("density outside residual curve")
        i = min(max(int(np.searchsorted(self.log_rho, x)) - 1, 0), len(self.log_rho) - 2)
        width = float(self.log_rho[i + 1] - self.log_rho[i])
        s = (x - float(self.log_rho[i])) / width
        return float(
            (2 * s**3 - 3 * s**2 + 1) * self.values[i]
            + (s**3 - 2 * s**2 + s) * width * self.slopes[i]
            + (-2 * s**3 + 3 * s**2) * self.values[i + 1]
            + (s**3 - s**2) * width * self.slopes[i + 1]
        )


@dataclass(frozen=True)
class Response:
    isothermal_s: float
    isoenergetic_s: float
    isothermal_sensitivities: tuple[float, ...]
    isoenergetic_sensitivities: tuple[float, ...]
    thermal_metric_strength: float
    susceptibility_ratio: float


def _relaxation(matrix: Vec, directions: Vec, fluxes: Vec, observable: Vec) -> tuple[float, Vec]:
    scale = float(np.max(np.abs(matrix)))
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("invalid relaxation matrix")
    rates, modes = np.linalg.eigh((matrix + matrix.T) / (2 * scale))
    if rates[0] <= 0 or rates[-1] / rates[0] > 1e12:
        raise ValueError("rank-deficient or unresolved relaxation spectrum")
    solved = modes @ ((modes.T @ observable) / rates) / scale
    numerator = float(observable @ solved)
    denominator = float(observable @ observable)
    if not math.isfinite(numerator) or min(numerator, denominator) <= 0:
        raise ValueError("unresolved energy response")
    sensitivities = fluxes * (directions.T @ solved) ** 2 / numerator
    if abs(float(np.sum(sensitivities)) - 1) > 1e-7:
        raise ValueError("unresolved rate-sensitivity sum")
    return numerator / denominator, sensitivities


def response(
    temp: float,
    densities: Vec,
    internal_energies: Vec,
    cv_j_m3_k: float,
    rate_scales: Vec | None = None,
) -> Response:
    """Matched iso-T and iso-E bond-energy autocorrelation relaxation times.

    For y=Q^T delta(n)/sqrt(n), energy conservation gives
    delta(T)=-a^T y/Cv, with a=Q^T sqrt(n) u. Van't Hoff gives metric
    H=I+a*a^T/(kB*T²*Cv) and dy/dt=-B*H*y. Symmetrize as
    C=sqrt(H)*B*sqrt(H); fixed-energy covariance is H^-1. The observable
    is therefore g=H^-1/2 Q^T sqrt(n) e0, not the unchanged iso-T weights.
    Rate sensitivities are -d ln(tau)/d ln(q_r), nonnegative and summing to one.
    """
    if not math.isfinite(cv_j_m3_k) or cv_j_m3_k <= 0:
        raise ValueError("need positive finite fixed-composition heat capacity")
    if not 1000 <= temp <= 2e6 or densities.shape != (6,) or not np.all(densities > 0):
        raise ValueError("invalid temperature or populations")
    if (
        internal_energies.shape != (6,)
        or not np.all(np.isfinite(densities))
        or not np.all(np.isfinite(internal_energies))
    ):
        raise ValueError("need finite species inputs")
    scales = np.ones(len(kn.REACTIONS)) if rate_scales is None else rate_scales
    if scales.shape != (10,) or not np.all(np.isfinite(scales)) or np.any(scales < 0):
        raise ValueError("need ten finite nonnegative rate scales")
    root_n = np.sqrt(densities)
    nu = np.column_stack([r.stoichiometry for r in kn.REACTIONS])
    if np.max(np.abs(nu.T @ (np.log(densities) - kn.equilibrium_log_reference(temp)))) > 1e-7:
        raise ValueError("not an equilibrium state")
    basis, _ = np.linalg.qr(nu[:, [0, 1, 4, 7]] / root_n[:, None], mode="reduced")
    directions = basis.T @ (nu / root_n[:, None])
    fluxes = scales * np.array([r.forward_flux(temp, densities) for r in kn.REACTIONS])
    matrix = (directions * fluxes) @ directions.T
    b = basis.T @ (kn.GROUND_ENERGIES * root_n)
    a = basis.T @ (internal_energies * root_n)
    capacity = ew.K_B * temp**2 * cv_j_m3_k
    strength = float(a @ a) / capacity
    root = math.sqrt(1 + strength)
    hroot = np.eye(4) + np.outer(a, a) / (capacity * (root + 1))
    hinvroot = np.eye(4) - np.outer(a, a) / (capacity * root * (root + 1))
    g = hinvroot @ b
    iso, si = _relaxation(matrix, directions, fluxes, b)
    thermal, st = _relaxation(hroot @ matrix @ hroot, hroot @ directions, fluxes, g)
    return Response(
        iso,
        thermal,
        tuple(float(s) for s in si),
        tuple(float(s) for s in st),
        strength,
        float(g @ g) / float(b @ b),
    )
