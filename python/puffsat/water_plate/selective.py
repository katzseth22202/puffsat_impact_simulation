"""Molecular abundances fixed, atomic ionization equilibrated (no reaction rates).

This constrained ideal-mixture equilibrium retains each parcel's molecular
inventory and bond store. It is the intermediate closure for an ordered,
conditional decomposition, not a unique molecular/ionization attribution.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from puffsat import eos_water as ew
from puffsat.water_plate.freeze import FrozenGas, gas_parameters


@dataclass(frozen=True)
class IonLadder:
    abundance: float
    energies: list[float]
    log_prefactors: list[float]


@dataclass(frozen=True)
class IonizingGas:
    base: FrozenGas
    ladders: list[IonLadder]
    boltzmann: float = ew.K_B
    formula_mass: float = ew.M_H2O

    def populations(self, rho: float, temp: float) -> list[list[float]]:
        """Charge-neutral Saha populations per formula unit, including neutrals."""
        if not (rho > 0 and temp >= 1000 and math.isfinite(rho) and math.isfinite(temp)):
            raise ValueError("selective gas requires finite rho > 0 and T >= 1000 K")
        ln_nf = math.log(rho / self.formula_mass)
        maximum = sum(a.abundance * (len(a.energies) - 1) for a in self.ladders)
        if maximum == 0:
            return [[0.0] * len(a.energies) for a in self.ladders]

        def fractions(ln_ne: float) -> list[list[float]]:
            result = []
            for a in self.ladders:
                logs = [
                    c + 1.5 * z * math.log(temp) - e / (self.boltzmann * temp) - z * ln_ne
                    for z, (c, e) in enumerate(zip(a.log_prefactors, a.energies, strict=True))
                ]
                peak = max(logs)
                weights = [math.exp(x - peak) for x in logs]
                total = sum(weights)
                result.append([a.abundance * x / total for x in weights])
            return result

        def residual(ln_ne: float) -> float:
            q = sum(z * p for row in fractions(ln_ne) for z, p in enumerate(row))
            return ln_ne - ln_nf - math.log(max(q, 1e-300))

        hi = ln_nf + math.log(maximum)
        lo = hi - 700
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if residual(mid) > 0:
                hi = mid
            else:
                lo = mid
        return fractions(0.5 * (lo + hi))

    def pressure_energy(self, rho: float, temp: float) -> tuple[float, float]:
        populations = self.populations(rho, temp)
        electrons = sum(z * p for row in populations for z, p in enumerate(row))
        ion = (
            sum(
                p * e
                for a, row in zip(self.ladders, populations, strict=True)
                for p, e in zip(row, a.energies, strict=True)
            )
            / self.formula_mass
        )
        electron_r = self.boltzmann / self.formula_mass * electrons
        return (
            rho * (self.base.gas_constant + electron_r) * temp,
            self.base.thermal(temp) + 1.5 * electron_r * temp + ion,
        )


def ionizing_parameters(rho: float, temp: float) -> IonizingGas:
    y = ew.frozen_composition(rho, temp)
    original = gas_parameters(rho, temp)
    electron_r = ew.K_B / ew.M_H2O * y.y_e
    base = FrozenGas(
        original.gas_constant - electron_r,
        original.linear_cv - 1.5 * electron_r,
        original.bond_store,
        original.vibrations,
        original.bond_store,
        0.0,
    )

    def ladder(abundance: float, energies: list[float], degeneracies: list[float]) -> IonLadder:
        # ln K_z(T) = c_z + 1.5 ln T - IP_z/(kT); accumulate stage products.
        constants = [0.0]
        for i in range(1, len(energies)):
            ip = energies[i] - energies[i - 1]
            c = ew.ln_k_saha(ip, degeneracies[i], degeneracies[i - 1], 10000.0)
            constants.append(constants[-1] + c - 1.5 * math.log(10000) + ip / (ew.K_B * 10000))
        return IonLadder(abundance, energies, constants)

    return IonizingGas(
        base,
        [
            ladder(y.y_h + y.y_hp, [0.0, ew.IP_H], [ew.G_H, ew.G_HP]),
            ladder(y.y_o + sum(y.y_o_ions), [0.0, *ew.E_O_CUM], list(ew.G_O_LADDER)),
        ],
    )
