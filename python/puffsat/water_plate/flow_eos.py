"""Screening EOS joining real-fluid water to the equilibrium gas on one energy zero.

Above 1000 K continue the real-fluid residual with heat capacity proportional
to T^-2, matching pressure, energy, entropy and heat capacity at the join.
All fields derive from this one Helmholtz potential. This is a provisional
dense/hot closure, not validated warm-dense water. No loss/opacity claim follows
from the transparent placeholder fields required by the loader.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import CoolProp.CoolProp as CP
import numpy as np

from puffsat import eos_water
from puffsat.water_plate.thermodynamics import energy_reference

ENERGY_SHIFT = 4e6  # Add to EVERY phase; positive storage for the existing log-energy loader.
COLD_MAX_T = 1000.0
TABLE_PATH = Path("data/tables/water_plate_flow.json")
RESIDUAL_POWER = 2.0


def gas_free_energy(rho: float, temp: float) -> float:
    """Helmholtz energy from the SAME partition functions/species as eos_water.

    Sum n*[epsilon0 + kT*(ln(n/q)-1)] / rho. Choose H's elemental energy zero
    as zero and O's as D_AT, so H2O's bound state is zero. Elemental gauge
    choices cancel in every reaction and reproduce the existing energy sum.
    """
    ew = eos_water
    comp = ew.composition(rho, temp)
    species = [
        (comp.n_h2o, ew.M_H2O, ew.G_H2O * ew._z_rot(temp) * ew._z_vib(temp), 0.0),
        (comp.n_h, ew.M_H, ew.G_H, 0.0),
        (comp.n_o, ew.M_O, ew.G_O, ew.D_AT),
        (comp.n_hp, ew.M_H, ew.G_HP, ew.IP_H),
        (comp.n_e, ew.M_E, ew.G_E, 0.0),
    ]
    species.extend(
        (n, ew.M_O, ew.G_O_LADDER[i + 1], ew.D_AT + ew.E_O_CUM[i])
        for i, n in enumerate(comp.n_o_ions)
    )
    for n, d, elemental in zip(
        (comp.n_oh, comp.n_h2, comp.n_o2),
        ew.INTERMEDIATES,
        (ew.D_AT, 0.0, 2 * ew.D_AT),
        strict=True,
    ):
        z = d.degeneracy * ew._z_rot_linear(d.theta_rot, d.symmetry, temp)
        z *= ew._z_vib_mode(d.theta_vib, temp)
        species.append((n, d.mass, z, elemental - d.d0))
    return (
        sum(
            n * (e0 + ew.K_B * temp * (math.log(n) - ew._ln_trans(mass, temp) - math.log(z) - 1))
            for n, mass, z, e0 in species
            if n > 0
        )
        / rho
    )


@cache
def entropy_reference_shift() -> float:
    """Align the additive entropy constant at the existing 400 K, 0.01 kg/m3 bridge."""
    rho, temp = 0.01, 400.0
    e = float(CP.PropsSI("U", "D", rho, "T", temp, "Water")) + energy_reference().offset_j_kg
    s = float(CP.PropsSI("S", "D", rho, "T", temp, "Water"))
    return (e - temp * s - gas_free_energy(rho, temp)) / temp


def _cold_free_energy(rho: float, temp: float) -> float:
    energy = float(CP.PropsSI("U", "D", rho, "T", temp, "Water"))
    entropy = float(CP.PropsSI("S", "D", rho, "T", temp, "Water"))
    return energy + energy_reference().offset_j_kg - temp * (entropy + entropy_reference_shift())


@dataclass(frozen=True)
class ResidualReference:
    """Real-fluid minus ideal equilibrium-gas derivatives at the matching state."""

    energy: float
    entropy: float
    cv: float
    pressure: float
    pressure_t: float
    pressure_tt: float


@cache
def residual_reference(rho: float) -> ResidualReference:
    temp = COLD_MAX_T
    p, e = eos_water.pressure_energy(rho, temp)
    de = float(CP.PropsSI("U", "D", rho, "T", temp, "Water")) + energy_reference().offset_j_kg - e
    df = _cold_free_energy(rho, temp) - gas_free_energy(rho, temp)
    dt = 1.0
    pm, pp = (eos_water.pressure_energy(rho, temp + d)[0] for d in (-dt, dt))
    pt = float(CP.PropsSI("d(P)/d(T)|D", "D", rho, "T", temp, "Water")) - (pp - pm) / (2 * dt)
    ptt = float(CP.PropsSI("d(d(P)/d(T)|D)/d(T)|D", "D", rho, "T", temp, "Water"))
    ptt -= (pp - 2 * p + pm) / dt**2
    em, ep = (eos_water.pressure_energy(rho, temp + d)[1] for d in (-0.1, 0.1))
    cv = float(CP.PropsSI("CVMASS", "D", rho, "T", temp, "Water")) - (ep - em) / 0.2
    return ResidualReference(
        de,
        (de - df) / temp,
        cv,
        float(CP.PropsSI("P", "D", rho, "T", temp, "Water")) - p,
        pt,
        ptt,
    )


def dense_residual(rho: float, temp: float) -> tuple[float, float, float]:
    """Return residual (p,e,f), integrating cv_res=cv_ref*(T_ref/T)^n.

    Integrating de=cv*dT and ds=cv*dT/T fixes the residual potential uniquely
    from the matching state. Its density derivative uses the Maxwell identity
    rho²*d(cv_res)/d(rho)=-T_ref*d²(p_res)/dT². Unlike switching potentials
    with a temperature weight, this adds no uncontrolled heat-capacity term.
    n is an explicit sensitivity parameter, not a fitted physical rate.
    """
    ref = residual_reference(rho)
    n, t0 = RESIDUAL_POWER, COLD_MAX_T
    q = t0 / temp
    energy_basis = t0 / (n - 1) * (1 - q ** (n - 1))
    entropy_basis = (1 - q**n) / n
    potential_basis = energy_basis - temp * entropy_basis
    energy = ref.energy + ref.cv * energy_basis
    potential = ref.energy - temp * ref.entropy + ref.cv * potential_basis
    pressure = ref.pressure + (temp - t0) * ref.pressure_t - t0 * ref.pressure_tt * potential_basis
    return pressure, energy, potential


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """Return p, e+shift; the same collision-paid liquid/gas zero as the mixing screen."""
    if temp <= COLD_MAX_T:
        pressure = float(CP.PropsSI("P", "D", rho, "T", temp, "Water"))
        energy = float(CP.PropsSI("U", "D", rho, "T", temp, "Water"))
        return pressure, energy + energy_reference().offset_j_kg + ENERGY_SHIFT
    hot_p, hot_e = eos_water.pressure_energy(rho, temp)
    residual_p, residual_e, _ = dense_residual(rho, temp)
    return hot_p + residual_p, hot_e + residual_e + ENERGY_SHIFT


def derivatives(rho: float, temp: float) -> tuple[float, float, float, float]:
    """p_rho, e_rho, p_T, e_T at fixed other coordinate, by centered differences."""
    dr, dt = 1e-4 * rho, 1e-4 * temp
    rlo, rhi = pressure_energy(rho - dr, temp), pressure_energy(rho + dr, temp)
    tlo, thi = pressure_energy(rho, temp - dt), pressure_energy(rho, temp + dt)
    return (
        (rhi[0] - rlo[0]) / (2 * dr),
        (rhi[1] - rlo[1]) / (2 * dr),
        (thi[0] - tlo[0]) / (2 * dt),
        (thi[1] - tlo[1]) / (2 * dt),
    )


def sound_speed(rho: float, temp: float) -> float:
    """c² = p_rho + p_T*(p/rho²-e_rho)/e_T along de=p/rho² d(rho).

    Fail on an unstable closure rather than flooring a negative acoustic speed.
    """
    p = pressure_energy(rho, temp)[0]
    return _checked_sound_speed(rho, temp, p, derivatives(rho, temp))


def _checked_sound_speed(
    rho: float, temp: float, p: float, partials: tuple[float, float, float, float]
) -> float:
    pr, er, pt, et = partials
    cs2 = pr + pt * (p / rho**2 - er) / et
    if et <= 0 or cs2 <= 0 or not math.isfinite(cs2):
        raise ValueError(f"unstable EOS at rho={rho:g}, T={temp:g}: cv={et:g}, c²={cs2:g}")
    return math.sqrt(cs2)


def liquid_source_temperature(rho: float) -> float:
    """Small equilibrium flash of stored liquid, at conserved source internal energy.

    The source bulk density is an input; pressure/T follow from its liquid energy.
    A density requiring cooling below the table's liquid range is rejected.
    """
    target = energy_reference().liquid_energy_j_kg + ENERGY_SHIFT
    lo, hi = 273.2, 400.0
    if pressure_energy(rho, lo)[1] > target or pressure_energy(rho, hi)[1] < target:
        raise ValueError("source energy requires a phase state outside the liquid/vapor screen")
    for _ in range(50):
        mid = (lo + hi) / 2
        if pressure_energy(rho, mid)[1] < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main() -> None:
    """Generate and audit the provisional table, without hiding unstable states."""
    global RESIDUAL_POWER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=TABLE_PATH)
    parser.add_argument("--residual-power", type=float, default=RESIDUAL_POWER)
    args = parser.parse_args()
    if not math.isfinite(args.residual_power) or args.residual_power <= 1:
        parser.error("residual power must be finite and greater than 1")
    RESIDUAL_POWER = args.residual_power
    rho_grid = np.geomspace(1e-5, 2000, 72)
    t_grid = np.unique(
        np.concatenate(
            (
                np.linspace(273.2, 400, 18),
                np.geomspace(400, 2000, 32),
                np.geomspace(2000, 2e6, 72),
                np.array([COLD_MAX_T]),
            )
        )
    )
    p: list[float] = []
    e: list[float] = []
    cs: list[float] = []
    worst = (0.0, 0.0, 0.0)
    min_cv = math.inf
    for rho in rho_grid:
        for temp in t_grid:
            r, t = float(rho), float(temp)
            pressure, energy = pressure_energy(r, t)
            if not all(math.isfinite(v) and v > 0 for v in (pressure, energy)):
                raise ValueError(f"nonpositive or nonfinite EOS field at rho={r}, T={t}")
            partials = derivatives(r, t)
            p.append(pressure)
            e.append(energy)
            cs.append(_checked_sound_speed(r, t, pressure, partials))
            min_cv = min(min_cv, partials[3])
            if t > COLD_MAX_T:
                _, er, pt, _ = partials
                residual = abs(r * r * er - pressure + t * pt) / pressure
                if residual > worst[0]:
                    worst = (residual, r, t)
        print(f"EOS density {float(rho):.5g} kg/m3 complete.", flush=True)
    if not np.all(np.diff(np.array(e).reshape(len(rho_grid), -1), axis=1) > 0):
        raise ValueError("nonmonotone tabulated energy prevents unique temperature inversion")
    if worst[0] > 1e-3:
        raise ValueError(f"Maxwell residual exceeds 0.1%: {worst}")
    transparent = [1e-10] * len(p)
    table = {
        "rho_grid": rho_grid.tolist(),
        "T_grid": t_grid.tolist(),
        "shape": [len(rho_grid), len(t_grid)],
        "fields": {
            "p": p,
            "e": e,
            "c_s": cs,
            "kappa_rosseland": transparent,
            "kappa_planck": transparent,
        },
        "provenance": {
            "status": (
                "IAPWS-95 with a cv-matched residual extension of ideal equilibrium water; "
                "dense/hot model sensitivity required"
            ),
            "coolprop_version": energy_reference().coolprop_version,
            "liquid_energy_on_gas_zero": energy_reference().liquid_energy_j_kg,
            "energy_storage_shift_j_kg": ENERGY_SHIFT,
            "real_fluid_max_k": COLD_MAX_T,
            "dense_hot_cv_residual_power": RESIDUAL_POWER,
            "dense_hot_model": "cv_res=cv_ref*(T_ref/T)^n; provisional thermodynamic continuation",
            "max_relative_maxwell_residual": worst,
            "minimum_cv_j_kg_k": min_cv,
            "radiation": "disabled; opacity fields are loader placeholders",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(table) + "\n")
    print(f"Wrote {args.output}; worst Maxwell residual (fraction,rho,T): {worst}")


if __name__ == "__main__":
    main()
