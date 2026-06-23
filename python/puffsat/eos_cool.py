"""Cool-gas two-phase equilibrium water EOS for the low-v (3.2 km/s) anchor (Rung C, C1).

At 3.2 km/s the stagnated gas is cool (~1700 K), neutral, and optically thin, so the bounce is
gas-dynamic re-expansion minus condensation (design §3, ADR-0004). Unlike the high-v package's
from-scratch plasma EOS ([`puffsat.eos_water`]), the low-v range is real-fluid water near and across
the saturation dome — exactly CoolProp's IAPWS-95 home. This module is the typed boundary around
CoolProp (untyped; isolated per-module in mypy, CLAUDE.md): it returns `(p, e, c_s, liquid_frac,
k_gas)` on a `(rho, T)` grid for the ADR-0007 table.

The **bulk vapor-pressure-collapse** channel (ADR-0004) is captured automatically — inside the dome
the equilibrium EOS returns the lower saturation pressure and folds latent heat into `e(rho, T)`.
The **wall-deposition** channel needs the condensed mass fraction, returned here as `liquid_frac`
for the kernel's wall-sticking sink (C3).

Two CoolProp specifics handled here: its speed of sound is undefined in two-phase, so `c_s` is
`sqrt((dp/drho)_s)` by a `(D, S)` finite difference (matches the analytic value in single phase, and
gives the low Wood sound speed in the dome); and the condensed fraction is bracketed by the
saturation densities rather than read from the `Q` sentinel, which is only valid inside the dome.
"""

from __future__ import annotations

import CoolProp.CoolProp as CP
import numpy as np
from numpy.typing import NDArray

Vec = NDArray[np.float64]

FLUID = "Water"
_E_FLOOR = 1.0  # J/kg — positivity safety for the loader's log(e); CoolProp U is already > 0 here.
_DS_FRAC = 1e-4  # relative density step for the (D, S) sound-speed finite difference.
_T_CRIT = float(CP.PropsSI("Tcrit", FLUID))  # 647.096 K — above this water is single-phase.


def _props(prop: str, rho: float, temp: float) -> float:
    """One CoolProp `(rho, T)` query as a float (the untyped-library boundary)."""
    return float(CP.PropsSI(prop, "D", rho, "T", temp, FLUID))


def pressure_energy(rho: float, temp: float) -> tuple[float, float]:
    """Equilibrium `(p [Pa], e [J/kg])` at `(rho, temp)` from CoolProp/IAPWS-95.

    `e` is the specific internal energy `U` (referenced to IAPWS triple-point liquid); it is already
    positive across the low-v box, floored at `_E_FLOOR` for the loader's log interpolation.
    """
    p = _props("P", rho, temp)
    e = max(_props("U", rho, temp), _E_FLOOR)
    return p, e


def sound_speed(rho: float, temp: float) -> float:
    """Equilibrium sound speed `c_s = sqrt((dp/drho)_s)` [m/s], by a `(D, S)` central difference.

    CoolProp's analytic speed of sound is undefined in the two-phase region, so we differentiate `p`
    at constant entropy instead. This matches CoolProp's `A` in single phase and yields the physical
    (low) Wood sound speed inside the dome.
    """
    s = _props("S", rho, temp)
    d = _DS_FRAC * rho
    p_hi = float(CP.PropsSI("P", "D", rho + d, "S", s, FLUID))
    p_lo = float(CP.PropsSI("P", "D", rho - d, "S", s, FLUID))
    cs2 = (p_hi - p_lo) / (2.0 * d)
    return float(np.sqrt(max(cs2, 0.0)))


def liquid_fraction(rho: float, temp: float) -> float:
    """Condensed mass fraction in `[0, 1]`: 0 in vapor/supercritical, 1 in compressed liquid, and
    the lever-rule `1 - Q` inside the two-phase dome.

    The phase is bracketed by the saturation densities `rho_g <= rho <= rho_f` at `temp` (CoolProp's
    `Q` is only meaningful there), so the result is robust without relying on a single-phase `Q`
    sentinel.
    """
    if temp >= _T_CRIT:
        return 0.0  # supercritical: single phase, no distinct condensate
    rho_g = float(CP.PropsSI("D", "T", temp, "Q", 1, FLUID))  # saturated vapor
    rho_f = float(CP.PropsSI("D", "T", temp, "Q", 0, FLUID))  # saturated liquid
    if rho <= rho_g:
        return 0.0  # superheated / unsaturated vapor
    if rho >= rho_f:
        return 1.0  # compressed liquid
    q = _props("Q", rho, temp)  # vapor quality, well-defined between the saturation densities
    return 1.0 - q


def conductivity(rho: float, temp: float) -> float:
    """Gas thermal conductivity `k_gas [W/m/K]` from CoolProp/IAPWS transport (the B-flux conduction
    operator's gas-side property; ADR-0005).

    CoolProp's transport conductivity is undefined for a two-phase `(rho, T)` input, so inside the
    saturation dome we evaluate it at the **saturated-vapor** density `rho_g(temp)` — the near-wall
    gas that the cold plate cools to drive condensation is vapor, not the mixture. Single-phase
    states (vapor, compressed liquid, supercritical) are queried directly. Strictly positive, so it
    rides the loader's log-interpolation path like the opacities.
    """
    if temp < _T_CRIT:
        rho_g = float(CP.PropsSI("D", "T", temp, "Q", 1, FLUID))  # saturated vapor
        rho_f = float(CP.PropsSI("D", "T", temp, "Q", 0, FLUID))  # saturated liquid
        if rho_g < rho < rho_f:  # two-phase: transport undefined → use the saturated vapor
            return _props("conductivity", rho_g, temp)
    return _props("conductivity", rho, temp)


def eos_grid_lowv(rho_grid: Vec, t_grid: Vec) -> tuple[Vec, Vec, Vec, Vec, Vec]:
    """Evaluate `(p, e, c_s, liquid_frac, k_gas)` on the `(rho, T)` grid, row-major over `(rho, T)`.

    Returns five `(n_rho, n_T)` arrays. Each `(rho, T)` node is independent (CoolProp is stateless),
    matching `eos_water.eos_grid` so the table assembly (`tables.py`) treats the two EOS the same.
    `k_gas` is the gas thermal conductivity for the B-flux conduction operator (ADR-0005).
    """
    n_rho, n_t = len(rho_grid), len(t_grid)
    p = np.empty((n_rho, n_t))
    e = np.empty((n_rho, n_t))
    cs = np.empty((n_rho, n_t))
    lf = np.empty((n_rho, n_t))
    kg = np.empty((n_rho, n_t))
    for i, rho in enumerate(rho_grid):
        for j, temp in enumerate(t_grid):
            r, t = float(rho), float(temp)
            p[i, j], e[i, j] = pressure_energy(r, t)
            cs[i, j] = sound_speed(r, t)
            lf[i, j] = liquid_fraction(r, t)
            kg[i, j] = conductivity(r, t)
    return p, e, cs, lf, kg
