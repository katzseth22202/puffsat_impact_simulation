"""Documented water PLOG rate: conditional extrapolation, not dense-water validation.

The Cantera example implements Singal et al. (2024), with Jasper (2022)
collider efficiencies. Its finite pressure table is NOT a high-pressure limit.
Temperature applicability at the study states is unestablished. Never silently
clamp or extend the pressure table. No mechanism thermodynamics are imported.
"""

from __future__ import annotations

import math

import numpy as np

from puffsat import eos_water as ew
from puffsat.water_plate.kinetics import Vec

SOURCE_URL = (
    "https://raw.githubusercontent.com/Cantera/cantera-example-data/"
    "1a5d27e508a38b1791543e9fded80ffd5c5b8d75/ammonia-CO-H2-Alzueta-2023.yaml"
)
SOURCE_SHA256 = "9eed1a29952ec41d154185798ecef7127a6329501322ee474be0b638011733f7"
ATM_PA = 101325.0
# H + OH (+M) <=> H2O (+M), reference collider N2, m³/kmol/s, cal/mol.
PLOG_NODES = (
    (1e-4, 5.30514e12, -2.80725, 499.267),
    (1e-3, 5.25581e13, -2.80630, 499.946),
    (1e-2, 5.18795e14, -2.80495, 501.765),
    (1e-1, 5.13043e15, -2.80388, 508.801),
    (1, 5.47458e16, -2.81214, 550.629),
    (10, 1.04665e18, -2.89077, 827.164),
    (100, 6.24786e18, -2.80241, 1433.20),
    (1000, 4.28006e15, -1.57172, 980.056),
    (10000, 6.79586e12, -0.577830, 456.911),
)
WATER_EFFICIENCY = (0.104529, 0.550787, -232.675)


def log_arrhenius(temp: float, a: float, b: float, ea: float) -> float:
    if not math.isfinite(temp) or temp <= 0:
        raise ValueError("temperature must be finite and positive")
    return math.log(a) + b * math.log(temp) - ea * 4.184 / (ew.N_A * ew.K_B * temp)


def water_efficiency(temp: float) -> float:
    """Relative to N2; extrapolated beyond the supplied 300/1000/2000 K data."""
    return math.exp(log_arrhenius(temp, *WATER_EFFICIENCY))


def effective_pressure(temp: float, densities: Vec) -> float:
    """Neutral ideal collider pressure times mixture efficiency, in Pa.

    Only H2O has a specific efficiency among the represented six neutrals;
    H/O/H2/OH/O2 inherit the source's default N2 value of one. This assumption
    is not evidence for hot radical collision efficiencies. Ions are absent.
    Because every collider uses the SAME reference PLOG curve, the linear-
    Burke mixture sum reduces exactly to k_ref(P_neutral * epsilon_mix).
    """
    if densities.shape != (6,) or not np.all(np.isfinite(densities)) or np.any(densities < 0):
        raise ValueError("need six finite nonnegative neutral densities")
    return (
        ew.K_B
        * temp
        * (float(np.sum(densities[:-1])) + water_efficiency(temp) * float(densities[-1]))
    )


def water_coefficient(temp: float, effective_pressure_pa: float) -> float:
    """Log(k)/log(P) interpolation, returning m³/molecule/s.

    No extra third-body density multiplies this bimolecular association rate.
    The upper pressure node is not k_infinity; reject outside the supplied grid.
    """
    p = effective_pressure_pa / ATM_PA
    if not math.isfinite(p) or not PLOG_NODES[0][0] <= p <= PLOG_NODES[-1][0]:
        raise ValueError("effective pressure outside source PLOG table")
    logs = [log_arrhenius(temp, a, b, ea) for _, a, b, ea in PLOG_NODES]
    value = float(np.interp(math.log(p), np.log([r[0] for r in PLOG_NODES]), logs))
    return math.exp(value) / (1000 * ew.N_A)
