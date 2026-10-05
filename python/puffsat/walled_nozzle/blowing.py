"""Carbon blowing parameter B'(T, p) of a pitch/graphite wall in hot hydrogen (N18 item 1).

**This module computes no thermochemistry of its own.** The table is the companion repository's
full-species solve, copied verbatim:

- source: `katzseth22202/aim_is_all_you_need`, `docs/companion_replies_2026-10-04.md`, section H1,
  @ `acfe859` (reproduce there with `make carbon-equilibrium`);
- method: Cantera 3.2.0, NASA Glenn data, every C/H species in `nasa_gas.yaml` plus `C(gr)`, ideal
  gas, pure hydrogen edge gas at equilibrium with unit-activity graphite;
- unit: kg of carbon held in the gas per kg of hydrogen.

Use this rather than `near_term.carbon_blowing_parameter`: that function lacks C4H2, C6H2 and
C2H4, so it reads `ln(1 + B')` about 15% low at 3900 K (the cross-check in H1).

**Edges of the table.** It spans 2500-4100 K and 10-818 bar. Beyond them the lookup clamps in
pressure and in the hot direction, and extrapolates **below 2500 K** along `ln B'` linear in `1/T`
from the two coolest columns (an Arrhenius tail). That tail is an assumption, not data; it only
matters for the first milliseconds, while the surface is still cold and B' is small anyway.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

Vec = NDArray[np.float64]

#: Table temperatures [K] (rows) and pressures [bar] (columns).
TEMPS = np.array([2500.0, 2800.0, 3100.0, 3400.0, 3700.0, 3900.0, 4100.0])
PRESSURES_BAR = np.array([10.0, 30.0, 62.0, 100.0, 300.0, 818.0])

#: B' [kg C / kg H] at `(TEMPS[i], PRESSURES_BAR[j])`.
B_PRIME = np.array(
    [
        [0.14, 0.16, 0.19, 0.22, 0.38, 0.70],
        [0.40, 0.42, 0.44, 0.46, 0.57, 0.81],
        [0.96, 0.98, 1.00, 1.02, 1.10, 1.27],
        [2.00, 2.02, 2.03, 2.05, 2.11, 2.23],
        [3.90, 3.74, 3.72, 3.72, 3.75, 3.81],
        [6.50, 5.55, 5.36, 5.31, 5.27, 5.28],
        [13.51, 8.61, 7.73, 7.45, 7.19, 7.10],
    ]
)

_LN_B = np.log(B_PRIME)
_LN_P = np.log(PRESSURES_BAR)


def _ln_b_at_pressure(row: int, ln_p: float) -> float:
    """`ln B'` on one temperature row, linear in `ln p`, clamped to the table's pressures."""
    return float(np.interp(ln_p, _LN_P, _LN_B[row]))


def blowing_parameter(temp: float, pressure: float) -> float:
    """B' at wall temperature `temp` [K] and pressure `pressure` [Pa].

    Bilinear in `(1/T, ln p)` on `ln B'`: B' is Arrhenius-like in temperature and nearly a power
    of pressure, so this keeps it positive and smooth across two decades.
    """
    ln_p = math.log(max(pressure, 1.0) / 1.0e5)
    col = np.array([_ln_b_at_pressure(i, ln_p) for i in range(len(TEMPS))])
    inv_t = 1.0 / np.asarray(TEMPS)
    x = 1.0 / temp
    if temp >= TEMPS[-1]:
        return float(math.exp(col[-1]))
    if temp < TEMPS[0]:
        slope = (col[1] - col[0]) / (inv_t[1] - inv_t[0])
        return float(math.exp(col[0] + slope * (x - inv_t[0])))
    # `np.interp` wants an increasing abscissa; 1/T decreases with T, so reverse both.
    return float(math.exp(np.interp(x, inv_t[::-1], col[::-1])))


def mass_flux_factor(b_prime: float) -> float:
    """`ln(1 + B')`: the blown mass flux in units of the transfer coefficient `g = h_g / c_p`."""
    return math.log1p(b_prime)


def blowing_correction(b_prime: float) -> float:
    """`ln(1 + B') / B'`: the factor by which the blown film cuts the convective heat flux."""
    return math.log1p(b_prime) / b_prime
