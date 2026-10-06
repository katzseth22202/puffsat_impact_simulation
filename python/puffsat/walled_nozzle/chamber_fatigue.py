"""Fatigue life of the near-term chambers' autofrettaged Cr-Mo shell (`sec:carbon_overwrap`).

The parent sizes a layered wall for the 20 m^3 chamber of the 2.5 kg rod. From the inside out it
has a GRCop-84 liner (hydrogen only), a 25 mm Cr-Mo shell, 5 mm of aluminium, an IM7/8552
overwrap and a 3 mm glass skin. Autofrettage leaves the steel in compression at rest. The design
rule lets each pulse swing the steel from -0.75 Y to +0.75 Y, with Y = 750 MPa and the dynamic
overshoot `D = 1.7` included. This module asks whether that steel survives its pulses. The parent
pairs the rod with hydrogen at 5500 K (818 bar) and methane at 7000 K (496 bar).

# The method is defect-tolerant, as hydrogen vessels are qualified

Hydrogen vessels are qualified by ASME Section VIII Div. 3 KD-10: assume a crack, grow it cycle
by cycle, and fail it at the fracture limit. The parent cites the same practice
(`sanmarchi2017_thickwall`). Here a semicircular surface crack on the shell's inner face starts
at 1 mm deep (0.5 and 2 mm are the sensitivity cases), with `K = 0.713 sigma sqrt(pi a)`. It fails
when `K_max` reaches 45 MPa m^0.5. That is the parent's hydrogen threshold for steels under 950
MPa, applied to both gases as a floor. It also fails when the crack reaches 80% of the 25 mm
wall.

Crack-growth curves:

- **Hydrogen:** the ASME Code Case 2938 master curve for quenched-and-tempered vessel steels
  (Ronevich et al., DOE Hydrogen Program AMR 2023, SCS005, slide 8). The rate is the smaller of
    - `3.5e-14 (1 + 0.4286 R)/(1 - R) dK^6.5 f(P)`, with `f(P) = 0.19 + 0.000763 P[bar]`, and
    - `1.5e-11 (1 + 2R)/(1 - R) dK^3.66`.

  The slide gives no units. m/cycle with MPa m^0.5 is the only reading under which hydrogen
  grows cracks faster than air does (5-10x, as observed), so that reading is used.
- **Methane:** the wall waits between pulses in methane, not hydrogen (parent
  `sec:methane_7000_near_term`), so the base case is the martensitic-steel air curve
  `1.36e-10 dK^2.25` (Barsom and Rolfe). The hydrogen curve at 496 bar is its bound.

Cycles with a compressive minimum are taken at `R = 0` with `dK = K_max`, the usual code
treatment.

# The load: the breathing mode rings

The pressure ramps over `t_r`, then blows down with the chamber's e-folding time
(`wall_heat.csv`). The shell's breathing mode (single mode, damping ratio zeta) is driven by it.
Steel stress is `-0.75 Y + S x(t)`, where `x` is the response normalised so a slow load gives
`x = p/p_peak`. `S = 1.5 Y / 1.7` is the parent's sizing. The overshoot leaves the wall ringing
about the slowly falling static stress, and every ring is a tensile cycle. Ramps run from instant
to one acoustic crossing (`2r/c`), and zeta from 0.005 to 0.05.

**Not modelled.**
- Residual pressure between pulses: at 4 Hz the chamber has not emptied, which narrows the swing.
- Thermal stress.
- Thick-wall gradients.
- The liner's own plasticity (GRCop-84 yields each pulse in the parent's design).
- Fragment or shock spikes.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

OUTPUT = Path("data/results/walled_nozzle/near_term/chamber_fatigue.csv")

YIELD = 750.0e6
"""Cr-Mo shell yield [Pa]: the parent's steel held under 950 MPa tensile strength."""
REST_STRESS = -0.75 * YIELD
STRESS_PER_STATIC = 1.5 * YIELD / 1.7
"""Steel stress [Pa] per unit static pressure ratio: the parent's sizing with `D = 1.7`."""
WALL = 0.025
K_LIMIT = 45.0
"""Failure `K_max` [MPa m^0.5]."""
CRACK_STARTS = (0.5e-3, 1.0e-3, 2.0e-3)
PULSE_PERIOD = 0.25
"""4 Hz [s] (parent: 630-780 pulses in 159-194 s)."""
PULSES_PER_DEPARTURE = 780
PULSES_FIVE_FLIGHTS = 5 * PULSES_PER_DEPARTURE
RADIUS = 1.684
"""Radius [m] of the 20 m^3 sphere."""


@dataclass(frozen=True)
class Layer:
    name: str
    modulus: float
    density: float
    thickness: float


#: Moduli and densities: steel 200 GPa/7830 (MIL-HDBK-5J, as the parent); aluminium 70/2700;
#: IM7/8552 quasi-isotropic 57.8 GPa/1579 (parent); glass-epoxy ~25 GPa (assumed)/1850; GRCop-84
#: 8620 kg/m^3 (parent) with an assumed copper modulus of 115 GPa. Carbon thickness is 18 cm for
#: methane (parent) and ~33 cm for hydrogen (the parent's 36 cm at 900 bar scaled to 818 bar).
STEEL = Layer("Cr-Mo", 200e9, 7830.0, 0.025)
COMMON = (Layer("aluminium", 70e9, 2700.0, 0.005), Layer("glass", 25e9, 1850.0, 0.003))


@dataclass(frozen=True)
class Chamber:
    name: str
    peak_bar: float
    sound_speed: float
    efold_s: tuple[float, ...]
    layers: tuple[Layer, ...]

    def breathing_period(self) -> float:
        """`T_b = 2 pi r / sqrt(2 sum(E t) / ((1 - nu) sum(rho t)))`, nu = 0.3."""
        stiff = sum(layer.modulus * layer.thickness for layer in self.layers)
        mass = sum(layer.density * layer.thickness for layer in self.layers)
        speed = math.sqrt(2.0 * stiff / (0.7 * mass))
        return 2.0 * math.pi * RADIUS / speed

    def crossing(self) -> float:
        """One acoustic crossing `2r/c` [s]: the slowest plausible pressure rise."""
        return 2.0 * RADIUS / self.sound_speed


HYDROGEN = Chamber(
    "hydrogen 5500 K",
    818.1,
    5471.0,
    (0.03378, 0.04488),
    (Layer("GRCop-84", 115e9, 8620.0, 0.005), STEEL, *COMMON, Layer("IM7", 57.8e9, 1579.0, 0.33)),
)
METHANE = Chamber(
    "methane 7000 K",
    495.9,
    4042.0,
    (0.04626, 0.06147),
    (STEEL, *COMMON, Layer("IM7", 57.8e9, 1579.0, 0.18)),
)


# ---- Crack growth ---------------------------------------------------------------------------


def rate_hydrogen(dk: float, r: float, pressure_bar: float) -> float:
    """ASME CC2938 master curve [m/cycle], the smaller of its two branches."""
    f_p = 0.19 + 0.000763 * pressure_bar
    low = 3.5e-14 * (1 + 0.4286 * r) / (1 - r) * dk**6.5 * f_p
    high = 1.5e-11 * (1 + 2 * r) / (1 - r) * dk**3.66
    return float(min(low, high))


def rate_air(dk: float, r: float, pressure_bar: float) -> float:
    """Barsom and Rolfe martensitic steels [m/cycle], no R correction."""
    del r, pressure_bar
    return float(1.36e-10 * dk**2.25)


RATES: dict[str, Callable[[float, float, float], float]] = {
    "hydrogen": rate_hydrogen,
    "air": rate_air,
}


def surface_k(stress: float, depth: float) -> float:
    """Semicircular surface crack, deepest point [MPa m^0.5]."""
    return 0.713 * stress / 1e6 * math.sqrt(math.pi * depth)


# ---- Load and response ----------------------------------------------------------------------


def pressure_ratio(t: NDArray[np.float64], ramp: float, efold: float) -> NDArray[np.float64]:
    """`p/p_peak`: linear ramp over `ramp`, then exponential blowdown."""
    rise = np.clip(t / ramp, 0.0, 1.0) if ramp > 0.0 else np.ones_like(t)
    fall = np.exp(-np.clip(t - ramp, 0.0, None) / efold)
    return np.asarray(np.where(t < ramp, rise, fall), dtype=np.float64)


def breathing_response(
    period: float, zeta: float, ramp: float, efold: float, duration: float = PULSE_PERIOD
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """`x(t)` under the ramp-and-blowdown load, for the period's own time grid."""
    t = np.arange(0.0, duration, period / 400.0, dtype=np.float64)
    return t, breathing_response_to(period, zeta, t, pressure_ratio(t, ramp, efold))


def breathing_response_to(
    period: float, zeta: float, t: NDArray[np.float64], f: NDArray[np.float64]
) -> NDArray[np.float64]:
    """`x(t)` for `x'' + 2 zeta w x' + w^2 x = w^2 f(t)` on a uniform grid `t`, by Newmark
    average acceleration (unconditionally stable, second order)."""
    omega = 2.0 * math.pi / period
    dt = float(t[1] - t[0])
    x = np.zeros_like(f)
    v = 0.0
    a = omega**2 * (f[0] - x[0])
    c = 2.0 * zeta * omega
    k_eff = omega**2 + 2.0 * c / dt + 4.0 / dt**2
    for n in range(1, t.size):
        rhs = omega**2 * f[n] + (4.0 / dt**2 + 2.0 * c / dt) * x[n - 1] + (4.0 / dt + c) * v + a
        x[n] = rhs / k_eff
        v_new = 2.0 / dt * (x[n] - x[n - 1]) - v
        a = 4.0 / dt**2 * (x[n] - x[n - 1]) - 4.0 / dt * v - a
        v = v_new
    return x


# ---- Rainflow (ASTM E1049) ------------------------------------------------------------------


def rainflow(series: NDArray[np.float64]) -> list[tuple[float, float, float]]:
    """`(minimum, maximum, count)` cycles."""
    d = np.diff(series)
    turning = np.flatnonzero(np.sign(d[1:]) * np.sign(d[:-1]) < 0) + 1
    points = [float(v) for v in series[np.concatenate(([0], turning, [len(series) - 1]))]]
    stack: list[float] = []
    cycles: list[tuple[float, float, float]] = []
    for p in points:
        stack.append(p)
        while len(stack) >= 3:
            x = abs(stack[-1] - stack[-2])
            y = abs(stack[-2] - stack[-3])
            if x < y:
                break
            lo, hi = sorted((stack[-2], stack[-3]))
            if len(stack) == 3:
                cycles.append((lo, hi, 0.5))
                stack.pop(0)
            else:
                cycles.append((lo, hi, 1.0))
                last = stack.pop()
                del stack[-2:]
                stack.append(last)
    for a, b in pairwise(stack):
        cycles.append((min(a, b), max(a, b), 0.5))
    return cycles


def pulses_to_failure(
    cycles: list[tuple[float, float, float]],
    a0: float,
    rate: Callable[[float, float, float], float],
    pressure_bar: float,
    cap: int = 1_000_000,
) -> int:
    """Pulses until `K_max` reaches the limit or the crack reaches 80% of the wall."""
    tensile = [(lo, hi, c) for lo, hi, c in cycles if hi > 0.0]
    if not tensile:
        return cap
    peak = max(hi for _, hi, _ in tensile)
    a = a0
    for pulse in range(cap):
        if surface_k(peak, a) >= K_LIMIT or a >= 0.8 * WALL:
            return pulse
        for lo, hi, count in tensile:
            r = max(lo, 0.0) / hi
            dk = surface_k(hi, a) * (1.0 - r)
            a += count * rate(dk, r, pressure_bar)
    return cap


@dataclass(frozen=True)
class Row:
    chamber: str
    growth: str
    breathing_ms: float
    ramp_ms: float
    efold_ms: float
    zeta: float
    overshoot: float
    peak_mpa: float
    rings_over_100mpa: int
    pulses_a05: int
    pulses_a1: int
    pulses_a2: int


def evaluate(chamber: Chamber, growth: str, ramp: float, efold: float, zeta: float) -> Row:
    """One chamber, crack-growth curve, ramp, blowdown and damping."""
    period = chamber.breathing_period()
    _, x = breathing_response(period, zeta, ramp, efold)
    stress = REST_STRESS + STRESS_PER_STATIC * x
    cycles = rainflow(stress)
    rate = RATES[growth]
    lives = [pulses_to_failure(cycles, a0, rate, chamber.peak_bar) for a0 in CRACK_STARTS]
    return Row(
        chamber=chamber.name,
        growth=growth,
        breathing_ms=period * 1e3,
        ramp_ms=ramp * 1e3,
        efold_ms=efold * 1e3,
        zeta=zeta,
        overshoot=float(x.max()),
        peak_mpa=float(stress.max()) / 1e6,
        rings_over_100mpa=sum(1 for lo, hi, c in cycles if hi - lo > 100e6 and c == 1.0),
        pulses_a05=lives[0],
        pulses_a1=lives[1],
        pulses_a2=lives[2],
    )


def cases() -> list[tuple[Chamber, str, float, float, float]]:
    out = []
    for chamber, growths in ((HYDROGEN, ("hydrogen",)), (METHANE, ("air", "hydrogen"))):
        for growth in growths:
            for ramp in (0.0, 0.4e-3, chamber.crossing()):
                for efold in chamber.efold_s:
                    for zeta in (0.005, 0.02, 0.05):
                        out.append((chamber, growth, ramp, efold, zeta))
    return out


def main() -> None:
    rows = [evaluate(*c) for c in cases()]
    print(
        f"steel at rest {REST_STRESS / 1e6:.0f} MPa, {STRESS_PER_STATIC / 1e6:.0f} MPa per static "
        f"unit; need {PULSES_PER_DEPARTURE} pulses (one departure), {PULSES_FIVE_FLIGHTS} (five)\n"
    )
    print(
        f"{'chamber':16s} {'curve':8s} {'T_b':>5} {'ramp':>5} {'efold':>5} {'zeta':>5} "
        f"{'D':>5} {'peak':>5} {'rings':>5} {'N(0.5mm)':>9} {'N(1mm)':>9} {'N(2mm)':>9}"
    )
    for r in rows:
        print(
            f"{r.chamber:16s} {r.growth:8s} {r.breathing_ms:5.2f} {r.ramp_ms:5.2f} "
            f"{r.efold_ms:5.1f} {r.zeta:5.3f} {r.overshoot:5.2f} {r.peak_mpa:5.0f} "
            f"{r.rings_over_100mpa:5d} {r.pulses_a05:9d} {r.pulses_a1:9d} {r.pulses_a2:9d}"
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Row.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([getattr(r, f) for f in fields])
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
