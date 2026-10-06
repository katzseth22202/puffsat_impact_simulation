"""The skirt's hoop response and mass (ADR-0055, Q56).

The skirt is a thin ring of radius `R`: a steel liner wrapped in fibre. Gas pressure on its inner
face drives the ring's breathing mode, a single oscillator per unit area,

    m'' x'' + k'' x = p(t),   m'' = sum(rho_i t_i),   k'' = sum(E_i t_i) / R^2,

with hoop strain `x / R`. The layers strain together (the wrap is wound tight on the liner). The
pressure history is taken as constant over each sample interval, and each interval is stepped
exactly, so a short pulse reproduces the impulsive strain `I / (m'' c)` and a held pressure twice
the static `p R / (E t)`. The peak includes the free ring's ringing after the last sample.

`main` sizes the wrap for each wall row of the recorded 2-D skirt histories
(`make water-plate-argon-2d-skirt`), scaled to the real pulse.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class Material:
    name: str
    modulus: float
    density: float


@dataclass(frozen=True)
class Layer:
    material: Material
    thickness: float


def peak_hoop_strain(
    time: NDArray[np.float64],
    pressure: NDArray[np.float64],
    radius: float,
    layers: list[Layer],
) -> float:
    """Peak hoop strain of the ring under `pressure[i]` held over `[time[i], time[i+1])`."""
    mass = sum(layer.material.density * layer.thickness for layer in layers)
    stiffness = sum(layer.material.modulus * layer.thickness for layer in layers) / radius**2
    omega = math.sqrt(stiffness / mass)
    x, v = 0.0, 0.0
    peak = 0.0
    for i in range(len(time) - 1):
        dt = float(time[i + 1] - time[i])
        x_eq = float(pressure[i]) / stiffness
        # Within the interval the motion is x_eq + amp cos(omega s + phase). Its top (cos = 1) and
        # bottom (cos = -1) count if they fall inside the interval.
        amp = math.hypot(x - x_eq, v / omega)
        phase = math.atan2(-v / omega, x - x_eq)
        for target, sign in ((0.0, 1.0), (math.pi, -1.0)):
            s_hit = ((target - phase) % (2.0 * math.pi)) / omega
            if 0.0 < s_hit <= dt:
                peak = max(peak, abs(x_eq + sign * amp))
        c, s = math.cos(omega * dt), math.sin(omega * dt)
        x, v = x_eq + (x - x_eq) * c + v / omega * s, -(x - x_eq) * omega * s + v * c
        peak = max(peak, abs(x))
    # After the last sample the ring rings freely about zero.
    peak = max(peak, math.hypot(x, v / omega))
    return peak / radius


# ---- Sizing the skirt from the recorded 2-D histories ------------------------------------------

HISTORIES = Path("data/results/water_plate/spray_2d_skirt.jsonl")
OUTPUT = Path("data/results/water_plate/skirt_hoop.csv")

#: The real pulse the 2-D floor impulse is scaled to [N s] (the 100 t plate's gas spring, P14),
#: and the closing speed that sets the time unit [m/s].
PULSE_IMPULSE = 8.0e6
CLOSING_SPEED = 45580.0

#: Maraging 300 liner: E ~190 GPa, rho 8000 (as `plate_thermal.MARAGING`). Aramid (Kevlar 49
#: class): E 112 GPa, rho 1440 (manufacturer's datasheet values, not verified here). Carbon
#: (T700 class) as the fallback wrap: E 230 GPa, rho 1800 (datasheet class values).
MARAGING = Material("maraging 300", modulus=190.0e9, density=8000.0)
ARAMID = Material("aramid", modulus=112.0e9, density=1440.0)
CARBON = Material("carbon", modulus=230.0e9, density=1800.0)

#: Hoop strain allowables: maraging yields near 1.05% (2.0 GPa / 190 GPa), so 0.7% keeps a 1.5
#: margin on the liner. 0.5% and 1.0% bracket it.
STRAIN_ALLOWABLES = (0.005, 0.007, 0.010)
LINER_THICKNESSES = (0.002, 0.004)


def wrap_thickness(
    time: NDArray[np.float64],
    pressure: NDArray[np.float64],
    radius: float,
    liner: Layer,
    wrap: Material,
    allowable: float,
) -> float:
    """The thinnest wrap over `liner` that keeps the peak hoop strain at `allowable` [m]."""
    if peak_hoop_strain(time, pressure, radius, [liner]) <= allowable:
        return 0.0
    lo, hi = 0.0, 0.01
    while peak_hoop_strain(time, pressure, radius, [liner, Layer(wrap, hi)]) > allowable:
        hi *= 2.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if peak_hoop_strain(time, pressure, radius, [liner, Layer(wrap, mid)]) > allowable:
            lo = mid
        else:
            hi = mid
    return hi


def main() -> None:
    rows_out: list[dict[str, object]] = []
    with HISTORIES.open() as fh:
        records = [json.loads(line) for line in fh]
    for rec in records:
        unit = float(rec["unit_m"])
        radius = float(rec["r_plate"]) * unit
        tau = unit / CLOSING_SPEED
        # The 2-D floor impulse, 2 pi restored, maps onto the real pulse.
        scale_impulse = PULSE_IMPULSE / (2.0 * math.pi * float(rec["wall_impulse"]))
        scale_pressure = scale_impulse / (unit * unit * tau)
        time = np.asarray(rec["time"], dtype=np.float64) * tau
        pressure = np.asarray(rec["pressure"], dtype=np.float64) * scale_pressure
        z = np.asarray(rec["z"], dtype=np.float64) * unit
        dz = float(z[1] - z[0])
        impulse = np.sum(pressure[:-1] * np.diff(time)[:, None], axis=0)
        peak = pressure.max(axis=0)
        worst = int(np.argmax(impulse))
        period = 2.0 * math.pi * radius / math.sqrt(ARAMID.modulus / ARAMID.density)
        for wrap in (ARAMID, CARBON):
            for t_liner in LINER_THICKNESSES:
                liner = Layer(MARAGING, t_liner)
                for allowable in STRAIN_ALLOWABLES:
                    wraps = [
                        wrap_thickness(time, pressure[:, j], radius, liner, wrap, allowable)
                        for j in range(len(z))
                    ]
                    band = 2.0 * math.pi * radius * dz
                    steel = band * len(z) * MARAGING.density * t_liner
                    fibre = band * wrap.density * float(np.sum(wraps))
                    rows_out.append(
                        {
                            "case": rec["case"],
                            "wrap": wrap.name,
                            "liner_mm": t_liner * 1e3,
                            "allowable": allowable,
                            "max_wrap_mm": max(wraps) * 1e3,
                            "steel_t": steel / 1e3,
                            "fibre_t": fibre / 1e3,
                            "skirt_t": (steel + fibre) / 1e3,
                            "worst_row_m": float(z[worst] - z[0]),
                            "worst_impulse_kPa_s": float(impulse[worst]) / 1e3,
                            "worst_peak_MPa": float(peak[worst]) / 1e6,
                            "worst_duration_ms": float(impulse[worst] / peak[worst]) * 1e3,
                            "ring_period_ms": period * 1e3,
                        }
                    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
        writer.writeheader()
        writer.writerows(rows_out)
    for r in rows_out:
        print(
            f"{r['case']:<22} {r['wrap']:<6} liner {r['liner_mm']:.0f} mm "
            f"eps {r['allowable']:.3f}: "
            f"wrap max {r['max_wrap_mm']:5.1f} mm, skirt {r['skirt_t']:5.1f} t "
            f"(steel {r['steel_t']:4.1f} + fibre {r['fibre_t']:4.1f}); worst row "
            f"{r['worst_impulse_kPa_s']:.2f} kPa s, peak {r['worst_peak_MPa']:.1f} MPa, "
            f"~{r['worst_duration_ms']:.2f} ms vs ring {r['ring_period_ms']:.1f} ms"
        )


if __name__ == "__main__":
    main()
