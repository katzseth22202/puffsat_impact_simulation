"""Spall and delamination of the layered chamber wall under a reflected blast spike.

`chamber_fatigue.py` treats the wall as a breathing shell under the chamber's mean pressure. This
module covers the other failure mode: a short, sharp pressure spike on the inner face. The
parent's `sec:carbon_overwrap` names the risk. A spike crosses the steel in ~4 µs. At each drop
in acoustic impedance part of it reflects as tension. When it reaches the free outer skin it
returns entirely as tension. Steel is strong in tension. The resin between carbon plies is not:
IM7/8552 holds 64 MPa across its fibres.

# The model

The wave is radial, plane and linear elastic, through the stack in the parent's order:
liner (hydrogen only), Cr-Mo, aluminium (optional), IM7/8552, glass skin, free outside. Each
layer carries right- and left-going stress waves on its own grid, `dx_i = c_i dt`, so advection
inside a layer is exact. At each interface the stress waves transmit with
`2 Z_b / (Z_a + Z_b)` and reflect with `(Z_b - Z_a)/(Z_a + Z_b)`. The inner face carries the
traction `sigma = -p(t)`.

A radial bias is added, the at-rest compression that autofrettage leaves at the steel/overwrap
interface: `2 sigma_rest t / r` ≈ 16.7 MPa. It falls linearly to zero across the steel and the
carbon. Because the spike arrives first, before the mean pressure builds, this bias is the only
compression helping.

Being linear, the model reports **tension per GPa of spike**, in each layer. From that it gives
the spike amplitude at which each material fails:

- steel spall: ~2.9 GPa. That is `0.5 rho c_b W` from 38KhN3MFA, a Q&T Cr-Ni-Mo-V steel
  (Mescheryakov et al. 2017, Table 3).
- carbon and glass laminates, and the bond lines: 64 MPa across the fibres (Hexcel 8552, as the
  parent).

# The spike

The amplitude is not solved here; `chamber_blast.py` and `vessel_blast.py` solve it (see
`docs/walled_nozzle_chamber_blast.md`). `sedov_spike` below is the first estimate those solves
replaced. It is a strong-shock Sedov arrival with normal reflection, using xi = 1.033 at
gamma 1.4. An earlier version used 1.15, the gamma = 5/3 value, and overstated the shock speed
by ~30%.

Impedances follow the parent (`sec:carbon_overwrap`):
- steel 5,900 m/s and 7,830 kg/m^3;
- aluminium 6,320 and 2,700;
- IM7/8552 through-thickness 4.8 MPa s/m at 1,579 kg/m^3;
- glass-epoxy 2,740 and 1,850.

The GRCop-84 speed of 4,700 m/s is assumed.

# Unbonded interfaces

A layer with `bonded=False` meets the next layer in contact only, as the shells of a
shrink-fitted multilayer wall or a dry wrap do. The contact passes compression exactly as a bond
does. It opens when the interface stress would exceed the shrink-fit `preload` in tension. While
open, both faces are free (`sigma = 0`), and the gap grows at the faces' velocity difference,
`u = (s - r)/Z`. It recloses when the gap shuts. A short pulse is then trapped in the outer layer,
the momentum trap of spall experiments, and no tension returns inward.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

OUTPUT = Path("data/results/walled_nozzle/near_term/wall_waves.csv")

RADIUS = 1.684
REST_HOOP = -562.5e6
"""Autofrettaged steel hoop stress at rest [Pa] (-0.75 Y, Y = 750 MPa)."""
STEEL_SPALL = 2.9e9
"""Q&T Cr-Ni-Mo-V steel, `0.5 rho c_b W` at W = 160 m/s, the lowest of [M17] Table 3."""
RESIN_TENSION = 64.0e6
"""IM7/8552 transverse tensile strength [Pa] (Hexcel datasheet, as the parent)."""


Load = Callable[[float], float]
"""Inner-face pressure [Pa] against time [s]."""


@dataclass(frozen=True)
class Layer:
    name: str
    thickness: float
    speed: float
    density: float
    strength: float
    bonded: bool = True
    """Whether this layer is bonded to the next one out (False: contact only)."""

    @property
    def impedance(self) -> float:
        return self.density * self.speed


def liner() -> Layer:
    return Layer("GRCop-84", 0.005, 4700.0, 8620.0, math.inf)


def steel() -> Layer:
    return Layer("Cr-Mo", 0.025, 5900.0, 7830.0, STEEL_SPALL)


def aluminium() -> Layer:
    return Layer("aluminium", 0.005, 6320.0, 2700.0, math.inf)


def carbon(thickness: float) -> Layer:
    return Layer("IM7/8552", thickness, 4.8e6 / 1579.0, 1579.0, RESIN_TENSION)


def glass() -> Layer:
    return Layer("glass skin", 0.003, 2740.0, 1850.0, RESIN_TENSION)


def interface_bias() -> float:
    """At-rest radial compression [Pa] at the steel/overwrap interface (thin-sphere balance)."""
    return 2.0 * REST_HOOP * steel().thickness / RADIUS


@dataclass(frozen=True)
class Stack:
    name: str
    layers: tuple[Layer, ...]


def stacks() -> list[Stack]:
    out = []
    for gas, carbon_m, lined in (("hydrogen", 0.33, True), ("methane", 0.18, False)):
        head = (liner(),) if lined else ()
        for with_al in (True, False):
            mid = (aluminium(),) if with_al else ()
            label = f"{gas}, {'with' if with_al else 'no'} Al"
            out.append(Stack(label, (*head, steel(), *mid, carbon(carbon_m), glass())))
    return out


def response(
    stack: Stack, amplitude: float, decay: float, cells_thin: int = 8
) -> list[tuple[Layer, NDArray[np.float64]]]:
    """Peak tension per layer and cell for `p(t) = amplitude exp(-t/decay)`, instant rise.

    Runs for eight decay times plus six round trips of the whole stack.
    """
    transit = sum(layer.thickness / layer.speed for layer in stack.layers)
    return response_to(
        stack, lambda t: amplitude * math.exp(-t / decay), 8.0 * decay + 6.0 * transit, cells_thin
    )


def response_to(
    stack: Stack, load: Load, duration: float, cells_thin: int = 8
) -> list[tuple[Layer, NDArray[np.float64]]]:
    """Peak tension per layer and cell [Pa], bias included, under inner-face pressure `load(t)`."""
    return simulate(stack, load, duration, cells_thin).peak


@dataclass(frozen=True)
class Run:
    """One wave run: peak tension per layer and cell [Pa, bias included], each layer's final
    momentum per unit area [N s/m^2, outward positive], and each interface's final gap [m]."""

    peak: list[tuple[Layer, NDArray[np.float64]]]
    momentum: list[float]
    gaps: list[float]


def simulate(
    stack: Stack, load: Load, duration: float, cells_thin: int = 8, preload: float = 0.0
) -> Run:
    """Waves through `stack` under inner-face pressure `load(t)`; contacts open past `preload`."""
    dt = min(layer.thickness / layer.speed for layer in stack.layers) / cells_thin
    cells = [max(2, round(layer.thickness / (layer.speed * dt))) for layer in stack.layers]
    steps = int(duration / dt)
    z = [layer.impedance for layer in stack.layers]
    r = [np.zeros(n) for n in cells]
    s = [np.zeros(n) for n in cells]
    bias = _bias_profiles(stack, cells)
    peak = [np.full(n, -math.inf) for n in cells]
    last = len(cells) - 1
    gap = [0.0] * last
    for k in range(steps):
        p = load(k * dt)
        out_r = [float(a[-1]) for a in r]
        out_s = [float(a[0]) for a in s]
        for i in range(len(cells)):
            r[i][1:] = r[i][:-1]
            s[i][:-1] = s[i][1:]
        r[0][0] = -p - out_s[0]
        s[last][-1] = -out_r[last]
        for i in range(last):
            za, zb = z[i], z[i + 1]
            joined = 2 * zb / (za + zb) * out_r[i] + (za - zb) / (za + zb) * out_s[i + 1]
            if stack.layers[i].bonded or (gap[i] <= 0.0 and joined + out_s[i + 1] <= preload):
                gap[i] = 0.0
                r[i + 1][0] = joined
                s[i][-1] = (zb - za) / (za + zb) * out_r[i] + 2 * za / (za + zb) * out_s[i + 1]
                continue
            # Open (or opening): both faces free; the gap grows at u_right - u_left.
            r[i + 1][0] = -out_s[i + 1]
            s[i][-1] = -out_r[i]
            gap[i] = max(0.0, gap[i] + (2 * out_s[i + 1] / zb + 2 * out_r[i] / za) * dt)
        for i in range(len(cells)):
            np.maximum(peak[i], r[i] + s[i], out=peak[i])
    # rho u dx = rho (s - r)/(rho c) c dt = (s - r) dt
    momentum = [float(np.sum(s[i] - r[i]) * dt) for i in range(len(cells))]
    return Run([(layer, peak[i] + bias[i]) for i, layer in enumerate(stack.layers)], momentum, gap)


def _bias_profiles(stack: Stack, cells: list[int]) -> list[NDArray[np.float64]]:
    """Radial bias: 0 at the steel's inner face, `interface_bias` at its back, through any
    aluminium, then falling linearly to 0 across the carbon. The liner and glass carry none."""
    b = interface_bias()
    out = []
    for layer, n in zip(stack.layers, cells, strict=True):
        x = (np.arange(n) + 0.5) / n
        if layer.name == "Cr-Mo":
            out.append(b * x)
        elif layer.name == "aluminium":
            out.append(np.full(n, b))
        elif layer.name == "IM7/8552":
            out.append(b * (1.0 - x))
        else:
            out.append(np.zeros(n))
    return out


@dataclass(frozen=True)
class Row:
    stack: str
    decay_us: float
    layer: str
    tension_per_gpa_mpa: float
    bias_mpa: float
    fail_spike_gpa: float


def evaluate(stack: Stack, decay: float) -> list[Row]:
    """Tension per GPa of spike in each layer, and the spike that breaks it.

    Linear in amplitude apart from the bias, so it is run at 1 GPa and at 2 GPa: the difference
    gives tension per GPa, and the intercept gives the bias at the worst point.
    """
    one = response(stack, 1e9, decay)
    two = response(stack, 2e9, decay)
    rows = []
    for (layer, p1), (_, p2) in zip(one, two, strict=True):
        worst = int(np.argmax(p2))
        per = float(p2[worst] - p1[worst])
        bias = float(2 * p1[worst] - p2[worst])
        fail = (
            (layer.strength - bias) / per * 1e9
            if per > 0 and math.isfinite(layer.strength)
            else math.inf
        )
        rows.append(Row(stack.name, decay * 1e6, layer.name, per / 1e6, bias / 1e6, fail / 1e9))
    return rows


def sedov_spike(energy: float, density: float, radius: float, gamma: float) -> tuple[float, float]:
    """(incident, reflected) strong-shock pressure [Pa] at the wall, Sedov arrival speed.

    `U = (2/5) xi^(5/2) (E/rho)^(1/2) R^(-3/2)` with xi = 1.033 at gamma 1.4 (alpha = 0.851,
    Kamm and Timmes) and ~1.0 lower; the
    reflected pressure is `(3 gamma - 1)/(gamma - 1)` times the incident.
    """
    xi = 1.033 if gamma >= 1.4 else 1.0
    speed = 0.4 * xi**2.5 * math.sqrt(energy / density) * radius**-1.5
    incident = 2.0 * density * speed**2 / (gamma + 1.0)
    return incident, incident * (3.0 * gamma - 1.0) / (gamma - 1.0)


def main() -> None:
    energy = 0.5 * 2.5 * 75_000.0**2
    print("Reflected-spike estimate (Sedov arrival, strong normal reflection):")
    for gas, rho in (("hydrogen", 3.34), ("methane", 3.78)):
        for gamma in (1.4, 1.2):
            inc, refl = sedov_spike(energy, rho, RADIUS, gamma)
            print(
                f"  {gas:8s} gamma {gamma}: incident {inc / 1e5:6.0f} bar, reflected "
                f"{refl / 1e9:.1f} GPa"
            )
    print(
        f"\nat-rest radial bias at the steel/overwrap interface: {interface_bias() / 1e6:.1f} MPa\n"
    )
    rows = [
        row for st in stacks() for decay in (10e-6, 30e-6, 100e-6) for row in evaluate(st, decay)
    ]
    print(f"{'stack':22s} {'decay':>6} {'layer':12s} {'MPa/GPa':>8} {'bias':>6} {'fails at':>9}")
    for r in rows:
        print(
            f"{r.stack:22s} {r.decay_us:5.0f}us {r.layer:12s} {r.tension_per_gpa_mpa:8.0f} "
            f"{r.bias_mpa:6.1f} {r.fail_spike_gpa:8.2f} GPa"
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
