"""The chamber wall's hoop overshoot under the solved blast, with neighbouring wall connected.

`dry_wrap.py` and `vessel_blast.py` load each band of wall as a free ring, unconnected to the
wall beside it. Under the blast's short, local spikes that puts the steel past yield in hoop, an
effective overshoot D of 3-6 against the sizing's 1.7. A real shell is connected. A band struck
harder than its neighbours is held back by bending to them, over the shell's characteristic
length `sqrt(R h) / (3(1 - nu^2))^(1/4)`. This module adds that coupling.

# The model

The axisymmetric shell as a beam on an elastic foundation along the meridian `s`, the standard
reduction for axisymmetric shell response:

    m w_tt + c w_t + D w_ssss + (sum(E t) / R(s)^2) w = p(s, t)

- `w`: radial displacement;
- `m`: areal mass;
- `D`: bending stiffness;
- `R(s)`: the local radius.

The meridian runs from the domed port head, along the cylinder, to the elliptical nose. Both ends
are clamped (`w = w_s = 0`), at the port ring and the throat insert. The meridional curvature of
the head and nose is neglected. The pressure `p(s, t)` is the solved wall history of the long plug
(10 kg, 2.5 m, back-weighted; 0.5 cm cells), scaled to the real methane chamber. For a larger
volume it is scaled as the blast scales: pressure as `20/V`, lengths and times as `(V/20)^(1/3)`.

Each wall gets the free-ring answer (`D = 0`) too, so the coupling's effect shows directly. The
reported **effective overshoot** is the peak hoop strain `w/R` over the static strain at the
quasi-static pressure, `p_QSP R / sum(E t)`. The walls are sized for 1.7 of that.

Walls:
- **monolithic steel, 20 m^3:** sized `t = 1.7 p_QSP R_c / (0.75 Y)`.
- **multilayer steel, 20 m^3:** the same steel in 8 unbonded shrink-fitted shells. It has the
  same hoop stiffness and 1/64 of the bending stiffness.
- **dry Kevlar 49 over 10 mm steel, 80 m^3:** sized as in `dry_wrap.py`. Only the liner bends;
  the wound layers slip.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import dry_wrap as dw

HISTORIES = Path("data/results/walled_nozzle/near_term/vessel_blast_long_plug_fine.jsonl")
OUTPUT = Path("data/results/walled_nozzle/near_term/shell_response.csv")
NU = 0.3
E_STEEL = 200e9
RHO_STEEL = 7830.0
Y = 750e6
QSP_20 = 49.59e6
MODEL_QSP = 69.7e6
ZETA = 0.02
DS = 0.01


@dataclass(frozen=True)
class WallSpec:
    name: str
    volume: float
    hoop_stiffness: float  # sum(E t) [N/m]
    bending: float  # D [N m]
    areal_mass: float  # [kg/m^2]


def plate_bending(e: float, t: float) -> float:
    return e * t**3 / (12.0 * (1.0 - NU**2))


def walls() -> list[WallSpec]:
    t_steel = 1.7 * QSP_20 * 1.4 / (0.75 * Y)
    mono = WallSpec(
        "monolithic steel",
        20.0,
        E_STEEL * t_steel,
        plate_bending(E_STEEL, t_steel),
        RHO_STEEL * t_steel,
    )
    multi = WallSpec(
        "multilayer steel x8",
        20.0,
        E_STEEL * t_steel,
        8 * plate_bending(E_STEEL, t_steel / 8),
        RHO_STEEL * t_steel,
    )
    v = 80.0
    r = 1.4 * (v / 20.0) ** (1.0 / 3.0)
    t_w = dw.sized_wrap(QSP_20 * 20.0 / v, r, dw.KEVLAR)
    wrap = WallSpec(
        "Kevlar dry wrap + 10 mm steel",
        v,
        E_STEEL * dw.STEEL_T + dw.KEVLAR.hoop_modulus * t_w,
        plate_bending(E_STEEL, dw.STEEL_T),
        RHO_STEEL * dw.STEEL_T + dw.KEVLAR.density * t_w,
    )
    return [mono, multi, wrap]


def load_field() -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """`(s, R(s), p[t, s] in real Pa at 20 m^3, dt)` on a uniform meridian grid."""
    rec = json.loads(HISTORIES.read_text().splitlines()[-1])
    z = np.asarray(rec["station_z"], dtype=np.float64)
    r = np.asarray(rec["station_r"], dtype=np.float64)
    p = np.asarray(rec["pressure_pa"], dtype=np.float64)
    t = np.asarray(rec["time_s"], dtype=np.float64)
    side = z > 0.0
    z, r, p = z[side], r[side], p[:, side]
    s_station = np.r_[0.0, np.cumsum(np.hypot(np.diff(z), np.diff(r)))]
    s = np.arange(0.0, s_station[-1], DS)
    radius = np.interp(s, s_station, r)
    field = np.array([np.interp(s, s_station, row) for row in p]) * QSP_20 / MODEL_QSP
    return s, radius, field, float(t[1] - t[0])


def respond(
    spec: WallSpec,
    s: NDArray[np.float64],
    radius: NDArray[np.float64],
    field: NDArray[np.float64],
    dt_field: float,
    t_end: float,
    bending: bool = True,
) -> float:
    """Peak hoop strain `max w/R` over the run."""
    n = s.size
    profile = respond_profile(
        np.full(n, spec.hoop_stiffness),
        np.full(n, spec.bending),
        np.full(n, spec.areal_mass),
        s,
        radius,
        field,
        dt_field,
        t_end,
        bending,
    )
    return float(profile.max())


def respond_profile(
    hoop_stiffness: NDArray[np.float64],
    bending_stiffness: NDArray[np.float64],
    areal_mass: NDArray[np.float64],
    s: NDArray[np.float64],
    radius: NDArray[np.float64],
    field: NDArray[np.float64],
    dt_field: float,
    t_end: float,
    bending: bool = True,
) -> NDArray[np.float64]:
    """Peak hoop strain `w/R` at each meridian node, for a wall that varies along `s`.

    `hoop_stiffness` is `sum(E t)` [N/m], `bending_stiffness` is `D` [N m] and `areal_mass` is in
    kg/m^2, one value per node. The bending term is the second derivative of `D w_ss`, which is
    the uniform fourth-difference stencil exactly when `D` is constant. Damping is
    `2 zeta omega_max m`, as before.
    """
    n = s.size
    ds = float(s[1] - s[0])
    k = hoop_stiffness / radius**2
    m = areal_mass
    d = bending_stiffness if bending else np.zeros(n)
    omega = math.sqrt(float((k / m).max()))
    dt = 0.1 / omega
    if float(d.max()) > 0.0:
        dt = min(dt, 0.2 * ds**2 / math.sqrt(float((d / m).max())))
    c = 2.0 * ZETA * omega * m
    # D at the n + 2 nodes where w_ss is formed (the ends extend the end values).
    d_pad = np.r_[d[0], d, d[-1]]
    w = np.zeros(n)
    v = np.zeros(n)
    tail = field[-max(1, field.shape[0] // 10) :].mean(axis=0)
    peak = np.zeros(n)
    steps = int(t_end / dt)
    for step in range(steps):
        t = step * dt
        x = t / dt_field
        i = int(x)
        inside = i + 1 < field.shape[0]
        p = field[i] + (x - i) * (field[i + 1] - field[i]) if inside else tail
        # Clamped ends: w_0 = 0, with mirrored ghosts w_-1 = w_1 and w_-2 = w_2.
        wp = np.r_[w[2], w[1], w, w[-2], w[-3]]
        moment = d_pad * (wp[:-2] - 2 * wp[1:-1] + wp[2:]) / ds**2
        b4 = (moment[:-2] - 2 * moment[1:-1] + moment[2:]) / ds**2
        acc = (p - k * w - b4 - c * v) / m
        acc[0] = acc[-1] = 0.0
        v += dt * acc
        w += dt * v
        w[0] = w[-1] = 0.0
        np.maximum(peak, w / radius, out=peak)
    return peak


@dataclass(frozen=True)
class Row:
    wall: str
    volume_m3: float
    char_length_m: float
    static_strain_pct: float
    peak_ring_pct: float
    peak_shell_pct: float
    overshoot_ring: float
    overshoot_shell: float


def evaluate(spec: WallSpec) -> Row:
    s, radius, field, dt_field = load_field()
    stretch = (spec.volume / 20.0) ** (1.0 / 3.0)
    s = s * stretch
    radius = radius * stretch
    field = field * 20.0 / spec.volume
    dt_field *= stretch
    t_end = field.shape[0] * dt_field + 3.0e-3 * stretch
    qsp = QSP_20 * 20.0 / spec.volume
    r_c = float(radius.max())
    static = qsp * r_c / spec.hoop_stiffness
    ring = respond(spec, s, radius, field, dt_field, t_end, bending=False)
    shell = respond(spec, s, radius, field, dt_field, t_end, bending=True)
    h_eq = math.sqrt(12.0 * (1 - NU**2) * spec.bending / spec.hoop_stiffness)
    return Row(
        wall=spec.name,
        volume_m3=spec.volume,
        char_length_m=math.sqrt(r_c * h_eq) / (3.0 * (1.0 - NU**2)) ** 0.25,
        static_strain_pct=static * 100,
        peak_ring_pct=ring * 100,
        peak_shell_pct=shell * 100,
        overshoot_ring=ring / static,
        overshoot_shell=shell / static,
    )


def main() -> None:
    rows = [evaluate(spec) for spec in walls()]
    print(
        f"{'wall':32s} {'V':>4} {'l_c m':>6} {'static%':>8} {'ring%':>7} {'shell%':>7} "
        f"{'D ring':>7} {'D shell':>7}"
    )
    for r in rows:
        print(
            f"{r.wall:32s} {r.volume_m3:4.0f} {r.char_length_m:6.3f} {r.static_strain_pct:8.3f} "
            f"{r.peak_ring_pct:7.3f} {r.peak_shell_pct:7.3f} {r.overshoot_ring:7.2f} "
            f"{r.overshoot_shell:7.2f}"
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
