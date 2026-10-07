"""A dry (unimpregnated) fibre wrap over a steel liner, under the solved blast spike.

`wall_waves.py` showed that a bonded overwrap fails in its resin: the reflected spike pulls the
plies apart at 0.12-0.22 GPa. A dry wrap has no resin to fail. Its wound layers press on each
other but cannot pull, so a reflected tension pulse lifts them off instead. Each lifted layer is
a ring: its hoop tension stops it and draws it back. What limits the wrap is then the **hoop
strain** of the flung layers, against the fibre's breaking strain.

# The model

A radial chain of thin layers, each a ring of the chamber's local radius `r` (the cylinder):

- **Steel liner:** `n_s` bonded cells, through-thickness stiffness `M/dt` (uniaxial-strain
  modulus `M`), carrying tension.
- **Dry wrap:** `n_w` wound layers in contact only. The force between neighbours is
  `max(0, b + k_c (u_i - u_{i+1}))`, with the through-thickness stiffness `k_c = E_z / t` of
  packed dry fibre, which is not well known; it is swept.
- **Hoop:** every layer feels `-(E_h t / r^2) u` per unit area, plus the constant force that
  holds the rest state in equilibrium.
- **Preload `b`:** the autofrettage contact compression, 16.7 MPa at the steel/wrap interface
  (as `wall_waves.interface_bias`), rising linearly from zero across the steel and falling
  linearly to zero across the wrap.

The inner face carries the solved wall pressure: the 10 kg, 2.5 m long plug at the cylinder's
hardest-hit station (`vessel_blast_long_plug.jsonl`). It is scaled to each chamber volume as the
blast scales:
- pressure as `20/V` against the solved quasi-static pressure;
- time as `(V/20)^(1/3)`;
- times a grid factor of 1 to 2, because the 1 cm spike is not converged.

**Sizing is by strain compatibility, as the parent's autofrettage rule requires.** The 10 mm steel
liner is bonded in strain to the wrap, and may swing only +/-0.75 Y, a total strain of
`1.5 Y / E` = 0.42%. The wrap must stiffen the wall so the design peak `D p_QSP` (D = 1.7) strains
it by no more than that: `sum(E_h t) = D p_QSP r / 0.0042`. A strength allowable alone would let a
compliant wrap strain the wall far past the steel's yield. The preload is the liner's autofrettage
compression, carried by the wrap: `0.75 Y t_s / r`, about 2.5 MPa.

**Material values are datasheet-class, not sourced here.**
- Fibre moduli and breaking strains: Kevlar 49, 112 GPa and 2.4%; Vectran HT, ~65 GPa and ~3.3%.
- Packing fraction 0.7, so the wrap is ~1000 kg/m^3 with hoop moduli 78 and 45 GPa.
- Through-thickness stiffness of packed dry fibre, `E_z`, is swept over 1, 3 and 10 GPa.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

HISTORIES = Path("data/results/walled_nozzle/near_term/vessel_blast_long_plug.jsonl")
OUTPUT = Path("data/results/walled_nozzle/near_term/dry_wrap.csv")
LOAD_CASE = "10 kg, 2.5 m, back"
MODEL_QSP = 69.7e6
"""The domed chamber's quasi-static pressure in the gamma-law runs [Pa] (`vessel_blast_domed`)."""
REAL_QSP_20 = 49.59e6
"""The methane 7000 K chamber's solved pressure at 20 m^3 [Pa]."""
R_C_20 = 1.4
"""Cylinder radius of the 20 m^3 domed chamber [m]."""

STEEL_T = 0.010
STEEL_M = 255.8e9
STEEL_E = 200.0e9
STEEL_RHO = 7830.0
STEEL_ALLOW = 562.5e6
DESIGN_D = 1.7
ZETA = 0.05


@dataclass(frozen=True)
class Fibre:
    name: str
    hoop_modulus: float  # packed wrap [Pa]
    density: float  # packed wrap [kg/m^3]
    breaking_strain: float


KEVLAR = Fibre("Kevlar 49", 0.7 * 112e9, 1000.0, 0.024)
VECTRAN = Fibre("Vectran HT", 0.7 * 65e9, 1000.0, 0.033)
CARBON_TOW = Fibre("dry carbon", 0.7 * 230e9, 1260.0, 0.021)
"""T700-class tow wound dry: stiff, but brittle and abrasion-prone without a matrix."""
SWING = 1.5 * 750e6 / STEEL_E
"""The liner's allowed strain swing, -0.75 Y to +0.75 Y (Y = 750 MPa)."""


@dataclass(frozen=True)
class Wall:
    radius: float
    wrap: Fibre
    wrap_t: float
    e_z: float
    n_s: int = 10
    n_w: int = 40
    bonded: bool = False

    def layers(self) -> tuple[NDArray[np.float64], ...]:
        """Per-layer `(thickness, density, hoop modulus)` and per-interface contact stiffness."""
        ts = np.r_[np.full(self.n_s, STEEL_T / self.n_s), np.full(self.n_w, self.wrap_t / self.n_w)]
        rho = np.r_[np.full(self.n_s, STEEL_RHO), np.full(self.n_w, self.wrap.density)]
        eh = np.r_[np.full(self.n_s, STEEL_E), np.full(self.n_w, self.wrap.hoop_modulus)]
        # Interface j joins layer j and j+1 (n-1 interfaces); stiffness from the two half-cells.
        mod = np.r_[np.full(self.n_s, STEEL_M), np.full(self.n_w, self.e_z)]
        half = 0.5 * ts / mod
        kc = 1.0 / (half[:-1] + half[1:])
        return ts, rho, eh, kc


def sized_wrap(qsp: float, radius: float, fibre: Fibre) -> float:
    """Wrap thickness [m] that holds the design peak to the liner's strain swing."""
    needed = DESIGN_D * qsp * radius / SWING
    return max(0.0, (needed - STEEL_E * STEEL_T) / fibre.hoop_modulus)


def preload(radius: float) -> float:
    """The liner's autofrettage compression carried by the wrap [Pa]: `0.75 Y t_s / r`."""
    return 0.75 * 750e6 * STEEL_T / radius


@dataclass(frozen=True)
class Response:
    peak_wrap_strain: float
    rest_wrap_strain: float
    max_gap_mm: float
    outer_speed: float
    steel_radial_tension_mpa: float
    peak_steel_hoop_mpa: float


def simulate(wall: Wall, load: NDArray[np.float64], dt_load: float, t_end: float) -> Response:
    """Integrate the layer chain under the inner-face pressure history `load` (sampled at
    `dt_load`, held at its last value beyond)."""
    ts, rho, eh, kc = wall.layers()
    n = ts.size
    m = rho * ts
    kh = eh * ts / wall.radius**2
    # Preload: interface bias rises linearly across the steel to b0, then falls across the wrap.
    n_s = wall.n_s
    b0 = preload(wall.radius)
    bias = np.r_[
        b0 * np.arange(1, n_s) / n_s,
        b0 * (1.0 - np.arange(0, wall.n_w) / wall.n_w),
    ]
    # Hoop force at rest balances the interface biases: H_i = F_i - F_{i-1}.
    f_in = np.r_[0.0, bias]
    f_out = np.r_[bias, 0.0]
    h_rest = f_in - f_out  # inward hoop force per area needed at rest (positive = inward)
    rest_strain = h_rest * wall.radius / (eh * ts)
    bonded = np.r_[np.ones(n_s - 1, dtype=bool), np.full(wall.n_w, wall.bonded)]
    damp = 2.0 * ZETA * np.sqrt(kc * 0.5 * (m[:-1] + m[1:]))
    dt = 0.2 * float(np.min(np.sqrt(np.minimum(m[:-1], m[1:]) / kc)))
    u = np.zeros(n)
    v = np.zeros(n)
    peak_strain = -math.inf
    max_gap = 0.0
    outer_speed = 0.0
    steel_tension = 0.0
    steel_hoop = -math.inf
    steps = int(t_end / dt)
    for k in range(steps):
        t = k * dt
        p = float(np.interp(t, np.arange(load.size) * dt_load, load))
        rel = u[:-1] - u[1:]
        force = bias + kc * rel + damp * (v[:-1] - v[1:])
        force = np.where(bonded, force, np.maximum(force, 0.0))
        acc_force = np.r_[p, force] - np.r_[force, 0.0] - h_rest - kh * u
        v += dt * acc_force / m
        u += dt * v
        strain = rest_strain + u / wall.radius
        peak_strain = max(peak_strain, float(strain[n_s:].max()))
        steel_hoop = max(steel_hoop, float(STEEL_E * strain[:n_s].max()))
        gaps = np.where(bonded, 0.0, -(rel + bias / kc))
        max_gap = max(max_gap, float(gaps.max()))
        outer_speed = max(outer_speed, float(v[-1]))
        steel_tension = min(steel_tension, float(force[: n_s - 1].min()))
    return Response(
        peak_wrap_strain=peak_strain,
        rest_wrap_strain=float(rest_strain[n_s:].max()),
        max_gap_mm=max_gap * 1e3,
        outer_speed=outer_speed,
        steel_radial_tension_mpa=-steel_tension / 1e6,
        peak_steel_hoop_mpa=steel_hoop / 1e6,
    )


def solved_load() -> tuple[NDArray[np.float64], float]:
    """The long plug's hardest-hit cylinder station [Pa, model units] and its sample spacing."""
    for line in HISTORIES.open():
        rec = json.loads(line)
        if rec["deposition"] != LOAD_CASE:
            continue
        z = np.asarray(rec["station_z"], dtype=np.float64)
        p = np.asarray(rec["pressure_pa"], dtype=np.float64)
        t = np.asarray(rec["time_s"], dtype=np.float64)
        cyl = (z >= 0.7) & (z < float(rec["z_c"]))
        j = int(np.flatnonzero(cyl)[np.argmax(p[:, cyl].max(axis=0))])
        return p[:, j], float(t[1] - t[0])
    raise ValueError(f"{LOAD_CASE} not in {HISTORIES}")


@dataclass(frozen=True)
class Row:
    volume_m3: float
    fibre: str
    e_z_gpa: float
    grid_factor: float
    bonded: bool
    radius_m: float
    wrap_cm: float
    wall_t_per_m2: float
    spike_gpa: float
    peak_wrap_strain_pct: float
    rest_wrap_strain_pct: float
    breaking_strain_pct: float
    max_gap_mm: float
    outer_speed_m_s: float
    steel_radial_tension_mpa: float
    peak_steel_hoop_mpa: float


def evaluate(volume: float, fibre: Fibre, e_z: float, grid: float, bonded: bool = False) -> Row:
    raw, dt_raw = solved_load()
    stretch = (volume / 20.0) ** (1.0 / 3.0)
    scale = REAL_QSP_20 / MODEL_QSP * 20.0 / volume * grid
    load = raw * scale
    qsp = REAL_QSP_20 * 20.0 / volume
    radius = R_C_20 * stretch
    wall = Wall(radius, fibre, sized_wrap(qsp, radius, fibre), e_z, bonded=bonded)
    dt_load = dt_raw * stretch
    t_end = load.size * dt_load + 4.0e-3 * stretch
    res = simulate(wall, load, dt_load, t_end)
    return Row(
        volume_m3=volume,
        fibre=fibre.name,
        e_z_gpa=e_z / 1e9,
        grid_factor=grid,
        bonded=bonded,
        radius_m=radius,
        wrap_cm=wall.wrap_t * 100,
        wall_t_per_m2=(STEEL_RHO * STEEL_T + fibre.density * wall.wrap_t) / 1e3,
        spike_gpa=float(load.max()) / 1e9,
        peak_wrap_strain_pct=res.peak_wrap_strain * 100,
        rest_wrap_strain_pct=res.rest_wrap_strain * 100,
        breaking_strain_pct=fibre.breaking_strain * 100,
        max_gap_mm=res.max_gap_mm,
        outer_speed_m_s=res.outer_speed,
        steel_radial_tension_mpa=res.steel_radial_tension_mpa,
        peak_steel_hoop_mpa=res.peak_steel_hoop_mpa,
    )


def main() -> None:
    rows = []
    for volume in (40.0, 60.0, 80.0, 120.0, 160.0):
        for fibre in (VECTRAN, KEVLAR, CARBON_TOW):
            for grid in (1.0, 2.0):
                rows.append(evaluate(volume, fibre, 3e9, grid))
    for e_z in (1e9, 10e9):
        rows.append(evaluate(80.0, VECTRAN, e_z, 2.0))
    rows.append(evaluate(80.0, VECTRAN, 3e9, 2.0, bonded=True))
    print(
        f"{'V':>4} {'fibre':10s} {'Ez':>3} {'grid':>4} {'bond':>4} {'r':>4} {'wrap':>5} "
        f"{'t/m2':>5} "
        f"{'spike':>6} {'strain%':>7} {'rest%':>5} {'break%':>6} {'gap mm':>6} {'v_out':>6} "
        f"{'steel r':>7} {'steel h':>7}"
    )
    for r in rows:
        print(
            f"{r.volume_m3:4.0f} {r.fibre:10s} {r.e_z_gpa:3.0f} {r.grid_factor:4.1f} "
            f"{int(r.bonded):4d} {r.radius_m:4.2f} {r.wrap_cm:4.1f}c {r.wall_t_per_m2:5.3f} "
            f"{r.spike_gpa:5.2f}G "
            f"{r.peak_wrap_strain_pct:7.2f} {r.rest_wrap_strain_pct:5.2f} "
            f"{r.breaking_strain_pct:6.1f} {r.max_gap_mm:6.2f} {r.outer_speed_m_s:6.0f} "
            f"{r.steel_radial_tension_mpa:6.0f}M {r.peak_steel_hoop_mpa:6.0f}M"
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
