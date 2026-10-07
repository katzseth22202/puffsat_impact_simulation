"""Which wall, where: zone-by-zone sizing of the chosen chamber against hoop, spall and cracks.

The chosen shape (`docs/walled_nozzle_chamber_blast.md` §3c) is the 40 m^3 long cone with an
r 2.2 m bulge around the plug, closed by a 2 m taper (~63 m^3). Its solved wall history (rod and
10 kg plug as material, 0.5 cm cells, `--shape-40b`) is scaled to the real methane chamber. The
wall is split into zones along the meridian: port head, the bulge band (ramp, plateau, taper),
cone and throat. Each zone takes one of three walls:

- **steel:** monolithic autofrettaged Cr-Mo;
- **multilayer:** the same steel in 8 shrink-fitted shells in contact only: the same hoop
  stiffness and mass, 1/64 of the bending stiffness;
- **wrap:** a 10 mm Cr-Mo liner under a dry Kevlar 49 wrap (`dry_wrap.KEVLAR`), unbonded.
- **maraging:** monolithic maraging 300 behind the alumina hydrogen barrier. Its spall strength
  (4.1 GPa, the low end of [M17]) and toughness (`K_Ic` 66.5 MPa m^0.5, NI76) are higher than
  Cr-Mo's.
  - It is rated "extreme" in hydrogen, so it is scored on the air crack curve only: the
    barrier must work.
  - Its allowed hoop swing is taken as 1.2 GPa (amplitude 0.6 GPa, about NI76's notched life
    at 10^4 cycles), not +/-0.75 Y; an assumption.

# Hoop

Every zone is first sized for the design peak `D p_QSP` with D = 1.7 at its own radius. The
allowed peak is the autofrettaged swing `dry_wrap.SWING` (-0.75 Y to +0.75 Y, 0.5625%). The parent's
`sec:carbon_overwrap` rule is used for every wall. `shell_response.respond_profile` then solves the
coupled shell under the real history. Any zone that strains past the swing is thickened by the
ratio, and the solve is repeated until every zone fits. The wrap zones thicken their wrap; the
steel zones their steel.

# Spall and cracks

For each zone the hardest-hit stations' histories run through `wall_waves.simulate`, with
contact interfaces where the wall is unbonded:
- steel: one bonded layer;
- multilayer: 8 shells in contact;
- wrap: the liner in contact with a wrap of through-thickness impedance ~2.8 MPa s/m.

The worst tension in any steel layer grows an embedded penny crack parallel to the face (the
plane the reflected tension opens). It starts at 0.5 mm radius, `dK = (2/pi) sigma sqrt(pi a)`,
R = 0, on the air curve and on the ASME CC2938 hydrogen curve at the quasi-static pressure. It
fails at `K = 45` MPa m^0.5 (the parent's hydrogen floor) or at half the layer's thickness. The
requirement is ~4,200 pulses: one departure of a 500-600 t stack at 5.43 km/s.

# Fling

An unbonded wall cannot spall, because it separates instead (the momentum trap), and the
separated outer layers fly outward until their hoop stops them. `dry_wrap.simulate`, the radial
contact chain of rings, gives the peak hoop strain of the flung layers at each zone's hardest-hit
station:
- wrap: the Kevlar wrap, against its 2.4% breaking strain;
- multilayer: the 8 steel shells, approximated as the 10 mm liner and 7 shells, against the
  autofrettaged swing.

Each station is a free ring, with no coupling along the wall, so the fling is overstated. In the
shell model, coupling cut the overshoot by 25-40%.

# Limits

- 1-D through the thickness, elastic. The spike's ~5 cm footprint would spread in 2-D, so the
  tensions are upper bounds.
- The spikes are 0.5 cm values sampled every 2 µs. Peaks rose 15-25% per halving of the cell, so
  the hydrogen-curve lives (exponent 3.7-6.5) may be ~2-3x shorter.
- The meridional curvature of the head, bulge and nose is neglected in the shell, as before.
- **The wave model runs over the spike only** (`SPIKE_WINDOW`), with the pressure measured from
  its value at the window's start. It has no hoop restoring force. Over the whole history the
  gas would push separated shells apart without limit and slam them back together, an
  artifact. Spall happens within one or two round trips (36 µs through 105 mm), and in 100 µs
  the ~1.5 ms ring barely moves. Dropping the earlier quasi-static compression is conservative.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections.abc import Callable
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import chamber_fatigue as cf
from puffsat.walled_nozzle import dry_wrap as dw
from puffsat.walled_nozzle import shell_response as sr
from puffsat.walled_nozzle import wall_waves as ww

HISTORIES = Path("data/results/walled_nozzle/near_term/vessel_blast_shape_40b.jsonl")
OUTPUT = Path("data/results/walled_nozzle/near_term/wall_zones.csv")
CASE = "bulge r 2.2 m z 1.4-2.8, 2 m taper, plug 1.5 m"
MODEL_QSP = 69.7e6
REAL_QSP_20 = 49.59e6
REQUIRED_PULSES = 4200
SHELLS = 8
LINER = dw.STEEL_T
WRAP_SPEED, WRAP_DENSITY = 2000.0, 1380.0
"""Dry aramid through the thickness (as `water_plate.skirt_delamination`), assumed."""
K_FAIL = 45.0
"""Cr-Mo: the parent's hydrogen floor for steels under 950 MPa [MPa m^0.5]."""


@dataclass(frozen=True)
class Steel:
    name: str
    modulus: float
    density: float
    speed: float
    spall: float
    k_fail: float
    swing: float
    hydrogen_qualified: bool


CRMO = Steel("Cr-Mo", dw.STEEL_E, dw.STEEL_RHO, 5900.0, ww.STEEL_SPALL, K_FAIL, dw.SWING, True)
MARAGING = Steel("maraging 300", 190e9, 8000.0, 5654.0, 4.1e9, 66.5, 1.2e9 / 190e9, False)


def steel_of(kind: str) -> Steel:
    return MARAGING if kind == "maraging" else CRMO


CRACK_START = 0.5e-3
ZONES: tuple[tuple[str, float, float], ...] = (
    ("head", 0.0, 0.7),
    ("bulge ramp", 0.7, 1.4),
    ("bulge plateau", 1.4, 2.8),
    ("taper", 2.8, 4.8),
    ("cone", 4.8, math.inf),
)
"""(name, z from, z to) [m]; the throat is the cone's last 0.15 m, reported with the cone."""
BAND = ("bulge ramp", "bulge plateau", "taper")
SPALL_STATIONS = 6
SPIKE_WINDOW = (30e-6, 100e-6)
"""The wave model runs from this long before each station's peak to this long after [s]."""


@dataclass(frozen=True)
class History:
    """The chosen chamber's wall: station `z`, `r` [m], time [s], real pressure [t, station]."""

    z: NDArray[np.float64]
    r: NDArray[np.float64]
    time: NDArray[np.float64]
    pressure: NDArray[np.float64]
    volume: float
    energy: float = 1.0
    """Pulse energy relative to the 2.5 kg rod."""

    @property
    def qsp(self) -> float:
        return REAL_QSP_20 * 20.0 / self.volume * self.energy

    @property
    def length_scale(self) -> float:
        """`energy^(1/3)`: zone boundaries are set for the 2.5 kg chamber and scale with it."""
        return float(self.energy ** (1.0 / 3.0))

    def scaled(self, energy: float) -> History:
        """The same chamber for `energy` times the pulse, by exact cube-root scaling. Inviscid
        flow has no length scale, so lengths and times grow as `energy^(1/3)` and pressures are
        unchanged."""
        k = energy ** (1.0 / 3.0)
        return History(
            self.z * k, self.r * k, self.time * k, self.pressure, self.volume * energy, energy
        )


def load_history(path: Path = HISTORIES, case: str = CASE) -> History:
    for line in path.open():
        rec = json.loads(line)
        if rec["deposition"] != case:
            continue
        z = np.asarray(rec["station_z"], dtype=np.float64)
        r = np.asarray(rec["station_r"], dtype=np.float64)
        keep = z > 0.0  # the port face is the head's flat closure; the meridian starts at the wall
        # Model and real QSP both scale as 1/V at fixed energy, so one ratio serves any volume.
        scale = REAL_QSP_20 / MODEL_QSP
        p = np.asarray(rec["pressure_pa"], dtype=np.float64)[:, keep] * scale
        t = np.asarray(rec["time_s"], dtype=np.float64)
        return History(z[keep], r[keep], t, p, float(rec["volume_m3"]))
    raise ValueError(f"{case!r} not in {path}")


def zone_of(z: float) -> str:
    return next(name for name, lo, hi in ZONES if lo <= z < hi)


@dataclass(frozen=True)
class Wall:
    """One wall concept: the kind of wall in each zone."""

    name: str
    kinds: dict[str, str]  # zone -> "steel" | "multilayer" | "wrap"


def walls() -> list[Wall]:
    every = [name for name, _, _ in ZONES]
    return [
        Wall("all steel", dict.fromkeys(every, "steel")),
        Wall("all multilayer", dict.fromkeys(every, "multilayer")),
        Wall("all dry wrap", dict.fromkeys(every, "wrap")),
        Wall(
            "multilayer band, wrap elsewhere",
            {z: ("multilayer" if z in BAND else "wrap") for z in every},
        ),
        Wall("steel band, wrap elsewhere", {z: ("steel" if z in BAND else "wrap") for z in every}),
        Wall("all maraging (barrier)", dict.fromkeys(every, "maraging")),
        Wall(
            "maraging band (barrier), wrap elsewhere",
            {z: ("maraging" if z in BAND else "wrap") for z in every},
        ),
    ]


def _initial_thickness(kind: str, qsp: float, radius: float) -> float:
    """Steel [m] for steel/multilayer, wrap [m] for wrap: sized for D = 1.7 at the swing."""
    if kind == "wrap":
        return dw.sized_wrap(qsp, radius, dw.KEVLAR)
    mat = steel_of(kind)
    return dw.DESIGN_D * qsp * radius / (mat.modulus * mat.swing)


def _section(kind: str, t: float) -> tuple[float, float, float]:
    """(hoop stiffness sum(E t) [N/m], bending D [N m], areal mass [kg/m^2]) for one zone."""
    if kind == "wrap":
        hoop = dw.STEEL_E * LINER + dw.KEVLAR.hoop_modulus * t
        return (
            hoop,
            sr.plate_bending(dw.STEEL_E, LINER),
            dw.STEEL_RHO * LINER + dw.KEVLAR.density * t,
        )
    mat = steel_of(kind)
    bend = sr.plate_bending(mat.modulus, t)
    if kind == "multilayer":
        bend = SHELLS * sr.plate_bending(mat.modulus, t / SHELLS)
    return mat.modulus * t, bend, mat.density * t


@dataclass(frozen=True)
class ZoneResult:
    wall: str
    zone: str
    kind: str
    radius_m: float
    thickness_mm: float  # steel, or wrap over the 10 mm liner
    tonnes: float
    peak_strain_pct: float
    overshoot: float  # zone peak strain over the zone's largest static strain
    tension_over_spall: float
    pulses_air: float
    pulses_h2: float
    fling_strain_pct: float  # flung layers' peak hoop strain (0 for a bonded wall)
    fling_allow_pct: float
    swing_pct: float
    needs_barrier: bool  # scored on the air curve: only valid behind a working barrier

    @property
    def passes(self) -> bool:
        cracks = self.pulses_air if self.needs_barrier else self.pulses_h2
        return (
            self.peak_strain_pct <= self.swing_pct * 1.001
            and cracks >= REQUIRED_PULSES
            and self.fling_strain_pct <= self.fling_allow_pct
        )


def _meridian(h: History) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Uniform meridian `s`, its radius and its `z`."""
    s_st = np.r_[0.0, np.cumsum(np.hypot(np.diff(h.z), np.diff(h.r)))]
    s = np.arange(0.0, s_st[-1], sr.DS)
    return s, np.interp(s, s_st, h.r), np.interp(s, s_st, h.z)


def size_hoop(
    wall: Wall, h: History, iterations: int = 6
) -> tuple[dict[str, float], NDArray[np.float64], NDArray[np.int64], NDArray[np.float64]]:
    """Zone thicknesses after the hoop iteration, with each node's peak strain, zone index and
    static strain."""
    s, radius, z = _meridian(h)
    s_st = np.r_[0.0, np.cumsum(np.hypot(np.diff(h.z), np.diff(h.r)))]
    field = np.array([np.interp(s, s_st, row) for row in h.pressure])
    dt = float(h.time[1] - h.time[0])
    names = [name for name, _, _ in ZONES]
    zone_ix = np.array(
        [names.index(zone_of(float(zz) / h.length_scale)) for zz in z], dtype=np.int64
    )
    t = {
        name: _initial_thickness(
            wall.kinds[name], h.qsp, float(radius[zone_ix == k].max(initial=1.4))
        )
        for k, name in enumerate(names)
    }
    strain = np.zeros(s.size)
    static = np.zeros(s.size)
    for _ in range(iterations):
        hoop, bend, mass = (np.zeros(s.size) for _ in range(3))
        for k, name in enumerate(names):
            sec = _section(wall.kinds[name], t[name])
            hoop[zone_ix == k], bend[zone_ix == k], mass[zone_ix == k] = sec
        strain = sr.respond_profile(
            hoop, bend, mass, s, radius, field, dt, field.shape[0] * dt + 3e-3
        )
        static = h.qsp * radius / hoop
        grown = False
        for k, name in enumerate(names):
            mask = zone_ix == k
            if not mask.any():
                continue
            kind = wall.kinds[name]
            swing = dw.SWING if kind == "wrap" else steel_of(kind).swing
            over = float(strain[mask].max()) / swing
            if over > 1.001:
                grown = True
                hoop_now = _section(wall.kinds[name], t[name])[0]
                target = hoop_now * over
                if wall.kinds[name] == "wrap":
                    t[name] = (target - dw.STEEL_E * LINER) / dw.KEVLAR.hoop_modulus
                else:
                    t[name] = target / steel_of(wall.kinds[name]).modulus
        if not grown:
            break
    return t, strain, zone_ix, static


def _stack(kind: str, t: float) -> ww.Stack:
    if kind in ("steel", "maraging"):
        mat = steel_of(kind)
        return ww.Stack(kind, (ww.Layer(mat.name, t, mat.speed, mat.density, mat.spall),))
    if kind == "multilayer":
        shell = ww.Layer("Cr-Mo", t / SHELLS, 5900.0, 7830.0, ww.STEEL_SPALL, bonded=False)
        last = ww.Layer("Cr-Mo", t / SHELLS, 5900.0, 7830.0, ww.STEEL_SPALL)
        return ww.Stack("multilayer", (*([shell] * (SHELLS - 1)), last))
    liner = ww.Layer("Cr-Mo", LINER, 5900.0, 7830.0, ww.STEEL_SPALL, bonded=False)
    return ww.Stack(
        "wrap", (liner, ww.Layer("wrap", max(t, 0.01), WRAP_SPEED, WRAP_DENSITY, math.inf))
    )


def crack_pulses(
    sigma_pa: float,
    layer_t: float,
    rate: Callable[[float, float, float], float],
    p_bar: float,
    k_fail: float = K_FAIL,
) -> float:
    """Pulses for an embedded penny crack parallel to the face to reach failure."""
    sigma = sigma_pa / 1e6
    if sigma <= 0.0:
        return math.inf
    a_c = min((k_fail * math.pi / (2.0 * sigma)) ** 2 / math.pi, 0.5 * layer_t)
    if a_c <= CRACK_START:
        return 0.0
    edges = np.geomspace(CRACK_START, a_c, 400)
    n = 0.0
    for a0, a1 in pairwise(edges):
        dk = 2.0 / math.pi * sigma * math.sqrt(math.pi * a0)
        n += (a1 - a0) / rate(dk, 0.0, p_bar)
    return n


def spall(kind: str, t: float, h: History, zone: str) -> tuple[float, float]:
    """(worst steel tension [Pa], the steel layer's thickness [m]) over the zone's hardest hits."""
    mask = np.array([zone_of(float(zz) / h.length_scale) == zone for zz in h.z])
    if not mask.any():
        return 0.0, t
    idx = np.flatnonzero(mask)
    top = idx[np.argsort(h.pressure[:, idx].max(axis=0))[::-1][:SPALL_STATIONS]]
    stack = _stack(kind, t)
    worst = 0.0
    for j in top:
        series = h.pressure[:, j]
        t_peak = float(h.time[int(np.argmax(series))])
        t0 = max(0.0, t_peak - SPIKE_WINDOW[0])
        base = float(np.interp(t0, h.time, series))

        def load(
            tt: float, series: NDArray[np.float64] = series, t0: float = t0, base: float = base
        ) -> float:
            return float(np.interp(t0 + tt, h.time, series)) - base

        run = ww.simulate(stack, load, t_peak + SPIKE_WINDOW[1] - t0, cells_thin=16)
        steels = (CRMO.name, MARAGING.name)
        worst = max(worst, max(float(a.max()) for la, a in run.peak if la.name in steels))
    layer_t = {"steel": t, "maraging": t, "multilayer": t / SHELLS, "wrap": LINER}[kind]
    return worst, layer_t


SHELL_RINGS = dw.Fibre("Cr-Mo shells", dw.STEEL_E, dw.STEEL_RHO, dw.SWING)
"""Steel shells as `dry_wrap` layers, failing at the autofrettaged swing."""


def fling(kind: str, t: float, h: History, zone: str) -> tuple[float, float]:
    """(peak flung-layer hoop strain, its allowable) at the zone's hardest-hit station."""
    if kind in ("steel", "maraging"):
        return 0.0, math.inf
    idx = np.flatnonzero([zone_of(float(zz) / h.length_scale) == zone for zz in h.z])
    if idx.size == 0:
        return 0.0, math.inf
    j = idx[int(np.argmax(h.pressure[:, idx].max(axis=0)))]
    if kind == "wrap":
        chain = dw.Wall(float(h.r[j]), dw.KEVLAR, t, 3e9)
        allow = dw.KEVLAR.breaking_strain
    else:
        chain = dw.Wall(
            float(h.r[j]), SHELL_RINGS, t - dw.STEEL_T, dw.STEEL_M, n_s=4, n_w=SHELLS - 1
        )
        allow = dw.SWING
    dt = float(h.time[1] - h.time[0])
    res = dw.simulate(chain, h.pressure[:, j], dt, float(h.time[-1]) + 4e-3)
    return res.peak_wrap_strain, allow


def evaluate(wall: Wall, h: History) -> list[ZoneResult]:
    t, strain, zone_ix, static = size_hoop(wall, h)
    _, radius, _ = _meridian(h)
    names = [name for name, _, _ in ZONES]
    out = []
    for k, name in enumerate(names):
        mask = zone_ix == k
        if not mask.any():
            continue
        kind = wall.kinds[name]
        _, _, mass = _section(kind, t[name])
        area = float(np.sum(2.0 * math.pi * radius[mask] * sr.DS))
        tension, layer_t = spall(kind, t[name], h, name)
        flung, allow = fling(kind, t[name], h, name)
        p_bar = h.qsp / 1e5
        mat = steel_of(kind)
        swing = dw.SWING if kind == "wrap" else mat.swing
        out.append(
            ZoneResult(
                wall=wall.name,
                zone=name,
                kind=kind,
                radius_m=float(radius[mask].max()),
                thickness_mm=t[name] * 1e3,
                tonnes=mass * area / 1e3,
                peak_strain_pct=float(strain[mask].max()) * 100,
                overshoot=float(strain[mask].max() / static[mask].max()),
                tension_over_spall=tension / mat.spall,
                pulses_air=crack_pulses(tension, layer_t, cf.rate_air, p_bar, mat.k_fail),
                pulses_h2=(
                    crack_pulses(tension, layer_t, cf.rate_hydrogen, p_bar, mat.k_fail)
                    if mat.hydrogen_qualified
                    else 0.0
                ),
                fling_strain_pct=flung * 100,
                fling_allow_pct=allow * 100,
                swing_pct=swing * 100,
                needs_barrier=not mat.hydrogen_qualified,
            )
        )
    return out


def main(argv: list[str] | None = None) -> None:
    """`[--histories PATH] [--case LABEL]`: score another solved chamber."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--histories", type=Path, default=HISTORIES)
    parser.add_argument("--case", default=CASE)
    parser.add_argument("--energy", type=float, default=1.0, help="pulse energy / 2.5 kg rod's")
    parser.add_argument("--suffix", default="", help="appended to the output name, e.g. fine")
    args = parser.parse_args(argv)
    h = load_history(args.histories, args.case)
    if args.energy != 1.0:
        h = h.scaled(args.energy)
    print(f"{args.case}, pulse energy x{args.energy:g}")
    print(f"chosen chamber {h.volume:.1f} m^3, QSP {h.qsp / 1e6:.1f} MPa, swing {dw.SWING:.4%}")
    rows = []
    for wall in walls():
        res = evaluate(wall, h)
        rows.extend(res)
        total = sum(r.tonnes for r in res)
        ok = all(r.passes for r in res)
        print(f"\n{wall.name}: {total:.1f} t, {'passes' if ok else 'FAILS'}")
        print(
            f"  {'zone':14s} {'kind':10s} {'r m':>5} {'t mm':>6} {'tonnes':>6} {'strain%':>7} "
            f"{'D':>5} {'x spall':>7} {'N air':>9} {'N H2':>9} {'fling%':>7} {'allow':>6}"
        )
        for r in res:
            print(
                f"  {r.zone:14s} {r.kind:10s} {r.radius_m:5.2f} {r.thickness_mm:6.0f} "
                f"{r.tonnes:6.1f} {r.peak_strain_pct:7.3f} {r.overshoot:5.2f} "
                f"{r.tension_over_spall:7.2f} {min(r.pulses_air, 9.99e8):9.3g} "
                f"{min(r.pulses_h2, 9.99e8):9.3g} {r.fling_strain_pct:7.2f} "
                f"{min(r.fling_allow_pct, 99.0):6.2f}"
            )
    tag = args.case.split(",")[0].replace(" ", "_")
    if args.energy != 1.0:
        tag += f"_energy{args.energy:g}"
    if args.suffix:
        tag += f"_{args.suffix}"
    output = OUTPUT.with_name(f"wall_zones_{tag}.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(ZoneResult.__dataclass_fields__)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([getattr(r, f) for f in fields])
    print(f"\nwrote {len(rows)} rows -> {output}")


if __name__ == "__main__":
    main()
