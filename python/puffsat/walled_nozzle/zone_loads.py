"""Wall forces per pulse, zone by zone: hydrocode ask B3 (2026-10-12), and the hydrogen fill (B5).

The bulged 212 m^3 chamber against unbulged 106 and 53 m^3 chambers (and the 80 m^3 long cone
already solved), all at the 5 kg rod, the same 10 kg plug, charge and throat. Each solved 2.5 kg
history (rod and plug as material, 0.5 cm cells) is read at 5 kg by cube-root scaling: pressures
unchanged, times and lengths x2^(1/3). Each is scaled to its own real quasi-static pressure at the
fixed charge (`blowdown.solve`), not the fixed-temperature `1/V` law.

Per zone (`wall_zones.ZONES`, the bulge's z-bands; in the unbulged chambers they are bands of
the cylinder and cone at the same z):

- **spike:** the worst station's peak, its overpressure impulse above the quasi-static pressure
  over `-30/+100 µs` of the peak (scaled), and the width at half its overpressure;
- **quasi-static history:** the zone's mean pressure over the QSP at the end of the 1.2 ms solve,
  then the first 10 ms of the 2-D blowdown (`vessel_blowdown_quick.jsonl`, sampled every 1 ms);
  past that the real-gas 0-D decay (`blowdown`);
- **hoop overshoot `D`:** `wall_zones.evaluate`, the coupled shell under the whole history, with
  each zone's wall sized to fit (dry wrap; all steel), and the wall's mass.
"""

from __future__ import annotations

import csv
import dataclasses
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import near_term as nt
from puffsat.walled_nozzle import wall_zones as wz

DIR = nt.OUTPUT_DIR
OUTPUT = DIR / "zone_loads.csv"
QS_OUTPUT = DIR / "zone_quasi_static.csv"
ENERGY = 2.0
"""The 5 kg rod over the 2.5 kg rod the histories were solved with."""
K5 = ENERGY ** (1.0 / 3.0)
SPIKE_WINDOW = (30e-6, 100e-6)


@dataclass(frozen=True)
class Case:
    """One solved chamber: where its 2.5 kg history is and how it is read at 5 kg."""

    name: str
    path: Path
    deposition: str
    #: Volume at 5 kg [m^3] and its fixed-charge chamber, for the real QSP.
    volume_5kg: float
    spec: bd.Spec
    #: The whole-blowdown 2-D run of this chamber (`vessel_blowdown.jsonl` case), if any.
    blowdown_case: str | None = None
    #: Pulse energy over the 2.5 kg rod's: 2 for the 5 kg methane chamber, 1 for hydrogen's.
    energy: float = ENERGY


def cases() -> list[Case]:
    a25 = bd.CH4_THROAT_5KG / K5**2
    ch4 = {s.volume: s for s in bd.methane_specs()}
    v80 = dataclasses.replace(
        ch4[106.0],
        name="CH4 80 m^3 unbulged",
        volume=80.0,
        contour=bd.long_cone(40.0, a25).scaled(K5),
    )
    return [
        Case(
            "CH4 212 m^3 bulged",
            DIR / "vessel_blast_shape_40c.jsonl",
            "bulge r 3.0 m z 1.4-2.8, 3 m taper, plug 1.5 m",
            212.0,
            ch4[212.0],
            "bulge r 3.0 m, 106 m^3 (212 at 5 kg)",
        ),
        Case(
            "CH4 212 m^3 bulged, 0.25 cm",
            DIR / "vessel_blast_shape_40c_fine.jsonl",
            "bulge r 3.0 m z 1.4-2.8, 3 m taper, plug 1.5 m",
            212.0,
            ch4[212.0],
            "bulge r 3.0 m, 106 m^3 (212 at 5 kg)",
        ),
        Case(
            "CH4 106 m^3 unbulged",
            DIR / "vessel_blast_asks_unbulged.jsonl",
            "unbulged 53 m^3 (106 at 5 kg), plug 0.8-2.3 m",
            106.0,
            ch4[106.0],
            "unbulged 53 m^3 (106 at 5 kg)",
        ),
        Case(
            "CH4 80 m^3 unbulged",
            DIR / "vessel_blast_nose_study_40.jsonl",
            "40 m^3 long cone, 10 kg 0.8-2.3 m",
            80.0,
            v80,
        ),
        Case(
            "CH4 53 m^3 unbulged",
            DIR / "vessel_blast_asks_unbulged.jsonl",
            "unbulged 26.5 m^3 (53 at 5 kg), plug 0.8-2.3 m",
            53.0,
            ch4[53.0],
            "unbulged 26.5 m^3 (53 at 5 kg)",
        ),
        Case(
            "H2 106 m^3 bulged (2.5 kg)",
            DIR / "vessel_blast_asks_hydrogen.jsonl",
            "hydrogen, bulge r 3.0 m z 1.4-2.8, 3 m taper, plug 1.5 m",
            106.0,
            bd.hydrogen_specs()[1],
            energy=1.0,
        ),
        *(
            Case(spec.name, DIR / "vessel_blast_compact.jsonl", label, spec.volume, spec)
            for spec, (label, *_) in zip(bd.compact_specs(), bd.COMPACT, strict=True)
        ),
        *(
            Case(
                spec.name,
                DIR / "vessel_blast_compact_hydrogen.jsonl",
                label,
                spec.volume,
                spec,
                energy=1.0,
            )
            for spec, (label, *_) in zip(bd.compact_specs("hydrogen_5500"), bd.COMPACT, strict=True)
            if "throat x2" in label
        ),
    ]


def history(case: Case) -> wz.History:
    """The case's history at its rod (`case.energy`), scaled to its real fixed-charge QSP."""
    h = wz.load_history(case.path, case.deposition, case.spec.proxy_volume)
    law = wz.REAL_QSP_20 * 20.0 / h.volume  # what load_history scaled to
    real = bd.solve(case.spec).pressure
    h = dataclasses.replace(h, pressure=h.pressure * real / law, qsp_real=real)
    return h.scaled(case.energy) if case.energy != 1.0 else h


@dataclass(frozen=True)
class ZoneLoad:
    case: str
    zone: str
    z_from: float
    z_to: float
    radius_max: float
    qsp_bar: float
    peak_gpa: float
    peak_z: float
    spike_impulse_kpas: float
    spike_width_us: float
    #: Zone-mean pressure over the QSP at the end of the 1.2 ms solve (x2^(1/3) at 5 kg).
    end_ratio: float
    overshoot_wrap: float
    tonnes_wrap: float
    passes_wrap: bool
    #: The dry wrap's flung-layer hoop strain against its 2.4% [%].
    fling_wrap_pct: float
    overshoot_steel: float
    tonnes_steel: float
    #: Solid Cr-Mo: back-face tension over spall, and crack life in hydrogen [pulses].
    spall_steel: float
    pulses_h2_steel: float


def _spike(h: wz.History, j: int) -> tuple[float, float]:
    """(overpressure impulse [Pa s], width at half overpressure [s]) of station `j`'s spike."""
    series = h.pressure[:, j] - h.qsp
    k = int(np.argmax(series))
    t_pk = float(h.time[k])
    ls = h.length_scale
    window = (h.time >= t_pk - SPIKE_WINDOW[0] * ls) & (h.time <= t_pk + SPIKE_WINDOW[1] * ls)
    dt = float(h.time[1] - h.time[0])
    impulse = float(np.sum(np.clip(series[window], 0.0, None)) * dt)
    # The contiguous run above half the overpressure, around the peak.
    lo = k
    while lo - 1 >= 0 and series[lo - 1] >= 0.5 * series[k]:
        lo -= 1
    hi = k
    while hi + 1 < series.size and series[hi + 1] >= 0.5 * series[k]:
        hi += 1
    return impulse, (hi - lo + 1) * dt


def zone_mask(h: wz.History, zone: str) -> NDArray[np.bool_]:
    return np.array([wz.zone_of(float(z) / h.length_scale) == zone for z in h.z])


def loads(case: Case) -> list[ZoneLoad]:
    h = history(case)
    walls = {w.name: w for w in wz.walls()}
    res = {
        name: {r.zone: r for r in wz.evaluate(walls[name], h, iterations=12)}
        for name in ("all dry wrap", "all steel")
    }
    wrap_pass = all(r.passes for r in res["all dry wrap"].values())
    out = []
    for name, lo, hi in wz.ZONES:
        mask = zone_mask(h, name)
        if not mask.any():
            continue
        idx = np.flatnonzero(mask)
        peaks = h.pressure[:, idx].max(axis=0)
        j = int(idx[int(np.argmax(peaks))])
        spikes = [_spike(h, int(i)) for i in idx]
        worst = max(range(len(idx)), key=lambda i: spikes[i][0])
        tail = h.time >= h.time[-1] - 0.1e-3 * h.length_scale
        out.append(
            ZoneLoad(
                case=case.name,
                zone=name,
                z_from=lo * h.length_scale,
                z_to=min(hi * h.length_scale, float(h.z.max())),
                radius_max=float(h.r[idx].max()),
                qsp_bar=h.qsp / 1e5,
                peak_gpa=float(peaks.max()) / 1e9,
                peak_z=float(h.z[j]),
                spike_impulse_kpas=spikes[worst][0] / 1e3,
                spike_width_us=spikes[worst][1] * 1e6,
                end_ratio=float(h.pressure[np.ix_(tail, idx)].mean()) / h.qsp,
                overshoot_wrap=res["all dry wrap"][name].overshoot,
                tonnes_wrap=res["all dry wrap"][name].tonnes,
                passes_wrap=wrap_pass,
                fling_wrap_pct=res["all dry wrap"][name].fling_strain_pct,
                overshoot_steel=res["all steel"][name].overshoot,
                tonnes_steel=res["all steel"][name].tonnes,
                spall_steel=res["all steel"][name].tension_over_spall,
                pulses_h2_steel=res["all steel"][name].pulses_h2,
            )
        )
    return out


QUICK_2D = bd.BLOWDOWN_2D.with_name("vessel_blowdown_quick.jsonl")
QS_UNTIL = 10e-3
"""The 2-D blowdown is read only this far (2.5 kg time): its immersed wall leaks mass beyond."""


def quasi_static(case: Case, path: Path = QUICK_2D) -> list[tuple[str, float, float]]:
    """`(zone, time [s] at 5 kg, zone-mean pressure / initial QSP)` from the 2-D blowdown.

    The 2-D run is gamma-law; its pressure is read against its own uniform QSP
    `(gamma - 1) E0 / V`, so the ratio carries over to the real chamber's QSP. Only its first
    `QS_UNTIL` is used (`blowdown.early_outflow`: the wall leaks a few percent by 20 ms); the
    real gas then empties on the 0-D clock.
    """
    if case.blowdown_case is None or not path.exists():
        return []
    rec = next(
        (json.loads(x) for x in path.open() if json.loads(x)["case"] == case.blowdown_case), None
    )
    if rec is None:
        return []
    qsp = (rec["gamma"] - 1.0) * rec["energy0_j"] / rec["volume_m3"]
    z = np.asarray(rec["station_z"])
    p = np.asarray(rec["wall_pressure_pa"])
    out = []
    for name, _, _ in wz.ZONES:
        mask = np.array([wz.zone_of(float(x)) == name for x in z]) & (z > 0.0)
        if not mask.any():
            continue
        for t, row in zip(rec["wall_time_s"], p, strict=True):
            if t <= QS_UNTIL:
                scale = case.energy ** (1.0 / 3.0)
                out.append((name, float(t) * scale, float(row[mask].mean()) / qsp))
    return out


def write(rows: Sequence[ZoneLoad], qs: Sequence[tuple[str, str, float, float]]) -> None:
    DIR.mkdir(parents=True, exist_ok=True)
    fields = list(ZoneLoad.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in rows:
            w.writerow([getattr(r, f) for f in fields])
    with QS_OUTPUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "zone", "time_ms_5kg", "p_over_qsp"])
        for c, zone, t, ratio in qs:
            w.writerow([c, zone, f"{1e3 * t:.2f}", f"{ratio:.4f}"])


def main(argv: list[str] | None = None) -> None:
    rows: list[ZoneLoad] = []
    qs: list[tuple[str, str, float, float]] = []
    for case in cases():
        if not case.path.exists():
            print(f"{case.name}: {case.path} not solved yet, skipped")
            continue
        try:
            got = loads(case)
        except ValueError as err:
            print(f"{case.name}: {err}")
            continue
        rows.extend(got)
        qs.extend((case.name, zone, t, ratio) for zone, t, ratio in quasi_static(case))
        print(
            f"\n{case.name}: QSP {got[0].qsp_bar:.0f} bar; dry wrap "
            f"{sum(r.tonnes_wrap for r in got):.1f} t "
            f"({'passes' if got[0].passes_wrap else 'FAILS'}),"
            f" all steel {sum(r.tonnes_steel for r in got):.1f} t"
        )
        print(
            f"  {'zone':14s} {'z m':>9} {'r m':>5} {'peak GPa':>8} {'I kPa s':>7} {'w us':>5} "
            f"{'end/QSP':>7} {'D wrap':>6} {'fling%':>6} {'D steel':>7} {'x spall':>7} "
            f"{'N H2':>8}"
        )
        for r in got:
            print(
                f"  {r.zone:14s} {r.z_from:4.1f}-{r.z_to:4.1f} {r.radius_max:5.2f} "
                f"{r.peak_gpa:8.2f} {r.spike_impulse_kpas:7.1f} {r.spike_width_us:5.0f} "
                f"{r.end_ratio:7.2f} {r.overshoot_wrap:6.2f} {r.fling_wrap_pct:6.2f} "
                f"{r.overshoot_steel:7.2f} {r.spall_steel:7.2f} "
                f"{min(r.pulses_h2_steel, 9.9e8):8.3g}"
            )
    write(rows, qs)
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
