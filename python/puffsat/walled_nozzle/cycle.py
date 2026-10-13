"""The chamber's cycle between rods: hydrocode asks B2 and B5's rate (2026-10-12).

Each pulse, the chamber must blow down, have films laid over both openings, take its next plug,
fill, and lock 20 ms before the next rod (parent ADR 0053, `templateArxiv.tex` "Thin membranes
could seal both openings"). The parent times the cycle to the **fill pressure**: the next charge,
cold, in the chamber. That gate is not physical. Two things are:

- **The throat film is laid over a residue, which it then traps.** Whatever is left when the
  throat is sealed joins the next charge. It is ~4,000 K gas, mostly atomic hydrogen, so it heats
  the cold fill, and the mixed charge's pressure is what the film must hold.
- **A film cannot be laid across a sonic jet.** The seal waits until the throat's flow is weak.
  Here that is posed as a seal pressure: the parent's fill pressure, 0.3 bar or 0.1 bar.

Per chamber and seal pressure this module reports the seal time, the trapped residue, the mixed
charge's temperature and pressure (two bounds: the residue's thermal energy only, or all of its
energy once it recombines towards the cold feed), and the throat film it takes to hold it:
`t = p r / (2 sigma)` over a hemispherical cap of the throat's radius (the parent's rule).

Mechanism times are **assumptions**, stated as such: film laying, injector flow, lock.

The fill gas's heat capacity is the EOS's own: pure CH4 (or H2) thermal energy
(`eos_methane.thermal_energy`), the cold feed being the zero of `e`.
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Sequence
from dataclasses import dataclass

from puffsat import eos_methane
from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import near_term as nt

OUTPUT = nt.OUTPUT_DIR / "blowdown" / "cycle.csv"
SEAL_PRESSURES = (math.nan, 0.3e5, 0.1e5)
"""NaN stands for the chamber's own fill pressure (the parent's gate)."""
KAPTON = (69e6, 1420.0)
"""Kapton HN yield [Pa] and density [kg/m^3] (the parent's `dupont_kapton_hn_datasheet`)."""
ARAMID = (250e6, 1440.0)
"""Woven aramid at the parent's working stress [Pa], density [kg/m^3]; plus a 12.7 µm Kapton
liner."""
KAPTON_LINER = 12.7e-6
FILM_LAY = 0.030
"""Assumed: laying both films once the throat flow is weak [s]."""
LOCK = 0.020
"""The port insert locks this long before the rod (parent) [s]."""
INJECTOR = 1000.0
"""Assumed fill flow per 5 kg of rod [kg/s]: 108 kg in ~0.11 s. The parent sizes 800 kg/s."""
LEAK_C = 0.67
"""Choked-leak coefficient of a cold gas, `mdot = C A p / sqrt(R T)` (parent eq:throat_leak)."""


def _pure(feed_name: str, rho: float) -> eos_methane.Composition:
    """The cold feed as one molecule: CH4 or H2 at density `rho`."""
    mol = eos_methane.H2 if feed_name.startswith("hydrogen") else eos_methane.CH4
    mass = 2.0 * eos_methane.M_H if mol is eos_methane.H2 else eos_methane.M_CH4
    n_mol = tuple(rho / mass if m is mol else 0.0 for m in eos_methane.MOLECULES)
    return eos_methane.Composition(0.0, 0.0, 0.0, (), 0.0, n_mol)


def fill_energy(feed_name: str, rho: float, temp: float) -> float:
    """Thermal energy per kg [J/kg] of the cold feed gas at `temp` (zero at 0 K, bound)."""
    return eos_methane.thermal_energy(_pure(feed_name, rho), temp, "janaf") / rho


def temperature_for(feed_name: str, rho: float, e: float) -> float:
    """Invert `fill_energy` by bisection [K]."""
    lo, hi = 20.0, 6000.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if fill_energy(feed_name, rho, mid) < e:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class Seal:
    chamber: str
    seal_bar: float
    t_seal_ms: float
    residue_kg: float
    residue_temp: float
    #: Residue energy [MJ]: thermal part, and all of it (thermal + chemical, cold feed = 0).
    residue_thermal_mj: float
    residue_total_mj: float
    #: The mixed next charge: temperature [K] and pressure [bar], low and high bound.
    mix_temp_low: float
    mix_temp_high: float
    mix_bar_low: float
    mix_bar_high: float
    #: Throat film per shot holding the high-bound mix [kg]: Kapton alone; aramid + Kapton.
    film_kapton_kg: float
    film_aramid_kg: float
    #: Fill time at `INJECTOR` and the period it leaves [ms]: seal + films + fill + lock.
    fill_ms: float
    period_ms: float
    max_rate_hz: float
    #: Charge lost per 10 ms that the fill overlaps an open throat, at the mixed pressure [%].
    open_leak_pct_per_10ms: float


def seals(r: bd.Result, report: nt.Report = print) -> list[Seal]:
    spec, sol = r.spec, r.sol
    rows = r.runs["none"]
    br = nt.branch("equilibrium", sol)
    hydrogen = spec.pairing.startswith("hydrogen")
    gas_r = bd.R_H2 if hydrogen else bd.R_CH4
    m0 = sol.rho * spec.volume
    charge = sol.charge_mass * spec.energy_ratio
    radius = math.sqrt(spec.throat_area / math.pi)
    cap = 2.0 * math.pi * radius**2
    out = []
    for target in SEAL_PRESSURES:
        p_seal = r.p_fill if math.isnan(target) else target
        t = bd.time_to(rows, lambda x: x.pressure, p_seal)
        if not math.isfinite(t):
            continue
        row = bd.state_at(rows, t)
        m_res = row.mass_fraction * m0
        rho_res = m_res / spec.volume
        comp, _ = br.state(rho_res, row.temp)
        e_th = eos_methane.thermal_energy(comp, row.temp, sol.feed.thermo) / rho_res
        mass = charge + m_res
        rho_mix = mass / spec.volume
        e_fill = fill_energy(spec.pairing, rho_mix, spec.fill_temp)
        temps = [
            temperature_for(spec.pairing, rho_mix, (charge * e_fill + m_res * e_res) / mass)
            for e_res in (e_th, row.energy)
        ]
        bars = [mass * gas_r * tm / spec.volume / 1e5 for tm in temps]
        p_hold = bars[1] * 1e5
        t_kapton = p_hold * radius / (2.0 * KAPTON[0])
        t_aramid = p_hold * radius / (2.0 * ARAMID[0])
        fill = charge / (INJECTOR * spec.rod / 5.0)
        period = t + FILM_LAY + fill + LOCK
        leak = LEAK_C * spec.throat_area * p_hold / math.sqrt(gas_r * temps[1]) * 0.010
        s = Seal(
            chamber=spec.name,
            seal_bar=p_seal / 1e5,
            t_seal_ms=1e3 * t,
            residue_kg=m_res,
            residue_temp=row.temp,
            residue_thermal_mj=m_res * e_th / 1e6,
            residue_total_mj=m_res * row.energy / 1e6,
            mix_temp_low=temps[0],
            mix_temp_high=temps[1],
            mix_bar_low=bars[0],
            mix_bar_high=bars[1],
            film_kapton_kg=KAPTON[1] * t_kapton * cap,
            film_aramid_kg=(ARAMID[1] * t_aramid + KAPTON[1] * KAPTON_LINER) * cap,
            fill_ms=1e3 * fill,
            period_ms=1e3 * period,
            max_rate_hz=1.0 / period,
            open_leak_pct_per_10ms=100.0 * leak / charge,
        )
        out.append(s)
        report(
            f"{spec.name:30s} seal at {s.seal_bar:5.2f} bar: {s.t_seal_ms:4.0f} ms, residue "
            f"{s.residue_kg:5.2f} kg at {s.residue_temp:4.0f} K ({s.residue_total_mj:4.0f} MJ); "
            f"mix {s.mix_temp_low:4.0f}-{s.mix_temp_high:4.0f} K, {s.mix_bar_low:4.1f}-"
            f"{s.mix_bar_high:4.1f} bar; film {s.film_kapton_kg:5.2f} kg Kapton / "
            f"{s.film_aramid_kg:4.2f} aramid; period {s.period_ms:4.0f} ms -> {s.max_rate_hz:4.1f}"
            f" Hz; open-throat leak {s.open_leak_pct_per_10ms:.1f}%/10 ms"
        )
    return out


def hydrogen_throats() -> list[bd.Spec]:
    """B5: the bulged hydrogen chamber on wider throats, for its rate."""
    base = bd.hydrogen_specs()[1]
    return [
        dataclasses.replace(
            base,
            name=f"H2 106 m^3 bulged, A* {a:.2f}",
            throat_area=a,
            contour=bd.long_cone(40.0, a, bulge=True),
        )
        for a in (0.5, 0.75, 1.0)
    ]


def write(rows: Sequence[Seal]) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Seal.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for s in rows:
            w.writerow([getattr(s, f) for f in fields])


def main() -> None:
    rows: list[Seal] = []
    specs = bd.methane_specs() + bd.hydrogen_specs() + hydrogen_throats()
    for spec in specs:
        res = bd.run_spec(spec, area_ratios=(100.0, 300.0))
        for i in res.impulses:
            print(
                f"    AR {i.area_ratio:3.0f}: I {i.impulse / 1e3:.0f} kN s "
                f"(wall heat {i.impulse_low / 1e3:.0f}-{i.impulse_high / 1e3:.0f}), "
                f"exit {math.sqrt(4.0 * spec.throat_area * i.area_ratio / math.pi):.1f} m"
            )
        rows.extend(seals(res))
    write(rows)
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
