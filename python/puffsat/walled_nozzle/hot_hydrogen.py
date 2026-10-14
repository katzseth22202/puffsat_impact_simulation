"""The bulged hydrogen chamber held hotter: the calcs repo's hot-hydrogen asks H1-H5 (2026-10-14).

`aim_is_all_you_need` (ADR 0063, `docs/hydrocode_asks_2026-10-14_hot_hydrogen.md`) asks whether
the bulged hydrogen chamber (106 m^3 at the 2.5 kg rod, 0.84 m^2 throat, 37.0 t dry Kevlar wrap,
today 4884 K on 59.9 kg) should run at 7000 K or 10 000 K on a smaller charge. A hotter charge
carries more energy per kilogram, so less of it holds the temperature; but the pulse energy is
fixed by the rod and momentum goes as `sqrt(2 m E)`, so each pulse gives less impulse unless the
chamber cycles faster. Every piece is an existing tool of this study, pointed at the hot charge:

- **H1, peak state.** `near_term.solve_charge` at the bulge's volume: the equilibrium charge
  (rod and plug included as polyethylene) that holds each temperature, on `eos_methane`'s
  partition functions with Saha ionisation. No fitted polynomial is extrapolated.
- **H2, recombination and impulse.** `finite_rate.pulse`: H + H + M at its measured rate through
  the AR 300 cone, every parcel of the blowdown, carbon bracketed by its two models (ADR-0057).
  Net impulse is B6's `m_rod ((1 - gate)(1 + k + P) u - w)`.
- **H3, the wall.** The gas is optically thick (`tops_optical_depth`, LANL TOPS gray Rosseland
  means for these H:C mixtures), so radiation to the wall is close to `sigma T^4` from the core,
  not the H- floor. The blowdown is integrated with that loss (`TopsHeat`) and the copper liner's
  surface rise is Duhamel's on the solved flux history. The wrap's fling comes from 2-D
  rod-and-plug solves with the lighter fills (`sweep --vessel-blast --hot-hydrogen`).
- **H4, the cycle.** `cycle.seals` on the blowdown: seal, films, fill, lock.
- **H5, the sweet spot.** Impulse per kilogram carried (charge with its tank, plug and film) and
  thrust at the H4 rate, from 4884 K to 10 000 K.
"""

from __future__ import annotations

import csv
import dataclasses
import math
from collections.abc import Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import cache
from itertools import pairwise
from pathlib import Path

import numpy as np

from puffsat import eos_methane, tops
from puffsat.walled_nozzle import blowdown as bd
from puffsat.walled_nozzle import cycle, hydrogen
from puffsat.walled_nozzle import finite_rate as fr
from puffsat.walled_nozzle import near_term as nt
from puffsat.walled_nozzle import zone_loads as zl

OUTPUT_DIR = nt.OUTPUT_DIR / "hot_hydrogen"

#: The bulged hydrogen chamber (B5, `walled_nozzle_hydrocode_asks.md` §6a), 2.5 kg rod.
VOLUME = 106.0
THROAT = 0.8378
AREA_RATIO = 300.0
CONTOUR = bd.long_cone(40.0, THROAT, bulge=True)
#: Today's charge: k = 24, which the bulge holds at 4884 K.
TODAY_CHARGE = bd.H2_CHARGE_25
#: The sweep for H5; the asks name 7000 and 10 000 K.
TEMPERATURES = (5500.0, 6000.0, 7000.0, 8000.0, 9000.0, 10000.0)
ASKED = (7000.0, 10000.0)
#: The ends of the burn and its reference [m/s] (the ask: w runs 74-81 km/s, k follows it).
SPEEDS = (74.0e3, 75.0e3, 81.0e3)
REFERENCE_SPEED = 75.0e3
#: Tank per kilogram of hydrogen (parent `chamber_isp`), and the burn's pulse count (the ask).
TANK_PER_KG = 0.205
BURN_PULSES = 4500

#: Copper liner: conductivity [W/m/K], density, heat capacity, melting point, and the hot-wall
#: temperature convection is quoted against. Ordinary OFHC copper; GRCop-84 conducts ~25% less.
COPPER = (390.0, 8960.0, 385.0)
COPPER_MELT = 1358.0
COPPER_WALL_TEMP = 800.0
#: Regenerative coolant's outlet temperature [K]: the copper's hot face caps it.
COOLANT_OUT = 700.0

#: The two TOPS pulls (H:C by number), for each charge's mixture.
TOPS_PULLS = (
    (50.6, Path("data/tables/tops/tops_hydrogen_7000_gray.html")),
    (39.7, Path("data/tables/tops/tops_hydrogen_10000_gray.html")),
)

#: The 2-D rod-and-plug runs of the hotter fills (`sweep --vessel-blast --hot-hydrogen`).
HOT_2D = nt.OUTPUT_DIR / "vessel_blast_hot_hydrogen.jsonl"


# ---- H1: the chamber ------------------------------------------------------------------------


def _pairing(temp: float) -> nt.Pairing:
    base = bd.pairing_named("hydrogen_5500")
    return dataclasses.replace(base, name=f"hydrogen_{temp:.0f}", temp=temp, paper_charge=25.0)


def today_spec() -> bd.Spec:
    """Today's bulged chamber: 59.9 kg, its temperature (4884 K) following from the charge."""
    return bd.Spec(
        "H2 bulged, today (k = 24)",
        "hydrogen_5500",
        nt.ROD_MASS,
        VOLUME,
        THROAT,
        CONTOUR,
        TODAY_CHARGE,
        bd.FILL_TEMP_H2,
    )


@cache
def chamber(temp: float | None, speed: float = REFERENCE_SPEED) -> nt.ChamberSolution:
    """The peak state: `temp` held by its charge at closing speed `speed`; `None` is today's."""
    if temp is None:
        if speed != REFERENCE_SPEED:
            raise ValueError("today's chamber is quoted at the reference speed only")
        return bd.solve(today_spec())
    return nt.solve_charge(_pairing(temp), "janaf", VOLUME, energy=0.5 * nt.ROD_MASS * speed**2)


def spec_for(sol: nt.ChamberSolution) -> bd.Spec:
    """The blowdown spec of a solved chamber (its fixed charge reproduces its temperature)."""
    return dataclasses.replace(
        today_spec(), name=f"H2 bulged {sol.pairing.temp:.0f} K", charge_25=sol.charge_mass
    )


@cache
def _tops_tables() -> tuple[tuple[float, tops.TopsGray], ...]:
    return tuple((hc, tops.parse_tops_gray(path.read_text())) for hc, path in TOPS_PULLS)


def hc_ratio(sol: nt.ChamberSolution) -> float:
    """Hydrogen atoms per carbon atom in the chamber (charge plus polyethylene)."""
    poly = nt.ROD_MASS + nt.PLUG_MASS
    m_ch2 = eos_methane.M_C + 2.0 * eos_methane.M_H
    n_c = poly / m_ch2
    n_h = sol.charge_mass / eos_methane.M_H + 2.0 * n_c
    return n_h / n_c


def tops_rosseland(hc: float, rho: float, temp: float) -> float:
    """TOPS gray Rosseland mean [m^2/kg], log-log bilinear, on the pull nearest `hc`; NaN below
    the OPLIB floor (5802 K), where TOPS has no data (and no molecules)."""
    _, g = min(_tops_tables(), key=lambda t: abs(math.log(t[0] / hc)))
    if temp < float(g.t_grid[0]):
        return math.nan
    lr, lt = math.log(rho), math.log(temp)
    xr, xt = np.log(g.rho_grid), np.log(g.t_grid)
    i = int(np.clip(np.searchsorted(xr, lr) - 1, 0, len(xr) - 2))
    j = int(np.clip(np.searchsorted(xt, lt) - 1, 0, len(xt) - 2))
    fr_ = (lr - xr[i]) / (xr[i + 1] - xr[i])
    ft = (lt - xt[j]) / (xt[j + 1] - xt[j])
    k = np.log(g.kappa_rosseland[i : i + 2, j : j + 2])
    val = (
        (1 - fr_) * (1 - ft) * k[0, 0]
        + fr_ * (1 - ft) * k[1, 0]
        + (1 - fr_) * ft * k[0, 1]
        + fr_ * ft * k[1, 1]
    )
    return float(math.exp(val))


def beam_length() -> float:
    """Mean beam length of the bulge, `4 V / A` [m]."""
    return 4.0 * VOLUME / CONTOUR.area()


@dataclass(frozen=True)
class PeakState:
    """H1's row."""

    label: str
    speed: float
    temp: float
    charge: float
    slug_ratio: float
    rho: float
    pressure_bar: float
    #: Energy per kilogram of everything in the chamber [MJ/kg], and its share in broken bonds.
    energy_mj: float
    bond_share: float
    hydrogen_dissociated: float
    carbon_free: float
    ionisation: float
    #: Pure hydrogen at the parent's 0.57 kg/m^3, for its Cantera column [MJ/kg above H2 at 0 K].
    pure_h2_energy_mj: float
    gamma_eff: float
    #: Electron-ion coupling `e^2 / (4 pi eps0 a kT)`: the ideal-gas EOS needs it << 1.
    coupling: float
    hc_ratio: float
    #: Rosseland optical depth across the beam length (TOPS); NaN below the OPLIB floor.
    tau_rosseland: float
    #: In-band H- bound-free floor (`near_term.optical_depth_floor`), for comparison.
    tau_h_minus: float


def peak_state(temp: float | None, speed: float = REFERENCE_SPEED) -> PeakState:
    sol = chamber(temp, speed)
    t = sol.pairing.temp
    comp = eos_methane.composition(sol.rho, t, sol.feed)
    n_e = max(comp.n_e, 1.0)
    a = (3.0 / (4.0 * math.pi * n_e)) ** (1.0 / 3.0)
    coupling = 1.602176634e-19**2 / (4.0 * math.pi * 8.8541878128e-12 * a * eos_methane.K_B * t)
    hc = hc_ratio(sol)
    kappa = tops_rosseland(hc, sol.rho, t)
    return PeakState(
        label="today (k = 24)" if temp is None else f"{t:.0f} K",
        speed=speed,
        temp=t,
        charge=sol.charge_mass,
        slug_ratio=sol.slug_ratio,
        rho=sol.rho,
        pressure_bar=sol.pressure / 1e5,
        energy_mj=sol.energy / 1e6,
        bond_share=sol.bond_share,
        hydrogen_dissociated=sol.hydrogen_dissociated,
        carbon_free=sol.carbon_free,
        ionisation=sol.ionisation_fraction,
        pure_h2_energy_mj=hydrogen.state(0.57, t).energy / 1e6,
        gamma_eff=1.0 + sol.pressure / (sol.rho * sol.energy),
        coupling=coupling,
        hc_ratio=hc,
        tau_rosseland=kappa * sol.rho * beam_length(),
        tau_h_minus=nt.optical_depth_floor(comp, t, beam_length()),
    )


# ---- H2: finite-rate recombination and the impulse ------------------------------------------


@dataclass(frozen=True)
class PulseJob:
    temp: float | None
    speed: float
    model: str
    rate_scale: float = 1.0


def pulse_jobs() -> list[PulseJob]:
    out = [PulseJob(t, REFERENCE_SPEED, m) for t in (None, *TEMPERATURES) for m in fr.MODELS]
    for t in ASKED:
        for w in (SPEEDS[0], SPEEDS[2]):
            out += [PulseJob(t, w, m) for m in fr.MODELS]
        out += [
            PulseJob(t, REFERENCE_SPEED, "exchange_carbon", s)
            for s in (fr.RATE_SCALES[0], fr.RATE_SCALES[2])
        ]
    return out


@cache
def _model_and_tab(temp: float | None, speed: float, model: str) -> tuple[fr.Model, nt.Isentrope]:
    sol = chamber(temp, speed)
    mod = fr.model(model, sol)
    tab = nt.isentrope(mod.branch, sol.rho, sol.pairing.temp, rho_ratio_end=1e-8, steps=900)
    return mod, tab


def run_pulse(job: PulseJob) -> tuple[PulseJob, fr.Pulse]:
    mod, tab = _model_and_tab(job.temp, job.speed, job.model)
    return job, fr.pulse(mod, tab, THROAT, AREA_RATIO, rate_scale=job.rate_scale)


def net_impulse(sol: nt.ChamberSolution, speed: float, closing: float) -> float:
    """B6's net impulse per pulse [N s]: gate charged, head-on debit `m_rod w` taken."""
    expelled = 1.0 + sol.slug_ratio + nt.PLUG_RATIO
    return nt.ROD_MASS * ((1.0 - bd.GATE_THRUST_COST) * expelled * speed - closing)


# ---- H3: the wall ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TopsHeat:
    """`q(rho, T, p)` [W/m^2]: optically thick radiation from the core plus Bartz convection.

    Radiation is `sigma T^4 (1 - exp(-tau))`, `tau` the TOPS Rosseland depth across the beam
    length. Below TOPS's 5802 K floor the in-band H- floor is used instead, the low edge's own
    expression; by then the chamber has emptied most of its energy. Convection is
    `blowdown.HeatModel`'s `p^0.8 (T - T_w)` scaling from its peak, on the copper wall.
    """

    conv_peak: float
    p0: float
    t0: float
    hc: float
    beam: float
    wall_temp: float = COPPER_WALL_TEMP

    def radiation(self, br: nt.Branch, rho: float, temp: float) -> float:
        kappa = tops_rosseland(self.hc, rho, temp)
        if math.isnan(kappa):
            comp, _ = br.state(rho, temp)
            tau = nt.optical_depth_floor(comp, temp, self.beam)
            return float(bd.SIGMA_SB * temp**4 * nt.band_fraction(temp) * -math.expm1(-tau))
        return float(bd.SIGMA_SB * temp**4 * -math.expm1(-kappa * rho * self.beam))

    def convection(self, temp: float, p: float) -> float:
        return float(
            self.conv_peak
            * (p / self.p0) ** 0.8
            * max(0.0, temp - self.wall_temp)
            / (self.t0 - self.wall_temp)
        )

    def __call__(self, br: nt.Branch, rho: float, temp: float, p: float) -> float:
        return self.radiation(br, rho, temp) + self.convection(temp, p)


@dataclass(frozen=True)
class Wall:
    """H3's row: per pulse, then at the H4 rate over the burn."""

    label: str
    temp: float
    #: Peak fluxes [MW/m^2]: radiation (TOPS-thick), liner convection low-high, throat
    #: convection low-high.
    q_rad_peak: float
    q_conv_low: float
    q_conv_high: float
    q_throat_low: float
    q_throat_high: float
    #: Per pulse [MJ/m^2] on the liner: radiation and convection (low and high edge).
    fluence_rad: float
    fluence_conv_low: float
    fluence_conv_high: float
    #: Share of the pulse energy the wall takes (high convection edge).
    wall_share: float
    #: Wall heat per kilogram of charge [MJ/kg] against what the charge can take as regen
    #: coolant, 77 K -> `COOLANT_OUT`.
    heat_per_kg_charge: float
    coolant_capacity: float
    #: Copper surface rise over one pulse [K] (Duhamel, semi-infinite, from its base).
    copper_rise: float
    #: Share of the exhaust momentum the wall heat removes (`blowdown.heat_debit`).
    impulse_debit: float


def copper_rise(times: Sequence[float], fluxes: Sequence[float]) -> float:
    """Peak surface rise [K] of a semi-infinite copper face under `q(t)` (Duhamel).

    `dT_s(t) = 2 / (e sqrt(pi)) sum_i q_i (sqrt(t - t_i) - sqrt(t - t_{i+1}))` for `q` held
    piecewise constant on `[t_i, t_{i+1}]`, `e = sqrt(k rho c)` the effusivity.
    """
    t = np.asarray(times)
    q = 0.5 * (np.asarray(fluxes)[1:] + np.asarray(fluxes)[:-1])
    eff = math.sqrt(COPPER[0] * COPPER[1] * COPPER[2])
    best = 0.0
    for n in range(1, len(t)):
        a = np.sqrt(np.maximum(t[n] - t[:n], 0.0))
        b = np.sqrt(np.maximum(t[n] - t[1 : n + 1], 0.0))
        best = max(best, float(np.sum(q[:n] * (a - b))) * 2.0 / (eff * math.sqrt(math.pi)))
    return best


def coolant_capacity() -> float:
    """Enthalpy [J/kg] hydrogen takes as regen coolant from the 77 K fill to `COOLANT_OUT`."""
    rho = 1.0
    de = cycle.fill_energy("hydrogen", rho, COOLANT_OUT) - cycle.fill_energy(
        "hydrogen", rho, bd.FILL_TEMP_H2
    )
    return de + bd.R_H2 * (COOLANT_OUT - bd.FILL_TEMP_H2)


@dataclass(frozen=True)
class Blowdown:
    """One chamber's blowdown, adiabatic and with the thick-gas wall heat, and its cycle."""

    label: str
    temp: float
    result: bd.Result
    lossy: bd.Result
    wall: Wall
    seals: list[cycle.Seal]
    seals_lossy: list[cycle.Seal]


def _result(spec: bd.Spec, sol: nt.ChamberSolution, rows: list[bd.Row]) -> bd.Result:
    """A `blowdown.Result` carrying one run, as `cycle.seals` reads it."""
    return bd.Result(
        spec, sol, spec.contour.area(), bd.fill_pressure(spec, sol), math.nan,
        (math.nan, math.nan), {"none": rows}, [], math.nan,
    )  # fmt: skip


def run_blowdown(temp: float | None) -> Blowdown:
    sol = chamber(temp)
    spec = today_spec() if temp is None else spec_for(sol)
    br = nt.branch("equilibrium", sol)
    tab = nt.isentrope(br, sol.rho, sol.pairing.temp, rho_ratio_end=1.0e-8, steps=900)
    cal = bd.calibration(br, tab)
    none = bd.integrate(spec, sol, br, bd.HeatModel("none", 0.0, 1.0, 1.0, 1.0), cal)
    lo, hi = nt.bartz_liner(sol, tab, THROAT, COPPER_WALL_TEMP)
    th_lo, th_hi = nt.bartz_liner(sol, tab, THROAT, COPPER_WALL_TEMP, area_ratio=1.0)
    hc = hc_ratio(sol)
    runs = {}
    for edge, conv in (("low", lo), ("high", hi)):
        heat = TopsHeat(conv, sol.pressure, sol.pairing.temp, hc, beam_length())
        runs[edge] = bd.integrate(spec, sol, br, heat, cal)
    heat_hi = TopsHeat(hi, sol.pressure, sol.pairing.temp, hc, beam_length())

    def fluence(rows: list[bd.Row], part: str, heat: TopsHeat) -> float:
        total = 0.0
        for a, b in pairwise(rows):
            qa = (
                heat.radiation(br, a.mass_fraction * sol.rho, a.temp)
                if part == "rad"
                else heat.convection(a.temp, a.pressure)
            )
            qb = (
                heat.radiation(br, b.mass_fraction * sol.rho, b.temp)
                if part == "rad"
                else heat.convection(b.temp, b.pressure)
            )
            total += 0.5 * (qa + qb) * (b.time - a.time)
        return total

    heat_lo = dataclasses.replace(heat_hi, conv_peak=lo)
    speeds = bd.exhaust_speeds(tab, AREA_RATIO)
    hi_rows = runs["high"]
    area = spec.contour.area()
    per_pulse = hi_rows[-1].wall_share * nt.PULSE_ENERGY
    wall = Wall(
        label=spec.name,
        temp=sol.pairing.temp,
        q_rad_peak=heat_hi.radiation(br, sol.rho, sol.pairing.temp) / 1e6,
        q_conv_low=lo / 1e6,
        q_conv_high=hi / 1e6,
        q_throat_low=th_lo / 1e6,
        q_throat_high=th_hi / 1e6,
        fluence_rad=fluence(hi_rows, "rad", heat_hi) / 1e6,
        fluence_conv_low=fluence(runs["low"], "conv", heat_lo) / 1e6,
        fluence_conv_high=fluence(hi_rows, "conv", heat_hi) / 1e6,
        wall_share=hi_rows[-1].wall_share,
        heat_per_kg_charge=per_pulse / sol.charge_mass / 1e6,
        coolant_capacity=coolant_capacity() / 1e6,
        copper_rise=copper_rise([r.time for r in hi_rows], [r.q_wall for r in hi_rows]),
        impulse_debit=bd.heat_debit(hi_rows, speeds),
    )
    del area
    res = _result(spec, sol, none)
    lossy = _result(spec, sol, hi_rows)
    quiet: nt.Report = lambda _: None  # noqa: E731
    return Blowdown(
        spec.name, sol.pairing.temp, res, lossy, wall, cycle.seals(res, quiet),
        cycle.seals(lossy, quiet),
    )  # fmt: skip


# ---- H4 + H5: the rate, and impulse per kilogram carried ------------------------------------

#: Kapton's working rating (parent ADR 0053), the film's limit on the mixed charge [K].
FILM_RATING = 673.0
#: Fill rates to price [kg/s at the 2.5 kg rod]: `cycle.INJECTOR` (1 t/s per 5 kg) and twice it.
INJECTORS = (cycle.INJECTOR / 2.0, cycle.INJECTOR)


def period(s: cycle.Seal, charge: float, injector: float) -> float:
    """Seal + films + fill at `injector` + lock [s]."""
    return 1e-3 * s.t_seal_ms + cycle.FILM_LAY + charge / injector + cycle.LOCK


def best_seal(seals: Sequence[cycle.Seal]) -> cycle.Seal:
    """The earliest seal whose mixed charge (high bound) the polymer film can hold."""
    ok = [s for s in seals if s.mix_temp_high <= FILM_RATING]
    return min(ok or seals, key=lambda s: s.t_seal_ms)


@dataclass(frozen=True)
class SweetSpot:
    label: str
    temp: float
    charge: float
    impulse: float
    impulse_wall: float
    rate_hz: float
    rate_fast_hz: float
    film_kg: float
    carried_kg: float
    #: Impulse per kilogram carried [m/s]: the departure's effective exhaust speed.
    per_kg: float
    #: Mean thrust [kN] at `rate_hz` and at the faster fill.
    thrust: float
    thrust_fast: float


# ---- Output ---------------------------------------------------------------------------------


def _write(path: Path, rows: Sequence[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(type(rows[0]).__dataclass_fields__)  # type: ignore[attr-defined]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in rows:
            w.writerow(
                [f"{v:.6g}" if isinstance(v, float) else v for v in (getattr(r, f) for f in fields)]
            )


def hot_wall_cases() -> list[zl.Case]:
    """The 2-D runs of the hotter fills, read at their own fixed-charge QSP."""
    out = []
    for temp in ASKED:
        sol = chamber(temp)
        out.append(
            zl.Case(
                f"H2 106 m^3 bulged, {temp:.0f} K",
                HOT_2D,
                f"hydrogen {temp:.0f} K, bulge r 3.0 m z 1.4-2.8, 3 m taper, plug 1.5 m",
                VOLUME,
                spec_for(sol),
                energy=1.0,
            )
        )
    return out


def main() -> None:
    temps: list[float | None] = [None, *TEMPERATURES]
    with ProcessPoolExecutor() as pool:
        pulses_f = pool.map(run_pulse, sorted(pulse_jobs(), key=lambda j: j.model))
        blow_f = pool.map(run_blowdown, temps)
        pulses = list(pulses_f)
        blows = list(blow_f)

    peaks = [peak_state(t) for t in temps] + [
        peak_state(t, w) for t in ASKED for w in (SPEEDS[0], SPEEDS[2])
    ]
    print("## H1: peak state (106 m^3, 2.5 kg rod)")
    for p in peaks:
        print(
            f"{p.label:15s} w {p.speed / 1e3:4.0f} | charge {p.charge:5.1f} kg "
            f"k {p.slug_ratio:5.2f}"
            f" | rho {p.rho:.3f} p {p.pressure_bar:5.0f} bar | u {p.energy_mj:5.0f} MJ/kg "
            f"bonds {p.bond_share:.2f} | H free {p.hydrogen_dissociated:.2f} C free "
            f"{p.carbon_free:.2f} ion {p.ionisation:.1e} | pure H2 @0.57 "
            f"{p.pure_h2_energy_mj:4.0f} "
            f"| G_eff {p.gamma_eff:.3f} coupling {p.coupling:.1e} | H:C {p.hc_ratio:.0f} "
            f"tau_R {p.tau_rosseland:.1f} tau_H- {p.tau_h_minus:.2f}"
        )

    by = {(j.temp, j.speed, j.model, round(j.rate_scale, 4)): p for j, p in pulses}
    pulse_rows = []
    print("\n## H2: finite-rate AR 300 (0.84 m^2, 15 deg), whole blowdown")
    for (t, w, m, scale), pu in sorted(
        by.items(), key=lambda kv: (kv[0][0] or 0, kv[0][1], kv[0][2])
    ):
        sol = chamber(t, w)
        v = pu.speed_shifting * pu.thrust_ratio_finite
        imp = net_impulse(sol, v, w)
        imp_eq = net_impulse(sol, pu.speed_shifting, w)
        pulse_rows.append(
            {
                "temp": sol.pairing.temp, "speed": w, "model": m, "rate_scale": scale,
                "charge": sol.charge_mass, "eta_shifting": pu.eta_shifting,
                "eta_finite": pu.eta_finite, "eta_frozen": pu.eta_frozen,
                "bond_shifting": pu.bond_shifting, "bond_finite": pu.bond_finite,
                "bond_frozen": pu.bond_frozen, "thrust_ratio_finite": pu.thrust_ratio_finite,
                "speed_effective": v, "impulse_eq": imp_eq, "impulse": imp,
                "peak_freeze_area_ratio": pu.peak.freeze_area_ratio,
                "peak_throat_damkohler": pu.peak.throat_damkohler,
            }
        )  # fmt: skip
        print(
            f"{sol.pairing.temp:6.0f} K w {w / 1e3:3.0f} {m:15s} x{scale:.3g} | charge "
            f"{sol.charge_mass:5.1f} | eta eq {pu.eta_shifting:.3f} fin {pu.eta_finite:.3f} frz "
            f"{pu.eta_frozen:.3f} | bond eq {pu.bond_shifting:.2f} fin {pu.bond_finite:.2f} frz "
            f"{pu.bond_frozen:.2f} | thrust {pu.thrust_ratio_finite:.4f} | v {v:5.0f} | I eq "
            f"{imp_eq / 1e3:4.0f} fin {imp / 1e3:4.0f} kN s | freeze AR "
            f"{pu.peak.freeze_area_ratio:.0f} Da* {pu.peak.throat_damkohler:.3g}"
        )

    print("\n## H3: the wall (copper at 800 K; radiation TOPS-thick)")
    for b in blows:
        w_ = b.wall
        print(
            f"{w_.label:28s} q_rad {w_.q_rad_peak:5.0f} conv {w_.q_conv_low:4.0f}-"
            f"{w_.q_conv_high:4.0f} throat {w_.q_throat_low:4.0f}-{w_.q_throat_high:4.0f} MW/m2 | "
            f"per pulse rad {w_.fluence_rad:5.1f} conv {w_.fluence_conv_low:4.1f}-"
            f"{w_.fluence_conv_high:4.1f} MJ/m2 | share {w_.wall_share:.3f} | "
            f"{w_.heat_per_kg_charge:5.1f} MJ/kg charge vs coolant {w_.coolant_capacity:.1f} | "
            f"Cu rise {w_.copper_rise:5.0f} K | I debit {w_.impulse_debit:.3f}"
        )

    print("\n## H4: the cycle (seal / films 30 ms / fill / lock 20 ms)")
    for b in blows:
        for tag, ss in (("adiabatic", b.seals), ("wall heat", b.seals_lossy)):
            for s in ss:
                charge = b.result.sol.charge_mass
                print(
                    f"{b.label:24s} {tag:9s} seal {s.seal_bar:5.2f} bar at {s.t_seal_ms:4.0f} ms"
                    f" | residue {s.residue_kg:4.2f} kg {s.residue_temp:5.0f} K | mix "
                    f"{s.mix_temp_low:4.0f}-{s.mix_temp_high:4.0f} K {s.mix_bar_low:4.2f}-"
                    f"{s.mix_bar_high:4.2f} bar | film {s.film_aramid_kg:.2f} kg | "
                    + " / ".join(f"{1.0 / period(s, charge, inj):.2f}" for inj in INJECTORS)
                    + " Hz (0.5 / 1 t/s)"
                )

    # H5: per temperature, the exchange-carbon finite-rate impulse with the wall debit, at the
    # rate the lossy cycle allows with the film's rating respected.
    spots = []
    for b, t in zip(blows, temps, strict=True):
        sol = b.result.sol
        pu = by[(t, REFERENCE_SPEED, "exchange_carbon", 1.0)]
        v = pu.speed_shifting * pu.thrust_ratio_finite
        imp = net_impulse(sol, v, REFERENCE_SPEED)
        imp_wall = net_impulse(sol, v * (1.0 - b.wall.impulse_debit), REFERENCE_SPEED)
        seal = best_seal(b.seals_lossy)
        rates = [1.0 / period(seal, sol.charge_mass, inj) for inj in INJECTORS]
        carried = sol.charge_mass * (1.0 + TANK_PER_KG) + nt.PLUG_MASS + seal.film_aramid_kg
        spots.append(
            SweetSpot(
                b.label, sol.pairing.temp, sol.charge_mass, imp, imp_wall, rates[0], rates[1],
                seal.film_aramid_kg, carried, imp_wall / carried, imp_wall * rates[0] / 1e3,
                imp_wall * rates[1] / 1e3,
            )
        )  # fmt: skip
    print("\n## H5: impulse per kg carried and thrust (exchange carbon, wall heat debited)")
    for sp in spots:
        print(
            f"{sp.label:28s} charge {sp.charge:5.1f} | I {sp.impulse / 1e3:4.0f} -> "
            f"{sp.impulse_wall / 1e3:4.0f} kN s with wall | {sp.rate_hz:.2f} / "
            f"{sp.rate_fast_hz:.2f} Hz | carried {sp.carried_kg:5.1f} kg -> {sp.per_kg:5.0f} m/s "
            f"| thrust {sp.thrust:4.0f} / {sp.thrust_fast:4.0f} kN"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _write(OUTPUT_DIR / "peak_states.csv", peaks)
    with (OUTPUT_DIR / "pulses.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(pulse_rows[0]))
        wr.writeheader()
        wr.writerows(
            {k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()}
            for r in pulse_rows
        )
    _write(OUTPUT_DIR / "wall.csv", [b.wall for b in blows])
    seal_rows = [
        dataclasses.replace(s, chamber=f"{s.chamber} ({tag})")
        for b in blows
        for tag, ss in (("adiabatic", b.seals), ("wall heat", b.seals_lossy))
        for s in ss
    ]
    _write(OUTPUT_DIR / "cycle.csv", seal_rows)
    _write(OUTPUT_DIR / "sweet_spot.csv", spots)

    wrap()


def wrap() -> None:
    """H3's dynamic half: the dry wrap, zone by zone, under today's and the hotter fills."""
    if not HOT_2D.exists():
        print(f"\n{HOT_2D} not solved yet: run `sweep --vessel-blast --hot-hydrogen`")
        return
    today = [c for c in zl.cases() if c.name == "H2 106 m^3 bulged (2.5 kg)"]
    rows: list[zl.ZoneLoad] = []
    print("\n## H3: the wrap under the hotter fills (2-D rod and plug, 0.5 cm)")
    for case in today + hot_wall_cases():
        got = zl.loads(case)
        rows.extend(dataclasses.replace(r, case=case.name) for r in got)
        print(
            f"{case.name}: QSP {got[0].qsp_bar:.0f} bar; dry wrap "
            f"{sum(r.tonnes_wrap for r in got):.1f} t "
            f"({'passes' if got[0].passes_wrap else 'FAILS'})"
        )
        for r in got:
            print(
                f"  {r.zone:14s} peak {r.peak_gpa:5.2f} GPa  I {r.spike_impulse_kpas:5.1f} "
                f"kPa s  D wrap {r.overshoot_wrap:5.2f}  fling {r.fling_wrap_pct:5.2f}%"
            )
    _write(OUTPUT_DIR / "wrap.csv", rows)


if __name__ == "__main__":
    main()
