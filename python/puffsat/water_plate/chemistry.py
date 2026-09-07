"""Conditional recombination-location screen, not injected plate hydrodynamics.

Follow the collision-paid mixed states along de = -p d(1/rho), using the existing
equilibrium water EOS. Locate separate bond/ion store returns. A prescribed
uniform hemisphere, r = r0 + U*t, supplies an illustrative clock and size only;
neither its velocity nor its wall pressure follows from a momentum equation.
The actual flow must supply density histories and cumulative wall impulse.

Rate coefficients come from puffsat.recombination. The H and equilibrium-OH
partner clocks are sensitivity proxies, not rigorous bounds on a reacting
network. High-temperature extrapolation, steam third-body efficiencies and
pressure falloff remain unvalidated. Radiative capture is reported separately
from collisional capture: a photon is not automatically useful gas heat.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from itertools import pairwise
from pathlib import Path

from puffsat import eos_water, expansion, recombination
from puffsat.water_plate.ledger import _check
from puffsat.water_plate.thermodynamics import OUTPUT_DIR, MixedState, mixed_state

EXPANSION_RATIO = 1e5
STEPS = 256
EXPANSION_SPEED_M_S = 10_000.0  # Prescribed illustration; clocks scale as 1/U.


@dataclass(frozen=True)
class HemisphereClock:
    """Equal-volume radius and prescribed clock; not a resolved plume boundary."""

    radius_m: float
    time_s: float
    expansion_time_s: float


def hemisphere_clock(mass_kg: float, rho0: float, rho: float, speed_m_s: float) -> HemisphereClock:
    """M = (2*pi/3)*rho*r^3; r = r0 + U*t; |dln(rho)/dt|^-1 = r/(3U)."""
    for name, value in (("mass", mass_kg), ("rho0", rho0), ("rho", rho), ("speed", speed_m_s)):
        _check(name, value, positive=True)
    if rho > rho0 * (1 + 1e-12):
        raise ValueError("clock requires expansion, rho <= rho0")
    radius = (3 * mass_kg / (2 * math.pi * rho)) ** (1 / 3)
    radius0 = (3 * mass_kg / (2 * math.pi * rho0)) ** (1 / 3)
    return HemisphereClock(
        radius, max(0.0, (radius - radius0) / speed_m_s), radius / (3 * speed_m_s)
    )


def _blend_state(
    lo: expansion.ExpansionState, hi: expansion.ExpansionState, weight: float
) -> expansion.ExpansionState:
    """Log-interpolate rho and T, then re-evaluate pressure and energy consistently."""
    rho = math.exp(math.log(lo.rho) + weight * math.log(hi.rho / lo.rho))
    temp = math.exp(math.log(lo.temp) + weight * math.log(hi.temp / lo.temp))
    pressure, energy = eos_water.pressure_energy(rho, temp)
    return expansion.ExpansionState(rho, temp, pressure, energy)


def release_crossing(
    states: Sequence[expansion.ExpansionState], stores: Sequence[float], returned: float
) -> expansion.ExpansionState | None:
    """First crossing of a fraction returned from the INITIAL store; None if absent.

    The sampled store defines the crossing; the interpolated EOS store can differ
    slightly. Output includes that actual value and a resolution comparison.
    """
    if len(states) != len(stores) or len(states) < 2 or stores[0] <= 0:
        raise ValueError("need matching histories and a positive initial store")
    if not 0 < returned < 1:
        raise ValueError("returned fraction must be between zero and one")
    target = stores[0] * (1 - returned)
    for i in range(1, len(states)):
        if stores[i] <= target < stores[i - 1]:
            weight = (stores[i - 1] - target) / (stores[i - 1] - stores[i])
            return _blend_state(states[i - 1], states[i], weight)
    return None


def _ion_store(rho: float, comp: eos_water.Composition) -> float:
    return (
        comp.n_hp * eos_water.IP_H
        + sum(n * e for n, e in zip(comp.n_o_ions, eos_water.E_O_CUM, strict=True))
    ) / rho


def _inverse(rate: float) -> float:
    return 1 / rate if rate > 0 else math.inf


@dataclass(frozen=True)
class Diagnostic:
    """Stores and reaction clocks evaluated on an assumed equilibrium history."""

    density_kg_m3: float
    temperature_k: float
    pressure_pa: float
    volume_ratio: float
    bond_energy_j_kg: float
    ion_energy_j_kg: float
    bond_returned_fraction: float
    ion_returned_fraction: float
    internal_work_j_kg: float
    hemisphere_radius_m: float
    hemisphere_time_s: float
    hemisphere_expansion_time_s: float
    tau_ion_three_body_s: float
    tau_ion_radiative_s: float
    tau_bond_h_proxy_s: float
    tau_bond_equilibrium_oh_s: float
    da_ion_three_body: float
    da_bond_h_proxy: float
    da_bond_equilibrium_oh: float


def diagnose(state: expansion.ExpansionState, mixed: MixedState) -> Diagnostic:
    """Evaluate gas energy and clocks without assigning any of them to wall impulse."""
    rho, temp = state.rho, state.temp
    comp = eos_water.composition(rho, temp)
    bond = eos_water.bond_energy_held(rho, temp)
    ion = _ion_store(rho, comp)
    clock = hemisphere_clock(
        mixed.projectile_kg * (1 + mixed.k), mixed.mixed_density_kg_m3, rho, EXPANSION_SPEED_M_S
    )
    tau_ion = _inverse(recombination.three_body_coefficient(temp) * comp.n_e**2)
    tau_rad = _inverse(recombination.radiative_recombination(temp) * comp.n_e)
    tau_h = recombination.atom_recombination_time(temp, comp.n_neutral_heavy, comp.n_h)
    tau_oh = recombination.atom_recombination_time(temp, comp.n_neutral_heavy, comp.n_oh)
    return Diagnostic(
        rho,
        temp,
        state.pressure,
        mixed.mixed_density_kg_m3 / rho,
        bond,
        ion,
        1 - bond / mixed.bond_energy_j_kg,
        1 - ion / mixed.ionization_energy_j_kg,
        mixed.mixed_internal_energy_j_kg - state.energy,
        clock.radius_m,
        clock.time_s,
        clock.expansion_time_s,
        tau_ion,
        tau_rad,
        tau_h,
        tau_oh,
        clock.expansion_time_s / tau_ion,
        clock.expansion_time_s / tau_h,
        clock.expansion_time_s / tau_oh,
    )


def _at_density(
    states: Sequence[expansion.ExpansionState], target: float
) -> expansion.ExpansionState | None:
    for lo, hi in pairwise(states):
        if hi.rho <= target <= lo.rho:
            return _blend_state(lo, hi, math.log(target / lo.rho) / math.log(hi.rho / lo.rho))
    return None


def trajectory(mixed: MixedState, steps: int = STEPS) -> list[expansion.ExpansionState]:
    """Loss-free equilibrium adiabat; no magnetic-nozzle geometry/constants imported."""
    return expansion.expand(
        mixed.mixed_density_kg_m3,
        mixed.temperature_k,
        mixed.mixed_density_kg_m3 / EXPANSION_RATIO,
        eos_water.pressure_energy,
        steps=steps,
    )


CsvValue = str | int | float


def _write(path: Path, rows: list[dict[str, CsvValue]]) -> None:
    # Missing milestones have explicit status and blank state fields, never zeros.
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _milestones(
    mixed: MixedState, states: list[expansion.ExpansionState], steps: int
) -> list[dict[str, CsvValue]]:
    bond = [eos_water.bond_energy_held(s.rho, s.temp) for s in states]
    ion = [_ion_store(s.rho, eos_water.composition(s.rho, s.temp)) for s in states]
    rows: list[dict[str, CsvValue]] = []
    for name, stores in (("ionization", ion), ("molecular_bonds", bond)):
        for fraction in (0.1, 0.5, 0.9):
            state = release_crossing(states, stores, fraction)
            row: dict[str, CsvValue] = {
                "case": mixed.case,
                "closing_speed_m_s": mixed.closing_speed_m_s,
                "k": mixed.k,
                "initial_density_kg_m3": mixed.mixed_density_kg_m3,
                "steps": steps,
                "expansion_ratio_end": EXPANSION_RATIO,
                "prescribed_expansion_speed_m_s": EXPANSION_SPEED_M_S,
                "store": name,
                "target_returned_fraction": fraction,
                "status": "crossed" if state is not None else "not_reached",
            }
            if state is not None:
                row.update(asdict(diagnose(state, mixed)))
            rows.append(row)
    return rows


def write_figure() -> None:
    """Plot the generated paths; the plotting library stays behind this typed boundary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, layout="constrained")
    for rho0, color in zip((0.1, 1.0, 10.0), ("#0072B2", "#D55E00", "#009E73"), strict=True):
        with (OUTPUT_DIR / f"chemistry_history_rho_{rho0:g}.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        radius = [float(row["hemisphere_radius_m"]) for row in rows]
        for axis, field in zip(
            axes[:2], ("ion_returned_fraction", "bond_returned_fraction"), strict=True
        ):
            axis.plot(
                radius,
                [100 * float(row[field]) for row in rows],
                color=color,
                label=f"Initial density {rho0:g} kg/m³",
            )
        for field, style in (("da_bond_equilibrium_oh", "-"), ("da_bond_h_proxy", ":")):
            axes[2].plot(radius, [float(row[field]) for row in rows], style, color=color)
    for axis in axes:
        axis.set_xscale("log")
        axis.set_xlim(2, 500)
        axis.axvline(15, color="#666666", linestyle="--", linewidth=1)
        axis.grid(True, alpha=0.2)
    for axis in axes[:2]:
        axis.set_ylim(0, 100)
    axes[0].set_ylabel("Ionization store returned (%)")
    axes[0].legend(loc="lower right")
    axes[1].set_ylabel("Molecular bond store returned (%)")
    axes[1].text(16, 85, "15 m: maximum plate radius\n(size reference only)", color="#555555")
    axes[2].set_yscale("log")
    axes[2].set_ylim(1e-7, 1e5)
    axes[2].axhline(1, color="#555555", linewidth=1)
    axes[2].axhline(10, color="#555555", linewidth=1, linestyle="--")
    axes[2].set_ylabel("Molecular reaction Da = τ expansion / τ reaction")
    axes[2].text(
        30,
        1e3,
        "Solid: equilibrium OH partner\nDotted: H-density proxy\nDa ≳ 10: fast; Da ≲ 0.1: slow",
    )
    axes[2].set_xlabel("Equal-volume hemisphere radius (m) — assumed shape, not a flow solution")
    fig.suptitle(
        "Water plate chemistry screen: 45.58 km/s, k = 8.5\n"
        "Equilibrium adiabats; prescribed radial speed 10 km/s for reaction clocks"
    )
    fig.savefig(OUTPUT_DIR / "chemistry.png", dpi=160)
    plt.close(fig)


def main() -> None:
    """Screen the first overtake anchor at three densities; refine the central path."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    milestones: list[dict[str, CsvValue]] = []
    sizes: list[dict[str, CsvValue]] = []
    convergence: list[dict[str, CsvValue]] = []
    for rho0 in (0.1, 1.0, 10.0):
        mixed = mixed_state(45_580.0, 8.5, rho0)
        states = trajectory(mixed)
        rows = _milestones(mixed, states, STEPS)
        milestones.extend(rows)
        # Detailed paths remain generated scratch; compact cited tables are versioned.
        _write(
            OUTPUT_DIR / f"chemistry_history_rho_{rho0:g}.csv",
            [{"initial_density_kg_m3": rho0, **asdict(diagnose(s, mixed))} for s in states],
        )
        for radius in (5.0, 10.0, 15.0):
            target = mixed.projectile_kg * (1 + mixed.k) / ((2 * math.pi / 3) * radius**3)
            state = _at_density(states, target)
            row: dict[str, CsvValue] = {
                "case": mixed.case,
                "closing_speed_m_s": mixed.closing_speed_m_s,
                "k": mixed.k,
                "initial_density_kg_m3": rho0,
                "test_radius_m": radius,
                "steps": STEPS,
                "prescribed_expansion_speed_m_s": EXPANSION_SPEED_M_S,
                "status": "crossed" if state is not None else "outside_trajectory",
            }
            if state is not None:
                row.update(asdict(diagnose(state, mixed)))
            sizes.append(row)
        print(f"Finished rho0={rho0:g} kg/m3 ({STEPS} steps).", flush=True)
        if rho0 == 1.0:
            convergence.extend(rows)
            for steps in (STEPS // 2, 2 * STEPS):
                convergence.extend(_milestones(mixed, trajectory(mixed, steps), steps))
                print(f"Central trajectory refinement: {steps} steps.", flush=True)
    _write(OUTPUT_DIR / "chemistry_milestones.csv", milestones)
    _write(OUTPUT_DIR / "chemistry_sizes.csv", sizes)
    _write(OUTPUT_DIR / "chemistry_convergence.csv", convergence)
    write_figure()
    print("Equilibrium locations and prescribed clocks only; no recovered wall impulse.")


if __name__ == "__main__":
    main()
