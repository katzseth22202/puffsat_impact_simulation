"""LTE validity of the sweep's plasma model (Q5 of the 16-63 km/s extension).

Two independent pieces of the high-v model assume **local thermodynamic equilibrium**: the EOS is
an equilibrium Saha ladder (`eos_water`), and the opacity is LTE TOPS/OPLIB gray data. Both take
the level populations to be collisionally controlled. If the plasma is instead
collisional-radiative — radiative decay outrunning collisional (de)excitation — the ionization
balance is wrong, and *both* the EOS sink and the opacity err in the same direction rather than
bracketing each other.

The standard test is the **McWhirter criterion** (McWhirter 1965): LTE requires

    n_e  >=  1.6e12 * sqrt(T) * dE^3     [cm^-3, T in K, dE in eV]

where `dE` is the largest energy gap in the level structure — conventionally the resonance
transition, ground to first excited state. It is the density at which collisional excitation and
de-excitation run at least 10x faster than the competing radiative decay.

**It is a necessary, not a sufficient, condition**, and it is derived for a *homogeneous,
stationary* plasma. A bounce is neither: the slug is stratified and the state changes on the
hydrodynamic timescale. So a comfortable pass here does not prove LTE, it only fails to disprove
it; a failure, by contrast, is decisive. Any margin this module reports should be read that way.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_water

# McWhirter coefficient in the literature's CGS form [cm^-3 K^-1/2 eV^-3], converted to SI by the
# 1e6 cm^-3 -> m^-3 factor. Kept separate so the published number stays recognizable.
MCWHIRTER_COEFF_CGS = 1.6e12
CM3_PER_M3 = 1.0e6


def mcwhirter_threshold(temp_k: float, delta_e_ev: float) -> float:
    """Minimum electron density [m^-3] for LTE at `temp_k` given a largest gap `delta_e_ev`.

    The criterion's SI form: `1.6e18 * sqrt(T) * dE^3` m^-3.
    """
    return MCWHIRTER_COEFF_CGS * CM3_PER_M3 * math.sqrt(temp_k) * delta_e_ev**3


HC_EV_NM = 1239.842  # hc in eV*nm, so dE[eV] = HC_EV_NM / lambda[nm]


@dataclass(frozen=True)
class SpeciesGap:
    """A species' resonance transition: the strongest dipole-allowed line out of the ground term,
    which sets its McWhirter gap."""

    label: str  # spectroscopic notation, e.g. "O VI"
    line_nm: float  # vacuum wavelength [nm], NIST ASD
    transition: str

    @property
    def gap_ev(self) -> float:
        return HC_EV_NM / self.line_nm


# Resonance transitions keyed by the `eos_water.Composition` species label, retrieved from the
# **NIST Atomic Spectra Database** (2026-08-17) as the longest-wavelength *strong* line out of the
# ground term — see `todos/nist_resonance/` for the raw query responses. Gaps are derived from the
# wavelength rather than stored in eV so the provenance stays visible.
#
# Selecting these correctly needs care; three ways to get a plausible wrong number:
#   - a wavelength window is not a selector. O II has 4p-5s lines at 834 A, right on top of its
#     834.47 A resonance line but from levels ~30 eV above ground.
#   - the longest ground-connected allowed line is often a weak *intercombination* line (O V
#     1S-3P* at 121.8 nm, gA 7e3 vs the resonance line's 8.6e9) at a smaller gap, which would make
#     the criterion look easier to satisfy than it is.
#   - O II and O III resonance lines nearly coincide (83.45 / 83.53 nm) — the solar "834 A
#     multiplet" — so they are easy to conflate.
#
# `H+` is deliberately absent: a bare proton has no bound levels, hence no transition to hold in
# LTE, even though it dominates the high-v turnaround by number.
RESONANCE_GAP_EV: dict[str, SpeciesGap] = {
    "H": SpeciesGap("H I", 121.5670, "1s 2S - 2p 2P* (Lyman-alpha)"),
    "O": SpeciesGap("O I", 130.6029, "2s2.2p4 3P - 2s2.2p3(4S*).3s 3S*"),
    "O1+": SpeciesGap("O II", 83.4466, "2s2.2p3 4S* - 2s.2p4 4P"),
    "O2+": SpeciesGap("O III", 83.5292, "2s2.2p2 3P - 2s.2p3 3D*"),
    "O3+": SpeciesGap("O IV", 79.0199, "2s2.2p 2P* - 2s.2p2 2D"),
    "O4+": SpeciesGap("O V", 62.9730, "1s2.2s2 1S - 1s2.2s.2p 1P*"),
    "O5+": SpeciesGap("O VI", 103.7615, "1s2.2s 2S - 1s2.2p 2P* (lithium-like)"),
    "O6+": SpeciesGap("O VII", 2.1601, "1s2 1S - 1s.2p 1P* (helium-like, K-shell)"),
    "O7+": SpeciesGap("O VIII", 1.8973, "1s 2S - 2p 2P* (hydrogenic Lyman-alpha)"),
}

# Species below this fraction of the heavy-particle population are excluded from the verdict: a
# trace ion cannot be allowed to condemn a state it barely populates.
ABUNDANCE_FLOOR = 0.01


@dataclass(frozen=True)
class StateVerdict:
    """Where one `(T, n_e)` state sits relative to the McWhirter criterion."""

    temp_k: float
    n_e: float
    critical_gap_ev: float
    governing: str
    governing_gap_ev: float
    margin: float
    abundant: tuple[str, ...]
    failing: tuple[str, ...]

    @property
    def lte_valid(self) -> bool:
        """True when every abundant species' resonance gap is collisionally controlled here."""
        return not self.failing


def evaluate_state(temp_k: float, n_e: float, abundances: dict[str, float]) -> StateVerdict:
    """Apply the criterion at one state, given species fractions of the heavy-particle population.

    The governing species is the one with the **largest** resonance gap among those above
    `ABUNDANCE_FLOOR` — the gap enters cubed, so the most demanding species decides, not the most
    abundant one.
    """
    d_crit = critical_gap(temp_k, n_e)
    present = sorted(
        (s for s, f in abundances.items() if f >= ABUNDANCE_FLOOR and s in RESONANCE_GAP_EV),
        key=lambda s: RESONANCE_GAP_EV[s].gap_ev,
        reverse=True,
    )
    if not present:  # e.g. a fully-stripped state: no bound levels anywhere
        return StateVerdict(temp_k, n_e, d_crit, "-", 0.0, math.inf, (), ())
    governing = present[0]
    gap = RESONANCE_GAP_EV[governing].gap_ev
    return StateVerdict(
        temp_k=temp_k,
        n_e=n_e,
        critical_gap_ev=d_crit,
        governing=governing,
        governing_gap_ev=gap,
        margin=n_e / mcwhirter_threshold(temp_k, gap),
        abundant=tuple(present),
        failing=tuple(s for s in present if RESONANCE_GAP_EV[s].gap_ev > d_crit),
    )


def critical_gap(temp_k: float, n_e: float) -> float:
    """Largest level gap [eV] still collisionally controlled at `(temp_k, n_e)`.

    The criterion inverted: `dE_crit = (n_e / (1.6e18 * sqrt(T)))^(1/3)`. Reporting the state this
    way needs **no atomic data** — a state is characterized by one number, and the species'
    resonance gaps enter only afterwards as a comparison. That keeps the primary result
    independent of line-list accuracy, which is worth having: an earlier draft of
    `RESONANCE_GAP_EV` had O III wrong by 16% (see the selection traps documented there).
    """
    return float((n_e / (MCWHIRTER_COEFF_CGS * CM3_PER_M3 * math.sqrt(temp_k))) ** (1.0 / 3.0))


@dataclass(frozen=True)
class ProbeState:
    """One turnaround state recorded by a `--frozen-probe-*` run."""

    v: float
    rho_impact: float
    rho_star: float
    t_star: float


def heavy_particle_fractions(comp: eos_water.Composition) -> dict[str, float]:
    """Species fractions of the heavy-particle (non-electron) population, keyed as in
    `RESONANCE_GAP_EV`. Electrons are excluded: they are the perturber, not a level system."""
    pops: dict[str, float] = {
        "H2O": comp.n_h2o,
        "H": comp.n_h,
        "O": comp.n_o,
        "H+": comp.n_hp,
    }
    for k, n in enumerate(comp.n_o_ions):
        pops[f"O{k + 1}+"] = n
    total = sum(pops.values())
    if total <= 0.0:
        return dict.fromkeys(pops, 0.0)
    return {s: n / total for s, n in pops.items()}


def read_probe(path: Path) -> list[ProbeState]:
    """Read a `--frozen-probe-*` JSONL, tolerating blank lines."""
    return [
        ProbeState(
            v=float(d["v"]),
            rho_impact=float(d["rho_impact"]),
            rho_star=float(d["rho_star"]),
            t_star=float(d["t_star"]),
        )
        for line in path.read_text().splitlines()
        if line.strip()
        for d in (json.loads(line),)
    ]


def evaluate_probe(path: Path) -> list[tuple[ProbeState, StateVerdict]]:
    """Evaluate the McWhirter criterion at every turnaround state in a probe file.

    The composition comes from `eos_water.composition` — the *same* equilibrium Saha ladder the
    sweep's EOS uses — so the criterion is applied to the populations actually assumed rather than
    to an independent estimate of them.
    """
    out = []
    for st in read_probe(path):
        comp = eos_water.composition(st.rho_star, st.t_star)
        out.append((st, evaluate_state(st.t_star, comp.n_e, heavy_particle_fractions(comp))))
    return out


DEFAULT_PROBE_PATHS = (
    Path("data/results/frozen_probe_heavyplate.jsonl"),
    Path("data/results/frozen_probe_jupiter.jsonl"),
)
DEFAULT_SUMMARY_PATH = Path("data/results/lte_validity.csv")

CSV_HEADER = (
    "probe,v,rho_impact,rho_star,t_star,n_e,critical_gap_ev,"
    "governing,governing_label,governing_gap_ev,margin,lte_valid,failing"
)


def write_summary(
    probes: tuple[Path, ...] = DEFAULT_PROBE_PATHS, path: Path = DEFAULT_SUMMARY_PATH
) -> None:
    """Write the per-state LTE verdict table (the Q5 deliverable)."""
    lines = [CSV_HEADER]
    for probe in probes:
        if not probe.exists():
            continue
        for st, v in evaluate_probe(probe):
            label = RESONANCE_GAP_EV[v.governing].label if v.governing in RESONANCE_GAP_EV else "-"
            lines.append(
                f"{probe.stem},{st.v},{st.rho_impact},{st.rho_star},{st.t_star},{v.n_e:.6e},"
                f"{v.critical_gap_ev:.4f},{v.governing},{label},{v.governing_gap_ev:.4f},"
                f"{v.margin:.6e},{v.lte_valid},{'|'.join(v.failing)}"
            )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    """CLI: `python -m puffsat.lte` -> data/results/lte_validity.csv + a stdout summary."""
    write_summary()
    worst: tuple[float, str] | None = None
    n_fail = 0
    for probe in DEFAULT_PROBE_PATHS:
        if not probe.exists():
            continue
        for st, v in evaluate_probe(probe):
            if not v.lte_valid:
                n_fail += 1
                continue
            tag = f"v={st.v / 1000:.0f} km/s rho={st.rho_impact:.3f} ({v.governing})"
            if worst is None or v.margin < worst[0]:
                worst = (v.margin, tag)
    if worst is not None:
        print(f"python: tightest passing McWhirter margin {worst[0]:.2f}x at {worst[1]}")
    print(
        f"python: {n_fail} state(s) fail on a trace high-charge stage; wrote {DEFAULT_SUMMARY_PATH}"
    )


if __name__ == "__main__":
    main()
