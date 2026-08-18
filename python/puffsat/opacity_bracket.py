"""The reported opacity uncertainty band on `f` (Q18/Q10 of the 16-63 km/s extension).

The sweep inherited a **0.1x-10x** opacity bracket. That range was sized when the table carried a
*placeholder* Kramers opacity that CONCLUSION.md records as ~2000x low — it expressed ignorance
about a stand-in, not uncertainty in a measurement. The table now carries real TOPS/OPLIB gray
means, so the bracket should express what is actually uncertain about *those*, which is far less.

Three inputs, each measured or cited rather than assumed:

1. **Published accuracy of the opacity data.** Farag et al. 2024 (ApJ 968:16, arXiv:2406.02845)
   section 2.1, citing Huebner & Barfield (2014), gives the OPLIB Rosseland-mean uncertainty by
   dominant process — see `OPLIB_UNCERTAINTY`. The paper also notes the uncertainty falls as ionic
   charge rises toward the hydrogenic limit, which is why the hot, highly-ionized states here take
   the tighter number.

2. **Non-LTE contribution — measured and negligible.** Q5's McWhirter check (`lte.py`) finds every
   abundant species collisionally controlled across the band; the only failure is helium-like
   O VII at ~1% abundance in the most dilute 69 km/s state. Propagating even a 100% error in that
   stage's population gives `dZbar/Zbar = 0.71%`, and with the table's measured `dln(kappa)/
   dln(Zbar)` of order unity that is `dkappa/kappa <= 0.93%` — an order below the data's own
   accuracy, and identically zero at every other probed state. **The bracket is set by the opacity
   data, not by the LTE question.**

3. **Sensitivity of the answer.** `de_eff/dln(kappa)` is measured from the sweep's own
   opacity-scale rows (`measure_slopes`), which span 0.0009 to 0.071 across the probed states —
   growing steeply with velocity and cloud length, falling steeply with impact density.

The band is then `d(f) = eta_capture * de_eff/dln(kappa) * ln(1 + fraction) / 2`.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# Farag et al. 2024 (arXiv:2406.02845) sec. 2.1, quoting Huebner & Barfield (2014), verbatim:
#   "estimates of the uncertainty in the opacity are ~5% when electron scattering dominates at
#    high T and low rho. As the rho increases and free-free processes become more important, the
#    uncertainty is less than ~10%. As T decreases and bound-free processes become important, the
#    uncertainty increases to ~20%. As T decreases further, bound-bound processes can contribute,
#    and the uncertainty rises to ~30%."
OPLIB_UNCERTAINTY: dict[str, float] = {
    "electron_scattering": 0.05,
    "free_free": 0.10,
    "bound_free": 0.20,
    "bound_bound": 0.30,
}

# Code-to-code spread for cross-check (same paper, abstract): OPLIB vs OPAL differ by "~20-80%
# across individual chemical mixtures", narrowing to ~7-15% in the well-constrained solar interior.
# Reported alongside the bracket rather than used as it: those mixtures are stellar, and the 80%
# tail is driven by compositions unlike water.
OPLIB_VS_OPAL = (0.20, 0.80)

# Mean oxygen charge above which the plasma is lithium-like or better, where bound-bound line
# forests have thinned out and the cited uncertainty falls toward the bound-free value.
ZBAR_IONIZED = 2.0


def bracket_fraction(temp_k: float, zbar: float) -> float:
    """Fractional opacity uncertainty applicable at a state of temperature `temp_k` and mean
    oxygen charge `zbar`.

    Cool and barely ionized means a dense forest of bound-bound lines and the widest bracket;
    hot and stripped means bound-free/free-free and the tighter one. `temp_k` is accepted (and
    unused below the charge test) because the cited classification is stated in terms of falling
    temperature — `zbar` is simply the sharper proxy for the same thing in this mixture.
    """
    if zbar >= ZBAR_IONIZED:
        return OPLIB_UNCERTAINTY["bound_free"]
    return OPLIB_UNCERTAINTY["bound_bound"]


# A bounce whose e_eff falls below this fraction of its own opacity-row median did not merely
# respond to opacity — it stalled. Legitimate opacity response across 0.1x-10x never exceeds a
# factor of ~1.3 in the measured data, so this threshold cannot fire on physical variation.
STALL_FRACTION = 0.5


def measure_slopes_checked(
    sweep_path: Path,
) -> tuple[dict[tuple[float, float, float], float], list[tuple[float, float, float, float]]]:
    """`de_eff/dln(kappa)` per `(v, rho_impact, length)`, plus the stalled rows rejected.

    Uses the widest symmetric pair available around `kappa = 1` so the slope is centered on the
    real-opacity row rather than extrapolated from one side, skipping any pair that includes a
    stalled row.

    A stalled run reports a truncated impulse integral, so its `e_eff` is not a physical opacity
    response, and fitting through it can flip the slope's sign. Such rows are excluded and returned
    so the caller can report them rather than silently dropping evidence of a solver problem.

    Rows carrying the kernel's `converged` flag are judged by it. The magnitude heuristic is only
    the fallback for older files written before the flag existed, and it can catch a stall solely
    when the artifact happens to look implausible — a stall that lands on a believable number would
    slip past it.
    """
    grouped: dict[tuple[float, float, float], dict[float, float]] = defaultdict(dict)
    flagged: set[tuple[float, float, float, float]] = set()
    for line in sweep_path.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        key = (float(r["v"]), float(r["rho_impact"]), float(r["length"]))
        scale = float(r["opacity_scale"])
        grouped[key][scale] = float(r["e_eff"])
        if r.get("converged") is False:
            flagged.add((*key, scale))

    out: dict[tuple[float, float, float], float] = {}
    stalled: list[tuple[float, float, float, float]] = []
    for key, by_scale in sorted(grouped.items()):
        explicit = {s for s in by_scale if (*key, s) in flagged}
        if explicit:
            bad = explicit
        else:
            values = sorted(by_scale.values())
            median = values[len(values) // 2]
            bad = {s for s, e in by_scale.items() if e < STALL_FRACTION * median}
        stalled.extend(sorted((*key, s) for s in bad))
        for lo, hi in ((0.1, 10.0), (0.3, 3.0)):
            if lo in by_scale and hi in by_scale and lo not in bad and hi not in bad:
                out[key] = (by_scale[hi] - by_scale[lo]) / (math.log(hi) - math.log(lo))
                break
    return out, stalled


def measure_slopes(sweep_path: Path) -> dict[tuple[float, float, float], float]:
    """`de_eff/dln(kappa)` per `(v, rho_impact, length)`, stalled rows excluded."""
    return measure_slopes_checked(sweep_path)[0]


@dataclass(frozen=True)
class FBand:
    """The opacity uncertainty band on one `f` value."""

    fraction: float
    delta_e_eff: float
    delta_f: float


def f_band(slope: float, fraction: float, eta_capture: float) -> FBand:
    """Translate a fractional opacity uncertainty into a band on `f`.

    `f = eta_capture * (1 + e_eff) / 2`, so `d(f) = eta_capture * d(e_eff) / 2`, and
    `d(e_eff) = slope * ln(1 + fraction)` since the slope is per natural log of kappa.
    """
    d_e = slope * math.log1p(fraction)
    return FBand(fraction=fraction, delta_e_eff=d_e, delta_f=eta_capture * d_e / 2.0)


# Each sweep is paired with *every* probe file that can supply a turnaround state for it. The
# heavy-plate arm needs two: the Q4 freeze-bracket probe covers 16-28 km/s at the full density
# grid, while the diagnostic probe covers the 45-63 km/s tau-check states. Those live in a separate
# file because the frozen-table builder emits one Saha table per probe row, and those tables are
# only wanted at the three freeze anchors.
DEFAULT_SWEEPS: tuple[tuple[Path, tuple[Path, ...]], ...] = (
    (
        Path("data/results/sweep_heavyplate.jsonl"),
        (
            Path("data/results/frozen_probe_heavyplate.jsonl"),
            Path("data/results/frozen_probe_heavyplate_diag.jsonl"),
        ),
    ),
    (
        Path("data/results/sweep_jupiter.jsonl"),
        (Path("data/results/frozen_probe_jupiter.jsonl"),),
    ),
)
DEFAULT_SUMMARY_PATH = Path("data/results/opacity_bracket.csv")

# eta_capture along the survivable cloud schedule (design §7; the geometry sweep gives 0.967-0.992
# there). Held fixed here because this module reports the *opacity* contribution to the f band.
ETA_NOMINAL = 0.98

CSV_HEADER = (
    "sweep,v,rho_impact,length,zbar,t_star,de_eff_dlnkappa,bracket_fraction,"
    "delta_e_eff,delta_f,delta_f_inherited"
)


def _zbar_by_state(
    *probe_paths: Path,
) -> dict[tuple[float, float], tuple[float, float]]:
    """`(t_star, mean oxygen charge)` keyed by `(v, rho_impact)`, merged over probe files.

    Missing files are skipped rather than raising: a probe that has not been generated yet costs
    coverage, and the caller already drops states it cannot find a composition for.
    """
    from puffsat import eos_water

    out: dict[tuple[float, float], tuple[float, float]] = {}
    for probe_path in probe_paths:
        if not probe_path.exists():
            continue
        for line in probe_path.read_text().splitlines():
            if not line.strip():
                continue
            d = json.loads(line)
            comp = eos_water.composition(float(d["rho_star"]), float(d["t_star"]))
            o_tot = comp.n_o + sum(comp.n_o_ions)
            zbar = (
                sum((k + 1) * n for k, n in enumerate(comp.n_o_ions)) / o_tot if o_tot > 0 else 0.0
            )
            out[(float(d["v"]), float(d["rho_impact"]))] = (float(d["t_star"]), zbar)
    return out


def write_summary(path: Path = DEFAULT_SUMMARY_PATH) -> list[tuple[float, float, float, float]]:
    """Write the per-state opacity band on `f`; return any stalled rows encountered."""
    lines = [CSV_HEADER]
    all_stalled: list[tuple[float, float, float, float]] = []
    for sweep_path, probe_paths in DEFAULT_SWEEPS:
        if not sweep_path.exists() or not any(p.exists() for p in probe_paths):
            continue
        slopes, stalled = measure_slopes_checked(sweep_path)
        all_stalled.extend(stalled)
        states = _zbar_by_state(*probe_paths)
        for (v, rho, length), slope in sorted(slopes.items()):
            if (v, rho) not in states:
                continue
            t_star, zbar = states[(v, rho)]
            frac = bracket_fraction(t_star, zbar)
            band = f_band(slope, frac, ETA_NOMINAL)
            inherited = f_band(slope, 9.0, ETA_NOMINAL)  # the old 1x -> 10x arm
            lines.append(
                f"{sweep_path.stem},{v},{rho},{length},{zbar:.4f},{t_star:.1f},{slope:.6e},"
                f"{frac:.2f},{band.delta_e_eff:.6e},{band.delta_f:.6e},{inherited.delta_f:.6e}"
            )
    path.write_text("\n".join(lines) + "\n")
    return all_stalled


def main() -> None:
    """CLI: `python -m puffsat.opacity_bracket` -> data/results/opacity_bracket.csv."""
    stalled = write_summary()
    for v, rho, length, scale in stalled:
        print(
            f"python: WARNING stalled row excluded — v={v / 1000:.0f} km/s rho={rho:.3f} "
            f"L={length:.0f} kappa={scale}x (table rho ceiling, Q6/ADR-0035)"
        )
    print(f"python: wrote {DEFAULT_SUMMARY_PATH}")


if __name__ == "__main__":
    main()
