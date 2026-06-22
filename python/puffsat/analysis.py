"""Frontier extraction + plots for the 16 km/s sweep (B5d-2).

Reads the sweep's JSONL (`crates/sweep` output, ADR-0019), builds the `e_eff(rho_impact)` frontier
and the per-density loss decomposition (each loss channel as a fraction of the total loss), writes a
summary CSV, and renders the figures. The extraction is **stdlib-only and unit-tested** (no science
deps); plotting imports matplotlib lazily (the `sci` extra), so this module imports and tests fine
without it.

This pass plots `e_eff` vs `rho` and the stacked loss decomposition; the opacity-sensitivity overlay
is added by B5d-3 (which generates the 0.1x/10x-opacity comparison data).
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, fields
from itertools import pairwise
from pathlib import Path

DEFAULT_SWEEP_PATH = Path("data/results/sweep.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/frontier.csv")
DEFAULT_PLOT_DIR = Path("data/results")


@dataclass(frozen=True)
class SweepRow:
    """One sweep result row (the fields the analysis needs from the JSONL schema)."""

    rho_impact: float
    v: float
    e_eff: float
    loss_radiative_wall: float
    loss_escape_space: float
    loss_conductive: float


@dataclass(frozen=True)
class FrontierPoint:
    """One point on the `e_eff(rho)` frontier with its loss decomposition (fractions sum to 1 when
    there is any loss, else all 0)."""

    rho_impact: float
    e_eff: float
    total_loss: float
    frac_radiative_wall: float
    frac_escape_space: float
    frac_conductive: float


def read_sweep(path: Path = DEFAULT_SWEEP_PATH) -> list[SweepRow]:
    """Parse the JSONL sweep results (one JSON object per line; blank lines tolerated)."""
    rows: list[SweepRow] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        rows.append(
            SweepRow(
                rho_impact=float(d["rho_impact"]),
                v=float(d["v"]),
                e_eff=float(d["e_eff"]),
                loss_radiative_wall=float(d["loss_radiative_wall"]),
                loss_escape_space=float(d["loss_escape_space"]),
                loss_conductive=float(d["loss_conductive"]),
            )
        )
    return rows


def frontier(rows: list[SweepRow]) -> list[FrontierPoint]:
    """Build the `e_eff(rho)` frontier (ascending in rho) with its loss decomposition."""
    points: list[FrontierPoint] = []
    for r in sorted(rows, key=lambda row: row.rho_impact):
        total = r.loss_radiative_wall + r.loss_escape_space + r.loss_conductive
        if total > 0.0:
            fa, fb, fc = (
                r.loss_radiative_wall / total,
                r.loss_escape_space / total,
                r.loss_conductive / total,
            )
        else:
            fa = fb = fc = 0.0
        points.append(
            FrontierPoint(
                rho_impact=r.rho_impact,
                e_eff=r.e_eff,
                total_loss=total,
                frac_radiative_wall=fa,
                frac_escape_space=fb,
                frac_conductive=fc,
            )
        )
    return points


def write_summary(points: list[FrontierPoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write the frontier to a CSV (header = the `FrontierPoint` field names), ascending in rho."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [f.name for f in fields(FrontierPoint)]
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for p in points:
            writer.writerow([getattr(p, name) for name in header])


def plot_frontier(points: list[FrontierPoint], out_dir: Path = DEFAULT_PLOT_DIR) -> list[Path]:
    """Render `e_eff` vs rho and the stacked loss decomposition. Returns the saved figure paths.

    matplotlib is imported lazily (the `sci` extra) and forced onto the headless `Agg` backend so a
    `make analysis` run needs no display. The numeric inputs are plain Python floats built above, so
    no `Any` from the untyped matplotlib API escapes this function.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    rho = [p.rho_impact for p in points]
    e_eff = [p.e_eff for p in points]
    saved: list[Path] = []

    # --- e_eff(rho) frontier ---
    fig, ax = plt.subplots()
    ax.plot(rho, e_eff, "o-")
    ax.set_xlabel(r"$\rho_\mathrm{impact}$ [kg/m$^3$]")
    ax.set_ylabel(r"$e_\mathrm{eff}$")
    ax.set_title(r"Restitution frontier $e_\mathrm{eff}(\rho)$ at 16 km/s")
    ax.grid(True, alpha=0.3)
    frontier_path = out_dir / "e_eff_frontier.png"
    fig.savefig(frontier_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(frontier_path)

    # --- stacked loss decomposition (fractions of total loss) ---
    fa = [p.frac_radiative_wall for p in points]
    fb = [p.frac_escape_space for p in points]
    fc = [p.frac_conductive for p in points]
    base_b = list(fa)
    base_c = [a + b for a, b in zip(fa, fb, strict=True)]
    fig, ax = plt.subplots()
    width = 0.6 * min((b - a for a, b in pairwise(rho)), default=1.0)
    ax.bar(rho, fa, width=width, label="1a radiative→wall")
    ax.bar(rho, fb, width=width, bottom=base_b, label="1b escape→space")
    ax.bar(rho, fc, width=width, bottom=base_c, label="2 conductive (deferred)")
    ax.set_xlabel(r"$\rho_\mathrm{impact}$ [kg/m$^3$]")
    ax.set_ylabel("fraction of total loss")
    ax.set_title("Loss-channel decomposition (ADR-0016)")
    ax.legend(loc="upper right", fontsize="small")
    losses_path = out_dir / "loss_decomposition.png"
    fig.savefig(losses_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved.append(losses_path)

    return saved


def main() -> None:
    """Read the sweep, extract the frontier, write the CSV summary, and render the figures."""
    rows = read_sweep()
    points = frontier(rows)
    write_summary(points)
    figs = plot_frontier(points)
    print(
        f"python: wrote {DEFAULT_SUMMARY_PATH} and {len(figs)} figure(s): "
        + ", ".join(map(str, figs))
    )


if __name__ == "__main__":
    main()
