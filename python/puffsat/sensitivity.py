"""Opacity-insensitivity scan (B5d-3) — the methodological gate licensing the interim opacity.

Design §3 predicts that at 16 km/s the stagnated gas is `tau >> 1`, so `e_eff` is EOS- and
gas-dynamics-dominated and **insensitive to the opacity**; the real ADR-0007 opacity table is
therefore deferrable without compromising the `e_eff(rho)` number. This module *demonstrates* that:
it regenerates the water table at several opacity scales (0.1x / 1x / 10x), re-runs the sweep on
each, and compares the frontiers. If `e_eff` barely moves while the radiative loss channels (1a/1b)
scale with the opacity, the interim bracket (B5c-2) is justified.

The comparison core (`compare`) is pure and unit-tested; the scan (`run_scan`/`main`) is the
experiment — it shells out to the Rust sweep binary and is driven by `make sensitivity`, not by the
test suite.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from puffsat import tables
from puffsat.analysis import SweepRow, read_sweep

SCALES: tuple[float, ...] = (0.1, 1.0, 10.0)
REFERENCE_SCALE = 1.0
SCAN_DIR = Path("data/results/opacity_scan")
DEFAULT_REPORT_PATH = SCAN_DIR / "sensitivity.json"
# `e_eff` is deemed opacity-insensitive if it moves by less than this (relative) across the scales.
PASS_THRESHOLD = 0.05

# The Rust sweep binary, invoked as `<cmd...> <table_path> <result_path>`.
SWEEP_CMD: tuple[str, ...] = ("cargo", "run", "--release", "--quiet", "-p", "sweep", "--")


@dataclass(frozen=True)
class SensitivityResult:
    """Outcome of the opacity scan: per-scale frontiers and the `e_eff` spread across them."""

    scales: list[float]
    rho: list[float]
    e_eff_by_scale: dict[float, list[float]]
    total_loss_by_scale: dict[float, list[float]]
    max_abs_de_eff: float
    max_rel_de_eff: float

    @property
    def passes(self) -> bool:
        """True if `e_eff` is opacity-insensitive at the [`PASS_THRESHOLD`] relative tolerance."""
        return self.max_rel_de_eff < PASS_THRESHOLD


def compare(rows_by_scale: dict[float, list[SweepRow]]) -> SensitivityResult:
    """Compare per-scale sweep results, measuring how far `e_eff` moves from the 1x reference.

    All scales must share the same (sorted) `rho` grid, and the reference scale must be present.
    Returns the max absolute and relative `e_eff` excursion over every `(scale, rho)`.
    """
    if REFERENCE_SCALE not in rows_by_scale:
        raise ValueError(f"reference scale {REFERENCE_SCALE} missing from the scan")

    scales = sorted(rows_by_scale)
    by_scale_sorted = {s: sorted(rows_by_scale[s], key=lambda r: r.rho_impact) for s in scales}
    rho = [r.rho_impact for r in by_scale_sorted[REFERENCE_SCALE]]
    for s in scales:
        grid = [r.rho_impact for r in by_scale_sorted[s]]
        if grid != rho:
            raise ValueError(f"scale {s} rho grid {grid} != reference {rho}")

    e_eff_by_scale = {s: [r.e_eff for r in by_scale_sorted[s]] for s in scales}
    total_loss_by_scale = {
        s: [
            r.loss_radiative_wall + r.loss_escape_space + r.loss_conductive
            for r in by_scale_sorted[s]
        ]
        for s in scales
    }

    ref = e_eff_by_scale[REFERENCE_SCALE]
    max_abs = 0.0
    max_rel = 0.0
    for s in scales:
        for got, base in zip(e_eff_by_scale[s], ref, strict=True):
            abs_d = abs(got - base)
            max_abs = max(max_abs, abs_d)
            if base != 0.0:
                max_rel = max(max_rel, abs_d / abs(base))
    return SensitivityResult(
        scales=scales,
        rho=rho,
        e_eff_by_scale=e_eff_by_scale,
        total_loss_by_scale=total_loss_by_scale,
        max_abs_de_eff=max_abs,
        max_rel_de_eff=max_rel,
    )


def run_scan(
    scales: tuple[float, ...] = SCALES,
    scan_dir: Path = SCAN_DIR,
    sweep_cmd: tuple[str, ...] = SWEEP_CMD,
) -> dict[float, list[SweepRow]]:
    """Generate an opacity-scaled table and run the sweep for each scale; return the parsed rows.

    Each scale gets its own `water_{scale}.json` table and `sweep_{scale}.jsonl` result, so the runs
    are independent and the artifacts are inspectable.
    """
    scan_dir.mkdir(parents=True, exist_ok=True)
    rows_by_scale: dict[float, list[SweepRow]] = {}
    for scale in scales:
        tag = f"{scale:g}"
        table_path = scan_dir / f"water_{tag}.json"
        result_path = scan_dir / f"sweep_{tag}.jsonl"

        table = tables.build_table(kappa_scale=scale)
        with table_path.open("w") as fh:
            json.dump(table, fh)

        subprocess.run([*sweep_cmd, str(table_path), str(result_path)], check=True)
        rows_by_scale[scale] = read_sweep(result_path)
    return rows_by_scale


def write_report(result: SensitivityResult, path: Path = DEFAULT_REPORT_PATH) -> None:
    """Write the scan verdict + per-scale frontiers as JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "pass_threshold_rel": PASS_THRESHOLD,
        "max_abs_de_eff": result.max_abs_de_eff,
        "max_rel_de_eff": result.max_rel_de_eff,
        "passes": result.passes,
        "scales": result.scales,
        "rho": result.rho,
        "e_eff_by_scale": {f"{s:g}": result.e_eff_by_scale[s] for s in result.scales},
        "total_loss_by_scale": {f"{s:g}": result.total_loss_by_scale[s] for s in result.scales},
    }
    with path.open("w") as fh:
        json.dump(report, fh, indent=2)


def plot_overlay(result: SensitivityResult, out_dir: Path = SCAN_DIR) -> Path:
    """Overlay `e_eff(rho)` for every opacity scale (lazy matplotlib, headless `Agg`)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots()
    for s in result.scales:
        ax.plot(result.rho, result.e_eff_by_scale[s], "o-", label=rf"$\kappa \times {s:g}$")
    ax.set_xlabel(r"$\rho_\mathrm{impact}$ [kg/m$^3$]")
    ax.set_ylabel(r"$e_\mathrm{eff}$")
    ax.set_title(
        rf"Opacity insensitivity: max $\Delta e_\mathrm{{eff}}$ = {result.max_rel_de_eff:.2%}"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = out_dir / "opacity_sensitivity.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    """Run the opacity scan, compare, write the report + overlay, and print the verdict."""
    result = compare(run_scan())
    write_report(result)
    overlay = plot_overlay(result)
    verdict = "PASS — e_eff is opacity-insensitive" if result.passes else "FAIL — e_eff moved"
    print(
        f"python: opacity scan {result.scales}: "
        f"max |Δe_eff| = {result.max_abs_de_eff:.2e} ({result.max_rel_de_eff:.2%} rel) -> {verdict}"
    )
    print(f"python: wrote {DEFAULT_REPORT_PATH} and {overlay}")


if __name__ == "__main__":
    main()
