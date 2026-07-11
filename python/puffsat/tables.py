"""Assemble the unified EOS/opacity table (ADR-0007) the Rust kernel loads (B5c-2).

The equilibrium water EOS (`eos_water`) supplies the load-bearing `(p, e, c_s)` that governs
`e_eff`. The opacity attached here is an explicit **interim bracket**: a Kramers-shaped model whose
coefficient is calibrated so the column optical depth `tau = rho * kappa_R * L` lands in the
design's `[1e2, 1e5]` band at the nominal stagnation (design §3). At `tau >> 1` `e_eff` is
opacity-insensitive
(demonstrated by the B5d-3 sensitivity sweep), so this bracket is adequate for the `e_eff(rho)`
deliverable; the **real** opacity table (CEA/TOPS/OPLIB/ExoMol — gated external data) is a separate
later rung feeding the survivability flux. Provenance records exactly this provisional status.

Output is the ADR-0007 JSON the Rust loader (`crates/tables`) validates: ascending positive grids,
`shape == [n_rho, n_T]`, and five positive-definite fields flattened **row-major over (rho, T)**
(rho outer, T inner — matching the loader's `i*n_T + j` indexing).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from puffsat import eos_water as ew
from puffsat import tops
from puffsat.eos_water import Vec

# --- (rho, T) grid (SI). rho spans the impact range plus strong-shock compression; T runs from
#     cold vapor to a hot plasma. Both log-spaced — the loader interpolates in log-log. ---
RHO_RANGE = (0.01, 20.0)  # kg/m^3
N_RHO = 40
T_RANGE = (300.0, 60_000.0)  # K
N_T = 60

# --- Interim bracketing opacity (PROVISIONAL). Kramers free-free/bound-free *shape* only:
#     kappa_R = KAPPA0_R (rho/rho_ref) (T/T_ref)^-3.5. The coefficient is fixed by requiring the
#     column optical depth tau = rho*kappa_R*L equal TAU_TARGET at the nominal stagnation point:
#     tau = rho_ref*KAPPA0_R*L_ref  =>  KAPPA0_R = TAU_TARGET/(rho_ref*L_ref). ---
KAPPA_RHO_REF = 1.0  # kg/m^3 — representative stagnated (shock-compressed) density
KAPPA_T_REF = 2.0e4  # K       — representative stagnation temperature
KAPPA_L_REF = 1.0  # m       — nominal slug column length
TAU_TARGET = 1.0e3  # mid-band of the design's [1e2, 1e5] optical-depth range
KAPPA0_R = TAU_TARGET / (KAPPA_RHO_REF * KAPPA_L_REF)  # [m^2/kg]
PLANCK_OVER_ROSSELAND = 3.0  # kappa_P >= kappa_R (Planck mean >= Rosseland mean)
TAU_BAND = (1.0e2, 1.0e5)

DEFAULT_TABLE_PATH = Path("data/tables/water.json")

# --- Jupiter-retrograde 69 km/s scenario table (special-scenario rung): the stagnated plasma sits
#     at ~1.3-1.8e5 K with oxygen at Z_bar ~ 4-4.5 (multi-stage Saha ladder), so the T grid must
#     reach well past it; the survivable impact densities are ~0.02-0.16 kg/m^3 and the rebound
#     rarefies far below the impact density, so the rho grid extends much lower too. ---
RHO_RANGE_JUPITER = (1.0e-4, 30.0)  # kg/m^3
N_RHO_JUPITER = 48
T_RANGE_JUPITER = (300.0, 1.2e6)  # K
N_T_JUPITER = 88
# Thomson electron-scattering floor on the Rosseland mean [m^2/kg]: sigma_T * (10 e-/molecule)
# / m_H2O ~ 0.022 for fully-ionized water. The Kramers T^-3.5 shape underflows at 1e5+ K; free
# electrons still scatter, so radiative *transport* never gets cheaper than Thomson. Applied to
# kappa_R only (scattering does not emit -> kappa_P keeps the Kramers shape).
KAPPA_THOMSON_FLOOR = 0.022
DEFAULT_TABLE_PATH_JUPITER = Path("data/tables/water_jupiter.json")

# --- Low-v (Rung C) cool-gas two-phase table: real-fluid water across the saturation dome
#     (CoolProp), a generous (rho, T) box (log-spaced packs more points at low T, near the dome edge
#     T_crit = 647 K). Radiation is off at 3.2 km/s (design §3), so the opacities are transparent.
RHO_RANGE_LOWV = (0.01, 100.0)  # kg/m^3 — re-expanded vapor up through compressed/condensing gas
N_RHO_LOWV = 48
T_RANGE_LOWV = (280.0, 1800.0)  # K — cold vapor through the ~1700 K stagnation (supercritical top)
N_T_LOWV = 72
KAPPA_TRANSPARENT = 1.0e-10  # m^2/kg — radiation off at low-v; kept positive for the loader
DEFAULT_TABLE_PATH_LOWV = Path("data/tables/water_lowv.json")


def opacity_grid(rho_grid: Vec, t_grid: Vec, kappa_scale: float = 1.0) -> tuple[Vec, Vec]:
    """Bracketing `(kappa_R, kappa_P)` [m^2/kg] on the `(rho, T)` grid, shaped `(n_rho, n_T)`.

    `kappa_R = kappa_scale * KAPPA0_R (rho/rho_ref) (T/T_ref)^-3.5`;
    `kappa_P = PLANCK_OVER_ROSSELAND * kappa_R`. Positive-definite (the loader interpolates
    `ln kappa`). PROVISIONAL — see the module docstring. `kappa_scale` multiplies both means by a
    constant factor: the B5d-3 opacity-insensitivity scan generates 0.1x / 1x / 10x tables to show
    `e_eff` does not move with it (which licenses this interim bracket).
    """
    kappa_r = (
        kappa_scale
        * KAPPA0_R
        * (rho_grid[:, None] / KAPPA_RHO_REF)
        * (t_grid[None, :] / KAPPA_T_REF) ** -3.5
    )
    kappa_p = PLANCK_OVER_ROSSELAND * kappa_r
    return kappa_r, kappa_p


def _flatten(a: Vec) -> list[float]:
    """Row-major (C-order) flatten to a plain `list[float]` for JSON (no numpy scalars leak)."""
    return [float(x) for x in a.reshape(-1)]


def _provenance(
    rho_range: tuple[float, float],
    n_rho: int,
    t_range: tuple[float, float],
    n_t: int,
    kappa_scale: float,
) -> dict[str, object]:
    """Human-readable provenance (free-form to the Rust loader). Flags the opacity PROVISIONAL."""
    return {
        "schema": "ADR-0007",
        "generated_by": "puffsat.tables.build_table",
        "grid": {
            "rho_range": list(rho_range),
            "n_rho": n_rho,
            "T_range": list(t_range),
            "n_T": n_t,
            "spacing": "geometric (log-spaced) in both axes",
            "units": "rho kg/m^3, T K, p Pa, e J/kg, c_s m/s, kappa m^2/kg",
        },
        "eos": {
            "model": "chemical-equilibrium water (dissociation H2O<=>2H+O + Saha ionization)",
            "species": ["H2O", "H", "O", "H+", "O+", "e-"],
            "energy_reference": "bound molecular H2O at T->0 = 0 (all e > 0)",
            "reference": (
                "Zel'dovich & Raizer, Physics of Shock Waves and High-Temperature "
                "Hydrodynamic Phenomena, Ch. III"
            ),
            "status": "analytic equilibrium EOS — the load-bearing physics for e_eff",
        },
        "opacity": {
            "status": "PROVISIONAL — bracketing placeholder, NOT a measured/computed opacity table",
            "model": (
                "Kramers shape kappa_R = KAPPA0_R (rho/rho_ref)(T/T_ref)^-3.5; "
                "kappa_P = PLANCK_OVER_ROSSELAND * kappa_R"
            ),
            "calibration": (
                "KAPPA0_R fixed so tau = rho*kappa_R*L = TAU_TARGET at the nominal stagnation "
                "(rho_ref, T_ref, L_ref)"
            ),
            "kappa0_rosseland": KAPPA0_R,
            "kappa_scale": kappa_scale,
            "rho_ref": KAPPA_RHO_REF,
            "T_ref": KAPPA_T_REF,
            "L_ref": KAPPA_L_REF,
            "tau_target": TAU_TARGET,
            "tau_band_design": list(TAU_BAND),
            "planck_over_rosseland": PLANCK_OVER_ROSSELAND,
            "real_table_rung": (
                "deferred — CEA/TOPS/OPLIB/ExoMol survivability-flux rung (design §9); "
                "e_eff is opacity-insensitive at tau>>1 (verified by the B5d-3 sweep)"
            ),
        },
    }


def build_table(
    rho_range: tuple[float, float] = RHO_RANGE,
    n_rho: int = N_RHO,
    t_range: tuple[float, float] = T_RANGE,
    n_t: int = N_T,
    kappa_scale: float = 1.0,
) -> dict[str, object]:
    """Build the ADR-0007 table dict: log-spaced grids, the equilibrium EOS fields, the bracketing
    opacity (scaled by `kappa_scale`), and nested provenance. Fields are flattened row-major over
    `(rho, T)`. `kappa_scale != 1` rescales only the opacity (the EOS fields are untouched) — the
    knob the B5d-3 insensitivity scan turns."""
    rho_grid = np.geomspace(rho_range[0], rho_range[1], n_rho)
    t_grid = np.geomspace(t_range[0], t_range[1], n_t)

    p, e, cs = ew.eos_grid(rho_grid, t_grid)
    kappa_r, kappa_p = opacity_grid(rho_grid, t_grid, kappa_scale)

    return {
        "rho_grid": [float(x) for x in rho_grid],
        "T_grid": [float(x) for x in t_grid],
        "shape": [n_rho, n_t],
        "fields": {
            "p": _flatten(p),
            "e": _flatten(e),
            "c_s": _flatten(cs),
            "kappa_rosseland": _flatten(kappa_r),
            "kappa_planck": _flatten(kappa_p),
        },
        "provenance": _provenance(rho_range, n_rho, t_range, n_t, kappa_scale),
    }


def build_table_jupiter(
    kappa_scale: float = 1.0, tops_path: Path | None = None
) -> dict[str, object]:
    """Build the Jupiter-retrograde (69 km/s) scenario table: the same equilibrium EOS on the
    extended `(rho, T)` grid (multi-stage O ladder engaged), with the interim Kramers opacity
    floored at Thomson scattering on the Rosseland mean.

    With `tops_path` (a saved TOPS gray pull, `puffsat.fetch_tops`), the REAL ATOMIC/OPLIB water
    means replace the interim opacity for `T >=` the OPLIB floor (~5802 K) — everywhere the 69 km/s
    bounce is radiatively active; the interim shape survives only in the cold molecular tail below.
    `kappa_scale` still multiplies the final field either way, so the sweep's 0.1x/10x bracket
    turns into a sensitivity band *around the real table* instead of around the placeholder."""
    rho_grid = np.geomspace(RHO_RANGE_JUPITER[0], RHO_RANGE_JUPITER[1], N_RHO_JUPITER)
    t_grid = np.geomspace(T_RANGE_JUPITER[0], T_RANGE_JUPITER[1], N_T_JUPITER)

    p, e, cs = ew.eos_grid(rho_grid, t_grid)
    kappa_r, kappa_p = opacity_grid(rho_grid, t_grid, kappa_scale)
    kappa_r = np.maximum(kappa_r, kappa_scale * KAPPA_THOMSON_FLOOR)
    if tops_path is not None:
        tops_pull = tops.load_tops_gray(tops_path)
        kappa_r_unscaled, kappa_p_unscaled = opacity_grid(rho_grid, t_grid, kappa_scale=1.0)
        kappa_r_unscaled = np.maximum(kappa_r_unscaled, KAPPA_THOMSON_FLOOR)
        kappa_r, kappa_p = tops.stitch_opacity(
            rho_grid, t_grid, kappa_r_unscaled, kappa_p_unscaled, tops_pull
        )
        kappa_r = kappa_scale * kappa_r
        kappa_p = kappa_scale * kappa_p

    prov = _provenance(RHO_RANGE_JUPITER, N_RHO_JUPITER, T_RANGE_JUPITER, N_T_JUPITER, kappa_scale)
    eos_prov = prov["eos"]
    assert isinstance(eos_prov, dict)
    eos_prov["species"] = ["H2O", "H", "O", "H+"] + [f"O{k}+" for k in range(1, 9)] + ["e-"]
    eos_prov["model"] = (
        "chemical-equilibrium water (dissociation H2O<=>2H+O + full multi-stage O Saha ladder)"
    )
    opac_prov = prov["opacity"]
    assert isinstance(opac_prov, dict)
    opac_prov["thomson_floor_rosseland"] = KAPPA_THOMSON_FLOOR
    if tops_path is not None:
        opac_prov["status"] = (
            "REAL for T >= 5802 K (the OPLIB floor): TOPS gray Rosseland/Planck means for "
            "water (2 H : 1 O atomic), LANL ATOMIC/OPLIB elemental opacities, "
            "aphysics2.lanl.gov; interim Kramers+Thomson survives only below the floor "
            "(cold molecular tail, radiatively inactive in this scenario)"
        )
        opac_prov["tops_pull"] = str(tops_path)
        opac_prov["real_table_rung"] = "landed for this scenario (ADR-0007 amendment)"
    prov["scenario"] = (
        "Jupiter-retrograde 69 km/s special scenario (2026-07): extended T grid past the "
        "~1.5e5 K stagnation, extended rho grid for the dilute survivable clouds"
    )

    return {
        "rho_grid": [float(x) for x in rho_grid],
        "T_grid": [float(x) for x in t_grid],
        "shape": [N_RHO_JUPITER, N_T_JUPITER],
        "fields": {
            "p": _flatten(p),
            "e": _flatten(e),
            "c_s": _flatten(cs),
            "kappa_rosseland": _flatten(kappa_r),
            "kappa_planck": _flatten(kappa_p),
        },
        "provenance": prov,
    }


DEFAULT_TABLE_DIR_FROZEN = Path("data/tables/frozen")


def frozen_table_name(v: float, rho_impact: float) -> str:
    """Per-case frozen-table file name — the contract shared with the Rust `--frozen` sweep."""
    return f"v{round(v):05d}_rho{rho_impact:.2f}.json"


def _build_table_frozen_common(
    y: ew.FrozenComposition,
    eos_provenance: dict[str, object],
    rho_range: tuple[float, float],
    n_rho: int,
    t_range: tuple[float, float],
    n_t: int,
) -> dict[str, object]:
    """Shared frozen-table assembly: frozen EOS fields + the interim bracketing opacity (unused by
    the EOS-only frozen bounces, kept for the loader contract) + provenance."""
    rho_grid = np.geomspace(rho_range[0], rho_range[1], n_rho)
    t_grid = np.geomspace(t_range[0], t_range[1], n_t)

    p, e, cs = ew.eos_grid_frozen(rho_grid, t_grid, y)
    kappa_r, kappa_p = opacity_grid(rho_grid, t_grid)

    prov = _provenance(rho_range, n_rho, t_range, n_t, kappa_scale=1.0)
    prov["generated_by"] = "puffsat.tables.build_table_frozen"
    prov["eos"] = eos_provenance

    return {
        "rho_grid": [float(x) for x in rho_grid],
        "T_grid": [float(x) for x in t_grid],
        "shape": [n_rho, n_t],
        "fields": {
            "p": _flatten(p),
            "e": _flatten(e),
            "c_s": _flatten(cs),
            "kappa_rosseland": _flatten(kappa_r),
            "kappa_planck": _flatten(kappa_p),
        },
        "provenance": prov,
    }


def build_table_frozen(
    rho_star: float,
    t_star: float,
    rho_range: tuple[float, float] = RHO_RANGE,
    n_rho: int = N_RHO,
    t_range: tuple[float, float] = T_RANGE,
    n_t: int = N_T,
) -> dict[str, object]:
    """Build a **frozen-composition** ADR-0007 table: the equilibrium composition at the freeze
    reference state `(rho_star, t_star)` (the bounce's mass-weighted turnaround state) held fixed
    across the whole grid — the sudden-freeze rebound EOS of the frozen-recombination bounding run.
    Agrees exactly with the equilibrium table at the reference state (splice continuity)."""
    y = ew.frozen_composition(rho_star, t_star)
    eos_prov: dict[str, object] = {
        "model": "FROZEN-composition water: equilibrium composition at the freeze state held "
        "fixed (no chemistry); constant chemical energy offset, ideal mixture pressure",
        "species": ["H2O", "H", "O", "H+", "O+", "e-"],
        "freeze_state": {"rho_star": rho_star, "T_star": t_star},
        "fractions_per_formula_unit": {
            "y_h2o": y.y_h2o,
            "y_h": y.y_h,
            "y_o": y.y_o,
            "y_hp": y.y_hp,
            "y_op": y.y_op,
            "y_e": y.y_e,
        },
        "energy_reference": "bound molecular H2O at T->0 = 0 (all e > 0)",
        "status": "sudden-freeze pessimistic bound (frozen-recombination check)",
    }
    return _build_table_frozen_common(y, eos_prov, rho_range, n_rho, t_range, n_t)


def build_table_frozen_h2o(
    rho_range: tuple[float, float] = RHO_RANGE,
    n_rho: int = N_RHO,
    t_range: tuple[float, float] = T_RANGE,
    n_t: int = N_T,
) -> dict[str, object]:
    """Build the **chemistry-free** pure-H2O frozen table (freeze *before* the plate): molecular
    water with dissociation/ionization switched off — the optimistic bracket of the
    frozen-recombination check (no chemical sink at all)."""
    eos_prov: dict[str, object] = {
        "model": "FROZEN pure molecular H2O (no dissociation/ionization): thermal-only water "
        "vapor — the no-chemical-sink bracket",
        "species": ["H2O"],
        "energy_reference": "bound molecular H2O at T->0 = 0 (all e > 0)",
        "status": "freeze-before-the-plate optimistic bound (frozen-recombination check)",
    }
    return _build_table_frozen_common(ew.PURE_H2O_FROZEN, eos_prov, rho_range, n_rho, t_range, n_t)


def build_frozen_tables_from_probe(probe_path: Path, outdir: Path) -> list[Path]:
    """Read the Rust `--frozen-probe` JSONL (one `{v, rho_impact, rho_star, t_star}` row per case)
    and emit one frozen-composition table per case into `outdir`, plus the shared pure-H2O
    `h2o.json`. Returns the written paths (per-case tables in probe order, then `h2o.json`)."""
    outdir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    with probe_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            table = build_table_frozen(float(row["rho_star"]), float(row["t_star"]))
            out = outdir / frozen_table_name(float(row["v"]), float(row["rho_impact"]))
            with out.open("w") as out_fh:
                json.dump(table, out_fh)
            written.append(out)

    h2o_out = outdir / "h2o.json"
    with h2o_out.open("w") as out_fh:
        json.dump(build_table_frozen_h2o(), out_fh)
    written.append(h2o_out)
    return written


def build_table_lowv(
    rho_range: tuple[float, float] = RHO_RANGE_LOWV,
    n_rho: int = N_RHO_LOWV,
    t_range: tuple[float, float] = T_RANGE_LOWV,
    n_t: int = N_T_LOWV,
    k_gas_scale: float = 1.0,
) -> dict[str, object]:
    """Build the low-v (Rung C / B-flux) ADR-0007 table from the CoolProp two-phase water EOS.

    Same JSON shape as `build_table` plus a `liquid_frac` field (the condensed mass fraction the
    wall-sticking sink reads, C3) and a `k_gas` field (the gas thermal conductivity the B-flux
    conduction operator reads, ADR-0005). `k_gas_scale != 1` rescales only `k_gas` (the 0.1x/10x
    conduction sensitivity scan). Opacities are transparent (radiation off at 3.2 km/s, design §3).
    CoolProp is imported lazily so this module still imports without the `sci` extra (the high-v
    `build_table` path is CoolProp-free)."""
    from puffsat import eos_cool as ec

    rho_grid = np.geomspace(rho_range[0], rho_range[1], n_rho)
    t_grid = np.geomspace(t_range[0], t_range[1], n_t)

    p, e, cs, liquid_frac, k_gas = ec.eos_grid_lowv(rho_grid, t_grid)
    transparent = [KAPPA_TRANSPARENT] * (n_rho * n_t)

    return {
        "rho_grid": [float(x) for x in rho_grid],
        "T_grid": [float(x) for x in t_grid],
        "shape": [n_rho, n_t],
        "fields": {
            "p": _flatten(p),
            "e": _flatten(e),
            "c_s": _flatten(cs),
            "kappa_rosseland": list(transparent),
            "kappa_planck": list(transparent),
            "liquid_frac": _flatten(liquid_frac),
            "k_gas": _flatten(k_gas * k_gas_scale),
        },
        "provenance": {
            "schema": "ADR-0007",
            "generated_by": "puffsat.tables.build_table_lowv",
            "grid": {
                "rho_range": list(rho_range),
                "n_rho": n_rho,
                "T_range": list(t_range),
                "n_T": n_t,
                "spacing": "geometric (log-spaced) in both axes",
                "units": "rho kg/m^3, T K, p Pa, e J/kg, c_s m/s, liquid_frac [0,1], k_gas W/m/K",
            },
            "eos": {
                "model": "real-fluid equilibrium water across the dome (CoolProp/IAPWS95)",
                "two_phase": "p -> p_sat(T); latent heat folded into e (bulk channel, ADR-0004)",
                "liquid_frac": "condensed mass fraction for the C3 wall-sticking sink (channel 3)",
                "c_s": "sqrt((dp/drho)_s) via a (D,S) finite difference (two-phase-safe)",
                "k_gas": "gas thermal conductivity for the B-flux conduction operator (ADR-0005); "
                f"CoolProp/IAPWS transport, dome -> saturated vapor; k_gas_scale={k_gas_scale}",
            },
            "opacity": {
                "status": "TRANSPARENT placeholder — radiation is off at 3.2 km/s (design §3); "
                "kept positive only for the loader",
            },
        },
    }


def main() -> None:
    """Generate the water table (Makefile `tables`/`tables-lowv`). `--lowv` builds the Rung C
    cool-gas two-phase table; `--out PATH`/`--kappa-scale S` let the B5d-3 scan emit variants."""
    parser = argparse.ArgumentParser(description="Generate the ADR-0007 water EOS/opacity table.")
    parser.add_argument("--out", type=Path, default=None, help="output JSON path")
    parser.add_argument(
        "--kappa-scale", type=float, default=1.0, help="multiply the bracketing opacity by this"
    )
    parser.add_argument(
        "--k-gas-scale", type=float, default=1.0, help="multiply the low-v gas conductivity by this"
    )
    parser.add_argument(
        "--lowv", action="store_true", help="build the low-v (Rung C) cool-gas two-phase table"
    )
    parser.add_argument(
        "--jupiter",
        action="store_true",
        help="build the Jupiter-retrograde 69 km/s scenario table (extended rho/T grid)",
    )
    parser.add_argument(
        "--tops",
        type=Path,
        default=None,
        metavar="TOPS_HTML",
        help="with --jupiter: overlay the real TOPS/OPLIB gray means from this saved pull "
        "(puffsat.fetch_tops) for T >= the OPLIB floor",
    )
    parser.add_argument(
        "--frozen-from-probe",
        type=Path,
        default=None,
        metavar="PROBE_JSONL",
        help="build the per-case frozen-composition tables (+ pure-H2O h2o.json) from the Rust "
        "--frozen-probe output, into --outdir",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_TABLE_DIR_FROZEN,
        help="output directory for --frozen-from-probe",
    )
    args = parser.parse_args()

    if args.frozen_from_probe is not None:
        written = build_frozen_tables_from_probe(args.frozen_from_probe, args.outdir)
        print(f"python: wrote {len(written)} frozen tables -> {args.outdir}")
        return

    if args.lowv:
        table = build_table_lowv(k_gas_scale=args.k_gas_scale)
        out: Path = args.out or DEFAULT_TABLE_PATH_LOWV
        label = (
            f"cool-gas two-phase table -> {out} "
            f"({N_RHO_LOWV}x{N_T_LOWV} nodes, k_gas_scale={args.k_gas_scale})"
        )
    elif args.jupiter:
        table = build_table_jupiter(kappa_scale=args.kappa_scale, tops_path=args.tops)
        out = args.out or DEFAULT_TABLE_PATH_JUPITER
        opacity_label = "TOPS/OPLIB real opacity" if args.tops else "interim opacity"
        label = (
            f"Jupiter 69 km/s scenario table -> {out} "
            f"({N_RHO_JUPITER}x{N_T_JUPITER} nodes, {opacity_label}, "
            f"kappa_scale={args.kappa_scale})"
        )
    else:
        table = build_table(kappa_scale=args.kappa_scale)
        out = args.out or DEFAULT_TABLE_PATH
        label = f"EOS/opacity table -> {out} ({N_RHO}x{N_T} nodes, kappa_scale={args.kappa_scale})"

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        json.dump(table, fh)
    print(f"python: wrote water {label}")


if __name__ == "__main__":
    main()
