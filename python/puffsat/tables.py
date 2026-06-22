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

import json
from pathlib import Path

import numpy as np

from puffsat import eos_water as ew
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


def opacity_grid(rho_grid: Vec, t_grid: Vec) -> tuple[Vec, Vec]:
    """Bracketing `(kappa_R, kappa_P)` [m^2/kg] on the `(rho, T)` grid, shaped `(n_rho, n_T)`.

    `kappa_R = KAPPA0_R (rho/rho_ref) (T/T_ref)^-3.5`; `kappa_P = PLANCK_OVER_ROSSELAND * kappa_R`.
    Positive-definite (the loader interpolates `ln kappa`). PROVISIONAL — see the module docstring.
    """
    kappa_r = (
        KAPPA0_R * (rho_grid[:, None] / KAPPA_RHO_REF) * (t_grid[None, :] / KAPPA_T_REF) ** -3.5
    )
    kappa_p = PLANCK_OVER_ROSSELAND * kappa_r
    return kappa_r, kappa_p


def _flatten(a: Vec) -> list[float]:
    """Row-major (C-order) flatten to a plain `list[float]` for JSON (no numpy scalars leak)."""
    return [float(x) for x in a.reshape(-1)]


def _provenance(
    rho_range: tuple[float, float], n_rho: int, t_range: tuple[float, float], n_t: int
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
) -> dict[str, object]:
    """Build the ADR-0007 table dict: log-spaced grids, the equilibrium EOS fields, the bracketing
    opacity, and nested provenance. Fields are flattened row-major over `(rho, T)`."""
    rho_grid = np.geomspace(rho_range[0], rho_range[1], n_rho)
    t_grid = np.geomspace(t_range[0], t_range[1], n_t)

    p, e, cs = ew.eos_grid(rho_grid, t_grid)
    kappa_r, kappa_p = opacity_grid(rho_grid, t_grid)

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
        "provenance": _provenance(rho_range, n_rho, t_range, n_t),
    }


def main() -> None:
    """Generate the production water table at `DEFAULT_TABLE_PATH` (Makefile `tables` target)."""
    table = build_table()
    DEFAULT_TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with DEFAULT_TABLE_PATH.open("w") as fh:
        json.dump(table, fh)
    print(f"python: wrote water EOS/opacity table -> {DEFAULT_TABLE_PATH} ({N_RHO}x{N_T} nodes)")


if __name__ == "__main__":
    main()
