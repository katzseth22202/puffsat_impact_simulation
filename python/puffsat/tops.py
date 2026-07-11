"""Parse a TOPS gray-opacity pull and stitch it into a table's `(rho, T)` opacity grid.

TOPS (aphysics2.lanl.gov) computes gray Rosseland/Planck means for a mixture from the LANL
ATOMIC/OPLIB elemental opacities — the *real* plasma opacity ADR-0007 defers to. The raw pull is
the results HTML of the two-stage web-form submission (see `puffsat.fetch_tops`), kept verbatim as
the provenance artifact; this module parses it and interpolates it (log-log bilinear) onto a
target table grid.

TOPS covers `T >= ~5802 K` (0.0005 keV, the OPLIB floor) — atomic/ionic physics only, no molecular
bands. Below the floor the caller's interim opacity is kept: for the 69 km/s scenario that regime
is the late re-expansion tail where radiation no longer moves the bounce, so the stitch seam is
not radiatively active (the sweep's `kappa_scale` bracket still covers residual model error).

Units in the HTML: T in keV, rho in g/cc, kappa in cm^2/g (SI here: K, kg/m^3, m^2/kg).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from puffsat.eos_water import Vec

KEV_TO_KELVIN = 1.16045e7
GCC_TO_SI = 1.0e3  # g/cc -> kg/m^3
CM2G_TO_SI = 0.1  # cm^2/g -> m^2/kg

_BLOCK_HEADER = re.compile(r"T=\s*([0-9.E+-]+)")
_ROW = re.compile(r"^\s*([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)\s+")


@dataclass(frozen=True)
class TopsGray:
    """One TOPS gray pull on its native grid: `kappa[i_rho, j_T]` [m^2/kg], grids in SI."""

    rho_grid: Vec  # [kg/m^3], ascending
    t_grid: Vec  # [K], ascending
    kappa_rosseland: Vec  # (n_rho, n_T)
    kappa_planck: Vec  # (n_rho, n_T)


def parse_tops_gray(html: str) -> TopsGray:
    """Parse the TOPS results page (datype=gray) into SI arrays.

    The page lists, per temperature block (`... T= <keV>`), rows of
    `rho  ross  planck  n_free  av_sq_free`. Tags and `&nbsp;` entities are stripped first.
    """
    text = re.sub(r"<[^>]+>", "", html).replace("&nbsp;", " ")
    lines = text.splitlines()

    t_kev: list[float] = []
    per_t_rho: list[list[float]] = []
    per_t_ross: list[list[float]] = []
    per_t_planck: list[list[float]] = []
    in_block = False
    for line in lines:
        header = _BLOCK_HEADER.search(line)
        if header is not None and "Density" in line:
            t_kev.append(float(header.group(1)))
            per_t_rho.append([])
            per_t_ross.append([])
            per_t_planck.append([])
            in_block = True
            continue
        if in_block:
            row = _ROW.match(line)
            if row is None:
                in_block = False
                continue
            per_t_rho[-1].append(float(row.group(1)))
            per_t_ross[-1].append(float(row.group(2)))
            per_t_planck[-1].append(float(row.group(3)))

    if not t_kev:
        raise ValueError("no TOPS temperature blocks found — not a gray results page?")
    rho_gcc = per_t_rho[0]
    if any(block != rho_gcc for block in per_t_rho[1:]):
        raise ValueError("TOPS density grid differs between temperature blocks")

    # native layout is per-T rows of rho -> transpose to kappa[i_rho, j_T]
    ross = np.array(per_t_ross, dtype=np.float64).T * CM2G_TO_SI
    planck = np.array(per_t_planck, dtype=np.float64).T * CM2G_TO_SI
    t_kelvin = np.array(t_kev, dtype=np.float64) * KEV_TO_KELVIN
    rho_si = np.array(rho_gcc, dtype=np.float64) * GCC_TO_SI
    if not (np.all(np.diff(t_kelvin) > 0) and np.all(np.diff(rho_si) > 0)):
        raise ValueError("TOPS grids are not strictly ascending")
    if not (np.all(ross > 0) and np.all(planck > 0)):
        raise ValueError("non-positive TOPS opacity — the loader interpolates ln kappa")
    return TopsGray(rho_grid=rho_si, t_grid=t_kelvin, kappa_rosseland=ross, kappa_planck=planck)


def load_tops_gray(path: Path) -> TopsGray:
    """Read and parse a saved TOPS results page."""
    return parse_tops_gray(path.read_text())


def _interp_loglog(tops_grid: Vec, target: Vec, ln_field: Vec, axis: int) -> Vec:
    """Linear interpolation of `ln_field` along `axis` in `ln(grid)`, clamped at the grid ends."""
    x = np.log(tops_grid)
    xq = np.clip(np.log(target), x[0], x[-1])
    idx = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    w = (xq - x[idx]) / (x[idx + 1] - x[idx])
    lo = np.take(ln_field, idx, axis=axis)
    hi = np.take(ln_field, idx + 1, axis=axis)
    shape = [1, 1]
    shape[axis] = len(target)
    w = w.reshape(shape)
    return np.asarray(lo * (1.0 - w) + hi * w, dtype=np.float64)


def stitch_opacity(
    rho_grid: Vec,
    t_grid: Vec,
    kappa_r_interim: Vec,
    kappa_p_interim: Vec,
    tops: TopsGray,
) -> tuple[Vec, Vec]:
    """Overlay TOPS means onto the interim `(n_rho, n_T)` opacity for `T >=` the TOPS floor.

    Log-log bilinear interpolation (matching the Rust loader's convention), clamped at the TOPS
    grid edges in rho. Below the TOPS temperature floor the interim values pass through unchanged.
    Returns new arrays; the inputs are not modified.
    """
    hot = t_grid >= tops.t_grid[0]
    kappa_r = kappa_r_interim.copy()
    kappa_p = kappa_p_interim.copy()
    for interim, field in ((kappa_r, tops.kappa_rosseland), (kappa_p, tops.kappa_planck)):
        ln_on_rho = _interp_loglog(tops.rho_grid, rho_grid, np.log(field), axis=0)
        ln_full = _interp_loglog(tops.t_grid, t_grid[hot], ln_on_rho, axis=1)
        interim[:, hot] = np.exp(ln_full)
    return kappa_r, kappa_p
