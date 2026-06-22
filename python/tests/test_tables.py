"""Contract tests for the assembled ADR-0007 table (B5c-2).

The Rust loader (`crates/tables`) requires: `shape == [len(rho_grid), len(T_grid)]`, strictly
ascending positive grids, each field of length `n_rho * n_T`, and every field value `> 0` (it
interpolates the logs). These tests replicate that contract on the Python side, and check the
interim bracketing opacity puts the column optical depth `tau = rho*kappa_R*L` in the design's
`[1e2, 1e5]` band at the nominal stagnation, with the right scaling. (The *real* Rust `Table::load`
round-trip is exercised by the B5d-1 sweep, which loads the generated file.)
"""

from __future__ import annotations

from itertools import pairwise
from typing import cast

import numpy as np

from puffsat import tables


def test_table_satisfies_loader_contract() -> None:
    """The emitted dict obeys the Rust loader's structural + positivity invariants (ADR-0007)."""
    # Small grid: the EOS is a per-node Newton solve, so keep the test fast.
    table = tables.build_table(rho_range=(0.05, 5.0), n_rho=6, t_range=(300.0, 50_000.0), n_t=8)

    shape = cast("list[int]", table["shape"])
    n_rho, n_t = shape
    rho_grid = cast("list[float]", table["rho_grid"])
    t_grid = cast("list[float]", table["T_grid"])
    fields = cast("dict[str, list[float]]", table["fields"])

    # shape agrees with the grids, and the grids are strictly ascending and positive.
    assert (n_rho, n_t) == (6, 8)
    assert len(rho_grid) == n_rho
    assert len(t_grid) == n_t
    assert rho_grid[0] > 0.0
    assert t_grid[0] > 0.0
    assert all(b > a for a, b in pairwise(rho_grid))
    assert all(b > a for a, b in pairwise(t_grid))

    # every field is fully populated and strictly positive (log interpolation requires > 0).
    assert set(fields) == {"p", "e", "c_s", "kappa_rosseland", "kappa_planck"}
    for name, vals in fields.items():
        assert len(vals) == n_rho * n_t, name
        assert all(v > 0.0 for v in vals), name


def test_interim_opacity_tau_in_design_band() -> None:
    """At the nominal stagnation the column optical depth hits TAU_TARGET, inside the band."""
    rho = np.array([tables.KAPPA_RHO_REF])
    temp = np.array([tables.KAPPA_T_REF])
    kappa_r, kappa_p = tables.opacity_grid(rho, temp)

    tau = float(rho[0] * kappa_r[0, 0] * tables.KAPPA_L_REF)
    lo, hi = tables.TAU_BAND
    np.testing.assert_allclose(tau, tables.TAU_TARGET, rtol=1e-12)
    assert lo <= tau <= hi

    # Planck mean exceeds Rosseland mean (re-emission), and both are positive.
    assert kappa_p[0, 0] == tables.PLANCK_OVER_ROSSELAND * kappa_r[0, 0]
    assert kappa_p[0, 0] > kappa_r[0, 0] > 0.0


def test_opacity_scales_as_kramers() -> None:
    """`kappa_R` is linear in rho and follows the Kramers `T^-3.5` temperature law."""
    rho = np.array([1.0, 2.0])
    temp = np.array([tables.KAPPA_T_REF, 2.0 * tables.KAPPA_T_REF])
    kappa_r, _ = tables.opacity_grid(rho, temp)

    # linear in rho at fixed T: doubling rho doubles kappa_R.
    np.testing.assert_allclose(kappa_r[1, 0], 2.0 * kappa_r[0, 0], rtol=1e-12)
    # T^-3.5 at fixed rho: doubling T scales kappa_R by 2^-3.5.
    np.testing.assert_allclose(kappa_r[0, 1], kappa_r[0, 0] * 2.0**-3.5, rtol=1e-12)
