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
import pytest

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


def test_lowv_table_has_liquid_frac_and_loader_contract() -> None:
    """The Rung C low-v table (CoolProp) carries a `liquid_frac` field; the five log-fields stay
    positive and `liquid_frac` is in `[0,1]`. (Skipped without the `sci` extra / CoolProp.)"""
    pytest.importorskip("CoolProp")
    # Small grid for speed (CoolProp is a per-node query).
    table = tables.build_table_lowv(rho_range=(0.05, 50.0), n_rho=5, t_range=(300.0, 1700.0), n_t=8)

    n_rho, n_t = cast("list[int]", table["shape"])
    fields = cast("dict[str, list[float]]", table["fields"])
    assert set(fields) == {"p", "e", "c_s", "kappa_rosseland", "kappa_planck", "liquid_frac"}
    for name in ("p", "e", "c_s", "kappa_rosseland", "kappa_planck"):
        assert len(fields[name]) == n_rho * n_t, name
        assert all(v > 0.0 for v in fields[name]), name
    lf = fields["liquid_frac"]
    assert len(lf) == n_rho * n_t
    assert all(0.0 <= v <= 1.0 for v in lf)


def test_kappa_scale_rescales_opacity_only() -> None:
    """`kappa_scale` multiplies both opacity means by the factor and leaves p, e, c_s untouched.

    This is the knob the B5d-3 insensitivity scan turns: only the opacity moves, so any change in
    `e_eff` it produces is attributable to opacity alone.
    """
    rho = np.array([0.2, 1.0])
    temp = np.array([1.0e4, 2.0e4])
    kr1, kp1 = tables.opacity_grid(rho, temp, 1.0)
    kr10, kp10 = tables.opacity_grid(rho, temp, 10.0)
    np.testing.assert_allclose(kr10, 10.0 * kr1, rtol=1e-12)
    np.testing.assert_allclose(kp10, 10.0 * kp1, rtol=1e-12)

    # The EOS fields in a scaled table are identical to the unscaled ones.
    g = ((0.05, 5.0), 5, (300.0, 50_000.0), 6)  # rho_range, n_rho, t_range, n_t
    base = tables.build_table(g[0], g[1], g[2], g[3], kappa_scale=1.0)
    scaled = tables.build_table(g[0], g[1], g[2], g[3], kappa_scale=10.0)
    base_fields = base["fields"]
    scaled_fields = scaled["fields"]
    assert isinstance(base_fields, dict)
    assert isinstance(scaled_fields, dict)
    for name in ("p", "e", "c_s"):
        assert scaled_fields[name] == base_fields[name]
    np.testing.assert_allclose(
        scaled_fields["kappa_rosseland"],
        10.0 * np.array(base_fields["kappa_rosseland"]),
        rtol=1e-12,
    )
