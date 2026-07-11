"""Tests for the TOPS gray-pull parser and the opacity stitch (puffsat.tops)."""

from __future__ import annotations

import numpy as np
import pytest

from puffsat import tops

# A minimal synthetic results page in the real TOPS layout: two temperature blocks over a
# three-point density grid, tags + &nbsp; entities included as the site emits them.
_HTML = """
<html><body><pre>
Rosseland&nbsp;and&nbsp;Planck&nbsp;opacities&nbsp;and&nbsp;free&nbsp;electrons
&nbsp;Density&nbsp;&nbsp;Ross&nbsp;opa&nbsp;&nbsp;Planck&nbsp;opa&nbsp;&nbsp;No.&nbsp;Free&nbsp;&nbsp;Av&nbsp;Sq&nbsp;Free&nbsp;&nbsp;T=&nbsp;&nbsp;1.0000E-03
&nbsp;&nbsp;1.0000E-07&nbsp;&nbsp;2.0000E-02&nbsp;&nbsp;3.0000E+01&nbsp;&nbsp;4.0000E-04&nbsp;&nbsp;4.0000E-04
&nbsp;&nbsp;1.0000E-05&nbsp;&nbsp;4.0000E-02&nbsp;&nbsp;5.0000E+01&nbsp;&nbsp;3.0000E-04&nbsp;&nbsp;3.0000E-04
&nbsp;&nbsp;1.0000E-03&nbsp;&nbsp;8.0000E-02&nbsp;&nbsp;7.0000E+01&nbsp;&nbsp;2.0000E-04&nbsp;&nbsp;2.0000E-04
&nbsp;Density&nbsp;&nbsp;Ross&nbsp;opa&nbsp;&nbsp;Planck&nbsp;opa&nbsp;&nbsp;No.&nbsp;Free&nbsp;&nbsp;Av&nbsp;Sq&nbsp;Free&nbsp;&nbsp;T=&nbsp;&nbsp;1.0000E-02
&nbsp;&nbsp;1.0000E-07&nbsp;&nbsp;1.0000E-03&nbsp;&nbsp;2.0000E+00&nbsp;&nbsp;1.0000E+00&nbsp;&nbsp;1.0000E+00
&nbsp;&nbsp;1.0000E-05&nbsp;&nbsp;3.0000E-03&nbsp;&nbsp;4.0000E+00&nbsp;&nbsp;1.0000E+00&nbsp;&nbsp;1.0000E+00
&nbsp;&nbsp;1.0000E-03&nbsp;&nbsp;9.0000E-03&nbsp;&nbsp;8.0000E+00&nbsp;&nbsp;1.0000E+00&nbsp;&nbsp;1.0000E+00
</pre></body></html>
"""


def test_parse_shapes_and_units() -> None:
    """Grids convert keV->K and g/cc->kg/m^3; kappa converts cm^2/g->m^2/kg, laid out (rho, T)."""
    pull = tops.parse_tops_gray(_HTML)
    assert pull.t_grid == pytest.approx([1.16045e4, 1.16045e5])
    assert pull.rho_grid == pytest.approx([1.0e-4, 1.0e-2, 1.0])
    assert pull.kappa_rosseland.shape == (3, 2)
    # first block, first row: 2e-2 cm^2/g = 2e-3 m^2/kg at (rho0, T0)
    assert pull.kappa_rosseland[0, 0] == pytest.approx(2.0e-3)
    assert pull.kappa_planck[0, 0] == pytest.approx(3.0)
    # second block, last row lands at (rho2, T1)
    assert pull.kappa_rosseland[2, 1] == pytest.approx(9.0e-4)


def test_parse_rejects_non_gray_page() -> None:
    """A page without temperature blocks (e.g. the request form) is a loud error."""
    with pytest.raises(ValueError, match="no TOPS temperature blocks"):
        tops.parse_tops_gray("<html>TOPS Opacities form</html>")


def test_stitch_overlays_hot_keeps_cold() -> None:
    """Above the TOPS floor the table takes TOPS values (exact at TOPS nodes); below, interim."""
    pull = tops.parse_tops_gray(_HTML)
    # target grid: one T below the TOPS floor, two exactly at TOPS nodes; rho at TOPS nodes
    rho_grid = np.array([1.0e-4, 1.0e-2])
    t_grid = np.array([300.0, 1.16045e4, 1.16045e5])
    interim = np.full((2, 3), 7.0)
    kappa_r, kappa_p = tops.stitch_opacity(rho_grid, t_grid, interim, interim, pull)
    assert kappa_r[:, 0] == pytest.approx([7.0, 7.0])  # cold column untouched
    assert kappa_r[0, 1] == pytest.approx(2.0e-3)  # exact TOPS node reproduced
    assert kappa_r[1, 2] == pytest.approx(3.0e-4)
    assert kappa_p[1, 1] == pytest.approx(5.0)
    assert interim[0, 1] == pytest.approx(7.0)  # inputs not mutated


def test_stitch_interpolates_loglog() -> None:
    """Between TOPS nodes the stitch is linear in (ln rho, ln T, ln kappa)."""
    pull = tops.parse_tops_gray(_HTML)
    rho_grid = np.array([1.0e-3])  # geometric midpoint of 1e-4 and 1e-2
    t_grid = np.array([1.16045e4])
    interim = np.full((1, 1), 7.0)
    kappa_r, _ = tops.stitch_opacity(rho_grid, t_grid, interim, interim, pull)
    assert kappa_r[0, 0] == pytest.approx(np.sqrt(2.0e-3 * 4.0e-3))


def test_stitch_clamps_rho_outside_tops_grid() -> None:
    """rho beyond the TOPS grid clamps to the edge value instead of extrapolating."""
    pull = tops.parse_tops_gray(_HTML)
    rho_grid = np.array([1.0e-6, 30.0])  # below/above the TOPS rho range
    t_grid = np.array([1.16045e4])
    interim = np.full((2, 1), 7.0)
    kappa_r, _ = tops.stitch_opacity(rho_grid, t_grid, interim, interim, pull)
    assert kappa_r[0, 0] == pytest.approx(2.0e-3)
    assert kappa_r[1, 0] == pytest.approx(8.0e-3)
