"""Tests for the skirt's hoop response (`skirt_hoop`), ADR-0055 Q56.

The skirt is a thin ring loaded by the gas pressure on its inner face. Its breathing mode is a
single oscillator, `m'' x'' + (E t / R^2) x = p(t)`, so it has two closed-form limits (the seam
agreed before these tests were written):

* a pulse much shorter than the ring's period delivers an impulse `I`, and the peak hoop strain is
  `I / (m'' c)` with `c = sqrt(E / rho)` for one material;
* a pressure applied suddenly and held is a step, and the peak strain is twice the static
  `p R / (E t)`.
"""

from __future__ import annotations

import numpy as np
import pytest

from puffsat.water_plate import skirt_hoop as sh

ARAMID = sh.Material("aramid", modulus=112.0e9, density=1440.0)


def test_a_short_pulse_gives_the_impulsive_strain() -> None:
    # 1,000 Pa s in 1 us on 1 cm of aramid at R = 10 m. Hand-worked: c = sqrt(112e9 / 1440) =
    # 8819.2 m/s, m'' = 14.4 kg/m^2, so the strain is 1000 / (14.4 * 8819.2) = 7.874e-3.
    time = np.array([0.0, 1.0e-6, 2.0e-6])
    pressure = np.array([1.0e9, 0.0, 0.0])
    strain = sh.peak_hoop_strain(time, pressure, 10.0, [sh.Layer(ARAMID, 0.01)])
    assert strain == pytest.approx(7.874e-3, rel=1e-2)


def test_a_held_pressure_gives_twice_the_static_strain() -> None:
    # 1 MPa held for many periods (the period is 2 pi R / c ~ 7 ms) on 1 cm of aramid at R = 10 m.
    # Hand-worked: 2 p R / (E t) = 2 * 1e6 * 10 / (112e9 * 0.01) = 1.786e-2.
    time = np.linspace(0.0, 0.2, 20001)
    pressure = np.full_like(time, 1.0e6)
    strain = sh.peak_hoop_strain(time, pressure, 10.0, [sh.Layer(ARAMID, 0.01)])
    assert strain == pytest.approx(1.786e-2, rel=5e-3)
