"""Tests for the seeded-plasma conductivity `sigma(T, rho, x_K)` (Study 1).

The paper's `tab:seed_window` reports a magnetic Reynolds number `Rm` against temperature, and the
leak schedule rests on where `Rm` falls through 1. `Rm = mu0 sigma v L`, so the whole table is a
conductivity calculation, and conductivity is set by the electron density -- which has two sources
that swap dominance. Below ~5000 K the potassium seed supplies essentially all the electrons; above
it, water's own ionisation does, by ~38x at 15 000 K. Getting only one of them is how the first hand
estimate went wrong.
"""

from __future__ import annotations

from itertools import pairwise

import pytest

from puffsat import conductivity


def test_electron_density_follows_the_saha_law_where_the_seed_dominates() -> None:
    """Slice 1: the seed-dominated limit, against the closed-form Saha solution.

    At 2500 K water is effectively un-ionised (`eos_water` gives `n_e ~ 6e10 m^-3`), so the electron
    density is potassium's alone and has an exact closed form. For `K = n_e n_K+ / n_K0` with
    `n_K+ = n_e` and `n_K0 = n_K - n_e`, `n_e = (-K + sqrt(K^2 + 4 K n_K))/2`.

    Worked independently of the implementation, at `rho = 0.32 kg/m^3`, `x_K = 0.01` (mass
    fraction), `chi_K = 4.34 eV`, `g+/g0 = 1/2`:
        n_K = x_K rho / m_K                        = 4.92882e22 m^-3
        K(2500 K) = (2 g+/g0)(2 pi m_e k T/h^2)^1.5 exp(-chi/kT)
                                                   = 5.3792e17 m^-3
        n_e                                        = 1.62559e20 m^-3
    """
    n_e = conductivity.electron_density(2500.0, 0.32, 0.01)
    assert n_e == pytest.approx(1.62559e20, rel=1e-3)


def test_seed_saturates_and_water_takes_over_at_high_temperature() -> None:
    """Slice 1b: the other limit, and the error the first hand estimate made.

    At 15 000 K the seed is **fully ionised** -- the weak-ionisation `sqrt(n_K K)` law would give
    2.76e24, well above the 4.93e22 potassium atoms available, so the seed saturates at `n_K`.
    Water's own ionisation then dominates, supplying ~38x more electrons than the seed. A model
    that counts only the seed is low by that factor, which is what the audit caught.
    """
    n_e = conductivity.electron_density(15000.0, 0.32, 0.01)
    seed_ceiling = 4.92882e22
    assert n_e > 20.0 * seed_ceiling, "water's own ionisation must dominate here"
    assert n_e == pytest.approx(1.8624e24 + 4.9272e22, rel=0.02)


def test_conductivity_reaches_the_spitzer_limit_when_neutrals_are_switched_off() -> None:
    """Slice 2: the analytic acceptance test. With electron-neutral collisions off, the conductivity
    must be Spitzer's.

    Spitzer conductivity is famously **independent of electron density** except through `ln Lambda`
    (more electrons carry more current but also scatter it more, and the two cancel). The standard
    result, from the NRL Plasma Formulary's `tau_e = 3.44e5 T[eV]^1.5/(n_e[cm^-3] ln Lambda)` with
    the `1.96` electron-electron correction for `Z = 1`:

        sigma_Spitzer = 1.899e4 * T_e[eV]^1.5 / ln Lambda   [S/m]

    Evaluated independently at `T = 15 000 K`, `rho = 0.32`, `x_K = 0.01`:
        n_e        = 1.9117e24 m^-3        (water 1.8624e24 + seed 4.927e22)
        lambda_D   = sqrt(eps0 k T / n_e e^2)      = 6.1127e-9 m
        Lambda     = 12 pi n_e lambda_D^3          = 16.460
        ln Lambda  = 2.8010
        T_e        = 1.29260 eV,  T_e^1.5          = 1.46955
        sigma      = 1.899e4 * 1.46955 / 2.8010    = 9963 S/m

    `ln Lambda = 2.8` is **small**, and Spitzer's derivation assumes it is large. That is a real
    caveat on this whole model, recorded rather than hidden: see `sigma`'s docstring.
    """
    s = conductivity.sigma(15000.0, 0.32, 0.01, q_en=0.0)
    assert s == pytest.approx(9963.0, rel=0.02)


def test_blended_conductivity_reproduces_the_audit_hand_calculation() -> None:
    """Slice 3: both scattering channels together, against an independent hand calculation.

    The 2026-08-21 audit re-did the 15 000 K point by hand after finding that the first estimate had
    omitted water's own ionisation. Its corrected blend is **~6950 S/m**, against the ~15 900 S/m
    that `tab:seed_window`'s `Rm = 361` implies -- a gap of 2.3x, not the 28x first reported.

    This is a genuinely independent source of truth: a different person, working the same physics by
    hand, with the electron density, both collision channels and the Spitzer correction assembled
    separately from this module. Reproducing it to a few percent is the check that the assembly
    here is not quietly missing a term.
    """
    s = conductivity.sigma(15000.0, 0.32, 0.01)
    assert s == pytest.approx(6950.0, rel=0.05)


def test_magnetic_reynolds_is_the_diffusion_to_expansion_time_ratio() -> None:
    """Slice 4: `Rm = mu0 sigma v L`, and what it means.

    Item 10's leak argument is that the ratio of field-diffusion time to expansion time is
    `tau_d/t_exp = (mu0 sigma L^2)/(L/v) = mu0 sigma v L`, which *is* the magnetic Reynolds
    number -- so the leak fraction goes as `1/Rm` and no MHD code is needed to bracket it.

    Worked by hand at `sigma = 6950 S/m`, `v L = 100 m^2/s`:
        Rm = 4 pi x 1e-7 * 6950 * 100 = 0.8734
    """
    rm = conductivity.magnetic_reynolds(6950.0, v_l=100.0)
    assert rm == pytest.approx(0.8734, rel=1e-3)
    # Linear in the product, which is why `v` and `L` only ever enter together.
    assert conductivity.magnetic_reynolds(6950.0, v_l=200.0) == pytest.approx(2.0 * rm, rel=1e-12)


def test_the_cliff_temperature_depends_on_a_product_the_paper_never_states() -> None:
    """Slice 4b: where `Rm` falls through 1 -- and why that number cannot be quoted alone.

    `Rm = mu0 sigma(T) v L`, so the cliff temperature is set jointly by the conductivity *and* the
    product `v L`. The paper's `tab:seed_window` gives an `Rm` column with **no stated `v` and no
    stated `L`** (audit Q-G), so the cliff it implies is not reproducible from the paper alone. It
    is reported here as a function of `v L` instead of collapsed to one number.

    The only claim made is monotonicity: a larger `v L` holds the field longer, so the cliff moves
    to lower temperature."""
    cold = conductivity.cliff_temperature(0.32, 0.01, v_l=1.0e3)
    warm = conductivity.cliff_temperature(0.32, 0.01, v_l=2.0e4)
    assert warm < cold, "a larger v*L must push the cliff to lower temperature"
    for t in (cold, warm):
        assert 1500.0 < t < 15000.0


def test_no_cliff_is_reported_rather_than_extrapolated_when_rm_never_reaches_one() -> None:
    """A small `v L` never holds the field at any temperature in the window, so there is no cliff to
    report. Returning the range edge would look like an answer; raising says what is true.

    At `v L = 10 m^2/s`, `Rm` tops out near 0.09 even at 15 000 K."""
    with pytest.raises(ValueError, match="not bracketed"):
        conductivity.cliff_temperature(0.32, 0.01, v_l=10.0)


def test_seed_window_caps_the_leak_where_the_field_is_not_held() -> None:
    """Item 10's bracket is `leak ~ 1/Rm`, which is only meaningful while `Rm > 1`.

    Below the cliff `1/Rm` exceeds 1 and would read as "a leak of 40x" -- reporting the breakdown of
    an approximation as though it were a measurement. The table caps it at 1 and says the field is
    not held. The ionised fraction is capped the same way: the seed cannot supply more electrons
    than it has atoms, which is exactly the saturation the first hand estimate missed.
    """
    rows = conductivity.seed_window(0.32, 0.01, 1.81e4)
    assert [r.temp for r in rows] == list(conductivity.SEED_WINDOW_TEMPS)
    for r in rows:
        assert 0.0 <= r.leak_fraction <= 1.0
        assert 0.0 <= r.ionised_fraction <= 1.0
    # Conductivity rises monotonically with temperature across the window.
    assert all(a.sigma < b.sigma for a, b in pairwise(rows))
    # The seed is fully ionised well before the top of the window.
    assert rows[-1].ionised_fraction == pytest.approx(1.0)
    assert rows[0].leak_fraction == pytest.approx(1.0), "at 2000 K the field is not held at all"
