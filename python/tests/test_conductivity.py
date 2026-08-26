"""Tests for the seeded-plasma conductivity `sigma(T, rho, x_K)` (Study 1).

The paper's `tab:seed_window` reports a magnetic Reynolds number `Rm` against temperature, and the
leak schedule rests on where `Rm` falls through 1. `Rm = mu0 sigma v L`, so the whole table is a
conductivity calculation, and conductivity is set by the electron density -- which has two sources
that swap dominance. Below ~5000 K the potassium seed supplies essentially all the electrons; above
it, water's own ionisation does, by ~38x at 15 000 K. Getting only one of them is how the first hand
estimate went wrong.
"""

from __future__ import annotations

import math
from itertools import pairwise

import pytest

from puffsat import conductivity
from puffsat.eos_water import K_B


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


def test_the_reference_v_l_is_the_solved_expansion_exit_not_a_back_solved_table() -> None:
    """Deferred item D9: `REF_V_L` moved off the retired 1.81e4 back-solve.

    `v L` is an **output** -- `expansion.HistoryRow.v_l` is the local flow speed times the local
    flux-tube radius -- so the value this module reports the seed window at has to be a station of
    the solved expansion, not a number reverse-engineered from the paper's own `Rm` column. The
    stated leg's nozzle exit is 7.44e4, quoted 7.4e4, and the four flown legs span 5.5e4-9.7e4
    between them.

    The retired value is kept, and this test is what stops it drifting back in.
    """
    assert conductivity.REF_V_L == 7.4e4
    assert conductivity.RETIRED_V_L == 1.81e4
    lo, hi = conductivity.V_L_BAND
    assert lo < conductivity.REF_V_L < hi, "the stated leg must sit inside the flown band"
    assert lo > conductivity.RETIRED_V_L, "the retired value is below every flown leg's exit"


def test_the_cliff_at_the_stated_expansion_is_2450_k_and_the_band_barely_moves_it() -> None:
    """Deferred item D9/P1: the number the paper prints, and its sensitivity.

    Two claims, and the second is the one that matters. At the stated `v L` the cliff is 2450 K
    (2449 K at the paper's rounded `rho` = 0.32, 2450 K at the flown 0.323 -- 1 K apart, well below
    anything this model resolves). Across the *whole* flown band it moves only 2386-2524 K, so the
    leak limit at 3800 K binds by ~1300 K no matter which leg is flown.
    """
    stated = conductivity.cliff_temperature(0.32, 0.01, conductivity.REF_V_L)
    assert stated == pytest.approx(2449.0, abs=2.0)
    assert conductivity.cliff_temperature(0.323, 0.01, conductivity.REF_V_L) == pytest.approx(
        2450.0, abs=2.0
    )
    lo, hi = conductivity.V_L_BAND
    cold = conductivity.cliff_temperature(0.32, 0.01, hi)
    warm = conductivity.cliff_temperature(0.32, 0.01, lo)
    assert cold == pytest.approx(2386.0, abs=3.0)
    assert warm == pytest.approx(2524.0, abs=3.0)
    assert warm - cold < 150.0, "the whole band spans less than 150 K of cliff"
    assert 3800.0 - warm > 1250.0, "the leak limit still binds well before the field lets go"


def test_interpolating_the_six_tabulated_conductivities_overstates_the_cliff() -> None:
    """Deferred item D9: *why* 2570 K was wrong, kept executable so it is not re-derived.

    `tab:seed_window` samples every 1000 K while `sigma` climbs 60x between its first two rows, and
    the crossing lies inside that first interval. Log-interpolating the table and solving `Rm = 1`
    on the interpolant therefore guesses high -- by about 120 K at the stated `v L`, which is the
    whole of the 2570-versus-2450 discrepancy.

    The tabulated values are not the problem: they are the same `sigma` the solver uses. Asking six
    points for a crossing inside their first interval is.
    """
    for v_l in (conductivity.RETIRED_V_L, conductivity.REF_V_L, *conductivity.V_L_BAND):
        solved = conductivity.cliff_temperature(0.32, 0.01, v_l)
        guessed = conductivity.interpolated_cliff_temperature(0.32, 0.01, v_l)
        assert guessed > solved, "interpolation must err high, never low"
        assert 40.0 < guessed - solved < 130.0

    stated = conductivity.REF_V_L
    assert conductivity.interpolated_cliff_temperature(0.32, 0.01, stated) == pytest.approx(
        2566.0, abs=3.0
    )
    # The crossing really is inside the table's first interval, which is the reason for all of it.
    assert conductivity.SEED_WINDOW_TEMPS[0] < 2449.0 < conductivity.SEED_WINDOW_TEMPS[1]


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
    rows = conductivity.seed_window(0.32, 0.01, conductivity.REF_V_L)
    assert [r.temp for r in rows] == list(conductivity.SEED_WINDOW_TEMPS)
    for r in rows:
        assert 0.0 <= r.leak_fraction <= 1.0
        assert 0.0 <= r.ionised_fraction <= 1.0
    # Conductivity rises monotonically with temperature across the window.
    assert all(a.sigma < b.sigma for a, b in pairwise(rows))
    # The seed is fully ionised well before the top of the window.
    assert rows[-1].ionised_fraction == pytest.approx(1.0)
    assert rows[0].leak_fraction == pytest.approx(1.0), "at 2000 K the field is not held at all"


def test_hall_parameter_agrees_with_the_conductivity_identity() -> None:
    """Slice 5: `beta = omega_c tau = e B / (m_e nu)`.

    Checked against a route that shares no arithmetic with the implementation. Since
    `sigma = n_e e^2/(m_e nu)`, it follows that `beta = sigma B / (n_e e)` -- the same number
    assembled from three already-tested public quantities instead of from the collision frequency.
    If the two disagree, either `sigma` or `beta` has the collision frequency wrong.
    """
    t, rho, x_k, b = 3000.0, 0.32, 0.01, 2.0
    beta = conductivity.hall_parameter(t, rho, x_k, b)
    identity = (
        conductivity.sigma(t, rho, x_k)
        * b
        / (conductivity.electron_density(t, rho, x_k) * conductivity.E_CHARGE)
    )
    assert beta == pytest.approx(identity, rel=1e-12)
    # Linear in B, by construction.
    assert conductivity.hall_parameter(t, rho, x_k, 4.0) == pytest.approx(2.0 * beta, rel=1e-12)


def test_ionisation_sensitivity_matches_the_saha_exponent_where_the_seed_is_weak() -> None:
    """Slice 6: `S = d ln n_e / d ln T_e`, the gain in the instability's feedback loop.

    The electrothermal runaway is: a local rise in `n_e` concentrates current, which raises local
    heating, which raises `T_e`, which raises `n_e`. `S` is the gain of that last link, so it
    decides whether the loop can close at all.

    In the weakly-ionised seed limit `n_e = sqrt(n_K K)` with `K ~ T^1.5 exp(-chi/kT)`, so
    `n_e ~ T^0.75 exp(-chi/2kT)` and the exponent follows analytically:

        S = 3/4 + chi / (2 k T_e)

    At 2500 K with `chi = 4.34 eV`: `kT = 0.2154333 eV`, so `S = 0.75 + 10.0727 = 10.823`.
    A gain of ~11 is very large -- this is why seeded plasmas are prone to the instability at all.
    """
    s = conductivity.ionisation_sensitivity(2500.0, 0.32, 0.01)
    assert s == pytest.approx(10.823, rel=0.02)


def test_ionisation_sensitivity_collapses_once_the_seed_is_fully_ionised() -> None:
    """Slice 6b: the stabilising limit, and it is a hard one.

    Once every potassium atom is ionised, raising `T_e` cannot liberate more seed electrons -- the
    feedback loop's gain collapses and the ionisation runaway is choked off at source. Water's own
    ionisation takes over, but at 12.6-13.6 eV it is far less sensitive than the 4.34 eV seed at
    these temperatures.

    This is a *computable* stabilisation boundary, not a literature criterion, which is why it is
    the part of the instability story this repository can actually settle."""
    weak = conductivity.ionisation_sensitivity(2500.0, 0.32, 0.01)
    saturated = conductivity.ionisation_sensitivity(8000.0, 0.32, 0.01)
    assert saturated < 0.25 * weak, "the gain must collapse once the seed saturates"


def test_electrothermal_screen_needs_both_hall_drive_and_ionisation_gain() -> None:
    """Slice 7: the screen, and what it is honestly for.

    The electrothermal runaway needs *both* links of its loop intact: a Hall parameter large enough
    for a conductivity perturbation to redirect current and change the local heating (`beta` above
    ~2), and an ionisation gain large enough for the resulting `T_e` rise to make more electrons
    (`S` above ~1). Killing either breaks the loop.

    The regime map is the opposite of the intuitive one, and the numbers below are measured rather
    than assumed. At the bag's 0.32 kg/m^3 the plasma is **strongly collisional**: `beta` is only
    0.42 at 3000 K in a 1 T field, and it *falls* with temperature (0.52 at 2500 K to 0.023 at
    15 000 K) because Coulomb collisions grow faster than the mobility.

    **The hot end is safe on the Hall link alone, not on both** (ADR-0038). An earlier version of
    this test asserted the gain was gone up there too. It is not -- that reading came from
    differentiating only the seed, and once water's ionisation is allowed to respond the gain at
    11 000 K is ~6.5, not ~0. The verdict is unchanged; the reason is not, and the reason is what
    decides whether a design change that raises the temperature actually helps.
    """
    cool_strong_b = conductivity.electrothermal_screen(3000.0, 0.32, 0.01, b_field=6.0)
    cool_weak_b = conductivity.electrothermal_screen(3000.0, 0.32, 0.01, b_field=0.005)
    hot = conductivity.electrothermal_screen(11000.0, 0.32, 0.01, b_field=6.0)

    assert cool_strong_b.at_risk, "cool + strongly magnetised is where the instability lives"
    assert not cool_weak_b.at_risk, "no Hall drive, no feedback"
    assert not hot.at_risk, "collisional enough that the Hall link is broken"
    # The hot end is carried by the Hall link only. The gain is *large* there, which is exactly
    # the correction ADR-0038 records: water takes the seed's place as the sensitive species.
    assert hot.hall_parameter < conductivity.BETA_CRIT < cool_strong_b.hall_parameter
    assert hot.ionisation_sensitivity > 1.0, "water's ionisation is sensitive, not inert"


def test_the_bag_is_collisional_enough_that_a_strong_field_is_needed_to_close_the_loop() -> None:
    """The quantitative form of the same finding, and the one number worth carrying to the paper.

    `beta = eB/(m_e nu)` and `nu ~ 3.7e11 /s` at the bag's density and 3000 K, so reaching the
    screening threshold takes several tesla. Below that the Hall link is broken and the uniform
    two-temperature description behind Q-F is safe from filamentation regardless of the gain.

    Diluting the plasma reverses this fast: `nu` falls with the neutral density, so a tenth of the
    density needs roughly a tenth of the field.

    **Re-baselined 2026-08-25, 4.73 -> 5.00 T (+5.6%),** when OH, H2 and O2 entered the species
    set. At 3000 K they are ~5% of the heavy particles, and they are collision partners: counting
    them raises `nu` and so raises the field the Hall link needs. The direction is the check --
    leaving them out of `_nu_electron_neutral` while the solve produced them made `sigma` rise and
    the required field *fall* to 4.18 T, which is the artifact, not the physics."""
    b_needed = 2.0 / conductivity.hall_parameter(3000.0, 0.32, 0.01, 1.0)
    assert b_needed == pytest.approx(5.00, rel=0.05)
    b_needed_dilute = 2.0 / conductivity.hall_parameter(3000.0, 0.032, 0.01, 1.0)
    assert b_needed_dilute < 0.2 * b_needed


def test_seed_window_carries_the_instability_data_without_assuming_a_field() -> None:
    """The table reports the field needed to reach the Hall threshold, not a Hall parameter at some
    assumed field. That keeps it free of a number this repository does not own, and it is the more
    useful form anyway: it is what the bag field has to be compared against.

    The field needed to close the Hall link rises monotonically, because Coulomb collisions make
    the plasma more collisional as it ionises.

    **The gain is not monotone** (ADR-0038). It falls as the seed saturates, bottoms out where
    potassium is spent and water has not yet started, then climbs again as water ionises -- and
    water at 13.6 eV against a hotter gas is just as sensitive as the seed at 4.34 eV against a
    cool one. An earlier version asserted a monotone fall to <0.05, which was the seed-only
    artifact rather than the physics."""
    rows = conductivity.seed_window(0.32, 0.01, conductivity.REF_V_L)
    assert all(a.b_field_for_beta_crit < b.b_field_for_beta_crit for a, b in pairwise(rows))

    gains = [r.ionisation_sensitivity for r in rows]
    trough = gains.index(min(gains))
    assert 0 < trough < len(gains) - 1, "the gain must have an interior minimum, not an endpoint"
    assert all(a > b for a, b in pairwise(gains[: trough + 1])), "falls as the seed saturates"
    assert gains[-1] > gains[trough], "and recovers once water takes over"
    # At the bag's density the trough never breaks the loop: the gain link is live throughout.
    assert min(gains) > conductivity.SENSITIVITY_CRIT
    assert rows[0].ionisation_sensitivity > 10.0


def test_critical_hall_parameter_diverges_as_the_plasma_approaches_equilibrium() -> None:
    """The analytic limit that decides this whole question (ADR-0038).

    Petit and Geffray put `(T_e - T_gas)` in the denominator of `s`, so `beta_cr -> infinity` as
    the electron temperature falls to the gas temperature. Physically: with no electron heating
    there is no drive, and no Hall parameter can make the runaway go.

    This is why `BETA_CRIT = 2` does not transfer to this plume. Two is the `s -> 0` limit of
    their fully-ionised form -- a *strongly* two-temperature plasma, which an MHD generator is and
    a near-equilibrium nozzle plume is not.
    """
    e_i, f = conductivity.IP_K, 1.0
    assert conductivity.critical_hall_parameter(4000.0, 4000.0, e_i, f) == math.inf
    assert conductivity.critical_hall_parameter(3900.0, 4000.0, e_i, f) == math.inf

    # Monotone: the closer to equilibrium, the harder it is to destabilise.
    near = conductivity.critical_hall_parameter(4100.0, 4000.0, e_i, f)
    far = conductivity.critical_hall_parameter(10000.0, 4000.0, e_i, f)
    assert near > far > 2.0


def test_critical_hall_parameter_reproduces_the_published_fully_ionised_limit() -> None:
    """`beta_cr = 1.935 f + 0.065 + s` must collapse to their stated `beta_cr ~ 2 + s` at `f = 1`.

    Worked independently at `T_e = 10 000 K`, `T_gas = 4000 K`, `E_i = 4.34 eV`:
        k T_e            = 0.861733 eV
        1.5 k T_e / E_i  = 0.297833
        s = 2(0.861733)(10000/6000)/4.34 / 1.297833 = 0.510044
        beta_cr = 1.935 + 0.065 + 0.510044 = 2.510044
    """
    beta_cr = conductivity.critical_hall_parameter(10000.0, 4000.0, conductivity.IP_K, 1.0)
    assert beta_cr == pytest.approx(2.5100, rel=1e-3)
    assert beta_cr == pytest.approx(2.0 + 0.5100, rel=1e-3)


def test_electron_energy_balance_actually_balances() -> None:
    """The residual test: at the returned `T_e`, Joule heating must equal elastic **plus**
    inelastic loss.

    This is the acceptance test for the solve itself, independent of whether the physics feeding
    it is right -- if the fixed point has not converged, everything downstream is decoration.
    """
    temp, rho, x_k = 4596.0, 0.02522, 0.01
    bal = conductivity.electron_energy_balance(temp, rho, x_k, 5.0, 6.0, 2.71e-3)
    assert bal.converged
    assert bal.t_e > bal.t_gas, "a driven plasma must sit above the gas temperature"

    elastic = (
        1.5
        * conductivity.electron_density(temp, rho, x_k, bal.t_e)
        * K_B
        * conductivity.energy_relaxation_rate(temp, rho, x_k, bal.t_e)
        * bal.elevation
    )
    inelastic = conductivity.inelastic_loss(temp, rho, x_k, bal.t_e)
    assert elastic + inelastic == pytest.approx(bal.joule_heating, rel=1e-4)
    # Both channels are live here: neither may be quietly zero at the station that decides the leg.
    assert inelastic > 0.2 * elastic


def test_a_vanishing_field_leaves_the_plume_in_equilibrium_and_unconditionally_stable() -> None:
    """The off-limit. No field means no driving current, so no electron heating, so no instability
    at any Hall parameter -- the `beta_cr -> infinity` branch reached through the physics rather
    than by asserting it directly."""
    loop = conductivity.electrothermal_loop(4596.0, 0.02522, 0.01, 1.0e-9, 6.0, 2.71e-3)
    assert loop.balance.elevation == pytest.approx(0.0, abs=1e-6)
    assert loop.critical_hall_parameter == math.inf
    assert not loop.unstable
    assert loop.e_folding_time == math.inf


def test_the_cold_leg_exit_is_unstable_and_the_hot_legs_are_not() -> None:
    """The verdict this module exists to deliver, at the nozzle exit of three closing speeds.

    The separation is not marginal in either direction, which is what makes it reportable. On the
    hot legs the plume is in equilibrium to within a few kelvin, `s` runs to 10^2-10^4 and
    `beta_cr` with it, against a Hall parameter of order 1. On the cold leg the plume falls far
    enough that Joule heating elevates the electrons ~1000 K, `beta_cr` drops to ~2, and `beta`
    is above it -- with an e-folding time of microseconds against a 2.7 ms transit.
    """
    hot = conductivity.electrothermal_loop(16224.0, 0.02512, 0.01, 5.0, 6.0, 1.69e-3)
    mid = conductivity.electrothermal_loop(11682.0, 0.02456, 0.01, 5.0, 6.0, 2.22e-3)
    cold = conductivity.electrothermal_loop(4596.0, 0.02522, 0.01, 5.0, 6.0, 2.71e-3)

    assert not hot.unstable and not mid.unstable
    assert hot.critical_hall_parameter > 1000.0, "equilibrium plume: beta_cr runs away"
    assert mid.critical_hall_parameter > 100.0

    assert cold.unstable, "the cold leg is where the criterion actually bites"
    assert cold.critical_hall_parameter < 3.0
    assert cold.screen.hall_parameter > cold.critical_hall_parameter
    # Microseconds against a millisecond transit: being there briefly is not a defence.
    assert cold.e_folding_time < 1.0e-4
    assert cold.e_folding_time * 100 < 2.71e-3


def test_inelastic_loss_vanishes_at_equilibrium_and_reverses_below_it() -> None:
    """Detailed balance, and it is the reason inelastic channels cannot rescue this plume.

    The heavies' internal states are held at `T_gas` by heavy-heavy collisions, so electron-impact
    excitation and super-elastic de-excitation cancel exactly when `T_e == T_gas`. Below it the
    electrons *gain* energy from the excited population, so the sign must flip.

    The consequence (ADR-0038): expanding the detailed-balance bracket for a small elevation gives
    `(dE/k) dT / T^2`, which is **linear in the elevation exactly as the elastic channel is**. So
    the inelastic-to-elastic ratio is a property of the state, not something that grows as the
    plasma approaches equilibrium -- inelastic losses cannot be a rescue mechanism for a
    near-equilibrium plasma.
    """
    temp, rho, x_k = 4596.0, 0.02522, 0.01
    assert conductivity.inelastic_loss(temp, rho, x_k, temp) == pytest.approx(0.0, abs=1e-9)
    assert conductivity.inelastic_loss(temp, rho, x_k, temp + 800.0) > 0.0
    assert conductivity.inelastic_loss(temp, rho, x_k, temp - 400.0) < 0.0


def test_the_alkali_resonance_dominates_and_water_has_no_target_left() -> None:
    """Channel ranking, and it is set by the composition rather than by the cross-sections.

    K 4s-4p is a resonance transition with oscillator strength ~1, so its cross-section is two to
    three orders above the molecular channels. The water channels are carried for completeness but
    the plume is 99.7-99.9998% dissociated wherever any of this matters (Q-M), so at the crossing
    station `n_H2O` is ~2e-6 of the heavies and they contribute nothing measurable.
    """
    original = conductivity.INELASTIC_CHANNELS
    try:
        per_channel = {}
        for channel in original:
            conductivity.INELASTIC_CHANNELS = (channel,)
            per_channel[channel.name] = conductivity.inelastic_loss(7569.0, 0.05810, 0.01, 8025.0)
    finally:
        conductivity.INELASTIC_CHANNELS = original

    alkali = per_channel["K 4s-4p resonance"]
    water = sum(v for k, v in per_channel.items() if k.startswith("H2O"))
    assert alkali > 100.0 * water, "the alkali channel must dominate at the crossing station"


def test_inelastic_channels_cool_the_electrons_but_do_not_save_the_cold_leg() -> None:
    """The verdict on ADR-0038's cheapest open item: real, helpful, and not enough.

    Including the channels drops the cold-leg exit elevation ~1055 K -> ~814 K and lifts `beta_cr`
    from ~1.84 to ~2.07. `beta` is 4.63, so the station remains unstable by a factor >2, and the
    unstable stretch of the leg is unchanged.

    **Re-baselined 2026-08-25** for the OH/H2/O2 species set: elevation 814 -> 922 K, `beta_cr`
    2.07 -> 1.99, `beta` 4.63 -> 4.81. Every structural claim is unchanged and the margin is
    slightly *worse*, so the finding survives its own refinement. `beta_cr` crossing back under
    the round 2.0 is why the bound below is stated against the pre-inelastic elevation instead of
    a constant: what the test is for is "helps, not enough", not the second digit.
    """
    exit_state = conductivity.electrothermal_loop(4596.0, 0.02522, 0.01, 5.0, 6.0, 2.71e-3)
    assert exit_state.balance.elevation < 1000.0, "inelastic cooling must reduce the elevation"
    assert exit_state.critical_hall_parameter > 1.9, "and must lift beta_cr"
    assert exit_state.unstable, "but not far enough: beta is still well above beta_cr"
    assert exit_state.screen.hall_parameter > 2.0 * exit_state.critical_hall_parameter

    # The hot legs stay in equilibrium and stay stable by orders of magnitude.
    hot = conductivity.electrothermal_loop(16224.0, 0.02512, 0.01, 5.0, 6.0, 1.69e-3)
    assert not hot.unstable
    assert hot.critical_hall_parameter > 1000.0
