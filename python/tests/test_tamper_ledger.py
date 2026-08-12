"""Acceptance tests for Rung 0 of the **tamped-nozzle study** (`puffsat_tamper_isp_prd.md`).

*Different study from `f(v)`.* These pin the analytic reference ledger — the single calculator
that owns every closed-form number the PRD quotes — so the PRD, the ADRs, and the eventual
analysis cannot drift apart (PRD §10, Rung 0).

**Where the targets come from.** Every expected value here is quoted from the PRD or from prior
work in this repository, *not* read back off the implementation, so these are acceptance tests in
the project's sense (CLAUDE.md: write the known answer first). Three families:

1. **PRD §8 analytic anchors** — the ceiling, the bare ballistic limit, the `k <= 1` zero-capture
   floor, the free-plate elastic bound, and the `k -> 0` degeneracy. These are the rung's exit
   criteria.
2. **Prior-work anchors** — `beta_bare(7.06) = 0.9087` and 1014 s at `w_bar = 77.28 km/s`, which
   confirm the inherited ballistic model before anything is built on it (PRD §3.2).
3. **Model invariants with teeth** — no model may exceed its own ejecta-mass ceiling, the
   ballistic fireball must conserve energy exactly, and `beta` must be invariant under `w`.

The *external-literature* benchmarks this study also needs (ice Hugoniots and release paths,
Sod/Noh/Sedov, Marshak, a two-material shock tube with an exact interface solution, and a
published hypervelocity penetration case) belong to Rung 1 and later, where kernel logic starts:
Rung 0 is closed-form algebra with no solver in it, so its only honest oracles are analytic
limits and the prior work it must reproduce.
"""

from __future__ import annotations

import math

import pytest

from puffsat.tamper import ledger

# ---------------------------------------------------------------------------------------------
# PRD §8 — analytic anchors specific to this device. These are the rung's exit criteria.
# ---------------------------------------------------------------------------------------------


def test_ceiling_is_the_ideal_collimation_limit() -> None:
    """§8: in the perfectly-collimated limit, `j = J/m_i -> w(sqrt(1+K_ej) - 1)`.

    The ceiling is a reformulation of prior work's ideal-collimation coefficient, so at
    `K_ej = k` it must reproduce `beta_ideal = sqrt(1+k) - 1` exactly (PRD §3.2).
    """
    for k_ej in (0.5, 1.0, 7.06, 14.12, 32.0):
        assert ledger.beta_ideal(k_ej) == pytest.approx(math.sqrt(1.0 + k_ej) - 1.0)
        assert ledger.j_max(k_ej) == pytest.approx(ledger.W_CLOSING_MS * ledger.beta_ideal(k_ej))

    # The two values §3.4 quotes for the reference cases. Its 2.889 was a rounding slip for
    # sqrt(15.12) - 1 = 2.8884; corrected in the PRD by this rung (see RECONCILIATIONS below).
    assert ledger.beta_ideal(7.06) == pytest.approx(1.839, abs=5e-4)
    assert ledger.beta_ideal(14.12) == pytest.approx(2.8884, abs=5e-4)


def test_bare_ballistic_limit() -> None:
    """§8: at collisionless expansion, `beta_bare(7.06) -> 0.9087`, capture fraction -> 0.3118."""
    assert ledger.beta_bare(7.06) == pytest.approx(0.90871, rel=1e-4)
    assert ledger.ballistic_capture_fraction(7.06) == pytest.approx(0.3118, abs=5e-5)
    # §3.4's second bare row, the same mass spent as extra slug.
    assert ledger.beta_bare(14.12) == pytest.approx(1.68353, rel=1e-4)


def test_ballistic_zero_capture_floor_below_k_equals_one() -> None:
    """§8: ballistically nothing reaches the plate below `k = 1`.

    The floor is a limit of the pressure-free ballistic model, not a general zero-thrust claim —
    measuring the small nonzero hydrodynamic thrust there is a Rung 1 job (PRD §2.2, §8). What
    Rung 0 owns is that the closed form is exactly zero and does not go negative.
    """
    for k in (0.0, 0.1, 0.5, 0.9, 1.0):
        assert ledger.ballistic_capture_fraction(k) == 0.0
        assert ledger.beta_bare(k) == 0.0
        assert ledger.beta_flat(k) == 0.0
    assert ledger.ballistic_capture_fraction(1.0 + 1e-9) == pytest.approx(0.0, abs=1e-9)
    assert ledger.ballistic_capture_fraction(4.0) == pytest.approx(0.25)


def test_free_plate_elastic_bound() -> None:
    """§8/§6.3: reflected/incident -> 0.274 at `tau_t = 1` in the ballistic limit.

    §6.3 states the cross-check as reproducing the prior analytic lower bound "to two figures",
    which is the precision this anchor is asserted at; the exact 1-D elastic value from the
    §6.3 mass pairing is 0.2732 (see `test_free_plate_bound_reconciliation`).
    """
    assert ledger.free_plate_reflected_fraction(tau_t=1.0, k=7.06) == pytest.approx(0.274, abs=5e-3)
    # Monotone in tamper mass, and it must bracket the mirror limit rather than reach it.
    ratios = [ledger.free_plate_reflected_fraction(tau_t=t, k=7.06) for t in (0.5, 1.0, 2.0, 4.0)]
    assert ratios == sorted(ratios)
    assert all(0.0 <= r < 1.0 for r in ratios)


def test_k_to_zero_degeneracy() -> None:
    """§8: zero net impulse as slug mass vanishes. A code producing thrust at `k = 0` is wrong."""
    assert ledger.beta_bare(0.0) == 0.0
    assert ledger.beta_mirror(0.0) == pytest.approx(0.0, abs=1e-12)
    assert ledger.beta_ideal(0.0) == pytest.approx(0.0)
    # Approached continuously from above, not merely defined at the endpoint.
    assert ledger.beta_mirror(1e-6) == pytest.approx(0.0, abs=1e-5)


# ---------------------------------------------------------------------------------------------
# Prior-work anchors (PRD §3.2, Rung 0's "confirm the inherited model" item).
# ---------------------------------------------------------------------------------------------


def test_reproduces_prior_work_1014_seconds() -> None:
    """§3.2: reproduce prior work's 1014 s at `w_bar = 77.28 km/s` from the same `beta_bare`."""
    isp = ledger.isp_eff(
        beta=ledger.beta_bare(ledger.K_BARE_OPTIMUM),
        c_ratio=ledger.K_BARE_OPTIMUM,
        w=ledger.W_BAR_PRIOR_MS,
    )
    assert isp == pytest.approx(1014.0, abs=1.0)


def test_bare_pass1_optimum_is_7_06() -> None:
    """§3.2/§3.5: the ballistic bare-plate Isp optimum sits at `k* = 7.060`, and is flat.

    Prior work's quoted 7.057 is the same optimum to optimiser tolerance. The flatness is
    load-bearing for the study (it is why `K` must be swept 6-32, PRD §3.5), so it is pinned too.
    """
    k_star = ledger.optimal_k_bare()
    assert k_star == pytest.approx(7.060, abs=5e-4)

    peak = ledger.isp_bare(k_star)
    assert peak == pytest.approx(984.0, abs=1.0)
    for k in (6.0, 8.0):
        assert ledger.isp_bare(k) / peak > 1.0 - 0.006  # flat to +-0.6% over k = 6-8


# ---------------------------------------------------------------------------------------------
# Model invariants — the checks that would catch a wrong ledger even with every anchor removed.
# ---------------------------------------------------------------------------------------------


def test_no_model_exceeds_its_own_ejecta_mass_ceiling() -> None:
    """Cauchy-Schwarz: `beta <= sqrt(1+K_ej) - 1` for the mass that actually moves (§3.2).

    Non-trivial for the perfect mirror, which reaches 96.9% of the `K_ej = k` ceiling: a sign or
    normalisation error in the reflected hemisphere breaks this immediately.
    """
    for k in (1.5, 3.0, 6.0, 7.06, 14.12, 32.0):
        ceiling = ledger.beta_ideal(k)
        assert ledger.beta_bare(k) < ceiling
        assert ledger.beta_flat(k) < ledger.beta_bare(k)
        assert ledger.beta_mirror(k) < ceiling
        # A finite plate can only lose material relative to R -> infinity.
        assert ledger.beta_finite(k, r_over_d=1.5) < ledger.beta_bare(k)


def test_ballistic_fireball_conserves_energy_exactly() -> None:
    """The isotropic ballistic model converts *all* merge-dissipated energy back to directed KE.

    Total ejecta KE per projectile kg is `(1/2)(1+k)(u^2 + V^2) = (1/2)w^2` identically — the
    incoming projectile KE. This is the energy audit the hydrocode must later close to <1% (§7.6),
    stated where it is exact.
    """
    w = ledger.W_CLOSING_MS
    for k in (0.5, 1.0, 7.06, 32.0):
        fb = ledger.Fireball.from_k(k, w=w)
        ke = 0.5 * (1.0 + k) * (fb.u**2 + fb.v_cm**2)
        assert ke == pytest.approx(0.5 * w**2, rel=1e-12)
        # ...and the blob's total momentum is exactly the incoming debit, so J = 0 with no plate.
        assert (1.0 + k) * fb.v_cm == pytest.approx(w, rel=1e-12)


def test_beta_is_dimensionless_and_w_invariant() -> None:
    """`beta` depends only on the mass ratios and geometry, never on `w` (PRD §0.2)."""
    for k in (2.0, 7.06, 20.0):
        assert ledger.beta_bare(k) == pytest.approx(ledger.beta_bare(k))
        for w in (40_000.0, 75_000.0, 81_000.0):
            fb = ledger.Fireball.from_k(k, w=w)
            assert fb.u / w == pytest.approx(math.sqrt(k) / (1.0 + k))
            assert fb.v_cm / w == pytest.approx(1.0 / (1.0 + k))
    # Isp does scale with w, linearly.
    assert ledger.isp_eff(beta=1.0, c_ratio=7.06, w=150_000.0) == pytest.approx(
        2.0 * ledger.isp_eff(beta=1.0, c_ratio=7.06, w=75_000.0)
    )


def test_mass_ledger_closes() -> None:
    """§0.1: `m_enc = m_i + m_hydro`, `K = k(1+tau_t+mu)`, `C = K + a_abl + a_other`."""
    led = ledger.MassLedger.from_ratios(k=7.06, tau_t=1.0, mu=0.35, a_abl=0.4, a_other=0.1)
    assert led.k_hydro == pytest.approx(7.06 * (1.0 + 1.0 + 0.35))
    assert led.c_charged == pytest.approx(led.k_hydro + 0.4 + 0.1)
    assert led.k_ej == pytest.approx(led.k_hydro)  # Pass 1 convention: closed, no ablator ejecta
    assert led.m_i + led.m_hydro == pytest.approx(led.m_enc)
    assert led.m_s + led.m_t + led.m_int == pytest.approx(led.m_hydro)
    assert led.m_charged == pytest.approx(led.c_charged * led.m_i)


# ---------------------------------------------------------------------------------------------
# PRD §4 — the analytic no-interlayer reference case.
# ---------------------------------------------------------------------------------------------


def test_section_4_reference_masses() -> None:
    """§4: the 200 kg reference encounter at `k = 7.06`, `mu = 0`, `tau_t` = 0 and 1."""
    bare = ledger.MassLedger.from_ratios(k=7.06, tau_t=0.0)
    assert bare.m_i == pytest.approx(24.81, abs=5e-3)
    assert bare.m_s == pytest.approx(175.2, abs=5e-2)
    assert bare.m_t == 0.0
    assert ledger.projectile_energy(bare.m_i) == pytest.approx(69.8e9, rel=1e-3)

    tamped = ledger.MassLedger.from_ratios(k=7.06, tau_t=1.0)
    assert tamped.m_i == pytest.approx(13.23, abs=5e-3)
    assert tamped.m_s == pytest.approx(93.4, abs=5e-2)
    assert tamped.m_t == pytest.approx(93.4, abs=5e-2)
    assert ledger.projectile_energy(tamped.m_i) == pytest.approx(37.2e9, rel=1e-3)


def test_section_4_fireball_state() -> None:
    """§4: `V = 9.31 km/s`, `u = 24.72 km/s`, `e = 305.7 MJ/kg = 57.1 eV` per H2O molecule."""
    fb = ledger.Fireball.from_k(7.06)
    assert fb.v_cm == pytest.approx(9.31e3, abs=5.0)
    assert fb.u == pytest.approx(24.72e3, abs=5.0)
    assert fb.specific_internal_energy == pytest.approx(305.7e6, rel=1e-3)
    assert ledger.ev_per_water_molecule(fb.specific_internal_energy) == pytest.approx(
        57.1, abs=0.05
    )


def test_section_4_1_vehicle_context() -> None:
    """§4.1: 1.69e6 N.s and 1.69 m/s per pulse; the departure burn ledger."""
    ctx = ledger.vehicle_context(k=7.06, tau_t=0.0)
    assert ctx.impulse_per_pulse == pytest.approx(1.69e6, rel=5e-3)
    assert ctx.delta_v_per_pulse == pytest.approx(1.69, rel=5e-3)
    assert ctx.acceleration_g(1.0) == pytest.approx(0.17, abs=5e-3)
    assert ctx.acceleration_g(4.0) == pytest.approx(0.69, abs=5e-3)

    burn = ledger.departure_burn(isp_s=ctx.isp_eff, carried_per_pulse_kg=ctx.m_charged)
    assert burn.mass_ratio == pytest.approx(2.08, abs=5e-3)
    assert burn.charged_mass_t == pytest.approx(519.0, abs=1.0)
    assert burn.pulses == pytest.approx(2960.0, rel=5e-3)
    assert burn.duration_s(cadence_hz=4.0) == pytest.approx(740.0, rel=5e-3)


# ---------------------------------------------------------------------------------------------
# PRD §3.4 — the reference comparison, row by row.
# ---------------------------------------------------------------------------------------------


def test_section_3_4_reference_comparison() -> None:
    """§3.4: bare / extra-slug / perfect-mirror / break-even / ceiling, at `w = 75 km/s`."""
    rows = {row.label: row for row in ledger.reference_comparison()}

    bare = rows["bare plate, k = 7.06"]
    assert bare.beta == pytest.approx(0.9087, abs=5e-4)
    assert bare.realization == pytest.approx(0.494, abs=5e-4)
    assert bare.isp_s == pytest.approx(984.0, abs=1.0)

    slug = rows["bare plate, k = 14.12"]
    assert slug.beta == pytest.approx(1.6835, abs=5e-4)
    assert slug.realization == pytest.approx(0.583, abs=5e-4)
    assert slug.isp_s == pytest.approx(912.0, abs=1.0)

    mirror = rows["perfect-mirror tamper, tau_t = 1"]
    assert mirror.beta == pytest.approx(1.7825, abs=5e-4)
    assert mirror.realization == pytest.approx(0.617, abs=5e-4)
    assert mirror.isp_s == pytest.approx(965.0, abs=1.0)

    breakeven = rows["break-even against the bare plate at k = 7.06"]
    assert breakeven.beta == pytest.approx(1.817, abs=5e-4)
    assert breakeven.realization == pytest.approx(0.629, abs=5e-4)
    assert breakeven.isp_s == pytest.approx(984.0, abs=1.0)

    ceiling = rows["ceiling at K = 14.12"]
    assert ceiling.beta == pytest.approx(2.8884, abs=5e-4)
    assert ceiling.realization == pytest.approx(1.0)
    assert ceiling.isp_s == pytest.approx(1565.0, abs=1.0)
    assert rows["ceiling at K = 7.06"].isp_s == pytest.approx(1992.0, abs=1.0)

    # The headline reading: the perfect mirror loses to not spending the mass at all.
    assert mirror.realization < breakeven.realization
    assert mirror.realization > slug.realization


def test_break_even_is_a_function_not_a_constant() -> None:
    """§3.4: 62.9% is one evaluation of the break-even rule, not a universal gate."""
    beta_ref, c_ref = ledger.beta_bare(7.06), 7.06

    required = ledger.required_realization(
        beta_ref=beta_ref, c_ref=c_ref, c_candidate=14.12, k_ej_candidate=14.12
    )
    assert required == pytest.approx(0.629, abs=5e-4)

    # A candidate that charges *less* mass needs less of its own ceiling; the number moves.
    cheaper = ledger.required_realization(
        beta_ref=beta_ref, c_ref=c_ref, c_candidate=10.0, k_ej_candidate=10.0
    )
    assert cheaper != pytest.approx(0.629, abs=1e-2)

    # Self-consistency: the reference compared with itself needs exactly its own realization.
    same = ledger.required_realization(
        beta_ref=beta_ref, c_ref=c_ref, c_candidate=c_ref, k_ej_candidate=c_ref
    )
    assert same == pytest.approx(ledger.realization_fraction(beta_ref, 7.06))

    # And the break-even beta scales with charged mass, as §3.4's 1.817 = 2 x 0.90871 does.
    assert ledger.break_even_beta(
        beta_ref=beta_ref, c_ref=c_ref, c_candidate=14.12
    ) == pytest.approx(2.0 * beta_ref)


# ---------------------------------------------------------------------------------------------
# PRD §3.5 — the mass-ratio trade table.
# ---------------------------------------------------------------------------------------------


def test_section_3_5_mass_ratio_table() -> None:
    """§3.5: Isp degrades gently with `K` while energy per unit impulse improves steeply."""
    expected = {
        6.0: (0.768, 979.0, 48.8e3),
        7.06: (0.909, 984.0, 41.3e3),
        8.0: (1.027, 981.0, 36.5e3),
        10.0: (1.260, 963.0, 29.8e3),
        14.12: (1.684, 912.0, 22.3e3),
        16.0: (1.858, 888.0, 20.2e3),
        32.0: (3.068, 733.0, 12.2e3),
    }
    for k, (beta, isp, e_per_j) in expected.items():
        assert ledger.beta_bare(k) == pytest.approx(beta, abs=5e-4)
        assert ledger.isp_bare(k) == pytest.approx(isp, abs=1.0)
        assert ledger.energy_per_impulse(ledger.beta_bare(k)) == pytest.approx(e_per_j, rel=5e-3)

    # Plate heat load and projectile consumption are the same function of beta (§3.5).
    for k in (6.0, 10.0, 32.0):
        beta = ledger.beta_bare(k)
        assert ledger.projectiles_per_impulse(beta) == pytest.approx(
            1.0 / (beta * ledger.W_CLOSING_MS)
        )
        assert ledger.energy_per_impulse(beta) * ledger.projectiles_per_impulse(beta) > 0.0

    # The ceiling has no interior optimum in K_ej: v_e,max falls monotonically (§3.5).
    v_e = [ledger.v_e_max(k) for k in (6.0, 7.06, 10.0, 14.12, 32.0)]
    assert v_e == sorted(v_e, reverse=True)


# ---------------------------------------------------------------------------------------------
# PRD §3.6 — the worked example, line by line. "Rung 0's calculator must reproduce every line."
# ---------------------------------------------------------------------------------------------


def test_section_3_6_worked_example_at_k_6() -> None:
    """§3.6: the full ledger behind one row of §3.5, at `k = 6`."""
    ex = ledger.worked_example(k=6.0)

    assert ex.v_cm == pytest.approx(10.71e3, abs=5.0)
    assert ex.u == pytest.approx(26.24e3, abs=5.0)
    assert ex.u_over_v == pytest.approx(2.45, abs=5e-3)
    assert ex.mu_capture == pytest.approx(0.408, abs=5e-4)
    assert ex.capture_fraction == pytest.approx(0.296, abs=5e-4)
    assert ex.captured_mass_per_projectile_kg == pytest.approx(2.07, abs=5e-3)

    assert ex.best_aimed_axial == pytest.approx(15.53e3, abs=5.0)
    assert ex.mean_axial_approach == pytest.approx(7.77e3, abs=5.0)
    assert ex.mean_speed_leaving == pytest.approx(20.04e3, abs=5.0)
    assert ex.delta_p_per_captured_kg == pytest.approx(27.81e3, abs=5.0)

    assert ex.v_e == pytest.approx(9.60e3, abs=5.0)
    assert ex.isp_s == pytest.approx(979.0, abs=1.0)

    # The flat-plate counterpart: the collimation prize is inside the 979 s, not on top of it.
    assert ledger.beta_flat(6.0) == pytest.approx(0.429, abs=5e-4)
    assert ex.isp_flat_s == pytest.approx(547.0, abs=1.0)

    # Why it is ~half the ceiling.
    assert ex.coherent_speed == pytest.approx(28.35e3, abs=5.0)
    assert ex.v_e_max == pytest.approx(20.57e3, abs=5.0)
    assert ex.isp_max_s == pytest.approx(2098.0, abs=1.0)
    assert ex.realization == pytest.approx(0.467, abs=5e-4)

    # Nothing but the plate contributes: the blob alone carries exactly the incoming debit.
    assert ex.blob_momentum_per_projectile_kg == pytest.approx(ledger.W_CLOSING_MS, rel=1e-12)


def test_flat_plate_closed_form_matches_the_integral() -> None:
    """§3.6: `beta_flat = (sqrt(k)-1)^2 / (2 sqrt(k))` must equal the integrated `2 v_z` capture."""
    for k in (1.5, 6.0, 7.06, 20.0):
        closed = (math.sqrt(k) - 1.0) ** 2 / (2.0 * math.sqrt(k))
        assert ledger.beta_flat(k) == pytest.approx(closed, rel=1e-12)
        assert ledger.beta_finite(k, r_over_d=math.inf, plate="flat") == pytest.approx(
            closed, rel=1e-9
        )


# ---------------------------------------------------------------------------------------------
# PRD §13.13 — the three capture-fraction conventions Rung 0 must reconcile.
# ---------------------------------------------------------------------------------------------


def test_three_capture_conventions_reconcile() -> None:
    """§3.6/§13.13: 31.2% (`R -> inf`), 22.3% (blob-frame rim angle), 10.6% (finite-plate ray)."""
    k = 7.06
    assert ledger.ballistic_capture_fraction(k) == pytest.approx(0.312, abs=5e-4)
    assert ledger.capture_fraction_rim_angle(k, r_over_d=1.5) == pytest.approx(0.223, abs=5e-4)
    assert ledger.capture_fraction_finite(k, r_over_d=1.5) == pytest.approx(0.106, abs=5e-4)
    assert ledger.capture_fraction_finite(k, r_over_d=2.5) == pytest.approx(0.164, abs=5e-4)

    # The finite-plate ray value is the lower edge of the bracket and rises to the R -> inf limit.
    assert ledger.capture_fraction_finite(k, r_over_d=1e6) == pytest.approx(
        ledger.ballistic_capture_fraction(k), abs=1e-5
    )
    for r_over_d in (1.0, 1.5, 2.5, 10.0):
        assert ledger.capture_fraction_finite(k, r_over_d) <= ledger.ballistic_capture_fraction(k)


def test_section_3_6_finite_plate_table() -> None:
    """§3.6: the `R/d` table — parabolic and flat `beta` and Isp, and the parabola/flat ratio."""
    k = 7.06
    expected = {
        1.5: (0.339, 368.0, 0.292, 317.0, 1.16),
        2.5: (0.511, 553.0, 0.400, 434.0, 1.28),
    }
    for r_over_d, (b_par, isp_par, b_flat, isp_flat, ratio) in expected.items():
        par = ledger.beta_finite(k, r_over_d, plate="parabolic")
        flat = ledger.beta_finite(k, r_over_d, plate="flat")
        assert par == pytest.approx(b_par, abs=5e-4)
        assert flat == pytest.approx(b_flat, abs=5e-4)
        assert ledger.isp_eff(beta=par, c_ratio=k) == pytest.approx(isp_par, abs=1.0)
        assert ledger.isp_eff(beta=flat, c_ratio=k) == pytest.approx(isp_flat, abs=1.0)
        assert par / flat == pytest.approx(ratio, abs=5e-3)

    # R -> infinity row: 0.909 (984 s) parabolic against 0.517 (560 s) flat. The PRD's quoted
    # 1.79x ratio disagreed with its own beta pair (0.9087/0.5167 = 1.759); corrected by this rung.
    assert ledger.beta_flat(k) == pytest.approx(0.517, abs=5e-4)
    assert ledger.isp_eff(beta=ledger.beta_flat(k), c_ratio=k) == pytest.approx(560.0, abs=1.0)
    assert ledger.beta_bare(k) / ledger.beta_flat(k) == pytest.approx(1.759, abs=5e-3)


# ---------------------------------------------------------------------------------------------
# Reconciliations Rung 0 owes the PRD, pinned so they cannot silently drift back.
# ---------------------------------------------------------------------------------------------


def test_free_plate_bound_reconciliation() -> None:
    """The exact 1-D elastic value from §6.3's own mass pairing is 0.2732, not 0.274.

    §6.3 pairs a tamper of `tau_t * k` against an away-going plume of `(1+k)/2` per projectile kg
    and quotes the ratio to two figures (0.27). Pinned here so the reconciliation stands.
    """
    assert ledger.free_plate_reflected_fraction(tau_t=1.0, k=7.06) == pytest.approx(
        0.2732, abs=5e-5
    )
    assert ledger.free_plate_areal_ratio(tau_t=1.0, k=7.06) == pytest.approx(1.752, abs=5e-4)


def test_superseded_lambda_framing_is_reported_not_used() -> None:
    """§1's 591-965 s bracket is the *superseded* `Lambda` framing (ADR-0030, CONTEXT.md).

    `Lambda = beta_tamped/beta_bare` hides the ceiling and miscounts recoil as loss, so the ledger
    exposes it only as provenance for the quoted bracket, never as a metric.
    """
    lam = ledger.legacy_lambda(
        beta_tamped=ledger.beta_mirror(7.06), beta_reference=ledger.beta_bare(7.06)
    )
    assert lam == pytest.approx(1.962, abs=5e-4)  # handoff's Lambda_ideal
    # The bracket endpoints §1 quotes, reproduced from Lambda and the bare Isp.
    isp_bare = ledger.isp_bare(7.06)
    assert isp_bare * 1.962 / 2.0 == pytest.approx(965.0, abs=1.0)
    assert isp_bare * 1.2 / 2.0 == pytest.approx(591.0, abs=1.0)


def _beta_by_quadrature(k: float, r_over_d: float, plate: str, panels: int = 20_000) -> float:
    """Independent midpoint-rule evaluation of the same cone integral `beta_finite` closed-forms."""
    u, v = math.sqrt(k) / (1.0 + k), 1.0 / (1.0 + k)
    mu_lo = ledger.mu_capture_finite(k, r_over_d)
    if mu_lo >= 1.0:
        return 0.0
    total = 0.0
    for i in range(panels):
        mu = mu_lo + (i + 0.5) * (1.0 - mu_lo) / panels
        v_z = u * mu - v
        speed = math.sqrt(u * u + v * v - 2.0 * u * v * mu)
        total += 2.0 * v_z if plate == "flat" else speed + v_z
    return (1.0 + k) * 0.5 * total * (1.0 - mu_lo) / panels


def test_closed_forms_match_independent_quadrature() -> None:
    """The analytic cone integrals must agree with brute-force numerical integration.

    The closed forms exist so a "reference" number can never carry a silent convergence error, but
    that only helps if the algebra is right. This checks the algebra by a genuinely different
    method, which is the same discipline the kernels get from their order-of-accuracy tests.
    """
    plates: tuple[ledger.PlateShape, ...] = ("parabolic", "flat")
    for k in (2.0, 6.0, 7.06, 14.12, 32.0):
        for r_over_d in (1.5, 2.5, math.inf):
            for plate in plates:
                closed = ledger.beta_finite(k, r_over_d, plate=plate)
                assert closed == pytest.approx(
                    _beta_by_quadrature(k, r_over_d, plate), rel=1e-6, abs=1e-9
                )


def test_plate_soak_chain_resolves_the_alpha_th_discrepancy() -> None:
    """PRD §6.5 quoted a 173 um soak depth, which implies `alpha_th = 1.0e-5`, not §0.6's 1.2e-5.

    Rung 0's checklist names this as a known open reconciliation and requires it be resolved in one
    place, because the depth propagates into an areal mass, a capacity, a whole-plate basis, four
    regenerative shares, and the regenerative cap. Pinning both ends here shows the discrepancy was
    exactly the diffusivity and nothing else.
    """
    # The stale figure is reproduced *only* by the stale diffusivity — that identifies the drift.
    assert ledger.plate_soak_chain(alpha_th=1.0e-5).depth_m == pytest.approx(173e-6, abs=1e-6)

    chain = ledger.plate_soak_chain()
    assert chain.alpha_th == 1.2e-5
    assert chain.depth_m == pytest.approx(190e-6, abs=1e-6)
    assert chain.areal_mass_kg_m2 == pytest.approx(1.49, abs=5e-3)
    assert chain.capacity_j_m2 == pytest.approx(1.04e6, rel=5e-3)
    assert chain.share_of_fluence == pytest.approx(0.083, abs=5e-4)
    assert chain.plate_area_m2 == pytest.approx(706.9, abs=0.1)
    assert chain.basis_j == pytest.approx(737e6, rel=5e-3)


def test_regenerative_budget_and_cap() -> None:
    """§6.5.3: the free (mass-free) regenerative budget, restated on the corrected basis."""
    rows = {row.water_state: row for row in ledger.regenerative_budget()}
    assert rows["ice warmed to just under melt"].share_of_basis == pytest.approx(0.090, abs=5e-4)
    assert rows["melted to liquid at 273 K"].share_of_basis == pytest.approx(0.179, abs=5e-4)
    assert rows["liquid at 373 K"].share_of_basis == pytest.approx(0.293, abs=5e-4)
    assert rows["saturated steam"].share_of_basis == pytest.approx(0.906, abs=5e-4)

    # Phase, not enthalpy, sets the cap: boiled coolant cannot go back to being a slug.
    assert not rows["saturated steam"].admissible_as_slug
    assert ledger.regenerative_cap() == pytest.approx(0.293, abs=5e-4)


def test_conducted_soak_is_two_orders_below_the_capacity() -> None:
    """§6.5.2: the physical soak while the ablator pins the surface, against the capacity bound.

    This is the lower edge of the "plausible soak range", restated from its own derivation: it is
    the conducted-in energy at the two `T_abl` design points, not a rounded inheritance.
    """
    lo = ledger.conducted_per_pulse_j_m2(delta_t=800.0 - 753.0)
    hi = ledger.conducted_per_pulse_j_m2(delta_t=1000.0 - 905.0)
    assert lo / ledger.INCIDENT_FLUENCE_J_M2 == pytest.approx(0.00075, abs=5e-6)
    assert hi / ledger.INCIDENT_FLUENCE_J_M2 == pytest.approx(0.00151, abs=5e-6)
    # Two orders below the penetration-depth capacity, which is why the plate is not the constraint.
    assert hi < ledger.plate_soak_chain().capacity_j_m2 / 50.0


def test_assumption_audit_covers_every_closed_form_model() -> None:
    """Every model function must name the assumptions it inherits, and each must be adjudicated.

    Seth's standing instruction: sanity-check that the inherited assumptions are still relevant
    and appropriate at each step, rather than carrying them forward silently.
    """
    audit = ledger.assumption_audit()
    assert audit, "the audit must not be empty"
    for entry in audit:
        assert entry.verdict in {"holds", "holds-with-caveat", "departs"}
        assert entry.why and entry.retired_by
        # A departure must say what replaces it, not merely that it fails.
        if entry.verdict == "departs":
            assert entry.retired_by != "not yet"

    covered = {name for entry in audit for name in entry.applies_to}
    for model in ("beta_bare", "beta_flat", "beta_mirror", "beta_ideal", "beta_finite"):
        assert model in covered, f"{model} has no recorded assumption"
