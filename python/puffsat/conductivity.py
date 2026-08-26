"""Electrical conductivity of the potassium-seeded water plasma (Study 1).

The paper's `tab:seed_window` reports a magnetic Reynolds number against temperature, and the field
leak that `tab:bag_state` depends on is set by where `Rm = mu0 sigma v L` falls through 1. So the
table is a conductivity calculation, and conductivity is set by the electron density.

**The electron density has two sources and they swap dominance.** Below ~5000 K the potassium seed
(`chi = 4.34 eV`) supplies essentially all the electrons; above it water's own ionisation does, by
~38x at 15 000 K. The first hand estimate counted only the seed and was low by that factor above
5000 K. Both are carried here.

**Electron temperature is explicit.** `t_e` defaults to the gas temperature, which reproduces the
equilibrium result exactly. It is a parameter rather than an assumption because Kerrebrock's
non-equilibrium MHD result -- that the electron temperature can decouple from the gas in a seeded
plasma -- would move the cliff, and possibly remove it. Carrying `t_e` from the start makes that a
sweep rather than a rewrite.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from puffsat.eos_water import AMU, K_B, M_E, ln_k_saha
from puffsat.eos_water import composition as water_composition

#: Potassium first ionisation potential [J]. The alkali seed's whole purpose: 4.34 eV against
#: water's 12.6-13.6 eV, so it ionises where the water cannot.
IP_K = 4.34 * 1.602176634e-19
#: Ground-state degeneracies: neutral K is `2S_{1/2}` (g = 2), K+ is closed-shell Ar-like (g = 1).
G_K_NEUTRAL = 2.0
G_K_ION = 1.0
#: Potassium atomic mass [kg].
M_K = 39.0983 * AMU


def seed_number_density(rho: float, x_k: float) -> float:
    """Total potassium number density [m^-3] from the seed **mass** fraction `x_k`."""
    return x_k * rho / M_K


def electron_density(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """Free-electron density [m^-3] from water's own ionisation plus the potassium seed.

    Water is solved by `eos_water.composition` on the water mass alone; the seed's ionisation is
    then solved in that electron field. Because the seed's ions and its electrons are the same
    population, the closure is exact rather than iterative: with `K = n_e n_K+ / n_K0` and
    `n_K0 + n_K+ = n_K`,

        n_e^2 + n_e (K - n_w) - K (n_w + n_K) = 0

    whose positive root is returned. That reduces to `sqrt(n_K K)` when the seed is weakly ionised
    and to `n_w + n_K` when it is fully ionised, which are the two regimes the table spans.

    **Known approximation.** Water's ionisation is computed without the seed's electrons in the
    field, so the seed's suppression of water's own Saha balance (Le Chatelier) is not fed back.
    The seed adds ~2.6% to the electron density at 15 000 K, which would push water's contribution
    down by ~1.3% -- below the uncertainty in `Q_en` and `ln Lambda` that dominate `sigma`.
    """
    t_e = temp if t_e is None else t_e
    n_w = water_composition((1.0 - x_k) * rho, temp).n_e
    n_k = seed_number_density(rho, x_k)
    k_saha = math.exp(ln_k_saha(IP_K, G_K_ION, G_K_NEUTRAL, t_e))
    b = k_saha - n_w
    return 0.5 * (-b + math.sqrt(b * b + 4.0 * k_saha * (n_w + n_k)))


#: Electron-neutral momentum-transfer cross-section [m^2]. A single generic value stands in for the
#: real, species- and energy-dependent set (H2O, H, O, K each differ, and each varies with electron
#: energy). It is the model's weakest input: the audit flagged it as hand-picked, and `sigma` is
#: directly inverse in it wherever neutrals dominate the scattering. Exposed as a parameter of
#: `sigma` so it can be swept rather than believed.
Q_EN = 1.0e-19
#: Vacuum permittivity [F/m] and elementary charge [C].
EPS0 = 8.8541878128e-12
E_CHARGE = 1.602176634e-19
#: Spitzer's electron-electron correction for `Z = 1`: the Lorentz-gas result is 1.96x too small.
SPITZER_EE = 1.96


def coulomb_logarithm(n_e: float, t_e: float) -> float:
    """`ln Lambda` with `Lambda = 12 pi n_e lambda_D^3`, the number of electrons in a Debye sphere.

    **This runs small here -- about 2.5-2.8 across the seed window -- and Spitzer's transport theory
    assumes it is large.** At `ln Lambda ~ 2.5` the weak-coupling expansion it comes from is at the
    edge of its validity, so the Coulomb collision frequency (and hence `sigma` wherever Coulomb
    scattering dominates) carries an uncertainty of order tens of percent that no amount of care
    inside this function removes. Reported rather than smoothed: `Rm` claims that depend on a factor
    of ~2 in `sigma` are not decidable by this model.
    """
    lambda_d = math.sqrt(EPS0 * K_B * t_e / (n_e * E_CHARGE * E_CHARGE))
    return float(math.log(12.0 * math.pi * n_e * lambda_d**3))


def _nu_coulomb(n_e: float, t_e: float) -> float:
    """Electron-ion Coulomb collision frequency [s^-1], NRL formulary form for `Z = 1`.

    `nu_ei = 2.91e-6 n_e[cm^-3] ln Lambda T_e[eV]^-1.5`, divided by Spitzer's `1.96` so that the
    pure-Coulomb limit of `sigma` reproduces Spitzer's conductivity rather than the Lorentz-gas
    value. Mixing that correction with the neutral channel below is approximate -- the `1.96` is
    derived for a fully-ionised gas -- but it is right in the limit that matters for the test.
    """
    t_ev = K_B * t_e / E_CHARGE
    n_e_cgs = n_e * 1.0e-6
    return float(2.91e-6 * n_e_cgs * coulomb_logarithm(n_e, t_e) * t_ev**-1.5 / SPITZER_EE)


def _nu_electron_neutral(temp: float, rho: float, x_k: float, t_e: float, q_en: float) -> float:
    """Electron-neutral collision frequency [s^-1], `nu = n_n q_en <v_th>`.

    `<v_th> = sqrt(8 k T_e / pi m_e)` is the mean electron thermal speed. Neutrals are counted from
    the water solve -- H2O, OH, H2, O2, H, O -- plus un-ionised potassium. The three molecular
    intermediates matter here: through the seed window they hold up to ~5% of the heavy particles,
    and leaving them out would count that fraction as no collision partner at all, overstating
    `sigma` exactly where the conductivity cliff is being located.
    """
    comp = water_composition((1.0 - x_k) * rho, temp)
    n_neutral = comp.n_neutral_heavy
    v_th = math.sqrt(8.0 * K_B * t_e / (math.pi * M_E))
    return n_neutral * q_en * v_th


def sigma(
    temp: float,
    rho: float,
    x_k: float,
    t_e: float | None = None,
    q_en: float = Q_EN,
) -> float:
    """Electrical conductivity [S/m] of the seeded plasma.

    `sigma = n_e e^2 / (m_e (nu_ei + nu_en))` -- a Drude/Lorentz form with the two scattering
    channels added as rates (Matthiessen). Coulomb scattering dominates once the gas is appreciably
    ionised; electron-neutral scattering dominates at the cool end, where the seed is the only
    electron source and the gas is almost all neutral.

    **Two inputs carry most of the uncertainty, and neither is settled here.** `q_en` is a single
    generic cross-section standing in for four species with energy-dependent values, and
    `ln Lambda ~ 2.5-2.8` sits at the edge of where Spitzer's theory applies. Together they make a
    factor-of-two claim on `sigma` undecidable by this model. See `coulomb_logarithm`.
    """
    t_e = temp if t_e is None else t_e
    n_e = electron_density(temp, rho, x_k, t_e)
    nu = _nu_coulomb(n_e, t_e) + _nu_electron_neutral(temp, rho, x_k, t_e, q_en)
    return n_e * E_CHARGE * E_CHARGE / (M_E * nu)


#: The six temperatures the paper's `tab:seed_window` reports.
SEED_WINDOW_TEMPS: tuple[float, ...] = (2000.0, 3000.0, 5000.0, 8000.0, 11000.0, 15000.0)
#: Reference conditions. `rho` and `x_K` are the paper's slug bag values.
REF_RHO = 0.32
REF_X_K = 0.01

#: The expansion product `v L` this table is reported at [m^2/s] -- an **output of the solved
#: expansion**, not a back-solve of the paper. `expansion.HistoryRow.v_l` is `speed x flux-tube
#: radius` at every station; this is its **nozzle-exit** value on the 56.53 km/s equilibrium leg
#: (7.44e4, quoted 7.4e4), which is the leg the paper states. Run `make analysis-expansion` and read
#: the last station of `data/results/cooling_history.csv` to reproduce it.
REF_V_L = 7.4e4
#: The band the four flown legs actually span at the nozzle exit, low to high: 5.54e4 at 45.58 km/s
#: on the frozen branch, 9.72e4 at 75 km/s on the equilibrium branch. Quoted to two digits because
#: the cliff moves only 2386-2524 K across the whole of it, so nothing rests on the middle.
V_L_BAND: tuple[float, float] = (5.5e4, 9.7e4)
#: Retired (2026-08-26, deferred item D9). It came from back-solving the paper's own `Rm = 361` row
#: at 15 000 K, and it encodes a conductivity this module contradicts: `361 / (mu0 x 1.81e4)` is
#: `15 872 S/m`, against the `6993 S/m` computed here. Those are the same 2.3x gap that
#: `test_blended_conductivity_reproduces_the_audit_hand_calculation` records against the audit's
#: corrected hand blend of ~6950 S/m -- so at this `v L` the module gives `Rm = 159`, not 361, and
#: the constant and the row it was fitted to cannot both be right. Superseded because `v L` is now
#: an output of the solved expansion instead of a fit to a table. Kept named so `main` prints the
#: cliff it implied (2859 K) beside the one that holds.
#:
#: (It happens to equal the *throat-end* `v L` of the 56.53 km/s equilibrium leg to three digits.
#: That is a coincidence -- a fit to a hand estimate has no reason to land on a station -- and
#: nothing is built on it.)
RETIRED_V_L = 1.81e4
DEFAULT_SEED_WINDOW_PATH = Path("data/results/seed_window.csv")

#: Vacuum permeability [H/m].
MU_0 = 4.0e-7 * math.pi


def magnetic_reynolds(sigma_value: float, v_l: float) -> float:
    """`Rm = mu0 sigma (v L)` -- the field-diffusion time over the expansion time.

    **`v` and `L` enter only as their product**, which is why this takes `v_l` rather than the two
    separately: `tau_d/t_exp = (mu0 sigma L^2)/(L/v) = mu0 sigma v L`. Item 10's leak bracket is
    `~1/Rm`, so this is the whole of what that argument needs -- no MHD solve.

    There is deliberately **no default** for `v_l`. The paper's `tab:seed_window` reports an `Rm`
    column without stating either quantity, so any default here would be inventing the input that
    makes the column reproducible. The caller must say what expansion it means. `REF_V_L` is what
    `main` says it means, and it is an output of the solved expansion rather than an invention --
    but it is a module constant the caller can see, not a default hidden in this signature.
    """
    return MU_0 * sigma_value * v_l


def cliff_temperature(
    rho: float,
    x_k: float,
    v_l: float,
    t_lo: float = 1500.0,
    t_hi: float = 15000.0,
    t_e_offset: float = 0.0,
) -> float:
    """Temperature [K] at which `Rm` falls through 1 -- the conductivity cliff, as an **output**.

    The paper sets the seed window's floor by where the field stops being held. Rather than assert a
    floor, this solves for it: bisect on `Rm(T) = 1`, with `sigma` rising steeply through the seed's
    ionisation so the crossing is sharp and the bisection well-conditioned.

    `t_e_offset` raises the electron temperature above the gas temperature by a fixed amount, the
    crude form of Kerrebrock's non-equilibrium effect. It is here so the question "does the cliff
    survive electron-temperature decoupling?" can be *asked*; a real answer needs an electron energy
    balance, not an offset.
    """

    def excess(t: float) -> float:
        t_e = t + t_e_offset
        return magnetic_reynolds(sigma(t, rho, x_k, t_e=t_e), v_l) - 1.0

    lo, hi = t_lo, t_hi
    if excess(lo) > 0.0 or excess(hi) < 0.0:
        raise ValueError(
            f"Rm = 1 is not bracketed on [{t_lo:.0f}, {t_hi:.0f}] K at v*L = {v_l:g}: "
            f"Rm spans {magnetic_reynolds(sigma(lo, rho, x_k), v_l):.3g} to "
            f"{magnetic_reynolds(sigma(hi, rho, x_k), v_l):.3g}"
        )
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if excess(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def interpolated_cliff_temperature(
    rho: float,
    x_k: float,
    v_l: float,
    temps: tuple[float, ...] = SEED_WINDOW_TEMPS,
) -> float:
    """The cliff a reader gets from **the six tabulated `sigma` values alone** [K].

    Kept deliberately, and only for contrast: this is the wrong answer, and it is the wrong answer
    that `tab:seed_window` invites. Log-interpolating `sigma` between the table's rows and solving
    `Rm = 1` on the interpolant overstates the crossing by +51 to +117 K across the flown `v L`
    band, because the table samples every 1000 K while `sigma` climbs 60x between its first two
    rows -- the shape of the steepest part of the curve is simply not in it.

    The tabulated values themselves are not at fault: they *are* the `sigma` the solver uses, at the
    six temperatures the table reports. What fails is asking a six-point table for a crossing that
    lies inside its first interval. Use
    `cliff_temperature`, which bisects the continuous `sigma`; this exists so `main` can print the
    two side by side and nobody re-derives the interpolated one by hand.
    """
    log_sigma = [math.log(sigma(t, rho, x_k)) for t in temps]

    def interpolated(t: float) -> float:
        for i in range(len(temps) - 1):
            if temps[i] <= t <= temps[i + 1]:
                frac = (t - temps[i]) / (temps[i + 1] - temps[i])
                return math.exp(log_sigma[i] + frac * (log_sigma[i + 1] - log_sigma[i]))
        raise ValueError(f"{t:.0f} K is outside the tabulated window")

    lo, hi = temps[0], temps[-1]
    if magnetic_reynolds(interpolated(lo), v_l) > 1.0 or (
        magnetic_reynolds(interpolated(hi), v_l) < 1.0
    ):
        raise ValueError(f"Rm = 1 is not bracketed by the tabulated window at v*L = {v_l:g}")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if magnetic_reynolds(interpolated(mid), v_l) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def hall_parameter(
    temp: float,
    rho: float,
    x_k: float,
    b_field: float,
    t_e: float | None = None,
    q_en: float = Q_EN,
) -> float:
    """`beta = omega_c tau = e B / (m_e nu)` -- how far an electron gyrates between collisions.

    It is the control parameter for the **electrothermal (Velikhov) instability**: the Hall effect
    is what lets a local conductivity perturbation redirect the current and change the local Joule
    heating, which is the feedback that drives the runaway. Below `beta ~ 1` the electrons are
    collision-dominated, the feedback loop is broken, and a uniform two-temperature description is
    safe. Well above it, a seeded plasma can filament into hot streamers instead of heating evenly.

    There is no default `B`: this repository does not own the bag field, and inventing one would be
    the same error as inventing the `v L` product in `magnetic_reynolds`.
    """
    t_e = temp if t_e is None else t_e
    n_e = electron_density(temp, rho, x_k, t_e)
    nu = _nu_coulomb(n_e, t_e) + _nu_electron_neutral(temp, rho, x_k, t_e, q_en)
    return E_CHARGE * b_field / (M_E * nu)


def ionisation_sensitivity(
    temp: float, rho: float, x_k: float, t_e: float | None = None, rel_step: float = 1.0e-4
) -> float:
    """`S = d ln n_e / d ln T_e` -- the gain in the electrothermal feedback loop.

    The runaway is: a local rise in `n_e` concentrates the current (via the Hall effect), which
    raises the local Joule heating, which raises `T_e`, which raises `n_e`. `S` is the gain of the
    last link, and it is the term this repository can compute without any literature.

    **Every electron source must respond, not just the seed** (ADR-0038). The gain asks how the
    *whole* electron population answers a temperature perturbation, so it differentiates water's
    ionisation as well as potassium's. Differentiating only the seed -- which is what
    `electron_density` does, deliberately and correctly for `sigma` -- reports `S ~ 0` wherever
    water supplies the electrons, i.e. above ~7000 K, where the true gain is 5-7. That artifact
    read as "the loop is broken at the hot end" and it is not: the hot end is protected by a low
    Hall parameter instead.

    Three regimes, all falling out of the Saha balance already implemented:

    - **Weak seed ionisation:** `n_e ~ T^0.75 exp(-chi_K/2kT)`, so `S -> 3/4 + chi_K/(2 k T_e)`,
      which is ~11 at 2500 K. A gain that large is why seeded plasmas are unstable at all.
    - **Saturated seed, water not yet started (~5000-6500 K):** potassium cannot liberate more and
      water is still too cold, so `S` dips below 1. A genuine, narrow stabilising window.
    - **Water ionising:** `chi_H2O = 12.6-13.6 eV` against a hotter gas gives
      `chi/(2 k T_e) ~ 10` again, so the gain comes *back*, peaking near 7 around 10 000 K.

    Central difference in `ln T_e` rather than the analytic exponent, so it stays correct through
    both knees where the closed form does not apply.
    """
    t_e = temp if t_e is None else t_e
    hi = _n_e_all_sources_at(rho, x_k, t_e * (1.0 + rel_step))
    lo = _n_e_all_sources_at(rho, x_k, t_e * (1.0 - rel_step))
    return (math.log(hi) - math.log(lo)) / (2.0 * rel_step)


def _n_e_all_sources_at(rho: float, x_k: float, t_ion: float) -> float:
    """`n_e` with **every** ionising balance -- water's as well as the seed's -- driven by `t_ion`.

    `electron_density` holds water at the gas temperature and lets only the seed respond to `t_e`.
    That is the right contract for `sigma`, whose non-equilibrium use in this repository is the
    2000-5000 K band where the seed supplies the electrons. It is the wrong contract for the
    electrothermal gain, which is why this exists separately rather than as a flag on that one.

    **Approximation, and the reason this is private.** Perturbing `eos_water.composition` moves
    dissociation as well as ionisation, and dissociation is a heavy-particle process that should
    follow the gas rather than the electrons. It is second order for the gain -- dissociation
    changes the neutral count, not the electron count -- but it means the two-temperature value
    is an estimate. At `t_ion == temp` (the equilibrium case, which Q-M established is the one
    the nozzle actually runs) there is no approximation at all.
    """
    n_w = water_composition((1.0 - x_k) * rho, t_ion).n_e
    n_k = seed_number_density(rho, x_k)
    k_saha = math.exp(ln_k_saha(IP_K, G_K_ION, G_K_NEUTRAL, t_ion))
    b = k_saha - n_w
    return 0.5 * (-b + math.sqrt(b * b + 4.0 * k_saha * (n_w + n_k)))


#: Hall parameter above which the Hall-driven feedback is taken to be active. Engineering practice
#: for seeded MHD plasmas puts the critical value around 2; the sign change in the feedback occurs
#: near 1. **A screening threshold, not a dispersion relation** -- see `electrothermal_screen`.
BETA_CRIT = 2.0
#: Ionisation gain above which the `T_e -> n_e` link is taken to be live. Order-unity by
#: construction: at `S < 1` a fractional rise in `T_e` makes a smaller fractional rise in `n_e`, so
#: the loop attenuates.
SENSITIVITY_CRIT = 1.0


@dataclass(frozen=True)
class ElectrothermalScreen:
    """Screening result for the electrothermal (Velikhov) instability at one state."""

    temp: float
    b_field: float
    hall_parameter: float
    ionisation_sensitivity: float
    at_risk: bool


def electrothermal_screen(
    temp: float,
    rho: float,
    x_k: float,
    b_field: float,
    t_e: float | None = None,
) -> ElectrothermalScreen:
    """Screen a state for electrothermal-instability risk. **A screen, not a stability analysis.**

    The runaway needs both links of its loop: the Hall effect must be strong enough to turn a
    conductivity perturbation into a heating perturbation (`beta > BETA_CRIT`), and the ionisation
    must be sensitive enough to turn that heating into more electrons
    (`S > SENSITIVITY_CRIT`). Breaking either breaks the loop, so requiring both is a *necessary*
    condition and this function can only rule states **out**, never in.

    **What this repository can and cannot settle.** `beta` and `S` are computed exactly from the
    conductivity model. The threshold `BETA_CRIT = 2` is engineering practice for seeded MHD
    plasmas, taken on authority: the real criterion comes from a linearised dispersion relation
    (Velikhov 1962; Kerrebrock 1964) whose sources are not available here, and it depends on the
    degree of non-equilibrium, the seed fraction and the geometry. So a state flagged `at_risk`
    means "the loop is closed and the literature criterion must be checked", not "this filaments".

    Why it matters: if the instability triggers, the plasma breaks into hot streamers, the uniform
    two-temperature description behind Q-F fails, and the *effective* conductivity is **lower** than
    a smooth calculation predicts -- the opposite direction from the decoupling itself.
    """
    beta = hall_parameter(temp, rho, x_k, b_field, t_e)
    gain = ionisation_sensitivity(temp, rho, x_k, t_e)
    return ElectrothermalScreen(
        temp=temp,
        b_field=b_field,
        hall_parameter=beta,
        ionisation_sensitivity=gain,
        at_risk=beta > BETA_CRIT and gain > SENSITIVITY_CRIT,
    )


#: Vacuum permeability [H/m].
MU0 = 4.0e-7 * math.pi


def energy_relaxation_rate(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """`delta nu` [s^-1] -- the rate at which electrons hand their excess energy to the heavies.

    `delta = 2 m_e / M` is the fraction of an electron's energy lost per elastic collision, so the
    electron energy relaxes as `1/(delta nu)`. The Coulomb and neutral channels see different
    heavy masses, so each is weighted by its own collision rate rather than by a single mean `M`.

    Inelastic channels (water's rotation and vibration) are **not** included, and they are faster.
    Omitting them makes the relaxation slower, hence the electron temperature elevation *larger*
    and the loop gain in `electrothermal_loop` *larger* -- so the omission is conservative.
    """
    t_e = temp if t_e is None else t_e
    comp = water_composition((1.0 - x_k) * rho, temp)
    n_e = electron_density(temp, rho, x_k, t_e)
    nu_c = _nu_coulomb(n_e, t_e)
    nu_n = _nu_electron_neutral(temp, rho, x_k, t_e, Q_EN)

    n_o_ions = sum(comp.n_o_ions)
    n_kp = max(0.0, n_e - comp.n_hp - sum((j + 1) * n for j, n in enumerate(comp.n_o_ions)))
    ions = ((comp.n_hp, 1.008), (n_o_ions, 15.999), (n_kp, 39.0983))
    neutrals = (
        (comp.n_h2o, 18.015),
        (comp.n_oh, 17.007),
        (comp.n_h2, 2.016),
        (comp.n_o2, 31.998),
        (comp.n_h, 1.008),
        (comp.n_o, 15.999),
    )

    def weighted(pairs: tuple[tuple[float, float], ...]) -> float:
        total = sum(n for n, _ in pairs)
        if total <= 0.0:
            return 0.0
        return sum(n * (2.0 * M_E / (m * AMU)) for n, m in pairs) / total

    return weighted(ions) * nu_c + weighted(neutrals) * nu_n


@dataclass(frozen=True)
class InelasticChannel:
    """One electron-impact excitation channel: a threshold, a cross-section, and its target."""

    name: str
    #: Excitation threshold [J].
    threshold: float
    #: Near-threshold momentum-independent cross-section [m^2]. **The weak input**, like `Q_EN`:
    #: a single value stands in for an energy-dependent curve. Exposed so it can be swept.
    cross_section: float
    #: Which population the channel excites -- `"k_neutral"` or `"h2o"`.
    target: str


#: Electron-impact excitation channels, in the order they matter here.
#:
#: **K 4s-4p dominates**, which is the standard result for alkali-seeded plasmas: it is a resonance
#: transition with an oscillator strength near 1, so its cross-section is 2-3 orders above the
#: molecular channels, and its 1.61 eV threshold sits close enough to `k T_e` at 5000-8000 K to be
#: reachable without being negligible in energy.
#:
#: **The water channels are carried but have almost no target.** The plume is 99.7-99.9998%
#: dissociated across the band where any of this matters (Q-M), so `n_H2O` is 1e-6 to 3e-3 of the
#: heavies. Rotation is additionally suppressed by its tiny quantum: the net transfer scales as
#: `(dE/kT_e)^2` near equilibrium, and `dE ~ 0.005 eV` against `k T_e ~ 0.5 eV` costs four orders.
INELASTIC_CHANNELS: tuple[InelasticChannel, ...] = (
    InelasticChannel("K 4s-4p resonance", 1.61 * E_CHARGE, 3.0e-19, "k_neutral"),
    InelasticChannel("H2O stretch nu1/nu3", 0.453 * E_CHARGE, 3.0e-21, "h2o"),
    InelasticChannel("H2O bend nu2", 0.198 * E_CHARGE, 5.0e-21, "h2o"),
    InelasticChannel("H2O rotation", 0.005 * E_CHARGE, 1.0e-18, "h2o"),
)


def neutral_seed_density(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """Un-ionised potassium [m^-3] -- the target for the dominant inelastic channel.

    The seed's ions are whatever electrons water's own ionisation did not supply, so the neutral
    remainder follows without a second Saha solve.
    """
    t_e = temp if t_e is None else t_e
    n_k = seed_number_density(rho, x_k)
    n_water_e = water_composition((1.0 - x_k) * rho, temp).n_e
    n_kp = max(0.0, min(n_k, electron_density(temp, rho, x_k, t_e) - n_water_e))
    return n_k - n_kp


def inelastic_loss(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """Net electron energy loss to internal excitation [W/m^3]. Positive when `T_e > T_gas`.

    For each channel, with a step cross-section above threshold and a Maxwellian at `T_e`:

        k_exc = sigma_0 sqrt(8 k T_e / pi m_e) (1 + dE/k T_e) exp(-dE / k T_e)
        Q     = n_e n_target dE k_exc [1 - exp((dE/k)(1/T_e - 1/T_gas))]

    **The bracket is detailed balance and it is the whole story here.** The heavy species' internal
    states are kept at `T_gas` by heavy-heavy collisions, so electron-impact excitation and
    super-elastic de-excitation nearly cancel when `T_e` is close to `T_gas`; the bracket vanishes
    identically at `T_e == T_gas` and goes negative below it (electrons *gain*).

    The consequence is that inelastic losses **cannot rescue a near-equilibrium plasma from the
    Velikhov criterion**: expanding the bracket for small elevations gives `(dE/k) dT / T^2`, so
    the net inelastic loss is linear in the elevation exactly as the elastic channel is, and their
    ratio is a property of the state rather than something that grows as the plasma equilibrates.
    See ADR-0038.
    """
    t_e = temp if t_e is None else t_e
    n_e = electron_density(temp, rho, x_k, t_e)
    v_th = math.sqrt(8.0 * K_B * t_e / (math.pi * M_E))
    targets = {
        "k_neutral": neutral_seed_density(temp, rho, x_k, t_e),
        "h2o": water_composition((1.0 - x_k) * rho, temp).n_h2o,
    }
    total = 0.0
    for channel in INELASTIC_CHANNELS:
        n_target = targets[channel.target]
        if n_target <= 0.0:
            continue
        ratio = channel.threshold / (K_B * t_e)
        k_exc = channel.cross_section * v_th * (1.0 + ratio) * math.exp(-ratio)
        exponent = (channel.threshold / K_B) * (1.0 / t_e - 1.0 / temp)
        total += n_e * n_target * channel.threshold * k_exc * (1.0 - math.exp(exponent))
    return total


@dataclass(frozen=True)
class ElectronEnergyBalance:
    """Steady two-temperature state: Joule heating against elastic transfer to the heavies."""

    t_gas: float
    t_e: float
    sigma: float
    #: Thickness the driving current actually occupies [m] -- see `electron_energy_balance`.
    current_length_scale: float
    joule_heating: float
    converged: bool

    @property
    def elevation(self) -> float:
        """`T_e - T_gas` [K]."""
        return self.t_e - self.t_gas


def _current_length_scale(
    length_scale: float, transit_time: float, sigma_value: float, use_skin_depth: bool
) -> float:
    """Thickness the driving current occupies [m] -- the dominant lever on the verdict.

    `use_skin_depth=True` (shipped) takes the smaller of the field-gradient scale and the
    magnetic diffusion length `sqrt(t/(mu0 sigma))`, on the reasoning that where the field has
    not diffused across the flow the current sits in a skin and `j` is correspondingly larger.

    `use_skin_depth=False` puts the current across the full `length_scale` instead. That is the
    *physical* reading at this plume's low plasma beta (`2 mu0 p / B^2 ~ 0.016`): a flow that
    cannot appreciably distort the field does not concentrate its current into a resistive skin,
    and the diamagnetic current is distributed over the pressure-gradient scale.

    **This is not a tuning knob, it is ADR-0038 Addendum 3's open question.** A factor 3 in
    thickness is 9x in `Q_joule`, hence 9x in the elevation, hence 9x in `s` -- enough to flip
    every unstable station in this study. The skin choice is retained as the default because it
    is the conservative one, not because it is known to be right.
    """
    if not use_skin_depth:
        return length_scale
    return min(length_scale, math.sqrt(transit_time / (MU0 * sigma_value)))


def electron_energy_balance(
    temp: float,
    rho: float,
    x_k: float,
    b_field: float,
    length_scale: float,
    transit_time: float,
    use_skin_depth: bool = True,
    iterations: int = 200,
    relaxation: float = 0.3,
) -> ElectronEnergyBalance:
    """Solve `Q_joule = Q_loss` for the electron temperature (ADR-0038).

    Q-F(b) established that this balance is **algebraic, not an ODE**: the electron energy
    relaxes in ~1e-8 s against a ~1e-3 s transit, so `T_e` tracks the local balance instantly.

        Q_joule = j^2 / sigma,  j = |curl B| / mu0 ~ B / (mu0 L_eff)
        Q_loss  = (3/2) n_e k_B (delta nu) (T_e - T_gas)

    **The driving current is the field-gradient current, not `sigma u B`.** At the magnetic
    Reynolds numbers here (~40-900) the field is largely frozen to the plasma, so the plasma-frame
    field is small and the current is only what sustains the gradient. Taking `j = sigma u B`
    instead -- the low-`Rm` generator form -- overstates the heating by ~3 orders of magnitude and
    is the wrong regime; see ADR-0038.

    **`L_eff` is the smaller of the field-gradient scale and how far the field has diffused.**
    `sqrt(transit_time / (mu0 sigma))` is the magnetic diffusion length over the transit; where it
    is shorter than `length_scale` the current sits in a skin and `j` is correspondingly larger.
    Taking the minimum is conservative in the direction that matters.

    Solved by damped fixed point because the feedback is negative -- hotter electrons raise
    `sigma`, which both thickens the skin and cuts `j^2/sigma` -- so it contracts.

    `length_scale` and `transit_time` have no defaults for the same reason `hall_parameter` has no
    default `B`: this repository does not own the nozzle geometry.
    """
    if length_scale <= 0.0 or transit_time <= 0.0:
        raise ValueError("length_scale and transit_time must be positive")

    t_e = temp
    sigma_value = sigma(temp, rho, x_k, t_e)
    length_eff = _current_length_scale(length_scale, transit_time, sigma_value, use_skin_depth)
    joule = 0.0
    converged = False
    for _ in range(iterations):
        sigma_value = sigma(temp, rho, x_k, t_e)
        length_eff = _current_length_scale(length_scale, transit_time, sigma_value, use_skin_depth)
        joule = (b_field / (MU0 * length_eff)) ** 2 / sigma_value
        elastic_per_kelvin = (
            1.5
            * electron_density(temp, rho, x_k, t_e)
            * K_B
            * energy_relaxation_rate(temp, rho, x_k, t_e)
        )
        # Inelastic is not linear in the elevation, so it cannot be folded into a per-kelvin
        # coefficient; it is subtracted from the heating instead and the elastic channel closes.
        inelastic = inelastic_loss(temp, rho, x_k, t_e)
        target = temp + max(0.0, joule - inelastic) / elastic_per_kelvin
        if abs(target - t_e) < 1.0e-6 * t_e:
            t_e = target
            converged = True
            break
        t_e += relaxation * (target - t_e)
    return ElectronEnergyBalance(
        t_gas=temp,
        t_e=t_e,
        sigma=sigma_value,
        current_length_scale=length_eff,
        joule_heating=joule,
        converged=converged,
    )


#: Ionisation energy of the water-derived electron sources (H 13.598 eV, O 13.618 eV) [J].
IP_WATER = 13.6 * 1.602176634e-19


def mobility_sensitivity(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """`f = -(d mu/mu)/(d n_e/n_e)` -- how electron mobility answers a density perturbation.

    Petit and Geffray's `f`. It has two clean limits and this interpolates between them by which
    collision channel carries the momentum transfer:

    - **Coulomb-dominated:** `nu ~ n_e`, so `mu ~ 1/n_e` and `f -> 1`. Substituting into
      `beta_cr = 1.935 f + 0.065 + s` recovers their stated fully-ionised form `beta_cr ~ 2 + s`.
    - **Neutral-dominated:** `nu` is set by the neutral density and does not see `n_e`, so
      `f -> 0`.
    """
    t_e = temp if t_e is None else t_e
    n_e = electron_density(temp, rho, x_k, t_e)
    nu_c = _nu_coulomb(n_e, t_e)
    nu_n = _nu_electron_neutral(temp, rho, x_k, t_e, Q_EN)
    total = nu_c + nu_n
    return nu_c / total if total > 0.0 else 0.0


def ionisation_energy(temp: float, rho: float, x_k: float, t_e: float | None = None) -> float:
    """`E_i` [J] -- the ionisation energy of whichever species is actually supplying electrons.

    The seed at 4.34 eV below ~6000 K, water's 13.6 eV above it, blended by electron share. `E_i`
    enters `critical_hall_parameter` inversely, so using water's value where the seed is supplying
    the electrons would understate `s` by ~3x and wrongly call a stable state unstable.
    """
    t_e = temp if t_e is None else t_e
    n_e = electron_density(temp, rho, x_k, t_e)
    if n_e <= 0.0:
        return IP_K
    water_share = max(0.0, 1.0 - min(seed_number_density(rho, x_k), n_e) / n_e)
    return (1.0 - water_share) * IP_K + water_share * IP_WATER


def critical_hall_parameter(t_e: float, t_gas: float, e_i: float, f: float) -> float:
    """`beta_cr` for the Velikhov instability -- Petit and Geffray (2009), eqns p. 1170.

        s       = 2 k T_e^2 / [E_i (T_e - T_gas)] * 1 / (1 + 1.5 k T_e / E_i)
        beta_cr = 1.935 f + 0.065 + s

    **`(T_e - T_gas)` is in the denominator, so `beta_cr` diverges as the plasma approaches
    thermal equilibrium.** That is the whole reason `BETA_CRIT = 2` does not transfer here: 2 is
    the `s -> 0` limit of their fully-ionised form, i.e. a *strongly* two-temperature plasma. The
    published criterion is explicitly for `T_e > T_gas` regimes.

    Returns infinity at `T_e <= T_gas`: with no electron heating there is no drive, and the
    instability cannot run at any Hall parameter.
    """
    if t_e <= t_gas:
        return math.inf
    s = (2.0 * K_B * t_e * t_e) / (e_i * (t_e - t_gas)) / (1.0 + 1.5 * K_B * t_e / e_i)
    return 1.935 * f + 0.065 + s


@dataclass(frozen=True)
class ElectrothermalLoop:
    """Velikhov stability at one state, against the published criterion rather than a constant."""

    screen: ElectrothermalScreen
    balance: ElectronEnergyBalance
    ionisation_energy: float
    mobility_sensitivity: float
    critical_hall_parameter: float
    #: Linear growth rate [s^-1], zero when `beta <= beta_cr`. Petit and Geffray's `g`.
    growth_rate: float
    unstable: bool

    @property
    def e_folding_time(self) -> float:
        """Time for the perturbation to grow by `e` [s]; infinite when stable."""
        return math.inf if self.growth_rate <= 0.0 else 1.0 / self.growth_rate


def electrothermal_loop(
    temp: float,
    rho: float,
    x_k: float,
    b_field: float,
    length_scale: float,
    transit_time: float,
    use_skin_depth: bool = True,
) -> ElectrothermalLoop:
    """Velikhov stability at one state, against Petit and Geffray (2009) (ADR-0038).

        g = sigma E*^2 / [n_e (E_i + 1.5 k T_e)(1 + beta^2)] * (beta - beta_cr)

    Unlike `electrothermal_screen`, this **can rule states in as well as out**, because it is the
    published linearised criterion rather than a two-link necessary condition against a constant.
    Three ingredients the screen does not have:

    1. `T_e` from the self-consistent `electron_energy_balance`, which is what sets `beta_cr`;
    2. `beta_cr` from the criterion itself, which diverges near equilibrium;
    3. `E*` -- the field actually driving the current, `j/sigma`, taken from the same
       field-gradient current the balance uses, **not** from `u x B`.

    Returns `growth_rate = 0` when `beta <= beta_cr`. A nonzero rate is a linear growth rate, so
    the thing to compare against is the *residence time*: this plasma's relaxation is microseconds
    against a millisecond transit, so a state that is unstable at all is unstable many times over.
    """
    balance = electron_energy_balance(
        temp, rho, x_k, b_field, length_scale, transit_time, use_skin_depth
    )
    screen = electrothermal_screen(temp, rho, x_k, b_field, balance.t_e)
    e_i = ionisation_energy(temp, rho, x_k, balance.t_e)
    f = mobility_sensitivity(temp, rho, x_k, balance.t_e)
    beta_cr = critical_hall_parameter(balance.t_e, temp, e_i, f)

    beta = screen.hall_parameter
    growth = 0.0
    if math.isfinite(beta_cr) and beta > beta_cr:
        n_e = electron_density(temp, rho, x_k, balance.t_e)
        e_star = b_field / (MU0 * balance.current_length_scale) / balance.sigma
        denom = n_e * (e_i + 1.5 * K_B * balance.t_e) * (1.0 + beta * beta)
        growth = balance.sigma * e_star * e_star / denom * (beta - beta_cr)
    return ElectrothermalLoop(
        screen=screen,
        balance=balance,
        ionisation_energy=e_i,
        mobility_sensitivity=f,
        critical_hall_parameter=beta_cr,
        growth_rate=growth,
        unstable=growth > 0.0,
    )


@dataclass(frozen=True)
class SeedWindowRow:
    """One row of the regenerated `tab:seed_window`."""

    temp: float
    n_e: float
    ionised_fraction: float
    sigma: float
    rm: float
    leak_fraction: float
    #: `S = d ln n_e / d ln T_e`, the electrothermal feedback gain.
    ionisation_sensitivity: float
    #: Field [T] needed to reach `BETA_CRIT`. Reported instead of a Hall parameter at an assumed
    #: `B`, because this repository does not own the bag field: it is the number the field must be
    #: compared against, and it needs no assumption to compute.
    b_field_for_beta_crit: float


def seed_window(
    rho: float,
    x_k: float,
    v_l: float,
    temps: Sequence[float] = SEED_WINDOW_TEMPS,
    t_e_offset: float = 0.0,
) -> list[SeedWindowRow]:
    """Regenerate `tab:seed_window` as an output rather than a hand table.

    `leak_fraction ~ 1/Rm` is item 10's bracket: the fraction of the field that diffuses out during
    one expansion time. It is capped at 1 -- below `Rm = 1` the field is simply not held, and
    quoting "a leak of 12x" would be reporting the failure of the approximation as a number.
    """
    out: list[SeedWindowRow] = []
    for t in temps:
        t_e = t + t_e_offset
        n_e = electron_density(t, rho, x_k, t_e)
        n_k = seed_number_density(rho, x_k)
        s = sigma(t, rho, x_k, t_e=t_e)
        rm = magnetic_reynolds(s, v_l)
        out.append(
            SeedWindowRow(
                temp=t,
                n_e=n_e,
                ionised_fraction=min(n_e / n_k, 1.0),
                sigma=s,
                rm=rm,
                leak_fraction=min(1.0 / rm, 1.0),
                ionisation_sensitivity=ionisation_sensitivity(t, rho, x_k, t_e),
                b_field_for_beta_crit=BETA_CRIT / hall_parameter(t, rho, x_k, 1.0, t_e),
            )
        )
    return out


def write_seed_window(rows: list[SeedWindowRow], path: Path = DEFAULT_SEED_WINDOW_PATH) -> None:
    """Write the regenerated table as CSV."""
    lines = [
        "T_K,n_e_m3,ionised_fraction,sigma_S_per_m,Rm,leak_fraction,"
        "ionisation_sensitivity,B_T_for_beta_crit"
    ]
    lines += [
        f"{r.temp:.0f},{r.n_e:.6e},{r.ionised_fraction:.6f},"
        f"{r.sigma:.6e},{r.rm:.6e},{r.leak_fraction:.6f},"
        f"{r.ionisation_sensitivity:.4f},{r.b_field_for_beta_crit:.4f}"
        for r in rows
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _report_cliff_band() -> None:
    """Print the conductivity cliff across the flown `v L` band, solved beside interpolated.

    Both rows are printed on purpose (deferred item D9). The paper quoted 2570 K for years because
    it came from the interpolated row; printing the two together is what stops that recurring.
    """
    lo, hi = V_L_BAND
    products = ((RETIRED_V_L, "retired"), (lo, "band lo"), (REF_V_L, "stated"), (hi, "band hi"))
    print(f"  cliff (Rm = 1) at rho = {REF_RHO} kg/m^3, x_K = {REF_X_K}:")
    print(f"    {'v L [m^2/s]':>12} {'solved':>9} {'interpolated':>13} {'error':>7}  leg")
    for v_l, label in products:
        try:
            solved = cliff_temperature(REF_RHO, REF_X_K, v_l)
            guess = interpolated_cliff_temperature(REF_RHO, REF_X_K, v_l)
        except ValueError as exc:
            print(f"    {v_l:12.3g} no cliff in range: {exc}")
            continue
        print(f"    {v_l:12.3g} {solved:8.0f}K {guess:12.0f}K {guess - solved:+6.0f}K  {label}")
    print(
        "    solved bisects the continuous sigma; interpolated log-interpolates the six rows "
        "above and is wrong:\n    the crossing lies inside the table's first interval. "
        "Do not quote the interpolated column."
    )


def main() -> None:
    """Regenerate the seed window and report the cliff, at the solved `v L` and across its band."""
    rows = seed_window(REF_RHO, REF_X_K, REF_V_L)
    write_seed_window(rows)
    print(
        f"python: seed window at rho = {REF_RHO} kg/m^3, x_K = {REF_X_K}, v*L = {REF_V_L:g} m^2/s"
    )
    header = f"  {'T [K]':>7} {'n_e [m^-3]':>12} {'ionised':>8} {'sigma [S/m]':>12}"
    print(f"{header} {'Rm':>10} {'leak':>7} {'gain S':>8} {'B@beta=2':>9}")
    for r in rows:
        print(
            f"  {r.temp:7.0f} {r.n_e:12.3e} {r.ionised_fraction:8.4f} "
            f"{r.sigma:12.4g} {r.rm:10.4g} {r.leak_fraction:7.3f} "
            f"{r.ionisation_sensitivity:8.2f} {r.b_field_for_beta_crit:8.2f}T"
        )
    _report_cliff_band()
    print(
        "  electrothermal screen (superseded, ADR-0038): the two-link form needs beta > "
        f"{BETA_CRIT:g} AND S > {SENSITIVITY_CRIT:g}, but {BETA_CRIT:g} is the strongly "
        "two-temperature limit and does not transfer to an equilibrium plume."
    )
    print(
        "  the criterion that decides it is `electrothermal_loop` -- see `make "
        "analysis-electrothermal`, which walks it along the cooling history."
    )
    print(f"python: wrote {DEFAULT_SEED_WINDOW_PATH}")


if __name__ == "__main__":
    main()
