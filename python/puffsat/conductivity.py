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
    the water solve (H2O, H, O) plus un-ionised potassium.
    """
    comp = water_composition((1.0 - x_k) * rho, temp)
    n_neutral = comp.n_h2o + comp.n_h + comp.n_o
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
#: Reference conditions. `rho` and `x_K` are the paper's slug bag values; `v_l` is **not** the
#: paper's -- it never states one (audit Q-G) -- and is the product implied by back-solving its own
#: `Rm = 361` row at 15 000 K against the audit's hand-computed conductivity. Stated, not hidden.
REF_RHO = 0.32
REF_X_K = 0.01
REF_V_L = 1.81e4
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
    makes the column reproducible. The caller must say what expansion it means.
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

    Two limits matter and both fall out of the Saha balance already implemented:

    - **Weak seed ionisation:** `n_e ~ T^0.75 exp(-chi/2kT)`, so `S -> 3/4 + chi/(2 k T_e)`, which
      is ~11 at 2500 K. A gain that large is why seeded plasmas are unstable at all.
    - **Saturated seed:** every potassium atom is ionised, so `T_e` cannot liberate more and `S`
      collapses. The runaway is choked at source, independently of the Hall parameter.

    Central difference in `ln T_e` rather than the analytic exponent, so it stays correct through
    the saturation knee where the closed form does not apply.
    """
    t_e = temp if t_e is None else t_e
    hi = electron_density(temp, rho, x_k, t_e * (1.0 + rel_step))
    lo = electron_density(temp, rho, x_k, t_e * (1.0 - rel_step))
    return (math.log(hi) - math.log(lo)) / (2.0 * rel_step)


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


def main() -> None:
    """Regenerate the seed window and report the cliff, for a stated `v L`."""
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
    try:
        print(f"  cliff (Rm = 1) at {cliff_temperature(REF_RHO, REF_X_K, REF_V_L):.0f} K")
    except ValueError as exc:
        print(f"  no cliff in range: {exc}")
    print(
        "  electrothermal screen: the loop needs beta > "
        f"{BETA_CRIT:g} AND S > {SENSITIVITY_CRIT:g}; the gain collapses above ~5000 K, and at "
        "this density the Hall link needs several tesla."
    )
    print(f"python: wrote {DEFAULT_SEED_WINDOW_PATH}")


if __name__ == "__main__":
    main()
