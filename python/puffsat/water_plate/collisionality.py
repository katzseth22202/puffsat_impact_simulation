"""Does the arriving PuffSat gas collide with the argon spray, or stream through it? (R3)

**It collides, by about six decades.** The spray cup's `eta_jet` comes from fluid solves, which
assume the two gases meet at a contact and shock. That needs the arriving atoms to lose their
momentum in a short distance compared with the ~4 m spray cloud. `continuum.py` checks the slug
on thermal cross-sections. Those do not apply here: the gases close at 46-68 km/s, i.e. 10-480 eV
in the centre-of-mass frame, where cross-sections shrink.

# Cross-sections at the relative velocity

The PuffSat arrives as cold neutral water and the spray as cold neutral argon, so the first
collisions are neutral-neutral at the closing speed. A 100-300 eV collision breaks the 5 eV O-H
bond, so the arriving water is treated as its atoms, O and H, each meeting Ar. A whole molecule
is larger than its atoms, so splitting it is generous.

The interaction is the **ZBL universal screened-Coulomb potential** (Ziegler, Biersack and
Littmark, *The Stopping and Range of Ions in Solids*, 1985), the standard repulsive potential
for atom-atom collisions at tens of eV and up. The classical deflection angle is integrated
exactly at each impact parameter. The momentum-transfer cross-section is
`sigma_mt = 2 pi int (1 - cos chi) b db`. Two values are reported:

- `sigma_mt`, the cross-section that sets stopping.
- `pi r0^2`, the head-on distance of closest approach squared. That is the **hard-sphere value at
  the relative velocity** and the smaller of the two, so it is the generous case. This is the
  same choice `continuum.py` makes.

ZBL is fitted to within a few per cent at short range and is less sure at the 1-2 Angstrom
separations reached here. A factor of 2 either way does not touch a six-decade verdict.

# Ion channels, alongside

Once the gas ionizes, two channels open:

- **Ion-neutral**, Ar+ on Ar at ~80 eV: `sigma_m` 0.7e-18 and `sigma_CT` 0.3e-18 m^2 (Phelps
  1991, as used by Merritt et al. 2014). That is resonant charge exchange. O+ on Ar is not
  resonant, so this is an upper reference, not our pair.
- **Ion-ion Coulomb drag** at the relative velocity, from the NRL fast-test-particle slowing
  rate. This is the rate Merritt et al. 2014 use (their Eqs. 5-7), with their factor of 4 for
  slowing to rest.

Ion-neutral is ~50x shorter than the neutral path. Coulomb drag is **longer**: its cross-section
falls as `v^-4`, so at 50-70 km/s and spray densities a fully ionized spray stops an O+ in 4-23
um, against 2-4 um for neutrals. That is still five decades under the cloud depth. Ionization
does not rescue interpenetration, and it does not help the verdict either.

# The leading edge

The arriving gas's own low-density front is where a fluid model might fail. Kn = 0.1 against the
cloud depth fixes a threshold density. How much PuffSat mass arrives below it depends on the
delivered shape, which is a design input (necklace charges, ADR-0055). So it is bracketed:

- an exponential tail, the fattest one an isothermal expansion makes;
- a 1-D Gaussian and a 3-D Gaussian, from ballistic spreading.

That threshold is the PuffSat's **own** collisionality. A leading-edge atom entering the spray
still meets argon at the spray's density, and stops within millimetres. It is deposited, not
streamed through. Only the spray's column matters for interpenetration, and `k` fixes that column.

# Droplets

All of the above assumes the spray is vapour. Liquid droplets of radius `r` are spaced so that
an atom crosses `4 r rho_liquid / (3 rho_spray)` before it hits one. That is the one way the
spray can stay porous. The droplet radius at which it reaches Kn = 0.1 is reported.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

OUTPUT = Path("data/results/water_plate/collisionality.csv")

AMU = 1.66053906660e-27
EV = 1.602176634e-19
BOHR = 0.529177210903e-10
COULOMB_J_M = 1.439964547e-9 * EV
"""`e^2 / (4 pi eps0)` [J m]."""

M_AR = 39.948
M_O = 15.999
M_H = 1.008
M_H2O = 18.015
Z_AR, Z_O, Z_H = 18, 8, 1

# ---- Design point (ADR-0055 amended, Q59-Q62) -------------------------------------------------

PUFFSAT_KG_AT_45 = 92.0
"""PuffSat per 12 MN s pulse at 45.58 km/s."""
K_SPRAY = 8.52
FOOTPRINT_RADIUS_M = 5.0
CLOUD_DEPTH_M = 4.0
"""Spray-cloud depth and arriving-pulse length, both ~4 m in the 2-D shape runs."""
RHO_LIQUID_ARGON = 1395.0
"""Liquid argon near its boiling point [kg/m^3]."""

SPEEDS_M_S = (45_580.0, 56_530.0, 61_830.0, 65_130.0, 68_000.0)
"""The four overtake speeds, plus 68 km/s as the top of the request's closing range."""

KN_THRESHOLD = 0.1
"""Kn above which a fluid model is in doubt, as in `continuum.py`."""

# Ion-neutral reference (Phelps 1991 via Merritt et al. 2014, Ar+ on Ar at ~80 eV) [m^2].
SIGMA_ARPLUS_AR_MOMENTUM = 0.7e-18
SIGMA_ARPLUS_AR_CHARGE_EXCHANGE = 0.3e-18


# ---- The ZBL potential and classical scattering ----------------------------------------------


def zbl_screening(x: NDArray[np.float64] | float) -> NDArray[np.float64] | float:
    """ZBL universal screening function `phi(r / a_U)`."""
    return (
        0.1818 * np.exp(-3.2 * x)
        + 0.5099 * np.exp(-0.9423 * x)
        + 0.2802 * np.exp(-0.4029 * x)
        + 0.02817 * np.exp(-0.2016 * x)
    )


@dataclass(frozen=True)
class Potential:
    """A central repulsive potential `V(r)` [J]: ZBL for a pair `(z1, z2)`, or bare Coulomb.

    `screened = False` gives `z1 z2 e^2 / r`, which has a closed-form deflection angle and is
    the test case for the integrator.
    """

    z1: int
    z2: int
    screened: bool = True

    @property
    def screening_length(self) -> float:
        """ZBL universal screening length `a_U` [m]."""
        return 0.8854 * BOHR / (math.pow(self.z1, 0.23) + math.pow(self.z2, 0.23))

    def __call__(self, r: NDArray[np.float64] | float) -> NDArray[np.float64] | float:
        coulomb = self.z1 * self.z2 * COULOMB_J_M / r
        if not self.screened:
            return coulomb
        return coulomb * zbl_screening(r / self.screening_length)


_GL_X, _GL_W = np.polynomial.legendre.leggauss(256)
_S = 0.5 * (_GL_X + 1.0)
_W = 0.5 * _GL_W


def closest_approach(b: float, energy: float, pot: Potential) -> float:
    """Root `r0` of `1 - V(r)/E - b^2/r^2 = 0` [m], by bisection in `ln r`."""

    def f(r: float) -> float:
        return 1.0 - float(pot(r)) / energy - (b / r) ** 2

    lo, hi = 1e-16, 2.0 * b + 1e-9
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return hi


def deflection(b: float, energy: float, pot: Potential) -> float:
    """Centre-of-mass deflection angle `chi(b)` [rad] at CM energy `energy` [J].

    `chi = pi - 2 (b/r0) int_0^1 du / sqrt(1 - V(r0/u)/E - (b u/r0)^2)`. The integrable
    singularity at `u = 1` is removed by `u = 1 - s^2`, then Gauss-Legendre in `s`.
    """
    r0 = closest_approach(b, energy, pot)
    u = 1.0 - _S**2
    g = 1.0 - np.asarray(pot(r0 / u)) / energy - (b * u / r0) ** 2
    integral = float(np.sum(_W * 2.0 * _S / np.sqrt(np.maximum(g, 1e-300))))
    return math.pi - 2.0 * b / r0 * integral


def momentum_transfer_cross_section(energy: float, pot: Potential) -> float:
    """`sigma_mt = 2 pi int (1 - cos chi) b db` [m^2], on a log grid in `b`."""
    bs = np.logspace(-14.0, -8.5, 800)
    one_minus_cos = np.array([1.0 - math.cos(deflection(float(b), energy, pot)) for b in bs])
    return float(2.0 * math.pi * np.trapezoid(one_minus_cos * bs**2, np.log(bs)))


def hard_sphere_at_speed(energy: float, pot: Potential) -> float:
    """`pi r0(b = 0)^2` [m^2]: the hard-sphere size at this collision energy."""
    return math.pi * closest_approach(0.0, energy, pot) ** 2


# ---- Coulomb drag (NRL fast test particle; Merritt et al. 2014 Eqs. 5-7) --------------------


def coulomb_stopping_length(
    v_rel: float, m_test: float, m_field: float, n_field: float, n_e: float, t_e_ev: float
) -> float:
    """Penetration length [m] of a singly charged ion into a singly ionized field plasma.

    `nu = 9.0e-8 n' lnL (1/mu + 1/mu') mu^0.5 / eps^1.5` with `n'` in cm^-3, masses in amu and
    the test particle's energy `eps` in eV. `lnL = 43 - ln[(mu+mu')/(mu mu' beta^2) (ne/Te)^0.5]`
    for counter-streaming ions in warm electrons. The length is `v_rel / (4 nu)`.
    """
    n_cm = n_field * 1e-6
    ne_cm = n_e * 1e-6
    beta2 = (v_rel / 2.99792458e8) ** 2
    ln_lambda = 43.0 - math.log(
        (m_test + m_field) / (m_test * m_field * beta2) * math.sqrt(ne_cm / t_e_ev)
    )
    eps_ev = 0.5 * m_test * AMU * v_rel**2 / EV
    nu = (
        9.0e-8 * n_cm * ln_lambda * (1.0 / m_test + 1.0 / m_field) * math.sqrt(m_test) / eps_ev**1.5
    )
    return float(v_rel / (4.0 * nu))


# ---- Leading-edge mass fractions ------------------------------------------------------------


def fraction_below_exponential(density_ratio: float) -> float:
    """Mass fraction below `rho*/rho0` in an exponential tail `rho0 exp(-x/h)`: just the ratio."""
    return min(1.0, density_ratio)


def fraction_below_gaussian_1d(density_ratio: float) -> float:
    """Mass fraction below `rho*/rho0` in a 1-D Gaussian: `erfc(sqrt(ln(rho0/rho*)))`."""
    if density_ratio >= 1.0:
        return 1.0
    return math.erfc(math.sqrt(math.log(1.0 / density_ratio)))


def fraction_below_gaussian_3d(density_ratio: float) -> float:
    """Mass fraction below `rho*/rho0` in a spherical Gaussian: the chi-squared(3) tail."""
    if density_ratio >= 1.0:
        return 1.0
    x = 2.0 * math.log(1.0 / density_ratio)
    return math.erfc(math.sqrt(x / 2.0)) + math.sqrt(2.0 * x / math.pi) * math.exp(-x / 2.0)


# ---- The design point -----------------------------------------------------------------------


def reduced_mass(m1: float, m2: float) -> float:
    return m1 * m2 / (m1 + m2)


def cm_energy(m1: float, m2: float, v_rel: float) -> float:
    """Centre-of-mass collision energy [J] for masses in amu."""
    return 0.5 * reduced_mass(m1, m2) * AMU * v_rel**2


@dataclass(frozen=True)
class SpeedRow:
    """Collisionality of PuffSat water meeting the argon spray, at one closing speed."""

    speed_km_s: float
    puffsat_kg: float
    ecm_o_ar_ev: float
    ecm_h_ar_ev: float
    sigma_mt_o_ar: float
    sigma_mt_h_ar: float
    sigma_hs_o_ar: float
    sigma_hs_h_ar: float
    n_argon: float
    argon_column: float
    n_puffsat: float
    stop_o_m: float
    stop_h_m: float
    kn_neutral_generous: float
    relaxation_lengths_in_column: float
    stop_ion_neutral_m: float
    stop_coulomb_m: float
    n_star: float
    rho_star_ratio: float
    frac_exponential: float
    frac_gauss_1d: float
    frac_gauss_3d: float
    droplet_radius_kn01_m: float


def evaluate(speed: float, t_e_ev: float = 3.0) -> SpeedRow:
    """Every number the request asks for, at one closing speed."""
    puffsat_kg = PUFFSAT_KG_AT_45 * 45_580.0 / speed  # fixed 12 MN s at fixed eta
    argon_kg = K_SPRAY * puffsat_kg
    area = math.pi * FOOTPRINT_RADIUS_M**2
    n_ar = argon_kg / (area * CLOUD_DEPTH_M) / (M_AR * AMU)
    column = argon_kg / area / (M_AR * AMU)
    rho_puff = puffsat_kg / (area * CLOUD_DEPTH_M)
    n_puff = rho_puff / (M_H2O * AMU)

    pot_o, pot_h = Potential(Z_O, Z_AR), Potential(Z_H, Z_AR)
    e_o, e_h = cm_energy(M_O, M_AR, speed), cm_energy(M_H, M_AR, speed)
    smt_o, smt_h = (
        momentum_transfer_cross_section(e_o, pot_o),
        momentum_transfer_cross_section(e_h, pot_h),
    )
    shs_o, shs_h = hard_sphere_at_speed(e_o, pot_o), hard_sphere_at_speed(e_h, pot_h)

    # Momentum relaxation length of a test atom: (m1/mu) / (n sigma_mt). The cross-section grows
    # as the atom slows, so taking it at the arrival speed overstates the length.
    stop_o = (M_O / reduced_mass(M_O, M_AR)) / (n_ar * smt_o)
    stop_h = (M_H / reduced_mass(M_H, M_AR)) / (n_ar * smt_h)
    sigma_generous = min(smt_o, smt_h, shs_o, shs_h)
    kn = max(stop_o, stop_h, 1.0 / (n_ar * sigma_generous)) / CLOUD_DEPTH_M
    relax = column * smt_o * reduced_mass(M_O, M_AR) / M_O

    stop_in = 1.0 / (n_ar * (SIGMA_ARPLUS_AR_MOMENTUM + SIGMA_ARPLUS_AR_CHARGE_EXCHANGE))
    stop_c = coulomb_stopping_length(speed, M_O, M_AR, n_ar, n_ar, t_e_ev)

    n_star = 1.0 / (KN_THRESHOLD * CLOUD_DEPTH_M * sigma_generous)
    ratio = n_star / n_puff
    rho_ar = argon_kg / (area * CLOUD_DEPTH_M)
    r_drop = KN_THRESHOLD * CLOUD_DEPTH_M * 3.0 * rho_ar / (4.0 * RHO_LIQUID_ARGON)

    return SpeedRow(
        speed_km_s=speed / 1e3,
        puffsat_kg=puffsat_kg,
        ecm_o_ar_ev=e_o / EV,
        ecm_h_ar_ev=e_h / EV,
        sigma_mt_o_ar=smt_o,
        sigma_mt_h_ar=smt_h,
        sigma_hs_o_ar=shs_o,
        sigma_hs_h_ar=shs_h,
        n_argon=n_ar,
        argon_column=column,
        n_puffsat=n_puff,
        stop_o_m=stop_o,
        stop_h_m=stop_h,
        kn_neutral_generous=kn,
        relaxation_lengths_in_column=relax,
        stop_ion_neutral_m=stop_in,
        stop_coulomb_m=stop_c,
        n_star=n_star,
        rho_star_ratio=ratio,
        frac_exponential=fraction_below_exponential(ratio),
        frac_gauss_1d=fraction_below_gaussian_1d(ratio),
        frac_gauss_3d=fraction_below_gaussian_3d(ratio),
        droplet_radius_kn01_m=r_drop,
    )


# ---- PLX benchmark --------------------------------------------------------------------------

PLX_CASES = (
    # (label, stopping length [m], interaction scale [m], observed outcome)
    ("Merritt 2014, oblique, 100% Ar", 0.018, 0.03, "collisional stagnation layer"),
    ("Merritt 2014, oblique, Ar/O/Al mix", 0.008, 0.03, "collisional stagnation layer"),
    ("Moser & Hsu 2015, head-on, t = 35 us", 10.0, 0.25, "interpenetration"),
    ("Moser & Hsu 2015, head-on, t = 40 us", 1.4, 0.29, "transition to stagnation"),
)
"""Inter-jet stopping length against the observed interaction width (Merritt et al., Phys.
Plasmas 21, 055703 (2014), Table II, few-cm layer; Moser and Hsu, Phys. Plasmas 22, 055707
(2015), Table I, 21-30 cm width). The experimental boundary between interpenetration and
stagnation sits at stopping length / width of order 1-5."""


def main() -> None:
    rows = [evaluate(v) for v in SPEEDS_M_S]

    print("== R3: does PuffSat water collide with the argon spray? ==")
    print("   (k = 8.52, 12 MN s, 4 m cloud)\n")
    print(
        f"{'v km/s':>7} {'Ecm O-Ar':>9} {'s_mt O-Ar':>10} {'s_mt H-Ar':>10} {'piR0^2 O':>10} "
        f"{'n_Ar m-3':>10} {'stop O':>9} {'Kn(gen)':>9} {'relax lens':>10}"
    )
    for r in rows:
        print(
            f"{r.speed_km_s:7.2f} {r.ecm_o_ar_ev:8.0f}eV {r.sigma_mt_o_ar:10.2e} "
            f"{r.sigma_mt_h_ar:10.2e} {r.sigma_hs_o_ar:10.2e} {r.n_argon:10.2e} "
            f"{r.stop_o_m * 1e6:7.2f}um {r.kn_neutral_generous:9.1e} "
            f"{r.relaxation_lengths_in_column:10.1e}"
        )
    print("\nIon channels at the same density (ionized spray, Te = 3 eV):")
    for r in rows:
        print(
            f"  {r.speed_km_s:6.2f} km/s  ion-neutral {r.stop_ion_neutral_m * 1e6:6.3f} um   "
            f"Coulomb drag {r.stop_coulomb_m * 1e6:8.3f} um"
        )
    print("\nLeading edge of the arriving PuffSat (Kn = 0.1 against 4 m):")
    for r in rows:
        print(
            f"  {r.speed_km_s:6.2f} km/s  n* = {r.n_star:.2e} m^-3"
            f" = {r.rho_star_ratio:.1e} of mean;"
            f" mass below: exp {r.frac_exponential:.1e}, 1-D Gauss {r.frac_gauss_1d:.1e},"
            f" 3-D Gauss {r.frac_gauss_3d:.1e}"
        )
    print("\nDroplets: liquid-argon radius at which the spray is porous at Kn = 0.1:")
    for r in rows:
        print(f"  {r.speed_km_s:6.2f} km/s  r = {r.droplet_radius_kn01_m * 1e3:.2f} mm")
    print("\nPLX benchmark (stopping length / interaction width -> outcome):")
    for label, stop, width, outcome in PLX_CASES:
        print(f"  {label:38s} {stop / width:6.2f}  {outcome}")
    worst = max(r.kn_neutral_generous for r in rows)
    decades = math.log10(1.0 / worst)
    print(f"\nThis study, generous neutral case: {worst:.1e}, about {decades:.1f} decades")
    print("below the PLX stagnation boundary.")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(SpeedRow.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([f"{getattr(r, f):.6g}" for f in fields])
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
