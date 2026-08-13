"""Rung 0 — the analytic reference ledger for the tamped-nozzle study (PRD §10).

**One calculator owns every closed-form number the PRD quotes**, so the PRD, the ADRs, and the
eventual analysis cannot drift apart. Downstream rungs are permitted to quote closed-form figures
only from here. There is no solver in this module: it is exact algebra over the ballistic model
this study exists to replace, plus the thermodynamic ceiling it scores everything against.

**This is the tamped-nozzle study, not the `f(v)` study.** See `puffsat.tamper.__init__`.

---

## The model, in one pass

A projectile of mass `m_i` arrives head-on at closing speed `w`, passes through the plate's vertex
hole, and merges completely inelastically with a slug of `k` kg per projectile kg. Per projectile
kg the merged blob is `(1+k)`, its centre of mass recedes from the plate at `V = w/(1+k)`, and the
merge dissipates `(1/2)w^2 k/(1+k)` of specific energy — which the ballistic model returns *in
full* as isotropic expansion at `u = w sqrt(k)/(1+k)` in the blob frame. So the ejecta's kinetic
energy is exactly the projectile's incoming kinetic energy, and the blob's momentum is exactly the
incoming debit: with no plate, `J = 0` identically (PRD §3.6).

Everything downstream is bookkeeping on that fireball:

- an element emitted at blob-frame polar angle `theta` (from the plate-directed `+z` axis) has
  vehicle-frame velocity `v_z = u cos(theta) - V`, `v_r = u sin(theta)`;
- it reaches the plate only if it outruns the recoil, `cos(theta) > 1/sqrt(k)` — the **ballistic
  capture fraction** `max[0, (1 - 1/sqrt(k))/2]`, zero at `k <= 1`;
- the plate turns what it catches, and `J = sum over captured mass of Delta_p`, with
  `Delta_p = |v| + v_z` for a focus-matched paraboloid and `Delta_p = 2 v_z` for a flat plate.

Impulse is reported dimensionlessly as `beta = J/(m_i w)`, **including the projectile momentum
debit**, and the two deliverables are effective Isp `beta w/(g0 C)` and projectile economy
`beta w` (PRD §3.1). Nothing here is a physical velocity except `u`, `V`, and `w` themselves.

## Sign and frame conventions (PRD §0)

`+z` is the thrust direction. The projectile arrives travelling `-z`, so its momentum is a debit;
useful ejecta also leaves travelling `-z`. This is the origin of the `-1` in every impulse law.
Angles are measured from `+z`, so `mu = cos(theta) = 1` is the best-aimed, plate-directed element.

## What this module does *not* do

It is the ballistic model plus the ceiling. It has no pressure term, and the whole point of the
study is that the ballistic model omits pressure in the one regime where that is least defensible
(at `k ~ 7` the fireball's sound speed and its recoil velocity agree to within 5%, PRD §3.4). Every
number here is therefore a *reference* to score a hydrocode against, never an answer. The
assumptions carried in from prior work are enumerated and adjudicated in `assumption_audit()`
rather than left implicit.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Literal

# ---- Fixed physical and reference constants (PRD §0, §4) ----------------------------------------

G0 = 9.80665
"""Standard gravity [m/s^2] (PRD §0.1)."""

W_CLOSING_MS = 75_000.0
"""Closing speed `w` used by this study [m/s]; the envelope is 74-81 km/s (PRD §0.2)."""

W_BAR_PRIOR_MS = 77_280.0
"""Mass-weighted mean closing speed `w_bar`, used *only* to reproduce prior work's 1014 s figure."""

K_BARE_OPTIMUM = 7.06
"""The bare-plate ballistic Isp optimum `k*_bare` quoted throughout the PRD (§3.2, §3.5)."""

M_ENC_REF_KG = 200.0
"""Reference encounter mass per pulse `m_enc` = projectile + near-field carried mass (PRD §4)."""

M_VEHICLE_KG = 1.0e6
"""Vehicle *initial* (wet) mass, pinned at 1000 t (PRD §4.1). Never charged in Isp."""

DEPARTURE_DELTA_V_MS = 7060.0
"""Departure-burn `Delta_v` used in the §4.1 mission ledger [m/s]."""

M_H2O_KG = 0.018015 / 6.02214076e23
"""Mass of one water molecule [kg] — for the eV/molecule statement of blob internal energy."""

J_PER_EV = 1.602176634e-19

DEFAULT_OUTPUT_DIR = Path("data/results/tamper")

PlateShape = Literal["parabolic", "flat"]
"""`parabolic` = focus-matched paraboloid (`Delta_p = |v| + v_z`); `flat` = specular (`2 v_z`).

`beta_bare` is the *parabolic* case: the collimation prize is already inside the headline 984 s
rather than upside on top of it (PRD §3.6).
"""


# ---- The isotropic ballistic fireball -----------------------------------------------------------


def _u_hat(k: float) -> float:
    """Blob-frame expansion speed in units of `w`: `u/w = sqrt(k)/(1+k)`."""
    return math.sqrt(k) / (1.0 + k)


def _v_hat(k: float) -> float:
    """Blob centre-of-mass recoil speed in units of `w`: `V/w = 1/(1+k)`."""
    return 1.0 / (1.0 + k)


def _mean_speed_over_cone(k: float, mu_lo: float) -> float:
    """Mass-weighted mean vehicle-frame *speed* `|v|` over the cone `mu in [mu_lo, 1]`, in `w`.

    For an isotropic blob-frame emission the mass element is uniform in `mu`, and
    `|v|^2 = u^2 + V^2 - 2 u V mu`, so the mean has the closed form

        (2/(3b)) [ (a - b mu_lo)^{3/2} - (a - b)^{3/2} ] / (1 - mu_lo),
        a = u^2 + V^2,  b = 2 u V.

    Exact, not quadrature: this is one of the few places where an integral would otherwise invite
    a silent convergence error into a "closed-form" reference.
    """
    if mu_lo >= 1.0:
        return 0.0
    u, v = _u_hat(k), _v_hat(k)
    a = u * u + v * v
    b = 2.0 * u * v
    if b <= 1e-300:  # k = 0 (no expansion): every element moves at V.
        return math.sqrt(max(a, 0.0))
    lo = max(a - b * mu_lo, 0.0)
    hi = max(a - b, 0.0)
    lo_32, hi_32 = lo * math.sqrt(lo), hi * math.sqrt(hi)
    return (2.0 / (3.0 * b)) * (lo_32 - hi_32) / (1.0 - mu_lo)


def _mean_axial_over_cone(k: float, mu_lo: float) -> float:
    """Mass-weighted mean axial approach speed `v_z = u mu - V` over the cone, in units of `w`."""
    if mu_lo >= 1.0:
        return 0.0
    return _u_hat(k) * (1.0 + mu_lo) / 2.0 - _v_hat(k)


@dataclass(frozen=True)
class Fireball:
    """The merged blob's ballistic state at one slug ratio (PRD §2.2, §4)."""

    k: float
    w: float
    v_cm: float
    """Centre-of-mass recoil speed `V = w/(1+k)` [m/s], directed *away* from the plate."""
    u: float
    """Blob-frame expansion speed `u = w sqrt(k)/(1+k)` [m/s]."""
    specific_internal_energy: float
    """`e = (1/2) w^2 k/(1+k)^2` [J/kg] — the merge-dissipated energy per kg of blob."""

    @classmethod
    def from_k(cls, k: float, w: float = W_CLOSING_MS) -> Fireball:
        return cls(
            k=k,
            w=w,
            v_cm=w * _v_hat(k),
            u=w * _u_hat(k),
            specific_internal_energy=0.5 * w * w * k / (1.0 + k) ** 2,
        )

    @property
    def u_over_v(self) -> float:
        """`u/V = sqrt(k)` — the ratio that sets the ballistic capture fraction."""
        return self.u / self.v_cm

    @property
    def sound_speed(self) -> float:
        """`c_s = sqrt(gamma (gamma-1) e)` at `gamma_eff = 1.25` — ~9.8 km/s here (PRD §4)."""
        gamma = 1.25
        return math.sqrt(gamma * (gamma - 1.0) * self.specific_internal_energy)


def projectile_energy(m_i: float, w: float = W_CLOSING_MS) -> float:
    """Incoming projectile kinetic energy per pulse, `E = (1/2) m_i w^2` [J] (PRD §0.2)."""
    return 0.5 * m_i * w * w


def ev_per_water_molecule(specific_energy: float) -> float:
    """Specific energy [J/kg] expressed as eV per H2O molecule (PRD §4: 305.7 MJ/kg = 57.1 eV)."""
    return specific_energy * M_H2O_KG / J_PER_EV


# ---- Capture geometry, in its three conventions (PRD §3.6, §13.13) -------------------------------


def mu_capture_infinite(k: float) -> float:
    """`cos(theta)` threshold for outrunning the recoil: `1/sqrt(k)`, clamped at 1 for `k <= 1`."""
    if k <= 1.0:
        return 1.0
    return 1.0 / math.sqrt(k)


def ballistic_capture_fraction(k: float) -> float:
    """`max[0, (1 - 1/sqrt(k))/2]` — the `R -> infinity` convention (PRD §0.1).

    Share of the isotropic, pressure-free ballistic fireball that outruns its own recoil. **It
    imposes no plate radius.** Its zero below `k = 1` is a limit of the ballistic model, not a
    claim that hydrodynamic thrust is exactly zero there (PRD §2.2, §8).
    """
    return max(0.0, (1.0 - mu_capture_infinite(k)) / 2.0)


def _rim_cutoff_roots(k: float, cos_theta_max: float) -> tuple[float, float]:
    """The two `mu` roots where a ray leaves at exactly the rim half-angle `theta_max`.

    A ray is caught only if its own direction lies inside the cone the rim subtends *at the
    source*: `cos(theta_ray) = v_z/|v| > cos(theta_max)`. Squaring gives a quadratic in `mu`
    whose roots are

        mu_pm = [ V(1-c^2) +- c sqrt(u^2 - V^2 (1-c^2)) ] / u,   c = cos(theta_max),

    with `mu_+` the *forward* capture threshold and `mu_-` the threshold for material that is
    turned around first (§6.3a). Both are needed: the tamper's whole job lives below `mu_-`.
    """
    u, v = _u_hat(k), _v_hat(k)
    s2 = 1.0 - cos_theta_max * cos_theta_max
    disc = u * u - v * v * s2
    if disc <= 0.0:
        return (-1.0, 1.0)
    root = cos_theta_max * math.sqrt(disc)
    return ((v * s2 - root) / u, (v * s2 + root) / u)


def rim_cos_cutoff_flat(*, r_rim: float, standoff: float) -> float:
    """`cos(theta_max)` for a **flat** plate: its rim sits in the plane a distance `d` away."""
    return math.cos(math.atan2(r_rim, standoff))


def rim_cos_cutoff_paraboloid(*, r_rim: float, focal_length: float) -> float:
    """`cos(theta_max)` for a **focus-matched paraboloid** — and it is *not* the flat value.

    With the source at the focus the surface is `rho(theta) = 2F/(1+cos theta)`, whose cylindrical
    radius is `2F tan(theta/2)`. The rim at radius `R` therefore sits at

        tan(theta_max/2) = R/(2F),

    and **the rim stands off the vertex toward the source** by the dish depth `delta = R^2/(4F)`.
    So a dish of the same rim radius subtends a much wider cone than a flat plate at the same
    vertex standoff, and captures far more. Past `delta/D = 0.25` (`F = R/2`) the rim reaches the
    source plane, `theta_max` exceeds 90 degrees, and `cos(theta_max)` goes negative — the dish
    begins catching *away-going* material, which is the tamper's job done by geometry.
    """
    return math.cos(2.0 * math.atan2(r_rim, 2.0 * focal_length))


def paraboloid_depth(*, r_rim: float, focal_length: float) -> float:
    """Dish depth `delta = R^2/(4F)` [m] — the **wall height**, vertex to rim, along the axis."""
    return r_rim * r_rim / (4.0 * focal_length)


def paraboloid_depth_ratio(*, r_rim: float, focal_length: float) -> float:
    """The shape parameter `delta/D = R/(8F)` — dish depth over dish *diameter* (`D = 2R`).

    0 is a flat plate; 0.25 is a bowl a quarter as deep as it is wide; 0.5 is a hemisphere.
    """
    return r_rim / (8.0 * focal_length)


def paraboloid_area(*, r_rim: float, focal_length: float) -> float:
    """Surface area of the dish [m^2], `(8 pi F^2/3)[(1 + R^2/(4F^2))^{3/2} - 1]`.

    This is the ablator's coating area and the plate's mass basis, so it is the honest cost side
    of dish depth. It grows far more slowly than capture does.
    """
    x = 1.0 + r_rim**2 / (4.0 * focal_length**2)
    return (8.0 * math.pi * focal_length**2 / 3.0) * (x * math.sqrt(x) - 1.0)


def mu_capture_cutoff(k: float, cos_theta_max: float) -> float:
    """Forward-capture `mu` threshold for a rim half-angle `theta_max` (see `_rim_cutoff_roots`)."""
    if k <= 1.0 and cos_theta_max >= 0.0:
        return 1.0
    mu_c = mu_capture_infinite(k) if cos_theta_max >= 0.0 else -1.0
    if cos_theta_max <= -1.0:
        return -1.0
    _, mu_plus = _rim_cutoff_roots(k, cos_theta_max)
    return min(1.0, max(mu_c, mu_plus))


def mu_capture_finite(k: float, r_over_d: float) -> float:
    """`cos(theta)` threshold for a *finite flat* plate of radius `R` at standoff `d` (PRD §3.6).

    A ray is caught only if it both outruns the recoil and lands inside the rim,
    `v_r/v_z <= R/d`. Rays near the capture threshold arrive nearly grazing and land at large
    radius, so a finite plate loses them — which is why the geometry is far more sensitive to
    `R/d` than the `R -> infinity` figures suggest.

    **This is the flat-plate rim position.** A dish of the same rim radius has its rim standing
    off toward the source, so it must be scored with `rim_cos_cutoff_paraboloid` instead.
    """
    if k <= 1.0:
        return 1.0
    if math.isinf(r_over_d):
        return mu_capture_infinite(k)
    return mu_capture_cutoff(k, math.cos(math.atan(r_over_d)))


def capture_fraction_finite(k: float, r_over_d: float) -> float:
    """Ray-consistent capture fraction for a finite plate — 10.6% at `R/d = 1.5`, `k = 7.06`."""
    return max(0.0, (1.0 - mu_capture_finite(k, r_over_d)) / 2.0)


def capture_fraction_rim_angle(k: float, r_over_d: float) -> float:
    """The §5.3 `Sigma` convention: rim angle applied to the *blob-frame emission* angle.

    Geometrically inconsistent with `capture_fraction_finite` — it asks which elements are emitted
    into the rim's solid angle rather than which *rays* land inside the rim — and gives 22.3% at
    `R/d = 1.5` against the ray value's 10.6%. Retained because `Sigma`, `Phi`, and `tau_opt` in
    the PRD inherit it; see `capture_convention_bracket()` for the reconciliation (PRD §13.13).
    """
    mu_rim = 1.0 / math.sqrt(1.0 + r_over_d * r_over_d)
    return max(0.0, (1.0 - max(mu_capture_infinite(k), mu_rim)) / 2.0)


# ---- Impulse coefficients ------------------------------------------------------------------------


def _cone_impulse(k: float, mu_lo: float, mu_hi: float, plate: PlateShape) -> float:
    """Impulse `sum m (|v| + v_z)` (or `2 v_z` for a flat plate) over the band `[mu_lo, mu_hi]`."""
    if mu_hi <= mu_lo:
        return 0.0
    u, v = _u_hat(k), _v_hat(k)
    mass = (1.0 + k) * (mu_hi - mu_lo) / 2.0
    mean_axial = u * (mu_lo + mu_hi) / 2.0 - v
    if plate == "flat":
        return mass * 2.0 * mean_axial
    a, b = u * u + v * v, 2.0 * u * v
    lo, hi = max(a - b * mu_lo, 0.0), max(a - b * mu_hi, 0.0)
    mean_speed = (2.0 / (3.0 * b)) * (lo * math.sqrt(lo) - hi * math.sqrt(hi)) / (mu_hi - mu_lo)
    return mass * (mean_speed + mean_axial)


def beta_cutoff(k: float, cos_theta_max: float, plate: PlateShape = "parabolic") -> float:
    """`beta = J/(m_i w)` for a plate whose rim subtends half-angle `theta_max` at the source.

    `J = sum over captured mass of Delta_p`, because the blob's own momentum is exactly the
    incoming debit and cancels it (PRD §3.6): with no plate, `beta = 0`.
    """
    mu_lo = mu_capture_cutoff(k, cos_theta_max)
    return _cone_impulse(k, mu_lo, 1.0, plate) if mu_lo < 1.0 else 0.0


def beta_finite(k: float, r_over_d: float, plate: PlateShape = "parabolic") -> float:
    """`beta` for a **flat** plate of radius ratio `R/d`. Kept as the flat-plate entry point."""
    mu_lo = mu_capture_finite(k, r_over_d)
    return _cone_impulse(k, mu_lo, 1.0, plate) if mu_lo < 1.0 else 0.0


def beta_mirror_cutoff(k: float, cos_theta_max: float, plate: PlateShape = "parabolic") -> float:
    """A perfect-mirror tamper in front of a plate of *finite* rim angle — and its new failure mode.

    At `R -> infinity` a tamper can only help: everything it turns around is eventually caught. A
    real plate breaks that. Splitting the fireball by the sign of `v_z`:

    * `v_z > 0` — already plate-bound; caught iff `v_z/|v| > c`, exactly as in the bare case;
    * `v_z < 0` — the tamper reverses it, so it is caught iff `-v_z/|v| > c`, i.e. `mu < mu_-`;
    * `mu_- < mu < mu_c` — **turned around and then missed.** This mass was moving `-z`, which is
      the *useful* direction, and the tamper has flipped it to `+z` where it flies past the rim.
      Its contribution changes sign, costing `2 m |v_z|`.

    The plate's shape decides whether the *caught* half of that is worth anything. A collimating
    dish sends it to `-z` at full speed, converting it into impulse. A **flat** plate merely
    reverses the axial component, handing the gas back the velocity the tamper had just taken from
    it — the pair cancels exactly, and the tamper's only remaining effect is the debit. That is why
    a tamper in front of a flat plate produces *negative* net impulse.

    So `beta = sum_captured m(|v| + v_z) + 2 sum_turned-and-missed m v_z`, the second term being
    negative. **A tamper is only worth carrying if the plate can catch what it turns around** —
    which is why dish geometry is prerequisite to the tamper question, not parallel to it.
    """
    mu_c = mu_capture_infinite(k)
    if cos_theta_max < 0.0:
        # The rim already reaches past the source, so nothing the tamper turns can miss: every
        # element ends up moving -z at its own speed, which is the tamper's ideal. Note this is
        # the *same* value as at the knee — a tamper's benefit saturates at delta/D = 0.25 and a
        # deeper dish adds nothing to the tamped case (it only helps the bare one).
        mu_plus = mu_minus = mu_c
    else:
        mu_minus, mu_plus = _rim_cutoff_roots(k, cos_theta_max)
        mu_plus = min(1.0, max(mu_c, mu_plus))
        mu_minus = max(-1.0, min(mu_c, mu_minus))
    caught = _cone_impulse(k, mu_plus, 1.0, plate)
    if plate == "flat":
        # A flat plate reverses only the axial component, so gas the tamper turned around and
        # the plate then turns back leaves with its *original* momentum: the pair cancels and
        # contributes exactly nothing. Only the collimating dish converts it into impulse.
        pass
    else:
        caught += _cone_impulse(k, -1.0, mu_minus, plate)
    # Turned around, then missed: a credit becomes a debit.
    turned_missed_mass = (1.0 + k) * (mu_c - mu_minus) / 2.0
    mean_axial_missed = _u_hat(k) * (mu_minus + mu_c) / 2.0 - _v_hat(k)
    return caught + 2.0 * turned_missed_mass * mean_axial_missed


def beta_bare(k: float) -> float:
    """The bare-plate ballistic coefficient: `R -> infinity`, focus-matched paraboloid.

    `beta_bare(7.06) = 0.9087`. Two idealisations are already inside it (PRD §3.6): it collimates
    each captured element to `-z` at its *full* arrival speed, and it captures everything with
    `v_z > 0` at any radius. Both are named in `assumption_audit()`.
    """
    return beta_finite(k, math.inf, "parabolic")


def beta_flat(k: float) -> float:
    """Flat-plate specular counterpart, `(sqrt(k)-1)^2 / (2 sqrt(k))` (PRD §0.2, §3.6).

    Closed form; the integral in `beta_finite(k, inf, "flat")` reduces to exactly this.
    """
    if k <= 1.0:
        return 0.0
    root = math.sqrt(k)
    return (root - 1.0) ** 2 / (2.0 * root)


def beta_mirror(k: float) -> float:
    """Perfect-mirror tamper: an infinitely massive wall reversing each element at unchanged speed.

    Then *every* ejecta element ends up moving `-z` at its ballistic speed `|v|`, so
    `beta = (1+k) <|v|>/w - 1` over the full sphere. Gives 1.7825 at `k = 7.06` — 96.9% of the
    `K_ej = k` ceiling, and still short of the 1.817 needed to beat the bare plate at equal
    charged mass (PRD §3.4).

    **This is an upper bound, not a model of a tamper.** "Reflects all energy, absorbs none, still
    recoils" is self-contradictory at finite mass, which is why ADR-0030 frames the real tamper as
    an isentropic piston and credits its recoil instead.
    """
    return (1.0 + k) * _mean_speed_over_cone(k, -1.0) - 1.0


def beta_ideal(k_ej: float) -> float:
    """The thermodynamic ceiling `sqrt(1+K_ej) - 1` (PRD §3.2).

    From Cauchy-Schwarz on `|P_ejecta| <= sqrt(2 M_ej E)`: equality needs every ejecta element
    moving `-z` at one common speed. It depends only on `K_ej`, never on how that mass splits
    between slug, interlayer, tamper, and ablator — which is why **the tamper can never be
    justified as a momentum multiplier**, only as a realizability device (PRD §3.3a).
    """
    return math.sqrt(1.0 + k_ej) - 1.0


def j_max(k_ej: float, w: float = W_CLOSING_MS) -> float:
    """Ceiling impulse per projectile mass, `j_max = w (sqrt(1+K_ej) - 1)` [N.s/kg]."""
    return w * beta_ideal(k_ej)


def v_e_max(k_ej: float, w: float = W_CLOSING_MS) -> float:
    """Ceiling effective exhaust velocity `j_max/K_ej` [m/s]. Falls monotonically with `K_ej`."""
    if k_ej <= 0.0:
        return 0.0
    return j_max(k_ej, w) / k_ej


# ---- The deliverables and their comparison metrics (PRD §3.1, §3.3, §3.5) ------------------------


def isp_eff(*, beta: float, c_ratio: float, w: float = W_CLOSING_MS) -> float:
    """Effective specific impulse `beta w/(g0 C)` [s] — the study's deliverable (PRD §3.1).

    `C` is the *charged*-mass ratio: all expended carried mass per projectile kg. In Pass 1 the
    ablator is excluded, so every Pass-1 figure is an **upper bound** and must be labelled one.
    """
    if c_ratio <= 0.0:
        return 0.0
    return beta * w / (G0 * c_ratio)


def isp_bare(k: float, w: float = W_CLOSING_MS) -> float:
    """Bare-plate Pass-1 Isp at slug ratio `k` (no tamper, no interlayer, no ablator)."""
    return isp_eff(beta=beta_bare(k), c_ratio=k, w=w)


def realization_fraction(beta: float, k_ej: float) -> float:
    """`r_real = beta/(sqrt(1+K_ej) - 1)` — share of the same-ejecta-mass ceiling (PRD §3.3).

    A hydrodynamic comparison metric feeding Isp, not a replacement for it, and not the superseded
    tamper multiplier `Lambda` (ADR-0030).
    """
    ceiling = beta_ideal(k_ej)
    if ceiling <= 0.0:
        return 0.0
    return beta / ceiling


def projectile_economy(beta: float, w: float = W_CLOSING_MS) -> float:
    """`beta w = J/m_projectile` [N.s/kg] — the second reported metric (PRD §3.1).

    A velocity-equivalent impulse normalisation, *not* a physical velocity, and deliberately never
    combined with Isp: converting projectiles to payload-equivalent needs program economics this
    study does not model (PRD §3.1, D2b).
    """
    return beta * w


def energy_per_impulse(beta: float, w: float = W_CLOSING_MS) -> float:
    """`E/J = w/(2 beta)` [J per N.s] — plate heat load per unit impulse (PRD §3.5).

    Contains neither encounter mass nor cadence, so it is invariant under that trade (PRD §4.1).
    """
    if beta <= 0.0:
        return math.inf
    return w / (2.0 * beta)


def projectiles_per_impulse(beta: float, w: float = W_CLOSING_MS) -> float:
    """`1/(beta w)` [kg of projectile per N.s]. The same function of `beta` as heat load."""
    if beta <= 0.0:
        return math.inf
    return 1.0 / (beta * w)


def break_even_beta(*, beta_ref: float, c_ref: float, c_candidate: float) -> float:
    """`beta` a candidate needs to match the reference's Isp at its own charged mass (PRD §3.4).

    From `beta_candidate/C_candidate > beta_reference/C_reference`. Configuration-specific: 62.9%
    of ceiling is *one evaluation* of this rule, not a universal gate.
    """
    if c_ref <= 0.0:
        return math.inf
    return beta_ref * c_candidate / c_ref


def required_realization(
    *, beta_ref: float, c_ref: float, c_candidate: float, k_ej_candidate: float
) -> float:
    """The candidate's required realization fraction — the §3.4 rule in `r_real` form."""
    return realization_fraction(
        break_even_beta(beta_ref=beta_ref, c_ref=c_ref, c_candidate=c_candidate), k_ej_candidate
    )


# ---- Free-plate elastic lower bound, and the superseded `Lambda` framing ------------------------


def free_plate_areal_ratio(*, tau_t: float, k: float) -> float:
    """Tamper mass over away-going plume mass, per projectile kg: `tau_t k / ((1+k)/2)` (§6.3)."""
    plume = (1.0 + k) / 2.0
    if plume <= 0.0:
        return math.inf
    return tau_t * k / plume


def free_plate_reflected_fraction(*, tau_t: float, k: float) -> float:
    """1-D elastic free-plate reflected/incident speed ratio, `|r-1|/(r+1)` (PRD §6.3, §8).

    The tamper is *not* immovable: it masses `tau_t k` against an away-going plume of `(1+k)/2`,
    and a 1-D elastic encounter between two free bodies returns only 27% of the incident speed at
    `tau_t = 1` — the honest lower bound against the mirror's 100%. Magnitude is taken because
    below `r = 1` the lighter tamper is pushed along rather than reflecting anything back.
    """
    r = free_plate_areal_ratio(tau_t=tau_t, k=k)
    if math.isinf(r):
        return 1.0
    return abs(r - 1.0) / (r + 1.0)


def legacy_lambda(*, beta_tamped: float, beta_reference: float) -> float:
    """The **superseded** tamper multiplier `Lambda = impulse(tamped)/impulse(bare)` (ADR-0030).

    Exposed only as provenance for the 591-965 s bracket §1 quotes from prior work. Do not use it
    as a metric: it is a ratio of two realized numbers, so it hides the ceiling and miscounts the
    tamper's credited recoil as loss. `realization_fraction` replaces it (CONTEXT.md).
    """
    if beta_reference <= 0.0:
        return math.inf
    return beta_tamped / beta_reference


# ---- Mass ledger (PRD §0.1) ----------------------------------------------------------------------


@dataclass(frozen=True)
class MassLedger:
    """Per-pulse mass accounting at a stated encounter mass, with §0.1's conventions applied.

    The three mass ratios differ only when an interlayer or ablator is present, which is exactly
    when confusing them changes an answer: `K` inventories near-field hydrodynamic mass, `K_ej`
    the mass in the ceiling, and `C` the carried mass Isp charges.
    """

    k: float
    tau_t: float
    mu: float
    a_abl: float
    a_other: float
    m_enc: float

    @classmethod
    def from_ratios(
        cls,
        *,
        k: float,
        tau_t: float = 0.0,
        mu: float = 0.0,
        a_abl: float = 0.0,
        a_other: float = 0.0,
        m_enc: float = M_ENC_REF_KG,
    ) -> MassLedger:
        return cls(k=k, tau_t=tau_t, mu=mu, a_abl=a_abl, a_other=a_other, m_enc=m_enc)

    @property
    def k_hydro(self) -> float:
        """`K = k(1 + tau_t + mu)` — near-field hydrodynamic carried mass per projectile kg."""
        return self.k * (1.0 + self.tau_t + self.mu)

    @property
    def c_charged(self) -> float:
        """`C = K + a_abl + a_other` — the Isp denominator per projectile kg."""
        return self.k_hydro + self.a_abl + self.a_other

    @property
    def k_ej(self) -> float:
        """`K_ej` — nonprojectile ejecta mass per projectile kg.

        Pass 1 is a closed no-ablator calculation, so `K_ej = K`. Pass 2 must take the ablator
        mass *actually ejected* from the wall model rather than assuming it equals `a_abl`.
        """
        return self.k_hydro

    @property
    def m_i(self) -> float:
        """Projectile mass per pulse [kg]: `m_enc/(1+K)`. Never charged in Isp."""
        return self.m_enc / (1.0 + self.k_hydro)

    @property
    def m_s(self) -> float:
        return self.k * self.m_i

    @property
    def m_t(self) -> float:
        return self.tau_t * self.m_s

    @property
    def m_int(self) -> float:
        return self.mu * self.m_s

    @property
    def m_hydro(self) -> float:
        return self.k_hydro * self.m_i

    @property
    def m_abl(self) -> float:
        return self.a_abl * self.m_i

    @property
    def m_charged(self) -> float:
        """All expended carried mass per pulse [kg] — the Isp denominator."""
        return self.c_charged * self.m_i


# ---- Mission ledger (PRD §4.1) -------------------------------------------------------------------


@dataclass(frozen=True)
class VehicleContext:
    """Per-pulse vehicle-scale figures at one configuration (PRD §4.1)."""

    ledger: MassLedger
    beta: float
    w: float

    @property
    def impulse_per_pulse(self) -> float:
        """`J = beta m_i w` [N.s]."""
        return self.beta * self.ledger.m_i * self.w

    @property
    def delta_v_per_pulse(self) -> float:
        """`Delta_v_vehicle = J/M_vehicle` [m/s] on the 1000 t *initial* mass."""
        return self.impulse_per_pulse / M_VEHICLE_KG

    def acceleration_g(self, cadence_hz: float) -> float:
        """Vehicle acceleration at a stated cadence, in `g`."""
        return self.delta_v_per_pulse * cadence_hz / G0

    def carried_mass_flow(self, cadence_hz: float) -> float:
        """Charged carried-mass flow [kg/s]. Excludes the externally supplied projectile."""
        return self.ledger.m_charged * cadence_hz

    def encounter_mass_flow(self, cadence_hz: float) -> float:
        """Total encounter mass flow [kg/s], *including* the projectile."""
        return self.ledger.m_enc * cadence_hz

    @property
    def isp_eff(self) -> float:
        return isp_eff(beta=self.beta, c_ratio=self.ledger.c_charged, w=self.w)

    @property
    def m_charged(self) -> float:
        return self.ledger.m_charged


def vehicle_context(
    *,
    k: float,
    tau_t: float = 0.0,
    mu: float = 0.0,
    m_enc: float = M_ENC_REF_KG,
    w: float = W_CLOSING_MS,
    beta: float | None = None,
) -> VehicleContext:
    """§4.1's per-pulse ledger. `beta` defaults to the bare ballistic value at this `k`."""
    led = MassLedger.from_ratios(k=k, tau_t=tau_t, mu=mu, m_enc=m_enc)
    return VehicleContext(ledger=led, beta=beta_bare(k) if beta is None else beta, w=w)


@dataclass(frozen=True)
class DepartureBurn:
    """The §4.1 departure-burn row: mass ratio, charged mass, pulse count, duration."""

    isp_s: float
    delta_v: float
    m_initial: float
    carried_per_pulse_kg: float

    @property
    def mass_ratio(self) -> float:
        return math.exp(self.delta_v / (self.isp_s * G0))

    @property
    def charged_mass_t(self) -> float:
        return self.m_initial * (1.0 - 1.0 / self.mass_ratio) / 1000.0

    @property
    def pulses(self) -> float:
        return self.charged_mass_t * 1000.0 / self.carried_per_pulse_kg

    def duration_s(self, *, cadence_hz: float) -> float:
        return self.pulses / cadence_hz


def departure_burn(
    *,
    isp_s: float,
    carried_per_pulse_kg: float,
    delta_v: float = DEPARTURE_DELTA_V_MS,
    m_initial: float = M_VEHICLE_KG,
) -> DepartureBurn:
    return DepartureBurn(
        isp_s=isp_s,
        delta_v=delta_v,
        m_initial=m_initial,
        carried_per_pulse_kg=carried_per_pulse_kg,
    )


# ---- Plate soak chain (PRD §6.5) -----------------------------------------------------------------
#
# Rung 0 owns this chain because it had already drifted once: §6.5 quoted a 173 um soak depth,
# which implies alpha_th = 1.0e-5 m^2/s against the 1.2e-5 stated in §0.6. The depth propagates
# into an areal mass, a heat capacity, the 4-row regenerative table, and the regenerative cap, so
# patching it row by row is how the next inconsistency gets introduced. One function, one source.

STEEL_K_COND = 45.0
"""Steel thermal conductivity [W/m/K] (PRD §0.6)."""

STEEL_C_P = 500.0
"""Steel specific heat [J/kg/K] (PRD §0.6)."""

STEEL_RHO = 7850.0
"""Steel density [kg/m^3] — implied by §6.5's own depth-to-areal-mass step."""

STEEL_ALPHA_TH = 1.2e-5
"""Steel thermal diffusivity [m^2/s] (PRD §0.6) — the single stated value this chain uses.

`k_cond/(rho c_p)` with the quoted constants gives 1.15e-5, a 4% difference that moves no
conclusion here; the stated `alpha_th` governs so the three constants cannot drift apart again.
"""

DELTA_T_TO_MELT_K = 1400.0
"""Surface-to-melt temperature interval used in the soak *capacity* (PRD §6.5)."""

RESIDENCE_TIME_S = 750.0e-6
"""Arrival window at the plate — the sustained-feed residence time (PRD §5.2, §6.4)."""

INCIDENT_FLUENCE_J_M2 = 12.6e6
"""Plate-area-averaged incident fluence at `R` = 15 m, `d` = 10 m (PRD §0.6). Order of magnitude."""

PLATE_RADIUS_REF_M = 15.0


@dataclass(frozen=True)
class SoakChain:
    """The §6.5 soak *capacity* chain, from one thermal diffusivity.

    Every row is a capacity, not a prediction: it assumes the full thermal penetration depth
    reaches melting, which happens only if the vapour curtain fails completely. §6.5.2's
    self-limiting argument puts the physical soak two orders below it.
    """

    alpha_th: float
    residence_s: float
    depth_m: float
    """`sqrt(4 alpha_th t)` — the thermal penetration depth over one pulse."""
    areal_mass_kg_m2: float
    capacity_j_m2: float
    share_of_fluence: float
    plate_area_m2: float
    basis_j: float
    """Whole-plate soak capacity — the basis the §6.5.3 regenerative shares are taken against."""


def plate_soak_chain(
    *,
    alpha_th: float = STEEL_ALPHA_TH,
    residence_s: float = RESIDENCE_TIME_S,
    plate_radius_m: float = PLATE_RADIUS_REF_M,
    fluence_j_m2: float = INCIDENT_FLUENCE_J_M2,
) -> SoakChain:
    """Soak depth -> areal mass -> capacity -> share of fluence -> whole-plate basis (PRD §6.5)."""
    depth = math.sqrt(4.0 * alpha_th * residence_s)
    areal_mass = depth * STEEL_RHO
    capacity = areal_mass * STEEL_C_P * DELTA_T_TO_MELT_K
    area = math.pi * plate_radius_m**2
    return SoakChain(
        alpha_th=alpha_th,
        residence_s=residence_s,
        depth_m=depth,
        areal_mass_kg_m2=areal_mass,
        capacity_j_m2=capacity,
        share_of_fluence=capacity / fluence_j_m2,
        plate_area_m2=area,
        basis_j=capacity * area,
    )


def conducted_per_pulse_j_m2(
    *,
    delta_t: float,
    residence_s: float = RESIDENCE_TIME_S,
    alpha_th: float = STEEL_ALPHA_TH,
    k_cond: float = STEEL_K_COND,
) -> float:
    """Semi-infinite conducted-in energy per pulse, `k_cond dT sqrt(t/(pi alpha_th))` (PRD §6.5.2).

    This is the *physical* soak while an unbroken ablating layer pins the surface at `T_abl`, as
    opposed to the penetration-depth capacity above. `delta_t` is ablator surface to plate bulk.
    """
    return k_cond * delta_t * math.sqrt(residence_s / (math.pi * alpha_th))


@dataclass(frozen=True)
class RegenerativeRow:
    """One row of §6.5.3's free-regenerative-cooling budget."""

    water_state: str
    delta_h_j_kg: float
    per_pulse_j: float
    share_of_basis: float
    admissible_as_slug: bool
    """Boiled coolant cannot double as propellant: a steam slug is ~300x too dilute (§6.2)."""


def regenerative_budget(
    *, basis_j: float | None = None, coolant_per_pulse_kg: float = M_ENC_REF_KG
) -> list[RegenerativeRow]:
    """§6.5.3: water is already the propellant, so this coolant flows regardless — it is mass-free.

    The cap is set by *phase*, not by enthalpy: melting is admissible because liquid water at
    1000 kg/m^3 still has ample column density to stop a projectile, but anything boiled must
    reject its latent heat again before it can become a slug.
    """
    basis = plate_soak_chain().basis_j if basis_j is None else basis_j
    states = [
        ("ice warmed to just under melt", 0.33e6, True),
        ("melted to liquid at 273 K", 0.66e6, True),
        ("liquid at 373 K", 1.08e6, True),
        ("saturated steam", 3.34e6, False),
    ]
    return [
        RegenerativeRow(
            water_state=name,
            delta_h_j_kg=delta_h,
            per_pulse_j=delta_h * coolant_per_pulse_kg,
            share_of_basis=delta_h * coolant_per_pulse_kg / basis,
            admissible_as_slug=admissible,
        )
        for name, delta_h, admissible in states
    ]


def regenerative_cap() -> float:
    """The free regenerative budget's cap: the best share reachable without boiling (§6.5.3)."""
    return max(row.share_of_basis for row in regenerative_budget() if row.admissible_as_slug)


# ---- Optimisation over the mass ratio (PRD §3.5) -------------------------------------------------


def _argmax_golden(f: Callable[[float], float], lo: float, hi: float, tol: float = 1e-9) -> float:
    """Golden-section maximisation of a unimodal `f` on `[lo, hi]`. No SciPy dependency."""
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - inv_phi * (b - a), a + inv_phi * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = f(d)
    return (a + b) / 2.0


def optimal_k_bare(lo: float = 2.0, hi: float = 30.0, w: float = W_CLOSING_MS) -> float:
    """`k*_bare` — the bare-plate Pass-1 Isp optimum, 7.060 (PRD §3.2).

    Prior work's quoted 7.057 is the same optimum to optimiser tolerance: the curve is flat to
    ±0.6% over `k` = 6-8, which is *why* `K` must be swept 6-32 rather than clustered here.
    """
    return _argmax_golden(lambda k: isp_bare(k, w=w), lo, hi)


# ---- Assembled reference tables ------------------------------------------------------------------


@dataclass(frozen=True)
class ComparisonRow:
    """One row of the §3.4 reference comparison."""

    label: str
    k_ratio: float
    """`K` = `K_ej` = `C` for these no-ablator reference rows."""
    beta: float
    realization: float
    isp_s: float


def _row(label: str, k_ratio: float, beta: float, w: float = W_CLOSING_MS) -> ComparisonRow:
    return ComparisonRow(
        label=label,
        k_ratio=k_ratio,
        beta=beta,
        realization=realization_fraction(beta, k_ratio),
        isp_s=isp_eff(beta=beta, c_ratio=k_ratio, w=w),
    )


def reference_comparison(
    *, k: float = K_BARE_OPTIMUM, tau_t: float = 1.0, w: float = W_CLOSING_MS
) -> list[ComparisonRow]:
    """PRD §3.4's reference comparison, regenerated (`mu = a_abl = 0`).

    Read row by row: at `tau_t = 1` the tamper beats spending the same mass as extra slug, but
    loses to not spending it at all — *and that is with a perfect mirror*. So the reference tamper
    pays only if pressure coupling, omitted from every model here, closes the 1.2-point gap.
    """
    k_tamped = k * (1.0 + tau_t)
    beta_ref = beta_bare(k)
    return [
        _row(f"bare plate, k = {k:g}", k, beta_ref, w),
        _row(f"bare plate, k = {k_tamped:g}", k_tamped, beta_bare(k_tamped), w),
        _row(f"perfect-mirror tamper, tau_t = {tau_t:g}", k_tamped, beta_mirror(k), w),
        _row(
            f"break-even against the bare plate at k = {k:g}",
            k_tamped,
            break_even_beta(beta_ref=beta_ref, c_ref=k, c_candidate=k_tamped),
            w,
        ),
        _row(f"ceiling at K = {k:g}", k, beta_ideal(k), w),
        _row(f"ceiling at K = {k_tamped:g}", k_tamped, beta_ideal(k_tamped), w),
    ]


@dataclass(frozen=True)
class MassRatioRow:
    """One row of the §3.5 mass-ratio trade."""

    k: float
    beta: float
    isp_s: float
    isp_vs_peak: float
    energy_per_impulse_j: float
    energy_vs_peak: float
    realization: float
    projectile_economy_ms: float


def mass_ratio_table(
    ks: Sequence[float] = (6.0, 7.06, 8.0, 10.0, 14.12, 16.0, 32.0), w: float = W_CLOSING_MS
) -> list[MassRatioRow]:
    """PRD §3.5: Isp degrades gently with `K` while energy per unit impulse improves steeply.

    Plate heat load and projectile consumption are *the same function* (both `1/beta`), so they
    are one quantity rather than two currencies needing an exchange rate (PRD §3.5, D2b).
    """
    k_star = optimal_k_bare(w=w)
    peak_isp = isp_bare(k_star, w=w)
    peak_energy = energy_per_impulse(beta_bare(k_star), w=w)
    rows: list[MassRatioRow] = []
    for k in ks:
        beta = beta_bare(k)
        rows.append(
            MassRatioRow(
                k=k,
                beta=beta,
                isp_s=isp_bare(k, w=w),
                isp_vs_peak=isp_bare(k, w=w) / peak_isp - 1.0,
                energy_per_impulse_j=energy_per_impulse(beta, w=w),
                energy_vs_peak=energy_per_impulse(beta, w=w) / peak_energy - 1.0,
                realization=realization_fraction(beta, k),
                projectile_economy_ms=projectile_economy(beta, w=w),
            )
        )
    return rows


@dataclass(frozen=True)
class WorkedExample:
    """PRD §3.6 line by line — where the ballistic numbers come from, and what they assume."""

    k: float
    w: float
    v_cm: float
    u: float
    u_over_v: float
    mu_capture: float
    capture_fraction: float
    captured_mass_per_projectile_kg: float
    blob_momentum_per_projectile_kg: float
    best_aimed_axial: float
    mean_axial_approach: float
    mean_speed_leaving: float
    delta_p_per_captured_kg: float
    beta: float
    v_e: float
    isp_s: float
    beta_flat: float
    isp_flat_s: float
    coherent_speed: float
    v_e_max: float
    isp_max_s: float
    realization: float


def worked_example(*, k: float = 6.0, w: float = W_CLOSING_MS) -> WorkedExample:
    """Regenerate every line of §3.6's ledger at one `k`.

    Two of its steps are not obvious. **The mean axial approach is half the best-aimed element**,
    because the cone's rim barely outruns the recoil and its axial speed goes to zero there. And
    **`beta_bare` already assumes a perfectly collimating plate** — the concave-plate prize is
    inside the 979 s, not upside on top of it.
    """
    fb = Fireball.from_k(k, w=w)
    mu_c = mu_capture_infinite(k)
    capture = ballistic_capture_fraction(k)
    beta = beta_bare(k)
    captured_mass = (1.0 + k) * capture
    mean_axial = _mean_axial_over_cone(k, mu_c) * w
    mean_speed = _mean_speed_over_cone(k, mu_c) * w
    return WorkedExample(
        k=k,
        w=w,
        v_cm=fb.v_cm,
        u=fb.u,
        u_over_v=fb.u_over_v,
        mu_capture=mu_c,
        capture_fraction=capture,
        captured_mass_per_projectile_kg=captured_mass,
        blob_momentum_per_projectile_kg=(1.0 + k) * fb.v_cm,
        best_aimed_axial=fb.u - fb.v_cm,
        mean_axial_approach=mean_axial,
        mean_speed_leaving=mean_speed,
        delta_p_per_captured_kg=mean_axial + mean_speed,
        beta=beta,
        v_e=beta * w / k,
        isp_s=isp_bare(k, w=w),
        beta_flat=beta_flat(k),
        isp_flat_s=isp_eff(beta=beta_flat(k), c_ratio=k, w=w),
        coherent_speed=w / math.sqrt(1.0 + k),
        v_e_max=v_e_max(k, w=w),
        isp_max_s=v_e_max(k, w=w) / G0,
        realization=realization_fraction(beta, k),
    )


@dataclass(frozen=True)
class CaptureRow:
    """One `R/d` row of §3.6's finite-plate table."""

    r_over_d: float
    capture_ray: float
    capture_rim_angle: float
    beta_parabolic: float
    isp_parabolic_s: float
    beta_flat: float
    isp_flat_s: float
    parabola_over_flat: float


def capture_convention_bracket(
    *, k: float = K_BARE_OPTIMUM, r_over_ds: Sequence[float] = (1.5, 2.5), w: float = W_CLOSING_MS
) -> list[CaptureRow]:
    """PRD §3.6 / §13.13 — the finite-plate table and the three capture conventions.

    Three mutually inconsistent capture fractions are in circulation at `k = 7.06`: **31.2%**
    (`R -> infinity`, which every headline Isp uses and which imposes no plate radius), **22.3%**
    (§5.3's `Sigma`, applying the rim angle to the blob-frame emission angle), and **10.6%** (the
    geometrically consistent finite-plate ray value at the reference `R/d` = 1.5).

    **Recommended convention** (Rung 0's reconciliation): quote the ray-consistent finite-plate
    value as the lower edge and the `R -> infinity` value as the upper edge of an explicit
    bracket, and retire the rim-angle form — it is not a capture fraction of anything, since it
    asks which elements are *emitted* into the rim's solid angle rather than which rays land
    inside the rim. Pure ray-tracing *understates* capture, because the flow is pressure-bearing
    at Mach ~2.5 and steers inward — this study's central thesis — so the truth lies between the
    two edges and only a simulation places it (Rung 4).
    """
    rows: list[CaptureRow] = []
    for q in (*r_over_ds, math.inf):
        par = beta_finite(k, q, "parabolic")
        flat = beta_finite(k, q, "flat")
        rows.append(
            CaptureRow(
                r_over_d=q,
                capture_ray=capture_fraction_finite(k, q),
                capture_rim_angle=(
                    ballistic_capture_fraction(k)
                    if math.isinf(q)
                    else capture_fraction_rim_angle(k, q)
                ),
                beta_parabolic=par,
                isp_parabolic_s=isp_eff(beta=par, c_ratio=k, w=w),
                beta_flat=flat,
                isp_flat_s=isp_eff(beta=flat, c_ratio=k, w=w),
                parabola_over_flat=par / flat if flat > 0.0 else math.inf,
            )
        )
    return rows


# ---- Inherited-assumption audit ------------------------------------------------------------------


@dataclass(frozen=True)
class DishRow:
    """One dish-depth row of §6.6: geometry, cost, and what it captures."""

    label: str
    depth_ratio: float
    """`delta/D` — dish depth over dish diameter. 0 is flat, 0.25 the knee, 0.5 a hemisphere."""
    focal_length_m: float
    wall_height_m: float
    """`delta` — the physical depth of the bowl, vertex to rim. The buildability constraint."""
    rim_vs_source_m: float
    """Axial position of the rim relative to the source; negative means it wraps behind it."""
    theta_max_deg: float
    area_m2: float
    plate_mass_t: float
    capture: float
    beta_bare: float
    isp_bare_s: float
    beta_tamped: float
    """Perfect-mirror tamper at `tau_t = 1` in front of *this* rim (`beta_mirror_cutoff`)."""
    isp_tamped_s: float
    turned_and_missed: float
    """Ejecta fraction the tamper turns around into a rim that cannot catch it — a pure debit."""


PLATE_AREAL_DENSITY_KG_M2 = 25.0
"""Plate areal density [kg/m^2] implied by §6.5.4's "50 t plate, 3.2 mm steel at R = 25 m"."""

PLATE_MASS_CEILING_T = 50.0
"""Plate mass ceiling: 5% of the 1000 t vehicle (PRD §4.1)."""


def dish_table(
    *,
    k: float = K_BARE_OPTIMUM,
    r_rim: float = PLATE_RADIUS_REF_M,
    standoff_flat: float = 10.0,
    focal_lengths: Sequence[float] = (12.0, 10.0, 9.0, 8.0, 7.5, 6.5, 5.4),
    w: float = W_CLOSING_MS,
) -> list[DishRow]:
    """PRD §6.6 — capture, impulse, and cost against dish depth, with the flat plate as row 0."""

    def row(
        label: str,
        cos_cut: float,
        focal: float,
        depth: float,
        area: float,
        plate: PlateShape = "parabolic",
    ) -> DishRow:
        # Each plate is scored with the momentum transfer it actually performs: a flat plate
        # reverses only the axial component, a focus-matched dish returns the full speed along
        # the axis. Pairing one plate's rim position with the other's momentum transfer is the
        # error this table exists to correct (§6.6.2).
        beta_b = beta_cutoff(k, cos_cut, plate)
        beta_t = beta_mirror_cutoff(k, cos_cut, plate)
        mu_minus, _ = _rim_cutoff_roots(k, cos_cut)
        mu_c = mu_capture_infinite(k)
        missed = max(0.0, (mu_c - max(-1.0, min(mu_c, mu_minus))) / 2.0) if cos_cut >= 0.0 else 0.0
        return DishRow(
            label=label,
            depth_ratio=depth / (4.0 * r_rim),
            focal_length_m=focal,
            wall_height_m=depth,
            rim_vs_source_m=focal - depth,
            theta_max_deg=math.degrees(math.acos(max(-1.0, min(1.0, cos_cut)))),
            area_m2=area,
            plate_mass_t=area * PLATE_AREAL_DENSITY_KG_M2 / 1000.0,
            capture=max(0.0, (1.0 - mu_capture_cutoff(k, cos_cut)) / 2.0),
            beta_bare=beta_b,
            isp_bare_s=isp_eff(beta=beta_b, c_ratio=k, w=w),
            beta_tamped=beta_t,
            isp_tamped_s=isp_eff(beta=beta_t, c_ratio=2.0 * k, w=w),
            turned_and_missed=missed,
        )

    rows = [
        row(
            "flat plate",
            rim_cos_cutoff_flat(r_rim=r_rim, standoff=standoff_flat),
            standoff_flat,
            0.0,
            math.pi * r_rim**2,
            "flat",
        )
    ]
    for focal in focal_lengths:
        depth = paraboloid_depth(r_rim=r_rim, focal_length=focal)
        rows.append(
            row(
                f"dish, delta/D = {paraboloid_depth_ratio(r_rim=r_rim, focal_length=focal):.3f}",
                rim_cos_cutoff_paraboloid(r_rim=r_rim, focal_length=focal),
                focal,
                depth,
                paraboloid_area(r_rim=r_rim, focal_length=focal),
            )
        )
    return rows


@dataclass(frozen=True)
class Assumption:
    """One inherited modelling assumption, adjudicated for *this* study at *this* rung."""

    ident: str
    statement: str
    source: str
    applies_to: tuple[str, ...]
    verdict: Literal["holds", "holds-with-caveat", "departs"]
    why: str
    retired_by: str
    """Which rung or mechanism removes the assumption — never "not yet" for a departure."""


def assumption_audit() -> list[Assumption]:
    """Every assumption the Rung 0 closed forms inherit, and whether it is still appropriate.

    A standing instruction on this project: assumptions are re-checked at each step against the
    step that uses them, not carried forward silently. At Rung 0 there is no solver, so what is
    audited is the *ballistic model* and the *ceiling* — and the honest summary is that the
    ballistic model's central assumption (pressure-free flight) is the one the whole study exists
    to remove, while the ceiling's assumptions are thermodynamic identities that hold everywhere.
    """
    return [
        Assumption(
            ident="pressure-free-ballistic-flight",
            statement=(
                "Fireball elements fly in straight lines at fixed speed after the merge; only "
                "rho v^2 acts, never P."
            ),
            source="prior work's ballistic model, PRD §2.2",
            applies_to=("beta_bare", "beta_flat", "beta_finite", "beta_mirror"),
            verdict="departs",
            why=(
                "At k ~ 7 the fireball's sound speed (~9.8 km/s) and its recoil velocity "
                "(9.31 km/s) agree to within 5%, and terminal Mach at the transport hand-off is "
                "~2.5 — the worst possible regime in which to omit pressure. Closing the 1.2-point "
                "gap this omission leaves is the study's whole purpose (PRD §3.4)."
            ),
            retired_by="Rung 1 (1-D spherical screen) and Rung 3/4 (2-D resolved transport)",
        ),
        Assumption(
            ident="isotropic-emission",
            statement="The merged blob expands isotropically at a single speed u in its own frame.",
            source="prior work's ballistic model, PRD §2.2",
            applies_to=("beta_bare", "beta_flat", "beta_finite", "beta_mirror"),
            verdict="departs",
            why=(
                "An axial projectile burying itself in a slug deposits energy anisotropically, and "
                "§6.2's snowplow screen shows deposited-energy fraction spanning 34-75% depending "
                "on projectile geometry alone. A single-speed shell also has zero velocity "
                "variance, whereas the real spread is 8-20 km/s (§6.4)."
            ),
            retired_by="Rung 1 (prescribed-deposition bracket), Rung 3 (resolved penetration)",
        ),
        Assumption(
            ident="complete-inelastic-merge-and-vaporisation",
            statement="Projectile and slug merge completely and everything vaporises.",
            source="PRD §2.2",
            applies_to=("beta_bare", "beta_mirror", "beta_ideal"),
            verdict="holds-with-caveat",
            why=(
                "The energy budget supports it — 305.7 MJ/kg is 57.1 eV per H2O molecule, far past "
                "dissociation plus first ionisation — but 'merges completely' is exactly what the "
                "§6.2 penetration screen contradicts for a slender projectile, which punches "
                "through and leaves the slug a spectator."
            ),
            retired_by="Rung 1 material qualification; Rung 3 resolved penetration",
        ),
        Assumption(
            ident="perfect-collimation-by-the-plate",
            statement=(
                "The plate redirects each captured element to -z at its full arrival speed "
                "(Delta_p = |v| + v_z), i.e. a focus-matched paraboloid."
            ),
            source="prior work's beta_bare, PRD §3.6",
            applies_to=("beta_bare", "beta_finite"),
            verdict="departs",
            why=(
                "This is an idealisation *already inside* the headline 984 s, not upside on top of "
                "it: the flat-plate counterpart is beta_flat = 0.517 and 560 s at the same k. A "
                "real stagnating plenum returns some fraction of that prize, and how much is what "
                "the shape sweep asks."
            ),
            retired_by="Rung 4 plate sweep (flat plus the paraboloid family, both tapers)",
        ),
        Assumption(
            ident="infinite-plate-radius",
            statement="Everything with v_z > 0 is caught, at any radius.",
            source="prior work's beta_bare, PRD §3.6",
            applies_to=("beta_bare", "beta_flat"),
            verdict="departs",
            why=(
                "Rays near the capture threshold arrive nearly grazing and land at large "
                "radius, so "
                "a finite plate loses them: capture falls 31.2% -> 10.6% at the reference "
                "R/d = 1.5. The headline number is an R -> infinity idealisation and must be "
                "quoted as the upper edge of a bracket (§13.13)."
            ),
            retired_by=(
                "Rung 4 plate-radius sweep; `capture_convention_bracket()` states the bracket"
            ),
        ),
        Assumption(
            ident="rigid-non-recoiling-plate",
            statement="The plate is rigid during the pulse and does not recoil appreciably.",
            source="inherited from the f(v) study, PRD §12",
            applies_to=("beta_bare", "beta_flat", "beta_finite"),
            verdict="holds",
            why=(
                "Plate recoil is 34 m/s per pulse against ~20 km/s gas — a factor of ~600 in hand "
                "— and a 30 m plate's first flexural mode is 10-100 ms against a 750 us pulse "
                "(§4.1, §6.4). This one transfers cleanly and passes trivially."
            ),
            retired_by="n/a — verified analytically at this study's loads",
        ),
        Assumption(
            ident="no-radiative-loss",
            statement="No energy leaves the system as radiation before the impulse is delivered.",
            source="ballistic model; PRD §7.2 Tier 1",
            applies_to=("beta_bare", "beta_mirror", "beta_ideal"),
            verdict="holds-with-caveat",
            why=(
                "Defensible for the *near field*, which is optically thick, and radiative loss can "
                "only lower the realized impulse, so the ceiling stays a ceiling. But at the plate "
                "tau_opt spans 0.63-63 and straddles 1, where flux-limited diffusion is weakest — "
                "and this project has already been burned by a 2000x opacity error there (§5.3)."
            ),
            retired_by="Rung 4's coupled rad-hydro spot-check (D6 gate)",
        ),
        Assumption(
            ident="single-closing-speed",
            statement="One closing speed w = 75 km/s, rather than the 74-81 km/s envelope.",
            source="PRD §0.2",
            applies_to=("beta_bare", "beta_finite", "beta_mirror"),
            verdict="holds",
            why=(
                "beta is dimensionless and independent of w, so w enters Isp and projectile "
                "economy "
                "linearly and only as a stated scale. The one place it must not be conflated is "
                "prior work's w_bar = 77.28 km/s, kept separate here."
            ),
            retired_by="n/a — parameterised, not assumed",
        ),
        Assumption(
            ident="cauchy-schwarz-ceiling",
            statement=(
                "All carried material starts at rest in the vehicle frame and ejecta energy is at "
                "most the projectile's incoming kinetic energy, giving j_max = w(sqrt(1+K_ej)-1)."
            ),
            source="PRD §3.2",
            applies_to=("beta_ideal", "v_e_max", "realization_fraction"),
            verdict="holds",
            why=(
                "A thermodynamic identity, not a model: it bounds any hydrodynamics whatsoever, "
                "which is precisely why it is the study's scoring reference. Its one real "
                "assumption is that no *other* energy source enters — true here, since the vehicle "
                "supplies none."
            ),
            retired_by="n/a — the invariant every rung is checked against",
        ),
        Assumption(
            ident="infinite-mass-mirror",
            statement="The tamper reverses each element at unchanged speed without recoiling.",
            source="prior work's upper bound, PRD §2 of the handoff",
            applies_to=("beta_mirror",),
            verdict="departs",
            why=(
                "Self-contradictory at finite mass — 'reflects all energy, absorbs none, still "
                "recoils'. Retained *only* as an upper bound (it reaches 96.9% of the K_ej = k "
                "ceiling), never as a tamper model. ADR-0030 replaces it with the isentropic "
                "piston framing and credits the tamper's recoil instead of counting it as loss."
            ),
            retired_by="ADR-0030; measured by Rungs 1-3",
        ),
        Assumption(
            ident="pass1-excludes-ablator",
            statement="K_ej = K = C: no ablator mass in the ceiling or the Isp denominator.",
            source="PRD §3.1, D4",
            applies_to=("beta_ideal", "isp_eff", "realization_fraction"),
            verdict="holds-with-caveat",
            why=(
                "Deliberate and labelled: it isolates the tamper question from the plate question. "
                "But the ablator is uncertain by 27x and could be 2-60% of the mass budget, so "
                "*every* Pass-1 Isp is an upper bound and must be quoted as one."
            ),
            retired_by="Rung 6 (Pass 2), which re-optimises every family with its measured ablator",
        ),
    ]


# ---- Report generation ---------------------------------------------------------------------------


def _write_csv(path: Path, rows: Iterable[object]) -> int:
    """Write dataclass rows to CSV, taking the header from the dataclass fields."""
    materialised = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialised:
        path.write_text("")
        return 0
    names = [f.name for f in fields(materialised[0])]  # type: ignore[arg-type]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in materialised:
            writer.writerow([_cell(getattr(row, name)) for name in names])
    return len(materialised)


def _cell(value: object) -> object:
    """Render one field for CSV. Sequences become `a; b; c` rather than a Python repr."""
    if isinstance(value, tuple | list):
        return "; ".join(str(item) for item in value)
    return value


@dataclass(frozen=True)
class AnchorRow:
    """One PRD §8 analytic anchor: what the document quotes against what the ledger computes."""

    anchor: str
    quoted: str
    computed: float
    agrees: bool
    note: str


def anchor_table() -> list[AnchorRow]:
    """The §8 analytic anchors recomputed, with residual disagreements stated, not hidden."""
    k = K_BARE_OPTIMUM
    reflected = free_plate_reflected_fraction(tau_t=1.0, k=k)
    return [
        AnchorRow(
            anchor="ceiling, j = w(sqrt(1+K_ej)-1) at K_ej = 7.06",
            quoted="1.839",
            computed=beta_ideal(k),
            agrees=abs(beta_ideal(k) - 1.839) < 5e-4,
            note="identity; equals prior work's beta_ideal exactly",
        ),
        AnchorRow(
            anchor="bare ballistic limit, beta_bare(7.06)",
            quoted="0.9087",
            computed=beta_bare(k),
            agrees=abs(beta_bare(k) - 0.90871) < 1e-4,
            note="confirms the inherited ballistic model",
        ),
        AnchorRow(
            anchor="ballistic capture fraction at k = 7.06",
            quoted="0.3118",
            computed=ballistic_capture_fraction(k),
            agrees=abs(ballistic_capture_fraction(k) - 0.3118) < 5e-5,
            note="R -> infinity convention; imposes no plate radius",
        ),
        AnchorRow(
            anchor="k <= 1 ballistic zero-capture floor",
            quoted="0",
            computed=ballistic_capture_fraction(1.0),
            agrees=ballistic_capture_fraction(1.0) == 0.0,
            note="ballistic-model limit; hydrodynamic thrust there is a Rung 1 measurement",
        ),
        AnchorRow(
            anchor="free-plate elastic bound at tau_t = 1",
            quoted="0.274",
            computed=reflected,
            agrees=abs(reflected - 0.274) < 5e-3,
            note=(
                "exact 1-D elastic value from §6.3's own mass pairing is 0.2732; §6.3 quotes the "
                "agreement to two figures (0.27), so 0.274 is a rounding artifact, not a model "
                "difference"
            ),
        ),
        AnchorRow(
            anchor="k -> 0 degeneracy (no thrust as slug mass vanishes)",
            quoted="0",
            computed=beta_bare(1e-9),
            agrees=beta_bare(1e-9) == 0.0,
            note="a code that produces thrust at k = 0 is wrong",
        ),
        AnchorRow(
            anchor="prior work's Isp at w_bar = 77.28 km/s",
            quoted="1014 s",
            computed=isp_eff(beta=beta_bare(k), c_ratio=k, w=W_BAR_PRIOR_MS),
            agrees=abs(isp_eff(beta=beta_bare(k), c_ratio=k, w=W_BAR_PRIOR_MS) - 1014.0) < 1.0,
            note="inherited ballistic model reproduced before anything is built on it",
        ),
        AnchorRow(
            anchor="bare Pass-1 optimum k*_bare",
            quoted="7.060",
            computed=optimal_k_bare(),
            agrees=abs(optimal_k_bare() - 7.060) < 5e-4,
            note="prior work's 7.057 is the same optimum to optimiser tolerance",
        ),
    ]


def main(argv: Sequence[str] | None = None) -> int:
    """Regenerate every closed-form table this study quotes, as CSV plus a stdout summary."""
    parser = argparse.ArgumentParser(
        prog="python -m puffsat.tamper.ledger",
        description=(
            "Rung 0 of the tamped-nozzle effective-Isp study (puffsat_tamper_isp_prd.md) — the "
            "analytic reference ledger. A different study from the f(v) per-collision work."
        ),
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT_DIR, help="output directory")
    parser.add_argument("--w", type=float, default=W_CLOSING_MS, help="closing speed [m/s]")
    args = parser.parse_args(argv)
    out: Path = args.out
    w: float = args.w

    comparison = reference_comparison(w=w)
    ratios = mass_ratio_table(w=w)
    captures = capture_convention_bracket(w=w)
    example = worked_example(w=w)
    anchors = anchor_table()
    audit = assumption_audit()
    reference_ledgers = [MassLedger.from_ratios(k=K_BARE_OPTIMUM, tau_t=tau) for tau in (0.0, 1.0)]

    _write_csv(out / "ledger_reference_comparison.csv", comparison)
    _write_csv(out / "ledger_mass_ratio.csv", ratios)
    _write_csv(out / "ledger_capture_conventions.csv", captures)
    _write_csv(out / "ledger_worked_example.csv", [example])
    _write_csv(out / "ledger_anchors.csv", anchors)
    _write_csv(out / "ledger_assumptions.csv", audit)
    _write_csv(out / "ledger_dish_depth.csv", dish_table(w=w))
    _write_csv(out / "ledger_plate_soak.csv", [plate_soak_chain()])
    _write_csv(out / "ledger_regenerative.csv", regenerative_budget())

    print("Rung 0 — tamped-nozzle analytic reference ledger  (PRD §10)")
    print(f"  w = {w / 1000.0:.2f} km/s, g0 = {G0} m/s^2, reference encounter {M_ENC_REF_KG:g} kg")
    print()

    print("§4  reference masses per pulse")
    for led in reference_ledgers:
        print(
            f"  tau_t={led.tau_t:g}: m_i={led.m_i:6.2f} kg  slug={led.m_s:6.1f} kg  "
            f"tamper={led.m_t:6.1f} kg  K={led.k_hydro:5.2f}  "
            f"E={projectile_energy(led.m_i, w) / 1e9:5.1f} GJ"
        )
    fireball = Fireball.from_k(K_BARE_OPTIMUM, w=w)
    print(
        f"  fireball: V={fireball.v_cm / 1000:.2f} km/s  u={fireball.u / 1000:.2f} km/s  "
        f"e={fireball.specific_internal_energy / 1e6:.1f} MJ/kg "
        f"({ev_per_water_molecule(fireball.specific_internal_energy):.1f} eV/molecule)"
    )
    print()

    print("§3.4  reference comparison (Pass 1 — ablator excluded, so every Isp is an UPPER BOUND)")
    print(f"  {'configuration':<48}{'K':>7}{'beta':>9}{'% ceiling':>11}{'Isp [s]':>10}")
    for row in comparison:
        print(
            f"  {row.label:<48}{row.k_ratio:>7.2f}{row.beta:>9.4f}"
            f"{row.realization * 100:>10.1f}%{row.isp_s:>10.0f}"
        )
    print()

    print("§3.5  mass-ratio trade (bare ballistic family)")
    print(f"  {'K':>7}{'beta':>9}{'Isp [s]':>10}{'vs peak':>10}{'E/J [kJ/N.s]':>15}{'vs peak':>10}")
    for ratio_row in ratios:
        print(
            f"  {ratio_row.k:>7.2f}{ratio_row.beta:>9.3f}{ratio_row.isp_s:>10.0f}"
            f"{ratio_row.isp_vs_peak * 100:>9.1f}%{ratio_row.energy_per_impulse_j / 1000:>15.1f}"
            f"{ratio_row.energy_vs_peak * 100:>9.0f}%"
        )
    print(f"  k*_bare = {optimal_k_bare(w=w):.4f}")
    print()

    print("§3.6/§13.13  capture conventions at k = 7.06 (the three figures in circulation)")
    print(
        f"  {'R/d':>8}{'ray capture':>14}{'rim-angle':>12}{'beta (par)':>12}{'Isp':>8}"
        f"{'beta (flat)':>13}{'Isp':>8}"
    )
    for cap in captures:
        label = "inf" if math.isinf(cap.r_over_d) else f"{cap.r_over_d:.1f}"
        print(
            f"  {label:>8}{cap.capture_ray * 100:>13.1f}%{cap.capture_rim_angle * 100:>11.1f}%"
            f"{cap.beta_parabolic:>12.3f}{cap.isp_parabolic_s:>8.0f}"
            f"{cap.beta_flat:>13.3f}{cap.isp_flat_s:>8.0f}"
        )
    print("  convention: quote the ray value as the lower edge and R -> inf as the upper edge of")
    print("  an explicit bracket; the rim-angle form is retired (see capture_convention_bracket).")
    print()

    print("§8  analytic anchors")
    for anchor in anchors:
        mark = "ok " if anchor.agrees else "XX "
        print(
            f"  {mark}{anchor.anchor:<46}quoted {anchor.quoted:<10}computed {anchor.computed:.5g}"
        )
    failures = [a for a in anchors if not a.agrees]
    print()

    print("§6.6  dish depth: capture is set by where the RIM is, not by the vertex standoff")
    print(
        f"  {'configuration':<22}{'wall':>7}{'rim-src':>9}{'th_max':>8}{'capture':>9}"
        f"{'beta':>7}{'Isp':>6}  |{'+tamper':>9}{'Isp':>6}{'lost':>7}{'plate':>8}"
    )
    for d in dish_table(w=w):
        print(
            f"  {d.label:<22}{d.wall_height_m:>6.1f}m{d.rim_vs_source_m:>8.1f}m"
            f"{d.theta_max_deg:>7.1f}d{d.capture * 100:>8.1f}%{d.beta_bare:>7.3f}"
            f"{d.isp_bare_s:>6.0f}  |{d.beta_tamped:>9.3f}{d.isp_tamped_s:>6.0f}"
            f"{d.turned_and_missed * 100:>6.1f}%{d.plate_mass_t:>7.1f}t"
        )
    print(f"  plate mass ceiling {PLATE_MASS_CEILING_T:.0f} t — not binding anywhere in this range")
    print()

    soak = plate_soak_chain()
    print("§6.5  plate soak chain (capacity, not prediction — reached only if the curtain fails)")
    print(
        f"  alpha_th={soak.alpha_th:.2g} m^2/s -> depth={soak.depth_m * 1e6:.0f} um  "
        f"areal={soak.areal_mass_kg_m2:.2f} kg/m^2  capacity={soak.capacity_j_m2 / 1e6:.2f} MJ/m^2 "
        f"({soak.share_of_fluence * 100:.1f}% of fluence)  basis={soak.basis_j / 1e6:.0f} MJ"
    )
    conducted_hi = conducted_per_pulse_j_m2(delta_t=1000.0 - 905.0)
    print(
        f"  physical conducted-in soak while the ablator pins the surface: "
        f"{conducted_per_pulse_j_m2(delta_t=800.0 - 753.0) / INCIDENT_FLUENCE_J_M2 * 100:.3f}%"
        f"-{conducted_hi / INCIDENT_FLUENCE_J_M2 * 100:.3f}% of fluence "
        f"(T_abl 800-1000 K) — two orders below the capacity"
    )
    print(f"  free regenerative cap (no boiling): {regenerative_cap() * 100:.0f}% of the basis")
    print()

    print("inherited-assumption audit")
    for entry in audit:
        print(f"  [{entry.verdict:<18}] {entry.ident:<38} retired by: {entry.retired_by}")
    print()

    print(f"wrote 9 CSVs to {out}/")
    if failures:
        print(f"WARNING: {len(failures)} anchor(s) disagree with the PRD — reconcile before use.")
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
