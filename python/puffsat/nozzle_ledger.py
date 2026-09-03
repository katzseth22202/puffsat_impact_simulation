"""Rung 0 -- the analytic reference ledger for the nozzle asks N1-N7.

**One calculator owns every closed-form number the N-asks quote**, so the paper, the answer
document, and the later rungs cannot drift apart. There is no solver here: it is exact algebra
over the paper's own stated model, plus the arithmetic that reconciles two of its numbers.
Downstream rungs quote closed-form figures only from here.

Source of the asks: `katzseth22202/Balloon-Pulse-Propulsion` @ `36080e1`,
`docs/nozzle_asks_for_impact_sim.md`, raised 2026-09-03. Working plan in
`todos/nozzle_asks_plan.md`.

---

## What this module settles, and why each is closed-form

Three of the seven asks turn out to be **partly tautological**, and saying so precisely is the
answer rather than a dodge:

- **`eps_b` (N3)** is `1/(gamma-1)` identically, for every pulse, whatever the magnet. Standoff
  sizing *defines* the bore field energy to be `pV`, and a monatomic plume holds `1.5 pV`. No
  retuning reaches Zakharov's 0.4, and no simulation moves it.
- **max beta (N3)** is 1 at every station by construction, because the paper's graded field is
  *derived from* the standoff condition. The answerable question is `p_actual/p_design`, which
  needs the solved expansion and belongs to Rung 5.
- **the reflection baseline (N1)** is a two-line average over an isotropic cone, and the whole
  N1 sensitivity table follows from `sqrt(2 alpha/pi)`.

The genuinely new closed-form result here is the **signed drift** of `f_d`: the paper applies one
formula to two legs that face opposite ways, and the sign it drops is worth four times the
correction it keeps.

## Sign and frame conventions

`+z` is the thrust direction: exhaust leaves along `-z`. Two legs, and they are not symmetric.

- **Leg 2, the Jupiter departure -- head-on.** The PuffSat arrives against the craft, so the
  merged centre of mass moves *retrograde*: already pointing out the back. A nozzle doing
  nothing still yields `+sqrt(f_d)`, but that drift delivers **zero net thrust**, because its
  momentum is exactly the debit the arriving PuffSat imposed.
- **Leg 1, the Earth growth push -- overtake.** The returning PuffSat catches the craft from
  behind, so the merged centre of mass moves *prograde*: the wrong way for the exhaust. A nozzle
  doing nothing yields `-sqrt(f_d)`, and the field must reverse the drift before it can profit
  from it.

`toll.py` already carries this asymmetry (the `-1` head-on against the `+1` overtake, Seth
2026-08-25). The paper does not: `eq:reflection_baseline` is derived once in a head-on framing and
`sec:two_leg_nozzle` never cites it.

## What this module does *not* do

No radial dimension, no field topology, no solver. `alpha` is an *input* here, not an output --
computing it is Rung 4's merge run, and this module only prices what a given `alpha` costs. The
all-space field-energy correction is likewise carried as the paper's stated bracket rather than
recomputed, because recomputing it needs a winding the paper never specifies (Rung 1).
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

# ---- Fixed constants: the flown Jupiter-only chain (sec:needle_through_fog, sec:minimum_nozzle)

MU0 = 4.0e-7 * math.pi
"""Vacuum permeability [H/m]."""

SIGMA_SB = 5.670374419e-8
"""Stefan-Boltzmann constant [W m^-2 K^-4]. Same value as `expansion.SIGMA_SB`."""

M_PROJECTILE_KG = 25.0
"""Ice projectile mass [kg]. Common to both legs (paper tex:986, tex:1570)."""

M_SLUG_KG = 213.0
"""Water slug mass [kg]. Common to both legs."""

K_FLOWN = 8.5
"""Flown slug ratio. **Both cadences and both legs** -- stated at tex:986, tex:2254, tex:1138.

Not to be confused with the pusher-plate payload boost ratio (2S 8.21-8.69, 3S 7.43-7.56), which
is a leg-1 plate quantity and not a slug ratio at all. Conflating them cost an hour on 2026-09-03.
"""

CADENCE_HZ = 2.0
"""Pulse repetition rate [Hz] (paper tex:2652, tex:2673 -- the only rate stated for this cycle).

Independently capped near 250 Hz by the ~4 ms blob round trip, so 2 Hz is a chosen design point
with ~125x of headroom, not a derived ceiling.
"""

W_HEAD_ON_MS = 75_000.0
"""Leg-2 closing speed [m/s]. The paper gives only "about 75 km/s", unsplit by cadence."""

BORE_RADIUS_M = 3.0
"""Throat/bore *radius* [m]. Confirmed a radius, not a diameter: the paper's 28 m^2 cross-section
and `r = sqrt(660/(pi*23)) = 3.02` both agree. Matches `expansion.THROAT_RADIUS`."""

COLUMN_LENGTH_M = 23.0
"""Column length [m] as the *paper* states it. Note `expansion.FIELD_LENGTH` is 23.8 m here; the
3.5% discrepancy is unreconciled and is reported rather than silently resolved."""

GRADED_FIELD_T: tuple[tuple[float, float], ...] = (
    (1.0, 20.0),
    (3.0, 12.0),
    (6.0, 9.0),
    (COLUMN_LENGTH_M, 5.0),
)
"""The paper's graded profile: `(z [m], B [T])` at four stations (tex:1129).

**This is the entire magnet specification in the paper.** No coil count, no radii, no axial
positions, no currents, no turn count. `12 T` and `9 T` appear nowhere in this repo, whose only
field model is `electrothermal.local_field`'s two-point `20 T -> 5 T` flux-conservation relation.
"""

GAMMA_MONATOMIC = 5.0 / 3.0

ZAKHAROV_THRESHOLD = 0.4
"""Rupture threshold `eps_b` inside a solenoid (Zakharov, via Schilling 2022)."""

EPS_B_ALLSPACE_BRACKET = (0.23, 0.52)
"""Fractional extra field energy outside the bore, from the paper's own all-space integral.

Carried as a stated input, not recomputed: the integral models the magnet as a current sheet at
the 3.0 m bore, and doing better needs a winding the paper never gives. Rung 1 replaces this.
"""

# ---- Radiative / thermal constants (paper sec:minimum_nozzle, tab:bag_sizing) --------------------

STRUCTURE_TEMP_K = 1500.0
"""Structure temperature the passive-radiator capacity is quoted at [K]."""

RADIATING_CAPACITY_W = 1.0e9
"""Passive shedding capacity [W]: "a few thousand square meters of structure at 1500 K sheds about
a gigawatt". Reproduced by `radiating_area_for_capacity()`."""

SKY_FRACTION = 0.1
"""Share of the sky filled by coils and structure -- "of order a tenth"."""

NEAR_SUN_PULSE_ENERGY_J = 1.0e16
"""Near-Sun collision energy over one periapsis burn [J] (tex:1578)."""

NEAR_SUN_BURN_S = 1000.0
"""Near-Sun periapsis burn duration [s] (tex:1578)."""

BOOKED_FLASH_W = 42.6e6
"""Intercepted flash the paper books for the Jupiter case [W] (tex:2673, tab:bag_state)."""

BAG_RADIATED_SHARE_BRACKET = (0.010, 0.036)
"""`tab:bag_sizing`'s radiated share at the flown 5.4 m bag, cold to hot end of the burn."""

Leg = Literal["head-on", "overtake"]
"""`head-on` = leg 2, the Jupiter departure. `overtake` = leg 1, the Earth growth push."""

DEFAULT_OUTPUT_DIR = Path("data/results")


# ---- The reflection baseline (N1, N7) -- eq:reflection_baseline ----------------------------------


def drift_fraction(k: float = K_FLOWN) -> float:
    """`f_d = 1/(1+k)` -- share of pulse energy carried as bulk drift (tex:2267).

    Only the projectile arrives with kinetic energy while the slug supplies mass alone, so the
    merged centre of mass moves at `w/(1+k)`. **This is an energy share and is therefore positive
    on both legs**; which way it points is `drift_sign`, and the paper's formula cannot see it.
    """
    if k < 0.0:
        raise ValueError("k must be non-negative")
    return 1.0 / (1.0 + k)


def drift_sign(leg: Leg) -> int:
    """`+1` head-on, `-1` overtake -- the sign the paper's single `f_d` drops.

    Head-on the merged centre of mass recedes *retrograde*, already aimed out the back. On the
    overtake it moves *prograde*, the wrong way, and the field has to reverse it first.
    """
    return 1 if leg == "head-on" else -1


def reflection_baseline(f_d: float) -> float:
    """`eq:reflection_baseline` -- what a perfect one-sided mirror returns, as a share of `v_g`.

    A particle's ship-frame axial velocity is `V + u mu` with `mu = cos(theta)` uniform on
    `[-1, 1]` for an isotropic expansion. Averaging `|V + u mu|` gives `V` when `V >= u` and
    `(u^2 + V^2)/(2u)` when `V < u`; with `V = sqrt(f_d) v_g` and `u = sqrt(1-f_d) v_g`,

        eta_refl = 1/(2 sqrt(1-f_d))   (f_d <= 1/2),      sqrt(f_d)   (f_d >= 1/2)

    and the two halves agree at `f_d = 1/2`. **The isotropy is the load-bearing assumption** --
    it is the `mu` uniform on `[-1,1]` step, and it is exactly what N1 exists to test.

    The magnitude is the *same on both legs*: `|-V + u mu|` has the same distribution as
    `|V + u mu|` under `mu -> -mu`. What differs is the floor it is measured against.
    """
    if not 0.0 <= f_d <= 1.0:
        raise ValueError("f_d must lie in [0, 1]")
    if f_d >= 0.5:
        return math.sqrt(f_d)
    return 1.0 / (2.0 * math.sqrt(1.0 - f_d))


def maxwell_boltzmann_baseline() -> float:
    """The drift-free baseline with heat drawn from a Maxwellian rather than one speed: `0.4607`.

    For a 3D Maxwellian of 1D width `sigma`, `<|v_z|> = sigma sqrt(2/pi)` while
    `v_g = sqrt(2e) = sigma sqrt(3)`, so the ratio is `sqrt(2/(3 pi))`. The paper quotes 0.46; it
    is the `alpha = 1/3` entry of `anisotropy_baseline`, and the two must agree by construction.
    """
    return math.sqrt(2.0 / (3.0 * math.pi))


def anisotropy_baseline(alpha: float) -> float:
    """`sqrt(2 alpha/pi)` -- the drift-free baseline for a Gaussian plume of axial share `alpha`.

    `alpha = <v_z^2>/<v^2>` over the thermal remainder, so `alpha = 1/3` is isotropic and smaller
    is a pancake. For an axisymmetric Gaussian, `<|v_z|> = sigma_z sqrt(2/pi)` and
    `v_g = sqrt(<v^2>)`, so the ratio is `sqrt(alpha) sqrt(2/pi)`.

    **This is N1's whole sensitivity in one line**, and it is a square root, which is why a 40%
    shortfall in `alpha` costs only 23% of the baseline rather than 40%.
    """
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("alpha must lie in [0, 1]")
    return math.sqrt(2.0 * alpha / math.pi)


def passthrough_floor(k: float = K_FLOWN, leg: Leg = "head-on") -> float:
    """What a nozzle doing *nothing* returns: `+sqrt(f_d)` head-on, `-sqrt(f_d)` on the overtake.

    The paper names this floor for the head-on case (`eta_jet = sqrt(m_rp)`, 0.5 at the
    three-to-one mix and 0.2 at `k = 24`) and never states its sign flip. On the overtake the
    bulk drift is prograde, so an inert nozzle is actively *negative*: the gas leaves the way the
    ship is going.
    """
    return drift_sign(leg) * math.sqrt(drift_fraction(k))


@dataclass(frozen=True)
class LegBaseline:
    """The reflection baseline and its floor on one leg -- the N1/N7 comparison, made per leg."""

    leg: Leg
    k: float
    f_d: float
    baseline: float
    """`eq:reflection_baseline`. Identical on both legs."""
    floor: float
    """Signed pass-through floor: what an inert nozzle returns."""
    nozzle_work: float
    """`baseline - floor` -- how far the field has to lift the plume to reach a plain mirror."""


def leg_baseline(k: float = K_FLOWN, leg: Leg = "head-on") -> LegBaseline:
    """Assemble one leg's baseline, floor, and the gap between them.

    The gap is the finding: at `k = 8.5` the head-on nozzle must supply 0.205 of `v_g` to match a
    mirror, and the overtake nozzle must supply 0.853 -- **four times the work for the same
    baseline**, because it has to turn the drift around before it counts.
    """
    f_d = drift_fraction(k)
    base = reflection_baseline(f_d)
    floor = passthrough_floor(k, leg)
    return LegBaseline(leg=leg, k=k, f_d=f_d, baseline=base, floor=floor, nozzle_work=base - floor)


# ---- eps_b, and why no run can move it (N3) ------------------------------------------------------


def epsilon_b_bore(gamma: float = GAMMA_MONATOMIC) -> float:
    """`eps_b = 1/(gamma-1)` counting bore field energy only -- **1.5, identically**.

    Standoff sizing sets `B^2/2mu0 = p`, so the field energy inside the bore is `pV`. A monatomic
    plume holds internal energy `pV/(gamma-1)`. Divide, and the volume cancels along with the
    pressure, the field, the propellant and the pulse: `eps_b = 1/(gamma-1)` for every pulse this
    architecture flies.

    **This is a tautology of the sizing rule, not a property of the magnet.** It is why N3 cannot
    be answered by retuning, and why the honest reply to "you are 3x short of Zakharov" is a
    statement about what the criterion measures rather than a better number.
    """
    if gamma <= 1.0:
        raise ValueError("gamma must exceed 1")
    return 1.0 / (gamma - 1.0)


def epsilon_b_allspace(extra_outside: float, gamma: float = GAMMA_MONATOMIC) -> float:
    """`eps_b` once field energy stored *outside* the bore is credited.

    The only lever that lowers `eps_b` at all. `extra_outside` is the fractional addition to the
    bore-only field energy, which the paper's current-sheet integral puts at 23-52%, bringing
    `eps_b` to 0.99-1.22 -- still 2.5-3x Zakharov's 0.4.
    """
    if extra_outside < 0.0:
        raise ValueError("extra_outside must be non-negative")
    return epsilon_b_bore(gamma) / (1.0 + extra_outside)


# ---- The standoff field, which is where 20/12/9/5 T comes from (N3, N5, N6) ----------------------


def standoff_pressure(b_field: float) -> float:
    """`p = B^2/(2 mu0)` [Pa] -- the plume pressure a field of `B` stands off."""
    return b_field * b_field / (2.0 * MU0)


def standoff_field(pressure: float) -> float:
    """`B = sqrt(2 mu0 p)` [T] -- the field required to stand off a plume pressure `p`.

    **This is the relation the paper's graded profile is an output of**: "working
    `eq:bore_from_length`'s geometry against the swept mass gives a field requirement near 20 T a
    meter in, 12 T at 3 m, 9 T at 6 m, and 5 T at the exit". So the four stated field values are
    a *pressure* profile in disguise, and `beta = p/(B^2/2mu0)` is 1 at every one of them by
    construction. N3's "max beta anywhere" is only meaningful as `p_actual/p_design`.
    """
    if pressure < 0.0:
        raise ValueError("pressure must be non-negative")
    return math.sqrt(2.0 * MU0 * pressure)


@dataclass(frozen=True)
class StandoffStation:
    """One station of the paper's graded profile, with the pressure it implies."""

    z_m: float
    b_field_t: float
    design_pressure_pa: float


def standoff_profile(
    stations: tuple[tuple[float, float], ...] = GRADED_FIELD_T,
) -> list[StandoffStation]:
    """The paper's four stations, each turned back into the design pressure it stands off.

    Reported so Rung 5 has a `p_design(z)` to divide the solved expansion's `p_actual(z)` by.
    """
    return [
        StandoffStation(z_m=z, b_field_t=b, design_pressure_pa=standoff_pressure(b))
        for z, b in stations
    ]


# ---- The radiative power balance (N4), and the reconciliation the paper never makes --------------


def pulse_energy(
    m_projectile: float = M_PROJECTILE_KG,
    m_slug: float = M_SLUG_KG,
    w: float = W_HEAD_ON_MS,
) -> float:
    """`E = (1/2) mu w^2` [J] on the reduced mass `mu = m_p m_s/(m_p + m_s)`.

    At the flown 25/213 kg and 75 km/s this is `mu = 22.37 kg` and `E = 62.9 GJ`, reproducing
    tex:1570.
    """
    mu = m_projectile * m_slug / (m_projectile + m_slug)
    return 0.5 * mu * w * w


def burn_power(pulse_j: float, cadence_hz: float = CADENCE_HZ) -> float:
    """Average burn power [W] -- pulse energy times repetition rate.

    The quantity the passive-structure gate is actually about. The paper computes it for the
    near-Sun case (10 TW) and never for the Jupiter case, which is why a near-Sun gate ended up
    quoted against Jupiter numbers in the ask.
    """
    return pulse_j * cadence_hz


def radiating_area_for_capacity(
    capacity_w: float = RADIATING_CAPACITY_W, temp: float = STRUCTURE_TEMP_K
) -> float:
    """Structure area [m^2] needed to shed `capacity_w` at `temp` -- reproduces "a few thousand".

    `A = P/(sigma T^4)`; at 1 GW and 1500 K this is 3483 m^2.
    """
    return capacity_w / (SIGMA_SB * temp**4)


def passive_structure_gate(burn_w: float, capacity_w: float = RADIATING_CAPACITY_W) -> float:
    """The allowed `(radiated share) x (sky fraction)` product for a burn of a given power.

    `gate = capacity/burn_power`. **The gate is a power balance and scales inversely with burn
    power**, which is the whole reconciliation: the paper's `1e-4` belongs to the near-Sun 10 TW
    burn, and the Jupiter-only chain runs 79x below that.
    """
    if burn_w <= 0.0:
        raise ValueError("burn power must be positive")
    return capacity_w / burn_w


def allowed_radiated_share(burn_w: float, sky_fraction: float = SKY_FRACTION) -> float:
    """Radiated share of pulse energy the gate permits at a stated sky fraction."""
    if not 0.0 < sky_fraction <= 1.0:
        raise ValueError("sky_fraction must lie in (0, 1]")
    return passive_structure_gate(burn_w) / sky_fraction


def implied_sky_fraction(absorbed_w: float, burn_w: float, radiated_share: float) -> float:
    """Sky fraction implied by a booked absorbed power, given a radiated share.

    Used to reconcile the paper's two routes to the same Jupiter heat load: `42.6 MW` booked
    against `tab:bag_sizing`'s radiated share. If the answer is far from the stated "of order a
    tenth", one of the three numbers is wrong.
    """
    if radiated_share <= 0.0:
        raise ValueError("radiated_share must be positive")
    return absorbed_w / (burn_w * radiated_share)


@dataclass(frozen=True)
class ThermalBalance:
    """One regime's passive-structure balance -- the N4 reconciliation, per case."""

    label: str
    pulse_energy_j: float
    cadence_hz: float
    burn_power_w: float
    gate_product: float
    """Allowed `(radiated share) x (sky fraction)`."""
    allowed_share: float
    """Allowed radiated share at the stated sky fraction."""


def thermal_balances() -> list[ThermalBalance]:
    """The near-Sun case beside the Jupiter case, which the paper never puts side by side.

    Near-Sun reproduces the paper's `1e-4` and its "below about a tenth of a percent". Jupiter,
    computed the same way from the cadence the paper states twice, gives a gate 79x looser. The
    ask carries the near-Sun figure across to the Jupiter numbers; the paper does not.
    """
    jupiter_pulse = pulse_energy()
    jupiter_power = burn_power(jupiter_pulse)
    near_sun_power = NEAR_SUN_PULSE_ENERGY_J / NEAR_SUN_BURN_S
    return [
        ThermalBalance(
            label="near-Sun, 4 R_sun periapsis burn",
            pulse_energy_j=NEAR_SUN_PULSE_ENERGY_J,
            cadence_hz=1.0 / NEAR_SUN_BURN_S,
            burn_power_w=near_sun_power,
            gate_product=passive_structure_gate(near_sun_power),
            allowed_share=allowed_radiated_share(near_sun_power),
        ),
        ThermalBalance(
            label="Jupiter-only chain, head-on departure",
            pulse_energy_j=jupiter_pulse,
            cadence_hz=CADENCE_HZ,
            burn_power_w=jupiter_power,
            gate_product=passive_structure_gate(jupiter_power),
            allowed_share=allowed_radiated_share(jupiter_power),
        ),
    ]


@dataclass(frozen=True)
class FlashDiscrepancy:
    """The paper's two routes to the same Jupiter heat load, and how far apart they are."""

    booked_w: float
    """`tab:bag_state`'s intercepted flash."""
    implied_w: float
    """What `tab:bag_sizing`'s radiated share at the stated sky fraction would give."""
    ratio: float
    booked_product: float
    implied_sky_fraction_at_hot: float
    implied_sky_fraction_at_cold: float
    margin_against_gate: float
    """How far the *booked* figure sits inside the Jupiter gate."""


def flash_discrepancy(
    share_bracket: tuple[float, float] = BAG_RADIATED_SHARE_BRACKET,
) -> FlashDiscrepancy:
    """Reconcile the booked 42.6 MW against `tab:bag_sizing`'s 1.0-3.6% radiated share.

    These are the same physical quantity by two routes and they disagree by an order of
    magnitude. Backing the sky fraction out of the booked figure gives ~1% at the hot end and
    ~3% at the cold end, against the "of order a tenth" the gate is stated with. So either the
    sky fraction is much smaller than a tenth for this nozzle, or the booked flash is low.

    **Either way the Jupiter case clears its own gate**, which is the load-bearing conclusion;
    the discrepancy matters because a reader can otherwise reach the opposite one.
    """
    cold_share, hot_share = share_bracket
    jupiter_power = burn_power(pulse_energy())
    implied = hot_share * SKY_FRACTION * jupiter_power
    booked_product = BOOKED_FLASH_W / jupiter_power
    return FlashDiscrepancy(
        booked_w=BOOKED_FLASH_W,
        implied_w=implied,
        ratio=implied / BOOKED_FLASH_W,
        booked_product=booked_product,
        implied_sky_fraction_at_hot=implied_sky_fraction(BOOKED_FLASH_W, jupiter_power, hot_share),
        implied_sky_fraction_at_cold=implied_sky_fraction(
            BOOKED_FLASH_W, jupiter_power, cold_share
        ),
        margin_against_gate=passive_structure_gate(jupiter_power) / booked_product,
    )


# ---- Assembled reference tables ------------------------------------------------------------------


def anisotropy_table(
    alphas: tuple[float, ...] = (1.0 / 3.0, 0.300, 0.250, 0.200, 0.100),
) -> list[tuple[float, float]]:
    """N1's sensitivity table: `(alpha, drift-free baseline)`.

    Reproduces the ask's 0.461 / 0.437 / 0.399 / 0.357 / 0.252.
    """
    return [(a, anisotropy_baseline(a)) for a in alphas]


def baseline_table(
    ks: tuple[float, ...] = (3.0, K_FLOWN, 24.0),
) -> list[LegBaseline]:
    """The reflection baseline on both legs across the slug ratios the paper quotes.

    At `k = 3`, 8.5 and 24 the baseline is 0.577, 0.529 and 0.510 -- reproducing tex:2267 -- while
    the *floor* moves 0.500, 0.324, 0.200 head-on and the negative of each on the overtake.
    """
    legs: tuple[Leg, ...] = ("head-on", "overtake")
    return [leg_baseline(k, leg) for k in ks for leg in legs]


def write_baselines(rows: list[LegBaseline], path: Path) -> None:
    """Write the per-leg baseline table as the deliverable's evidence artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["leg", "k", "f_d", "baseline", "floor", "nozzle_work"])
        for r in rows:
            writer.writerow(
                [
                    r.leg,
                    f"{r.k:g}",
                    f"{r.f_d:.6f}",
                    f"{r.baseline:.6f}",
                    f"{r.floor:+.6f}",
                    f"{r.nozzle_work:.6f}",
                ]
            )


def main() -> None:
    """Print every closed-form number the N-asks quote, beside the paper's stated value."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    print("== N1: the anisotropy sensitivity, sqrt(2 alpha/pi) ==")
    print(f"{'alpha':>8} {'baseline':>10}  shape")
    shapes = ["isotropic", "mildly oblate", "oblate", "pancake", "strong pancake"]
    for (alpha, base), shape in zip(anisotropy_table(), shapes, strict=True):
        print(f"{alpha:8.3f} {base:10.4f}  {shape}")
    print(
        f"\nMaxwell-Boltzmann drift-free baseline: {maxwell_boltzmann_baseline():.4f} "
        f"(paper prints 0.46)"
    )

    print("\n== N1/N7: the reflection baseline, per leg -- the sign the paper drops ==")
    print(f"{'k':>6} {'leg':>10} {'f_d':>8} {'baseline':>10} {'floor':>9} {'nozzle work':>12}")
    rows = baseline_table()
    for r in rows:
        print(
            f"{r.k:6.1f} {r.leg:>10} {r.f_d:8.4f} {r.baseline:10.4f} "
            f"{r.floor:+9.4f} {r.nozzle_work:12.4f}"
        )
    flown = [r for r in rows if r.k == K_FLOWN]
    ratio = flown[1].nozzle_work / flown[0].nozzle_work
    print(
        f"\nAt the flown k = {K_FLOWN}, the overtake nozzle does {ratio:.2f}x the work of the "
        f"head-on one\nto reach the same baseline. The paper applies one formula to both."
    )

    print("\n== N3: eps_b is a tautology of the sizing rule ==")
    lo, hi = EPS_B_ALLSPACE_BRACKET
    print(f"bore field energy only          : {epsilon_b_bore():.3f}  = 1/(gamma-1), exactly")
    print(f"crediting {hi:.0%} outside the bore  : {epsilon_b_allspace(hi):.3f}")
    print(f"crediting {lo:.0%} outside the bore  : {epsilon_b_allspace(lo):.3f}")
    print(
        f"Zakharov's threshold            : {ZAKHAROV_THRESHOLD:.3f}  "
        f"-> short by {epsilon_b_allspace(hi) / ZAKHAROV_THRESHOLD:.1f}"
        f"-{epsilon_b_allspace(lo) / ZAKHAROV_THRESHOLD:.1f}x"
    )

    print("\n== N3/N5/N6: the graded field is a pressure profile in disguise ==")
    print(f"{'z [m]':>7} {'B [T]':>7} {'p_design [MPa]':>16}")
    for st in standoff_profile():
        print(f"{st.z_m:7.1f} {st.b_field_t:7.1f} {st.design_pressure_pa / 1e6:16.2f}")
    print(
        "beta = p/(B^2/2mu0) is 1 at every station by construction; Rung 5 reports "
        "p_actual/p_design."
    )

    print("\n== N4: the passive-structure gate is a power balance ==")
    print(f"{'case':>38} {'burn power':>12} {'gate product':>14} {'allowed share':>14}")
    balances = thermal_balances()
    for b in balances:
        print(
            f"{b.label:>38} {b.burn_power_w / 1e9:9.1f} GW {b.gate_product:14.2e} "
            f"{b.allowed_share:13.2%}"
        )
    print(
        f"\nradiating area behind the 1 GW capacity: {radiating_area_for_capacity():.0f} m^2 "
        f'("a few thousand")'
    )
    print(
        f"Jupiter runs {balances[0].burn_power_w / balances[1].burn_power_w:.0f}x below the "
        f"near-Sun burn, so its gate is that much looser."
    )

    fd = flash_discrepancy()
    print("\n== N4: the paper's two routes to the same Jupiter heat load ==")
    print(f"booked intercepted flash (tab:bag_state) : {fd.booked_w / 1e6:8.1f} MW")
    print(
        f"implied by tab:bag_sizing at 3.6% x 0.1  : {fd.implied_w / 1e6:8.1f} MW "
        f"({fd.ratio:.1f}x the booked figure)"
    )
    print(
        f"sky fraction the booked figure implies   : "
        f"{fd.implied_sky_fraction_at_hot:.1%} (hot) to "
        f"{fd.implied_sky_fraction_at_cold:.1%} (cold), against a stated 10%"
    )
    print(
        f"booked figure sits {fd.margin_against_gate:.0f}x inside the Jupiter gate -- "
        f"the flown case clears it either way."
    )

    out = args.output_dir / "nozzle_baselines.csv"
    write_baselines(rows, out)
    print(f"\nwrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
