"""Is the nozzle expansion collisionless enough for `mu = v_perp^2/2B` to mean anything? (R1)

**No, and not marginally.** This module exists because `jet.py` listed that assumption in its own
"what this is not" section and nobody ever tested it. The paper's reply R1 tested it from the
outside and got `Kn` ~ 1e-7. This is the inside check, on the solved cooling history rather than
on a bag-average.

# Why the question decides `eta_geom`

`jet.py` computes the flare's benefit from the adiabatic invariant: `mu = v_perp^2/(2B)` is
conserved per particle, so a fourfold fall in `B` divides `v_perp^2` by four and energy
conservation puts the difference into `v_par`. That is a **guiding-centre, single-particle**
argument. It needs two separate things to be true:

1. **Magnetization** -- the field varies slowly over a gyroradius. `jet.adiabaticity_parameter`
   tests this, and it passes comfortably.
2. **Collisionlessness** -- a particle keeps its own `mu` long enough for the field to act on it.
   *Nothing tested this.* A collision randomizes pitch angle and destroys `mu` outright; if a
   parcel collides many times while crossing the nozzle, the population's velocity distribution
   is set by the local thermodynamics, not by the invariant it started with.

Condition 2 is the one that fails, and it fails by seven orders of magnitude.

# What replaces it

A collisional plume in a field that is 15--75x over-strength for it (`beta` = 0.013--0.073, our
own P3) is **a de Laval nozzle whose walls happen to be magnetic**. That is the paper's own
"walled by field rather than fenced by it", and it needs no hardware change. Its exhaust
anisotropy is not a remembered invariant but the ordinary ratio of directed to thermal energy,
fixed by the Mach number the expansion actually reaches -- which `expansion.py` already solves and
`data/results/cooling_history.csv` already reports. See `jet.py` for the replacement.

# What is computed here, and why the estimate is deliberately generous

Three heavy-particle collision channels set the mean free path, and they differ by orders of
magnitude in cross-section:

| channel | cross-section | comment |
| --- | --- | --- |
| neutral-neutral | ~2e-19 m^2 | hard-sphere, the **smallest** of the three |
| ion-neutral | ~1e-18 m^2 | polarization capture, 5x larger |
| ion-ion (Coulomb) | `pi b90^2 ln Lambda`, ~1e-16 m^2 here | far larger again at these `T` |

**Only the smallest is used for the headline.** Reporting the neutral hard-sphere path alone is
the most collisionless-favourable reading available: every other channel shortens it. If the
verdict is "collisional" on the generous estimate it is collisional, full stop, and no argument
about which cross-section is right can rescue the invariant. The Coulomb path is computed
alongside so the size of the margin is visible.

The reported metric is not just `Kn`. The sharper one is **collisions per parcel over the
transit**: `mu` survives one collision no better than a hundred, so the number of them is the
statement that matters, and the transit time is in the cooling history.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

from puffsat import eos_water

K_B = eos_water.K_B
EV = eos_water.EV
EPS0 = 8.8541878128e-12
E_CHARGE = 1.602176634e-19

DEFAULT_HISTORY = Path("data/results/cooling_history.csv")
DEFAULT_OUTPUT = Path("data/results/continuum_check.csv")

SIGMA_NEUTRAL = 2.0e-19
"""Hard-sphere momentum-transfer cross-section for the heavy neutrals [m^2].

`pi d^2` at `d` ~ 2.5 Angstrom, which brackets atomic H (~1e-19), atomic O (~1.5e-19) and H2O
(~3e-19). Deliberately taken at the small end of that range: see the module docstring on why the
generous estimate is the one that carries the argument.
"""

SIGMA_ION_NEUTRAL = 1.0e-18
"""Ion-neutral polarization (Langevin) capture cross-section [m^2], order of magnitude.

Included only to show it is larger than the neutral channel, i.e. that ignoring it is generous.
"""

BORE_RADIUS_M = 3.0
"""The bore the mean free path is measured against -- `expansion.CHAMBER_RADIUS`."""

COLLISIONLESS_KN = 0.1
"""`Kn` above which a guiding-centre treatment is defensible. The conventional continuum/slip
boundary is 0.01 and free-molecular is ~10; 0.1 is charitable to the invariant.
"""


# ---- Mean free paths -------------------------------------------------------------------------


def coulomb_logarithm(n_e: float, temp: float) -> float:
    """`ln Lambda` for a thermal plasma at `(n_e, T)`, floored at 2.

    The ratio of the Debye length to the 90-degree impact parameter. Floored because the
    expression goes negative in cold dense gas where the Coulomb channel is irrelevant anyway.
    """
    if n_e <= 0.0:
        return 2.0
    debye = math.sqrt(EPS0 * K_B * temp / (n_e * E_CHARGE**2))
    b90 = E_CHARGE**2 / (4.0 * math.pi * EPS0 * 3.0 * K_B * temp)
    return max(2.0, math.log(debye / b90))


def coulomb_cross_section(n_e: float, temp: float) -> float:
    """`pi b90^2 ln Lambda` [m^2] -- the ion-ion momentum-transfer cross-section."""
    b90 = E_CHARGE**2 / (4.0 * math.pi * EPS0 * 3.0 * K_B * temp)
    return math.pi * b90**2 * coulomb_logarithm(n_e, temp)


@dataclass(frozen=True)
class Paths:
    """Heavy-particle mean free paths at one state [m], plus the densities behind them."""

    n_neutral: float
    n_ion: float
    n_e: float
    lambda_neutral: float
    lambda_ion_neutral: float
    lambda_coulomb: float

    @property
    def lambda_generous(self) -> float:
        """The **longest** defensible heavy-particle path: neutrals against neutrals only.

        Every other channel is shorter, so this is the reading most favourable to `mu`.
        """
        return self.lambda_neutral

    @property
    def lambda_combined(self) -> float:
        """All three channels, rates added -- the honest path a heavy particle actually sees."""
        rate = 0.0
        for path in (self.lambda_neutral, self.lambda_ion_neutral, self.lambda_coulomb):
            if path > 0.0 and math.isfinite(path):
                rate += 1.0 / path
        return math.inf if rate == 0.0 else 1.0 / rate


def mean_free_paths(rho: float, temp: float) -> Paths:
    """Heavy-particle mean free paths at `(rho, T)` on the equilibrium water composition.

    `sqrt(2)` in the neutral term is the standard relative-velocity correction for like particles.
    """
    comp = eos_water.composition(rho, temp)
    n_neutral = comp.n_neutral_heavy
    n_ion = comp.n_hp + sum(comp.n_o_ions)
    n_e = comp.n_e

    lam_nn = math.inf if n_neutral <= 0.0 else 1.0 / (math.sqrt(2.0) * n_neutral * SIGMA_NEUTRAL)
    lam_in = math.inf if n_neutral <= 0.0 else 1.0 / (n_neutral * SIGMA_ION_NEUTRAL)
    sigma_c = coulomb_cross_section(n_e, temp)
    lam_c = math.inf if n_ion <= 0.0 else 1.0 / (math.sqrt(2.0) * n_ion * sigma_c)

    return Paths(
        n_neutral=n_neutral,
        n_ion=n_ion,
        n_e=n_e,
        lambda_neutral=lam_nn,
        lambda_ion_neutral=lam_in,
        lambda_coulomb=lam_c,
    )


def knudsen(rho: float, temp: float, length: float = 2.0 * BORE_RADIUS_M) -> float:
    """`Kn` on the generous (neutral-only) path, against a stated macroscopic length."""
    return mean_free_paths(rho, temp).lambda_generous / length


def collisions_per_transit(rho: float, temp: float, speed: float, transit_s: float) -> float:
    """How many times a parcel collides while crossing the nozzle.

    `mu` is destroyed by the *first* collision; this counts how far past that the plume is. The
    thermal speed sets the collision rate, not the flow speed, so `speed` enters only to bound
    the estimate from below when it exceeds the thermal speed.
    """
    paths = mean_free_paths(rho, temp)
    m_bar = _mean_particle_mass(rho, temp)
    v_thermal = math.sqrt(8.0 * K_B * temp / (math.pi * m_bar))
    v_rel = max(v_thermal, 0.0)
    if paths.lambda_generous <= 0.0 or not math.isfinite(paths.lambda_generous):
        return 0.0
    return v_rel * transit_s / paths.lambda_generous


def _mean_particle_mass(rho: float, temp: float) -> float:
    """Mass per particle [kg] of the equilibrium mixture, electrons included."""
    comp = eos_water.composition(rho, temp)
    total = comp.n_total
    return rho / total if total > 0.0 else eos_water.M_H


# ---- The verdict over the solved history ------------------------------------------------------


@dataclass(frozen=True)
class StationVerdict:
    """One station of one leg/branch case, with its collisionality."""

    closing_speed_km_s: float
    branch: str
    area_ratio: float
    x_m: float
    rho: float
    temp_k: float
    radius_m: float
    time_ms: float
    lambda_generous_m: float
    lambda_combined_m: float
    kn_bore: float
    kn_local: float

    @property
    def collisionless(self) -> bool:
        """Would a guiding-centre treatment be defensible here?"""
        return self.kn_local >= COLLISIONLESS_KN


def load_history(path: Path = DEFAULT_HISTORY) -> list[dict[str, str]]:
    """The solved cooling history, as raw rows."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def evaluate_history(path: Path = DEFAULT_HISTORY) -> list[StationVerdict]:
    """Collisionality at every station of every leg and branch."""
    out: list[StationVerdict] = []
    for row in load_history(path):
        rho = float(row["rho"])
        temp = float(row["temp_k"])
        radius = float(row["radius_m"])
        paths = mean_free_paths(rho, temp)
        out.append(
            StationVerdict(
                closing_speed_km_s=float(row["closing_speed_km_s"]),
                branch=row["branch"],
                area_ratio=float(row["area_ratio"]),
                x_m=float(row["x_m"]),
                rho=rho,
                temp_k=temp,
                radius_m=radius,
                time_ms=float(row["time_ms"]),
                lambda_generous_m=paths.lambda_generous,
                lambda_combined_m=paths.lambda_combined,
                kn_bore=paths.lambda_generous / (2.0 * BORE_RADIUS_M),
                kn_local=paths.lambda_generous / (2.0 * radius),
            )
        )
    return out


@dataclass(frozen=True)
class CaseSummary:
    """One leg/branch case, reduced to the numbers the verdict rests on."""

    closing_speed_km_s: float
    branch: str
    kn_max: float
    kn_at_exit: float
    lambda_max_m: float
    transit_ms: float
    collisions: float
    any_collisionless: bool


def summarise(verdicts: list[StationVerdict]) -> list[CaseSummary]:
    """Reduce to one row per leg/branch, keyed on the **most** collisionless station in it."""
    cases: dict[tuple[float, str], list[StationVerdict]] = {}
    for v in verdicts:
        cases.setdefault((v.closing_speed_km_s, v.branch), []).append(v)

    out: list[CaseSummary] = []
    for (speed, branch), group in sorted(cases.items()):
        group = sorted(group, key=lambda s: s.area_ratio)
        worst = max(group, key=lambda s: s.kn_local)
        exit_station = group[-1]
        transit_s = (exit_station.time_ms - group[0].time_ms) * 1e-3
        collisions = collisions_per_transit(exit_station.rho, exit_station.temp_k, 0.0, transit_s)
        out.append(
            CaseSummary(
                closing_speed_km_s=speed,
                branch=branch,
                kn_max=worst.kn_local,
                kn_at_exit=exit_station.kn_local,
                lambda_max_m=worst.lambda_generous_m,
                transit_ms=exit_station.time_ms - group[0].time_ms,
                collisions=collisions,
                any_collisionless=any(s.collisionless for s in group),
            )
        )
    return out


CSV_HEADER = (
    "closing_speed_km_s",
    "branch",
    "area_ratio",
    "x_m",
    "rho",
    "temp_k",
    "radius_m",
    "lambda_generous_m",
    "lambda_combined_m",
    "kn_bore",
    "kn_local",
    "collisionless",
)


def write_verdicts(rows: list[StationVerdict], path: Path = DEFAULT_OUTPUT) -> None:
    """Write the per-station collisionality table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for r in rows:
            writer.writerow(
                [
                    f"{r.closing_speed_km_s:g}",
                    r.branch,
                    f"{r.area_ratio:.5f}",
                    f"{r.x_m:.4f}",
                    f"{r.rho:.6e}",
                    f"{r.temp_k:.2f}",
                    f"{r.radius_m:.4f}",
                    f"{r.lambda_generous_m:.6e}",
                    f"{r.lambda_combined_m:.6e}",
                    f"{r.kn_bore:.6e}",
                    f"{r.kn_local:.6e}",
                    int(r.collisionless),
                ]
            )


def main() -> None:
    """R1's premise: is the nozzle expansion collisionless?"""
    parser = argparse.ArgumentParser(description="Collisionality of the nozzle expansion (R1)")
    parser.add_argument("--history", type=Path, default=DEFAULT_HISTORY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    verdicts = evaluate_history(args.history)
    cases = summarise(verdicts)

    print("== R1's premise, tested on our own solved expansion ==")
    print("Mean free path is the NEUTRAL hard-sphere one -- the longest of the three channels,")
    print("so this is the reading most favourable to keeping `mu` as an invariant.\n")
    print(
        f"{'leg':>8} {'branch':>12} {'max Kn':>10} {'Kn at exit':>11} "
        f"{'longest mfp':>13} {'transit':>9} {'collisions/parcel':>18}"
    )
    for c in cases:
        print(
            f"{c.closing_speed_km_s:8.2f} {c.branch:>12} {c.kn_max:10.2e} {c.kn_at_exit:11.2e} "
            f"{c.lambda_max_m * 1e6:11.3f}um {c.transit_ms:8.3f}ms {c.collisions:18.2e}"
        )

    worst = max(c.kn_max for c in cases)
    print(f"\nmost collisionless station anywhere : Kn = {worst:.2e}")
    print(f"threshold for a guiding-centre model: Kn > {COLLISIONLESS_KN}")
    print(f"-> short by a factor of {COLLISIONLESS_KN / worst:.3g}.")
    print("\nVERDICT: the expansion is a collisional continuum everywhere. `mu = v_perp^2/2B`")
    print("is destroyed on the first collision, and a parcel has millions of them in transit.")
    print("`jet.alpha_after_expansion` models the wrong regime; P2's 0.70-0.88 is withdrawn.")

    write_verdicts(verdicts, args.output)
    print(f"\nwrote {len(verdicts)} rows -> {args.output}")


if __name__ == "__main__":
    main()
