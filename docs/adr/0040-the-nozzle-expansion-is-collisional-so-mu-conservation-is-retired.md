# The nozzle expansion is a collisional continuum, so `mu` conservation is retired and `eta_geom` is the de Laval decomposition

Status: accepted (2026-09-05, working reply R1 of
`katzseth22202/Balloon-Pulse-Propulsion` `docs/nozzle_replies_to_impact_sim.md` @ `b6dba6f`).
Supersedes the physical model of `jet.py`, which is retired in place. Answers P11 and P12 of
`docs/nozzle_replies_answered.md`.

## What was wrong

`jet.py` computed the flare's benefit from the adiabatic invariant `mu = v_perp^2/(2B)`: a
fourfold fall in `B` divides `v_perp^2` by four, and energy conservation puts the difference into
`v_par`, giving `alpha_exit = 1 - (1 - alpha_0)/(A/A*)`.

That is a **guiding-centre, single-particle** argument, and it needs two independent things:

1. **Magnetization** — the field varies slowly over a gyroradius. `adiabaticity_parameter` tested
   this, and it passes comfortably.
2. **Collisionlessness** — a particle keeps its own `mu` long enough for the field to act on it.
   **This was never tested.** The module's own "what this is not" section listed it as an
   assumption and no code checked it.

The paper found it from the outside first. Tested here on the solved cooling history, on the
**longest** of the three heavy-particle collision channels (neutral hard-sphere, which is the
reading most favourable to keeping `mu`), the most collisionless station anywhere in the flown
expansion has

    Kn = 2.5e-7,   against the ~0.1 a guiding-centre treatment would want

and a parcel collides between **two and ten million times** crossing the nozzle. A collision
randomizes pitch angle and destroys `mu` outright, so the first one ends the argument.

## The decision

**`eta_geom` is computed as a continuum nozzle quantity**, in `extension.py`:

    eta_geom = <cos theta> / sqrt(1 + <v_th^2>/u^2),    <v_th^2> = 3 k T / m_bar

- `<cos theta>` is **mass-weighted over flux tubes traced through the real Biot–Savart field**
  (`fluxtube.py`), not through a paraxial expansion.
- `<v_th^2>` uses the mean **heavy-particle** mass from `eos_water.composition`, because
  electrons carry thermal energy but essentially no mass and `v_g` is energy per unit mass.
- `m_bar` rather than an assumed `gamma`, because the mixture's effective `gamma` varies through
  the expansion. This is algebraically identical to the paper's `<cos theta>/sqrt(1+3/(gamma M^2))`
  on a `gamma`-law gas, pinned by test.

**`jet.py` is retired in place**, not deleted: its docstring carries the retirement, its `ETA_CHEM`
and `TARGET_ETA_JET` remain current data, and the mirror half (leg 1's reflection) is flagged as
resting on the same collisionless assumption and needing re-derivation before it is quoted again.

## What the flux tubes keep, and why

Retiring `mu` does **not** retire flux tubes. `mu` is a statement about a *particle*; a flux tube
is a statement about the *fluid*, and its warrant is the magnetic Reynolds number, which
`cooling_history.csv` reports in the hundreds through the column. A fluid element stays on its
tube even though no particle remembers its `mu`.

This splits two things the previous framework conflated, and the distinction is load-bearing:

| | station-weighted? | why |
| --- | --- | --- |
| **which tube a parcel rides** (geometry) | **yes** | fixed by where it starts; it never sees the field upstream of that |
| **how much thermal energy it converts** (energetics) | **no** | fixed by the total pressure drop it falls through |

`mu` was the thing that tied them together. Without it, the paper's station weighting survives in
`fluxtube.py` and does not survive in its own R4.

## What changes numerically

`eta_geom` = **0.58–0.86** on the capped winding as flown, against `jet.py`'s withdrawn 0.70–0.88
and against the paper's own corrected estimate of 0.48–0.64. The thermal term is 0.86–0.97 — the
plume is already well converted at the exit — so the whole of the remaining loss is **pointing**,
not conversion, and it happens entirely downstream of the last coil.

## Consequences

- **The "deeper flare is always better" claim is inverted.** `<cos theta>` falls as the flare
  opens while the thermal term rises, so a magnetic nozzle has an optimum. `jet.py`'s docstring
  said the opposite.
- **ADR-0012's cap becomes a gain rather than a compromise** (P13): the smaller area ratio leaves
  the plume denser at the same exit field, which raises `M_A`, which detaches sooner, which fans
  less. `eta_geom` moves +0.06 to +0.15 across the cap.
- **The loss-cone criterion goes with it** — see ADR-0041.
- **The detachment window is worth more than P3 implied**: `<cos theta>` runs 0.95 at 1.0 exit
  radii and 0.69 at 2.0, so narrowing `M_A` at the exit is a 27% lever on `eta_geom`.

## Provenance

`make analysis-continuum`, `make analysis-nozzle-fluxtube`, `make analysis-nozzle-extension`.
Validated by `extension.continued_history` reproducing `data/results/cooling_history.csv` exactly
at the flown geometry before being run anywhere new.
