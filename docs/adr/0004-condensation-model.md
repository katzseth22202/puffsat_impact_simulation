# Condensation: two channels, equilibrium baseline, kinetics gated

The low-v worst case (water at 3.2 km/s) is condensation-dominated. The model splits the loss
into two physically distinct channels and treats the kinetic refinement as deferred-and-gated.

**Two channels** (they have *opposite* latent-heat bookkeeping, so they must be kept separate):

- **Bulk vapor-pressure collapse** — vapor condensing in the gas removes the pressure that drives
  re-expansion (weaker bounce), but the latent heat (~2.26 MJ/kg, ~44% of the 5.1 MJ/kg impact KE
  at 3.2 km/s) returns *to the gas* and the condensate co-moves and can still rebound as a mist.
  Handled by a **two-phase equilibrium EOS** (saturation curve + latent heat in the energy equation).
- **Wall deposition** — condensate that reaches the wall and sticks loses its mass, momentum, *and*
  latent heat to the wall. The hard momentum sink. Handled by a **wall sticking boundary condition**,
  sticking coefficient `α` (baseline `α = 1`, most pessimistic; sensitivity reported).

**Equilibrium condensation is the conservative baseline.** Instantaneous condensation to saturation
condenses the most and removes the most vapor pressure → lowest `e_eff`. If water at 3.2 km/s clears
a useful `f` under this pessimistic assumption, the low-v case is settled (every other case is
easier) and the kinetic model is never built. The **kinetic nucleation/growth model is deferred
behind that useful-`f` gate**, because kinetic inhibition (condensation slower than the ~µs
residence time) can only *raise* `e_eff`. A condensation-timescale-vs-residence estimate up front
indicates which regime applies.

**Why wall deposition stays active despite a nominally hot face.** During the ~µs pulse the
SiC/ablator surface sits within ~25 K of its *cold* initial temperature, not the 1700 K gas: the
contact-interface temperature is governed by the effusivity ratio `e_solid/e_gas ~ 10³`
(`√(kρc)`: ablator ~550, water vapor ~10), so the solid dominates and the surface stays below
water's saturation temperature (`T_sat ≈ 520–620 K` at 4–16 MPa stagnation). Net wall condensation
is therefore favored — confirming condensation as the dominant low-v sink. Between-pulse heat soak
would change this but is out of scope (§11).

**Amendment (Rung C, 2026-06): the wall-deposition channel is conduction-gated.** The two channels
turn out to be *sequenced*, not independent. Rung C implemented both — a CoolProp two-phase EOS (bulk
channel: `p → p_sat`, latent heat in `e`) and a true mass-sink wall-sticking BC (deposition channel) —
and swept water at 3.2 km/s with conduction deferred (the B-flux gas-side-resistance gap, ADR-0005).
Result: `e_eff(ρ) ≈ 0.78→0.735`, with **wall deposition ≈ 0**. The diagnostic shows why and confirms
this ADR's own reasoning: wall condensation requires the near-wall gas to be cooled below `T_sat`,
which (per the effusivity argument above) is the *cold wall doing the cooling* — i.e. **conduction**.
With conduction off, the adiabatic wall cell stays superheated and never condenses. Bulk condensation
*does* occur (~25 % in the cool re-expansion tail, captured by the EOS) but sits in the low-pressure
tail, so it barely moves `e_eff` — which is why 3.2 km/s bounces *better* (0.74) than 16 km/s (0.63)
in the adiabatic model, *contrary* to the worst-case framing. **Consequence:** the condensation-
dominated low-v worst case only materializes once **B-flux** adds the wall-cooling conduction that
drives deposition; until then Rung C's `e_eff` is an adiabatic upper bound. The deposition machinery
(`liquid_frac` table field + `CondensingBounce` sink) is built and verified, dormant until then.

## Considered Options

- **Single lumped "condensation loss."** Rejected: blurs two mechanisms whose latent-heat
  bookkeeping is opposite (bulk returns heat to the gas; wall deposition loses it).
- **Kinetic condensation from the start.** Rejected as baseline: equilibrium is the conservative
  bound and kinetics only helps, so it is built only if the equilibrium baseline fails the gate.
