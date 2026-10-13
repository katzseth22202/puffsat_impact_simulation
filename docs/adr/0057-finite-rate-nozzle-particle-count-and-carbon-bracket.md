# ADR-0057: The finite-rate nozzle marches the particle count, and brackets carbon between two limits

Status: accepted, 2026-10-13. Answers the paper repo's S26
(`Balloon-Pulse-Propulsion/docs/throat_kinetics_asks_for_impact_sim.md`); results in
[`walled_nozzle_throat_kinetics.md`](../walled_nozzle_throat_kinetics.md). Extends
[ADR-0050](0050-walled-thermal-nozzle-methane-chamber.md),
[ADR-0052](0052-near-term-chamber-feeds-janaf-and-carbon-branches.md) and
[ADR-0056](0056-hydrocode-asks-blowdown-clock-and-fixed-charge.md).

## Context

B6 (ADR-0056 item 6) debits the shifting-equilibrium impulse with a Bray sudden freeze of
`H + H + M`. The atomic H held at the first station where the Damkoehler number drops below one
is written off, and recombination after that point is ignored. That gives an upper bound on the
debit, not a value. S26 asks for a finite-rate quasi-1D expansion with real carbon and hydrogen
rates for three nozzles and the sphere.

Two facts constrain the answer:
- carbon has no evaluated three-body rate, and the provenance rule forbids inventing one
  (ADR-0050);
- the paper's `eta` is a blowdown figure, while the ask names only the peak state.

## Decision

1. **One kinetic variable: particles per kilogram, `nu`.** Bimolecular exchanges conserve it,
   and only three-body association and dissociation change it. Its only sink is `H + H + M` at
   `freeze.mixture_coefficient` (Jacobs, Giedt & Cohen's within-study efficiencies), with the
   reverse rate from detailed balance on the EOS's own `K_H2`. No other three-body channel is
   given a rate.
2. **Composition is held off equilibrium by one multiplier `phi`.** Every mass-action constant
   becomes `K phi^(dnu)`. This is Keck's rate-controlled constrained equilibrium with the
   particle count as its single constraint (`eos_methane.constrained_composition`).
   `n_H^2 / n_H2 = K phi` exactly, so the net rate is `-k n_M K n_H2 (phi - 1) / rho`, with no
   cancellation near equilibrium.
3. **Carbon is bracketed by two composition laws, both limits:**
   - `frozen_carbon`: carbon frozen at the peak (ADR-0052's branch), with `phi` acting on the
     free H pool only. This is the floor.
   - `exchange_carbon`: constrained equilibrium. Carbon re-forms molecules at no kinetic cost,
     but only out of particles that `H + H + M` has removed. This is the upper realistic case.
   Full shifting equilibrium is the ceiling above both. `required_carbon_rate` reports the Bray
   threshold for carbon (the coefficient that keeps up to the half-release station) beside
   Slack's single `C + C + Ar` measurement. It is a threshold and is used for no verdict.
4. **The march.** It steps in `ln rho` (0.01) from `A/A* = 1.05` on the branch isentrope. It
   solves `(T, zeta = ln n_H^2/n_H2)` by Newton, with backward-Euler kinetics. Area comes from
   continuity, and the clock from a straight cone, `x = r* (sqrt(A/A*) - 1) / tan(theta)`, as
   in `near_term.freeze_check`.
   - The start is supersonic because the sonic point's area is flat in density, so `dx/u` is
     ill-posed there. The throat's Damkoehler number is reported.
   - Newton runs in `zeta`, not `ln phi`, because far from equilibrium the composition depends
     on `T` and `ln phi` almost only through `zeta`.
5. **Every parcel of the blowdown is marched** (14, geometric in chamber density down to 1e-3).
   The coarse quadrature is used only as a ratio against the same quadrature on the shifting
   exits, and that ratio is applied to the fine `efficiencies` value.
6. **Impulse is debited by the vacuum-thrust ratio `u + p / (rho u)`, not by `eta`.** At a fixed
   area ratio a partly frozen gas keeps a higher gamma, so its bare exit speed can exceed the
   shifting one while it leaves colder and at lower pressure. `eta` (kinetic energy only, the
   study's convention) is still reported, as the ask requests. The thrust ratio is the honest
   finite-to-shifting comparison, and it is ordered in the rate (tested).
7. **Wall angle: 10°, 15° and 20° rows.** The cone sets only the clock.

## Consequences

- B6's Bray debit is superseded where this module has rows. `blowdown.run_spec` still carries the
  Bray figure, and the ledger owed to the parent quotes this one.
- `near_term.FrozenCarbon.state` takes an optional `ln_phi`. At its default the branch is
  unchanged.
- Found while building this: `eos_methane.janaf_correction` holds `delta` below 1000 K but drops
  its slope, so `e(T)` steps by ~2e-4 of `h0` at 1000 K for the corrected molecules. The march
  brackets the step by bisection (`finite_rate._bracket`), so the error is spent once per parcel.
  The EOS is not changed here, because every earlier result rests on it.
