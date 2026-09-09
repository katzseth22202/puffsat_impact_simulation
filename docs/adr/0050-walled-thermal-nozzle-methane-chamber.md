# The walled thermal nozzle is its own study, and methane gets its own EOS

Status: **accepted**, 2026-09-09. Answers asks N9 and N10 of
`Balloon-Pulse-Propulsion` @ `6faf171`, `docs/nozzle_asks_for_impact_sim.md`,
raised 2026-09-08 against that repo's ADR-0016.

## Scope and what it is not

The walled thermal nozzle is a **third study** in this repository, beside `f(v)`
and the head-on tamper study, and separate again from the water-injected overtake
plate. It shares the EOS machinery, the Bray freeze criterion and the expansion
integrator; it shares **none** of the plate-side conventions, none of the magnetic
nozzle's geometry, and none of the `f(v)` bag.

Its geometry is ADR-0016's, not N1-N8's: a 3 m bore, a column of 200, 400 or
673 m^3, **no field**, 25 kg of impactor at 75 km/s into 474 kg of methane
pre-charged at 1.4 bar and 111 K. Nothing here inherits `expansion.BAG_RHO`,
`expansion.AREA_RATIO_EXIT` or the graded 20-5 T field, and a run that reads one
of those constants is reading the wrong study.

Its Python is `python/puffsat/walled_nozzle/`, its EOS is the top-level peer
`python/puffsat/eos_methane.py`, its outputs are under
`data/results/walled_nozzle/`, and its targets are prefixed `walled-nozzle-`.

## The chamber equilibrium is the first deliverable, and it is shared

N9 asks what the wall sees; N10 asks what the chemistry returns. Ask item N10.4b
makes them one question: both depend on how much of methane's atomisation store
is **charged** at chamber conditions, and that is an equilibrium composition at
`(rho, T)` and nothing else. So the composition solve lands first, and both asks
read it rather than each assuming a charge state.

`eos_methane.py` is built in `eos_water.py`'s frame and to the same standard:
partition functions from spectroscopic constants, dissociation energies from
NIST-JANAF **0 K** heats of formation, law of mass action for every molecule, a
Saha ladder for every ionisation stage, closed by element conservation and charge
neutrality. Three log-unknowns (`ln n_C`, `ln n_H`, `ln n_e`); every other species
follows.

**0 K heats of formation, deliberately.** ADR-0016's paper-side arithmetic used
298 K values (its 103.7 MJ/kg atomisation and 4.66 MJ/kg formation enthalpy are
the 298 K numbers). This side uses 0 K throughout, matching `eos_water`, and gets
102.35 and 4.15. The 1.3% gap is a reference-state convention, not a
disagreement, and it is stated rather than reconciled because the 0 K convention
is the one the partition functions are referenced to.

## Species set, and what it omits

Carried: `CH4`, `CH3`, `CH2`, `CH`, `C2H2`, `C2`, `C3`, `H2`, atomic `C` and `H`,
`H+`, the full `C+ .. C6+` Saha ladder, and `e-`.

Omitted, each for a reason that is stated rather than assumed away:

- **Condensed carbon.** Graphite is not a gas-phase equilibrium species and it
  does not form by a three-body collision, so it is not in the mass-action set.
  It is N10 item 5 and gets its own nucleation treatment.
- **Excited electronic states.** Ground terms only, as in `eos_water`. `C2`'s
  low-lying `a 3Pi_u` at 716 cm^-1 is the one place this bites; it understates
  `C2`, which never exceeds a few percent of the carbon here.
- **`C2H`, `C4`, `C5`, larger clusters.** They matter on the soot pathway rather
  than in the chamber, and they belong with the nucleation work.

## What a run may and may not claim

The composition solve reports equilibrium at a stated `(rho, T)`. It does **not**
decide whether that equilibrium is reached: that is the Damkoehler comparison of
N10 items 1-3, run separately, and a composition quoted without it is a charge
ceiling rather than a charge.

`k` reported from the wall-cap energy density is a **ledger identity** at the
adopted `eta_geom`, not a solved nozzle. It inherits ADR-0016's scoring
conventions so the two repositories' tables can be compared line for line, and it
must not be read as an independent confirmation of them.

## Amendment 2026-09-09: rate coefficients enter, under a provenance rule

N10 items 1-3 are answered, which required this study to carry three-body
recombination rate coefficients for the first time. They live in
`python/puffsat/walled_nozzle/rates.py` under the rule `recombination.py`
already states and this module inherits: **it computes no rate coefficient of
its own.** Every constant is a published value carrying its bath gas, its fitted
temperature window, the evaluation's own stated uncertainty, and the NIST record
identifier it was read from, so any number can be traced to its source without
this repository in the loop.

Three consequences are load-bearing and are recorded here rather than only in
docstrings.

**A rate belongs to a third body.** Published `H + H + M` coefficients differ by
an order of magnitude across `M = Ar`, `H2` and atomic `H`, and at the flown
10 kK this charge is 93-98% dissociated, so the third body is atomic hydrogen.
`freeze.mixture_coefficient` weights the species actually present rather than
adopting one coefficient; anything unmeasured is charged at argon's efficiency,
the slowest measured, so an unknown collider cannot flatter the result.

**The answer is a margin, not a point value.** Verdicts are reported as decades
of headroom in the rate coefficient, because that is the form in which a result
resting on measured constants can be honestly compared with their uncertainty.
A verdict whose margin does not exceed the coefficient's stated uncertainty is
not a verdict.

**Carbon gets no coefficient, and that is a finding.** The evaluated literature
holds one `C + C + M` measurement and no `C + H + M` at all. This is not a gap to
be filled by a harder search: the carbon store condenses rather than
associating, which is why nobody measured it. ADR-0050's exclusion of condensed
carbon therefore stands, and N10 item 5 remains a nucleation problem outside this
machinery. Adding an invented carbon rate to close the ledger is specifically
forbidden.

No no-field integrator was needed. `expansion.cooling_history` was already
parametrised by the EOS and the sound speed and knows nothing about a magnetic
nozzle; running it on `eos_methane` with ADR-0016's geometry is the walled case.
