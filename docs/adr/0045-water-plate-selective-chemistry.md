# Water-plate chemistry: molecular-frozen, ionization-active intermediate closure

Scope: an ordered conditional decomposition of ADR-0044's same-state reaction
switches. All geometry, source-energy and provisional dense-water restrictions
from ADR-0043 remain. This is not a finite-rate network or a physical freeze-out
model, and does not implement a distributed injector or finite plate.

## Constrained equilibrium

Hold each parcel's H2O, OH, H2 and O2 abundances fixed at its switch state.
Also fix the number of H and O nuclei outside those molecules. Only their
atomic charge states and the electron population may equilibrate. Use the
existing H/H+ and complete O/O+..O8+ Saha constants and species partition
functions. The fixed molecules retain their translational, rotational and
vibrational energies. Their bond store stays constant; the ionization store
and electron thermal energy remain active.

At a specified density and temperature, solve one monotone charge-neutrality
equation for electron density. For an atomic ladder the stage weights are
proportional to `product(K_1..K_z) / n_e^z`; normalizing them conserves that
element's available atomic inventory. Solve in log electron density with stable
log weights. Python independently bisects the scalar equation; Rust uses a
safeguarded Newton solve on the hydrodynamic hot path.

Differentiate the normalized populations and implicit charge constraint
analytically. In particular, with electron abundance Q and the inventory-weighted
charge variance V, `d ln(n_e)/d ln(rho) = Q/(Q+V)`. Temperature derivatives use
covariances with `d ln(product K)/dT = 1.5*z/T + E_z/(k*T^2)`.
Include changing electron particle count in p, e and their derivatives;
ionization-energy derivatives alone would violate thermodynamics.
Use the first-law acoustic identity, not a fixed gamma.

Add the same dense Helmholtz residual and per-parcel parent-table matching
corrections as ADR-0044. The selective ideal mixture agrees with the original
equilibrium mixture at the switch. No species, conserved energy or flow state
jump is introduced. Retain the 1000 K–2 MK domain and all conservation,
positive-cv/acoustics and source cross-language acceptance gates.

## Ordered differences, not a unique additive attribution

Let E mean full equilibrium, M fixed molecules with atomic ionization active,
and F all species fixed. For the matched water-layer gain G, report

```
Delta G_total = G_E - G_F
Delta G_molecular | ions active = G_E - G_M
Delta G_ionization | molecules fixed = G_M - G_F
```

The two conditional differences telescope exactly to the old total. They are
**ordered**: a reverse-order constraint would in general allocate interactions
differently. They include forward and reverse reactions after the switch, and
the ensuing hydrodynamic response. Do not call them a uniquely measured split
of recombination work or extrapolate them to all chemistry before the switch.

Subtract the unsprayed response in every difference. A negative ionization
effect on G means ionization helps the unsprayed control more than the combined
configuration under this constraint; it need not mean negative wall impulse or
an endothermic recombination reaction. A small total can conceal sizable,
opposing conditional effects. Ratios of one effect to their small net sum are
not reaction-energy fractions or bounded percentages.

The ordered report requires matching parent metadata, actual staggered restart
arrays, storage offsets, area, residual knots, switch times, timestep fraction
and complete output grids, in addition to the independent branch gates.
Equilibrium final impulses must agree between the two experiments. Missing or
invalid branches cannot produce an ordered difference.

## Acceptance and outstanding work

Tests cover switch continuity against the full source EOS, element-inventory
and Saha/charge conservation, the analytic single-stage quadratic Saha limit,
multi-stage derivatives and Maxwell consistency, positive heat capacity and
acoustics, no available atomic inventory, and mismatched-report rejection.
Actual restarts independently cross-check Rust source p/e against Python before
applying parent interpolation corrections.

The study's Step 6 records the accepted finite-window flow comparison and its
resolution check. Molecular reaction rates, pressure falloff and dense-fluid
activity corrections remain unvalidated. Equilibrium is a closure assumption,
not a rate measurement. Before interpreting the selective switch as real-water
nozzle performance, validate finite-rate chemistry and its coupling to the
actual injected flow, with table refinement and wall/geometry losses retained.
