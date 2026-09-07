# Water-plate chemistry: same-state parcel-frozen restarts

Scope: a conditional reaction-switch experiment on ADR-0043's planar,
stratified, prepositioned-water reference. This does not implement a distributed
injector, finite plate, reaction network, wall losses or recoil.

## Counterfactual and bookkeeping

Restart the equilibrium and parcel-frozen branches from identical staggered
positions, nodal velocities, cell masses, specific energies and accumulated wall
impulse. Never reconstruct nodal velocities from cell averages. Each parcel
retains its own switch-state equilibrium species mix; a globally averaged mix
would not preserve the stratified state.

The frozen branch disables **all subsequent molecular and ionization reactions**,
including forward reactions during further compression. It keeps the fixed
chemical store, translational/rotational/vibrational thermal energy, and the same
provisional dense residual. This is not physical freeze-out, not a molecular-only
measurement, not a strict bound, and not the total benefit of chemistry. Earlier
chemistry, dissociation avoidance and the switch-state flow remain common to
both branches. Some combined-case parcels are still compressing at both switches.

Compare both combined and unsprayed configurations:

```
G_branch = J_combined,branch - J_unsprayed,branch - J_cold-only
Delta G = (J_combined,eq - J_combined,frozen)
        - (J_unsprayed,eq - J_unsprayed,frozen)
```

The analytic cold-only control is unchanged and cancels from Delta G. Require
four completed branches at the same resolution, switch and endpoint, identical
starting energy/impulse within each pair, finite diagnostics, <=0.5% restart
energy drift and momentum closure. A domain exit produces no matched difference.
Finite-window impulse is not settled impulse.

## Frozen thermodynamics and switch matching

Python exports species coefficients from the existing water partition functions.
Rust evaluates their thermal energies and heat capacities analytically, retaining
the molecular and ion stores as constant energy. The n=2 residual from ADR-0043
is represented by cubic Hermite interpolation of reference e, s, cv and their
log-density derivatives on 512 knots. Derive residual pressure and acoustics from
that potential, rather than independently interpolating p and e. The first-law
sound-speed identity and the Maxwell relation remain explicit tests.

The parent equilibrium table independently interpolates fields. At each
switch-state density and temperature, match it using two disclosed corrections:

```
Delta e = e_parent - e_source
Delta R = (p_parent - p_source)/(rho*T)
f_correction = Delta e - Delta R*T*ln(rho/rho_reference)
```

The constant energy offset adds no heat capacity. The entropy-pressure correction
adds `rho*Delta R*T` to pressure but no internal energy or heat capacity. The
actual restart state is unchanged. These numerical matching terms are not energy
injected at the switch and must not be counted as chemical release. Report their
magnitudes; reject summed absolute energy offsets above 0.5% of physical restart
energy or pressure matching corrections above 1%. These gates are diagnostics,
not a convergence proof or physical error bound. Parent table/energy-offset
refinement remains necessary, especially before interpreting small Delta G.

Cross-check the independent Rust source against Python before matching. Require
relative p/T jumps below 1e-8 on **both** restart branches; a wrong parent table
must not pass merely because the frozen branch was matched to its checkpoint.
Frozen parcels must stay at 1000 K–2 MK and within the table density domain.
Reject nonpositive pressure, heat capacity or squared sound speed. Energy
inversion may safeguard intermediate queries but cannot turn an out-of-domain
trajectory into a valid result. No cold frozen-mixture continuation is invented.

The new indexed EOS dispatch in `hydro1d` is limited to bare hydrodynamic
advancement/history with an explicit staggered restart. Radiation, ablation and
unported bounce diagnostics are not enabled for parcel-indexed closures;
unindexed calls fail explicitly. Uniform-EOS callers retain their old behavior.

## Accepted comparison and limits

At a common 0.5 ms endpoint, switches at 153.5 and 193 microseconds complete for
N=40 and N=80. Early-switch Delta G is 86.10 and 88.04 kN s respectively,
about 12.7% of the equilibrium water-layer gain. The late-switch values are
0.679 and 0.579 kN s: small differences of much larger responses, not a
precision-resolved chemistry benefit. See the study's Step 5 for artifacts.

The early N=40 frozen combined trajectory leaves its thermal domain at
0.6635 ms, so its 1 ms probe has **no** matched gain. The late N=40 1 ms
comparison remains valid but does not repair the early comparison. The result
supports a conditional contribution during wall loading in this closure, not
real-water nozzle performance. Separating molecular and ionization effects,
finite-rate validation, table convergence and the actual injection/plate geometry
remain open.
