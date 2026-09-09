# Water-plate kinetics: local neutral-network and applicability screen

Scope: test the rate margin assumed by Steps 4–6 on their **saved equilibrium
planar flow**. This is not a finite-rate hydrodynamic calculation or a validated
rate set for warm-dense water. No new impulse prediction is made.

## Rate source and thermodynamic matching

Use the forward coefficients and collider efficiencies in the
[Cantera v3.2.0 h2o2.yaml](https://github.com/Cantera/cantera/blob/v3.2.0/data/h2o2.yaml)
hydrogen–oxygen submechanism of GRI-Mech 3.0. The SHA-256 of the inspected source
is `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
The result provenance lists every retained coefficient and its source reaction
number. Convert cm/mol/s to m/molecule/s: multiply a reaction of order m by
`(1e-6/N_A)^(m-1)`. Activation energy is cal/mol and converts with 4.184 J/cal.

Retain reactions 1, 2, 3, 11–15, 21 and 23 among H, O, H2, OH, O2 and H2O.
The three-body water-formation rate is `2.2e22*T^-2 cm6 mol^-2 s^-1`, about
`6.066e-38*T^-2 m6 molecule^-2 s^-1`, close to the older water-study OH-partner
proxy. Its H2 and H2O collider efficiencies are 0.73 and 3.65. The separate
H2/H2O collider channels for H2 formation must not be double-counted in the
generic third body. Rates for `2 H` use the mechanism's event convention
`q=k[H]^2[M]`, with two H consumed per event, not an extra factor of 1/2.

Reverse coefficients are **not** taken from GRI's NASA7 thermodynamics. Match
them by detailed balance to this study's existing molecular partition functions,
so every retained reaction equilibrates at the same composition and on the
same chemical energy zero as the flow. Reject a local state that fails these
mass-action identities. This hybrid rate/thermodynamic diagnostic must be named
as such; it is not an unchanged GRI-Mech simulation.

Ions and electrons are fixed spectators. Their reaction channels and collider
effects are absent. HO2/H2O2 pathways are omitted because these species are not
in the current water EOS. Adding them requires matching their thermodynamics
and energy bookkeeping as well as rates. The dense residual is composition
independent at fixed mass, density and temperature; it does not change the
matched ideal-mixture equilibrium constants, nor supply dense-fluid rate
corrections.

## Local linear relaxation, not one-way equilibrium heat release

At equilibrium each reversible reaction has equal forward and reverse flux q.
Let n be the six neutral number densities and nu the product-minus-reactant
stoichiometric vector. The negative linearized kinetic Jacobian, transformed
by sqrt(n), is the symmetric positive semidefinite matrix

```
B = sum_r q_r (nu_r/sqrt(n)) (nu_r/sqrt(n))^T.
```

Two elemental conservation modes are null. Project onto four independent
stoichiometric directions, diagonalize there, and weight each decaying mode
by its projection on chemical ground-state energy `e0*sqrt(n)`. The diagnostic
time is `tau_E = sum(w_i/lambda_i)/sum(w_i)`, the integrated normalized local
bond-energy relaxation time. It is independent of additive conserved elemental
energy gauges. This is an **isothermal, fixed-density, fixed-ion** relaxation;
reaction/temperature feedback and moving equilibrium are not integrated.
Equilibrium one-way fluxes do not imply net chemical heating or useful work.

Nonpositive eigenvalues or a reduced spectral condition number above 1e12
are reported unresolved, never clipped into a fast or slow verdict. Tests
check detailed balance, elemental conservation, event/unit conventions,
finite-difference reversible Jacobians, elemental energy-gauge invariance,
and the exact inverse scaling of relaxation time when all rates are scaled.

## Flow clocks and weighted applicability

Audit 100–250 microseconds, retaining exact boundaries at 153.5 and 193
microseconds. In each fixed-mass parcel use positive consecutive bond-store
decreases as weights. Repeated decreases after reheating count again: these
are neither unique reaction energy, decreases from one global peak, nor causal
impulse weights. Evaluate rates and the actual velocity-gradient expansion
clock at the interval's right endpoint. Thin outputs while retaining all
window boundaries to check this assignment.

Report `Da_E = 1/(expansion_rate*tau_E)` only when a measured positive expansion
rate exists. Also report the sampled mean-held-store / return-rate clock, which
can capture rapid chemistry changes not described by density alone. The
tracking screen uses the smaller of those two clocks. Compression gets no
invented expansion clock. This sampled clock is not a resolved continuous
composition forcing rate or a nonlinear kinetics integration.

All positive decreases stay in the summary denominator, including unresolved
states, missing expansion clocks and tiny interval changes. Rate evaluation
is skipped below `1e-10 * parcel full-atomization energy`, with that weight
reported separately. Cases below `1e-6 * case full-atomization energy` are
flagged negligible; trace numerical decreases in the unsprayed control must
not drive a physical freeze-out conclusion. Quantiles of clocks refer only
to the explicitly resolved subset; weight coverage is reported beside them.

## What this can and cannot validate

The source file's six neutral-species NASA7 fits end at **3500 K**. This is a
thermodynamic-fit limit, **not** a validated temperature range for every forward
rate. Our study's thermodynamics is used instead, but the range marker shows
why an off-the-shelf h2o2.yaml reactor is not a drop-in validation.

The retained association reactions have three-body rate laws with **no
high-pressure falloff limit** in that source. The retained network also lacks
a demonstrated dense-fluid rate/activity correction and uncertainty envelope
at the flow's temperatures and pressures. Pressure bins of 100 MPa and 1 GPa,
the difference between gas and flow pressure, and ionized-nuclei fractions
are applicability diagnostics, not certified rate limits.
The file's optional Redlich–Kwong phase does not establish that its kinetic
rate laws were validated at these pressures or match our dense-fluid closure.

Even a large rate margin therefore remains conditional. A uniform slowdown
by 1000 is a sensitivity calculation (Da divided by 1000), not a credible
uncertainty interval for missing physics. All summaries retain
`kinetics_validated=false`. A finite-rate nozzle prediction still requires
an evaluated reaction model in this regime, energy-consistent coupling,
and the actual injected-flow/finite-plate geometry and losses.
