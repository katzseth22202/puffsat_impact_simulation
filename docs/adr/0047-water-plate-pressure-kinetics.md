# Water-association pressure-source coverage and channel-removal sensitivity

## Scope

Continue ADR-0046 on its exact saved equilibrium histories. Change only the
direct H + OH <=> H2O channel's forward and reverse flux together, keeping
the other nine reactions, equilibrium composition, chemical energy reference,
isothermal energy susceptibility and sampled flow clocks fixed. No new
hydrodynamic trajectory or wall impulse is calculated.

## Source investigation

The inspected [Burke–Dryer 2011 mechanism](https://github.com/Pele-Suite/PelePhysics/blob/ab2c54a4183cc97ff6a390b1ddc0da5367535e1f/Mechanisms/BurkeDryer/mechanism.inp)
(published as Burke et al.'s
*Comprehensive H2/O2 Kinetic Model for High-Pressure Combustion*) still uses
uncapped three-body water dissociation. Its water rate cites Srinivasan and
Michael (2006), with generic and explicit steam-collider channels. The
[PelePhysics FFCM1_Red mechanism](https://github.com/Pele-Suite/PelePhysics/blob/ab2c54a4183cc97ff6a390b1ddc0da5367535e1f/Mechanisms/FFCM1_Red/mechanism.inp)
likewise has those two uncapped water channels.
Having pressure falloff for HO2 does not establish water falloff. These are
mechanism-file observations, not inspection of the cited 2006 article or
proof that no applicable evaluation exists.

A pressure-dependent **association** law is available in the Cantera
example-data implementation of Singal, Lee, Lei, Speth and Burke (2024),
*Implementation of New Mixture Rules Has a Substantial Impact on Combustion
Predictions for H2 and NH3*. It modifies Alzueta et al. (2023), and cites
Jasper (2022), *Predicting third-body collision efficiencies for water and
other polyatomic baths*, for collision efficiencies.

- [Pinned mechanism file](https://github.com/Cantera/cantera-example-data/blob/1a5d27e508a38b1791543e9fded80ffd5c5b8d75/ammonia-CO-H2-Alzueta-2023.yaml#L966):
  SHA-256 `9eed1a29952ec41d154185798ecef7127a6329501322ee474be0b638011733f7`.
- [Author-group collider database](https://github.com/TheBurkeLab/LMRRfactory/blob/ab812a1fe29038fd48a16ba402bc34c5c3ece1b3/src/LMRRfactory/thirdbodydatabase.yaml#L219):
  SHA-256 `4ff0b6af69f8a56390732dd55d2401983ec9392c35f1bf23524f7826c574ca75`.
  Water-system relative-to-Ar efficiencies are supplied at **300, 1000,
  2000 K**, for N2 (2.51, 1.62, 1.20) and H2O (8.97, 8.55, 8.75).
- [Cantera 3.2.0 mixture-rule implementation](https://github.com/Cantera/cantera/blob/v3.2.0/src/kinetics/LinearBurkeRate.cpp).

The pressure grid runs from 1e-4 to 1e4 atm (upper node **1.01325 GPa**).
This is a finite interpolation table, **not an explicit high-pressure limit**.
Neither these supplied pressure nodes nor the cooler collision-efficiency
data establish an uncertainty envelope at this study's 5–9 kK loading
states. Publisher/evaluation access was incomplete; do not label the
primary articles as fully inspected or their kinetics as validated here.

## Conditional pressure-source variant

Extract only the water association PLOG curve and H2O collision efficiency.
The source specifies m3/kmol/s for this reaction, overriding its global
cm/mol units. Divide by `1000*N_A` for m3/molecule/s. Interpolate log(k)
linearly in log(P); **do not multiply by another collider density**.

N2 is the reference collider. Among our six neutrals only H2O has a
specific source efficiency: `epsilon_H2O=0.104529*T^0.550787*exp(232.675/Rcal/T)`.
All five others inherit its default efficiency of one. This is an explicit
source assumption, not a hot-radical collider evaluation. Ar is absent.
With no collider-specific alternative rate curve, the linear-Burke mixture
sum reduces exactly to the reference curve at

```
P_eff = k_B T * (n_H + n_O + n_H2 + n_OH + n_O2 + epsilon_H2O n_H2O).
```

This is **neutral ideal collider pressure**, not the dense EOS pressure;
ions/electrons are not silently made colliders. Reject evaluation outside
the supplied effective-pressure table, rather than Cantera's default PLOG
endpoint clamping. Retain that weight as no supplied rate. Temperatures
above the supplied collider-data points remain explicitly extrapolated even
when effective pressure is tabulated.

Use the extracted forward association coefficient with the **study's**
equilibrium constant for the reverse, not the source NASA thermodynamics.
Replace reaction 15's equilibrium flux with that value; no equilibrium or
energy jump is introduced. Derivatives of the composition-dependent forward
coefficient cancel from the equilibrium Jacobian because they multiply zero
net progress; the same detailed-balanced outer-product construction applies.
This is a single-channel hybrid variant, **not**
the full Alzueta/Singal mechanism. In particular its separate common
H2O + H2O dissociation reaction and other water-producing channels are not
added on top of the extracted collider-dependent rate.

## Direct-channel removal as a local slow endpoint

Also set reaction 15's forward and reverse flux to **zero**. Water can still
form through H2 + OH <=> H + H2O and 2 OH <=> O + H2O. This does not remove
all molecular chemistry or all other three-body associations.

In the four-dimensional reactive subspace of ADR-0046, let
`B(s)=B_other+s*v*v^T`, where s is any nonnegative direct-channel equilibrium
flux. For the fixed projected bond-energy vector b,

```
tau_E(s) = b^T B(s)^(-1) b / (b^T b)
d tau_E / ds = -(v^T B(s)^(-1) b)^2 / (b^T b) <= 0.
```

The other nine reactions retain full reactive rank. Thus removal gives a
**slow endpoint for this fixed local linear network**, encompassing any
nonnegative replacement of this one channel. It is not an uncertainty
bound on all water chemistry, a nonlinear/adiabatic theorem, or a bound on
impulse. Other retained association rates still lack high-pressure limits.

## Verification and evidence contract

Preserve the original Step 7 artifact and default calculations. A separate
`make water-plate-pressure-kinetics-audit` target reuses the parent gates,
exact windows, fixed-parcel positive store-decrease weights and two spatial
grids/output strides. Every summary keeps `kinetics_validated=false`.
Untabulated pressure, tiny changes, missing clocks and unresolved spectra
remain in denominators. Quantiles refer only to explicitly covered subsets.

Tests check pressure nodes, log interpolation, molar units, out-of-table
rejection, mixture pressure, agreement with independent Cantera 3.2.0
forward-rate evaluations, monotonic relaxation versus direct-channel
strength, and agreement of the zero-channel eigenanalysis with an independent
element-conservation-nullspace linear solve. The flux change scales both
directions together, preserving the original detailed-balance construction.
