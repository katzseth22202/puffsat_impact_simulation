# Water-plate flow reference: cold source energy and a matched dense-fluid continuation

Scope: the **water-injected overtake plate study**, not a replacement for the
original `f(v)` or head-on tamper EOS. This decision enables a diagnostic planar
prepositioned-water reference, not a validated injector or warm-dense-water model.

## Source and hydrodynamics

Use the existing staggered Lagrangian kernel (ADR-0022), a rigid adiabatic wall,
and a vacuum outer boundary. Start with stationary cold water against the wall
and the incoming pulse immediately beyond it. There is no interpenetration,
thermal diffusion across the material contact, droplet model, feed boundary,
lateral escape, radiation or wall loss. Consequently this is **not** the study's
Stage 1 immediate-local-mixing implementation. The initial water bulk density
sets its extent, not a prescribed shocked or mixed-gas density.

The cold source carries stored-liquid internal energy from the study's existing
CoolProp/gas bridge. Flashing, heating and dissociation are collision-paid. Add a
common 4 MJ/kg storage offset to every phase solely for the positive/log-energy
table interface; subtract it in physical energy reports. A mass-weighted nodal
velocity projection deposits its kinetic-energy deficit locally as heat. Wall
pinning has a separately recorded initial impulse. Both effects shrink with
spatial refinement and must not be hidden as physical mixing or wall chemistry.

Integrate the wall force using the actual two Verlet force evaluations, including
artificial viscosity, and reconcile with nodal momentum. Record energy drift
and table-domain exits; invalid runs cannot provide a matched gain. A completed
finite time window does not mean settled impulse or completed exhaust ejection.

Spatial refinement alone left about 1.5% shock energy drift in the unsprayed
control at the original CFL 0.4. The reference therefore uses a quarter of the
usual timestep (CFL 0.1), with a further CFL 0.05 check; existing callers are
unchanged. A new ideal-gas impact test verifies that shock energy drift falls
with smaller timesteps. The matched-gain diagnostic requires <=0.5% energy drift
in each simulated collision control, as well as momentum closure and matching
geometry/time windows. This tolerance does not establish table convergence.

The cold-only free surface enters an unmodelled sub-freezing phase on refined
grids. Do not clamp that energy or claim a complete spatial history. For a
uniform stationary layer the wall remains at its initial pressure until the
vacuum rarefaction head arrives at `t_head = L/c_initial`. The 1 m anchor has
`t_head = 0.14665 s`, much longer than the 0.001 s comparison. Its cold-only wall
control is therefore **analytic**, `J = p_initial*A*t = 91.94 N s`, with no
invented numerical energy/momentum residuals or downstream cell history. Coarse
numerical controls, while still within the table, reproduce this wall impulse.
The runner allows this analytic control only before the head-arrival time;
longer windows still require an adequate cold-phase EOS.

## Thermodynamic continuation

Below and at `T0 = 1000 K`, use CoolProp/IAPWS-95 on the study's common energy
zero. The existing ideal equilibrium-water partition functions define the gas
Helmholtz potential, including OH/H2/O2 and the oxygen ion ladder. Align the
additive entropy constant at the same 400 K, 0.01 kg/m3 reference as the energy
bridge. This liquid/vapor model omits high-pressure ice; above 1 GPa the IAPWS
values themselves are extrapolations, including some cold table corners. Stability
does not validate those states. At each density define residual quantities `Delta = real fluid - gas`
at T0. Above T0 choose the explicit provisional model

```
cv_res(rho,T) = Delta_cv(rho) * (T0/T)^n,      n = 2 nominal
e_res = Delta_e + Delta_cv * T0/(n-1) * [1-(T0/T)^(n-1)]
s_res = Delta_s + Delta_cv/n * [1-(T0/T)^n]
f_res = e_res - T*s_res
p_res = rho^2 * (d f_res/d rho)_T
```

The pressure is evaluated using reference pressure derivatives and the identity
`rho^2*d(Delta_cv)/d rho = -T0*d^2(Delta_p)/dT^2`. The total potential is
`f_gas + f_res`. Pressure, energy, entropy and heat capacity match at T0;
first derivatives of pressure and energy, and acoustic speed, are continuous.
The cold EOS is never evaluated above T0. This continuation retains a dense-fluid
pressure contribution without extrapolating liquid water straight into an ideal
gas. It is **not** a calibrated high-temperature interaction/composition model.
Changing n from 2 to 3 is a limited closure-sensitivity check, not a credible
uncertainty bound on warm-dense water.

Compute acoustics from the same first-law derivatives:
`c^2 = p_rho + p_T*(p/rho^2-e_rho)/e_T`. Reject nonpositive heat capacity,
nonpositive acoustic speed squared, nonpositive fields, nonmonotone tabulated
energy or a Maxwell-identity residual above 0.1%. The source model passes a
dense temperature/density transition grid, not just sparse endpoints. Table
interpolation is still the existing log-bilinear field interpolation; its
thermodynamic error and storage-offset dependence require table refinement
before precision performance claims. Source consistency does not make separate
interpolated fields an exactly consistent potential.

## Rejected attempts and interpretation

A temperature-weighted pressure/energy blend violated the Maxwell identity.
Blending Helmholtz potentials fixed that identity, but the provisional
temperature-independent dense residual gave negative heat capacity near
303.822 kg/m3 and 1624.95 K. Widening its transition shifted, rather than removed,
the unstable band. Matching only the residual entropy/energy tangent still
failed near 2000 kg/m3. None of these rejected closures is a flow result.

The current heat-capacity construction removes that uncontrolled switching term;
it does not establish the physical accuracy of the dense/hot model. Distinguish
cold compacted water from hot dense gas, report the mass in the provisional
dense/hot region, and retain the unsprayed and cold-water-only controls.

Chemical stores decompose the ideal equilibrium-gas part of the potential, not
the real-fluid residual. Per-parcel sampled peak decreases and actual expansion
clocks can establish temporal overlap with wall loading **within this model**.
They do not isolate a causal chemistry contribution to impulse. Finite-rate
chemistry and a controlled frozen-composition comparison remain separate work.
