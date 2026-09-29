# The near-term chamber: mixed feeds, a JANAF partition correction, and three carbon branches

Status: **accepted**, 2026-09-29. Answers asks **A1-A4** (A5 and A6 lightly) of
`Balloon-Pulse-Propulsion` @ `fe04471`, `todos/impact_sim_asks_methane_7000K.md`
(paper `sec:methane_7000_near_term`). Extends
[ADR-0050](0050-walled-thermal-nozzle-methane-chamber.md) and
[ADR-0051](0051-walled-nozzle-grid-geometry-and-the-load-case.md); every scope
statement there still holds, including that this study shares only the EOS
machinery with `f(v)` and the tamper study.

## Context

The paper proposes methane at 7000 K on a pitch-lined steel wall as the near-term
walled chamber: a 2.5 kg polyethylene rod at 75 km/s into a 20 m^3 sphere, with a
4.4 kg polyethylene plug that leaves with the exhaust. It prices this with
efficiencies borrowed from ADR-0016's 10 kK, 200 m^3, 25 kg cells. The ask wants
them solved on the paper's own geometry, with the carbon question bracketed.

Four things were not possible on the existing machinery. `eos_methane` hard-codes
H:C = 4. The Isp ledger has no carried mass besides the slug. Nothing condenses
carbon. And, found while answering A1, the harmonic partition functions are
badly wrong for C3 and C2 in exactly the 5-7 kK band this chamber sits in.

## Decision 1: a feed is an H:C ratio and a reference energy, and nothing else

`eos_methane.Feed` carries the H:C ratio and the cold feed's atomisation energy
per carbon nucleus, which is the zero of `e`. `mixture_feed` builds one from
`(mass, component)` pairs (CH4 gas, H2 gas, solid polyethylene). Once the feed is
hot, the equilibrium cannot know which molecules it came as, so a methane charge
with a polyethylene rod and plug is the same solve at `r = 3.77`, and a hydrogen
charge with that polyethylene is `r ~ 121`.

`METHANE` is the default everywhere. Its residual keeps the exact
CH4-divided-out form, and every earlier result is **bit-identical** (checked by
direct comparison, and by the existing tests). A non-4:1 feed takes the plain
total-ratio row, an off-stoichiometric molecular initial guess (CH4 + H2 or
CH4 + C2H2), and a warm start from the previous solution. Pure methane never
takes that path, so its results cannot depend on call history.

The polyethylene heat of formation, -27.0 +- 2 kJ/mol per CH2, comes from its
heat of combustion. The WebBook carries no value. It is a 298 K number against
the module's 0 K set, which is a <0.3% effect on a component that is 10% of the
charge.

## Decision 2: JANAF partition functions are opt-in, and only the partition part is imported

Against NIST-JANAF, with the heat-of-formation difference removed, the harmonic
C3 is over-bound by `e^3.45` (~30x) at 6000 K and C2 under-bound by `e^2.2`. At
7000 K that moves the methane chamber's charge from 72.5 to 68.8 kg, and C3's
share of the store from 16% to 2%. It is not a small correction there. At 10 kK,
the flown ADR-0016 point, it moves the charge by 1%.

`Feed.thermo = "janaf"` shifts each molecule's atomisation constant by a fitted
`delta(T)` (`a + b ln T + c T + d/T` over JANAF's 1000-6000 K; residual <0.006),
and adds `k T^2 d(delta)/dT` to that molecule's energy. That is the same shift
written on `ln z`, so `e` and `K` remain one thermodynamics, which a van't Hoff
test checks. **Only the partition-function part is imported**: JANAF's C2H heat
of formation is ~92 kJ/mol below the modern value this module carries, and
taking JANAF whole would import that error.

It is opt-in rather than a default change because every published W-result is on
the harmonic set. Changing it under them would silently move numbers the paper
already cites. The near-term study runs JANAF as primary and harmonic as the
sensitivity. Above 6000 K the fit is extrapolated: 1000 K at the near-term
chamber, and 4000 K at 10 kK.

Tables are transcribed verbatim with their JANAF ids into `janaf_carbon.py`.

## Decision 3: three carbon branches, each a single constrained isentrope

- **equilibrium**: the existing treatment.
- **frozen_carbon**: every carbon-bearing species frozen per kilogram at the
  *peak chamber* composition, for the whole pulse. Hydrogen alone re-equilibrates
  (H + H + M is the one measured channel, and the freeze clock says it keeps up).
  Freezing at peak rather than at a throat means chamber blowdown and nozzle lie
  on one isentrope, and it is the conservative floor the reactor comparison
  needs.
- **graphite_ceiling**: equilibrium plus condensed carbon wherever monatomic C
  exceeds its saturation density over graphite. Graphite is NASA-Glenn C(gr)
  against the EOS's own C(g). It reproduces JANAF's `log Kf` for C(g) to 0.01 dex
  at 3000-5000 K. This is **a thermodynamic ceiling**: soot at no kinetic cost, in
  velocity and thermal equilibrium with the gas. It is not a nucleation result,
  and N10 item 5 stays open.

Carbon association still has no rate (`rates.CARBON_CHANNEL_NOTE`, and the
provenance rule forbids inventing one), so the branches bracket it rather than
decide it.

## Decision 4: eta is reported on the blowdown, beside the companion's peak convention

The companion's conversion `v^2 / 2u` expands the *peak* state from its
stagnation enthalpy `h0 = u + p/rho`, which credits every kilogram with the flow
work `p/rho`. A closed vessel emptying adiabatically hands out exactly `U0`
(`integral h d rho = rho0 u0` along the isentrope, which is tested). So
`eta_blowdown` integrates each parcel's exit over the chamber's own emptying
history. Isp uses the impulse-averaged speed, which sits below the paper's
`w sqrt(eta/(k+1+P))` by Jensen's inequality. Both conventions are written, and
the paper's `eta` means the blowdown one.

`chamber.IspLedger` gains `plug_ratio` (default 0, bit-identical), so the
head-on sign is still written in one place: `[(1+k+P) v - w] / ((k+P) g0)`.

## Decision 5: the wall load is integrated over the real blowdown

A4 runs the chamber's state against time on the real-gas choked flux, so there
is no `gamma = 1.2` coefficient. Radiation is bracketed:

- the lower edge is the **H- bound-free floor**, in-band (0.25-1.55 um, ~90% of
  a 7000 K blackbody), with the optical depth taken over the mean beam length;
- the upper edge is the opaque blackbody.

Convection is Bartz at the liner (`wall.bartz_coefficient`, split out of
`bartz_flux` unchanged) with its bracket spread to the corners. The paper's own
cost model is reused for A6 (leak `2f`, wall `r/2`), with the gas state and wall
load re-solved at each volume rather than scaled.

## Consequences

- The shared EOS gains a feed and an optional partition correction. Neither
  changes a methane result unless asked for.
- The ADR-0016 10 kK cells should eventually be re-run on JANAF. The effect there
  is ~1% on the charge, but C2 and C3 matter more on the cooler stations of the
  expansion.
- `make walled-nozzle-near-term` produces
  `data/results/walled_nozzle/near_term/*.csv`. The results are written up in
  `docs/walled_nozzle_near_term.md`.

## Amendment (2026-09-29): area ratios to 1000, and eta net of hydrogen freezing

**Area ratios extend to 1000.** `AREA_RATIOS` is 4, 14, 30, 100, 300, 500 and 1000. The ask
named 4-100. The higher values are in family for vacuum upper stages: the RL10B-2 flies 285:1 on
a carbon-carbon extension, and high-expansion studies run 500-1000:1. The ratio is conventional;
the absolute size is not. On the paper's 0.149-0.198 m^2 throats, 1000:1 is a 14-16 m exit. The
throat cannot shrink to compensate: the entrance leak (share = port area / throat area) and the
4 Hz emptying time both pin it. So the results report exit diameter beside the ratio. The default
isentrope runs to 1e-9 of the chamber density, so late-blowdown parcels still resolve their exit
at 1000:1.

**`eta_net` is the eta to quote.** `FreezeCheck.eta_net` is the equilibrium `eta_blowdown` less
the atomic-hydrogen energy still held at the Bray sudden-freeze station (first `Da < 1` on
H + H + M) and not returned by the exit. This is an upper-bound debit: after freezing, the flow
still converts its thermal share. It matters because recombination freezes inside these small
nozzles, at A/A* ~ 100 for methane and ~ 23 for hydrogen. Past that point, equilibrium eta
overstates what a real nozzle returns.
