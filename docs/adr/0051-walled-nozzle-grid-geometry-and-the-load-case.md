# The N16 grid runs on a cone, and the N9 load case is a snowplow with a clock

Status: **accepted**, 2026-09-10. Answers asks **N16** and **N9 items 1-7** of
`Balloon-Pulse-Propulsion` @ `dd6ad2a`, `docs/nozzle_asks_for_impact_sim.md`.
Extends [ADR-0050](0050-walled-thermal-nozzle-methane-chamber.md); every scope
statement there still holds.

## Context

N16 replaces N14 and N15 and asks for one surface over chamber temperature,
chamber volume and throat area, for methane **and** water, plus the density
dependence of the held store, which it calls "the most valuable thing this ask
can return". N9 items 1-7 ask what the arriving front does to the wall, and the
ask itself names the gate: below about 170 m^3 the front may leave the column
before its cone reaches the liner, which would close items 1-4 rather than
answer them.

Two things about the geometry moved on the paper side between ADR-0050 and this
ask, and they are what force new decisions rather than a wider sweep of the old
machinery: chamber volumes fell to 50-400 m^3, and throat areas fell to
0.05-0.2 m^2 -- area ratios of 140 to 600 against the flown 4.

## Decision 1: the throat is not a third dimension of the grid

A chamber at `(rho_c, T_c)` has exactly one isentrope. The throat area decides
only **where along it** the gas stops, because the exit area is the 3 m bore in
every case. So `surface.py` runs one deep expansion per `(fluid, T_c, V)` out to
the widest area ratio in the grid and reads every throat off it.

That is 32 expansions rather than 256, and -- the reason it is a decision and not
an optimisation -- it makes the throat comparison exactly like-for-like: the same
gas, the same adiabat, one cut point moved. `test_walled_nozzle_surface.py` pins a
point read off the deep isentrope against an expansion run directly to that ratio.

## Decision 2: the nozzle clock is rebuilt on a cone, not inherited from a length

`expansion.cooling_history` opens the area **linearly in `A`** over a length the
caller supplies. That was right for `freeze.py`, where the geometry was a field
profile of fixed length. It is wrong across two decades of area ratio: at fixed
length a larger area ratio means a faster opening, so the Damkoehler number would
fall with area ratio for a reason that is pure bookkeeping.

So the grid puts a **conical nozzle at 15 degrees** from the throat radius out to
the 3 m bore and integrates `dt = dx/u` along it:

    x(A) = (sqrt(A/pi) - r*) / tan(theta),   L = (R_bore - r*) / tan(theta)

Consequences, both reported rather than hidden: across the throats the ask
recommends (0.5 m^2 and below) the length is nearly constant at 9.7-10.7 m, so
those columns differ only in where the cone starts; the flown 7 m^2 throat is a
genuinely *shorter* nozzle at 5.6 m, because a 1.5 m throat radius is already half
the bore. `Da` scales linearly in the length, so the half-angle is a stated
sensitivity and not a hidden one.

## Decision 3: the equilibrium branch is reported with a freeze cap beside it

The deep throats this ask wants **do** freeze -- which is new, and is the opposite
of what W5 found inside the flown bore. The equilibrium conversion past a freeze
station is not a physical number, so every grid row carries a second column,
`conversion_capped`, which charges the equilibrium branch for the chemical energy
it releases *downstream* of the freeze station.

It is a **lower bound**, and deliberately: the frozen gas is also colder, so part
of that energy would not have reached the jet as thermal enthalpy either. The
truth is between the two columns, and the two coincide wherever nothing freezes.
Re-solving a genuinely frozen methane expansion needs a
`pressure_energy_frozen` that `eos_methane` does not have (`eos_water` does), and
that is recorded as deferred rather than approximated.

## Decision 4: the mixture N16 names is assembled, not solved, and says so

N16 asks to price "90% water with 10% liquid hydrogen" first. There is no joint
H/O equilibrium EOS here, so `surface.mixture` composes `eos_water` and
`hydrogen` as **two independent subsystems** at a common temperature and their own
partial densities. That is exact for the thermal terms and each component's
internal chemistry, and wrong about exactly one thing: the free hydrogen does not
shift water's dissociation equilibrium. Le Chatelier puts the real mixture's held
store *above* this, so the mixture rung is an optimistic edge and carries
`estimated=True` on every row -- the same discipline `propellants.ammonia_rung`
already uses.

## Decision 5: the front's wall load is priced on the layer's clock, not its temperature

`wall.py` integrates `coupling.snowplow`'s system with a clock and a general EOS,
pinned in the tests against `front.integrate` on the magnetic geometry so the two
models cannot drift apart.

The load case is then split three ways and only one of them is a real delivery:

- **Radiative** is `held x (1 - exp(-dt/tau_rad))`, with `tau_rad` the binding
  branch of `expansion.radiation_check`, capped at `sigma T^4 dt`. **That cap
  alone settles the item**, which is why the verdict does not depend on the
  opacity model.
- **Convective** is reported as an *incident* stagnation-enthalpy flux times a
  Stanton bracket of 0.001-0.01, because gas arriving at a wall does not deposit
  its energy -- the boundary layer decides what crosses. This is the crudest
  number in the module and it is the one that binds, which is exactly what the
  ask says of its own item 7.
- **The ceiling** -- everything the shocked column holds -- is printed beside both
  so a reader can see how far from it the answers are.

## Decision 6: the opacity is the free-free floor, bracketed, and the floor is the conservative choice

The shocked layer is a partly ionised C/H plasma at tens of thousands of kelvin.
`plasma_opacity` computes Thomson plus Kramers free-free on the **solved**
composition and nothing else. Bound-free and line opacity are omitted, and every
omitted term makes the layer *thicker*; a thicker layer in the diffusion regime
delivers *less*. So the floor is conservative for a question about how hard the
wall is hit, and `OPACITY_BRACKET` carries 0.1x-10x to show the verdict surviving
it. This is the same posture `opacity_bracket.py` takes for the magnetic nozzle
and the opposite of quoting one number from a table this study does not have.

## Decision 7: throat carbon is decided at the wall-adjacent state, not the free stream

Item 5 asks whether the throat self-cleans. The first attempt compared the **free stream's**
monatomic-carbon partial pressure against graphite's vapour pressure at the **wall's**
temperature. That mixes two states and is wrong: pressure is constant across a boundary layer and
the gas at a no-slip wall sits at the wall temperature, so the state that decides whether carbon
sticks is a **third** one -- the free stream's pressure at the wall's temperature, denser than
either end.

`wall_gas_density` builds that state and `equilibrium_saturation` evaluates `S = p_C / p_vap` on
it. The verdict reverses: `S` = 0.03 at graphite's working temperature rather than the
supersaturation the mixed-state comparison reported. The mechanism is acetylene, which is the same
sink W9 found from the recovery side, and it means the C3 partition function that is
`eos_methane`'s worst weakness does not reach this answer -- C3 never carries more than 8% of the
carbon at any wall temperature in range.

**Both edges of the fork are kept**, because the answer depends entirely on which one the boundary
layer sits in and they differ by six orders of magnitude: `deposition_threshold` is the frozen
edge (renamed and re-documented rather than deleted, since it is what the earlier reading
compared), `equilibrium_saturation` is the other.

**The fork is closed on a bound, not on a rate.** `boundary_layer` gives the layer its thickness
(`k/h`, which is free of Bartz's `c_p` uncertainty because `h` carries the same factor) and a
diffusion clock, and `decades_below_gas_kinetic` reports how far under a **collision frequency**
the carbon chemistry could run and still finish inside it. That is 5.0-5.3 decades. No rate
coefficient is asserted -- `rates.py`'s provenance rule forbids inventing one -- and the number is
stated as what would have to be true to flip the verdict.

`hertz_knudsen_thickness` is kept as the free-flight ceiling on a plating rate. It is deliberately
useless as a prediction (centimetres per pulse) and is retained only to show that the frozen edge,
if it applied, would not be rate-limited either.

## Decision 8: the contraction is priced from the energy budget, not from a flux

Decision 5's split answers the *side wall*. Below the gate volume the front never reaches it, and
the load moves to the convergent section, which at these throat areas is 96-99% of the bore --
effectively a flat end wall taking a normal-incidence strike.

`end_wall_strike` prices it as two terms with different scalings, and keeping them apart is the
decision:

- the **arrival**, the front's own remaining kinetic energy over the patch its cone has opened to.
  It grows as the chamber shortens, because a shorter column decelerates the front less.
- the **bulk**, the whole pulse energy over the contraction area. It is **volume-independent** by
  construction -- the same 70.3 GJ turns the same corner -- so chamber volume is not a lever on it
  at all.

The bulk term is an **energy-budget ceiling**, reached the way `strike`'s `fluence_ceiling` is:
everything the impactor brought must cross that surface on its way out, so dividing by the area
gives what would land if none of it turned. What crosses into the wall is that times the same
Stanton bracket, with the same caveat. A flux-based treatment would need a boundary-layer solve on
a curved converging surface, which is recorded as deferred rather than approximated.

The consequence is a recommendation the grid alone would not have produced: shrinking the chamber
has a **floor** near 150 m^3, below which the arrival load climbs faster than anything else
improves.

## Decision 9: effective Isp carries both mass-ledger corrections, and it lives in one dataclass

Two corrections separate the real Isp `u_e/g0` from the number a mission plan should use, and
**they run in opposite directions**:

- a **credit** of `(1+k)/k`, for the impactor mass the vehicle never lifted. This one was already
  in the study.
- a **debit** of `w/(k g0)`, for the momentum that same impactor brings in **head-on**. It arrives
  against the ship's motion, so the exhaust must cancel `m_p w` before any of it is thrust. This
  one was missing everywhere, and adding it is the decision.

`isp_effective` now means both, which is the companion's own head-on form (`templateArxiv.tex`:
`I/(mw) = eta_jet sqrt(1+k) - 1` head-on, against `+ 1` on an overtake):

    Isp_effective = [ (1+k) u_e - w ] / (k g0) = w (eta_jet sqrt(1+k) - 1) / (k g0)

**The debit is the larger correction, by exactly `w/u_e`.** Both go as `1/k`, so the credit does
not dominate by being "the mass one"; since `u_e` runs four to nine times below `w` here, the
corrected figure lands **below** the real Isp rather than above it. That reverses the sign of a
whole column's worth of intuition: the old inequality `isp_effective > isp_true` was asserted in a
test, and the corrected one asserts the opposite.

**The naming is deliberate and it changes the meaning of a published column.** The alternative --
keeping `isp_effective` credit-only and adding `isp_net` beside it -- was implemented first and
rejected: it leaves the study's most quotable name attached to a quantity no reader should quote,
and every downstream table then has to remember which of two nearly-identical names it wanted. So
`isp_effective` is the corrected figure of merit and the credit-only intermediate is renamed
`isp_carried_gross`, explicitly "not a figure of merit". The cost is that `isp_effective_s` in the
CSVs changes meaning; it is paid down by shipping `isp_carried_gross_s` alongside, so the header
itself tells a reader which vintage of the file they hold, and by saying so in both answer
documents.

**The credit-only column is kept, not deleted.** It is what W5 and W8 published, `isp_effective`
is built from it, and the `Isp ~ 1/sqrt(mean atomised particle mass)` scaling law is a property of
*that* convention -- the debit compresses the ladder precisely because it punishes the small-`k`
fluids the law rewards, so asserting the law on the corrected column would be asserting the debit
away.

**The algebra lives in `chamber.IspLedger` and nowhere else.** The first implementation put it in
`surface.py` for the N16 grid and left `propellants.py` -- the W8 ladder the paper lifts directly
-- on the uncorrected convention. That is the failure mode the sign is most exposed to, so the
ledger is now a single dataclass upstream of both, and `chamber.isp_scaling` was corrected to
`(sqrt(1+k) - 1)/k` at the same time. Three things are pinned by test rather than assumed:

- **`eta_jet = sqrt(conversion)`**, since `(1+k) u = w^2/2` makes `w^2/(1+k) = 2u`. A conversion
  fraction quoted *as* an `eta_jet` would be wrong by a square root.
- **the thrust floor is `1/sqrt(1+k)`** -- not a separate assertion but the root of
  `isp_effective`; below it the nozzle pushes the ship backwards. Every cell clears it by 3.5-5.6x.
- **agreement with the tamper study.** At `eta_jet = 1` the form collapses to `sqrt(1+k) - 1`,
  which is `tamper.ledger.beta_ideal`, derived independently from Cauchy-Schwarz under a different
  PRD. A test asserts the two match; if they ever diverge, one study has the sign wrong.

**The sign is the hazard, not the algebra.** The water-plate study in this same repository records
its incoming momentum as a **credit**, because that one is an overtake. Inheriting that `+1` here
would flip a ~30% penalty into a ~30% bonus of the same size.
`test_the_head_on_burn_is_a_debit_and_never_a_credit` asserts the corrected figure is below both
the credit-only one and the real one, for every cell, so the mistake cannot land silently.

## Decision 10: a narrow throat is priced from the blowdown dwell, not from the peak flux

W6 recommended narrowing the throat and W19 found a wall past `A/A*` ~ 140, but neither priced
the cost, and the paper's only stated cost is the pulse stretching from 8 ms to 28 ms. Crossing
the two needed a decision about *which* thermal quantity does the pricing, and the intuitive
answer is wrong.

Bartz gives the **throat** flux as `D*^-0.2 ~ A*^-0.1`, so halving the throat area raises it by
about 7% -- a factor 1.6 across a factor 140 in area. The **blowdown** goes as `1/A*`. Since
fluence is flux times time, the dwell carries essentially the whole cost: a narrow throat does
not heat the throat harder, it heats it for twice as long.

`throat_life` therefore reports a **pulse count** rather than a flux or a temperature, which is
only meaningful because of Decision 7: W15 established that the throat plates nothing back at any
survivable wall temperature, so recession is a consumable rate rather than a survival question.
Recession is radial and the throat *opens* as it erodes, so the metric is fractional area growth
-- the nozzle drifts back up the trade curve it was narrowed to climb.

**The exponent is derived before the solver runs and is asserted, not fitted.** Recession goes as
`A*^-1.1`; a fixed fractional area growth needs a recession proportional to `r* ~ A*^0.5`; so
throat life goes as **`A*^1.6`**, a factor 3.03 per halving. The joined Bartz/blowdown/ablation
chain reproduces it to 2% across the ladder, and `THROAT_LIFE_EXPONENT` plus
`test_throat_life_follows_the_analytic_area_exponent` pin it -- the same analytic-first discipline
the kernels are built under, applied to a correlation whose absolute calibration is not
trustworthy. It is what lets W24 quote an exchange rate while disowning the millimetres.

The trade table is joined in `wall.py` rather than `surface.py` because the cost side is the
wall's -- flux, dwell, ablation -- while the gain side is one column read off the N16 grid.

## Consequences

- New module `python/puffsat/walled_nozzle/surface.py` (N16) and
  `python/puffsat/walled_nozzle/wall.py` (N9 items 1-7).
- Item 5's verdict is that the throat is a **consumable**: it is chemically eroded on top of
  being thermally ablated, and no redeposition credit is available at any survivable wall
  temperature. The companion's "0.34 to 3.7% redeposition suffices" line prices a mechanism that
  does not operate.
- New targets `make walled-nozzle-surface` and `make walled-nozzle-wall`, both
  added to `make walled-nozzle-test`.
- Six new committed CSVs under `data/results/walled_nozzle/`, un-ignored for the
  cross-repo reason ADR-0025 and the existing walled-nozzle entries already record.
- **Every specific impulse in both answer documents is corrected.** No ranking changes: the N16
  grid's temperature, volume, throat and methane-over-water verdicts all survive, with the
  temperature margin narrowing from 1.74x to 1.56x, and W8's hydrogen > methane > water ordering
  survives with hydrogen's lead compressed from 1.88x to 1.55x.
- **Three W-series conclusions did move**, all in the ladder rather than the grid, and all
  recorded in place: ammonia's crossover with methane falls from 0.601 to 0.477 conversion (a
  factor 1.18 above methane, no longer "a large ask"); ammonia at the bottom of its assumed
  conversion range now falls *below* water rather than clearing it; and W5's two ADR-0016 anchor
  rows, previously 2.9% apart with the gap flagged as unexplained, close to 1.0% -- the head-on
  debit was most of what they disagreed about, which is independent evidence that the companion's
  own column already carries the `-1`.
- **`chamber.isp_scaling` moved with it**: the 200 -> 673 m^3 gain reads 1.3% rather than 1.8%,
  strengthening the "chamber volume is not an Isp dial" verdict.
- The answer document is `docs/walled_nozzle_grid_and_wall.md`, numbered W10
  onward so it continues the W series the paper repository already cites.
- **W24 asks the paper for something it does not currently state: a throat-replacement
  interval.** The recommendation is 2 m^2 where the throat is hard to service and 1 m^2 where it
  is a scheduled consumable, and total impulse per throat (pulses x Isp) falls 50x between them --
  so the choice is an economics one and the study supplies only the exchange rate.
