# Water-injected overtake plate study

Status: scope, analytic ledger, kinematic profiles, local equilibrium mixing,
conditional recombination-location screen, and a provisional planar stratified-water
flow reference implemented. The reference shows temporal overlap between molecular
store return and wall impulse. No distributed-injection performance, isolated causal
chemistry contribution to impulse, or survivability result yet.

Source: [companion proposal at 9440b27](https://github.com/katzseth22202/Balloon-Pulse-Propulsion/blob/9440b2789195d14d18409b4eb619c80dba44f71d/docs/water_injected_overtake_plate_for_impact_sim.md).
User refinements recorded 2026-09-05 are listed below. This is distinct from the
original per-collision `f(v)` and head-on tamped-nozzle studies.

## Question and agreed design freedom

Can a material plate receiving a dispersed 25 kg water PuffSat and 125–250 kg of
sprayed water deliver more useful impulse per carried consumable kg than the
magnetic nozzle during the overtake growth push? A gain from accelerating added
mass counts as scientific success even without recombination recovery. Recovering
chemical energy while the gas still pushes would be an additional advantage.
Engineering simplicity is a potential benefit, not an established conclusion.

- Sweep injected-water ratio `k = m_water/m_projectile` from 5 to 10, retaining
  8.5 and the unsprayed k = 0 control. Refine near useful cases and extend the
  search if the best result sits at an endpoint.
- Optimize incoming footprint and duration; the projectile arrives dispersed.
  For a uniform pulse, `rho_in = m_projectile/(pi*r_foot^2*w*duration)`.
  This specifies incoming density, not shocked mixture density. Numerical search
  bounds and the initial gas temperature/profile remain to be specified.
- Compare pre-spray, during-pulse injection, and splits at fixed total water mass.
- Projected plate diameter `D <= 30 m` (`R <= 15 m`), mass no more than 50 tons.
  Working mass interpretation is 50,000 kg of permanent plate assembly;
  inclusion of injection hardware/supports must be explicit before structural
  screening. Consumable inventory is separate. No structure has been sized yet.
- Prefer smaller/lighter designs when enlargement offers only marginal benefit.
  Report the performance/mass/diameter trade and diminishing returns. No numeric
  penalty weight or acceptable performance sacrifice has been chosen.
- Reopen concavity and taper for this encounter. Flat and parabolic
  `delta/D = 0.10, 0.15` are initial geometry anchors, not a ban on other shapes.
  Thickness taper versus conical surface taper remains to be clarified. Central
  peak stress is a hypothesis; loading may be concentrated elsewhere.
- Prioritize thin protection: screen 0, 1, 3, 10, 30, and 100 micrometers, with
  100 micrometers retained as a conservative comparison, not the default design.
  Zero means no separate sacrificial coating; bulk reaction-water injection still
  operates. Adequacy follows wall loading, not thickness alone. These are design
  trials, not verified Orion performance figures.
- Ablator elements must be normally present in the atmosphere: carbon allowed,
  lead excluded. A C/H/O-only formulation is a provisional candidate family;
  specific ingredients and additives have not been selected. The old study's
  silicone/material prescriptions are not imported automatically.

The four speed anchors are 45.58 and 56.53 km/s (3-synodic overtake), and 61.83
and 65.13 km/s (2-synodic overtake). The proposed first hydrodynamic anchor is
45.58 km/s, k = 8.5, with controls. The 75 km/s head-on ledger does not apply.

## Stages and accounting

Stage 1 assumes gasification and immediate local mixing/thermal equilibration
where flows meet; it does not redistribute mass throughout the bowl instantly.
The collision pays for heating and gasification. Injection profiles must specify
position, velocity, temperature, density, and timing with consistent source energy.
Stage 2 follows reported Stage 1 results and derives droplet sizes and timing
that retain the useful gain. Injector engineering is not a Stage 1 prerequisite.

Use the initial vehicle frame, +z prograde. The projectile supplies incoming
momentum `+m_projectile*w`. At settled completion,
`J_vehicle = m_projectile*w - P_out,z`, including radiation. At a finite stop,
also subtract residual stored momentum. Reconcile this with plate and feed-system
forces; do not add incoming momentum a second time to a complete wall impulse.

At each geometry, compare matched incoming mass/speed/footprint/duration in the
unsprayed, injection-only, and combined cases. The collision-induced gain is
`G = J_combined - J_injection_only - J_unsprayed`. Report net J and G separately.
Future droplet retention contours compare G only when ideal G is positive.

Report beta, impulse per projectile kg, and impulse per expended carried kg.
Effective Isp is `J/(g0*m_charged)` only for nonzero charged mass. Charge missed
water, lost coating, spray waste, and other expendables. If protective water is
part of the reaction-water allocation, count it once. Ablator mass leaving as
exhaust affects both the ejecta-mass energy ceiling and the consumable denominator.

Compare water and argon nozzles with all composition differences named. Nozzle
geometry remains a performance range/input unless computed. Report the nozzle
efficiency each simulated plate would beat rather than inventing one benchmark.

Track energy never spent on dissociation separately from energy returned before
pressure coupling ends; separate molecular and ionization stores. Carbon-bearing
ablator products may bind oxygen and change that chemistry. Equilibrium alone
does not establish recovery. Follow density/species/energy against cumulative
impulse and check reaction times through the cooling trajectory.

## Step 0: reproducible analytic reference

Run `make water-plate-ledger`; verify with `make water-plate-test`.
Implementation: [`ledger.py`](../python/puffsat/water_plate/ledger.py).
The two generated, versioned tables are:

- [`ledger_pulses.csv`](../data/results/water_plate/ledger_pulses.csv): four speeds
  crossed with k = 0, 5, 6, 7, 8, 8.5, 9, 10 (32 rows).
- [`ledger_coatings.csv`](../data/results/water_plate/ledger_coatings.csv): 10/20/30 m
  diameters, flat/0.10/0.15 parabolic depths, and 0/1/3/10/30/100 micrometer layers
  (54 rows).
  These diameters are inventory examples, not the selected optimization grid.

The ideal no-ablator ceiling is `beta_max = 1 + sqrt(1+k)`. It assumes negligible
vehicle recoil, initially cold stationary water, no added injection energy,
complete ejection, and loss-free single-speed backward exhaust. This is a common
energy/momentum upper reference for material and magnetic redirection, not an
achievable plate prediction. The free inelastic merge identity partitions input
KE into bulk drift and heat; further plate stagnation can also use the bulk part.

At the 45.58 km/s anchor, incoming KE is 25.969 GJ:

| k | Injected water (kg) | Gasification scale (GJ) | Share of incoming KE | Ideal impulse (MN s) | Water-only ideal Isp (s) |
| --- | --- | --- | --- | --- | --- |
| 5 | 125 | 0.325 | 1.25% | 3.931 | 3207 |
| 8.5 | 212.5 | 0.5525 | 2.13% | 4.652 | 2232 |
| 10 | 250 | 0.650 | 2.50% | 4.919 | 2006 |

These Isp ceilings exclude ablator and all other expendables. Their decline with
k while total impulse rises demonstrates why the two objectives are distinct;
it does not locate the optimum of the real device.

Gasification is a rounded 2.6 MJ/kg engineering scale for liquid near 293 K heated
to boiling and vaporized. It is not a vacuum phase trajectory or an EOS zero-point
offset. Molecular dissociation uses `eos_water.FULL_ATOMIZATION_ENERGY`, about
50.94 MJ/kg of water vapor, excluding ionization. At k = 8.5 full dissociation of
all 237.5 kg represents 12.098 GJ, or 46.59% of incoming KE. That is a possible
store size, not a measured stranded fraction or a claim that all water dissociates.
Actual source enthalpies and final retained stores require consistent EOS states.
Do not permanently subtract every intermediate store from the ideal ceiling;
recoverable energy must not be charged twice.

Coating inventory is `rho_layer*h*A_actual`, with thickness normal to the surface.
At illustrative density 1000 kg/m3 and 100 micrometers, full-face inventories for
10/20/30 m flat plates are 7.85/31.42/70.69 kg. At 30 m, parabolic depths 0.10/0.15
raise inventory to 73.44/76.71 kg. A full fresh flat layer each pulse would add
33.26% to the k = 8.5 water bill. Consumption could be smaller if material remains;
replenish the computed loss and account for deposition waste. Restricting coverage
requires resolving the wall heating footprint, including radiation outside the
incoming footprint. No coating survival or thermal properties have been assumed.

At the maximum flat diameter, the thinner 1/3/10 micrometer layers hold
0.707/2.121/7.069 kg, respectively, at the same illustrative density. Full renewal
would add 0.33/1.00/3.33% to the k = 8.5 water bill. These inventory figures support
exploring thin protection; they do not establish that any of these layers survives.

## Progress and findings to take back to the paper

This is the living receiving-repository summary. Update this checklist and attach
evidence as each step completes. Paper prose is edited separately; no claim of
architecture superiority is ready to transfer.

- [x] Read and pin the companion proposal; record the user-refined scope.
- [x] Implement the overtake mass/energy/impulse reference and its conservation checks.
- [x] Quantify gasification versus full molecular-dissociation energy scales.
- [x] Quantify coating inventory versus diameter, curvature, and thinner layers.
- [x] Specify and mass-normalize bounded incoming pulse and injection profiles.
- [x] Align stored liquid and vapor energy references; test collision-paid local mixing.
- [x] Screen local equilibrium states for initial dissociation avoidance and energy stores.
- [x] Locate equilibrium store return along expansion paths and screen reaction clocks
  at the first speed (conditional density/geometry assumptions; Step 3).
- [x] Run a planar prepositioned-water reference with conservative source bookkeeping,
  matched wall controls, spatial/time-step checks and actual parcel chemistry clocks.
- [ ] Implement the flow source/boundary with the same energy reference and feed work.
- [ ] Validate the injected-flow calculation against analytic limits and convergence.
- [ ] Measure combined, unsprayed, and injection-only impulse at the first speed.
- [ ] Resolve chemistry and wall loading through the impulse-producing interval.
- [ ] Optimize k, footprint, duration, geometry, and timing; check search boundaries.
- [ ] Determine minimum adequate protection and assess the permanent mass limit.
- [ ] Compare against water/argon nozzle ranges at matched consumable accounting.
- [ ] Report Stage 1, then derive Stage 2 droplet and timing requirements.

| Finding | Evidence | Paper disposition |
| --- | --- | --- |
| Bulk-momentum cancellation alone does not bound hot-mixture plate impulse at 2f; the common ideal overtake ceiling is 1 + sqrt(1+k). | Conservation identities and `ledger_pulses.csv`; already identified by the proposal. | Ready as an accounting qualification, not a plate-performance claim. |
| At the 45.58 km/s, k = 8.5 reference, gasification costs about 2.13% of incoming KE; full molecular dissociation represents a 46.59% store. | `ledger_pulses.csv`, with the reference-state qualifications above. | Useful motivation; neither actual dissociation nor recovered energy has been measured. |
| A 100 micrometer full-face layer need not be a small consumable; 1–10 micrometers is much cheaper in inventory. | `ledger_coatings.csv`; actual area and illustrative density explicit. | Design trade only; do not prescribe a protective thickness yet. |
| Distributed mixing depends on radial overlap: at r_foot/R = 0.75, uniform full-plate injection places 56.25% of its mass inside the incoming footprint; an r-squared source places 31.64% there. | `profile_injection.csv` and independently integrated source normalization tests. | Input geometry only. The outer water can still interact with outward flow; these fractions are not capture efficiency. |
| In local equilibrium mixing at 45.58 km/s, k = 5–10 and assumed mixed densities 0.1–10 kg/m3 retain 97.9–100% of the full molecular dissociation store immediately after the merge. | `thermo_mixing.csv`; common liquid/vapor energy reference and source-state sensitivity below. | Conditional finding: simple uniform dilution in this box does not avoid most initial dissociation. Recombination during later pressure coupling remains untested. |
| At 45.58 km/s and k = 8.5, dense initial mixtures allow some molecular return much earlier than dilute ones: 10% equilibrium return occurs at equal-volume hemisphere radii 37.4/13.8/5.02 m for initial densities 0.1/1/10 kg/m3. | `chemistry_milestones.csv`; Step 3 states the geometry and rate limitations. | Conditional screening result, not a plate impulse prediction. The dense case's first 10% is kinetically plausible under both rate proxies; most of the bond store remains a later, slower process. |
| The planar stratified-water reference has molecular store return within centimetres of the wall during the main impulse interval; its 1 ms impulse exceeds the unsprayed control. | `flow_reference_milestones.csv` and `flow_reference_convergence.csv`; Step 4 states the EOS, sampling and control limitations. | Provisional reference result only: no distributed spray, finite plate, causal chemistry increment or survivability claim. |
| Water injection beats a bare plate or a magnetic nozzle, recovers recombination energy, or simplifies survivable hardware. | No simulation/engineering evidence yet. | No result to transfer; corresponding checklist items remain open. |

Step 3 checkpoint verification: 39 study acceptance tests pass; `make lint` passes
(including strict Python typing, clippy, and Rust formatting). Step 3 includes
the resolution comparison below. Full Rust tests pass in release mode: 175 passed,
11 repository-marked ignored. A pre-existing float-to-integer cast in the Rust
test's text histogram was replaced with a bounded integer count to clear the
repository's denied cast lints; its report test also passes after that edit.
Full Python verification passes: **477 passed** with the `sci` extra, in 23m40s.
Cargo, rustfmt, clippy, and the system linker have now been installed; the old
missing-Cargo limitation is resolved. Source `/home/agent/.cargo/env` in this
environment to expose Cargo. Rust builds use a fresh temporary `CARGO_TARGET_DIR`
because the existing `target/` produced directory-creation errors. The slow debug
run was stopped and replaced by the completed `cargo test --release --locked` run.

## Restart checkpoint

- Branch: `main`; find the latest study checkpoint with
  `git log --oneline --grep='water.*plate'` and inspect `git status --short` first.
- Completed: Step 0 (08dc8f3), kinematic profiles (3485abd), and the local equilibrium
  mixing screen (f9c2df8), followed by the conditional chemistry-location screen
  (2026-09-06 working-tree checkpoint; inspect `git status` before editing).
  Reproduce with `make water-plate-ledger`, `make water-plate-profiles`,
  `make water-plate-thermo`, `make water-plate-chemistry`, and `make water-plate-test`
  (last three use the `sci` extra).
- Active next step: implement and validate the first conservative injected-flow
  calculation, carrying the liquid/gas energy offset and the feed-work ledger.
  Extract parcel density/temperature and separate chemical stores against cumulative
  wall impulse. Step 3 prioritizes whether a dense cooling layer actually forms and
  stays pressure-coupled; it does not supply that layer as an initial condition.
  No hydrodynamic result exists yet; the local mixed density is still an input.
- Scratch: `todos/water_plate_study_scope.md` points here; it is not the canonical
  task state. The companion checkout is disposable and pinned at 9440b27.
- Keep this checkpoint and the progress checklist synchronized with completed
  work. Record paper-ready findings with their evidence and limits in the table.

## Step 1: bounded kinematic input profiles

`make water-plate-profiles` writes [`profile_pulses.csv`](../data/results/water_plate/profile_pulses.csv)
(45 cases) and [`profile_injection.csv`](../data/results/water_plate/profile_injection.csv)
(18 cases). These are screening inputs, not a parameter optimization or flow solve.
Their spatial and temporal mass integrals are checked independently by quadrature.

The incoming screen at 45.58 km/s uses diameters 10/20/30 m, r_foot/R = 0.5/0.75/1,
and durations 10/30/100/300/1000 microseconds. These are provisional search bounds;
extend any binding edge. The initial shape is a uniform cylindrical pulse with
zero radial divergence. Density is derived from 25 kg and never held fixed while
changing area or duration. Incoming ram pressure is an input scale, not wall pressure.

The first illustrative reference is D = 10 m, r_foot/R = 0.75, duration = 100
microseconds: length 4.558 m, incoming density 0.12415 kg/m3, ram pressure 257.93 MPa.
This is not a chosen optimum or a structurally qualified plate.

Injection examples use k = 8.5 and compare uniform footprint-matched delivery,
uniform full-plate delivery, and full-plate delivery weighted by r-squared. Flux
is per projected area. For a curved plate convert to actual face flux using its
local normal; source momentum is initially axial, not implicitly surface-normal.
The cumulative radial mass fraction is (r/R_source)^(power+2), clipped at one.

Temporal examples are all pre-spray, all during-pulse, or a half-and-half split.
The pre-spray occupies [-2, -1] ms and the during part [0, 0.1] ms; t = 0 is
unperturbed front arrival at the vertex plane. A solver must start early enough
to represent the pre-spray and its earlier encounter with the incoming pulse,
not initialize it at t = 0 as though no collision could already have occurred.
Injection speeds are 10 and 100 m/s toward the pulse. Lead times, source duration,
and velocities are illustrative and will need a wider scan once the solver exists.

For the footprint-matched example, spreading 212.5 kg over the 1 ms pre-spray needs
212,500 kg/s, versus 2,125,000 kg/s during the 100 microsecond pulse. At 100 m/s the
area-mean equivalent stream densities are 48.1 and 481.0 kg/m3. These are continuity
requirements m_dot/(A*u), not physical equilibrium vapor densities or shocked
mixture densities. The r-squared case's local rim flux is twice its area mean.
Such rates are recorded without making injector engineering a prerequisite to
the ideal performance screen. Actual prespray expansion changes density/overlap.

At 100 m/s the injected water carries 21,250 N s of backward momentum and 1.0625 MJ
of KE (0.0041% of incoming KE). That is a prescribed kinetic input; actual feed
pressure work is still unspecified. Injection-only control forces, not this scale
alone, determine the recoil subtraction in a flow calculation.

## Step 2: local equilibrium mixing with collision-paid gasification

`make water-plate-thermo` produces [`thermo_mixing.csv`](../data/results/water_plate/thermo_mixing.csv)
(84 states: four speeds, seven nonzero k values, three mixed densities) and
[`thermo_sensitivity.csv`](../data/results/water_plate/thermo_sensitivity.csv)
(24 source/reference variations at the cold k = 8.5 anchor).

This is a closed-parcel calculation. It conserves signed momentum and total energy
when projectile vapor and carried liquid mix inelastically, then inverts the
existing equilibrium EOS at a chosen final gas density. It does not impose a
global mixture on the bowl, compute where mixing occurs, or equate the local k to
the global k without qualification. It represents a parcel where that ratio is
realized. No wall stagnation, radiation, expansion, or species evolution is solved.

Stored water is liquid at 293.15 K and 101325 Pa. Match CoolProp/IAPWS-95 internal
energy to `eos_water` at dilute molecular vapor (400 K, 0.01 kg/m3). Add the resulting
constant offset to the liquid's internal energy. This yields about -1.914 MJ/kg
on the bound-molecular-gas zero, with CoolProp 7.2.0; negative is physically correct
on this reference. Heating/gasifying it to the 400 K vapor reference costs about
2.471 MJ/kg of internal energy. The earlier 2.6 MJ/kg remains a rounded enthalpy
scale, not a second debit. No dissociation cost is subtracted separately because
the gas EOS already includes it. The bridge transfers an internal-energy
difference; it is not a globally stitched phase/plasma EOS for a flow solver.

Incoming vapor is provisionally 400 K at 0.12415 kg/m3, the input reference density.
This cold molecular energy is held fixed for the local screen across speeds; the
screen does not assert that the resulting mixture follows the input geometry.
Mixed gas density is sampled at 0.1, 1, and 10 kg/m3, independently of stream
density. It must ultimately come from the flow, not be fitted to aid recombination.
Baseline injection velocity is zero; a 100 m/s opposing velocity is included in
sensitivity. For local ratio k, the merge velocity is `(w-k*u)/(1+k)` and the
heat generated per projectile kg is `k*(w+u)^2/[2*(1+k)]`. The remaining bulk KE
is retained explicitly. The plate can later dissipate more of it.

At 45.58 km/s, k = 8.5, with stationary carried liquid:

| Assumed mixed density (kg/m3) | Equilibrium temperature (K) | Molecular bond store held | Ionization store (MJ/kg mixture) |
| --- | --- | --- | --- |
| 0.1 | 13,947 | 99.974% | 14.38 |
| 1 | 15,998 | 99.827% | 10.50 |
| 10 | 18,242 | 98.775% | 6.73 |

The molecular fraction here is **stored bond energy divided by full atomization
energy**, including OH/H2/O2 bonds, not the fraction of molecules destroyed. Across
all k = 5–10 at this speed it spans 97.897–99.991%. At these sampled local states,
adding this much water does not avoid most initial dissociation. This is evidence
for prioritizing the later recombination/pressure-coupling question, not evidence
that water injection cannot help or that equilibrium recovery will occur.

Sensitivity rows vary the dilute matching temperature 400/600 K, incoming vapor
temperature 400/800 K, and injection speed 0/100 m/s. The energy-offset change
between the two matching states is 1.295 kJ/kg (small compared with the phase
energy); the table records the corresponding solved states rather than silently
treating the two EOS models as identical. No extra composition or ablator has yet
been added. Carbon effects remain a separate requirement.

## Step 3: where the chemical stores could return

`make water-plate-chemistry` follows the collision-paid k = 8.5 mixing states at
45.58 km/s (3-synodic overtake), for the same assumed initial mixed densities
0.1, 1, and 10 kg/m3. It writes:

- [`chemistry_milestones.csv`](../data/results/water_plate/chemistry_milestones.csv):
  first 10/50/90% return of each initial store; absent crossings are explicit.
- [`chemistry_sizes.csv`](../data/results/water_plate/chemistry_sizes.csv):
  states at equal-volume hemisphere radii 5/10/15 m, when on the trajectory.
- [`chemistry_convergence.csv`](../data/results/water_plate/chemistry_convergence.csv):
  central-density milestones at 128/256/512 integration steps.
- [`chemistry.png`](../data/results/water_plate/chemistry.png): stores and reaction
  clocks against the illustrative size. Full sampled histories are generated but
  gitignored (`chemistry_history_rho_*.csv`).

**Method and limits.** Reuse `expansion.expand` to integrate the reversible,
adiabatic relation `de = -p d(1/rho)` through a 100,000-fold volume increase.
This imports the integrator and water EOS, not the magnetic nozzle's density,
field geometry, chamber size, or residence time. Radiation, conduction, additional
wall stagnation, feed sources, and finite-rate species evolution are absent.
The two stores are measured as energy per kg, including intermediate molecular
bonds and all EOS oxygen ion stages. Return fractions reference the **initial
store**, not the full-atomization reference used in Step 2. Small negative molecular
returns near the initial state mean a little further dissociation, not a loss
of normalization. Internal-energy decrease is parcel expansion work, not wall
impulse, nor work uniquely attributable to chemistry.

Density locates a state without requiring a geometry. For a size/time illustration
only, put all 237.5 kg in a uniform hemisphere:
`r = [3m/(2*pi*rho)]^(1/3)`, prescribe `r = r0 + U*t`, and use `U = 10 km/s`.
Then `tau_exp = |dlnrho/dt|^-1 = r/(3U)`. Initial radii are 10.43, 4.84, and
2.25 m respectively. These are **not nozzle stations, actual source-to-wall
distances, or simulated flow boundaries**. A real distributed spray need not form
a hemisphere. The initial kinetic energy and subsequent acceleration of this
prescribed geometry are not solved. Changing U rescales times and Damkohler
numbers as 1/U, without changing equilibrium states or sizes. The 15 m comparison
is the allowed maximum plate radius, **not an assumed cutoff of wall force**:
gas outside the rim can still exert pressure on the plate, and gas inside it may
already be escaping. Neither size nor elapsed time replaces cumulative impulse.

Reaction rates reuse the literature coefficients documented in `recombination.py`.
Ionization has separate collisional and radiative capture clocks; radiated photons
are not credited as gas heat. The ion clock is a hydrogenic screening proxy, not
an oxygen charge-state network. Molecular rates use `H + OH + M -> H2O + M`, with
the actual equilibrium neutral third-body density and either equilibrium OH or
the older H-density substitution as partner. These are **sensitivity proxies,
not mathematically proven bounds** on the full OH/H2/O2/water network. The
low-pressure coefficient, steam third-body efficiencies, pressure falloff, and
extrapolation to these temperatures still require validation. `Da = tau_exp/tau_rec`
above 10 supports fast chemistry; below 0.1 indicates slow chemistry for that
proxy. Small Da where the relevant store is already empty is not a freeze verdict.
Beyond an inconsistent reaction clock, an equilibrium curve is a hypothetical
path, not a prediction of recovered energy.

**Results at the first speed.** Molecular milestones, with all sizes and times
conditional on the geometry above:

| Initial mixed density (kg/m3) | T at 10% molecular return (K) | Density there (kg/m3) | Radius there (m) | Time there (ms) | Da, H / equilibrium-OH proxies | Radius at 50% return (m) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.1 | 3783 | 0.00217 | 37.41 | 2.698 | 0.135 / 0.00635 | Not reached by 484 m |
| 1 | 4808 | 0.0434 | 13.77 | 0.893 | 12.2 / 0.648 | 131.6 |
| 10 | 6577 | 0.896 | 5.02 | 0.277 | 987 / 62.9 | 34.28 |

At 15 m the equilibrium molecular returns are **0.0028%, 12.15%, and 36.44%**
respectively. The high-density path returns 9.88% by 5 m and 27.99% by 10 m;
its equilibrium-OH Da falls from 63.3 to 7.35 to 1.30 across 5/10/15 m. Thus the
first roughly 10% is kinetically plausible in this dense case, while the full
36% cannot be assumed to track equilibrium. At its 50% milestone even the
H-density proxy has Da = 0.0923 (OH: 0.0272), so the late equilibrium recovery
already fails this screen. None of the three paths reaches 90% molecular return
within the chosen integration domain; that is not proof it never returns.

Ionization is earlier and faster. The 90% return points are 9380/10468/11581 K,
at radii 16.84/7.21/3.17 m and elapsed times 0.642/0.237/0.093 ms for the three
initial densities. Collisional-ion Da there is about 3.1e3/7.1e4/1.1e6, with
collisional capture faster than radiative capture. This supports early ion-store
return along these assumed paths; a pusher plate's useful share remains unmeasured.
The magnetic-nozzle study already credits early ionization recovery, so this
alone is not an advantage unique to a material plate.

**Interpretation.** A dense cooling layer is the promising configuration to test
for **partial molecular recovery during pressure coupling**. A dilute freely
expanding mixture does not gain access to most of the bond store on plate-sized
scales in this illustration. The dense state must arise from the pulse/spray/wall
solution, not be imposed to obtain this result. The next decisive output is
chemical energy versus cumulative wall impulse, including the actual local
expansion rate. This screen establishes neither a thrust gain nor the ranking
against an unsprayed plate or either magnetic-nozzle composition.

**Numerics.** Doubling 256 to 512 steps changes central-density molecular crossing
densities by 0.091% (10% return) and 0.368% (50% return); successive 128/256/512
density differences shrink by factors 4.08/4.09, consistent with second order.
The three ionization crossing densities change by at most 0.021%; their
interpolation errors are not uniformly second order. The reused integrator has
independent ideal-gas analytic and order tests. New acceptance tests cover the
mass-conserving clock and its density derivative, missing/interpolated crossings,
source-store consistency, and separation of radiative from collisional capture.
This is numerical convergence within the chosen EOS, not validation of its
equilibrium or rate assumptions.

## Step 4: planar stratified-water flow reference

**Result within this model:** molecular store return begins close enough to the
wall and early enough to overlap the main pressure impulse. This is a stronger
result than the prescribed hemisphere in Step 3, but not yet an isolated thrust
increment caused by recombination. The water layer also changes mass loading,
shock transmission and the amount of energy ever spent on dissociation.

Evidence:

- [`flow_reference_convergence.csv`](../data/results/water_plate/flow_reference_convergence.csv)
  records matched controls, conservation errors, resolution/sampling checks and
  sampled chemistry peaks.
- [`flow_reference_milestones.csv`](../data/results/water_plate/flow_reference_milestones.csv)
  puts parcel stores and local reaction clocks beside cumulative wall impulse.
- [`flow_reference_model.json`](../data/results/water_plate/flow_reference_model.json)
  pins the reference table provenance and exact source/output-time inputs.
- [`flow_reference.png`](../data/results/water_plate/flow_reference.png) shows
  wall impulse, held stores and sampled per-parcel peak decreases separately.

### What is solved, and what is not

Use 25 kg incoming water at 45.58 km/s (3-synodic overtake), 400 K and a
100 microsecond pulse, with footprint radius 3.75 m. The constant planar area is
44.179 m2, pulse length 4.558 m and incoming density 0.12415 kg/m3. Place
212.5 kg of stationary water in a 1 m layer touching a rigid, adiabatic wall;
its initial bulk density is 4.8100 kg/m3. It carries stored-liquid internal
energy, not the larger energy of cold ideal steam. Equilibrium flashing at
that density gives about 291.36 K. The collision pays subsequent heating and
phase changes. No feed energy or momentum has been added.

The existing 1D Lagrangian kernel tracks the material contact and its pressure
waves. There is **no interpenetration or heat exchange across that contact**,
distributed source, droplet physics, radial escape, finite plate boundary,
radiation, wall heat loss, ablation or plate recoil. This is a stratified
prepositioned-water reference, **not Stage 1's immediate local spray mixing**.
The footprint is an area normalization, not a selected finite plate design.
The 50,000 kg permanent assembly has not been modeled or sized.

The cold/dense EOS and rejected attempts are documented in
[ADR-0043](adr/0043-water-plate-cold-source-and-provisional-dense-eos.md).
It uses CoolProp/IAPWS below 1000 K and a thermodynamically matched residual
extension of the existing equilibrium gas above it. A common 4 MJ/kg storage
offset keeps the log-energy table positive and is subtracted from physical
energy reports. This dense/hot continuation is provisional; high-pressure ice
is absent and IAPWS itself is extrapolated above 1 GPa. Positive heat capacity
and sound speed are necessary checks, not physical validation of warm-dense water.
The full nominal table spans 1e-5–2000 kg/m3 and 273.2 K–2 MK; its worst checked
Maxwell residual is 2.52e-5 relative. The alternative residual exponent gives
1.95e-5. No negative heat capacity or out-of-domain flow is silently accepted.

The cold-only control needs special care: refining its free surface eventually
cools it below the liquid/vapor table. The correct short-window wall control
is available analytically. The rarefaction head reaches the wall only after
`L/c_initial = 0.14665 s`, so for the 1 ms comparison its wall pressure stays
2081.09 Pa and `J_water_only = 91.94 N s`. This agrees with the coarse numerical
controls before their domain exits. It is explicitly marked analytic and has
no invented downstream flow history or measured numerical conservation error.

### Wall impulse and numerical checks

All numbers below stop at **1 ms**, not a settled-outflow criterion. N is the
incoming cell count; the combined case has 8.5N water cells of equal mass.
The default reference reduces the existing CFL 0.4 timestep to CFL 0.1.

| N | Combined J (MN s) | Unsprayed J (MN s) | Matched G (MN s) |
| --- | --- | --- | --- |
| 20 | 2.60558 | 1.92563 | 0.67986 |
| 40 | 2.62539 | 1.92026 | 0.70504 |
| 80 | 2.63440 | 1.91640 | 0.71790 |

From N=40 to 80 the individual impulses change by 0.34% and 0.20%, but the
subtracted gain still changes by 1.8%. Do not attach a sub-percent convergence
claim to G. At N=80 the maximum total-energy drift is 0.136% combined and
0.315% unsprayed; momentum closes to about 1e-8 N s. At N=20, halving the
timestep again to CFL 0.05 changes combined J by 0.037% and G by 0.059%.
Changing the dense residual exponent from 2 to 3 at N=40 changes combined J
by 0.075% and G by 0.28%. That is a **limited closure sensitivity**, not a
validated model uncertainty interval. EOS-table refinement, including the
log-interpolation dependence on storage offset, remains open.

Earlier CFL 0.4 runs had up to about 1.5% shock energy drift; spatial refinement
alone did not remove it. Those runs fail the 0.5% energy gate and provide no
matched gain. The new timestep-control test checks decreasing energy drift
against an energy-conserving ideal-gas impact. Mass-weighted nodal initialization
also has explicit local projection heat and an initial wall impulse; these are
numerical interface effects, recorded and reduced with spatial refinement.

The fine spatial reference has approximately **5.79 GPa peak physical wall
pressure**, versus **0.316 GPa unsprayed**. Thus the extra impulse comes with a
much sharper wall load in this geometry. Neither coating adequacy nor structural
survival follows from this pressure-only calculation. At 1 ms the combined wall
still has about 1.85 MPa pressure; the reported impulse is not its asymptote.

### Chemistry during the push

The timeline uses N=40 with additional sub-microsecond output near the sharp
wall-loading interval. Chemical stores decompose the ideal-gas part of the
extended EOS; they do not include its dense-fluid residual. Each fixed-mass
parcel retains its own sampled maximum store. The reported decrease is
`sum m_j*(sampled_peak_j - current_store_j)`, not a net drop from one global
peak, a reaction-energy flux into the plate, or work uniquely caused by chemistry.
The midpoint locations weight **positive store decreases within the preceding
sampling interval**; they approximate where those changes occur, not the
location of a separately tracked reaction front.

Coarse logarithmic output under-sampled those parcel peaks: 64 and 128 output
times gave 2.317 and 2.488 GJ of molecular peak decrease by 1 ms. Adding 1.0
microsecond samples across 0.10–0.25 ms raises this to 2.704 GJ; halving that
spacing to 0.5 microseconds gives **2.717 GJ**, a further 0.50% change. The
corresponding ion value changes by 0.056%, to 5.471 GJ. Wall impulse is effectively
unchanged by these output refinements. The milestone artifact uses the final
0.5 microsecond grid; spatial convergence of these parcel maxima still requires
equally dense time sampling on finer spatial grids.

The main wall impulse rises around 0.14–0.20 ms. Molecular store return is
already present within centimetres of the wall during this interval. By about
0.16 ms the summed parcel molecular stores have fallen by roughly 0.6 GJ from
their sampled peaks, while a substantial fraction of the 1 ms wall impulse
is still to come. By the sampled 90% impulse point near 0.19 ms, roughly
1.5 GJ of molecular store has returned. These quantities are equilibrium
store changes, **not 1.5 GJ measured as useful thrust work**. Detailed current
values and sampling brackets are in the milestone CSV.

More specifically, the first sampled 10% molecular peak decrease is at
**0.1575 ms**, with 0.603 GJ returned and **70.9%** of the 1 ms wall impulse
already delivered. The store-decrease-weighted interval midpoint is **1.69 cm**
from the wall. At **0.1925 ms**, 90.0% of the wall impulse is delivered and
the molecular decrease is **1.555 GJ**; that interval's midpoint is **11.6 cm**
from the wall. These are planar material positions, not hemisphere radii or
finite-plate capture distances. They establish partial temporal/spatial overlap,
not recovery of the entire chemical store.

Ion-store return starts earlier, before the transmitted pressure wave loads the
wall, but its late tail is still important: reaching half the summed parcel
ion peaks takes about 0.56 ms, after about 98% of the 1 ms wall impulse.
Late return should not be confused with the early, strongly pressure-coupled
part. Much of the incoming material remains hot; the stratified solution is
not the single uniform adiabat used in Step 3.

Reaction clocks now use the actual cell velocity gradient:
`tau_exp = width/(u_right-u_left)` only for expanding cells. Compression,
stationary cells, missing gradients and the cold phase get no invented expansion
clock. The CSV reports the **held store** in expanding gas with OH-partner or
collisional-ion Da above 10 or below 0.1, rather than the minimum Da of the
emptiest cell. The near-wall molecular return occurs in a regime where these
proxies support rapid chemistry. Their pressure falloff, third-body efficiencies
and extrapolation remain unvalidated; this is not a finite-rate network.

### Reproduction and verification

`make water-plate-flow-table`, then `make water-plate-flow WATER_FLOW_CELLS=80`
and `make water-plate-flow-analysis WATER_FLOW_CELLS=80` regenerate the fine
spatial reference. `WATER_FLOW_LAYER`, `WATER_FLOW_SAMPLES`,
`WATER_FLOW_TIMESTEP_FRACTION`, `WATER_FLOW_COUPLING_SAMPLES`, `WATER_FLOW_TABLE` and `WATER_FLOW_OUTPUT`
expose the checks; give variants distinct output filenames. The finer timeline
uses `python -m puffsat.water_plate.flow --cells 40 --layer 1 --coupling-samples 301`
under `PYTHONPATH=python uv run --extra sci`; its metadata pins the full time grid.
Generate the exponent-3 table with `flow_eos --residual-power 3 --output
data/tables/water_plate_flow_n3.json`. Analyze each raw file with `flow --analyze`,
then assemble the eight comparison runs using `make water-plate-flow-report`
(or `flow_report PATHS... --reference REFERENCE` for an explicit selection).
Bulk tables, JSONL and per-snapshot CSVs remain ignored; the four cited artifacts
are retained. Cargo is installed; in this environment source
`/home/agent/.cargo/env` and use `CARGO_TARGET_DIR=/tmp/puffsat-rust-check.ihcX6w`
because the inherited repository `target/` directory is not usable here.

Full checkpoint verification passed 512 Python tests (22m35s) and 180 Rust
tests in release mode (11 intentionally ignored). The final focused run passes
78 Python study tests plus 5 Rust history/driver tests. These checks
cover the updated history analysis, sampling inputs, analytic cold-only control,
thermodynamic transition grid and shock time refinement; `make lint` passes.

## Step 5: same-state chemistry-switch comparison

The planar reference now has an energy-preserving, parcel-specific reaction
switch, specified in [ADR-0044](adr/0044-water-plate-parcel-frozen-restarts.md).
Both branches restart the **identical** positions, nodal velocities, masses,
energies and accumulated impulse. One continues equilibrium chemistry; the
other holds every parcel's species mix fixed. The latter disables molecular
**and ionization** reactions in both directions, retaining the chemical store
and the same provisional dense residual. It is not a finite-rate calculation,
physical freeze-out or a molecular-only attribution.

The early switch is at **153.5 microseconds**, near peak molecular storage;
the late switch is at **193 microseconds**, just after about 90% of the
reference's 1 ms impulse has arrived. Some water-layer parcels still compress:
3.125 kg and 5.938 kg respectively at N=80. Thus freezing does not solely
disable recombination in an everywhere-expanding gas.

### Conditional impulse effect

Use a common **0.5 ms endpoint**. Subtract the chemistry response of the
unsprayed control as well as comparing the combined branches:
`Delta G = (J_combined,eq - J_combined,frozen)
- (J_unsprayed,eq - J_unsprayed,frozen)`.
The unchanged cold-only impulse is 45.97 N s; it is subtracted once in each
gain and cancels from Delta G.

| Switch (microseconds) | N | Combined response (kN s) | Unsprayed response (kN s) | Delta G (kN s) | Share of equilibrium G |
| --- | --- | --- | --- | --- | --- |
| 153.5 | 40 | 126.367 | 40.267 | 86.099 | 12.67% |
| 153.5 | 80 | 127.229 | 39.191 | 88.037 | 12.73% |
| 193 | 40 | 22.453 | 21.775 | 0.679 | 0.100% |
| 193 | 80 | 22.072 | 21.493 | 0.579 | 0.084% |

At N=80, early-switch equilibrium and frozen water-layer gains are
**0.691812 and 0.603775 MN s**. The extra 0.088037 MN s establishes that
continued chemistry after this switch contributes to the water-layer advantage
**within the provisional planar closure**. Refining N=40 to 80 changes that
effect by 2.25%, so it is not sub-percent converged. The late-switch difference
is much smaller and changes by about 15% under refinement; do not interpret its
digits as a precision measurement. Late chemistry helps both configurations,
with little additional water-specific advantage in this comparison.

The remaining 87% is **not** established as chemistry-independent: both branches
inherit all chemistry, dissociation avoidance and flow development before the
switch. Likewise Delta G is not additive molecular recombination work or a
strict bound on real finite-rate chemistry.

The early N=40 frozen combined run reaches its 1000 K lower domain limit at
**0.663512 ms**. Its attempted 1 ms comparison is explicitly invalid and reports
no matched gain. The late N=40 switch completes to 1 ms, giving Delta G =
1.470 kN s (0.209% of equilibrium G). This separate late result cannot extend
the early-switch result beyond its valid common window.

### Numerical checks and evidence

Actual switch p/T jumps are below 2e-11 relative, with no change of conserved
energy or momentum at the restart. The independently evaluated source EOS is
cross-checked before matching. At N=80 the largest restart energy drift is
1.22e-6 of physical energy. Small per-parcel constant-energy and entropy-pressure
corrections match the analytic frozen closure to the interpolated parent table:
summed absolute energy offsets are at most **30.98 MJ** (about 0.12% of physical
energy), and the largest pressure correction is **0.459%**. These are disclosed
closure adjustments, not an actual state jump or released chemical energy.
Table interpolation, storage-offset and residual-representation sensitivity
remain open; the corrections do not supply an uncertainty bound, especially
for the small late-switch difference.

Evidence retained in the repository:

- [Matched comparisons](../data/results/water_plate/freeze_comparison.csv), including the invalid probe.
- [Per-branch conservation and switch diagnostics](../data/results/water_plate/freeze_runs.csv).
- [Reference provenance](../data/results/water_plate/freeze_model.json).
- [Impulse history comparison](../data/results/water_plate/freeze_comparison.png).

Reproduce with `make water-plate-freeze-parent`, then
`make water-plate-freeze` (defaults: N=80, switches 153.5/193 microseconds,
endpoint 0.5 ms). `make water-plate-freeze-prepare` only prepares the restart
JSON. For the coarse evidence, use `WATER_FREEZE_CELLS=40` for the parent,
then also set `WATER_FREEZE_INPUT=data/results/water_plate/freeze_inputs_05ms.json`
and `WATER_FREEZE_OUTPUT=data/results/water_plate/freeze_05ms.jsonl` for the run.
For its 1 ms probe set `WATER_FREEZE_END=0.001`,
`WATER_FREEZE_INPUT=data/results/water_plate/freeze_inputs.json` and
`WATER_FREEZE_OUTPUT=data/results/water_plate/freeze_probe.jsonl`.
`make water-plate-freeze-report` assembles those three files, using N=80 as the
reference; override `WATER_FREEZE_REPORT_RUNS` and `WATER_FREEZE_REFERENCE` to
report another explicit selection. `WATER_FREEZE_TIMES_US` permits other switches
only when their exact checkpoints exist and all parcels satisfy the hot domain.
Bulk parent/restart files remain ignored. The flow table and Cargo environment
requirements are unchanged from Step 4.

## Next calculation

Separate molecular and ionization responses and validate finite-rate chemistry
on this controlled flow, with EOS table convergence and the dense-water closure
limits kept visible. The same-state all-reaction switch above is complete as a
conditional test, not a full chemistry attribution. Then
implement the actual locally mixing injection prescription and finite plate
escape. The positive reference G and early store return are reasons to continue,
not permission to transfer an injected-nozzle performance claim.

Carry the local energy reference into the distributed source/boundary prescription,
then run the locally mixing injected-flow comparison. Do not impose a cold ideal-gas state at the equivalent
stream densities above and silently grant it free vaporization energy. Stored
liquid and gaseous products must share a consistent energy zero, with local
conversion drawing energy from the collision. Initial incoming temperature is
provisionally 400 K, bracketed to 800 K in the local screen. Existing `euler2d` is lossless and fixed-gamma; it cannot
resolve the chemistry advantage. Injection, geometry, EOS, and wall losses must
be checked together before a performance claim. Report required coating thickness,
local pressure/heat flux, and permanent structure alongside consumable economy.
The original ADR-0014/0015/0021/0024/0026 results inform the checks but their thin
ablator, temperature, and incidence assumptions do not automatically transfer.
