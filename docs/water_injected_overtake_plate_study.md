# Water-injected overtake plate study

Status: Steps 0–9 and the **nonlinear prescribed-history parcel calculation**
(Step 11) implemented. Tested documented-rate scenarios return roughly 23–25%
of the molecular store held at 153.5 microseconds by 193 microseconds, and about
30–31% by 250 microseconds. Extreme slow-rate scenarios return less. These are
conditional amounts, not likelihood intervals or useful-work fractions.
Step 10 establishes source-applicability limits but does not validate the hot/dense
rates. No coupled finite-rate flow, distributed-injection performance, or
survivability result yet.

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
- [x] Isolate the conditional post-switch chemistry effect using same-state parcel-frozen
  controls at two switch times and two spatial resolutions (Step 5).
- [x] Separate an ordered molecular-reaction effect from atomic ionization using a
  fixed-molecule/active-ion intermediate closure (Step 6; not finite-rate kinetics).
- [x] Audit a reduced neutral network's local relaxation and rate applicability on
  the saved equilibrium flow (Step 7; kinetic validation remains open).
- [x] Check direct-water pressure-source coverage and channel removal (Step 8).
- [x] Resolve local fixed-energy feedback and association rate control (Step 9).
- [x] Trace evaluated H2/collider sources and compare a documented H2 replacement
  (Steps 10–11; hot/dense physical validation not passed).
- [ ] Establish an applicable hot/dense rate and thermodynamic uncertainty envelope.
- [x] Verify nonlinear fixed-energy neutral chemistry and integrate prescribed
  parcel density/energy-forcing histories (Step 11; ions fixed, flow not recomputed).
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

Historical Step 3 verification (superseded by the Step 9 verification below):
39 study acceptance tests pass; `make lint` passes
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
- Completed: Steps 0–9 and nonlinear conditional parcel scenarios (Step 11), with
  substantial uncommitted code, tests and evidence. Do not restart at Step 3 or
  overwrite the working tree. Current verification is recorded at Step 11.
- Physical-data gate still open: hot/dense H2 association/collider and other
  association rates (Step 10 below).
  NIST's kinetics database and Princeton's Li manuscript became accessible after
  the network retry. The original Cohen/Tsang PDFs remain unread; the NIST records
  are separately identified as database transcriptions of those evaluations.
  The Li H2 replacement and nonlinear parcel integration are now implemented;
  source agreement does not turn the missing physical envelope into a passed gate.
- Reproduce local diagnostics with `make water-plate-kinetics-audit`,
  `make water-plate-pressure-kinetics-audit`, and
  `make water-plate-thermal-kinetics-audit`. They use the saved N40/N80
  `freeze_parent_h1_*.jsonl` histories, not new injected flow.
- Reproduce the new amounts with `make water-plate-parcels` and
  `make water-plate-parcel-report`; Step 11 names grid/sampling/solver overrides.
- Next numerical gate: conservative chemistry/flow coupling, including pressure
  feedback, dense-EOS convergence and the missing ion/loss physics. The parcel
  results cannot supply a finite-rate wall impulse by postprocessing alone.
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

Verification: the full checkpoint suite passed **519 Python tests** (21m29s)
and **183 Rust release tests** (11 intentionally ignored). Subsequent additions
and final edits pass **98 focused Python tests and 9 Rust study/history tests**,
including the additional vibrational-heat-capacity test. `make lint` passes.
All three reported restart files were rerun with the final equilibrium-splice
guard; invalid-domain and mismatched comparisons have explicit regression tests.

## Step 6: separating molecular reactions from atomic ionization

An intermediate closure now holds each parcel's H2O, OH, H2 and O2 abundances
fixed while allowing its remaining atomic H/O inventory and electrons to
equilibrate through the existing Saha ion ladder. The fixed molecular bond store
does not return, but ionization can still exchange energy with the thermal pool.
The thermal degrees of freedom, dense-water residual and exact conserved restart
state remain unchanged. Details and thermodynamic acceptance tests are in
[ADR-0045](adr/0045-water-plate-selective-chemistry.md).

This compares three closures: full equilibrium E; molecules fixed, ionization
active M; and all species fixed F. The **ordered** differences are
`G_E-G_M` (allow molecular reactions with ions active) and `G_M-G_F`
(allow ionization with molecules fixed). They sum to Step 5's `G_E-G_F`.
Changing the order could allocate coupled chemistry/flow effects differently;
this is not a unique decomposition of chemical work. Forward and reverse
reactions are both included, not recombination alone.

### Conditional molecular contribution

All gains below use the same 0.5 ms endpoint and subtract the unsprayed and
analytic cold-only controls. Positive contributions increase the water-layer
advantage; a negative contribution need not mean negative thrust.

| Switch (microseconds) | N | Molecular effect, ions active (kN s) | Ionization effect, molecules fixed (kN s) | All-chemistry effect (kN s) |
| --- | --- | --- | --- | --- |
| 153.5 | 40 | 126.456 | -40.356 | 86.099 |
| 153.5 | 80 | 127.219 | -39.182 | 88.037 |
| 193 | 40 | 22.516 | -21.837 | 0.679 |
| 193 | 80 | 22.143 | -21.564 | 0.579 |

The early-switch molecular contribution is positive and material: at N=80 it
is **18.4%** of the equilibrium water-layer gain. The molecular effect changes
by 0.60% from N=40 to 80, while the late-switch molecular effect changes by
1.66%. These are limited grid checks, not model uncertainty bounds.
Halving the N=40 restart timestep changes either molecular contribution by
less than **0.21 N s** (0.00017% for the early switch); the late net
all-chemistry difference changes by 0.10 N s. This is a post-switch timestep
check only: the saved pre-switch shock history is identical, not independently
time-refined by this experiment.
The unsprayed pulse benefits much
more from enabling atomic ionization under the fixed-molecule constraint,
reducing the **relative** water-layer advantage. That is why the ordered
molecular contribution can exceed the net all-chemistry contribution without
violating energy conservation.

The late-switch result corrects a possible over-reading of Step 5: its small
net difference does **not** mean that both chemical effects are negligible.
At N=80 roughly +22.1 kN s of molecular effect and -21.6 kN s of relative
ionization effect largely cancel. Neither conditional effect is an independently
usable engineering gain or a unique recombination-energy fraction.

This strengthens the evidence that molecular chemistry is early enough to
affect wall loading in the provisional planar reference. It does not establish
the real finite-rate recovery fraction, total pre-switch chemistry benefit,
or performance of locally mixed injected water on a finite pusher plate.

At N=80 the largest actual switch p/T jump is below 4e-13 relative; the
independent Rust/Python source checks agree within 3.01e-7 in pressure and
2.89e-10 in energy. Maximum restart energy drift is 3.27e-6 of physical energy.
Matching the analytic constrained EOS to the parent table still requires up to
30.98 MJ of summed absolute constant-energy offsets (about 0.12% of physical
energy) and a maximum 0.459% pressure correction, essentially the same as
Step 5. Small conservation errors and source agreement do not validate the
provisional dense residual or remove its table-interpolation uncertainty.

### Evidence and reproduction

- [Ordered differences](../data/results/water_plate/selective_comparison.csv).
- [Molecular-switch comparisons](../data/results/water_plate/molecular_comparison.csv).
- [Branch conservation and splice diagnostics](../data/results/water_plate/molecular_runs.csv).
- [Reference provenance](../data/results/water_plate/molecular_model.json).
- [Molecular-switch impulse histories](../data/results/water_plate/molecular_comparison.png).

Use the exact Step 5 parents, then `make water-plate-molecular` (N=80 by
default), and repeat with `WATER_FREEZE_CELLS=40`. The molecular preparation
is separately available as `make water-plate-molecular-prepare`.
`WATER_MOLECULAR_INPUT` and `WATER_MOLECULAR_OUTPUT` name distinct variant files;
the endpoint and switches use `WATER_FREEZE_END` and `WATER_FREEZE_TIMES_US`.
The optional `WATER_FREEZE_TIMESTEP_FRACTION=0.125` halves the reference
restart timestep, without changing the pre-switch parent flow; give both
all-frozen and molecular-frozen variants distinct filenames when using it.
The default report includes this N=40 half-timestep pair. Reproduce it with:

```
make water-plate-molecular WATER_FREEZE_CELLS=40 WATER_FREEZE_TIMESTEP_FRACTION=0.125 WATER_MOLECULAR_INPUT=data/results/water_plate/molecular_inputs_n40_halfdt_05ms.json WATER_MOLECULAR_OUTPUT=data/results/water_plate/molecular_n40_halfdt_05ms.jsonl
make water-plate-freeze WATER_FREEZE_CELLS=40 WATER_FREEZE_TIMESTEP_FRACTION=0.125 WATER_FREEZE_INPUT=data/results/water_plate/freeze_inputs_n40_halfdt_05ms.json WATER_FREEZE_OUTPUT=data/results/water_plate/freeze_n40_halfdt_05ms.jsonl
```

`make water-plate-molecular-report` reports the molecular runs and
`make water-plate-selective-report` pairs them with the all-frozen controls.
The latter checks actual restart arrays and integration grids as well as
metadata and matched equilibrium outcomes. Override `WATER_MOLECULAR_RUNS`,
`WATER_MOLECULAR_REFERENCE` and `WATER_SELECTIVE_ALL_RUNS` for other selections;
paired lists must use the same order. Bulk JSON inputs/trajectories remain ignored.

Verification: the full checkpoint passed **546 Python tests** (22m59s) and
**186 Rust release tests** (11 intentionally ignored). Final additions and
report changes pass **123 focused Python tests and 11 Rust study/history tests**;
`make lint` and `git diff --check` pass. The final reports include both spatial
resolutions and the matched N=40 half-timestep pair, with all branches completing
the common window and passing conservation, thermodynamic-domain and splice gates.

## Step 7: reaction-rate margin and applicability audit

The rate check now includes ten reversible reactions among H, O, H2, OH, O2
and H2O, using documented forward coefficients and collision-partner efficiencies
from the [Cantera v3.2.0 GRI-Mech hydrogen–oxygen submechanism](https://github.com/Cantera/cantera/blob/v3.2.0/data/h2o2.yaml).
Reverse rates are matched to this study's molecular partition functions, not
silently substituted from the mechanism's NASA7 thermodynamics. The new
[ADR-0046](adr/0046-water-plate-kinetics-applicability.md) records the source hash,
unit/event conventions, linear-response construction and applicability limits.

This is **not** finite-rate hydrodynamics. It measures a local isothermal,
fixed-density, fixed-ion, bond-energy-weighted relaxation time from the four
reactive eigenmodes of the detailed-balanced neutral network. It includes both
directions of each reaction and multiple routes; it is not the former one-way
OH-partner collision time. HO2/H2O2 and ionic reaction/collider channels remain
absent. Equal forward and reverse equilibrium fluxes do not count as heat release.

### Conditions during modeled molecular return

Audit the existing 100–250 microsecond loading history, weighted by **positive
consecutive parcel bond-store decreases**. These weights count repeated returns
after reheating again; they are not unique released energy, parcel-peak decreases
or causal impulse weights. Rates, temperatures, pressures and measured velocity-
gradient clocks use interval right endpoints. Every quoted temperature range
below contains the central 80% of the corresponding store-decrease weight.

| Interval (microseconds) | Positive store-decrease weight (GJ) | Temperature p10–p90 (K) | Median pressure (MPa) | Median tracking Da |
| --- | --- | --- | --- | --- |
| 100–153.5 | 0.236 | 5,135–9,449 | 2,537 | 3.19e6 |
| 153.5–193 | 1.406 | 5,653–8,270 | 373 | 1.48e5 |
| 193–250 | 0.369 | 4,861–5,368 | 39.8 | 3.28e4 |

These are N=80, fine-output results. Tracking Da compares the network relaxation
time with the smaller of the measured expansion time and the sampled
mean-held-store / store-return-rate time. The latter checks rapid chemical-store
changes that a density-only clock might miss; it is not a continuous nonlinear
composition forcing calculation. Compression receives no invented expansion clock.
In the main 153.5–193 microsecond interval, the store-decrease-weighted median
density is 71 kg/m3 and the median parcel position is **2.1 cm from the wall**.
These are sampled material positions, not a separately resolved reaction front.

Across the whole interval, **99.74%** of the store-decrease weight has a resolved
expansion clock and tracking Da above 10. About **0.26%** has no positive expansion
clock and remains unclassified, not frozen. Numerically tiny decreases remain in
the denominator; no spectrum failure is relabeled a fast reaction. The tracking
Da p10/median/p90 are approximately **3.12e4 / 1.28e5 / 2.40e6** for the resolved
subset. Uniformly reducing all retained forward and reverse rates by 1000 would
still leave **99.62%** of the total weight above Da=10 **on this saved equilibrium
history**. That is a sensitivity calculation, not a validated uncertainty bound
or a recomputed wall impulse.

The weighted median of the old OH-partner time divided by the new energy-weighted
network time is about 4.17. Neither clock is a rigorous upper or lower
bound on real kinetics. Steam's relative collider abundance evolves: the effective
water-formation collider multiplier has a weighted median near 1.50, not a fixed
3.65 for the whole mixture. Almost none of the molecular-return weight is in
material with more than 1% ionized nuclei. This does not erase the separate ion
impulse effect in the hotter unsprayed configuration.

### What remains unvalidated

**All** of the positive molecular store-decrease weight in this loading interval
is above the source mechanism's **3500 K thermodynamic-fit ceiling**. This is
not a demonstrated upper validity limit for each forward rate; it shows that
the supplied mechanism cannot be used unchanged as a thermodynamic validation.
The present screen uses the study's own extended thermodynamics and extrapolates
the listed forward rate laws.

The retained association reactions, including H + OH + M -> H2O + M, have **no
high-pressure falloff law in that source**. About 24% of the fine-grid return
weight occurs above 1 GPa; roughly 17–19% has more than a 10% difference between
the flow pressure and the ideal-mixture pressure across the sampling checks.
The pressure thresholds are diagnostics, not certified rate limits. No evaluated
uncertainty envelope or dense-fluid activity/rate correction at these conditions
has been established here. Thus every audit summary explicitly retains
`kinetics_validated = false`.

The numerical rate margin is substantial, but it does not yet validate the
127 kN s conditional molecular impulse contribution from Step 6 as real finite-
rate performance. Missing pressure-dependent kinetics, nonlinear thermal feedback,
and the provisional dense EOS still matter; actual injection, finite-plate escape
and wall losses remain separate work.

### Numerical checks and evidence

The total positive store-decrease weight is 2.0018 GJ at N=40 and 2.0108 GJ at
N=80, a 0.45% change. Thinning the N=80 output from approximately 0.5 to 1.0
microsecond spacing changes it by 0.53% and the median tracking Da by about 3.4%.
The fraction retaining Da above 10 after a uniform 1000-fold rate slowdown stays
between 99.55% and 99.77% across the two spatial grids and output samplings.
Threshold-based nonideality fractions are less precisely converged; quote their
sampling spread, not fine-grid digits as physical accuracy. The unsprayed control's
trace interval decreases are only around 10–16 J and are explicitly flagged
negligible, not used to declare a molecular freeze-out result.

Evidence: [weighted rate/applicability audit](../data/results/water_plate/kinetics_audit.csv)
and [pinned coefficients and parent provenance](../data/results/water_plate/kinetics_model.json).
`make water-plate-kinetics-audit` regenerates both spatial grids and both output
samplings from the exact Step 5 parents. `WATER_KINETICS_PARENTS` and
`WATER_KINETICS_STRIDES` select alternatives. The audit requires finite, accepted
parent conservation diagnostics, fixed parcel masses and exact window boundaries.
It makes no changes to the Rust flow or to the chemistry-switch impulse artifacts.

Verification: the full checkpoint passed **569 Python tests** (23m37s) and
**186 Rust release tests** (11 intentionally ignored). Final additions pass
**138 focused Python tests and 11 Rust study/history tests**. `make lint` and
`git diff --check` pass. Tests cover units, reaction-event conventions, detailed
balance, reversible Jacobians, elemental energy-gauge invariance, rate scaling,
sampling boundaries, parent conservation/mass gates and unresolved-weight handling.

## Step 8: pressure-source coverage and removing direct water association

**The fast local chemistry result does not depend solely on the uncapped
H + OH + M <=> H2O + M reaction.** Removing this reaction in both directions
slows the local energy relaxation by a weighted median factor of **2.54**,
not orders of magnitude. This narrows one rate uncertainty; it does not
validate all the other rates or the Step 6 impulse prediction.

### A documented pressure-dependent source, with incomplete coverage

The [Cantera example-data implementation of Singal et al. (2024)](https://github.com/Cantera/cantera-example-data/blob/1a5d27e508a38b1791543e9fded80ffd5c5b8d75/ammonia-CO-H2-Alzueta-2023.yaml#L966)
supplies a pressure-dependent direct water-association curve and a
temperature-dependent steam collision efficiency. Its reference-collider
pressure table reaches **10,000 atm = 1.01325 GPa**. The final table point
is **not** an explicit high-pressure limit. The associated author-group
collider database supplies values at 300, 1000 and 2000 K; all of our
positive store-decrease weight is hotter than those points. This marker
is not a certified temperature boundary for the whole reaction model.
[ADR-0047](adr/0047-water-plate-pressure-kinetics.md) records the pinned
sources, extracted coefficients, matching rules and limitations.

Importantly, the pressure argument is the **neutral ideal collider
pressure multiplied by mixture efficiency**, not the measured dense-flow
pressure. Its weighted median during 153.5–193 microseconds is **1.265 GPa**,
although the flow-pressure median is only 373 MPa. Consequently only
**42.65%** of the main-window store-decrease weight is inside the supplied
effective-pressure table. Whole-window coverage is **48.98%** at N=80 with
fine sampling, ranging **48.98–50.30%** across both grids and samplings.
Do not silently clamp the remaining states to the highest pressure node.

On the pressure-covered subset, extrapolating the source's temperature
fits gives a weighted median direct association coefficient **7.44 times**
the baseline value. Replacing that channel alone makes the network
relaxation time a median **0.195 times** baseline. About **48.84% of the
total weight** has a supplied-pressure rate, a resolved clock and Da>10;
the other weight is not declared slow. The reference's five other neutral
colliders inherit its default N2 efficiency, accounting for a whole-window
median 79.3% of neutral particles. These defaults are not evaluated
collision efficiencies for a hot radical-rich mixture.

This is an extracted **single-channel hybrid**, not the full source
mechanism. Reverse rates still use the study's equilibrium constants, and
the other nine reactions remain unchanged. The source's separate explicit
steam dissociation channel and other pathways are not added. The forward
rate implementation agrees with independent Cantera 3.2.0 evaluations to
better than 1e-13 relative at three checked in-table states. Agreement
verifies implementation, not applicability at 5–9 kK or dense-fluid conditions.

### Does losing this one channel make chemistry too late?

Set the direct association/dissociation flux to zero, keeping the other
nine rates and the saved equilibrium history fixed. Water still forms
through **H2 + OH <=> H + H2O** and **2 OH <=> O + H2O**. The remaining
network retains all four reactive modes. With its other rates fixed,
removing this one channel is mathematically a **slow endpoint for the
isothermal local energy-relaxation time**: adding any nonnegative
detailed-balanced direct-channel flux can only shorten that time
(ADR-0047). It is not a bound on nonlinear chemistry or thrust.

| Interval (microseconds) | Median local slowing factor after removal | Median tracking Da after removal | Total store-decrease weight with Da>10 |
| --- | --- | --- | --- |
| 100–153.5 | 1.96 | 1.63e6 | 99.28% |
| 153.5–193 | 2.47 | 5.96e4 | 99.75% |
| 193–250 | 3.17 | 1.04e4 | 99.97% |

These are N=80, fine-sampling, right-endpoint results with the same
positive consecutive parcel store-decrease weights as Step 7. They are
neither unique released energy nor causal impulse weights. Whole-window
tracking Da p10/median/p90 after removal are **1.01e4 / 4.83e4 / 1.19e6**,
and **99.74%** of the total weight remains above Da=10. The unclassified
fraction is almost entirely the same missing positive expansion clocks;
no unresolved spectrum is recast as a fast rate. Numerically tiny changes
remain in the denominator.

The weighted median slowing factor stays **2.54–2.56** across N=40/80
and both samplings. The fraction above Da=10 stays **99.55–99.79%**.
If the remaining nine rates are also uniformly reduced by 1000, **90.06%**
of fine-grid weight still has Da>10 (89.77–90.43% across samplings/grids).
This additional slowdown is only an algebraic sensitivity on the saved
history, not an uncertainty distribution or a new flow calculation.

The result supports early molecular relaxation even if direct water
association is strongly suppressed. **Other retained three-body
associations still have no high-pressure limits**, and the exchange rates,
nonlinear thermal response and dense EOS remain unvalidated here. Thus
the 127 kN s conditional molecular impulse effect is still not an
established finite-rate nozzle result. The modeled return positions and
loading times from Step 7 have not been recalculated or moved by this audit.

Evidence: [pressure-source and removal audit](../data/results/water_plate/pressure_kinetics_audit.csv)
and [rate/parent provenance](../data/results/water_plate/pressure_kinetics_audit.json).
`make water-plate-pressure-kinetics-audit` regenerates both grids and samplings;
the Step 7 parent/stride overrides also apply. All 32 baseline summary rows
agree with the original Step 7 artifact, which is preserved. Tests cover
units, pressure interpolation/rejection, composition-dependent detailed
balance, the independent Cantera anchors, the removal bound and a separate
conservation-projected linear solve. The full checkpoint passes **596 Python
tests** (11m27s with eight parallel test workers) and **186 Rust release
tests** (11 intentionally ignored). The focused new suite passes 20 tests.
`make lint` and `git diff --check` pass. The ordinary debug `make test` was
started, then its slow debug run was stopped in favor of the complete release
Rust suite and parallel Python suite above; no test failures were bypassed.

## Step 9: remaining rate controls and local temperature feedback

The next check permits temperature to respond to small composition
disturbances while conserving **local total energy at fixed volume**.
It also measures which of the ten rates control the relaxation time and
slows all six association channels together. This is still a local linear
diagnostic on the saved equilibrium history, **not nonlinear finite-rate
chemistry or a new impulse prediction**. No newly evaluated physical rate
set has been substituted. [ADR-0048](adr/0048-water-plate-local-thermal-kinetics.md)
defines the constraint, observable, rate controls and acceptance tests.

### What temperature feedback means in this check

Use per-species total internal energies and the **fixed-composition** heat
capacity, including neutral internal modes, fixed ion/electron sensible
heat and the same provisional dense residual. Using equilibrium heat
capacity here would count the composition response twice. Both the
relaxation operator and the bond-energy susceptibility are transformed
for the fixed-energy constraint; changing only the eigenvalues and keeping
the isothermal energy weights would mix two different diagnostics.

The infinite-heat-capacity limit reproduces Step 7's isothermal time.
An independent finite-difference check solves temperature from conserved
energy after each composition perturbation and verifies the local kinetic
Jacobian. This is not a claim that the modeled cooling history reverses
or that a shorter local relaxation time releases more usable energy.

At N=80 with fine output sampling, the weighted median fixed-energy time
is **0.269 times** the fixed-temperature time: roughly a 3.7-fold faster
decay of these small disturbances. The corresponding median chemical
susceptibility is only **0.272 times** its isothermal value. Whole-window
median tracking Da is **4.87e5**, with 99.74% of the total store-decrease
weight above Da=10. Thus this local temperature feedback does not erase
the rate margin, but it also does not imply an extra energy or thrust gain.
The whole-window median frozen heat capacity is 2440 J/kg/K; its dense
residual contributes a median 0.154%. The equilibrium compositions and
flow pressures remain subject to the earlier, separate dense-EOS limits.

### Why the remaining association rates matter

Rate control is measured by `S_r = -d ln(tau) / d ln(k_r)`, with both
directions of a reaction scaled together and equilibrium fixed. These
nonnegative sensitivities sum to one. Their positive-store-decrease-
weighted means are **not reaction heat fractions or impulse fractions**.
In particular, a fast alternate water-forming exchange reaction can
carry flux without being the slow step that controls complete relaxation.

The grouped slowdown changes source reactions 1, 2, 12, 13, 14 and 15:
O2 formation, OH formation, all three H2-formation collider channels,
and direct H2O formation. The four exchange rates stay unchanged. A
1,000-fold or million-fold slowdown is a parameterized sensitivity, not
an evaluated uncertainty interval. Exactly deleting every association
creates an extra conserved neutral-particle count: the exchange-only
network has rank three rather than four. It cannot be labeled fully
equilibrating by clipping the resulting zero eigenvalue.

The whole-window fixed-energy rate-control means are:

| Rate group | Baseline control weight | With direct water association removed |
| --- | --- | --- |
| Direct H + OH association to H2O | 59.74% | 0% |
| H2 formation, all three collider channels | 30.70% | 76.66% |
| H + O association to OH | 7.28% | 17.89% |
| O + O association to O2 | 1.67% | 4.04% |
| All four exchange reactions | 0.62% | 1.41% |

These means condition on resolved rate pairs, covering effectively all
of the non-tiny store-decrease weight. Association channels carry **99.38%**
of baseline rate control and **98.59%** when direct water association is
absent. In the latter case, the explicit **H2O-collider H2-formation channel
alone carries 41.32%**. That makes H2 recombination with realistic collision
partners a concrete next rate-validation priority. Redundant routes explain
why direct water association can be important in the baseline sensitivity
yet its removal need not prevent rapid equilibration.

### Stress test of all association channels together

| Interval (microseconds) | Median fixed-energy Da, baseline | With associations 1,000x slower | With associations 1,000,000x slower |
| --- | --- | --- | --- |
| 100–153.5 | 8.13e6 | 8.18e3 | 8.21 |
| 153.5–193 | 5.49e5 | 552 | 0.553 |
| 193–250 | 1.62e5 | 162 | 0.163 |

At a 1,000-fold slowdown, **99.73%** of whole-window weight still has a
resolved clock and Da>10. At a million-fold slowdown that fraction is
only **5.52%**; the main and late loading windows no longer have even a
median Da above one. This extreme parameter choice would undermine the
local equilibrium-tracking screen. It is **not evidence that actual water
rates are that slow**, nor a recomputed physical freeze-out or impulse.
All of these diagnostics use the unchanged positive consecutive parcel
store decreases and right-endpoint clocks from Step 7, not unique energy
release or causal impulse weights.

Across N=40/80 and both output samplings, the median thermal/isothermal
time ratio is **0.267–0.270**. The Da>10 fraction is **99.55–99.79%** for
the 1,000-fold slowdown and **5.37–5.96%** for the million-fold slowdown.
The baseline/zero-direct/1,000-fold variants have no additional unresolved
rate weight beyond the tiny-change exclusion. The extreme million-fold
case has about **0.176%** of fine-grid weight without a resolved rate
pair, retained in the denominator and not declared fast or frozen.

The 512-to-1024-knot residual check changes the reported non-extreme
thermal time/Da quantiles by less than **8e-11 relative**, with unchanged
Da>10 fractions. The extreme million-fold case has marginal numerical
acceptance: its lower Da quantile shifts by up to **0.19%**, and per-window
resolved coverage by up to **0.020 percentage points**. Both constraints
must pass the spectrum and sensitivity-sum gates, so this affects their
common accepted subset; it is not evidence for a large physical heat-
capacity interpolation error. Do not overinterpret extreme-case tail digits.

Evidence: [thermal/rate-control audit](../data/results/water_plate/thermal_kinetics_audit.csv),
[model and parent provenance](../data/results/water_plate/thermal_kinetics_audit.json),
and the [1024-knot residual check](../data/results/water_plate/thermal_kinetics_refined.csv).
`make water-plate-thermal-kinetics-audit` runs both spatial grids and both
output samplings; the existing kinetics parent/stride overrides apply.
All 32 original baseline summaries match Step 7 exactly. Verification:
**612 full-suite Python tests** passed (11m55s, eight workers), plus the
final **53 focused rate/thermal tests**, including two tests added after
the full run collected its cases. The full Rust release suite passed
**186 tests** (11 intentionally ignored). `make lint` and `git diff --check`
pass. The full Rust/Python checks use the release/parallel equivalents of
the repository's `make test` commands; no production Rust flow was changed.

## Step 10 in progress: evaluated-rate source access and limits

Network retry on **2026-09-08** obtained the author-hosted Li manuscript and
NIST Chemical Kinetics Database records. This materially improves the source
audit, but does not validate the rates at the saved flow's dense/hot states.

### Directly inspected sources

[Li et al. (2004), author-hosted manuscript](<https://www.princeton.edu/~combust/research/publications/H2_O2%20Mech%20(Li%20et%20al%20Int%20J%20Chem%20Kin%2036%20566%20575,%202004).pdf>),
*An updated comprehensive kinetic model of hydrogen combustion*,
[DOI 10.1002/kin.20026](https://doi.org/10.1002/kin.20026):

- PDF page 9 reports mechanism comparisons over **298–3000 K, 0.3–87 atm**
  (upper pressure 8.82 MPa). This is the paper's collection of experimental
  conditions, not proof of validation at every point of that rectangle or a
  reaction-specific steam uncertainty envelope.
- Table I, PDF page 13, gives H2 + M -> 2 H + M with
  **A = 4.58e19, n = -1.40, E = 104.38 kcal/mol**, in cm/mol/s units;
  footnote a gives H2O efficiency **12.0**, H2 efficiency **2.5**.
  No high-pressure limit is supplied for this reaction.
- The 300–3000 K collider-efficiency averaging statement on PDF page 6 concerns
  **H + O2 (+M) -> HO2 (+M)**, not H2 formation. It must not be reassigned.
- PDF page 7 explains that direct H + OH + M water formation was adjusted to
  improve flame speeds, and that transport/rate uncertainty permits nonunique
  fits. Whole-mechanism agreement is not independent validation of each rate.
- Download SHA-256:
  `4749d34d818dafbb4fa0feb5f76f2e908ae7f2158b21e6be4b526618e803da3b`.

The following are **NIST database records**, not a claim to have read the
underlying Cohen–Westberg or Tsang–Hampson articles. The displayed rate convention
is `k(T) = A*(T/298 K)^n*exp(-Ea/RT)`, with **molecules**, not moles.

| NIST record | Reaction / bath | Listed temperature range | Displayed A, n, Ea (J/mol) | Reported uncertainty |
| --- | --- | --- | --- | --- |
| [1986TSA/HAM1087:91](https://kinetics.nist.gov/kinetics/Detail?id=1986TSA%2FHAM1087%3A91) | H2 dissociation / N2 | 600–2000 K | 2.61e-8 cm3/molecule/s, -1.40, 436510 | 3.0 |
| [1983COH/WES531:40](https://kinetics.nist.gov/kinetics/Detail?id=1983COH%2FWES531%3A40) | 2 H association / H2O | 300–2000 K | 9.26e-32 cm6/molecule2/s, -1.00, 0 | 5.0 |
| [1983COH/WES531:42](https://kinetics.nist.gov/kinetics/Detail?id=1983COH%2FWES531%3A42) | 2 H association / H2 | 50–5000 K | 9.04e-33 cm6/molecule2/s, -0.60, 0 | about 2.51 |
| [1983COH/WES531:18](https://kinetics.nist.gov/kinetics/Detail?id=1983COH%2FWES531%3A18) | H2 dissociation / Ar | 600–5000 K | 1.88e-8 cm3/molecule/s, -1.10, 436510 | 2.0 |

These detail records do not provide a pressure-validity interval or dense-fluid
correction. In particular the factor 3 inherited by the FFCM2 trial table now has
an identifiable **600–2000 K nitrogen-bath** evaluation behind it. It is not a
5–9 kK steam uncertainty bound. The older steam-specific record does supply an
uncertainty, unlike FFCM2's `N/A` for its adopted steam efficiency; these are
different prescriptions and neither validates the hot/dense extrapolation.
The Ar record cannot validate a steam collider, and Ar is absent from this study.

For the source-coefficient check, the NIST N2 expression converts to approximately
`4.58e19*T^-1.4*exp(-436510/RT) cm3/mol/s` by multiplying A by
`N_A*298^1.4`. The result is `4.574e19`, within 0.2% of the Li executable
coefficient using the rounded database inputs. Thus the **4.58e16** prefactor in the
[pinned FFCM2 trial web table](https://github.com/dongwd016/FFCM2_Website/blob/64c78825158ac2b240c2ea32e911ed7e9cf9973f/assets/data/trialmodel/trialmodel.json)
is not used. The cause of that web-table discrepancy remains unresolved; it is
not silently interpreted as a unit conversion or a physical uncertainty range.

Crossref verified the original evaluation identifiers:
[Cohen and Westberg (1983), 10.1063/1.555692](https://doi.org/10.1063/1.555692)
and [Tsang and Hampson (1986), 10.1063/1.555759](https://doi.org/10.1063/1.555759).
Their publisher PDFs returned HTTP 403; NIST's PDF archive returned 503.
Wayback availability queries returned no snapshots for those exact current
publisher PDF URLs or the Li manuscript URL. CDX queries also encountered
503/timeouts. This is a retrieval outcome, **not evidence that no archived
copies exist**. Princeton supplied Li directly; OpenAlex/Semantic Scholar
did not identify a usable open Cohen/Tsang full-text copy in the checked records.

### Consequence and exit criteria

Step 7's main loading interval has a 10th–90th percentile temperature range
of **5653–8270 K** and median flow pressure **373 MPa**. These new source
limits do not support treating its large extrapolated Da as physically validated.
The early/near-wall equilibrium finding survives as a conditional result; no
finite-rate impulse result is added by this source retrieval.

The conditional Li H2 replacement is implemented in Step 11 without double-counting
the three original H2 channels. It retains study-matched reverse rates and explicitly
extrapolated status. A useful rate-validation conclusion requires either an
applicable temperature/collider/pressure uncertainty envelope or an explicit
unresolved physical-data gate; numerical agreement among inherited mechanisms
does not pass that gate. Nonlinear parcel tests may proceed as conditional
numerical work, not as a substitute for missing physical rate data.

The remaining performance gates are: (1) conservative finite-rate flow and EOS
convergence; (2) actual local injection plus finite-plate escape with matched
combined/unsprayed/injection-only controls; (3) losses, material-specific
protection, permanent assembly mass, and matched nozzle comparison. Each needs
analytic/conservation checks and resolution tests where applicable. A missing
material or rate input must be named, not assigned a fictitious validated value.

## Step 11: nonlinear conditional amounts during the loading window

This is the first **nonlinear amount calculation**, not another conversion of a
local Da into a presumed release fraction. The Rust integrator evolves the four
independent molecular populations and temperature together. Hydrogen and oxygen
budgets determine the two atomic populations. Reverse rates and species energy
use the study's common partition functions. The algorithm and limitations are
specified in [ADR-0049](adr/0049-water-plate-nonlinear-parcels.md).

### Ensemble and energy ledger

All **237.5 kg** of the N40/N80 combined planar reference start from their exact
153.5 microsecond states. There are 380/760 fixed-mass parcels; no spatial
representative-subset weighting is substituted. At N80 the initial molecular
dissociation store is **5.627720 GJ**. For each endpoint the amount is

```
net return = sum_parcels mass * (bond_store_at_153.5us - bond_store_at_endpoint).
```

Negative parcel returns are retained. Successive positive decreases are **not**
added; reheating cannot earn the same bond energy twice. These fractions are
relative to the remaining store at the stated start, not incoming kinetic energy,
the maximum full-atomization capacity, water-molecule fraction, or total release
since impact. The 193/250 microsecond endpoints refer to the existing reference
loading windows, not independently solved pressure-decoupling times.

Density is replayed log-linearly between parent samples. Temperature is solved,
not replayed. The implicit update obeys
`delta(e) = -p_new*delta(1/rho) + prescribed_input`, with the existing dense
Helmholtz residual. Two forcing cases are compared:

- Reversible expansion/compression work only, with no external input.
- The same work plus the parent's sampled non-mechanical residual,
  `delta(e_parent) + mean(p_parent)*delta(1/rho_parent)`.

The latter includes shock heating **and** sampling/EOS/discretization defects;
it is not measured chemical heating. At 193 microseconds its ensemble input is
about **0.156/0.198 GJ** with fine/coarse parent sampling. Its sensitivity is kept
visible instead of silently forcing the parcel to follow the equilibrium temperature.
The initial energy offset includes constant frozen-ion chemical energy and the
elemental/storage energy gauge as well as matching the parent table; do not
interpret the full offset as interpolation error.

Ions/electrons are fixed populations per unit mass, but carry sensible heat.
Initially the mass-mean ionized-nuclei fraction is about **10.12%**, concentrated
largely in the incoming hot material; 10.66% of mass has more than 1% ionized
nuclei. This is not a complete ion/molecular recombination calculation.

### Conditional result

The rounded ranges below include the tested grid, history-sampling, solver and
prescribed-forcing variations. They are **not statistical confidence intervals**
or physical bounds on untested rate laws. The first row includes baseline rates,
the documented Li H2 replacement, and removal of both direct water association
and the explicit steam-assisted H2 channel. Their differences are smaller than
the checked numerical/forcing spread; that is not independent rate validation.

| Rate scenario | Net return by 193 microseconds | Net return by 250 microseconds |
| --- | --- | --- |
| Documented-rate alternatives / named channel removals | about **1.3–1.4 GJ**, **23–25%** | about **1.7–1.8 GJ**, **30–31%** |
| All six association channels slowed 1,000-fold, exchanges unchanged | about **1.3–1.4 GJ**, **23–25%** | about **1.7–1.8 GJ**, **30–31%** |
| All six association channels slowed 1,000,000-fold, exchanges unchanged | about **0.55–0.59 GJ**, **10%** | about **0.67–0.71 GJ**, **12%** |
| All chemistry frozen, numerical control | zero to roundoff | zero to roundoff |

The two slowdown factors are deliberate stress parameters inherited from Step 9,
not measured uncertainty factors. In particular **10% is not a defensible physical
lower bound**. Slower untested chemistry, different flow, losses, additional
species or different dense-water thermodynamics can move the result outside this
scenario set. Conversely this test does not assert that molecular return ceases
after 250 microseconds; later return may be less useful to the plate.

The saved parent equilibrium trajectory gives net returns of **1.369796 GJ**
and **1.735269 GJ** at the same two endpoints. That comparison is informative,
but not a convergence target for fixed-ion, independently energy-solved parcels.
The reported amounts are not wall work, and none has been multiplied by an
invented conversion efficiency to produce thrust.

### Numerical checks and retained evidence

Kernel acceptance checks cover analytic reversible chemistry, the frozen ideal-gas
adiabat, first-order temporal convergence, independent Python/Rust partition
functions and source terms, analytic composition/temperature Jacobians, finite
water perturbations returning to a known fixed-energy equilibrium, equilibrium
stationarity, positivity, elemental conservation and the energy ledger. Invalid
rates/settings are rejected. Missing endpoints or failed energy checks stay in
the ensemble denominator in the report; signed-return/missing-coverage tests
exercise that contract.

The initial N40/N80 runs use every second saved output, relative tolerance 1e-4,
and maximum step 0.25 microseconds. Separate N80 checks use (a) tolerance 1e-5
and maximum step 0.125 microseconds, and (b) every saved output with the original
solver settings. Among reacting cases, maximum relative net-return changes are
**2.17%** for N40/N80, **0.73%** for the solver refinement, and **2.02%** for
history sampling. Thus fine digits are not warranted. A combined finer check
uses every saved output, tolerance 1e-6 and maximum step 0.0625 microseconds for
baseline and the two slowdown cases. Its results are retained separately.
Relative to the fine-history/original-solver run, this last solver check changes
net return by at most **0.93%**. Its baseline values are **1.361–1.380 GJ**
at 193 microseconds and **1.725–1.743 GJ** at 250 microseconds across the two
forcing prescriptions. The rounded scenario table deliberately retains the
wider observed numerical spread rather than presenting these as precise bounds.

Every tested parcel reaches both endpoints: **100% mass/start-store coverage**
to summation roundoff, with no inferred contribution assigned to unresolved
parcels. The largest accumulated energy-ledger error across all runs is less
than **1.14e-8 of the parcel's initial internal energy**. The frozen control's
ensemble bond return is below **3e-9 J** in magnitude. Conservation is much
tighter than the trajectory discretization uncertainty; the two are not the
same accuracy test.

Evidence: [net-return scenarios and convergence](../data/results/water_plate/parcel_return.csv)
and [source, input and output provenance](../data/results/water_plate/parcel_return.json).
Large replay inputs/JSONL remain ignored but are hashed in the provenance. The
independent thermochemistry fixture is retained under
`crates/water_plate/tests/fixtures/neutral_parcel.json` and checked against its
Python generator. `make water-plate-parcel-check` runs the focused acceptance tests.

Reproduce each base run with `make water-plate-parcels`, overriding
`WATER_PARCEL_PARENT`, `WATER_PARCEL_INPUT`, and `WATER_PARCEL_OUTPUT` for N40/N80.
Use `WATER_PARCEL_STRIDE=1` for the history check; use
`WATER_PARCEL_TOLERANCE=1e-5 WATER_PARCEL_MAX_STEP=1.25e-7` for the separate solver
check. The combined fine input is generated with `python -m
puffsat.water_plate.parcel --stride 1 --tolerance 1e-6 --maximum-step 6.25e-8
--variants baseline association_1e3_slower association_1e6_slower`, with explicit
`--parent` and `--output`, under `PYTHONPATH=python uv run --extra sci`.
`make water-plate-parcel-report` combines all five outputs.

Verification: **614 full-suite Python tests** passed (11m02s, eight workers),
followed by **66 focused parcel/rate/thermal tests**, including the 13 new parcel
tests added after full-suite collection. The final full Rust release suite passed
**193 tests**, with 11 intentionally ignored. The numerical results remain
conditional despite those software checks.
`make water-plate-parcel-check`, `make lint`, and `git diff --check` pass.

### Useful conclusion at this gate

For this reference density history, documented extrapolated rates permit an
appreciable amount of the remaining molecular store to return during the principal
loading interval. Thousand-fold association slowing barely changes that broad
amount; million-fold slowing reduces it substantially but does not eliminate it.
This supports pursuing conservative chemistry/flow coupling. It does **not**
establish a likely real-world range, a finite-rate thrust increment, or the
survivability/superiority of a sprayed-water plate.

## Final cooling check and bounded conclusion (2026-09-08)

The high-pressure objection is well founded: the nominal calculations do **not**
show molecular recombination being substantially prevented by slow reactions.
The fine neutral parcel replay returns about 99.4% of the parent equilibrium
net molecular return at both 193 and 250 microseconds. These are not identical
thermodynamic controls (the replay freezes ions), but they rule out describing
the nominal result as a large neutral kinetic shortfall on this history.

The remaining store is principally a temperature/history issue. The reference
is stratified, without the proposed local injection/mixing. Its incoming 25 kg
remains much hotter than its 212.5 kg water layer. Of the 5.628 GJ molecular
store at the common 153.5 microsecond start, 4.354 GJ is in the water layer and
1.274 GJ in the incoming material. The incoming material returns less than
0.001 MJ net through 1 ms in the equilibrium parent. A conventional H2/O2
rocket's cooler combustion products and nozzle expansion are therefore not a
matched thermodynamic trajectory, even if its density is lower. High pressure
helps association but does not force complete molecular binding at high T.

### Cooling versus loading through the saved tail

Recomputed from all fixed-mass cells of `freeze_parent_h1_n80.jsonl`, using
signed molecular-store differences from 153.5 microseconds, not summed positive
release or parcel peak drops:

| Time (microseconds) | Water mass-mean T (K) | Equilibrium net molecular return (GJ) | Cumulative wall impulse (MN s) | Fraction of impulse accumulated by 1 ms |
|---:|---:|---:|---:|---:|
| 153.5 | 8830 | 0 | 1.5604 | 59.23% |
| 193 | 5457 | 1.370 | 2.3778 | 90.26% |
| 250 | 4858 | 1.735 | 2.4737 | 93.90% |
| 481.437 | 4205 | 2.200 | 2.5686 | 97.50% |
| 1000 | 3800 | 2.527 | 2.6344 | 100% |

The last two molecular amounts are **equilibrium parent diagnostics**, not an
extension of the finite-rate replay beyond 250 microseconds. The N40 parent
gives 1.350, 1.716, 2.181 and 2.507 GJ at the four nonzero-return checkpoints;
the largest N40/N80 change is about 1.5% of the N80 return. Water-only return
is about 40% of its starting molecular store by 250 microseconds and 58% by
1 ms; the corresponding whole-ensemble fractions are about 31% and 45%.

Between 250 microseconds and 1 ms, another 0.792 GJ returns while wall impulse
increases by 0.1607 MN s, 6.10% of its 1 ms value. This establishes temporal
overlap and diminishing loading, **not** a causal energy-to-impulse efficiency.
At 1 ms wall pressure remains 1.85 MPa: the impulse is not demonstrably settled,
and this table does not identify a physical finite-plate decoupling time.

### Source-EOS consistency check

At each of the five N80 checkpoints, invert the underlying source EOS at each
cell's saved density and physical internal energy, then compare temperature and
pressure with the saved interpolated fields. The correct inversion is
`pressure_energy(rho,T)[1] - ENERGY_SHIFT = cell.energy_j_kg`: the table includes
the common 4 MJ/kg storage offset whereas snapshot energy excludes it.
Brent root solves used temperature tolerance 1e-5 K, with all 760 cells checked
per checkpoint (3800 checks). Maximum absolute relative discrepancies over this
sample are **0.489% in T** and **0.522% in pressure**. All inversions completed.
An initial diagnostic omitted the offset and was discarded; its larger
discrepancies are not evidence of an EOS defect.

This sampled check finds no large interpolation temperature error explaining
the incomplete return. It is not a table-refined flow integration, an audit of
every intervening shock state, or independent physical validation of the dense
EOS. The previously small n=2/n=3 impulse sensitivity is likewise only a limited
closure test. No new flow or chemistry evolution was run for this final check.

Parent SHA-256 identifiers:

- N80: `8ae0d4d1f74d8d8dfce817d57ddc45a469c7184cdc69851f437eaf18a25b9991`
- N40: `8d44881cc663c02a6e868314c5d12d53f8800863c40209e689385bde26ddf1be`

### Study conclusion

**Water recombination is not ruled out as an early contributor to pusher-plate
impulse.** On this lossless planar reference, finite-rate neutral chemistry
returns roughly 1.3–1.4 GJ by 193 microseconds and 1.7 GJ by 250 microseconds,
during appreciable wall loading. The ordered same-state closure test also gives
a positive molecular impulse increment (127 kN s through 0.5 ms), conditional
on its equilibrium-versus-frozen closures. Neither is a validated finite-rate
thrust prediction for the proposed sprayed plate.

**A useful chemistry feasibility conclusion is reached; an engineering
performance conclusion is not.** The remaining decisive uncertainties are the
physical hot/dense EOS and rate applicability, actual local mixing and finite
plate escape with chemistry/flow feedback, and radiation, wall heat transfer,
ablation and structural survival. More arbitrary rate multipliers cannot resolve
those uncertainties. Do not quote a likely real-world recombination range,
claimed nozzle superiority, or a survivable plate design from this reference.
Further work would be a new physical-model/design stage, not another necessary
iteration of the present nominal-rate screen.

## Deferred work if the design study resumes

Prioritize evaluated **H2 association and collider rates**, including
steam as collision partner, alongside direct water and OH association.
The local sensitivity work identifies where rate validation matters; it
does not supply a high-temperature/dense-fluid uncertainty envelope.
The pressure-dependent water curve found in Step 8 still leaves about
half the weighted states beyond its effective-pressure table, with
temperature applicability unestablished. Step 11 now supplies nonlinear,
energy-conserving neutral parcel amounts on prescribed histories. Next check
**chemistry/flow pressure feedback**, rather than treating that replay as a new
hydrodynamic solution. Match additional
species and reverse rates to the same energy reference. Do not treat a
standard mechanism reactor, pressure-table endpoint or large extrapolated
local Da as performance validation.
Continue EOS table convergence
and the dense-water closure limits kept visible. The all-reaction and ordered
molecular/ionization switches above are conditional tests, not a full chemistry
attribution or measured reaction rates. Then
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
