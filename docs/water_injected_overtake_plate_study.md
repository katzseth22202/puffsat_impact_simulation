# Water-injected overtake plate study

Status: scope agreed in discussion; Step 0 analytic ledger implemented. No injected
plate hydrodynamic performance, chemistry recovery, or survivability result yet.

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
- [ ] Specify bounded incoming pulse and injection profiles and source-energy accounting.
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
| Water injection beats a bare plate or a magnetic nozzle, recovers recombination energy, or simplifies survivable hardware. | No simulation/engineering evidence yet. | No result to transfer; corresponding checklist items remain open. |

Checkpoint verification: 16 new acceptance tests pass; all Python lint, format,
and strict type checks pass. A broader Python run was interrupted after 164
passes. Full repository tests/Rust lint could not execute because this environment
has no `cargo`; this is not a full-suite pass.

## Restart checkpoint

- Branch: `main`; find the latest study checkpoint with
  `git log --oneline --grep='water.*plate'` and inspect `git status --short` first.
- Completed: Step 0 and the checked items above. Reproduce with
  `make water-plate-ledger && make water-plate-test`.
- Active next step: define and verify bounded incoming-pulse and injection-source
  profiles for the first flow comparison. No hydrodynamic result exists yet.
- Scratch: `todos/water_plate_study_scope.md` points here; it is not the canonical
  task state. The companion checkout is disposable and pinned at 9440b27.
- Keep this checkpoint and the progress checklist synchronized with completed
  work. Record paper-ready findings with their evidence and limits in the table.

## Next calculation

Specify bounded incoming and injected gas profiles, then run the first controlled
hydrodynamic comparison. Existing `euler2d` is lossless and fixed-gamma; it cannot
resolve the chemistry advantage. Injection, geometry, EOS, and wall losses must
be checked together before a performance claim. Report required coating thickness,
local pressure/heat flux, and permanent structure alongside consumable economy.
The original ADR-0014/0015/0021/0024/0026 results inform the checks but their thin
ablator, temperature, and incidence assumptions do not automatically transfer.
