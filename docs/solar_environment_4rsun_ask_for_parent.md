# Asks for the companion repos: the near-Sun environment, and the charge pusher's growth cycle

From `katzseth22202/puffsat_impact_simulation`, 2026-10-15. This is written to be copied into
`katzseth22202/Balloon-Pulse-Propulsion` and worked there, as `docs/solar_environment_ask_from_impact_sim.md`.
It is self-contained: nothing in it needs this repository's code.

Two groups:
- E1-E6 (the near-Sun environment) are for the parent.
- G1-G10 (the charge pusher's growth cycle, Results 10-12 and ADR-0060 to ADR-0063) are split
  between the parent and the calcs repo (`aim_is_all_you_need`). Each says who owes it. The
  calcs-repo items G1-G4 and G10 are gathered into one re-evaluation, D1-D7, in
  `solar_dive_cycle_ask_for_calcs.md` (2026-10-16).

## Why it is asked

The charge pusher boosts its Sun ships at a solar periapsis.
- **Baseline:** Parker Solar Probe's 9.86 R_sun, which has flown.
- **Upgrade:** the paper's 4 R_sun. It is the route to a 150 km/s Earth return, and to a shorter
  boost (13.6 km/s instead of 20.8 km/s for a 100 km/s crossing).
- **Context:** this repository's Results 10-12 and ADR-0060/0062. Seth wants to know whether
  4 R_sun is materially harder than 9.86 R_sun beyond its steady heat, before the upgrade is
  argued to investors.

This repository has only rough, unsourced estimates. They are not recorded as results.
- **Heat:** steady sunlight (~3.9 MW/m^2, ~2890 suns) dominates the heat by orders of magnitude.
  The solar wind gives ~1e-5 of it, a dense fast CME ~5e-4, and a large flare's light ~1e-4.
- **Push:** a CME's pressure on a ~300 m^2 plate is ~1 N, ~0.1 m/s on a 100 t ship over a
  pass, which is negligible.
- **Open hazard:** energetic-particle (SEP) dose, which may scale as 1/r^2 to 1/r^3 toward the
  Sun.

## The asks

**E1. Coronal structure at 4 R_sun across the solar cycle.** Are closed loops, helmet-streamer
cusps and pseudo-streamers below 4 R_sun at solar maximum, or can a trajectory meet closed field
there? How often would a pass cross the streamer belt or the heliospheric current sheet, and with
what density and field contrast? Cite observations or models (PFSS source-surface heights, coronagraph
streamer statistics, Parker encounter data).

**E2. CME encounters.** At solar maximum, what is the chance that a pass of a few hours at 4 R_sun
meets a CME? What are the typical and extreme density, speed and field inside one there? Confirm
or correct the estimates above that the push and the heat are negligible beside sunlight. Parker's
reported CME passages are the anchor.

**E3. Energetic-particle dose (the one most likely to matter).** What proton and heavy-ion
fluence and total ionizing dose would a large SEP event deliver at 4 R_sun? Give it behind a few
light shieldings (e.g. 1, 3 and 10 mm Al-equivalent) over a periapsis pass of a few hours, against
the same event at 1 AU and at 9.86 R_sun. Include the radial scaling used and its source, and the
rate of such events at solar maximum. This sets the guidance electronics' hardening on every Sun
ship and every Sun-wave projectile.

**E4. Flares.** Peak irradiance and duration at 4 R_sun from an X10-class flare, including hard
X-rays and the particle burst, as a fraction of the steady load.

**E5. Dust.** The dust density and size distribution near 4 R_sun (the dust-free or sublimation
zone), and the impact rate at ~300 km/s relative speed on a ~300 m^2 plate and on a small
projectile.

**E6. Charging and plasma interaction.** Spacecraft charging levels at 4 R_sun (the photoelectron
sheath, SEP-driven deep dielectric charging), and anything that would disturb a Kapton sheet or a
polyethylene charge held in the ship's shadow.

## What would close it

A short table for 4 R_sun against 9.86 R_sun, at solar minimum and maximum:
- the steady heat, CME push and heat, flare heat, SEP dose behind shielding, dust impacts and
  charging;
- each as a number, with its source;
- a one-line verdict: "the steady heat is the design load", or the item that overtakes it.

Return the answer to this repository as `docs/solar_environment_answers_for_impact_sim.md` in the
parent, and this repository will fold any dose or heat figure into the Sun ship and Sun-wave
projectile budgets.

---

## The charge pusher's growth cycle (G1-G10)

The context, for a reader in the parent:
- A **Sun wave** of guided projectiles (graphite disks, $1000/kg) boosts **Sun ships** at the
  solar periapsis. The ships' cheap polyethylene charges are the reaction mass.
- Each Sun ship's whole remainder returns at ~100 km/s as the **Earth wave**: Kapton sheets that
  push the next departure.
- On the 3S Jovian cycle (Earth -> Jupiter -> Sun -> Earth) the departure burn's thrust meets the
  returning stream at ~101 degrees. So the plate is a flat **canted plate** beside which the
  merged fireball explodes as it slides past.
- With a square polyethylene charge, this repository's growth optimizer (`growth.py`, ADR-0063)
  finds 7-8x per 3S at 9.86 R_sun and 9-11x at 4 R_sun. Those values count Sun-wave and
  Earth-wave impactor alike, with polyethylene's carbon left atomic. Its inputs from the calcs
  repo are copied, not computed here.

**G1 (calcs repo): confirm the 3S closures this repository copied.** From commit 868e88a,
`solve_synodic_closure(3, dive_solar_radii=R, return_excess=X)` and `departure_nozzle_ledger(...,
30.0, 0.6)` gave:
- at 9.86 R_sun, X = 90 km/s: a 10.43 km/s departure excess aimed -13.5 degrees, a 4.15 km/s burn
  from the 10-day parking orbit, the stream at 100-102 degrees to the thrust axis and 103 km/s,
  and a 20.45 km/s periapsis boost;
- at 4 R_sun, X = 150 km/s: 10.54 km/s, 4.23 km/s, 94-96 degrees at 158 km/s, and 35.04 km/s.

Please confirm them against the current ledger, and say whether X = 90 km/s at 9.86 R_sun is the
right stand-in for a 100 km/s heliocentric Earth crossing.

**G2 (calcs repo): replace the departure's impulse law with the canted plate's.** The calcs repo
prices the departure with a magnetic nozzle's impact-angle law. This repository's flat canted plate
gets 29-37% of the canted ceiling `m w (cos theta + sqrt(1 + k - sin^2 theta))` with a square
polyethylene charge, and up to 43% for an ideal point fireball (`canted_2d.jsonl`). Re-run the cycle
growth ledger with that law, and with the Sun stage priced from this repository's
`charge_ledger.csv` scaled to the dive's closing speed. Compare against this repository's 7-11x
per 3S.

**G3 (calcs repo): charge the departure from the pad.** The 4.15 km/s burn starts at the 10-day
parking orbit's periapsis. The calcs repo's own `split_push_ledger` found growth falls ~6.6x when the
departure is charged from the pad. Price the canted plate's growth that way too.

**G4 (calcs repo): the 5.3 km/s 3S outbound.** `growth_cost_for_parent.md` quotes ~5.3 km/s for a
3S outbound on the retrograde-return cycle. Confirm the departure excess it implies (~12 km/s from
near-escape, by `sqrt(v_inf^2 + v_esc^2) - v_esc`), so the paper can say why the solar-dive cycle's
4.15 km/s is cheaper: Jupiter only has to stop the ship, not reverse it.

**G5 (parent): the Sun ship's plate as impactor.** The breeding factor assumes the Sun ship releases
its whole remainder as Earth-wave sheets. Its plate (the δ/D 0.5 paraboloid, ~0.16 of its final
mass) cannot obviously be released. Is it cut up, released whole as an impactor, or kept? Keeping
it costs ~16% of growth (ADR-0062).

**G6 (parent): the Earth-wave sheet's guidance.** A 35-55 µm Kapton sheet at 0.05-0.07 kg/m^2 is a
solar sail. It drifts 30,000-100,000 km on its 0.5-1 AU coast (`overtake_sheet.csv`). It must
arrive face-on within ~1 degree and tens of cm. Size the spin, attitude control and tracking, and
the guidance mass per sheet.

**G7 (parent): how the Sun-wave projectiles reach their orbit.** The cost basis is $1000/kg to
manufacture. Delivering a projectile onto a retrograde orbit with a 9.86 or 4 R_sun periapsis is
likely to dominate its real cost. State how they are placed and what that costs per kg.

**G8 (parent): the Sun-wave projectile at 4 R_sun.** There the closing speed is ~617 km/s and the
sunlight ~2890 suns. Say whether the 0.5 m graphite disk, its guidance and its shade survive the
approach, and what peak plate pressure limits the projectile per pulse (this repository scales it
as energy over plate size cubed).

**G9 (parent): the charge material.** Polyethylene's 84 MJ/kg of bonds caps k at the Earth wave
(the gas has only `w^2/2/(1+k)` per kg). Whether the carbon condenses before the plate moves the
optimum k by 2-3x (`growth_cycle.csv`, "bonds returned" rows). A bond-free charge, such as argon
ice or a metal, would remove the cap. Is there a reason the ship must carry polyethylene, or can
the Earth-departure charge be argon, as in the calcs repo's spray plate?

**G10 (calcs repo): split the impactor stream from the ship at Jupiter.** The Sun-wave impactors
leave Earth with the ship. At Jupiter the ship takes the prograde dive (arrival excess 11.70 km/s
against its 11.33 km/s floor, so it rides unpowered). The impactors take the mirror-image
retrograde dive, so both reach periapsis together, head-on. A flyby cannot change the size of the
excess, so the impactors' carrier must burn at perijove to raise it to the 14.78 km/s retrograde
floor (12.17 -> 14.16 km/s at 4 R_sun). This repository prices the burn with Oberth alone, without
bend or timing: 0.70 km/s at 1.1 R_J and 1.37 km/s at 5 R_J, ≲1% of the ship in propellant
(`mission_jupiter_split.csv`). With the flyby code (`opposing_stream_placement`,
`solve_synodic_closure(3, ...)`):
1. Check that the flyby at 14.78 km/s still turns the impactors onto the retrograde dive at each
   perijove from 1.1 to 10 R_J. The ship's own flyby has 135° of bend available against 77° needed;
   the impactors' is not checked.
2. Solve the two departures from Jupiter as true mirror images that meet at the dive periapsis, and
   give the timing tolerance on the perijove burn.
3. Give the lowest perijove the disks' guidance electronics survive in Jupiter's radiation belts.
   Above it the burn roughly doubles.
4. Compare with the alternative: a separate impactor departure from Earth at ~11.3 km/s excess,
   ~0.6 km/s more than the ship's burn, with its own Jupiter timing.
