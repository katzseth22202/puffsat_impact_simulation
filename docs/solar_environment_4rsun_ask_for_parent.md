# Ask for the parent: the near-Sun environment at 4 solar radii during solar maximum

From `katzseth22202/puffsat_impact_simulation`, 2026-10-15. This is written to be copied into
`katzseth22202/Balloon-Pulse-Propulsion` and worked there, as `docs/solar_environment_ask_from_impact_sim.md`.
It is self-contained: nothing in it needs this repository's code.

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
