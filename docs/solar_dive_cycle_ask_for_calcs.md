# Ask for the calcs repo: re-evaluate growth on the solar-dive charge-pusher cycle

> **Already answered, 2026-10-16.** The calcs repo's ADR 0072 (`aim_is_all_you_need` 127162f,
> 2026-10-15) priced this cycle with plates, dry mass, the carrier, launch, the seed and the
> child's comparison before this ask was sent: 3.60x per cycle at 9.86 R_sun and 4.92x at 4 R_sun,
> dear-seed break-even $104 / $67 per kg. Its corrections are taken into ADR-0063's amendment of
> 2026-10-16. This ask is closed; its 4.15 km/s, 14.78 km/s and 0.16 below are superseded.

From `katzseth22202/puffsat_impact_simulation`, 2026-10-16. This is written to be copied into
`katzseth22202/aim_is_all_you_need` and worked there. It is self-contained: nothing in it needs
this repository's code. The asks are D1-D7. They gather, and take over, the calcs-repo items
G1-G4 and G10 of `solar_environment_4rsun_ask_for_parent.md`.

## Why it is asked

The child paper (`puffsats_for_datacenters`, `templateArxiv.tex` §`seed_finance`) flies the
Earth-Jupiter-Earth cycle with a pulsed-chamber departure:
- **Doubling:** 2.0-2.8 yr behind the spray cup, and 7.9 yr for the plate-only methalox cycle.
- **Seed:** $337-$9293 per kg sent toward Jupiter. The seed's cost per fleet kilogram falls only
  as one over the growth multiple: 10.5x in ten years for methane, 2.3x for methalox.
- **Largest term:** the growth charge (the `1/G` of every return reinvested).

So the seed and the growth charge, the two costs that dominate, are both set by the doubling time.
ADR 0068's shared split makes it slower still (methane 4.46 yr, hydrogen 2.96 yr). Your ADR 0063
found the departure is where growth is lost. Its chamber gives only ~7 km/s of effective exhaust
against a 7.3 km/s burn, and it expends the chamber.

This repository has a different architecture that doubles in ~1-1.2 years in its own model. The
author wants it re-evaluated with the ledger's realistic assumptions before the child decides
whether to lead with it. Beyond cost, it is likely to appeal more to the space community:
- it sets up streams of very high-energy impactors crossing Earth's orbit at 100-150 km/s;
- those streams open fast travel to the outer solar system;
- it leads into the parent paper's more exotic ideas.

## The architecture (Results 10-12, ADR-0060 to ADR-0063 here)

1. **Departure.** A ship in the 10-day parking orbit burns 4.15 km/s (9.86 R_sun closure) onto
   the 3S cycle: Earth -> Jupiter -> solar dive -> Earth, 3.28 yr. The **Earth wave**, a stream
   of Kapton sheets returning at ~103 km/s, pushes it. The stream meets the thrust axis at
   ~101 degrees, so it sweeps past a **flat canted plate**. Each sheet strikes a square
   polyethylene **charge** held 1 m off the plate, and the merged fireball explodes beside it.
   - No chamber, rods, plug or pitch are expended. The reaction mass is polyethylene.
   - The burn is small because Jupiter only has to stop the ship, not reverse it.
2. **Jupiter.** Unpowered for the ship, which falls prograde into the dive. Its **Sun-wave**
   projectiles (graphite disks, ~1/60 of the ship) travel with it from Earth. At perijove a
   carrier bus burns to put them on the mirror-image retrograde dive, so they meet the ship
   head-on at periapsis.
3. **Solar periapsis (9.86 R_sun).** The disks strike polyethylene charges in a δ/D 0.5
   paraboloid at ~392 km/s closing, head-on. This stage is efficient: ~3000 s.
4. **Return.** The Sun ship's whole remainder is released as the next Earth wave (ADR-0062).

## This repository's numbers (inputs for D1)

From `growth_cycle.csv` (square charge, polyethylene's carbon left atomic, bonds charged):

| closure | `k_E` | `k_S` | `v_e` Earth / Sun | mass ratio Earth / Sun | departs per kg of Earth-wave impactor | Sun-wave impactor per kg of Earth-wave impactor | growth per 3S (per all impactor) | doubling |
|---|---|---|---|---|---|---|---|---|
| 9.86 R_sun, 103 km/s stream | 24 | 30 | 5.74 / 29.9 km/s | 2.06 / 1.98 | 22.7 kg | 0.37 kg | 7.0-8.4 | 1.07-1.17 yr |
| 4 R_sun, 158 km/s stream | 54 | 30 | 5.65 / 41.5 km/s | 2.11 / 2.33 | 48.7 kg | 0.92 kg | 9.2-10.9 | 0.95-1.03 yr |

The ranges run from the Sun ship's plate kept as dead mass (0.16 of its final mass) to the plate
released as impactor.

- **Departure law** (`canted_2d.jsonl`, planar 2-D Euler with polyethylene's rebound swapped in):
  `J_n = share(k) * m w (cos theta + sqrt(1 + k - sin^2 theta))`. The share is 0.36 at k = 10,
  0.32 at k = 30 and 0.29 at k = 60. The rebound is e ≈ 0.75, not elastic. The bond toll is ~13% of
  `v_e` at k = 24. A methane-faced half-cylinder charge raises the share by ~7-10%
  (`canted_2d_layers.jsonl`); it is not yet the design charge.
- **Jupiter split** (`mission_jupiter_split.csv`): 11.70 -> 14.78 km/s excess (12.17 -> 14.16 at
  4 R_sun). The burn is 1.37 km/s at 5 R_J and 1.77 at 10 R_J, which is 0.88-1.21% of the ship
  in propellant at 330 s. Your ADR 0070 pins the perijove by the turn at 5.7-11 R_J, so the lower
  perijoves in that file are probably not reachable. This is not ADR 0068's dear division. There,
  ~47% of a wave's mass went to methalox. Here the divided stream is ~1/60 of the ship.

## What this repository's growth figure leaves out

`G` is Earth-wave impactor returned per kg of all impactor spent, on ideal closures. It does not
count:
- **The plates.** The departure's 20 m flat plate has no mass in the model, and its fate after the
  burn (dropped, carried, or released into the Earth wave) is not set. The Sun ship's paraboloid
  enters only as the 0.16 dead fraction at this repository's survival floor. The bore collar
  (~0.5-1 t) and the charge mounts are not counted.
- **The ship's dry mass**: avionics, attitude control, the disks' carrier bus and its stage.
- **The lofted mass.** The 22.7 kg departing per kg of Earth-wave impactor has to be launched and
  assembled in the 10-day orbit every cycle. That is the ledger's launch unit, and `G` does not
  price it.
- **The impactors' cost.** The disks' guidance and manufacture are taken at $1000/kg. Converting
  the Sun ship's remainder into guided sheets, and the share of sheets that miss, are not modelled.
- **The real ephemeris.** The closures are your `solve_synodic_closure(3, ...)` at 868e88a, so
  their timing, launch windows and the two streams' mirror-image timing are idealized.

## The asks

**D1: fly this architecture through the growth ledger.** Use your closures (confirming G1), the
canted-plate law above for the departure, and the Sun stage at the `v_e(k)` above (`k_S` ≈ 30).
Report the doubling time and the growth per cycle on the same footing as the child paper's
chamber and methalox rows.

**D2: charge the plates.** Size the departure's flat plate and the Sun ship's paraboloid by their
peak pressures. This repository gives ~0.1-0.3 GPa per kg of projectile on the paraboloid away from
the bore, and its canted runs record peak plate pressure per `m w^2`. Then choose each plate's
fate:
- **The departure plate:** dropped after the burn, carried to the Sun, or cut up into the Earth
  wave.
- **The Sun ship's paraboloid:** kept (G5), or released.

Carry the collar and the charge mounts as dead mass.

**D3: charge the ship's dry mass and the Jupiter split's carrier.** Add avionics and attitude
control, and the carrier bus with its dry mass and stage. Use the perijove the turn allows (ADR
0070), and check the impactors' bend at 14.78 km/s (G10).

**D4: price the launch unit.** Launch and assemble the departing mass in the 10-day orbit on the
child's price books. Also price the departure charged from the pad (G3), which your
`split_push_ledger` found costs ~6.6x.

**D5: price the impactors.** Use the child's assumptions for disk manufacture and guidance. Size
the sheets' guidance mass per sheet (G6), and charge a miss fraction.

**D6: the seed for this architecture.** The first cycle has no Earth wave. The first ships leave
on a chemical or chamber burn of 4.15 km/s from the 10-day orbit, with their disks, and the first
growth return is at 3S. Price the seed on the child's seed model (`seed_finance`, Earth loop
included). Report the seed per fleet kilogram after ten years and the break-even at L1.

**D7: the comparison the child needs.** Give one table:
- the doubling time, steady cost, seed per fleet kilogram and break-even of this cycle;
- next to the methane chamber, the hydrogen chamber and methalox, behind the spray cup;
- at both the 9.86 R_sun and the 4 R_sun closures.

Say which of D2-D5 moves the doubling time most. The child decides from that whether this
becomes its lead architecture.

## What would close it

- A doubling time with plates, dry mass, launch and impactors all charged.
- It beats the chamber cycle the child leads with now if it doubles faster than 2.0-2.8 yr, or
  2.96-4.46 yr once ADR 0068's shared split is charged.
- If it fails, say which of D2-D5 broke it, so this repository knows where to work: the canted
  plate's share, the charge, or the plate masses.
