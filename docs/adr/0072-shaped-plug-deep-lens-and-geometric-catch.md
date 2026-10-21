# ADR-0072: The ice plug is a lens twice the bag's depth, and the gun fires past its geometric catch

Status: proposed, 2026-10-21. Result 21 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --shaped-plug [overlap [fine]|backup [far]] [coarse]`, the dished and open-core
plug of `euler2d::bounce::BumperCharge` (`charge_dish`, `charge_hole`), and `shaped_plug.py`
(`make xray-collision-shaped-plug`). Answers the parent's S52 (its
`docs/shaped_plug_asks_for_impact_sim.md`, commit ec5995d); answers doc
[`shaped_plug_answers_for_parent.md`](../shaped_plug_answers_for_parent.md). Builds on ADR-0071.

## Context

Seth (2026-10-21) took Result 20's ice charge and graded paraboloid, and turned the gun's ice
bullet from the spreader into a backup that fires only on a predicted miss. The plug the bag
strikes should itself spread the bag into a near-point burst. The parent sized that plug as a
squat lens as deep as the bag (Birkhoff's penetration depth, ~0.38 m) and asked for its leverage
at k = 1-15, the thinnest plug that stops the bag, the share of the bag that may miss the plug
before the face passes 2.3 GPa (which sets the gun's trigger and replaces the ledger's 5%
placeholder backup rate), the backup pulse, the gun's firing window and the plate.

## Decision

1. **The plug is an ice cylinder twice the bag's depth (~0.73 m) at 2 F.** The parent's 1x lens
   punches through at every k (2.8-2.9 GPa on the axis, from the flattened bag's jet). At 2x,
   k = 4-15, the lens gives the point fireball's leverage to within 3% (8.18x, 8.81x, 8.39x;
   frozen 7.32x, 7.65x, 7.37x) at a face peak of 0.69-0.90 GPa. The 1.5x lens is the thinnest
   that stops the bag (1.1-2.0 GPa) but leaves too little margin for the grid's +23% on
   refinement.
2. **Keep 2 F; no dish.** 1.5 F gains 1.5% at k = 8.5 for twice the peak; 2.5 F halves the peak
   for 1-3%. A dished face cups the jet onto the axis (5.9-7.5 GPa at 1x) and gains only at k = 4.
3. **The lens pays only above k ~ 4.** Below, the compact column is already short and the three
   shapes agree within 2%. At k = 1-2 no depth reliably holds the plate.
4. **The trigger is the geometric catch.** The gun fires when the predicted offset of the bag's
   centre from the plug's exceeds the plug's radius less the bag's: 0.04 m at k = 4, 0.14 m at
   k = 8.5, 0.24 m at k = 15. The 2-D overlap brackets (an open core, a narrow plug) do not
   converge (the safe share falls 56% -> 19% -> 8% as the cell halves twice), since axisymmetry
   focuses any passing ring onto the axis. They are no basis for credit, so none is taken until a
   3-D run.
5. **The backup rate is the Rayleigh tail past that catch.** At k = 8.5 the parent's 5% needs a
   track error at the plug of 5.6 cm per axis (1%: 4.5 cm). At k = 4, 1.5 cm. Result 20's
   0.1-0.15 m was the gun's aim, not a known track.
6. **A faster gun.** The bullet's meeting distance is set by the cloud's spread: 30.6 m for
   100 g with no plug in the cloud, 45.9 m with the k = 8.5 lens in it (decision 7). At 230 m/s the
   gun then fires 6.15 km or 9.43 km out, before the plug's last 100 ms. At >= 500 m/s it fires
   after (2.2 km at 1 km/s for the lens case), so the decision sees the steered plug. A 50 g
   bullet does not help.
7. **The backup pulse spends the plug, and at k = 8.5 the bullet meets the bag ~46 m out.** With
   the 2x lens in the cloud's middle, the pulse gives 2.4-2.6 MN s (4.16x at k = 8.5). It holds
   1.80 GPa at 46 m, against 2.27 GPa (no margin) at 31 m. The plug is kicked ~98 m/s off
   station. At k = 4 no meeting holds (2.41-3.04 GPa): the narrow plug turns the spread cloud
   into an axial jet. The gun fires 9.4 km out at 230 m/s, or 2.2 km at 1 km/s.
8. **The Earth bag plate is 27-39 t graded** (57-70 mm of steel at the centre, 0.3-4.7% ringing),
   and its graphite liner takes 0.4-18 kg a pulse (1-35 t a push), least at high k.

## Consequences

- The parent's ledger carries the lens's leverage as its simulated base (the point fireball's
  former "optimistic ceiling"). The backup rate is a function of k and of the track error at the
  plug, not a constant. Lower k now costs catch as well as leverage.
- The catch radius is a lower bound. A 3-D off-centre strike is the one run that could widen it.
- `BumperCharge` gains `charge_dish` and `charge_hole` (both 0 for every earlier run). Its vacuum
  stand-in's pressure is capped so its heat stays under 0.8x the bag's. That binds only on the
  46 m far-bullet domain.
