# The thin 1 kg disk on a wider charge, and 4 R_sun: ask S62 (for `puffsat_impact_simulation`)

Raised 2026-10-27 against sim `3581fbd` (Result 23, ADR-0074) and aim `d37a7ec` (ADR 0086).
**Written to be copied verbatim into `katzseth22202/puffsat_impact_simulation`.** Register entry:
S62 in `docs/deferred_to_companion_repos.md`.

Same ground rules as your own handbacks. If a premise here is wrong, say so rather than working
around it.

## Why

aim's ADR 0086 turned Result 23 into an aim tolerance and re-priced the Sun wave. Its cheapest
train is the headline's thin disk (1 kg, 40 cm, 8 kg/m^2, 2 Hz) on a charge 1.3x its radius,
which breaks even at $129 / $155 against the $131 / $158 the paper prints. That train flown on its
matched charge at the real aim costs $142 / $181, almost all of it the unstruck lens.

The parent has **not** adopted $155. It rests on two inputs Result 23 did not measure for this
disk:

- the thin disk on x1.3 ran only on the 0.025 m grid (1.08x; the grids agree to 0.4%);
- its bore clearance and side impulse were not traced. aim took a 5 cm clearance (0.25 m bore,
  0.20 m disk radius) and the thick disk's x1.3 side-impulse slope, 13.5% per metre of offset.

The paper prints $142 / $181 as the honest figure for today's design and $129 / $155 as pending
this ask.

## Ask

1. **Thin disk, x1.3 charge, finest grid.** Run the 1 kg, 8 kg/m^2 disk on the x1.3 charge at
   0.0125 m. Report the coupling over the thin disk on its matched charge.
2. **Its bore clearance and side impulse.** Trace the offset at which the rim clips the bowl's
   bore, and the side impulse against offset (at least two offsets, so linearity can be checked).
   Same table shape as `thick_disk_offset.csv`.
3. **If x1.3 is not the thin disk's best,** say which width is, with its clearance and side
   impulse. A x1.5 charge through a wider bore is welcome as a second point.
4. **4 R_sun.** Result 23 ran 392 km/s only. Repeat items 1-2 at the 4 R_sun closing speed so the
   $85 / $89 headline there can be restated under the same aim.

## Outcome

If items 1-2 land near aim's assumptions, the parent adopts the 1 kg disk on the wider charge as
its headline everywhere (8 places at Parker's depth plus the sensitivities that hang off $158),
and aim reruns `make disk-train-aim` with the measured inputs.
