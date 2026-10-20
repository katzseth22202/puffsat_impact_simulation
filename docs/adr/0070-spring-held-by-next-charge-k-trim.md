# ADR-0070: The next charge's k holds the lossless spring, thrusters take the net, and theta = pi is unstable

Status: proposed, 2026-10-20. Result 19 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --spring-trim-2d`, `spring_trim.py` (`make xray-collision-spring-trim`). Answers the
parent's S48 (its `docs/spring_trim_asks_for_impact_sim.md`, commit 94f8808).

## Context

ADR-0068 put the 63 t push plate on a nearly lossless nitrogen spring that loses 0.14-0.37% of
its swing a cycle. A lossless spring never forgets an impulse error. The parent estimated that a
2-10% pulse scatter walks the plate off its 8 m stroke in 25-625 pulses. It ruled out burst and
coast for the plate (0.85 TJ a push), and asked which per-pulse trim holds the spring. Seth named
two levers (2026-10-20): conventional thrusters for the gap between the push expected and the
push achieved (a small budget, not yet charged), and modest changes to the next charge's k.

## Decision

1. **Trim k on every pulse, not only when needed.** Untrimmed, the flown spring (2 pi/3) leaves
   its stroke by 4 m in a median 85 pulses at 1% scatter, and 40 at 2%. Every seed is lost within
   the push at every loss Result 17 found. That is 5-30x faster than the parent's estimate,
   because an error rings at the free frequency (`dv / omega` = 0.24 s x dv) while the stroke is
   only a 2 pi/3 arc.
2. **Thrusters take only the summed error.** A discrete LQR on the one-period map sets the next
   charge's k from the swing measured over the 0.5 s interval. The summed impulse error is one of
   its states, so the k trim also holds the craft's delta-v, and the thrusters clear what is left
   at the push's end: 0.16-0.45 m/s p95, or 61-173 kg at 300 s. Thrusters cannot hold the spring
   themselves: cancelling each pulse's error at 1% would take 117 MN s and 40 t.
3. **The trim is modest.** At the budgeted 0.85% scatter the next charge's k varies by 0.24 rms
   (6.5-7.7 at 1-99%), or 6.3 kg of a 185 kg charge. The 2-D slope is 3.3% of J per unit of k at
   fixed sheet mass. The plate then strays ~1 m (p95) past each end of its steady stroke, so build
   the bags for ~10 m, not 8. The trim holds to 2% scatter (2.4 m, k 5.8-8.5). At 5% the
   ladder's k = 5-9 saturates and the spring is lost.
4. **Do not fly theta = pi. Prefer 5 pi/6.** On the nonlinear column, theta = pi is unstable (a
   real eigenvalue of -2.006: an error doubles every pulse). The k trim holds it only at <= 1%
   scatter and a lag of <= 1 pulse. The cycle stays stable up to theta = 0.850 pi (an 11.5 m
   column). At 5 pi/6 the column is 11.9 m, not 18.6, and the trim holds it better than the flown
   spring (0.77 m p95, 130 kg of thruster propellant, and still ~1 m at a lag of 2 pulses).
5. **Standoff, timing and an orifice are not needed.** The standoff is a lever of about the same
   size (3.6% per metre), but it means placing the charge against a plate that swings +-4 m. At
   2 pi/3 the next charge reaches the swing's phase as well as its amplitude, so timing
   adds nothing: a 0.1 ms arrival error leaves the result unchanged. Result 17's loss barely
   matters, so a switched orifice would have to be far lossier than the bags to help, and it puts
   a valve in a sealless bag.
6. **The aim, not the heave, sets the tolerance.** Rocking shares the heave's map, and k cannot
   reach it. A centroid 10 cm off the axis kicks the floor's edge at 4.2% of the heave kick, five
   times the heave budget. With the charge's areal density tilted (a wedge of up to 30%), the
   rocking holds at an aim of ~3 cm per axis or better (edge 1.4-2.0 m p95). The alternative is
   the parent's "alignment guides": guides that hold the tilt hand the torque to the craft, which
   has no spring. Thrusters then take only its sum (6.9 MN s, 2.4 t at 10 cm aim). The guides are
   named as an option and deliberately not sized (Seth, 2026-10-20).

## Consequences

- The paper can say the spring is held by trimming each charge's k from the swing measured
  between pulses, with thrusters clearing the residual delta-v. The stroke grows ~25%, and theta
  should be kept below ~0.85 pi.
- The heave budget is dominated by how well each sheet's mass is known on arrival (1% assumed
  gives 0.77% of J). Weighing and tracking the sheets is the cheap way to shrink it.
- Open: the taper of the pulse as propellant burns (a feed-forward; not modelled). k outside
  5-9. The loss as thermal hysteresis rather than an equivalent drag. Rocking on the 4-cell bag layout, whose rocking stiffness differs from its
  heave stiffness.
