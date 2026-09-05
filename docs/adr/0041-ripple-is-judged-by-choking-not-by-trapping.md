# A field ripple is judged by choking, not by trapping, and P8's coil-count requirement is withdrawn

Status: accepted (2026-09-05, working reply R10). Withdraws P8's ">= 36 coils" from
`docs/nozzle_asks_answered.md`. Follows from ADR-0040.

## What was wrong, on both sides

`field.ripple_sweep` reported whether an off-axis `|B|` minimum **exists**, and P8 concluded that
the paper's residence claim needs at least 36 coils. The paper's reply R10 objects, correctly,
that existence is the wrong test: a magnetic trap only holds what its mirror ratio can hold, so
the criterion should be `R > 1/sin^2(theta)`.

**R10 is right about the defect and wrong about the fix**, and the correct fix retires the
criterion rather than sharpening it.

Mirror trapping *is* `mu` conservation seen from the other side: a particle turns around because
its own `mu` is fixed, and the loss cone is the set of pitch angles for which that never happens.
ADR-0040 shows there is no `mu`. The fluid statement is cleaner and needs no kinetic theory at
all:

    (J x B) . B = 0   identically

so **the Lorentz force has no component along a field line**. In a collisional plasma the pressure
is a scalar, the parallel momentum equation is `rho Du/Dt = -dp/ds` with no magnetic term, and a
`|B|` variation exerts no force along the flow whatsoever. The mirror force that appears in
kinetic theory, `-(p_perp - p_par) d ln B/ds`, is proportional to an anisotropy that collisions
erase on the mean-free-path timescale — about 1e-10 s here, against a 2 ms transit.

R10's threshold also takes `alpha` = 0.088 from P1, which is the anisotropy of a **free expansion
with no nozzle**, not of plasma in the bore.

## The decision

**The criterion is choking.** A flux tube's cross-section is `A ~ 1/|B|`, so a local `|B|`
*maximum* is a local area *minimum* — a throat. Supersonic flow through a contraction decelerates,
and a deep enough contraction drives it to `M` = 1 and stands a shock there. Isentropic flow at
Mach `M` survives a contraction up to `A/A*(M)`, and a ripple's contraction ratio is exactly the
mirror ratio `field.scan_axial` already computes. So:

    a ripple is admissible iff   R < A/A*(M_local)

**Each minimum is judged at its own station**, with `M(z)` read from the solved cooling history
(`residence.mach_profile`), and **only minima in supersonic flow are tested** — upstream of the
sonic point a contraction accelerates the flow and the worst it can do is move the sonic point.

## Where it binds, and it is not where P8 looked

`A/A*` goes to 1 at the sonic point, so the margin *vanishes* at the chamber and grows fast
downstream: 1.03 at `M` = 1.2, 1.83 at `M` = 2.0, 9.67 at `M` = 3.4. **Every binding minimum in
the sweep lands at `z` = 0.6–2.8 m**, on ADR-0012's flat 12 T shelf, where there is no background
gradient to swamp ripple and the flow has almost no Mach margin. That is exactly the interaction
the paper's R10 predicted from the other direction.

## What changes numerically

| criterion | straight winding, flown | flared winding, capped |
| --- | ---: | ---: |
| **choking (adopted)** | **>= 12 coils** | **>= 18 coils** |
| R10's loss cone at `R` = 1.096 | >= 24 | >= 24 |
| P8's existence test (withdrawn) | >= 36 | >= 72 |

P8's requirement falls by a factor of two to four. R10's own proposal is much closer than P8's but
still about twice as strict as the physics demands.

## Consequences

- `sec:jet_efficiency`'s residence claim survives at a far lower coil count, and the tension the
  paper notes between "continuous winding" and "graded profile" is much weaker than P8 made it.
- The flared winding needs **more** coils than the straight one, because its binding minima sit on
  the flat shelf. That is a real cost of ADR-0012's cap and the first one found.
- `field.RippleResult.trapped_fraction` is now a historical quantity. It is kept because it is
  what P8 reported and the answer document has to be able to say what changed.

## Provenance

`make analysis-nozzle-residence`. The area-Mach relation is pinned against the standard
`gamma` = 1.4 table (`A/A*` = 1.6875 at `M` = 2, 4.2346 at `M` = 3) so the criterion rests on an
external anchor rather than on this repository's own arithmetic.
