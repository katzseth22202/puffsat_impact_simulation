# The standoff wall is the liner, not the bag bore, and the bag volume stays 660 m³

Status: accepted (2026-09-05, working replies R14 and R15). Unblocks the paper's ADR-0013, which
is proposed pending this. Answers P19 and P20 of `docs/nozzle_replies_answered.md`.

Two decisions, recorded together because they are the same measurement asked twice — *which
surface, and how big is it*.

## 1. The standoff requirement is written against the liner

The paper's ADR-0012 says the field's job at the chamber is "wall standoff — the field is there to
keep plasma off the liner", and then measures its cap table to the 3.02 m **bag bore**. Its
ADR-0013 proposes measuring to the 3.50 m flared liner instead and asks this repository which
surface the requirement means.

**The liner.** The mist is inside the 3.0 m bag and the clearance gap between bag and liner is
**vacuum**. A front that expands past the bag bore therefore sweeps no *less* mass — it already
spans the bag's entire cross-section, which is all there is to sweep — and no *more*, because
there is nothing out there. Coupling is indifferent to the excursion. The only thing the field
must prevent is plasma reaching the graphite.

`front.integrate` caps swept area at the bag radius rather than at the wall for exactly this
reason, and confirms the swept mass is identical to nine digits whether the wall is placed at
3.02 m or 3.50 m.

**Two independent effects make the liner reading the softer one**, not one: the contact station
moves out, *and* expanding from 3.00 m to 3.50 m dilutes the front's pressure by
`(3.5/3.0)^{-2 gamma}` = 0.60 at `gamma` = 5/3.

*Caveat, and it is ours.* The lateral excursion into the gap is momentum going sideways that a
hard wall would have reflected. The gap is 0.5 m against a 3 m radius, so the volume involved is
1.36x and the effect is second-order. This is an argument, not a resolved run.

**Consequence for P9.** The graded profile 20/12/9/5 T is derived by setting the snowplow's
pressure against `eq:bore_from_length`'s bore area *at every station*, so **every station of it,
the 5 T exit included, is a standoff number evaluated at 3.02 m**. The paper asked whether that is
so; it is, which reopens the `1/r`-graded profile it had set aside.

## 2. The bag volume stays 660 m³; the bore follows from it

The paper's R14 adopts this repository's 23.8 m and infers 672.9 m³ from it, on the strength of
three of its own round figures matching. **Do not adopt it.** Our two numbers came from different
places and were never asserted as a pair:

- `expansion.FIELD_LENGTH = 23.8` came from the paper's own 660 m³ at its own aspect ratio 4:
  `8 pi r^3 = 660` gives `r` = 2.9722 m and `l = 8r` = 23.777 m.
- `expansion.CHAMBER_RADIUS = 3.0` came separately from the paper's quoted "3 m bore".
- `BAG_RHO = 0.323` = 213/660, i.e. on the 660 basis.

So our own radius and density sit on bases 2% apart in area — a defect on this side, recorded
here. And 672.9 m³ contradicts the 660 m³ that `PV = nR_gT` fixes (`field.V_STANDOFF_M3`).

**The volume is the physics; the bore and the aspect ratio are consequences.** The four-way
consistent solution is 660 m³, `l` = 23.78 m, `r` = **2.972 m**, cross-section 27.75 m², aspect
4.000 — which still rounds to every figure the paper quotes ("3.0 m", "28 m²", "4") *and* keeps
the constraint.

**Not yet applied here.** Moving `CHAMBER_RADIUS` to 2.972 m is a 1.9% change in area that would
invalidate every artifact `docs/nozzle_replies_answered.md` cites, for an effect below the width
of every bracket in it. It is recorded as owed and should be made in one pass once the paper
settles which geometry it flies.

## Provenance

`make analysis-nozzle-front` for the contact stations and the swept-mass invariance;
`field.V_STANDOFF_M3` and `expansion.py`'s constants for the geometry provenance.
