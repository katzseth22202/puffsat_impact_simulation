# Plate shape: flat floor + shallow concave; the deep dish is foreclosed at high-v, not a contingency

The plate is one built object (§5); its shape is a 2D-sweep variable (flat / 0.10 / 0.15
depth-to-diameter), not a free knob. This ADR pins which shapes are in play and — the non-obvious
part — *why the deep dish is foreclosed precisely where its rebound re-collimation would help most*,
so a reviewer's "why not just deepen the dish to recover the sideways loss Orion suffered?" is
answered up front. Scope is mechanism #2 only: gas that lands, stagnates, and squirts radially
*across* the plate (the rebound-axiality loss in `eta_capture`, channel 4 of ADR-0016). Interception
(gas missing the plate) is the cloud schedule's job, not the plate's.

**Rebound axiality is primarily a cloud-shape race; curvature is the recovery lever.** For a flat
plate, how axially the stagnated gas rebounds is a race between **radial relief** (`~r_foot/c_s`, the
squirt out to the free edge) and **axial drain** (`~L/c_s`, re-expansion back along the incoming
axis). A wide, short cloud (disk) is closest to the 1D plane-wave ideal — slow radial relief, fast
axial drain → axial rebound → high `eta_capture` — but it stagnates all at once → high peak pressure
→ fails the survivability frontier. A thin, long cloud (cylinder) is pressure-friendly but splats
radially → low `eta_capture`. So `L/D` trades `eta_capture` directly against survivability, and the
pressure limit — worst at **16 km/s** (~393 MPa disk-like) — *forces* the cloud toward the splat-prone
elongated shape. Curvature catches that forced radial outflow and bends it back axial, so its
`eta_capture` value **peaks at the 16 km/s anchor**, the same place survivability forces the stretch.

**The deep dish is unusable exactly where it would help most.** Maximum re-collimation would come from
a deep dish, but its focal hot spot concentrates the rebound in the **strongest-radiating** (`flux ~
v⁸`), **`τ≫1`** gas and couples radiation into the otherwise radiation-free 2D geometry track
(ADR-0003/0008), breaking the factorization. Both penalties are worst at 16 km/s — so the
maximum-re-collimation geometry is foreclosed at the one anchor that most wants it. It stays an
**unrun upper-bound note, not a high-v contingency.**

## Amendment (2026-06): the shallow-concave `eta_capture` gain is measured, and it rises with depth

The flat / 0.10 / 0.15 depth-to-diameter shapes are now run through the 2D kernel's ghost-cell IBM
(ADR-0023 amended), confirming curvature *is* the rebound-axiality recovery lever this ADR posits:

- **Shallow concave lifts `eta_capture` above the flat floor, monotone in depth.** At a short
  disk-like cloud (`L/D = 0.3`, `r_foot/R = 0.3`, M ≈ 5) `eta_capture` rises **0.88 → 1.02 → 1.03**
  over `d/D = 0 / 0.10 / 0.15`; the gain holds across the swept `L/D` and footprint grid. Concave can
  **over-collimate past 1** — it bends the rebound *more* axial than a flat plane wave (the 1D limit),
  so `eta_capture > 1` is physical here, not a numerical artifact (and not a contradiction of ADR-0003's
  "1D ⇒ `eta_capture = 1`", which fixes the *flat* plane-wave denominator, not a ceiling on shaped plates).
- **It clears the useful gate where the cloud need not be stretched.** Combined with the 1D `e_eff`
  (ADR-0003 `f`-reconciliation): at the transitional dip `e_eff = 0.57` the flat floor is `f ≈ 0.696`
  and shallow concave reaches **`f ≈ 0.83` (> the 0.8 useful gate, ADR-0009)**; at 16 km/s
  (`e_eff = 0.63`) concave reaches `f ≈ 0.86`. The recovery lever is real at the low/mid-v anchors
  where there is no focal penalty — exactly the regime this ADR reserved it for.

This is a shallow-plate result at two Mach anchors; the deep dish stays foreclosed (above), and the
focal-radiation penalty that forecloses it at 16 km/s is untouched by these radiation-free runs.

*Corrected magnitudes (2026-07):* the numbers above were from the unconverged 56×40 grid at M ≈ 5.
On the converged 112×80 grid at the physical Mach anchors (10/20), the same short-disk case gives
`eta_capture` **0.915 → 0.993 → 1.009** over `d/D = 0 / 0.10 / 0.15` at M = 10 (0.88 → 0.92 → 0.94
at M = 20): the depth-monotone concave lift and the marginal over-collimation past 1 both survive,
with smaller magnitude. The corresponding concave `f` maxima are ≈ 0.792 (dip) / 0.822 (16 km/s);
the decision (shallow concave is the lever, deep dish foreclosed) is unchanged.

*Second correction (2026-07-10, ADR-0023 kernel fix — the dish rim's side face was missing from the
immersed boundary, inflating the M = 10 concave etas and depressing M = 20):* the same case now gives
**0.915 → 0.977 → 0.994** at M = 10 and 0.92 → 0.97 → 0.99 at M = 20 — the two Mach anchors agree to
< 1 % (the earlier M = 10/M = 20 gap was the bug), the depth-monotone concave lift stands, but the
**over-collimation past 1 is gone** (sweep-wide max `eta_capture` = 0.994): the `eta > 1` readings
were the rim-corner artifact, not physics. Concave `f` maxima become ≈ 0.780 (dip) / 0.810 (16 km/s);
the decision is again unchanged.

## Amendment (2026-06): Orion's geometric loss is *interception* (collimation `C₀ ≈ 0.5`), not sideways rebound — a conditional PuffSat advantage to claim later

Sanity-checking the §3.2 Orion comparison turned up a **distinct, larger geometric difference than the
rebound-axiality lever this ADR scopes** — but it lives in mechanism #1 (interception), the cloud
schedule's job, so it is recorded here as a *deferred finding, not yet a paper claim*.

- **What Orion concluded (verified against the Orion literature).** Orion's plasma **reflects and
  transfers momentum** to the pusher plate — it does not "escape sideways"; reflection *is* the
  propulsion principle (~6000 s Iₛₚ). Its dominant *geometric* loss is the **collimation factor
  `C₀ ≈ 0.5`** — only about half the bomb's debris hits the plate, the rest overspilling the edge from
  a quasi-point-source expansion — improvable by shaping the propellant (a disk reaction mass → a
  focused "cigar" jet) and matching plate diameter to fireball diameter. This is **interception**
  (debris *missing* the plate), not rebound axiality (`eta_capture`), and not the radiative §3.2 thesis.
- **The PuffSat difference (conditional).** A *controlled, shaped gas cloud sized inside the plate
  footprint* (`r_foot/R < 1`, "room before the edge") can drive PuffSat's collimation factor toward
  **~1 vs Orion's ~0.5** — a near-2× geometric advantage, and a *separate multiplicative factor* from
  the `eta_capture` rebound term (the 2D sweep already assumes the gas lands, i.e. implicitly credits
  collimation ≈ 1). **Gate:** this is only a paper claim *if the deferred cloud-schedule / delivery
  study confirms PuffSat can actually deliver the cloud onto the footprint*; until then it stays a note.
- **Citation caveat.** `C₀ ≈ 0.5` is from secondary Orion references (standard Orion literature, ultimately
  Dyson / General Atomics). The paper's "Balcomb 1970" = **LA-4541-MS, "Nuclear Pulsed Space Propulsion
  Systems," J.D. Balcomb et al., LASL, Oct 1970** — the *laser-driven* micro-pulse variant using the
  Orion pusher plate, adjacent to bomb-Orion. The primary PDF was not read (OSTI firewalled), so confirm
  the exact `C₀` attribution against the original Orion sources before quoting it.

## Amendment (2026-06, Rung S): the `L/D ↔ survivability` trade is now measured

This ADR asserts that `L/D` "trades `eta_capture` directly against survivability" and that the pressure
limit "*forces* the cloud toward the splat-prone elongated shape." The survivability frontier (design
§7, ADR-0010/0011) now measures it. Via the `Σ = m/(π r_foot²) = ρL` contract, a short disk (small
`L/D`) or tight footprint (small `r_foot/R`) packs the mass into a dense column → high peak pressure;
so the **`f`-maximizing corner (short disk + tight footprint) is exactly the densest, least survivable
case** — it peaks at ~2.3 GPa at 16 km/s (with the concave local-peak focusing folded in, below) and is
foreclosed. The survivable optimum is the intermediate, elongated, wider-footprint shape this ADR
predicted: best survivable `f ≈ 0.80` at the dip, ≈ 0.78 / 0.84 at 16 km/s (baseline / relaxed).

**The concave focusing penalty is also now quantified** (Rung S `euler2d` `max_plate_pressure`): a
shallow dish focuses the rebound to a peak *local* facesheet pressure **1.0–2.4×** the flat plane-wave
value (the survivability cost of the re-collimation gain this ADR credits concavity with — and the
shallow-plate cousin of the focal hot spot that forecloses the deep dish). It is folded into the
frontier, so the survivable-`f` numbers above already carry it (e.g. the best dip survivor peaks at
370 MPa *after* a 1.78× focusing penalty). Curvature thus cuts both ways — it lifts `eta_capture` but
raises the local peak — and the shallow `d/D ≤ 0.15` band stays net-favorable; a deeper dish would not.

## Consequence

`eta_capture(v)` is a measured sweep output, **allowed to come in below 1 at high-v** — it is never
propped up by assuming shallow concavity closes the gap. If shallow concave + the best cloud miss the
target at 16 km/s, the recovery order is: (1) push `r_foot/R → 1` and trim divergence (free, no focal
penalty); then (2) **accept and report** the resulting `f(16)` via the dual curve (ADR-0013) — a
sub-0.8 anchor is a reported result, not a failure. The deep dish is **never** reopened, because it is
worst exactly where the shortfall would occur.

## Considered Options

- **Deep dish** (highest `f` ceiling, strongest rebound re-collimation). Rejected: focal hot spot
  couples radiation into the radiation-free 2D track and lands in the strongest-radiating gas — both
  worst at the high-v anchor where its recovery would be most wanted.
- **Flat only** (simplest build, conservative floor). Rejected as the *sole* shape: it leaves the
  cheap rebound-fanning recovery on the table at the low/mid-v anchors where there is no focal penalty
  and the cloud need not be stretched as hard.
