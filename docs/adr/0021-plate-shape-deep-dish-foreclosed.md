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
