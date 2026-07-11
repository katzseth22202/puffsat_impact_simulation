# Frozen-recombination bracket: equilibrium `e_eff` stays the headline; the sudden-freeze bound is recorded, not adopted

The 2026-07 physics audit flagged one approximation stacked *optimistically* in an otherwise
conservative-floor study: the equilibrium water EOS **returns the dissociation/ionization energy
during the rebound**. On compression the gas banks a large chemical store (at the ~11 km/s dip
the EOS sink is exactly why `e_eff` bottoms out); on re-expansion the equilibrium assumption
hands that store back as pressure work. If recombination instead *freezes* on the fast rebound
(the classic nozzle-flow effect), the chemical energy leaves as inert enthalpy in the escaping
gas and `e_eff` falls **below** the quoted floor. Crucially, a FLASH cross-check run with the
same equilibrium EOS **shares this assumption and cannot catch it** — so it had to be bounded
in-house. This ADR records the bounding study (the `frozen` pipeline: `make sweep-frozen-probe
→ tables-frozen → sweep-frozen → analysis-frozen`) and the decision on how to quote it.

## Method: bracket the freeze timing both ways

Finite-rate chemistry is out of scope, but the two extremes of freeze *timing* bound it:

- **Freeze after the plate (pessimistic): sudden-freeze splice at global turnaround.** Each
  transitional-grid case runs on the equilibrium table to the global momentum zero (where the
  chemical store is maximal), records the mass-weighted `(rho*, T*)`, and swaps to a
  **frozen-composition table** built at that state: equilibrium fractions locked → ideal
  mixture with constant mean molecular weight plus a constant, *inert* chemical-energy offset.
  The splice is temperature-continuous and zero-point-consistent (it reproduces the equilibrium
  `p, e` exactly at the freeze reference state); out-of-table cells carry their energy
  deficit/excess across the zero-point shift. The bookkeeping jump is a logged diagnostic:
  ≤ 7.3 % of incident KE everywhere, ≤ 1.7 % at the dip, and inert (chemical offset only —
  the thermal state driving the dynamics is unchanged).
- **Freeze before the plate (optimistic): frozen throughout as pure H₂O.** The whole bounce
  runs on a chemistry-free water-vapor table — no dissociation/ionization sink at all.
- **Equilibrium** (the study curve) is re-run alongside as the reference.

## Result: the bracket is wide, and ordered correctly everywhere

ρ-mean `e_eff(v)`, EOS-only transitional grid (`data/results/frontier_frozen.csv`, figure
`data/results/frozen_e_eff_v.png`); the invariant *frozen-throughout ≥ equilibrium ≥
sudden-freeze* holds at every `v`:

| v (km/s) | frozen throughout | equilibrium | sudden freeze |
|---|---|---|---|
| 5  | 0.708 | 0.697 | 0.605 |
| 8  | 0.674 | 0.596 | 0.467 |
| **11 (dip)** | **0.661** | **0.570** | **0.398** |
| 16 | 0.653 | 0.640 | 0.509 |

Translated through `f = eta_capture·(1+e_eff)/2` at the survivable geometries:

- **Dip:** `Δe_eff = 0.172` → `Δf ≈ 0.085`; best survivable `f` 0.777 → **≈ 0.69** under
  sudden freeze (frozen-throughout: ≈ 0.82).
- **16 km/s:** `Δe_eff = 0.131` → `Δf ≈ 0.064`; `f` 0.805 → **≈ 0.74** (frozen-throughout:
  ≈ 0.81).

*(2026-07-10, ADR-0023 kernel correction: the equilibrium anchors become 0.768 dip / 0.798 at
16 km/s, so the bracket translates to ≈ 0.68 / ≈ 0.74 under sudden freeze and ≈ 0.81 at the
frozen-throughout end — the conclusion below is unchanged.)*

So the pessimistic end of the bracket **crosses the 0.8 gate** — this is ~3× the ±0.03
numerics band, the largest single physics uncertainty surfaced by the study. It cannot be
waved into the noise; it needs a physical argument, which exists:

**Sudden freeze at turnaround is a maximal bound, not an expectation.** It freezes at the
instant of *maximum* chemical storage and forbids recombination for the entire rebound. At the
probed turnaround states (dip: `rho* ≈ 2–8 kg/m³`, `T* ≈ 7 kK`, i.e. `n ~ 10²⁰–10²¹ cm⁻³`),
three-body recombination runs at `k·n² ~ 10⁷–10⁹ s⁻¹` per atom (`k ~ 10⁻³² cm⁶/s`) — a ~1–100 ns
timescale against the ~µs rebound. The gas therefore tracks equilibrium through the dense phase
where essentially all of the momentum is exchanged; composition freezes only late in the
expansion, at ~100× lower density, when the push is already delivered. The physical curve sits
near the equilibrium end of the bracket.

## Decision

1. **The equilibrium curve remains the headline** (`f ≈ 0.77–0.82` across the envelope, as
   corrected 2026-07-10 — ADR-0023), on the kinetics argument above.
2. **The bracket is quoted as a named caveat wherever `f` is quoted** (CONCLUSION.md, white
   paper): "under a maximal sudden-freeze-at-turnaround bound, best survivable `f` drops to
   ≈ 0.68 (dip) / 0.74 (16 km/s); recombination kinetics at the probed densities argue the
   physical answer is near the equilibrium end."
3. **The FLASH cross-check does not discharge this caveat** if run with an equilibrium EOS —
   it is flagged so the cross-check is not over-credited. The refinement that would close it is
   a **finite-rate (or partial-equilibrium) chemistry rung**: e.g. freeze each case at the
   density where `τ_rec` first exceeds the local expansion time, instead of at turnaround.
   Non-blocking: it interpolates a bracket whose physical end is already argued.

## Amendment (2026-07-11): the 69 km/s Jupiter scenario has its own, wider bracket — and the kinetics defence is weaker

The Jupiter-retrograde 69 km/s special scenario (ADR-0007 Jupiter table, `sweep-jupiter`) needed
its own freeze bracket: the turnaround is a **multi-charge oxygen plasma at ~140–170 kK**, which
the transitional grid (T ≤ 60 kK, single O⁺ stage) cannot even represent. A parallel pipeline was
built on the extended grid — `sweep-frozen-probe-jupiter → tables-frozen-jupiter →
sweep-frozen-jupiter → analysis-frozen-jupiter` — with `FrozenComposition` carrying the **full
8-stage O Saha ladder** (each stage freezing its cumulative ionization energy) so the splice stays
exact at the hot reference state. Anchored at the realistic `L = 12 m` cloud over the survivable
`JUP_RHO` grid (`data/results/frontier_frozen_jupiter.csv`):

| ρ_impact (kg/m³) | frozen throughout (pure H₂O) | equilibrium (EOS-only) | sudden freeze |
|---|---|---|---|
| 0.01 | 0.647 | 0.695 | 0.502 |
| 0.07 | 0.647 | 0.711 | 0.533 |
| 0.16 | 0.647 | 0.716 | 0.548 |

Two things differ from the baseline anchors:

1. **The ordering inverts: pure-H₂O is *not* the upper bound here.** At this ionization the
   chemistry-free and equilibrium tables are genuinely different EOS (not a sink on/off toggle),
   so the transitional invariant *frozen-throughout ≥ equilibrium* fails — equilibrium sits above
   both freeze timings. The load-bearing bracket is therefore **[sudden-freeze, equilibrium]**.
   Translated through `f = eta·(1+e_eff)/2` onto the coupled headline at the `R = 22 m` concave
   design point (`ρ ≈ 0.04`, `Δe_eff = 0.182 → Δf ≈ 0.089`): equilibrium headline **`f = 0.781`**
   → **sudden freeze `f = 0.692`** (pure-H₂O bound `f = 0.752`). Splice energy-jump diagnostic
   ≤ 3.5e-3 of incident KE. So the 69 km/s bracket on `f` is **[0.69, 0.78]**.
2. **The kinetics argument that keeps equilibrium the headline is *weaker* at 69 km/s.** The
   defence above rests on the turnaround being dense (`n ~ 10²⁰–10²¹ cm⁻³`) so three-body
   recombination outruns the rebound. The 69 km/s turnaround is **dilute** — `ρ* ≈ 0.09–1.2 kg/m³`
   (`n ~ 10¹⁹ cm⁻³`, ~10–100× below the dip) — **and** multi-charge, so there is an 8-stage
   ionization ladder to unwind against a fast, low-density expansion. Equilibrium is still the
   better central estimate (the dense stagnation phase, where the momentum is exchanged, does
   track it), but the sudden-freeze end is **less firmly foreclosed** here than for the baseline
   anchors. The `[0.69, 0.78]` bracket is correspondingly more load-bearing and must accompany any
   quoted 69 km/s `f`. The named refinement (finite-rate / freeze-out-density chemistry) matters
   *more* for this scenario than for the envelope.

## Considered Options

- **Adopt the sudden-freeze curve as the new conservative floor.** Rejected: it is not a floor
  estimate but a physically-foreclosed extreme (it needs recombination ~10²–10³× slower than
  the three-body rate at the probed densities); adopting it would misstate the study's central
  value while the honest instrument — the quoted bracket + kinetics argument — already exists.
- **Ignore the frozen bound (equilibrium is obviously fine).** Rejected: the bracket is 3× the
  numerics band and the pessimistic end crosses the gate; unquantified, it was the audit's top
  open finding, and no external cross-check with an equilibrium EOS would ever surface it.
- **Full finite-rate chemistry in the kernel.** Rejected for now: real kinetics networks are a
  scope step-change for a bounding question; the two-sided bracket plus the timescale argument
  answers "is `f` drastically influenced?" without it. Recorded as the named refinement.
