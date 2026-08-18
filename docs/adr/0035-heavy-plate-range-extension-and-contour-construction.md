# The heavy plate runs to 63 km/s, and `f(v)` is read off a contour rather than a shape grid

The heavy-plate scenario (design §12.1) was swept 16–28 km/s. The question that opened this pass —
*what is the likely fudge factor at 45–63 km/s, and is the physics we built adequate there?* — is
answered by extending the same scenario on unchanged knobs (100 kg pulse, `R = 15 m`, ≤ 40 t plate,
water) rather than standing up a separate high-velocity study. One continuous curve, not an island.

**`f` stays near 0.8 across the whole range.** On the contour construction below it holds
**0.79–0.82** from 16 to 63 km/s, and the 16 km/s point reproduces the published 0.811, so the
extension is consistent with what it extends.

| v [km/s] | 16 | 22 | 28 | 34 | 40 | 45 | 50 | 55 | 63 |
|---|---|---|---|---|---|---|---|---|---|
| ρ contour [kg/m³] | 0.582 | 0.582 | 0.427 | 0.287 | 0.206 | 0.160 | 0.130 | 0.106 | 0.081 |
| binds | box | box | surv | surv | surv | surv | surv | surv | surv |
| **f** | 0.794 | 0.816 | 0.813 | 0.813 | 0.816 | 0.818 | 0.819 | 0.817 | **0.809** |

## The velocity grid has two legs, because the halves are known to different degrees

3 km/s over 16–43 (10 anchors) and 1 km/s over 45–63 (19), for 29. The lower half was already swept
and varies smoothly there, so it is reconnaissance; the upper half has never been sampled and its
`τ` structure was unknown, so it gets the resolution. The 43 → 45 seam is a deliberate 2 km/s
join — finer than the coarse leg, coarser than the fine one — so nothing is lost across it.

Density likewise runs 12 non-uniform points, dense through the steep on-contour band 0.06–0.20 where
`ρ_ceiling(v)` actually lives, coarse across the flat top where `e_eff` barely moves.

## `f(v)` is read off the contour, not off the shape grid

The original construction evaluated 27 fixed cloud shapes per velocity and kept the best survivor.
That let the **shape grid** decide where the answer steps rather than the physics: across 45 → 46
km/s the optimum jumped `r_foot/R` 0.5 → 0.7 and `f` fell 0.817 → 0.791, purely because no
intermediate footprint existed to be chosen.

The contour is the curve the schedule actually flies:

```
ρ(v) = min(ρ_ceiling(v), ρ_max),    ρ_ceiling = P_limit / (c_stag v²)
```

solved back through the Σ contract (ADR-0003) for the shape that lands on it. The same crossing
becomes smooth — `Δ(r_foot/R)` 0.200 → 0.007, and the step in the reported quantity 0.0259 → 0.0003.
The 27 discrete shapes become validation points *on* the contour rather than the only places it may
be sampled.

**Which limit binds is reported, not smoothed over.** Below ~25 km/s the shape box binds — the cloud
cannot be built dense enough to reach the survivable ceiling. Above it survivability binds. Those are
different engineering situations with different levers, and collapsing them into one number hides
which one you are up against.

Fixing `ρ` is one equation in two shape parameters, so the point is **chosen**, not determined:
`e_eff` is already pinned by `ρ`, so the choice maximizes `eta_capture` along the remaining
one-parameter family. It lands at the shortest feasible `L/D` every time (0.300–0.302 across the
whole range), which is why the discrete construction kept selecting `L/D = 0.3`.

**The `eta_capture` interpolation is claimed only along this contour.** Across the full shape box
`eta` spans 0.79–0.99; along the schedule it spans 0.972–0.978. A bilinear fit over a 3×3×3 sweep is
adequate for the second and would be advertising unearned precision in the first.

## The cloud schedule is worth 0.017 in `f`

Design §7 treats the cloud as a per-shot schedule `shape(v)`, available because the plate is one
built object and the cloud is not. The pinned-shape companion asks what happens without that freedom:
one shape flown everywhere, forced by the most demanding velocity since it must survive there.

The worst gap across 16–63 km/s is **0.017 in `f`**, against a freeze-timing band of ~0.05. The
schedule buys less than the dominant uncertainty. Worth knowing before anyone builds per-shot cloud
shaping to recover it.

## What the sweep may report about itself

Two contracts settled here because they are the same decision — what a row is allowed to claim.

**The Jupiter table's ρ grid was extended** node-preservingly to 1000 kg/m³ (48 → 62 nodes, ceiling
30 → 1284 kg/m³). `build_table_jupiter` never received `build_table`'s `extended_rho_grid` fix, so a
radiatively-cooled wall cell compressed past the table ceiling, where the clamped `p(ρ)` stops
arresting Lagrangian compression — the radiative collapse of `4ddaed5`, still live. The first N nodes
are the original array elements bit-for-bit, so the preserved block is unchanged. The extension is a
*numerical* provision for the compressed wall cell, not a claim the scenario samples those densities,
and stitched opacity clamps above the TOPS coverage ceiling.

**A stalled bounce now says so.** `BounceResult.converged` records whether the run reached the
10⁻³-of-peak tail guard or exhausted its step budget. The contract for consumers is that
**`!converged` means no result, not a low one** — a stalled run reports a truncated impulse integral,
so reading its `e_eff` as a low restitution is reading a solver artifact as physics. The analysis
rejects such rows rather than averaging them in.

This was not hypothetical. One row (28 km/s, ρ = 0.6, κ = 10×) was stalled at `e_eff = 0.0562`
against neighbours at 0.6749–0.6783, and it **passed the regression gate** because the baseline
carried the identical corruption. Old-vs-new comparison cannot see physical implausibility. With the
grid extended it returns 0.6799 and fits the monotone trend.
