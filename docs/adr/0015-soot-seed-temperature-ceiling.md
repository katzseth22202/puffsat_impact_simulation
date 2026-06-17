# The dark-oil seed works only below soot's sublimation ceiling — its niche is the cool-to-low-transition window, not the hot/transitional core

The dark-oil opacity seed (§7) pyrolyzes to soot, and **soot is condensed carbon with a
temperature ceiling.** That ceiling, not the radiative loss budget, is what bounds where the lever
can act. This inverts the §7 draft's expectation that the "bigger payoff" was at the
transitional/high-v end.

**The ceiling.** Carbon sublimes at ~3915 K at 1 atm; its triple point is ~10 MPa / ~4600 K, so at
the stagnation pressures here (tens–hundreds of MPa) condensed carbon survives to ~5000–5500 K and
no further. Above that the soot is vaporized and then ionized — it becomes a minor low-Z additive to
the water plasma's *own* opacity, not a Rayleigh particle absorber.

**Consequence — where the lever is and is not real:**
- **Hot end (16 km/s) and transitional core (10–20 kK):** no soot exists. And the high-v gas is
  already `τ≫1` on its own (Orion). The seed adds essentially nothing here. The §7 "bigger payoff at
  the transitional/high-v end" framing is **dropped.**
- **Cool anchor (3.2 km/s, ~1700 K):** soot survives everywhere in the stagnated gas, but this is
  exactly where total radiative loss is tiny (flux ~`v⁸`, ~5×10⁻⁵ of the 11 km/s level). Making an
  already-negligible channel optically thick changes a small number to another small number.
- **Real niche:** the window where soot still survives *and* radiation is climbing —
  roughly **`v ≲ 5 km/s` (T ≲ ~5000 K), i.e. the cool anchor up to the low edge of the transition
  (~4–5 km/s),** where radiative loss is ~10⁻³ of the high-v level. Plus a secondary role trapping
  **outward/sideways** radiation as it passes through the cooler peripheral/upstream gas surrounding
  a hotter core, where the particles survive.

**Ties.** This is the bulk-gas cousin of the ablating wall's near-wall carbon-bearing absorber
(ADR-0014): same Rayleigh-carbon physics, sourced in the cloud rather than at the wall, and bound by
the same condensed-carbon temperature ceiling in the hot near-wall layer. The seed's `f`-recovery
role remains gated on the loss decomposition (ADR-0009/0013); this ADR narrows *where* that gate can
ever open.

**Oxidation narrows the window further, from the top.** The sublimation ceiling (~5000-5500 K) is an
upper bound; soot consumption by O/OH reaches it first. As `v` climbs toward the ceiling, water
dissociates and floods the gas with O/OH, so the soot oxidizes within the pulse — likely pulling the
top of the *usable* window down to ~4 km/s. The seed model therefore carries a two-clock lifetime
check (inception time < pulse; oxidation + sublimation consumption time > pulse) as a precondition;
see §7. Feedstock stays dark oil (dispersion + biocompatibility), with a pre-sooted/heavy-aromatic
fallback if the inception clock is too slow.

## Considered Options

- **Seed targets the transitional/high-v end (original §7).** Rejected: soot cannot exist above
  ~5000–5500 K, and the high-v gas is already optically thick without it.
- **Carbon-plasma line/continuum opacity at the hot end.** Rejected as a meaningful lever: at 0.5%
  mass it is a minor perturbation to the water plasma's own ionization opacity, which already gives
  `τ≫1`.
