# A jet-delivered slug separates column density from bulk density: Arm M3

The magnetic arms want two things from the slug that have been treated as one quantity, and read as
opposed. §6.2 requires **column density** — mass per unit area along the projectile's path — or the
projectile punches through and the slug becomes a spectator. §6.9 wants **low bulk density**, because
the temperature reached at a given radius follows the expansion adiabat from the initial density, and
that temperature is what decides `r_σ`, the outer edge of the nozzle window. Arm M2 resolves the
tension radially, with a dense anvil backed by snow, and pays the entropy cost of the anvil driving a
strong shock into it (§4.2).

**They are not one quantity.** Column density depends only on the slug's *radius*; bulk density
depends on radius *and length*. A slug elongated along the collision axis at fixed radius lowers bulk
density while leaving column density untouched. At 60 kg:

| jet diameter | column density | retained speed (0.2 m rod, 421 kg/m²) | energy deposited |
|---|---|---|---|
| 0.20 m | 1910 kg/m² | 0.181 | 96.7% |
| 0.42 m | 433 kg/m² | 0.493 | 75.7% |
| 1.0 m | 76 kg/m² | 0.846 | 28.4% |
| 2.0 m | 19 kg/m² | 0.957 | 8.5% |

| length at 0.42 m diameter | bulk density | vs solid ice | column density |
|---|---|---|---|
| 1 m | 433 kg/m³ | 2.1× below | **433 kg/m² in every row** |
| 10 m | 43 kg/m³ | 21× below | " |
| 30 m | 14 kg/m³ | 64× below | " |

The reference solid ice slug is 531.7 kg/m² of column at 917 kg/m³ bulk, so a 0.42 m × 10 m jet keeps
81% of the coupling at 1/21 the density. **Phase is irrelevant to column density** — 60 kg inside a
0.42 m column is 433 kg/m² as steam, as crystals, or as solid ice. What killed the first version of
this idea was geometry, not phase: the same 60 kg as a free-expanding 3 m sphere is 3.2 kg/m², which
is 167× short and leaves the projectile 99.3% of its speed.

**Mean free path is not a third lever and does not enter.** At every density in this study's range it
is ≲1 µm — 0.28 µm in the 3 m cloud, `Kn` ≈ 10⁻⁷ — and reaching `Kn` = 1 would need densities 10⁷×
lower. The binding constraint is the §6.2 snowplow, which is momentum conservation and is indifferent
to collisionality.

**What the elongation buys, and what it does not.** It does *not* move `r_β`: once expanded past its
initial size the cloud has ρ(r) = M/(4/3·πr³), set by mass and radius, so §4.2's claim that a larger
initial volume moves `r_β` inward is a near-field effect and weaker than it reads. It moves `r_σ`
outward, via `T = T₀(ρ/ρ₀)^(γ-1)` — 1.84× at γ = 1.2, 3.40× at γ = 1.4 — and `r_σ` is the side that
binds (§6.9: the window "closes hard below ~2500 K").

**It nearly, but probably not quite, removes the alkali seed.** Saha at ρ = 10⁻² kg/m³, calibrated so
the seeded row reproduces §6.9's table:

| T | unseeded `Rm` | seeded `Rm` (1 wt% K) |
|---|---|---|
| 2500 K | 7×10⁻⁸ | 5.8 |
| 4000 K | 0.014 | 210 |
| 5000 K | 0.84 | 303 |
| 6000 K | 13.4 | 326 |

Unseeded water needs ~5000–6000 K where seed works at 2500 K — a 2.2–2.4× requirement against the
1.84–3.40× elongation supplies. **It lands on the knife edge and the deciding parameter is the
effective γ**, which for a dissociating, ionising plasma is low, so the pessimistic branch is the
likely one. This is ADR-0026's frozen-vs-equilibrium bracket, which §6.9 already notes cuts
favourably for a magnetic nozzle. The seed is retained; what changes is that it may be margin rather
than load-bearing.

**The delivery nozzle is load-bearing for three separate reasons, and only their intersection works.**
The magnetic arms carry no ablator (§4.2), so the coil has no mass-transfer sink at all and the
parasitic eddy and unreflected-radiation load must leave with the propellant. But holding a
0.42 m × 10 m column requires 2.1 mrad of divergence, and a warm jet destroys that on its own vapour
pressure — radial blowoff `a ≈ (p_sat/2)/(ρ_bulk·r)` gives Δv_r = 561 m/s at 373 K over a 100 ms
flight, against a 0.21 m/s budget, and only ≲230 K is quiescent. A jet cold enough to stay collimated
carries 14–38 MW of sensible heat, against the 259 MW the load needs. **Cold cannot carry the heat;
warm cannot hold the jet. Only a de Laval expansion that takes the water in hot and delivers it out
cold satisfies both** — converting the thermal load into directed axial kinetic energy that departs
with the mass. This also voids §6.5.3's prohibition on boiling ("phase change is one-way without a
radiator"), which is scoped to an *assembled-body* slug that must recondense on board; a jet that
recondenses during nozzle expansion throws the latent heat overboard instead, raising the usable sink
from 259 MW to ~800 MW.

## Consequences

- **Arm M3 is added** (§4.2, §14): seeded water boiled against the coil's eddy shield, expanded
  through a de Laval nozzle, delivered as a collimated cold-crystal column. **M1 is its
  aspect-ratio-1 limit**, so M3 is a generalisation of M1 rather than a competitor to it.
- **Rung 1A-i still leads with M1**, not M3. The gate asks whether the field couples to this plasma at
  all, and M3 would add two-phase nucleation kinetics, a nozzle contour, millisecond collimation, and
  a burst feed system to a screen whose value is that it is clean. §4.2 selects M1 precisely because
  "initial mixing is the least uncertain of any arm."
- **But the 1A-i gate must report quantitatively, not as a binary.** §6.9's "a closed window kills
  both" predates M3. A window that misses by less than the ~1.8× elongation supplies leaves M3 live;
  one that misses by 10× kills the architecture. The rung must therefore report how far `r_β` and
  `r_σ` are from overlapping, in temperature at radius — otherwise a marginal M1 result is a false
  negative for the magnetic arms as a class.
- **EOS/opacity tables are generated with `(ρ, T)` margin for the elongated case.** An elongated slug
  reaches ~1.8× higher T at the same ρ. Table coverage is the one thing here that is not a cheap
  re-run (§7.3 already lists the dense-ice → vapour-plasma handoff as missing).
- **The eddy shield and burst feed system are dry mass, not charged mass** — reported beside Isp, not
  inside it, per §3.1 and the existing "charged in Isp; the coil is not" convention.
- **Achievable nozzle exit divergence now prices this whole branch**: it caps aspect ratio → caps the
  bulk-density reduction → caps the `r_σ` gain → decides whether the seed is margin or load-bearing.
  At a realistic 5–10 mrad the usable length is ~2–4 m and the temperature gain falls to ~1.35–1.55×.
  2.1 mrad is a stretch target, not an assumption.
- **Condensation is assumed to occur inside the nozzle and this is not yet justified.** Condensation
  in rapidly expanding flow is kinetically limited; late condensation releases latent heat into
  already-diverging flow and spoils the collimation being bought. Equilibrium condensation is the
  optimistic branch and requires a nucleation calculation, not an assumption.

## Considered Options

- **Deliver the slug pre-vaporised as a free-expanding cloud** (the original proposal). Rejected on
  §6.2: a 3 m cloud is 3.2 kg/m², 167× short, and deposits 1.5% of the projectile's energy. The slug
  becomes the spectator §6.2 warns about, and in a *pure* magnetic arm there is no plate to catch what
  it misses, so the event produces no fireball at all.
- **Deliver a warm liquid or slush jet.** Rejected on collimation: self-driven vapour-pressure blowoff
  exceeds the budget by 100× at 100 ms of flight. This also reverses the intuition that a slow, gentle
  jet is the simpler engineering choice — it is the one configuration that cannot exist.
- **Take M2's route and buy initial volume radially, with snow backing.** Retained as Arm M2, not
  superseded. It costs column density where elongation does not, but it is an assembled body with no
  feed system, no nozzle contour, and no condensation kinetics — a materially lower-risk way to buy
  some of the same effect.
- **Cool the coil by running water over the copper stabiliser.** Not available: that layer is part of
  the cold mass, bonded to the superconductor at 4–77 K, and its function is to conduct heat *to* the
  cryogen during a local quench. Water at 273 K on it is a quench, and the water freezes. The loop
  belongs on a warm normal-conducting eddy shield standing off the cryostat — which §6.9 already
  requires as "a radiation shield with its own thermal path" and which, with an L/R longer than the
  750 µs pulse, also excludes the plasma's AC perturbation from the winding pack deliberately rather
  than by accident.
