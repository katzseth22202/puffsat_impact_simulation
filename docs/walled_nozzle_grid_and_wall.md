# N16's conversion surface and N9's load case — second return

Answers to `katzseth22202/Balloon-Pulse-Propulsion` @ `dd6ad2a`,
`docs/nozzle_asks_for_impact_sim.md`, worked in
`katzseth22202/puffsat_impact_simulation` on 2026-09-10. Scope is
[ADR-0051](adr/0051-walled-nozzle-grid-geometry-and-the-load-case.md), extending
[ADR-0050](adr/0050-walled-thermal-nozzle-methane-chamber.md).

**This continues the W series**, so the items below are **W10 onward** and the earlier document
[`walled_nozzle_asks_answered.md`](walled_nozzle_asks_answered.md) (W1–W9) still stands except
where an item here says otherwise.

**This document is written to be copied into the paper repository and worked there.** Every
number is stated in full, because an agent in the paper repository cannot run the code that
produced it. Every figure came from a run; the [provenance table](#provenance) says which.

---

## Status

| item | verdict | does the paper change? |
| --- | --- | --- |
| **N9 items 1–2** (contact station, state there) | **Answered, and the temperature premise does not survive.** The front reaches the wall at 5.6–6.2 m — the ask's own "6 m" is right — but by then it has ploughed through 41–178 kg and slowed to 9–28 km/s, so the layer is at **4,900–30,700 K, not 200,000 K**. | **yes** (W10, W11) |
| **the gate** (N9's own precondition) | **Answered, and the paper side's guess is right to 5%.** The front never reaches the wall below **178 m³**, so items 1–4 **close for the 50 and 100 m³ chambers**. On the fast end of the spreading bracket the gate falls to 92 m³. | **yes** (W10) |
| **the contraction** (raised here, not asked) | **Priced, and it does not shrink with the chamber.** The bulk load on the convergent section is **volume-independent** — the same 70.3 GJ turns the same corner — and the arrival load gets *worse* as the chamber shortens. | **yes** (W23) |
| **N9 item 3** (fluence, split) | **Answered, and it reverses which channel matters.** Radiative delivery is **0.001–0.009 MJ/m²** against the 1.6 MJ/m² the section already clears. The load is **convective**, and it *grows* with chamber volume. | **yes** (W12) |
| **N9 items 4 and 6** (the film) | **Answered, and generously.** Even 0.02 kg/m² carries 2.8 MJ/m² of capacity — more than the whole radiative transient at any volume. A sprayed film is not ruled out here the way it is for the magnetic liner. | **yes** (W13) |
| **N9 item 7** (convective flux) | **Answered.** The paper's Bartz-like 123–281 MW/m² is a good *liner* number (we get 93–174 at 10 kK and 200 m³) and roughly 2–3× low at the *throat*. | **yes** (W14) |
| **N9 item 5** (throat carbon) | **Answered, and there is no self-healing to be had.** The exhaust is **undersaturated** in carbon at any wall temperature graphite survives, so it *erodes* the throat rather than plating it. The boundary layer settles the fork with 5 decades to spare. | **yes** (W15) |
| **N16** — `held(T, ρ)` | **Answered, and it retires the fit.** The exponent is **not −0.21 and not constant**: it runs 0 to −0.25 for methane and 0 to −0.42 for water, vanishing at both ends. | **yes** (W16) |
| **N16** — the temperature dial | **Answered, and it reverses.** **Hotter is better, monotonically.** Across the grid 6,000 K returns 497–743 s effective where 12,000 K returns 800–1,156 s effective. | **yes** (W17) |
| **N16** — the volume dial | **Answered, and "free" is wrong.** Volume buys slug ratio and costs freeze margin; the optimum is interior and moves with temperature. | **yes** (W18) |
| **N16** — the throat dial | **Answered, and there is a wall.** The chemistry **does** freeze past `A/A*` ≈ 140. The store returned **peaks near `A/A*` = 56 and then falls**. | **yes** (W19, W20) |
| **the throat's price** (raised here, not asked) | **Priced, and the paper states no replacement interval.** Narrowing from 7 m² buys **+22% at 2 m²** and +31% at 1 m², but throat life goes as `A*^1.6` — **a factor of three per halving** — because the cost is blowdown *dwell*, not peak flux. | **yes** (W24) |
| **N16** — water on the same grid | **Answered, and it loses everywhere.** Water returns **52–83%** of methane's effective impulse at every one of the 128 grid cells, including at 400 m³ and at the cold end. | **yes** (W21) |
| **N16** — the named configuration | **Priced, with plain methane beside it, and the hydrogen fraction swept.** Methane beats the 90% water + 10% H₂ mix by 4–12% at the ask's own point. Hydrogen is a real dial worth +18–23% over pure water. | **yes** (W22) |
| **N13** (exit-plane velocity distribution) | **Not started, and confirmed unreachable from here.** It needs a 2-D solve; see [Deferred](#deferred-with-cost). | no |
| **N11** (radiative escape, carbon-bearing plume) | **Not started.** | no |

**Fifteen items fall out, W10–W24.** **If you read only four:** **W17** (the temperature dial
runs the other way, which unwinds the whole "Going hotter was considered and declined" section),
**W24** (narrowing the throat is worth +22 to +35%, priced in throat life at a factor of three per
halving — and the paper's stated 8 → 28 ms cost corresponds to a 2 m² throat and no deeper),
**W10** (items 1–4 close for the chamber you were going to build), and **W16** (the `ρ^-0.21`
hack, replaced by the computed surface). **W19** remains the one to read before quoting any deep
throat.

**Eight of the ask's stated premises did not survive**, three of them paper-side beliefs that
reverse outright; they are collected in
[What the ask got wrong](#what-the-ask-got-wrong), with what the ask got right beside them.

---

## What the paper should change

### W10. The front never reaches the wall below 178 m³, so N9 items 1–4 close for the chamber you want to build

**Locate:** N9's "Before designing any mitigation, price the strike", and its "if the front exits
first below about 170 m³, say so and items 1 to 4 close for the short chamber".

**Now:** the section carries the strike as an open condition at every volume.

**Should be:** the gate is real and the paper side's estimate is right to 5%.

| spreading rate | gate volume | gate column length |
| --- | ---: | ---: |
| **1.0 × the shocked layer's sound speed** | **178 m³** | 6.28 m |
| 1.6 × | 109 m³ | 3.87 m |
| 1.9 × | 92 m³ | 3.25 m |

At the conservative (sound-speed) reading the **50 m³ and 100 m³ chambers are both safe**: the
front leaves the column at a radius of 1.19 m and 2.02 m against a 3 m wall. At the fast end of
the bracket the 100 m³ chamber becomes marginal (gate 92 m³), so the honest statement is
**"safe at 50 m³; safe at 100 m³ unless the front spreads faster than its own sound speed"**.

**The bracket is the same one `sec:needle_through_fog` already uses** (1.0/1.6/1.9 × the shocked
layer's sound speed), so this is not a new uncertainty — it is the existing one, evaluated.

**And it is not a free pass.** The front leaving the *column* is not the front leaving the
*hardware*. At a 0.05–1 m² throat the convergent section is 96–99% of the bore area — very nearly
a flat end wall — so a short chamber trades a **grazing side-wall strike for a normal-incidence
strike on the contraction**. That is priced in [W23](#w23-shrinking-the-chamber-moves-the-strike-to-the-contraction-and-the-arrival-term-gets-worse-not-better),
and it does not go the way the gate does.

---

### W11. Where the front does reach the wall, the layer is at 4,900–30,700 K — not 200,000 K

**Locate:** N9's "a first estimate puts it near 200,000 K", and `sec:needle_through_fog`'s
94,600 K arrival.

**Now:** the wall's load case is stated against a layer near 200,000 K.

**Should be:** 200,000 K is right for the *entering* front and wrong for the *arriving* one. At
75 km/s into methane the freshly shocked layer really is at **166,000–205,000 K** across the
candidate pre-charge densities (50 to 673 m³) — the strong-shock piston solution `e = v²/2` on `eos_methane`
reproduces the ask's estimate closely. But the front decelerates hard, because the chamber is
0.7 to 11 kg/m³ rather than the magnetic bag's 0.32.

| chamber | column | pre-charge ρ | cone half-angle | contact `x` | arrival `t` | speed there | swept | local `k` | **layer T** |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 m³ | 1.77 m | 11.04 kg/m³ | 30.3° | **never** (leaves at r = 1.19 m) | — | — | — | — | — |
| 100 m³ | 3.54 m | 5.12 kg/m³ | 29.7° | **never** (leaves at r = 2.02 m) | — | — | — | — | — |
| 200 m³ | 7.07 m | 2.45 kg/m³ | 29.3° | **6.15 m** | 249 µs | 9.2 km/s | 178 kg | 7.13 | **4,867 K** |
| 400 m³ | 14.15 m | 1.19 kg/m³ | 28.9° | **5.64 m** | 135 µs | 19.2 km/s | 72.5 kg | 2.90 | **17,186 K** |
| 673 m³ | 23.80 m | 0.70 kg/m³ | 28.7° | **5.57 m** | 108 µs | 28.2 km/s | 41.5 kg | 1.66 | **30,726 K** |

**Two corrections inside this.** The cone is **28.7–30.3°**, not the 24.9° ADR-0016 carries —
that angle belongs to water at 45.58 km/s, and methane at 75 km/s opens wider. And the contact
station is **essentially independent of chamber volume** at 5.6–6.2 m, because a denser
pre-charge decelerates the front (which widens the cone) exactly as fast as it shortens the
column. ADR-0016's own "6 m of the column" is therefore right, and right for a reason it did
not state.

**The temperature is what moves.** In the 200 m³ chamber the front has swept seven times its own
mass before it touches, so the layer at contact is 37× cooler than the layer at entry. A wall
argument made against 200,000 K is being made against a state that never sees the wall.

---

### W12. The strike is radiatively negligible; the load is convective, and it gets worse with a bigger chamber

**Locate:** N9 item 3, "fluence delivered to the wall over the contact transient, in MJ/m²,
separated into radiative and convective parts".

**Now:** unpriced, with the radiative channel implicitly the worry (it is what `tau ≈ 12.5` and
`54.8 MW/m²` are about).

**Should be:** the radiative channel is two to three orders of magnitude below the equilibrium
load, and the convective channel is the whole answer.

| chamber | Rosseland `τ` | regime | `σT⁴` | **radiative** | **convective** (St = 0.001–0.01) | ceiling |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 200 m³ | 0.30 | emission | 3.2e7 W/m² | **0.0032** | **0.11 – 1.09** | 625 |
| 400 m³ | 6,560 | diffusion | 4.9e9 W/m² | **0.0009** | **3.15 – 31.5** | 7,514 |
| 673 m³ | 9,920 | diffusion | 5.1e10 W/m² | **0.0088** | **8.79 – 87.9** | 20,387 |

All fluences in MJ/m²; the reference the section already clears the wall against is **1.6 MJ/m²
per pulse**.

**Why radiation loses.** At 673 m³ the layer is genuinely hot, but it is also optically thick by
four orders of magnitude, so photons random-walk out and the flux available is `4σT⁴/3τ`, not
`σT⁴`. At 200 m³ the layer is thin but only 4,867 K, so even the *blackbody* cap — 3.2e7 W/m²
over a 100 µs transient, 0.0032 MJ/m² — is 500× below the equilibrium load. **The verdict does
not depend on the opacity model**, because at one end the blackbody cap settles it and at the
other end more opacity means less delivery:

| opacity scale at 673 m³ | `τ` | radiative fluence |
| --- | ---: | ---: |
| 0.1 × | 992 | 0.0878 MJ/m² |
| 1 × | 9,921 | 0.0088 MJ/m² |
| 10 × | 99,210 | 0.0009 MJ/m² |

**The convective number is the crude one and it is the one that binds.** What is computed is the
stagnation enthalpy the layer *brings* to the wall — 109, 3,152 and 8,793 MJ/m² of incident flux
at the three volumes — times a Stanton number. The bracket 0.001–0.01 is the Reynolds analogy on
a turbulent boundary layer at Re = 1e6–1e8. **A tighter number needs a boundary-layer solve** and
is the largest single uncertainty in this document.

**The direction is unambiguous even so.** Convective delivery scales as ρ_shock × u_radial × e ×
Δt, and every one of those grows as the chamber gets larger and thinner. **200 m³ is 0.1–1.1
MJ/m², comparable to the equilibrium load; 673 m³ is 8.8–88 MJ/m², 5–55× worse.** So the wall
argument and the Isp argument now point the same way: **make the chamber small.**

---

### W13. The film has an order of magnitude more capacity than the transient needs

**Locate:** N9 item 4 ("whether a continuously injected cold film at 0.02 to 0.2 kg/m²
measurably reduces (3)"), and item 6 (whether a sprayed film survives the front).

**Now:** the film is proposed and unpriced, and `sec:watering_it_down` rules a sprayed film out
for the magnetic liner.

**Should be:** a cold methane film's capacity is its areal density times the specific energy a
kilogram of methane absorbs at chamber conditions — atomisation plus sensible heat, 140.5 MJ/kg:

| film | capacity |
| ---: | ---: |
| 0.02 kg/m² | 2.8 MJ/m² |
| 0.05 kg/m² | 7.0 MJ/m² |
| 0.10 kg/m² | 14.1 MJ/m² |
| 0.20 kg/m² | 28.1 MJ/m² |

**Against the radiative transient (≤ 0.009 MJ/m²) even the thinnest film is 300× oversized**, so
the answer to item 4 as posed is *the film is not the lever, because radiation is not the load*.
Against the convective transient the film is sized correctly at 200 m³ (2.8–28 MJ/m² of capacity
against 0.11–1.09 MJ/m² of load) and undersized at 673 m³ on the pessimistic Stanton edge
(28 MJ/m² against 88).

**Item 6 — the sprayed film — is not ruled out here, and `sec:watering_it_down`'s rejection does
not carry over.** That rejection is about a liner facing a plume at 45.58 km/s; this liner is
grazed by a layer whose radial arrival speed is 2.6–12.2 km/s and whose radiative load is
negligible. The film's problem in the walled case is mechanical scouring, not radiative
vaporisation, and nothing here prices scouring.

---

### W14. Bartz at the liner matches the paper; the throat is 3.5× the liner and is the new item

**Locate:** N9 item 7, "the paper-side estimate is a Bartz-like scaling off a single anchor,
giving 123 to 281 MW/m² between 10,000 and 18,000 K against a radiative 16 to 75".

**Now:** the estimate carries an unbounded uncertainty and is what sets the chamber temperature
ceiling.

**Should be:** the paper's number is a **good liner number** and should be labelled as one.

| chamber | `T_c` | `p_c` | liner (bore) | throat | throat surface `T` at radiative balance |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 200 m³ | 10,000 K | 642 bar | **93 – 174** | 328 – 610 | 10,183 K |
| 200 m³ | 12,000 K | 723 bar | 137 – 197 | 482 – 692 | 10,511 K |
| 200 m³ | 15,000 K | 815 bar | 198 – 351 | 696 – 1,234 | 12,147 K |
| 400 m³ | 10,000 K | 321 bar | 55 – 81 | 193 – 283 | 8,406 K |
| 400 m³ | 15,000 K | 401 bar | 112 – 228 | 393 – 801 | 10,900 K |
| 673 m³ | 10,000 K | 190 bar | 36 – 49 | 128 – 173 | 7,436 K |
| 673 m³ | 15,000 K | 235 bar | 73 – 168 | 256 – 590 | 10,099 K |

All fluxes MW/m². The two numbers in each cell are the **frozen** and **equilibrium** `c_p`
edges; the factor between them (1.9–2.2×) is the honest width of a Bartz answer on a dissociating
gas, and it is the largest term in the correlation's uncertainty.

**The paper's 123–281 MW/m² overlaps the liner column across 10–15 kK** — it is a good match at
200 m³ and runs high at 673 m³. So the
conclusion that convection is 80–90% of the wall load, and that it rather than radiation sets the
temperature ceiling, **survives** — with the qualification that the underlying Bartz correlation
was fitted on chemical rockets at tens of bar and a few thousand kelvin, and is here extrapolated
to a thousand bar and a partly ionised boundary layer. **Read the size, not the digits.**

**What is new is the throat, and it is the item.** At the flown 7 m² throat the flux is 2–4× the
liner's. At the throats N9 and N16 now recommend it becomes the design's binding component:

| throat | flux (frozen – equilibrium) | blowdown | fluence per pulse | graphite ablated |
| ---: | ---: | ---: | ---: | ---: |
| 7.00 m² | 302 – 1,445 MW/m² | 5 ms | 8 MJ/m² | 0.06 mm |
| 2.00 m² | 342 – 1,638 | 19 ms | 30 | 0.23 mm |
| 0.50 m² | 393 – 1,882 | 74 ms | 140 | 1.07 mm |
| **0.20 m²** | 431 – 2,062 | 185 ms | 382 | **2.94 mm** |
| **0.10 m²** | 462 – 2,210 | 371 ms | 820 | **6.29 mm** |
| 0.05 m² | 495 – 2,369 | 742 ms | 1,757 | 13.49 mm |

(100 m³ at 8,000 K; ablation at the equilibrium-`c_p` edge with no transpiration credit, so read
it as the pessimistic end of a factor-of-two band.)

**The flux barely moves with throat area — the fluence moves by 200×**, because the blowdown
lengthens in proportion. That is the cost the paper has not booked for narrowing the throat: the
throat is exposed for 100× longer per pulse.

**And a throat under this load is not a passive structure.** The surface temperature at radiative
balance is 7,400–12,100 K, far above graphite's 3,900 K working ceiling. The throat ablates,
and the ablation is millimetres per pulse at the recommended throats.

---

### W15. The throat does not self-heal — and it is not plating at all, it is being chemically eroded

**Locate:** N9 item 5, "the throat is the hottest and fastest station so it should self-clean,
but nothing shows it", and "0.34 to 3.7% redeposition suffices" for the liner.

**Now:** self-cleaning is assumed on the grounds that the throat is hot, and the liner is assumed
to be self-healing if a small share of the exhaust carbon redeposits.

**Should be:** neither. The exhaust does not give carbon to a graphite surface at any temperature
that surface can survive, so there is nothing to self-heal *with*.

**Three states have to be kept apart, and conflating any two gets the wrong answer:**

| state | where | condition at the flown throat |
| --- | --- | --- |
| the free stream | mid-channel | 7,244 K, 623 bar, 23% of carbon as free atoms |
| **the wall-adjacent gas** | at the surface | **623 bar** (pressure is constant across a boundary layer) at **the wall's temperature** |
| the surface | graphite | 3,900 K, its working ceiling |

Deposition is driven by the **second against the third**, not the first against the third.
Pressure does not vary across a boundary layer and the gas at a no-slip wall is at the wall
temperature, so the state that decides whether carbon sticks is the free stream's *pressure* at
the wall's *temperature* — a third state, denser than either end (15.1 kg/m³ against the free
stream's 4.4).

**And at that state the gas is undersaturated, because the carbon is in acetylene.** Saturation
ratio `S = p_C / p_vap(graphite)`, at 623 bar:

| wall `T` | 1,500 K | 1,800 K | **1,912 K** | 2,000 K | 2,500 K | 3,000 K | **3,900 K** | 4,500 K |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `S` | 6.13 | 1.52 | **1.00** | 0.750 | 0.210 | 0.0886 | **0.0316** | 0.0182 |

`S > 1` plates; `S < 1` means the gas can dissolve *more* carbon than it is carrying. **At
graphite's own working temperature `S` = 0.032 — the gas has 32× spare carbon capacity, so it
takes graphite off the surface rather than putting it on.** The band where nothing plates runs
from **1,912 K upward**, with no upper edge out to 8,000 K; below the cold edge everything is
supersaturated, which is just "solid carbon is the stable phase down there".

**The mechanism is the same one W9 found from the recovery side.** In a hydrogen-rich carbon gas
between about 1,900 and 5,000 K the carbon parks itself in **C₂H₂** — 68% of carbon nuclei at
2,500 K (the rest as CH₄) and 92% at 3,900 K — and that sink is deep enough to hold the
monatomic-carbon activity
below graphite's own. It is why hydrogen-rich flames do not soot, and it is why a graphite throat
inside this window cannot collect a deposit. **Nearly all the throat's carbon is in acetylene,
not C₃**, so the C₃ partition function that is `eos_methane`'s worst admitted weakness (W9,
ADR-0050 weakness 1) does **not** reach this answer: C₃ carries under 3% of the carbon at
graphite's working temperature and under 9% anywhere below 4,300 K. (It does take over above
4,700 K — 20% at 4,700 and 46% at 5,500 — so a *non-carbon* refractory running hotter than
graphite would put the question back inside the EOS's weak spot.)

**This turns on the boundary layer being in chemical equilibrium, and it is, with five decades to
spare.** If the layer were frozen at the free stream's composition, free carbon atoms would arrive
at the wall and `S` would be **41** instead of 0.032 — a six-order-of-magnitude fork, and the only
thing that decides it:

| throat | boundary layer | diffusion clock | gas-kinetic clock | margin |
| ---: | ---: | ---: | ---: | ---: |
| 2.0 m² | 18.7 µm | 4.05 µs | ~24 ps | **5.3 decades** |
| 0.5 m² | 16.3 µm | 3.07 µs | ~24 ps | **5.2 decades** |
| 0.2 m² | 14.9 µm | 2.56 µs | ~24 ps | **5.1 decades** |
| 0.1 m² | 13.9 µm | 2.23 µs | ~24 ps | **5.0 decades** |

The carbon-to-acetylene chemistry would have to run **more than five decades below a gas-kinetic
collision frequency** to freeze inside that layer. These are radical reactions with small
barriers at 3,000–4,000 K; nothing plausible is five decades slow. **No rate coefficient is
asserted here** — `rates.py`'s provenance rule forbids inventing one — and the margin is stated
as "how far below the collision frequency the chemistry could run and still finish", which is a
bound on what would have to be true to flip the verdict rather than a claim about any reaction.

**The boundary-layer thickness is free of the largest Bartz uncertainty.** `δ = k/h` with
`k ∝ μ c_p` and Bartz's `h ∝ c_p`, so `c_p` cancels: the frozen and equilibrium branches, which
differ by a factor of two in *flux*, give identical thicknesses. What is left depends on `μ` and
`Pr`, and the margin is 5 decades, so it survives an order of magnitude on either.

**So the throat's ledger has no credit side:**

| throat | ablated per pulse (convective, no transpiration credit) | replaced by deposition |
| ---: | ---: | ---: |
| 2.0 m² | 0.23 mm | **none** |
| 0.5 m² | 1.07 mm | **none** |
| 0.2 m² | 2.94 mm | **none** |
| 0.1 m² | 6.29 mm | **none** |

**and the erosion is worse than the thermal number**, because the chemical attack is additive to
it and is not counted in the table. **Delete the "0.34 to 3.7% redeposition suffices" line.** It
prices a mechanism that does not operate: the redeposition is not rate-limited, it is
thermodynamically absent.

**What could be done about it, and it is not encouraging.** Plating needs a wall below ~1,900 K.
The same surface is taking 400–2,400 MW/m², whose radiative-equilibrium temperature alone is
7,400–12,100 K, so holding it under 1,900 K means removing essentially the whole throat load
actively. **A self-healing carbon throat and a passively cooled one are mutually exclusive**, and
the honest conclusion is that the throat is a consumable. The remaining levers are transpiration
(worth about a factor of two on the flux), a refractory that is not carbon and therefore has no
deposition question at all, and **water**, which has no carbon to deposit or erode — which is
W21's point arriving from the liner side, exactly as N16 anticipated it might.

---

### W16. `held(T, ρ)` computed directly — the `ρ^-0.21` fit should be retired

**Locate:** N16's "That single curve is the most valuable thing this ask can return", and
`held(T, ρ) = held_W9(T) · ρ^-0.21`.

**Now:** a three-point fit extrapolated thirtyfold in density, load-bearing for every cold-end
number in ADR-0016.

**Should be:** the surface itself. **Methane**, fraction of full atomisation still held:

| `T` \ ρ [kg/m³] | 0.001 | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2500 | 0.162 | 0.138 | 0.124 | 0.118 | 0.114 | 0.112 | 0.109 | 0.104 | 0.093 |
| 3000 | 0.381 | 0.290 | 0.214 | 0.168 | 0.139 | 0.125 | 0.118 | 0.113 | 0.107 |
| 3500 | 0.625 | 0.521 | 0.402 | 0.307 | 0.227 | 0.176 | 0.144 | 0.128 | 0.117 |
| 4000 | 0.773 | 0.698 | 0.593 | 0.481 | 0.364 | 0.277 | 0.206 | 0.164 | 0.137 |
| 4500 | 0.922 | 0.833 | 0.737 | 0.637 | 0.511 | 0.399 | 0.296 | 0.226 | 0.173 |
| 5000 | 0.991 | 0.956 | 0.864 | 0.764 | 0.644 | 0.523 | 0.397 | 0.303 | 0.225 |
| 6000 | 0.999 | 0.998 | 0.991 | 0.959 | 0.862 | 0.742 | 0.599 | 0.470 | 0.350 |
| 8000 | 1.000 | 1.000 | 0.999 | 0.998 | 0.993 | 0.976 | 0.904 | 0.773 | 0.608 |

**Water**, same axes:

| `T` \ ρ [kg/m³] | 0.001 | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 | 1 | 3 | 10 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2500 | 0.103 | 0.065 | 0.040 | 0.027 | 0.017 | 0.011 | 0.007 | 0.005 | 0.003 |
| 3000 | 0.469 | 0.312 | 0.193 | 0.124 | 0.077 | 0.050 | 0.032 | 0.021 | 0.014 |
| 3500 | 0.866 | 0.721 | 0.521 | 0.358 | 0.227 | 0.148 | 0.093 | 0.061 | 0.039 |
| 4000 | 0.975 | 0.931 | 0.821 | 0.661 | 0.467 | 0.319 | 0.204 | 0.134 | 0.086 |
| 4500 | 0.994 | 0.982 | 0.945 | 0.864 | 0.706 | 0.531 | 0.359 | 0.243 | 0.156 |
| 5000 | 0.998 | 0.994 | 0.981 | 0.948 | 0.860 | 0.719 | 0.530 | 0.375 | 0.246 |
| 6000 | 1.000 | 0.999 | 0.997 | 0.990 | 0.968 | 0.915 | 0.796 | 0.636 | 0.455 |
| 8000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.996 | 0.988 | 0.963 | 0.903 | 0.776 |

**The exponent is not a constant and cannot be.** `d ln(held)/d ln(ρ)` runs from 0 to −0.25 for
methane and 0 to −0.42 for water. It **vanishes at both ends** — a fully atomised gas has nothing
left to shift and a fully recombined one has nothing left to break — and reaches its magnitude
where `held` ≈ 0.4. That is mass action, not a fitting artefact, and no single power law can
reproduce it.

**Where the paper's fit is wrong, and in which direction.** Against the solved surface at 1 kg/m³
(the fit's own anchor) it is exact by construction. Away from it:

| `T` | ρ | solved `held` | `held_W9 · ρ^-0.21` | fit / solved |
| ---: | ---: | ---: | ---: | ---: |
| 3000 K | 0.001 | 0.381 | 0.502 | **1.32** |
| 3000 K | 0.01 | 0.214 | 0.309 | **1.44** |
| 3000 K | 0.1 | 0.139 | 0.191 | **1.37** |
| 3500 K | 0.01 | 0.402 | 0.378 | 0.94 |
| 4000 K | 0.01 | 0.593 | 0.542 | 0.91 |
| 5000 K | 0.01 | 0.864 | 1.000 (clipped) | 1.16 |

Between 3,500 and 5,000 K the fit is good to ±10%. **At 3,000 K — which is exactly where the deep
expansions this ask wants actually exit — it overstates `held` by 32–44%**, and overstating
`held` *understates* recovery. So the paper's corrected cold-end numbers (707 s at the flown
case) are pessimistic there rather than optimistic, which is the opposite of the direction the
correction was made to guard against.

**Use the tables.** They are 72 rows each and they are committed as
`data/results/walled_nozzle/held_surface.csv` with the local exponent alongside.

---

### W17. Cooler is *not* better — the temperature dial runs the other way, and it is not close

**Locate:** N16's "Cooler is better and hotter is worse ... 8,000 K returns about 817 s against
10,000 K's 722", and ADR-0016's "Going hotter was considered and declined".

**Now:** the paper believes exit temperature is a fixed fraction of chamber temperature, so a
hotter chamber leaves more residual heat in the exhaust and returns less.

**Should be:** the first clause is right and the conclusion does not follow. Effective specific
impulse, 100 m³ methane, unnormalised (see the conventions note below):

| throat | 6,000 K | 8,000 K | 10,000 K | 12,000 K |
| ---: | ---: | ---: | ---: | ---: |
| 7.00 m² | 509 | 637 | 738 | 808 |
| 2.00 m² | 625 | 780 | 899 | 976 |
| 0.50 m² | 702 | 884 | 1,012 | 1,093 |
| 0.20 m² | 716 | 903 | 1,040 | 1,124 |
| 0.10 m² | 725 | 909 | 1,045 | 1,129 |

**Hotter wins at every throat, and by 56–59%, not by the 2.5%-per-2,000 K ADR-0016 books.**
(On the gross convention the same rows read 677 → 1,231 and 892 → 1,552, so the reversal is not a
convention artefact — the debit narrows the margin from 1.74× to 1.56× and leaves the ranking.)

**Why the paper's reasoning misses it.** Exit temperature *is* a near-fixed fraction of chamber
temperature (W20 confirms it), so the residual-heat argument is sound as far as it goes. What it
omits is that **the store charged collapses as the chamber cools**:

| chamber | 6,000 K | 8,000 K | 10,000 K | 12,000 K |
| --- | ---: | ---: | ---: | ---: |
| store charged (100 m³) | 0.336 | 0.657 | 0.887 | 0.961 |
| slug ratio `k` | 45.71 | 27.01 | 20.48 | 18.11 |
| chamber pressure | 1,008 bar | 1,127 bar | 1,284 bar | 1,448 bar |

A 6,000 K chamber does not atomise the methane — it parks two thirds of the store in C₃ and
acetylene, which is *already the recombined state*. So a kilogram holds much less, more kilograms
are needed to absorb the same 70.3 GJ, `k` rises by 2.5×, and the exhaust speed `w/√(1+k)` falls
with it. **The cold chamber is cheap on the wall and expensive on the ledger, and the ledger term
is much the larger.**

**What this does not overturn.** Going hotter still costs wall load, and W14 says how much: the
liner flux rises from 93–174 MW/m² at 10 kK to 198–351 at 15 kK. **The temperature ceiling is a
wall argument and should be stated as one** — not, as ADR-0016 currently has it, as a case where
the gain is small enough not to bother. The gain is not small.

---

### W18. Chamber volume is not a free variable — it buys slug ratio and sells freeze margin

**Locate:** N16's "Smaller is better and costs nothing. Vessel mass runs as `nRT·ρ/σ` with `nRT`
fixed by the pulse, so 100 m³ at 1,284 bar weighs what 200 m³ at 642 does."

**Now:** volume is treated as free, decided on blowdown alone.

**Should be:** the vessel-mass argument is correct and the *performance* is not flat. Effective
Isp, methane at a 0.5 m² throat:

| chamber | 6,000 K | 8,000 K | 10,000 K | 12,000 K |
| ---: | ---: | ---: | ---: | ---: |
| 50 m³ | 692 | 856 | 995 | **1,091** |
| 100 m³ | 702 | 884 | **1,012** | **1,093** |
| 200 m³ | **707** | **893** | 1,008 | 1,081 |
| 400 m³ | 706 | 879 | 971 | 1,039 |

**The optimum is interior and it moves with temperature** — 200 m³ at the cold end, 100 m³ at the
hot end, and the whole row is flat to within 1% at 12,000 K. Two effects cross:

- **A bigger chamber is thinner, so it dissociates more** (store charged 0.548 → 0.864 from 50 to
  400 m³ at 8,000 K), so `k` falls (30.27 → 22.34) and the exhaust speed rises. **This is the
  same physics as W17** and it is why "smaller is better" fails.
- **A thinner chamber freezes sooner**, because a three-body rate goes as density squared. The
  freeze margin at 8,000 K and a 0.5 m² throat is +0.64 decades at 50 m³ and −0.77 at 400 m³.

At 10,000 K and above the second effect catches the first and the curve turns over. **100 m³ at
12,000 K is the best cell in the grid that also empties inside 400 ms and keeps its chemistry**, at
**1,093 s effective** — so the paper's 100 m³ recommendation is right, for a reason it did not give, and
only in company with a *hotter* chamber rather than a cooler one.

---

### W19. The deep throat's impulse does not arrive — the chemistry freezes past `A/A*` ≈ 140

**Locate:** N16's "The chemistry finishes well past the flown throat. 707 s at `A/A*` 4.04 rising
to about 1,074 at 200, on the corrected `held`", and W6's own throat recommendation.

**Now:** the area ratio is treated as a free lever limited only by blowdown.

**Should be:** **there is a chemical wall and it sits at `A/A*` ≈ 140.** Freeze margin in decades
against the equilibrium threshold, methane at 8,000 K (negative = frozen):

| chamber \ throat | 7.00 | 2.00 | 1.00 | 0.50 | 0.20 | 0.10 | 0.05 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 m³ | 3.48 | 2.10 | 1.37 | 0.64 | **−0.32** | −1.05 | −1.77 |
| 100 m³ | 3.00 | 1.65 | 0.93 | 0.20 | **−0.74** | −1.46 | −2.17 |
| 200 m³ | 2.48 | 1.16 | 0.45 | **−0.27** | −1.19 | −1.90 | −2.60 |
| 400 m³ | 0.76 | 0.49 | **−0.05** | −0.77 | −1.68 | −2.37 | −3.07 |

**This does not contradict W5** — it extends it. W5 asked whether the flow freezes inside
ADR-0016's *bore at the flown throat* and the answer was no, with 0.4–2.2 decades of margin. This
grid asks the same question at 35–140× the area ratio, and there the answer is yes.

**What it costs.** The store actually returned peaks and then *falls*, methane at 100 m³ / 8,000 K:

| throat | `A/A*` | exit `T` | equilibrium conversion | **conversion, freeze-capped** | store returned |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2.00 m² | 14.1 | 4,234 K | 0.504 | 0.504 | 0.380 |
| 1.00 m² | 28.3 | 3,897 K | 0.558 | 0.558 | 0.428 |
| **0.50 m²** | **56.5** | 3,615 K | 0.606 | 0.606 | **0.472** ← peak |
| 0.20 m² | 141.4 | 3,306 K | 0.661 | **0.626** | 0.471 |
| 0.10 m² | 282.7 | 3,109 K | 0.697 | **0.633** | 0.461 |
| 0.05 m² | 565.5 | 2,935 K | 0.730 | **0.637** | 0.451 |

The equilibrium column keeps rising because the equilibrium branch keeps releasing chemistry the
real flow no longer has time to release. **The capped column is the honest one and it is flat
past `A/A*` ≈ 140.**

**So the recommendation W6 made, and N16 inherited, should be trimmed.** The useful range is
`A/A*` = 30–140 — a **0.2 to 1.0 m² throat**, not 0.05. Below 0.2 m² the impulse gained is
under 1%, the blowdown eats the whole 400 ms pulse period (371 ms at 0.1 m², 742 ms at 0.05 m² —
over it) and the throat fluence roughly doubles for every halving of throat area (W14). **Three independent constraints agree on the same
place**, which is the strongest thing in this document.

**The capped column is a lower bound**, deliberately: it charges the flow for the chemical energy
the equilibrium branch releases downstream of the freeze station, and a genuinely frozen gas is
also colder, so part of that energy would not have reached the jet anyway. The truth is between
the two columns. Re-solving a frozen methane expansion needs machinery `eos_methane` does not yet
have (see [Deferred](#deferred-with-cost)).

---

### W20. Exit temperature never steepens on the equilibrium branch — a 2,000 K exit is not reachable, and it is also the wrong target

**Locate:** N16's "The single most valuable curve in the whole grid: where does exit temperature
stop falling slowly with area ratio?", and its "whether a 2,000 K exit needs an area ratio of 600
or of 6,200".

**Now:** `T_e/T_c = 0.667 (A/A*)^-0.138`, fitted to four points and extrapolated a hundredfold,
with the expectation that the curve steepens back toward the ideal −0.364 once recombination
finishes.

**Should be:** **it does not steepen — it keeps flattening.** Methane, 100 m³, 10,000 K chamber:

| `A/A*` | `T_e` | `T_e/T_c` | exit ρ | local slope |
| ---: | ---: | ---: | ---: | ---: |
| 3.9 | 5,737 K | 0.574 | 4.4e-1 | **−0.189** |
| 10.1 | 4,892 K | 0.489 | 1.5e-1 | −0.152 |
| 30.6 | 4,204 K | 0.420 | 4.7e-2 | −0.124 |
| 100.9 | 3,671 K | 0.367 | 1.3e-2 | −0.104 |
| 297.5 | 3,305 K | 0.331 | 4.3e-3 | −0.091 |
| 564.5 | 3,124 K | 0.312 | 2.2e-3 | −0.085 |
| 943.7 | 2,994 K | 0.299 | 1.3e-3 | **−0.080** |

The paper's −0.138 is a fair *average* over `A/A*` = 4–140 and a bad local slope at either end.
Extrapolating the far-field slope, **a 2,000 K exit needs `A/A*` ≈ 1.4e5** — a 0.8 cm throat
radius. Not 600, not 6,200: 144,000, and unreachable.

**But the equilibrium branch is not valid out there**, because W19 says the chemistry froze at
`A/A*` ≈ 140. A genuinely frozen gas *does* steepen toward the ideal, and on the paper's own
−0.364 a 2,000 K exit would then arrive near `A/A*` ≈ 700 — close to the paper's optimistic
guess. **So the answer to "where does it steepen" is: at the freeze station, and the equilibrium
curve is structurally incapable of showing it.**

**And a 2,000 K exit reached that way is worthless.** Cold-because-frozen is not the same as
cold-because-finished: the gas is cold precisely because it stopped handing energy back.
**Exit temperature was only ever a proxy for "the chemistry finished", and once the flow freezes
the proxy inverts.** The variable to optimise is the store returned (W19), which peaks at
`A/A*` ≈ 56 while the exhaust is still at **3,615 K**. The paper should retire the 2,000 K
target rather than chase it.

Water behaves the same way with slightly shallower slopes (−0.177 to −0.071 over the same range).
Both curves are committed as `data/results/walled_nozzle/exit_temperature.csv`.

---

### W21. Water loses to methane in all 128 grid cells, including at 400 m³ and at the cold end

**Locate:** N16's "Paper-side work puts water within 5 to 9% of methane at the cold end and
*ahead* of it at 400 m³", and "The fluid choice may come down to the liner rather than the
impulse".

**Now:** water is carried as roughly competitive and better at large volume.

**Should be:** water's effective impulse is **52–83% of methane's** at every cell of the grid.
The ratio, water ÷ methane:

| `T_c` \ chamber | 50 m³ | 100 m³ | 200 m³ | 400 m³ |
| ---: | ---: | ---: | ---: | ---: |
| 6,000 K | 0.772 | 0.780 | 0.790 | **0.800** |
| 8,000 K | 0.798 | 0.800 | 0.797 | 0.763 |
| 10,000 K | 0.792 | 0.787 | 0.786 | 0.570 |
| 12,000 K | 0.786 | 0.787 | 0.606 | 0.615 |

(at a 2 m² throat; the full 128-cell table is in `surface.csv`. The low cells at 400 m³ are where
water freezes and methane has not yet.)

**A ratio, deliberately.** Both fluids are run through identical machinery with identical
conventions, so every normalisation the paper applies — vessel mass, `eta_geom`, the drift term —
cancels. **The ratio survives the conventions; an absolute number would not.** The best it ever
gets is 0.800, at 6,000 K and 400 m³ — which is also the worst corner of the grid in absolute
terms (706 s effective for methane, 565 for water). **The debit narrows the gap** by ~6 points, because
water's larger `k` means a smaller per-kilogram debit (85–208 s against methane's 167–422 s), but
it never closes it.

**Why.** Water's mean atomised particle mass is 6.0 amu against methane's 3.2, so a kilogram of
water holds fewer particles and less energy at the same temperature. It needs `k` = 40–50 where
methane needs 20–27, and both the exhaust speed `w/√(1+k)` and the free-impactor bonus `(1+k)/k`
punish that twice. This is W8's finding, unchanged, now confirmed across the whole grid rather
than at one point.

**Water's real argument is the one N16 states last, and it stands.** It has no condensed-carbon
problem, so it is exactly the fluid whose cold end this repository *can* model. If the liner or
the throat turns out to decide the design (W14, W15), water's 26–30% impulse penalty may still be
the right trade. **But it should be argued as a liner choice, not sold as a near-tie on impulse.**

---

### W22. The configuration to price first, priced — plain methane still beats it, and hydrogen is a bigger dial than 10%

**Locate:** N16's "100 m³ at 1,284 bar, a 0.1 m² throat, and a bulk slug of 90% water with 10%
liquid hydrogen ... If only one point can be run, run that one and plain methane beside it."

**Should be:**

| fluid | `k` | `p_c` | exit `T` | conv raw / capped | freeze at `A/A*` | **Isp effective** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| **90% water + 10% H₂*** | 34.42 | 1,289 bar | 3,108 K | 0.763 / 0.610 | 11.2 | **811 – 933 s** |
| water | 50.26 | 1,163 bar | 3,228 K | 0.703 / 0.596 | 21.3 | 689 – 761 s |
| **plain methane** | 27.01 | 1,127 bar | 3,109 K | 0.697 / 0.633 | 47.6 | **909 – 968 s** |

(100 m³, 8,000 K, 0.1 m² throat. \* **ESTIMATE** — see below. All three freeze at this throat,
so each is a bracket: capped conversion at the low end, equilibrium at the high end.)

**Three readings.** The hydrogen dilution is worth **+18 to +23%** over plain water, which is the
ask's own point and it is real. **Plain methane still beats the mix by 4–12%**, at both ends of
the bracket. And all three freeze here (W19), so none of these is an equilibrium figure — methane's
lead is partly that it freezes *latest* (`A/A*` 47.6 against 11.2), because its return channel is
`H + H + M` with both partners the same species while the water fluids need a scarce OH.

**Hydrogen is a real dial and 10% is not the end of it.** Swept at the best methane cell
(12,000 K, 100 m³, 0.5 m² throat), effective Isp against hydrogen mass fraction:

| H₂ by mass | `k` | `u` | **net, upper** | **net, lower** | min `Da` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0% (pure water) | 36.81 | 74.4 MJ/kg | 846 | 816 | 2.89 |
| 5% | 30.44 | 89.5 | 933 | 885 | 1.68 |
| **10%** | 26.05 | 104.0 | **1,009** | 737 † | 1.04 |
| 20% | 20.25 | 132.4 | 1,140 | 827 † | 0.445 |
| 30% | 16.55 | 160.3 | 1,251 | 900 † | 0.206 |
| *methane, same cell* | *18.11* | *147.2* | *1,093* | *1,093* | *23* |

**Hydrogen buys energy density and sells freeze margin** — `u` rises 74 → 160 MJ/kg, which is the
whole benefit and is real, but diluting the water slows its own three-body return (∝ ρ_water²) and
`Da` falls monotonically. You would need ~15% hydrogen to match methane on the upper bound and ~20%
to pass it, and by then the bracket is 313 s wide against methane's single hard number.

**† Do not quote those three lower bounds.** They are a **step-function artefact**, not physics.
Between 5% and 10% the freeze station jumps to `A/A*` = 1.00 — the throat itself, at 9,509 K with
**95% of the store still held** — because `Da` at the throat crosses the threshold of 10, and the
cap treats first-crossing-of-10 as *total* quench. The real freeze band is 10 down to 0.1, so
`Da` ≈ 1 is mid-band, not frozen. The raw column is smooth and monotone throughout
(0.681 → 0.706 → 0.728 → 0.760 → 0.783), which is what identifies the cliff as bookkeeping. Scored
on the *optimistic* edge of the OH bracket instead — which is the fair comparison, since methane's
single channel carries no scarcity bracket at all — the 10% mixture returns **983 s**, still under
methane.

**The mixture is an estimate and must be labelled one wherever it is quoted.** There is no joint
H/O equilibrium EOS here, so it is composed of `eos_water` and a hydrogen EOS as two independent
subsystems at a common temperature. That is exact for the thermal terms and each component's
internal chemistry and wrong about one thing: the free hydrogen does not shift water's own
dissociation equilibrium. Le Chatelier puts the real mixture's held store **above** this, so
**every mixture number here is an optimistic edge**. Since the mixture already loses to methane,
the conclusion is safe in the direction it matters.

**At a throat that does not freeze** (0.5 m², `A/A*` = 56.5, the W19 recommendation), methane at
100 m³ and 8,000 K returns **884 s effective** and at 12,000 K **1,093 s effective**, both as hard numbers
rather than brackets.

---

### W23. Shrinking the chamber moves the strike to the contraction, and the arrival term gets worse, not better

**Locate:** nothing in the ask — this is the successor to W10, and it is the load case for the
configuration W18 and W19 jointly recommend, so leaving it open would mean recommending a chamber
whose wall load nobody had computed.

**Now:** the section treats the side wall as the only wall.

**Should be:** below the gate the front leaves the *column* but not the *hardware*. With a
0.05–1 m² throat against a 3 m bore, the convergent section is 96–99% of the bore area, so the
front stagnates on it at **normal incidence** rather than grazing past. Two loads land there and
they scale in opposite directions:

| chamber | front leaves at `r` | speed there | share of pulse still in the front | **arrival** (St = 0.001–0.01) | **bulk** (St = 0.001–0.01) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 50 m³ | 1.19 m | 31.4 km/s | 41.8% | **6.67 – 66.7** | 2.53 – 25.3 |
| 100 m³ | 2.02 m | 16.1 km/s | 21.4% | **1.18 – 11.8** | 2.53 – 25.3 |
| 150 m³ | 2.67 m | 9.9 km/s | 13.2% | **0.41 – 4.14** | 2.53 – 25.3 |
| 178 m³ (the gate) | 3.00 m | 8.0 km/s | 10.6% | 0.26 – 2.65 | 2.53 – 25.3 |

All MJ/m², against the same 1.6 MJ/m² the section already clears the wall against.

**The bulk column is identical at every volume**, and that is the structural point: the same
70.3 GJ has to turn the corner and leave through the same throat however long the column in front
of it was. Chamber volume is simply not a lever on it. At 2.53–25.3 MJ/m² it is 1.6× to 16× the
side wall's equilibrium reference, and it sits on the same surface that W14 already has ablating
millimetres per pulse and W15 has being chemically eroded.

**The arrival column runs the wrong way.** A shorter column decelerates the front less, so it
arrives faster and on a smaller patch: at 50 m³ the front is still doing 31.4 km/s and still holds
42% of the pulse energy when it hits, concentrated on 4.4 m² of the 28.3 m² bore. That is
6.7–67 MJ/m² — **4× to 42× the side-wall reference**, and worse than any side-wall number in W12
except the 673 m³ corner.

**So "make the chamber small" has a floor, and it is around 150 m³.** By then the front has slowed
to 9.9 km/s and the arrival is down to 0.41–4.14 MJ/m², a factor of 16 below the 50 m³ case, while
the side wall is still untouched (the gate is 178 m³). Below that the arrival climbs faster than
anything else improves.

**And 150 m³ costs almost nothing on the ledger**, which is what makes it a real recommendation
rather than a compromise. Methane at 12,000 K through a 0.5 m² throat:

| chamber | `p_c` | exit `T` | conversion | freeze margin | **Isp effective** | blowdown |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100 m³ | 1,448 bar | 4,061 K | 0.674 | **+0.36** | **1,093 s** | 50 ms |
| **150 m³** | 965 bar | 3,998 K | 0.668 | +0.07 | 1,091 s | 75 ms |
| 200 m³ | 723 bar | 3,951 K | 0.657 | −0.15 (freezing) | 1,081 s | 99 ms |

**The impulse is flat across all three** — 12 s out of 1,093, well inside every uncertainty in
this document, and flatter on the corrected convention than on the credit-only one. So the choice is not an impulse choice at all:

- **150 m³** is 2.9× softer on the contraction's arrival load, and 0.29 decades thinner on freeze
  margin. Its exit at 3,998 K is 2 K below the condensed-carbon floor, so it is quoted as a bound.
- **100 m³** has a comfortable freeze margin and an exit at 4,061 K, **above** the carbon floor —
  the only cell in the whole grid that is both near-optimal and free of ADR-0050's weakness 4.

**The paper should carry both and pick on the wall**, because that is the axis where they actually
differ. If the Stanton bracket tightens toward its low end, 100 m³ is fine and the chemistry
argument wins. If it tightens toward the high end, 150 m³ is the only one of the three whose
contraction load is under 5 MJ/m².

**What this does not price.** The contraction is a curved, converging surface and the front meets
it obliquely near the bore and normally near the axis; treating it as a flat end wall overstates
the normal component near the rim. And the bulk term is an energy-budget ceiling times a Stanton
number, not a boundary-layer solve — the same weakness W12 carries, on the same footing.

---

### W24. Narrowing the throat is worth +22 to +35%, and it is priced in throat life at a factor of three per halving

**Locate:** W6 recommended a narrower throat and W19 found a wall past `A/A*` ≈ 140, but neither
priced the cost. This item crosses the two, measured from the ask's **own 7 m² throat** at its own
200 m³ chamber and 10,000 K.

**Now:** the ask treats the throat as free to choose, and the paper's only stated cost of
narrowing is the pulse stretching from 8 ms to 28 ms.

**Should be:** the pulse stretch is real but it is not the binding cost. **The throat is a
consumable** (W15: it plates nothing back at any survivable wall temperature), and narrowing it
shortens its life as `A*^1.6` — a factor of three for every halving of throat area.

#### What it buys

| `A*` | `A/A*` | exit `T` | conv raw / capped | **Isp effective** | vs 7 m² | per halving |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **7.00 m²** | 4.0 | 5,583 K | 0.406 / 0.406 | **738 s** | — | — |
| 4.00 m² | 7.1 | 5,067 K | 0.468 / 0.468 | 822 s | +11.3% | +11.3% |
| **2.00 m²** | 14.1 | 4,561 K | 0.532 / 0.532 | **902 s** | **+22.1%** | +9.7% |
| **1.00 m²** | 28.3 | 4,159 K | 0.585 / 0.585 | **965 s** | **+30.7%** | +7.0% |
| 0.50 m² | 56.5 | 3,831 K | 0.631 / 0.623 | 1,008 s | +36.5% | +4.4% |
| 0.20 m² | 141.4 | 3,481 K | 0.683 / **0.629** | 1,015 s | +37.5% | +0.7% |
| 0.10 m² | 282.7 | 3,262 K | 0.717 / **0.633** | 1,020 s | +38.1% | +0.4% |

Hydrogen runs the same curve a little steeper — 1,143 s at 7 m², 1,431 s at 2 m², **1,544 s at
1 m²** (+35.1%), then flat at about 1,630 s.

**Why it flattens is in the two conversion columns.** Raw conversion climbs all the way to 0.748,
but *capped* conversion saturates at **0.633**. Past the freeze station the equilibrium branch
keeps releasing chemistry the real flow no longer has time to release, so **below about 0.5 m² the
extra nozzle delivers nothing.** That is W19's wall, located: it starts to bite between 1 and
0.5 m², not at `A/A*` = 140.

#### What it costs — and the cost is dwell, not flux

**This is the part that changes the recommendation.** Bartz gives the *throat* flux as `D*^-0.2`,
so halving the throat area raises it by only about **7%** — 610 MW/m² at 7 m² against 932 at
0.10 m², a factor 1.5 across a factor **70** in area. What doubles is the **blowdown time**,
because the same chamber empties through a smaller hole: **8 ms at 7 m², 56 ms at 1 m², 560 ms at
0.10 m²** — a factor 70, matching the area exactly. Fluence is flux × time, so the dwell carries
essentially the whole cost.

> **A narrower throat does not heat the throat harder. It heats it for twice as long.**

Combined with W15 — the throat is eroded, not plated, at every wall temperature graphite survives
— that dwell converts directly into a pulse count. Recession is **radial**, so an eroding throat
*opens*, drifting back up the very trade curve it was narrowed to climb:

| `A*` | `r*` | recession / pulse | `dA/A` per pulse | pulses to **+10%** area | pulses to **2×** area |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 7.00 m² | 1,493 mm | 0.04 mm | 0.01% | **1,947** | 16,524 |
| 4.00 m² | 1,128 mm | 0.07 mm | 0.01% | 795 | 6,749 |
| **2.00 m²** | 798 mm | 0.15 mm | 0.04% | **262** | 2,226 |
| **1.00 m²** | 564 mm | 0.32 mm | 0.11% | **87** | 734 |
| 0.50 m² | 399 mm | 0.68 mm | 0.34% | 29 | 242 |
| 0.20 m² | 252 mm | 1.87 mm | 1.48% | 7 | 56 |
| 0.10 m² | 178 mm | 4.01 mm | 4.49% | 2 | 18 |

**The exponent is derived, not fitted, which is why it is the part to trust.** Flux goes as
`A*^-0.1` and blowdown as `A*^-1`, so recession goes as `A*^-1.1`; a fixed *fractional* area
growth needs a recession proportional to `r* ~ A*^0.5`; so **life ~ `A*^1.6`**, a factor 3.03 per
halving. The solved columns reproduce it to 2% across the whole ladder, and a test pins it.

> **The exchange rate: each halving of throat area buys 4–10% effective Isp and costs a factor of
> three in throat life** — and the rate worsens as you go. 7 → 2 m² pays 7× the life for +22%;
> 1 → 0.5 m² pays 3× for +4.4%.

#### The three limits, in the order they arrive

1. **Diminishing returns.** Marginal gain falls under 5% per halving at 1 → 0.5 m².
2. **The carbon floor — methane only.** Exit temperature drops under 4,000 K at 0.5 m², where
   `eos_methane` has no condensed-carbon phase to offer. Below that the EOS quotes a gas-phase
   equilibrium reality would partly replace with soot, so **methane's sub-0.5 m² numbers are
   model-limited independently of everything else.** Hydrogen has no such floor — one more point
   in its favour, and it is why hydrogen's curve can be trusted a step deeper.
3. **Blowdown against the 400 ms pulse period.** Methane **fails at 0.10 m²** (560 ms against
   the 400 ms budget); hydrogen still fits there (311 ms) and crosses 400 ms just below, near
   0.08 m², since blowdown goes as `1/A*`. Either way this limit arrives well after the returns
   have gone flat, so it never binds in practice — it is the backstop, not the constraint.

#### The recommendation

**2 m² if the throat is hard to service; 1 m² if it is a scheduled consumable.**

- **2 m²** — +22% (methane) / +25% (hydrogen), **262 pulses** to +10% area, equilibrium branch
  everywhere, above the carbon floor, 28 ms blowdown. This is W6's existing recommendation, and
  the 7 → 2 m² step is exactly the 8 ms → 28 ms pulse stretch the paper already quotes: **the
  paper's stated cost of narrowing corresponds to this throat and no deeper.**
- **1 m²** — +31% / +35%, but **87 pulses**. Three times the erosion for a third more impulse.
- **Below 0.5 m²** — under 1% per halving, bought at double the recession each time. Don't.

**If the throat cannot be serviced at all the calculus inverts.** Total impulse delivered per
throat goes as pulses × Isp, and that falls **50×** from 7 m² to 1 m² (1,947 × 738 against
87 × 965). So this is not a physics optimum — it is an economics one, and the tables above are
the exchange rate to price it with. **The paper should state a throat-replacement interval
alongside its throat area**; at present it states neither.

**What this does not price.** Bartz is far outside its fit — it was calibrated on chemical rockets
at tens of bar and a few thousand kelvin, and this is a thousand bar at ten thousand kelvin with a
partly ionised boundary layer — so read the exponent and not the millimetres. `ablation_depth`
takes **no transpiration credit**, and a blowing boundary layer typically halves the net flux at
these rates, so **every pulse count above is the pessimistic edge of a factor-of-two band**;
the optimistic edge doubles all of them. The freeze cap is the step function known weakness 1
records, but here it is benign: the freeze lands at `A/A*` ≈ 56–67, deep in the nozzle where
little store remains, so the cap is nearly tight. It would not be quotable if the freeze sat at
the throat. Nothing here prices the *structural* consequence of an opening throat — a drifting
`A*` moves the chamber pressure and the blowdown with it, which is a coupled problem this study
does not solve.

---

## What the ask got wrong

- **"A first estimate puts it near 200,000 K" (N9).** Right for the entering front, wrong by 6×
  to 37× for the arriving one. The layer that touches the wall has been decelerated by seven
  times its own mass of pre-charge. (W11)
- **"The 24.9° cone" (ADR-0016, carried into N9).** That is water at 45.58 km/s. Methane at
  75 km/s opens at **28.7–30.3°**. The contact station survives the correction because the two
  errors offset. (W11)
- **"Cooler is better and hotter is worse" (N16).** It reverses. The residual-heat argument is
  sound; the charged-store argument is larger and points the other way. (W17)
- **"Smaller is better and costs nothing" (N16).** Vessel mass is flat, performance is not.
  Halving the chamber raises `k`. The optimum is interior. (W18)
- **"Water within 5 to 9% of methane at the cold end and ahead of it at 400 m³" (N16).** Water is
  17–24% behind on the cells where both fluids stay in equilibrium, and 39–48% behind on the
  cells where water freezes first. It is never ahead. (W21)
- **"The chemistry finishes well past the flown throat ... about 1,074 s at `A/A*` = 200"
  (N16).** It freezes at `A/A*` ≈ 140 and the extra impulse does not arrive. (W19)
- **`held(T, ρ) = held_W9(T) · ρ^-0.21` (N16).** Good to ±10% between 3,500 and 5,000 K, wrong by
  32–44% at 3,000 K, and structurally incapable of being right at both ends. (W16)
- **"The throat is the hottest and fastest station so it should self-clean" (N9 item 5), and
  "0.34 to 3.7% redeposition suffices".** There is nothing to clean off and nothing to redeposit:
  the exhaust is *undersaturated* in carbon at every wall temperature graphite survives, because
  the carbon is in acetylene. The throat is chemically eroded, not plated. (W15)

**One thing the ask got right and should be credited for:** the 170 m³ gate estimate is correct
to 5% (W10), and the "throat heat flux stops being a footnote and becomes the item" call is
exactly right (W14).

---

## Deferred, with cost

- **N13 — the exit-plane velocity distribution.** Confirmed unreachable from here, which is the
  item's own point. Everything in this repository's walled-nozzle track is quasi-1-D; a
  divergence factor needs `⟨v_z²⟩/⟨v²⟩` over the expelled mass, which needs the 2-D axisymmetric
  Euler kernel pointed at this geometry. **Cost:** the kernel exists (`crates/`) but has never
  been run on a walled nozzle, so this is days rather than hours. **Consequence of leaving it:**
  the assumed 0.98 divergence factor stands, worth about 30 s out of 1,500 if it is wrong by a
  realistic amount — so it still does not gate anything.
- **A genuinely frozen methane expansion.** `eos_water` has `pressure_energy_frozen`;
  `eos_methane` does not. Without it, the freeze-capped conversion in W19 is a lower bound rather
  than the answer. **Cost:** moderate — the species set is large but the frozen EOS is the same
  partition functions at fixed composition. **Consequence:** the capped column understates the
  deep throats by an unknown amount that cannot exceed the gap to the equilibrium column
  (0.03–0.09 in conversion).
- **The Stanton number in W12.** The convective strike is bracketed over a factor of ten and it is
  the binding channel. **Cost:** a boundary-layer solve on the contact geometry. **Consequence:**
  the 200 m³ verdict ("comparable to the equilibrium load") is safe at either edge; the 673 m³
  verdict spans "5×" to "55× worse", which is the difference between a design margin and a
  redesign.
- **The carbon-to-acetylene kinetics themselves (W15).** The equilibrium verdict is settled and
  the boundary layer has five decades of margin against a *collision-frequency* bound, but no
  evaluated rate coefficient for the C → C₂H₂ path exists in this repository — W9 already named
  building one as the next thing. **Cost:** a literature pass under `rates.py`'s provenance rule.
  **Consequence:** the verdict would have to be wrong by five decades to flip, so this is
  confirmation rather than a gate.
- **The chemical erosion rate (W15).** That the gas erodes graphite is settled; *how fast* is
  not — it needs the same boundary-layer mass-transfer solve, now with the sign reversed.
  **Consequence:** the ablation column in W14 is a floor, not the total.
- **Softening the freeze cap from a step to a ramp** (weakness 1). **Cost:** small — it is one
  function. **Consequence:** it would tighten every bracket in W19 and W22 and remove the cliff,
  but it moves every number in this document, so it is the paper's call rather than ours.
- **A boundary-layer solve on the contraction (W23).** The end-wall loads are energy-budget
  ceilings times the same Stanton bracket as W12, and the surface is curved rather than flat.
  **Consequence:** the 150-versus-100 m³ choice is decided by where inside that bracket the truth
  sits, so it is the same deferral as W12 and closing one closes both.
- **N11 entirely.** Untouched.

---

## Conventions, so nothing here is compared with the wrong number

- **Every specific impulse in this document is `Isp effective`**, and that is the figure of merit.
  It carries **both** mass-ledger corrections:

      Isp_effective = [ (1+k) u_e − w ] / (k g₀) = w (η_jet √(1+k) − 1) / (k g₀)

  — a **credit** of `(1+k)/k` for the impactor mass the vehicle never lifted, and a **debit** of
  `w/(k g₀)` for the momentum that same impactor brings in against the ship's motion. This is the
  companion's own **head-on** form (`templateArxiv.tex`: `I/(mw) = η_jet√(1+k) − 1` head-on,
  against `+ 1` on an overtake). **The head-on burn is a momentum *debit*, not a credit.** The
  exhaust has to cancel `m_p w` before any of it becomes thrust. **Do not inherit the overtake
  study's `+1`.**
- **Both corrections scale as `1/k` and the debit is the larger**, by exactly `w/u_e`. So
  `Isp effective` sits **below** `Isp true` everywhere in this grid, and a reader who applies the
  credit but forgets the debit reads every number too high by `(w − u_e)/(k g₀)`.
- The algebra lives in one place, `chamber.IspLedger`, shared with the propellant ladder of W8 so
  the sign cannot be right in one study and wrong in the other. At `η_jet = 1` it reduces to
  `√(1+k) − 1`, which is the tamper study's independently-derived `beta_ideal`.
- **`η_jet` is exactly `√(conversion)`.** Because `(1+k) u = w²/2` makes `w²/(1+k) = 2u`, the ratio
  of squares is `u_e²/(2u)`, which is the conversion fraction. So the companion can lift `η_jet`
  straight out of the grid — and a conversion fraction quoted *as* an `η_jet` would be wrong by a
  square root.
- **The thrust floor is `1/√(1+k)`.** Below it the nozzle pushes the ship backwards. It runs
  0.146–0.233 across this grid against η_jet of 0.76–0.86, so **every cell clears it by 3.5–5.6×**
  and nothing here reverses.
- **The credit-only figure is still in the CSV** as `isp_carried_gross_s`, because it is what W5
  and W8 first reported and because `isp_effective` is built from it. It is **not** a figure of
  merit on this burn. **`isp_effective_s` in this CSV is the corrected column** — if you are
  holding a copy of `surface.csv` whose header has no `isp_carried_gross_s`, its
  `isp_effective_s` is the *old credit-only* quantity and is too high by `momentum_debit_s`.
  No column carries the launch-ledger normalisations this repository does not own — no `eta_geom`,
  no vessel mass — so **do not compare any of them directly with ADR-0016's 1,080 s or 709 s.**
- **`conversion`** is `u_e²/(2u)`, the share of the chamber's energy density that left as directed
  motion. It is convention-free and is the safest quantity to lift.
- **`held`** is the fraction of *full atomisation* still locked in broken bonds, net of the bonds
  the intermediates hold — not `1 − n_intact/n_formula`, which would charge full atomisation for
  every C₃ and C₂H₂ (see W9's note and `eos_methane.bond_energy_held`).
- **The freeze margin** is `log₁₀(Da_min / 10)`. Positive means the equilibrium branch is
  legitimate with that many decades of rate-coefficient error to spare.

---

## Provenance

| item | number | target | artifact |
| --- | --- | --- | --- |
| W10 | gate volumes 178 / 109 / 92 m³ | `make walled-nozzle-wall` | `data/results/walled_nozzle/wall_front.csv` |
| W11 | contact station, time, layer state | `make walled-nozzle-wall` | `wall_front.csv` |
| W12 | radiative and convective fluence, opacity bracket | `make walled-nozzle-wall` | `wall_strike.csv` |
| W13 | film capacity | `make walled-nozzle-wall` | printed |
| W14 | Bartz fluxes, throat fluence, ablation | `make walled-nozzle-wall` | `wall_bartz.csv` |
| W15 | saturation ratios, the no-plating band, the boundary-layer margin | `make walled-nozzle-wall` | printed |
| W16 | `held(T, ρ)` and its exponent | `make walled-nozzle-surface` | `held_surface.csv` |
| W17–W19, W21 | the 256-row grid (2 fluids x 128 cells), `isp_effective_s` column | `make walled-nozzle-surface` | `surface.csv` |
| W22 | the mixture and the hydrogen-fraction sweep | `surface.mixture()`, scored inline | printed |
| W20 | `T_e(A/A*)` and its local slope | `make walled-nozzle-surface` | `exit_temperature.csv` |
| W23 | contraction arrival and bulk loads | `make walled-nozzle-wall` | printed |
| W24 | the throat trade: Isp gain, recession per pulse, pulses to +10% and 2x area | `make walled-nozzle-wall` | `data/results/walled_nozzle/wall_throat_life.csv` |

**W24's hydrogen rows are in `wall_throat_life.csv`, not `surface.csv`.** The 256-row grid carries
methane and water only; hydrogen is run as a third fluid on the same machinery for this item, so
that CSV is the place to look for it.

Tests: `make walled-nozzle-test` (39 new tests across
`test_walled_nozzle_surface.py` and `test_walled_nozzle_wall.py`). The two that matter most are
`test_one_deep_isentrope_reproduces_a_dedicated_run_at_each_throat`, which pins the shortcut the
whole grid rests on, and `test_the_front_integrator_reproduces_front_integrate`, which pins the
N9 integrator against the established magnetic-nozzle one.

---

## Known weaknesses of this answer

1. **The freeze cap is a step function where the physics is a ramp.** `conversion_capped` charges
   the whole store still held at the *first* station where `Da` falls under 10, but real
   sudden-freezing takes about a decade in `Da` to complete. Where the crossing happens deep in
   the nozzle (all the methane rows) little is held there and the cap is nearly tight. Where it
   happens at the throat — which it does for the hydrogen-diluted mixtures above 5% — the cap
   charges 95% of the store and produces a cliff that is not physical. **Read `freeze_area_ratio`
   in the CSV: a freeze at `A/A*` = 1.00 means the cap is biting at the throat and should not be
   quoted.** Softening it to a ramp across the 10→0.1 band is deferred, because it would move
   every number in this document and is the paper's call.
2. **Water and the mixture are scored on the conservative edge of the OH scarcity bracket while
   methane carries no bracket at all**, since `H + H + M` has both partners as the same species.
   That asymmetry favours methane. Re-scoring the water fluids on their optimistic edge raises
   `Da` by about 3.3× and the mixture's effective Isp from 737 to 983 s — **still under methane**, so
   the W21 and W22 verdicts survive it, but a like-for-like comparison needs the evolved OH
   network `recombination.py` already names.
3. **The Stanton bracket (W12) is the weakest number here and it is the one that binds.** A factor
   of ten on the only channel that delivers meaningful energy to the wall.
4. **Bartz is extrapolated far outside its calibration (W14)** — from tens of bar and a few
   thousand kelvin to a thousand bar and a partly ionised boundary layer — and its equilibrium/
   frozen `c_p` edges differ by 2×. Read the size, not the digits.
5. **The freeze-capped conversion (W19) is a lower bound, not the answer** (see Deferred), and
   weakness 1 says when that bound is loose.
6. **The mixture rung (W22) omits cross-chemistry** and is an optimistic edge.
7. **`eos_methane` has no condensed carbon below ~4,000 K** (ADR-0050 weakness 4), and **most of
   this grid's interesting rows exit below it** — exit temperatures run 2,900–4,300 K. The one
   place this does *not* bite is the cell W18 recommends: 100 m³ at 12,000 K with a 0.5 m² throat
   exits at **4,061 K**, just above the floor. Every deeper throat is below it, which is a further
   reason to stop at `A/A*` ≈ 56.
8. **The graphite vapour-pressure model is Clausius–Clapeyron with a constant sublimation
   enthalpy**, anchored at the measured 1 atm sublimation point. It reproduces the JANAF value at
   2,500 and 3,000 K to within about a factor of two. At graphite's working temperature `S` is
   0.032, so it would have to be 30× wrong to flip W15's verdict; the *cold edge* of the band
   moves about ±200 K per factor of two, so quote it as 1,900 K ± 200.
9. **The nozzle is assumed conical at 15°** (ADR-0051). `Da` scales linearly in the resulting
   length, so a 10° cone would improve every freeze margin by 0.18 decades and a 20° cone would
   worsen it by 0.13. That does not move W19's conclusion but it moves the exact station.
10. **The front integrator is a snowplow**, so it assumes the swept gas moves with the front and
   the shocked layer is uniform. Both are the standard treatment and both are what
   `sec:needle_through_fog` already assumes.
