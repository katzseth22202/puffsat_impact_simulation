# Calculations owed to `aim_is_all_you_need`

Everything new from **2026-08-19 to 2026-08-21** (commits `60c7334` through `4787376`).
None of it is in either repo. All of it was derived by hand in grill sessions and lives
only in `CONTEXT.md` and in the paper prose.

> **Read the audit section first (2026-08-21).** This file was written without opening
> `puffsat_impact_simulation`. That repo was cloned and read on 2026-08-21: item 1's solve
> already exists there, item 1b is already done, item 1a was routed backwards, and its
> `f(v)` extension to 63 km/s **contradicts the premise of ADR 0015**, which was committed
> here the same day. Nothing below the routing table should be started before that section
> is read.

## The three rules this list was built under

1. **Inclusion filter: solves and sweeps only.** A calculation earns a module if it is an
   iteration, a sweep or a root-find, or if a whole paper table's cells move when one
   upstream constant changes. Printed algebra and datasheet lookups stay in the paper with
   the cites they already have. That took ~43 new numbers down to the 13 below.
2. **The chain is coded in the paper's order, with its hand-cuts intact**, so every published
   digit reproduces and the citation is honest. A *separate* fixed-point check iterates the
   loop to convergence and reports the gap as a number.
3. **Citation mechanism: a reproduction line at the end of each computed table's caption**,
   naming a `make` target. This is the convention ADR 0015 already uses internally
   ("Reproducing: `make two-leg`, plus per-cell ..."), extended into the paper. It requires
   no change to the 26 existing bare `\cite{Katz_aim_is_all_you_need_2025}` calls.

## The loop, and where the paper cuts it

```
bag volume V
   -> slug density rho = m/V
      -> Saha(rho, T) + energy conservation  ==>  plume T, ionisation fraction f
         -> pressure P = rho Rg T / M_eff(f)
            -> field B = sqrt(2 mu0 P)
               -> stored energy E_B = P V = n Rg T
                  -> field leak = (leak fraction) x E_B
                     -> waste heat  (+ blackbody intercept)
                        -> vapour fraction x = (waste - melt) / latent
                           -> saturation curve  ==>  mist T, mist P
                              -> film mass = F rho_f x Rg T / (M sigma)
   ^                                                                    |
   +---- V was chosen from a radiative-loss and opacity trade ----------+
         that depends on the plume T at the bottom of this chain
```

**Cut 1.** The plume was held at 15 000 K while the bag was sized; the real value is
14 700-26 200 K. *(Was wrong. Fixed 2026-08-21.)*
**Cut 2.** The leak fraction was held at 4.4% while the `E_B` it multiplies moved 4.43 -> 12.2 GJ.
*(Still open. This is the 3.7 kg vs 31 kg bag.)*
**Cut 3.** The vapour fraction's effect on the ignition budget was held at the paper's
81.3 MJ/kg rather than recomputed from the bag state. *(Small, probably harmless, unverified.)*

---

# Which repo: routing all thirteen

Three destinations, and the criterion is what physics the number needs.

- **`aim_is_all_you_need`** -- closed-form models, sweeps, root-finds. Inputs are known constants
  or another repo's output.
- **`puffsat_impact_simulation`** -- anything that needs the collision or the expansion actually
  simulated: densities, cooling histories, opacities, mixing.
- **Neither, yet** -- needs a capability neither repo has.

| # | calculation | repo | why |
| ---: | --- | --- | --- |
| 1 | plume state (Saha + energy) | **split** | the Saha solve is `eos_water.py` and already exists; the burn sweep and bag consequence are `aim`. **Revised 2026-08-21** |
| 1a | *the density it is solved at* | **aim** | it is `213/659.6 = 0.3229`, the bag's own sizing. The *expanded fireball* density is impact sim; the initial one is not. **Revised 2026-08-21** |
| 1b | *is the plume even in LTE?* | **already done** | `lte.py` + `data/results/lte_validity.csv`, checked directly at 45-63 km/s on 2026-08-17. **Revised 2026-08-21** |
| 2 | conductivity `sigma(T, rho, x_K)` | **impact sim** | transport coefficient; shares its Saha solve with opacity. **Rerouted 2026-08-21** |
| 3 | `tab:seed_window` | **split** | still split, but the K Saha should ride `eos_water`'s solver rather than a second one. That repo has **no alkali species today**, so it arrives with Study 1. **Revised 2026-08-21** |
| 4 | `tab:bag_sizing` | **aim** | ideal-gas pressure, closed form |
| 5 | `E_B = n Rg T` invariance | **aim** | algebra, worth a regression test |
| 6 | saturation-curve inversion | **aim** | root-find against steam tables |
| 7 | `tab:bag_state` cascade | **aim** | closed form once the leak fraction is an input |
| 8 | film mass with shape factor | **aim** | closed form |
| 9 | `tab:axial_bag` | **aim** | closed form |
| **10** | **leak fraction** | **split** | see below; the quadrature is `aim`, the cooling history is impact sim |
| 11 | snowplow field profile | **split** | 1-D snowplow is `aim`; whether a real fireball snowplows rather than jets is impact sim |
| 12 | mirror stagnation pressure | **split** | pressure balance is `aim`; the blob's arrival state is impact sim |
| 13 | two-term nozzle mass | **aim** | closed form once `E_B` is fixed |
| 14 | ice sublimation equilibrium | **aim** | root-find on a radiative balance |

---

# Audit against the real `puffsat_impact_simulation`, 2026-08-21

Everything above this line was written without opening that repo. It was cloned and read
today: HEAD `a13996a`, dated **2026-08-17**. Three of the thirteen routings change, one
checklist item turns out to be already done, and one finding contradicts ADR 0015 outright
(see the next section).

## What that repo already has

`python/puffsat/eos_water.py` -- a chemical-equilibrium water EOS, cold vapour to plasma:

- species `H2O, H, O, H+`, the **full `O+ .. O8+` Saha ladder**, and `e-`
- dissociation `H2O <=> 2H + O` by law of mass action, every ionisation stage by Saha,
  closed by **H:O = 2:1 element conservation and charge neutrality**
- real degeneracies (`g_O = 9`, `g_O+ = 4`, ...) and real potentials (`IP_H` = 13.598 eV,
  `IP_O` = 13.618 eV, then the eight-stage oxygen ladder)
- API: `composition(rho, T)`, `pressure_energy(rho, T) -> (P, e)`, `sound_speed`, `eos_grid`
- energy zero is bound molecular `H2O` at `T -> 0`, so its `e` already carries dissociation
  (917.7 kJ/mol) **and** ionisation. Item 1's "subtract ~54 MJ/kg first" step is not needed
  against it -- solve `e(rho, T) = dissipated` directly.
- documented first-pass simplification: `OH`, `H2`, `O2` omitted, which reshapes the
  2000-6000 K *transition* but preserves both endpoints.
- **no potassium.** The species set carries no alkali, so item 3's seed Saha is genuinely
  absent and should arrive as part of Study 1 rather than be built twice.

`python/puffsat/lte.py` -- the McWhirter criterion, with verdicts already written to
`data/results/lte_validity.csv`. Commit `4e89105` (2026-08-17), "Check LTE at 45-63 km/s
directly instead of inferring it."

## Item 1's three "Known gaps" are exactly what `eos_water` already closes

As written they are: ignores oxygen's second ionisation (35.1 eV), treats O's 13.618 eV as
H's 13.598, and uses a single-species Saha rather than a proper mixture. That is a list of
`eos_water.py`'s features. Coding item 1 here as specified means hand-rolling a strictly
cruder duplicate of tested code that already exists.

## The hand calculation is right, which is the argument for citing it, not redoing it

Solving `pressure_energy(rho, T)[1] = dissipated` at `rho = 213/659.6 = 0.3229 kg/m^3`,
with `dissipated = (1/2) k w^2 / (1+k)^2` at `k = 8.5`:

| closing speed | dissipated | T (item 1) | T (`eos_water`) | f (item 1) | f (`eos_water`) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 75 km/s | 265 MJ/kg | 26 200 K | 26 514 K | 0.573 | 0.580 |
| 65 | 199 | 22 400 K | 22 683 K | 0.371 | 0.380 |
| 56.53 | 150 | 19 400 K | 19 708 K | 0.217 | 0.226 |
| 45.58 | 98 | 14 700 K | 15 165 K | 0.053 | 0.062 |

1-3% in temperature, under 0.01 in ionisation fraction. **The published rows survive.** What
changes is where the number should come from. Reproduce with the impact-sim repo's `python/`
on `PYTHONPATH`; `f` there is `n_e` divided by `3 rho / m_H2O`, item 1's per-atom basis.

Worth carrying across: `eos_water`'s own docstring warns that without the oxygen ladder a
table "would overshoot the stagnation temperature severely" at high specific energy. At
26 500 K the ladder is only marginally engaged, which is why the crude version survives here.
It is not a guarantee at hotter pulses.

## Revised division of labour for item 1

| piece | repo | why |
| --- | --- | --- |
| the equilibrium solve `(rho, e) -> (T, f, P)` | **impact sim** | `eos_water.py`, already built and tested |
| the LTE validity of that solve | **impact sim, done** | `lte.py`, checked at 45-63 km/s |
| the *expanded fireball* density during the ~200 us expansion | **impact sim** | hydro; this is what the recombination-freeze check below 0.01 kg/m^3 needs |
| the *initial* bag density `rho = m_slug / V` | **aim** | a design choice, item 4's bag sizing |
| the burn-envelope sweep `(1/2) k w^2 / (1+k)^2` | **aim** | `aim` owns the closing speeds (the two-wave chain) and `k` |
| `P/P0`, `B/B0`, `E_B = n Rg T` and everything downstream | **aim** | items 4-13 |

So `src/plume_state.py` still exists, but as a **thin consumer of an impact-sim-produced
plume-state table**, not as a Saha implementation. That is the doc's own criterion applied
honestly: "`aim` -- inputs are known constants **or another repo's output**."

---

# The 63 km/s extension contradicts ADR 0015

**This is the most important thing in this file and it is not a calculation-routing matter.**

`puffsat_impact_simulation` commit `957c63d`, **2026-08-17**: "Extend `f(v)` to 63 km/s: it
stays near 0.8." Its ADR 0035 reports, on the heavy-plate scenario read off a survivability
contour:

| v [km/s] | 16 | 22 | 28 | 34 | 40 | **45** | **50** | **55** | 63 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| rho contour [kg/m^3] | 0.582 | 0.582 | 0.427 | 0.287 | 0.206 | 0.160 | 0.130 | 0.106 | 0.081 |
| binds | box | box | surv | surv | surv | surv | surv | surv | surv |
| **f** | 0.794 | 0.816 | 0.813 | 0.813 | 0.816 | **0.818** | **0.819** | **0.817** | 0.809 |

ADR 0015 (written here 2026-08-20, committed `4236e88`) rests on this premise:

> `STD_FUDGE_FACTOR = 0.8` comes from a hydrodynamic sweep across **3.2-16 km/s**. The growth
> push runs from **56.53 km/s down to 45.58 km/s**. That is roughly four times the top of the
> validated envelope and about nineteen times the specific energy. [...] **nothing rebounds
> elastically at 46 km/s.**

The growth push's 45.58-56.53 km/s is now **inside** the swept envelope, at `f ~ 0.817-0.818`.
The extension predates the ADR 0015 draft by three days.

## What keeps this from being a flat refutation

Two things, and both are load-bearing:

1. **It is a different plate.** Design section 12.1 is "Heavy-plate (16-63 km/s, 100 kg pulse,
   30 m / <= 40 t plate)". The paper's reference plate is `R = 5.0 m`, 25 kg pulse. Four times
   the pulse mass and six times the radius.
2. **It holds only on the survivability contour.** `rho(v) = min(rho_ceiling(v), rho_max)` with
   `rho_ceiling = P_limit / (c_stag v^2)`, so the cloud must be *diluted* as `1/v^2`:
   0.106 kg/m^3 at 55 km/s against 0.582 at 16, a ~5.5x dilution. Above 28 km/s it is
   survivability that binds, not the box. That dilution has an operational cost the ADR
   comparison never priced.

It is also still a single-code result; CONCLUSION.md names the FLASH cross-check as the one
open validation gate.

## The question that decides whether ADR 0015 stands

**Does `f ~ 0.818` transfer from the 100 kg / 30 m heavy plate on a diluted contour to the
paper's 25 kg / 5 m plate at the density the growth push actually delivers?**

If yes, ADR 0015's verdict should be retracted and ADR 0014's "the plate stays" reinstated,
because the matched-recovery comparison was built on the claim that 0.8 is indefensible up
there. If no -- if the paper's plate cannot fly that contour -- then ADR 0015 stands but its
Context paragraph is wrong as written and must be rewritten to say *the measured 0.8 does not
transfer to this plate*, which is a different and weaker argument than *nothing rebounds
elastically at 46 km/s*.

Do not treat ADR 0015 as settled until this is answered. Its Consequences section already
says "the decision now rests on one unmeasured number, the plate's true restitution at
45-65 km/s." **That number is no longer unmeasured.**

Related: does the PuffSat cloud density at intercept resemble the contour's 0.106-0.160
kg/m^3? Do not confuse this with the slug bag's 0.323 kg/m^3, which is a different body.

## Item 10 in detail: the leak does **not** need an MHD code

`tau_d/t_exp = mu0 sigma L^2 / (L/v) = mu0 sigma v L`, which **is** `Rm`. So

```
leak fraction  ~  1/Rm
```

and `tab:seed_window` already tabulates `Rm`. The paper has been carrying the answer two pages
earlier and quoting 4.4% as though it were unrelated. Reading `1/Rm` off that table, with the bag
consequence at `E_B` = 12.2 GJ:

| plume T | `Rm` | leak | `x` | |
| ---: | ---: | ---: | ---: | --- |
| 15 000 K | 361 | 0.28% | -0.21 | slug never finishes melting |
| 6 000 K | 400 | 0.25% | -0.22 | slug never finishes melting |
| 5 000 K | 238 | 0.42% | -0.17 | slug never finishes melting |
| 4 000 K | 76.5 | 1.31% | 0.05 | bag nearly unnecessary |
| 3 000 K | 9.2 | 10.9% | **2.48** | **exceeds the whole slug** |
| 2 000 K | 0.1 | 1000% | -- | field is gone |

**There is a cliff between 4000 K and 3000 K and the design sits on it.** Above it there is barely a
bag problem. Below it the field soaks in faster than the slug can absorb the heat, which is not a
heavier bag, it is no confinement. The published 4.4% implies `Rm ~ 22.7`, between those two rows.
Above ~5% leak, `x` passes 1 and the question stops being about bag mass at all.

So the calculation is **a quadrature, not a simulation**: weight `1/Rm(T)` by how long the plume
spends at each `T` while the field is doing work, and integrate. That splits as

- **`puffsat_impact_simulation` owes `T(t)`** -- the cooling history through the expansion. This is
  what a hydrocode with radiation transport produces natively.
- **`aim_is_all_you_need` does the integral** -- it already owns the `Rm` column.
- **No new MHD capability is required.** The `1/Rm` identity is what buys that.

**Reconciliation owed in the paper:** `tab:seed_window`'s `Rm` column and `tab:bag_state`'s 4.4%
leak line are the same physical quantity and are currently presented as independent. They should be
joined, and the cliff stated.

**Unreconciled input, revised 2026-08-21.** The factor of 28 first reported here was mostly my own
error (see item 2). Corrected, the gap is 2.3x at 15 000 K but **6x at 3000 K**, which is where the
cliff is. Settling it is Study 1 in the impact-sim repo, and it is what fixes the cliff temperature.

## What the paper calls one thing and is really six

`templateArxiv.tex` defers to "the radiation-hydrodynamic calculation" **eight times**, as if it were
a single deliverable. It is not, and the pieces need different physics:

| deferral | what it actually needs | repo |
| --- | --- | --- |
| near-Sun radiative loss, iron dopant fraction (`sec:solid_PuffSats`) | radiation transport in a dense plasma | impact sim |
| nozzle thermal limit, radiated share x sky fraction (`sec:minimum_nozzle`) | same, plus geometry | impact sim |
| fireball density, does recombination freeze below 0.01 kg/m^3 (`sec:watering_it_down`) | hydro plus recombination kinetics | impact sim |
| plume state, O's second ionisation (`sec:watering_it_down`) | equilibrium thermo | **aim** (mostly done, item 1) |
| leak fraction (`sec:watering_it_down`) | a cooling history, then a quadrature | **split**, item 10 |
| RT growth and whether axisymmetry survives (`sec:two_leg_nozzle`) | ideal-MHD stability analysis | **neither repo today** |

Only the last has no home. It is also the least load-bearing of the six, since `eq:rt_efolds` already
bounds it and minimum-B geometry is a known fix.

---

# The thirteen calculations

## `src/plume_state.py` -- `make plume-state`

### 1. Self-consistent plume state (Saha + energy conservation)
**Paper:** `sec:watering_it_down`, the paragraph beginning "Dividing 180 by 219 would suggest".
**Why coded:** it is a root-find, and it reversed a published claim within hours of my making it.

*How I did it.* Bisection on `T` over `[3e3, 2e5]` K, 90 iterations. At each `T`, Saha gives
the ionisation fraction, treating all three atoms of water as hydrogen-like at `chi` = 13.6 eV:

```
S(T) = 2 (g_i/g_0) (2 pi m_e k T / h^2)^{3/2} exp(-chi/kT) / n_a      g_i/g_0 = 1/2
f^2/(1-f) = S    ->    f = (-S + sqrt(S^2 + 4S)) / 2
n_a = 3 rho / (18.015 amu)     (atom density; rho = 0.323 kg/m^3 at the flown bag)
```

Energy available is the dissipated collision energy less ~54 MJ/kg for vaporising and
dissociating. It is spent on thermal motion plus stripped electrons:

```
dissipated per kg of blob = (1/2) k w^2 / (1+k)^2        (k = 8.5)
u(T,f) = 1.5 N_a (1+f) k T  +  f N_a chi                 N_a = 3/0.018015 x N_A per kg
```

*Results.* **The last row is the finding**: at the coldest pulse the fleet flies, the paper's
original 15 000 K assumption is the answer rather than a guess.

| closing speed | dissipated | T | f | P/P0 | B/B0 | E_B |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 75 km/s | 265 MJ/kg | 26 200 K | 0.573 | 2.75 | 1.66 | 12.2 GJ |
| 65 | 199 | 22 400 K | 0.371 | 2.05 | 1.43 | 9.1 GJ |
| 56.53 | 150 | 19 400 K | 0.217 | 1.58 | 1.26 | 7.0 GJ |
| 45.58 | 98 | 14 700 K | 0.053 | 1.03 | 1.02 | 4.6 GJ |

*Known gaps.* Ignores oxygen's second ionisation (35.1 eV), treats O's 13.618 eV as H's
13.598, and uses a single-species Saha rather than a proper mixture with the potassium seed.

**Rescoped 2026-08-21.** All three of those gaps are closed already by
`puffsat_impact_simulation`'s `python/puffsat/eos_water.py`, which nobody had looked at when
this was written. Its ladder answer agrees with the table above to 1-3% in `T`, so the rows
stand -- but the solve moves there and `aim` cites it. Only the burn-envelope sweep and the
`P/P0`, `B/B0`, `E_B` consequence stay here. See the audit section.

### 2. Conductivity `sigma(T, rho, x_K)` -- MOVED to `puffsat_impact_simulation`
**Paper:** `sec:watering_it_down`, the field-leak discussion after `tab:bag_state`.
**Status: rerouted 2026-08-21.** Full spec now in
[`impact_sim_conductivity_and_bag.md`](impact_sim_conductivity_and_bag.md), Study 1. It is a
**transport coefficient**, it shares its Saha solve with the **opacity** that repo needs anyway,
and it consumes `T(t)` and `rho(t)` which that repo owns. Summarised here so the list stays whole.

*What I did, and the error in it.* First estimate was 569 S/m at 15 000 K, electron-neutral
limited with `n_e` from **1% potassium alone**. That omitted **water's own ionisation, 5.9% at
15 000 K by Saha, which supplies 38x more electrons than the seed does.** The seed is not the
dominant electron source above ~5000 K. Corrected blend gives ~6 950 S/m at 15 000 K against the
~15 900 that `tab:seed_window`'s `Rm` = 361 implies, so the gap is **2.3x, not the 28x first
reported**.

*Where the gap still bites.* 68 S/m modelled against ~405 implied at **3000 K**, a factor of 6,
sitting exactly on the cliff of item 10 and running in the direction that makes it worse.

*Consequence already applied to the paper.* "Four hundred times more conductive" conflated the
**electron-density** ratio with the **conductivity** ratio. Conductivity saturates once Coulomb
collisions dominate; the real factor against the cool end is ~100x. Fixed in `515aca3`.

*Known gaps.* `ln Lambda` comes out **2.5**, marginal for a Spitzer formula that assumes it is
large. `Q_en` = 1e-19 m^2 is a hand-picked generic atomic value. Both sides of the comparison are
weakly sourced: the paper's `Rm` column has no published `sigma` and no stated `v` or `L`.
**Validation data exists and is already in `references.bib`**: `kerrebrock1964nonequilibrium`,
`rosa1968mhd`, `messerle1995mhd` measured potassium-seeded conductivity at 2000-3000 K, the cliff
regime. Kerrebrock also raises an unasked question -- **if electron temperature decouples from gas
temperature, the cliff may not exist at all.**

### 3. `tab:seed_window` -- potassium ionisation and `Rm` versus temperature
**Paper:** `tab:seed_window`, six rows from 2000 K to 15 000 K.
**Why coded:** Saha again, and every cell moves with the seed fraction and the density.

*How I did it.* Saha for potassium at `chi` = 4.34 eV, then
`Rm = mu0 sigma v L` with the plume's own scale length. Currently hand-tabulated.

**Revised again 2026-08-21 (audit).** `eos_water.py` carries **no alkali species at all**, so
the K Saha is not sitting there waiting to be reused. It should be *added* to that solver, as
part of Study 1, rather than written a second time here: `sigma` needs the same electron
density that the seed Saha produces, and the two must not disagree.

**Superseded in part 2026-08-21.** The `Rm` column is now known to be the **leak schedule** (item
10), and its `sigma` should come from Study 1 in the impact-sim repo rather than be asserted here.
The paper's stated window floor also moved from 2500 K (potassium condensation) to **~3300 K**,
because the leak binds ~800 K earlier. This table should be *regenerated*, not just reproduced.

---

## `src/bag_state.py` -- `make bag-state`

### 4. `tab:bag_sizing` -- field, density and radiative loss versus enclosed volume
**Paper:** `tab:bag_sizing`, six rows, now with two field columns.
**Why coded:** all 36 cells move with slug mass, plume temperature, or the opacity constant.

*How I did it.* `B = sqrt(2 mu0 rho Rg T / M_eff)`, with `M_eff` = 6 g/mol for the
dissociated-neutral column (three particles per 18 g) and `18.015/(3(1+f))` for the ionised
one. **This reproduces every published row to three figures**, which is how I confirmed the
table was ideal-gas pressure and not `(gamma-1)E` as I had first assumed.

*Results.* Cold-pulse column 21.3 / 7.9 / 4.1 / 2.0 / 0.70 / 0.35 T. Hot-pulse column is
1.66x those.

### 5. `E_B = n Rg T` invariance
**Paper:** `sec:watering_it_down`, "One thing the bag does not buy is a lighter nozzle".
**Why coded:** it is a one-line assertion that pins the nozzle mass, and it should be
regression-tested rather than trusted.

*How I did it.* `E_B = B^2 V / 2 mu0 = P V = n Rg T`. 213 kg of dissociated water is
35 500 mol, so at 15 000 K, `E_B` = **4.427 GJ at every row of `tab:bag_sizing`**, verified
numerically across a 15x range in radius.

### 6. Saturation-curve inversion
**Paper:** `tab:bag_state`, the "Temperature, from the saturation curve" row.
**Why coded:** a root-find against steam tables, currently done by interpolating by hand.

*How I did it.* Given vapour mass and bag volume, `rho_vap = x m / V`, then search the
saturation curve for the `T` where `rho_sat(T) = rho_vap`. Should use `nist_webbook_water`
or CoolProp rather than my hand interpolation.

*Results.* 0.0355 kg/m^3 -> 306 K -> 4.9 kPa (Jupiter). With the plug, 0.0371 -> 307 K.

### 7. `tab:bag_state` -- the waste-heat cascade
**Paper:** `tab:bag_state`, two columns (Jupiter 122 K storage, Earth 278 K).
**Why coded:** it is a chain of six dependent steps and it moved twice on 2026-08-21.

*How I did it.* Blackbody intercept (1.2% radiated x one tenth of the sky) plus field leak
(4.4% of 20.8 MJ/kg) gives waste heat. Subtract warming and melting (0.73 MJ/kg from 122 K:
0.26 warm ice + 0.33 melt + 0.14 warm liquid; 0.21 MJ/kg from 278 K). Divide the remainder
by 2.26 MJ/kg latent heat for `x`. Then item 6, then item 8.

### 8. Film mass with shape factor
**Paper:** `eq:bag_film_mass` and `eq:bag_shape_factor`.
**Why coded:** the coefficient folds in a derating choice (half the quoted fibre strength for
weave and seams) that is an upstream constant nobody can see in the printed number.

*How I did it.* `m_bag/m_slug = F rho_f x Rg T / (M sigma)`, numerically `2.6e-4 F x T` for
polyethylene at `sigma` = 1.75 GPa (half of Toyobo's 3.5) and `rho_f` = 970 kg/m^3.
`F = (2L+2r)/(L + 4r/3)`, running 1.5 (sphere) to 2.0 (long tube).

*Results.* 2.8 kg sphere -> 3.55 kg capsule at aspect 4 -> 3.72 kg once the plug's thermal
credit shifts the vapour fraction.

### 9. `tab:axial_bag` -- bore and conductor versus column length
**Paper:** `tab:axial_bag`, five rows.
**Why coded:** every cell moves with the standoff volume, which item 1 just changed.

*How I did it.* `r = sqrt(V/pi l)` at fixed `V` = 659.6 m^3; conductor as `B r l` normalised
to the sphere; `F` from item 8; film as `2.8 x F/1.5`.

*Results.* 5.40 m bore at 10.8 m length down to 2.05 m bore at 50 m, conductor 1.00 -> 1.76.

### 10. The leak bracket -- 3.7 kg against 31 kg
**Paper:** `sec:watering_it_down`, the paragraph after `tab:bag_state`.
**Why coded:** **this is the highest-value item on the list.** It is the difference between a
bag that is a rounding error and one that costs 15% of the slug every pulse.

*How I did it.* Hold the leak fraction at 4.4% and scale only `E_B` from 4.43 to 12.2 GJ.
Leak goes 0.92 -> 2.52 MJ/kg, waste heat 1.02 -> 2.62, `x` 0.126 -> 0.836, mist 306 -> 352 K,
film 3.7 -> **31.0 kg**. Polyethylene survives with 71 K of margin.

---

## `src/nozzle_geometry.py` -- `make nozzle-geom`

### 11. Snowplow field profile down the bore
**Paper:** `sec:needle_through_fog`, closing paragraph.
**Why coded:** it is an integration over swept mass, quoted as four spot values.

*How I did it.* Swept mass `m(x) = m_0 + lambda x` with `lambda = rho A`; dissipated energy
to position `x` is `(1/2) m_imp v_0^2 (1 - m_imp/m(x))`; pressure `0.2 E(x)/(A x)`;
`B = sqrt(2 mu0 P)`.

*Results.* ~20 T at 1 m, 12 T at 3 m, 9 T at 6 m, 5 T at the exit.

### 12. Mirror stagnation pressure versus plug position
**Paper:** `sec:two_leg_nozzle`, closing run.
**Why coded:** it is the number that decides where the plug sits on the overtake leg, and the
two answers differ by 7x in field.

*How I did it.* Ram-to-static ratio `M v^2 / ((gamma-1) E)`. Plug at the mirror end: 62.5 kg
at 22.4 km/s with 23.5 GJ in one metre of bore gives ratio 6.67 and **56 T**. Plug at the
throat end after sweeping the full column: 238 kg at 5.88 km/s through 660 m^3 gives ratio
1.17 and **7.6 T**.

### 13. Two-term nozzle mass at the revised `E_B`
**Paper:** `sec:minimum_nozzle`; the "on the order of a tonne" claim in `sec:space_mortgages`.
**Why coded:** the model exists in `todos/nozzle_rewrite_plan_2026-07-14.md` but was never
coded, and item 1 just moved its dominant input from 4.43 to 12.2 GJ.

*How I did it.* Virial floor `M >= (rho/sigma_eff) E_B` with `sigma_eff` = 0.4-1.2 MJ/kg,
plus REBCO tape as `2 pi r B l / (mu0 I_tape)` at 1.5 g/m thin-substrate tape.

**Open discrepancy:** the paper's Jupiter nozzle is 1.5 t from linear Mini-Mag scaling, while
the virial floor at 4.43 GJ gives 3.7-11 t and at 12.2 GJ gives 10-30 t. These have never been
reconciled and the paper quotes the small one.

---

## `src/cruise_thermal.py` -- `make cruise-thermal`

### 14. Ice sublimation equilibrium for the projectile
**Paper:** `sec:needle_through_fog`, the ice-versus-polyethylene paragraph.
**Why coded:** a root-find on a radiative-plus-latent energy balance. I got it wrong by six
orders of magnitude on the first attempt.

*How I did it.* Balance absorbed `340 W/m^2` (sphere-averaged at 1 AU) against
`sigma T^4 + L x flux(T)`, with `flux` from Hertz-Knudsen `P sqrt(M/2 pi R T)` and `P` from
Marti and Mauersberger, `log10(P/Pa) = -2663.5/T + 12.537`, `L` = 2.83 MJ/kg.

*Results.* Bare equilibrium **194 K**, 0.064 Pa, 7.4 kg/m^2/day, against a rod's 55 kg/m^2,
so gone in a week. Shaded at 150 K: 0.58 kg/m^2 over two years, **1.4% of a 25 kg rod**.

---

## Extension, not a new module

### `src/plume_thermal.py`
It already charges 50.4 MJ/kg to atomise the water, 60% of the 84.4 MJ/kg bill. It has no
ionisation line, and item 1 says ionisation is where most of the collision energy actually
goes: 125 MJ/kg at the hottest pulse against 86 MJ/kg of thermal motion. Add it, and make the
15 000 K assumption an argument rather than a constant.

---

# Checklist

## Modules
- [ ] `src/plume_state.py` + `make plume-state` -- items 1, 3. **Rescoped 2026-08-21:** a thin
      consumer of an impact-sim plume-state table, *not* a Saha implementation. `eos_water.py`
      already solves `(rho, e) -> (T, f, P)` better than the spec in item 1 does. What stays
      here is the burn-envelope sweep and `P/P0`, `B/B0`, `E_B`. See the audit section above.
- [ ] `src/bag_state.py` + `make bag-state` -- items 4-10
- [ ] `src/nozzle_geometry.py` + `make nozzle-geom` -- items 11, 12, 13
- [ ] `src/cruise_thermal.py` + `make cruise-thermal` -- item 14
- [ ] Extend `src/plume_thermal.py` with the ionisation term

## Owed to `puffsat_impact_simulation` instead
See [`impact_sim_conductivity_and_bag.md`](impact_sim_conductivity_and_bag.md) for both specs.
- [ ] **Study 1**: `sigma(T, rho, x_K)` across all three collision regimes, validated against the
      1960s MHD-generator literature already in `references.bib`. Reports `Rm(T)` and `1/Rm(T)`
      directly, so `tab:seed_window` is regenerated and the cliff temperature falls out as an output.
- [ ] **Study 2**: does a projectile couple to a **droplet cloud** at 0.32 kg/m^3 the way it couples
      to a vapour? Decides whether `k = 8.5` survives if the leak turns out small and the bag becomes
      an unpressurised container. Note the bag does *not* become unnecessary: spreading 213 kg over
      660 m^3 is a field requirement independent of heat.
- [ ] The cooling history `T(t)`, which item 10's quadrature consumes.
- [ ] Fireball density, for the recombination-freeze check below 0.01 kg/m^3.
- [x] ~~LTE check: item 1's Saha solve assumes it and nobody has verified it.~~ **Already done.**
      `python/puffsat/lte.py` (McWhirter) + `data/results/lte_validity.csv`; commit `4e89105`,
      2026-08-17, checked directly at 45-63 km/s rather than inferred by bracketing.

**Added 2026-08-21 by the audit:**
- [ ] **Publish the plume-state table `aim` will cite.** `(w, rho) -> (T, f, P)` across the
      burn envelope from `eos_water.pressure_energy`, so `src/plume_state.py` consumes a
      number instead of re-deriving one. `aim` supplies `w` and `rho = m_slug / V`; that repo
      supplies the solve. Cross-check already run: agrees with the hand table to 1-3% in `T`.
- [ ] **Add an alkali species to `eos_water`.** There is no potassium in the species set
      today, and item 3's seed Saha plus Study 1's `sigma` both need one. Build it once, in
      the solver that already exists, not a second time here.
- [x] **Answer whether `f ~ 0.818` transfers to the paper's plate.** **Answered 2026-08-21.**
      **Short version: the *feasibility* transfers, the *number* does not.** The paper's plate
      flies the whole growth push, but at `f ~ 0.75`, not `0.818`.

      *Method.* Ran this repo's own `puffsat/contour.py` construction with
      `(mass, plate_radius) = (25, 5)` instead of `(100, 15)`. Nothing else touched -- same
      geometry sweep, same `e_eff` table, same `c_stag`, same `P_limit`, same shape box. The
      heavy-plate column reproduces `data/results/frontier_contour_heavyplate.csv` exactly
      (16 km/s: `eta = 0.971714`, `f = 0.794298`), so this is the published harness, not a
      re-implementation. Reproduce: `uv run python audit/plate_transfer_check.py`.

      | v [km/s] | 16 | 28 | 34 | 40 | **45** | **50** | **55** | 60 | 63 |
      | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
      | heavy rho contour | 0.582 | 0.427 | 0.287 | 0.206 | 0.160 | 0.130 | 0.106 | 0.089 | 0.081 |
      | paper rho contour | 1.270 | 0.427 | 0.287 | 0.206 | 0.160 | 0.130 | 0.106 | -- | -- |
      | heavy `rf/R` | 0.300 | 0.332 | 0.380 | 0.424 | 0.461 | 0.495 | 0.529 | 0.561 | 0.580 |
      | heavy `L/D` | 0.300 | 0.300 | 0.300 | 0.300 | 0.300 | 0.300 | 0.300 | 0.300 | 0.300 |
      | **heavy `f`** | 0.794 | 0.813 | 0.813 | 0.816 | **0.818** | **0.820** | **0.817** | 0.812 | 0.809 |
      | paper `rf/R` | 0.437 | 0.628 | 0.700 | 0.700 | 0.700 | 0.700 | 0.700 | -- | -- |
      | paper `L/D` | 0.300 | 0.300 | 0.323 | 0.450 | 0.579 | 0.716 | 0.873 | -- | -- |
      | **paper `f`** | 0.798 | 0.812 | 0.804 | 0.785 | **0.767** | **0.757** | **0.745** | infeas | infeas |
      | `delta f` | +0.004 | -0.001 | -0.010 | -0.031 | **-0.052** | **-0.063** | **-0.071** | -- | -- |

      **1. The paper's plate does fly the contour across the growth push.** `rho_ceiling =
      P_limit/(c_stag v^2)` is a *facesheet* property -- 400 MPa SiC+Ti and a gas-dynamic
      coefficient -- so the contour density is **identical** for both plates. The paper's shape
      box delivers down to `rho = 0.0928 kg/m^3`, which is `rho_ceiling` at **58.9 km/s**. The
      growth push's 45.58-56.53 km/s is **inside** that. And `e_eff ~ 0.66-0.67` there: the
      rebound is not close to inelastic.

      **=> ADR 0015's stated premise -- "nothing rebounds elastically at 46 km/s" -- is false.**
      *(Read directly, 2026-08-21: `aim_is_all_you_need` is not mounted in the sandbox but is
      reachable over HTTPS git; cloned shallow to scratch. ADR 0015 is
      "matched-recovery-reverses-the-plate-verdict", 2026-08-20, status accepted.)*

      **ADR 0015 names its own decision criterion, and the measurement lands in one branch.** Its
      Consequences say: *"The decision now rests on one unmeasured number, the plate's true
      restitution at 45-65 km/s. If it is near 0.8 the plate holds; if it is near 0.4-0.5, which
      is what the shock regime suggests, the nozzle wins outright."* Measured: **0.745-0.767** on
      the paper's plate, **0.817-0.820** on the heavy plate. **That is the "plate holds" branch.
      The 0.4-0.5 branch is closed** -- the shock regime does not suggest what 0015 says it does.

      **What survives, and what has to change:**
      - **The Decision survives on its own terms.** "Compare at matched quality (`f = e1`)" is a
        *methodological* choice about not benchmarking a candidate against an incumbent's
        unmeasured number. It does not depend on `f` being low, and the nozzle still wins the
        matched diagonal at every `e >= 0.3`, including at 0.8 (36.2x).
      - **The Context paragraph must be rewritten.** Both of its factual claims are stale. `f` is
        no longer from a 3.2-16 km/s sweep -- the envelope reached 63 km/s in commit `957c63d`,
        **three days before 0015 was written** -- so "roughly four times the top of the validated
        envelope" is wrong, and "nothing rebounds elastically at 46 km/s" is wrong (`e_eff ~
        0.66-0.67` there).
      - **The crossover moves, slightly.** Log-interpolating 0015's own matched table, the nozzle
        needs `e1 ~ 0.676` to match a plate at `f = 0.80`, against `e1 ~ 0.640-0.655` at the
        measured 0.745-0.767. **A ~0.03 shift in the required `e1`** -- real, in the nozzle's
        favour, and nowhere near the `e1 ~ 0.471` that `f = 0.5` would have bought it.

      **The audit's own prediction at "Amend ADR 0015" (below) is wrong.** It says today's finding
      *strengthens* 0015. It does not: it removes the collapse that 0015's Context is built on.
      The slug-cancellation argument (`2 f m w` regardless of `k`) may be perfectly correct on its
      own terms, but it cannot "close one fewer escape route" from an objection that the
      measurement has already retired.

      **If 0015 wants a plate-side objection that survives contact with this repo, it is ablation,
      not restitution.** 0015 gestures at it in one clause -- *"A field touches nothing and cannot
      ablate"* -- and never prices it. This repo does, on the same contour
      (`frontier_contour_heavyplate.csv`, Q7/ADR-0014):

      | v [km/s] | 45 | 50 | 55 | 63 |
      | --- | ---: | ---: | ---: | ---: |
      | recession [um/pulse] | 2717-13583 | 3625-18125 | 5104-25519 | 8286-41431 |

      **Millimetres to centimetres of sacrificial layer per pulse** across the growth push. It is
      an upper bound (nothing credited to re-radiation or ADR-0014's vapor curtain) and a
      diagnostic rather than a gate -- but even the low end is mm/pulse, and over an 11-cycle
      chain it is a consumable-mass problem. **It is a measured number where `f`'s collapse was a
      supposed one.** That is the argument 0015 should be making.

      **Question back to the paper, which I cannot answer from here:** which plate does the growth
      push actually fly? On the 15 m heavy plate `f = 0.817-0.820` and `STD_FUDGE_FACTOR = 0.8`
      is defensible **as-is, unchanged**. Only on the paper's 5 m plate does it fall to ~0.75. A
      chain described as *growth* plausibly scales the vehicle, and if it does, 0015's premise
      fails even harder. If it does not, 0.75 is the number.

      **2. The entire gap is `eta_capture`, not `e_eff`.** *Above 28 km/s* the two plates fly the
      **same** contour density -- survivability binds for both -- so `e_eff` is identical by
      construction and every bit of the `delta f` is geometry. (At 16-22 km/s they diverge: the
      heavy plate's *shape box* binds at 0.582 while the paper's plate, which can build denser
      clouds, is still on the survivability ceiling at 1.270/0.686. That is why the paper plate is
      marginally *better* at 16 km/s. It does not affect the growth-push conclusion, which lives
      entirely in the identical-`rho` regime.) The smaller plate hits the footprint box edge `rf/R = 0.7` at ~34 km/s and
      from there can only keep diluting by **stretching the cloud** (`L/D` 0.32 -> 0.87). Long
      clouds splat with more radial relief, so `eta` falls 0.976 -> 0.891. The heavy plate never
      has to: it has the radius to stay at the optimal `L/D = 0.300` the whole way.

      *Ruled out as explanations (each checked, not assumed):*
      - **The `L = 10 m` anchor.** `e_eff` is interpolated on `rho` alone at `LENGTH_ANCHOR = 10.0`,
        while the contour actually delivers `L ~ 1.3-6.1 m` for both plates. I expected this to
        bias the comparison. It does not: at contour densities (0.10-0.28) the `e_eff` spread over
        `L in {6,10,14}` is only **0.002-0.006**. Design SS12.1's quoted 0.052 is at `rho = 0.01`,
        which neither plate's contour ever visits. Effect on `delta f` < 0.002.
      - **`P_limit`, `c_stag`** -- material and gas-dynamic, plate-independent by inspection.
      - **Focusing** -- does not enter the contour construction for either plate (only the
        discrete frontier), so it cannot be a differentiator.
      - **Interpolation error.** Both plates walk a *node line* of the 3x3 geometry grid (heavy
        along `L/D = 0.3`; paper along `rf/R = 0.7`), so this is 1D interpolation, not 2D
        extrapolation. Quadratic-vs-linear through the three nodes moves `eta` by <= 0.006
        (`f` by <= 0.003), and in the direction that makes the paper plate slightly *worse*.

      **3. What I do not know -- and it is the whole ballgame.** *Is `rf/R > 0.7` usable?* If the
      paper's plate could use its full radius it would want `rf/R = 0.87` at 45 km/s and **1.00 at
      55 km/s** to hold the optimal `L/D = 0.3`. A 5 m plate is almost exactly the right size to
      fly the growth push at heavy-plate efficiency -- **the box edge 0.7, not the physics, is what
      costs the 0.05.** But as `rf/R -> 1` gas spills past the plate edge and ADR-0003 counts none
      of it, so there is a real competing loss and the true optimum is probably an interior peak
      around `rf/R ~ 0.8-0.9`. Nobody has run it.

      **RESOLVED 2026-08-21 -- the run is done.** `make sweep-geometry-wide`
      (`crates/sweep --geometry-wide`, 36 cases, **13 s**). Reproduce the reduction with
      `uv run python audit/wide_footprint_check.py`.

      **First, the 0.7 edge was never a design decision.** `GEO_RFOOT_OVER_R`'s own doc comment
      reads *"design SSsweep 0.3-1.0"* and design SS7's sweep-grid table (line 167) specifies
      **`r_foot/R` | shared | 0.3-1.0 | 4 pts**. The implemented grid stopped at 0.3/0.5/0.7 and
      `contour.R_FOOT_BOX` inherited that as its upper edge. So this run **completes design SS7's
      own specified grid**; it does not extend past it. The paper should not describe 0.7 as a
      constraint -- it was an unfinished sweep.

      **`eta_capture` at d/D = 0.10 (new nodes starred):**

      | `rf/R` | 0.3 | 0.5 | 0.7 | 0.8* | 0.9* | 1.0* |
      | --- | ---: | ---: | ---: | ---: | ---: | ---: |
      | `L/D` = 0.3 | 0.9717 | **0.9783** | 0.9671 | 0.9535 | 0.9490 | 0.8719 |
      | `L/D` = 0.6 | 0.9339 | 0.9474 | 0.9115 | 0.9206 | 0.8745 | 0.8391 |
      | `L/D` = 1.0 | 0.9227 | 0.8869 | 0.8820 | 0.8511 | 0.8104 | 0.7513 |

      The predicted interior peak is real: `eta` tops out at `rf/R = 0.5` and the edge-spill
      collapse arrives at 1.0 (0.9671 -> 0.8719 at `L/D = 0.3`). Widening the box is not free.

      **Paper plate, box `<= 0.7` vs the completed box `<= 1.0`:**

      | v [km/s] | 34 | 40 | **45** | **50** | **55** | 58 | 63 |
      | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
      | clamped `f` | 0.804 | 0.785 | 0.767 | 0.757 | 0.745 | 0.738 | infeas |
      | completed `f` | 0.805 | 0.797 | **0.796** | **0.787** | **0.771** | 0.761 | 0.740 |
      | gain | +0.002 | +0.012 | **+0.029** | **+0.031** | **+0.026** | +0.023 | -- |
      | chosen `rf/R` | 0.718 | 0.801 | 0.871 | 0.900 | 0.800 | 0.800 | 0.821 |
      | still below heavy | -0.008 | -0.019 | -0.022 | -0.032 | -0.046 | -0.053 | -0.068 |

      **Answer: it recovers about half the gap, not all of it.** Over the growth push the paper's
      plate goes from **0.745-0.767 to 0.771-0.796** -- so `f = 0.8` is *reached* at 45 km/s
      (0.796) and closely approached at 50 (0.787), but the plate still trails the heavy plate by
      0.02-0.05 and the deficit widens with `v`. **`f ~ 0.78-0.80` is now the defensible number
      for the paper's plate over 45-55 km/s**, against 0.75 before the run.

      **The reach limit also moves.** Clamped at 0.7 the paper's plate fell off the survivable
      contour at 58.9 km/s; with the completed grid it flies the whole swept range including
      63 km/s. *I am not quoting the ~100 km/s the box arithmetic alone suggests* -- that answer
      leans on the one unconverged corner below **and** on extrapolating both `rho_ceiling` and
      `e_eff` past the 63 km/s top of the sweep. The honest statement is that the binding
      constraint is no longer the shape box; it is the `e_eff` table's own 63 km/s ceiling.

      **Domain convergence was measured, not assumed.** `r_max = 1.4 r_plate` measures escape room
      from the *plate* rim, so room past the *cloud* edge shrinks from `3.7 r_foot` at `rf/R = 0.3`
      to `0.4 r_foot` at 1.0 -- a too-tight domain would confine the rebound and inflate `eta`
      exactly where the new rows are. Re-ran every new node at `r_max/r_plate` = 1.4 / 2.0 / 2.8:

      - `rf/R = 0.8, L/D = 0.3`: 0.9535 / 0.9521 / 0.9530 -- spread **0.0014**, converged.
      - `rf/R = 0.9, L/D = 0.3`: 0.9490 / 0.9564 / 0.9505 -- spread 0.0074, non-monotone (noise).
      - `rf/R = 0.8, L/D = 0.6`: 0.9206 / 0.9242 / 0.9259 -- spread 0.0053, drifting *up*, so the
        reported value is if anything conservative.
      - `rf/R = 1.0, L/D = 1.0`: 0.7513 / 0.7416 / 0.7200 -- spread **0.031, NOT converged.**

      **The one unconverged corner is `rf/R = 1.0` with long clouds, and the optimizer never
      selects it** (chosen points sit at `rf/R` 0.72-0.90, `L/D` 0.30-0.71). Domain uncertainty on
      the reported `f` is **<= 0.006**, well inside the +-0.04-0.05 freeze band. If anyone later
      pushes the contour into that corner, re-run it at `r_max/r_plate >= 2.8` first.

      **What this does to ADR 0015.** Log-interpolating 0015's matched table again, the nozzle must
      reach `e1 = 0.666` against a paper plate at the recovered `f = 0.785`, versus **0.686** at
      `f = 0.80` and 0.648 at the clamped 0.757. **0015's shift away from ADR 0014 has now almost
      entirely closed** -- from a 0.038 relaxation in the required `e1` down to 0.020, against
      0015's supposed 0.471. Its verdict is not overturned (the matched-quality *method* still
      stands on its own), but the quantitative concession it won from the `f` number is now small
      enough to state as such.

      **Regression + hygiene.** `run_eta_case` was refactored to expose the domain margin; re-ran
      `--geometry-m40` and diffed **bit-identical** against the pre-change file, so the five
      frontier CSVs that consume it are untouched. New rows go to their own
      `sweep_geometry_wide.jsonl`. `cargo test -p sweep` 29/29 pass, `cargo fmt --check` clean,
      clippy's denied `all` group clean (the remaining warnings are pre-existing `pedantic`).

      **Still not done:** I did not fold the new nodes into `contour.R_FOOT_BOX` or regenerate
      `frontier_contour_heavyplate.csv`. That is a real deliverable change -- it moves published
      numbers and wants an ADR -- and it is your call, not mine. The heavy plate barely benefits
      (its optimum already sits near the `eta` peak at `rf/R ~ 0.5`); the change is worth ~0.03 in
      `f` to the *core* 5 m study, which is the one the paper cites.

      *Also unresolved, smaller:* where `R_FOOT_BOX = (0.3, 0.7)` comes from. It is cited as
      "design SS13" but I did not find its derivation. If it is a *delivery-dispersion* assumption
      rather than a plate-geometry one, it is a PuffSat design parameter -- exactly the kind of
      thing a first mover optimizes, and worth naming as such in the paper.

      *Inherited caveat:* still a single code. The FLASH cross-check `CONCLUSION.md` names as the
      one open validation gate covers all of the above.

      *Not patched:* I did not touch ADR 0015 (other repo) and did not adjust any number in this
      repo to make the two agree. Note also that `contour.contour_point()` **cannot** be run on
      another plate as shipped -- it takes no `mass`/`plate_radius` and picks them up from
      `rho_max_achievable()` / `l_over_d_for()` defaults. Worked around in the scratch script
      rather than changing the shipped API; worth a small refactor if this comparison recurs.

## The convergence check (rule 2)
- [ ] `make bag-converge` -- iterate volume, plume state, leak and bag state to a fixed point
      and **report the gap against the published tables as a number**, not as a replacement.
      Expect it to disagree; the interesting output is by how much.

## Regression tests
- [ ] `tab:bag_sizing` reproduces to three figures, both columns
- [ ] `tab:bag_state` reproduces both columns
- [ ] `tab:axial_bag` reproduces all five rows
- [ ] `tab:seed_window` reproduces all six rows
- [ ] `E_B = n Rg T` holds across a 15x radius sweep
- [ ] The 2026-08-21 burn sweep reproduces (26 200 K / 0.573 down to 14 700 K / 0.053)

## Paper edits (rule 3)
- [ ] Add "Reproduce with `make <target>`" to the caption of `tab:bag_sizing`,
      `tab:bag_state`, `tab:axial_bag`, `tab:seed_window`
- [ ] Leave the 26 existing bare cites alone

## Documentation in the companion repo
- [ ] **ADR 0016**: the plume state is solved, not assumed. Records that Saha allows only
      5.9% ionisation at 15 000 K, that the temperature must therefore run to 26 200 K at the
      hottest pulse, and that the original 15 000 K is the *coldest-pulse* answer.
- [ ] **ADR 0017**: the slug bag is an axial capsule with an ice plug. Records the
      bore-for-conductor trade and the aperture argument.
- [ ] **Amend ADR 0015**: today's finding *strengthens* it. A slugged plate delivers `2 f m w`
      exactly as an unslugged one does, because `(1+k)` in the mass cancels `1/(1+k)` in the
      drift speed. So the plate genuinely cannot buy its way out of the 46 km/s objection, and
      0015's argument holds with one fewer escape route. Add to Consequences.
- [ ] `CONTEXT.md`: note that `f = 0.8` now has a second unmeasured companion, the leak
      fraction, and that the leak decides a 3.7 kg versus 31 kg bag.

## Ranked by value if only some get done

**Revised 2026-08-21.** A new item takes the top slot, and item 1 drops off the `aim` list.

-1. **Does the measured `f(v)` transfer to the paper's plate?** ADR 0015 was committed here on
   2026-08-21 (`4236e88`) asserting the plate's restitution at 45-65 km/s is unmeasured. The
   impact sim measured it on 2026-08-17 and got ~0.818. Until the plate-transfer question is
   answered, an accepted ADR in this repo rests on a premise its companion has already
   addressed. Nothing else on either list changes a published verdict this directly.

0. **Impact-sim Study 1, conductivity.** It sets where the cliff of item 10 sits, the validation
   data already exists and is unused, and one literature check (Kerrebrock's electron-temperature
   decoupling) could dissolve the problem entirely. Cheapest high-value item on either list.
1. **Item 10, the leak bracket.** Decides the vehicle, not a digit.
2. ~~**Item 1, the plume state solve.**~~ **Moved to the impact sim 2026-08-21** -- the solve
   exists there as `eos_water.py`. What is left here is the burn sweep that reads its output,
   which is bookkeeping rather than physics and carries no rank of its own.
3. **Item 13, the nozzle mass.** There is an unreconciled 1.5 t versus 10-30 t sitting in the
   paper right now.
4. **Item 7, the bag-state cascade.** Six chained steps, moved twice in one day.
5. Everything else.

---

# New questions raised by the impact-sim work

**Opened 2026-08-21 while working the routed items.** These are *not* on anyone's list yet.
Nothing here is a blocker; they are recorded so we can decide what to attack rather than
rediscovering them later. Each says what it would take to settle and what it is worth.

**Framing (Seth, 2026-08-21):** structural optimization of the vehicle will move these numbers in
the real world once there are empirical results from an actual Jupiter round trip. That is the
point, not a caveat -- the first mover who flies and tunes gets a compounding advantage over a
later competitor, because the benefit rises steeply with efficiency. So the bar for this repo is
**"is the physics obviously wrong, and does everything the paper says have plausible backing?"**
-- not "is this number final." Where a question is honestly open, say so; where simulation can
settle it cheaply, say that too, and say how cheaply.

## Q-A. Which plate does the growth push fly? *(blocks a paper number)*
`STD_FUDGE_FACTOR = 0.8` is defensible **unchanged** on the 15 m heavy plate (`f = 0.817-0.820`)
and needs to become ~0.78-0.80 on the paper's 5 m plate. The two are a real 0.02-0.05 apart and
the gap widens with `v`. **Cost to settle: zero simulation** -- it is a statement about what the
paper means. **Worth: it decides whether an ADR 0015 rewrite is needed at all.**

## Q-B. Should `contour.R_FOOT_BOX` be widened to design SS7's specified `0.3-1.0`? *(deliverable change)*
The data now exists (`sweep_geometry_wide.jsonl`). Widening moves published numbers -- the core
5 m study gains ~0.03 in `f` over 45-55 km/s; the heavy plate barely moves. **Cost: ~1 h** (box
constant, regenerate `frontier_contour_heavyplate.csv`, an ADR, refresh the figures).
**Worth: high for the core study, near-zero for the heavy plate.** Deliberately not done
unilaterally.

## Q-C. Is the `rf/R = 1.0` / long-cloud corner ever entered? *(numerics hygiene)*
That corner is **not** domain-converged (`eta` drifts 0.031 over `r_max/r_plate` 1.4 -> 2.8) and a
too-tight domain inflates `eta`. Today's optimizer never selects it, so nothing published is
affected. **If Q-B is taken, add a guard** so a future contour cannot silently optimize into it.
**Cost: minutes** (assert, or re-run those nodes at `r_max/r_plate >= 2.8`).

## Q-D. The contour drops the focusing factor the discrete frontier applies. **CONFIRMED BUG.**
**Worked 2026-08-21.** `heavyplate_frontier` (and `analysis.py`) compute
`peak = c_stag rho v^2 * focusing`, where `focusing = peak_local_pressure_concave /
peak_local_pressure_flat` from the geometry sweep. `contour.rho_ceiling = P_limit/(c_stag v^2)`
**omits it.** Focusing runs **1.15 to 2.20** across the shape box and rises steeply with
footprint, so the omission is not small. Reproduce: `uv run python audit/focusing_check.py`.

**Cross-validation says the contour is wrong, not my correction.** Re-solving the contour with
focusing as a fixed point (focusing depends on shape, shape on `rho`, `rho` on focusing; damped,
converges in ~8 iterations) puts the heavy plate at 45 km/s at **`rho = 0.1261`, focusing 1.27**.
The discrete frontier's own best surviving shape there is **`rho = 0.1258`, focusing 1.27,
peak = 399 MPa** against the 400 MPa baseline -- essentially the same point. The **published**
contour sits at `rho = 0.1604`, which at focusing 1.27 is **~509 MPa, over the baseline it claims
to respect.** The discrete construction always had this right; the ADR-0035 contour refactor
dropped it.

**Consequence -- the headline deliverable is fine, the core study is not:**

| v [km/s] | 28 | 34 | 40 | 45 | 50 | 55 | 63 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **heavy** published | 0.813 | 0.813 | 0.816 | 0.818 | 0.820 | 0.817 | 0.809 |
| **heavy** corrected | 0.819 | 0.812 | 0.815 | 0.817 | 0.813 | 0.805 | 0.783 |
| shift | +0.006 | -0.001 | -0.001 | -0.001 | -0.007 | -0.012 | **-0.026** |
| **paper** (box 1.0) | 0.812 | 0.805 | 0.797 | 0.796 | 0.787 | 0.771 | 0.740 |
| **paper** corrected | 0.797 | 0.789 | 0.767 | 0.744 | 0.734 | 0.716 | 0.684 |
| shift | -0.014 | -0.017 | -0.030 | **-0.051** | **-0.054** | **-0.055** | **-0.057** |

**The 15 m heavy plate's published `f = 0.79-0.82` survives** (worst case -0.026 at 63 km/s, well
inside the +-0.04-0.05 freeze band). **The 5 m core plate does not:** it loses ~0.05 across the
growth push.

**And the plate should be flat, not concave.** Once local pressure concentration is priced, the
shallow dish's `eta` gain no longer pays for its focusing penalty on a *small* plate at high `v`.
Best curvature flips to `d/D = 0` at 45 km/s and stays there: at 55 km/s flat gives 0.716 against
concave's 0.679. ADR-0021 foreclosed the *deep* dish; this says shallow concave also loses its
edge on the small plate up here. On the heavy plate `d/D = 0.10` still wins throughout.

**Net for the paper's plate, all corrections applied: `f ~ 0.72-0.74` over 45-55 km/s, on a flat
plate.** Note this lands close to the pre-wide-sweep 0.745-0.767 -- **the wide-footprint gain
(+0.03) and the focusing correction (-0.05) very nearly cancel.** My interim "0.78-0.80" was
computed on the repo's contour before this bug was found; it should not be quoted.

**Correction to Q-B.** Widening the box is no longer merely an optimization worth +0.03. With
focusing applied, the paper's plate clamped at `rf/R <= 0.7` is **infeasible above ~34 km/s** --
the required density falls below what that box can deliver. Widening is what lets the core plate
fly the growth push at all.

**ADR 0015 crossover, final:** nozzle needs `e1 = 0.631` against the corrected paper plate
(`f = 0.731`), against 0.676 at `f = 0.80` and 0.683 against the corrected heavy plate. Still
nowhere near the 0.471 that 0015's "0.4-0.5" branch assumed.

**FIXED 2026-08-21**, on sign-off. `contour.focusing_at` + `contour.survivable`, and
`contour_point` reformulated as a **constrained maximum over the shape box** rather than
"pick `rho` first, then maximize `eta` on the iso-density curve". Making `rho` shape-dependent made
the old formulation circular; stating the search over the 2D box removes the circularity instead of
iterating it away, and lets the optimizer trade `eta` against `focusing`. The inner search is a
bisection on `L/D`, since `peak`, `eta` and `rho` all fall monotonically with cloud length
(pinned by `test_peak_and_score_fall_with_cloud_length`), so the contour lands **exactly** on the
pressure limit rather than one grid step inside it.

| v [km/s] | 16 | 22 | 28 | 34 | 40 | 45 | 50 | 55 | 63 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| published `f` | 0.794 | 0.816 | 0.813 | 0.813 | 0.816 | 0.818 | 0.819 | 0.817 | 0.809 |
| corrected `f` | 0.794 | 0.816 | 0.812 | 0.812 | 0.815 | 0.817 | 0.813 | 0.805 | **0.783** |
| corrected `rho` | 0.582 | 0.582 | 0.359 | 0.233 | 0.165 | 0.126 | 0.093 | 0.070 | 0.047 |
| focusing | 1.15 | 1.15 | 1.18 | 1.22 | 1.25 | 1.27 | 1.39 | 1.52 | 1.70 |

**Headline moves 0.79-0.82 -> 0.78-0.82.** Every point now flies at or under 400 MPa; before the
fix they ran **470-565 MPa** across 28-60 km/s. Also **the schedule value moved 0.017 -> 0.040**:
the pinned shape must survive the worst velocity, and focusing makes that bite harder.

*Done test-first, after a false start.* The fix was written before its test, which is backwards for
this repo. The red step was supplied afterwards by reverting `contour.py` to the committed version
and running the new regression check against it: **9/9 velocities over the limit, 470-565 MPa.**
A test that has never been seen to fail is not known to be a test.

Updated: ADR-0035 (correction block + the 0.040 schedule figure), design SS12.1, `CONCLUSION.md`,
`frontier_contour_heavyplate.csv` (now carries `focusing` and `peak_mpa` columns so the flown
pressure is auditable from the deliverable itself).

## Q-I. Is the focusing model itself right? **ANSWERED 2026-08-22: yes, to 1.4%.**
`peak = c_stag rho v^2 * (P_local_concave / P_local_flat)` treats the **flat 2D** peak as equal to
the **1D plane-wave** peak, then scales by the concave/flat ratio. Q-D's fix is a *consistency* fix
regardless, but this premise underneath it had never been checked, and it scales every
survivability verdict in the repo.

*Correction to my own scoping note:* I said the sweep "already records `peak_local_pressure` for
the confined run". It did not -- the kernel **computes** it for both runs and the sweep recorded
only the free one. Exposing it was still a one-field change (`GeoRecord.peak_local_pressure_confined`),
so "nearly free" held, but the reason was wrong.

**Seam** (agreed before writing the test): the free **flat** run against the **confined**
(plane-wave) run of the *same kernel*, so scheme error is common-mode and divides out -- the
construction ADR-0003 already uses for `eta_capture`. **Tolerance: 10%** on peak pressure, which
sits well inside the SiC+Ti margin (400 MPa is the conservative floor of a 400/700/900 band).

**Result: the flat plate runs 1.0011-1.0142x the plane-wave peak across the shape box** -- worst
departure **1.4%**, inside even the stricter 5% bar that ADR-0003's cross-kernel check uses.

Two things worth naming:
- **The sign is the opposite of the obvious guess.** I predicted a finite footprint would relieve
  sideways and push the local peak *below* plane-wave, making the model conservative. It does not:
  flat sits slightly **above** plane-wave, rising monotonically with both `L/D` and `r_foot/R`, so
  the survivability model is marginally **optimistic**. At <= 1.4% that is negligible -- but it is
  not zero, and it is not the direction anyone would assume.
- **The plane-wave denominator is genuinely shape-independent.** All nine confined runs agree to a
  part in 10^4 despite differing in slab length and domain radius. That is what licenses treating
  `focusing` as pure geometry, and it is an independent check on the confined boundary condition.
  Kept as its own test.

**So Q-D's fix rests on a validated premise.** Done red-green this time: the test was written
first and failed on `KeyError: 'peak_local_pressure_confined'` before the field existed.

**Still open (Seam B, not taken):** the *cross-kernel* leg -- whether the confined-2D peak matches
the 1D Lagrangian `hydro1d` `peak_wall_pressure` that `c_stag` is actually measured from. Seam A
cannot see a normalization error between the two kernels. ADR-0003 already cross-checks `1+e_eff`
between them to ~5%, so this is a narrowing rather than an unknown. **Cost: ~half a day** (matched
configs across two kernels, different units).

## Q-E. Does the ~5 m cloud length the contour delivers matter, given `e_eff` is read at `L = 10 m`?
Measured and **answered no** at contour densities (spread 0.002-0.006), but the *reason* is that
`tau >> 1` there. The design's own 2026-08-18 correction says `tau ~ 1` at the dilute end. If any
future scenario flies `rho < 0.04`, the `L = 10 m` anchor stops being harmless (spread 0.0155 at
0.04, 0.0517 at 0.01). **Cost: an assert. Worth: cheap insurance, no present error.**

## Q-F. Kerrebrock's decoupling -- does the conductivity cliff exist at all? *(physics, open)*
Item 2 already flags it: if electron temperature decouples from gas temperature in the seeded
plume, the ~3300 K cliff the leak schedule rests on may not be there. This is a **real** effect in
MHD generators (it is what non-equilibrium MHD *is*), and it runs in the direction that helps the
paper. **Cost: it falls out of Study 1 if the two-temperature option is built in from the start;
expensive to retrofit.** **Worth: high -- it is load-bearing for item 10's leak bracket.**

## Q-G. Both sides of the `sigma` comparison are weakly sourced. *(paper backing)*
The paper's `Rm` column has no published `sigma`, no stated `v`, and no stated `L`; the model side
has `ln Lambda = 2.5` (marginal for a Spitzer formula that assumes it is large) and a hand-picked
`Q_en = 1e-19 m^2`. **Until Study 1 lands, no factor-of-a-few claim in either direction is
defensible.** Validation data is already in `references.bib` (`kerrebrock1964nonequilibrium`,
`rosa1968mhd`, `messerle1995mhd`), measured in exactly the 2000-3000 K cliff regime.

## Q-H. Ablation is the plate's real exposure, and it is a diagnostic rather than a gate.
2.7-13.6 mm/pulse at 45 km/s rising to 8.3-41.4 mm at 63 (upper bound, nothing credited to
re-radiation or ADR-0014's vapor curtain). Over an 11-cycle chain this is a consumable-mass
problem, and it is **measured** where the `f` collapse was supposed. **Cost to tighten: crediting
re-radiation + the vapor curtain is a real modelling task, not a re-run.** **Worth: high -- this,
not restitution, is the strongest plate-side objection, and the paper currently under-argues it.**
