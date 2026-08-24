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
- [x] **Study 1**: `sigma(T, rho, x_K)`. **DONE 2026-08-22.** `python/puffsat/conductivity.py`,
      `make analysis-conductivity`, 7 tests. Built test-first with the seams agreed in advance:
      `electron_density`, `sigma`, `magnetic_reynolds` public; collision frequencies private.

      **It reproduces every number the audit computed by hand, independently.** That is the check
      that the assembly is not quietly missing a term:

      | quantity | audit (by hand) | this module |
      | --- | ---: | ---: |
      | water/seed electron ratio at 15 000 K | 38x | **37.8x** |
      | water ionised fraction at 15 000 K | 5.9% | **5.8%** |
      | blended `sigma` at 15 000 K | ~6950 S/m | **6993 S/m** |
      | `sigma` at 3000 K | 68 S/m | **67.1 S/m** |

      Plus the analytic acceptance test the repo's standards ask for: with electron-neutral
      collisions switched off, `sigma` reproduces **Spitzer's** conductivity
      (`1.899e4 T_e[eV]^1.5 / ln Lambda`) to 2%.

      **Regenerated `tab:seed_window`** (rho = 0.32, x_K = 0.01, `v L = 1.81e4`):

      | T [K] | 2000 | 3000 | 5000 | 8000 | 11000 | 15000 |
      | --- | ---: | ---: | ---: | ---: | ---: | ---: |
      | seed ionised | 0.0002 | 0.020 | 0.565 | 1.000 | 1.000 | 1.000 |
      | `sigma` [S/m] | 1.06 | 67.1 | 452 | 731 | 2184 | 6993 |
      | `Rm` | 0.024 | 1.53 | 10.3 | 16.6 | 49.7 | 159 |
      | leak `~1/Rm` | 1.000 | 0.655 | 0.097 | 0.060 | 0.020 | 0.006 |

      **The cliff is an output, not an assertion: `Rm = 1` at 2845 K** for that `v L`, against the
      paper's revised floor of ~3300 K.

- [x] **Item 3 (`tab:seed_window`)** falls out of Study 1 above -- regenerated, not reproduced.
- [x] **Add an alkali species.** Done as a *trace layer* over `eos_water` rather than inside its
      Newton solve: potassium's Saha is closed exactly against water's electron field
      (`n_e^2 + n_e(K - n_w) - K(n_w + n_K) = 0`), so a verified solver every EOS table depends on
      is left untouched. The one approximation is named in the docstring: water's own ionisation is
      not re-solved with the seed's electrons present, worth ~1.3% at 15 000 K.
- [x] **The Velikhov criterion, and which legs it actually exposes. DONE 2026-08-24.** ADR-0038,
      `python/puffsat/electrothermal.py` + `make analysis-electrothermal` ->
      `data/results/electrothermal_scan.csv`, plus `electron_energy_balance`,
      `critical_hall_parameter` and `electrothermal_loop` in `conductivity.py`. 13 tests.

      **Not on anyone's list -- it was `BETA_CRIT = 2`, which Q-F(b) recorded as "engineering
      practice taken on authority" and Q-K said was worth obtaining rather than designing around.**
      Obtained (Petit and Geffray 2009, restating Velikhov 1962 via Petit and Valensi 1969), and
      **both halves of the old screen were wrong, in opposite directions**: the ionisation gain was
      differentiating the seed alone and read ~0 above 7000 K when the true value is 2-6, and
      `beta_cr = 2` is the strongly two-temperature limit of a quantity that **diverges** as a
      plasma approaches equilibrium -- which Q-M had just established this plume does.

      | leg | `T` at exit | `T_e - T_g` | `beta_cr` | `beta` | unstable | e-folding |
      | --- | ---: | ---: | ---: | ---: | ---: | ---: |
      | 75 km/s | 16 224 K | 0.4 K | 6 535 | 0.38 | none | -- |
      | 65 | 14 151 | 1.2 | 1 851 | 0.58 | none | -- |
      | 56.53 | 11 681 | 7.9 | 203.9 | 1.20 | none | -- |
      | 45.58 | 4 597 | **813.5** | **2.07** | **4.63** | **1.70 of 2.71 ms** | **3.2 us** |

      **What is settled: the hot legs are stable under every variant tested**, including with the
      near-equilibrium term deleted entirely, so no remaining uncertainty reaches them. **Q-F closes
      as a by-product** -- the electron elevation is now computed rather than swept (+0 K hot,
      +814 K at the cold exit), and every leg exits above its own conductivity cliff.

      **What is NOT settled, and it is the cold leg's whole verdict: see Q-O.** The gap is driven
      entirely by the field-gradient current, whose layer thickness is estimated by a scaling. The
      other defensible choice -- and the low plasma beta (0.016) argues for it -- makes **every**
      station stable. Both readings ship; the conservative one is the default. **I don't know which
      is right, and settling it is a 2D resistive-MHD solve this repo does not have.**
- [ ] **Study 2**: does a projectile couple to a **droplet cloud** at 0.32 kg/m^3 the way it couples
      to a vapour? Decides whether `k = 8.5` survives if the leak turns out small and the bag becomes
      an unpressurised container. Note the bag does *not* become unnecessary: spreading 213 kg over
      660 m^3 is a field requirement independent of heat.
- [x] **The cooling history `T(t)`, which item 10's quadrature consumes. BUILT 2026-08-22.**
      `python/puffsat/expansion.py` + `make analysis-expansion` -> `data/results/cooling_history.csv`.
      Quasi-1D steady isentropic expansion on `eos_water`, run on **both** branches of
      ADR-0026. Full findings in Q-L below; the three headline consequences are that the
      residence time is **1.7-2.8 ms, not the 5.5 ms Q-K assumed**, that `v L` is now an
      **output** (5.5e4-9.7e4) rather than the guess Q-G could not source, and that the
      equilibrium/frozen bracket on exit temperature is **a factor 3**, which spans the whole
      cliff. So the answer to "is the plume still hot at the exit" is *yes on every leg in
      equilibrium, and no on the cold leg if recombination freezes* -- and **Q-M then settled
      which of those holds: equilibrium, by 2 to 5 decades. Yes on every leg.**
- [x] **Fireball density, for the recombination-freeze check below 0.01 kg/m^3. ANSWERED
      2026-08-24: it freezes, at the lip, with the store full.** `python/puffsat/fireball.py` +
      `make analysis-fireball` -> `data/results/fireball_freeze.csv`, 6 tests. **Q-P below carries
      the finding and the consequence; this is the summary.**

      Q-M found the dissociation store never returns inside the nozzle and said its marginal
      Damkohler number would bind downstream instead. It does, **immediately** downstream:

      | leg | `rho` at freeze | `T` | `f_diss` still held | stranded | of the dissipated budget |
      | --- | ---: | ---: | ---: | ---: | ---: |
      | 75 km/s | 2.35e-2 | 16 063 K | **1.0000** | 50.9 MJ/kg | **19.2%** |
      | 65 | 2.29e-2 | 13 974 | **1.0000** | 50.9 | **25.6%** |
      | 56.53 | 2.10e-2 | 11 271 | **1.0000** | 50.9 | **33.9%** |
      | 45.58 | 1.01e-2 | 3 908 | 0.9170 | 46.7 | **47.7%** |

      The three hot legs cross `Da = 1` at the **first station past the nozzle lip** (exit `rho` is
      ~2.5e-2), at 100% of the store held. The cold leg gets ~0.4 decades further out because it is
      cold enough for the rate to be fast. **The paper's own 0.01 kg/m^3 threshold is where the
      cold leg lands, to within 1%** -- which reads as confirmation of where the paper drew it, not
      as a coincidence worth leaning on.
- [x] ~~LTE check: item 1's Saha solve assumes it and nobody has verified it.~~ **Already done.**
      `python/puffsat/lte.py` (McWhirter) + `data/results/lte_validity.csv`; commit `4e89105`,
      2026-08-17, checked directly at 45-63 km/s rather than inferred by bracketing.

**Added 2026-08-21 by the audit:**
- [x] **Publish the plume-state table `aim` will cite. DONE 2026-08-24.**
      `python/puffsat/plume.py` + `make analysis-plume` -> `data/results/plume_state.csv`,
      6 tests. **171 rows, `w` = 44-76 km/s on a 2 km/s grid with the four quoted anchors
      inserted exactly, `rho` = 0.05-2.0 kg/m^3.** The CSV is **committed to the tree**, not
      just gitignored output: a consumer in another repository cannot run `make`.

      **The audit's 1-3% cross-check reproduces, and the sign of it is now explained.** The solve
      runs *warmer* at every anchor, because the audit charged 54 MJ/kg for vaporisation plus
      dissociation while `eos_water`'s bond energy is **50.9 MJ/kg** -- the ~3 MJ/kg difference
      stays in the thermal pool. The gap widens toward the cold end because there is least energy
      there for it to hide in.

      | `w` [km/s] | dissipated | `T` solved | `T` hand | `f` solved | `f` hand | `P` [MPa] |
      | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
      | 75 | 264.9 MJ/kg | 26 514 K | 26 200 | 0.5805 | 0.573 | 18.74 |
      | 65 | 199.0 | 22 684 | 22 400 | 0.3798 | 0.371 | 14.00 |
      | 56.53 | 150.5 | 19 708 | 19 400 | 0.2259 | 0.217 | 10.80 |
      | 45.58 | 97.8 | 15 165 | 14 700 | 0.0616 | 0.053 | 7.20 |

      **The table has to be two-dimensional, which the original item did not assume but `aim`
      needs.** Dissipated energy depends only on `w` and `k`, never on density -- but Saha does,
      so the same budget lands at a different temperature in a different bag. At 56.53 km/s,
      `rho` 0.05 -> 1.0 moves the plume **16 857 -> 21 795 K** (and `f` 0.254 -> 0.206): denser
      pushes recombination, spends less of the budget stripping electrons, leaves more as heat.
      Since `aim` sets `rho = m_slug / V` and the enclosed volume is a live design variable
      (item 4), a single row would not have served.

      **One trap, and it is load-bearing.** `eos_water` references `e` to bound molecular H2O at
      `T -> 0`, so the bond energy is **already inside `e`** -- the balance is `e(rho, T) =
      e_dissipated` with nothing subtracted. Porting the audit's formula, which subtracts 54 MJ/kg
      *first*, double-charges the bond: 11% low at the hot anchor (plausible-looking) and **4 672 K
      instead of 15 165 K at the cold one.** Pinned by a test for exactly that reason.

      **What did NOT change: `expansion.PLUME_STATES` still carries the hand temperatures**, so
      the cooling history, Q-M and ADR-0038 are all anchored on them rather than on this solve.
      Re-basing them would move `T(t)` by 1-3% and is a deliberate call, not a cleanup -- there is
      a test asserting the two agree to 3.5%, which is what would start failing if they drift.
- [x] **Add an alkali species to `eos_water`.** **Resolved 2026-08-22 -- differently from how it
      is written here, and deliberately.** This asks for potassium *inside* `eos_water`'s Newton
      solve. It was built as a **trace layer over** it instead (the checked item above), because
      potassium's Saha closes exactly against water's electron field
      (`n_e^2 + n_e(K - n_w) - K(n_w + n_K) = 0`) and doing it that way leaves a verified solver
      that every EOS table depends on untouched. Same consumers served, smaller blast radius.
      **The one approximation the shortcut buys is named:** water's own ionisation is not re-solved
      with the seed's electrons present, worth ~1.3% at 15 000 K. If a consumer ever needs better
      than that, this item comes back as originally written.
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

## Q-I(b). Seam B: the absolute level. **DONE 2026-08-22 -- and it is analytic, not cross-kernel.**
`crates/euler2d/tests/wall_pressure.rs`.

**The seam was easier and sharper than billed.** `bounce.rs` normalizes to `rho0 = 1`, `v = 1`,
`p0 = 1/(gamma M^2)`, so `peak_local_pressure` is *already* in units of `rho0 v^2` -- it **is** a
`c_stag`. So Seam B is not a fuzzy cross-kernel comparison at all: the plane-wave limit has a
closed form (the piston problem), and at `gamma = 1.4, M = 40` it is **1.20097**, approaching the
strong-shock `(gamma+1)/2 = 1.2`.

**Measured: 1.33479, i.e. +11.1%.** And it does **not converge away** -- identical to five decimals
from `nz = 80` to `640`, an 8x refinement.

**What it is.** Not wall heating, and not a startup spike. The wall pressure rises *smoothly* to
1.335 at `t ~ 0.08`, decays back **through** the exact 1.201, and settles near 1.164 as the rear
rarefaction arrives. A resolved, converged **impact overshoot** of a finite slug arriving through a
low-density ambient. (Distinct from the 1D kernel's `peak_wall_force` spike at ~`2.0 rho v^2`, which
ADR-0010 correctly diagnosed as an artificial-viscosity artifact.)

**A suspicion that did not pan out.** The kernel's own Noh acceptance test measures "in a band that
avoids ... the axis (the classic wall-heating anomaly)", while `max_plate_pressure` samples exactly
that wall cell. I expected contamination. There is none worth the name: the post-shock plateau sits
at `1.0057x` exact **in the wall cell itself** and `0.999x` a few cells in.

**Does it corrupt `focusing`? No -- and this is the result that matters.** `focusing` is a ratio of
two peak-over-time pressures, so a common overshoot divides out. Measured at `r_foot/R = 0.5,
L/D = 0.3`: **1.2714 from peak-over-time against 1.2627 from the sustained load, 0.7% apart.** Same
cancellation ADR-0003 relies on for `eta_capture`. Pinned by
`focusing_does_not_depend_on_measuring_the_overshoot_or_the_steady_load`.

**Correction to how I reported Seam A.** I wrote that the flat 2D plate "sees the plane-wave load"
to 1.4%. The **number stands** -- flat and confined peaks agree to 1.4% -- but the gloss was
incomplete: *both sides carry the +11% overshoot*, because both are peak-over-time. So Seam A shows
the two are **measured consistently**, not that the flat plate sees the *steady* plane-wave load.
Seam B is what supplies the absolute level, and the two must be read together.

**The one residual optimism, and it is real.** The model computes
`peak = c_stag * rho * v^2 * focusing`, where `c_stag` is the **steady** reflected-shock coefficient
from the 1D kernel. So it prices the steady stagnation load times a geometric concentration, and
**omits the ~11% impact overshoot** the 2D kernel says is physically there. If that overshoot is
real for a delivered pulse, true facesheet peaks run ~11% above what the survivability frontier
classifies. That sits inside the 400 -> 700/900 MPa SiC+Ti margin, and it is the *same direction* as
Seam A's 1.4% -- the model is consistently a little optimistic, never conservative.

**Worth deciding (new, small):** should `c_stag` be redefined to include the overshoot, or should
`P_limit` carry it as an explicit margin? Right now it is carried by neither and is simply absent.
**Cost: ~1 h either way. Worth: an 11% uniform shift in every survivability verdict.**

*Process notes.* Two mistakes cost time and are worth recording. `Grid2D::run_to(t)` resets `t = 0`
internally, so it takes a **duration**, not an absolute time -- two diagnostics mislabelled their
time axis before I checked, which is what made the first "peak at t=0.044" reading look impossible
and prompted the check. And I briefly broke `bounce.rs` by adding a doc comment to
`init_slug_grid`, which was **already public** with one; reverted.

## Q-E. Does the ~5 m cloud length the contour delivers matter, given `e_eff` is read at `L = 10 m`?
Measured and **answered no** at contour densities (spread 0.002-0.006), but the *reason* is that
`tau >> 1` there. The design's own 2026-08-18 correction says `tau ~ 1` at the dilute end. If any
future scenario flies `rho < 0.04`, the `L = 10 m` anchor stops being harmless (spread 0.0155 at
0.04, 0.0517 at 0.01). **Cost: an assert. Worth: cheap insurance, no present error.**

## Q-F. Kerrebrock's decoupling. **CLOSED 2026-08-24: the elevation is computed, and the cliff is not reached.**
> **Closed by ADR-0038.** The section below asks the right question and leaves it open -- "deciding
> needs an electron energy balance (Joule heating against elastic loss), which needs an E-field and
> is a real piece of work". It is now written (`conductivity.electron_energy_balance`), so the
> elevation is an **output** rather than a sweep parameter: **+0 K on the hot legs, +814 K at the
> cold exit.** Read against the table below -- which puts the cliff at 1841-2378 K for +1000 K and
> 2341-2884 K for +500 K -- **every leg exits above its own cliff temperature, so the conductivity
> cliff is not reached anywhere.** The factor-6 discrepancy against the paper's `Rm` column is
> *not* explained by non-equilibrium electrons after all: the elevation that would explain it is
> ~+1000 K and the computed one at bag conditions is ~0 K. **That discrepancy is now unexplained
> and routed back** rather than attributed. The sweep below stands as written; only its status
> changes, from an open possibility to a closed calculation.

**The original entry, 2026-08-22, unchanged:**
`sigma` carries an explicit electron temperature (default `T_e = T_gas`, reproducing equilibrium
exactly), so the question could be asked rather than retrofitted. Raising `T_e` above the gas:

| `T_e - T_gas` | +0 K | +500 K | +1000 K | +2000 K |
| --- | ---: | ---: | ---: | ---: |
| cliff at `v L = 5e3` | 3422 K | 2884 K | 2378 K | **none in range** |
| cliff at `v L = 1.81e4` | 2845 K | 2341 K | 1841 K | **none in range** |

**A ~2000 K electron elevation removes the cliff entirely** from the 1500-15 000 K window. So the
seed window's floor is not a property of the gas temperature alone, and the leak schedule that
rests on it is conditional on electron-ion equilibrium that nobody has established.

**And it may explain the audit's own discrepancy.** The gap the audit found at 3000 K -- 68 S/m
modelled against ~405 S/m implied by the paper's `Rm` column, a factor of 6 -- is almost exactly
what a **+1000 K** electron temperature produces: `sigma` rises 67 -> **390 S/m**. That is close
enough to be worth stating and *not* close enough to be a demonstration. A fixed offset is not a
model. **The honest reading: the factor of 6 is the size a plausible non-equilibrium electron
temperature would produce, so it is not necessarily an error in either calculation -- it may be a
missing term.** Deciding needs an electron energy balance (Joule heating against elastic loss),
which needs an E-field and is a real piece of work. `cliff_temperature(..., t_e_offset=)` exists to
ask the question, not to answer it.

### Q-F(b). Hall parameter and instability screen. **SUPERSEDED 2026-08-24 by ADR-0038.**
> **Both halves of this screen were wrong, in opposite directions.** Kept in full because the
> steady-state finding stands and because the two errors are instructive.
>
> **1. The gain `S` was differentiating the seed alone.** `ionisation_sensitivity` claimed to be
> `d ln n_e / d ln T_e`, but `electron_density` evaluates *water's* contribution at the **gas**
> temperature -- so above ~7000 K, where water supplies most of the electrons, the derivative saw
> almost nothing move. The stabilising headline below, "the gain collapses above ~5000 K ... the
> hot end is safe on **both** counts", **is an artifact and is not real**:
>
> | `T` [K] | 2000 | 3000 | 5000 | 8000 | 11 000 | 15 000 |
> | --- | ---: | ---: | ---: | ---: | ---: | ---: |
> | `S`, seed only (the table below) | 13.3 | 9.05 | 3.51 | **0.11** | **0.01** | **0.00** |
> | `S`, all sources (corrected) | 13.3 | 9.05 | 3.52 | **2.18** | **6.48** | **5.68** |
>
> Water at 13.6 eV against a hotter gas is as sensitive as potassium at 4.34 eV against a cooler
> one -- `chi/(2kT) ~ 10` either way. The hot end is **not** stabilised by gain collapse.
>
> **2. `BETA_CRIT = 2` is the wrong constant for this plume, and it is wrong the other way.**
> Two is the `s -> 0` limit of Velikhov's fully-ionised form, i.e. a *strongly* two-temperature
> plasma. The real `beta_cr` carries `(T_e - T_gas)` in a **denominator**, so it **diverges as the
> plasma approaches equilibrium** -- and Q-M established this plume *is* in equilibrium. Against
> the published criterion the hot legs clear by 2-4 decades, not by a hair.
>
> **Net: the verdict is narrower and firmer than either error allowed.** The screen's own caveat
> -- "it can only rule states *out*, never in" -- was the right instinct and is now discharged:
> `conductivity.electrothermal_loop` is the published criterion and rules states in as well as out.
> **What survives unchanged: the steady-state justification below.** `T_e` relaxing in ~1e-8 s
> against a ~1e-3 s transit is what makes the balance algebraic, and ADR-0038 is built on it.

**Q-F(b). Hall parameter and instability screen. *(original entry, 2026-08-22)***
Before the two-temperature balance is worth writing, two things had to be checked: whether a
steady-state treatment is even valid, and whether a *uniform* one is -- because the electrothermal
(Velikhov) instability filaments a non-equilibrium seeded plasma into hot streamers, and if it
triggers the whole uniform picture fails.

**Steady state is comfortably valid.** Electrons relax by elastic collisions at rate `delta nu`
with `delta = 2 m_e/M = 6.1e-5`:

| T [K] | 2500 | 3000 | 5000 | 10000 |
| --- | ---: | ---: | ---: | ---: |
| `tau_relax` [s] | 5.0e-8 | 4.4e-8 | 1.3e-8 | 8.3e-9 |
| vs the ~200 us expansion | 4000x | 4500x | 15000x | 24000x |

So `T_e` tracks the local balance instantaneously. **No ODE and no transient** -- the balance is
algebraic at every point, which is most of what would have made it hard.

**The instability loop needs two links, and this repo can compute both.** A local rise in `n_e`
redirects current via the Hall effect (link 1, needs `beta = eB/(m_e nu)` large), which raises local
heating and `T_e`, which makes more electrons (link 2, needs the gain
`S = d ln n_e / d ln T_e` large). Break either and the loop opens.

| T [K] | 2000 | 3000 | 5000 | 8000 | 11000 | 15000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| gain `S` | 13.3 | 9.05 | 3.51 | 0.11 | 0.01 | 0.00 |
| `B` [T] for `beta = 2` | 3.36 | 4.73 | 19.7 | 26.4 | 39.8 | 87.1 |

**Both results run opposite to the intuitive guess, and both are reassuring.**
- **The gain collapses above ~5000 K** as the seed saturates -- once every potassium atom is
  ionised, hotter electrons cannot make more. This is a *computable* stabilisation boundary, not a
  literature criterion, which is why it is the part this repository can actually settle.
- **`beta` *falls* with temperature** (0.52 at 2500 K to 0.023 at 15 000 K in 1 T), because Coulomb
  collisions grow faster than the mobility. So the hot end is safe on **both** counts. I had assumed
  the opposite and wrote a test asserting it; the data corrected me.

**At the bag's 0.32 kg/m^3 the plasma is strongly collisional and the Hall link is probably open:**
`beta = 0.42` at 3000 K in 1 T, so closing it needs **~4.7 T**. Diluting reverses this fast -- `nu`
tracks the neutral density, so a tenth of the density needs about a tenth of the field.

**The verdict, stated at the strength it deserves.** The risk sits at the cool end (2500-5000 K),
which is exactly where the cliff and the Q-F question live -- so it is not a comfortable
separation. But it needs several tesla at bag density. **`electrothermal_screen` can only rule
states *out*, never in:** `beta` and `S` are computed exactly, but `BETA_CRIT = 2` is engineering
practice taken on authority. The real criterion is a linearised dispersion relation (Velikhov 1962;
Kerrebrock 1964) whose sources are not available here -- the same `references.bib` gap as Q-G. A
state flagged `at_risk` means "check the literature criterion", not "this filaments".

**Consequence for Q-F.** If the instability does trigger, the effective conductivity is **lower**
than a smooth calculation gives -- the *opposite* direction from the decoupling itself. So the two
effects fight, and quoting the decoupling gain without the instability check would be one-sided.

### **The field is now known, and the loop is CLOSED. 2026-08-22.** **SUPERSEDED 2026-08-24.**
> **This section's verdict is reversed, and the reversal is the main result of ADR-0038.** It
> concluded "the electrothermal instability is a live concern in this design ... live at exactly
> 2000-5000 K", using the bag's own density and `BETA_CRIT = 2`. Both inputs were wrong for the
> question: at bag density the plume is in **equilibrium**, where `beta_cr` is thousands rather
> than 2, so **the bag is not exposed at all**. The table below is arithmetically correct and
> physically inapplicable.
>
> **What replaces it** (`make analysis-electrothermal` -> `data/results/electrothermal_scan.csv`,
> the criterion walked along the cooling history rather than evaluated at assumed temperatures):
>
> | leg | `T` at exit | `T_e - T_g` | `beta_cr` | `beta` | unstable | e-folding |
> | --- | ---: | ---: | ---: | ---: | ---: | ---: |
> | 75 km/s | 16 224 K | 0.4 K | 6 535 | 0.38 | none | -- |
> | 65 | 14 151 | 1.2 | 1 851 | 0.58 | none | -- |
> | 56.53 | 11 681 | 7.9 | 203.9 | 1.20 | none | -- |
> | 45.58 | 4 597 | **813.5** | **2.07** | **4.63** | **1.70 of 2.71 ms** | **3.2 us** |
>
> **The exposure moved from "the whole cool end of every leg" to "the cold leg's second half".**
> That is a real narrowing. But see Q-O: the cold-leg entry is **not robust** -- it rests on the
> current-layer thickness, and the other defensible choice makes every station stable.

**The original entry, 2026-08-22:**
I wrote above that closing the Hall link "needs several tesla" and left it there. The paper states
the field: **~20 T at 1 m down the bore, 12 T at 3 m, 9 T at 6 m, 5 T at exit.** That is above the
threshold everywhere.

| T [K] | `beta` at 5 T | at 9 T | at 12 T | at 20 T | gain `S` | loop |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 2000 | 2.95 | 5.31 | 7.07 | 11.8 | 13.3 | **CLOSED** |
| 3000 | 2.10 | 3.77 | 5.03 | 8.39 | 9.05 | **CLOSED** |
| 4000 | 1.03 | 1.85 | 2.47 | 4.11 | 6.33 | **CLOSED** |
| 5000 | 0.50 | 0.91 | 1.21 | 2.01 | 3.52 | **CLOSED** (at 20 T) |
| 8000 | 0.38 | 0.68 | 0.90 | 1.50 | 0.12 | open (gain gone) |

**So the electrothermal instability is a live concern in this design, not a hypothetical** -- and it
is live at exactly 2000-5000 K, the band where the seed window's floor is decided. Above 5000 K the
seed saturates and the loop opens on its own.

This does **not** mean the plasma filaments: `BETA_CRIT = 2` is engineering practice and the real
criterion is a dispersion relation whose sources are still unavailable (Q-G). The screen rules
states *out*, never in, and it no longer rules these out. **This is now the most valuable open item
on the seed side**, because it bears on the floor the whole leak schedule rests on.

**What is still needed to close Q-F itself: the expansion `u`.** `B` is now known; `u` can be taken
as `v_in/(1+k)` (7.89 km/s at 75 km/s, 5.26 at 50) per Q-G(b), so the balance is now writable. It
was not before.

## Q-G. Both sides of the `sigma` comparison are weakly sourced. **PARTLY CLOSED, honestly.**
**The model side is now explicit rather than hand-picked.** `Q_en = 1e-19 m^2` is a parameter of
`sigma`, not a buried constant, so it can be swept; `sigma` is directly inverse in it wherever
neutrals dominate. `ln Lambda` is computed, and comes out **2.5-2.8** across the window -- which
`coulomb_logarithm`'s docstring says plainly is at the edge of where Spitzer's weak-coupling theory
applies. Together these two make a **factor-of-two claim on `sigma` undecidable by this model**, and
the code says so rather than implying a precision it does not have.

**The paper side is unchanged, and one part of it is now sharper.** `v` and `L` never enter
separately -- only as the product `v L`, since `Rm = mu0 sigma v L`. The paper states neither, so
`magnetic_reynolds` deliberately takes **no default**: a default would be inventing the input that
makes the column reproducible. Back-solving the paper's own `Rm = 361` at 15 000 K against the
audit's conductivity gives `v L ~ 1.81e4 m^2/s`, which is what the regenerated table uses and
labels as an inference, not a citation.

**`references.bib` located 2026-08-22** in `katzseth22202/Balloon-Pulse-Propulsion` (HEAD
`515aca3`, the same commit the audit cites for the "four hundred times" fix). Two things follow, one
useful and one much sharper than expected.

**1. The literature is a weaker match than its filenames suggest.** The bib gives entries, not data --
the two MHD sources are textbooks behind ISBNs (`rosa1968mhd`, `messerle1995mhd`), and
`kerrebrock1964nonequilibrium` is **"I. Theory"**, whose own note records that its worked example is
**argon** with a potassium seed. Argon is a monatomic noble gas with a Ramsauer minimum; water is a
polar molecule whose permanent dipole makes it a far stronger electron scatterer, and `Q_en` is the
input `sigma` is directly inverse in wherever neutrals dominate. So these sources **cannot validate
`Q_en` for steam** -- they validate the seeding practice and the two-temperature law, not the
cross-section. Still open, and now more sharply: what is `Q_en` for **H2O**, not for argon.

**2. The `Rm` column cannot be converted to a conductivity at a single `v L` -- taken that way, four
of its six rows are physically impossible.** The paper's table (`CONTEXT.md`) reads:

| plume T | 15 000 | 6000 | 5000 | 4000 | 3000 | 2000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| paper `Rm` | 361 | 400 | 238 | 76.5 | 9.2 | 0.1 |
| implied `sigma` at `v L = 1.81e4` | 15872 | 17586 | 10464 | 3363 | 404 | 4.4 |
| **Spitzer ceiling** [S/m] | 9954 | 2121 | 1644 | 1101 | 599 | 234 |
| implied / ceiling | **1.6x** | **8.3x** | **6.4x** | **3.1x** | 0.68 | 0.02 |

Spitzer conductivity is the value with **Coulomb collisions only**; adding neutrals can only lower
it, and it is independent of `n_e` except through `ln Lambda`, so **no seed fraction can beat it**.
An implied `sigma` above the ceiling is not a number that can be made right by a better model.

**The likely resolution is benign, and it is still a finding.** Different rows are different *times*
in one expansion, so `v` and `L` plausibly differ per row -- the table is a schedule, not an
isotherm sweep. But then the `Rm` column is **not** a conductivity column, and converting it to one
(which the 2026-08-21 audit did, and which I repeated) is invalid except at whichever single row the
chosen `v L` belongs to. To stay under the ceiling, `v L` must vary by **at least 441x** across the
six rows; against this module's `sigma`, by **13.5x**.

**So the honest status of the "2.3x gap" is: not established.** It was computed by pinning `v L`
from the 15 000 K row and applying it to the others. The comparison is only defensible row by row,
with each row's own `v` and `L`.

**What the paper needs to state, and it is small:** the `v` and `L` behind each row of
`tab:seed_window`. Until it does, no gap between this model and that column is measurable in either
direction -- including the ones I reported earlier in this file, which should be read as
provisional.

### Q-G(b). `v` and `L` derived from the paper's own geometry. **2026-08-22.**
They do not have to be guessed. `Balloon-Pulse-Propulsion`'s `CONTEXT.md` gives the bag as a
**capsule 23 m long with a 3.0 m bore** (`pi r^2 l` = 650 m^3 against the stated 660), the slug as
213 kg (0.323 kg/m^3), the slug ratio `k = 8.5`, and a snowplow transit of ~2.3 ms.

- **Axial `v`** follows from `k`: momentum gives `v_out = v_in/(1+k)`. At 75 km/s that is
  **7.89 km/s**; at 50 km/s, **5.26 km/s**. A 23 m transit at ~10 km/s mean is 2.3 ms, which is the
  paper's own stated figure -- so the geometry closes on itself.
- **Radial `v`** is the plume's sound speed, 1.7-3.8 km/s over 3000-15 000 K.
- **`L`** is the field-gradient scale: 3.0 m (bore radius, the shortest escape path), ~6 m (plume at
  exit -- the field falls 20 T to 5 T, so flux conservation puts the radius at 2x), or 23 m (column).

**So `v L` lies in `[6.0e3, 1.8e5] m^2/s`**, and the audit's back-solved 1.81e4 sits inside it.

**The cliff across that whole range is 2251-3317 K** -- a tight band, because the Saha exponential
makes the cliff only logarithmically sensitive to `v L`.

**The two legs barely differ.** 75 km/s gives 2251-3317 K; 50 km/s gives 2336-3317 K. `v_out`
changes by only 1.5x between them and the cliff moves ~85 K. **Field retention is not a reason to
prefer either leg.**

**Two corners of the range each reproduce one of the paper's own numbers -- and they disagree.**

| what it reproduces | `v L` needed | implied `L`, `v` |
| --- | ---: | --- |
| window floor ~3300 K | **6.0e3** | bore radius 3 m x sound speed 2 km/s |
| `tab:bag_state` leak 4.4% at ~3500 K | **1.0e5** | column 23 m x 4.35 km/s |

At `v L` = 1.0e5 the model gives **4.45% at 3500 K against the paper's 4.4%** -- essentially exact.
But that same `v L` puts the cliff at 2377 K, not 3300 K. **The paper's floor and its leak imply
`v L` values 17x apart.** Both cannot hold with one `v`, one `L`.

**And a single `v L` still does not fit the whole column** (ratios 0.18-2.43 at `v L` = 1.0e5). The
model is high at 15 000 K, low at 5000-6000 K, and **agrees to 10-30% at 2000-3000 K**. That the
disagreement is worst in the middle and vanishes at the cool end supports the trajectory reading:
each row is a different *time*, with its own density as well as its own `v` and `L`. It is also the
reassuring shape -- the rows that agree are the ones that decide the floor.

## Q-H. Ablation is the plate's real exposure, and it is a diagnostic rather than a gate.
2.7-13.6 mm/pulse at 45 km/s rising to 8.3-41.4 mm at 63 (upper bound, nothing credited to
re-radiation or ADR-0014's vapor curtain). Over an 11-cycle chain this is a consumable-mass
problem, and it is **measured** where the `f` collapse was supposed. **Cost to tighten: crediting
re-radiation + the vapor curtain is a real modelling task, not a re-run.** **Worth: high -- this,
not restitution, is the strongest plate-side objection, and the paper currently under-argues it.**

## Q-J. The magnetic nozzle's field retention is coupled to the frozen-recombination bracket.
**Opened 2026-08-22, from a design question of Seth's:** *should the nozzle push while the gas is
still reasonably hot, by keeping the initial steam pressure down?* The direction is right, the
mechanism named is not, and the steam may already solve it without the lever.

**The goal is right.** The leak falls by a factor of 6 between 3000 K and 5000 K (65% -> 10%), so
being on the hot side of the cliff during the push is exactly the design target.

**But initial pressure does not set the exit temperature.** In a nozzle `T_exit/T_0` is fixed by how
much you expand -- the area ratio and `gamma` -- not by the absolute starting pressure. Halving
`P_0` at the same `T_0` and area ratio gives the *same* exit temperature and less thrust. Pushing
**is** cooling: the gas cools because work is extracted from it, and no pressure setting avoids that.

**What lower density does do, and what it costs.** At fixed temperature, diluting genuinely helps
conduction -- and genuinely hurts stability:

| `rho` [kg/m^3] | 0.32 | 0.032 | 0.0032 |
| --- | ---: | ---: | ---: |
| `sigma` at 3000 K [S/m] | 67 | 153 | **235** |
| `beta` at 1 T | 0.42 | 3.1 | **16.3** |

**3.5x better conduction against 39x more Hall drive** -- straight from "instability impossible"
into the at-risk band of Q-F(b). Diluting trades field retention against filamentation; it is not a
free win.

**The good news: the dissociation buffer probably already does the job.** Along an adiabatic
expansion from 8000 K through a 100x density drop, with `v L ~ rho^(-1/3)` at fixed mass:

| `gamma_eff` | `Rm` along the expansion | leak at the end |
| --- | --- | ---: |
| **1.2** | *rises* 16.6 -> 29 | **3.5%** |
| 1.3 | rises, then collapses to 1.1 | 89% |

And `eos_water` says the real value is **`gamma_eff` = 1.10-1.25**, sitting at 1.11-1.16 through the
4000-8000 K band. So the expansion lands in the good column: water dissociation absorbs the energy
that would otherwise appear as cooling, the gas coasts down slowly, and `Rm` actually **increases**
during the push because the cloud grows faster than `sigma` falls.

**So the real lever is not pressure -- it is whether the recombination energy comes back.** That is
**ADR-0026's frozen-recombination bracket**, which this repo already owns for the plate side. In
equilibrium, `gamma_eff` stays ~1.15 and the field is held throughout. If recombination **freezes**,
the buffering stops, `gamma_eff` climbs toward the frozen value, and the expansion moves toward the
1.3 column where the field collapses at the end of the push.

**This is a coupling neither list has**, and it reaches the frozen bracket from a direction nobody
has used: the same question that brackets `e_eff` on the plate also decides whether the magnetic
nozzle keeps its field. Worth stating in the paper as one question, not two.

**Do not over-trust the numbers above.** The expansion is a constant-`gamma` adiabat, not
`eos_water`'s actual isentrope, and `v L ~ rho^(-1/3)` assumes fixed mass and constant expansion
speed. The 1.2-vs-1.3 split is far enough from the real 1.10-1.25 range that the *direction* is
solid, but the **margin is not trustworthy at this level of care**. **Cost to settle: ~1 h** (run it
on the real isentrope through `eos_water`, and take `v(t)` and `L(t)` from the companion repo's
expansion rather than assuming them). **Worth: high -- it converts a design intuition about nozzle
pressure into a statement about a bracket that is already funded and half-built.**

## Q-K. Design levers that keep the electrothermal loop open. **2026-08-22.**
Given Q-F(b)'s finding that the loop is **closed** at 2000-5000 K with the paper's stated field,
the follow-on question is whether geometry or a weaker magnet reopens it. Three levers, with very
different economics. All at `T` = 3000 K, 213 kg slug, unless stated.

**`B` is not a free parameter.** The paper sets it by pressure balance -- "static 10.6 MPa,
ram-to-static 1.17, total ~23 MPa, **~7.6 T**". The field must stand off the plasma at the *hot*
end. What *is* adjustable is the field the plume sees once it is **cool**, since the bore already
expands 20 T -> 5 T. (Sanity check on the geometry: aspect 4 means `l = 8r` and `V = 8 pi r^3`, so
660 m^3 gives `r` = 2.97 m and `l` = 23.8 m -- the paper's 3.0 m bore and 23 m length.)

**Lever 1 -- the exit field. The cheapest, and the margin is 5%.**

| T [K] | 2000 | 2500 | **3000** | 4000 | 5000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `B` below which the loop is open [T] | 3.39 | 3.92 | **4.77** | 9.72 | 19.88 |

The paper's exit field is **5 T** against a 4.77 T threshold at 3000 K -- **closed by ~5%**. A
slightly longer or more divergent nozzle (20 T -> 4.5 T rather than 20 -> 5) opens it. Note how
steeply the threshold climbs with temperature: **the plume is safe while hot**, and the exposure is
only the last stretch as it cools through 2000-3000 K.

**Lever 2 -- bag geometry. It works, and it is a bad trade.** At a fixed 5 T exit field:

| V [m^3] | 200 | 400 | **660** | 1200 | 2500 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `beta` | 0.68 | 1.32 | **2.10** | 3.62 | 6.91 |
| loop | open | open | **CLOSED** | CLOSED | CLOSED |
| leak | 28.9% | 16.8% | **11.5%** | 7.4% | 4.4% |

A **smaller** bag opens the loop but takes the leak from 11.5% to 29% -- trading a possible
instability for a certain confinement failure, and running against the paper's own reason for the
big bag. Reject.

**Lever 3 -- more seed. The only lever that helps both.**

| `x_K` | 0.5% | **1%** | 2% | 5% |
| --- | ---: | ---: | ---: | ---: |
| leak | 15.9% | **11.5%** | 8.4% | **5.5%** |
| `beta` at 5 T | 2.15 | **2.10** | 2.03 | **1.94** |

More seed means more electrons, hence more **Coulomb** collisions, hence a *lower* `beta`. So it
improves retention **and** opens the loop. **The paper's "seed fraction is not a lever" is correct
for retention alone** (`Rm ~ sqrt(seed)`, needing 69x to buy a temperature step) **but it has not
noticed that seed is also a stability lever**, and there it pushes the right way for free.

**Recommendation: more seed plus a little more expansion.** They act in the same direction and
neither is expensive. A weaker magnet on its own does not work, because `B` is set at the hot end.

**The caveat dominates the recommendation.** A 5% margin is far inside the uncertainty on
`BETA_CRIT` itself, which is taken on authority (Q-G). If the true critical value is 1, everything
here is closed and no geometry saves it; if it is 3, everything is already open and none of this
matters. **This says the design sits on a boundary whose location is unknown** -- an argument for
obtaining the Velikhov criterion, not for redesigning around an unsourced number.

> **Resolved 2026-08-24 by ADR-0038 -- and the levers are re-scoped, not confirmed.** The
> criterion was obtained (Petit and Geffray 2009, restating Velikhov 1962 via Petit and Valensi
> 1969). The paragraph above was right that the design sat on a boundary whose location was
> unknown, and right to refuse to redesign around it. But the boundary is not near 1-3 at all:
> `beta_cr` carries `(T_e - T_gas)` in a denominator and runs to **200-6500** on the hot legs.
>
> **The consequence for this section is that all three levers act on the wrong variable.** Exit
> field, bag geometry and seed fraction move `beta`; what binds is `(T_e - T_gas)`, i.e. **how
> cold the plume gets before it leaves the field**. On the exposed leg `beta` is 4.63 against
> `beta_cr` 2.07, and no achievable move in `beta` closes a factor 2.2 -- the seed lever's own
> table above buys 10% in `beta` for 4x the seed.
>
> **The lever that does act on it is `k`.** Sharing momentum with less carried mass leaves the
> plume hotter: `k ~ 4-5` on the cold tail puts the exit above ~11 000 K, where `beta_cr` is ~200
> against `beta ~ 1.2`. Cost ~2-7% Isp on the bare-plate ballistic model -- **which is the tamper
> study's model, not this one's, so that number needs checking in `aim` before it is quoted.**
> Sizing it should wait on Q-O, which decides whether the cold leg is exposed at all.

**It also reduces to a residence time, which is `T(t)`.** The threshold is 20 T at 5000 K and 4.8 T
at 3000 K, so a plume that **leaves the field region before cooling below ~4000 K** never gives the
instability time to grow. Field residence is `l/v` = 23.8/4350 ~ **5.5 ms**. Whether the plume is
still above 4000 K at that point is exactly the cooling history the companion repo asked this repo
for, and which is still unbuilt. **`T(t)` is now the item the most threads are waiting on:** Q-F,
Q-K and item 10's quadrature all reduce to it.

> **Superseded 2026-08-22 by Q-L.** The 5.5 ms above is wrong: it divides the bore length by an
> assumed constant 4.35 km/s, and the plume **accelerates** through the nozzle. The solved
> residence is **1.69-2.71 ms**. The error ran against the design, not for it -- less residence
> means less time for the instability to grow -- and Q-L narrows this section's exposure to a
> single corner, the cold leg on the frozen branch. The three levers and the `BETA_CRIT` caveat
> below are unaffected; only the clock changes.

## Q-L. The cooling history `T(t)`. **BUILT 2026-08-22.** The bracket is the answer.
`python/puffsat/expansion.py`, `make analysis-expansion`, `data/results/cooling_history.csv`.

**The model, and the one assumption in it.** Quasi-1D steady isentropic expansion of the
post-impact plume, parametrised by density so nothing is root-found on the nozzle: the adiabat
`de = -p d(1/rho)` is integrated on the real `eos_water` EOS (which publishes no entropy
function, so the adiabat *is* the path), `u` follows from `h + u^2/2 = h0`, and the area from
`rho u A = const`, referenced to the sonic throat. **The paper's `20 T -> 5 T` is `A/A* = 4` by
flux conservation with no further assumption** -- so `T` against area ratio is assumption-free.
Only the *clock* needs a shape for `A(x)` along the bore, and a linear opening over 23.8 m is
assumed there. Consumers that can work in area ratio are free of it.

**It runs on both branches of ADR-0026, and that turns out to be the finding.**

| closing speed | `T` exit, equilibrium | `T` exit, frozen | residence | `v L` |
| ---: | ---: | ---: | ---: | ---: |
| 75 km/s | 16 224 K | 5 297 K | 1.69 ms | 9.7e4 |
| 65 | 14 151 K | 4 528 K | 1.95 ms | 8.5e4 |
| 56.53 | 11 681 K | 3 922 K | 2.22 ms | 7.4e4 |
| 45.58 | 4 597 K | 2 972 K | 2.71 ms | 6.0e4 |

**1. The residence time is 1.7-2.8 ms, not 5.5 ms.** Q-K divided 23.8 m by an assumed 4.35 km/s.
The plume does not travel at a constant 4.35 km/s -- it **accelerates**, 8.0 -> 16.2 km/s on the
hot leg and 5.0 -> 9.9 km/s on the cold one, because that is what a nozzle does. Halving the
residence runs *in the design's favour*: the instability has less time to grow. Independent
corroboration: the cold leg's 2.71 ms is the paper's own stated ~2.3 ms snowplow transit, which
was not an input to this calculation.

**2. Q-K's question, answered.** The threshold field is 9.7 T at 4000 K against the 5 T flown, so
"is the plume still above 4000 K when it leaves the field" decides it.

- **75 and 65 km/s: yes on both branches.** No exposure.
- **56.53: yes in equilibrium (11 681 K); marginal frozen (3 922 K).**
- **45.58: yes in equilibrium (4 597 K); no if recombination freezes (2 972 K).**

So **the exposure is one corner -- the cold leg, on the frozen branch -- and nowhere else.** That
is a much narrower claim than Q-K could make, and it is narrower in the safe direction.

**3. `v L` is now an output.** Q-G's finding was that no gap between this model and
`tab:seed_window`'s `Rm` column is measurable, because the paper never states the `v` and `L`
behind it. The solve supplies both: local flow speed times local flux-tube radius, giving
`v L` = **5.5e4 to 9.7e4 m^2/s**. Against Q-G(b)'s two disagreeing corners:

| what it reproduces | `v L` needed | against the computed 5.5e4-9.7e4 |
| --- | ---: | --- |
| `tab:bag_state` leak 4.4% at ~3500 K | 1.0e5 | **at the top of the range** -- corroborated |
| window floor ~3300 K | 6.0e3 | **10x below the range** -- not reachable |
| the 2026-08-21 audit's back-solve | 1.8e4 | **3x below the range** |

**The paper's leak number survives and its window floor does not.** The `Rm` column cannot be an
isotherm sweep at 6.0e3; the geometry does not produce that product anywhere in the expansion.

**4. Item 10's quadrature, run here rather than routed.** Residence-weighted `1/Rm`, the exact
thing item 10 asks for ("weight `1/Rm(T)` by how long the plume spends at each `T`"):

| closing speed | 75 | 65 | 56.53 | 45.58 |
| --- | ---: | ---: | ---: | ---: |
| leak, equilibrium | 0.11% | 0.17% | 0.29% | **2.54%** |
| leak, frozen | 1.21% | 1.80% | 2.58% | **5.38%** |

**The paper's published 4.4% sits inside this bracket**, near the cold-and-frozen corner. It is
not reproduced by any single row, which is the point: it is a trajectory average, and averaging
over the trajectory is what produces a number of that size.

**5. Radiation was checked, not assumed away.** The isentrope is adiabatic, which is a claim. The
plume is diffusion-limited while hot (Rosseland depth up to 420 across the bore) and the radiated
fraction of internal energy over the whole transit is **1.3% to 13.2%**, worst at 56.53 km/s in
equilibrium. So adiabatic holds -- but it is an **underestimate of the cooling by up to 13%**,
always toward a colder exit, so the equilibrium column above is an upper bound.

The resolution is set by this integral rather than by `T(t)`: the loss is concentrated in a narrow
window around the opacity crossover, and a coarse sampling that resolves the temperature history
perfectly well understates the radiated fraction by ~3x. At the shipped settings it is converged
to three digits (0.1324 against 0.1326 at double the stations). The exit temperatures and the leak
column do not move at all between those resolutions.

**6. A trap worth recording, because I fell into it.** The first version of the radiation check
branched on the **Rosseland** optical depth and then computed the emission from the **Planck**
mean. TOPS puts `kappa_P/kappa_R ~ 100` for water at 9000-15 000 K, so states that read
"optically thin" by Rosseland are still deeply thick to their own lines, and the loss came out
two to three orders of magnitude too high. It reported the plume radiating away **five times its
own internal energy**, and the tell was that the number **grew without bound under grid
refinement** rather than converging. Fixed by taking the flux-limited minimum of three rates --
emission, free-streaming, diffusion -- which is ADR-0006's own convention. The regression test is
the invariant it violated: nothing radiates faster than a blackbody at its own temperature
through its own surface.

**What this does not settle.** Which branch of the bracket is real. That is a factor 3 in exit
temperature and it spans the entire cliff, so it is now the largest single uncertainty on the
nozzle side -- exactly as ADR-0026 found it to be on the plate side.

> **Settled the same day, in Q-M: the equilibrium branch.** Recombination beats the expansion by
> 2.2 to 5.2 decades wherever the ionisation store still holds anything, so the equilibrium column
> above is the answer and the frozen column is not reached in the nozzle. Two consequences flow
> back into this section: the instability exposure Q-L narrowed to one corner is **closed**, and
> item 10's leak is the equilibrium row alone, **0.11%-2.5%** -- which puts the paper's published
> 4.4% above the model's whole range rather than inside it. See Q-M.

## Q-M. Which branch of the frozen bracket is real? **ANSWERED 2026-08-22: equilibrium.**
Everything above brackets rather than answers, and the bracket is a factor 3 in exit temperature,
2.5x in leak, and the difference between "no exposure anywhere" and "the cold leg is exposed".
**The criterion is a rate comparison and it is computable**: recombination rate against expansion
rate, `tau_rec` vs `rho/(d rho/dt)`, along the history this repo now owns. Neither repo has done
it, and both are currently quoting brackets because of it. **Cost: moderate** -- it needs
three-body and radiative recombination coefficients for H and O, which are tabulated. **Worth:
the highest on the list**, because a single calculation collapses a bracket that two separate
studies are each carrying at full width.

### Q-M answered 2026-08-22. **The bracket collapses to equilibrium. The nozzle does not freeze.**
`python/puffsat/recombination.py`, `make analysis-recombination`. No new artifact -- the answer is
a verdict, not a table.

**Method.** Race the two clocks. `Da = tau_expansion / tau_recombination`, the classical Bray
criterion; `Da >> 1` means the chemistry keeps up and the equilibrium branch is right. Run the
*equilibrium* history and ask at every station whether the chemistry was fast enough to have
produced it -- a self-consistency test, which is why no finite-rate solver was needed.

| closing speed | binding `Da` (ionisation) | margin | verdict | `f_diss` at exit |
| ---: | ---: | ---: | --- | ---: |
| 75 km/s | 1.6e6 | 5.2 decades | **equilibrium** | 1.0000 |
| 65 | 8.2e5 | 4.9 | **equilibrium** | 1.0000 |
| 56.53 | 1.8e5 | 4.3 | **equilibrium** | 1.0000 |
| 45.58 | 1.5e3 | 2.2 | **equilibrium** | 0.9931 |

**The store returns, and it returns with orders of magnitude to spare.** Three-body recombination
(`X+ + e + e -> X + e`) dominates the radiative channel by ~1e5 here, because the plume is dense
-- 1e25 electrons/m^3, about a fifth of atmospheric number density. At 75 km/s the ionisation
store falls 83% -> 54% of the reservoir value across the transit, every bit of it at `Da > 2e6`.

**The statistic matters, and the obvious one inverts the answer.** `min Da` over the history is
**0.0093** on the cold leg, which reads "frozen". That minimum sits at 5070 K, at a station
holding **0.01% of the ionisation store**. A store that has already been returned cannot freeze.
Weighting by what is actually still held gives 1.5e3, and the two differ by five orders of
magnitude. Recorded because it is an easy way to get this exactly backwards.

**A separate finding that fell out: the dissociation store never returns in the nozzle at all.**
Equilibrium water is **fully dissociated at every station** (`f_diss` = 1.0000 to 0.9931), so
there is no dissociation energy to give back on either branch, and the factor-3 bracket was
purely the ionisation store all along. Its Damkohler number is marginal (`min Da_diss` = 1.3 to
10.8) -- but it does not bind *here*, because the reaction equilibrium is not asking for it. It
binds **downstream, in the fireball**, which is the still-open "does recombination freeze below
0.01 kg/m^3" item. **That item now has its rate comparison already built.**

> **Settled 2026-08-24, and it binds sooner than "downstream" suggested: at the lip.** See Q-P.
> The three hot legs cross `Da = 1` at the **first station past the nozzle exit**, with **100% of
> the dissociation store still held**, stranding 50.9 MJ/kg -- 19-34% of the dissipated budget.
> This is the mirror image of the ionisation case above: there the smallest `Da` sat on a spent
> store and the freeze was an artifact; here the store is full when it freezes.

**Consequences, and one of them cuts against the paper.**

1. **Q-K's instability exposure is closed.** The only exposed corner was the cold leg on the
   frozen branch (2 972 K). The equilibrium answer is 4 692 K, above the 4000 K line. **No leg is
   exposed on either the hot or cold end.**

   > **Partly overturned 2026-08-24 by ADR-0038, against a criterion this answer did not have.**
   > The "4000 K line" is Q-K's own `beta > 2` screen, and `beta_cr = 2` does not transfer to an
   > equilibrium plume. Under the published criterion **the hot legs clear by 2-4 decades**, which
   > is a much stronger statement than clearing a line by 700 K -- so this consequence *understated*
   > its own case there. **The cold leg does not clear it.** Its exit sits at `beta/beta_cr = 2.2`
   > and it is unstable over 1.70 ms of its 2.71 ms transit. Equilibrium made the cold leg *hotter*
   > and it is still the exposed one, because what binds is the electron-gas temperature gap rather
   > than the gas temperature itself. **Whether it is genuinely exposed is undecided -- see Q-O.**
2. **Item 10's leak is the equilibrium column: 0.11% to 2.5%**, not the 0.11-5.4% bracket.
3. **The paper's published 4.4% is now *above* the model's whole range.** Q-L had it sitting
   inside the bracket; with the bracket collapsed it does not. The model says the leak is roughly
   **2x smaller** than published. That is a discrepancy in the *conservative* direction -- the
   paper is under-claiming -- but it is a discrepancy and it should not be quietly absorbed.
   Two readings: either the paper's 4.4% carries margin it does not state, or `v L` is smaller
   than the nozzle solve says. **Routing back rather than resolving.**

**What this does NOT settle: ADR-0026's plate bracket.** Everything above is the *nozzle*, whose
transit is ~2 ms. The plate bounce is ~**1 microsecond**, a thousand times shorter, at a different
density history. Both `Da` channels scale as density squared, so the plate is *not* simply worse,
and it cannot be assumed either way -- but nothing here transfers to it. **The module now exists
to check it**, which is a much smaller job than building it was, and it is the obvious next use.

**Honest weak point.** The three-body atomic coefficient (`6.1e-26 T^-2 cm^6/s`, Baulch et al.)
carries a factor 2-3, and steam's third-body efficiency differs from the argon/nitrogen the fit is
anchored on. That uncertainty is ~0.5 decades against the 2.2-5.2 decades of margin, so it does
not threaten the ionisation verdict -- but the *dissociation* numbers above sit inside it, which
is exactly why they are reported as marginal rather than as an answer. The ionic coefficient's
published form is quoted with `T` in **eV**; reading it as kelvin changes the rate by ~1e11 and
inverts the verdict, so it is pinned in the tests by an independent Thomson estimate as well as by
the literal.

## Q-N. The nozzle exhaust speed is not the momentum-sharing speed, and the paper uses one number.
The solve gives an exhaust of **9.9-16.2 km/s** relative to the ship. The paper's `v_out =
w/(1+k)` gives **5.26-7.89 km/s**. These are different quantities -- one is the merged slug's bulk
speed from momentum sharing, the other is what the nozzle converts thermal energy back into -- and
both are correct for what they describe. But **effective Isp depends on which one leaves the
vehicle**, and the paper does not distinguish them. Flagging rather than resolving: this is the
tamper study's regime (`puffsat_tamper_isp_prd.md`), not this one's. **Cost: low** to state.
**Worth: high** if the Isp claim rests on the smaller number, since the nozzle solve says the
larger one is available.

> **2026-08-24: this now matters more, because Q-P changed what is at stake.** The nozzle exhaust
> leaves **chemically frozen**, carrying 19-48% of the dissipated energy away as inert bond
> enthalpy. If the paper's Isp rests on `w/(1+k)`, that costs it nothing -- it never claimed the
> thermal energy. If it rests on the nozzle exhaust speed, the frozen-flow loss is charged against
> exactly that claim. **The question is unchanged and the answer is now load-bearing.**

## Q-O. What thickness does the driving current occupy? *(decides the cold leg, undecided here)*
**Opened 2026-08-24 by ADR-0038 Addendum 3, on Seth's question about where a `T_e > T_g` gap comes
from at all in a plume that is neither driven nor strongly two-temperature.** The question is well
founded, and working it downgraded the cold-leg finding **from a result to a model-dependent
possibility.**

Addenda 1 and 2 of ADR-0038 eliminated every other candidate source of the gap. The expansion
does not outrun electron-heavy equilibration (it loses by ~1e6). Three-body recombination
depositing its released energy into the third electron is real -- 4.7e9 W/m^3 at the throat -- but
spread over milliseconds against a nanosecond redistribution time it is worth **fractions of a
kelvin**. Both would have *raised* `T_e` and so run against the design, and both are simply too
small. So the gap is set entirely by the field-gradient current, `Q_joule = j^2/sigma` with
`j = B/(mu0 L_eff)` -- and **`L_eff` is estimated by a scaling, not solved.**

| station | `L_eff` reading | `L` [m] | `T_e - T_g` | `beta_cr` | `beta` | verdict |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| crossing (7569 K) | skin, `sqrt(t/(mu0 sigma))` -- shipped | 1.72 | 456 K | 3.86 | 3.93 | **unstable** |
| crossing | full flux-tube radius | 4.21 | **77 K** | **18.1** | 3.93 | **stable** |
| exit (4597 K) | skin -- shipped | 2.12 | 814 K | 2.07 | 4.63 | **unstable** |
| exit | full flux-tube radius | 6.00 | **79 K** | **10.6** | 4.74 | **stable** |

**A factor 3 in thickness flips every unstable station**, because it is 9x in `Q_joule`, hence 9x
in the elevation, hence 9x in `s`. Reproduce both columns: `make analysis-electrothermal`.

**A physical argument points the stable way, and it is not weak.** Plasma beta at the exit is
`2 mu0 p / B^2` ~ **0.016** -- magnetic pressure exceeds gas pressure ~60x. A plume that cannot
appreciably distort the field does not concentrate its current into a thin resistive skin; the
diamagnetic current is small and spread over the pressure-gradient scale, which is the flux-tube
radius. The skin estimate is retained as shipped because it is the *conservative* one, **not
because it is believed.**

**Cost to settle: this is the one genuinely expensive item on the list.** The current distribution
in a low-plasma-beta magnetic nozzle is an MHD problem neither `expansion.py` nor this repository's
kernels solve -- quasi-1D cannot see it, because the whole question is radial. It is a 2D
axisymmetric resistive-MHD solve on a prescribed field, which is a new kernel rather than a new
consumer of an existing one.

**Worth: it decides whether a 2-7% Isp cost is owed.** If the cold leg is exposed, Q-K's `k` lever
is the fix and has to be sized; if it is not, nothing is owed and ADR-0038's design consequence
falls away entirely. Note the asymmetry -- **the hot legs are stable under every variant tested**,
clearing even the `s = 0` floor, so no answer here reaches them. This bears on one leg only.

**Honest statement of where it stands: I don't know.** The cold leg is the only exposed candidate;
whether it is actually exposed is undecided, and the cheaper physical argument says it is not.

## Q-P. The exhaust leaves chemically frozen, and that is a frozen-flow loss nobody has charged.
**Answered 2026-08-24** while closing the fireball item. `python/puffsat/fireball.py`,
`make analysis-fireball`. This is the downstream half of Q-M: that answer established the nozzle
stays in equilibrium and noted the dissociation store's Damkohler number was marginal and would
bind further out. It binds **at the lip.**

**Method, and it is the same self-consistency test Q-M used.** Continue the equilibrium isentrope
past `A/A* = 4` into the free jet, and ask at every station whether atomic three-body recombination
was fast enough to have produced the composition equilibrium claims. Inside the field the answer is
yes (Q-M). One station outside, it is no.

**Two things happen at the lip, and the first is the trigger.** The paper's bore opens `A/A*` 1 to
4 over 23.8 m -- about **7.9 m of travel per unit area ratio.** A 45-degree free jet covers the same
unit in **0.75 m.** So the local expansion time steps down ~8x at the boundary (2.66 ms just inside
to 0.33 ms just outside) and `Da` falls from ~35 to ~1.8 in a single step. **The magnetic nozzle
was holding the plume in a slow expansion; the freeze begins where the field lets go.** Then the
rates finish it: three-body recombination is a density-*cubed* process, so each further decade of
density costs two decades of rate while `tau_exp` grows only slowly.

| leg | `rho` at freeze | `T` | `f_diss` held | stranded | of dissipated budget |
| --- | ---: | ---: | ---: | ---: | ---: |
| 75 km/s | 2.35e-2 kg/m^3 | 16 063 K | **1.0000** | 50.9 MJ/kg | **19.2%** |
| 65 | 2.29e-2 | 13 974 | **1.0000** | 50.9 | **25.6%** |
| 56.53 | 2.10e-2 | 11 271 | **1.0000** | 50.9 | **33.9%** |
| 45.58 | 1.01e-2 | 3 908 | 0.9170 | 46.7 | **47.7%** |

**This is the case `binding_damkohler` was built to tell apart, and it comes out the other way.**
In Q-M the smallest `Da` sat at a station holding 0.01% of the store, and the freeze was an
artifact. Here the crossing station holds **100%** on three legs of four. The freeze is real and
it strands the bond energy.

**The verdict is robust to the one free parameter.** `tan theta` scales the clock linearly, so
`Da ~ 1/tan theta`. Over **5 to 75 degrees -- a 43x range in `tan`** -- the cold leg's freeze
density moves only 12x (1.9e-3 to 2.2e-2) and the stranded energy moves **40.3 to 50.0 MJ/kg**, a
24% spread. On the hot legs `theta` does not move the answer at all, because the crossing is
already at the first station past the lip. So the freeze *density* is conditional on the jet
geometry and the *freeze* is not.

### What this means, stated at the strength it deserves
**It is not, by itself, a contradiction of the paper.** It is a constraint on one specific claim:
that the nozzle converts the dissipated thermal energy into directed exhaust. **19% to 48% of that
energy leaves as inert chemical enthalpy instead** -- classic frozen-flow loss, and it is not
charged anywhere in the ledger I have seen.

**Whether it costs the paper anything depends on Q-N**, and the two questions are now coupled. If
the paper's Isp rests on `v_out = w/(1+k)` -- the momentum-sharing speed -- then it never claimed
the thermal energy at all and frozen flow costs it nothing. If it rests on the nozzle exhaust
speed, this is a 19-48% haircut on the energy available to produce it. **Q-N asked which number
the paper uses and did not resolve it. It now matters more than it did.**

### Honest weak points, and one of them points the wrong way for the design
- **`eos_water` carries no OH.** The real dominant path is `H + OH + M -> H2O + M`, and the model
  has only `H` to work with. The code comments say plainly that using `n_H` **overestimates** the
  rate where OH is scarce -- so the true freeze is *earlier* than modelled, not later. The finding
  is conservative in the direction that matters, but the number is soft.
- **The three-body atomic coefficient carries a factor 2-3** (Q-M's own recorded weak point).
  A factor 3 in rate moves the freeze density by ~1.7x, well inside the `theta` bracket.
- **Past the crossing the equilibrium history is not the real one.** Composition holds near where
  it froze and the temperature stops being buffered, so the real plume runs *colder* than the
  curve past that point. That does not move where the freeze starts, which is all that is claimed.
- **What would settle it properly: a finite-rate calculation with OH, H2 and O2 in the species
  set.** That is a real piece of work -- a new species set in `eos_water` plus a reaction network,
  not a consumer of an existing kernel. **Cost: high. Worth: it converts a 19-48% bracket into a
  number, and the bracket currently spans "a correction" to "half the budget".**
