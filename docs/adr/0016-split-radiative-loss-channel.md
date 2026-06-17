# Split the radiative loss into to-plate and escape-to-space — the two halves have different recovery levers

The §6 loss decomposition itemized one radiative channel ("radiative-to-plate"). That is too coarse:
the two recovery levers in the study — the ablating-wall vapor shield (ADR-0014) and the bulk soot
seed (ADR-0015) — act on **two physically distinct halves** of the radiative loss, and a single
channel hides which lever can help and corrupts the seed gate.

**The 1D geometry forces the split.** The hot stagnation layer sits against the wall; temperature
falls outward. Radiation emitted by that near-wall gas divides:
- **(1a) radiative-to-plate** — the half heading at the cold absorbing wall (ADR-0005). It traverses
  the *hottest* gas, which holds no condensed soot (ADR-0015's ceiling), so the **seed cannot touch
  it**. The **ablating-wall vapor shield can** — it sits exactly at the wall and absorbs incoming
  radiation in the near-wall layer (ADR-0014).
- **(1b) radiative-escape-to-space** — the half heading outward through the cooler outer gas to the
  open (non-wall) far boundary. Optically thin, it escapes and cools the gas, lowering rebound
  pressure (an `e_eff` loss). The cooler outer gas *can* hold soot, so the **seed can recover it**
  (τ≫1 there → reabsorbed → retained as rebound pressure). The wall vapor shield cannot reach it.

**Consequence — the two levers are complementary, not redundant.** Ablation recovers (1a); the seed
recovers (1b). Each covers the half the other physically cannot. Both gates sharpen: the **seed gate
keys on (1b)**, not on the old lumped channel (which was dominated by (1a), the part the seed can't
help); the ablation gate keys on (1a).

**Decomposition becomes five channels** (per velocity anchor): (1a) radiative-to-plate, (1b)
radiative-escape-to-space, (2) conductive-to-plate, (3) condensation/recombination, (4) sideways
escape / non-capture (geometric, the 2D `eta_capture` track). The 1D track already needs a
far-boundary BC to expand into; (1b) just itemizes the energy leaving there.

## Considered Options

- **Keep one lumped "radiative" channel.** Rejected: it hides which of the two levers can act, and
  points the seed gate at the to-plate fraction the seed cannot recover.
- **Track only radiative-to-plate, treat far-boundary escape as negligible.** Rejected: the outward
  escape is precisely the seed's recovery target; dropping it pre-judges the seed study to null.
