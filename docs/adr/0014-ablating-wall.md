# The ablating wall is a Phase-2 e_eff lever, not just protection; the rigid wall is the conservative floor

The survivability frontier runs on a passive ablating/transpiring wall (Phase 2). Counter to the
intuition that ablation only "costs" energy, it tends to *raise* `e_eff` — so the rigid wall
(ADR-0005) is the conservative floor and the ablating wall is the best-estimate refinement.

**Why ablation raises `e_eff`.** The vapor acts on two loss channels: **blowing** (injected vapor
thickens the boundary layer and cuts the conductive flux into the solid) and **vapor shielding**
(the vapor absorbs incoming radiation in the near-wall layer before it reaches the cold solid).
Both keep energy *in the gas* instead of losing it to the wall. The only true energy sink is the
ablation enthalpy `Q*·ṁ`, and `ṁ` is small — a few µm of ablator ≈ 0.4 kg vs the 25 kg pulse,
~1.5%. The recovery outweighs the cost, so `e_eff` rises relative to the rigid wall.

**Ties to the dual-curve and the transition.** Rigid wall = conservative floor for *both* `e_eff`
and survivability (full wall losses, full flux to the SiC); ablating wall = best-estimate. Per
ADR-0013, ablation is a gated refinement (rigid → lower-bound curve, ablating → best-estimate curve),
run where the floor dips below `f = 0.8` (ADR-0009). Its `e_eff` leverage should *peak at the
transitional anchor*: at 3.2 km/s the gas barely radiates, at 16 km/s the radiation is already
trapped (`τ≫1`), but at `τ~1` radiation reaches the wall and vapor shielding has the most to recover
— so the ablating wall is the prime candidate to fill the transitional `e_eff` dip (ADR-0012).

**Model.** Quasi-steady surface energy balance: incoming flux → ablation rate via an effective heat
of ablation `Q*`, with vapor injected as a **mass/momentum/energy source at the 1D wall boundary**.
Three coupled effects, all self-consistent in the existing 1D rad-hydro: blowing correction on the
conductive flux; vapor opacity in the FLD for shielding; injected mass in the bounce. **Passive
ablation** of the MEMS-renewed sacrificial layer — not pumped coolant; renewal is between pulses,
passive during the pulse.

**Vapor approximated, not speciated (baseline).** The vapor is treated as a `Q*` enthalpy sink + a
gray near-wall absorber + a mass source, without full silicone/ionic-liquid pyrolysis chemistry —
even though the ablating wall injects ablator species the water-only table pipeline (ADR-0007) does
not cover. Full ablator-vapor speciation is gated on whether spot-checks show it moves `e_eff`. `Q*`
is parameterized (silicone ~2–10 MJ/kg literature), with sensitivity reported.

**RT-coherence caveat.** The vapor-shield and blowing recovery assume the near-wall layer stays
coherent. Rayleigh–Taylor mixing (ADR-0020) is the mechanism that erodes it — hot gas mixed down
through the shield reaches the solid early. So this `e_eff` recovery is an *upper estimate* to the
extent the shield stays RT-coherent; the ablating wall is where RT bites hardest, not the rigid floor.

## Considered Options

- **Transient / full-pyrolysis ablation** (char layer, decomposition kinetics, gas-phase products)
  from the start. Deferred: quasi-steady `Q*` captures the first-order energy interception; full
  pyrolysis only if spot-checks show vapor chemistry matters.
- **Treat ablation mass-loading as a first-order perturbation needing full vapor chemistry.**
  Rejected as baseline: `ṁ ~ 1.5%` of the pulse; the dominant effect is on the loss channels
  (blowing, shielding), captured with an effective `Q*` + gray absorber.

## Outcome (Rung E, 2026-06) — landed, with three corrections to this ADR's framing

The `AblatingBounce` kernel landed (1D Lagrangian, `wall = None`): a quasi-steady `Q*` surface
energy balance with a wall **mass source** `ṁ = q_in/Q*`, a **blowing factor** `φ = 1/(1+B)` on the
conductive flux, and a **vapor shield**. Each is verified against its off-limit (`Q* → ∞` recovers
the rigid floor with mass/energy closure; `φ → 1` recovers the bare conductive flux; `κ_vapor → 0`
recovers the bare wall; shielding monotonically cuts `loss_radiative_wall` and the recovery tracks
the recoverable loss — the τ-leverage). The `--ablating` sweep × `analysis.py --axis ablating`
deliver the recovery as a **τ-bracket** (opacity-scale knob) against the rigid floor, folded through
the Σ→ρ→peak survivability frontier (Rung S). Three findings amend the framing above:

1. **The dip is *not* radiatively fillable — this ADR's central "leverage peaks at the transitional
   dip" thesis does not hold for the dip's worst case.** The thesis assumed the dip is a radiative
   `τ~1` leak the shield can recover. But the measured 0.57 dip (ADR-0012) is an **EOS dissociation
   specific-heat sink**, not a wall radiative loss — energy goes into breaking bonds, not into
   radiation headed at the plate — so the shield has almost nothing to claw back there: ablating
   lifts the 11 km/s `e_eff` from ≈ 0.566 only to **[0.570, 0.580]** (recovery ≤ +0.014). The shield
   *can* recover a real radiative loss (its τ-leverage is confirmed where one exists), but the dip
   is not where that loss lives.

2. **The recovery is mass-injection-dominated, not shielding-dominated** (refining the "two channels
   = blowing + shielding" model). **Blowing is null at both science anchors** — the high-v table
   carries no `k_gas` and low-v conduction is negligible (Rung C), so it is a verified-and-bounded
   channel, not a live lever. Of what remains, mass injection leads: at 16 km/s the transparent
   (no-shield) floor is already `+0.017…+0.028`, with the shielding τ-leverage adding only `~+0.017`
   on top (~38 % of the recovery). So the live channels at the anchors are **injection (primary) +
   shielding (secondary)**.

3. **The shield is a throttled-conductance BC, and `ṁ` is larger than the ~1.5 % assumed.** The
   shield is implemented as `RadBc::MarshakAttenuated` scaling the Marshak surface conductance by
   `1/(1+τ_v)` and *retaining the intercepted radiation in the near-wall field* (energy-conserving,
   self-consistent) — chosen after a χ-augmentation attempt failed (the Marshak flux is a surface
   conductance independent of cell-0 opacity) and superseding the manual energy-return first
   committed. The measured ablated fraction at 16 km/s is **~3.7 % (`Q* = 10 MJ/kg`) to ~8.9 %
   (`Q* = 2 MJ/kg`)** — well above the ~1.5 % estimate, and at the aggressive end it strains the
   thin-curtain assumption *and* plate durability (the back-propagated MEMS-replenishment
   requirement, design §7).

**The recovery-lever decision (was deferred to this rung):** the ablating wall does **not** robustly
lift 16 km/s over the `f = 0.8` gate. Best survivable `f` 0.784 → **[0.788, 0.807]**: the gate
**straddles** the bracket, cleared only at the optimistic low-`Q*`/high-τ corner (`f` 0.807, at
~8.9 % ablation), with the conservative `Q* = 10 MJ/kg` end at 0.788, just under. It clears
comfortably (0.859) only if the structure tolerates the relaxed 900 MPa `P_limit`.

*Corrected (2026-07):* with the physical stagnation coefficient (`c_stag ≈ 1.24`, not the AV-artifact
2.0; ADR-0010 correction) and the converged 112×80 geometry grid, the 16 km/s rigid baseline itself
clears at `f ≈ 0.805`, and the ablating bracket becomes **[0.809, 0.829] — the gate clears across the
whole `Q*`/τ bracket**, no longer straddling; the relaxed 900 MPa limit adds almost nothing (0.830)
because the corrected peak pressure barely binds. The physics conclusion of this ADR (recovery real
but small and injection-dominated; the dip not radiatively fillable) is unchanged. **E5** (real
per-regime opacity) was **not pulled**: the decision is `Q*`/EOS-limited, not τ-limited (real opacity
sharpens only the ~38 % shielding sub-component and cannot fill an EOS sink), and the data
(HITEMP/ExoMol/TOPS/OPLIB) is firewall-blocked by default-deny — see ADR-0007. The **RT-coherence
caveat above still stands**: the shielding component of the recovery is an upper estimate.

*Further correction (2026-07-10, ADR-0023 kernel fix — concave `eta_capture` corners were ~1 % high):*
the 16 km/s rigid baseline becomes `f ≈ 0.798` and the ablating bracket **[0.802, 0.821]** — the gate
still clears across the `Q*`/τ bracket, though only barely at the conservative `Q* = 10 MJ/kg` end.
The physics conclusions are unchanged.
