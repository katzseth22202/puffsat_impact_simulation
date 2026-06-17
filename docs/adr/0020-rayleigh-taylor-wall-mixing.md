# Rayleigh–Taylor mixing at the wall boundary layer is out of scope, gated on the conductive channel

A reviewer who knows Project Orion's follow-on hydrocode work will ask whether Rayleigh–Taylor (RT)
instability at the pusher-plate plasma boundary degrades `f`. It does not, to first order, for a
momentum-conservation reason — so RT is not modeled, with one gated exception tied to the
conductive-to-plate loss channel. §3.2 leans on Orion, so this is answered explicitly rather than
left implicit.

**Naming.** The unstable interface is specifically the **near-wall boundary layer / ablation front** —
the hot, low-density gas at the wall supporting the denser stagnated bulk against the stagnation
deceleration (∇ρ·∇P < 0) — not the bulk contact surface. The shock-impact phase is technically
**Richtmyer–Meshkov** (impulsive) seeding the **Rayleigh–Taylor** (sustained-deceleration) growth.

**RT cannot move the impulse; it can only move a loss channel.** `e_eff` is a pure momentum ratio
(ADR-0001): `J_wall = ∫P_wall·A dt`, fixed by the global momentum balance — the gas arrives with
`p_in` and is turned around, `J_wall = p_in(1+e_eff)`, whether the stagnation layer is laminar or
RT-turbulent. RT corrugates the interface (more area, faster mixing) but neither creates nor destroys
axial momentum. So RT reaches `f` only by lowering `e_eff` through a loss channel, never through the
bulk impulse. The same argument clears `eta_capture`: it is a momentum-conserving 2D/1D ratio
(ADR-0003), and the 2D track is inviscid Euler where grid-scale RT would be a numerical artifact
anyway.

**This is not a timescale dismissal — RT does grow.** The stagnation deceleration is
`g ~ v/τ_pulse ≈ 10⁴ / 10⁻⁴ = 10⁸ m/s²` (~10⁷ g). With `γ ≈ √(A·k·g)` (Atwood `A ~ 0.5`, mm-scale
`k ~ 10³ m⁻¹`) the growth rate is `~10⁵ s⁻¹`, so over a few-hundred-µs pulse RT runs ~10–20
e-foldings and goes **fully nonlinear within the pulse**. We do not claim RT is too slow; we claim it
is confined to a channel that does not carry the deliverable.

**The one channel it reaches: conductive-to-plate (channel 2, ADR-0016).** RT turbulent mixing brings
hot gas to the wall and roughens the boundary layer, enhancing wall heat transfer **above** the laminar
effusivity solve (ADR-0005). This is one-sided — RT can only *lower* `e_eff` — and, critically, the
enhanced conduction can exceed the floor's laminar value, so the conservative floor is **not**
automatically RT-safe on this channel.

**Gated trigger (the watchdog).** RT stays unmodeled *unless* the loss decomposition (ADR-0016) shows
conductive-to-plate is a non-trivial slice of `(1−f)` at some anchor. Where it is, that anchor gets a
**bounding turbulent-conduction estimate** — a mixing-enhancement multiplier on the laminar effusivity
flux — applied as a one-sided correction to the floor *before* `f` is quoted there. Conduction is
expected negligible at the high-v end (`τ≫1`, energy stays trapped in the gas) and at the low-v end
(condensation-dominated loss, weak RT drive in cold low-Atwood gas); the transitional anchor (ADR-0012)
is where to watch.

**Best-estimate caveat (ADR-0014).** RT bites hardest on the *ablating-wall* curve, not the rigid
floor: ADR-0014 credits a clean blowing/vapor-shield layer that raises `e_eff` by intercepting flux,
and RT mixing hot gas down through that shield is the mechanism that degrades it. The vapor-shield
recovery is therefore an **upper estimate** to the extent the shield stays RT-coherent.

## Consequence

RT (and its impulsive cousin Richtmyer–Meshkov) is not modeled. The conservative floor stands as a
lower bound on `f` **subject to the channel-2 watchdog**; the best-estimate ablating-wall vapor shield
is RT-optimistic; and RT concentrates peak local flux/ablation, so the ablation-per-pulse number
back-propagated to the MEMS replenishment system (§7) is an **underestimate** — RT loosens that
logistics requirement rather than tightening it, so it is conservative to omit.
