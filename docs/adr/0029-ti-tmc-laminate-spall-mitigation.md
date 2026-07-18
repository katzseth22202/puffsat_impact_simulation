# A Ti / titanium-matrix-composite laminate spall mitigation is deferred: it scatters the pulse but hardens a non-binding layer, and multi-pulse durability is the gate

A reviewer question in the spall family (with ADR-0011 and its Ti back-face amendment): *"intersperse
the solid titanium backing with layers of titanium-matrix composite (TMC), so the composite's
heterogeneous microstructure breaks up the high-frequency content of the compression pulse and reduces
spall."* Paired with the honest worry that motivated it: *would the composite itself crack under the
reflected tension and lose its effectiveness?* This ADR records the evaluation so the idea is logged as
considered — the mechanism is sound, but for this stack the leverage is at the wrong interface and the
make-or-break is a durability question the repo already defers.

**The scattering mechanism is real — and a low-fidelity version is already in the design.** Acoustic
scattering off impedance heterogeneities is strongly frequency-dependent (Rayleigh-like, rising steeply
with frequency while the scatterers are small vs the wavelength), so a heterogeneous TMC preferentially
scatters and attenuates the *high-frequency* content of a shock front. The result is a longer rise time
and a lower peak. Because spall is governed by the *peak* tensile stress at the reflection plane, a
dispersed pulse spalls less — this is exactly the principle behind graded/laminated "spall liners,"
functionally-graded impedance stacks, and architected mitigators. The design already exploits a coarse
form of it: ADR-0011's graded interlayer / compliant braze "doubles as impedance smoothing," and design
§5 / ADR-0011 name pulse rise-time (the stretched high-v cloud softening the stress gradient) as a spall
lever. So the proposal is a higher-fidelity extension of a lever already in play, not a new principle.

**Cracking is a *feature* on a single pulse and a *risk* over many.** The worry splits by regime. On one
pulse, distributed microcracking — fiber/matrix debonding, particle decohesion, crack deflection at
interfaces — is precisely the composite toughening mechanism: it converts one catastrophic
through-thickness scab into many small, benign, energy-absorbing debonds. A composite designed for this
"spalls" gracefully and distributed, which is the point of a sacrificial heterogeneous layer at the
tension plane. Over the mission's *many* pulses, though, the question is whether that microdamage
**shakes down** (forms once, then arrests, stable thereafter) or **ratchets** (accumulates
pulse-over-pulse to loss of function) — set by whether the per-pulse tension sits below the composite's
fatigue/shakedown threshold. A composite that ratchets loses both its scattering and its
impedance-termination, and a progressively delaminated layer is a growing *near-free* surface, which
ADR-0011's corollary shows drives the reflection coefficient back toward `R = −1` (prompt spall). So the
mitigation's own failure mode is real, it is not analyzable in closed form, and it is the same
cyclic-fatigue-over-pulse-count axis the design parks for Phase-2 (design §11; the Vectran back-face
carries the identical caveat).

**The binding tensile mode is the brittle SiC, not the titanium — so TMC *behind the Ti* hardens the
wrong layer.** From the two spall checks now operationalized (ADR-0011 + its amendment, in
`analysis.classify_survivability`):

| Tensile mode | Stress at the plane | Dynamic spall strength | Governs? |
|---|---|---|---|
| SiC–Ti interface (brittle SiC) | `0.15·peak` | `0.3–1.0 GPa` | **yes** — reaches its limit at the foreclosed ~2 GPa f-max corner |
| Ti back-face free surface (ductile Ti) | `0.85·peak` | `2.5–4.5 GPa` | no — never binds (peaks at ~335 MPa on the heavy plate) |

The titanium sees ~5.7× more tension yet has ~8× the strength, so the tension/strength ratio at the SiC
interface (`0.15/0.3 = 0.50 GPa⁻¹`) exceeds the Ti back face's (`0.85/2.5 = 0.34 GPa⁻¹`) at *every*
load: the SiC spalls first, always. Interspersing scattering layers *behind* the Ti therefore stiffens a
layer that already carries a factor of ~7 of margin — it solves a problem the stack does not have. The
tensile risk that actually binds lives between the gas-facing SiC and the Ti, which is exactly where the
graded interlayer / compliant braze already sits, and where the design table already flags the C/SiC
composite hot-face's "fibers vulnerable if matrix cracks" — the same worry, already noted for the layer
that matters.

**If pursued, the disciplined form is particulate TMC at the SiC→Ti transition — not fiber, not behind
the Ti.** Two material constraints fix the design. (1) **Particulate, not continuous-fiber.** Fiber TMCs
(e.g. SCS-6 SiC/Ti) are weak in *transverse / through-thickness* tension — and spall tension is exactly
through-thickness, so in-plane fibers build a designed cleavage plane in the worst possible orientation
(plus brittle SiC/Ti reaction layers). A particulate TMC (TiB / TiC-reinforced Ti) scatters waves
isotropically with no preferential weak plane. (2) **Thermal placement.** TMC interfaces and any SiC
reinforcement degrade at temperature and the SiC/Ti reaction layer is thermally unstable, so the layer
must stay cool (below the Ti's ~700 K ceiling, away from the hot face) — compatible only near the
SiC→Ti transition once that region is under temperature, which is already a design constraint on the Ti.

**What the current kernels can and cannot say.** The 1D and 2D kernels model the stack as *impedance
layers only* — no microstructure, no fracture, no fiber/matrix — so they cannot evaluate cracking at
all. The peak/rise-time *benefit* of a laminated or graded impedance column is modelable and in scope: a
1D layered stack with reflection/transmission at each interface would put a number on how much the
interspersing lowers the peak tension. The cracking half — shakedown vs ratcheting, interfacial fatigue
— is irreducibly an FEA-with-damage-model plus coupon-test question, out of scope here (two hydro
kernels, no structural/fracture solver; design §5/§11, ADR-0011, ADR-0027).

## Consequence

The interspersed Ti/TMC laminate is **not adopted now and not relitigated**; it is deferred to the
Phase-2 structural/durability track. The current stack stands unchanged — ADR-0011's graded interlayer
plus a continuous, solid, void-free, impedance-terminating Ti layer — and the explicit Ti back-face
spall check (ADR-0011 amendment) confirms the titanium carries ample margin, so *no* spall mitigation
behind the Ti is warranted. Should a scattering/spall-liner layer ever be added, this ADR fixes its
shape: a **particulate** TMC at the **SiC→Ti** transition (not fiber, not behind the Ti), fed by the
solver's peak-load / rise-time / footprint outputs, not changing `f(v)` (a parallel Phase-2 task like
the truss/Vectran body of ADR-0027), and **gated by the multi-pulse shakedown-vs-ratcheting question**,
which is the same deferred cyclic-fatigue axis (§11) — the make-or-break, not the single-pulse behavior.
If a quantitative estimate is wanted before Phase-2, the in-scope step is a 1D layered/graded impedance
column to bound the peak-tension reduction; the damage question waits for FEA + coupons.

## Considered Options

- **Intersperse solid Ti with TMC layers behind the Ti back face (as proposed).** Deferred / not
  adopted: it hardens a layer that already has ~7× margin, because the binding tensile mode is the
  brittle SiC interface, not the ductile Ti back face (ADR-0011 amendment). The scattering mechanism is
  sound but aimed at the wrong plane.
- **Continuous-fiber (SiC/Ti) TMC at the tension plane.** Rejected on materials grounds: fiber TMCs are
  weak in transverse / through-thickness tension — the spall direction — so in-plane fibers introduce a
  designed weak plane in the worst orientation, plus brittle interfacial reaction layers.
- **Particulate (TiB/TiC-Ti) TMC scattering layer at the SiC→Ti transition.** Noted as the disciplined
  form *if* a spall liner is ever pursued: isotropic scattering, no weak plane, thermally placeable — but
  gated by multi-pulse durability (shakedown vs ratcheting), the deferred fatigue axis, and unnecessary
  while the Ti back-face check passes with margin.
- **Model the idea in the current kernels.** Partially in scope: the impedance / rise-time peak-reduction
  is modelable as a 1D layered column; the cracking / damage accumulation is not (needs an FEA damage
  model and coupon tests), so a full evaluation cannot be produced from the hydro kernels alone.
