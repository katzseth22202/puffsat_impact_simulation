# Shock reflection at the SiC–Ti interface: tensile spall is a Phase-2 structural check, decoupled from the shock absorber

The compressive shock that loads the facesheet partially reflects at the SiC–Ti bond as a *tensile*
wave back into the brittle SiC, because Ti has lower acoustic impedance than SiC.

Impedances `Z = ρc`: SiC ~`3.6×10⁷`, Ti ~`2.75×10⁷` Pa·s/m, so
`R = (Z_Ti − Z_SiC)/(Z_Ti + Z_SiC) ≈ −0.15` — about 15% of the incident compressive amplitude
returns as tension, the direction ceramics are weak in.

**The Ti backing is protective, not the cause.** A bare SiC layer with a *free* back surface would
reflect the pulse as *full* tension (`R = −1`) and spall promptly; the Ti backing converts that to
`R ≈ −0.15`. This is the quantitative mechanism behind §5's "the ductile backing keeps the brittle
ceramic from spalling catastrophically." A higher-impedance backing would give `R ≥ 0` (compression
reflects as compression, zero tension) at a mass cost; the graded interlayer / compliant braze already
required for the SiC–Ti CTE mismatch doubles as impedance smoothing that further cuts the reflection.

**Sub-dominant at baseline loads.** At the 400 MPa baseline limit the reflected tension is ~60 MPa,
below SiC dynamic spall strength (~0.3–1 GPa), so the compressive limit binds first; incident stress
would need ~2 GPa for the reflected tension to reach spall. Cyclic fatigue could erode this over many
pulses, but multi-pulse accumulation is out of scope (§11), so single-pulse spall is the check. The
survivability frontier therefore tests **both** the compressive facesheet limit and the
reflected-tensile (spall) limit at the interface.

**The shock absorber cannot and need not mitigate it.** Two problems ~10³× apart in time: the absorber
acts on the plate's *bulk recoil* (~230 m/s/pulse) on the vehicle/cadence scale (ms–s), while the spall
reflection is an *internal stress wave* in the stack on the µs scale (SiC transit ~0.1 µs/mm, Ti ~µs).
The wave reverberates and decays long before the absorber responds. Spall mitigation lives in the stack
— impedance grading, pulse rise-time (the stretched high-v cloud softens the stress gradient and peak
tension), Ti toughness — not the absorber.

**Scope.** The tensile/spall analysis is a Phase-2 structural task (with the SiC–Ti bond, §5), fed by
the 1D solver's peak-load history. The sim supplies the load; the spall check is structural analysis
on top.

## Amendment (2026-06, Rung S): the reflected-tensile check, operationalized

The survivability frontier (design §7) now tests both limits. The reflected tensile stress at the
SiC–Ti interface is taken as `|R|·peak ≈ 0.15·peak` (this ADR's impedance result), against SiC dynamic
spall strength 0.3–1 GPa. With the binding compressive peaks landing at the 400 MPa baseline, the
reflected tension is ~60 MPa — far below spall — so **the compressive facesheet limit binds first**,
exactly as this ADR's "sub-dominant at baseline" analysis predicted, now confirmed numerically across
the swept shapes. The reflected-tensile check is carried in `classify_survivability` as the second
gate but never controls the frontier at these loads; it would only bind if incident stress reached
~2 GPa (which the foreclosed `f`-max corner does reach, but that corner already fails on compression).

## Amendment (2026-07, Rung S): the Ti back-face free-surface spall check, added explicitly

The interface check above covers the *brittle SiC*. It does not cover the *Ti layer itself*, whose
back face is a second spall site: the compression **transmitted** through the SiC–Ti step travels
through the solid Ti and reflects off its back surface as tension. That back surface is free (or, per
the corollary below, terminated by a low-impedance truss / tensioned-fiber back-face — a *near*-free
surface), so `R ≈ −1` and the full transmitted amplitude returns as tension.

Using the same impedance step, the transmitted stress fraction is `T = 1 + R = 1 − |R| ≈ 0.85` (stress
convention), so the peak Ti back-face tension is `≈ 0.85·peak` — much larger than the SiC interface's
`0.15·peak`. But Ti is **ductile**, with dynamic spall strength ~2.5–4.5 GPa (Ti-6Al-4V / CP-Ti),
~8× the SiC's 0.3–1 GPa. Comparing tension-to-strength ratios, the SiC interface is closer to failure
at every load — `0.15/0.3 GPa⁻¹ > 0.85/2.5 GPa⁻¹` — so **the brittle SiC interface spalls first and the
Ti back-face check never controls** (it would only bind if Ti spall strength fell below ~1.7 GPa, far
under any titanium). At the foreclosed `f`-max corner (~2 GPa incident) the SiC interface reaches its
0.3 GPa limit while the Ti back-face tension (~1.7 GPa) is still under Ti spall.

**Operationalized.** `classify_survivability` now carries a third gate (`survives_back_spall`, via
`back_face_tensile = 0.85·peak` against `TI_SPALL_LO = 2.5 GPa`), ANDed into the baseline/relaxed
survivability verdicts across the Rung S, Jupiter, and heavy-plate frontiers; `structure.py` reports
it as sub-check (3b) alongside the SiC interface (3a). It is confirmatory — a defensible, explicit
guard that the ductile Ti backing does not itself spall — and passes at every swept load with margin,
consistent with this ADR's "the ductile backing keeps the ceramic from spalling" narrative now applied
to the backing's own back face. Conservatism is doubled in: full `R = −1` free-surface reflection
(a solid-terminated back face reflects less) and the low end of the Ti spall band. The same
multi-pulse-fatigue caveat (§11, out of scope) applies to the Ti back face as to the SiC interface.

## Considered Options

- **Delegate spall mitigation to the shock absorber.** Rejected: timescale mismatch ~10³× — the
  absorber handles bulk recoil (ms–s), not the µs internal stress wave.
- **Test only the compressive limit at the interface.** Rejected: the impedance mismatch reflects
  tension into the brittle layer, so the frontier must test the reflected-tensile/spall limit too
  (even though it is sub-dominant at baseline).
- **Higher-impedance (e.g. tungsten) backing to zero out the tensile reflection.** Noted, not chosen:
  mass penalty; the graded interlayer plus Ti toughness keep the baseline reflection well below spall.

## Corollary — plate structural architecture

The reflection result constrains the plate's structure: **no voids may sit directly behind the SiC**
(truss/corrugated core, membrane gaps) — a void is a near-free surface that drives `R` back toward
`−1` (prompt spall). A continuous, solid, impedance-terminating Ti layer is therefore required
immediately behind the SiC. Any lightweight structural body — a Ti truss/corrugated core with a
tensioned high-strength-fiber (e.g. Vectran) tension back-face, a Phase-2 mass-optimization far more
efficient in bending than solid Ti — may begin only *behind* that solid layer. The structural choice
does not change `f(v)` provided the body's first vibrational-mode period stays ≫ the pulse duration
(rigid-wall). Vectran is thermally limited (polymer; keep cool, far from the hot face) and needs a
cyclic-fatigue check over the mission pulse count.

A higher-fidelity variant of this — interspersing the solid Ti with titanium-matrix-composite layers
to scatter the pulse's high-frequency content and mitigate spall — is evaluated and **deferred** in
ADR-0029: the scattering mechanism is sound, but the binding tension is at the SiC interface (not the
Ti back face, which the amendment above shows has ample margin), and the make-or-break is multi-pulse
durability, the same deferred cyclic-fatigue axis.
