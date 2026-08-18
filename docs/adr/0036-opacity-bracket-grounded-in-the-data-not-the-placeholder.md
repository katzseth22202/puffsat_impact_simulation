# The opacity bracket expresses what is uncertain about the real table, not ignorance of a placeholder

The sweep inherited a **0.1×–10×** opacity bracket. That range was sized when the table carried a
*placeholder* Kramers opacity that CONCLUSION.md records as ~2000× low: it expressed ignorance about
a stand-in, not uncertainty in a measurement. The table now carries real TOPS/OPLIB gray means, so
the bracket should express what is actually uncertain about *those*.

**Result: the band on `f` narrows 12.6×**, from `Δf = 0.069` to `Δf ≤ 0.0055`. It is no longer the
dominant uncertainty; the freeze-timing bracket (~0.05) is.

## Three inputs, each measured or cited rather than assumed

**1 — Published accuracy of the data.** Farag et al. 2024 (ApJ 968:16, arXiv:2406.02845) §2.1, citing
Huebner & Barfield (2014), gives the OPLIB Rosseland-mean uncertainty by dominant process: ~5%
electron scattering, ~10% free-free, ~20% bound-free, ~30% bound-bound. The uncertainty falls as
ionic charge rises toward the hydrogenic limit, which is why hot, highly-ionized states here take the
tighter number. Mean charge `Z̄` is the sharper proxy for the classification in this mixture.

**2 — The non-LTE contribution, measured and negligible.** Both the Saha EOS and the LTE TOPS opacity
assume collisionally controlled level populations, so that assumption was checked before either was
trusted (below). Propagating even a 100% error in the one failing stage's population gives
`dZ̄/Z̄ = 0.71%`, and with the table's measured `dln κ/dln Z̄` of order unity that is `dκ/κ ≤ 0.93%` —
an order below the data's own accuracy, and identically zero at every other probed state. **The
bracket is set by the opacity data, not by the LTE question.**

**3 — Sensitivity of the answer, measured from the sweep itself.** `de_eff/dln κ` is read off the
sweep's own opacity-scale rows rather than assumed, spanning 0.0009–0.071 across the probed states:
growing steeply with velocity and cloud length, falling steeply with impact density.

The band is then `Δf = eta_capture · (de_eff/dln κ) · ln(1 + fraction) / 2`.

## Explicitly not the Rosseland/Planck ratio

That ratio measures **spectral structure**, not model error. This matters because the Sₙ transport
check (ADR-0037) finds the two single-mean tallies disagreeing by up to 164% — more than either
disagrees with FLD. That spread is the gray model's inability to represent a spectrum (ADR-0006), and
reporting it as a band on `f` would convert a known modelling limitation into a fictitious
measurement uncertainty. Only multigroup resolves it.

## LTE validity, checked rather than assumed

Griem/McWhirter at the probe turnaround states: `n_e ≥ 1.6e12 √T ΔE³ cm⁻³`. The primary reported
quantity is the criterion **inverted** — the critical energy gap at the state's own `(n_e, T)` — so
the headline needs no atomic data and does not rest on line-list accuracy.

| band | governs | margin | verdict |
|---|---|---|---|
| 16 km/s | O II | 2.5× – 70× | PASS |
| 22 km/s | O II | 7.4× – 206× | PASS |
| 28 km/s | O II | 9.3× – 295× | PASS |
| 45 km/s | O IV | 7.5× – 147× | PASS |
| 55 km/s | O V | 3.6× – 71× | PASS |
| 63 km/s | O V | 3.9× – 74× | PASS |
| 69 km/s, ρ ≥ 0.02 | O V | 7.8× – 49× | PASS |
| 69 km/s, ρ = 0.01 | **O VII** (574 eV, ~1%) | **1.7e-4** | **FAIL** |

**The tightest margin is at the coldest end, not the hottest** — `n_e` falls faster than `ΔE` on the
way down. That inverts the intuition the range extension was worried about: if any part of 16–63 is
LTE-marginal it is the bottom, which the *existing* study already published. The single failure is a
trace helium-like stage at 69 km/s, outside the extended range.

45–63 km/s is **checked directly**, not inferred by bracketing between 28 and 69. Bracketing would
have assumed a monotonicity the data had already contradicted.

**Caveat carried forward:** McWhirter is **necessary, not sufficient**, and is derived for
homogeneous stationary plasmas. A passing margin means the criterion does not rule LTE out; it is not
a proof of LTE in a transient expanding flow.

## The atomic data is NIST-verified, and that mattered

Three documented selection traps, each of which yields a plausible wrong number:

1. **A wavelength window is not a selector.** O II's 834 Å multiplet has its lower level ~30 eV above
   ground (4p–5s) — it is not a resonance line.
2. **`unit=1` is nanometres, not Ångström.**
3. **The longest ground-connected allowed line is usually a weak *intercombination* line** (O V
   ¹S–³P° at 121.8 nm, `gA` 7e3 against the resonance line's 8.6e9) at a *smaller* gap, which would
   flatter the criterion. Selection is therefore by longest wavelength whose `gA` is within 1e-3 of
   the strongest ground-connected line.

This caught **O III at 17.64 eV where the true resonance gap is 14.84 eV** (702.90 Å was not the
resonance line; 835.29 Å is). The error had wrongly made O III govern at 28 km/s and understated that
margin 1.7×. Pinned by a regression test.
