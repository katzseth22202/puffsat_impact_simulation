# EOS/opacity tables: a regime-stitched pipeline with split, per-regime provenance

No single source spans water from ~300 K to ~43 kK across the EOS and opacity needs (two-phase
condensation, dissociation/ionization, molecular bands, plasma opacity). The Python cold path
therefore generates **one unified table** `(ρ, T) → (p, e, c_s, κ_Rosseland, κ_Planck)` by stitching
regime-specific sources, with provenance cited **per regime**.

**EOS seam:**
- **CoolProp** — low-T real fluid, the saturation curve and two-phase region (the condensation
  physics of ADR-0004), valid to ~1000–1300 K.
- **CEA / Saha** — high-T equilibrium dissociation + ionization (the §4 specific-heat buffering);
  ideal-gas-mixture based, no condensation.
- Seam ~1–3 kK, enforced **C¹-continuous in `p(ρ,T)` and sound speed** so the hydro does not ring
  at the join.

**Opacity seam — the provenance gap this ADR closes.** TOPS and OPLIB are *plasma* opacity codes;
they do **not** cover cool molecular water vapor — the rovibrational bands that set the Planck mean
at 3.2 km/s (ADR-0006). So opacity provenance splits the same way as the EOS:
- **HITEMP / ExoMol** — hot H₂O vapor bands (low-v Rosseland + Planck means).
- **soot Mie/Rayleigh model** — the dark-oil seed (§7).
- **TOPS or OPLIB** — the ionized high-v end.

Each regime's source is cited, so the low-v optical fraction is no longer un-sourced.

**Amendment (ADR-0019): serialization format.** The unified table is serialized as **JSON**,
not the HDF5 originally implied by §8 — it is small and loaded whole into RAM, so a binary
container is unjustified. The per-regime provenance this ADR mandates is carried as a nested
`provenance` object inside the table JSON. See ADR-0019.

**Amendment (Rung B, high-v pass): what landed, and what is deferred.** The 16 km/s `e_eff(ρ)`
deliverable needs only the high-T half of this pipeline, so Rung B implemented:
- **EOS — a self-contained analytic equilibrium model** (`python/puffsat/eos_water.py`): H₂O⇌2H+O
  by mass action + H/O Saha ionization, closed by element conservation + charge neutrality
  (Zel'dovich & Raizer Ch. III). This *replaces*, for the high-v pass, the CoolProp+CEA seam above —
  CoolProp only reaches ~1300 K, and the low-T condensation seam is not exercised at 16 km/s. The
  CoolProp two-phase seam (ADR-0004) is deferred to the low-v (3.2 km/s) anchor that needs it.
- **Opacity — an explicitly PROVISIONAL Kramers-shaped bracket** (`python/puffsat/tables.py`):
  `κ_R = κ₀(ρ/ρ_ref)(T/T_ref)^-3.5`, `κ_P = 3κ_R`, calibrated so `τ = ρκ_R L` sits mid-band of the
  design's `[10², 10⁵]` at the nominal stagnation. The **real per-regime opacity** (HITEMP/ExoMol +
  soot + TOPS/OPLIB) is deferred to its own *survivability-flux* rung.

  This deferral is **empirically licensed**, not assumed. The B5d-3 scan re-ran the sweep at
  **0.1× / 1× / 10× opacity** (`python/puffsat/sensitivity.py`, `make sensitivity`): `e_eff` moved
  by at most **1.0e-2 (1.6% relative)** across that 100× range, while the radiative loss channels
  (1a/1b) *did* move with κ — exactly design §3 (`e_eff` is EOS/gas-dynamics-dominated and
  opacity-insensitive at `τ≫1`; the opacity is load-bearing only for the survivability flux). The
  loader's JSON boundary lets the real opacity table hot-swap in later with no kernel change.

**Amendment (Rung C / B-flux): two optional low-v fields.** The low-v (3.2 km/s) CoolProp table
(`build_table_lowv`) carries two fields beyond the five standard ones, both optional and `#[serde(
default)]` so the high-v table and the loader stay backward-compatible:
- **`liquid_frac` ∈ [0,1]** — the condensed mass fraction for the wall-sticking sink (ADR-0004,
  channel 3). Interpolated **linearly** (it is legitimately `0`, so it does *not* ride the positive
  log-interp path the other fields use).
- **`k_gas` > 0** — the gas thermal conductivity for the B-flux gas-side conduction operator
  (ADR-0005), from CoolProp/IAPWS transport (`PropsSI("conductivity", …)`; the two-phase dome, where
  transport is undefined, uses the saturated-vapor value). Strictly positive, so it rides the same
  **log-interp** path as the opacities. A `k_gas_scale` knob mirrors `kappa_scale` for a sensitivity
  scan. The **high-v plasma `k_gas`** (Spitzer-like transport) is *not* tabulated — the high-v table
  has no `k_gas`, so high-v conduction stays off (the deferred B-flux high-v sibling).

## Considered Options

- **Single-source EOS or opacity across the whole range.** Rejected: none exists — CoolProp lacks
  ionization, CEA lacks condensation, TOPS/OPLIB lack molecular bands.
- **Citing TOPS/OPLIB alone** (as the original §9 draft did). Rejected: covers only the plasma half,
  leaving the low-v optical fraction — which gates the seed study — without a citable source.
