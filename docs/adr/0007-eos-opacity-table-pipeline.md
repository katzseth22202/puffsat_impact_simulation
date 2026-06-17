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

## Considered Options

- **Single-source EOS or opacity across the whole range.** Rejected: none exists — CoolProp lacks
  ionization, CEA lacks condensation, TOPS/OPLIB lack molecular bands.
- **Citing TOPS/OPLIB alone** (as the original §9 draft did). Rejected: covers only the plasma half,
  leaving the low-v optical fraction — which gates the seed study — without a citable source.
