# Spray-plate vapor curtain: a sourced gray `kappa_vapor`

Research note for the argon arm of the spray plate ([ADR-0055](adr/0055-argon-arm-of-the-spray-plate.md),
[`argon_plate_handoff.md`](argon_plate_handoff.md) step 1b). The ablating wall (ADR-0014,
`crates/hydro1d/src/kernel.rs`, `Ablation::with_vapor_opacity`) sets the curtain's optical depth to
`tau = kappa_vapor * ablated_mass` (kg/m^2) and passes `1/(1+tau)` of the incoming flux to the steel.
The sweep uses 200, 1500 and 1e4 m^2/kg. This note asks what one gray number the curtain really has
over the plasma's 20-60 kK spectrum. No code was changed.

## Answer

- **Recommended `kappa_vapor` = 5e3 m^2/kg** in place of the unsourced 1e4.
- **Range 2e3-3e4 m^2/kg** for a cold curtain that still carries its soot. The low end is a 20 kK
  spectrum with graphite-sphere soot. The high end is a 60 kK spectrum with an aromatic vapor.
  Uncertainty is about a factor of 2 at fixed temperature, from the soot's UV absorption.
  Across 20-60 kK the value changes by a factor of about 5.
- **1e4 is about right for a 40 kK spectrum and too high below ~30 kK.** For argon's 31-61 kK
  starting gas it is defensible. For water's 14-30 kK it is about 2-3x too high (~4e3 fits better).
- **The Planck mean is much larger (2e4-7e4) and is the wrong number for this model.** The
  transmission `1/(1+tau)` at the film's 0.08-0.3 kg/m^2 is set by the most transparent part of
  the spectrum, not the average one.
- **The soot carries the curtain.** Without it, the oil vapor has a window below ~7.5 eV. That window
  holds 40-50% of a 20 kK Planck spectrum, and the effective gray value collapses to 5-100 m^2/kg.
- **Pessimistic bound: 10-100 m^2/kg** if the vapor is atomized (soot sublimed) but not ionized.
  Then everything below the carbon edge at 11.26 eV leaks through. See "Physical state".

## Data used

1. **Atomic photoabsorption (C, H, O).** NIST FFAST (Chantler 1995, 2000), from the site's own
   `FFAST_data.sqlite` on physics.nist.gov, column `mu_rho_pe`. Cross-checked against the Henke,
   Gullikson & Davis (1993) `.nff` files (f2) from henke.lbl.gov. Henke
   `mu/rho = 2 r_e lambda f2 N_A / A`.
2. **Molecular photoabsorption.** The Leiden photodissociation and photoionisation database
   (Heays, Bosman & van Dishoeck 2017), file set `all_cross_sections_text_continuum_0.1nm.zip`
   from home.strw.leidenuniv.nl. The calculation first used the copy in the VULCAN repository
   (Tsai et al. 2017). **The two agree exactly** for CH4, C2H6, C2H2, CO, H2O, C2H4, H2, C and H at
   7, 10, 15, 20, 30, 50 and 100 eV, both as point values and as 1% band averages. VULCAN adds its own
   long-wavelength tails (CH4 above 145 nm from PhiDrates; C2H6 above 150 nm from Au et al. 1993).
   Benzene is not in the Leiden set. Its VULCAN file (Boechat-Roberty et al. 2004; Capalbo et al.
   2016) is **not verified**.
3. **Soot UV absorption (5-11 eV), three models.** All are anchored or checked at 550 nm, where
   Bond & Bergstrom (2006) give 7.5 m^2/g.
   - **Graphite spheres.** Draine (2003b) graphite dielectric functions (`callindex.out_CpeD03_0.01`
     and `_CpaD03_0.01`). Below 30 eV these equal Draine & Lee (1984), built from laboratory data.
     Mie absorption for spheres of radius 5, 15 and 25 nm, density 1.8 g/cm^3, orientation averaged
     2/3 perpendicular + 1/3 parallel. This gives 3.3e3 m^2/kg at 550 nm, 2.3x below Bond &
     Bergstrom. It is used both as computed and scaled up by 2.3x.
   - **Measured flame soot.** Gavilan et al. (2016), SOLEIL synchrotron transmission of a ~50 nm
     ethylene-flame soot film down to 50 nm wavelength. Their Table 1 fit (a Lorentzian at
     4.86 eV, Gaussians at 6.7, 12.0 and 23.8 eV, plus a linear term) gives the spectral shape. It
     is scaled to 7.5 m^2/g at 550 nm. Widths are read as FWHM, and the linear term as
     `0.01 x` with zero intercept. Both readings are assumptions. The paper gives optical
     constants only as figures.
   - **The old power law** (`lambda^-1` from 550 nm) is kept for comparison.

   The measured shapes rise 4-13x from 550 nm to the 4.9-6 eV pi-band. They dip at 8-9 eV and rise
   again past 10 eV. The power law has no band, but its 5-11 eV average lands close to theirs.
   Above 11 eV soot is taken as atomic carbon (FFAST), except in the graphite runs, which use
   Draine's dielectric function throughout.

## 1. Cross sections and mass absorption

Conversion: `kappa = sigma * N_A / M`, with sigma in m^2 (1 Mb = 1e-22 m^2), `N_A = 6.022e23`/mol
and M in kg/mol. Example: CH4 at 15 eV, sigma = 46.3 Mb = 4.63e-21 m^2.
`4.63e-21 * 6.022e23 / 0.016043 = 1.74e5 m^2/kg`.

Leiden values, interpolated at the listed point, verified against the Leiden originals (C6H6 excepted).

| species (M, g/mol) | 7 eV | 10 eV | 15 eV | 20 eV | 30 eV | 50 eV | 100 eV |
|---|---|---|---|---|---|---|---|
| CH4 (16.04) sigma, Mb | ~0 | 18.2 | 46.3 | 33.2 | 14.3 | 3.5 | 0.52 |
| CH4 kappa, m^2/kg | ~0 | 6.8e4 | 1.7e5 | 1.2e5 | 5.4e4 | 1.3e4 | 1.9e3 |
| C2H6 (30.07) sigma | 0 | 24.1 | 68.0 | 52.5 | 20.8 | 5.5 | n/t |
| C2H6 kappa | 0 | 4.8e4 | 1.4e5 | 1.1e5 | 4.2e4 | 1.1e4 | n/t |
| C2H4 (28.05) sigma | 5.6 | 26.0 | n/t | n/t | n/t | n/t | n/t |
| C2H4 kappa | 1.2e4 | 5.6e4 | n/t | n/t | n/t | n/t | n/t |
| C2H2 (26.04) sigma | 1.0 | 6.1 (32) | 52.0 | 35.6 | 16.8 | 6.0 | 0.89 |
| C2H2 kappa | 2.4e3 | 1.4e4 | 1.2e5 | 8.2e4 | 3.9e4 | 1.4e4 | 2.1e3 |
| C6H6 (78.11) sigma, unverified | 213 | 34.3 | 106 | 119 | 61.1 | n/t | n/t |
| C6H6 kappa, unverified | 1.6e5 | 2.6e4 | 8.2e4 | 9.2e4 | 4.7e4 | n/t | n/t |
| CO (28.01) sigma | 0 | 0.05 | 27.4 (23) | 22.1 | 17.8 | 9.5 | 2.3 |
| CO kappa | 0 | 110 | 5.9e4 | 4.8e4 | 3.8e4 | 2.0e4 | 5.0e3 |
| H2O (18.02) sigma | 2.7 | 19.1 (7.6) | 17.1 | 22.6 | 15.6 | 7.1 | 1.8 |
| H2O kappa | 9.1e3 | 6.4e4 | 5.7e4 | 7.6e4 | 5.2e4 | 2.4e4 | 5.9e3 |
| H2 (2.016) sigma | 0 | 0 | 34.4 (11) | 7.6 | 2.7 | 0.53 | n/t |
| H2 kappa | 0 | 0 | 1.0e6 | 2.3e5 | 8.1e4 | 1.6e4 | n/t |
| C atom (12.01) sigma | ~0 | ~0 | 14.0 | 10.0 | 5.5 | 2.1 | 0.50 |
| C atom kappa | ~0 | ~0 | 7.0e4 | 5.0e4 | 2.8e4 | 1.1e4 | 2.5e3 |
| H atom (1.008) sigma | ~0 | ~0 | 4.9 | 2.2 | 0.71 | 0.16 | 0.019 |
| H atom kappa | ~0 | ~0 | 2.9e5 | 1.3e5 | 4.2e4 | 9.5e3 | 1.2e3 |

n/t = outside the tabulated range. In the spectral means, energies above a file's range use the
FFAST atomic sum for that formula, and energies below it are taken as transparent. A value in
brackets is the mean over +/-1% in wavelength, shown where it differs from the point value by more
than 10%. Those points sit on sharp bands. The spectral means use the full-resolution data, so this
does not affect them. CO lines below its ionization edge are not in the continuum file.

Atomic mass absorption, FFAST / Henke, m^2/kg:

| element | 10 eV | 15 eV | 20 eV | 30 eV | 50 eV | 100 eV |
|---|---|---|---|---|---|---|
| C | 6.4e4 / 2.8e4 | 3.7e4 / 6.7e4 | 2.7e4 / 6.5e4 | 1.7e4 / 3.1e4 | 8.1e3 / 9.9e3 | 2.3e3 / 2.4e3 |
| H | n/t / 0 | 2.3e5 / 0 | 1.1e5 / 1.5e5 | 3.6e4 / 4.2e4 | 8.6e3 / 9.5e3 | 1.1e3 / 1.2e3 |
| O | 7.2e4 / 1.9e4 | 5.4e4 / 1.5e4 | 4.2e4 / 5.0e4 | 3.0e4 / 4.1e4 | 1.7e4 / 2.3e4 | 5.2e3 / 6.8e3 |

The Henke files from henke.lbl.gov are byte-identical to the copy first used. The Henke page says
the grid runs from 10 eV, with f1 set to -9999 below 29 eV.

Graphite spheres (Draine 2003b, Mie, radius 15 nm), m^2/kg: 3.3e3 at 2.25 eV, 2.3e4 at 5 eV, 4.3e4
at 6 eV, 1.1e4 at 8 eV, 1.6e4 at 10 eV, 4.9e4 at 15 eV, 3.2e4 at 20 eV, 3.1e3 at 100 eV.

- **Above ~50 eV every C/H/O species agrees within ~30%:** 1e4 m^2/kg at 50 eV, 2e3 at 100 eV.
  Molecular binding does not matter there. Graphite agrees too.
- **At 10-30 eV everything absorbs at 1e4-2e5 m^2/kg.** The two atomic tables disagree there by up
  to 2x. Henke's title range starts at 50 eV, and FFAST puts carbon's edge at 6.4 eV against the
  measured 11.26 eV ionization energy. Neither is reliable below ~30 eV. The molecular data rule there.
- **Below ~10 eV the species split.** Alkanes and CO are near zero below ~8 eV. C2H2, C2H4 and
  especially benzene have strong bands at 6-8 eV. Atoms are transparent below their edges. Soot
  and graphite absorb at 1e4-5e4 per kg of carbon.

## 2. Gray means over 20, 40 and 60 kK

Planck-spectrum energy below a given photon energy (computed):

| T | < 6 eV | < 8 eV | < 11.26 eV (C edge) | < 13.6 eV (H edge) |
|---|---|---|---|---|
| 20 kK | 50% | 70% | 90% | 96% |
| 40 kK | 13% | 25% | 45% | 59% |
| 60 kK | 5% | 10% | 22% | 31% |

So at 20 kK most of the energy sits below the photoionization edges, in the band where cold vapor
is selective. "EUV" describes only the 40-60 kK cases well.

Three means, on a 0.3 eV-3 keV grid. `kappa_P` weights by B_nu. `kappa_R` is the Rosseland
harmonic mean weighted by dB_nu/dT. `kappa_T` is the gray value that reproduces the actual
Planck-weighted transmission `<1/(1+kappa_nu Sigma)>` at Sigma = 0.1 kg/m^2. That is the number
this model needs. It changes by under 15% between 0.01 and 0.3 kg/m^2 whenever soot is present.

Mixtures are per kg of film: 80% oil vapor + 20% soot, unless stated.

| curtain (soot model) | kappa_P 20/40/60 kK | kappa_R 20/40/60 | **kappa_T 20/40/60** |
|---|---|---|---|
| alkane + soot (Gavilan shape) | 1.8e4 / 5.7e4 / 6.1e4 | 5.3e3 / 1.7e4 / 2.4e4 | **3.8e3 / 1.0e4 / 1.8e4** |
| alkane + soot (graphite, 5-25 nm) | 1.7e4 / 5.6e4 / 6.0e4 | 3.4e3 / 1.5e4 / 2.3e4 | 1.9e3 / 7.5e3 / 1.6e4 |
| alkane + soot (graphite x2.3) | 2.3e4 / 6.6e4 / 6.9e4 | 7.7e3 / 2.7e4 / 3.2e4 | 4.5e3 / 1.5e4 / 2.6e4 |
| alkane + soot (old power law) | 1.8e4 / 5.8e4 / 6.1e4 | 5.3e3 / 1.9e4 / 2.5e4 | 3.5e3 / 1.1e4 / 1.9e4 |
| aromatic (C6H6) + soot (Gavilan) | 2.0e4 / 4.7e4 / 5.2e4 | 7.8e3 / 2.5e4 / 2.8e4 | 4.9e3 / 1.5e4 / 2.4e4 |
| aromatic + soot (graphite x2.3) | 2.6e4 / 5.6e4 / 6.0e4 | 8.8e3 / 3.0e4 / 3.4e4 | 4.8e3 / 1.7e4 / 2.9e4 |
| pyrolysis gas (C2H2 + H2) + soot (Gavilan) | 2.3e4 / 6.2e4 / 6.4e4 | 5.6e3 / 1.8e4 / 2.5e4 | 3.9e3 / 1.1e4 / 1.9e4 |
| atomic C + H, soot particles kept (Gavilan) | 3.7e4 / 7.1e4 / 6.1e4 | 3.4e3 / 8.5e3 / 1.4e4 | 2.9e3 / 5.6e3 / 9.8e3 |
| alkane, no soot | 1.8e4 / 6.5e4 / 7.0e4 | ~0 | 5 / 36 / 100 |
| atomic C + H, neutral, soot gone | 3.7e4 / 7.4e4 / 6.5e4 | ~0 | 2 / 15 / 43 |
| soot alone, per kg of soot (Gavilan) | n/c | n/c | 1.4e4 / 1.8e4 / 1.9e4 |

n/c = not computed.

Reading the table:

- **With soot, `kappa_T` is 2e3-5e3 at 20 kK, 7e3-1.7e4 at 40 kK and 1.6e4-2.9e4 at 60 kK.** The
  soot model moves it by about 2x. The vapor species move it by less than 1.5x.
- **The measured soot curves confirm the old power law in the mean.** Gavilan's soot shape gives
  within 10% of the power law at all three temperatures. Absolute graphite spheres give ~0.6x at
  20 kK and ~0.85x at 60 kK. The pi-band near 5-6 eV is real, but the 8-9 eV dip offsets it.
- **The existing 1500 (20% x 7.5 m^2/g) is the visible-light floor of the same soot.** The UV and EUV
  lift it 1.3-20x, depending on temperature and soot model.
- **Below-threshold treatment.** Under a molecular edge the tables give the measured (often near
  zero) value, with no fill. Under the atomic edges, only the tabulated lines and continuum are
  kept. With no soot, those windows decide `kappa_T`, and the Rosseland mean goes to zero. This is
  why the pessimistic bound is so low. Rosseland with soot (3e3-3.4e4) sits just above `kappa_T`.

## 3. Physical state of the curtain

Rough reasoning only.

- **Energy available.** Argon at 1 m puts about 1.3 GJ on the 78.5 m^2 face, about 1.7e7 J/m^2.
  The film boiled per pulse is 0.08-0.3 kg/m^2. If the curtain kept what it intercepts, that is
  ~6e7-2e8 J/kg.
- **What that buys.** Breaking all bonds in (CH2)n costs ~4e7 J/kg (~6 eV per CH2, estimate, not
  sourced). Subliming graphite costs ~6e7 J/kg (~7 eV/atom, estimate, not sourced). Singly ionizing
  carbon costs 9.0e7 J/kg (11.26 eV / 12 g/mol). Ionizing hydrogen costs 1.3e9 J/kg (13.6 eV /
  1 g/mol). So the curtain can be fully atomized, soot included, and its carbon partly ionized.
  Its hydrogen cannot be ionized.
- **Saha check for carbon** (`n_e n_+ / n_0 ~ 2.4e21 T^1.5 exp(-11.26 eV/kT)` m^-3, statistical
  weights ~1). At n ~ 1e25 m^-3, carbon is ~7% ionized at 1 eV (11.6 kK) and ~80% at 2 eV (23 kK).
- **Stratification is likely.** The kernel keeps the intercepted energy in the near-wall gas, so
  not all of it heats the vapor. The outer, plasma-side vapor takes the radiation first. It should
  atomize, lose its soot and partly ionize within the pulse. The inner vapor, near the 400 K wall
  and behind the outer layer, should stay molecular and sooty. That inner layer is what the
  `kappa_T` values above describe.
- **Ionized carbon shifts the spectrum.** C+ has no continuum from 11.26 eV up to its own 24.38 eV
  edge. A hot ionized layer therefore opens an 11-24 eV window. Free-free absorption, bound-free
  absorption from excited states, and line broadening at these densities partly close it.
  None of that is quantified here.
- **Net.** A cold, sooty inner layer gives 2e3-3e4. Losing the soot without ionizing heavily gives
  the 10-100 bound. A hot outer layer also re-radiates part of what it absorbs toward the wall at
  its own lower temperature, which a single gray `1/(1+tau)` cannot represent.

## 4. Recommendation

| use | kappa_vapor (m^2/kg) | basis |
|---|---|---|
| replace the unsourced 1e4 | **5e3** | `kappa_T`, cold sooty curtain, 20-40 kK, middle of the soot models |
| sweep low (cold spectrum, weak soot) | 2e3 | 20 kK, graphite-sphere soot |
| sweep high (hot spectrum) | 3e4 | 60 kK, aromatic vapor or strong soot |
| pessimistic, soot destroyed | 30-100 | atomic neutral or soot-free vapor at 40-60 kK |
| keep as visible-light floor | 1.5e3 | Bond & Bergstrom, 550 nm |

- **Uncertainty:** factor of ~2 at fixed T from the soot's UV absorption, and a factor of ~5
  across 20-60 kK. A soot-free curtain would sit far below the range.
- **1e4 verdict:** about right for argon's hot gas (7e3-1.7e4 at 40 kK), but too high for 20-30 kK
  spectra and for water. Its stated basis, "EUV absorption", is misleading. EUV values are
  3e4-2e5, but they do not set the transmission. The soot's near-UV window does.
- **200 verdict:** below every cold-curtain value here, and above only the soot-free case. It is a
  hedge for partial soot loss, not a measurement.
- A temperature-dependent gray value, about `kappa_vapor ~ 3.8e3 * (T / 20 kK)^1.4` (fit to the
  Gavilan alkane row, 20-60 kK), would follow the table better than any constant.

## Limitations

- Room-temperature molecular cross sections. Hot molecules have broader, red-shifted bands, which
  would partly fill the 6-8 eV window. That is not included.
- Oil vapor is stood in for by C2H6 or C6H6. Real pyrolysis gas is a mixture. The table shows the
  choice matters < 1.5x while soot is present. The C6H6 data are not verified.
- Soot UV absorption is still the largest uncertainty. Graphite spheres sit 2.3x below Bond &
  Bergstrom at 550 nm. Gavilan's shape is read from fit parameters whose width convention the
  paper does not state. Draine (2003b) calls graphite's E-parallel-c function in the UV uncertain.
- Soot as a per-mass absorber assumes particles smaller than an absorption length. Mie runs at 5,
  15 and 25 nm radius change `kappa_T` by under 15%. Large aggregates would self-shield in the EUV.
- Planck weighting assumes the face sees a blackbody at one temperature. The kernel's gray
  flux-limited radiation does not carry a spectrum.
- Section 3 is order-of-magnitude only. The bond and sublimation energies there are not sourced.
- The Bond & Bergstrom 7.5 m^2/g comes from the repo's existing usage. The paper itself was not
  re-read (see below). The 1.8 g/cm^3 soot density in the Mie runs is an assumption, not sourced.

## Sources

- Chantler, C. T. (1995). J. Phys. Chem. Ref. Data 24, 71. doi:10.1063/1.555974. Chantler, C. T.
  (2000). J. Phys. Chem. Ref. Data 29, 597. doi:10.1063/1.1321055. Data:
  https://physics.nist.gov/PhysRefData/FFast/FFAST_data.sqlite (downloaded 2026-10-05).
- Henke, B. L., Gullikson, E. M. & Davis, J. C. (1993). X-ray interactions: photoabsorption,
  scattering, transmission, and reflection at E = 50-30,000 eV, Z = 1-92. At. Data Nucl. Data
  Tables 54, 181-342. doi:10.1006/adnd.1993.1013. Data:
  https://henke.lbl.gov/optical_constants/sf/c.nff (and h.nff, o.nff).
- Heays, A. N., Bosman, A. D. & van Dishoeck, E. F. (2017). A&A 602, A105.
  doi:10.1051/0004-6361/201628742. Data: https://home.strw.leidenuniv.nl/~ewine/photo/ ,
  `data/photo_data/all_cross_sections_text_continuum_0.1nm.zip`. The file headers cite the
  underlying measurements, for example CH4: Au et al. (1993), doi:10.1016/0301-0104(93)80142-V,
  and Kameta et al. (2002); CO: Chan et al. (1993); H2O: Chan et al. (1993) and others.
- Tsai, S.-M. et al. (2017). VULCAN. ApJS 228, 20. doi:10.3847/1538-4365/228/2/20 (verified).
  Copy: https://github.com/exoclime/VULCAN/tree/master/thermo/photo_cross (C6H6 file only).
- Draine, B. T. (2003b). Scattering by interstellar dust grains. II. X-rays. ApJ 598, 1026.
  doi:10.1086/379123. Open copy: https://arxiv.org/abs/astro-ph/0308251. Data:
  https://www.astro.princeton.edu/~draine/dust/dust.diel.html. Draine, B. T. & Lee, H. M.
  (1984). ApJ 285, 89. doi:10.1086/162480. Taft, E. A. & Philipp, H. R. (1965), Phys. Rev. 138,
  A197, doi:10.1103/PhysRev.138.A197, is one of DL84's graphite inputs. That link is **not
  verified** here; DL84 was not re-read.
- Gavilan, L. et al. (2016). VUV spectroscopy of carbon dust analogs. A&A 586, A106.
  doi:10.1051/0004-6361/201527098. Open copy: https://arxiv.org/abs/1601.08037 (Table 1 used).
- Bond, T. C. & Bergstrom, R. W. (2006). Aerosol Sci. Technol. 40, 27-67.
  doi:10.1080/02786820500421521 (DOI verified; text not reachable).
- NIST Atomic Spectra Database, ionization energies (verified 2026-10-05):
  H I 13.598434599702 eV, C I 11.2602880 eV, C II 24.383143 eV.
  https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
- Repo: `docs/argon_plate_handoff.md` (face energy, film masses, 78.5 m^2, gas temperatures);
  `crates/hydro1d/src/kernel.rs` (curtain model). Means computed by a scratch script, not checked in.

## Blocked sources

- www.tandfonline.com (Bond & Bergstrom 2006): Cloudflare "Just a moment" challenge. No open copy
  found (OpenAlex lists it as closed). The Wayback Machine returned HTTP 429 (rate limit).
- www.astro.uni-jena.de (Jena optical-constants database, for Zubko 1996 and Jäger 1998
  amorphous carbon): redirects to a university landing page. The archived copy returned HTTP 429.
- pubs.aip.org and www.sciencedirect.com: publisher bot walls. Not needed after the open copies.
- Berkowitz (2002) and Palik are books and were not consulted.
