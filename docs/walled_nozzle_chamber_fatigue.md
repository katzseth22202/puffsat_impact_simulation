# Walled nozzle: fatigue, spall and overwrap delamination of the autofrettaged chamber wall

Status: analysis note, 2026-10-06. Checks the parent's layered, autofrettaged wall
(`sec:carbon_overwrap`, `tab:layered_wall_mass`) for the 20 m^3 near-term chambers of the 2.5 kg
rod: hydrogen at 5500 K (818 bar) and methane at 7000 K (496 bar).

Reproduce:
- `make walled-nozzle-chamber-fatigue` (`chamber_fatigue.py`): crack growth under the breathing
  mode.
- `make walled-nozzle-wall-waves` (`wall_waves.py`): through-thickness spall and delamination.

The required life is 630-780 pulses per departure (parent `sec:methane_7000_near_term`), and
about 3,900 if the chamber is recovered and flown five times. The growth ledger expends it, so
one departure is the requirement there.

## 1. Bottom line

1. **As designed, the hydrogen chamber's steel does not last one departure.** The parent's own
   sizing has the steel swing -0.75 Y to +0.75 Y with an overshoot of D = 1.7. Grow a 1 mm crack
   on the ASME CC2938 hydrogen curve and the wall fails in **82-482 pulses**. The breathing mode
   rings after each pulse, 7-74 tensile cycles per pulse, and each ring grows the crack. At
   ~600 MPa peak in hydrogen, the critical crack is only ~3.5 mm deep.
2. **The methane chamber mostly survives one departure.** On the air curve, with the parent's
   D = 1.7, it lasts 718-4,565 pulses. With a rise of one acoustic crossing (D 1.3) it lasts
   10,800-39,000. If hydrogen from dissociation does embrittle its wall, it drops to the
   hydrogen chamber's numbers.
3. **The fix is a rise time, as on the spray plate.** If the chamber pressure rises over at
   least one breathing period of the wall (`T_b` ≈ 1.1 ms), D falls to 1.0-1.2. The steel then
   peaks at 100-234 MPa instead of ~600, and **both chambers last 8,000 to >10^6 pulses, even on
   the hydrogen curve.** The parent already names this lever ("unless the pressure rise is
   deliberately slowed"). This analysis turns it into a requirement.
4. **The reflected blast spike is the open hazard, and it breaks the overwrap, not the steel.**
   A sharp spike on the inner face crosses the steel, enters the carbon, and reflects off the
   free outer skin as tension. The resin's 64 MPa across the fibres fails at a spike of only
   **0.12-0.22 GPa**. Steel spall needs 3.2-4.1 GPa. A Sedov blast from the rod's 7 GJ reaches the
   1.68 m wall as a 2,100-3,900 bar shock. Strong normal reflection makes that ~2.8-3.2 GPa. That
   estimate is crude and owed a solved run, but even the incident shock alone exceeds the
   overwrap's threshold. **If it holds, the overwrap delaminates on the first pulse.**
5. **The aluminium layer does not help at these pulse lengths.** It takes 0.8 µs to cross, so
   10-100 µs spikes pass through it as though it were not there (carbon tension within 3%). The
   parent's "splits the step in two" holds only for sub-microsecond spikes. Likewise, the 81%
   reflection at the steel holds only for spikes shorter than the steel's ~4 µs transit. Longer
   spikes pass the full pressure into the carbon.

## 2. Breathing-mode fatigue

The shell's breathing mode, damped single mode with `T_b` 1.10-1.11 ms for the parent's
layups, is driven by a ramp-and-blowdown pressure. Ramps are instant, 0.4 ms (the parent's
D = 1.7) and one acoustic crossing (0.62 ms H2, 0.83 ms CH4). Blowdown e-folding times come from
`wall_heat.csv`. Steel stress is `-562 MPa + 662 MPa x(t)`, the parent's sizing. Rainflow-counted
cycles grow a semicircular surface crack (`K = 0.713 sigma sqrt(pi a)`) until `K_max` reaches 45
MPa m^0.5 or the crack reaches 80% of the wall. Pulses to failure from a 1 mm crack:

| chamber, curve | instant (D ~1.9) | 0.4 ms (D ~1.7) | one crossing (D 1.3-1.5) | rise ≥ T_b (D 1.0-1.2) |
|---|---|---|---|---|
| hydrogen, CC2938 H2 | 31-213 | 82-482 | 321-1,522 | 8,800 to >10^6 |
| methane, air | 313-2,356 | 718-4,565 | 10,800-39,000 | 27,500 to >10^5 |
| methane, CC2938 H2 (bound) | 27-203 | 71-458 | 2,150-8,200 | 7,900 to >10^6 |

The ranges span damping ratio 0.005-0.05 and the two throats. Starting cracks of 0.5 and 2 mm
roughly double and halve these lives (`chamber_fatigue.csv`).

## 3. Spall and delamination under a spike

The model is a linear elastic plane wave through the parent's stack, with exact transmission and
reflection at each impedance step and a free outer skin:
- steel 46 MPa s/m;
- aluminium 17;
- IM7/8552 4.8 (through thickness);
- glass 5.1;
- GRCop-84 ~40 (assumed).

Autofrettage leaves 16.7 MPa of radial compression at the steel/overwrap interface. It is added,
and it is the only compression present when the spike arrives, before the mean pressure builds.

| spike e-folding | steel tension per GPa | carbon tension per GPa | carbon fails at | glass fails at | steel spalls at |
|---|---|---|---|---|---|
| 10 µs | 0.72-0.80 | 0.36-0.39 | 0.20-0.22 GPa | 0.38-0.40 GPa | 3.6-4.0 GPa |
| 30 µs | 0.73-0.91 | 0.46-0.52 | 0.15-0.17 GPa | 0.35-0.39 GPa | 3.2-4.0 GPa |
| 100 µs | 0.71-0.89 | 0.54-0.60 | 0.12-0.14 GPa | 0.33-0.39 GPa | 3.3-4.1 GPa |

The steel spall strength is ~2.9 GPa, from a Q&T Cr-Ni-Mo-V steel (38KhN3MFA, Mescheryakov et
al. 2017, Table 3, `0.5 rho c_b W` at W = 160 m/s). The resin's is 64 MPa across the fibres
(Hexcel 8552, as the parent). The ranges span both chambers, with and without aluminium.

The slow mean-pressure ramp is harmless by comparison. Ramping over 0.4 ms against the carbon's
60-110 µs transit leaves about 7 MPa of radial tension, under the 16.7 MPa bias.

## 4. Limits

- **The fill is uniform gas here, and the design's is not.** The parent's methane chamber pumps
  part of its charge through the liner as coolant. That part fills the chamber as warm gas. The
  rest, from 14.1 kg of 19.6 at the low edge of the heat band to none at the top edge, is sprayed
  in as liquid in the last 0.1 s. 0.5-1 mm drops survive within ~0.8 m of the injectors (parent
  `tab:drop_survival`). The parent places that liquid on purpose, "so that it softens the local
  shocks that reach the wall", and leaves the deposition unmodeled. Neither analysis here
  includes it.
  - The spike thresholds in §3 are wall properties and stand.
  - The spike *estimate*, and the gas-crossing ramp in §2, are for a dry fill.
  - A near-wall droplet layer is the design's own answer to §3's hazard. It is also a candidate
    for slowing the rise in §2, since a two-phase layer carries pressure more slowly than gas.
  - Neither effect is quantified. Where no liquid is left (the top of the heat band), the dry
    numbers apply as they stand.
- **The spike amplitude is an estimate.** It is a strong-shock Sedov blast with normal reflection
  off a rigid flat wall. A decaying, dissociating, curved blast, energy deposited along the rod's
  track, and oblique incidence will all lower it. The thresholds in §3 are the result; the
  estimate only says the hazard is real. A solved blast-in-vessel run (2-D Euler kernel, real EOS)
  is owed.
- The 45 MPa m^0.5 failure `K` is the parent's lower bound for steels under 950 MPa in hydrogen,
  applied to methane too.
- CC2938 constants are from the Sandia slide, which states no units. m/cycle with MPa m^0.5 is the
  only reading in which hydrogen grows cracks faster than air (5-10x here).
- **Not modelled:**
  - Residual pressure between pulses: at 4 Hz the chamber is not empty, which narrows the swing
    and helps.
  - Thermal stress.
  - Thick-wall gradients.
  - The liner's own cyclic plasticity.
  - Fatigue knock-down of the resin under repeated through-thickness tension. No source here;
    composites usually lose a large part of their static strength over 10^3-10^6 cycles.

## 5. What the paper should change (owed to `Balloon-Pulse-Propulsion`)

- **`sec:steel_chamber_service`, after `eq:ramp_dlf`.** "Every shell mass grows by D unless the
  pressure rise is deliberately slowed" becomes a requirement. Add:
  - At the parent's own sizing (D = 1.7, steel -0.75 Y to +0.75 Y), a 1 mm flaw in the hydrogen
    chamber's Cr-Mo shell grows to failure in 82-482 pulses on the ASME CC2938 hydrogen curve.
    That is fewer than one 630-780-pulse departure, because the breathing mode rings 7-74
    tensile cycles per pulse.
  - A pressure rise of at least one breathing period (~1.1 ms) gives D ≈ 1.0-1.2 and >8,000
    pulses.
  - Methane on the air curve survives one departure at D = 1.7 (718-4,565 pulses).
- **`sec:carbon_overwrap`, the spike paragraph.**
  - The 81% reflection and the aluminium's halving apply only to spikes shorter than the steel
    (~4 µs) and the aluminium (~0.8 µs) transits.
  - For 10-100 µs spikes the full pressure reaches the carbon. It reflects off the free outer
    skin as 0.36-0.60 of the spike in tension, so the overwrap delaminates at a 0.12-0.22 GPa
    spike.
  - A first estimate of the reflected blast at the 20 m^3 wall is ~2-3 GPa. The section should
    say the overwrap's survival depends on keeping the reflected blast below ~0.1 GPa, and that a
    solved blast-in-vessel load is owed. Candidate mitigations, none sized: a crushable inner
    layer that attenuates the reflected shock; through-thickness reinforcement of the overwrap;
    a thicker tough outer skin that moves the reflection tension out of the carbon.
- **`sec:carbon_overwrap`, the liquid-share paragraph ("softens the local shocks").** That
  softening is load-bearing. Without it the overwrap's 0.12-0.22 GPa threshold sits at or under
  the paper's own ~0.2 GPa incident estimate, before any reflection. The paper should say the
  overwrap's survival depends on that liquid share, and that the share falls to zero at the top
  of the heat band.
- **The "not a high-cycle fatigue part" sentence** stays right on cycle count. It should add that
  the defect-tolerant check, not the count, is what sizes the hydrogen wall.

## 6. Sources

- ASME BPVC VIII-3 Code Case 2938 master curve as given by J. Ronevich et al., DOE Hydrogen
  Program Annual Merit Review 2023, SCS005, slide 8.
  https://www.hydrogen.energy.gov/pdfs/review23/scs005_ronevich_2023_p.pdf
- Yu. I. Mescheryakov et al., "Shock-induced structural instabilities and spall-strength of
  maraging steels", Mater. Phys. Mech. 32, 152-164 (2017), Table 3 (38KhN3MFA).
- J. M. Barsom and S. T. Rolfe, *Fracture and Fatigue Control in Structures*: martensitic-steel
  air curve. Not re-fetched.
- Parent `templateArxiv.tex`: `sec:steel_chamber_service`, `sec:carbon_overwrap`,
  `sec:methane_7000_near_term`, for the layups, impedances, 45 MPa m^0.5, 64 MPa and 950 MPa
  limits (with its citations `sanmarchi2017_thickwall`, `hexcel_8552_datasheet`,
  `alexander2012_cfrp_planar_impact`, `evident_sound_velocities`).
