# The Earth-departure overtake pusher: a 100 km/s impactor onto a charge held off the plate

Handoff to `Balloon-Pulse-Propulsion` (parent), for the departure wave. Evidence:
[`xray_collision_study.md`](xray_collision_study.md) Result 11 and
[ADR-0061](adr/0061-overtake-pusher-charge-off-the-plate.md);
`make xray-collision-overtake`, written to `data/results/xray_collision/overtake_ledger.csv`.

**Status: provisional.** The plate impulse comes from 2-D axisymmetric Euler with an effective
gamma of 1.4 on two grids, which agree to 2-8%. The finer grid gives more impulse, so the numbers
here are conservative. Polyethylene's radiative rebound is swapped in from the f(v) study's 1-D
bounce. Freeze-out comes from a one-zone model on a Saha/TOPS polyethylene table. The peak pressure
of the compact and sheet collisions is not converged.

## Conclusion for the paper

To give a departing ship 11 km/s with a third of its starting mass as propellant (mass ratio 1.5),
use **k ≈ 10-12 kg of polyethylene per kg of impactor**:
- k ≈ 10-11.5 on the parent's 20 m bowl;
- k ≈ 12-14 on a δ/D 0.5 paraboloid with the charge at 1.5 focal lengths.

The impactor overtakes at 100 km/s and strikes a charge a few metres behind the plate, never one
on it. The merged mass drifts onto the plate. The plate keeps 60-74% of the expansion's ideal
momentum, plus the impactor's own, which is better than the 400 km/s charge pusher. The fireball
is dark, ionization recombines, and H2 re-forms before the plate; only carbon's share of the bonds
(a 6% velocity loss at k = 10) is in doubt. k = 40 does not work: the mass ratio is 2.1-3.3.

## Ledger (mass ratio for an 11 km/s boost; carbon left atomic)

| plate | geometry | k = 10 | k = 20 | k = 40 | impactor per pulse (400 MPa) |
|---|---|---|---|---|---|
| parent bowl | compact charge at focus | 1.49-1.50 | 1.94-2.00 | 3.23 | 2.4-6 kg |
| δ/D 0.5 paraboloid at 1.5 F | compact charge | **1.42** | 1.79 | 2.76 | 3-5 kg |
| flat + 4 m skirt | 7 m impactor sheet onto a 7 m charge sheet | **1.46** | 1.86 | | ~11 kg |

- **Cost:** the impactor is 0.042-0.049 kg per kg of final mass. At $1000/kg that is ~$45 per kg of
  final mass, or $70-100 per kg of propellant.
- **The planar sheet buys loading, not impulse.** The peak scales with the impactor's mass per area,
  so the plate takes ~0.05-0.07 kg/m^2 of impactor at 400 MPa: ~35-50 µm of polyimide or carbon
  fibre. The ship's charge sheet is ~0.5-0.75 mm of polyethylene.
- **The sheet must last from 0.5 to 1 AU**, at ~470 K face-on at 0.5 AU.
  - Polyimide film (IKAROS-proven) is the first choice. Carbon fibre also works, in a polyimide or
    dry matrix. Kevlar only with a UV coating.
  - It is a solar sail (~0.3-0.5 mm/s^2 at 0.5 AU), so it needs spin for flatness, a controlled
    attitude and late guidance. Arrival tolerances are loose: ~1 µs of spread against a ~100 µs
    pulse.

## Suggested text for the parent

> The departure wave uses the same cheap charge. Each impactor overtakes the ship at 100 km/s and
> strikes a polyethylene charge, or a wide polyethylene sheet, held a few metres behind the
> plate. The merged mass drifts onto the plate as a dense, dark gas. Two-dimensional hydrodynamics
> with polyethylene's radiative rebound gives an 11 km/s departure at a mass ratio of 1.5 with
> 10-12 kg of polyethylene per kg of impactor \cite{Katz_xray_collision_2026}. A wide, thin
> impactor sheet spreads the load, so a 20 m plate takes ~11 kg per pulse at 400 MPa.

## Not modelled; open for the parent to flag

- Polyethylene's real EOS through dissociation, and carbon soot, in 2-D (the toll is bracketed).
- The sheet collision's resolved peak, and real sheets: sag, holes, a spin-flattened cone.
- Plate ablation at the 20-45 km/s arrival speed.
- Deploying and holding the ship's charge sheet, and the impactor sheet's guidance.
