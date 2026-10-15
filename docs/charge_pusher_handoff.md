# The charge pusher: a graphite disk into cheap polyethylene, onto the parent's own plate

Handoff to `Balloon-Pulse-Propulsion` (parent), §`pass_through_pusher`, as a simpler alternative
to the X-ray pusher for the near-Sun boost. It is written to be copied there. The evidence is
[`xray_collision_study.md`](xray_collision_study.md) Result 10:
- `make xray-collision-charge-2d` and `make xray-collision-charge`, written to
  `data/results/xray_collision/charge_ledger.csv`;
- the polyethylene rebound runs of Result 9, extended to 20-50 km/s.

**Status: provisional.** The plate impulse comes from 2-D axisymmetric Euler with an effective
gamma of 1.4 and an immersed plate. It converges to ~1% on a twice-finer grid. Polyethylene's
radiative rebound comes from the f(v) study's 1-D radiation-hydrodynamic bounce, and the
fireball's radiation from a one-zone model on a Saha/TOPS polyethylene table. Peak pressure is not
converged: a finer grid gives 18-30% less, so the values here are upper ones.

## Conclusion for the paper

The near-Sun boost does not need X-rays or a radiation case. A graphite disk at the 400 km/s
closing speed of a Parker-like periapsis (9.86 solar radii, 200 km/s met head-on) threads the
bore and stops in a compact **charge** of near-solid polyethylene. The charge is 20-100 times
the disk's mass, and is carried by the ship and set a few metres astern. The fireball is millions
of optical depths thick, so it radiates nothing and its whole heat becomes expansion onto the
plate. On the parent's own 20 m bowl (depth over diameter 0.30), unchanged, the boost from a
Jupiter fall to a 100 km/s Earth crossing needs a mass ratio of **2.0 at k = 30 and 2.9 at
k = 100**. A deeper, δ/D 0.5 paraboloid lowers this to 1.7-2.3. The projectile bill falls to
**0.013-0.034 kg of projectile per kilogram of final mass**, because the cheap charge, not the
projectile, is the reaction mass.

## The design

| element | choice |
|---|---|
| Projectile | A graphite disk that fits the ~0.5 m bore, 1-2.5 kg per pulse: 1 mm thick is 0.43 kg, 2.3 mm is 1 kg. The guidance and its foam shade sit behind it. |
| Charge | Near-solid polyethylene, `k` = 20-100 times the disk's mass, as a compact cylinder the disk's diameter (a column of ~k x 2.2 kg/m^2: ~7 cm at k = 30 for a 1 mm disk). It has a thin graphite skin and is kept in the ship's shadow until the pulse. |
| Plate | The parent's 20 m bowl, δ/D 0.30, charge at its focus 4.2 m from the vertex. Upgrade: a δ/D 0.5 paraboloid with the charge at 1.5x its focus, 3.75 m. |
| Mission | Fall from Jupiter to 9.86 R_sun, a 20.8 km/s boost there, Earth crossed at 100 km/s. A mass ratio of 3 needs Isp 1930 s. |

Three rules carry the design:
1. **The charge must be compact and in the disk's path.** A projectile couples only to the column
   it crosses. A loose cloud of foam balls is mostly missed. A long foam column becomes a line
   source and explodes sideways.
2. **Dense means dark.** The fireball cannot radiate, so the energy leaves as motion, which is
   what a plate wants. This reverses every X-ray design, whose difficulty was getting light out.
3. **Enclosure, not focusing.** Moving the charge off a paraboloid's focus costs 3% of capture
   and cuts peak pressure 2.8 times. Any deep bowl works, so the cheapest deep bowl to build
   should win.

## Ledger (per 1 kg projectile, impulsive boost; a 10 h arc adds ~1% to the mass ratio)

| plate | k | Isp | mass ratio | projectile per kg final | peak pressure |
|---|---|---|---|---|---|
| parent bowl | 30 | 2990 s | 2.03 | 0.034 | 133-163 MPa |
| parent bowl | 100 | 2015 s | 2.86 | 0.019 | 239 MPa |
| bowl + 4 m skirt | 30 | 3540 s | 1.82 | 0.027 | 133-162 MPa |
| bowl + 4 m skirt | 100 | 2330 s | 2.48 | 0.015 | 239 MPa |
| δ/D 0.5 paraboloid at 1.5 F | 30 | 3860 s | 1.73 | 0.024 | 195-280 MPa |
| δ/D 0.5 paraboloid at 1.5 F | 100 | 2520 s | 2.32 | 0.013 | 405 MPa |
| flat plate + 4 m skirt | 30 | 2940 s | 2.06 | 0.035 | 465 MPa |

For comparison, the parent prices its magnetic nozzle at 0.064-0.068 kg of projectile per kg of
final mass, for a 4 solar radii dive to a 150 km/s crossing. That is a different mission, so the
comparison is indicative only.

How the numbers move:
- **More charge per projectile lowers the projectile bill and the Isp.** From k = 30 to 100 the
  bill falls 45% while the Isp falls a third. The mass ratio stays under 3 throughout.
- **Polyethylene rebounds from the plate with `e_eff` 0.6-0.7** at the speeds here. The parent's
  pass-through ledger takes no rebound.
- **Peak pressure, not k, caps the projectile at ~1-2.5 kg per pulse** on a 20 m plate against a
  400 MPa facesheet. It goes as the projectile's energy over the plate's size cubed. A charge 1.5
  focal lengths off the bowl's vertex takes ~7 kg per pulse, at about 15% less capture. A 40 m
  plate takes eight times the projectile.
- At 1 kg per pulse, a 100 t final mass needs ~1500-2700 pulses, a few seconds to ~13 s apart
  across a ~10 h periapsis arc.

## Suggested text for the parent

For §`pass_through_pusher`, a closing paragraph or an alternative subsection:

> A simpler pusher needs neither X-rays nor a radiation case. Each graphite disk threads the bore
> and stops in a compact charge of polyethylene, tens of times its own mass, that the ship sets a
> few metres astern. Polyethylene is the cheapest mass the ship carries. The merged fireball is
> so dense that it cannot radiate, so its whole heat becomes expansion, and the parent's own 20 m
> bowl catches it as a pressure-borne flow. Two-dimensional hydrodynamics, with polyethylene's
> radiative rebound from the plate, gives 2000-3000 s from the unchanged bowl
> \cite{Katz_xray_collision_2026}. That is enough for a boost at Parker Solar Probe's periapsis,
> on a fall from Jupiter, to cross Earth's orbit at 100 km/s with a mass ratio of 2-2.9. The
> projectiles are the expensive part of the architecture, and here they are 1-3% of the final
> mass, because the cheap charge, not the projectile, is the reaction mass.

## Checked since (Result 10, ADR-0060)

- **Freeze-out.** Ionization (~18% of the heat) recombines before the plate at every k, and leaves
  ~0.03% frozen. Polyethylene's bonds do not re-form before the plate, because the gas is still
  0.4-1.4 eV there. They cost 1.6% of exhaust velocity at k = 30 and 5.4% at k = 100, charged in
  full.
- **Mixing.** The real disk-into-charge collision, in 2-D, gives -4% to +25% *more* plate impulse
  than the uniform fireball for a solid charge, so the ledger above is conservative. A 30 kg/m^3
  foam charge loses 15-39% as a line source, so the charge must be solid.
- **Uneven exit speeds.** Without a nozzle the exhaust's direction factor is 0.70-0.89 and its
  speed uniformity 0.64-0.83. In jet-efficiency terms that is 0.4-0.6 against a magnetic
  nozzle's assumed 0.8-1. The plate impulse already includes this, since it comes from the wall
  pressure of the full 2-D flow.

## Not modelled; open for the parent to flag

- **The bore collar.** The real collision throws a narrow back-splash jet down the axis. It lands
  just outside the bore (r 0.26-0.56 m) early in the pulse, at 5-20 GPa per kg of projectile on
  the production grid and 2-7 GPa on a grid twice as coarse. It is under-resolved, so the true
  value is at least this. The ring out to ~0.6-1 m needs a sacrificial or ablative collar, or a
  charge shaped to tamp its own splash. Elsewhere the plate sees ~0.1-0.3 GPa per kg.
- **The bore's flux toward the ship:** about five times the point fireball's share, not yet
  reported as energy.
- **Plate ablation** at 40-85 km/s: the `v/4` rule gives 10-21 kJ per N s.
- **The charge's thermal survival** in the ship's shadow until it is set out.
