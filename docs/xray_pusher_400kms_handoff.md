# Pass-through pusher at 400 km/s: a light foam on a dense anvil clears 2650 s

Handoff to `Balloon-Pulse-Propulsion` (parent), §`pass_through_pusher`. It is written to be copied
there. It asks whether the pass-through pusher works at a 400 km/s closing speed, a 200 km/s
periapsis at about Parker Solar Probe's distance met head-on, rather than at the 618 km/s of the 4
solar radii dive. The evidence is [`xray_collision_study.md`](xray_collision_study.md) Results 7-9:
- `make xray-collision` and `make xray-collision-anvil`, written to
  `data/results/xray_collision/near_sun_400.csv`, `anvil_*.csv` and `propellant_rebound.csv`;
- `make xray-collision-rebound`, the `f(v)` study's 1-D radiation-hydrodynamic bounce, run on a
  new polyethylene table.

**Status: provisional.** The collision is a one-zone fireball, using a Saha EOS and TOPS (LANL
OPLIB) gray opacities. The case is a cavity energy balance. The rebound is a 1-D Lagrangian
bounce off a rigid plate. Expect factor-of-order-unity errors in the collision and the case, and
tens of percent in the rebound.

## Conclusion for the paper

At 400 km/s the pass-through pusher still works, but with the density contrast reversed. The
parent sends a dense disk into a light lattice. Here a **light carbon-foam projectile** stops on a
**dense graphite anvil** inside the radiation case. In any collision the light partner takes the
heat, so the foam stagnates on the anvil at about 700 eV. Half of the projectile's kinetic energy
leaves as light. The case delivers 54-78% of that light to the channel filler, depending on the
entry hole. The polyethylene propellant plate then rebounds from the pusher plate with
`e_eff = 0.63-0.72`. The parent's ledger does not count that rebound. With it, specific impulse
is **2740-3300 s**, above the 2650 s bar. Without it the design fails at 1680-1920 s.

## Why not the parent's design at 400 km/s

A dense carbon disk into a carbon lattice (Result 7) still radiates about half its kinetic energy
at 400 km/s, but only from a lattice of about 1 kg/m^3. That is ten times thinner than the
parent's 10 kg/m^3, and at 10 kg/m^3 the share falls to 18%. A 1 kg/m^3 target of 2 kg is a
column 16 m long behind a 40 cm disk. The light is also softer: it is emitted near 110 eV, and
3% of the kinetic energy is above 1 keV, against 28% at 618 km/s. A steel lattice is no better.
Iron gets more energy per atom, but it shares it among about 20 free electrons and pays several
keV to free them, so it starts cooler than carbon.

## The design

| element | choice |
|---|---|
| Projectile | 10 kg in total. A carbon-foam cylinder at 10 kg/m^3, **0.5 m across and about 3.8 m long** (1 m across is the simpler fallback). Like Parker's heat shield, it is a sandwich: C-C sheets of 0.5 mm or less on the **sunward side only** (2.6 kg at 0.5 m, 1.5 kg at 1 m) and **no front cap**. |
| Anvil | **50 kg of graphite** (k = 5), as wide as the foam, faced with a thin high-Z liner (tungsten or lead) for albedo. It is spent on every pulse and leaves astern at 67 km/s. |
| Case | Short and wide, e.g. 1.6 m across and 0.5 m deep. Its walls and liner take the 1.5 kg/kg left in the parent's 6.5 kg/kg case budget. |
| Filler and propellant | A low-Z channel filler and a polyethylene propellant plate **of equal mass** (3-4 kg each per kg of projectile). The propellant strikes the plate at 82-94 km/s, under the 150 km/s cap. |

Choices that were tried and dropped:
- **A z-pinch or other magnetic conversion.** It adds two lossy stages, kinetic to electrical and
  back, ahead of the same radiating collision. The collision's self-generated field is about
  0.02 T, so a seed field and a ship-held circuit would be needed.
- **An iron anvil.** It is too thin per unit mass to outlast the pulse against Rayleigh-Taylor
  burn-through.
- **Spaced blanket foils in the case.** They explode within nanoseconds at a ~250 eV drive.
- **Heavy material behind a light propellant.** Eq. `pulse_split` is symmetric in the two
  masses. A light-forward, heavy-backward stack of the same total mass loses 13% of the impulse
  and sends the propellant at 172 km/s, over the cap.

## Per-kilogram-of-projectile ledger at 400 km/s

The projectile's kinetic energy is 80 GJ per kg. The table uses the parent's own pulse-split
optimizer, which reproduces `tab:pass_through_isp` to within 0.2%.

| | 1 m hole | **0.5 m hole** |
|---|---|---|
| light (sheets counted as dead mass) | 0.52 of KE | 0.46 |
| delivered to the filler (case 1.6 x 0.5 m, albedo 0.95) | 0.54 | **0.78** |
| propellant speed at the plate | 82 km/s | 94 km/s |
| Isp without rebound (the parent's 0.8 catch) | 1680 s | 1920 s |
| **Isp, rebound at a 1 kg/m^2 column (`e = 0.63`)** | **2740 s** | **3130 s** |
| Isp, rebound at a 3 kg/m^2 column (`e = 0.72`) | 2890 s | 3300 s |
| spent mass (anvil, walls, filler, propellant) | 13 kg | 13 kg |
| energy efficiency, `(1/2) M v_e^2 / KE`, with rebound | 6-6.5% | 7.6-8.5% |

For comparison, the parent's 618 km/s design with a 5 kg case gives 2930 s from 14 kg spent, with
no rebound credit. This design's light is as hard: the foam stagnates near 700 eV, and 37% of
the bare foam's kinetic energy leaves above 1 keV, against 28% at 618 km/s. That design needs the
4 solar radii dive; this one needs a 200 km/s periapsis.

## Three findings the parent's text should take

1. **Count the rebound.** The parent takes "the plate catches 80% of the forward momentum", with
   no bounce. The hydrocode gives a polyethylene propellant at ~95 km/s
   `e_eff = 0.63-0.72` for an arriving column of 1-3 kg/m^2. The plate then catches
   `0.8 (1 + e) = 1.30-1.38` times the forward momentum. This is the largest lever in the ledger:
   a factor of 1.6 in impulse. It holds only for a column of at least ~1 kg/m^2 and a propellant
   well under the cap. At 0.1 kg/m^2 the rebound falls to 0.44, and at 150 km/s it falls to 0.3-0.6.
   That puts a limit on plate area: 30-40 kg of propellant per pulse needs a plate under about
   30 m^2.
2. **The case's 80% delivery is reachable, and the entry hole sets it.** Light leaves a stagnated
   parcel within nanoseconds, so emission follows the projectile's arrival. While the foam's cold
   tail is still in the entry hole, it blocks it. This helps only if the foam is **longer than the
   case**. A 1 m foam in a short case delivers 53%, and a 0.5 m foam delivers 78%. Below ~0.35 m
   the returns are small.
3. **Spread the projectile, and make the anvil dense.** A dense projectile is the right choice
   when the ship's target is light: that is the parent's lattice design. When the ship holds a
   dense stopping mass, the projectile should be the light one. Every dense gram in the projectile
   collides dense-on-dense, stays opaque and is lost: the C-C sheets cost 7-13 points of light.

## Suggested text for the parent

For a closing paragraph of §`pass_through_pusher`:

> The pusher does not need the full 4 solar radii dive. At a 200 km/s periapsis met head-on, a
> 400 km/s closing speed, a dense disk would need a lattice ten times thinner than at 618 km/s.
> The roles can be reversed instead. A light carbon-foam projectile, sheathed like the Parker Solar
> Probe's heat shield on its sunward side, stops on a 50 kg graphite anvil at the far end of the
> case. The foam takes the heat, stagnating near 700 eV, and about half of its kinetic energy
> leaves as light \cite{Katz_xray_collision_2026}. The case delivers 54-78% of that light into the
> filler. Its tail blocks the entry hole while it arrives, so a narrower projectile loses less. The
> polyethylene propellant does not stick to the plate. A radiation-hydrodynamic bounce gives it a
> restitution of 0.63-0.72 at 95 km/s, so the plate takes 1.3-1.4 times the propellant's
> momentum. The pusher then reaches 2740-3300 s, above the 2650 s the boost requires. Without
> that rebound it would reach under 2000 s.

## Not modelled; open for the parent to flag

- **Case-wall mass.** This is the largest risk. A ~5 Mbar fireball held for microseconds needs
  about 20 kg/m^2 of wall by inertial tamping, ~100 kg for this case, against the ~15 kg the
  budget leaves. The 0.5 m foam lengthens the pulse to 9.4 us, and in that time the plasma spreads
  ~2.5 m sideways, which makes this worse.
- **Plate ablation.** The plate is rigid and inert. At ~95 km/s the `v/4` rule asks it to absorb
  about 24 kJ per N s.
- **The arriving propellant** is modelled as a uniform cold slug. A real one is a stretched jet,
  preheated by the drive.
- **The case balance** holds at each instant and has no shadowing. The in-flight foam's albedo
  (0.5-0.9) and the wall albedo (0.95) are estimates. The wall albedo should be computed from a
  Marshak wave in the liner material.
- **The collision** is one zone with an averaged stagnation. Each parcel is actually pinned
  against the anvil by about 24 Mbar of ram pressure as it radiates, which probably raises the
  radiated share.
- **The sheets' survival** in a ~10 solar radii dive (about 0.6 MW/m^2, a face near 1800 K) is a
  hand estimate.
