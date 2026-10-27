# Request to puffsat_impact_simulation: a thicker, heavier Sun disk

For the companion repo, then the parent, then a one-line copy down here. Raised 2026-10-27.
The aim questions this raises are in `todos/thick_disk_aim_request.md`, for
`aim_is_all_you_need`.

## Why

The Sun boost at the 9.86 R_sun baseline flies 1 kg graphite disks into charges in a bowl with
no hole, closing head-on at about 392 km/s. The disk inherited its shape from the parent's
4 R_sun X-ray pusher: 40 cm across, 3.6 mm thick, 8 kg/m².

A thicker disk is cheaper and tougher. Take 1 cm, or 22 kg/m² at the companion's
2200 kg/m³. Our rough reasons, for the companion to correct rather than adopt:

- Cratering depth scales roughly as d * v^(2/3). The dust grain that spalls the back face then
  grows from roughly 20-60 um to 55-170 um. Grains that large are 15-25 times rarer.
- Bending stress under launch and handling falls as 1/t^2. Diameter over thickness drops from
  about 110 to about 38.
- At 9.86 R_sun the flux is about 0.65 MW/m². A face-on disk settles near 1550 K and a
  sun-facing rim near 1840 K, so thickness costs nothing thermally.
- 1 cm plate cuts from molded block with less scrap, and may allow a coarser, cheaper grade.

The disk is also slower than the parent's. The kinetic energy per kilogram at 392 km/s is
(392/618)^2 = 0.40 of the 618 km/s value. A disk of 2.5 kg at 1 cm therefore carries the
same lab-frame energy as the parent's 1 kg disk, 192 GJ against 191 GJ. It is 38 cm across,
nearly the parent's 40 cm footprint. If the charge grows with the disk, aim precision stays
where it was.

## Cadence

The local paper gives 7.2 t of 1 kg disks (7,200 disks) and a cadence of 2 per second over
4.9 h (35,280 pulses). These do not match. The proposed cadence is 1 Hz. At 2.5 kg the same
7.2 t is 2,880 disks, a 48 minute burn. Each disk gives the ship about 3.4 m/s, if the
impulse per kilogram of disk stays at the baseline's. The average is 0.35 g. Run items 2 and
3 at that cadence.

## Ask

1. **Sweep or punch-through.** Does a 22 kg/m² disk still couple to the Sun bowl's charge at
   392 km/s, or does it tunnel? The parent's snowplow argument wants a target column of about
   k times the disk's areal density, which at 1 cm is 2.75 times the column the 8 kg/m² disk
   needs. That argument was made for the X-ray lattice, so say whether it carries over to the
   bowl charge at all. Report f (or eta) for 8 and 22 kg/m² at the same disk mass.
2. **Scaled pulse.** Run disks of 1, 2 and 2.5 kg at 1 cm, each with its charge scaled to keep
   the current charge-to-disk mass ratio. Report f, peak bowl load, and the energy per pulse
   reaching the bowl. The 2.5 kg case should match the parent's per-pulse energy. Say whether
   the 5 t bowl takes a 3.4 m/s kick to the ship each second, or what it needs instead.
3. **Off-center tolerance.** For the 2.5 kg case, the lateral offset (as a fraction of charge
   radius) at which f drops by 5% and 20%, and the side impulse at each. The aim repo turns
   this into a guidance requirement.
4. **A split disk.** Optional. A disk that a large grain has cracked in two arrives as two
   halves a few cm apart. Does the pulse care, or does it couple as one body?

## Loads we checked and do not ask about

Neither threatens the graphite. Both are rough, from textbook values, and are here so the
companion need not redo them. Graphite breaks at tens of MPa.

- **Coronal mass ejection.** Ambient wind at 10 R_sun is roughly 1e3-1e4 protons/cm³. A CME
  can run 10-100 times denser at up to about 2000 km/s. The pessimistic dynamic pressure is
  about 1e-3 Pa. Plate bending stress is about p * (radius/thickness)^2, or 0.4 Pa. Sunlight
  at 9.86 R_sun pushes harder, about 2e-3 Pa, and is also harmless.
- **Cold gas thrusters.** A CubeSat-class 10-100 mN thruster on 2.5 kg gives about
  0.04 m/s². Even 1 N through a 1 cm² mount is about 10 kPa. Spin at 1 rad/s gives hoop
  stress of tens of Pa, and at 10 revolutions per second about 0.1 MPa.

The loads that can break a disk are launch, the thermal mismatch at the pod and brow mounts
(graphite near 1500 K against an insulated pod), and a rare large grain.

## Alumina coat, probably not needed

Parker's white coat is alumina on carbon/carbon. Graphite takes it at least as well. Alumina
expands about 8 ppm/K. Molded isotropic graphites run about 4-8 ppm/K, and C/C is near zero
in-plane, so graphite is the closer match. Plasma-spraying ceramics onto graphite is routine
industrial practice. Hot carbon reduces alumina in vacuum at high temperature on either
substrate.

The disk does not need the coat at 9.86 R_sun. Bare graphite runs about 1550-1840 K and
sublimes almost nothing, and the parent already declines to count on a white coat. The coat
adds cost, would chip under dust at about 300 km/s, and would put its alumina into every pulse.
It stays available for the sunshade brow if the pod needs it.

5. **Optional, only if a coat flies.** Each 100 um of alumina at 3950 kg/m³ adds about
   0.4 kg/m², 1.8% of the 22 kg/m² disk. Does that much aluminum and oxygen in the pulse
   change f measurably?

## Outcome

If items 1 and 2 hold, the Sun disk becomes a 1 cm, 2.5 kg molded-graphite plate with the
parent's footprint and aim, and C/C drops out of consideration. If the thick disk tunnels,
the disk stays thin and the dust case gets a proper fluence estimate instead.
