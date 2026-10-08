"""X-ray pusher at rocket speeds: how much of a collision's heat leaves as light.

The parent's pass-through X-ray pusher (templateArxiv.tex §why_not_xray_drive, §pass_through_pusher)
was written for 618 km/s near-Sun dives. This study asks the same question at the rocket case: a
56 km/s impactor meeting a ship moving at it at 11 km/s, so 67 km/s closing in the ship's frame.

A ship-frame target of `k` times the impactor's mass turns `k/(1+k)` of the impactor's kinetic
energy into heat; the rest is centre-of-mass motion that leaves with the debris. Of that heat,
the question is what fraction radiates before expansion turns it into bulk motion, and at what
photon energy. See `docs/xray_collision_study.md`.

Entry point: `fireball.collide`.
"""
