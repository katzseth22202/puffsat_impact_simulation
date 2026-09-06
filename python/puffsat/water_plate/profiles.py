"""Prescribed mass-flux profiles and kinematic scales, not a hydrodynamic solver.

All radii/fluxes refer to projected area. Time zero is when the unperturbed incoming
front crosses the plate-vertex plane. Pre-spray can intercept it before that time.
Injection is axial toward the incoming pulse (-z in the study's prograde frame).
Temperature, phase, pressure work, and collision-powered gasification are not
assigned by this kinematic specification; see the study's open source-energy task.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from puffsat.water_plate.ledger import PROJECTILE_KG, REFERENCE_K, _check

OUTPUT_DIR = Path("data/results/water_plate")


@dataclass(frozen=True)
class IncomingPulse:
    """Uniform cylindrical incoming mass, purely axial with a rectangular time profile."""

    radius_m: float
    duration_s: float
    speed_m_s: float = 45_580.0
    mass_kg: float = PROJECTILE_KG

    def __post_init__(self) -> None:
        for name in ("radius_m", "duration_s", "speed_m_s", "mass_kg"):
            _check(name, getattr(self, name), positive=True)

    @property
    def length_m(self) -> float:
        return self.speed_m_s * self.duration_s

    @property
    def density_kg_m3(self) -> float:
        return self.mass_kg / (math.pi * self.radius_m**2 * self.length_m)


@dataclass(frozen=True)
class TimeBand:
    """Fraction of supplied mass uniformly injected over [start_s, end_s)."""

    start_s: float
    end_s: float
    fraction: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.start_s) or not math.isfinite(self.end_s):
            raise ValueError("band endpoints must be finite")
        if self.end_s <= self.start_s:
            raise ValueError("band must have positive duration")
        _check("fraction", self.fraction)
        if self.fraction > 1:
            raise ValueError("mass fraction exceeds one")


@dataclass(frozen=True)
class InjectionProfile:
    """Axisymmetric projected-area source proportional to r**radial_power.

    Power 0 is uniform, power 2 weights the outer radius. Its normalized cumulative
    radial mass fraction is (r/R)**(power+2), which preserves mass with either
    footprint-matched or full-plate injection. Bands must be disjoint.
    """

    radius_m: float
    radial_power: int
    bands: tuple[TimeBand, ...]
    water_kg: float = PROJECTILE_KG * REFERENCE_K

    def __post_init__(self) -> None:
        _check("radius_m", self.radius_m, positive=True)
        _check("water_kg", self.water_kg)
        if self.radial_power not in (0, 2):
            raise ValueError("radial_power must be 0 (uniform) or 2 (outer-weighted)")
        if not math.isclose(sum(b.fraction for b in self.bands), 1.0, abs_tol=1e-12):
            raise ValueError("time-band mass fractions must sum to one")
        ordered = sorted(self.bands, key=lambda b: b.start_s)
        if any(a.end_s > b.start_s for a, b in pairwise(ordered)):
            raise ValueError("time bands must be disjoint")

    def mass_flux_kg_m2_s(self, radius_m: float, time_s: float) -> float:
        """Prescribed source flux per projected square meter (not shocked gas density)."""
        _check("radius_m", radius_m)
        if not math.isfinite(time_s):
            raise ValueError("time_s must be finite")
        if radius_m >= self.radius_m:
            return 0.0
        temporal = sum(
            b.fraction / (b.end_s - b.start_s) for b in self.bands if b.start_s <= time_s < b.end_s
        )
        radial = (self.radial_power + 2) / (2 * math.pi * self.radius_m**2)
        radial *= (radius_m / self.radius_m) ** self.radial_power
        return self.water_kg * temporal * radial

    def mass_between(self, start_s: float, end_s: float, radius_m: float) -> float:
        """Exact supplied mass inside radius and time window; no interception prediction."""
        _check("radius_m", radius_m)
        if not math.isfinite(start_s) or not math.isfinite(end_s) or end_s < start_s:
            raise ValueError("integration endpoints must be finite and ordered")
        temporal = sum(
            b.fraction
            * max(0.0, min(end_s, b.end_s) - max(start_s, b.start_s))
            / (b.end_s - b.start_s)
            for b in self.bands
        )
        radial = min(1.0, radius_m / self.radius_m) ** (self.radial_power + 2)
        return self.water_kg * temporal * radial


def main() -> None:
    """Write bounded screening inputs. None of these rows is a simulated outcome."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_DIR / "profile_pulses.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "diameter_m",
                "footprint_ratio",
                "duration_s",
                "closing_speed_m_s",
                "projectile_kg",
                "length_m",
                "incoming_density_kg_m3",
                "incoming_ram_pa",
            )
        )
        for diameter in (10.0, 20.0, 30.0):
            for coverage in (0.5, 0.75, 1.0):
                for duration in (10e-6, 30e-6, 100e-6, 300e-6, 1e-3):
                    pulse = IncomingPulse(diameter * coverage / 2, duration)
                    writer.writerow(
                        (
                            diameter,
                            coverage,
                            duration,
                            pulse.speed_m_s,
                            pulse.mass_kg,
                            pulse.length_m,
                            pulse.density_kg_m3,
                            pulse.density_kg_m3 * pulse.speed_m_s**2,
                        )
                    )
    reference = IncomingPulse(3.75, 100e-6)
    with (OUTPUT_DIR / "profile_injection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "radial_mode",
                "timing_mode",
                "source_radius_m",
                "radial_power",
                "pre_fraction",
                "pre_start_s",
                "pre_end_s",
                "during_start_s",
                "during_end_s",
                "water_kg",
                "axial_injection_speed_m_s",
                "peak_total_mass_rate_kg_s",
                "peak_area_mean_stream_density_kg_m3",
                "direct_radial_overlap_fraction",
                "injection_kinetic_energy_j",
                "injection_momentum_scale_ns",
            )
        )
        for radial_mode, radius, power in (
            ("footprint-uniform", 3.75, 0),
            ("plate-uniform", 5.0, 0),
            ("plate-outer-weighted", 5.0, 2),
        ):
            for timing_mode, pre_fraction in (("pre", 1.0), ("during", 0.0), ("split", 0.5)):
                bands = (
                    TimeBand(-0.002, -0.001, pre_fraction),
                    TimeBand(0, reference.duration_s, 1 - pre_fraction),
                )
                profile = InjectionProfile(radius, power, bands)
                rate = profile.water_kg * max(b.fraction / (b.end_s - b.start_s) for b in bands)
                overlap = profile.mass_between(-1, 1, reference.radius_m) / profile.water_kg
                for speed in (10.0, 100.0):
                    writer.writerow(
                        (
                            radial_mode,
                            timing_mode,
                            radius,
                            power,
                            pre_fraction,
                            bands[0].start_s,
                            bands[0].end_s,
                            bands[1].start_s,
                            bands[1].end_s,
                            profile.water_kg,
                            speed,
                            rate,
                            rate / (math.pi * radius**2 * speed),
                            overlap,
                            0.5 * profile.water_kg * speed**2,
                            profile.water_kg * speed,
                        )
                    )
    print(f"Kinematic profiles written to {OUTPUT_DIR}; source thermodynamics remains open.")


if __name__ == "__main__":
    main()
