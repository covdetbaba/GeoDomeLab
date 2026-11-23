"""Circular hollow section (CHS / round HSS) catalogue.

Geodesic domes are almost universally built from round tube, so a single
section family is enough and it keeps the design checks clean: a CHS is
doubly symmetric, so the radius of gyration is the same about every axis and
there is only one slenderness ratio per strut.

Properties are computed from the outside diameter ``D`` and wall thickness
``t`` rather than being tabulated, which removes any chance of a transcription
error:

    A  = pi/4 * (D^2 - d^2)
    I  = pi/64 * (D^4 - d^4)
    r  = sqrt(I / A)
    S  = 2 I / D
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Section:
    """A circular hollow section.  ``D`` and ``t`` are stored in metres."""

    name: str
    D: float
    t: float

    # ----------------------------------------------------------- geometry --
    @property
    def d_inner(self):
        return self.D - 2.0 * self.t

    @property
    def A(self):
        """Cross-sectional area [m^2]."""
        return math.pi / 4.0 * (self.D ** 2 - self.d_inner ** 2)

    @property
    def I(self):
        """Second moment of area [m^4]."""
        return math.pi / 64.0 * (self.D ** 4 - self.d_inner ** 4)

    @property
    def r(self):
        """Radius of gyration [m]."""
        return math.sqrt(self.I / self.A)

    @property
    def S(self):
        """Elastic section modulus [m^3]."""
        return 2.0 * self.I / self.D

    @property
    def D_over_t(self):
        return self.D / self.t

    def mass_per_length(self, rho):
        """[kg/m]"""
        return self.A * rho

    def is_slender_compression(self, E, Fy):
        """AISC Table B4.1a case 9: round HSS is slender if D/t > 0.11 E / Fy."""
        return self.D_over_t > 0.11 * E / Fy

    def describe(self):
        return (f"{self.name}  A = {self.A * 1e4:.2f} cm2, "
                f"r = {self.r * 1e3:.1f} mm, D/t = {self.D_over_t:.0f}")


# Name, outside diameter [mm], wall thickness [mm] -- ordered by area so that
# the auto-sizing routine can walk the list monotonically.
_CATALOGUE_MM = [
    ("CHS 26.9 x 2.3", 26.9, 2.3),
    ("CHS 33.7 x 2.6", 33.7, 2.6),
    ("CHS 42.4 x 2.6", 42.4, 2.6),
    ("CHS 48.3 x 3.2", 48.3, 3.2),
    ("CHS 60.3 x 3.2", 60.3, 3.2),
    ("CHS 60.3 x 4.0", 60.3, 4.0),
    ("CHS 76.1 x 4.0", 76.1, 4.0),
    ("CHS 88.9 x 4.0", 88.9, 4.0),
    ("CHS 88.9 x 5.0", 88.9, 5.0),
    ("CHS 101.6 x 5.0", 101.6, 5.0),
    ("CHS 114.3 x 5.0", 114.3, 5.0),
    ("CHS 114.3 x 6.3", 114.3, 6.3),
    ("CHS 139.7 x 6.3", 139.7, 6.3),
    ("CHS 168.3 x 6.3", 168.3, 6.3),
    ("CHS 168.3 x 8.0", 168.3, 8.0),
    ("CHS 193.7 x 8.0", 193.7, 8.0),
    ("CHS 219.1 x 8.0", 219.1, 8.0),
    ("CHS 244.5 x 10.0", 244.5, 10.0),
    ("CHS 273.0 x 10.0", 273.0, 10.0),
    ("CHS 323.9 x 12.5", 323.9, 12.5),
]

CATALOGUE = [Section(n, D * 1e-3, t * 1e-3) for n, D, t in _CATALOGUE_MM]
CATALOGUE.sort(key=lambda s: s.A)

SECTION_NAMES = [s.name for s in CATALOGUE]
SECTIONS = {s.name: s for s in CATALOGUE}

DEFAULT_SECTION = "CHS 88.9 x 4.0"


def get_section(name):
    return SECTIONS[name]


def section_index(name):
    return SECTION_NAMES.index(name)


def next_larger(name):
    """Name of the next section up the catalogue, or the same name at the top."""
    i = section_index(name)
    return SECTION_NAMES[min(i + 1, len(SECTION_NAMES) - 1)]


def smallest_with(predicate, start=0):
    """First catalogue section at or after ``start`` satisfying ``predicate``."""
    for i in range(start, len(CATALOGUE)):
        if predicate(CATALOGUE[i]):
            return CATALOGUE[i]
    return None
