"""Structural material library (SI units).

Every material carries the elastic modulus, the specified yield and tensile
strengths, the mass density and the coefficient of linear thermal expansion
``alpha``.  ``alpha`` is used for the self-straining load T: a uniform
temperature change ``dT`` gives every strut the free (initial) strain

    eps0 = alpha * dT                      [-]

Sources for alpha
-----------------
* Carbon steel (ASTM A500, A53): 6.5e-6 per deg F = 11.7e-6 per deg C, the
  customary ambient-temperature value for structural carbon steel used in US
  practice (AISC Steel Construction Manual design tables).  AISC 360-16 itself
  gives a coefficient only for elevated temperatures (Appendix 4), and the
  EN 1993-1-1 value 12e-6 differs by less than 3 %.
* EN S235 / S355: EN 1993-1-1:2005, Section 3.2.6, alpha = 12e-6 per K.
* Aluminium 6061-T6: EN 1999-1-1:2007, Section 3.2.5, alpha = 23e-6 per K.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    name: str
    E: float          # Young's modulus              [Pa]
    Fy: float         # specified yield stress       [Pa]
    Fu: float         # specified tensile strength   [Pa]
    rho: float        # mass density                 [kg/m^3]
    # Coefficient of linear thermal expansion [1/deg C].  The default is the
    # AISC value for carbon steel so that a user-defined steel built with the
    # old five-argument call keeps a physically sensible value.
    alpha: float = 11.7e-6

    @property
    def unit_weight(self):
        """Weight density [N/m^3]."""
        return self.rho * 9.80665

    def describe(self):
        return (f"{self.name}: E = {self.E / 1e9:.0f} GPa, "
                f"Fy = {self.Fy / 1e6:.0f} MPa, rho = {self.rho:.0f} kg/m3, "
                f"alpha = {self.alpha * 1e6:.1f}e-6 /C")


ALPHA_CARBON_STEEL_AISC = 11.7e-6     # 6.5e-6 /F, AISC Manual
ALPHA_STEEL_EN = 12.0e-6              # EN 1993-1-1 3.2.6
ALPHA_ALUMINIUM_EN = 23.0e-6          # EN 1999-1-1 3.2.5

MATERIALS = {
    "ASTM A500 Gr. B (round)": Material("ASTM A500 Gr. B (round)",
                                        200e9, 317e6, 400e6, 7850.0,
                                        ALPHA_CARBON_STEEL_AISC),
    "ASTM A500 Gr. C (round)": Material("ASTM A500 Gr. C (round)",
                                        200e9, 345e6, 427e6, 7850.0,
                                        ALPHA_CARBON_STEEL_AISC),
    "ASTM A53 Gr. B pipe": Material("ASTM A53 Gr. B pipe",
                                    200e9, 240e6, 414e6, 7850.0,
                                    ALPHA_CARBON_STEEL_AISC),
    "EN S235": Material("EN S235", 210e9, 235e6, 360e6, 7850.0,
                        ALPHA_STEEL_EN),
    "EN S355": Material("EN S355", 210e9, 355e6, 490e6, 7850.0,
                        ALPHA_STEEL_EN),
    "Aluminium 6061-T6": Material("Aluminium 6061-T6",
                                  69e9, 240e6, 260e6, 2700.0,
                                  ALPHA_ALUMINIUM_EN),
}

DEFAULT_MATERIAL = "ASTM A500 Gr. C (round)"


def get_material(name):
    return MATERIALS[name]
