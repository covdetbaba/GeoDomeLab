"""Structural material library (SI units)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    name: str
    E: float          # Young's modulus              [Pa]
    Fy: float         # specified yield stress       [Pa]
    Fu: float         # specified tensile strength   [Pa]
    rho: float        # mass density                 [kg/m^3]

    @property
    def unit_weight(self):
        """Weight density [N/m^3]."""
        return self.rho * 9.80665

    def describe(self):
        return (f"{self.name}: E = {self.E / 1e9:.0f} GPa, "
                f"Fy = {self.Fy / 1e6:.0f} MPa, rho = {self.rho:.0f} kg/m3")


MATERIALS = {
    "ASTM A500 Gr. B (round)": Material("ASTM A500 Gr. B (round)",
                                        200e9, 317e6, 400e6, 7850.0),
    "ASTM A500 Gr. C (round)": Material("ASTM A500 Gr. C (round)",
                                        200e9, 345e6, 427e6, 7850.0),
    "ASTM A53 Gr. B pipe": Material("ASTM A53 Gr. B pipe",
                                    200e9, 240e6, 414e6, 7850.0),
    "EN S235": Material("EN S235", 210e9, 235e6, 360e6, 7850.0),
    "EN S355": Material("EN S355", 210e9, 355e6, 490e6, 7850.0),
    "Aluminium 6061-T6": Material("Aluminium 6061-T6",
                                  69e9, 240e6, 260e6, 2700.0),
}

DEFAULT_MATERIAL = "ASTM A500 Gr. C (round)"


def get_material(name):
    return MATERIALS[name]
