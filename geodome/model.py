"""The :class:`DomeModel` container -- geometry + sections + loads + supports.

The model is deliberately a plain data object: it knows how to describe
itself and how to serialise to JSON, but all of the mathematics lives in
``solver``, ``buckling`` and ``optimize``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict

import numpy as np

from . import geometry as geo
from .materials import MATERIALS, DEFAULT_MATERIAL, Material
from .sections import SECTIONS, DEFAULT_SECTION, Section


# ------------------------------------------------------------------ loads ----
@dataclass
class LoadSpec:
    """User-facing load definition.  Pressures are in Pa, forces in N."""

    include_self_weight: bool = True
    dead_pressure: float = 300.0        # cladding + services on surface  [Pa]
    snow_pressure: float = 1200.0       # ground-projected snow           [Pa]
    snow_unbalanced: bool = False       # 1.5 / 0.5 split across the dome
    live_pressure: float = 0.0          # projected live/maintenance      [Pa]

    wind_enabled: bool = True
    wind_speed: float = 45.0            # 3 s gust basic wind speed       [m/s]
    wind_Kz: float = 0.90               # velocity pressure exposure coeff
    wind_Kzt: float = 1.00              # topographic factor
    wind_Kd: float = 0.95               # directionality factor
    wind_G: float = 0.85                # gust effect factor
    wind_azimuth: float = 0.0           # wind comes from this bearing [deg]

    crown_point_load: float = 0.0       # downward point load at crown    [N]

    def copy(self):
        return LoadSpec(**asdict(self))


@dataclass
class LoadCombination:
    """An LRFD load combination: name and the factor on each load type."""

    name: str
    D: float = 0.0
    S: float = 0.0
    L: float = 0.0
    W: float = 0.0
    serviceability: bool = False        # deflection combos are unfactored

    def factors(self):
        return {"D": self.D, "S": self.S, "L": self.L, "W": self.W}

    def label(self):
        parts = []
        for k, v in self.factors().items():
            if abs(v) > 1e-12:
                parts.append(f"{v:g}{k}")
        return " + ".join(parts) if parts else "0"


def default_combinations():
    """ASCE 7 strength combinations (snow-governed roof) plus one service case."""
    return [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("LC3  1.2D + 1.0W + 0.5S", D=1.2, W=1.0, S=0.5),
        LoadCombination("LC4  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("SLS  1.0D + 1.0S", D=1.0, S=1.0, serviceability=True),
    ]


# ------------------------------------------------------------------ model ----
@dataclass
class DomeModel:
    """Complete analysis model for one geodesic dome."""

    geom: geo.DomeGeometry
    material_name: str = DEFAULT_MATERIAL
    # One section name per strut *family* (A, B, C, ...).  Assigning by family
    # rather than by individual member is how real geodesic domes are built and
    # it keeps the fabrication schedule short.
    group_sections: list = field(default_factory=list)
    support_type: str = geo.PINNED
    loads: LoadSpec = field(default_factory=LoadSpec)
    combinations: list = field(default_factory=default_combinations)

    # design parameters
    K_effective: float = 1.0            # effective length factor
    deflection_limit: float = 240.0     # span / this value
    phi_c: float = 0.90                 # LRFD resistance factor, compression
    phi_t: float = 0.90                 # LRFD resistance factor, tension

    def __post_init__(self):
        n = self.geom.n_strut_types()
        if not self.group_sections:
            self.group_sections = [DEFAULT_SECTION] * n
        elif len(self.group_sections) != n:
            fill = self.group_sections[-1] if self.group_sections else DEFAULT_SECTION
            self.group_sections = (self.group_sections + [fill] * n)[:n]

    # ---------------------------------------------------------- accessors --
    @property
    def material(self) -> Material:
        return MATERIALS[self.material_name]

    def section_of_group(self, g) -> Section:
        return SECTIONS[self.group_sections[g]]

    def member_sections(self):
        """List of :class:`Section` objects, one per member."""
        return [SECTIONS[self.group_sections[g]] for g in self.geom.strut_group]

    def member_areas(self):
        return np.array([s.A for s in self.member_sections()])

    def member_radii_gyration(self):
        return np.array([s.r for s in self.member_sections()])

    def EA(self):
        return self.material.E * self.member_areas()

    def fixed_dofs(self):
        return geo.support_dofs(self.geom, self.support_type)

    # ------------------------------------------------------------- weights --
    def steel_mass(self):
        """Total structural mass [kg]."""
        rho = self.material.rho
        L = self.geom.member_lengths()
        A = self.member_areas()
        return float(np.sum(A * L * rho))

    def mass_per_plan_area(self):
        return self.steel_mass() / max(self.geom.plan_area, 1e-12)

    # --------------------------------------------------------- description --
    def summary_lines(self):
        g = self.geom
        return [
            f"Frequency              {g.frequency}V  (Class I, icosahedral)",
            f"Sphere radius          {g.radius:.3f} m",
            f"Dome fraction          {g.dome_fraction:.3f} of diameter",
            f"Base diameter          {g.base_diameter:.3f} m",
            f"Crown height           {g.height:.3f} m",
            f"Nodes / struts / panels {g.n_nodes} / {g.n_members} / {g.n_faces}",
            f"Strut families         {g.n_strut_types()}",
            f"Supported base nodes   {len(g.base_nodes)}  ({g.frequency}V ring)",
            f"Support condition      {geo.SUPPORT_TYPES[self.support_type]}",
            f"Material               {self.material.describe()}",
            f"Faceted surface area   {g.surface_area:.1f} m2",
            f"Plan (footprint) area  {g.plan_area:.1f} m2",
            f"Structural steel mass  {self.steel_mass():.0f} kg "
            f"({self.mass_per_plan_area():.1f} kg/m2 of plan)",
        ]

    # ------------------------------------------------------ serialisation --
    def to_dict(self):
        return {
            "format": "GeoDomeLab model",
            "version": 1,
            "geometry": {
                "frequency": int(self.geom.frequency),
                "radius": float(self.geom.radius),
                "dome_fraction": float(self.geom.dome_fraction),
            },
            "material": self.material_name,
            "group_sections": list(self.group_sections),
            "support_type": self.support_type,
            "loads": asdict(self.loads),
            "combinations": [asdict(c) for c in self.combinations],
            "design": {
                "K_effective": self.K_effective,
                "deflection_limit": self.deflection_limit,
                "phi_c": self.phi_c,
                "phi_t": self.phi_t,
            },
        }

    def save(self, path):
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @staticmethod
    def from_dict(d):
        gd = d["geometry"]
        geom = geo.generate_dome(gd["frequency"], gd["radius"],
                                 gd["dome_fraction"])
        m = DomeModel(
            geom=geom,
            material_name=d.get("material", DEFAULT_MATERIAL),
            group_sections=list(d.get("group_sections", [])),
            support_type=d.get("support_type", geo.PINNED),
            loads=LoadSpec(**d.get("loads", {})),
            combinations=[LoadCombination(**c) for c in d.get("combinations", [])]
            or default_combinations(),
        )
        des = d.get("design", {})
        m.K_effective = des.get("K_effective", 1.0)
        m.deflection_limit = des.get("deflection_limit", 240.0)
        m.phi_c = des.get("phi_c", 0.90)
        m.phi_t = des.get("phi_t", 0.90)
        return m

    @staticmethod
    def load(path):
        with open(path, "r", encoding="utf-8") as fh:
            return DomeModel.from_dict(json.load(fh))


def default_model(frequency=4, radius=12.0, dome_fraction=0.5):
    """Convenience factory used by the GUI on start-up and by the tests."""
    return DomeModel(geom=geo.generate_dome(frequency, radius, dome_fraction))
