"""The :class:`DomeModel` container -- geometry + sections + loads + supports.

The model is deliberately a plain data object: it knows how to describe
itself and how to serialise to JSON, but all of the mathematics lives in
``solver``, ``buckling`` and ``optimize``.

Load types (ASCE 7 naming)
--------------------------
D dead, S snow, L live, W wind -- each a nodal load vector -- and T, the
**self-straining** load: a uniform temperature change ``dT`` of every strut
(initial strain alpha*dT) together with prescribed support settlements.  T is
not a load vector; ``solver.solve_model(..., T_factor=c.T)`` applies it with
the combination factor ``c.T``.

JSON format
-----------
Version 2 (written since GeoDomeLab 2.0)::

    {"format": "GeoDomeLab model", "version": 2,
     "geometry": {"kind": "geodesic", "frequency", "radius", "dome_fraction"}
              or {"kind": "custom", "nodes", "members", "faces", "base_nodes",
                  "groups", "title"},
     "material", "group_sections", "support_type",
     "spring_stiffness": [kx, ky, kz],            # N/m, ELASTIC supports
     "custom_supports": {"node": [cx, cy, cz]},   # CUSTOM supports
     "loads": {..., "temperature_change", "settlements", "nodal_loads"},
     "combinations": [{"name", "D", "S", "L", "W", "serviceability", "T"}],
     "design": {...}}

Version-1 files (no geometry kind, no T, no spring/custom supports) load
unchanged: every new field takes its default.
"""

from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass, field, asdict, fields

import numpy as np

from . import geometry as geo
from .materials import MATERIALS, DEFAULT_MATERIAL, Material
from .sections import SECTIONS, DEFAULT_SECTION, Section

FORMAT_VERSION = 2


# ------------------------------------------------------------------ loads ----
@dataclass
class LoadSpec:
    """User-facing load definition.  Pressures are in Pa, forces in N."""

    include_self_weight: bool = True
    dead_pressure: float = 300.0        # cladding + services on surface  [Pa]
    snow_pressure: float = 1200.0       # ground-projected snow           [Pa]
    snow_unbalanced: bool = False       # 1.5 p on x >= 0, 0.5 p on x < 0
    live_pressure: float = 0.0          # projected live/maintenance      [Pa]

    wind_enabled: bool = True
    wind_speed: float = 45.0            # 3 s gust basic wind speed       [m/s]
    wind_Kz: float = 0.90               # velocity pressure exposure coeff
    wind_Kzt: float = 1.00              # topographic factor
    wind_Kd: float = 0.95               # directionality factor
    wind_G: float = 0.85                # gust effect factor
    # direction the wind blows towards, degrees counter-clockwise from +x
    # (0 = wind from -x towards +x; see loads.wind_pressure_coefficients)
    wind_azimuth: float = 0.0

    crown_point_load: float = 0.0       # downward point load at crown    [N]

    # ---- self-straining load T (temperature + settlement) -----------------
    # Uniform temperature change of every strut [deg C], + = heating.
    temperature_change: float = 0.0
    # Prescribed movement of supported nodes [m]:
    #   [{"node": int, "ux": m, "uy": m, "uz": m}, ...]   (missing keys = 0)
    # Applies to restrained directions (and moves the ground end of a spring
    # support); a non-zero value in a free direction raises ValueError.
    settlements: list = field(default_factory=list)

    # ---- explicit nodal loads ------------------------------------------------
    #   [{"node": int, "Fx": N, "Fy": N, "Fz": N, "case": "D"|"S"|"L"|"W"}]
    # added to the load vector of their case ("L" when "case" is missing).
    nodal_loads: list = field(default_factory=list)

    def copy(self):
        """Deep copy (the settlement and nodal-load lists are copied too)."""
        return copy.deepcopy(self)

    @property
    def has_self_straining(self):
        return (abs(self.temperature_change) > 0.0
                or any(any(float(s.get(k, 0.0) or 0.0) != 0.0
                           for k in ("ux", "uy", "uz"))
                       for s in self.settlements))

    @staticmethod
    def from_dict(d):
        """Build from a dict, ignoring unknown keys (forward compatible)."""
        names = {f.name for f in fields(LoadSpec)}
        kw = {k: copy.deepcopy(v) for k, v in (d or {}).items() if k in names}
        return LoadSpec(**kw)


LOAD_KEYS = ("D", "S", "L", "W", "T")


@dataclass
class LoadCombination:
    """An LRFD load combination: name and the factor on each load type.

    ``T`` multiplies the self-straining load (temperature + settlement).
    It is the last field so that old positional calls keep working.
    """

    name: str
    D: float = 0.0
    S: float = 0.0
    L: float = 0.0
    W: float = 0.0
    serviceability: bool = False        # deflection combos are unfactored
    T: float = 0.0

    def factors(self):
        return {"D": self.D, "S": self.S, "L": self.L, "W": self.W,
                "T": self.T}

    def label(self):
        parts = []
        for k, v in self.factors().items():
            if abs(v) > 1e-12:
                parts.append(f"{v:g}{k}")
        return " + ".join(parts) if parts else "0"

    @staticmethod
    def from_dict(d):
        names = {f.name for f in fields(LoadCombination)}
        return LoadCombination(**{k: v for k, v in d.items() if k in names})


def default_combinations():
    """ASCE 7-16 Sec. 2.3.1 strength combinations plus one service case.

    Snow S stands for "Lr or S or R" (a roof).  The set covers basic
    combinations 1 (LC1), 2 (LC6, with the live load L), 3 (LC2 and LC5,
    the "L or 0.5W" term taken as 0.5W) and 4 and 5 (LC3, LC4) for D, S
    and W.  When L is not zero, ASCE 7 also adds L to combinations 3 and 4
    (1.2D + 1.6S + 1.0L, 1.2D + 1.0W + 1.0L + 0.5S): add those rows by hand.
    ``SLS`` is unfactored and used for the deflection check.
    """
    return [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("LC3  1.2D + 1.0W + 0.5S", D=1.2, W=1.0, S=0.5),
        LoadCombination("LC4  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("LC5  1.2D + 1.6S + 0.5W", D=1.2, S=1.6, W=0.5),
        LoadCombination("LC6  1.2D + 1.6L + 0.5S", D=1.2, L=1.6, S=0.5),
        LoadCombination("SLS  1.0D + 1.0S", D=1.0, S=1.0, serviceability=True),
    ]


# ------------------------------------------------------------- json help ----
def _jsonable(obj):
    """Convert numpy scalars/arrays and non-finite floats for strict JSON."""
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj)
    return obj


def _support_entry_to_json(v):
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if v is None:
        return False
    if isinstance(v, str):
        return v
    v = float(v)
    if math.isinf(v) and v > 0:
        return True                      # an infinitely stiff spring is fixed
    return v


def _support_entry_from_json(v):
    if isinstance(v, bool) or v is None or isinstance(v, str):
        return False if v is None else v
    return float(v)


# ------------------------------------------------------------------ model ----
@dataclass
class DomeModel:
    """Complete analysis model for one geodesic dome (or custom truss)."""

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
    # wall-thickness factor for the strength checks only (AISC 360-16 B4.2:
    # 0.93 for ERW HSS such as ASTM A500; 1.0 = nominal thickness)
    design_thickness_factor: float = 1.0

    # support data for the ELASTIC and CUSTOM support types
    # (kx, ky, kz) [N/m] of the spring at every base node (ELASTIC)
    spring_stiffness: tuple = geo.DEFAULT_SPRING_STIFFNESS
    # {node: [cx, cy, cz]}: True = fixed, False = free, k > 0 = spring [N/m]
    custom_supports: dict = field(default_factory=dict)

    def __post_init__(self):
        n = self.geom.n_strut_types()
        if not self.group_sections:
            self.group_sections = [DEFAULT_SECTION] * n
        elif len(self.group_sections) != n:
            fill = self.group_sections[-1] if self.group_sections else DEFAULT_SECTION
            self.group_sections = (self.group_sections + [fill] * n)[:n]
        self.spring_stiffness = tuple(float(k) for k in self.spring_stiffness)
        self.custom_supports = {int(k): list(v)
                                for k, v in (self.custom_supports or {}).items()}

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

    # ------------------------------------------------------------ supports --
    def support_conditions(self):
        """``(fixed_dofs, springs)`` -- rigid restraints and {dof: k} springs.

        See :func:`geometry.support_conditions` for the four support types.
        """
        return geo.support_conditions(self.geom, self.support_type,
                                      self.spring_stiffness,
                                      self.custom_supports)

    def fixed_dofs(self):
        """Rigidly restrained DOFs (the fixed part of support_conditions)."""
        return self.support_conditions()[0]

    def support_nodes(self):
        """Sorted node numbers that carry a restraint or a spring."""
        fixed, springs = self.support_conditions()
        dofs = np.concatenate([np.asarray(fixed, dtype=int),
                               np.asarray(sorted(springs), dtype=int)])
        return np.unique(dofs // 3)

    # ------------------------------------------------------------- weights --
    def steel_mass(self):
        """Total structural mass [kg]."""
        rho = self.material.rho
        L = self.geom.member_lengths()
        A = self.member_areas()
        return float(np.sum(A * L * rho))

    def mass_per_plan_area(self):
        """Mass per unit plan area [kg/m^2]; NaN when there are no panels."""
        pa = self.geom.plan_area
        if not pa > 0.0:
            return float("nan")
        return self.steel_mass() / max(pa, 1e-12)

    # --------------------------------------------------------- description --
    def _per_plan_text(self):
        per_area = self.mass_per_plan_area()
        if not math.isfinite(per_area):
            return "(no panels, so no plan area)"
        return f"({per_area:.1f} kg/m2 of plan)"

    def _support_line(self):
        text = geo.SUPPORT_TYPES.get(self.support_type, self.support_type)
        if self.support_type == geo.ELASTIC:
            k = ", ".join(f"{v / 1e6:g}" for v in self.spring_stiffness)
            text += f" (kx, ky, kz = {k} kN/mm)"
        elif self.support_type == geo.CUSTOM:
            text += f" ({len(self.custom_supports)} node(s) listed)"
        return text

    def summary_lines(self):
        g = self.geom
        if g.is_custom:
            R = g.radius
            sup = (f"{len(self.support_nodes())}"
                   if self.support_type != geo.CUSTOM or self.custom_supports
                   else "none defined")
            return [
                f"Structure              {g.title or 'Custom truss'} "
                f"(custom space truss)",
                f"Best-fit sphere radius "
                + (f"{R:.3f} m" if np.isfinite(R) else "n/a (not spherical)"),
                f"Span (support plan)    {g.base_diameter:.3f} m",
                f"Rise                   {g.height:.3f} m",
                f"Nodes / struts / panels {g.n_nodes} / {g.n_members} / {g.n_faces}",
                f"Strut families         {g.n_strut_types()}",
                f"Supported nodes        {sup}",
                f"Support condition      {self._support_line()}",
                f"Material               {self.material.describe()}",
                f"Faceted surface area   {g.surface_area:.1f} m2",
                f"Plan (footprint) area  {g.plan_area:.1f} m2",
                f"Structural steel mass  {self.steel_mass():.0f} kg "
                f"{self._per_plan_text()}",
            ]
        if self.support_type in (geo.PINNED, geo.ROLLER3):
            sup_line = (f"Supported base nodes   {len(g.base_nodes)}  "
                        f"({g.frequency}V ring)")
        else:
            try:
                ns = len(self.support_nodes())
            except ValueError:
                ns = 0
            sup_line = f"Supported nodes        {ns}"
        return [
            f"Frequency              {g.frequency}V  (Class I, icosahedral)",
            f"Sphere radius          {g.radius:.3f} m",
            f"Dome fraction          {g.dome_fraction:.3f} of diameter",
            f"Base diameter          {g.base_diameter:.3f} m",
            f"Crown height           {g.height:.3f} m",
            f"Nodes / struts / panels {g.n_nodes} / {g.n_members} / {g.n_faces}",
            f"Strut families         {g.n_strut_types()}",
            sup_line,
            f"Support condition      {self._support_line()}",
            f"Material               {self.material.describe()}",
            f"Faceted surface area   {g.surface_area:.1f} m2",
            f"Plan (footprint) area  {g.plan_area:.1f} m2",
            f"Structural steel mass  {self.steel_mass():.0f} kg "
            f"{self._per_plan_text()}",
        ]

    # ------------------------------------------------------ serialisation --
    def _geometry_dict(self):
        g = self.geom
        if g.is_custom:
            return {
                "kind": "custom",
                "title": g.title,
                "nodes": [[float(v) for v in row] for row in g.nodes],
                "members": [[int(v) for v in row] for row in g.members],
                "faces": [[int(v) for v in row] for row in g.faces],
                "base_nodes": [int(v) for v in g.base_nodes],
                "groups": [int(v) for v in g.strut_group],
            }
        return {
            "kind": "geodesic",
            "frequency": int(g.frequency),
            "radius": float(g.radius),
            "dome_fraction": float(g.dome_fraction),
        }

    def to_dict(self):
        return {
            "format": "GeoDomeLab model",
            "version": FORMAT_VERSION,
            "geometry": self._geometry_dict(),
            "material": self.material_name,
            "group_sections": list(self.group_sections),
            "support_type": self.support_type,
            "spring_stiffness": [float(k) for k in self.spring_stiffness],
            "custom_supports": {
                str(int(k)): [_support_entry_to_json(c) for c in v]
                for k, v in self.custom_supports.items()},
            "loads": _jsonable(asdict(self.loads)),
            "combinations": [_jsonable(asdict(c)) for c in self.combinations],
            "design": {
                "K_effective": self.K_effective,
                "deflection_limit": self.deflection_limit,
                "phi_c": self.phi_c,
                "phi_t": self.phi_t,
                "design_thickness_factor": self.design_thickness_factor,
            },
        }

    def save(self, path):
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @staticmethod
    def geometry_from_dict(gd):
        """Rebuild a :class:`DomeGeometry` from the JSON ``geometry`` block."""
        kind = gd.get("kind", "geodesic")
        if kind == "geodesic":
            return geo.generate_dome(gd["frequency"], gd["radius"],
                                     gd["dome_fraction"])
        if kind == "custom":
            faces = gd.get("faces") or None
            return geo.custom_geometry(
                gd["nodes"], gd["members"], faces=faces,
                base_nodes=gd.get("base_nodes"), groups=gd.get("groups"),
                title=gd.get("title", "Custom truss"))
        raise ValueError(f"Unknown geometry kind {kind!r}")

    @staticmethod
    def from_dict(d):
        version = int(d.get("version", 1))
        if version > FORMAT_VERSION:
            raise ValueError(f"Model file version {version} is newer than "
                             f"this GeoDomeLab ({FORMAT_VERSION}).")
        geom = DomeModel.geometry_from_dict(d["geometry"])
        m = DomeModel(
            geom=geom,
            material_name=d.get("material", DEFAULT_MATERIAL),
            group_sections=list(d.get("group_sections", [])),
            support_type=d.get("support_type", geo.PINNED),
            loads=LoadSpec.from_dict(d.get("loads", {})),
            combinations=[LoadCombination.from_dict(c)
                          for c in d.get("combinations", [])]
            or default_combinations(),
        )
        des = d.get("design", {})
        m.K_effective = des.get("K_effective", 1.0)
        m.deflection_limit = des.get("deflection_limit", 240.0)
        m.phi_c = des.get("phi_c", 0.90)
        m.phi_t = des.get("phi_t", 0.90)
        m.design_thickness_factor = float(
            des.get("design_thickness_factor", 1.0))
        if "spring_stiffness" in d:
            m.spring_stiffness = tuple(float(k) for k in d["spring_stiffness"])
        m.custom_supports = {
            int(k): [_support_entry_from_json(c) for c in v]
            for k, v in (d.get("custom_supports") or {}).items()}
        return m

    @staticmethod
    def load(path):
        with open(path, "r", encoding="utf-8") as fh:
            return DomeModel.from_dict(json.load(fh))


def default_model(frequency=4, radius=12.0, dome_fraction=0.5):
    """Convenience factory used by the GUI on start-up and by the tests."""
    return DomeModel(geom=geo.generate_dome(frequency, radius, dome_fraction))
