"""
GeoDomeLab -- Geodesic Dome Space-Truss Analyzer
================================================

A matrix-structural-analysis tool for single-layer Class I (icosahedral)
geodesic domes.

Internal unit system is strictly SI:

    length      m
    force       N
    stress      Pa
    density     kg/m^3
    mass        kg

The GUI converts to engineering display units (mm, kN, MPa, kN/m^2) only at
the user-interface boundary.

Modules
-------
geometry    Class I geodesic sphere generation, dome truncation, chord factors
sections    Circular hollow section (CHS) catalogue and section properties
materials   Structural material library
model       The DomeModel container that ties geometry + sections + loads
loads       Dead / snow / wind / point load builders and LRFD combinations
solver      3D space-truss direct stiffness solver
buckling    AISC member checks + linear (eigenvalue) buckling analysis
optimize    Fully-stressed-design auto-sizing
postprocess Result objects, utilisation and summary statistics
plotting    Matplotlib 3D renderers
report      CSV / JSON / HTML export
verification Built-in benchmark suite
nonlinear   Corotational geometric nonlinear analysis (Newton-Raphson,
            displacement control, Crisfield arc-length, imperfections)
modal       Natural frequencies, mode shapes, effective modal mass
fabrication Cut lists, hub and panel types, axial and dihedral angles
study       Parameter studies over frequency, radius and dome fraction
benchmarks  Published and closed-form benchmark models
opensees_bridge  Independent cross-check against OpenSees (optional)
gui         Qt front end
"""

__version__ = "2.0.0"
__author__ = "Jovdat Babazade"
__url__ = "https://github.com/covdetbaba/GeoDomeLab"
__all__ = [
    "geometry", "sections", "materials", "model", "loads", "solver",
    "buckling", "optimize", "postprocess", "plotting", "report",
    "verification", "nonlinear", "modal", "fabrication", "study",
    "benchmarks", "opensees_bridge",
]
