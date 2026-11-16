"""
Benchmark and validation presets -- ready-made :class:`DomeModel` objects.

Every preset is a complete model (geometry, sections, material, supports,
loads, combinations) that opens in the GUI, runs through the ordinary
analysis pipeline and has an *independent* reference answer, so anyone can
check the program on a structure whose solution is known before trusting it
on their own dome.  ``verify_preset(key)`` runs the preset and compares the
program with that reference; ``python -m geodome.benchmarks`` does it for all.

Registry
--------
``PRESETS = {key: (title, description, factory)}``; ``make_preset(key)``
returns a fresh model and ``write_preset_files()`` saves every preset to
``examples/models/<key>.json`` (File > Open in the GUI).

=================  ==========================================================
key                structure / purpose
=================  ==========================================================
validation_a       Validation Model A: 1V pentagonal cap (6 nodes, 10 struts)
validation_b       Validation Model B: 1V dome on its lower ring (11 / 25)
von_mises          von Mises two-bar truss (snap-through, closed form)
star_dome          24-member shallow star dome (snap-through benchmark)
three_bar          planar three-bar indeterminate truss (temperature and
                   settlement, force method)
case_study_6v      the 6V 15 m case study of ``examples/worked_example.py``
=================  ==========================================================

Validation Models A and B
-------------------------
Built exactly as ``validation/run_validation.py`` builds them (same
geometry ``generate_dome(1, 5 m, fraction)``, CHS 88.9 x 4.0 in ASTM A500
Gr. C, pinned base, K = 1.0, phi_c = phi_t = 0.90, cladding 300 Pa, snow
1200 Pa, wind off).  The hand side is ``validation/hand_calc.py``: Model A
by the direct stiffness method written out by hand (3 free DOF), Model B by
virtual work with the five-fold symmetry (18 DOF reduced to 3); see McGuire,
Gallagher & Ziemian (2000), *Matrix Structural Analysis*, 2nd ed., Ch. 3.
Each preset carries both validation load cases in one model:

* A1 / B1 -- the unfactored hand-check point load at the crown (A1:
  Fx = +20 kN, Fz = -60 kN as an explicit nodal load; B1: 50 kN downward as
  the crown point load), combination ``1.0L``.  It is flagged
  ``serviceability=True`` (it is an unfactored case), so the strength
  envelope, the governing case and the buckling reference stay on A2 / B2
  exactly as in the validation run.  The L case holds nothing else.
* A2 / B2 -- ``1.2D + 1.6S`` from self weight, cladding and snow.

Model A therefore reproduces every row of ``validation/output/comparison.csv``
including the AISC utilisation 0.27257 and lambda_cr = 5506.10 of A2.

von Mises two-bar truss (Pecknold, Ghaboussi & Healey 1985)
-----------------------------------------------------------
Two bars from the supports (-a, 0, 0), (+a, 0, 0) to the apex (0, 0, h),
same section as validation case F6 of ``validation/hand_calc_ext.py``
(a = 2.5 m, h = 0.125 m, CHS 88.9 x 4.0, reference load P_ref = 10 kN down
at the apex).  The truss lies in the x-z plane; CUSTOM supports pin both
supports and restrain the apex out of plane (uy), so the 3D program sees a
properly restrained planar structure with two free DOF (apex ux, uz).

With L = sqrt(a^2 + h^2), apex deflection w (down), y = h - w,
l = sqrt(a^2 + y^2) and engineering strain (l - L)/L, the exact equilibrium
path, its limit point and the linear results are::

    P(w)   = 2 EA y (1/l - 1/L)                    (downward apex load)
    l*     = (a^2 L)^(1/3),  y* = sqrt(l*^2 - a^2),  w* = h - y*
    P*     = 2 EA y* (1/l* - 1/L)                  (first limit point)
    K0     = dP/dw(0) = 2 EA h^2 / L^3,  w_lin = P / K0,  N_lin = -P L/(2h)
    lambda_cr P_ref = 2 EA h^3 / (L a^2)           (linearised buckling)

P = 0 at w = h (bars level, N most compressive) and at w = 2h (mirror
state, N = 0).  The linearised lambda_cr overestimates the true limit load
(by 5.2x for these proportions) -- the reason a nonlinear path is needed.

24-member shallow star dome (Hangai & Kawamata 1972)
-----------------------------------------------------
Geometry from Greco & Venturini (2006), *Lat. Am. J. Solids Struct.* 3(3),
Sec. 5.3 / Fig. 7, and Iu & Bradford (2012), *Adv. Steel Constr.* 8(2),
Sec. 5.5 / Fig. 8; the plan is also confirmed by Hrinda (2006, NASA NTRS
20060020675), Fig. 6 (all compiled in ``Final Project/dev/research_notes.md``,
Sec. 5c).  In cm: crown (0, 0, 8.216); inner ring of six nodes at radius 25
and z = 6.216 (one on the +x axis); six pinned supports at radius 50 and
z = 0, at 30 deg + k 60 deg.  The outer plan coordinate is the exact
50 cos 30 deg = 25 sqrt(3) = 43.30127 cm (the papers print 43.3, which moves
the first limit load by 2.5e-6, relative).  Coordinates are converted to
metres exactly (x 0.01), so a crown deflection of 0.01 m is 1 cm in the
papers.  All 24 bars have the same EA (any value: P/EA is dimensionless and
independent of the length unit).  The reference load is P_ref = 1e-4 EA
down at the crown, so a load factor lambda reads directly as 10^4 P/EA,
the ordinate of Iu & Bradford's plots.

Published and reference values (research notes, Sec. 5d):

* first limit point "near 0.75 cm" of crown deflection (Greco & Venturini);
  P/EA = 3.15655e-4 at w = 0.768 cm (corotational engineering strain;
  0.76844 cm when the flat maximum is refined);
* zero-load states: crown 19 mm with node 2 at -1.16 mm, i.e. moving up
  (Iu & Bradford: "approximately"), and crown 40 mm with node 2 at 0;
* exact mirror states (every bar back at its original length, N = 0):
  crown inverted through the ring plane w = 2 x 2.000 cm, whole dome
  inverted w = 2 x 8.216 cm;
* initial crown stiffness 0.090430 EA per m; linearised buckling
  lambda_cr P/EA = 1.52097e-3, about 4.8 times the true first limit load.

The linear symmetric response used by ``verify_preset`` is solved
independently of the global solver: under the crown load the six-fold
symmetric displacement field has three generalised coordinates
q = (w0 crown vertical, u_r ring radial, w1 ring vertical).  Each bar's
elongation is e = g.q with (spoke, ring, diagonal; d1 = crown height above
the ring, h_r = ring height above the supports, a, b = ring and support
radii)::

    g_spoke = (d1, a, -d1) / L1
    g_ring  = (0, 2 sin 30deg, 0)
    g_diag  = (0, a - b cos 30deg, h_r) / L3

and virtual work gives  [6 EA g_s g_s^T / L1 + 6 EA g_r g_r^T / L2 +
12 EA g_d g_d^T / L3] q = (-P, 0, 0),  N = EA (g.q) / L.

Planar three-bar indeterminate truss (force method)
---------------------------------------------------
Validation case F4 of ``validation/hand_calc_ext.py``: supports A (-4, 0, 0),
B (0, 0, 0), C (3, 0, 0), joint D (0, 0, 3) m; bars A-D and C-D CHS 88.9 x
4.0, bar B-D CHS 76.1 x 4.0; ASTM A500 Gr. C (alpha = 11.7e-6 /degC).
Load (15, 0, -60) kN at D; self-straining load T = uniform dT = +30 degC on
all three bars plus a 5 mm settlement of support C.  CUSTOM supports pin
A, B, C and restrain D out of plane.  (F4 heats bar B-D only; GeoDomeLab's
temperature load is uniform, so the preset heats every bar.  The force-
method reference below takes any per-bar dT and reproduces F4 as well.)

The reference is the force (flexibility) method -- unknowns are forces,
nothing is shared with the stiffness method of the solver.  Redundant
X1 = N2 (bar B-D); primary structure = bars 1 and 3; c_i = unit vector from
the support of bar i to D (x, z components)::

    sum_i N_i c_i = P                         (equilibrium of joint D)
    n0 = primary forces under P,   n1 = primary forces under X1 = 1
    f11 = sum n1_i^2 L_i/EA_i
    D1P = sum n1_i n0_i L_i/EA_i,  D1T = sum n1_i alpha dT_i L_i
    D1S = -sum_k Rbar_k . delta_k,  Rbar_k = -n1_i c_i  (support of bar i)
    X1  = -(D1P + D1T + D1S) / f11,   N_i = n0_i + X1 n1_i
    e_i = N_i L_i/EA_i + alpha dT_i L_i
    u_D . e = sum_{1,3} n_e,i e_i - sum_k Rtilde_e,k . delta_k   (unit load)
    R_k = -N_i c_i

6V case study
-------------
The 30 m span, 15 m rise exhibition-hall dome of ``examples/worked_example.py``
(trial section CHS 114.3 x 5.0 before auto-sizing).  It has no closed-form
answer; ``verify_preset`` checks the statics instead: the vertical reaction
of every gravity combination equals the dead and snow load computed directly
from the bar lengths and the true / projected panel areas (not from the
nodal tributary areas the load generator uses), plus the equilibrium and
energy identities.

References (verified, see ``Final Project/dev/research_notes.md``)
-------------------------------------------------------------------
* D. A. Pecknold, J. Ghaboussi, T. J. Healey (1985). Snap-through and
  bifurcation in a simple structure. J. Eng. Mech. 111(7), 909-922.
* Y. Hangai, S. Kawamata (1972). Nonlinear analysis of space frames and
  snap-through buckling of reticulated shell structures. Proc. IASS Pacific
  Symp. (1971), AIJ Tokyo, 803-816 (original star dome; cited, not seen).
* M. Greco, W. S. Venturini (2006). Stability analysis of three-dimensional
  trusses. Lat. Am. J. Solids Struct. 3(3), 325-344.
* C. K. Iu, M. A. Bradford (2012). Higher-order non-linear analysis of steel
  structures, Part I. Adv. Steel Constr. 8(2), 168-182.
* G. A. Hrinda (2006). Geometrically nonlinear static analysis of 3D trusses
  using the arc-length method. HPSM 2006 (NASA NTRS 20060020675).
* W. McGuire, R. H. Gallagher, R. D. Ziemian (2000). Matrix Structural
  Analysis, 2nd ed. Wiley (Ch. 3 direct stiffness, Ch. 9 critical load).
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

from . import geometry as geo
from .model import DomeModel, LoadSpec, LoadCombination
from .sections import SECTIONS
from .units import G_ACCEL

SQ5 = math.sqrt(5.0)
SQ3 = math.sqrt(3.0)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# examples/models of the source tree; when the package is installed without
# the examples folder, fall back to ./examples/models of the working folder.
MODELS_DIR = (os.path.join(_ROOT, "examples", "models")
              if os.path.isdir(os.path.join(_ROOT, "examples"))
              else os.path.join(os.getcwd(), "examples", "models"))


# =============================================================================
# shared data of the validation models  (= validation/hand_calc.py)
# =============================================================================
VAL_RADIUS = 5.0                          # sphere radius               [m]
VAL_MATERIAL = "ASTM A500 Gr. C (round)"
VAL_SECTION = "CHS 88.9 x 4.0"
VAL_K_EFF = 1.0
VAL_PHI_C = 0.90
VAL_PHI_T = 0.90
VAL_DEAD_PRESSURE = 300.0                 # cladding on the surface      [Pa]
VAL_SNOW_PRESSURE = 1200.0                # snow, projected              [Pa]
VAL_A1_LOAD = (20e3, 0.0, -60e3)          # crown load of case A1  (Fx, Fy, Fz) [N]
VAL_B1_LOAD = 50e3                        # crown load of case B1, downward     [N]
FRACTION_A = 0.5 * (1.0 - 1.0 / SQ5)      # cut at the upper ring of the 1V sphere
FRACTION_B = 0.5 * (1.0 + 1.0 / SQ5)      # cut at the lower ring

COMBO_A1 = "A1  1.0L"
COMBO_A2 = "A2  1.2D + 1.6S"
COMBO_B1 = "B1  1.0L"
COMBO_B2 = "B2  1.2D + 1.6S"

# Hand-calculation answers (validation/hand_calc.py, SI units).  The test
# suite recomputes them from hand_calc.py and checks these to 1e-12.
HAND_A = {
    "A1_ux": 0.00027239772602717297, "A1_uz": -0.0010697197577959816,
    "A1_N1": -32229.920427763263, "A1_sumRz": 60000.000000000015,
    "A2_uz": -0.0006937147165418648, "A2_N_spoke": -14802.274636332933,
    "A2_Rz_node": 23605.122314028544, "A2_total_Pz": -118025.61157014272,
    "A2_util": 0.27257000721147456, "A2_lambda_cr": 5506.096413150748,
}
HAND_B = {
    "B1": {"w0": -0.001645132269502631, "ur": 0.00032252353684526373,
           "w1": -0.00023184509321851396, "N_spoke": -19021.130325903076,
           "N_ring": 15388.417685876268, "N_zig": -5877.8525229247325,
           "Rz_node": 10000.000000000002, "Rr_node": -1909.830056250526,
           "tot_P": -50000.0},
    "B2": {"w0": -0.0018589407480980834, "ur": 0.00033597801610817655,
           "w1": -0.0006216021820204292, "N_spoke": -14802.274636332933,
           "N_ring": 16030.365088130498, "N_zig": -19245.547632195747,
           "Rz_node": 42397.98707733561, "Rr_node": -6253.257490488971,
           "tot_P": -211989.93538667803},
}

# F-case data shared with validation/hand_calc_ext.py
F4_SECTIONS = ("CHS 88.9 x 4.0", "CHS 76.1 x 4.0", "CHS 88.9 x 4.0")
F4_LOAD = (15e3, 0.0, -60e3)              # at joint D  (Fx, Fy, Fz)     [N]
F4_DT = 30.0                              # uniform heating            [degC]
F4_SETTLEMENT = -0.005                    # support C vertical         [m]
VM_HALF_SPAN = 2.5                        # von Mises a                [m]
VM_RISE = 0.125                           # von Mises h                [m]
VM_P_REF = 10e3                           # apex reference load        [N]

# Star dome reference numbers (research_notes.md Sec. 5d)
STAR_P_OVER_EA_LIMIT = 3.15655e-4         # first limit point, P/EA
# Crown deflection at that limit point [cm].  The research notes print
# 0.7685 (a rounded grid value); refining the flat maximum by golden-section
# search on the corotational path gives 0.76844 (tests/test_benchmarks.py).
STAR_W_LIMIT_CM = 0.76844
STAR_W_ZERO1_CM = 1.8838                  # first zero-load state       [cm]
STAR_NODE2_ZERO1_CM = 0.1162              # node 2 moves UP by this     [cm]
STAR_LAMBDA_CR_P_OVER_EA = 1.52097e-3     # linearised buckling, P/EA
STAR_K0_PER_EA = 0.090430                 # initial crown stiffness / EA [1/m]
STAR_LOAD_RATIO = 1e-4                    # P_ref / EA of the preset


# =============================================================================
# helpers
# =============================================================================
def _no_area_loads(**kw):
    """A LoadSpec with every area load, wind and self weight switched off."""
    base = dict(include_self_weight=False, dead_pressure=0.0,
                snow_pressure=0.0, snow_unbalanced=False, live_pressure=0.0,
                wind_enabled=False, crown_point_load=0.0)
    base.update(kw)
    return LoadSpec(**base)


def _validation_lc2_loads():
    """The loads of validation case 2 (run_validation.lc2_loads)."""
    return LoadSpec(include_self_weight=True, dead_pressure=VAL_DEAD_PRESSURE,
                    snow_pressure=VAL_SNOW_PRESSURE, snow_unbalanced=False,
                    live_pressure=0.0, wind_enabled=False,
                    crown_point_load=0.0)


def _validation_model(fraction, loads, combos):
    """Same construction as run_validation.program_model."""
    g = geo.generate_dome(1, VAL_RADIUS, fraction)
    m = DomeModel(geom=g, material_name=VAL_MATERIAL,
                  group_sections=[VAL_SECTION] * g.n_strut_types(),
                  support_type=geo.PINNED, loads=loads, combinations=combos)
    m.K_effective = VAL_K_EFF
    m.phi_c, m.phi_t = VAL_PHI_C, VAL_PHI_T
    return m


# =============================================================================
# factories
# =============================================================================
def validation_model_a():
    """Validation Model A: 1V pentagonal cap, cases A1 and A2.

    6 nodes, 10 struts, 5 pinned base nodes, 3 free DOF (the crown).
    A1 = explicit crown load (20, 0, -60) kN in case L, combination 1.0L
    (unfactored, flagged serviceability); A2 = 1.2D + 1.6S.
    """
    loads = _validation_lc2_loads()
    m = _validation_model(FRACTION_A, loads, [
        LoadCombination(COMBO_A1, L=1.0, serviceability=True),
        LoadCombination(COMBO_A2, D=1.2, S=1.6)])
    px, py, pz = VAL_A1_LOAD
    m.loads.nodal_loads = [{"node": int(m.geom.crown_node), "Fx": px,
                            "Fy": py, "Fz": pz, "case": "L"}]
    return m


def validation_model_b():
    """Validation Model B: 1V dome on its lower ring, cases B1 and B2.

    11 nodes, 25 struts, 5 pinned base nodes, 18 free DOF.  B1 = 50 kN
    crown point load (case L), combination 1.0L (unfactored, flagged
    serviceability); B2 = 1.2D + 1.6S.
    """
    loads = _validation_lc2_loads()
    loads.crown_point_load = VAL_B1_LOAD
    return _validation_model(FRACTION_B, loads, [
        LoadCombination(COMBO_B1, L=1.0, serviceability=True),
        LoadCombination(COMBO_B2, D=1.2, S=1.6)])


def von_mises_truss(half_span=VM_HALF_SPAN, rise=VM_RISE, P_ref=VM_P_REF,
                    section="CHS 88.9 x 4.0", material=VAL_MATERIAL):
    """von Mises two-bar truss in the x-z plane (validation case F6).

    Nodes 0 = (-a, 0, 0) and 1 = (+a, 0, 0) pinned, node 2 = apex (0, 0, h)
    with only uy restrained; P_ref downward at the apex in case L,
    combination ``P  1.0L``.
    """
    a, h = float(half_span), float(rise)
    if not (a > 0.0 and h > 0.0):
        raise ValueError("half_span and rise must be positive.")
    nodes = [[-a, 0.0, 0.0], [a, 0.0, 0.0], [0.0, 0.0, h]]
    members = [[0, 2], [1, 2]]
    g = geo.custom_geometry(nodes, members, base_nodes=[0, 1], groups=[0, 0],
                            title="von Mises two-bar truss")
    loads = _no_area_loads(nodal_loads=[
        {"node": 2, "Fx": 0.0, "Fy": 0.0, "Fz": -float(P_ref), "case": "L"}])
    pin = [True, True, True]
    return DomeModel(
        geom=g, material_name=material, group_sections=[section],
        support_type=geo.CUSTOM, loads=loads,
        combinations=[LoadCombination("P  1.0L  apex reference load", L=1.0)],
        custom_supports={0: list(pin), 1: list(pin),
                         2: [False, True, False]})


def star_dome_geometry(exact_hexagon=True):
    """Nodes [cm], members, faces and supports of the 24-member star dome.

    Returns ``(nodes_cm (13, 3), members (24, 2), faces (12, 3), supports)``
    with 0-based indices: 0 = crown, 1-6 = inner ring (node 1 on +x), 7-12 =
    supports.  Member order = research notes Sec. 5c (spokes 1-6, ring
    7-12, diagonals 13-24).  ``exact_hexagon=False`` uses the printed 43.3.
    """
    xo = 25.0 * SQ3 if exact_hexagon else 43.3          # 50 cos 30 deg
    yr = 12.5 * SQ3                                     # 25 sin 60 deg
    nodes = np.array([
        [0.0, 0.0, 8.216],                                       # 1 crown
        [25.0, 0.0, 6.216], [12.5, yr, 6.216], [-12.5, yr, 6.216],
        [-25.0, 0.0, 6.216], [-12.5, -yr, 6.216], [12.5, -yr, 6.216],
        [xo, 25.0, 0.0], [0.0, 50.0, 0.0], [-xo, 25.0, 0.0],
        [-xo, -25.0, 0.0], [0.0, -50.0, 0.0], [xo, -25.0, 0.0],
    ])
    one_based = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
                 (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 2),
                 (2, 8), (3, 8), (3, 9), (4, 9), (4, 10), (5, 10),
                 (5, 11), (6, 11), (6, 12), (7, 12), (7, 13), (2, 13)]
    members = np.array(one_based, dtype=int) - 1
    faces_1 = [(1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 6), (1, 6, 7),
               (1, 7, 2), (2, 3, 8), (3, 4, 9), (4, 5, 10), (5, 6, 11),
               (6, 7, 12), (7, 2, 13)]
    faces = np.array(faces_1, dtype=int) - 1
    supports = np.arange(7, 13)
    return nodes, members, faces, supports


def star_dome(section="CHS 26.9 x 2.3", material=VAL_MATERIAL,
              load_ratio=STAR_LOAD_RATIO, exact_hexagon=True):
    """24-member shallow star dome, published geometry converted cm -> m.

    Six pinned supports, one downward crown load P_ref = ``load_ratio`` x EA
    in case L (combination ``P  1.0L``), so the load factor equals
    P/EA / load_ratio (= 10^4 P/EA by default).  Strut families: A = the
    six crown spokes, B = the six inner-ring bars, C = the twelve diagonals.
    """
    nodes_cm, members, faces, supports = star_dome_geometry(exact_hexagon)
    groups = [0] * 6 + [1] * 6 + [2] * 12
    g = geo.custom_geometry(nodes_cm * 0.01, members, faces=faces,
                            base_nodes=supports, groups=groups,
                            title="24-member shallow star dome")
    m = DomeModel(geom=g, material_name=material,
                  group_sections=[section] * 3, support_type=geo.PINNED,
                  loads=_no_area_loads(),
                  combinations=[LoadCombination(
                      "P  1.0L  crown reference load", L=1.0)])
    EA = float(m.EA()[0])
    m.loads.nodal_loads = [{"node": 0, "Fx": 0.0, "Fy": 0.0,
                            "Fz": -float(load_ratio) * EA, "case": "L"}]
    return m


def three_bar_truss(temperature_change=F4_DT, settlement_uz=F4_SETTLEMENT,
                    load=F4_LOAD, material=VAL_MATERIAL,
                    sections=F4_SECTIONS):
    """Planar three-bar indeterminate truss (validation case F4 geometry).

    Nodes 0 = A (-4, 0, 0), 1 = B (0, 0, 0), 2 = C (3, 0, 0) pinned,
    3 = D (0, 0, 3) with uy restrained.  Bars A-D, B-D, C-D (families
    A, B, C in that order).  Load ``load`` at D (case L), self-straining T =
    uniform ``temperature_change`` on every bar + vertical settlement
    ``settlement_uz`` of support C.  Combinations ``P 1.0L``, ``T 1.0T``,
    ``P + T 1.0L + 1.0T``.
    """
    nodes = [[-4.0, 0.0, 0.0], [0.0, 0.0, 0.0], [3.0, 0.0, 0.0],
             [0.0, 0.0, 3.0]]
    members = [[0, 3], [1, 3], [2, 3]]
    g = geo.custom_geometry(nodes, members, base_nodes=[0, 1, 2],
                            groups=[0, 1, 2],
                            title="Planar three-bar indeterminate truss")
    fx, fy, fz = (float(v) for v in load)
    settlements = []
    if float(settlement_uz) != 0.0:
        settlements = [{"node": 2, "ux": 0.0, "uy": 0.0,
                        "uz": float(settlement_uz)}]
    loads = _no_area_loads(
        nodal_loads=[{"node": 3, "Fx": fx, "Fy": fy, "Fz": fz, "case": "L"}],
        temperature_change=float(temperature_change),
        settlements=settlements)
    pin = [True, True, True]
    return DomeModel(
        geom=g, material_name=material, group_sections=list(sections),
        support_type=geo.CUSTOM, loads=loads,
        combinations=[
            LoadCombination("P  1.0L", L=1.0),
            LoadCombination("T  1.0T  temperature + settlement", T=1.0),
            LoadCombination("P + T  1.0L + 1.0T", L=1.0, T=1.0)],
        custom_supports={0: list(pin), 1: list(pin), 2: list(pin),
                         3: [False, True, False]})


def case_study_6v():
    """The 6V, 15 m radius case study (= examples/worked_example.build_model)."""
    geometry = geo.generate_dome(frequency=6, radius=15.0, dome_fraction=0.50)
    model = DomeModel(geom=geometry, material_name="ASTM A500 Gr. C (round)",
                      support_type=geo.PINNED)
    model.group_sections = ["CHS 114.3 x 5.0"] * geometry.n_strut_types()
    model.loads = LoadSpec(
        include_self_weight=True, dead_pressure=350.0, snow_pressure=600.0,
        snow_unbalanced=False, live_pressure=0.0, wind_enabled=True,
        wind_speed=58.0, wind_Kz=0.95, wind_Kzt=1.0, wind_Kd=0.95,
        wind_G=0.85, crown_point_load=0.0)
    model.combinations = [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("LC3  1.2D + 1.0W + 0.5S", D=1.2, W=1.0, S=0.5),
        LoadCombination("LC4  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("SLS  1.0D + 1.0S", D=1.0, S=1.0, serviceability=True),
    ]
    model.K_effective = 1.0
    model.deflection_limit = 240.0
    return model


# =============================================================================
# registry
# =============================================================================
PRESETS = {
    "validation_a": (
        "Validation Model A - 1V pentagonal cap",
        "Hand-calculated validation model (validation/hand_calc.py): 1V cap "
        "of a 5 m sphere, 6 nodes, 10 CHS 88.9 x 4.0 struts, 5 pinned base "
        "nodes, 3 free DOF. A1 = unfactored crown load (+20 kN x, -60 kN z): "
        "crown ux = +0.2724 mm, uz = -1.0697 mm, spoke N1 = -32.23 kN. "
        "A2 = 1.2D + 1.6S: crown uz = -0.6937 mm, spokes -14.80 kN, "
        "utilisation 0.2726, lambda_cr = 5506.1. Every value matches the "
        "hand solution to round-off (validation/output/comparison.csv).",
        validation_model_a),
    "validation_b": (
        "Validation Model B - 1V dome on its lower ring",
        "Hand-calculated validation model solved by virtual work with "
        "five-fold symmetry: 11 nodes, 25 CHS 88.9 x 4.0 struts, 5 pinned "
        "nodes, 18 free DOF. B1 = 50 kN crown load: crown uz = -1.6451 mm, "
        "spokes -19.02 kN, upper-ring hoop +15.39 kN, zigzag -5.878 kN. "
        "B2 = 1.2D + 1.6S: crown uz = -1.8589 mm, total load 211.99 kN.",
        validation_model_b),
    "von_mises": (
        "von Mises two-bar truss (snap-through)",
        "Shallow two-bar truss, half span 2.5 m, rise 0.125 m, CHS 88.9 x 4.0, "
        "P_ref = 10 kN at the apex (validation case F6). Closed-form linear "
        "answer: apex w = P L^3/(2 EA h^2), N = -P L/(2h); linearised "
        "lambda_cr = 2 EA h^3/(L a^2 P_ref) = 5.328. The exact engineering-"
        "strain path P(w) = 2 EA y (1/l - 1/L) has its limit point at "
        "P* = 10.24 kN (w* = 52.9 mm), so lambda_cr overestimates the true "
        "snap-through load 5.2 times: run the nonlinear analysis to see it.",
        von_mises_truss),
    "star_dome": (
        "24-member shallow star dome (snap-through)",
        "The classic Hangai-Kawamata star dome benchmark, geometry of Greco & "
        "Venturini (2006) and Iu & Bradford (2012) converted from cm to m: "
        "crown 8.216 cm high, inner ring r = 25 cm at 6.216 cm, six pinned "
        "supports at r = 50 cm. Crown load P_ref = 1e-4 EA, so the load factor "
        "reads 10^4 P/EA. Reference: first limit point 10^4 P/EA = 3.1566 at a "
        "crown deflection of 0.768 cm (published: near 0.75 cm), zero-load "
        "state at 1.88 cm; linearised lambda_cr = 15.21 is 4.8 times too high.",
        star_dome),
    "three_bar": (
        "Planar three-bar truss - temperature and settlement",
        "Statically indeterminate planar truss (validation case F4 geometry): "
        "bars A-D, B-D, C-D meeting at D (0, 0, 3) m, out-of-plane restrained. "
        "Load (15, 0, -60) kN at D; self-straining load T = +30 degC on every "
        "bar plus a 5 mm settlement of support C. The reference is the force "
        "(flexibility) method with N2 as the redundant: self-straining forces "
        "are of the same order as the load effects in a stiff truss.",
        three_bar_truss),
    "case_study_6v": (
        "Case study - 6V dome, 30 m span",
        "The report's worked example: 6V Class I hemisphere, R = 15 m, ASTM "
        "A500 Gr. C, trial section CHS 114.3 x 5.0, cladding 0.35 kPa, snow "
        "0.60 kPa, wind 58 m/s, ASCE 7 combinations LC1-LC4 + SLS. Run "
        "Auto-size to reproduce examples/worked_example.py.",
        case_study_6v),
}


def preset_keys():
    """Preset keys in display order."""
    return list(PRESETS)


def preset_info(key):
    """``(title, description)`` of a preset."""
    _check_key(key)
    title, desc, _ = PRESETS[key]
    return title, desc


def make_preset(key):
    """A fresh :class:`DomeModel` for ``key`` (every call builds a new one)."""
    _check_key(key)
    return PRESETS[key][2]()


def _check_key(key):
    if key not in PRESETS:
        raise KeyError(f"Unknown preset {key!r}; choose one of "
                       f"{', '.join(PRESETS)}.")


def preset_path(key, directory=None):
    return os.path.join(directory or MODELS_DIR, f"{key}.json")


def write_preset_files(directory=None, keys=None):
    """Save every preset (or ``keys``) as ``<directory>/<key>.json``."""
    directory = directory or MODELS_DIR
    os.makedirs(directory, exist_ok=True)
    paths = []
    for key in (keys or preset_keys()):
        path = preset_path(key, directory)
        make_preset(key).save(path)
        paths.append(path)
    return paths


def load_preset_file(key, directory=None):
    """Load ``<directory>/<key>.json`` (written by :func:`write_preset_files`)."""
    return DomeModel.load(preset_path(key, directory))


# =============================================================================
# independent reference solutions
# =============================================================================
def von_mises_parameters(model):
    """``(a, h, EA, P_ref)`` of a von Mises preset (from its coordinates)."""
    X = model.geom.nodes
    a = 0.5 * float(np.linalg.norm(X[1] - X[0]))
    h = float(X[2, 2] - 0.5 * (X[0, 2] + X[1, 2]))
    EA = float(model.EA()[0])
    P_ref = -sum(float(p.get("Fz", 0.0) or 0.0)
                 for p in model.loads.nodal_loads if int(p["node"]) == 2)
    return a, h, EA, P_ref


def von_mises_load(w, a, h, EA):
    """Exact (engineering-strain) downward apex load at apex deflection ``w``.

    P = 2 EA (L - l)/L * y/l with y = h - w, l = sqrt(a^2 + y^2); L - l is
    evaluated as w (2h - w)/(L + l) so that no digits cancel.  ``w`` may be
    an array.
    """
    w = np.asarray(w, dtype=float)
    L = math.hypot(a, h)
    y = h - w
    l = np.hypot(a, y)
    return 2.0 * EA * (w * (2.0 * h - w) / (L + l)) / L * y / l


def von_mises_limit_point(a, h, EA):
    """First limit point of the exact path: ``{w, P, l, N, y}``.

    dP/dy = 2 EA (a^2/l^3 - 1/L) = 0  ->  l* = (a^2 L)^(1/3).
    """
    L = math.hypot(a, h)
    ls = (a * a * L) ** (1.0 / 3.0)
    ys = math.sqrt(ls * ls - a * a)
    P = 2.0 * EA * ys * (1.0 / ls - 1.0 / L)
    return {"w": h - ys, "P": P, "l": ls, "y": ys, "N": EA * (ls - L) / L}


def von_mises_linear(a, h, EA, P):
    """Linear answers: apex deflection (down), bar force, thrust, lambda_cr."""
    L = math.hypot(a, h)
    return {"w": P * L ** 3 / (2.0 * EA * h * h),
            "N": -P * L / (2.0 * h),
            "H": P * a / (2.0 * h),                  # horizontal support thrust
            "K0": 2.0 * EA * h * h / L ** 3,
            "lambda_cr": 2.0 * EA * h ** 3 / (L * a * a * P)}


def star_dome_symmetric_linear(model, P=None):
    """Linear crown-load response of the star dome by symmetric virtual work.

    Returns ``{"w0", "ur", "w1", "N_spoke", "N_ring", "N_diag", "K_crown"}``
    (SI; w0 and w1 up-positive, ur outward) for a downward crown load ``P``
    (default: the preset's own crown load).  See the module docstring.
    """
    X = model.geom.nodes
    EA = model.EA()
    ref_members = star_dome_geometry()[1]
    if (model.geom.n_nodes != 13
            or not np.array_equal(np.sort(model.geom.members, axis=1),
                                  np.sort(ref_members, axis=1))):
        raise ValueError("Not the 24-member star dome topology.")
    if P is None:
        P = -sum(float(p.get("Fz", 0.0) or 0.0)
                 for p in model.loads.nodal_loads if int(p["node"]) == 0)
    c = X[0]
    ring = X[1:7] - c
    sup = X[7:13] - c
    a = float(np.mean(np.hypot(ring[:, 0], ring[:, 1])))
    b = float(np.mean(np.hypot(sup[:, 0], sup[:, 1])))
    zc, zr, zs = 0.0, float(np.mean(ring[:, 2])), float(np.mean(sup[:, 2]))
    d1 = zc - zr                               # crown above the ring
    hr = zr - zs                               # ring above the supports
    half = math.pi / 6.0                       # 30 deg
    L1 = math.hypot(a, d1)
    L2 = 2.0 * a * math.sin(half)
    L3 = math.sqrt(a * a + b * b - 2.0 * a * b * math.cos(half) + hr * hr)
    g1 = np.array([d1, a, -d1]) / L1
    g2 = np.array([0.0, 2.0 * math.sin(half), 0.0])
    g3 = np.array([0.0, a - b * math.cos(half), hr]) / L3
    EA1, EA2, EA3 = (float(np.mean(EA[s])) for s in
                     (slice(0, 6), slice(6, 12), slice(12, 24)))
    K = (6.0 * EA1 / L1 * np.outer(g1, g1) + 6.0 * EA2 / L2 * np.outer(g2, g2)
         + 12.0 * EA3 / L3 * np.outer(g3, g3))
    q = np.linalg.solve(K, np.array([-float(P), 0.0, 0.0]))
    return {"w0": float(q[0]), "ur": float(q[1]), "w1": float(q[2]),
            "N_spoke": EA1 / L1 * float(g1 @ q),
            "N_ring": EA2 / L2 * float(g2 @ q),
            "N_diag": EA3 / L3 * float(g3 @ q),
            "K_crown": float(P) / abs(float(q[0])) if q[0] != 0.0 else math.inf,
            "L": (L1, L2, L3)}


def _cramer2(a11, a12, a21, a22, b1, b2):
    det = a11 * a22 - a12 * a21
    if abs(det) < 1e-14:
        raise ValueError("Singular 2 x 2 system in the force method.")
    return (b1 * a22 - a12 * b2) / det, (a11 * b2 - b1 * a21) / det


def three_bar_force_method(model, load=(0.0, 0.0), dT=(0.0, 0.0, 0.0),
                           settlement=None, redundant=1):
    """Force (flexibility) method for three bars meeting at one joint.

    ``load`` = (Px, Pz) at the joint [N]; ``dT`` = temperature change of
    each bar [degC]; ``settlement`` = {support node: (dx, dz)} [m].  Works
    in the x-z plane of the preset.  Returns ``{"X", "f11", "N" (3,),
    "u" (ux, uz), "R" {node: (Rx, Rz)}}``.  Equations: module docstring.
    """
    g = model.geom
    Mb = np.asarray(g.members, dtype=int)
    if len(Mb) != 3:
        raise ValueError("The force method here needs exactly three bars.")
    joint = set(Mb[0]) & set(Mb[1]) & set(Mb[2])
    if len(joint) != 1:
        raise ValueError("The three bars must meet at one joint.")
    D = joint.pop()
    sup = [int(i if j == D else j) for i, j in Mb]
    X = g.nodes
    EA = model.EA()
    alpha = model.material.alpha
    L = np.zeros(3)
    c = np.zeros((3, 2))
    for k in range(3):
        v = X[D] - X[sup[k]]
        L[k] = float(np.linalg.norm(v))
        c[k] = np.array([v[0], v[2]]) / L[k]
    fl = L / EA
    dT = np.asarray(dT, dtype=float)
    P = np.asarray(load, dtype=float)
    dl = {s: np.zeros(2) for s in sup}
    for node, (dx, dz) in (settlement or {}).items():
        dl[int(node)] = np.array([float(dx), float(dz)])

    r = int(redundant)
    p, q = [k for k in range(3) if k != r]

    def primary(rhs):
        n_p, n_q = _cramer2(c[p][0], c[q][0], c[p][1], c[q][1],
                            rhs[0], rhs[1])
        n = np.zeros(3)
        n[p], n[q] = n_p, n_q
        return n

    n0 = primary(P)
    n1 = primary(-c[r])
    n1[r] = 1.0
    f11 = float(np.sum(n1 ** 2 * fl))
    D1P = float(np.sum(n1 * n0 * fl))
    D1T = float(np.sum(n1 * alpha * dT * L))
    D1S = -sum(float((-n1[k] * c[k]) @ dl[sup[k]]) for k in range(3))
    Xr = -(D1P + D1T + D1S) / f11
    N = n0 + Xr * n1
    e = N * fl + alpha * dT * L
    u = np.zeros(2)
    for a in range(2):
        unit = np.zeros(2)
        unit[a] = 1.0
        ne = primary(unit)
        u[a] = (sum(ne[k] * e[k] for k in (p, q))
                - sum(float((-ne[k] * c[k]) @ dl[sup[k]]) for k in (p, q)))
    R = {sup[k]: tuple(-N[k] * c[k]) for k in range(3)}
    return {"X": Xr, "f11": f11, "N": N, "u": u, "R": R, "joint": D,
            "supports": sup}


# =============================================================================
# verification
# =============================================================================
def _row(rows, key, case, quantity, reference, program, unit, scale=1.0,
         tol=1e-9, source=""):
    r = float(reference) * scale
    p = float(program) * scale
    if abs(r) > 1e-12:
        diff, kind = abs(p - r) / abs(r), "rel"
    else:
        diff, kind = abs(p - r), "abs"
    rows.append({"preset": key, "case": case, "quantity": quantity,
                 "reference": r, "program": p, "unit": unit,
                 "difference": diff, "diff_kind": kind, "tolerance": tol,
                 "pass": bool(diff <= tol), "source": source})


def _identity_rows(rows, key, res):
    """Equilibrium, reaction balance and energy identity of every case."""
    for name in res.order:
        r = res.cases[name]
        scale = max(float(np.abs(r.F).sum()), float(np.abs(r.reactions).sum()),
                    1.0)
        _row(rows, key, name, "equilibrium residual (relative)", 0.0,
             r.equilibrium_residual, "-", tol=1e-9, source="identity")
        _row(rows, key, name, "|sum F + sum R| / scale", 0.0,
             float(np.linalg.norm(r.reaction_balance)) / scale, "-",
             tol=1e-9, source="identity")
        _row(rows, key, name, "energy identity error", 0.0, r.energy_error,
             "-", tol=1e-9, source="identity")


def _check_validation_a(model, res, rows):
    key = "validation_a"
    src = "validation/hand_calc.py (direct stiffness by hand)"
    g = model.geom
    c = g.crown_node
    hn2 = int(g.base_nodes[np.argmax(g.nodes[g.base_nodes, 0])])   # hand node 2
    spoke1 = next(k for k, (i, j) in enumerate(g.members)
                  if {int(i), int(j)} == {c, hn2})
    spokes = [k for k, (i, j) in enumerate(g.members) if c in (i, j)]
    r1 = res.cases[COMBO_A1]
    r2 = res.cases[COMBO_A2]
    _row(rows, key, COMBO_A1, "crown ux", HAND_A["A1_ux"], r1.u[3 * c],
         "mm", 1e3, source=src)
    _row(rows, key, COMBO_A1, "crown uz", HAND_A["A1_uz"], r1.u[3 * c + 2],
         "mm", 1e3, source=src)
    _row(rows, key, COMBO_A1, "N member 1 (spoke on +x)", HAND_A["A1_N1"],
         r1.N[spoke1], "kN", 1e-3, source=src)
    _row(rows, key, COMBO_A1, "sum Rz", HAND_A["A1_sumRz"],
         r1.reaction_vectors().sum(axis=0)[2], "kN", 1e-3, source=src)
    _row(rows, key, COMBO_A2, "crown uz", HAND_A["A2_uz"], r2.u[3 * c + 2],
         "mm", 1e3, source=src)
    _row(rows, key, COMBO_A2, "total vertical load", HAND_A["A2_total_Pz"],
         r2.F[2::3].sum(), "kN", 1e-3, source=src)
    _row(rows, key, COMBO_A2, "N spoke (mean of 5)", HAND_A["A2_N_spoke"],
         np.mean(r2.N[spokes]), "kN", 1e-3, source=src)
    _row(rows, key, COMBO_A2, "Rz per support node", HAND_A["A2_Rz_node"],
         r2.reaction_vectors()[hn2, 2], "kN", 1e-3, source=src)
    _row(rows, key, "envelope", "AISC utilisation of a spoke",
         HAND_A["A2_util"], res.util_env[spoke1], "-", source=src)
    lam = res.buckling.load_factors[0] if res.buckling is not None else np.nan
    _row(rows, key, res.buckling.reference_case if res.buckling else "-",
         "lambda_cr", HAND_A["A2_lambda_cr"], lam, "-", source=src)


def _check_validation_b(model, res, rows):
    key = "validation_b"
    src = "validation/hand_calc.py (virtual work, five-fold symmetry)"
    g = model.geom
    X = g.nodes
    c = g.crown_node
    zmax = X[:, 2].max()
    tol_z = 1e-9 * VAL_RADIUS
    upper = [n for n in range(g.n_nodes)
             if n != c and X[n, 2] > tol_z and abs(X[n, 2] - zmax) > tol_z]
    lower = [n for n in range(g.n_nodes) if abs(X[n, 2]) <= tol_z]
    us, ls = set(upper), set(lower)
    spokes = [k for k, (i, j) in enumerate(g.members) if c in (i, j)]
    hoop = [k for k, (i, j) in enumerate(g.members) if {i, j} <= us]
    zig = [k for k, (i, j) in enumerate(g.members)
           if (i in us and j in ls) or (j in us and i in ls)]
    for combo, tag in ((COMBO_B1, "B1"), (COMBO_B2, "B2")):
        hb = HAND_B[tag]
        r = res.cases[combo]
        U = r.u.reshape(-1, 3)
        ur = [U[n, 0] * X[n, 0] / math.hypot(X[n, 0], X[n, 1])
              + U[n, 1] * X[n, 1] / math.hypot(X[n, 0], X[n, 1])
              for n in upper]
        Rv = r.reaction_vectors()
        rr = [Rv[n, 0] * X[n, 0] / math.hypot(X[n, 0], X[n, 1])
              + Rv[n, 1] * X[n, 1] / math.hypot(X[n, 0], X[n, 1])
              for n in lower]
        _row(rows, key, combo, "total vertical load", hb["tot_P"],
             r.F[2::3].sum(), "kN", 1e-3, source=src)
        _row(rows, key, combo, "crown uz (w0)", hb["w0"], U[c, 2], "mm", 1e3,
             source=src)
        _row(rows, key, combo, "upper ring radial (u_r)", hb["ur"],
             np.mean(ur), "mm", 1e3, source=src)
        _row(rows, key, combo, "upper ring vertical (w1)", hb["w1"],
             np.mean(U[upper, 2]), "mm", 1e3, source=src)
        _row(rows, key, combo, "N spokes", hb["N_spoke"], np.mean(r.N[spokes]),
             "kN", 1e-3, source=src)
        _row(rows, key, combo, "N upper ring (hoop)", hb["N_ring"],
             np.mean(r.N[hoop]), "kN", 1e-3, source=src)
        _row(rows, key, combo, "N zigzag", hb["N_zig"], np.mean(r.N[zig]),
             "kN", 1e-3, source=src)
        _row(rows, key, combo, "Rz per support node", hb["Rz_node"],
             np.mean(Rv[lower, 2]), "kN", 1e-3, source=src)
        _row(rows, key, combo, "R radial per support node", hb["Rr_node"],
             np.mean(rr), "kN", 1e-3, source=src)


def _check_von_mises(model, res, rows):
    key = "von_mises"
    src = "closed form (linear statics of the two-bar truss)"
    a, h, EA, P = von_mises_parameters(model)
    ref = von_mises_linear(a, h, EA, P)
    name = res.order[0]
    r = res.cases[name]
    R = r.reaction_vectors()
    _row(rows, key, name, "apex uz = -P L^3/(2 EA h^2)", -ref["w"],
         r.u[3 * 2 + 2], "mm", 1e3, source=src)
    _row(rows, key, name, "apex ux (symmetry)", 0.0, r.u[3 * 2], "mm", 1e3,
         source=src)
    for k in range(2):
        _row(rows, key, name, f"N bar {k + 1} = -P L/(2h)", ref["N"], r.N[k],
             "kN", 1e-3, source=src)
    _row(rows, key, name, "Rz at each support = P/2", 0.5 * P,
         0.5 * (R[0, 2] + R[1, 2]), "kN", 1e-3, source=src)
    _row(rows, key, name, "thrust at left support = +P a/(2h)", ref["H"],
         R[0, 0], "kN", 1e-3, source=src)
    _row(rows, key, name, "thrust at right support = -P a/(2h)", -ref["H"],
         R[1, 0], "kN", 1e-3, source=src)
    lam = res.buckling.load_factors[0] if res.buckling is not None else np.nan
    _row(rows, key, name, "lambda_cr = 2 EA h^3/(L a^2 P)", ref["lambda_cr"],
         lam, "-", source="closed form (2-DOF eigenproblem)")


def _check_star_dome(model, res, rows):
    key = "star_dome"
    src = "symmetric virtual work, 3 generalised DOF (independent)"
    sym = star_dome_symmetric_linear(model)
    name = res.order[0]
    r = res.cases[name]
    U = r.u.reshape(-1, 3)
    EA = float(model.EA()[0])
    P = -float(r.F[2])
    _row(rows, key, name, "crown uz (w0)", sym["w0"], U[0, 2], "mm", 1e3,
         source=src)
    _row(rows, key, name, "ring vertical (w1)", sym["w1"],
         np.mean(U[1:7, 2]), "mm", 1e3, source=src)
    rad = [float(U[n, :2] @ model.geom.nodes[n, :2])
           / math.hypot(*model.geom.nodes[n, :2]) for n in range(1, 7)]
    _row(rows, key, name, "ring radial (u_r)", sym["ur"], np.mean(rad), "mm",
         1e3, source=src)
    _row(rows, key, name, "N spokes", sym["N_spoke"], np.mean(r.N[0:6]), "kN",
         1e-3, source=src)
    _row(rows, key, name, "N inner ring", sym["N_ring"], np.mean(r.N[6:12]),
         "kN", 1e-3, source=src)
    _row(rows, key, name, "N diagonals", sym["N_diag"], np.mean(r.N[12:24]),
         "kN", 1e-3, source=src)
    _row(rows, key, name, "N spread within the 3 families", 0.0,
         max(np.ptp(r.N[0:6]), np.ptp(r.N[6:12]), np.ptp(r.N[12:24])),
         "kN", 1e-3, source="six-fold symmetry")
    _row(rows, key, name, "initial crown stiffness P/(w EA) [1/m]",
         STAR_K0_PER_EA, P / abs(U[0, 2]) / EA, "-", tol=1e-4,
         source="research_notes.md Sec. 5d (0.090430, 5 digits)")
    lam = res.buckling.load_factors[0] if res.buckling is not None else np.nan
    _row(rows, key, name, "lambda_cr P/EA", STAR_LAMBDA_CR_P_OVER_EA,
         lam * P / EA, "-", tol=1e-5,
         source="research_notes.md Sec. 5d (1.52097e-3, 6 digits)")


def _three_bar_case_inputs(model, combo):
    """(Px, Pz), per-bar dT and settlements of one combination (from data)."""
    f = combo.factors()
    P = np.zeros(2)
    for p in model.loads.nodal_loads:
        fac = f.get(str(p.get("case", "L")).upper(), 0.0)
        P += fac * np.array([float(p.get("Fx", 0.0) or 0.0),
                             float(p.get("Fz", 0.0) or 0.0)])
    T = f.get("T", 0.0)
    dT = np.full(3, T * float(model.loads.temperature_change))
    sett = {}
    for s in model.loads.settlements:
        sett[int(s["node"])] = (T * float(s.get("ux", 0.0) or 0.0),
                                T * float(s.get("uz", 0.0) or 0.0))
    return P, dT, sett


def _check_three_bar(model, res, rows):
    key = "three_bar"
    src = "force (flexibility) method, redundant N2"
    for combo in model.combinations:
        r = res.cases[combo.name]
        P, dT, sett = _three_bar_case_inputs(model, combo)
        fm = three_bar_force_method(model, P, dT, sett)
        D = fm["joint"]
        for k in range(3):
            _row(rows, key, combo.name, f"N bar {k + 1}", fm["N"][k], r.N[k],
                 "kN", 1e-3, source=src)
        _row(rows, key, combo.name, "ux of D", fm["u"][0], r.u[3 * D], "mm",
             1e3, source=src + " + unit-load method")
        _row(rows, key, combo.name, "uz of D", fm["u"][1], r.u[3 * D + 2],
             "mm", 1e3, source=src + " + unit-load method")
        R = r.reaction_vectors()
        for s, (rx, rz) in fm["R"].items():
            _row(rows, key, combo.name, f"Rx at node {s}", rx, R[s, 0], "kN",
                 1e-3, source=src)
            _row(rows, key, combo.name, f"Rz at node {s}", rz, R[s, 2], "kN",
                 1e-3, source=src)
        _row(rows, key, combo.name, "out-of-plane reaction at D", 0.0,
             R[D, 1], "kN", 1e-3, source="planar structure")


def _check_case_study(model, res, rows):
    key = "case_study_6v"
    src = "statics: bar lengths x A x rho g + true / projected panel areas"
    g = model.geom
    X = g.nodes
    Lm = np.linalg.norm(X[g.members[:, 1]] - X[g.members[:, 0]], axis=1)
    W_self = float(np.sum(model.member_areas() * Lm)) * model.material.rho \
        * G_ACCEL
    cr = np.cross(X[g.faces[:, 1]] - X[g.faces[:, 0]],
                  X[g.faces[:, 2]] - X[g.faces[:, 0]])
    area_true = 0.5 * float(np.linalg.norm(cr, axis=1).sum())
    area_plan = 0.5 * float(np.abs(cr[:, 2]).sum())
    ls = model.loads
    D_tot = (W_self if ls.include_self_weight else 0.0) \
        + ls.dead_pressure * area_true
    S_tot = 0.0 if ls.snow_unbalanced else ls.snow_pressure * area_plan
    for combo in model.combinations:
        if abs(combo.W) > 0.0 or abs(combo.L) > 0.0 or ls.snow_unbalanced:
            continue                      # wind / live: no independent total
        r = res.cases[combo.name]
        expect = combo.D * D_tot + combo.S * S_tot
        _row(rows, key, combo.name, "sum Rz = factored gravity load", expect,
             r.reaction_vectors().sum(axis=0)[2], "kN", 1e-3, source=src)


_CHECKS = {
    "validation_a": _check_validation_a,
    "validation_b": _check_validation_b,
    "von_mises": _check_von_mises,
    "star_dome": _check_star_dome,
    "three_bar": _check_three_bar,
    "case_study_6v": _check_case_study,
}


def verify_preset(key, model=None, progress=None):
    """Run a preset and compare the program with its independent reference.

    Returns a list of row dicts (preset, case, quantity, reference, program,
    unit, difference, diff_kind "rel" | "abs", tolerance, pass, source).
    Values are in the display unit given in ``unit``.  Every row passes when
    the program is right.
    """
    from .postprocess import run_analysis           # avoids an import cycle
    _check_key(key)
    model = model if model is not None else make_preset(key)
    res = run_analysis(model, n_buckling_modes=6, do_buckling=True,
                       progress=progress)
    rows = []
    _CHECKS[key](model, res, rows)
    _identity_rows(rows, key, res)
    return rows


def verify_all(keys=None):
    """``{key: rows}`` for every preset (see :func:`verify_preset`)."""
    return {k: verify_preset(k) for k in (keys or preset_keys())}


def format_rows(rows):
    """Plain-text table of :func:`verify_preset` rows."""
    out = [f"{'case':<28}{'quantity':<42}{'reference':>14}{'program':>14}"
           f"{'unit':>6}{'diff':>11}  ok"]
    for r in rows:
        d = (f"{r['difference']:.1e}" if r["diff_kind"] == "rel"
             else f"{r['difference']:.1e}a")
        out.append(f"{r['case'][:27]:<28}{r['quantity'][:41]:<42}"
                   f"{r['reference']:>14.6g}{r['program']:>14.6g}"
                   f"{r['unit']:>6}{d:>11}  {'yes' if r['pass'] else 'NO'}")
    return "\n".join(out)


# =============================================================================
# command line
# =============================================================================
def main(argv=None):
    """``python -m geodome.benchmarks [--write [DIR]] [keys ...]``."""
    ap = argparse.ArgumentParser(
        prog="python -m geodome.benchmarks",
        description="List the GeoDomeLab benchmark presets and check each "
                    "one against its independent reference solution.")
    ap.add_argument("keys", nargs="*", help="presets to check (default all)")
    ap.add_argument("--write", nargs="?", const=MODELS_DIR, default=None,
                    metavar="DIR", help="also save the presets as JSON "
                    f"(default folder {MODELS_DIR})")
    ap.add_argument("--quiet", action="store_true",
                    help="print only the summary line per preset")
    args = ap.parse_args(argv)
    keys = args.keys or preset_keys()
    for k in keys:
        _check_key(k)
    if args.write is not None:
        for p in write_preset_files(args.write, keys):
            print(f"wrote {p}")
    total = failed = 0
    for k in keys:
        title, _ = preset_info(k)
        rows = verify_preset(k)
        nbad = sum(not r["pass"] for r in rows)
        total += len(rows)
        failed += nbad
        print("=" * 118)
        print(f"{k}: {title}  --  {len(rows) - nbad} of {len(rows)} checks pass")
        if not args.quiet:
            print("-" * 118)
            print(format_rows(rows))
    print("=" * 118)
    print(f"{total - failed} of {total} benchmark checks pass "
          f"({len(keys)} presets).")
    return 0 if failed == 0 else 1


if __name__ == "__main__":                            # pragma: no cover
    sys.exit(main())
