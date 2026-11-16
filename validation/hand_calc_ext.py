"""
Extended hand calculations for GeoDomeLab 2.0 (validation cases F1 to F8).

This file is the "by hand" side of the extended validation.  It does NOT
import the geodome package, so it cannot inherit a mistake from the program.
It imports validation/hand_calc.py only to reuse the exact definition of
validation Model A (1V pentagonal cap, R = 5 m, CHS 88.9 x 4.0, ASTM A500
Gr. C), and hand_calc.py itself never imports geodome either.

Every case is worked the way a person does it on paper:

    structure and data -> hand method (equilibrium, compatibility, force or
    displacement method) -> every intermediate number -> result ->
    an independent check (joint or free-body equilibrium, energy, symmetry,
    a second method, or a limiting case)

numpy is used as a calculator (dot and cross products, small sums).  No
general matrix solver is used for any structure: 2 x 2 systems are solved by
Cramer's rule and the 3 x 3 Model A systems by the written-out Gaussian
elimination of hand_calc.gauss_3x3.

CASES
  F1  heated fixed-fixed bar (inclined 3D bar)  +  two-segment stepped bar
  F2  Model A under a uniform temperature change  +  one sun-side spoke
  F3  Model A with one support settling (vertical)  +  one spreading (radial)
  F4  planar three-bar indeterminate truss by the FORCE (flexibility) METHOD:
      a point load, a temperature change in one bar, a support settlement
  F5  bar + spring: in series, in parallel, and a heated bar on a spring
  F6  von Mises two-bar truss: closed-form engineering-strain path, limit
      point, snap-through target, linearised buckling factor (2 DOF)
  F7  Model A crown natural frequencies with lumped masses
  F8  2V Class I geodesic sphere: chord factors, axial angles, face angles,
      dihedral angles, counts, hub and panel types

CONVENTIONS (same as hand_calc.py and GeoDomeLab)
  z vertical (up), tension positive, SI units inside, kN / mm / MPa printed.
  Member i -> j has unit vector c = (x_j - x_i)/L.
  Elongation e = c . (u_j - u_i).  Free (initial) strain eps0 = alpha*dT.
  Member force N = EA (e/L - eps0).
  Force of the member on its i end node = +N c, on its j end node = -N c.
  Reaction R at a support = minus (sum of member forces on that node + load).
  Thermal equivalent nodal load F0 = EA*eps0*[-c ; +c].

OUTPUTS
  printed step-by-step derivation (python validation/hand_calc_ext.py)
  results()      -> dict {key: {"value", "unit", "formula", "description"}}
                    keyed like "F4.N2_temperature", values in SI units
  definitions()  -> the structures of every case as plain data (nodes,
                    members, sections, supports, loads), 0-based indices
  validation/output/hand_calc_ext.md  readable derivation for the report

Run:  python validation/hand_calc_ext.py
"""

from __future__ import annotations

import math
import os
import sys
import textwrap

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import hand_calc as hc                       # noqa: E402  (no geodome inside)

OUT_DIR = os.path.join(HERE, "output")
MD_PATH = os.path.join(OUT_DIR, "hand_calc_ext.md")

# =============================================================================
# PARAMETERS  (edit here only)
# =============================================================================
E = hc.E                     # 200 GPa, ASTM A500 Gr. C                  [Pa]
FY = hc.FY                   # 345 MPa                                   [Pa]
RHO = hc.RHO                 # 7850 kg/m3
G = hc.G                     # 9.80665 m/s2
ALPHA = 11.7e-6              # thermal expansion of carbon steel  [1/degC]
#                              (AISC 6.5e-6 /degF = 11.7e-6 /degC)
DT = 30.0                    # uniform heating above erection temp.    [degC]

SEC_MAIN = ("CHS 88.9 x 4.0", hc.D_OUT, hc.T_WALL)     # validation section
SEC_SMALL = ("CHS 76.1 x 4.0", 0.0761, 0.0040)         # second section

# F1  heated bars
F1_X1 = (0.0, 0.0, 0.0)      # inclined bar, length exactly 3 m       [m]
F1_X2 = (2.0, 1.0, 2.0)
F1B_L1 = 2.0                 # stepped bar: CHS 88.9 x 4.0 segment     [m]
F1B_L2 = 1.5                 #              CHS 76.1 x 4.0 segment     [m]

# F3  Model A support movements
F3_SETTLE = 0.005            # vertical settlement of hand node 2      [m]
F3_SPREAD = 0.005            # outward radial spread of hand node 3    [m]

# F4  planar three-bar truss in the x-z plane (out-of-plane y restrained)
F4_NODES = {"A": (-4.0, 0.0, 0.0), "B": (0.0, 0.0, 0.0),
            "C": (3.0, 0.0, 0.0), "D": (0.0, 0.0, 3.0)}
F4_BARS = {1: ("A", "D", SEC_MAIN), 2: ("B", "D", SEC_SMALL),
           3: ("C", "D", SEC_MAIN)}
F4_P = (15e3, 0.0, -60e3)    # load at D (Fx, Fy, Fz)                   [N]
F4_HEATED_BAR = 2            # bar heated by DT
F4_SETTLE_NODE = "C"
F4_SETTLE = (0.0, 0.0, -0.005)   # settlement of support C              [m]

# F5  bar + spring
F5_L = 3.0                   # bar length along x                       [m]
F5_K = 5.0e7                 # spring stiffness (GeoDomeLab default)  [N/m]
F5_P = 50e3                  # axial load                              [N]

# F6  von Mises two-bar truss
VM_A = 2.5                   # half span                                [m]
VM_H = 0.125                 # rise  (h/a = 0.05, shallow)              [m]
VM_PREF = 10e3               # reference load for the eigenproblem      [N]
VM_SAMPLES = (0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5)   # w/h for the P(w) table

# F8  2V geodesic sphere
F8_R = 5.0                   # radius for the metric lengths            [m]
KENNER_2V = {"A": 0.54653, "B": 0.61803}   # Kenner (1976) chord factors

DEG = math.pi / 180.0
SQ5 = math.sqrt(5.0)


# =============================================================================
# small helpers
# =============================================================================
def chs(D, t):
    """Circular hollow section properties (same formulas as hand_calc)."""
    d = D - 2.0 * t
    A = math.pi / 4.0 * (D * D - d * d)
    I = math.pi / 64.0 * (D ** 4 - d ** 4)
    return {"D": D, "t": t, "A": A, "I": I, "r": math.sqrt(I / A)}


def g6(x, nd=6):
    """Compact number format, exact zero printed as 0."""
    x = float(x)
    if x == 0.0:
        return "0"
    return f"{x:.{nd}g}"


def vec(v, scale=1.0, nd=6):
    return "(" + ", ".join(g6(float(x) * scale, nd) for x in v) + ")"


def kN(x, nd=6):
    return g6(float(x) * 1e-3, nd)


def mm(x, nd=6):
    return g6(float(x) * 1e3, nd)


def cramer2(a11, a12, a21, a22, b1, b2):
    """Solve [a11 a12; a21 a22] x = b by Cramer's rule (hand method)."""
    det = a11 * a22 - a12 * a21
    x1 = (b1 * a22 - a12 * b2) / det
    x2 = (a11 * b2 - b1 * a21) / det
    return x1, x2, det


class Doc:
    """Collects the derivation once and renders it as text or markdown."""

    def __init__(self):
        self.blocks = []

    def h1(self, t):
        self.blocks.append(("h1", t))

    def h2(self, t):
        self.blocks.append(("h2", t))

    def h3(self, t):
        self.blocks.append(("h3", t))

    def p(self, t):
        self.blocks.append(("p", " ".join(t.split())))

    def eq(self, *lines):
        self.blocks.append(("eq", list(lines)))

    def table(self, header, rows):
        self.blocks.append(("t", list(header),
                            [[str(x) for x in r] for r in rows]))

    def items(self, lst):
        self.blocks.append(("l", list(lst)))

    # -------------------------------------------------------------- render --
    def markdown(self):
        out = []
        for b in self.blocks:
            kind = b[0]
            if kind == "h1":
                out += ["# " + b[1], ""]
            elif kind == "h2":
                out += ["## " + b[1], ""]
            elif kind == "h3":
                out += ["### " + b[1], ""]
            elif kind == "p":
                out += [b[1], ""]
            elif kind == "eq":
                out += ["```text"] + b[1] + ["```", ""]
            elif kind == "t":
                hdr, rows = b[1], b[2]

                def cell(x):
                    return x.replace("|", "\\|")
                out.append("| " + " | ".join(cell(x) for x in hdr) + " |")
                out.append("|" + "|".join("---" for _ in hdr) + "|")
                for r in rows:
                    out.append("| " + " | ".join(cell(x) for x in r) + " |")
                out.append("")
            elif kind == "l":
                out += ["- " + x for x in b[1]] + [""]
        return "\n".join(out)

    def text(self):
        out = []
        for b in self.blocks:
            kind = b[0]
            if kind == "h1":
                out += ["=" * 78, b[1], "=" * 78]
            elif kind == "h2":
                out += ["", "=" * 78, b[1], "=" * 78]
            elif kind == "h3":
                out += ["", "--- " + b[1] + " " + "-" * max(0, 72 - len(b[1]))]
            elif kind == "p":
                out += textwrap.wrap(b[1], 78)
            elif kind == "eq":
                out += ["    " + x for x in b[1]]
            elif kind == "t":
                hdr, rows = b[1], b[2]
                w = [max(len(hdr[i]), *(len(r[i]) for r in rows))
                     for i in range(len(hdr))]
                out.append("  " + "  ".join(h.rjust(w[i])
                                            for i, h in enumerate(hdr)))
                out.append("  " + "  ".join("-" * x for x in w))
                for r in rows:
                    out.append("  " + "  ".join(c.rjust(w[i])
                                                for i, c in enumerate(r)))
            elif kind == "l":
                for x in b[1]:
                    out += textwrap.wrap(x, 76, initial_indent="  - ",
                                         subsequent_indent="    ")
        return "\n".join(out)


class Results(dict):
    """Reference values: key -> value (SI), unit, formula, description."""

    def add(self, key, value, unit, formula, description):
        self[key] = {"value": float(value), "unit": unit,
                     "formula": formula, "description": description}


CHECKS = []          # (case, description, residual, tolerance)


def check(case, what, residual, tol):
    CHECKS.append((case, what, float(abs(residual)), float(tol)))
    return abs(residual) <= tol


# =============================================================================
# MODEL A helper (members, reactions, energy for any crown/support movement)
# =============================================================================
def model_a_state(MA, u_crown, us=None, eps0=None, P_crown=None):
    """Member forces and support reactions of Model A by joint equilibrium.

    u_crown : crown displacement (3,)
    us      : {hand node: prescribed displacement} of supported nodes
    eps0    : {member: free strain}
    P_crown : load applied at the crown (3,)
    """
    nodes, members, geo = MA["nodes"], MA["members"], MA["geo"]
    EA, k = MA["section"]["EA"], MA["k"]
    d = {n: np.zeros(3) for n in nodes}
    d[1] = np.asarray(u_crown, float)
    for n, v in (us or {}).items():
        d[n] = np.asarray(v, float)
    eps0 = eps0 or {}
    Pc = np.zeros(3) if P_crown is None else np.asarray(P_crown, float)
    N, elong = {}, {}
    for m, (i, j) in members.items():
        c = geo[m][1]
        elong[m] = float(c @ (d[j] - d[i]))
        N[m] = k * elong[m] - EA * eps0.get(m, 0.0)
    fn = {n: np.zeros(3) for n in nodes}          # member forces on nodes
    for m, (i, j) in members.items():
        c = geo[m][1]
        fn[i] += N[m] * c
        fn[j] -= N[m] * c
    crown_res = fn[1] + Pc                        # crown equilibrium
    R = {n: -fn[n] for n in nodes if n != 1}      # no loads on supports here
    sumR = sum(R.values()) + Pc
    mom = sum(np.cross(nodes[n], R[n]) for n in R) + np.cross(nodes[1], Pc)
    U = sum(N[m] ** 2 / (2.0 * k) for m in members)       # N^2 L / (2 EA)
    W_sup = 0.5 * sum(float(R[n] @ d[n]) for n in R)
    W_load = 0.5 * float(Pc @ d[1])
    return {"d": d, "N": N, "elong": elong, "R": R, "crown_res": crown_res,
            "sumR": sumR, "mom": mom, "U": U, "W_sup": W_sup,
            "W_load": W_load}


def _model_a_basics():
    MA = hc.model_a()
    return MA


# =============================================================================
# F1  heated fixed-fixed bar
# =============================================================================
def case_f1(doc, res):
    sec = chs(SEC_MAIN[1], SEC_MAIN[2])
    secb = chs(SEC_SMALL[1], SEC_SMALL[2])
    A, I = sec["A"], sec["I"]
    EA = E * A

    doc.h2("F1  Heated fixed-fixed bar")
    doc.p("Purpose. The most basic self-straining case. A bar that is held "
          "at both ends cannot expand, so heating it produces a compressive "
          "force that does not depend on its length. This checks the sign "
          "and size of the thermal strain eps0 = alpha dT and the thermal "
          "reactions in GeoDomeLab.")
    doc.h3("F1a  Inclined bar pinned at both ends")
    x1 = np.array(F1_X1, float)
    x2 = np.array(F1_X2, float)
    v = x2 - x1
    L = math.sqrt(float(v @ v))
    c = v / L
    eT = ALPHA * DT * L
    N = -EA * ALPHA * DT
    sigma = N / A
    eN = N * L / EA
    R1 = -N * c
    R2 = N * c
    F0_1 = -EA * ALPHA * DT * c
    F0_2 = EA * ALPHA * DT * c
    Pe = math.pi ** 2 * E * I / L ** 2
    dT_cr = Pe / (EA * ALPHA)

    doc.p(f"Data. {SEC_MAIN[0]}: A = {g6(A * 1e6)} mm2, "
          f"I = {g6(I * 1e12)} mm4. E = {g6(E / 1e9)} GPa, "
          f"alpha = {g6(ALPHA)} 1/degC, dT = +{g6(DT)} degC. "
          f"Node 1 = {vec(x1)} m, node 2 = {vec(x2)} m, both pinned "
          f"(all three translations fixed).")
    doc.eq(f"L = |x2 - x1| = sqrt(2^2 + 1^2 + 2^2) = {g6(L)} m",
           f"c = (x2 - x1)/L = {vec(c)}",
           f"EA = {kN(EA, 8)} kN",
           "",
           "Compatibility: the ends cannot move, so the total elongation is 0",
           "  e = N L/(EA) + alpha dT L = 0",
           "  N = -EA alpha dT",
           f"  free thermal elongation alpha dT L = {mm(eT)} mm",
           f"  N = -{kN(EA, 8)} x {g6(ALPHA)} x {g6(DT)} = {kN(N)} kN"
           "  (compression)",
           f"  sigma = N/A = -E alpha dT = {g6(sigma / 1e6)} MPa"
           f"  ({g6(abs(sigma) / FY, 4)} Fy)",
           "",
           "Reactions (force of the bar on node 1 is +N c, on node 2 is -N c):",
           f"  R1 = -N c = {vec(R1, 1e-3)} kN",
           f"  R2 = +N c = {vec(R2, 1e-3)} kN")
    doc.p("Independent check by the stiffness method. With every degree of "
          "freedom fixed, u = 0 and the reactions equal minus the thermal "
          "equivalent nodal loads F0 = EA alpha dT [-c ; +c].")
    doc.eq(f"F0 at node 1 = -EA alpha dT c = {vec(F0_1, 1e-3)} kN"
           f"  ->  R1 = -F0 = {vec(-F0_1, 1e-3)} kN",
           f"F0 at node 2 = +EA alpha dT c = {vec(F0_2, 1e-3)} kN"
           f"  ->  R2 = -F0 = {vec(-F0_2, 1e-3)} kN")
    compat = eN + eT
    sumR = R1 + R2
    mom = np.cross(x1, R1) + np.cross(x2, R2)
    check("F1a", "total elongation N L/EA + alpha dT L [m]", compat, 1e-15)
    check("F1a", "R1 + R2 [N]", np.linalg.norm(sumR), 1e-9)
    check("F1a", "moment of reactions about the origin [N m]",
          np.linalg.norm(mom), 1e-9)
    check("F1a", "R1 (equilibrium) - R1 (stiffness method) [N]",
          np.linalg.norm(R1 + F0_1), 1e-9)
    doc.items([f"Total elongation e = {g6(compat, 3)} m (must be 0).",
               f"Sum of reactions = {vec(sumR, 1, 3)} N, moment about the "
               f"origin = {vec(mom, 1, 3)} N m.",
               f"Reactions by equilibrium and by -F0 differ by "
               f"{g6(np.linalg.norm(R1 + F0_1), 3)} N.",
               f"Engineering note: the Euler load of this 3 m pin-ended bar "
               f"is Pe = pi^2 EI/L^2 = {kN(Pe)} kN, so the bar would buckle "
               f"at dT = Pe/(EA alpha) = {g6(dT_cr, 4)} degC. The compressive "
               f"thermal force does not depend on L, but the buckling "
               f"temperature falls with L^2."])

    res.add("F1.L", L, "m", "L = |x2 - x1|", "length of the inclined bar")
    res.add("F1.N", N, "N", "N = -E*A*alpha*dT",
            "axial force in the heated fixed-fixed bar")
    res.add("F1.sigma", sigma, "Pa", "sigma = -E*alpha*dT",
            "axial stress in the heated fixed-fixed bar")
    for ax, lab in enumerate("xyz"):
        res.add(f"F1.R1_{lab}", R1[ax], "N", "R1 = -N*c",
                f"reaction {lab} at node 1 (0,0,0)")
        res.add(f"F1.R2_{lab}", R2[ax], "N", "R2 = +N*c",
                f"reaction {lab} at node 2 (2,1,2)")
    res.add("F1.free_elongation", eT, "m", "e_T = alpha*dT*L",
            "elongation if the bar were free to expand")
    res.add("F1.dT_euler", dT_cr, "degC", "dT_cr = pi^2*I/(alpha*A*L^2)",
            "temperature rise that buckles the 3 m bar (information)")

    # ---------------------------------------------------------- F1b -------
    doc.h3("F1b  Stepped bar between two fixed ends (one redundant)")
    L1, L2 = F1B_L1, F1B_L2
    EA1, EA2 = E * A, E * secb["A"]
    doc.p(f"Data. Collinear bar along x: node 1 (0,0,0) fixed, node 2 "
          f"({g6(L1)},0,0) free in x (y, z restrained), node 3 "
          f"({g6(L1 + L2)},0,0) fixed. Segment 1-2 is {SEC_MAIN[0]} "
          f"(EA1 = {kN(EA1)} kN), segment 2-3 is {SEC_SMALL[0]} "
          f"(A = {g6(secb['A'] * 1e6)} mm2, EA2 = {kN(EA2)} kN). "
          f"Both segments are heated by dT = +{g6(DT)} degC. Force method "
          f"with the redundant X = axial force N (equal in both segments by "
          f"equilibrium of node 2): release node 3 in x, let the bar expand "
          f"freely, then push it back.")
    dfree = ALPHA * DT * (L1 + L2)
    f = L1 / EA1 + L2 / EA2
    Nb = -dfree / f
    u2 = (Nb / EA1 + ALPHA * DT) * L1
    e2 = (Nb / EA2 + ALPHA * DT) * L2
    k_sum = EA1 / L1 + EA2 / L2
    F0_2 = ALPHA * DT * (EA1 - EA2)
    u2K = F0_2 / k_sum
    R1x = -Nb
    R3x = Nb
    doc.eq(f"free expansion at the release  d_T = alpha dT (L1 + L2) = "
           f"{mm(dfree)} mm",
           f"flexibility  f = L1/EA1 + L2/EA2 = {g6(f)} m/N",
           "compatibility  d_T + f N = 0",
           f"  N = -d_T / f = {kN(Nb)} kN",
           "node 2 moves by the elongation of segment 1",
           f"  u2 = (N/EA1 + alpha dT) L1 = {mm(u2)} mm",
           f"reactions  R1x = -N = {kN(R1x)} kN,  R3x = +N = {kN(R3x)} kN")
    doc.p("Independent check by the stiffness method (one free DOF): "
          "(EA1/L1 + EA2/L2) u2 = F0 with F0 = alpha dT (EA1 - EA2), the "
          "j end load of segment 1 plus the i end load of segment 2.")
    doc.eq(f"u2 = alpha dT (EA1 - EA2)/(EA1/L1 + EA2/L2) = {mm(u2K)} mm",
           f"segment 2 elongation (N/EA2 + alpha dT) L2 = {mm(e2)} mm"
           " = -u2 (node 3 fixed)")
    check("F1b", "u2 force method - u2 stiffness method [m]", u2 - u2K, 1e-15)
    check("F1b", "segment 2 elongation + u2 [m]", e2 + u2, 1e-15)
    doc.items([f"u2 by the two methods differs by {g6(u2 - u2K, 3)} m.",
               f"Compatibility of segment 2: e2 + u2 = {g6(e2 + u2, 3)} m."])
    res.add("F1.stepped_N", Nb, "N",
            "N = -alpha*dT*(L1+L2)/(L1/EA1 + L2/EA2)",
            "axial force in both segments of the heated stepped bar")
    res.add("F1.stepped_u2", u2, "m", "u2 = (N/EA1 + alpha*dT)*L1",
            "x displacement of the middle node of the stepped bar")
    res.add("F1.stepped_R1x", R1x, "N", "R1x = -N",
            "x reaction at node 1 of the stepped bar")


# =============================================================================
# F2  Model A under a temperature change
# =============================================================================
def case_f2(doc, res, MA):
    EA, k, s, t, L = MA["section"]["EA"], MA["k"], MA["s"], MA["t"], MA["L"]
    Kff = MA["Kff"]
    geo = MA["geo"]
    H = hc.R - hc.H_RING
    phiPn = MA["aisc"]["phiPn_c"]

    doc.h2("F2  Model A (1V cap) under a uniform temperature change")
    doc.p(f"Structure. Validation Model A of hand_calc.py: the 1V pentagonal "
          f"cap of a sphere of radius R = {g6(hc.R)} m. Node 1 is the crown "
          f"(free), nodes 2 to 6 are the pinned base ring at 0, 72, 144, 216 "
          f"and 288 deg. Members 1 to 5 are the spokes (ring node to crown), "
          f"members 6 to 10 the ring (2-3, 3-4, 4-5, 5-6, 6-2). All ten "
          f"struts have the 1V length L = {g6(L)} m and the section "
          f"{SEC_MAIN[0]}, EA = {kN(EA)} kN, k = EA/L = {kN(k)} kN/m. "
          f"Crown height H = R - R/sqrt5 = {g6(H)} m, ring plan radius "
          f"a = 2R/sqrt5 = {g6(hc.A_RING)} m.")
    doc.eq("spoke direction (ring node -> crown)  c_i = (-s cos th_i, "
           "-s sin th_i, t)",
           f"  s = a/L = {g6(s)},  t = H/L = {g6(t)}",
           f"free block  Kff = k diag(5/2 s^2, 5/2 s^2, 5 t^2) = "
           f"diag({kN(Kff[0, 0])}, {kN(Kff[1, 1])}, {kN(Kff[2, 2])}) kN/m")
    doc.h3(f"F2a  Every strut heated by dT = +{g6(DT)} degC")
    doc.p("Hand method (compatibility first). The five ring struts join two "
          "pinned nodes, so each one is the fixed-fixed bar of F1 and carries "
          "N = -EA alpha dT. The crown is free. If it rises by w, every spoke "
          "lengthens by c_i . (0,0,w) = t w, the same for all five. Choosing "
          "t w = alpha dT L makes each spoke exactly as long as it wants to "
          "be, so the spokes carry no force. Then the crown is in equilibrium "
          "with zero spoke forces, every member is compatible and every "
          "force obeys the material law. A linear elastic solution that "
          "satisfies all three is the unique solution.")
    eT = ALPHA * DT * L
    w = eT / t
    N_ring = -EA * ALPHA * DT
    eps0 = {m: ALPHA * DT for m in range(1, 11)}
    st = model_a_state(MA, np.array([0.0, 0.0, w]), eps0=eps0)
    Rr = {n: float(st["R"][n] @ np.array([math.cos(72.0 * (n - 2) * DEG),
                                          math.sin(72.0 * (n - 2) * DEG),
                                          0.0])) for n in range(2, 7)}
    Rr_formula = -2.0 * math.sin(36.0 * DEG) * EA * ALPHA * DT
    doc.eq(f"alpha dT L = {mm(eT)} mm",
           "crown rise  w = alpha dT L / t = alpha dT L^2 / H",
           f"  w = {mm(eT)} / {g6(t)} = {mm(w)} mm",
           "spokes  N = k (t w) - EA alpha dT = 0",
           f"ring    N = -EA alpha dT = {kN(N_ring)} kN",
           "",
           "Reaction at a ring node (two ring bars in compression push the",
           "node outwards, the unit vectors to the two neighbours add to",
           "-2 sin36 e_r):",
           "  R = -2 sin36 EA alpha dT e_r   (inward)",
           f"  R_radial = {kN(Rr_formula)} kN,  R_z = 0")
    doc.table(["member", "type", "elongation [mm]", "N [kN]"],
              [[m, "spoke" if m <= 5 else "ring", mm(st["elong"][m]),
                kN(st["N"][m])] for m in range(1, 11)])
    doc.p("Check by the stiffness method. The thermal equivalent load at the "
          "crown (the j end of every spoke) is F0 = sum EA alpha dT c_i = "
          "EA alpha dT (0, 0, 5t). The ring loads act on supported DOFs only.")
    F0 = sum(EA * ALPHA * DT * geo[m][1] for m in range(1, 6))
    uK = F0 / np.diag(Kff)
    doc.eq(f"F0 = {vec(F0, 1e-3)} kN",
           f"u = Kff^-1 F0 = {vec(uK, 1e3)} mm")
    check("F2a", "crown uz compatibility - stiffness [m]", uK[2] - w, 1e-15)
    check("F2a", "crown ux, uy by stiffness [m]", np.linalg.norm(uK[:2]),
          1e-18)
    check("F2a", "crown equilibrium [N]", np.linalg.norm(st["crown_res"]),
          1e-8)
    check("F2a", "sum of reactions [N]", np.linalg.norm(st["sumR"]), 1e-8)
    check("F2a", "moment of reactions [N m]", np.linalg.norm(st["mom"]), 1e-7)
    check("F2a", "radial reaction joint sum - formula [N]",
          max(abs(Rr[n] - Rr_formula) for n in Rr), 1e-8)
    doc.items([f"Crown rise by compatibility and by Kff^-1 F0 differ by "
               f"{g6(uK[2] - w, 3)} m.",
               f"Crown equilibrium residual {g6(np.linalg.norm(st['crown_res']), 3)} N, "
               f"sum of reactions {g6(np.linalg.norm(st['sumR']), 3)} N, "
               f"moment {g6(np.linalg.norm(st['mom']), 3)} N m "
               f"(a self-straining load has self-equilibrated reactions).",
               f"Radial reaction from the joint sums equals the formula to "
               f"{g6(max(abs(Rr[n] - Rr_formula) for n in Rr), 3)} N at all "
               f"five supports (five-fold symmetry).",
               f"Engineering note: the ring force {kN(N_ring, 4)} kN does not "
               f"depend on the strut length. The AISC E3 design strength of "
               f"this {g6(L, 4)} m strut is phi Pn = {kN(phiPn, 4)} kN "
               f"(hand_calc.py), so a {g6(DT)} degC rise alone gives a "
               f"utilisation of {g6(abs(N_ring) / phiPn, 3)}. A ring between "
               f"pinned "
               f"supports needs either a temperature check or a sliding "
               f"bearing. The spokes are stress free because the crown can "
               f"rise {mm(w, 4)} mm, which is L^2/H^2 = {g6(1 / t ** 2, 4)} "
               f"times the free rise alpha dT H of an unrestrained cap."])
    res.add("F2.crown_uz", w, "m", "w = alpha*dT*L^2/H",
            "crown rise of Model A, all struts heated")
    res.add("F2.crown_ux", 0.0, "m", "symmetry", "crown x movement (zero)")
    res.add("F2.N_spoke", 0.0, "N", "N = k*t*w - EA*alpha*dT = 0",
            "force in each spoke (stress free)")
    res.add("F2.N_ring", N_ring, "N", "N = -E*A*alpha*dT",
            "force in each ring strut (fixed-fixed)")
    res.add("F2.R_radial", Rr_formula, "N", "R_r = -2*sin(36deg)*E*A*alpha*dT",
            "radial reaction at every support node (+ outward)")
    res.add("F2.R_z", 0.0, "N", "spokes carry no force",
            "vertical reaction at every support node")
    noise = 1e-12 * EA * ALPHA * DT          # round-off on exact zeros
    for n in range(2, 7):
        for ax, lab in enumerate("xyz"):
            val = st["R"][n][ax]
            res.add(f"F2.R{n}_{lab}", 0.0 if abs(val) < noise else val, "N",
                    "R = -2*sin36*EA*alpha*dT*e_r",
                    f"reaction {lab} at hand node {n}, all struts heated")

    # ------------------------------------------------------------- F2b -----
    doc.h3(f"F2b  Only spoke 1 heated by dT = +{g6(DT)} degC (sun on one side)")
    doc.p("A non-uniform temperature field. Only member 1 (node 2 to the "
          "crown) gets eps0 = alpha dT. The crown now moves sideways and the "
          "five spokes share a self-equilibrated set of forces. Stiffness "
          "method with the thermal equivalent load of member 1 at the crown.")
    c1 = geo[1][1]
    F0b = EA * ALPHA * DT * c1
    ub, steps = hc.gauss_3x3(Kff, F0b)
    stb = model_a_state(MA, ub, eps0={1: ALPHA * DT})
    doc.eq(f"F0 = EA alpha dT c_1 = {kN(EA * ALPHA * DT)} x {vec(c1)}"
           f" = {vec(F0b, 1e-3)} kN",
           "Kff u = F0   (Gaussian elimination, multipliers "
           + ", ".join(g6(x[2], 3) for x in steps["multipliers"]) + ")",
           f"u = {vec(ub, 1e3)} mm",
           "closed form  u = alpha dT L ( -1/(2.5 s), 0, 1/(5 t) )",
           f"            = {vec(ALPHA * DT * L * np.array([-1 / (2.5 * s), 0.0, 1 / (5 * t)]), 1e3)} mm",
           "N_1 = EA alpha dT (0.4 + 0.2 - 1) = -0.4 EA alpha dT",
           "N_i = EA alpha dT (0.4 cos th_i + 0.2)   for i = 2..5")
    closed = {1: -0.4}
    for i in range(2, 6):
        closed[i] = 0.4 * math.cos(72.0 * (i - 1) * DEG) + 0.2
    doc.table(["member", "N [kN]", "N / (EA alpha dT)", "closed form"],
              [[m, kN(stb["N"][m]), g6(stb["N"][m] / (EA * ALPHA * DT)),
                g6(closed[m]) if m in closed else "0"]
               for m in range(1, 11)])
    err = max(abs(stb["N"][m] / (EA * ALPHA * DT) - closed[m])
              for m in range(1, 6))
    check("F2b", "spoke forces vs closed form [-]", err, 1e-13)
    check("F2b", "crown equilibrium [N]", np.linalg.norm(stb["crown_res"]),
          1e-8)
    check("F2b", "sum of reactions [N]", np.linalg.norm(stb["sumR"]), 1e-8)
    check("F2b", "sum of spoke forces (z equilibrium) [N]",
          sum(stb["N"][m] for m in range(1, 6)), 1e-8)
    doc.items([f"Spoke forces match the closed form to {g6(err, 3)} "
               f"(relative to EA alpha dT).",
               f"Crown equilibrium residual "
               f"{g6(np.linalg.norm(stb['crown_res']), 3)} N, sum of "
               f"reactions {g6(np.linalg.norm(stb['sumR']), 3)} N.",
               "The five spoke forces add to zero, as vertical equilibrium "
               "of the unloaded crown requires (all spokes have the same "
               "rise t)."])
    for ax, lab in enumerate("xyz"):
        res.add(f"F2.oneSpoke_crown_u{lab}", ub[ax], "m",
                "u = Kff^-1 * EA*alpha*dT*c_1",
                f"crown {lab} displacement, only spoke 1 heated")
    for m in range(1, 6):
        res.add(f"F2.oneSpoke_N{m}", stb["N"][m], "N",
                "N1 = -0.4*EA*alpha*dT, Ni = EA*alpha*dT*(0.4cos(th_i)+0.2)",
                f"force in spoke {m}, only spoke 1 heated")


# =============================================================================
# F3  Model A with a support movement
# =============================================================================
def case_f3(doc, res, MA):
    k, s, t = MA["k"], MA["s"], MA["t"]
    Kff, geo = MA["Kff"], MA["geo"]

    doc.h2("F3  Model A with one support movement")
    doc.p("Prescribed support displacements are the second self-straining "
          "load. The stiffness method with prescribed values us on the "
          "restrained DOFs solves Kff uf = Ff - Kfs us. The crown is linked "
          "to a ring node only through that node's spoke, so the coupling "
          "block is Kfs = -k c c^T for that spoke and the right-hand side is "
          "k c (c . us). Ring struts do not touch the crown.")

    def solve(label, node, us_vec, member):
        c = geo[member][1]
        rhs = k * c * float(c @ us_vec)
        u, steps = hc.gauss_3x3(Kff, rhs)
        st = model_a_state(MA, u, us={node: us_vec})
        doc.eq(f"us (node {node}) = {vec(us_vec, 1e3)} mm,   spoke {member}"
               f" c = {vec(c)}",
               f"Kfs us = -k c (c . us),  c . us = {mm(c @ us_vec)} mm",
               f"right-hand side  -Kfs us = {vec(rhs, 1e-3)} kN",
               "",
               "3 x 3 system (kN/m, mm -> kN):",
               f"  [{kN(Kff[0, 0], 7):>10} {kN(Kff[0, 1], 3):>6} "
               f"{kN(Kff[0, 2], 3):>6} ] {{D1}}   {{{kN(rhs[0], 7):>10}}}",
               f"  [{kN(Kff[1, 0], 3):>10} {kN(Kff[1, 1], 7):>6} "
               f"{kN(Kff[1, 2], 3):>6} ] {{D2}} = {{{kN(rhs[1], 7):>10}}}",
               f"  [{kN(Kff[2, 0], 3):>10} {kN(Kff[2, 1], 3):>6} "
               f"{kN(Kff[2, 2], 7):>6} ] {{D3}}   {{{kN(rhs[2], 7):>10}}}",
               "  elimination multipliers (row, pivot, m): "
               + ", ".join(f"({a},{b},{g6(mlt, 3)})"
                           for a, b, mlt in steps["multipliers"]),
               "  back substitution  D3 = b3/U33, D2 = b2/U22, D1 = b1/U11",
               f"  D = {vec(u, 1e3)} mm")
        doc.table(["member", "type", "elongation [mm]", "N [kN]"],
                  [[m, "spoke" if m <= 5 else "ring", mm(st["elong"][m]),
                    kN(st["N"][m])] for m in range(1, 11)])
        doc.table(["support", "Rx [kN]", "Ry [kN]", "Rz [kN]"],
                  [[n, kN(st["R"][n][0]), kN(st["R"][n][1]),
                    kN(st["R"][n][2])] for n in range(2, 7)])
        return u, st

    # ------------------------------------------------------------- F3a -----
    doc.h3(f"F3a  Hand node 2 settles vertically by {mm(F3_SETTLE)} mm")
    d = F3_SETTLE
    us = np.array([0.0, 0.0, -d])
    ua, sa = solve("vert", 2, us, 1)
    doc.p("Closed form. With c_1 = (-s, 0, t) and us = (0, 0, -delta) the "
          "right-hand side is k t delta (s, 0, -t), so")
    cf_u = np.array([t * d / (2.5 * s), 0.0, -d / 5.0])
    cf_N = {1: 0.4}
    for i in range(2, 6):
        cf_N[i] = -0.4 * math.cos(72.0 * (i - 1) * DEG) - 0.2
    doc.eq("D1 = t delta/(2.5 s),  D2 = 0,  D3 = -delta/5",
           f"   = {vec(cf_u, 1e3)} mm",
           "N_1 = 0.4 k t delta,   N_i = -k t delta (0.4 cos th_i + 0.2)",
           f"k t delta = {kN(k * t * d)} kN")
    err_u = np.linalg.norm(ua - cf_u)
    err_N = max(abs(sa["N"][m] / (k * t * d) - cf_N[m]) for m in range(1, 6))
    check("F3a", "crown displacement vs closed form [m]", err_u, 1e-17)
    check("F3a", "spoke forces vs closed form [-]", err_N, 1e-13)
    check("F3a", "crown equilibrium [N]", np.linalg.norm(sa["crown_res"]),
          1e-8)
    check("F3a", "sum of reactions [N]", np.linalg.norm(sa["sumR"]), 1e-8)
    check("F3a", "moment of reactions [N m]", np.linalg.norm(sa["mom"]), 1e-7)
    check("F3a", "energy: U - 1/2 R.us [J]", sa["U"] - sa["W_sup"], 1e-9)
    doc.p("Independent checks. Joint equilibrium of the crown, global "
          "equilibrium of the self-equilibrated reactions, and Clapeyron's "
          "theorem for a prescribed displacement: the strain energy stored "
          "in the struts equals the work of the support reaction on the "
          "settlement, U = sum N^2 L/(2EA) = 1/2 R2 . us.")
    doc.items([f"Crown displacement and spoke forces equal the closed form "
               f"to {g6(err_u, 3)} m and {g6(err_N, 3)} (relative).",
               f"Crown equilibrium {g6(np.linalg.norm(sa['crown_res']), 3)} N, "
               f"sum of reactions {g6(np.linalg.norm(sa['sumR']), 3)} N, "
               f"moment {g6(np.linalg.norm(sa['mom']), 3)} N m.",
               f"U = {g6(sa['U'])} J, 1/2 R2 . us = {g6(sa['W_sup'])} J, "
               f"difference {g6(sa['U'] - sa['W_sup'], 3)} J.",
               "The ring struts stay stress free because they are horizontal "
               "and a vertical settlement does not change their length "
               "(first-order theory).",
               f"Engineering note: a {mm(d, 3)} mm settlement puts "
               f"{kN(min(sa['N'][m] for m in range(1, 6)), 4)} kN of "
               f"compression into spokes 2 and 5. Settlement forces scale "
               f"with EA, so stiffer struts attract more of them."])
    for ax, lab in enumerate("xyz"):
        res.add(f"F3.vert_crown_u{lab}", ua[ax], "m",
                "Kff*u = k*c1*(c1.us); u = (t*d/(2.5s), 0, -d/5)",
                f"crown {lab} displacement, node 2 settles {mm(d, 3)} mm")
    for m in range(1, 11):
        res.add(f"F3.vert_N{m}", sa["N"][m], "N",
                "N1 = 0.4 k t d, Ni = -k t d (0.4 cos th_i + 0.2), ring 0",
                f"force in member {m}, node 2 settles")
    for n in range(2, 7):
        for ax, lab in enumerate("xyz"):
            res.add(f"F3.vert_R{n}_{lab}", sa["R"][n][ax], "N",
                    "R = -(sum of member forces on the node)",
                    f"reaction {lab} at hand node {n}, node 2 settles")
    res.add("F3.vert_energy", sa["U"], "J", "U = sum N^2 L/(2EA) = 1/2 R.us",
            "strain energy, node 2 settles")

    # ------------------------------------------------------------- F3b -----
    doc.h3(f"F3b  Hand node 3 spreads radially outwards by {mm(F3_SPREAD)} mm")
    D = F3_SPREAD
    th = 72.0 * DEG
    er = np.array([math.cos(th), math.sin(th), 0.0])
    us = D * er
    ub, sb = solve("spread", 3, us, 2)
    cf_ub = 0.4 * D * er + np.array([0.0, 0.0, -s * D / (5.0 * t)])
    cf_Nb = {2: 0.4}
    for i in (1, 3, 4, 5):
        cf_Nb[i] = -0.4 * math.cos(72.0 * (i - 2) * DEG) - 0.2
    N_ring_cf = k * D * math.sin(36.0 * DEG)
    doc.p("Closed form. c_2 . us = -s Delta, so the crown moves by 0.4 Delta "
          "along e_r(72 deg) and drops by s Delta/(5 t). The two ring struts "
          "at node 3 are stretched by Delta sin36.")
    doc.eq("D = 0.4 Delta e_r(72) + (0, 0, -s Delta/(5 t))",
           f"  = {vec(cf_ub, 1e3)} mm",
           "N_2 = 0.4 k s Delta,  N_i = -k s Delta (0.4 cos(th_i - 72) + 0.2)",
           f"ring 6 and 7: N = k Delta sin36 = {kN(N_ring_cf)} kN")
    err_ub = np.linalg.norm(ub - cf_ub)
    err_Nb = max(abs(sb["N"][m] / (k * s * D) - cf_Nb[m]) for m in range(1, 6))
    err_ring = max(abs(sb["N"][6] - N_ring_cf), abs(sb["N"][7] - N_ring_cf),
                   abs(sb["N"][8]), abs(sb["N"][9]), abs(sb["N"][10]))
    check("F3b", "crown displacement vs closed form [m]", err_ub, 1e-17)
    check("F3b", "spoke forces vs closed form [-]", err_Nb, 1e-13)
    check("F3b", "ring forces vs closed form [N]", err_ring, 1e-8)
    check("F3b", "crown equilibrium [N]", np.linalg.norm(sb["crown_res"]),
          1e-8)
    check("F3b", "sum of reactions [N]", np.linalg.norm(sb["sumR"]), 1e-8)
    check("F3b", "moment of reactions [N m]", np.linalg.norm(sb["mom"]), 1e-7)
    check("F3b", "energy: U - 1/2 R.us [J]", sb["U"] - sb["W_sup"], 1e-9)
    doc.items([f"Crown displacement, spoke and ring forces equal the closed "
               f"form to {g6(err_ub, 3)} m, {g6(err_Nb, 3)} (relative) and "
               f"{g6(err_ring, 3)} N.",
               f"Crown equilibrium {g6(np.linalg.norm(sb['crown_res']), 3)} N, "
               f"sum of reactions {g6(np.linalg.norm(sb['sumR']), 3)} N, "
               f"moment {g6(np.linalg.norm(sb['mom']), 3)} N m.",
               f"U = {g6(sb['U'])} J equals 1/2 R3 . us = {g6(sb['W_sup'])} J "
               f"(difference {g6(sb['U'] - sb['W_sup'], 3)} J).",
               "F3a and F2b have the same pattern 0.4 / (0.4 cos + 0.2): "
               "heating one spoke and moving its support are both a change "
               "of that spoke's free length."])
    for ax, lab in enumerate("xyz"):
        res.add(f"F3.spread_crown_u{lab}", ub[ax], "m",
                "u = 0.4*D*e_r(72) + (0,0,-s*D/(5t))",
                f"crown {lab} displacement, node 3 spreads {mm(D, 3)} mm")
    for m in range(1, 11):
        res.add(f"F3.spread_N{m}", sb["N"][m], "N",
                "N2 = 0.4 k s D, ring 6,7 = k D sin36",
                f"force in member {m}, node 3 spreads")
    for n in range(2, 7):
        for ax, lab in enumerate("xyz"):
            res.add(f"F3.spread_R{n}_{lab}", sb["R"][n][ax], "N",
                    "R = -(sum of member forces on the node)",
                    f"reaction {lab} at hand node {n}, node 3 spreads")
    res.add("F3.spread_energy", sb["U"], "J",
            "U = sum N^2 L/(2EA) = 1/2 R.us", "strain energy, node 3 spreads")


# =============================================================================
# F4  planar three-bar truss, force method
# =============================================================================
def case_f4(doc, res):
    nodes = {kk: np.array(v, float) for kk, v in F4_NODES.items()}
    sup = {1: "A", 2: "B", 3: "C"}
    L, c, EA, fl, A, name = {}, {}, {}, {}, {}, {}
    for i, (a, b, secdef) in F4_BARS.items():
        v = nodes[b] - nodes[a]
        L[i] = math.sqrt(float(v @ v))
        c[i] = np.array([v[0], v[2]]) / L[i]          # (x, z) components
        sec = chs(secdef[1], secdef[2])
        A[i] = sec["A"]
        EA[i] = E * sec["A"]
        fl[i] = L[i] / EA[i]
        name[i] = secdef[0]
    P = np.array([F4_P[0], F4_P[2]])
    dT = {i: (DT if i == F4_HEATED_BAR else 0.0) for i in F4_BARS}
    dlt = {kk: np.zeros(2) for kk in ("A", "B", "C")}
    dlt[F4_SETTLE_NODE] = np.array([F4_SETTLE[0], F4_SETTLE[2]])

    doc.h2("F4  Planar three-bar truss by the force (flexibility) method")
    doc.p("Purpose. A statically indeterminate truss solved by the force "
          "method, which shares nothing with the stiffness method used by "
          "GeoDomeLab: the unknowns are forces, compatibility is written "
          "with the principle of virtual forces, and displacements come "
          "last from the unit-load method. Three load types are treated: "
          "a point load, a temperature change in one bar and a support "
          "settlement.")
    doc.p(f"Structure. Bars lie in the x-z plane and meet at the free node "
          f"D = {vec(nodes['D'])} m. Supports (pinned): "
          f"A = {vec(nodes['A'])}, B = {vec(nodes['B'])}, "
          f"C = {vec(nodes['C'])} m. Node D is restrained out of plane "
          f"(uy = 0) so it has two free DOFs (ux, uz). E = {g6(E / 1e9)} GPa, "
          f"alpha = {g6(ALPHA)} 1/degC.")
    doc.table(["bar", "from-to", "section", "A [mm2]", "L [m]", "cx", "cz",
               "L/EA [m/MN]"],
              [[i, f"{sup[i]}-D", name[i], g6(A[i] * 1e6), g6(L[i]),
                g6(c[i][0]), g6(c[i][1]), g6(fl[i] * 1e6)] for i in F4_BARS])
    doc.p(f"Loads. Case P: Px = {kN(P[0])} kN, Pz = {kN(P[1])} kN at D. "
          f"Case T: bar {F4_HEATED_BAR} heated by dT = +{g6(DT)} degC. "
          f"Case S: support {F4_SETTLE_NODE} settles by "
          f"{vec(F4_SETTLE, 1e3)} mm.")

    doc.h3("Step 1  Degree of indeterminacy and the primary structure")
    doc.eq("plane truss: i = m + r - 2j = 3 + 6 - 2*4 = 1",
           "redundant X1 = N2 (force in bar 2), released by cutting bar 2",
           "primary structure = bars 1 and 3 (a determinate two-bar truss)",
           "joint D: forces on D are -N_i c_i, so sum N_i c_i = P")

    def primary(rhs):
        n1, n3, det = cramer2(c[1][0], c[3][0], c[1][1], c[3][1],
                              rhs[0], rhs[1])
        return {1: n1, 2: 0.0, 3: n3}, det

    n0, det = primary(P)
    n1v, _ = primary(-c[2])
    n1v[2] = 1.0
    nex, _ = primary(np.array([1.0, 0.0]))
    nez, _ = primary(np.array([0.0, 1.0]))
    doc.h3("Step 2  Primary structure under the load (n0) and under X1 = 1 (n1)")
    doc.eq(f"[c1 c3] n = rhs,  det = c1x c3z - c3x c1z = {g6(det)}",
           "n_1 = (rhs_x c3z - c3x rhs_z)/det,  n_3 = (c1x rhs_z - rhs_x c1z)/det",
           f"load:      rhs = P   = {vec(P, 1e-3)} kN"
           f"  ->  n0 = ({kN(n0[1])}, 0, {kN(n0[3])}) kN",
           f"X1 = 1:    rhs = -c2 = {vec(-c[2])}"
           f"  ->  n1 = ({g6(n1v[1])}, 1, {g6(n1v[3])})",
           f"unit x:    rhs = (1,0)  ->  n_x = ({g6(nex[1])}, 0, {g6(nex[3])})",
           f"unit z:    rhs = (0,1)  ->  n_z = ({g6(nez[1])}, 0, {g6(nez[3])})")

    doc.h3("Step 3  Flexibility coefficient and compatibility terms")
    f11 = sum(n1v[i] ** 2 * fl[i] for i in F4_BARS)
    D1P = sum(n1v[i] * n0[i] * fl[i] for i in F4_BARS)
    D1T = sum(n1v[i] * ALPHA * dT[i] * L[i] for i in F4_BARS)
    Rbar = {sup[i]: -n1v[i] * c[i] for i in F4_BARS}   # reactions, X1 = 1
    D1S = -sum(float(Rbar[kk] @ dlt[kk]) for kk in Rbar)
    doc.eq("virtual forces (X1 = 1 system) x real deformations:",
           "  sum n1_i e_i - sum Rbar_k . delta_k = 0   (gap at the cut = 0)",
           "  e_i = (n0_i + X1 n1_i) L_i/EA_i + alpha dT_i L_i",
           "  f11 X1 + D1P + D1T + D1S = 0",
           "  f11 = sum n1_i^2 L_i/EA_i",
           "  D1P = sum n1_i n0_i L_i/EA_i",
           "  D1T = sum n1_i alpha dT_i L_i",
           "  D1S = -sum Rbar_k . delta_k,  Rbar_k = -n1_i c_i at the support"
           " of bar i")
    doc.table(["bar", "n0 [kN]", "n1 [-]", "n1^2 L/EA [m/MN]",
               "n1 n0 L/EA [mm]", "n1 alpha dT L [mm]"],
              [[i, kN(n0[i]), g6(n1v[i]), g6(n1v[i] ** 2 * fl[i] * 1e6),
                mm(n1v[i] * n0[i] * fl[i]),
                mm(n1v[i] * ALPHA * dT[i] * L[i])] for i in F4_BARS]
              + [["sum", "", "", g6(f11 * 1e6), mm(D1P), mm(D1T)]])
    doc.eq(f"Rbar_{F4_SETTLE_NODE} = -n1_3 c3 = "
           f"{vec(Rbar[F4_SETTLE_NODE])},  delta_{F4_SETTLE_NODE} = "
           f"{vec(dlt[F4_SETTLE_NODE], 1e3)} mm",
           f"D1S = {mm(D1S)} mm",
           f"f11 = {g6(f11)} m/N")

    cases = {"load": (D1P, P, False, False),
             "temperature": (D1T, np.zeros(2), True, False),
             "settlement": (D1S, np.zeros(2), False, True),
             "combined": (D1P + D1T + D1S, P, True, True)}
    out = {}
    doc.h3("Step 4  Redundant, bar forces, displacements and reactions")
    rows = []
    for cname, (Dl, Pc, hasT, hasS) in cases.items():
        X = -Dl / f11
        N = {i: (n0[i] if Pc.any() else 0.0) + X * n1v[i] for i in F4_BARS}
        dTc = {i: (dT[i] if hasT else 0.0) for i in F4_BARS}
        dl = {kk: (dlt[kk] if hasS else np.zeros(2)) for kk in dlt}
        e = {i: N[i] * fl[i] + ALPHA * dTc[i] * L[i] for i in F4_BARS}
        # unit-load method on the primary structure (bars 1 and 3)
        ux = (nex[1] * e[1] + nex[3] * e[3]
              - (float((-nex[1] * c[1]) @ dl["A"])
                 + float((-nex[3] * c[3]) @ dl["C"])))
        uz = (nez[1] * e[1] + nez[3] * e[3]
              - (float((-nez[1] * c[1]) @ dl["A"])
                 + float((-nez[3] * c[3]) @ dl["C"])))
        u = np.array([ux, uz])
        R = {sup[i]: -N[i] * c[i] for i in F4_BARS}
        # checks
        compat = max(abs(float(c[i] @ (u - dl[sup[i]])) - e[i])
                     for i in F4_BARS)
        eqD = np.linalg.norm(sum(N[i] * c[i] for i in F4_BARS) - Pc)
        sumR = np.linalg.norm(sum(R.values()) + Pc)
        xz = {kk: np.array([nodes[kk][0], nodes[kk][2]]) for kk in nodes}
        mom = (sum(xz[kk][1] * R[kk][0] - xz[kk][0] * R[kk][1] for kk in R)
               + xz["D"][1] * Pc[0] - xz["D"][0] * Pc[1])
        # displacement-method cross-check (2 x 2, Cramer)
        K11 = sum(EA[i] / L[i] * c[i][0] * c[i][0] for i in F4_BARS)
        K12 = sum(EA[i] / L[i] * c[i][0] * c[i][1] for i in F4_BARS)
        K22 = sum(EA[i] / L[i] * c[i][1] * c[i][1] for i in F4_BARS)
        Fe = (Pc + sum(EA[i] * ALPHA * dTc[i] * c[i] for i in F4_BARS)
              + sum(EA[i] / L[i] * c[i] * float(c[i] @ dl[sup[i]])
                    for i in F4_BARS))
        uKx, uKz, _ = cramer2(K11, K12, K12, K22, Fe[0], Fe[1])
        uK = np.array([uKx, uKz])
        NK = {i: EA[i] / L[i] * float(c[i] @ (uK - dl[sup[i]]))
              - EA[i] * ALPHA * dTc[i] for i in F4_BARS}
        dN = max(abs(NK[i] - N[i]) for i in F4_BARS)
        du = np.linalg.norm(uK - u)
        tag = f"F4 {cname}"
        check(tag, "compatibility of all bars [m]", compat, 1e-15)
        check(tag, "equilibrium of joint D [N]", eqD, 1e-8)
        check(tag, "sum of reactions + load [N]", sumR, 1e-8)
        check(tag, "moment about B [N m]", mom, 1e-8)
        check(tag, "N force method - displacement method [N]", dN, 1e-8)
        check(tag, "u force method - displacement method [m]", du, 1e-16)
        out[cname] = {"X": X, "N": N, "u": u, "R": R, "e": e,
                      "compat": compat, "eqD": eqD, "sumR": sumR, "mom": mom,
                      "dN": dN, "du": du, "uK": uK}
        rows.append([cname, kN(X), kN(N[1]), kN(N[2]), kN(N[3]), mm(u[0]),
                     mm(u[1])])
    doc.eq("X1 = -(D1P + D1T + D1S)/f11,   N_i = n0_i + X1 n1_i",
           "u_D . e = sum_{1,3} n_e,i e_i - sum Rtilde_e,k . delta_k"
           "   (unit load e on the primary)",
           "R_k = -N_i c_i at the support of bar i")
    doc.table(["case", "X1 = N2 [kN]", "N1 [kN]", "N2 [kN]", "N3 [kN]",
               "ux_D [mm]", "uz_D [mm]"], rows)
    doc.table(["case", "RA_x", "RA_z", "RB_x", "RB_z", "RC_x", "RC_z"],
              [[cn] + [kN(o["R"][kk][a]) for kk in ("A", "B", "C")
                       for a in (0, 1)] for cn, o in out.items()])
    doc.p("Reactions in kN. The out-of-plane reaction at D and all y "
          "reactions are zero.")
    doc.h3("Step 5  Independent checks")
    doc.p("Compatibility of every bar (the elongation computed from the "
          "unit-load displacements equals N L/EA + alpha dT L, including the "
          "cut bar 2), equilibrium of joint D, global force and moment "
          "equilibrium, and a second solution by the displacement method "
          "(2 x 2 stiffness matrix of joint D solved by Cramer's rule).")
    doc.table(["case", "max compat. [m]", "joint D [N]", "sum R + P [N]",
               "moment [N m]", "|dN| vs stiffness [N]", "|du| [m]"],
              [[cn, g6(o["compat"], 2), g6(o["eqD"], 2), g6(o["sumR"], 2),
                g6(o["mom"], 2), g6(o["dN"], 2), g6(o["du"], 2)]
               for cn, o in out.items()])
    comb = sum(out[cn]["X"] for cn in ("load", "temperature", "settlement"))
    check("F4", "superposition of X1 [N]", comb - out["combined"]["X"], 1e-8)
    doc.items([f"Superposition: X1(P) + X1(T) + X1(S) = {kN(comb)} kN, "
               f"combined case {kN(out['combined']['X'])} kN.",
               f"Engineering note: heating bar 2 alone by {g6(DT)} degC "
               f"gives N2 = {kN(out['temperature']['N'][2], 4)} kN, and a "
               f"{mm(-F4_SETTLE[2], 2)} mm settlement of support "
               f"{F4_SETTLE_NODE} gives N2 = "
               f"{kN(out['settlement']['N'][2], 4)} kN, compared with "
               f"{kN(out['load']['N'][2], 4)} kN from the "
               f"{kN(-P[1], 3)} kN load. Self-straining loads are of the "
               f"same order as gravity loads in a stiff indeterminate "
               f"truss."])
    res.add("F4.f11", f11, "m/N", "f11 = sum n1_i^2 L_i/(E A_i)",
            "flexibility coefficient of the redundant N2")
    for cn, o in out.items():
        for i in F4_BARS:
            res.add(f"F4.N{i}_{cn}", o["N"][i], "N",
                    "N_i = n0_i + X1*n1_i, X1 = -(D1P+D1T+D1S)/f11",
                    f"force in bar {i} ({sup[i]}-D), case {cn}")
        res.add(f"F4.ux_{cn}", o["u"][0], "m",
                "unit-load method on the primary structure",
                f"x displacement of node D, case {cn}")
        res.add(f"F4.uz_{cn}", o["u"][1], "m",
                "unit-load method on the primary structure",
                f"z displacement of node D, case {cn}")
        for kk in ("A", "B", "C"):
            for a, lab in ((0, "x"), (1, "z")):
                res.add(f"F4.R{kk}{lab}_{cn}", o["R"][kk][a], "N",
                        "R_k = -N_i*c_i", f"reaction {lab} at support {kk}, "
                        f"case {cn}")


# =============================================================================
# F5  bar + spring
# =============================================================================
def case_f5(doc, res):
    sec = chs(SEC_MAIN[1], SEC_MAIN[2])
    EA = E * sec["A"]
    L, k, P = F5_L, F5_K, F5_P
    kb = EA / L
    doc.h2("F5  Bar and elastic spring support")
    doc.p(f"Purpose. Elastic (spring) supports. A {SEC_MAIN[0]} bar of length "
          f"L = {g6(L)} m lies on the x axis from node 1 (0,0,0) to node 2 "
          f"({g6(L)},0,0). Its axial stiffness is EA/L = {kN(kb)} kN/m. The "
          f"spring has k = {kN(k)} kN/m (GeoDomeLab's default support "
          f"spring) and acts in x. The y and z DOFs of both nodes are "
          f"restrained, so each case has one or two free DOFs. The spring "
          f"force on its node is -k u, reported as the reaction.")

    doc.h3("F5a  In series: spring at node 1, load at node 2")
    u1 = P / k
    u2 = P / k + P * L / EA
    Ns = P
    Rs = -k * u1
    keq = 1.0 / (1.0 / k + L / EA)
    doc.eq("node 2:  bar force N = P (statically determinate)",
           "node 1:  the spring carries the same force, k u1 = N",
           "flexibilities add:  u2 = P/k + P L/EA,  k_eq = 1/(1/k + L/EA)",
           f"  u1 = {mm(u1)} mm,  u2 = {mm(u2)} mm,  N = {kN(Ns)} kN",
           f"  spring reaction at node 1 = -k u1 = {kN(Rs)} kN,"
           f"  k_eq = {kN(keq)} kN/m")
    Wext = 0.5 * P * u2
    Uint = 0.5 * Ns ** 2 * L / EA + 0.5 * k * u1 ** 2
    check("F5a", "Clapeyron 1/2 P u2 - (bar + spring energy) [J]",
          Wext - Uint, 1e-12)
    check("F5a", "load + spring reaction [N]", P + Rs, 1e-9)
    doc.items([f"Energy: 1/2 P u2 = {g6(Wext)} J, 1/2 N^2 L/EA + 1/2 k u1^2 = "
               f"{g6(Uint)} J.",
               f"Equilibrium: P + spring reaction = {g6(P + Rs, 3)} N."])
    res.add("F5.series_u1", u1, "m", "u1 = P/k", "spring node displacement")
    res.add("F5.series_u2", u2, "m", "u2 = P/k + P*L/EA",
            "loaded node displacement, bar and spring in series")
    res.add("F5.series_N", Ns, "N", "N = P", "bar force, series")
    res.add("F5.series_spring_reaction", Rs, "N", "R = -k*u1 = -P",
            "spring reaction at node 1, series")

    doc.h3("F5b  In parallel: node 1 fixed, spring and load at node 2")
    u2p = P / (kb + k)
    Np = kb * u2p
    Rp = -k * u2p
    R1p = -Np
    doc.eq("node 2:  (EA/L + k) u2 = P",
           f"  u2 = {mm(u2p)} mm",
           f"  bar N = (EA/L) u2 = {kN(Np)} kN  "
           f"(share EA/L/(EA/L + k) = {g6(kb / (kb + k))})",
           f"  spring reaction = -k u2 = {kN(Rp)} kN,"
           f"  fixed reaction at node 1 = -N = {kN(R1p)} kN")
    check("F5b", "load + reactions [N]", P + Rp + R1p, 1e-9)
    check("F5b", "Clapeyron [J]", 0.5 * P * u2p
          - (0.5 * kb * u2p ** 2 + 0.5 * k * u2p ** 2), 1e-12)
    doc.items([f"Equilibrium: P + R1 + spring = {g6(P + Rp + R1p, 3)} N.",
               "Energy: 1/2 P u2 = 1/2 (EA/L) u2^2 + 1/2 k u2^2 holds to "
               f"{g6(0.5 * P * u2p - 0.5 * (kb + k) * u2p ** 2, 3)} J."])
    res.add("F5.parallel_u2", u2p, "m", "u2 = P/(EA/L + k)",
            "loaded node displacement, bar and spring in parallel")
    res.add("F5.parallel_N", Np, "N", "N = (EA/L)*u2", "bar force, parallel")
    res.add("F5.parallel_spring_reaction", Rp, "N", "R = -k*u2",
            "spring reaction at node 2, parallel")
    res.add("F5.parallel_R1x", R1p, "N", "R1 = -N",
            "fixed reaction at node 1, parallel")

    doc.h3(f"F5c  Heated bar on a spring: node 1 fixed, spring at node 2, "
           f"dT = +{g6(DT)} degC")
    u2t = ALPHA * DT * L / (1.0 + k * L / EA)
    Nt = -EA * ALPHA * DT * k / (k + kb)
    Rt = -k * u2t
    doc.eq("thermal equivalent load at node 2:  F0 = EA alpha dT",
           "(EA/L + k) u2 = EA alpha dT",
           "  u2 = alpha dT L / (1 + k L/EA)",
           "  N  = EA (u2/L - alpha dT) = -EA alpha dT k/(k + EA/L)",
           f"  u2 = {mm(u2t)} mm   (free expansion {mm(ALPHA * DT * L)} mm)",
           f"  N  = {kN(Nt)} kN    (fully fixed {kN(-EA * ALPHA * DT)} kN)",
           f"  spring reaction -k u2 = {kN(Rt)} kN = N")
    Nchk = EA * (u2t / L - ALPHA * DT)
    check("F5c", "N from u2 - closed form [N]", Nchk - Nt, 1e-8)
    check("F5c", "node 2 equilibrium -N + spring [N]", -Nt + Rt, 1e-8)
    big = 1e22                      # "rigid" spring, kb/big ~ 1e-15
    lim_inf = -EA * ALPHA * DT * big / (big + kb)
    check("F5c", "k -> infinity gives F1 [relative]",
          lim_inf / (-EA * ALPHA * DT) - 1.0, 1e-13)
    tiny = 1e-6                     # "no" spring
    check("F5c", "k -> 0 gives free expansion [relative]",
          (ALPHA * DT * L / (1.0 + tiny * L / EA)) / (ALPHA * DT * L) - 1.0,
          1e-13)
    doc.items([f"Node 2 equilibrium: the compressed bar pushes node 2 by "
               f"-N = {kN(-Nt)} kN and the spring pulls back by "
               f"{kN(Rt)} kN.",
               "Limits: k -> 0 gives the free expansion alpha dT L with "
               "N = 0, k -> infinity gives the fixed bar of F1 "
               f"(N -> -EA alpha dT = {kN(-EA * ALPHA * DT)} kN).",
               f"Engineering note: an elastic support with k comparable to "
               f"EA/L relieves {g6(100 * (1 - k / (k + kb)), 3)} % of the "
               f"fully restrained thermal force."])
    res.add("F5.thermal_u2", u2t, "m", "u2 = alpha*dT*L/(1 + k*L/EA)",
            "spring node displacement, heated bar")
    res.add("F5.thermal_N", Nt, "N", "N = -EA*alpha*dT*k/(k + EA/L)",
            "bar force, heated bar on a spring")
    res.add("F5.thermal_spring_reaction", Rt, "N", "R = -k*u2 = N",
            "spring reaction, heated bar on a spring")


# =============================================================================
# F6  von Mises two-bar truss
# =============================================================================
def vm_props():
    sec = chs(SEC_MAIN[1], SEC_MAIN[2])
    EA = E * sec["A"]
    a, h = VM_A, VM_H
    L = math.hypot(a, h)
    return a, h, L, EA, sec


def vm_state(w):
    """Exact state of the von Mises truss at a downward apex deflection w.

    y = h - w,  l = sqrt(a^2 + y^2),  eps = (l - L)/L,  N = EA eps.
    L - l is evaluated without cancellation as
        L - l = (L^2 - l^2)/(L + l) = w (2h - w)/(L + l)
    which is the same number, written so that no digits are lost when the
    bars are almost as long as L.
    """
    a, h, L, EA, _ = vm_props()
    y = h - w
    l = math.hypot(a, y)
    L_minus_l = w * (2.0 * h - w) / (L + l)
    eps = -L_minus_l / L
    N = EA * eps
    P = 2.0 * EA * L_minus_l / L * y / l
    return {"y": y, "l": l, "eps": eps, "N": N, "P": P}


def vm_P(w):
    """Downward apex load in equilibrium with a downward apex deflection w.

    P(w) = 2 EA (L - l)/L * (h - w)/l,   l = sqrt(a^2 + (h - w)^2)
    (equal to -F_z, F_z = 2 EA (l - L)/L (h - w)/l the upward apex force).
    """
    return vm_state(w)["P"]


def vm_dPdw(w):
    """Tangent dP/dw = 2 [(EA/L)(y/l)^2 + (N/l)(a/l)^2] = 2 EA (1/L - a^2/l^3)."""
    a, h, L, EA, _ = vm_props()
    s = vm_state(w)
    y, l, N = s["y"], s["l"], s["N"]
    return 2.0 * ((EA / L) * (y / l) ** 2 + (N / l) * (a / l) ** 2)


def case_f6(doc, res):
    a, h, L, EA, sec = vm_props()
    doc.h2("F6  von Mises two-bar truss: snap-through")
    doc.p(f"Purpose. The classic benchmark for geometrically nonlinear "
          f"(large displacement, small strain) truss analysis. Two "
          f"{SEC_MAIN[0]} bars run from pinned supports at (-a, 0, 0) and "
          f"(+a, 0, 0) to the apex (0, 0, h), with a = {g6(a)} m and "
          f"h = {g6(h)} m (h/a = {g6(h / a)}, a shallow truss). The apex is "
          f"restrained out of plane (uy = 0) and carries a downward load P. "
          f"Strain is the engineering strain (l - L)/L, the measure "
          f"GeoDomeLab's corotational bar uses. EA = {kN(EA)} kN.")
    doc.h3("Equilibrium path (exact, closed form)")
    K0 = 2.0 * EA * h * h / L ** 3
    doc.eq(f"L = sqrt(a^2 + h^2) = {g6(L, 10)} m",
           "apex deflection w (down), y = h - w, l = sqrt(a^2 + y^2)",
           "eps = (l - L)/L,  N = EA eps  (both bars, by symmetry)",
           "vertical equilibrium of the apex in the deformed position:",
           "  P = -2 N y/l",
           "  P(w) = 2 EA (L - l)/L * (h - w)/l",
           "  (the same as F_z = 2 EA (l - L)/L (h - w)/l with F_z positive "
           "up)",
           "  evaluated with L - l = (L^2 - l^2)/(L + l) = w (2h - w)/(L + l)",
           "  (identical algebra, no loss of digits when l is close to L)",
           "tangent  dP/dw = 2 EA (1/L - a^2/l^3)",
           "         = 2 [ (EA/L)(y/l)^2 + (N/l)(a/l)^2 ]   (material +"
           " geometric)",
           f"at w = 0  dP/dw = 2 EA h^2/L^3 = {kN(K0)} kN/m"
           "  (the linear stiffness)")
    K0_fd = (vm_P(1e-7) - vm_P(-1e-7)) / 2e-7
    check("F6", "dP/dw(0) analytic - finite difference [relative]",
          K0_fd / K0 - 1.0, 1e-6)
    check("F6", "dP/dw(0) - 2 EA h^2/L^3 [relative]", vm_dPdw(0.0) / K0 - 1.0,
          1e-12)

    doc.h3("Limit point from dP/dw = 0")
    # dP/dw = 0  <=>  l^3 = a^2 L.  With z = y^2/a^2 and r = h/a:
    #   (1 + z)^(3/2) = L/a = (1 + r^2)^(1/2)  <=>  (1 + z)^3 = 1 + r^2
    #   <=>  z^3 + 3 z^2 + 3 z - r^2 = 0     (no cancellation, p'(z) ~ 3)
    r2 = (h / a) ** 2
    zs = math.expm1(math.log1p(r2) / 3.0)          # (1 + r^2)^(1/3) - 1
    ys = a * math.sqrt(zs)
    ls = a * math.sqrt(1.0 + zs)
    ls_direct = (a * a * L) ** (1.0 / 3.0)
    ws = h - ys
    sts = vm_state(ws)
    Ps, Ns, eps_s = sts["P"], sts["N"], sts["eps"]
    # Newton on the cubic, started at z = r^2/3 (the shallow estimate)
    z = r2 / 3.0
    iters = []
    for it in range(50):
        p = ((z + 3.0) * z + 3.0) * z - r2
        dp = (3.0 * z + 6.0) * z + 3.0
        dz = -p / dp
        z += dz
        iters.append((it + 1, z, p))
        if abs(dz) <= 2.0 * sys.float_info.epsilon * abs(z):
            break
    yN = a * math.sqrt(z)
    d2 = -6.0 * EA * a * a * ys / ls ** 5
    Kuu = 2.0 * ((EA / L) * (a / ls) ** 2 + (Ns / ls) * (ys / ls) ** 2)
    Ps_shallow = 2.0 * EA * h ** 3 / (3.0 * math.sqrt(3.0) * a ** 3)
    doc.eq("dP/dw = 0  <=>  a^2/l^3 = 1/L  <=>  l* = (a^2 L)^(1/3)",
           "written without cancellation, z = y^2/a^2, r = h/a:",
           "  (1 + z)^3 = 1 + r^2   <=>   z^3 + 3 z^2 + 3 z - r^2 = 0",
           "  z* = (1 + r^2)^(1/3) - 1 = expm1(log1p(r^2)/3)",
           f"  z* = {g6(zs, 16)}",
           f"  y* = a sqrt(z*) = {g6(ys, 16)} m",
           f"  l* = a sqrt(1 + z*) = {g6(ls, 16)} m"
           f"   ((a^2 L)^(1/3) = {g6(ls_direct, 16)} m)",
           f"  w* = h - y* = {g6(ws, 16)} m  = {g6(ws / h, 10)} h",
           f"  P* = P(w*) = {g6(Ps, 13)} N = {kN(Ps)} kN",
           f"  eps* = (l* - L)/L = {g6(eps_s)},  N* = EA eps* = {kN(Ns)} kN"
           f",  sigma* = {g6(Ns / sec['A'] / 1e6, 4)} MPa",
           "Newton iteration on p(z) = z^3 + 3z^2 + 3z - r^2 from z = r^2/3:",
           *[f"  it {i:2d}: z = {zz:.17g}   p = {pp:+.3e}"
             for i, zz, pp in iters],
           f"  Newton y* - closed form y* = {g6(yN - ys, 3)} m"
           f"  (relative {g6((yN - ys) / ys, 3)})",
           "second derivative  d2P/dw2 = -6 EA a^2 y/l^5 < 0  (a maximum)",
           f"  d2P/dw2(w*) = {g6(d2)} N/m^2",
           "shallow-truss approximation  P* ~ 2 EA h^3/(3 sqrt3 a^3)"
           f" = {kN(Ps_shallow)} kN,  w* ~ h (1 - 1/sqrt3) = "
           f"{g6(h * (1 - 1 / math.sqrt(3)))} m")
    check("F6", "Newton y* - closed form y* [relative]", (yN - ys) / ys,
          4.0 * sys.float_info.epsilon)
    check("F6", "l* = a sqrt(1+z*) vs (a^2 L)^(1/3) [relative]",
          ls / ls_direct - 1.0, 1e-14)
    check("F6", "dP/dw at the limit point [relative to K0]",
          vm_dPdw(ws) / K0, 1e-12)
    check("F6", "2EA(1/L - a^2/l*^3) at the limit point [relative to K0]",
          2.0 * EA * (1.0 / L - a * a / ls ** 3) / K0, 1e-12)
    check("F6", "P* = -2 N* y*/l* [N]", Ps + 2.0 * Ns * ys / ls, 1e-7)

    doc.h3("Snap-through target and the full path")
    w2 = h + ys                     # second (lower) limit point, P = -P*
    lo, hi = 2.0 * h, 2.0 * h + h
    while vm_P(hi) < Ps:
        hi += h
    nb = 0
    while hi - lo > 1e-15 * hi and nb < 200:
        mid = 0.5 * (lo + hi)
        if vm_P(mid) < Ps:
            lo = mid
        else:
            hi = mid
        nb += 1
    wsn = 0.5 * (lo + hi)
    N_sn = vm_state(wsn)["N"]
    doc.eq("Under load control the apex jumps at P* from w* to the stable",
           "inverted branch (w > 2h) where P(w) = P* again. Bisection on",
           "P(w) - P* over [2h, w_hi]:",
           f"  w_snap = {g6(wsn, 12)} m  = {g6(wsn / h, 8)} h"
           f"  ({nb} bisections)",
           f"  bar force there N = {kN(N_sn)} kN (tension, the bars are",
           "  longer than L in the inverted position)",
           f"second limit point  w = h + y* = {g6(w2)} m, P = -P*",
           "zero-load points  w = 0, w = h (bars horizontal), w = 2h")
    check("F6", "P(w_snap) - P* [N]", vm_P(wsn) - Ps, 1e-6)
    check("F6", "P(h) [N]", vm_P(h), 1e-9)
    check("F6", "P(2h) [N]", vm_P(2 * h), 1e-6)
    check("F6", "P(h + y*) + P* [N]", vm_P(w2) + Ps, 1e-6)
    rows = []
    for r in VM_SAMPLES:
        w = r * h
        sw = vm_state(w)
        rows.append([g6(r), mm(w), g6(sw["eps"]), kN(sw["N"]),
                     kN(sw["P"]), g6(sw["P"] / Ps)])
        res.add(f"F6.P_at_w{r:g}h", vm_P(w), "N",
                "P = 2EA(L-l)/L*(h-w)/l", f"load on the path at w = {r:g} h")
    doc.table(["w/h", "w [mm]", "strain", "N [kN]", "P [kN]", "P/P*"], rows)

    doc.h3("Linearised buckling factor (2-DOF eigenproblem)")
    Pref = VM_PREF
    Nref = -Pref * L / (2.0 * h)
    Kuu0 = 2.0 * EA * a * a / L ** 3
    Kww0 = 2.0 * EA * h * h / L ** 3
    Gu = 2.0 * Nref * h * h / L ** 3
    Gw = 2.0 * Nref * a * a / L ** 3
    lam_u = -Kuu0 / Gu
    lam_w = -Kww0 / Gw
    lam = min(lam_u, lam_w)
    Pcr = lam * Pref
    Pcr_cf = 2.0 * EA * h ** 3 / (a * a * L)
    # characteristic polynomial of the diagonal pencil, checked at lambda
    det_w = (Kuu0 + lam * Gu) * (Kww0 + lam * Gw)
    doc.eq(f"reference load P_ref = {kN(Pref)} kN downward, linear member "
           f"force",
           f"  N_ref = -P_ref L/(2h) = {kN(Nref)} kN",
           "apex DOFs (u horizontal, w vertical), bars c = (+-a, h)/L:",
           "  K  = (EA/L) sum c c^T        = diag(2EA a^2/L^3, 2EA h^2/L^3)",
           "  Kg = (N_ref/L) sum (I - c c^T) = diag(2N_ref h^2/L^3,"
           " 2N_ref a^2/L^3)",
           f"  K  = diag({kN(Kuu0)}, {kN(Kww0)}) kN/m",
           f"  Kg = diag({kN(Gu)}, {kN(Gw)}) kN/m",
           "det(K + lambda Kg) = (Kuu + lambda Guu)(Kww + lambda Gww) = 0",
           "  lambda_u = -Kuu/Guu = EA a^2/(-N_ref h^2)   (sway mode (1,0))",
           "  lambda_w = -Kww/Gww = EA h^2/(-N_ref a^2)   (symmetric mode (0,1))",
           f"  lambda_u = {g6(lam_u, 10)},  lambda_w = {g6(lam_w, 10)}",
           f"  lambda_cr = {g6(lam, 10)}  ->  P_cr,lin = lambda_cr P_ref = "
           f"{kN(Pcr)} kN",
           f"  closed form P_cr,lin = 2 EA h^3/(a^2 L) = {kN(Pcr_cf)} kN",
           f"  det(K + lambda_cr Kg) = {g6(det_w, 3)}")
    check("F6", "lambda_cr P_ref - 2EA h^3/(a^2 L) [relative]",
          Pcr / Pcr_cf - 1.0, 1e-13)
    check("F6", "no sway bifurcation before the limit (Kuu(w*) > 0)",
          0.0 if Kuu > 0 else 1.0, 0.5)
    ratio = Pcr / Ps
    doc.items([f"The exact limit load is P* = {kN(Ps, 6)} kN. The linearised "
               f"eigenvalue gives P_cr,lin = {kN(Pcr, 6)} kN, "
               f"{g6(ratio, 5)} times larger. For a shallow truss the ratio "
               f"tends to 3 sqrt3 = {g6(3 * math.sqrt(3), 5)}. Linear "
               f"buckling ignores the change of geometry before the limit "
               f"point and is unconservative for shallow lattices, which is "
               f"why GeoDomeLab 2.0 adds the nonlinear path.",
               f"Sway stiffness at the limit point Kuu(w*) = "
               f"{kN(Kuu)} kN/m > 0, so the apex reaches the symmetric limit "
               f"point before any sway bifurcation.",
               f"At the limit point the bar stress is "
               f"{g6(Ns / sec['A'] / 1e6, 4)} MPa < Fy = {g6(FY / 1e6)} MPa "
               f"and |N*| = {kN(abs(Ns), 4)} kN is below the Euler load of "
               f"one bar, {kN(math.pi ** 2 * E * sec['I'] / L ** 2, 4)} kN, "
               f"so the elastic two-bar model is physically meaningful."])
    res.add("F6.L", L, "m", "L = sqrt(a^2 + h^2)", "initial bar length")
    res.add("F6.EA", EA, "N", "EA", "axial rigidity of each bar")
    res.add("F6.K0", K0, "N/m", "K0 = 2EA h^2/L^3",
            "initial vertical stiffness of the apex")
    res.add("F6.l_limit", ls, "m", "l* = (a^2 L)^(1/3)",
            "bar length at the limit point")
    res.add("F6.w_limit", ws, "m", "w* = h - sqrt(l*^2 - a^2)",
            "apex deflection at the limit point")
    res.add("F6.P_limit", Ps, "N", "P* = 2EA(L-l*)/L*(h-w*)/l*",
            "limit (snap-through) load")
    res.add("F6.N_limit", Ns, "N", "N* = EA(l*-L)/L",
            "bar force at the limit point")
    res.add("F6.strain_limit", eps_s, "-", "eps* = (l*-L)/L",
            "engineering strain at the limit point")
    res.add("F6.w_snap", wsn, "m", "P(w_snap) = P*, w_snap > 2h",
            "apex deflection on the stable inverted branch at P*")
    res.add("F6.N_snap", N_sn, "N", "N = EA(l-L)/L at w_snap",
            "bar force after snap-through at P*")
    res.add("F6.w_second_limit", w2, "m", "w = h + y*",
            "apex deflection at the second limit point (P = -P*)")
    res.add("F6.P_ref", Pref, "N", "stated", "reference load of the "
            "eigenproblem")
    res.add("F6.N_ref", Nref, "N", "N_ref = -P_ref*L/(2h)",
            "linear bar force under P_ref")
    res.add("F6.lambda_cr", lam, "-", "lambda = EA h^2/(-N_ref a^2)",
            "linearised buckling factor at P_ref (symmetric mode)")
    res.add("F6.lambda_sway", lam_u, "-", "lambda = EA a^2/(-N_ref h^2)",
            "linearised buckling factor of the sway mode")
    res.add("F6.P_cr_linear", Pcr, "N", "P_cr = 2EA h^3/(a^2 L)",
            "linearised critical load")
    res.add("F6.ratio_linear_to_limit", ratio, "-", "P_cr,lin / P*",
            "overestimate of linear buckling (-> 3 sqrt3 when shallow)")


# =============================================================================
# F7  Model A crown natural frequencies
# =============================================================================
def case_f7(doc, res, MA):
    s, t, L = MA["s"], MA["t"], MA["L"]
    A = MA["section"]["A"]
    Kff = MA["Kff"]
    doc.h2("F7  Model A crown natural frequencies (lumped mass)")
    doc.p("Purpose. Free vibration. With pinned ring nodes only the crown "
          "moves, so Model A has three dynamic DOFs and the mass matrix is "
          "m I (the same lumped mass in x, y and z). Lumping rule: each strut "
          "puts half of its mass rho A L at each end node. Optionally the "
          "superimposed dead load (cladding pressure q times the surface "
          "tributary area, divided by g) is added as mass.")
    m_sw = 5.0 * 0.5 * RHO * A * L
    face = math.sqrt(3.0) / 4.0 * L * L
    trib = 5.0 / 3.0 * face
    m_dead = hc.Q_DEAD * trib / G
    doc.eq(f"rho A L = {g6(RHO)} x {g6(A)} x {g6(L)} = {g6(RHO * A * L)} kg"
           " per strut",
           f"self weight at the crown  m_sw = 5 x 1/2 rho A L = {g6(m_sw)} kg",
           f"crown surface tributary area = 5/3 x (sqrt3/4) L^2 = {g6(trib)}"
           " m2",
           f"dead load mass  m_q = q A_trib / g = {g6(hc.Q_DEAD)} x "
           f"{g6(trib)} / {g6(G)} = {g6(m_dead)} kg",
           "stiffness  Kxx = Kyy = 5/2 k s^2,  Kzz = 5 k t^2",
           f"  Kxx = {kN(Kff[0, 0])} kN/m,  Kzz = {kN(Kff[2, 2])} kN/m",
           "(K - omega^2 m I) phi = 0 with diagonal K:",
           "  omega_x = omega_y = sqrt(Kxx/m),  omega_z = sqrt(Kzz/m)",
           "  f = omega/(2 pi),  T = 1/f,  mass-normalised phi = e/sqrt(m)",
           "  f_z/f_xy = sqrt(2) t/s, independent of the mass")
    rows = []
    for tag, m in (("sw", m_sw), ("dead", m_sw + m_dead)):
        wx = math.sqrt(Kff[0, 0] / m)
        wz = math.sqrt(Kff[2, 2] / m)
        fx, fz = wx / (2 * math.pi), wz / (2 * math.pi)
        dst = m * G / Kff[2, 2]
        fz_st = math.sqrt(G / dst) / (2 * math.pi)
        # Rayleigh quotient with the static deflected shape under 1 g in x
        ux = m * G / Kff[0, 0]
        ray = (0.5 * Kff[0, 0] * ux ** 2) / (0.5 * m * ux ** 2)
        check(f"F7 {tag}", "f_z vs static deflection sqrt(g/d_st) [rel]",
              fz_st / fz - 1.0, 1e-14)
        check(f"F7 {tag}", "Rayleigh quotient vs omega_x^2 [rel]",
              ray / wx ** 2 - 1.0, 1e-14)
        check(f"F7 {tag}", "f_z/f_xy - sqrt2 t/s [-]",
              fz / fx - math.sqrt(2.0) * t / s, 1e-13)
        lab = "self weight only" if tag == "sw" else "self weight + dead"
        rows.append([lab, g6(m), g6(fx), g6(fx), g6(fz), g6(1 / fz),
                     mm(dst)])
        res.add(f"F7.m_crown_{tag}", m, "kg",
                "m = 5*rho*A*L/2" + (" + q*A_trib/g" if tag == "dead" else ""),
                f"lumped crown mass, {lab}")
        res.add(f"F7.f_xy_{tag}", fx, "Hz", "f = sqrt(5/2 k s^2/m)/(2 pi)",
                f"horizontal crown frequency (repeated pair), {lab}")
        res.add(f"F7.f_z_{tag}", fz, "Hz", "f = sqrt(5 k t^2/m)/(2 pi)",
                f"vertical crown frequency, {lab}")
        res.add(f"F7.T_z_{tag}", 1 / fz, "s", "T = 1/f_z",
                f"vertical crown period, {lab}")
        res.add(f"F7.omega_xy_{tag}", wx, "rad/s", "omega = sqrt(Kxx/m)",
                f"horizontal circular frequency, {lab}")
        res.add(f"F7.omega_z_{tag}", wz, "rad/s", "omega = sqrt(Kzz/m)",
                f"vertical circular frequency, {lab}")
        res.add(f"F7.phi_{tag}", 1 / math.sqrt(m), "1/sqrt(kg)",
                "phi = 1/sqrt(m)",
                f"mass-normalised crown mode amplitude, {lab}")
    doc.table(["mass", "m [kg]", "f_x [Hz]", "f_y [Hz]", "f_z [Hz]",
               "T_z [s]", "d_st [mm]"], rows)
    doc.items(["The x and y modes form a repeated pair: the five-fold "
               "symmetric cap is equally stiff in every horizontal "
               "direction (Kxx = Kyy, Kxy = 0).",
               "Check 1: the vertical frequency equals the single-DOF "
               "static-deflection formula f = sqrt(g/d_st)/(2 pi) with "
               "d_st = m g/Kzz.",
               "Check 2: the Rayleigh quotient of the static deflected shape "
               "under a 1 g horizontal load equals omega_x^2 (exact here "
               "because the shape is the mode).",
               "Check 3: each mode is a pure crown translation, so its "
               "effective modal mass is the whole free mass m in its own "
               "direction and 0 in the others (100 % participation). The "
               "modes are orthogonal with respect to M and K because both "
               "matrices are diagonal."])
    for tag in ("sw", "dead"):
        for d in "xyz":
            res.add(f"F7.meff_{d}_{tag}",
                    res[f"F7.m_crown_{tag}"]["value"], "kg",
                    "m_eff = (phi^T M r)^2 = m",
                    f"effective modal mass of the {d} mode, "
                    f"{'self weight only' if tag == 'sw' else 'self weight + dead'}")


# =============================================================================
# F8  2V Class I geodesic sphere
# =============================================================================
def icosahedron_unit():
    """Unit icosahedron in pole orientation, faces found from adjacency."""
    zr, rr = 1.0 / SQ5, 2.0 / SQ5
    V = [np.array([0.0, 0.0, 1.0])]
    for i in range(5):
        th = 72.0 * i * DEG
        V.append(np.array([rr * math.cos(th), rr * math.sin(th), zr]))
    for i in range(5):
        th = (72.0 * i + 36.0) * DEG
        V.append(np.array([rr * math.cos(th), rr * math.sin(th), -zr]))
    V.append(np.array([0.0, 0.0, -1.0]))
    V = np.array(V)
    n = len(V)
    dist = [[float(np.linalg.norm(V[i] - V[j])) for j in range(n)]
            for i in range(n)]
    e = min(dist[i][j] for i in range(n) for j in range(n) if i != j)
    adj = [[abs(dist[i][j] - e) < 1e-9 and i != j for j in range(n)]
           for i in range(n)]
    faces = []
    for i in range(n):
        for j in range(i + 1, n):
            for kk in range(j + 1, n):
                if adj[i][j] and adj[j][kk] and adj[i][kk]:
                    f = [i, j, kk]
                    nrm = np.cross(V[j] - V[i], V[kk] - V[i])
                    if nrm @ (V[i] + V[j] + V[kk]) < 0:
                        f = [i, kk, j]
                    faces.append(f)
    return V, faces, e


def geodesic_2v():
    """2V Class I: split each icosahedron face in 4 at the edge midpoints."""
    V, F, e1v = icosahedron_unit()
    pts = [p for p in V]
    mid = {}

    def m(i, j):
        key = (min(i, j), max(i, j))
        if key not in mid:
            p = V[i] + V[j]
            pts.append(p / np.linalg.norm(p))
            mid[key] = len(pts) - 1
        return mid[key]

    faces = []
    for a, b, c in F:
        ab, bc, ca = m(a, b), m(b, c), m(c, a)
        faces += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
    return np.array(pts), faces, V, F, e1v


def face_normal(P, f):
    a, b, c = (P[i] for i in f)
    n = np.cross(b - a, c - a)
    return n / np.linalg.norm(n)


def dihedral_formula(R, Ledge, tri1, tri2):
    """Interior dihedral angle across an edge of length Ledge shared by two
    triangles inscribed in a sphere of radius R (sides tri1, tri2), using the
    face circumradius rho = abc/(4 Area):
        d_f   = sqrt(R^2 - rho^2)        distance centre -> face plane
        h_e   = sqrt(R^2 - (L/2)^2)      distance centre -> edge midpoint
        delta = asin(d_1/h_e) + asin(d_2/h_e)
    valid when each circumcentre lies on its own face's side of the edge
    (true for acute triangles)."""
    def rho(sides):
        a, b, c = sides
        s = 0.5 * (a + b + c)
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        return a * b * c / (4.0 * area), area
    r1, _ = rho(tri1)
    r2, _ = rho(tri2)
    he = math.sqrt(R * R - 0.25 * Ledge * Ledge)
    d1 = math.sqrt(R * R - r1 * r1)
    d2 = math.sqrt(R * R - r2 * r2)
    return math.asin(d1 / he) + math.asin(d2 / he), (r1, r2, d1, d2, he)


def case_f8(doc, res):
    doc.h2("F8  2V Class I geodesic geometry: chord factors and angles")
    doc.p("Purpose. Fabrication data. The 2V Class I sphere splits each of "
          "the 20 icosahedron faces into four triangles through the edge "
          "midpoints, which are pushed out radially to the sphere. Every "
          "value below is computed twice, once from closed-form spherical "
          "trigonometry and once from the vertex vectors of an independently "
          "built 2V sphere, and the chord factors are compared with the "
          "published 2V chord factors of Kenner (1976), Geodesic Math and "
          "How to Use It: A = 0.54653, B = 0.61803.")
    P, faces, V, F, e1 = geodesic_2v()
    R = 1.0
    doc.h3("Closed form")
    c_ico = 1.0 / SQ5
    th_ico = math.acos(c_ico)
    thA = 0.5 * th_ico
    cosB = (1.0 + 3.0 * c_ico) / (2.0 + 2.0 * c_ico)
    thB = math.acos(cosB)
    cfA = 2.0 * math.sin(0.5 * thA)
    cfB = 2.0 * math.sin(0.5 * thB)
    axA = math.acos(cfA / 2.0)
    axB = math.acos(cfB / 2.0)
    doc.eq("icosahedron: central angle of an edge  cos th_ico = 1/sqrt5",
           f"  th_ico = {g6(th_ico / DEG, 12)} deg,  1V chord = "
           f"2 sin(th_ico/2) = {g6(2 * math.sin(th_ico / 2), 10)}",
           "strut A (5-valent vertex -> edge midpoint):  th_A = th_ico/2",
           f"  th_A = {g6(thA / DEG, 12)} deg,  CF_A = 2 sin(th_A/2) = "
           f"{g6(cfA, 12)}",
           "strut B (midpoint -> midpoint in one face):",
           "  cos th_B = (v1+v2).(v1+v3)/|v1+v2|^2 = (1 + 3c)/(2 + 2c),"
           " c = 1/sqrt5",
           f"         = {g6(cosB, 15)} = cos 36 deg = (1 + sqrt5)/4",
           f"  CF_B = 2 sin 18 deg = (sqrt5 - 1)/2 = {g6(cfB, 12)}",
           "axial angle (strut to the radius at its hub):",
           "  (p_j - p_i).(-p_i) = R^2 - p_i.p_j = L^2/2  ->  cos alpha = L/(2R)",
           f"  alpha_A = acos(CF_A/2) = {g6(axA / DEG, 10)} deg,"
           f"  alpha_B = acos(CF_B/2) = {g6(axB / DEG, 10)} deg",
           "  (the end-cut or bevel angle to the tangent plane is"
           " 90 - alpha = th/2)")
    check("F8", "CF_B - (sqrt5 - 1)/2", cfB - (SQ5 - 1) / 2, 1e-15)
    check("F8", "cos th_B - cos 36", cosB - math.cos(36 * DEG), 1e-15)

    # face angles and areas
    apex = math.acos(1.0 - cfB ** 2 / (2.0 * cfA ** 2))
    base = 0.5 * (math.pi - apex)
    area_AAB = 0.25 * cfB * math.sqrt(4.0 * cfA ** 2 - cfB ** 2)
    area_BBB = math.sqrt(3.0) / 4.0 * cfB ** 2
    dA, infoA = dihedral_formula(R, cfA, (cfA, cfA, cfB), (cfA, cfA, cfB))
    dB, infoB = dihedral_formula(R, cfB, (cfA, cfA, cfB), (cfB, cfB, cfB))
    dI, _ = dihedral_formula(R, e1, (e1, e1, e1), (e1, e1, e1))
    dI_exact = math.acos(-SQ5 / 3.0)
    doc.eq("panels: AAB (isosceles, 60 of them) and BBB (equilateral, 20)",
           f"  AAB apex angle (at the 5-valent hub) = acos(1 - B^2/(2A^2))"
           f" = {g6(apex / DEG, 10)} deg",
           f"  AAB base angles = {g6(base / DEG, 10)} deg,  BBB = 60 deg",
           f"  areas (R = 1): AAB = (B/4) sqrt(4A^2 - B^2) = {g6(area_AAB)},"
           f"  BBB = (sqrt3/4) B^2 = {g6(area_BBB)}",
           "dihedral angle across an edge (face circumradius rho = abc/(4"
           " Area)):",
           "  d_f = sqrt(R^2 - rho^2),  h_e = sqrt(R^2 - (L/2)^2)",
           "  delta = asin(d_1/h_e) + asin(d_2/h_e)",
           f"  A edge (AAB | AAB): rho = {g6(infoA[0])}, d = {g6(infoA[2])},"
           f" h_e = {g6(infoA[4])}",
           f"    delta_A = {g6(dA / DEG, 10)} deg,  bend = 180 - delta = "
           f"{g6(180 - dA / DEG, 10)} deg",
           f"  B edge (AAB | BBB): rho = {g6(infoB[0])}, {g6(infoB[1])},"
           f" d = {g6(infoB[2])}, {g6(infoB[3])}, h_e = {g6(infoB[4])}",
           f"    delta_B = {g6(dB / DEG, 10)} deg,  bend = "
           f"{g6(180 - dB / DEG, 10)} deg",
           f"  anchor: the same formula for the icosahedron gives "
           f"{g6(dI / DEG, 12)} deg,",
           f"  the exact value acos(-sqrt5/3) = {g6(dI_exact / DEG, 12)} deg")
    check("F8", "dihedral formula on the icosahedron - acos(-sqrt5/3) [rad]",
          dI - dI_exact, 1e-13)
    check("F8", "AAB triangle acute (circumcentre inside)",
          0.0 if apex < math.pi / 2 and base < math.pi / 2 else 1.0, 0.5)

    doc.h3("From the vertex vectors of an independently built 2V sphere")
    edges = {}
    for fi, f in enumerate(faces):
        for i, j in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            edges.setdefault((min(i, j), max(i, j)), []).append(fi)
    outward = all(face_normal(P, f) @ (P[f[0]] + P[f[1]] + P[f[2]]) > 0
                  for f in faces)
    lens = {e: float(np.linalg.norm(P[e[1]] - P[e[0]])) for e in edges}
    fam = {}
    for e, Le in lens.items():
        fam[e] = "A" if abs(Le - cfA) < 1e-9 else (
            "B" if abs(Le - cfB) < 1e-9 else "?")
    nA = sum(1 for e in fam if fam[e] == "A")
    nB = sum(1 for e in fam if fam[e] == "B")
    nU = sum(1 for e in fam if fam[e] == "?")
    val = {}
    hubtype = {}
    for (i, j) in edges:
        for n_ in (i, j):
            val[n_] = val.get(n_, 0) + 1
            hubtype.setdefault(n_, []).append(fam[(i, j)])
    n5 = sum(1 for v in val.values() if v == 5)
    n6 = sum(1 for v in val.values() if v == 6)
    hub_kinds = {}
    for n_, lst in hubtype.items():
        key = f"{lst.count('A')}A+{lst.count('B')}B"
        hub_kinds[key] = hub_kinds.get(key, 0) + 1
    ftype = {}
    for f in faces:
        key = "".join(sorted(fam[(min(f[i], f[(i + 1) % 3]),
                                  max(f[i], f[(i + 1) % 3]))]
                             for i in range(3)))
        ftype[key] = ftype.get(key, 0) + 1
    # axial angles from vectors (both ends of every strut)
    ax_err = {"A": 0.0, "B": 0.0}
    for (i, j), fm in fam.items():
        for p, q in ((i, j), (j, i)):
            dvec = P[q] - P[p]
            cosang = float(dvec @ (-P[p])) / (np.linalg.norm(dvec)
                                              * np.linalg.norm(P[p]))
            ang = math.acos(cosang)
            ref = axA if fm == "A" else axB
            ax_err[fm] = max(ax_err[fm], abs(ang - ref))
    # dihedral and bend angles from normals
    dih = {"A": [], "B": []}
    for e, fl_ in edges.items():
        n1 = face_normal(P, faces[fl_[0]])
        n2 = face_normal(P, faces[fl_[1]])
        bend = math.acos(max(-1.0, min(1.0, float(n1 @ n2))))
        dih[fam[e]].append(math.pi - bend)
    dvA = float(np.mean(dih["A"]))
    dvB = float(np.mean(dih["B"]))
    sprA = max(dih["A"]) - min(dih["A"])
    sprB = max(dih["B"]) - min(dih["B"])
    two_faces = all(len(v) == 2 for v in edges.values())
    nV, nE, nF = len(P), len(edges), len(faces)
    on_sphere = max(abs(float(np.linalg.norm(p)) - 1.0) for p in P)
    tot_area = sum(0.5 * float(np.linalg.norm(np.cross(P[f[1]] - P[f[0]],
                                                       P[f[2]] - P[f[0]])))
                   for f in faces)
    tot_cf = 60 * area_AAB + 20 * area_BBB
    doc.table(["quantity", "vectors", "closed form", "difference"],
              [["nodes V", nV, "10 v^2 + 2 = 42", nV - 42],
               ["struts E", nE, "30 v^2 = 120", nE - 120],
               ["panels F", nF, "20 v^2 = 80", nF - 80],
               ["Euler V - E + F", nV - nE + nF, "2", nV - nE + nF - 2],
               ["5-valent hubs", n5, "12", n5 - 12],
               ["6-valent hubs", n6, "30", n6 - 30],
               ["A struts", nA, "60", nA - 60],
               ["B struts", nB, "60", nB - 60],
               ["CF_A", g6(float(np.mean([lens[e] for e in fam
                                          if fam[e] == 'A'])), 12),
                g6(cfA, 12),
                g6(float(np.mean([lens[e] for e in fam if fam[e] == 'A']))
                   - cfA, 3)],
               ["CF_B", g6(float(np.mean([lens[e] for e in fam
                                          if fam[e] == 'B'])), 12),
                g6(cfB, 12),
                g6(float(np.mean([lens[e] for e in fam if fam[e] == 'B']))
                   - cfB, 3)],
               ["axial A [deg] (max dev.)", g6(axA / DEG, 10), g6(axA / DEG, 10),
                g6(ax_err["A"] / DEG, 3)],
               ["axial B [deg] (max dev.)", g6(axB / DEG, 10), g6(axB / DEG, 10),
                g6(ax_err["B"] / DEG, 3)],
               ["dihedral A [deg]", g6(dvA / DEG, 10), g6(dA / DEG, 10),
                g6((dvA - dA) / DEG, 3)],
               ["dihedral B [deg]", g6(dvB / DEG, 10), g6(dB / DEG, 10),
                g6((dvB - dB) / DEG, 3)],
               ["panel area sum (R = 1)", g6(tot_area, 12), g6(tot_cf, 12),
                g6(tot_area - tot_cf, 3)]])
    doc.p(f"Hub types: {', '.join(f'{v} x {k}' for k, v in sorted(hub_kinds.items()))}. "
          f"Panel types: {', '.join(f'{v} x {k}' for k, v in sorted(ftype.items()))}. "
          f"Every strut has exactly two panels: {two_faces}. All panels are "
          f"oriented outwards: {outward}. Largest distance of a node from "
          f"the unit sphere: {g6(on_sphere, 3)}. Spread of the dihedral "
          f"angle inside each family: A {g6(sprA / DEG, 3)} deg, "
          f"B {g6(sprB / DEG, 3)} deg. Faceted area / sphere area = "
          f"{g6(tot_area / (4 * math.pi), 8)}.")
    check("F8", "V, E, F counts", abs(nV - 42) + abs(nE - 120) + abs(nF - 80),
          0)
    check("F8", "unclassified struts", nU, 0)
    check("F8", "valences 12 x 5 and 30 x 6", abs(n5 - 12) + abs(n6 - 30), 0)
    check("F8", "A and B strut counts 60/60", abs(nA - 60) + abs(nB - 60), 0)
    check("F8", "axial angle vectors - acos(L/2R) [rad]",
          max(ax_err.values()), 1e-12)
    check("F8", "dihedral A vectors - formula [rad]", dvA - dA, 1e-12)
    check("F8", "dihedral B vectors - formula [rad]", dvB - dB, 1e-12)
    check("F8", "dihedral spread inside a family [rad]", max(sprA, sprB),
          1e-12)
    check("F8", "panel area sum vectors - formula", tot_area - tot_cf, 1e-12)
    check("F8", "CF_A vs Kenner 0.54653 (5 decimals)",
          cfA - KENNER_2V["A"], 5e-6)
    check("F8", "CF_B vs Kenner 0.61803 (5 decimals)",
          cfB - KENNER_2V["B"], 5e-6)
    check("F8", "every strut shared by two panels",
          0.0 if two_faces else 1.0, 0.5)
    doc.items([f"Kenner (1976): A = {KENNER_2V['A']}, B = {KENNER_2V['B']}. "
               f"Computed: A = {cfA:.8f}, B = {cfB:.8f}. Both agree to the "
               f"five published decimals (differences "
               f"{g6(cfA - KENNER_2V['A'], 3)} and "
               f"{g6(cfB - KENNER_2V['B'], 3)}).",
               f"For a sphere of R = {g6(F8_R)} m the struts are "
               f"{g6(cfA * F8_R, 6)} m (A) and {g6(cfB * F8_R, 6)} m (B).",
               "The axial angle is the same at both ends of a strut because "
               "both hubs lie on the sphere. The dihedral angle is the angle "
               "between the two panels measured inside the dome. The bend "
               "angle (180 - dihedral) is the angle between the panel "
               "normals and is what a hub or panel bracket must turn."])

    res.add("F8.CF_A", cfA, "-", "CF_A = 2 sin(acos(1/sqrt5)/4)",
            "2V chord factor of strut A (Kenner 0.54653)")
    res.add("F8.CF_B", cfB, "-", "CF_B = 2 sin 18 deg = (sqrt5-1)/2",
            "2V chord factor of strut B (Kenner 0.61803)")
    res.add("F8.L_A_R5", cfA * F8_R, "m", "L = CF_A*R",
            f"strut A length for R = {g6(F8_R)} m")
    res.add("F8.L_B_R5", cfB * F8_R, "m", "L = CF_B*R",
            f"strut B length for R = {g6(F8_R)} m")
    res.add("F8.axial_A_deg", axA / DEG, "deg", "alpha = acos(L/(2R))",
            "axial angle of strut A (strut to hub radius)")
    res.add("F8.axial_B_deg", axB / DEG, "deg", "alpha = acos(L/(2R))",
            "axial angle of strut B (strut to hub radius)")
    res.add("F8.bevel_A_deg", 90 - axA / DEG, "deg", "90 - alpha = th_A/2",
            "strut A end angle to the tangent plane")
    res.add("F8.bevel_B_deg", 90 - axB / DEG, "deg", "90 - alpha = th_B/2",
            "strut B end angle to the tangent plane")
    res.add("F8.dihedral_A_deg", dA / DEG, "deg",
            "asin(d1/h)+asin(d2/h), d = sqrt(R^2-rho^2)",
            "interior dihedral angle across A struts (AAB | AAB)")
    res.add("F8.dihedral_B_deg", dB / DEG, "deg",
            "asin(d1/h)+asin(d2/h), d = sqrt(R^2-rho^2)",
            "interior dihedral angle across B struts (AAB | BBB)")
    res.add("F8.bend_A_deg", 180 - dA / DEG, "deg", "180 - dihedral",
            "angle between panel normals across A struts")
    res.add("F8.bend_B_deg", 180 - dB / DEG, "deg", "180 - dihedral",
            "angle between panel normals across B struts")
    res.add("F8.face_AAB_apex_deg", apex / DEG, "deg",
            "acos(1 - B^2/(2A^2))", "AAB panel angle at the 5-valent hub")
    res.add("F8.face_AAB_base_deg", base / DEG, "deg", "(180 - apex)/2",
            "AAB panel base angles")
    res.add("F8.area_AAB_R1", area_AAB, "-", "(B/4) sqrt(4A^2 - B^2)",
            "AAB panel area on the unit sphere")
    res.add("F8.area_BBB_R1", area_BBB, "-", "(sqrt3/4) B^2",
            "BBB panel area on the unit sphere")
    res.add("F8.area_total_R1", tot_cf, "-", "60 AAB + 20 BBB",
            "faceted area of the 2V sphere, R = 1 (x R^2 for any R)")
    for key, v, f, dsc in (("nodes", 42, "10 v^2 + 2", "nodes of the sphere"),
                           ("struts", 120, "30 v^2", "struts of the sphere"),
                           ("panels", 80, "20 v^2", "panels of the sphere"),
                           ("struts_A", 60, "count", "A struts"),
                           ("struts_B", 60, "count", "B struts"),
                           ("hubs_valence5", 12, "count",
                            "5-valent hubs (type 5A)"),
                           ("hubs_valence6", 30, "count",
                            "6-valent hubs (type 2A+4B)"),
                           ("panels_AAB", 60, "count", "AAB panels"),
                           ("panels_BBB", 20, "count", "BBB panels")):
        res.add(f"F8.{key}", v, "-", f, f"2V Class I full sphere: {dsc}")


# =============================================================================
# assemble
# =============================================================================
_CACHE = {}


def _run():
    if "res" in _CACHE:
        return _CACHE["doc"], _CACHE["res"]
    CHECKS.clear()
    doc, res = Doc(), Results()
    MA = _model_a_basics()
    sec = MA["section"]
    doc.h1("GeoDomeLab extended hand calculations (validation cases F1 to F8)")
    doc.p("This document is produced by validation/hand_calc_ext.py. The "
          "script does not import the GeoDomeLab package. Each case defines "
          "a structure, solves it with a hand method (equilibrium, "
          "compatibility, the force method or the stiffness method as "
          "taught in the course), prints every intermediate number and ends "
          "with an independent check. The reference values feed "
          "validation/run_validation_ext.py, which runs the same structures "
          "through GeoDomeLab and compares.")
    doc.p("Conventions. z is vertical (up), tension is positive. A member "
          "from node i to node j has the unit vector c = (x_j - x_i)/L, the "
          "elongation e = c . (u_j - u_i) and the force N = EA (e/L - eps0) "
          "with the free strain eps0 = alpha dT. The member pushes or pulls "
          "its i end with +N c and its j end with -N c. A reaction is minus "
          "the sum of member forces and loads at its node. Values are kept "
          "in SI units and printed in kN, mm and MPa.")
    doc.table(["parameter", "value", "source"],
              [["E (A500 Gr. C)", f"{g6(E / 1e9)} GPa", "AISC 360 / ASTM A500"],
               ["Fy", f"{g6(FY / 1e6)} MPa", "ASTM A500 Gr. C round"],
               ["alpha", f"{g6(ALPHA)} 1/degC",
                "AISC Manual 6.5e-6 /degF for carbon steel"],
               ["rho", f"{g6(RHO)} kg/m3", "steel"],
               ["dT", f"+{g6(DT)} degC", "uniform rise above erection"],
               [SEC_MAIN[0], f"A = {g6(sec['A'] * 1e6)} mm2, "
                f"I = {g6(sec['I'] * 1e12)} mm4", "pi/4 (D^2 - d^2)"],
               [SEC_SMALL[0], f"A = {g6(chs(SEC_SMALL[1], SEC_SMALL[2])['A'] * 1e6)}"
                " mm2", "pi/4 (D^2 - d^2)"],
               ["Model A", f"R = {g6(hc.R)} m, L = {g6(MA['L'])} m, "
                f"k = EA/L = {kN(MA['k'])} kN/m", "hand_calc.py"]])
    case_f1(doc, res)
    case_f2(doc, res, MA)
    case_f3(doc, res, MA)
    case_f4(doc, res)
    case_f5(doc, res)
    case_f6(doc, res)
    case_f7(doc, res, MA)
    case_f8(doc, res)

    doc.h2("Summary of independent checks")
    npass = sum(1 for _, _, r, tol in CHECKS if r <= tol)
    doc.p(f"{npass} of {len(CHECKS)} internal checks are within tolerance. "
          f"A check compares two independent routes to the same number "
          f"(equilibrium against compatibility, the force method against "
          f"the displacement method, energy against work, closed form "
          f"against vectors, or a published value).")
    doc.table(["case", "check", "residual", "tolerance", "status"],
              [[cs, w, g6(r, 3), g6(tol, 2), "ok" if r <= tol else "FAIL"]
               for cs, w, r, tol in CHECKS])
    doc.h2("Reference values (results() dictionary)")
    doc.p("All values in SI units. These are the numbers the program must "
          "reproduce.")
    doc.table(["key", "value", "unit", "formula", "description"],
              [[kk, g6(v["value"], 10), v["unit"], v["formula"],
                v["description"]] for kk, v in res.items()])
    doc.h2("References")
    doc.items([
        "Kassimali, A. Matrix Analysis of Structures, 2nd ed. Cengage, 2012 "
        "(direct stiffness method for plane and space trusses, member "
        "forces from end displacements, support settlements and "
        "temperature changes by equivalent joint loads).",
        "Hibbeler, R. C. Structural Analysis, 10th ed. Pearson, 2018 "
        "(force method of analysis for trusses, principle of virtual work, "
        "unit-load method, temperature and fabrication errors, support "
        "settlement).",
        "Kenner, H. Geodesic Math and How to Use It. University of "
        "California Press, 1976 (Class I breakdown, 2V chord factors "
        "0.54653 and 0.61803).",
        "Crisfield, M. A. Non-linear Finite Element Analysis of Solids and "
        "Structures, Vol. 1. Wiley, 1991 (shallow two-bar truss, limit "
        "point, snap-through).",
        "Bazant, Z. P., and Cedolin, L. Stability of Structures. Oxford "
        "University Press, 1991 (von Mises truss, linearised versus exact "
        "limit load).",
        "Chopra, A. K. Dynamics of Structures, 5th ed. Pearson, 2017 "
        "(free vibration, lumped mass, mass-normalised modes, effective "
        "modal mass).",
        "AISC 360-22 Specification for Structural Steel Buildings, Chapter E "
        "(compression strength used for the utilisation notes), and the "
        "AISC Steel Construction Manual, 16th ed. (coefficient of thermal "
        "expansion of structural steel, 6.5e-6 /degF = 11.7e-6 /degC).",
    ])
    _CACHE["doc"], _CACHE["res"] = doc, res
    return doc, res


def results():
    """Reference values of every F case.

    Returns {key: {"value": float (SI), "unit": str, "formula": str,
    "description": str}}, keys like "F4.N2_temperature".
    """
    _, res = _run()
    return {k: dict(v) for k, v in res.items()}


def checks():
    """List of (case, description, residual, tolerance) internal checks."""
    _run()
    return list(CHECKS)


def definitions():
    """The structures of the F cases as plain data (SI, 0-based indices).

    Each entry has "nodes" (list of [x, y, z]), "members" (list of [i, j]),
    "sections" (name per member), "A" (area per member), "E", "alpha",
    "supports" ({node: [cx, cy, cz]} with True = fixed, False = free,
    float = spring N/m), and the load data of each case.
    """
    MA = hc.model_a()
    na = [list(map(float, MA["nodes"][n])) for n in sorted(MA["nodes"])]
    ma = [[i - 1, j - 1] for _, (i, j) in sorted(MA["members"].items())]
    sec = chs(SEC_MAIN[1], SEC_MAIN[2])
    secb = chs(SEC_SMALL[1], SEC_SMALL[2])
    pin = [True, True, True]
    modelA = {"nodes": na, "members": ma, "sections": [SEC_MAIN[0]] * 10,
              "A": [sec["A"]] * 10, "E": E, "alpha": ALPHA, "rho": RHO,
              "supports": {n: pin for n in range(1, 6)},
              "hand_numbering": "program index = hand node - 1, "
                                "hand member m -> index m - 1",
              "radius": hc.R}
    f4n = ["A", "B", "C", "D"]
    out = {
        "F1a": {"nodes": [list(F1_X1), list(F1_X2)], "members": [[0, 1]],
                "sections": [SEC_MAIN[0]], "A": [sec["A"]], "E": E,
                "alpha": ALPHA, "dT": DT, "supports": {0: pin, 1: pin}},
        "F1b": {"nodes": [[0.0, 0.0, 0.0], [F1B_L1, 0.0, 0.0],
                          [F1B_L1 + F1B_L2, 0.0, 0.0]],
                "members": [[0, 1], [1, 2]],
                "sections": [SEC_MAIN[0], SEC_SMALL[0]],
                "A": [sec["A"], secb["A"]], "E": E, "alpha": ALPHA, "dT": DT,
                "supports": {0: pin, 1: [False, True, True], 2: pin}},
        "F2": dict(modelA, dT=DT, heated_members_F2b=[0]),
        "F3": dict(modelA, settlement_F3a={"node": 1, "ux": 0.0, "uy": 0.0,
                                           "uz": -F3_SETTLE},
                   spread_F3b={"node": 2,
                               "ux": F3_SPREAD * math.cos(72 * DEG),
                               "uy": F3_SPREAD * math.sin(72 * DEG),
                               "uz": 0.0}),
        "F4": {"nodes": [list(F4_NODES[n]) for n in f4n],
               "node_names": f4n,
               "members": [[f4n.index(F4_BARS[i][0]),
                            f4n.index(F4_BARS[i][1])] for i in (1, 2, 3)],
               "sections": [F4_BARS[i][2][0] for i in (1, 2, 3)],
               "A": [chs(F4_BARS[i][2][1], F4_BARS[i][2][2])["A"]
                     for i in (1, 2, 3)],
               "E": E, "alpha": ALPHA,
               "supports": {0: pin, 1: pin, 2: pin, 3: [False, True, False]},
               "load": {"node": 3, "Fx": F4_P[0], "Fy": F4_P[1],
                        "Fz": F4_P[2]},
               "temperature": {"member": F4_HEATED_BAR - 1, "dT": DT},
               "settlement": {"node": f4n.index(F4_SETTLE_NODE),
                              "ux": F4_SETTLE[0], "uy": F4_SETTLE[1],
                              "uz": F4_SETTLE[2]}},
        "F5a": {"nodes": [[0.0, 0.0, 0.0], [F5_L, 0.0, 0.0]],
                "members": [[0, 1]], "sections": [SEC_MAIN[0]],
                "A": [sec["A"]], "E": E, "alpha": ALPHA,
                "supports": {0: [F5_K, True, True], 1: [False, True, True]},
                "load": {"node": 1, "Fx": F5_P, "Fy": 0.0, "Fz": 0.0}},
        "F5b": {"nodes": [[0.0, 0.0, 0.0], [F5_L, 0.0, 0.0]],
                "members": [[0, 1]], "sections": [SEC_MAIN[0]],
                "A": [sec["A"]], "E": E, "alpha": ALPHA,
                "supports": {0: pin, 1: [F5_K, True, True]},
                "load": {"node": 1, "Fx": F5_P, "Fy": 0.0, "Fz": 0.0}},
        "F5c": {"nodes": [[0.0, 0.0, 0.0], [F5_L, 0.0, 0.0]],
                "members": [[0, 1]], "sections": [SEC_MAIN[0]],
                "A": [sec["A"]], "E": E, "alpha": ALPHA, "dT": DT,
                "supports": {0: pin, 1: [F5_K, True, True]}},
        "F6": {"nodes": [[-VM_A, 0.0, 0.0], [VM_A, 0.0, 0.0],
                         [0.0, 0.0, VM_H]],
               "members": [[0, 2], [1, 2]],
               "sections": [SEC_MAIN[0]] * 2, "A": [sec["A"]] * 2, "E": E,
               "supports": {0: pin, 1: pin, 2: [False, True, False]},
               "load": {"node": 2, "Fx": 0.0, "Fy": 0.0, "Fz": -VM_PREF},
               "P_of_w": "vm_P(w), w = downward apex deflection"},
        "F7": dict(modelA, dead_pressure=hc.Q_DEAD, g=G),
        "F8": {"frequency": 2, "radius": F8_R, "class": "I (alternate)",
               "orientation": "icosahedron vertex at the pole"},
    }
    return out


def write_markdown(path=MD_PATH):
    doc, _ = _run()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc.markdown())
    return path


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:                                  # pragma: no cover
        pass
    doc, res = _run()
    print(doc.text())
    path = write_markdown()
    npass = sum(1 for _, _, r, tol in CHECKS if r <= tol)
    print()
    print("=" * 78)
    print(f"{len(res)} reference values, {npass} of {len(CHECKS)} internal "
          f"checks within tolerance.")
    print(f"geodome imported: {'geodome' in sys.modules}")
    print(f"markdown written to {path}")
    return 0 if npass == len(CHECKS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
