"""
Hand calculation of the two GeoDomeLab validation models.

This file is the "by hand" side of the validation.  It does NOT import the
geodome package.  Every number is built from first principles, in the same
order a person would do it on paper:

    geometry -> direction cosines -> element stiffness -> assembly of the
    free-DOF block -> solve {P} = [K]{D} -> member forces -> reactions ->
    equilibrium check -> AISC check -> buckling load factor

Only small arithmetic is used (sums, products, a 3 x 3 elimination).  numpy is
used as a calculator, never as a solver for the full structure.

Notation follows the course:  {P} = [K]{D}, free DOF numbered first,
z is vertical (up), x and y are horizontal, tension is positive.

MODEL A  1V pentagonal cap (the top cap of the icosahedron)
         6 nodes, 10 struts, 5 pinned base nodes, 3 free DOF (the crown).
         Solved by the direct stiffness method, fully by hand.

MODEL B  1V dome supported on its lower ring
         11 nodes, 25 struts, 5 pinned base nodes, 18 free DOF.
         Solved by the principle of virtual work, using the five-fold
         symmetry of the dome to reduce 18 unknowns to 3.

Run:  python validation/hand_calc.py
"""

from __future__ import annotations

import math

import numpy as np

# =============================================================================
# PARAMETERS  (edit here only)
# =============================================================================
R = 5.0                  # sphere radius                         [m]

E = 200e9                # Young's modulus, ASTM A500 Gr. C       [Pa]
FY = 345e6               # yield stress                          [Pa]
RHO = 7850.0             # steel density                         [kg/m3]
G = 9.80665              # gravity                               [m/s2]

D_OUT = 0.0889           # CHS 88.9 x 4.0 outside diameter        [m]
T_WALL = 0.0040          # wall thickness                        [m]

K_EFF = 1.0              # effective length factor
PHI_C = 0.90             # LRFD resistance factor, compression
PHI_T = 0.90             # LRFD resistance factor, tension

# Model A, load case A1: a single crown load with a horizontal part, so the
# response is NOT symmetric and every direction cosine matters.
PX_A1 = 20e3             # [N]  +x
PY_A1 = 0.0              # [N]
PZ_A1 = -60e3            # [N]  downward

# Load case 2 for both models: ASCE 7 combination 1.2D + 1.6S built from
# area loads (these are GeoDomeLab's default load values).
Q_DEAD = 300.0           # cladding on the faceted surface        [Pa]
P_SNOW = 1200.0          # snow on the horizontal projection      [Pa]
F_DEAD = 1.2
F_SNOW = 1.6

# Model B, load case B1: one crown load.
PZ_B1 = -50e3            # [N]  downward

DEG = math.pi / 180.0


# =============================================================================
# SECTION PROPERTIES  (CHS 88.9 x 4.0)
# =============================================================================
def section():
    d_in = D_OUT - 2.0 * T_WALL
    A = math.pi / 4.0 * (D_OUT ** 2 - d_in ** 2)
    I = math.pi / 64.0 * (D_OUT ** 4 - d_in ** 4)
    r = math.sqrt(I / A)
    return {"d_in": d_in, "A": A, "I": I, "r": r, "EA": E * A}


# =============================================================================
# GEOMETRY OF THE 1V ICOSAHEDRON (pole orientation)
# =============================================================================
# Crown on the z axis.  Upper ring of 5 nodes at height R/sqrt5 and radius
# a = 2R/sqrt5, at 0, 72, 144, 216, 288 deg.  Lower ring at -R/sqrt5,
# rotated by 36 deg.  Every edge of the icosahedron has the same length
#     L = R * sqrt(2 - 2/sqrt5) = 1.051462 R      (the 1V chord factor)
SQ5 = math.sqrt(5.0)
A_RING = 2.0 * R / SQ5                 # plan radius of both rings
H_RING = R / SQ5                       # height of a ring above/below centre
L_EDGE = R * math.sqrt(2.0 - 2.0 / SQ5)


def model_a_geometry():
    """Nodes 1..6 (hand numbering), base ring at z = 0 like the program."""
    H = R - H_RING                                # crown height above base
    nodes = {1: np.array([0.0, 0.0, H])}
    for i in range(5):
        th = 72.0 * i * DEG
        nodes[2 + i] = np.array([A_RING * math.cos(th),
                                 A_RING * math.sin(th), 0.0])
    # Members 1-5: spokes, from ring node (2..6) to the crown (1).
    # Members 6-10: ring, node 2->3, 3->4, 4->5, 5->6, 6->2.
    members = {}
    for i in range(5):
        members[1 + i] = (2 + i, 1)
    for i in range(5):
        members[6 + i] = (2 + i, 2 + (i + 1) % 5)
    return nodes, members


def model_b_geometry():
    """Nodes 1..11 (hand numbering), lower ring at z = 0 like the program."""
    zc = R + H_RING                               # crown height
    zu = 2.0 * H_RING                             # upper ring height
    nodes = {1: np.array([0.0, 0.0, zc])}
    for i in range(5):                            # upper ring 2..6
        th = 72.0 * i * DEG
        nodes[2 + i] = np.array([A_RING * math.cos(th),
                                 A_RING * math.sin(th), zu])
    for i in range(5):                            # lower ring 7..11
        th = (36.0 + 72.0 * i) * DEG
        nodes[7 + i] = np.array([A_RING * math.cos(th),
                                 A_RING * math.sin(th), 0.0])
    members = {}
    m = 1
    for i in range(5):                            # 1-5 spokes U -> C
        members[m] = (2 + i, 1); m += 1
    for i in range(5):                            # 6-10 upper ring
        members[m] = (2 + i, 2 + (i + 1) % 5); m += 1
    for i in range(5):                            # 11-20 zigzag L -> U
        members[m] = (7 + i, 2 + i); m += 1                 # L_i to U_i
        members[m] = (7 + i, 2 + (i + 1) % 5); m += 1       # L_i to U_i+1
    for i in range(5):                            # 21-25 lower ring
        members[m] = (7 + i, 7 + (i + 1) % 5); m += 1
    return nodes, members


def length_and_cosines(nodes, members):
    out = {}
    for m, (i, j) in members.items():
        v = nodes[j] - nodes[i]
        L = math.sqrt(float(v @ v))
        out[m] = (L, v / L)
    return out


def element_matrix(EA, L, c):
    """6x6 global element stiffness  (EA/L) [cc^T -cc^T ; -cc^T cc^T]."""
    cc = np.outer(c, c)
    return (EA / L) * np.block([[cc, -cc], [-cc, cc]])


# =============================================================================
# MODEL A  -  direct stiffness method by hand
# =============================================================================
def model_a():
    sec = section()
    EA = sec["EA"]
    nodes, members = model_a_geometry()
    geo = length_and_cosines(nodes, members)
    L = geo[1][0]
    k = EA / L                                    # axial stiffness of one strut

    # ---- direction cosines of the spokes (ring node -> crown) --------------
    #   c = ( -a cos(th), -a sin(th), H ) / L
    s = A_RING / L                                # horizontal part
    t = (R - H_RING) / L                          # vertical part

    # ---- free-DOF block [Kff] (crown ux, uy, uz = D1, D2, D3) --------------
    # The crown is the "j" end of every spoke, so each spoke adds +k c c^T.
    # Summed over 5 spokes at 72 deg spacing:
    #   sum cos^2 = sum sin^2 = 5/2,  sum cos = sum sin = sum cos*sin = 0
    # so [Kff] = k * diag( 5/2 s^2, 5/2 s^2, 5 t^2 ).  The ring struts join
    # two supported nodes and add nothing to [Kff].
    kxx = 2.5 * s * s
    kzz = 5.0 * t * t
    Kff = k * np.diag([kxx, kxx, kzz])

    # Cross-check of the closed form: add the five 3x3 blocks one by one.
    Kff_sum = np.zeros((3, 3))
    for m in range(1, 6):
        c = geo[m][1]
        Kff_sum += k * np.outer(c, c)

    # ---- load case A1: asymmetric crown load ------------------------------
    P1 = np.array([PX_A1, PY_A1, PZ_A1])
    D1 = P1 / np.diag(Kff)                        # [Kff] is diagonal
    a1 = _member_results_a(nodes, members, geo, k, D1, Pfree=P1)

    # ---- load case A2: 1.2D + 1.6S ----------------------------------------
    w_member = RHO * sec["A"] * L * G             # weight of one strut [N]
    face_area = math.sqrt(3.0) / 4.0 * L * L      # equilateral face, true area
    face_proj = 0.5 * A_RING ** 2 * math.sin(72.0 * DEG)   # its plan area
    crown_surf = 5.0 / 3.0 * face_area            # crown tributary, surface
    crown_proj = 5.0 / 3.0 * face_proj            # crown tributary, plan
    crown_self = 5.0 * 0.5 * w_member             # half of 5 spokes
    Dz = crown_self + Q_DEAD * crown_surf         # crown dead load   [N]
    Sz = P_SNOW * crown_proj                      # crown snow load   [N]
    Pz2 = -(F_DEAD * Dz + F_SNOW * Sz)
    P2 = np.array([0.0, 0.0, Pz2])
    D2 = P2 / np.diag(Kff)
    # Total load on the whole cap (for the reaction check)
    tot_dead = 10.0 * w_member + Q_DEAD * 5.0 * face_area
    tot_snow = P_SNOW * 5.0 * face_proj
    tot_P2 = -(F_DEAD * tot_dead + F_SNOW * tot_snow)
    # Load that lands directly on each supported ring node (2 faces, and
    # half of 1 spoke + 2 ring struts).  It goes straight into the reaction.
    ring_surf = 2.0 / 3.0 * face_area
    ring_proj = 2.0 / 3.0 * face_proj
    ring_self = 1.5 * w_member
    P_ring2 = -(F_DEAD * (ring_self + Q_DEAD * ring_surf)
                + F_SNOW * P_SNOW * ring_proj)
    a2 = _member_results_a(nodes, members, geo, k, D2, Pfree=P2,
                           ring_load_z=P_ring2)

    # ---- AISC 360 member check (spoke, load case A2) ----------------------
    KLr = K_EFF * L / sec["r"]
    limit = 4.71 * math.sqrt(E / FY)
    Fe = math.pi ** 2 * E / KLr ** 2
    Fcr = (0.658 ** (FY / Fe)) * FY if KLr <= limit else 0.877 * Fe
    phiPn_c = PHI_C * Fcr * sec["A"]
    phiPn_t = PHI_T * FY * sec["A"]
    Dt = D_OUT / T_WALL
    Dt_limit = 0.11 * E / FY
    N_spoke2 = a2["N"][1]
    util2 = abs(N_spoke2) / phiPn_c

    # ---- linear buckling of the crown (load case A2) ----------------------
    # Geometric stiffness of one spoke at the crown: (N/L)(I - c c^T).
    # Summed over 5 spokes with equal N:
    #   [Kg] = (N/L) diag( 5 - 5/2 s^2, 5 - 5/2 s^2, 5 - 5 t^2 )
    # [K + lambda Kg]{phi} = 0 with both matrices diagonal gives, per direction,
    #   lambda_i = - EA * k_i / ( N * g_i )        (L cancels)
    gxx = 5.0 - 2.5 * s * s
    gzz = 5.0 - 5.0 * t * t
    lam_z = -EA * kzz / (N_spoke2 * gzz)
    lam_xy = -EA * kxx / (N_spoke2 * gxx)

    return {
        "section": sec, "nodes": nodes, "members": members, "geo": geo,
        "L": L, "k": k, "s": s, "t": t, "kxx": kxx, "kzz": kzz,
        "Kff": Kff, "Kff_sum": Kff_sum,
        "K1": element_matrix(EA, L, geo[1][1]),
        "A1": {"P": P1, "D": D1, **a1},
        "A2": {"P": P2, "D": D2, "w_member": w_member,
               "face_area": face_area, "face_proj": face_proj,
               "crown_surf": crown_surf, "crown_proj": crown_proj,
               "crown_self": crown_self, "Dz": Dz, "Sz": Sz,
               "ring_surf": ring_surf, "ring_proj": ring_proj,
               "ring_self": ring_self, "P_ring": P_ring2,
               "tot_dead": tot_dead, "tot_snow": tot_snow, "tot_P": tot_P2,
               **a2},
        "aisc": {"KLr": KLr, "limit": limit, "Fe": Fe, "Fcr": Fcr,
                 "phiPn_c": phiPn_c, "phiPn_t": phiPn_t, "Dt": Dt,
                 "Dt_limit": Dt_limit, "N": N_spoke2, "util": util2,
                 "branch": "elastic (0.877 Fe)" if KLr > limit
                 else "inelastic (0.658^(Fy/Fe) Fy)"},
        "buckling": {"gxx": gxx, "gzz": gzz, "lam_z": lam_z,
                     "lam_xy": lam_xy, "lam_cr": min(lam_z, lam_xy)},
    }


def _member_results_a(nodes, members, geo, k, Dcrown, Pfree, ring_load_z=0.0):
    """Member forces and reactions of Model A for a crown displacement."""
    N = {}
    for m, (i, j) in members.items():
        c = geo[m][1]
        d_i = Dcrown if i == 1 else np.zeros(3)   # ring nodes do not move
        d_j = Dcrown if j == 1 else np.zeros(3)
        N[m] = k * float(c @ (d_j - d_i))
    # Reaction at a ring node = minus the force the spoke puts on it.
    # A spoke in tension N pulls its ring end towards the crown (+c), so the
    # support must push back with -N c.  Ring struts carry no force.
    # A load applied straight onto a support is carried by that support:
    # {R} = [Ksf]{Df} - {Ps}, so the ring node load enters with a minus sign.
    Rn = {}
    for m in range(1, 6):
        ring_node = members[m][0]
        Rn[ring_node] = -N[m] * geo[m][1] - np.array([0.0, 0.0, ring_load_z])
    Rsum = sum(Rn.values())
    Ptot = Pfree + np.array([0.0, 0.0, 5.0 * ring_load_z])
    return {"N": N, "R": Rn, "Rsum": Rsum, "Ptot": Ptot,
            "balance": Rsum + Ptot}


# =============================================================================
# MODEL B  -  principle of virtual work with five-fold symmetry
# =============================================================================
def model_b():
    sec = section()
    EA = sec["EA"]
    nodes, members = model_b_geometry()
    geo = length_and_cosines(nodes, members)
    L = geo[1][0]
    k = EA / L

    # Symmetric displacement field (vertical, symmetric loads):
    #   crown       u = (0, 0, w0)
    #   upper node  u = u_r * e_r(theta) + w1 * e_z       (no tangential part)
    #   lower ring  fixed
    # Generalised coordinates  q = { w0, u_r, w1 }.
    #
    # Elongation of each strut type as a row vector b (delta = b . q):
    s = A_RING / L
    t = (R - H_RING) / L                     # spoke rise / L  (same as Model A)
    b_spoke = np.array([t, s, -t])           # t*w0 + s*u_r - t*w1
    b_ring = np.array([0.0, 2.0 * math.sin(36.0 * DEG), 0.0])
    zz_r = A_RING * (1.0 - math.cos(36.0 * DEG)) / L
    zz_z = 2.0 * H_RING / L
    b_zig = np.array([0.0, zz_r, zz_z])      # zigzag, lower node fixed

    # Virtual work:  delta_q^T (sum n k b b^T) q = delta_q^T f
    #  -> [Ks] {q} = {f},  n = number of struts of each type (5, 5, 10).
    Ks = k * (5.0 * np.outer(b_spoke, b_spoke)
              + 5.0 * np.outer(b_ring, b_ring)
              + 10.0 * np.outer(b_zig, b_zig))
    Ks_over_k = Ks / k

    # Tributary areas and self weight for load case B2 --------------------
    w_member = RHO * sec["A"] * L * G
    face_area = math.sqrt(3.0) / 4.0 * L * L
    proj_top = 0.5 * A_RING ** 2 * math.sin(72.0 * DEG)
    proj_band = 0.5 * A_RING ** 2 * (2.0 * math.sin(36.0 * DEG)
                                     - math.sin(72.0 * DEG))
    # crown: 5 top faces.  upper node: 2 top + 3 band faces.
    crown_surf = 5.0 / 3.0 * face_area
    crown_proj = 5.0 / 3.0 * proj_top
    up_surf = 5.0 / 3.0 * face_area
    up_proj = (2.0 * proj_top + 3.0 * proj_band) / 3.0
    crown_self = 2.5 * w_member               # 5 spokes, half each
    up_self = 2.5 * w_member                  # 1 spoke + 2 ring + 2 zigzag
    # lower (supported) node: 3 band faces, half of 2 zigzag + 2 ring struts
    low_surf = 3.0 / 3.0 * face_area
    low_proj = 3.0 * proj_band / 3.0
    low_self = 2.0 * w_member
    Pl2 = -(F_DEAD * (low_self + Q_DEAD * low_surf)
            + F_SNOW * P_SNOW * low_proj)
    Pc2 = -(F_DEAD * (crown_self + Q_DEAD * crown_surf)
            + F_SNOW * P_SNOW * crown_proj)
    Pu2 = -(F_DEAD * (up_self + Q_DEAD * up_surf)
            + F_SNOW * P_SNOW * up_proj)
    # Whole dome, for the reaction check (15 faces, 25 struts)
    tot_dead = 25.0 * w_member + Q_DEAD * 15.0 * face_area
    tot_snow = P_SNOW * 5.0 * (proj_top + 2.0 * proj_band)
    tot_P2 = -(F_DEAD * tot_dead + F_SNOW * tot_snow)

    cases = {}
    for name, Pc, Pu, Pl, totP in (("B1", PZ_B1, 0.0, 0.0, PZ_B1),
                                   ("B2", Pc2, Pu2, Pl2, tot_P2)):
        f = np.array([Pc, 0.0, 5.0 * Pu])
        q, steps = gauss_3x3(Ks, f)
        w0, ur, w1 = q
        N_spoke = k * float(b_spoke @ q)
        N_ring = k * float(b_ring @ q)
        N_zig = k * float(b_zig @ q)
        # Each lower node holds two zigzag struts.  Unit vector from the
        # lower node to the upper node has vertical part zz_z (= 2h/L) and a
        # horizontal part along the line to that upper node.
        Rz_strut = -N_zig * zz_z                 # vertical reaction per strut
        Rz_node = 2.0 * Rz_strut - Pl            # + load placed on the support
        # Horizontal: the two struts' horizontal parts add along the radius
        # of the lower node (their tangential parts cancel).  The plan vector
        # of a zigzag, lower node -> upper node, has a radial part
        # a(cos36 - 1) along the lower node's own radius.
        rad_part = A_RING * (math.cos(36.0 * DEG) - 1.0) / L
        Rr_node = -2.0 * N_zig * rad_part        # reaction, + = outward
        # Crown equilibrium (vertical): each spoke in tension N pulls the
        # crown down along -t, so  -5 N_spoke t + Pc = 0.
        crown_eq = -5.0 * N_spoke * t + Pc
        # Upper node equilibrium.  Forces on the node from tension N:
        #   spoke  -> towards the crown:  radial -N s,      vertical +N t
        #   2 ring -> towards neighbours: radial -2 N sin36
        #   2 zig  -> towards lower ring: radial -2 N zz_r, vertical -2 N zz_z
        up_rad = (-N_spoke * s - N_ring * 2.0 * math.sin(36.0 * DEG)
                  - 2.0 * N_zig * zz_r)
        up_vert = (N_spoke * t - 2.0 * N_zig * zz_z) + Pu
        cases[name] = {
            "f": f, "q": q, "steps": steps, "w0": w0, "ur": ur, "w1": w1,
            "N_spoke": N_spoke, "N_ring": N_ring, "N_zig": N_zig,
            "Rz_node": Rz_node, "Rr_node": Rr_node, "Pc": Pc, "Pu": Pu,
            "Pl": Pl, "N_lowring": 0.0,
            "tot_P": totP, "crown_eq": crown_eq,
            "up_rad": up_rad, "up_vert": up_vert,
            "Rz_total": 5.0 * Rz_node,
        }
    return {
        "section": sec, "nodes": nodes, "members": members, "geo": geo,
        "L": L, "k": k, "s": s, "t": t,
        "b_spoke": b_spoke, "b_ring": b_ring, "b_zig": b_zig,
        "Ks": Ks, "Ks_over_k": Ks_over_k,
        "loads": {"w_member": w_member, "face_area": face_area,
                  "proj_top": proj_top, "proj_band": proj_band,
                  "crown_surf": crown_surf, "crown_proj": crown_proj,
                  "up_surf": up_surf, "up_proj": up_proj,
                  "low_surf": low_surf, "low_proj": low_proj,
                  "tot_dead": tot_dead, "tot_snow": tot_snow},
        **cases,
    }


def gauss_3x3(K, f):
    """Gaussian elimination written out as on paper; returns q and the steps."""
    A = np.array(K, dtype=float).copy()
    b = np.array(f, dtype=float).copy()
    steps = []
    for p in range(3):
        for r in range(p + 1, 3):
            m = A[r, p] / A[p, p]
            A[r, :] -= m * A[p, :]
            b[r] -= m * b[p]
            steps.append((r + 1, p + 1, m))
    q = np.zeros(3)
    for r in (2, 1, 0):
        q[r] = (b[r] - A[r, r + 1:] @ q[r + 1:]) / A[r, r]
    return q, {"multipliers": steps, "U": A, "b": b}


# =============================================================================
# PRINT
# =============================================================================
def _v(x, unit="", nd=6):
    return f"{x:.{nd}g} {unit}".strip()


def report():
    A = model_a()
    B = model_b()
    sec = A["section"]
    kN = 1e-3
    out = []
    p = out.append
    p("=" * 78)
    p("HAND CALCULATION - GeoDomeLab validation models")
    p("=" * 78)
    p(f"Section CHS 88.9x4.0: d = {sec['d_in']*1e3:.1f} mm, "
      f"A = {sec['A']*1e6:.2f} mm2, I = {sec['I']*1e12:.4g} mm4, "
      f"r = {sec['r']*1e3:.3f} mm, EA = {sec['EA']*kN:,.1f} kN")
    p(f"1V edge L = {L_EDGE:.6f} m = {L_EDGE/R:.6f} R (chord factor)")
    p("")
    p("MODEL A  1V cap: s = %.6f  t = %.6f  k = EA/L = %.3f kN/m"
      % (A["s"], A["t"], A["k"] * kN))
    p("  [Kff] (kN/m) =\n" + np.array2string(A["Kff"] * kN, precision=3))
    p("  check by summing 5 blocks:\n"
      + np.array2string(A["Kff_sum"] * kN, precision=3))
    for lc in ("A1", "A2"):
        r = A[lc]
        p(f"  {lc}: P = {np.round(r['P']*kN, 4)} kN")
        p(f"      D = {np.round(r['D']*1e3, 6)} mm")
        p("      N = " + ", ".join(f"{m}:{r['N'][m]*kN:+.4f}" for m in r["N"]))
        p(f"      sum R = {np.round(r['Rsum']*kN, 4)} kN")
    a2 = A["A2"]
    p(f"  A2 loads: w_member = {a2['w_member']:.3f} N, crown D = "
      f"{a2['Dz']*kN:.4f} kN, S = {a2['Sz']*kN:.4f} kN, Pz = "
      f"{a2['P'][2]*kN:.4f} kN, total = {a2['tot_P']*kN:.4f} kN")
    ai = A["aisc"]
    p(f"  AISC: KL/r = {ai['KLr']:.2f} (limit {ai['limit']:.2f}), "
      f"Fe = {ai['Fe']/1e6:.3f} MPa, Fcr = {ai['Fcr']/1e6:.3f} MPa, "
      f"phiPn_c = {ai['phiPn_c']*kN:.3f} kN, util = {ai['util']:.5f}")
    bk = A["buckling"]
    p(f"  Buckling: lambda_z = {bk['lam_z']:.4f}, lambda_xy = "
      f"{bk['lam_xy']:.4f}")
    p("")
    p("MODEL B  1V dome, PVW: [Ks]/k =\n"
      + np.array2string(B["Ks_over_k"], precision=6))
    for lc in ("B1", "B2"):
        r = B[lc]
        p(f"  {lc}: f = {np.round(r['f']*kN, 4)} kN")
        p(f"      w0 = {r['w0']*1e3:.6f} mm, u_r = {r['ur']*1e3:.6f} mm, "
          f"w1 = {r['w1']*1e3:.6f} mm")
        p(f"      N spoke = {r['N_spoke']*kN:+.4f}, ring = "
          f"{r['N_ring']*kN:+.4f}, zigzag = {r['N_zig']*kN:+.4f} kN")
        p(f"      Rz/node = {r['Rz_node']*kN:.4f} kN, R radial/node = "
          f"{r['Rr_node']*kN:.4f} kN, sum Rz = {r['Rz_total']*kN:.4f}, "
          f"total load = {r['tot_P']*kN:.4f}")
        p(f"      checks: crown {r['crown_eq']:.2e} N, upper node radial "
          f"{r['up_rad']:.2e} N, vertical {r['up_vert']:.2e} N")
    return "\n".join(out)


if __name__ == "__main__":
    print(report())
