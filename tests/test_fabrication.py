"""
Tests for the fabrication schedule (spec task C, ``geodome/fabrication.py``).

Every check is against an independent answer:

* published chord factors and strut counts: 2V (Kenner 1976; Domerama),
  3V to 6V domes (Domerama chord-factor chart, "method 1" icosahedral,
  https://domerama.com/calculators/chord-charts, read 2026-11-16), and the
  2V closed forms 2 sin(acos(1/sqrt5)/4) and (sqrt5 - 1)/2;
* complete Class I spheres 1V-6V: 10 nu^2 + 2 hubs, 30 nu^2 struts,
  20 nu^2 panels, twelve 5-valent hubs and all others 6-valent, Euler 2;
* axial angle from raw vectors (computed here, not by the module) and the
  closed form arccos(L/2R) at both ends of every strut; the 2V tangent
  (bend) angles 15.86 / 18 deg of the Domerama chart;
* dihedral angles: the F8 hand-calculation formula for 2V, the icosahedron
  acos(-sqrt5/3), the regular octahedron acos(-1/3), a dented (non-convex)
  octahedron checked with an inside/outside oracle, and a flat grid (180);
* identities: sum of panel areas = faceted surface area (sum of the nodal
  tributary areas), sum of valences = 2 x struts, table totals = counts,
  tube mass with zero hub offset = DomeModel.steel_mass();
* cut length L - 2 x hub offset, input validation, CSV export round trip,
  stock plan (closed-form pieces per bar, first-fit decreasing validity and
  the lower bound).

Run with pytest::

    pytest -q tests/test_fabrication.py

or standalone, without pytest installed::

    python tests/test_fabrication.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import fabrication as fb                               # noqa: E402
from geodome import geometry as geo                                 # noqa: E402
from geodome.model import DomeModel, default_model                  # noqa: E402

SQ5 = math.sqrt(5.0)
DEG = 180.0 / math.pi

# Published chord factors (as printed) and strut counts, sorted by chord
# factor.  2V: Kenner (1976) and the Domerama chart; 3V-6V: Domerama chart
# ("3V 3/8 = 75 faces", "4V 1/2 method 1, 160 faces", "5V 7/15 = 225
# faces", "6V 360 faces").
PUBLISHED = {
    2: [("0.54653", 30), ("0.61803", 35)],
    3: [("0.34862", 30), ("0.40355", 40), ("0.41241", 50)],
    4: [("0.25318", 30), ("0.29453", 60), ("0.29524", 30), ("0.29859", 30),
        ("0.31287", 70), ("0.32492", 30)],
    5: [("0.19814743", 30), ("0.22568578", 60), ("0.23159760", 30),
        ("0.23179025", 30), ("0.24508578", 50), ("0.24534642", 10),
        ("0.24724291", 60), ("0.25516701", 50), ("0.26159810", 30)],
    6: [("0.16256722", 30), ("0.18190825", 60), ("0.1873834", 30),
        ("0.19047686", 30), ("0.19801258", 60), ("0.20281969", 90),
        ("0.20590774", 130), ("0.21535373", 65), ("0.21662821", 60)],
}
PUBLISHED_FACES = {2: 40, 3: 75, 4: 160, 5: 225, 6: 360}


# ================================================================ helpers ===
def full_sphere(nu, R=1.0):
    """Complete Class I sphere as a custom truss (faces = panels)."""
    v, f = geo.icosahedron()
    pts, tris = (v, f) if nu == 1 else geo.subdivide_class1(v, f, nu)
    return geo.custom_geometry(pts * R, geo.edges_from_faces(tris),
                               faces=tris, title=f"{nu}V sphere")


def circumradius(a, b, c):
    s = 0.5 * (a + b + c)
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return a * b * c / (4.0 * area)


def f8_dihedrals(R=1.0):
    """2V dihedral angles from the F8 hand calculation (asin form)."""
    A = 2.0 * math.sin(math.acos(1.0 / SQ5) / 4.0)
    B = (SQ5 - 1.0) / 2.0
    r_aab = circumradius(A, A, B) * R
    r_bbb = circumradius(B, B, B) * R
    A, B = A * R, B * R
    d_aab = math.sqrt(R * R - r_aab ** 2)
    d_bbb = math.sqrt(R * R - r_bbb ** 2)
    hA = math.sqrt(R * R - A * A / 4.0)
    hB = math.sqrt(R * R - B * B / 4.0)
    dA = 2.0 * math.asin(d_aab / hA)
    dB = math.asin(d_aab / hB) + math.asin(d_bbb / hB)
    return dA * DEG, dB * DEG


def row_by(rows, **kw):
    out = [r for r in rows if all(r[k] == v for k, v in kw.items())]
    assert len(out) == 1, (kw, out)
    return out[0]


# ======================================================== chord factors ===
def test_2v_chord_factors_kenner_and_closed_form():
    R = 5.0
    t = fb.fabrication_tables(geo.generate_dome(2, R, 0.5))
    st = sorted(t["struts"], key=lambda r: r["chord_factor"])
    cfA = 2.0 * math.sin(math.acos(1.0 / SQ5) / 4.0)
    cfB = (SQ5 - 1.0) / 2.0
    assert abs(st[0]["chord_factor"] - 0.54653) < 5e-6        # Kenner 1976
    assert abs(st[1]["chord_factor"] - 0.61803) < 5e-6
    assert abs(st[0]["chord_factor"] - cfA) < 1e-13
    assert abs(st[1]["chord_factor"] - cfB) < 1e-13
    assert abs(st[0]["length_mm"] - cfA * R * 1e3) < 1e-9
    assert (st[0]["count"], st[1]["count"]) == (30, 35)       # Domerama 2V
    assert len(t["faces"]) == 40


def test_published_chord_factors_3v_to_6v():
    for nu, pub in PUBLISHED.items():
        g = geo.generate_dome(nu, 1.0, 0.5)
        t = fb.fabrication_tables(g)
        st = sorted(t["struts"], key=lambda r: r["chord_factor"])
        assert len(st) == len(pub), nu
        assert len(t["faces"]) == PUBLISHED_FACES[nu], nu
        for r, (txt, count) in zip(st, pub):
            dec = len(txt.split(".")[1])
            # 5-decimal charts: within half a unit of the last decimal;
            # 7-8-decimal charts: within one unit of the last decimal.
            tol = 0.5 * 10.0 ** -dec if dec <= 5 else 10.0 ** -dec
            assert abs(r["chord_factor"] - float(txt)) <= tol, \
                (nu, r["chord_factor"], txt)
            assert r["count"] == count, (nu, txt, r["count"], count)


PUBLISHED_7V_8V = {
    # Domerama chart: "7V (10/21), 455 faces" and "8V, 640 faces"
    7: (10.0 / 21.0, 455, [
        ("0.13774", 30), ("0.15197", 60), ("0.15664", 30), ("0.16154", 30),
        ("0.16480", 60), ("0.17066", 30), ("0.17098", 60), ("0.17132", 60),
        ("0.17353", 50), ("0.17585", 70), ("0.18155", 50), ("0.18161", 30),
        ("0.18237", 50), ("0.18548", 60), ("0.1879", 30)]),
    8: (0.5, 640, [
        ("0.11946", 30), ("0.13033", 60), ("0.13424", 30), ("0.14018", 30),
        ("0.14056", 60), ("0.14548", 60), ("0.14628", 30), ("0.14803", 60),
        ("0.14862", 60), ("0.15267", 60), ("0.15296", 70), ("0.15315", 30),
        ("0.15477", 60), ("0.15636", 90), ("0.16033", 60), ("0.16036", 30),
        ("0.16088", 70), ("0.16300", 60), ("0.16465", 30)]),
}


def test_published_chord_factors_7v_8v():
    """7V / 8V: the chart merges some lengths that differ by < 5e-5 R.

    With a tiny cut tolerance every distinct strut length is its own cut
    group.  Each group is assigned to the nearest published chord factor:
    it must lie within 6e-5 of it, the largest group of every published row
    within 1e-5 (one unit of the 5th decimal: the chart sometimes truncates
    its last digit, e.g. 7V "0.1879" for 0.1879085, like 6V "0.16256722"
    for 0.162567229), and the counts assigned to each published row must
    add up to the published count.
    """
    for nu, (fr, nfaces, pub) in PUBLISHED_7V_8V.items():
        g = geo.generate_dome(nu, 1.0, fr)
        t = fb.fabrication_tables(g, length_tol=1e-9)
        assert len(t["faces"]) == nfaces
        vals = np.array([float(p) for p, _ in pub])
        got = {k: [] for k in range(len(pub))}
        for r in t["struts"]:
            k = int(np.argmin(np.abs(vals - r["chord_factor"])))
            assert abs(vals[k] - r["chord_factor"]) < 6e-5, (nu, r)
            got[k].append((r["count"], r["chord_factor"]))
        for k, (p, count) in enumerate(pub):
            assert sum(c for c, _ in got[k]) == count, (nu, p, got[k])
            main = max(got[k])[1]
            assert abs(main - float(p)) <= 1e-5, (nu, p, main)


# ========================================================== full sphere ===
def test_full_sphere_counts_and_valences():
    for nu in range(1, 7):
        g = full_sphere(nu, R=3.0)
        t = fb.fabrication_tables(g)
        n, m, f = 10 * nu ** 2 + 2, 30 * nu ** 2, 20 * nu ** 2
        assert len(t["nodes"]) == n and len(t["members"]) == m
        assert len(t["faces"]) == f
        val = np.array([r["valence"] for r in t["nodes"]])
        assert int(np.sum(val == 5)) == 12, nu
        assert int(np.sum(val == 6)) == n - 12, nu
        assert int(val.sum()) == 2 * m                     # handshake lemma
        assert sum(r["count"] for r in t["hubs"]) == n
        assert sum(r["count"] for r in t["struts"]) == m
        assert sum(r["count"] for r in t["panels"]) == f
        assert sum(r["count"] for r in t["dihedrals"]) == m   # closed: all
        summ = {r["quantity"]: r["value"] for r in t["summary"]}
        assert summ["Euler characteristic V - E + F"] == 2
        assert summ["struts with two panels / one panel / none"] == \
            f"{m} / 0 / 0"
        assert summ["max hub distance from the sphere"] < 1e-12


def test_2v_sphere_hub_and_panel_types():
    t = fb.fabrication_tables(full_sphere(2))
    hubs = {(r["struts"], r["arrangement"]): r["count"] for r in t["hubs"]}
    # Every hub lies on the sphere; the default support set of a custom
    # truss is the lowest node, so add the counts over locations.
    tot = {}
    for r in t["hubs"]:
        tot[(r["struts"], r["arrangement"])] = \
            tot.get((r["struts"], r["arrangement"]), 0) + r["count"]
    assert tot == {("5A", "A-A-A-A-A"): 12, ("2A+4B", "A-B-B-A-B-B"): 30}, hubs
    p = {r["panel_type"]: r for r in t["panels"]}
    assert p["AAB"]["count"] == 60 and p["BBB"]["count"] == 20
    A = 2.0 * math.sin(math.acos(1.0 / SQ5) / 4.0)
    B = (SQ5 - 1.0) / 2.0
    apex = math.acos(1.0 - B * B / (2.0 * A * A)) * DEG      # F8
    assert abs(p["AAB"]["angle_3_deg"] - apex) < 1e-10        # opposite B
    assert abs(p["AAB"]["angle_1_deg"] - (180.0 - apex) / 2.0) < 1e-10
    area_aab = 0.25 * B * math.sqrt(4.0 * A * A - B * B)
    area_bbb = math.sqrt(3.0) / 4.0 * B * B
    assert abs(p["AAB"]["area_m2"] - area_aab) < 1e-13
    assert abs(p["BBB"]["area_m2"] - area_bbb) < 1e-13
    total = sum(r["total_area_m2"] for r in t["panels"])
    assert abs(total - (60 * area_aab + 20 * area_bbb)) < 1e-12
    assert total < 4.0 * math.pi                         # inscribed facets


# =========================================================== axial angle ===
def test_axial_angle_vectors_vs_formula():
    R = 7.5
    g = geo.generate_dome(4, R, 0.5)
    t = fb.fabrication_tables(g)
    # independent: the crown is the pole (0, 0, R) of the generating sphere
    c = np.array([0.0, 0.0, g.nodes[g.crown_node, 2] - R])
    assert np.allclose(np.linalg.norm(g.nodes - c, axis=1), R, atol=1e-12)
    worst = 0.0
    for r in t["members"]:
        i, j = r["node_i"], r["node_j"]
        xi, xj = g.nodes[i], g.nodes[j]
        L = np.linalg.norm(xj - xi)
        ref = math.degrees(math.acos(L / (2.0 * R)))
        for x0, x1, got in ((xi, xj, r["axial_i_deg"]),
                            (xj, xi, r["axial_j_deg"])):
            dvec, rad = x1 - x0, c - x0
            cosang = dvec @ rad / (np.linalg.norm(dvec) * np.linalg.norm(rad))
            vec = math.degrees(math.acos(cosang))
            worst = max(worst, abs(got - vec), abs(got - ref))
    assert worst < 1e-9, worst
    for r in t["struts"]:
        assert r["axial_check_deg"] < 1e-11
        assert r["axial_spread_deg"] < 1e-11
        assert abs(r["axial_angle_deg"] + r["tangent_angle_deg"] - 90.0) \
            < 1e-12
        cf = r["chord_factor"]
        assert abs(r["tangent_angle_deg"]
                   - math.degrees(math.asin(cf / 2.0))) < 1e-10


def test_2v_axial_and_tangent_angles_published():
    t = fb.fabrication_tables(geo.generate_dome(2, 3.0, 0.5))
    st = sorted(t["struts"], key=lambda r: r["chord_factor"])
    # Clinton relation alpha = arccos(cf/2).  Closed forms: the A strut
    # subtends half the icosahedron edge angle acos(1/sqrt5), so
    # alpha_A = 90 - acos(1/sqrt5)/4 = 74.1413 deg; the B strut subtends
    # 36 deg, so alpha_B = 72 deg exactly.  (A rounded table prints
    # 74.142, a rounding slip of 74.14126.)
    alpha_A = 90.0 - math.degrees(math.acos(1.0 / SQ5)) / 4.0
    assert abs(st[0]["axial_angle_deg"] - alpha_A) < 1e-11
    assert abs(st[0]["axial_angle_deg"] - 74.1413) < 5e-5
    assert abs(st[1]["axial_angle_deg"] - 72.0) < 1e-11
    # Domerama 2V bend angles 15.86 deg and 18 deg (= 90 - axial)
    assert abs(st[0]["tangent_angle_deg"] - 15.86) < 5e-3
    assert abs(st[1]["tangent_angle_deg"] - 18.0) < 1e-11


# ============================================================== dihedral ===
def test_2v_dihedral_f8_formula():
    dA, dB = f8_dihedrals()
    for g in (full_sphere(2, 1.0), geo.generate_dome(2, 4.0, 0.5)):
        t = fb.fabrication_tables(g)
        rA = row_by(t["dihedrals"], strut_family="A")
        rB = row_by(t["dihedrals"], strut_family="B")
        assert rA["panels"] == "AAB | AAB" and rB["panels"] == "AAB | BBB"
        assert abs(rA["dihedral_deg"] - dA) < 1e-10
        assert abs(rB["dihedral_deg"] - dB) < 1e-10
        assert abs(rA["bend_deg"] - (180.0 - dA)) < 1e-10
        assert rA["dihedral_spread_deg"] < 1e-10
        assert rA["dihedral_check_deg"] < 1e-10
        assert rB["dihedral_check_deg"] < 1e-10
    # counts: sphere 60 / 60 interior; hemisphere: the 10 base B struts
    # carry one panel only
    t = fb.fabrication_tables(geo.generate_dome(2, 4.0, 0.5))
    assert row_by(t["dihedrals"], strut_family="A")["count"] == 30
    assert row_by(t["dihedrals"], strut_family="B")["count"] == 25


def test_platonic_dihedrals():
    # icosahedron (1V sphere): acos(-sqrt5/3) = 138.1896851 deg
    t = fb.fabrication_tables(full_sphere(1, 2.0))
    (r,) = t["dihedrals"]
    assert abs(r["dihedral_deg"] - math.degrees(math.acos(-SQ5 / 3.0))) < 1e-10
    (s,) = t["struts"]
    assert abs(s["chord_factor"] - 1.05146) < 5e-6           # Domerama 1V
    assert abs(s["tangent_angle_deg"] - 31.72) < 5e-3
    # regular octahedron: dihedral acos(-1/3), L = sqrt2 R, axial 45 deg
    P = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1],
                  [0, 0, -1]], dtype=float) * 2.0
    F = [[4, 0, 2], [4, 2, 1], [4, 1, 3], [4, 3, 0],
         [5, 2, 0], [5, 1, 2], [5, 3, 1], [5, 0, 3]]
    g = geo.custom_geometry(P, geo.edges_from_faces(F), faces=F)
    t = fb.fabrication_tables(g)
    (r,) = t["dihedrals"]
    assert r["count"] == 12
    assert abs(r["dihedral_deg"] - math.degrees(math.acos(-1.0 / 3.0))) < 1e-10
    assert r["dihedral_check_deg"] < 1e-10
    (s,) = t["struts"]
    assert abs(s["axial_angle_deg"] - 45.0) < 1e-10
    assert abs(s["chord_factor"] - math.sqrt(2.0)) < 1e-12
    # faces are listed counter-clockwise seen from outside
    for fr in t["faces"]:
        a, b, c = (P[fr[f"node_{k}"]] for k in (1, 2, 3))
        assert np.cross(b - a, c - a) @ ((a + b + c) / 3.0) > 0.0


def _wedge_dihedral(P, e, q1, q2, inside):
    """Interior dihedral across edge e from the oracle ``inside(point)``."""
    a, b = P[e[0]], P[e[1]]
    t = (b - a) / np.linalg.norm(b - a)

    def w(q):
        v = P[q] - a
        v = v - (v @ t) * t
        return v / np.linalg.norm(v)

    w1, w2 = w(q1), w(q2)
    phi = math.acos(max(-1.0, min(1.0, float(w1 @ w2))))
    mid = 0.5 * (a + b)
    probe = mid + 1e-4 * (w1 + w2) / np.linalg.norm(w1 + w2)
    return phi if inside(probe) else 2.0 * math.pi - phi


def test_non_convex_dented_octahedron():
    """Re-entrant (valley) folds get dihedral > 180 deg."""
    P = np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
                  [0, 0, -0.5],                     # 4: dented top vertex
                  [0, 0, -1.0]])                    # 5: bottom vertex
    F = [[4, 0, 1], [4, 1, 2], [4, 2, 3], [4, 3, 0],
         [5, 1, 0], [5, 2, 1], [5, 3, 2], [5, 0, 3]]
    g = geo.custom_geometry(P, geo.edges_from_faces(F), faces=F,
                            title="dented octahedron")

    def inside(x):          # between the dent surface and the lower pyramid
        s = abs(x[0]) + abs(x[1])
        return s < 1.0 and (-1.0 + s) < x[2] < (-0.5 + 0.5 * s)

    t = fb.fabrication_tables(g)
    reflex = 0
    for r in t["members"]:
        i, j = r["node_i"], r["node_j"]
        opp = [[v for v in f if v not in (i, j)][0] for f in F
               if i in f and j in f]
        assert len(opp) == 2
        ref = _wedge_dihedral(P, (i, j), opp[0], opp[1], inside) * DEG
        assert abs(r["dihedral_deg"] - ref) < 1e-9, (i, j, r["dihedral_deg"],
                                                     ref)
        reflex += r["dihedral_deg"] > 180.0
    assert reflex == 4                       # the four dent ridges T-Ek
    summ = {r["quantity"]: r["value"] for r in t["summary"]}
    assert summ["Euler characteristic V - E + F"] == 2


def test_flat_grid_custom_truss():
    """Coplanar truss: no sphere (NaN), flat panels meet at 180 deg."""
    xs = np.linspace(0.0, 3.0, 4)
    P = np.array([[x, y, 2.0] for y in xs for x in xs])

    def idx(i, j):
        return 4 * j + i

    F = []
    for j in range(3):
        for i in range(3):
            F += [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)],
                  [idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]]
    g = geo.custom_geometry(P, geo.edges_from_faces(F), faces=F,
                            base_nodes=[0, 3, 12, 15], title="flat grid")
    t = fb.fabrication_tables(g)
    for r in t["struts"]:
        assert math.isnan(r["chord_factor"])
        assert math.isnan(r["axial_formula_deg"])
        assert abs(r["axial_angle_deg"] - 90.0) < 1e-12   # axis = +-z normal
    for r in t["dihedrals"]:
        assert abs(r["dihedral_deg"] - 180.0) < 1e-12
        assert math.isnan(r["dihedral_formula_deg"])
    total = sum(r["total_area_m2"] for r in t["panels"])
    assert abs(total - 9.0) < 1e-12
    assert abs(total - g.surface_area) < 1e-12
    loc = {r["node"]: r["location"] for r in t["nodes"]}
    assert [loc[k] for k in (0, 3, 12, 15)] == ["base"] * 4
    assert loc[1] == "edge" and loc[5] == "interior"


def test_custom_truss_without_panels():
    """A bare space truss (star dome of docs/VALIDATION.md, no faces)."""
    cm = 0.01
    P = [[0, 0, 8.216]]
    for k in range(6):
        a = math.radians(60.0 * k)
        P.append([25 * math.cos(a), 25 * math.sin(a), 6.216])
    for k in range(6):
        a = math.radians(30.0 + 60.0 * k)
        P.append([50 * math.cos(a), 50 * math.sin(a), 0.0])
    P = np.array(P) * cm
    M = [[0, k] for k in range(1, 7)]
    M += [[k, k % 6 + 1] for k in range(1, 7)]
    for k in range(6):
        M += [[1 + k, 7 + k], [1 + (k + 1) % 6, 7 + k]]
    g = geo.custom_geometry(P, M, title="star dome")
    t = fb.fabrication_tables(g, hub_offset=0.01)
    assert t["panels"] == [] and t["dihedrals"] == [] and t["faces"] == []
    assert len(t["members"]) == 24 and len(t["nodes"]) == 13
    assert sum(r["count"] for r in t["struts"]) == 24
    # families by length: A = ring 25.000 cm, B = crown spokes 25.080 cm,
    # C = outer struts 31.600 cm (docs/VALIDATION.md)
    lens = [round(r["length_mm"], 3) for r in t["struts"]]
    assert lens == [250.0, round(250.79872, 3), round(316.00240, 3)], lens
    assert [r["count"] for r in t["struts"]] == [6, 6, 12]
    crown = row_by(t["nodes"], node=0)
    assert crown["valence"] == 6 and crown["struts"] == "6B"
    ring = row_by(t["nodes"], node=1)
    assert ring["valence"] == 5 and ring["struts"] == "2A+1B+2C"
    assert [row_by(t["nodes"], node=k)["location"] for k in range(7, 13)] \
        == ["base"] * 6
    for r in t["members"]:
        assert all(math.isfinite(r[k]) for k in ("axial_i_deg", "axial_j_deg"))
        assert math.isnan(r["dihedral_deg"]) and r["panels"] == 0
    for r in t["struts"]:
        assert abs(r["cut_length_mm"] - (r["length_mm"] - 20.0)) < 1e-9
    fb.fabrication_text(t)
    with tempfile.TemporaryDirectory() as tmp:
        files = fb.export_fabrication_csv(t, tmp, prefix="star")
        names = {os.path.basename(p) for p in files}
        assert "star_fab_struts.csv" in names
        assert "star_fab_panels.csv" not in names           # empty -> skipped


# ============================================================ identities ===
def test_panel_area_sum_equals_surface_area():
    cases = [geo.generate_dome(nu, R, fr) for nu, R, fr in
             ((2, 4.0, 0.5), (3, 6.0, 0.62), (4, 10.0, 0.5), (5, 9.0, 0.45),
              (6, 15.0, 0.4))]
    cases.append(full_sphere(3, 2.0))
    for g in cases:
        t = fb.fabrication_tables(g)
        total = sum(r["total_area_m2"] for r in t["panels"])
        direct = 0.0
        for a, b, c in g.faces:
            direct += 0.5 * np.linalg.norm(np.cross(g.nodes[b] - g.nodes[a],
                                                    g.nodes[c] - g.nodes[a]))
        assert abs(total - g.surface_area) <= 1e-12 * g.surface_area
        assert abs(total - direct) <= 1e-12 * direct
        assert abs(sum(r["area_m2"] for r in t["faces"]) - direct) \
            <= 1e-12 * direct
        # each panel type: area x count = total; corner angles sum to 180
        for r in t["panels"]:
            assert abs(r["area_m2"] * r["count"] - r["total_area_m2"]) \
                <= 1e-12 * r["total_area_m2"]
            s = r["angle_1_deg"] + r["angle_2_deg"] + r["angle_3_deg"]
            assert abs(s - 180.0) < 1e-10
            assert r["size_spread_mm"] < 1e-8
        summ = {r["quantity"]: r["value"] for r in t["summary"]}
        assert summ["max dihedral: vectors - closed form"] < 1e-10
        assert summ["max axial angle: vectors - arccos(L/2R)"] < 1e-10


def test_tables_are_consistent_with_geometry():
    m = default_model(4, 8.0, 0.5)
    g = m.geom
    t = fb.fabrication_tables(m)
    sched = g.strut_schedule()
    assert [r["family"] for r in t["struts"]] == [s["label"] for s in sched]
    for r, s in zip(t["struts"], sched):
        assert r["count"] == s["count"]
        assert abs(r["length_mm"] - s["length"] * 1e3) < 1e-9
        assert abs(r["chord_factor"] - s["chord_factor"]) < 1e-14
        assert r["section"] == m.group_sections[s["id"]]
    hub_of = {r["node"]: r["hub_type"] for r in t["nodes"]}
    for r in t["members"]:
        assert r["hub_i"] == hub_of[r["node_i"]]
        assert r["hub_j"] == hub_of[r["node_j"]]
        assert r["family"] == chr(ord("A") + int(g.strut_group[r["member"]]))
    base = set(int(v) for v in g.base_nodes)
    for r in t["nodes"]:
        assert (r["location"] == "base") == (r["node"] in base)
    n_base = sum(r["count"] for r in t["hubs"] if r["location"] == "base")
    assert n_base == len(base)
    # every hub type: gaps around the hub add up to 360 deg
    for r in t["hubs"]:
        gaps = [float(x) for x in r["gaps_deg"].split()]
        assert len(gaps) == r["valence"]
        assert abs(sum(gaps) - 360.0) < 0.05
        assert r["gap_spread_deg"] < 1e-8
        if r["location"] == "base":            # open side listed last
            assert gaps[-1] == max(gaps)
            assert gaps[-1] > 170.0               # hemisphere: the ground
    # interior hubs of a 4V dome are 5- or 6-valent
    assert {r["valence"] for r in t["hubs"] if r["location"] == "interior"} \
        == {5, 6}
    # geometry-only call gives the same geometric columns
    t2 = fb.fabrication_tables(g)
    for a, b in zip(t["struts"], t2["struts"]):
        assert a["length_mm"] == b["length_mm"]
        assert a["axial_angle_deg"] == b["axial_angle_deg"]
        assert math.isnan(b["mass_each_kg"])


def test_cut_length_and_mass():
    m = default_model(3, 6.0, 0.62)
    t0 = fb.fabrication_tables(m, hub_offset=0.0)
    mass0 = sum(r["mass_total_kg"] for r in t0["struts"])
    assert abs(mass0 - m.steel_mass()) <= 1e-12 * m.steel_mass()
    o = 0.045
    t = fb.fabrication_tables(m, hub_offset=o)
    rho = m.material.rho
    for r, g in zip(t["struts"], range(m.geom.n_strut_types())):
        assert abs(r["cut_length_mm"] - (r["length_mm"] - 2e3 * o)) < 1e-9
        A = m.section_of_group(g).A
        assert abs(r["mass_each_kg"] - rho * A * r["cut_length_mm"] / 1e3) \
            < 1e-12
        assert abs(r["total_cut_length_m"]
                   - r["count"] * r["cut_length_mm"] / 1e3) < 1e-9
    for r in t["members"]:
        assert abs(r["cut_length_mm"] - (r["length_mm"] - 90.0)) < 1e-9
    summ = {r["quantity"]: r["value"] for r in t["summary"]}
    n_m = m.geom.n_members
    assert abs(summ["total cut length"]
               - (m.geom.member_lengths().sum() - 2 * o * n_m)) < 1e-9
    # hub offsets that leave no tube (or are negative) are rejected
    Lmin = float(m.geom.member_lengths().min())
    for bad in (-0.001, 0.5 * Lmin, float("nan")):
        try:
            fb.fabrication_tables(m, hub_offset=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"hub offset {bad} was accepted")
    try:
        fb.fabrication_tables("not a model")
    except TypeError:
        pass
    else:
        raise AssertionError("a string was accepted as a model")


def test_support_locations_follow_the_model():
    m = default_model(2, 5.0, 0.5)
    base = [int(v) for v in m.geom.base_nodes]
    m.support_type = geo.CUSTOM
    m.custom_supports = {base[0]: [True, True, True],
                         base[3]: [True, True, 2e7],
                         base[6]: [False, False, True]}
    t = fb.fabrication_tables(m)
    loc = {r["node"]: r["location"] for r in t["nodes"]}
    assert sorted(k for k, v in loc.items() if v == "base") == \
        sorted([base[0], base[3], base[6]])
    assert sum(1 for v in loc.values() if v == "edge") == len(base) - 3


# ================================================================ export ===
def test_csv_export_round_trip():
    m = default_model(4, 10.0, 0.5)
    t = fb.fabrication_tables(m, hub_offset=0.05, stock_length=6.0,
                              kerf=0.003)
    with tempfile.TemporaryDirectory() as tmp:
        files = fb.export_fabrication_csv(t, tmp, prefix="d4")
        names = [os.path.basename(p) for p in files]
        expect = [f"d4_fab_{k}.csv" for k in fb.TABLE_ORDER if t.get(k)]
        assert names == expect, names
        with open(os.path.join(tmp, "d4_fab_struts.csv"), newline="",
                  encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == len(t["struts"]) == 6
        assert list(rows[0].keys()) == list(t["struts"][0].keys())
        for r, ref in zip(rows, t["struts"]):
            assert r["family"] == ref["family"]
            assert int(r["count"]) == ref["count"]
            assert abs(float(r["cut_length_mm"]) - ref["cut_length_mm"]) \
                <= 0.005
            assert abs(float(r["chord_factor"]) - ref["chord_factor"]) <= 5e-9
            assert len(r["cut_length_mm"].split(".")[1]) == 2      # tidy
        with open(os.path.join(tmp, "d4_fab_members.csv"), newline="",
                  encoding="utf-8") as fh:
            mrows = list(csv.DictReader(fh))
        assert len(mrows) == m.geom.n_members
        # base struts have no dihedral: NaN is written as an empty cell
        assert any(r["dihedral_deg"] == "" for r in mrows)
        assert "nan" not in open(os.path.join(tmp, "d4_fab_members.csv"),
                                 encoding="utf-8").read().lower()
        # a model can be passed directly
        files2 = fb.export_fabrication_csv(m, os.path.join(tmp, "x"),
                                           hub_offset=0.05)
        assert len(files2) >= 7
    txt = fb.fabrication_text(t)
    assert "Strut schedule" in txt and "Hub schedule" in txt


# ============================================================ stock plan ===
def test_stock_plan_closed_form_and_ffd():
    # single-family pieces per bar n = floor((S + k)/(c + k))
    fam, pat, sec, notes = fb.stock_plan([("A", "s", 1.997, 7)], 6.0, 0.003)
    assert fam[0]["pieces_per_bar"] == 3            # 3*1.997 + 2*0.003 = 5.997
    assert fam[0]["bars"] == 3
    fam, *_ = fb.stock_plan([("A", "s", 1.999, 7)], 6.0, 0.003)
    assert fam[0]["pieces_per_bar"] == 2            # 3*1.999 + 0.006 > 6
    assert fam[0]["bars"] == 4
    # mixed plan beats one family per bar and reaches the lower bound
    fam, pat, sec, notes = fb.stock_plan(
        [("A", "s", 2.0, 5), ("B", "s", 3.5, 3)], 6.0, 0.0)
    (s,) = sec
    assert s["bars_single_family"] == 2 + 3
    assert s["bars_mixed"] == 4 and s["lower_bound"] == 4 and s["bars"] == 4
    placed = {"A": 0, "B": 0}
    for p in pat:
        assert p["used_mm"] <= 6000.0 + 1e-6
        assert p["bars"] >= 1
        for part in p["pattern"].split(" + "):
            q, lb = part.split(" x ")
            placed[lb] += int(q) * p["bars"]
    assert placed == {"A": 5, "B": 3}
    # a strut longer than the stock is flagged, not packed
    fam, pat, sec, notes = fb.stock_plan([("Z", "s", 7.0, 2)], 6.0, 0.0)
    assert fam[0]["pieces_per_bar"] == 0 and math.isnan(fam[0]["bars"])
    assert notes and "longer than" in notes[0]
    for bad in ((0.0, 0.0), (6.0, -0.001)):
        try:
            fb.stock_plan([("A", "s", 1.0, 1)], *bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"stock {bad} accepted")


def test_stock_plan_on_a_dome_is_valid():
    m = default_model(4, 10.0, 0.5)
    t = fb.fabrication_tables(m, hub_offset=0.05, stock_length=6.0,
                              kerf=0.003)
    counts = {r["family"]: r["count"] for r in t["struts"]}
    cuts = {r["family"]: r["cut_length_mm"] for r in t["struts"]}
    placed = {k: 0 for k in counts}
    n_bars = 0
    for p in t["stock_patterns"]:
        used = 0.0
        npieces = 0
        for part in p["pattern"].split(" + "):
            q, lb = part.split(" x ")
            placed[lb] += int(q) * p["bars"]
            used += int(q) * cuts[lb]
            npieces += int(q)
        used += 3.0 * (npieces - 1)
        assert abs(used - p["used_mm"]) < 1e-6
        assert used <= 6000.0 + 1e-6
        n_bars += p["bars"]
    assert placed == counts
    total = sum((cuts[k] + 3.0) * counts[k] for k in counts)
    assert n_bars >= math.ceil(total / 6003.0 - 1e-9)
    summ = [r for r in t["summary"] if r["quantity"].startswith("bars of")]
    assert summ and summ[0]["value"] == n_bars


# ========================================================== JSON models ===
def test_round_tripped_model_gives_same_schedule():
    m = default_model(3, 6.0, 0.62)
    m2 = DomeModel.from_dict(m.to_dict())
    a = fb.fabrication_tables(m, hub_offset=0.03)
    b = fb.fabrication_tables(m2, hub_offset=0.03)
    for k in ("struts", "hubs", "panels", "dihedrals"):
        assert len(a[k]) == len(b[k])
        for ra, rb in zip(a[k], b[k]):
            for key, va in ra.items():
                vb = rb[key]
                if isinstance(va, float):
                    assert (math.isnan(va) and math.isnan(vb)) or \
                        abs(va - vb) <= 1e-9 * max(1.0, abs(va)), (k, key)
                else:
                    assert va == vb, (k, key, va, vb)


# =============================================================== driver ===
def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = 0
    for name, fn in fns:
        try:
            fn()
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        else:
            npass += 1
            print(f"  pass  {name}")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
