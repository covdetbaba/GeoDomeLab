"""
Extended validation: the F1 to F8 hand calculations against GeoDomeLab.

The hand side is ``validation/hand_calc_ext.py``.  It never imports the
geodome package, and this script asks it for every reference value BEFORE
geodome is imported (the run records ``geodome loaded before the hand
values: False``), so no program number can leak into a reference.

The program side runs each structure through the normal GeoDomeLab API, the
same calls the GUI makes:

    F1  heated fixed-fixed bars       custom_geometry + CUSTOM supports +
                                      LoadSpec.temperature_change, run_analysis
    F2  Model A, uniform dT           benchmarks.validation_model_a, run_analysis
        Model A, one spoke heated     solver.solve_static(eps0 = one member)
    F3  Model A, support movements    LoadSpec.settlements, run_analysis
    F4  three-bar truss               benchmarks.three_bar_truss, run_analysis
                                      and solve_static (bar 2 heated alone);
                                      f11 from a stiffness solve of the primary
                                      structure; the program's own force-method
                                      routine (benchmarks.three_bar_force_method)
                                      is compared as well
    F5  bar + spring                  CUSTOM supports with spring entries
    F6  von Mises truss               benchmarks.von_mises_truss,
                                      nonlinear.run_nonlinear (arc length),
                                      the corotational element (TrussSystem),
                                      buckling.linear_buckling
    F7  Model A frequencies           modal.modal_analysis (lumped mass)
    F8  2V geometry                   geometry.icosahedron/subdivide_class1 +
                                      custom_geometry (full sphere) and
                                      generate_dome (hemisphere) ->
                                      fabrication.fabrication_tables

Comparison rule
---------------
For a reference value h and a program value p the difference is

    d = |p - h| / |h|            when h != 0 (plain relative difference)
    d = |p| / s                  when the exact answer is h = 0

where s is the largest |reference| with the same unit inside the same
sub-case (the natural magnitude of that quantity, e.g. the largest force of
the case for a stress-free spoke or a zero reaction), because a relative
difference is not defined for an exact zero.  A row passes when d <= tol:

    linear statics, linear buckling, geometry   tol = 1e-9
    natural frequencies, modal quantities       tol = 1e-9
    nonlinear path (limit point, snap target)   tol = 1e-6
    counts (nodes, struts, hubs, panels)        exact
    Kenner (1976) published chord factors       |p - h| <= 5e-6 (5 decimals)
    program identities (equilibrium residual,   value <= 1e-10
      energy mismatch, sum of loads + reactions)

Outputs (validation/output/):
    comparison_ext.csv   every compared quantity
    comparison_ext.md    summary per case, the full table, notes

References (the same as the hand calculations)
    Kassimali, Matrix Analysis of Structures, 2nd ed. (2012): stiffness method,
      thermal and settlement loads by equivalent joint loads.
    Hibbeler, Structural Analysis, 10th ed. (2018): force method, unit load.
    McGuire, Gallagher & Ziemian, Matrix Structural Analysis, 2nd ed. (2000):
      geometric stiffness, nonlinear solution procedures.
    Crisfield, Non-linear FE Analysis of Solids and Structures, Vol. 1 (1991):
      corotational truss, arc-length method.
    Pecknold, Ghaboussi & Healey (1985), J. Eng. Mech. 111(7): von Mises
      truss closed form.
    Chopra, Dynamics of Structures: lumped mass, effective modal mass.
    Kenner, Geodesic Math and How to Use It (1976): 2V chord factors.
    Clinton (1971), NASA CR-1734: axial angle = arccos(L/2R).

Run:  python validation/run_validation_ext.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# ---- 1. the hand side, BEFORE geodome exists in this process ---------------
import hand_calc as hc                                          # noqa: E402
import hand_calc_ext as hx                                      # noqa: E402

HAND = hx.results()
HAND_CHECKS = hx.checks()
DEFS = hx.definitions()
GEODOME_LOADED_BEFORE_HAND = "geodome" in sys.modules

# ---- 2. the program ----------------------------------------------------------
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from geodome import benchmarks as bm                            # noqa: E402
from geodome import buckling as bk                              # noqa: E402
from geodome import fabrication as fab                          # noqa: E402
from geodome import geometry as geo                             # noqa: E402
from geodome import loads as ld                                 # noqa: E402
from geodome import modal as md                                 # noqa: E402
from geodome import nonlinear as nl                             # noqa: E402
from geodome import solver as slv                               # noqa: E402
from geodome.model import DomeModel, LoadCombination, LoadSpec  # noqa: E402
from geodome.postprocess import run_analysis                    # noqa: E402

OUT = os.path.join(HERE, "output")
CSV_PATH = os.path.join(OUT, "comparison_ext.csv")
MD_PATH = os.path.join(OUT, "comparison_ext.md")

MATERIAL = bm.VAL_MATERIAL                       # ASTM A500 Gr. C (round)
SEC_MAIN = hx.SEC_MAIN[0]                        # CHS 88.9 x 4.0
SEC_SMALL = hx.SEC_SMALL[0]                      # CHS 76.1 x 4.0
PIN = [True, True, True]

TOL = {"linear": 1e-9, "modal": 1e-9, "nonlinear": 1e-6, "count": 0.0,
       "published": 5e-6, "identity": 1e-10}
TOL_TEXT = {
    "linear": "linear statics, linear buckling, geometry: d <= 1e-9",
    "modal": "natural frequencies and modal quantities: d <= 1e-9",
    "nonlinear": "nonlinear path (limit point, snap-through target): "
                 "d <= 1e-6",
    "count": "counts: exact",
    "published": "published 5-decimal chord factors: |p - h| <= 5e-6",
    "identity": "program identity (residual, energy mismatch, load "
                "balance): value <= 1e-10",
}
CASE_TITLES = {
    "F1": "Heated fixed-fixed bars",
    "F2": "Model A under a temperature change",
    "F3": "Model A with a support movement",
    "F4": "Three-bar truss by the force method",
    "F5": "Bar and elastic spring support",
    "F6": "von Mises truss snap-through",
    "F7": "Model A natural frequencies",
    "F8": "2V geodesic geometry and fabrication angles",
}


# =============================================================================
# rows
# =============================================================================
class Comparison:
    """Collects comparison rows and evaluates them (see module docstring)."""

    def __init__(self, hand):
        self.hand = hand
        self.rows = []
        self.used = set()

    def add(self, sub, key, quantity, program, source, kind="linear",
            ref=None, unit=None, scale=None):
        """Compare the hand value ``key`` (or an explicit ``ref``).

        ``scale`` is the magnitude used when the reference is exactly zero
        (default: the largest |reference| of the same unit in ``sub``).
        """
        if ref is None:
            h = self.hand[key]
            ref, unit = h["value"], h["unit"] if unit is None else unit
            self.used.add(key)
        self.rows.append({
            "case": sub.split()[0][:2], "sub": sub, "key": key,
            "quantity": quantity, "hand": float(ref),
            "program": float(program), "unit": unit or "-", "kind": kind,
            "source": source, "scale": scale})

    def identity(self, sub, label, res):
        """Equilibrium, energy and load-balance identities of a solve."""
        F = np.asarray(res.F, float)
        R = np.asarray(res.reactions, float)
        scale = max(np.linalg.norm(F), np.linalg.norm(R), 1e-30)
        bal = float(np.linalg.norm(res.reaction_balance)) / scale
        for q, v in ((f"{label}: equilibrium residual", res.equilibrium_residual),
                     (f"{label}: energy identity mismatch", res.energy_error),
                     (f"{label}: |sum F + sum R| / max(|F|, |R|)", bal)):
            self.rows.append({"case": sub[:2], "sub": sub, "key": "identity",
                              "quantity": q, "hand": 0.0, "program": float(v),
                              "unit": "-", "kind": "identity",
                              "source": "solver.StaticResult"})

    def evaluate(self):
        groups = {}
        for r in self.rows:
            if r["kind"] in ("linear", "modal", "nonlinear"):
                g = (r["sub"], r["unit"])
                groups[g] = max(groups.get(g, 0.0), abs(r["hand"]))
        for r in self.rows:
            h, p, kind = r["hand"], r["program"], r["kind"]
            err = abs(p - h)
            if kind == "identity":
                d, scale = abs(p), 1.0
            elif kind in ("count", "published"):
                d, scale = err, 1.0
            elif h != 0.0:
                scale = abs(h)                       # plain relative
                d = err / scale
            else:                                    # exact answer is zero
                scale = r["scale"] or groups[(r["sub"], r["unit"])]
                if not scale > 0.0:
                    raise ValueError(f"No magnitude for the zero reference "
                                     f"{r['sub']} {r['quantity']}.")
                d = err / scale
            r["scale"] = scale
            r["difference"] = d
            r["rel_diff"] = err / abs(h) if h != 0.0 else float("nan")
            r["tolerance"] = TOL[kind]
            r["pass"] = bool(np.isfinite(d) and d <= TOL[kind])
        return self.rows


# =============================================================================
# helpers
# =============================================================================
def _loads(**kw):
    return bm._no_area_loads(**kw)


def custom_model(nodes, members, sections, supports, groups=None, loads=None,
                 combos=None, title="F case", material=MATERIAL):
    """A DomeModel of an arbitrary truss with CUSTOM supports."""
    groups = list(range(len(members))) if groups is None else groups
    g = geo.custom_geometry(nodes, members, base_nodes=sorted(supports),
                            groups=groups, title=title)
    return DomeModel(geom=g, material_name=material,
                     group_sections=list(sections), support_type=geo.CUSTOM,
                     loads=loads or _loads(),
                     combinations=combos or [LoadCombination("T  1.0T",
                                                             T=1.0)],
                     custom_supports={int(k): list(v)
                                      for k, v in supports.items()})


def analyse(model, combo_name=None):
    """run_analysis (no buckling) and return the result of one combination."""
    res = run_analysis(model, do_buckling=False)
    name = combo_name or model.combinations[0].name
    return res.cases[name]


def model_vectors(model):
    fixed, springs = model.support_conditions()
    return fixed, springs


def node_map(hand_nodes, prog_nodes):
    """hand node number -> program node index, matched by coordinates."""
    out = {}
    for h, xyz in hand_nodes.items():
        d = np.linalg.norm(prog_nodes - np.asarray(xyz, float), axis=1)
        k = int(np.argmin(d))
        if d[k] > 1e-9:
            raise RuntimeError(f"hand node {h} has no program twin "
                               f"({d[k]:.2e} m)")
        out[h] = k
    return out


def member_map(hand_members, prog_members, nmap):
    lookup = {frozenset((int(i), int(j))): k
              for k, (i, j) in enumerate(prog_members)}
    return {m: lookup[frozenset((nmap[i], nmap[j]))]
            for m, (i, j) in hand_members.items()}


def model_a_maps(model):
    nodes, members = hc.model_a_geometry()
    nmap = node_map(nodes, model.geom.nodes)
    mmap = member_map(members, model.geom.members, nmap)
    return nmap, mmap


def model_a_T_only(temperature=0.0, settlements=None):
    """Validation Model A with only a self-straining load (combination T)."""
    m = bm.validation_model_a()
    m.loads.temperature_change = float(temperature)
    m.loads.settlements = list(settlements or [])
    m.combinations = [LoadCombination("T  1.0T", T=1.0)]
    return m


# =============================================================================
# F1  heated fixed-fixed bars
# =============================================================================
def case_f1(C):
    d = DEFS["F1a"]
    m = custom_model(d["nodes"], d["members"], [SEC_MAIN],
                     {0: PIN, 1: PIN},
                     loads=_loads(temperature_change=d["dT"]),
                     title="F1a heated inclined bar")
    r = analyse(m)
    src = "custom_geometry + CUSTOM pins, LoadSpec.temperature_change, " \
          "run_analysis (T = 1.0)"
    C.add("F1a", "F1.L", "bar length L", m.geom.strut_lengths[0],
          "custom_geometry")
    C.add("F1a", "F1.N", "axial force N", r.N[0], src)
    C.add("F1a", "F1.sigma", "axial stress N/A", r.stress[0], src)
    Rv = r.reaction_vectors()
    for a, lab in enumerate("xyz"):
        C.add("F1a", f"F1.R1_{lab}", f"reaction R{lab} at node 1",
              Rv[0, a], src)
        C.add("F1a", f"F1.R2_{lab}", f"reaction R{lab} at node 2",
              Rv[1, a], src)
    C.add("F1a", "F1.free_elongation", "free thermal elongation eps0 L",
          r.eps0[0] * r.lengths[0], "StaticResult.eps0 x L")
    sec = m.member_sections()[0]
    Pe = bk.euler_load(sec, m.material, float(m.geom.strut_lengths[0]), 1.0)
    C.add("F1a", "F1.dT_euler", "temperature rise at the Euler load",
          Pe / (m.EA()[0] * m.material.alpha),
          "buckling.euler_load / (EA alpha)")
    C.identity("F1a", "F1a", r)

    d = DEFS["F1b"]
    m = custom_model(d["nodes"], d["members"], [SEC_MAIN, SEC_SMALL],
                     {0: PIN, 1: [False, True, True], 2: PIN},
                     loads=_loads(temperature_change=d["dT"]),
                     title="F1b heated stepped bar")
    r = analyse(m)
    src = "two sections, run_analysis (T = 1.0)"
    C.add("F1b", "F1.stepped_N", "force in segment 1", r.N[0], src)
    C.add("F1b", "F1.stepped_N", "force in segment 2", r.N[1], src)
    C.add("F1b", "F1.stepped_u2", "x displacement of the middle node",
          r.u[3], src)
    C.add("F1b", "F1.stepped_R1x", "x reaction at node 1", r.reactions[0],
          src)
    C.identity("F1b", "F1b", r)


# =============================================================================
# F2  Model A under a temperature change
# =============================================================================
def case_f2(C):
    dT = hx.DT
    m = model_a_T_only(temperature=dT)
    nmap, mmap = model_a_maps(m)
    c = nmap[1]
    r = analyse(m)
    src = "validation_model_a, LoadSpec.temperature_change, run_analysis"
    C.add("F2a", "F2.crown_uz", "crown uz", r.u[3 * c + 2], src)
    C.add("F2a", "F2.crown_ux", "crown ux", r.u[3 * c], src)
    C.add("F2a", "F2.crown_ux", "crown uy (symmetry, 0)", r.u[3 * c + 1], src)
    for k in range(1, 6):
        C.add("F2a", "F2.N_spoke", f"N spoke {k}", r.N[mmap[k]], src)
    for k in range(6, 11):
        C.add("F2a", "F2.N_ring", f"N ring member {k}", r.N[mmap[k]], src)
    Rv = r.reaction_vectors()
    for n in range(2, 7):
        th = math.radians(72.0 * (n - 2))
        Rn = Rv[nmap[n]]
        C.add("F2a", "F2.R_radial", f"radial reaction at hand node {n}",
              Rn[0] * math.cos(th) + Rn[1] * math.sin(th), src)
        C.add("F2a", "F2.R_z", f"vertical reaction at hand node {n}", Rn[2],
              src)
        for a, lab in enumerate("xyz"):
            C.add("F2a", f"F2.R{n}_{lab}", f"R{lab} at hand node {n}", Rn[a],
                  src)
    C.identity("F2a", "F2a", r)

    # one spoke heated: per-member initial strain through solve_static
    fixed, springs = model_vectors(m)
    eps0 = np.zeros(m.geom.n_members)
    eps0[mmap[1]] = m.material.alpha * dT
    r = slv.solve_static(m.geom.nodes, m.geom.members, m.EA(),
                         np.zeros(3 * m.geom.n_nodes), fixed, name="F2b",
                         eps0=eps0, springs=springs or None)
    src = "solver.solve_static, eps0 on spoke 1 only"
    for a, lab in enumerate("xyz"):
        C.add("F2b", f"F2.oneSpoke_crown_u{lab}", f"crown u{lab}",
              r.u[3 * c + a], src)
    for k in range(1, 6):
        C.add("F2b", f"F2.oneSpoke_N{k}", f"N spoke {k}", r.N[mmap[k]], src)
    C.identity("F2b", "F2b", r)


# =============================================================================
# F3  Model A support movements
# =============================================================================
def case_f3(C):
    for sub, tag, key in (("F3a", "vert", "settlement_F3a"),
                          ("F3b", "spread", "spread_F3b")):
        mv = DEFS["F3"][key]
        hand_node = int(mv["node"]) + 1               # definitions: 0-based
        m0 = bm.validation_model_a()
        nmap, mmap = model_a_maps(m0)
        settle = [{"node": nmap[hand_node], "ux": mv["ux"], "uy": mv["uy"],
                   "uz": mv["uz"]}]
        m = model_a_T_only(settlements=settle)
        r = analyse(m)
        src = "validation_model_a, LoadSpec.settlements, run_analysis"
        c = nmap[1]
        for a, lab in enumerate("xyz"):
            C.add(sub, f"F3.{tag}_crown_u{lab}", f"crown u{lab}",
                  r.u[3 * c + a], src)
        for k in range(1, 11):
            C.add(sub, f"F3.{tag}_N{k}", f"N member {k}", r.N[mmap[k]], src)
        Rv = r.reaction_vectors()
        for n in range(2, 7):
            for a, lab in enumerate("xyz"):
                C.add(sub, f"F3.{tag}_R{n}_{lab}", f"R{lab} at hand node {n}",
                      Rv[nmap[n], a], src)
        C.add(sub, f"F3.{tag}_energy", "strain energy sum N^2 L/(2EA)",
              r.strain_energy_int, "StaticResult.strain_energy_int")
        C.add(sub, f"F3.{tag}_energy", "work 1/2 R . us (Clapeyron)",
              r.strain_energy_ext, "StaticResult.strain_energy_ext")
        C.identity(sub, sub, r)


# =============================================================================
# F4  three-bar truss
# =============================================================================
def case_f4(C):
    d = DEFS["F4"]
    alpha = hx.ALPHA
    dT_bar = np.zeros(3)
    dT_bar[d["temperature"]["member"]] = d["temperature"]["dT"]
    s = d["settlement"]
    # the preset with no temperature, support C settling: "P", "T", "P + T"
    m = bm.three_bar_truss(temperature_change=0.0, settlement_uz=s["uz"])
    assert np.allclose(m.geom.nodes, np.asarray(d["nodes"])), "F4 geometry"
    res = run_analysis(m, do_buckling=False)
    fixed, springs = model_vectors(m)
    basic = ld.basic_load_vectors(m)
    FP = ld.combine(basic, m.combinations[0])
    eps_t = alpha * dT_bar
    pres = {3 * int(s["node"]) + 2: s["uz"]}
    runs = {
        "load": (res.cases["P  1.0L"], "three_bar_truss, run_analysis 'P'"),
        "settlement": (res.cases["T  1.0T  temperature + settlement"],
                       "three_bar_truss (dT = 0), run_analysis 'T'"),
        "temperature": (slv.solve_static(
            m.geom.nodes, m.geom.members, m.EA(), np.zeros_like(FP), fixed,
            name="T bar 2", eps0=eps_t, springs=springs or None),
            "solver.solve_static, eps0 on bar 2 only"),
        "combined": (slv.solve_static(
            m.geom.nodes, m.geom.members, m.EA(), FP, fixed, name="P+T+S",
            eps0=eps_t, prescribed=pres, springs=springs or None),
            "solver.solve_static, load + eps0 (bar 2) + settlement"),
    }
    D = 3
    for cn, (r, src) in runs.items():
        sub = f"F4 {cn}"
        for i in (1, 2, 3):
            C.add(sub, f"F4.N{i}_{cn}", f"N bar {i}", r.N[i - 1], src)
        C.add(sub, f"F4.ux_{cn}", "ux of joint D", r.u[3 * D], src)
        C.add(sub, f"F4.uz_{cn}", "uz of joint D", r.u[3 * D + 2], src)
        Rv = r.reaction_vectors()
        for k, nm in enumerate(("A", "B", "C")):
            C.add(sub, f"F4.R{nm}x_{cn}", f"Rx at support {nm}", Rv[k, 0], src)
            C.add(sub, f"F4.R{nm}z_{cn}", f"Rz at support {nm}", Rv[k, 2], src)
        C.identity(sub, sub, r)

    # flexibility coefficient by the stiffness method on the primary structure
    # (bars 1 and 3): f11 = L2/EA2 + c2 . u,  K13 u = c2 at joint D.
    X = m.geom.nodes
    EA = m.EA()
    keep = [0, 2]
    L2 = float(np.linalg.norm(X[3] - X[1]))
    c2 = (X[3] - X[1]) / L2
    F = np.zeros(3 * len(X))
    F[3 * D:3 * D + 3] = c2
    rp = slv.solve_static(X, m.geom.members[keep], EA[keep], F, fixed,
                          name="primary")
    C.add("F4 flexibility", "F4.f11",
          "f11 = L2/EA2 + c2 . K13^-1 c2 (primary structure)",
          L2 / EA[1] + float(c2 @ rp.u[3 * D:3 * D + 3]),
          "solver.solve_static on bars 1 and 3")

    # the program's own force-method reference (used by verify_preset)
    P = (d["load"]["Fx"], d["load"]["Fz"])
    sett = {int(s["node"]): (s["ux"], s["uz"])}
    fm_cases = {"load": dict(load=P),
                "temperature": dict(dT=dT_bar),
                "settlement": dict(settlement=sett),
                "combined": dict(load=P, dT=dT_bar, settlement=sett)}
    for cn, kw in fm_cases.items():
        fm = bm.three_bar_force_method(m, **kw)
        sub = f"F4 force-method {cn}"
        src = "benchmarks.three_bar_force_method"
        for i in (1, 2, 3):
            C.add(sub, f"F4.N{i}_{cn}", f"N bar {i}", fm["N"][i - 1], src)
        C.add(sub, f"F4.ux_{cn}", "ux of joint D", fm["u"][0], src)
        C.add(sub, f"F4.uz_{cn}", "uz of joint D", fm["u"][1], src)
    C.add("F4 force-method load", "F4.f11", "f11", fm["f11"],
          "benchmarks.three_bar_force_method")


# =============================================================================
# F5  bar + spring
# =============================================================================
def case_f5(C):
    P = hx.F5_P
    load = _loads(nodal_loads=[{"node": 1, "Fx": P, "Fy": 0.0, "Fz": 0.0,
                                "case": "L"}])
    cP = [LoadCombination("P  1.0L", L=1.0)]
    for sub, key in (("F5a", "F5a"), ("F5b", "F5b"), ("F5c", "F5c")):
        d = DEFS[key]
        sup = {int(k): list(v) for k, v in d["supports"].items()}
        if sub == "F5c":
            m = custom_model(d["nodes"], d["members"], [SEC_MAIN], sup,
                             loads=_loads(temperature_change=d["dT"]),
                             title="F5c heated bar on a spring")
        else:
            m = custom_model(d["nodes"], d["members"], [SEC_MAIN], sup,
                             loads=load, combos=cP, title=f"{sub} bar+spring")
        r = analyse(m)
        src = "CUSTOM supports with a spring entry, run_analysis"
        if sub == "F5a":
            C.add(sub, "F5.series_u1", "u1 (spring node)", r.u[0], src)
            C.add(sub, "F5.series_u2", "u2 (loaded node)", r.u[3], src)
            C.add(sub, "F5.series_N", "bar force", r.N[0], src)
            C.add(sub, "F5.series_spring_reaction", "spring reaction",
                  r.reactions[0], src)
            C.add(sub, "F5.series_spring_reaction",
                  "StaticResult.spring_forces", r.spring_forces[0], src)
        elif sub == "F5b":
            C.add(sub, "F5.parallel_u2", "u2", r.u[3], src)
            C.add(sub, "F5.parallel_N", "bar force", r.N[0], src)
            C.add(sub, "F5.parallel_spring_reaction", "spring reaction",
                  r.reactions[3], src)
            C.add(sub, "F5.parallel_R1x", "fixed reaction at node 1",
                  r.reactions[0], src)
        else:
            C.add(sub, "F5.thermal_u2", "u2", r.u[3], src)
            C.add(sub, "F5.thermal_N", "bar force", r.N[0], src)
            C.add(sub, "F5.thermal_spring_reaction", "spring reaction",
                  r.reactions[3], src)
        C.identity(sub, sub, r)


# =============================================================================
# F6  von Mises truss
# =============================================================================
def _newton_at_load(system, F, lam, u0, tol=1e-14, max_iter=50):
    """Equilibrium f(u) = lam F of the program's corotational system."""
    u = np.array(u0, dtype=float)
    Ff = F[system.free]
    for _ in range(max_iter):
        st = system.state(u)
        r = lam * Ff - st.f_free
        if np.linalg.norm(r) <= tol * abs(lam) * np.linalg.norm(Ff):
            return st
        K = system.tangent(st)
        K = K.toarray() if hasattr(K, "toarray") else np.asarray(K)
        u[system.free] += np.linalg.solve(K, r)
    raise RuntimeError("Newton did not converge on the program residual.")


def _brent_root(f, a, b, xtol=1e-15, maxit=200):
    """Bracketed root by bisection + secant safeguard (no SciPy needed)."""
    fa, fb = f(a), f(b)
    if fa * fb > 0.0:
        raise ValueError("root not bracketed")
    for _ in range(maxit):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or 0.5 * (b - a) < xtol:
            return m
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def case_f6(C):
    m = bm.von_mises_truss()
    a, h, EA, P_ref = bm.von_mises_parameters(m)
    apex = 2
    iz = 3 * apex + 2
    combo = m.combinations[0]
    F = ld.combine(ld.basic_load_vectors(m), combo)
    C.add("F6 linear", "F6.L", "bar length", m.geom.strut_lengths[0],
          "custom_geometry")
    C.add("F6 linear", "F6.EA", "axial rigidity EA", EA, "DomeModel.EA")
    C.add("F6 linear", "F6.P_ref", "reference load (downward)", -F[iz],
          "LoadSpec.nodal_loads -> loads.combine")
    r = analyse(m)
    C.add("F6 linear", "F6.N_ref", "linear bar force under P_ref", r.N[0],
          "run_analysis")
    C.add("F6 linear", "F6.K0", "initial stiffness P_ref / w_lin",
          P_ref / abs(r.u[iz]), "run_analysis")
    C.identity("F6 linear", "F6 linear", r)

    b = bk.linear_buckling(m, r.N, n_modes=6)
    lam_cr = float(b.load_factors[0])
    C.add("F6 linear", "F6.lambda_cr", "linear buckling factor",
          lam_cr, "buckling.linear_buckling")
    C.add("F6 linear", "F6.P_cr_linear", "lambda_cr P_ref", lam_cr * P_ref,
          "buckling.linear_buckling")
    # linear_buckling keeps n_free - 2 modes (one here); the sway factor
    # comes from the same program matrices K and K_g on the free DOFs.
    fixed, _ = m.support_conditions()
    fd = slv.free_dofs(3 * m.geom.n_nodes, fixed)
    K = slv.assemble_stiffness(m.geom.nodes, m.geom.members, m.EA())
    Kg = slv.assemble_geometric_stiffness(m.geom.nodes, m.geom.members, r.N)
    K = (K.toarray() if hasattr(K, "toarray") else K)[np.ix_(fd, fd)]
    Kg = (Kg.toarray() if hasattr(Kg, "toarray") else Kg)[np.ix_(fd, fd)]
    mu = np.linalg.eigvals(np.linalg.solve(K, -Kg)).real
    lams = np.sort(1.0 / mu[mu > 1e-14])
    C.add("F6 linear", "F6.lambda_cr", "lowest root of det(K + lam Kg)",
          lams[0], "solver.assemble_stiffness / assemble_geometric_stiffness")
    C.add("F6 linear", "F6.lambda_sway", "sway root of det(K + lam Kg)",
          lams[1], "solver.assemble_stiffness / assemble_geometric_stiffness")

    # ---- nonlinear path --------------------------------------------------
    t0 = time.perf_counter()
    # first increment 0.05 (lambda* = 1.024) gives a densely sampled path
    res = nl.run_nonlinear(m, combo.name, method="arc-length",
                           control="crown_uz", first_increment=0.05,
                           stop_after_limit=None, control_max=3.2 * h,
                           max_steps=600)
    src = (f"run_nonlinear arc-length ({res.n_points - 1} steps, "
           f"{time.perf_counter() - t0:.2f} s)")
    C.add("F6 nonlinear", "F6.P_limit", "limit load lambda_limit P_ref",
          res.lambda_limit * P_ref, src, kind="nonlinear")
    C.add("F6 nonlinear", "F6.w_limit", "apex deflection at the limit",
          -res.control_limit, src, kind="nonlinear")
    st = res.u_limit
    sys_ = nl.model_system(m)
    stl = sys_.state(st)
    C.add("F6 nonlinear", "F6.l_limit", "bar length at the limit",
          stl.l[0], src, kind="nonlinear")
    C.add("F6 nonlinear", "F6.N_limit", "bar force at the limit",
          res.N_limit[0], src, kind="nonlinear")
    C.add("F6 nonlinear", "F6.strain_limit", "strain at the limit",
          stl.eps[0], src, kind="nonlinear")
    C.add("F6 nonlinear", "F6.ratio_linear_to_limit",
          "lambda_cr / lambda_limit", res.lambda_cr / res.lambda_limit, src,
          kind="nonlinear")
    C.add("F6 nonlinear", "F6.lambda_cr", "lambda_cr reported by "
          "run_nonlinear", res.lambda_cr, "run_nonlinear (linear_buckling)")
    # the whole traced path against the closed form P(w) of the hand file
    w_path = -res.control
    P_exact = np.array([hx.vm_P(float(w)) for w in w_path])
    Ps = HAND["F6.P_limit"]["value"]
    dev = float(np.max(np.abs(res.lam * P_ref - P_exact))) / Ps
    w1h = HAND["F6.w_limit"]["value"]
    w2h = HAND["F6.w_second_limit"]["value"]
    n_unstable = int(np.sum((w_path > w1h) & (w_path < w2h)))
    C.add("F6 nonlinear", "path", f"max |P_program - P_exact(w)| / P* over "
          f"all {res.n_points} path points (w up to "
          f"{w_path.max() / h:.2f} h, {n_unstable} on the unstable branch)",
          dev, src, kind="nonlinear", ref=0.0, unit="-", scale=1.0)
    C.add("F6 nonlinear", "path", "path passes the snap-through target "
          "(max w / w_snap >= 1)",
          float(w_path.max() >= HAND["F6.w_snap"]["value"]), src,
          kind="count", ref=1.0, unit="-")

    # limit points from the program's consistent tangent: under symmetric
    # motion dP/dw = K_T,zz, so both limit points are its zeros.
    def ktzz(w):
        u = np.zeros(sys_.ndof)
        u[iz] = -w
        Kt = sys_.tangent(sys_.state(u))
        Kt = Kt.toarray() if hasattr(Kt, "toarray") else Kt
        return Kt[sys_.fmap[iz], sys_.fmap[iz]]
    w1 = _brent_root(ktzz, 1e-6 * h, h)
    w2 = _brent_root(ktzz, h, 2.0 * h - 1e-6 * h)
    src_t = "zero of the program's tangent K_T,zz (TrussSystem.tangent)"
    C.add("F6 tangent", "F6.w_limit", "first limit point w*", w1, src_t,
          kind="nonlinear")
    C.add("F6 tangent", "F6.w_second_limit", "second limit point h + y*",
          w2, src_t, kind="nonlinear")

    # P(w) at the sample deflections: program internal force at the
    # symmetric state (apex ux = 0 is an equilibrium of the x equation)
    for rr in hx.VM_SAMPLES:
        w = rr * h
        u = np.zeros(sys_.ndof)
        u[iz] = -w
        stw = sys_.state(u)
        C.add("F6 path", f"F6.P_at_w{rr:g}h", f"P at w = {rr:g} h",
              -stw.f_full[iz], "TrussSystem.state (corotational element)",
              kind="nonlinear")
    # snap-through target: Newton on the program residual at lambda*,
    # started from the traced point on the inverted branch nearest lambda*
    lam_s = Ps / P_ref
    far = np.where(w_path > 2.0 * h)[0]
    k0 = int(far[np.argmin(np.abs(res.lam[far] - lam_s))])
    stn = _newton_at_load(sys_, res.F_ref, lam_s, res.u[k0])
    src_s = "Newton on the program residual at lambda* from the traced path"
    C.add("F6 snap", "F6.w_snap", "apex deflection after snap-through",
          -stn.u[iz], src_s, kind="nonlinear")
    C.add("F6 snap", "F6.N_snap", "bar force after snap-through", stn.N[0],
          src_s, kind="nonlinear")
    return res


# =============================================================================
# F7  Model A natural frequencies
# =============================================================================
def case_f7(C):
    m = bm.validation_model_a()
    nmap, _ = model_a_maps(m)
    c = nmap[1]
    for tag, dead in (("sw", False), ("dead", True)):
        sub = f"F7 {tag}"
        r = md.modal_analysis(m, n_modes=3, include_dead=dead)
        src = f"modal.modal_analysis (lumped, include_dead={dead})"
        C.add(sub, f"F7.m_crown_{tag}", "lumped crown mass",
              r.nodal_mass[c, 0], src, kind="modal")
        C.add(sub, f"F7.m_crown_{tag}", "free mass in x (= crown mass)",
              r.free_mass[0], src, kind="modal")
        dirs = {r.dominant_direction(k): k for k in range(r.n_modes)}
        assert set(dirs) == {"x", "y", "z"}, dirs
        for ax in "xy":
            k = dirs[ax]
            C.add(sub, f"F7.f_xy_{tag}", f"f of the {ax} mode",
                  r.frequencies[k], src, kind="modal")
            C.add(sub, f"F7.omega_xy_{tag}", f"omega of the {ax} mode",
                  r.omega[k], src, kind="modal")
        kz = dirs["z"]
        C.add(sub, f"F7.f_z_{tag}", "f of the z mode", r.frequencies[kz],
              src, kind="modal")
        C.add(sub, f"F7.omega_z_{tag}", "omega of the z mode", r.omega[kz],
              src, kind="modal")
        C.add(sub, f"F7.T_z_{tag}", "period of the z mode", r.periods[kz],
              src, kind="modal")
        for ax in "xyz":
            k = dirs[ax]
            C.add(sub, f"F7.phi_{tag}", f"max |phi| of the {ax} mode",
                  float(np.max(np.abs(r.modes[:, k]))), src, kind="modal")
            C.add(sub, f"F7.meff_{ax}_{tag}",
                  f"effective mass of the {ax} mode in {ax}",
                  r.effective_mass[k, "xyz".index(ax)], src, kind="modal")
        for q, v in (("M-orthonormality error max|Phi^T M Phi - I|",
                      r.orthogonality_error),
                     ("eigen residual (normwise backward error)",
                      r.residual)):
            C.rows.append({"case": "F7", "sub": sub, "key": "identity",
                           "quantity": q, "hand": 0.0, "program": float(v),
                           "unit": "-", "kind": "identity",
                           "source": "modal.ModalResult"})


# =============================================================================
# F8  2V geometry
# =============================================================================
def full_sphere_2v(R):
    """Full 2V Class I sphere from the program's own geometry routines."""
    v, f = geo.icosahedron()
    pts, tris = geo.subdivide_class1(v, f, 2)
    members = geo.edges_from_faces(tris)
    return geo.custom_geometry(pts * R, members, faces=tris,
                               title="2V Class I full sphere")


def _f8_angles(C, sub, t, R, src, counts):
    st = {r["family"]: r for r in t["struts"]}
    dh = {r["strut_family"]: r for r in t["dihedrals"]}
    for fam in "AB":
        C.add(sub, f"F8.CF_{fam}", f"chord factor {fam}",
              st[fam]["chord_factor"], src)
        C.add(sub, f"F8.L_{fam}_R5", f"strut {fam} length (R = 5 m)",
              st[fam]["length_mm"] / 1000.0 * hx.F8_R / R, src)
        C.add(sub, f"F8.axial_{fam}_deg", f"axial angle {fam}",
              st[fam]["axial_angle_deg"], src)
        C.add(sub, f"F8.bevel_{fam}_deg", f"tangent (bevel) angle {fam}",
              st[fam]["tangent_angle_deg"], src)
        C.add(sub, f"F8.dihedral_{fam}_deg", f"dihedral angle across {fam}",
              dh[fam]["dihedral_deg"], src)
        C.add(sub, f"F8.bend_{fam}_deg", f"bend angle across {fam}",
              dh[fam]["bend_deg"], src)
        C.add(sub, f"F8.CF_{fam}", f"chord factor {fam} vs Kenner (1976)",
              st[fam]["chord_factor"], src, kind="published",
              ref=hx.KENNER_2V[fam], unit="-")
    pn = {r["panel_type"]: r for r in t["panels"]}
    aab = pn["AAB"]
    # sides are ordered by family (A, A, B); angle k is opposite side k
    C.add(sub, "F8.face_AAB_apex_deg", "AAB apex angle (opposite B)",
          aab["angle_3_deg"], src)
    C.add(sub, "F8.face_AAB_base_deg", "AAB base angle 1", aab["angle_1_deg"],
          src)
    C.add(sub, "F8.face_AAB_base_deg", "AAB base angle 2", aab["angle_2_deg"],
          src)
    C.add(sub, "F8.area_AAB_R1", "AAB area / R^2", aab["area_m2"] / R ** 2,
          src)
    C.add(sub, "F8.area_BBB_R1", "BBB area / R^2",
          pn["BBB"]["area_m2"] / R ** 2, src)
    if not counts:
        return
    C.add(sub, "F8.area_total_R1", "total panel area / R^2",
          sum(r["total_area_m2"] for r in t["panels"]) / R ** 2, src)
    cnt = {r["family"]: r["count"] for r in t["struts"]}
    C.add(sub, "F8.struts_A", "A struts", cnt["A"], src, kind="count")
    C.add(sub, "F8.struts_B", "B struts", cnt["B"], src, kind="count")
    C.add(sub, "F8.panels_AAB", "AAB panels", pn["AAB"]["count"], src,
          kind="count")
    C.add(sub, "F8.panels_BBB", "BBB panels", pn["BBB"]["count"], src,
          kind="count")
    hubs = t["hubs"]
    C.add(sub, "F8.hubs_valence5", "5-valent hubs",
          sum(r["count"] for r in hubs if r["valence"] == 5), src,
          kind="count")
    C.add(sub, "F8.hubs_valence6", "6-valent hubs",
          sum(r["count"] for r in hubs if r["valence"] == 6), src,
          kind="count")


def case_f8(C):
    R = hx.F8_R
    g = full_sphere_2v(R)
    C.add("F8 sphere", "F8.nodes", "nodes", g.n_nodes,
          "icosahedron + subdivide_class1 + custom_geometry", kind="count")
    C.add("F8 sphere", "F8.struts", "struts", g.n_members,
          "edges_from_faces", kind="count")
    C.add("F8 sphere", "F8.panels", "panels", g.n_faces,
          "subdivide_class1", kind="count")
    t = fab.fabrication_tables(g)
    _f8_angles(C, "F8 sphere", t, R, "fabrication_tables (full 2V sphere)",
               counts=True)
    gd = geo.generate_dome(2, R, 0.5)
    td = fab.fabrication_tables(gd)
    _f8_angles(C, "F8 dome", td, R,
               "generate_dome(2V, 0.5) -> fabrication_tables", counts=False)
    # the formula columns the program prints next to the vector values
    for r in t["struts"]:
        C.add("F8 sphere", f"F8.axial_{r['family']}_deg",
              f"axial angle {r['family']}, program formula column",
              r["axial_formula_deg"], "fabrication_tables axial_formula_deg")
    for r in t["dihedrals"]:
        C.add("F8 sphere", f"F8.dihedral_{r['strut_family']}_deg",
              f"dihedral {r['strut_family']}, program formula column",
              r["dihedral_formula_deg"],
              "fabrication_tables dihedral_formula_deg")


# =============================================================================
# driver
# =============================================================================
CASES = (case_f1, case_f2, case_f3, case_f4, case_f5, case_f6, case_f7,
         case_f8)


def run_comparison():
    """Run every case; return (rows, info dict)."""
    C = Comparison(HAND)
    t0 = time.perf_counter()
    for fn in CASES:
        fn(C)
    rows = C.evaluate()
    missing = sorted(set(HAND) - C.used)
    info = {"elapsed": time.perf_counter() - t0, "missing": missing,
            "n_hand": len(HAND), "n_used": len(C.used),
            "hand_checks": HAND_CHECKS,
            "geodome_before_hand": GEODOME_LOADED_BEFORE_HAND}
    return rows, info


def _fmt(x, nd=10):
    x = float(x)
    if not np.isfinite(x):
        return "-"
    if x == 0.0:
        return "0"
    return f"{x:.{nd}g}"


def write_csv(rows, path=CSV_PATH):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = ["case", "sub", "key", "quantity", "hand", "program", "unit",
            "difference", "rel_diff", "scale", "tolerance", "kind", "pass",
            "source"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            out = {k: r[k] for k in cols}
            for k in ("hand", "program", "difference", "rel_diff", "scale"):
                out[k] = repr(float(r[k]))
            out["pass"] = "yes" if r["pass"] else "NO"
            w.writerow(out)
    return path


def summary_by_case(rows):
    out = {}
    for r in rows:
        s = out.setdefault(r["case"], {"n": 0, "pass": 0, "worst": {}})
        s["n"] += 1
        s["pass"] += int(r["pass"])
        w = s["worst"].get(r["kind"], 0.0)
        s["worst"][r["kind"]] = max(w, r["difference"])
    return out


def write_markdown(rows, info, path=MD_PATH):
    npass = sum(r["pass"] for r in rows)
    summ = summary_by_case(rows)
    nh_ok = sum(1 for _, _, res, tol in info["hand_checks"] if res <= tol)
    L = ["# GeoDomeLab extended validation: hand calculations F1-F8 vs the "
         "program", "",
         f"{npass} of {len(rows)} comparisons pass. "
         f"{info['n_used']} of {info['n_hand']} hand reference values are "
         f"compared (every one of them"
         + ("" if not info["missing"] else
            f"; NOT compared: {', '.join(info['missing'])}") + "). "
         f"The hand file's own {len(info['hand_checks'])} internal checks: "
         f"{nh_ok} within tolerance.", "",
         "Produced by `validation/run_validation_ext.py`. The reference values "
         "come from `validation/hand_calc_ext.py`, which does not import the "
         "geodome package; they were computed in this run before geodome was "
         f"imported (geodome loaded before the hand values: "
         f"{info['geodome_before_hand']}). The program values come from the "
         "normal GeoDomeLab API (column `source`).", "",
         "## Comparison rule", "",
         "For a non-zero reference d = |program - hand| / |hand|, the plain "
         "relative difference. Where the exact answer is zero (a stress-free "
         "spoke, a zero reaction, a symmetric displacement) a relative "
         "difference is not defined, so d = |program| / s with s the largest "
         "|hand value| of the same unit in the same sub-case (column "
         "`scale` of the CSV). A row passes when d is within its tolerance:",
         ""]
    L += [f"- {TOL_TEXT[k]}" for k in TOL_TEXT]
    L += ["", "## Summary by case", "",
          "| Case | Structure | Rows | Pass | Worst d (linear) | Worst d "
          "(modal) | Worst d (nonlinear) | Worst identity |",
          "|---|---|---|---|---|---|---|---|"]
    for case in sorted(summ):
        s = summ[case]
        w = s["worst"]

        def g(k):
            return f"{w[k]:.1e}" if k in w else "-"
        L.append(f"| {case} | {CASE_TITLES.get(case, '')} | {s['n']} | "
                 f"{s['pass']} | {g('linear')} | {g('modal')} | "
                 f"{g('nonlinear')} | {g('identity')} |")
    L += ["", "## How each case is run in GeoDomeLab", ""]
    L += ["- **F1** custom 2-node and 3-node trusses (`custom_geometry`, "
          "CUSTOM supports) with `LoadSpec.temperature_change = 30`, "
          "combination T = 1.0, `run_analysis`.",
          "- **F2** `benchmarks.validation_model_a()` with a uniform "
          "temperature change through `run_analysis`; the one-spoke case uses "
          "`solve_static` with a per-member initial strain (the GUI load "
          "input heats every strut uniformly).",
          "- **F3** Model A with `LoadSpec.settlements` (vertical settlement "
          "of hand node 2, radial spread of hand node 3), `run_analysis`.",
          "- **F4** `benchmarks.three_bar_truss` (support C settles 5 mm): "
          "load case and settlement case through `run_analysis`, bar 2 heated "
          "alone and the combined case through `solve_static`. The "
          "flexibility coefficient f11 is obtained independently of the hand "
          "force method by a stiffness solve of the primary structure. The "
          "program's own force-method routine "
          "(`benchmarks.three_bar_force_method`, used by the built-in preset "
          "check) is compared too.",
          "- **F5** CUSTOM supports with a spring stiffness entry "
          "(5e7 N/m) in series and in parallel with the bar, and a heated "
          "bar on a spring.",
          "- **F6** `benchmarks.von_mises_truss()`: linear solve and "
          "`buckling.linear_buckling`; `nonlinear.run_nonlinear` (Crisfield "
          "arc length) traced through the limit point, the unstable branch "
          "and onto the inverted stable branch; every traced point is "
          "compared with the closed form P(w). The sample loads P(w) are the "
          "program's corotational internal force at the symmetric state; "
          "both limit points are also located as zeros of the program's "
          "consistent tangent; the snap-through target is the program's "
          "equilibrium at lambda* found by Newton iteration on its residual. "
          "`linear_buckling` returns n_free - 2 modes, one for this 2-DOF "
          "truss, so the sway root uses the same program matrices K and K_g.",
          "- **F7** `modal.modal_analysis` on Model A, lumped mass, with and "
          "without the cladding dead load as mass.",
          "- **F8** a full 2V Class I sphere built with the program's own "
          "`icosahedron`, `subdivide_class1` and `custom_geometry`, and the "
          "2V hemisphere from `generate_dome`, both through "
          "`fabrication.fabrication_tables`; chord factors also against "
          "Kenner (1976)."]
    fails = [r for r in rows if not r["pass"]]
    L += ["", "## Failures", ""]
    if fails:
        for r in fails:
            L.append(f"- {r['sub']} {r['quantity']} ({r['key']}): hand "
                     f"{_fmt(r['hand'])}, program {_fmt(r['program'])}, "
                     f"d = {r['difference']:.2e} > {r['tolerance']:.0e}")
    else:
        L.append("None.")
    L += ["", "## Full table", "",
          "| Sub-case | Quantity | Hand key | Hand | GeoDomeLab | Unit | d | "
          "Tol | Pass | Program source |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['sub']} | {r['quantity'].replace('|', chr(92) + '|')} | "
                 f"{r['key']} | {_fmt(r['hand'])} | {_fmt(r['program'])} | "
                 f"{r['unit']} | {r['difference']:.1e} | "
                 f"{r['tolerance']:.0e} | {'yes' if r['pass'] else 'NO'} | "
                 f"{r['source']} |")
    L += ["", f"Run time {info['elapsed']:.1f} s.", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    return path


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:                                  # pragma: no cover
        pass
    rows, info = run_comparison()
    write_csv(rows)
    write_markdown(rows, info)
    print("=" * 112)
    print(f"{'sub-case':<26}{'quantity':<42}{'hand':>16}{'program':>16}"
          f"{'d':>10}  ok")
    print("-" * 112)
    for r in rows:
        print(f"{r['sub'][:25]:<26}{r['quantity'][:41]:<42}"
              f"{_fmt(r['hand'], 9):>16}{_fmt(r['program'], 9):>16}"
              f"{r['difference']:>10.1e}  {'ok' if r['pass'] else 'FAIL'}")
    print("-" * 112)
    summ = summary_by_case(rows)
    for case in sorted(summ):
        s = summ[case]
        worst = ", ".join(f"{k} {v:.1e}" for k, v in sorted(s["worst"].items()))
        print(f"{case}  {CASE_TITLES[case]:<46} {s['pass']:>3}/{s['n']:<3} "
              f"worst d: {worst}")
    npass = sum(r["pass"] for r in rows)
    print(f"\n{npass} of {len(rows)} comparisons pass; {info['n_used']} of "
          f"{info['n_hand']} hand reference values compared"
          + (f"; NOT compared: {info['missing']}" if info["missing"] else "")
          + f".  geodome loaded before the hand values: "
          f"{info['geodome_before_hand']}.  {info['elapsed']:.1f} s.")
    print(f"written: {CSV_PATH}\n         {MD_PATH}")
    ok = npass == len(rows) and not info["missing"]
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
