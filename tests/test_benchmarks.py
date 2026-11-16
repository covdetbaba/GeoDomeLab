"""
Tests for the benchmark presets (spec task G, ``geodome/benchmarks.py``).

What is checked, and against what:

* Validation Models A and B: every row of ``validation/output/comparison.csv``
  is recomputed from the *preset* models and must equal both the hand answer
  and the program value stored there (round-off only); the presets extend
  the saved validation models only by the hand-check point-load case; the
  constants equal ``validation/hand_calc.py``.
* von Mises truss: closed-form linear answers and linearised lambda_cr
  written out in this file; the closed-form path helpers against the
  research-note example and against validation case F6.
* Star dome: the published coordinate/length table; the linear response
  against the independent symmetric virtual-work solution; lambda_cr from an
  eigen solve assembled in this file; and the published snap-through points
  (Greco & Venturini 2006, Iu & Bradford 2012) reproduced by a small
  corotational path-following solver written here, independent of the
  package.
* Three-bar truss: the force method against validation case F4
  (``validation/hand_calc_ext.py``, an independent implementation) and the
  program against the force method for load, temperature, settlement and
  their sum.
* The 6V case study equals ``examples/worked_example.build_model()``.
* ``examples/models/<key>.json`` are current, every preset is properly
  supported (positive definite free stiffness), runs the full pipeline and
  draws, and ``verify_preset`` passes -- and fails when the model is wrong.

Run with pytest::

    pytest -q tests/test_benchmarks.py

or standalone, without pytest installed::

    python tests/test_benchmarks.py
"""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VALID = os.path.join(ROOT, "validation")
sys.path.insert(0, ROOT)
if VALID not in sys.path:
    sys.path.insert(0, VALID)

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")
from matplotlib.figure import Figure                                # noqa: E402
import scipy.linalg as sla                                          # noqa: E402

from geodome import benchmarks as bm                                # noqa: E402
from geodome import geometry as geo                                 # noqa: E402
from geodome import plotting as plot                                # noqa: E402
from geodome.model import DomeModel                                 # noqa: E402
from geodome.postprocess import run_analysis, summary_text, \
    member_table, node_table, reaction_table                        # noqa: E402
from geodome.sections import SECTIONS                               # noqa: E402
from geodome.materials import MATERIALS                             # noqa: E402
from geodome.solver import assemble_stiffness, free_dofs, \
    solve_static                                                    # noqa: E402

import hand_calc as hc                                              # noqa: E402

COMPARISON_CSV = os.path.join(VALID, "output", "comparison.csv")


def _hand_ext():
    """validation/hand_calc_ext.py (independent F cases) or None."""
    try:
        import hand_calc_ext                                        # noqa: F401
        return hand_calc_ext
    except Exception as exc:                         # pragma: no cover
        print(f"    (hand_calc_ext not importable: {exc}; F-case "
              f"cross-checks skipped)")
        return None


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


# =============================================================== registry ====
def test_registry_is_complete_and_builds_fresh_models():
    expected = ["validation_a", "validation_b", "von_mises", "star_dome",
                "three_bar", "case_study_6v"]
    assert bm.preset_keys() == expected
    for key in expected:
        title, desc, factory = bm.PRESETS[key]
        assert isinstance(title, str) and title
        assert isinstance(desc, str) and len(desc) > 40
        assert callable(factory)
        m1, m2 = bm.make_preset(key), bm.make_preset(key)
        assert isinstance(m1, DomeModel)
        assert m1 is not m2 and m1.loads is not m2.loads
        m1.loads.nodal_loads.append({"node": 0, "Fz": -1.0})
        m1.combinations.pop()
        m3 = bm.make_preset(key)
        assert m3.to_dict() == m2.to_dict()          # factories share nothing
        assert bm.preset_info(key) == (title, desc)
    try:
        bm.make_preset("no_such_preset")
    except KeyError:
        pass
    else:
        raise AssertionError("unknown key accepted")


# ======================================================= validation A / B ====
def test_validation_constants_equal_hand_calc():
    import run_validation as rv
    assert bm.VAL_RADIUS == hc.R
    sec = SECTIONS[bm.VAL_SECTION]
    assert sec.D == hc.D_OUT and sec.t == hc.T_WALL
    assert bm.VAL_SECTION == rv.SECTION and bm.VAL_MATERIAL == rv.MATERIAL
    mat = MATERIALS[bm.VAL_MATERIAL]
    assert (mat.E, mat.Fy, mat.rho) == (hc.E, hc.FY, hc.RHO)
    assert (bm.VAL_K_EFF, bm.VAL_PHI_C, bm.VAL_PHI_T) == \
        (hc.K_EFF, hc.PHI_C, hc.PHI_T)
    assert bm.VAL_DEAD_PRESSURE == hc.Q_DEAD
    assert bm.VAL_SNOW_PRESSURE == hc.P_SNOW
    assert bm.VAL_A1_LOAD == (hc.PX_A1, hc.PY_A1, hc.PZ_A1)
    assert bm.VAL_B1_LOAD == -hc.PZ_B1
    assert bm.FRACTION_A == rv.FRACTION_A and bm.FRACTION_B == rv.FRACTION_B
    # the reference numbers carried by the module = the hand calculation
    A, B = hc.model_a(), hc.model_b()
    ref = {"A1_ux": A["A1"]["D"][0], "A1_uz": A["A1"]["D"][2],
           "A1_N1": A["A1"]["N"][1], "A1_sumRz": A["A1"]["Rsum"][2],
           "A2_uz": A["A2"]["D"][2], "A2_N_spoke": A["A2"]["N"][1],
           "A2_Rz_node": A["A2"]["R"][2][2], "A2_total_Pz": A["A2"]["tot_P"],
           "A2_util": A["aisc"]["util"],
           "A2_lambda_cr": A["buckling"]["lam_cr"]}
    for k, v in ref.items():
        assert _rel(bm.HAND_A[k], float(v)) < 1e-12, k
    for case in ("B1", "B2"):
        for k, v in bm.HAND_B[case].items():
            assert _rel(v, float(B[case][k])) < 1e-12, (case, k)


def _program_rows_a(model):
    """Every quantity of comparison.csv model A, computed from the preset."""
    import run_validation as rv
    H = hc.model_a()
    g = model.geom
    nmap = rv.node_map(H["nodes"], g.nodes)
    mmap = rv.member_map(H["members"], g.members, nmap)
    c = nmap[1]
    out = {}

    def put(case, q, v, s=1.0):
        out[("A", case, q)] = float(v) * s

    put("geometry", "nodes", g.n_nodes)
    put("geometry", "struts", g.n_members)
    put("geometry", "supported nodes", len(model.support_nodes()))
    put("geometry", "strut length L", np.max(g.strut_lengths))
    put("geometry", "chord factor L/R", g.strut_schedule()[0]["chord_factor"])
    put("geometry", "crown height", g.height)
    K = assemble_stiffness(g.nodes, g.members, model.EA())
    K = K.toarray() if hasattr(K, "toarray") else K
    Kc = K[3 * c:3 * c + 3, 3 * c:3 * c + 3]
    for i, lab in enumerate(("Kff[1,1]", "Kff[2,2]", "Kff[3,3]")):
        put("stiffness", lab, Kc[i, i], 1e-3)
    put("stiffness", "largest off-diagonal of Kff",
        np.max(np.abs(Kc - np.diag(np.diag(Kc)))), 1e-3)

    res = run_analysis(model, n_buckling_modes=6)
    r1 = res.cases[bm.COMBO_A1]
    put("A1", "crown ux (D1)", r1.u[3 * c], 1e3)
    put("A1", "crown uy (D2)", r1.u[3 * c + 1], 1e3)
    put("A1", "crown uz (D3)", r1.u[3 * c + 2], 1e3)
    for m in range(1, 6):
        put("A1", f"N member {m}", r1.N[mmap[m]], 1e-3)
    put("A1", "N ring members 6-10 (max |N|)",
        max(abs(r1.N[mmap[m]]) for m in range(6, 11)), 1e-3)
    Rp = r1.reaction_vectors()
    for nd in (2, 4):
        put("A1", f"Rx at node {nd}", Rp[nmap[nd], 0], 1e-3)
        put("A1", f"Rz at node {nd}", Rp[nmap[nd], 2], 1e-3)
    put("A1", "sum Rx", Rp.sum(axis=0)[0], 1e-3)
    put("A1", "sum Rz", Rp.sum(axis=0)[2], 1e-3)

    r2 = res.cases[bm.COMBO_A2]
    put("A2", "crown load Pz", r2.F[3 * c + 2], 1e-3)
    put("A2", "ring node load Pz", r2.F[3 * nmap[2] + 2], 1e-3)
    put("A2", "total vertical load", r2.F[2::3].sum(), 1e-3)
    put("A2", "crown uz", r2.u[3 * c + 2], 1e3)
    put("A2", "N spoke (all 5)",
        np.mean([r2.N[mmap[m]] for m in range(1, 6)]), 1e-3)
    put("A2", "Rz per support node", r2.reaction_vectors()[nmap[2], 2], 1e-3)
    k1 = mmap[1]
    put("AISC E3", "KL/r", res.caps["slenderness"][k1])
    put("AISC E3", "Fcr", res.caps["Fcr"][k1], 1e-6)
    put("AISC E3", "phi_c Pn", res.caps["phiPn_c"][k1], 1e-3)
    put("AISC D2", "phi_t Pn", res.caps["phiPn_t"][k1], 1e-3)
    put("AISC E3", "utilisation", res.util_env[k1])
    assert res.buckling.reference_case == bm.COMBO_A2
    put("buckling", "lambda_cr", res.buckling.load_factors[0])
    return out


def _program_rows_b(model):
    """Every quantity of comparison.csv model B, computed from the preset."""
    import run_validation as rv
    H = hc.model_b()
    g = model.geom
    nmap = rv.node_map(H["nodes"], g.nodes)
    mmap = rv.member_map(H["members"], g.members, nmap)
    c = nmap[1]
    out = {}

    def put(case, q, v, s=1.0):
        out[("B", case, q)] = float(v) * s

    put("geometry", "nodes", g.n_nodes)
    put("geometry", "struts", g.n_members)
    put("geometry", "supported nodes", len(model.support_nodes()))
    put("geometry", "crown height", g.height)
    K = assemble_stiffness(g.nodes, g.members, model.EA())
    K = K.toarray() if hasattr(K, "toarray") else K
    T = np.zeros((3 * g.n_nodes, 3))
    T[3 * c + 2, 0] = 1.0
    for i in range(5):
        nd = nmap[2 + i]
        th = 72.0 * i * hc.DEG
        T[3 * nd + 0, 1] = math.cos(th)
        T[3 * nd + 1, 1] = math.sin(th)
        T[3 * nd + 2, 2] = 1.0
    Kp = T.T @ K @ T
    for (i, j, lab) in ((0, 0, "Ks[1,1]"), (1, 1, "Ks[2,2]"),
                        (2, 2, "Ks[3,3]"), (0, 1, "Ks[1,2]"),
                        (1, 2, "Ks[2,3]"), (0, 2, "Ks[1,3]")):
        put("stiffness", lab + " (T^T K T)", Kp[i, j], 1e-3)

    res = run_analysis(model, do_buckling=False)
    for name, combo in (("B1", bm.COMBO_B1), ("B2", bm.COMBO_B2)):
        r = res.cases[combo]
        if name == "B2":
            put(name, "crown load Pz", r.F[3 * c + 2], 1e-3)
            put(name, "upper node load Pz", r.F[3 * nmap[2] + 2], 1e-3)
            put(name, "lower node load Pz", r.F[3 * nmap[7] + 2], 1e-3)
        put(name, "total vertical load", r.F[2::3].sum(), 1e-3)
        put(name, "crown uz (w0)", r.u[3 * c + 2], 1e3)
        urs, uts, w1s = [], [], []
        for i in range(5):
            nd = nmap[2 + i]
            th = 72.0 * i * hc.DEG
            ux, uy, uz = r.u[3 * nd:3 * nd + 3]
            urs.append(ux * math.cos(th) + uy * math.sin(th))
            uts.append(-ux * math.sin(th) + uy * math.cos(th))
            w1s.append(uz)
        put(name, "upper node radial (u_r)", np.mean(urs), 1e3)
        put(name, "upper node vertical (w1)", np.mean(w1s), 1e3)
        put(name, "upper node tangential (max |u_t|)",
            max(abs(x) for x in uts), 1e3)
        for lab, ms in (("N spokes 1-5", range(1, 6)),
                        ("N upper ring 6-10 (hoop)", range(6, 11)),
                        ("N zigzag 11-20", range(11, 21))):
            vals = [r.N[mmap[m]] for m in ms]
            put(name, lab, np.mean(vals), 1e-3)
            put(name, lab + " spread (max-min)", max(vals) - min(vals), 1e-3)
        put(name, "N lower ring 21-25 (max |N|)",
            max(abs(r.N[mmap[m]]) for m in range(21, 26)), 1e-3)
        Rv = r.reaction_vectors()
        nd = nmap[7]
        beta = 36.0 * hc.DEG
        put(name, "Rz per support node", Rv[nd, 2], 1e-3)
        put(name, "R radial per support node",
            Rv[nd, 0] * math.cos(beta) + Rv[nd, 1] * math.sin(beta), 1e-3)
        put(name, "sum Rz", Rv.sum(axis=0)[2], 1e-3)
    return out


def _compare_with_csv(model_letter, rows):
    with open(COMPARISON_CSV, newline="", encoding="utf-8") as fh:
        csv_rows = [r for r in csv.DictReader(fh)
                    if r["model"] == model_letter]
    assert len(csv_rows) >= 30, "comparison.csv looks incomplete"
    keys = {(r["model"], r["case"], r["quantity"]) for r in csv_rows}
    missing = keys - set(rows)
    assert not missing, f"quantities not reproduced: {sorted(missing)}"
    worst = 0.0
    for r in csv_rows:
        k = (r["model"], r["case"], r["quantity"])
        p = rows[k]
        hand = float(r["hand"])
        prog = float(r["program"])
        if abs(hand) > 1e-9:
            d_hand = _rel(p, hand)
            d_prog = _rel(p, prog)
            assert d_hand < 1e-12, f"{k}: preset {p!r} vs hand {hand!r}"
            assert d_prog < 1e-12, f"{k}: preset {p!r} vs csv {prog!r}"
            worst = max(worst, d_hand)
        else:
            assert abs(p) < 1e-9, f"{k}: preset {p!r} should be 0"
            assert abs(p - prog) < 1e-9, f"{k}: {p!r} vs csv {prog!r}"
    return len(csv_rows), worst


def test_run_validation_has_its_own_pass_rule():
    """Review V05: run_validation.py carries the tolerance, a pass column
    and an exit code; the GUI and run_all_checks.py use the same limits."""
    import run_validation as rv
    import importlib.util
    from geodome.gui import validation_view as vv
    assert rv.REL_TOL == vv.RV_REL_TOL == 1e-12
    assert rv.ABS_TOL == vv.RV_ABS_TOL == 1e-9
    spec = importlib.util.spec_from_file_location(
        "rac", os.path.join(os.path.dirname(VALID), "run_all_checks.py"))
    rac = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rac)
    assert abs(rac.RV_REL_TOL_PERCENT - 100.0 * rv.REL_TOL) < 1e-24
    assert rac.RV_ABS_TOL == rv.ABS_TOL
    saved = list(rv.ROWS)
    try:
        rv.ROWS.clear()
        rv.row("X", "t", "1 % off", 1.0, 1.01, "-")
        rv.row("X", "t", "exact", 2.0, 2.0, "-")
        rv.row("X", "t", "zero, small", 0.0, 5e-10, "kN")
        rv.row("X", "t", "zero, large", 0.0, 2e-9, "kN")
        assert [r["pass"] for r in rv.ROWS] == [False, True, True, False]
    finally:
        rv.ROWS[:] = saved
    with open(os.path.join(VALID, "run_validation.py"),
              encoding="utf-8") as fh:
        assert "raise SystemExit(main())" in fh.read()
    with open(COMPARISON_CSV, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows and all(r["pass"] == "True" for r in rows)


def test_validation_a_reproduces_every_row_of_comparison_csv():
    n, worst = _compare_with_csv("A", _program_rows_a(bm.validation_model_a()))
    assert n >= 37 and worst < 1e-12


def test_validation_b_reproduces_every_row_of_comparison_csv():
    n, worst = _compare_with_csv("B", _program_rows_b(bm.validation_model_b()))
    assert n >= 40 and worst < 1e-12


def test_validation_presets_extend_the_saved_validation_models():
    """Same model as validation/output/*.json + only the point-load case."""
    out = os.path.join(VALID, "output")

    def saved(name):
        with open(os.path.join(out, name), encoding="utf-8") as fh:
            return json.load(fh)

    def jd(model):
        return json.loads(json.dumps(model.to_dict()))

    pa = jd(bm.validation_model_a())
    sa = saved("model_A_LC2.json")
    for k in ("geometry", "material", "group_sections", "support_type",
              "design", "spring_stiffness", "custom_supports"):
        assert pa[k] == sa[k], k
    la, lsa = dict(pa["loads"]), dict(sa["loads"])
    assert len(la.pop("nodal_loads")) == 1 and lsa.pop("nodal_loads") == []
    assert la == lsa
    assert pa["combinations"][1] == sa["combinations"][0]        # A2
    a1 = pa["combinations"][0]
    assert (a1["D"], a1["S"], a1["L"], a1["W"], a1["T"]) == (0, 0, 1, 0, 0)

    pb = jd(bm.validation_model_b())
    for fname, idx in (("model_B_LC1.json", 0), ("model_B_LC2.json", 1)):
        sb = saved(fname)
        for k in ("geometry", "material", "group_sections", "support_type",
                  "design"):
            assert pb[k] == sb[k], (fname, k)
        cp, cs = dict(pb["combinations"][idx]), dict(sb["combinations"][0])
        assert cp.pop("serviceability") is (idx == 0)
        cs.pop("serviceability")
        assert cp == cs, fname
    lb2 = dict(saved("model_B_LC2.json")["loads"])
    lpb = dict(pb["loads"])
    assert lpb.pop("crown_point_load") == -hc.PZ_B1
    assert lb2.pop("crown_point_load") == 0.0
    assert lpb == lb2
    assert saved("model_B_LC1.json")["loads"]["crown_point_load"] == -hc.PZ_B1


# ================================================================ supports ===
def test_every_preset_is_properly_supported():
    """Free stiffness (with springs) positive definite; planar ones stay planar."""
    for key in bm.preset_keys():
        m = bm.make_preset(key)
        g = m.geom
        fixed, springs = m.support_conditions()
        K = assemble_stiffness(g.nodes, g.members, m.EA())
        K = K.toarray() if hasattr(K, "toarray") else np.asarray(K)
        for d, k in springs.items():
            K[d, d] += k
        fd = free_dofs(3 * g.n_nodes, fixed)
        ev = np.linalg.eigvalsh(K[np.ix_(fd, fd)])
        assert ev[0] > 1e-10 * ev[-1], (key, ev[0], ev[-1])
    for key in ("von_mises", "three_bar"):
        m = bm.make_preset(key)
        assert m.support_type == geo.CUSTOM
        assert np.all(m.geom.nodes[:, 1] == 0.0)
        fixed, springs = m.support_conditions()
        assert not springs
        # every uy is restrained, in-plane DOF of the loaded joint are free
        assert set(range(1, 3 * m.geom.n_nodes, 3)) <= set(fixed.tolist())
        top = m.geom.crown_node
        assert 3 * top not in fixed and 3 * top + 2 not in fixed
        res = run_analysis(m, do_buckling=False)
        for r in res.cases.values():
            assert np.all(r.u[1::3] == 0.0)
            assert np.max(np.abs(r.reactions[1::3])) < 1e-9


# =============================================================== von Mises ===
def test_von_mises_linear_closed_form():
    m = bm.von_mises_truss()
    X = m.geom.nodes
    assert np.array_equal(X, [[-2.5, 0, 0], [2.5, 0, 0], [0, 0, 0.125]])
    assert m.geom.members.tolist() == [[0, 2], [1, 2]]
    a, h, P = 2.5, 0.125, 10e3
    EA = 200e9 * SECTIONS["CHS 88.9 x 4.0"].A
    L = math.hypot(a, h)
    res = run_analysis(m, n_buckling_modes=2)
    r = res.cases[res.order[0]]
    assert _rel(r.u[8], -P * L ** 3 / (2 * EA * h * h)) < 1e-12
    assert abs(r.u[6]) < 1e-18
    assert np.allclose(r.N, -P * L / (2 * h), rtol=1e-12, atol=0)
    R = r.reaction_vectors()
    assert _rel(R[0, 2], P / 2) < 1e-12 and _rel(R[1, 2], P / 2) < 1e-12
    assert _rel(R[0, 0], P * a / (2 * h)) < 1e-12
    assert _rel(R[1, 0], -P * a / (2 * h)) < 1e-12
    lam = res.buckling.load_factors
    assert _rel(lam[0], 2 * EA * h ** 3 / (L * a * a * P)) < 1e-10
    # the critical mode is the symmetric snap (apex vertical), not the sway
    mode = res.buckling.mode(0)[2]
    assert abs(mode[0]) < 1e-12 * abs(mode[2])
    assert bm.von_mises_parameters(m) == (a, h, EA, P)


def test_von_mises_closed_form_helpers():
    # example: b = 1, h = 0.1, EA = 1
    lp = bm.von_mises_limit_point(1.0, 0.1, 1.0)
    assert abs(lp["w"] - 0.042360746517) < 1e-11
    assert _rel(lp["P"], 3.810871904181e-4) < 1e-11
    assert _rel(bm.von_mises_linear(1.0, 0.1, 1.0, 1.0)["lambda_cr"],
                1.990074e-3) < 1e-6
    a, h, EA = 2.5, 0.125, 2.0e8
    lp = bm.von_mises_limit_point(a, h, EA)
    w = np.linspace(0.0, 2.5 * h, 200001)
    P = bm.von_mises_load(w, a, h, EA)
    k = int(np.argmax(P[w < h]))
    assert _rel(P[k], lp["P"]) < 1e-9 and abs(w[k] - lp["w"]) < 2e-5
    d = 1e-6
    slope = (bm.von_mises_load(lp["w"] + d, a, h, EA)
             - bm.von_mises_load(lp["w"] - d, a, h, EA)) / (2 * d)
    assert abs(slope) < 1e-7 * bm.von_mises_linear(a, h, EA, 1.0)["K0"]
    assert abs(bm.von_mises_load(h, a, h, EA)) < 1e-12 * lp["P"]
    assert abs(bm.von_mises_load(2 * h, a, h, EA)) < 1e-12 * lp["P"]
    s = bm.von_mises_load(1e-7, a, h, EA) / 1e-7
    assert _rel(s, bm.von_mises_linear(a, h, EA, 1.0)["K0"]) < 1e-5
    # mid-path mirror symmetry P(h + y) = -P(h - y)
    assert _rel(bm.von_mises_load(h + 0.03, a, h, EA),
                -bm.von_mises_load(h - 0.03, a, h, EA)) < 1e-12


def test_von_mises_preset_is_validation_case_F6():
    ext = _hand_ext()
    if ext is None:                                  # pragma: no cover
        return
    d = ext.definitions()["F6"]
    m = bm.von_mises_truss()
    assert np.array_equal(m.geom.nodes, np.asarray(d["nodes"], float))
    assert m.geom.members.tolist() == d["members"]
    assert [m.group_sections[g] for g in m.geom.strut_group] == d["sections"]
    assert m.custom_supports == {int(k): list(v)
                                 for k, v in d["supports"].items()}
    nl = m.loads.nodal_loads[0]
    assert (nl["node"], nl["Fx"], nl["Fy"], nl["Fz"]) == \
        (d["load"]["node"], d["load"]["Fx"], d["load"]["Fy"], d["load"]["Fz"])
    res = ext.results()
    a, h, EA, P = bm.von_mises_parameters(m)
    assert _rel(EA, res["F6.EA"]["value"]) < 1e-12
    lp = bm.von_mises_limit_point(a, h, EA)
    assert _rel(lp["w"], res["F6.w_limit"]["value"]) < 1e-10
    assert _rel(lp["P"], res["F6.P_limit"]["value"]) < 1e-10
    for r in (0.1, 0.25, 0.5, 1.5, 2.5):
        ref = res[f"F6.P_at_w{r:g}h"]["value"]
        assert _rel(float(bm.von_mises_load(r * h, a, h, EA)), ref) < 1e-12
    rr = run_analysis(m, n_buckling_modes=1)
    assert _rel(rr.buckling.load_factors[0], res["F6.lambda_cr"]["value"]) \
        < 1e-10
    w_lin = -rr.cases[rr.order[0]].u[8]
    assert _rel(P / w_lin, res["F6.K0"]["value"]) < 1e-12


# =============================================================== star dome ===
STAR_TABLE_CM = [   # docs/VALIDATION.md, star dome, nodes 1..13
    (0.00000, 0.00000, 8.216), (25.00000, 0.00000, 6.216),
    (12.50000, 21.65064, 6.216), (-12.50000, 21.65064, 6.216),
    (-25.00000, 0.00000, 6.216), (-12.50000, -21.65064, 6.216),
    (12.50000, -21.65064, 6.216), (43.30127, 25.00000, 0.0),
    (0.00000, 50.00000, 0.0), (-43.30127, 25.00000, 0.0),
    (-43.30127, -25.00000, 0.0), (0.00000, -50.00000, 0.0),
    (43.30127, -25.00000, 0.0)]


def test_star_dome_geometry_matches_published_table():
    m = bm.star_dome()
    g = m.geom
    assert (g.n_nodes, g.n_members, g.n_faces) == (13, 24, 12)
    X_cm = g.nodes * 100.0
    assert np.max(np.abs(X_cm - np.asarray(STAR_TABLE_CM))) < 6e-6
    L_cm = g.member_lengths() * 100.0
    assert np.allclose(L_cm[0:6], 25.07987, atol=6e-6)
    assert np.allclose(L_cm[6:12], 25.00000, atol=1e-9)
    assert np.allclose(L_cm[12:24], 31.60024, atol=6e-6)
    assert abs((g.nodes[0, 2] - g.nodes[1, 2]) * 100 - 2.0) < 1e-12
    assert abs(g.height * 100 - 8.216) < 1e-12
    assert g.crown_node == 0
    assert m.support_type == geo.PINNED
    assert m.support_nodes().tolist() == list(range(7, 13))
    assert len(np.setdiff1d(np.arange(39), m.fixed_dofs())) == 21
    assert np.bincount(g.strut_group).tolist() == [6, 6, 12]
    # every member of the papers' list, nothing else
    one_based = {frozenset((int(i) + 1, int(j) + 1)) for i, j in g.members}
    assert one_based == {frozenset(p) for p in
                         [(1, k) for k in range(2, 8)]
                         + [(k, 2 + (k - 1) % 6) for k in range(2, 8)]
                         + [(2, 8), (3, 8), (3, 9), (4, 9), (4, 10), (5, 10),
                            (5, 11), (6, 11), (6, 12), (7, 12), (7, 13),
                            (2, 13)]}
    # printed 43.3 variant
    X3 = bm.star_dome(exact_hexagon=False).geom.nodes * 100
    assert abs(X3[7, 0] - 43.3) < 1e-12
    # load: P_ref = 1e-4 EA down at the crown, nothing else
    EA = float(m.EA()[0])
    assert np.allclose(m.EA(), EA)
    assert m.loads.nodal_loads == [{"node": 0, "Fx": 0.0, "Fy": 0.0,
                                    "Fz": -1e-4 * EA, "case": "L"}]
    assert not m.loads.include_self_weight and not m.loads.wind_enabled


def _assemble_linear(X, M, EA):
    """Own assembly: K = sum EA/L [cc' -cc'; -cc' cc'] (no package code)."""
    n = len(X)
    K = np.zeros((3 * n, 3 * n))
    for k, (i, j) in enumerate(M):
        d = X[j] - X[i]
        L = math.sqrt(float(d @ d))
        c = d / L
        kk = EA[k] / L * np.outer(c, c)
        for a, b, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
            K[3 * a:3 * a + 3, 3 * b:3 * b + 3] += s * kk
    return K


def test_star_dome_linear_response_and_buckling():
    m = bm.star_dome()
    g = m.geom
    EA = m.EA()
    res = run_analysis(m, n_buckling_modes=6)
    r = res.cases[res.order[0]]
    sym = bm.star_dome_symmetric_linear(m)
    assert _rel(r.u[2], sym["w0"]) < 1e-12
    assert np.allclose(r.u[5::3][:6], sym["w1"], rtol=1e-11)
    assert np.allclose(r.N[0:6], sym["N_spoke"], rtol=1e-11)
    assert np.allclose(r.N[6:12], sym["N_ring"], rtol=1e-11)
    assert np.allclose(r.N[12:24], sym["N_diag"], rtol=1e-11)
    P = 1e-4 * EA[0]
    assert abs(P / abs(r.u[2]) / EA[0] - bm.STAR_K0_PER_EA) < 5e-6
    # independent eigenproblem: own K and K_g from own linear solution
    fd = np.setdiff1d(np.arange(39), np.arange(21, 39))
    K = _assemble_linear(g.nodes, g.members, EA)
    F = np.zeros(39)
    F[2] = -P
    u = np.zeros(39)
    u[fd] = np.linalg.solve(K[np.ix_(fd, fd)], F[fd])
    Kg = np.zeros((39, 39))
    for k, (i, j) in enumerate(g.members):
        d = g.nodes[j] - g.nodes[i]
        L = math.sqrt(float(d @ d))
        c = d / L
        N = EA[k] / L * float(c @ (u[3 * j:3 * j + 3] - u[3 * i:3 * i + 3]))
        kg = N / L * (np.eye(3) - np.outer(c, c))
        for a, b, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
            Kg[3 * a:3 * a + 3, 3 * b:3 * b + 3] += s * kg
    mu = sla.eigh(-Kg[np.ix_(fd, fd)], K[np.ix_(fd, fd)], eigvals_only=True)
    lam_own = 1.0 / mu.max()
    lam = res.buckling.load_factors[0]
    assert _rel(lam, lam_own) < 1e-10
    assert abs(lam * 1e-4 - bm.STAR_LAMBDA_CR_P_OVER_EA) < 5e-9
    assert res.buckling.reference_case == res.order[0]


def _corot(X0, M, EA, u):
    """Corotational bar, engineering strain (own code): f, K_T, N."""
    n = len(X0)
    x = X0 + u.reshape(-1, 3)
    f = np.zeros(3 * n)
    K = np.zeros((3 * n, 3 * n))
    N = np.zeros(len(M))
    I3 = np.eye(3)
    for k, (i, j) in enumerate(M):
        D0 = X0[j] - X0[i]
        L0 = math.sqrt(float(D0 @ D0))
        d = x[j] - x[i]
        l = math.sqrt(float(d @ d))
        N[k] = EA * (l * l - L0 * L0) / (L0 * (l + L0))   # EA (l - L0)/L0
        c = d / l
        f[3 * i:3 * i + 3] -= N[k] * c
        f[3 * j:3 * j + 3] += N[k] * c
        kk = EA / L0 * np.outer(c, c) + N[k] / l * (I3 - np.outer(c, c))
        for a, b, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
            K[3 * a:3 * a + 3, 3 * b:3 * b + 3] += s * kk
    return f, K, N


def test_star_dome_reproduces_published_snap_through_points():
    """Crown displacement control on the preset geometry (cm, EA = 1).

    Greco & Venturini (2006): first limit point "near 0.75 cm".
    Iu & Bradford (2012): zero-load state with the crown at ~19 mm and
    node 2 at ~-1.16 mm (i.e. moving up); docs/VALIDATION.md: P/EA =
    3.15655e-4 at w = 0.76844 cm, zero load at 1.8838 cm with node 2 up
    0.1162 cm, minimum -2.76000e-4 at 3.0278 cm, exact mirror state
    (P = 0, ring back at z) at w = 4.000 cm.
    """
    m = bm.star_dome()
    X0 = m.geom.nodes * 100.0                       # cm
    M = m.geom.members
    cz = 2
    free = [d for d in range(21) if d != cz]        # nodes 0..6 minus crown z

    def solve_at(w, u0):
        # Newton on the 20 unloaded DOF with the crown uz prescribed.  With
        # EA = 1 the forces are O(1e-4) and round-off leaves |r| ~ 5e-16,
        # so stop on a tiny correction and require |r| < 1e-13 (1e-9 of P).
        u = u0.copy()
        u[cz] = -w
        for _ in range(40):
            f, K, _ = _corot(X0, M, 1.0, u)
            r = f[free]
            du = np.linalg.solve(K[np.ix_(free, free)], -r)
            u[free] += du
            if np.linalg.norm(du) < 1e-13:
                break
        f, _, N = _corot(X0, M, 1.0, u)
        assert np.linalg.norm(f[free]) < 1e-13, "Newton did not converge"
        return -f[cz], u, N                         # P/EA (down), state

    ws = np.round(np.arange(0.0, 4.0001, 0.01), 10)
    Ps, Us = [], []
    u = np.zeros(39)
    for w in ws:
        P, u, _ = solve_at(w, u)
        Ps.append(P)
        Us.append(u.copy())
    Ps = np.array(Ps)

    def golden(k, sign):
        a_, b_ = ws[k - 1], ws[k + 1]
        g = (math.sqrt(5) - 1) / 2
        for _ in range(60):
            c_ = b_ - g * (b_ - a_)
            d_ = a_ + g * (b_ - a_)
            if sign * solve_at(c_, Us[k])[0] > sign * solve_at(d_, Us[k])[0]:
                b_ = d_
            else:
                a_ = c_
        w_ = 0.5 * (a_ + b_)
        return w_, solve_at(w_, Us[k])

    def zero(k):                                    # sign change in [k, k+1]
        a_, b_ = ws[k], ws[k + 1]
        Pa = Ps[k]
        for _ in range(60):
            c_ = 0.5 * (a_ + b_)
            Pc = solve_at(c_, Us[k])[0]
            if (Pc > 0) == (Pa > 0):
                a_, Pa = c_, Pc
            else:
                b_ = c_
        w_ = 0.5 * (a_ + b_)
        return w_, solve_at(w_, Us[k])

    k1 = int(np.argmax(Ps[ws < 1.5]))
    w1, (P1, u1, _) = golden(k1, +1)
    assert _rel(P1, bm.STAR_P_OVER_EA_LIMIT) < 2e-5, P1
    assert abs(w1 - bm.STAR_W_LIMIT_CM) < 1e-4, w1
    assert abs(w1 - 0.76844) < 2e-4                 # reference value
    assert abs(w1 - 0.75) < 0.05                    # Greco & Venturini
    ring_uz = u1[5:21:3]
    assert np.ptp(ring_uz) < 1e-12                  # symmetric path

    kz = int(np.where((Ps[:-1] > 0) & (Ps[1:] <= 0))[0][0])
    wz, (Pz, uz_, _) = zero(kz)
    assert abs(wz - bm.STAR_W_ZERO1_CM) < 2e-4, wz
    assert abs(uz_[5] - bm.STAR_NODE2_ZERO1_CM) < 2e-4, uz_[5]
    assert abs(wz - 1.9) < 0.05                     # Iu & Bradford, ~19 mm
    assert abs(uz_[5] - 0.116) < 0.005              # ... and ~1.16 mm up

    k2 = int(np.argmin(Ps[(ws > 2.0) & (ws < 3.5)])) + int(np.sum(ws <= 2.0))
    w2, (P2, _, _) = golden(k2, -1)
    assert _rel(P2, -2.76000e-4) < 2e-5 and abs(w2 - 3.0278) < 2e-4

    # the traced path passes the exact mirror state at w = 2 x 2.000 cm
    i4 = int(np.argmin(np.abs(ws - 4.0)))
    assert abs(Ps[i4]) < 1e-12
    assert np.max(np.abs(Us[i4][5:21:3])) < 1e-9
    # closed-form mirror states: every bar at its original length, N = 0
    for uz in ({0: -4.0}, {0: -2 * 8.216, **{n: -2 * 6.216
                                            for n in range(1, 7)}}):
        u = np.zeros(39)
        for n, v in uz.items():
            u[3 * n + 2] = v
        f, _, N = _corot(X0, M, 1.0, u)
        assert np.max(np.abs(N)) < 1e-13 and np.max(np.abs(f)) < 1e-13


# =============================================================== three-bar ===
def test_three_bar_preset_is_validation_case_F4():
    ext = _hand_ext()
    if ext is None:                                  # pragma: no cover
        return
    d = ext.definitions()["F4"]
    m = bm.three_bar_truss()
    assert np.array_equal(m.geom.nodes, np.asarray(d["nodes"], float))
    assert m.geom.members.tolist() == d["members"]
    assert [m.group_sections[g] for g in m.geom.strut_group] == d["sections"]
    assert np.allclose(m.member_areas(), d["A"], rtol=1e-15, atol=0)
    assert m.material.E == d["E"] and m.material.alpha == d["alpha"]
    assert m.custom_supports == {int(k): list(v)
                                 for k, v in d["supports"].items()}
    nl = m.loads.nodal_loads[0]
    assert (nl["node"], nl["Fx"], nl["Fy"], nl["Fz"]) == \
        (d["load"]["node"], d["load"]["Fx"], d["load"]["Fy"], d["load"]["Fz"])
    s = m.loads.settlements[0]
    assert (s["node"], s["ux"], s["uy"], s["uz"]) == \
        tuple(d["settlement"][k] for k in ("node", "ux", "uy", "uz"))
    assert m.loads.temperature_change == d["temperature"]["dT"]

    # the module's force method = F4's independent force method
    res = ext.results()
    dT1 = np.zeros(3)
    dT1[d["temperature"]["member"]] = d["temperature"]["dT"]
    sett = {d["settlement"]["node"]: (d["settlement"]["ux"],
                                      d["settlement"]["uz"])}
    P = (d["load"]["Fx"], d["load"]["Fz"])
    cases = {"load": (P, np.zeros(3), None),
             "temperature": ((0, 0), dT1, None),
             "settlement": ((0, 0), np.zeros(3), sett),
             "combined": (P, dT1, sett)}
    names = d["node_names"]
    for cn, (Pc, dTc, sc) in cases.items():
        fm = bm.three_bar_force_method(m, Pc, dTc, sc)
        for i in range(3):
            ref = res[f"F4.N{i + 1}_{cn}"]["value"]
            assert abs(fm["N"][i] - ref) <= 1e-9 * max(abs(ref), 1.0), (cn, i)
        assert abs(fm["u"][0] - res[f"F4.ux_{cn}"]["value"]) < 1e-15
        assert abs(fm["u"][1] - res[f"F4.uz_{cn}"]["value"]) < 1e-15
        for s_node, (rx, rz) in fm["R"].items():
            nm = names[s_node]
            assert abs(rx - res[f"F4.R{nm}x_{cn}"]["value"]) < 1e-8
            assert abs(rz - res[f"F4.R{nm}z_{cn}"]["value"]) < 1e-8
        # ... and the program reproduces F4 itself (eps0 per bar via solver)
        F = np.zeros(12)
        F[9], F[11] = Pc
        fixed, _ = m.support_conditions()
        presc = ({3 * k + 2: v[1] for k, v in sc.items()} if sc else None)
        r = solve_static(m.geom.nodes, m.geom.members, m.EA(), F, fixed,
                         prescribed=presc,
                         eps0=m.material.alpha * dTc if dTc.any() else None)
        for i in range(3):
            ref = res[f"F4.N{i + 1}_{cn}"]["value"]
            assert abs(r.N[i] - ref) <= 1e-9 * max(abs(ref), 1.0), (cn, i)
        assert abs(r.u[9] - res[f"F4.ux_{cn}"]["value"]) < 1e-15
        assert abs(r.u[11] - res[f"F4.uz_{cn}"]["value"]) < 1e-15


def test_three_bar_program_equals_force_method():
    m = bm.three_bar_truss()
    res = run_analysis(m, n_buckling_modes=2)
    alpha = m.material.alpha
    P = (15e3, -60e3)
    sett = {2: (0.0, -0.005)}
    dT = np.full(3, 30.0)
    inputs = {"P  1.0L": (P, np.zeros(3), None),
              "T  1.0T  temperature + settlement": ((0, 0), dT, sett),
              "P + T  1.0L + 1.0T": (P, dT, sett)}
    expect = {}
    X = m.geom.nodes
    for name, (Pc, dTc, sc) in inputs.items():
        fm = bm.three_bar_force_method(m, Pc, dTc, sc)
        expect[name] = fm
        r = res.cases[name]
        assert np.allclose(r.N, fm["N"], rtol=1e-11, atol=1e-6), name
        assert abs(r.u[9] - fm["u"][0]) < 1e-15 and \
            abs(r.u[11] - fm["u"][1]) < 1e-15, name
        R = r.reaction_vectors()
        for s, (rx, rz) in fm["R"].items():
            assert abs(R[s, 0] - rx) < 1e-7 and abs(R[s, 2] - rz) < 1e-7
        # compatibility of every bar: the elongation from the program's
        # displacements (incl. the settled support) = N L/EA + alpha dT L
        U = r.u.reshape(-1, 3)
        for k, (i, j) in enumerate(m.geom.members):
            c = (X[j] - X[i]) / r.lengths[k]
            e_prog = float(c @ (U[j] - U[i]))
            e_force = (fm["N"][k] * r.lengths[k] / m.EA()[k]
                       + alpha * dTc[k] * r.lengths[k])
            assert abs(e_prog - e_force) < 1e-15, (name, k)
        if sc:
            assert U[2, 2] == -0.005                 # support C settled
    # superposition: (P) + (T) = (P + T), member by member
    NP_, NT, NPT = (res.cases[n].N for n in expect)
    assert np.allclose(NP_ + NT, NPT, rtol=1e-12, atol=1e-7)
    # uniform heating = the sum of heating each bar alone (linearity)
    parts = sum(bm.three_bar_force_method(m, (0, 0), np.eye(3)[i] * 30.0)["N"]
                for i in range(3))
    whole = bm.three_bar_force_method(m, (0, 0), dT)["N"]
    assert np.allclose(parts, whole, rtol=1e-12, atol=1e-7)
    # a 30 degC heating of a stiff indeterminate truss is not negligible
    assert np.max(np.abs(whole)) > 0.2 * np.max(np.abs(NP_))
    # a determinate check: no settlement/temperature -> T case is zero
    m0 = bm.three_bar_truss(temperature_change=0.0, settlement_uz=0.0)
    r0 = run_analysis(m0, do_buckling=False)
    assert np.max(np.abs(r0.cases["T  1.0T  temperature + settlement"].N)) \
        == 0.0


# ============================================================== case study ===
def test_case_study_equals_worked_example():
    path = os.path.join(ROOT, "examples", "worked_example.py")
    spec = importlib.util.spec_from_file_location("_worked_example", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert bm.case_study_6v().to_dict() == mod.build_model().to_dict()


# ============================================================ JSON files ====
def test_example_json_files_are_current_and_reload_identically():
    for key in bm.preset_keys():
        path = bm.preset_path(key)
        assert os.path.isfile(path), f"missing {path}; run " \
            f"python -m geodome.benchmarks --write"
        with open(path, encoding="utf-8") as fh:
            on_disk = json.load(fh)
        fresh = json.loads(json.dumps(bm.make_preset(key).to_dict()))
        assert on_disk == fresh, f"{path} is stale; rerun --write"
        if key == "case_study_6v":
            continue                                # covered by verify_preset
        a = run_analysis(bm.make_preset(key), do_buckling=False)
        b = run_analysis(bm.load_preset_file(key), do_buckling=False)
        for name in a.order:
            assert np.array_equal(a.cases[name].u, b.cases[name].u), key
            assert np.array_equal(a.cases[name].N, b.cases[name].N), key
    with tempfile.TemporaryDirectory() as tmp:
        paths = bm.write_preset_files(tmp, ["von_mises", "star_dome"])
        assert [os.path.basename(p) for p in paths] == \
            ["von_mises.json", "star_dome.json"]
        m = bm.load_preset_file("star_dome", tmp)
        assert m.geom.is_custom and m.geom.n_members == 24


# ======================================================== self-checking ====
def test_verify_preset_passes_for_every_preset():
    counts = {}
    for key, rows in bm.verify_all().items():
        bad = [r for r in rows if not r["pass"]]
        assert not bad, f"{key}: {bad[:3]}"
        non_identity = [r for r in rows if r["source"] != "identity"]
        assert len(non_identity) >= 3, key
        counts[key] = len(rows)
    assert sum(counts.values()) >= 100
    text = bm.format_rows(bm.verify_preset("von_mises"))
    assert "lambda_cr" in text and "NO" not in text


def test_verify_preset_detects_a_wrong_model():
    """Negative controls: the checks must fail when the model is wrong."""
    m = bm.validation_model_a()
    m.group_sections = ["CHS 76.1 x 4.0"]               # not the hand section
    assert any(not r["pass"] for r in bm.verify_preset("validation_a", m))
    m = bm.star_dome()
    m.geom.nodes[0, 2] += 0.001                          # crown 1 mm too high
    rows = bm.verify_preset("star_dome", m)
    assert any(not r["pass"] for r in rows
               if r["source"].startswith("computed in this work"))
    # Controls the other way: the force-method reference is built from the
    # model's own E, alpha and loads, so an edited (still valid) truss passes.
    m = bm.three_bar_truss()
    m.material_name = "EN S355"
    assert all(r["pass"] for r in bm.verify_preset("three_bar", m))
    m = bm.three_bar_truss(temperature_change=-20.0, settlement_uz=0.0)
    assert all(r["pass"] for r in bm.verify_preset("three_bar", m))
    m = bm.von_mises_truss()
    m.custom_supports[2] = [False, False, False]         # apex free out of plane
    try:
        bm.verify_preset("von_mises", m)
    except Exception:
        pass                                             # mechanism detected
    else:
        raise AssertionError("an unrestrained planar truss was accepted")


def test_presets_run_the_full_pipeline_and_draw():
    for key in bm.preset_keys():
        m = bm.make_preset(key)
        res = run_analysis(m, n_buckling_modes=3)
        txt = summary_text(res)
        assert "ANALYSIS SUMMARY" in txt
        if m.geom.is_custom:
            assert m.geom.title in txt
        assert len(member_table(res)) == m.geom.n_members
        assert len(node_table(res)) == m.geom.n_nodes
        assert len(reaction_table(res)) == len(m.support_nodes()) + 1  # TOTAL
        assert m.summary_lines()
        for r in res.cases.values():
            assert r.equilibrium_residual < 1e-10 and r.energy_error < 1e-10
        fig = Figure()
        ax = fig.add_subplot(111, projection="3d")
        plot.draw_model(ax, m)
        r = res.worst_case_result()
        ax = fig.add_subplot(111, projection="3d")
        plot.draw_axial(ax, m, r.N)
        ax = fig.add_subplot(111, projection="3d")
        plot.draw_deflected(ax, m, r)
        if res.buckling is not None and len(res.buckling.load_factors):
            ax = fig.add_subplot(111, projection="3d")
            plot.draw_buckling_mode(ax, m, res.buckling, 0)


def test_command_line_runs_and_writes():
    buf = io.StringIO()
    with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(buf):
        code = bm.main(["--quiet", "--write", tmp, "von_mises", "three_bar"])
        files = sorted(os.listdir(tmp))
    assert code == 0
    assert files == ["three_bar.json", "von_mises.json"]
    assert "benchmark checks pass" in buf.getvalue()


# =============================================================================
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
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:                                # pragma: no cover
        pass
    raise SystemExit(_run_standalone())
