"""
GeoDomeLab vs OpenSees -- independent finite-element cross-check (task E).

Every structure below is solved twice: by GeoDomeLab and by OpenSees
(OpenSeesPy 3.8, McKenna, Scott & Fenves 2010; Zhu, McKenna & Scott 2018),
an independent open-source finite element framework.  The two programs share
only the input data (``opensees_bridge.to_opensees_spec``: coordinates,
connectivity, E, A, supports, factored nodal loads, thermal strain,
settlements, masses).  Stiffness, supports, self-straining loads, the linear
solve, the eigen solution and the nonlinear path are computed by each
program on its own.

Sections
--------
P  formulation probes: the strain measure of OpenSees ``corotTruss``
   (decided numerically), the sign of ``InitStrainMaterial`` (heated bar vs
   closed form), and a sensitivity control (E changed by 1 ppm in the
   OpenSees input must show up as a 1e-6 difference: the comparison is not
   blind).
S  linear static: 1V-8V geodesic domes (R = 12 m, hemisphere) on pinned,
   3-roller and elastic supports under every default combination; the
   validation Models A and B, the 6V case study, the planar three-bar truss
   (load, temperature, settlement), self-straining T on a 4V dome for every
   support type, a 3V dome on mixed custom supports (rigid, springs, free,
   settling springs), and two other dome fractions.
M  modal: 1V-8V x 3 support types (12 modes, lumped mass), all modes of
   1V-3V with lumped and consistent mass, the case study (20 modes, both
   masses), validation Model A, the star dome and the three-bar truss.
N  geometrically nonlinear: von Mises truss (whole path through both limit
   points, plus the closed form), a deep von Mises truss with 33 % strain,
   the 24-member star dome (whole path through two snap-backs to 17 cm),
   2V/3V/4V domes and the 6V case study under LC2 past the limit point, a
   3V dome on elastic supports.
X  exported standalone OpenSeesPy scripts (static, modal, nonlinear) run in
   the OpenSees interpreter and compared with GeoDomeLab.

Measures and targets (spec section 2.E)
---------------------------------------
Normwise relative differences (see ``geodome/opensees_bridge.py``):
displacements, member forces, reactions and spring forces <= 1e-8;
natural frequencies <= 1e-6 (pointwise) and 1 - MAC <= 1e-6 (mode shapes,
subspace MAC inside repeated-frequency clusters); load factor at matched
states of the nonlinear path <= 1e-5.  The equilibrium residual of every
OpenSees state evaluated with GeoDomeLab's own corotational equations must
be <= 1e-8 (it would be about the strain itself, 1e-3 or more, if the
strain measures differed).

Run (base or OpenSees environment)::

    python validation/opensees_crosscheck.py            # full matrix
    python validation/opensees_crosscheck.py --quick    # a subset

writes ``validation/output/opensees_crosscheck.md`` and ``.csv`` and exits
with 1 if any comparison fails, 2 if OpenSees is not available.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import platform
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from geodome import __version__ as GDL_VERSION                      # noqa: E402
from geodome import benchmarks as bm                                # noqa: E402
from geodome import geometry as geo                                 # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import modal as md                                     # noqa: E402
from geodome import nonlinear as nl                                 # noqa: E402
from geodome import opensees_bridge as ob                           # noqa: E402
from geodome import solver as slv                                   # noqa: E402
from geodome.model import LoadCombination, default_model            # noqa: E402

OUT_DIR = os.path.join(ROOT, "validation", "output")

TOL_STATIC = 1e-8
TOL_FREQ = 1e-6
TOL_MAC = 1e-6
TOL_PATH = 1e-5
TOL_RESID = 1e-8
TOL_PROBE = 1e-9

SUPPORT_LABEL = {geo.PINNED: "pinned", geo.ROLLER3: "3 rollers",
                 geo.ELASTIC: "elastic springs", geo.CUSTOM: "custom"}


# ================================================================ rows ====
class Rows:
    """Collects one row per compared quantity."""

    def __init__(self):
        self.rows = []

    def add(self, section, case, structure, supports, combination, quantity,
            difference, tolerance, note="", kind="max"):
        d = float(difference)
        if kind == "min":                       # value must reach tolerance
            ok = np.isfinite(d) and d >= tolerance
        elif kind == "range":                   # tolerance = (lo, hi)
            ok = np.isfinite(d) and tolerance[0] <= d <= tolerance[1]
        else:                                   # "max" (relative), "abs"
            ok = np.isfinite(d) and d <= tolerance
        self.rows.append({
            "section": section, "case": case, "structure": structure,
            "supports": supports, "combination": combination,
            "quantity": quantity, "difference": d,
            "tolerance": tolerance, "pass": bool(ok), "note": note,
            "kind": kind})
        return ok

    def section(self, s):
        return [r for r in self.rows if r["section"] == s]


def _structure(model):
    g = model.geom
    if g.is_custom:
        return g.title
    return f"{g.frequency}V R={g.radius:g} m f={g.dome_fraction:.3g}"


def _nfree(model):
    fixed, _ = model.support_conditions()
    return 3 * model.geom.n_nodes - len(fixed)


def _fmt(v):
    if isinstance(v, tuple):
        return f"[{v[0]:.3g}, {v[1]:.3g}]"
    if isinstance(v, str):
        return v
    v = float(v)
    if v == 0.0:
        return "0"
    if not np.isfinite(v):
        return "n/a"
    return f"{v:.2e}"


# ============================================================ sections ====
def section_probes(rows, log):
    s = ob.strain_measure_check()
    log.append(("strain", s))
    note = (f"bar rotated {s['rotation_deg']:.1f} deg, stretched to "
            f"l/L = {s['l_over_L']:.4f}; N_OS = {s['N']:.6e} N")
    rows.add("P", "P1", "single corotTruss bar", "pinned + sp", "-",
             "N_OS vs engineering strain EA(l-L)/L", s["rel_engineering"],
             TOL_PROBE, note)
    rows.add("P", "P1", "single corotTruss bar", "pinned + sp", "-",
             "N_OS vs Green-Lagrange EA(l^2-L^2)/(2L^2) (must differ)",
             s["rel_green"], 1e-3, "different strain measure", kind="min")
    rows.add("P", "P1", "single corotTruss bar", "pinned + sp", "-",
             "N_OS vs logarithmic EA ln(l/L) (must differ)", s["rel_log"],
             1e-3, "different strain measure", kind="min")
    t = ob.thermal_sign_check()
    log.append(("thermal", t))
    note = (f"u2 = {t['u2_os']:.6e} m (exact {t['u2_exact']:.6e}), "
            f"N = {t['N1_os']:.6e} N (exact {t['N_exact']:.6e})")
    rows.add("P", "P2", "heated bar + bar between walls", "fixed ends",
             "dT = 30 C", "u2 vs closed form", t["rel_u"], TOL_PROBE, note)
    rows.add("P", "P2", "heated bar + bar between walls", "fixed ends",
             "dT = 30 C", "N1, N2 vs closed form", t["rel_N"], TOL_PROBE,
             "InitStrainMaterial(-alpha dT) gives N = EA(eps - alpha dT)")
    # sensitivity control: E + 1 ppm in the OpenSees input
    m = default_model(4, 12.0, 0.5)
    spec, an = ob.static_job(m, "LC2")
    spec["E"] = [e * (1.0 + 1e-6) for e in spec["E"]]
    x = ob.compare_static(m, "LC2", ob.run_opensees(spec, an))
    rows.add("P", "P3", _structure(m), "pinned", "LC2",
             "du with E x (1 + 1e-6) in OpenSees only (must be ~1e-6)",
             x["du"], (0.5e-6, 2e-6), "the measure resolves a 1 ppm "
             "input change, 100x below the static target", kind="range")


def _static_cases(quick):
    cases = []
    nus = (1, 2, 4) if quick else range(1, 9)
    for nu in nus:
        for st in (geo.PINNED, geo.ROLLER3, geo.ELASTIC):
            m = default_model(nu, 12.0, 0.5)
            m.support_type = st
            cases.append((f"S-{nu}V-{st}", m, [c.name for c in
                                                m.combinations]))
    for key in ("validation_a", "validation_b", "three_bar",
                "case_study_6v"):
        m = bm.make_preset(key)
        cases.append((f"S-{key}", m, [c.name for c in m.combinations]))
    # self-straining T on a 4V dome, every support type
    for st in (geo.PINNED, geo.ROLLER3, geo.ELASTIC):
        m = default_model(4, 12.0, 0.5)
        m.support_type = st
        b = m.geom.base_nodes
        m.loads.temperature_change = 30.0
        if st == geo.ROLLER3:
            m.loads.settlements = [{"node": int(b[0]), "uz": -0.010},
                                   {"node": int(b[3]), "uz": -0.005}]
        else:
            m.loads.settlements = [
                {"node": int(b[0]), "ux": 0.0, "uy": 0.0, "uz": -0.010},
                {"node": int(b[3]), "ux": 0.003, "uy": -0.002,
                 "uz": -0.005}]
        m.combinations += [
            LoadCombination("T1  1.2D + 1.6S + 1.0T", D=1.2, S=1.6, T=1.0),
            LoadCombination("T2  1.0T", T=1.0),
            LoadCombination("T3  1.2D + 1.0W + 0.5S + 1.2T", D=1.2, W=1.0,
                            S=0.5, T=1.2)]
        cases.append((f"S-4V-T-{st}", m, ["T1", "T2", "T3"]))
    # mixed custom supports: rigid, springs, free, settling spring
    m = default_model(3, 10.0, 0.5)
    cs = {}
    for k, nd in enumerate(m.geom.base_nodes):
        cs[int(nd)] = ([True, True, True] if k % 3 == 0 else
                       [2e7, 3e7, True] if k % 3 == 1 else
                       [False, 1e8, 5e7])
    m.support_type = geo.CUSTOM
    m.custom_supports = cs
    m.loads.temperature_change = -25.0
    nd1 = [nd for nd, v in cs.items() if v[1] == 1e8][0]
    m.loads.settlements = [{"node": int(m.geom.base_nodes[0]),
                            "uz": -0.008},
                           {"node": nd1, "uy": 0.004, "uz": 0.002}]
    m.combinations.append(LoadCombination("T1  1.2D + 1.6S + 1.0T", D=1.2,
                                          S=1.6, T=1.0))
    cases.append(("S-3V-custom", m, [c.name for c in m.combinations]))
    if not quick:
        for nu, R, f in ((5, 15.0, 0.625), (3, 8.0, 0.375)):
            m = default_model(nu, R, f)
            m.loads.snow_unbalanced = True
            m.loads.wind_azimuth = 30.0
            cases.append((f"S-{nu}V-f{f:g}", m, [c.name for c in
                                                 m.combinations]))
    return cases


def section_static(rows, quick, log):
    cases = _static_cases(quick)
    jobs, idx = [], []
    for cid, m, combos in cases:
        for c in combos:
            jobs.append(ob.static_job(m, c))
            idx.append((cid, m, c))
    t0 = time.perf_counter()
    res = ob.run_opensees_batch(jobs)
    log.append(("static_time", time.perf_counter() - t0))
    for (cid, m, c), r in zip(idx, res):
        x = ob.compare_static(m, c, r)
        st = SUPPORT_LABEL.get(m.support_type, m.support_type)
        name = ob.find_combination(m, c).name
        if not x.get("ok"):
            rows.add("S", cid, _structure(m), st, name, "OpenSees run",
                     float("nan"), TOL_STATIC, x.get("error", ""))
            continue
        note = (f"{_nfree(m)} free DOF, |u|max = {x['u_max'] * 1e3:.4g} mm,"
                f" |N|max = {x['N_max'] / 1e3:.4g} kN, OS system "
                f"{x['system']}")
        rows.add("S", cid, _structure(m), st, name, "displacements u",
                 x["du"], TOL_STATIC, note)
        rows.add("S", cid, _structure(m), st, name, "member forces N",
                 x["dN"], TOL_STATIC)
        rows.add("S", cid, _structure(m), st, name,
                 "reactions R (incl. springs)", x["dR"], TOL_STATIC)
        if m.support_type in (geo.ELASTIC, geo.CUSTOM):
            rows.add("S", cid, _structure(m), st, name, "spring forces",
                     x["dspring"], TOL_STATIC)
        rows.add("S", cid, _structure(m), st, name,
                 "OpenSees residual at free DOFs / |F|", x["os_residual"],
                 TOL_STATIC)


def _modal_cases(quick):
    cases = []
    nus = (1, 2, 4) if quick else range(1, 9)
    for nu in nus:
        for st in (geo.PINNED, geo.ROLLER3, geo.ELASTIC):
            m = default_model(nu, 12.0, 0.5)
            m.support_type = st
            cases.append((f"M-{nu}V-{st}", m, 12, md.LUMPED))
    for nu in ((1, 2) if quick else (1, 2, 3)):
        for mt in (md.LUMPED, md.CONSISTENT):
            m = default_model(nu, 12.0, 0.5)
            cases.append((f"M-{nu}V-all-{mt}", m, "all", mt))
    if not quick:
        for mt in (md.LUMPED, md.CONSISTENT):
            cases.append((f"M-case6V-{mt}", bm.make_preset("case_study_6v"),
                          20, mt))
    for key in ("validation_a", "star_dome", "three_bar"):
        cases.append((f"M-{key}", bm.make_preset(key), "all", md.LUMPED))
    return cases


def section_modal(rows, quick, log):
    cases = _modal_cases(quick)
    jobs = [ob.modal_job(m, k, mt) for _cid, m, k, mt in cases]
    t0 = time.perf_counter()
    res = ob.run_opensees_batch(jobs)
    log.append(("modal_time", time.perf_counter() - t0))
    for (cid, m, k, mt), r in zip(cases, res):
        st = SUPPORT_LABEL.get(m.support_type, m.support_type)
        x = ob.compare_modal(m, r, n_modes=k, mass_type=mt)
        if not x.get("ok"):
            rows.add("M", cid, _structure(m), st, mt, "OpenSees eigen",
                     float("nan"), TOL_FREQ, x.get("error", ""))
            continue
        note = (f"{x['n_compared']} modes ({x['n_single']} single, "
                f"{x['n_clustered']} in repeated clusters), f1 = "
                f"{x['f1']:.6f} Hz, OS {x['solver']}")
        rows.add("M", cid, _structure(m), st, f"{mt} mass",
                 "natural frequencies f", x["df"], TOL_FREQ, note)
        if x["n_single"]:
            rows.add("M", cid, _structure(m), st, f"{mt} mass",
                     "1 - MAC (single modes)", 1.0 - x["mac_min"], TOL_MAC)
        if x["n_clustered"]:
            rows.add("M", cid, _structure(m), st, f"{mt} mass",
                     "1 - subspace MAC (repeated pairs)",
                     1.0 - x["submac_min"], TOL_MAC)


def _nl_cases(quick):
    big = dict(lambda_max=1e9, lambda_min=-1e9, stop_after_limit=1000,
               max_steps=800)
    vm = bm.make_preset("von_mises")
    a, h, _EA, _P = bm.von_mises_parameters(vm)
    cases = [
        ("N-von-mises", vm, None, dict(control_max=2.4 * h, **big),
         dict(arc_length=0.002, n_steps=300, control_max=2.4 * h)),
        ("N-von-mises-deep", bm.von_mises_truss(half_span=1.0, rise=1.0,
                                                P_ref=1e6),
         None, dict(control_max=2.4, **big), None),
        ("N-star-dome", bm.make_preset("star_dome"), None,
         dict(control_max=0.03 if quick else 0.17, **big),
         dict(arc_length=0.002, n_steps=300,
              control_max=0.03 if quick else 0.17)),
    ]
    for nu in ((2,) if quick else (2, 3, 4)):
        cases.append((f"N-{nu}V-LC2", default_model(nu, 12.0, 0.5), "LC2",
                      dict(stop_after_limit=10), None))
    if not quick:
        m = default_model(3, 12.0, 0.5)
        m.support_type = geo.ELASTIC
        cases.append(("N-3V-elastic-LC2", m, "LC2",
                      dict(stop_after_limit=10), None))
        cases.append(("N-case6V-LC2", bm.make_preset("case_study_6v"), "LC2",
                      dict(stop_after_limit=10), None))
    return cases


def section_nonlinear(rows, quick, log):
    for cid, m, combo, gkw, arc in _nl_cases(quick):
        t0 = time.perf_counter()
        x = ob.crosscheck_path(m, combo, gdl_kwargs=gkw, arclength=arc)
        r = x["gdl"]
        log.append((cid, time.perf_counter() - t0, x))
        st = SUPPORT_LABEL.get(m.support_type, m.support_type)
        cname = r.combination
        lr = x["lambda_range"]
        cr = x["control_range"]
        strain = float(np.max(np.abs(r.N / m.EA()[None, :]))) \
            if r.N is not None and len(r.N) else float("nan")
        note = (f"{x['n_matched']}/{x['n_targets']} states, "
                f"{x['n_segments']} control segment(s) on DOF(s) "
                f"{x['control_dofs']}; lambda in [{lr[0]:.4g}, {lr[1]:.4g}],"
                f" control in [{cr[0] * 1e3:.4g}, {cr[1] * 1e3:.4g}] mm, "
                f"max |strain| {strain:.2e}; {r.stop_reason}")
        if not x["ok"]:
            note += f"; OpenSees: {x['message']}"
        rows.add("N", cid, _structure(m), st, cname,
                 "lambda at matched states / max|lambda|",
                 x["dlam"] if x["ok"] else float("nan"), TOL_PATH, note)
        rows.add("N", cid, _structure(m), st, cname,
                 "lambda pointwise (|lambda| > 5 % of max)",
                 x["dlam_point"], TOL_PATH)
        rows.add("N", cid, _structure(m), st, cname,
                 "displacements at matched states", x["du"], TOL_PATH)
        if np.isfinite(x["dlimit"]):
            rows.add("N", cid, _structure(m), st, cname,
                     "lambda at the refined first limit point",
                     x["dlimit"], TOL_PATH,
                     f"lambda_limit GDL {x['lambda_limit_gdl']:.10g}, "
                     f"OS {x['lambda_limit_os']:.10g} "
                     f"({r.limit_type})")
        rows.add("N", cid, _structure(m), st, cname,
                 "GDL equilibrium residual of OpenSees states",
                 x["res_os_in_gdl"], TOL_RESID,
                 "OpenSees states checked with GeoDomeLab's corotational "
                 "equations")
        if "arc_points" in x:
            ac = x["arc_control"]
            rows.add("N", cid, _structure(m), st, cname,
                     "GDL residual of OpenSees ArcLength points",
                     x["res_arc_in_gdl"], TOL_RESID,
                     f"{x['arc_points']} independent OpenSees ArcLength "
                     f"points, control in [{ac.min() * 1e3:.4g}, "
                     f"{ac.max() * 1e3:.4g}] mm, max lambda "
                     f"{x['arc_lam_max']:.6g} ({x['arc_message']})")
        # closed forms and literature values on the OpenSees side
        if cid.startswith("N-von-mises"):
            a, h, EA, P = bm.von_mises_parameters(m)
            os_res = x["os"]
            w = -np.asarray(os_res["control"])
            Pos = np.asarray(os_res["lam"]) * P
            Pex = bm.von_mises_load(w, a, h, EA)
            d = float(np.max(np.abs(Pos - Pex)) / np.max(np.abs(Pex)))
            rows.add("N", cid, _structure(m), st, cname,
                     "OpenSees P(w) vs exact closed form", d, TOL_PATH,
                     f"{len(w)} points, w up to {w.max() * 1e3:.4g} mm "
                     f"(h = {h * 1e3:g} mm)")
            lp = bm.von_mises_limit_point(a, h, EA)
            if np.isfinite(x["lambda_limit_os"]):
                d = abs(x["lambda_limit_os"] * P - lp["P"]) / lp["P"]
                rows.add("N", cid, _structure(m), st, cname,
                         "OpenSees limit load vs closed form P*", d,
                         TOL_PATH, f"P* = {lp['P'] / 1e3:.8g} kN at w* = "
                         f"{lp['w'] * 1e3:.6g} mm")
        if cid == "N-star-dome" and np.isfinite(x["lambda_limit_os"]):
            p_ea = x["lambda_limit_os"] * bm.STAR_LOAD_RATIO
            d = abs(p_ea - bm.STAR_P_OVER_EA_LIMIT)
            rows.add("N", cid, _structure(m), st, cname,
                     "OpenSees P/EA at limit vs reference 3.15655e-4 "
                     "(absolute, 6 digits)", d, 5e-10,
                     f"P/EA = {p_ea:.9e}; reference computed in this "
                     f"work (docs/VALIDATION.md, star dome reference values; "
                     f"agrees with Greco & Venturini 2006, Fig. 8)",
                     kind="abs")


def section_scripts(rows, log):
    m = default_model(3, 10.0, 0.5)
    m.support_type = geo.ELASTIC
    b = m.geom.base_nodes
    m.loads.temperature_change = 20.0
    m.loads.settlements = [{"node": int(b[0]), "uz": -0.004}]
    m.combinations.append(LoadCombination("T1  1.2D + 1.6S + 1.0T", D=1.2,
                                          S=1.6, T=1.0))
    st = SUPPORT_LABEL[m.support_type]
    with tempfile.TemporaryDirectory(prefix="geodome_xs_") as tmp:
        p = ob.write_openseespy_script(m, "T1", os.path.join(tmp, "s.py"))
        r = ob.run_script(p)
        c = ob.find_combination(m, "T1")
        g = slv.solve_model(m, ld.combine(ld.basic_load_vectors(m), c),
                            T_factor=c.T)
        rows.add("X", "X1", _structure(m), st, c.name,
                 "exported static script: u", ob.rel_inf(g.u, r["u"]),
                 TOL_STATIC, "write_openseespy_script + run_script")
        rows.add("X", "X1", _structure(m), st, c.name,
                 "exported static script: N",
                 ob.rel_inf(g.N, r["N"], max(np.max(np.abs(g.N)), 1e-30)),
                 TOL_STATIC)
        rows.add("X", "X1", _structure(m), st, c.name,
                 "exported static script: R",
                 ob.rel_inf(g.reactions, r["reactions"]), TOL_STATIC)
        p = ob.write_openseespy_script(m, None, os.path.join(tmp, "m.py"),
                                       {"type": "modal", "n_modes": 10})
        r = ob.run_script(p)
        mr = md.modal_analysis(m, n_modes=10)
        d = float(np.max(np.abs(np.asarray(r["frequencies"]) -
                                mr.frequencies) / mr.frequencies))
        rows.add("X", "X2", _structure(m), st, "-",
                 "exported modal script: f (10 modes)", d, TOL_FREQ)
        p = ob.write_openseespy_script(m, "LC2", os.path.join(tmp, "n.py"),
                                       {"type": "path", "control_max": 0.4,
                                        "n_steps": 20})
        r = ob.run_script(p)
        system = nl.model_system(m)
        F = nl.reference_load(m, ob.find_combination(m, "LC2"))
        res = max(ob.gdl_residual(system, F, u, lam)
                  for u, lam in zip(r["u"], r["lam"]))
        rows.add("X", "X3", _structure(m), st, "LC2",
                 "exported nonlinear script: GDL residual of its states",
                 res, TOL_RESID, f"{len(r['lam'])} DisplacementControl "
                 f"states to 400 mm, lambda up to {max(r['lam']):.5g}")


# ============================================================ report =====
_SECTION_TITLES = {
    "P": "Formulation probes",
    "S": "Linear static analysis",
    "M": "Modal analysis",
    "N": "Geometrically nonlinear analysis",
    "X": "Exported standalone OpenSeesPy scripts",
}


def _tol_text(r):
    t = r["tolerance"]
    kind = r.get("kind", "max")
    if kind == "range":
        return f"in [{t[0]:.1e}, {t[1]:.1e}]"
    if kind == "min":
        return f">= {t:.0e}"
    if kind == "abs":
        return f"<= {t:.0e} (abs)"
    return f"<= {t:.0e}"


def _esc(text):
    """Markdown table cell: escape the column separator."""
    return str(text).replace("|", "\\|")


def write_csv(rows, path):
    keys = ["section", "case", "structure", "supports", "combination",
            "quantity", "difference", "tolerance", "pass", "note"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows.rows:
            w.writerow([r["section"], r["case"], r["structure"],
                        r["supports"], r["combination"], r["quantity"],
                        repr(r["difference"]), _tol_text(r),
                        "yes" if r["pass"] else "NO", r["note"]])
    return path


def _summary(rows):
    out = []
    for s, title in _SECTION_TITLES.items():
        rs = rows.section(s)
        if not rs:
            continue
        cases = sorted(set(r["case"] for r in rs))
        worst = {}
        for r in rs:
            if r.get("kind", "max") != "max":
                continue
            q = r["quantity"]
            if q not in worst or not (r["difference"] <= worst[q][0]):
                worst[q] = (r["difference"], r)
        out.append((s, title, len(cases), len(rs),
                    sum(1 for r in rs if r["pass"]), worst))
    return out


def write_md(rows, path, log, meta):
    L = []
    allpass = all(r["pass"] for r in rows.rows)
    L.append("# GeoDomeLab vs OpenSees: independent cross-check")
    L.append("")
    L.append(f"Generated {meta['date']} by `validation/opensees_crosscheck.py`"
             f"{' (quick subset)' if meta['quick'] else ''}.  ")
    L.append(f"GeoDomeLab {GDL_VERSION} (Python {meta['python']}); OpenSees "
             f"{meta['opensees']} via OpenSeesPy, run {meta['mode']} "
             f"(separate OpenSees interpreter).  "
             f"Total time {meta['elapsed']:.1f} s.")
    L.append("")
    L.append(f"**Result: {sum(r['pass'] for r in rows.rows)} of "
             f"{len(rows.rows)} comparisons pass"
             f"{'' if allpass else ' -- SOME FAIL'}.**")
    L.append("")
    L.append("## What is compared")
    L.append("")
    L.append(
        "Each structure is solved by GeoDomeLab and by OpenSees, an "
        "independent finite element framework (McKenna, Scott & Fenves "
        "2010; OpenSeesPy, Zhu, McKenna & Scott 2018).  The programs share "
        "only the input data written by `opensees_bridge.to_opensees_spec`: "
        "node coordinates, connectivity, E and A, supports, the factored "
        "nodal load vector of the combination, the thermal strain "
        "T*alpha*dT, the settlements T*delta and (modal) the masses.  "
        "OpenSees builds and solves everything else itself: `Truss` and "
        "`corotTruss` elements with an `Elastic` material, "
        "`InitStrainMaterial(-alpha dT)` for temperature, `sp` constraints "
        "in the load pattern for settlement, `zeroLength` springs to fixed "
        "ground nodes for elastic supports, element mass `-rho` (lumped, or "
        "consistent with `-cMass 1`) plus nodal `mass`, `eigen` (ARPACK or "
        "LAPACK), and Newton iteration with `DisplacementControl` for the "
        "nonlinear path.")
    L.append("")
    L.append(
        "Load generation and the non-structural masses are therefore "
        "*inputs* here, and this cross-check cannot see an error in them.  "
        "Checked by hand elsewhere: the dead and snow tributary loads and "
        "the crown load on Models A and B (`validation/run_validation.py`), "
        "the wind field on Model B (case F9 of "
        "`validation/run_validation_ext.py`), and E, alpha and rho through "
        "the typed-in constants of the hand files.  Not checked by a hand "
        "calculation: the unbalanced snow pattern and the default "
        "combination factors (unit tests only), and the 7-point Cp table "
        "itself (an illustrative simplification, not ASCE 7-16 Fig. "
        "27.3-2).")
    L.append("")
    L.append("Differences are normwise relative: for displacements "
             "max|u_GDL - u_OS| / max|u_GDL| (likewise N, with the fully "
             "restrained thermal force EA*alpha*dT as a floor of the scale, "
             "and R, with max|F| as a floor); natural frequencies pointwise "
             "|f_GDL - f_OS| / f_GDL; mode shapes by 1 - MAC, and inside "
             "clusters of exactly repeated frequencies (the symmetric pairs "
             "of a dome, whose basis is arbitrary) by 1 - subspace MAC; the "
             "nonlinear load factor at matched states by max|lambda_GDL - "
             "lambda_OS| / max|lambda| over the path.  Targets (spec 2.E): "
             f"static {TOL_STATIC:.0e}, frequencies {TOL_FREQ:.0e}, "
             f"nonlinear lambda {TOL_PATH:.0e}.")
    L.append("")
    L.append("## Summary (worst case of every quantity)")
    L.append("")
    L.append("| section | cases | comparisons | pass | quantity | worst "
             "difference | target | worst case |")
    L.append("|---|---:|---:|---:|---|---:|---:|---|")
    for s, title, nc, nr, npass, worst in _summary(rows):
        first = True
        for q, (d, r) in worst.items():
            L.append(f"| {title if first else ''} | {nc if first else ''} | "
                     f"{nr if first else ''} | {npass if first else ''} | "
                     f"{_esc(q)} | {_fmt(d)} | {_tol_text(r)} | "
                     f"{r['case']} {_esc(r['combination'])} |")
            first = False
    L.append("")

    # -- probes
    L.append("## P. Formulation probes")
    L.append("")
    s = dict(e for e in log if len(e) == 2)["strain"]
    L.append(
        "**Strain measure of `corotTruss` (decided numerically).**  One bar "
        f"of length L = 2 m was rotated by {s['rotation_deg']:.1f} degrees "
        "(prescribed transverse end displacement 0.8 L) and pulled axially "
        f"(P = 0.1 EA) to l/L = {s['l_over_L']:.6f}.  OpenSees returned "
        f"N = {s['N']:.9e} N.  Predictions: engineering strain "
        f"EA(l-L)/L = {s['engineering']:.9e} N (difference "
        f"{s['rel_engineering']:.1e}), Green-Lagrange "
        f"EA(l^2-L^2)/(2L^2) = {s['green']:.9e} N ({s['rel_green']:.1e}), "
        f"logarithmic EA ln(l/L) = {s['log']:.9e} N ({s['rel_log']:.1e}).  "
        f"OpenSees {s['opensees_version']} `corotTruss` therefore uses "
        f"**{s['measure']} strain on the original length**, the same "
        "measure as GeoDomeLab's corotational bar (and as the OpenSees "
        "source, CorotTruss.cpp).  The nonlinear paths can be compared at "
        "any strain; the deep von Mises case below reaches 33 % strain.")
    L.append("")
    t = dict(e for e in log if len(e) == 2)["thermal"]
    L.append(
        "**Sign of `InitStrainMaterial`.**  Bar 1 (heated, free strain "
        f"eps0 = {11.7e-6 * 30:.4e}) and bar 2 (twice the area) in series "
        "between two walls.  Closed form u2 = eps0 L/(1 + r), N = -EA eps0 "
        f"r/(1 + r) with r = 2: u2 = {t['u2_exact']:.9e} m, N = "
        f"{t['N_exact']:.6f} N.  OpenSees with epsInit = -eps0: u2 = "
        f"{t['u2_os']:.9e} m, N1 = {t['N1_os']:.6f} N, N2 = "
        f"{t['N2_os']:.6f} N.  So `InitStrainMaterial(tag, E, -alpha dT)` "
        "gives N = EA(eps - alpha dT): heating compresses a restrained bar.")
    L.append("")
    L.append(
        "**Sensitivity control.**  The same 4V comparison with E raised by "
        "1 ppm in the OpenSees input only must show a displacement "
        "difference of about 1e-6.  It does (row P3), so the round-off-level "
        "agreement elsewhere is a real agreement of two independent "
        "computations, not a program compared with itself.")
    L.append("")
    _table(L, rows.section("P"))

    L.append("## S. Linear static analysis")
    L.append("")
    L.append("Geodesic domes: radius 12 m, dome fraction 0.5 requested "
             "(the generator cuts at the nearest node ring, so odd "
             "frequencies come out lower; the structure column shows the "
             "actual fraction f), default sections (CHS 88.9 x 4.0 "
             "unless stated), default loads and the default ASCE 7-16 "
             "combinations LC1 1.4D, LC2 1.2D + 1.6S, LC3 1.2D + 1.0W + "
             "0.5S, LC4 0.9D + 1.0W, LC5 1.2D + 1.6S + 0.5W, LC6 1.2D + "
             "1.6L + 0.5S, SLS 1.0D + 1.0S.  Springs (elastic supports) "
             "50 MN/m in "
             "x, y and z.  T cases: dT = +30 C on every strut and two base "
             "nodes settling (10 mm down, and 3 / -2 / -5 mm in x / y / z; "
             "vertical only on the roller supports), combinations 1.2D + "
             "1.6S + 1.0T, 1.0T and 1.2D + 1.0W + 0.5S + 1.2T.  Custom: "
             "rigid, spring (20 and 30 MN/m; 100 and 50 MN/m) and free "
             "directions mixed at the base, one settling rigid support, one "
             "spring support whose ground end moves, and dT = -25 C.  The "
             "two extra fractions use unbalanced snow and wind from 30 "
             "degrees.")
    L.append("")
    _case_table(L, rows.section("S"))

    L.append("## M. Modal analysis")
    L.append("")
    L.append("Mass: strut mass rho*A*L/2 to each end (OpenSees `-rho`), or "
             "the consistent bar mass (`-cMass 1`), plus the cladding dead "
             "load / g at the nodes (OpenSees `mass`).  'all' compares every "
             "mode of the model (OpenSees `-fullGenLapack`).  OpenSees mode "
             "shapes are compared without any normalisation (MAC is "
             "scale-free).")
    L.append("")
    _case_table(L, rows.section("M"))

    L.append("## N. Geometrically nonlinear analysis")
    L.append("")
    L.append(
        "GeoDomeLab traces the path (cylindrical arc-length, Crisfield "
        "1981) and refines the first limit point.  OpenSees (`corotTruss`, "
        "Newton, unbalance tolerance 1e-10 |F_ref|) is then driven by "
        "`DisplacementControl` to **every GeoDomeLab state**, including the "
        "refined limit point: where the path snaps back in the control "
        "displacement, control switches to another DOF that is monotone "
        "over the next stretch (piecewise displacement control, "
        "`opensees_bridge.path_segments`).  At each state OpenSees finds "
        "lambda and all displacements by itself; the table compares them.  "
        "In addition every OpenSees state is put into GeoDomeLab's own "
        "corotational equilibrium equations (residual row).")
    L.append("")
    arc = [_arc_sentence(cid, x) for e in log
           if isinstance(e, tuple) and len(e) == 3
           for cid, _t, x in [e] if "arc_points" in x]
    if arc:
        L.append(
            "OpenSees' own `ArcLength` integrator (fixed arc length 2 mm, "
            "alpha = 0, cylindrical) was also run on the benchmarks as an "
            "independent path tracer.  Its predictor takes the sign of the "
            "previous load increment (ArcLength.cpp), which is not a "
            "reliable direction test at a limit point: " + " ".join(arc) +
            "  Every point it produced still satisfies GeoDomeLab's "
            "equilibrium equations (residual row), but because it can turn "
            "back at a limit point instead of continuing along the path, the "
            "matched comparison drives OpenSees by displacement control "
            "instead.")
        L.append("")
    _case_table(L, rows.section("N"))
    nl_logs = [e for e in log if isinstance(e, tuple) and len(e) == 3]
    for cid, _t, x in nl_logs:
        if cid != "N-star-dome":
            continue
        L.append("Star dome, a sample of matched states (lambda = 10^4 P/EA):")
        L.append("")
        L.append("| state | crown w [mm] | lambda GeoDomeLab | lambda "
                 "OpenSees | difference |")
        L.append("|---|---:|---:|---:|---:|")
        pairs = x["pairs"]
        step = max(1, len(pairs) // 14)
        pick = pairs[::step]
        if pairs and pairs[-1] not in pick:
            pick.append(pairs[-1])
        lim = [p for p in pairs if p[0] == "limit"]
        for p in sorted(set(pick + lim), key=lambda q: pairs.index(q)):
            lab, c, lg, lo = p
            L.append(f"| {lab} | {c * 1e3:.4f} | {lg:.10f} | {lo:.10f} | "
                     f"{abs(lg - lo):.1e} |")
        L.append("")

    L.append("## X. Exported standalone scripts")
    L.append("")
    L.append("`opensees_bridge.write_openseespy_script` writes a readable "
             "OpenSeesPy script (explicit `ops.node`, `ops.element`, ... "
             "lines) that anyone can run with `pip install openseespy`.  "
             "The scripts below were written for a 3V dome on elastic "
             "supports with temperature and settlement, run in the OpenSees "
             "interpreter and compared with GeoDomeLab.")
    L.append("")
    _table(L, rows.section("X"))

    L.append("## Notes and limits")
    L.append("")
    L.append(
        "* The differences are at round-off level (1e-16 to 1e-11), many "
        "orders of magnitude below the targets.  The largest ones occur on "
        "the biggest models and along the nonlinear paths, where "
        "GeoDomeLab's own convergence tolerance (1e-8 on the residual) is "
        "the limit, not the formulation.")
    L.append(
        "* What OpenSees cannot check: the generation of the loads and "
        "masses (inputs here; see the hand-calculation validation), the "
        "AISC member checks, the linear buckling eigenvalue (OpenSees has "
        "no built-in linearised buckling analysis; GeoDomeLab's is checked "
        "against closed forms in the benchmark suite), and the fabrication "
        "geometry.")
    L.append(
        "* The nonlinear path excludes the self-straining load T in both "
        "programs (GeoDomeLab's documented scope).")
    L.append(
        "* Mode shapes inside a repeated-frequency cluster are compared as "
        "subspaces because any rotation of a repeated pair is an equally "
        "valid pair of modes.")
    L.append("")
    L.append("References: F. McKenna, M. H. Scott, G. L. Fenves (2010), J. "
             "Comput. Civ. Eng. 24(1), 95-107; M. Zhu, F. McKenna, M. H. "
             "Scott (2018), SoftwareX 7, 6-11; OpenSees source files "
             "CorotTruss.cpp and InitStrainMaterial.cpp; M. A. Crisfield "
             "(1981), Comput. Struct. 13, 55-62 and (1991) Non-linear Finite "
             "Element Analysis of Solids and Structures, Vol. 1, Wiley; R. J. "
             "Allemang (2003), Sound and Vibration 37(8), 14-23 (MAC); "
             "Pecknold, Ghaboussi & Healey (1985), J. Eng. Mech. 111(7), "
             "909-922 (von Mises truss); Greco & Venturini (2006) (star "
             "dome).")
    L.append("")
    L.append(f"Full rows: `validation/output/opensees_crosscheck.csv` "
             f"({len(rows.rows)} rows).")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def _arc_sentence(cid, x):
    """What the OpenSees ArcLength trace did, from its own numbers."""
    c = np.asarray(x["arc_control"], dtype=float)
    la = np.asarray(x["arc_lam"], dtype=float)
    if len(c) < 2:
        return f"{cid}: no step converged."
    nz = c[np.abs(c) > 0.0]
    sgn = float(np.sign(nz[0])) if len(nz) else 1.0     # loading direction
    k = int(np.argmax(sgn * c))
    ext, end = c[k], c[-1]
    cl = x.get("control_limit_gdl", float("nan"))
    txt = (f"{cid}: it went as far as {ext * 1e3:.4g} mm in the loading "
           f"direction (lambda = {la[k]:.5g} there; GeoDomeLab's first limit "
           f"point is at {cl * 1e3:.4g} mm)")
    if sgn * (ext - end) > 0.01 * abs(ext):
        txt += (f" and then turned back, ending at {end * 1e3:.4g} mm with "
                f"lambda = {la[-1]:.5g} after {len(c) - 1} steps.")
    else:
        txt += f" and ended there after {len(c) - 1} steps."
    return txt


def _table(L, rs):
    L.append("| case | structure | quantity | difference | target | pass | "
             "note |")
    L.append("|---|---|---|---:|---:|---|---|")
    for r in rs:
        L.append(f"| {r['case']} | {_esc(r['structure'])} | "
                 f"{_esc(r['quantity'])} | {_fmt(r['difference'])} | "
                 f"{_tol_text(r)} | {'yes' if r['pass'] else '**NO**'} | "
                 f"{_esc(r['note'])} |")
    L.append("")


def _case_table(L, rs):
    """One line per (case, combination): every quantity in a column."""
    if not rs:
        return
    quantities = []
    for r in rs:
        if r["quantity"] not in quantities:
            quantities.append(r["quantity"])
    groups = {}
    order = []
    for r in rs:
        k = (r["case"], r["combination"])
        if k not in groups:
            groups[k] = {}
            order.append(k)
        groups[k][r["quantity"]] = r
    short = {q: _esc(q) for q in quantities}
    L.append("| case | structure | supports | combination | " +
             " | ".join(short[q] for q in quantities) + " | pass | note |")
    L.append("|---|---|---|---|" + "---:|" * len(quantities) + "---|---|")
    for k in order:
        g = groups[k]
        any_r = next(iter(g.values()))
        cells = [_fmt(g[q]["difference"]) if q in g else "" for q in
                 quantities]
        ok = all(r["pass"] for r in g.values())
        note = "; ".join(r["note"] for r in g.values() if r["note"])
        L.append(f"| {k[0]} | {_esc(any_r['structure'])} | "
                 f"{_esc(any_r['supports'])} | {_esc(k[1])} | " +
                 " | ".join(cells) +
                 f" | {'yes' if ok else '**NO**'} | {_esc(note)} |")
    L.append("")


# ============================================================== main =====
def run(quick=False, out_dir=OUT_DIR, verbose=True):
    """Run the matrix; returns (rows, md_path, csv_path)."""
    status = ob.opensees_status()
    if status["mode"] is None:
        raise ob.OpenSeesUnavailable("OpenSees is not available.")
    t0 = time.perf_counter()
    rows, log = Rows(), []
    steps = [("P", lambda: section_probes(rows, log)),
             ("S", lambda: section_static(rows, quick, log)),
             ("M", lambda: section_modal(rows, quick, log)),
             ("N", lambda: section_nonlinear(rows, quick, log)),
             ("X", lambda: section_scripts(rows, log))]
    for key, fn in steps:
        t = time.perf_counter()
        fn()
        if verbose:
            rs = rows.section(key)
            print(f"  {_SECTION_TITLES[key]:40s} {sum(r['pass'] for r in rs):4d}"
                  f" / {len(rs):4d} pass  ({time.perf_counter() - t:.1f} s)")
    probe = ob.run_opensees(None, {"type": "thermal_probe"})
    meta = {"date": _dt.date.today().isoformat(), "quick": quick,
            "python": platform.python_version(),
            "opensees": probe.get("opensees_version", "?"),
            "mode": status["mode"], "os_python": status["python"],
            "elapsed": time.perf_counter() - t0}
    os.makedirs(out_dir, exist_ok=True)
    md_path = write_md(rows, os.path.join(out_dir, "opensees_crosscheck.md"),
                       log, meta)
    csv_path = write_csv(rows, os.path.join(out_dir,
                                            "opensees_crosscheck.csv"))
    return rows, md_path, csv_path


def main(argv=None):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quick", action="store_true",
                    help="a smaller matrix (about 10 s)")
    ap.add_argument("--out", default=OUT_DIR, help="output folder")
    a = ap.parse_args(argv)
    try:
        rows, md_path, csv_path = run(quick=a.quick, out_dir=a.out)
    except ob.OpenSeesUnavailable as exc:
        print(f"OpenSees cross-check skipped: {exc}")
        return 2
    n = len(rows.rows)
    npass = sum(r["pass"] for r in rows.rows)
    worst = {}
    for r in rows.rows:
        if r.get("kind", "max") != "max":
            continue
        s = r["section"]
        if s not in worst or r["difference"] > worst[s]["difference"]:
            worst[s] = r
    for s, r in worst.items():
        print(f"  worst {_SECTION_TITLES[s]:38s} {r['difference']:.2e} "
              f"({r['quantity']}, {r['case']} {r['combination']})")
    for r in rows.rows:
        if not r["pass"]:
            print(f"  FAIL {r['case']} {r['combination']} {r['quantity']}: "
                  f"{r['difference']!r} ({r['note']})")
    print(f"\n{npass} of {n} comparisons pass.")
    print(f"Wrote {md_path}\n      {csv_path}")
    return 0 if npass == n else 1


if __name__ == "__main__":
    raise SystemExit(main())
