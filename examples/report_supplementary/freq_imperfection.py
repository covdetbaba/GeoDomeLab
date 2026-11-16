"""Supplementary run for Chapter 6 of the GeoDomeLab 2.0 final report.

    python examples/report_supplementary/freq_imperfection.py [out.json]

Uses the geodome package of this repository (two folders up). Writes its
results to a JSON file next to the script unless a path is given."""
"""Chapter 6 supplementary run: perfect and imperfect limit loads of the auto-sized
4V, 6V and 8V hemispheres (R = 15 m, case-study loads, LC2), with span-based and
hub-scale (10 mm) amplitudes, and the member-strength load factor of each.
Does not write into the repo."""
import sys, json, time
import numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else os.path.splitext(os.path.abspath(__file__))[0] + ".json"
from geodome import benchmarks as bm, study as st
from geodome import nonlinear as nl
from geodome.optimize import auto_size
from geodome.postprocess import run_analysis

base = bm.make_preset("case_study_6v")
OUT = []
for nu in (4, 6, 8):
    t0 = time.perf_counter()
    m = st.study_model(base, nu, 15.0, 0.5, section="CHS 114.3 x 5.0")
    auto_size(m, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
    res = run_analysis(m, n_buckling_modes=6)
    span = float(m.geom.base_diameter)
    U_lc2 = float(np.max(res.util[[n for n in res.order if n.startswith("LC2")][0]]))
    r = nl.run_nonlinear(m, "LC2", method="arc-length", control="crown_uz")
    rows = nl.imperfection_study(m, "LC2", amplitudes=[0.010, span / 1000, span / 500, span / 300])
    OUT.append(dict(nu=nu, mass=float(m.steel_mass()), sections=list(m.group_sections),
                    U_gov=float(res.governing_util), gov=res.governing_case,
                    U_LC2=U_lc2, lam_member_LC2=1.0 / U_lc2,
                    lam_cr_lowest=float(res.lambda_cr), lam_cr_case=res.lambda_cr_case,
                    lam_limit=float(r.lambda_limit), ratio=float(r.lambda_limit / r.lambda_cr),
                    lam_cr_nl=float(r.lambda_cr),
                    lam_member_design=r.lambda_member_design, lam_member_euler=r.lambda_member_euler,
                    imperf=[dict(a=x["amplitude_m"], lam=x["lambda_limit"], ratio=x["ratio"],
                                 sign=x["sign"], lam_plus=x["lambda_plus"], lam_minus=x["lambda_minus"])
                            for x in rows],
                    time_s=time.perf_counter() - t0))
    print(json.dumps(OUT[-1], indent=1, default=str), flush=True)
json.dump(OUT, open(OUTFILE, "w"), indent=1, default=str)
