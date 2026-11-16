"""Supplementary run for Chapter 6 of the GeoDomeLab 2.0 final report.

    python examples/report_supplementary/b42_wall.py [out.json]

Uses the geodome package of this repository (two folders up). Writes its
results to a JSON file next to the script unless a path is given."""
"""Chapter 6 supplementary run: effect of the AISC 360-16 B4.2 design wall (0.93 t) on
the sized 6V case study, without and with re-sizing. Writes only its own (git-ignored) JSON file."""
import sys, copy, json
import numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else os.path.splitext(os.path.abspath(__file__))[0] + ".json"
from geodome import benchmarks as bm
from geodome.postprocess import run_analysis
from geodome.optimize import auto_size

OUT = {}
sized = bm.make_preset("case_study_6v")
auto_size(sized, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
r0 = run_analysis(sized, n_buckling_modes=1, do_buckling=False)
OUT["nominal"] = dict(U=float(r0.governing_util), case=r0.governing_case, mass=float(sized.steel_mass()),
                      sections=list(sized.group_sections))
m = copy.deepcopy(sized)
m.design_thickness_factor = 0.93
r1 = run_analysis(m, n_buckling_modes=1, do_buckling=False)
OUT["b42_same_sections"] = dict(U=float(r1.governing_util), case=r1.governing_case)
m2 = bm.make_preset("case_study_6v")
m2.design_thickness_factor = 0.93
sz = auto_size(m2, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
r2 = run_analysis(m2, n_buckling_modes=1, do_buckling=False)
OUT["b42_resized"] = dict(U=float(r2.governing_util), case=r2.governing_case, mass=float(m2.steel_mass()),
                          sections=list(m2.group_sections), converged=bool(sz.converged))
print(json.dumps(OUT, indent=1, default=str))
json.dump(OUT, open(OUTFILE, "w"), indent=1, default=str)
