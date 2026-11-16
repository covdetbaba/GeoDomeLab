"""Supplementary run for Chapter 6 of the GeoDomeLab 2.0 final report.

    python examples/report_supplementary/mode_shapes.py [out.json]

Uses the geodome package of this repository (two folders up). Writes its
results to a JSON file next to the script unless a path is given."""
"""Where is the first LC2 buckling mode of the sized 4V, 6V, 8V domes, and how much
of the crown rise does each amplitude remove? Writes only its own (git-ignored) JSON file."""
import sys, json
import numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else os.path.splitext(os.path.abspath(__file__))[0] + ".json"
from geodome import benchmarks as bm, study as st
from geodome import nonlinear as nl
from geodome.optimize import auto_size

base = bm.make_preset("case_study_6v")
OUT = []
for nu in (4, 6, 8):
    m = st.study_model(base, nu, 15.0, 0.5, section="CHS 114.3 x 5.0")
    auto_size(m, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
    g = m.geom
    P = g.nodes
    c = int(g.crown_node)
    nb = sorted({int(j) for a, b in g.members for (i, j) in ((a, b), (b, a)) if i == c})
    F = nl.reference_load(m, nl._find_combination(m, "LC2"))
    lam, _, b = nl.linear_buckling_reference(m, F)
    phi = b.modes[:, 0].copy()
    phi /= np.max(np.abs(phi))
    p3 = phi.reshape(-1, 3)
    imax = int(np.argmax(np.abs(p3).max(axis=1)))
    sq = np.sum(p3 ** 2, axis=1)
    share6 = float((sq[c] + sq[nb].sum()) / sq.sum())
    rise0 = float(P[c, 2] - P[nb, 2].mean())
    drop_per_a = float(-(p3[c, 2] - p3[nb, 2].mean()))
    span = float(g.base_diameter)
    row = dict(nu=nu, lam_cr=float(lam), crown=c, node_of_max=imax, crown_phi=p3[c].tolist(),
               nb_z_mean=float(p3[nb, 2].mean()), share_crown_and_nb=share6, rise0=rise0,
               rise_loss_per_a=drop_per_a,
               a_flat=rise0 / drop_per_a if drop_per_a > 0 else None,
               span_over_a_flat=span / (rise0 / drop_per_a) if drop_per_a > 0 else None)
    for d in (1000, 500, 300):
        a = span / d
        row[f"rise_left_L{d}"] = rise0 - drop_per_a * a
    OUT.append(row)
    print(row, flush=True)
json.dump(OUT, open(OUTFILE, "w"), indent=1)
