"""Supplementary run for Chapter 6 of the GeoDomeLab 2.0 final report.

    python examples/report_supplementary/hub_tolerance.py [out.json]

Uses the geodome package of this repository (two folders up). Writes its
results to a JSON file next to the script unless a path is given."""
"""Chapter 6 supplementary run: imperfection sensitivity of the sized 6V dome at
hub-scale amplitudes (mode-affine and crown hub alone). Writes only its own (git-ignored) JSON file."""
import sys, json, time
import numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else os.path.splitext(os.path.abspath(__file__))[0] + ".json"
import geodome
from geodome import benchmarks as bm
from geodome import nonlinear as nl
from geodome.optimize import auto_size

OUT = {"geodome": geodome.__version__}
t0 = time.perf_counter()
model = bm.make_preset("case_study_6v")
sz = auto_size(model, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
OUT["sized_mass"] = float(model.steel_mass())
g = model.geom
P = g.nodes
crown = int(g.crown_node)
nb = sorted({int(j) for a, b in g.members for (i, j) in ((a, b), (b, a)) if i == crown})
rise0 = float(P[crown, 2] - P[nb, 2].mean())
span = float(g.base_diameter)
strut = float(np.mean(g.member_lengths()))
combo = nl._find_combination(model, "LC2")
F = nl.reference_load(model, combo)
lam_cr, _, b = nl.linear_buckling_reference(model, F)
phi = b.modes[:, 0].copy()
phi /= np.max(np.abs(phi))
phi3 = phi.reshape(-1, 3)
OUT.update(crown=crown, neighbours=nb, rise0=rise0, span=span, strut_mean=strut,
           lam_cr=float(lam_cr), phi_crown=phi3[crown].tolist(),
           phi_nb_z=phi3[nb, 2].tolist())

# 1) mode-affine at hub-scale amplitudes
amps = [0.003, 0.005, 0.010, 0.015, 0.020]
rows = nl.imperfection_study(model, "LC2", amplitudes=amps)
for r in rows:
    a = r["amplitude_m"]
    Q = P + a * phi3
    r["crown_rise_plus"] = float(Q[crown, 2] - Q[nb, 2].mean())
OUT["mode_affine"] = rows

# 2) crown hub alone moved down (and up) by a
single = []
for a in [0.005, 0.010, 0.020, 0.030, 0.060]:
    row = dict(a=a)
    for sgn in (1.0, -1.0):
        off = np.zeros_like(P)
        off[crown, 2] = -sgn * a          # +1: crown pushed DOWN
        r = nl.run_nonlinear(model, combo, node_offset=off, load_vector=F, lambda_cr=lam_cr,
                             compare=False, refine_bifurcation=True, member_checks=False)
        row["lam_down" if sgn > 0 else "lam_up"] = float(r.lambda_limit) if r.limit_found else None
        row["type_down" if sgn > 0 else "type_up"] = r.limit_type
    row["crown_rise_down"] = rise0 - a
    single.append(row)
    print(row, flush=True)
OUT["crown_only"] = single
OUT["time_s"] = time.perf_counter() - t0


def conv(o):
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (np.bool_,)): return bool(o)
    return str(o)


json.dump(OUT, open(OUTFILE, "w", encoding="utf-8"), indent=1, default=conv)
print(json.dumps(OUT["mode_affine"], indent=1, default=conv))
print("time", OUT["time_s"])
