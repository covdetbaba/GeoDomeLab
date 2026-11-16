"""Supplementary run for Chapter 6 of the GeoDomeLab 2.0 final report.

    python examples/report_supplementary/ring_and_rise.py [out.json]

Uses the geodome package of this repository (two folders up). Writes its
results to a JSON file next to the script unless a path is given."""
"""Chapter 6 supplementary run: (1) crown rise over its five neighbours for 2V to 8V
hemispheres (R = 15 m) against span/1000; (2) base-ring force under +30 C on the sized
6V design for horizontal support stiffness from rigid down to the ring's own radial stiffness.
Does not write into the repo."""
import sys, copy, json, math
import numpy as np
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else os.path.splitext(os.path.abspath(__file__))[0] + ".json"
from geodome import benchmarks as bm, geometry as geo
from geodome.postprocess import run_analysis
from geodome.optimize import auto_size
from geodome.model import LoadCombination

OUT = {}
# (1) crown rise
rows = []
for nu in (2, 4, 6, 8):
    g = geo.generate_dome(nu, 15.0, 0.5)
    P = g.nodes
    c = int(g.crown_node)
    nb = sorted({int(j) for a, b in g.members for (i, j) in ((a, b), (b, a)) if i == c})
    L = [float(np.linalg.norm(P[c] - P[j])) for j in nb]
    rise = float(P[c, 2] - P[nb, 2].mean())
    span = float(g.base_diameter)
    rows.append(dict(nu=nu, n_nb=len(nb), L_crown=float(np.mean(L)), rise=rise,
                     L2_over_2R=float(np.mean(L)) ** 2 / 30.0, span=span,
                     span1000_over_rise=span / 1000 / rise, span300_over_rise=span / 300 / rise,
                     Lmean=float(np.mean(g.member_lengths()))))
OUT["rise"] = rows
for r in rows:
    print(r)

# (2) ring force vs horizontal support stiffness, sized 6V, +30 C, 1.0T alone
sized = bm.make_preset("case_study_6v")
auto_size(sized, target_util=0.95, max_iter=12, enforce_deflection=True, slenderness_limit=200.0)
g = sized.geom
base = set(int(i) for i in sized.support_nodes())
ring = np.array([int(a) in base and int(b) in base for a, b in g.members])
EA = np.asarray(sized.EA())
Lm = g.member_lengths()
kr = EA[ring] / Lm[ring]
n_ring = int(ring.sum())
k_radial = float(np.mean(kr) * 4 * math.sin(math.pi / n_ring) ** 2)
OUT["ring"] = dict(n=n_ring, EA_over_L_mean=float(np.mean(kr)), k_radial_per_node=k_radial)
print(OUT["ring"])
res = []
for kh in (None, 5e7, 1e7, 5e6, k_radial, 1e6, 1e5):
    m = copy.deepcopy(sized)
    if kh is not None:
        m.support_type = geo.ELASTIC
        m.spring_stiffness = (kh, kh, 1e12)
    m.loads.temperature_change = 30.0
    m.loads.settlements = []
    m.combinations = [LoadCombination("T 1.0T", T=1.0)]
    r = run_analysis(m, n_buckling_modes=1, do_buckling=False)
    c = r.cases[r.order[0]]
    N = np.asarray(c.N) / 1e3
    U = np.asarray(r.util[r.order[0]])
    row = dict(kh=kh, ring_N_min=float(N[ring].min()), ring_N_max=float(N[ring].max()),
               U_ring=float(U[ring].max()), U_other=float(U[~ring].max()),
               relief=1 - float(N[ring].min()) / (-40.297033038407676))
    res.append(row)
    print(row, flush=True)
OUT["ring_vs_k"] = res
json.dump(OUT, open(OUTFILE, "w"), indent=1, default=str)
