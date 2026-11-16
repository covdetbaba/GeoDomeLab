"""
Run GeoDomeLab on the two hand-calculated validation models and compare.

The program side uses only the normal GeoDomeLab API, unchanged:
    geometry.generate_dome  ->  model.DomeModel  ->  postprocess.run_analysis
    (or solver.solve_model for the point-load case)
The hand side is validation/hand_calc.py, which never imports geodome.

Outputs (validation/output/):
    comparison.csv          every compared quantity, hand vs program
    comparison.md           the same as a markdown table
    model_A_sketch.png      node / member / DOF numbering, plan + elevation
    model_B_sketch.png
    program_A1_axial.png    GeoDomeLab's own 3D result plots
    program_B1_axial.png
    program_B2_deflected.png
    model_A_LC2.json        the validation models, ready for File > Open
    model_B_LC1.json
    model_B_LC2.json

A row passes when |program - hand| / |hand| <= REL_TOL (1e-12, that is
1e-10 %), or, where the hand answer is zero, when |program| <= ABS_TOL
(1e-9 in the row unit).  The tables carry a pass column, the script prints
"N of M pass" and exits with code 1 when any row fails.

Run:  python validation/run_validation.py
"""

from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
from matplotlib.figure import Figure                           # noqa: E402

import hand_calc as hc                                         # noqa: E402
from geodome import geometry as geo                            # noqa: E402
from geodome import plotting as plot                           # noqa: E402
from geodome.model import DomeModel, LoadSpec, LoadCombination  # noqa: E402
from geodome.postprocess import run_analysis                   # noqa: E402
from geodome.solver import solve_model, assemble_stiffness     # noqa: E402
from geodome.verification import run_all                       # noqa: E402

OUT = os.path.join(HERE, "output")
SECTION = "CHS 88.9 x 4.0"
MATERIAL = "ASTM A500 Gr. C (round)"

FRACTION_A = 0.5 * (1.0 - 1.0 / math.sqrt(5.0))   # cut at the upper ring
FRACTION_B = 0.5 * (1.0 + 1.0 / math.sqrt(5.0))   # cut at the lower ring


# =============================================================================
# program models
# =============================================================================
def program_model(fraction, loads, combos):
    g = geo.generate_dome(1, hc.R, fraction)
    m = DomeModel(geom=g, material_name=MATERIAL,
                  group_sections=[SECTION] * g.n_strut_types(),
                  support_type=geo.PINNED, loads=loads, combinations=combos)
    m.K_effective = hc.K_EFF
    m.phi_c, m.phi_t = hc.PHI_C, hc.PHI_T
    return m


def lc2_loads():
    return LoadSpec(include_self_weight=True, dead_pressure=hc.Q_DEAD,
                    snow_pressure=hc.P_SNOW, snow_unbalanced=False,
                    live_pressure=0.0, wind_enabled=False,
                    crown_point_load=0.0)


def b1_loads():
    return LoadSpec(include_self_weight=False, dead_pressure=0.0,
                    snow_pressure=0.0, live_pressure=0.0, wind_enabled=False,
                    crown_point_load=-hc.PZ_B1)


def node_map(hand_nodes, prog_nodes):
    """hand node number -> program node index, matched by coordinates."""
    out = {}
    for h, xyz in hand_nodes.items():
        d = np.linalg.norm(prog_nodes - xyz, axis=1)
        k = int(np.argmin(d))
        if d[k] > 1e-9:
            raise RuntimeError(f"hand node {h} has no program twin ({d[k]:.2e} m)")
        out[h] = k
    return out


def member_map(hand_members, prog_members, nmap):
    lookup = {frozenset((int(i), int(j))): k
              for k, (i, j) in enumerate(prog_members)}
    return {m: lookup[frozenset((nmap[i], nmap[j]))]
            for m, (i, j) in hand_members.items()}


# =============================================================================
# comparison rows
# =============================================================================
ROWS = []

# Pass rule (also used by the GUI Validation menu and run_all_checks.py).
REL_TOL = 1e-12          # |p - h| / |h|                   (1e-10 %)
ABS_TOL = 1e-9           # |p| in the row unit when the hand answer is 0


def row(model, case, quantity, hand, prog, unit, scale=1.0):
    h = float(hand) * scale
    p = float(prog) * scale
    if abs(h) > 1e-9:
        diff = abs(p - h) / abs(h) * 100.0
        kind = "%"
        ok = diff <= REL_TOL * 100.0
    else:
        diff = abs(p - h)
        kind = unit
        ok = diff <= ABS_TOL
    ROWS.append({"model": model, "case": case, "quantity": quantity,
                 "hand": h, "program": p, "unit": unit,
                 "difference": diff, "diff_kind": kind,
                 "pass": bool(ok)})


# =============================================================================
# MODEL A
# =============================================================================
def validate_a(H):
    # ---------- A1: asymmetric crown point load (solver API) ---------------
    mA = program_model(FRACTION_A, lc2_loads(),
                       [LoadCombination("A2  1.2D + 1.6S", D=1.2, S=1.6)])
    g = mA.geom
    nmap = node_map(H["nodes"], g.nodes)
    mmap = member_map(H["members"], g.members, nmap)
    c = nmap[1]

    row("A", "geometry", "nodes", 6, g.n_nodes, "-")
    row("A", "geometry", "struts", 10, g.n_members, "-")
    row("A", "geometry", "supported nodes", 5, len(g.base_nodes), "-")
    row("A", "geometry", "strut length L", H["L"],
        float(np.max(g.strut_lengths)), "m")
    row("A", "geometry", "chord factor L/R", H["L"] / hc.R,
        g.strut_schedule()[0]["chord_factor"], "-")
    row("A", "geometry", "crown height", H["nodes"][1][2], g.height, "m")

    K = assemble_stiffness(g.nodes, g.members, mA.EA())
    K = K.toarray() if hasattr(K, "toarray") else K
    Kc = K[3 * c:3 * c + 3, 3 * c:3 * c + 3]
    for (i, j, lab) in ((0, 0, "Kff[1,1]"), (1, 1, "Kff[2,2]"),
                        (2, 2, "Kff[3,3]")):
        row("A", "stiffness", lab, H["Kff"][i, j], Kc[i, j], "kN/m", 1e-3)
    off = float(np.max(np.abs(Kc - np.diag(np.diag(Kc)))))
    row("A", "stiffness", "largest off-diagonal of Kff", 0.0, off, "kN/m", 1e-3)

    F = np.zeros(3 * g.n_nodes)
    F[3 * c + 0] = hc.PX_A1
    F[3 * c + 1] = hc.PY_A1
    F[3 * c + 2] = hc.PZ_A1
    r1 = solve_model(mA, F, name="A1")
    a1 = H["A1"]
    row("A", "A1", "crown ux (D1)", a1["D"][0], r1.u[3 * c], "mm", 1e3)
    row("A", "A1", "crown uy (D2)", a1["D"][1], r1.u[3 * c + 1], "mm", 1e3)
    row("A", "A1", "crown uz (D3)", a1["D"][2], r1.u[3 * c + 2], "mm", 1e3)
    for m in range(1, 6):
        row("A", "A1", f"N member {m}", a1["N"][m], r1.N[mmap[m]], "kN", 1e-3)
    ring_max = max(abs(r1.N[mmap[m]]) for m in range(6, 11))
    row("A", "A1", "N ring members 6-10 (max |N|)", 0.0, ring_max, "kN", 1e-3)
    Rp = r1.reaction_vectors()
    for nd in (2, 4):
        for ax, lab in ((0, "Rx"), (2, "Rz")):
            row("A", "A1", f"{lab} at node {nd}", a1["R"][nd][ax],
                Rp[nmap[nd], ax], "kN", 1e-3)
    row("A", "A1", "sum Rx", a1["Rsum"][0], Rp.sum(axis=0)[0], "kN", 1e-3)
    row("A", "A1", "sum Rz", a1["Rsum"][2], Rp.sum(axis=0)[2], "kN", 1e-3)

    # ---------- A2: 1.2D + 1.6S through the full load generator ------------
    res = run_analysis(mA, n_buckling_modes=6)
    r2 = res.cases["A2  1.2D + 1.6S"]
    a2 = H["A2"]
    row("A", "A2", "crown load Pz", a2["P"][2], r2.F[3 * c + 2], "kN", 1e-3)
    row("A", "A2", "ring node load Pz", a2["P_ring"],
        r2.F[3 * nmap[2] + 2], "kN", 1e-3)
    row("A", "A2", "total vertical load", a2["tot_P"],
        r2.F[2::3].sum(), "kN", 1e-3)
    row("A", "A2", "crown uz", a2["D"][2], r2.u[3 * c + 2], "mm", 1e3)
    row("A", "A2", "N spoke (all 5)", a2["N"][1],
        float(np.mean([r2.N[mmap[m]] for m in range(1, 6)])), "kN", 1e-3)
    row("A", "A2", "Rz per support node", a2["R"][2][2],
        r2.reaction_vectors()[nmap[2], 2], "kN", 1e-3)

    ai = H["aisc"]
    k1 = mmap[1]
    row("A", "AISC E3", "KL/r", ai["KLr"], res.caps["slenderness"][k1], "-")
    row("A", "AISC E3", "Fcr", ai["Fcr"], res.caps["Fcr"][k1], "MPa", 1e-6)
    row("A", "AISC E3", "phi_c Pn", ai["phiPn_c"], res.caps["phiPn_c"][k1],
        "kN", 1e-3)
    row("A", "AISC D2", "phi_t Pn", ai["phiPn_t"], res.caps["phiPn_t"][k1],
        "kN", 1e-3)
    row("A", "AISC E3", "utilisation", ai["util"], res.util_env[k1], "-")
    row("A", "buckling", "lambda_cr", H["buckling"]["lam_cr"],
        res.lambda_cr_by_case["A2  1.2D + 1.6S"], "-")
    return mA, res, r1


# =============================================================================
# MODEL B
# =============================================================================
def validate_b(H):
    mB1 = program_model(FRACTION_B, b1_loads(),
                        [LoadCombination("B1  1.0L", L=1.0)])
    mB2 = program_model(FRACTION_B, lc2_loads(),
                        [LoadCombination("B2  1.2D + 1.6S", D=1.2, S=1.6)])
    g = mB1.geom
    nmap = node_map(H["nodes"], g.nodes)
    mmap = member_map(H["members"], g.members, nmap)
    c = nmap[1]

    row("B", "geometry", "nodes", 11, g.n_nodes, "-")
    row("B", "geometry", "struts", 25, g.n_members, "-")
    row("B", "geometry", "supported nodes", 5, len(g.base_nodes), "-")
    row("B", "geometry", "crown height", H["nodes"][1][2], g.height, "m")

    # Projection of the program's 18 x 18 [Kff] onto the symmetric pattern
    # must reproduce the virtual-work matrix [Ks].
    K = assemble_stiffness(g.nodes, g.members, mB1.EA())
    K = K.toarray() if hasattr(K, "toarray") else K
    T = np.zeros((3 * g.n_nodes, 3))
    T[3 * c + 2, 0] = 1.0                                   # w0
    for i in range(5):
        nd = nmap[2 + i]
        th = 72.0 * i * hc.DEG
        T[3 * nd + 0, 1] = math.cos(th)                     # u_r
        T[3 * nd + 1, 1] = math.sin(th)
        T[3 * nd + 2, 2] = 1.0                              # w1
    Kproj = T.T @ K @ T
    for (i, j, lab) in ((0, 0, "Ks[1,1]"), (1, 1, "Ks[2,2]"),
                        (2, 2, "Ks[3,3]"), (0, 1, "Ks[1,2]"),
                        (1, 2, "Ks[2,3]"), (0, 2, "Ks[1,3]")):
        row("B", "stiffness", lab + " (T^T K T)", H["Ks"][i, j], Kproj[i, j],
            "kN/m", 1e-3)

    results = {}
    for name, model in (("B1", mB1), ("B2", mB2)):
        res = run_analysis(model, do_buckling=False)
        r = res.cases[model.combinations[0].name]
        hb = H[name]
        if name == "B2":
            row("B", name, "crown load Pz", hb["Pc"], r.F[3 * c + 2], "kN", 1e-3)
            row("B", name, "upper node load Pz", hb["Pu"],
                r.F[3 * nmap[2] + 2], "kN", 1e-3)
            row("B", name, "lower node load Pz", hb["Pl"],
                r.F[3 * nmap[7] + 2], "kN", 1e-3)
        row("B", name, "total vertical load", hb["tot_P"], r.F[2::3].sum(),
            "kN", 1e-3)
        row("B", name, "crown uz (w0)", hb["w0"], r.u[3 * c + 2], "mm", 1e3)
        urs, uts, w1s = [], [], []
        for i in range(5):
            nd = nmap[2 + i]
            th = 72.0 * i * hc.DEG
            ux, uy, uz = r.u[3 * nd:3 * nd + 3]
            urs.append(ux * math.cos(th) + uy * math.sin(th))
            uts.append(-ux * math.sin(th) + uy * math.cos(th))
            w1s.append(uz)
        row("B", name, "upper node radial (u_r)", hb["ur"], np.mean(urs),
            "mm", 1e3)
        row("B", name, "upper node vertical (w1)", hb["w1"], np.mean(w1s),
            "mm", 1e3)
        row("B", name, "upper node tangential (max |u_t|)", 0.0,
            max(abs(x) for x in uts), "mm", 1e3)
        for lab, key, ms in (("N spokes 1-5", "N_spoke", range(1, 6)),
                             ("N upper ring 6-10 (hoop)", "N_ring", range(6, 11)),
                             ("N zigzag 11-20", "N_zig", range(11, 21))):
            vals = [r.N[mmap[m]] for m in ms]
            row("B", name, lab, hb[key], np.mean(vals), "kN", 1e-3)
            spread = max(vals) - min(vals)
            row("B", name, lab + " spread (max-min)", 0.0, spread, "kN", 1e-3)
        low = max(abs(r.N[mmap[m]]) for m in range(21, 26))
        row("B", name, "N lower ring 21-25 (max |N|)", 0.0, low, "kN", 1e-3)
        Rv = r.reaction_vectors()
        nd = nmap[7]
        beta = 36.0 * hc.DEG
        row("B", name, "Rz per support node", hb["Rz_node"], Rv[nd, 2],
            "kN", 1e-3)
        row("B", name, "R radial per support node", hb["Rr_node"],
            Rv[nd, 0] * math.cos(beta) + Rv[nd, 1] * math.sin(beta),
            "kN", 1e-3)
        row("B", name, "sum Rz", -hb["tot_P"], Rv.sum(axis=0)[2], "kN", 1e-3)
        results[name] = (res, r)
    return mB1, mB2, results


# =============================================================================
# figures
# =============================================================================
def sketch(H, path, title, free_nodes, dof_note, dof_arrows=False):
    nodes, members = H["nodes"], H["members"]
    fig = Figure(figsize=(12.0, 5.4), dpi=150)
    fig.patch.set_facecolor("white")
    axp = fig.add_subplot(1, 2, 1)
    axe = fig.add_subplot(1, 2, 2)
    for ax, (i0, i1), name in ((axp, (0, 1), "Plan (x-y)"),
                               (axe, (0, 2), "Elevation (x-z)")):
        for m, (i, j) in members.items():
            p, q = nodes[i], nodes[j]
            ax.plot([p[i0], q[i0]], [p[i1], q[i1]], color="#2b4c6f", lw=1.6,
                    zorder=1)
            if ax is axp:
                mid = 0.5 * (p + q)
                ax.text(mid[i0], mid[i1], str(m), fontsize=8, color="#b4161b",
                        ha="center", va="center", zorder=4,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                  ec="none", alpha=0.85))
        for n, xyz in nodes.items():
            sup = n not in free_nodes
            ax.scatter(xyz[i0], xyz[i1], s=70 if sup else 55,
                       marker="^" if sup else "o",
                       c="#111111" if sup else "#1f6fb4", zorder=3)
            if ax is axp:
                ax.text(xyz[i0] + 0.18, xyz[i1] + 0.18, str(n), fontsize=10,
                        fontweight="bold", color="#111111", zorder=5)
        ax.set_aspect("equal")
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.grid(alpha=0.25)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]" if ax is axp else "z [m]")
    # elevation: label nodes on the x-z silhouette that are visible
    for n, xyz in nodes.items():
        if abs(xyz[1]) < 1e-6 or n == 1:
            axe.text(xyz[0] + 0.15, xyz[2] + 0.15, str(n), fontsize=10,
                     fontweight="bold")
    # DOF arrows at the crown (free DOF numbered first, as in the course)
    cz = nodes[1]
    arr = dict(arrowstyle="-|>", color="#e07b00", lw=1.8, mutation_scale=14)
    if dof_arrows:
        axp.annotate("", xy=(cz[0] + 1.3, cz[1]), xytext=(cz[0], cz[1]),
                     arrowprops=arr, zorder=6)
        axp.annotate("", xy=(cz[0], cz[1] + 1.3), xytext=(cz[0], cz[1]),
                     arrowprops=arr, zorder=6)
        axp.text(cz[0] + 1.35, cz[1] - 0.45, "D1", color="#e07b00",
                 fontsize=10, fontweight="bold")
        axp.text(cz[0] - 0.55, cz[1] + 1.35, "D2", color="#e07b00",
                 fontsize=10, fontweight="bold")
        axe.annotate("", xy=(cz[0] + 1.3, cz[2]), xytext=(cz[0], cz[2]),
                     arrowprops=arr, zorder=6)
        axe.annotate("", xy=(cz[0], cz[2] + 1.3), xytext=(cz[0], cz[2]),
                     arrowprops=arr, zorder=6)
        axe.text(cz[0] + 1.35, cz[2] - 0.4, "D1", color="#e07b00",
                 fontsize=10, fontweight="bold")
        axe.text(cz[0] + 0.15, cz[2] + 1.35, "D3", color="#e07b00",
                 fontsize=10, fontweight="bold")
    # headroom so the note never sits on a member
    zmax = max(v[2] for v in nodes.values())
    axe.set_ylim(-0.6, zmax + (3.6 if dof_arrows else 3.2))
    axe.text(0.02, 0.98, dof_note, transform=axe.transAxes, fontsize=9,
             va="top", family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", fc="#f5f7f9", ec="#c8d2db"))
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.text(0.5, 0.01, "Black numbers = nodes, red = members, triangles = "
             "pinned supports, blue dots = free nodes.", ha="center",
             fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(path, facecolor="white")


def program_view(draw, model, *args, path, cbar=None, **kw):
    fig = Figure(figsize=(7.4, 5.8), dpi=150)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection="3d")
    art = draw(ax, model, *args, **kw)
    if cbar and hasattr(art, "set_array"):
        cb = fig.colorbar(art, ax=ax, shrink=0.6, pad=0.10)
        cb.set_label(cbar, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, facecolor="white")


# =============================================================================
def main():
    os.makedirs(OUT, exist_ok=True)
    HA = hc.model_a()
    HB = hc.model_b()
    ROWS.clear()

    mA, resA, rA1 = validate_a(HA)
    mB1, mB2, resB = validate_b(HB)

    # ---- tables ------------------------------------------------------------
    with open(os.path.join(OUT, "comparison.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ROWS[0].keys()))
        w.writeheader()
        w.writerows(ROWS)
    with open(os.path.join(OUT, "comparison.md"), "w", encoding="utf-8") as fh:
        fh.write("| Model | Case | Quantity | Hand | GeoDomeLab | Unit | "
                 "Difference | Pass |\n|---|---|---|---|---|---|---|---|\n")
        for r in ROWS:
            d = (f"{r['difference']:.1e} %" if r["diff_kind"] == "%"
                 else f"{r['difference']:.1e} {r['unit']}")
            fh.write(f"| {r['model']} | {r['case']} | {r['quantity']} | "
                     f"{r['hand']:.6g} | {r['program']:.6g} | {r['unit']} | "
                     f"{d} | {'yes' if r['pass'] else 'NO'} |\n")
        n_ok = sum(r["pass"] for r in ROWS)
        fh.write(f"\n{n_ok} of {len(ROWS)} pass (|p - h| / |h| <= "
                 f"{REL_TOL:g}, or |p| <= {ABS_TOL:g} in the row unit where "
                 f"the hand answer is zero).\n")

    # ---- figures ------------------------------------------------------------
    sketch(HA, os.path.join(OUT, "model_A_sketch.png"),
           "Model A - 1V pentagonal cap (6 nodes, 10 struts, 3 free DOF)",
           free_nodes={1},
           dof_note="Free DOF (numbered first):\n D1 = ux(1)\n D2 = uy(1)\n"
                    " D3 = uz(1)\nRestrained: D4-D18\n(nodes 2-6, pinned)",
           dof_arrows=True)
    sketch(HB, os.path.join(OUT, "model_B_sketch.png"),
           "Model B - 1V dome on its lower ring (11 nodes, 25 struts, 18 "
           "free DOF)",
           free_nodes={1, 2, 3, 4, 5, 6},
           dof_note="Free: nodes 1-6 (18 DOF)\nPinned: nodes 7-11\n"
                    "Symmetry -> 3 unknowns:\n w0  crown vertical\n"
                    " u_r ring radial\n w1  ring vertical")
    program_view(plot.draw_axial, mA, rA1.N,
                 path=os.path.join(OUT, "program_A1_axial.png"),
                 title="GeoDomeLab - Model A, case A1 axial force",
                 cbar="axial force [kN]  (+ tension)")
    resB1, rB1 = resB["B1"]
    program_view(plot.draw_axial, mB1, rB1.N,
                 path=os.path.join(OUT, "program_B1_axial.png"),
                 title="GeoDomeLab - Model B, case B1 axial force",
                 cbar="axial force [kN]  (+ tension)")
    resB2, rB2 = resB["B2"]
    program_view(plot.draw_deflected, mB2, rB2,
                 path=os.path.join(OUT, "program_B2_deflected.png"),
                 title="GeoDomeLab - Model B, case B2 deflected shape",
                 cbar="nodal displacement [mm]")

    # ---- GUI-ready model files ----------------------------------------------
    mA.save(os.path.join(OUT, "model_A_LC2.json"))
    mB1.save(os.path.join(OUT, "model_B_LC1.json"))
    mB2.save(os.path.join(OUT, "model_B_LC2.json"))

    # ---- print ------------------------------------------------------------
    print(hc.report())
    print()
    print("=" * 104)
    print(f"{'model':<6}{'case':<10}{'quantity':<40}{'hand':>14}"
          f"{'program':>14}{'unit':>7}{'diff':>13}")
    print("-" * 104)
    worst_pct = 0.0
    worst_abs = 0.0
    for r in ROWS:
        if r["diff_kind"] == "%":
            worst_pct = max(worst_pct, r["difference"])
            d = f"{r['difference']:.1e} %"
        else:
            worst_abs = max(worst_abs, r["difference"])
            d = f"{r['difference']:.1e}"
        print(f"{r['model']:<6}{r['case']:<10}{r['quantity']:<40}"
              f"{r['hand']:>14.6g}{r['program']:>14.6g}{r['unit']:>7}{d:>13}")
    print("-" * 104)
    print(f"{len(ROWS)} quantities compared.  Largest relative difference "
          f"{worst_pct:.2e} %.  Largest absolute value where the hand answer "
          f"is zero {worst_abs:.2e}.")
    n_ok = sum(r["pass"] for r in ROWS)
    print(f"{n_ok} of {len(ROWS)} pass (relative tolerance {REL_TOL:g}, "
          f"absolute {ABS_TOL:g} where the hand answer is zero).")
    for r in ROWS:
        if not r["pass"]:
            print(f"  FAIL  {r['model']} {r['case']} {r['quantity']}")
    suite = run_all()
    print(f"Built-in verification suite on the same code: "
          f"{sum(r['pass'] for r in suite)} of {len(suite)} checks pass.")
    return 0 if n_ok == len(ROWS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
