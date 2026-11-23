"""
GeoDomeLab in the browser
=========================

A small Streamlit front end for the GeoDomeLab core.  It imports the same
``geodome`` package as the desktop program (no copy of the analysis code and
no Qt), so every number here is the number the desktop program gives.

Run it locally from the repository root:

    streamlit run web/streamlit_app.py

The same file runs on Streamlit Community Cloud (web/requirements.txt) and in
the browser through stlite and Pyodide (web/build_pages.py builds that site).

What it does
------------
* Two one-click examples: the 6V case-study dome (30 m span) and validation
  Model A, checked against the independent hand calculation.
* A custom dome: frequency, radius, dome fraction, supports, section, loads.
* Outputs: key numbers with units, an interactive 3D view coloured by axial
  force, member and node tables, CSV downloads and a verification table.
"""

from __future__ import annotations

import copy
import io
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for _p in (ROOT, ROOT / "validation"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import geodome                                              # noqa: E402
from geodome import benchmarks as bm                        # noqa: E402
from geodome import geometry as geo                         # noqa: E402
from geodome import optimize as opt                         # noqa: E402
from geodome import postprocess as pp                       # noqa: E402
from geodome.model import DomeModel, LoadSpec               # noqa: E402
from geodome.sections import SECTION_NAMES, DEFAULT_SECTION  # noqa: E402

# Python plotly is optional.  Without it (the stlite build leaves it out to
# keep the first load short) the same figure is drawn by plotly.js inside a
# components.html frame.  GEODOMELAB_NO_PLOTLY=1 forces that path for tests.
try:
    if os.environ.get("GEODOMELAB_NO_PLOTLY"):
        raise ImportError("plotly switched off by GEODOMELAB_NO_PLOTLY")
    import plotly.graph_objects as go                       # noqa: E402
    HAVE_PLOTLY = True
except ImportError:
    go = None
    HAVE_PLOTLY = False

PLOTLY_JS = "https://cdn.plot.ly/plotly-gl3d-2.35.2.min.js"   # pinned
REPO_URL = geodome.__url__
IN_BROWSER = sys.platform == "emscripten"                  # stlite / Pyodide

# matplotlib "coolwarm", the colour map of the desktop axial-force view:
# blue = compression, red = tension, symmetric about zero.
COOLWARM = [[0.0, "#3b4cc0"], [0.125, "#6282ea"], [0.25, "#8db0fe"],
            [0.375, "#b9d0f9"], [0.5, "#dddcdc"], [0.625, "#f5c4ac"],
            [0.75, "#f4987a"], [0.875, "#dd5f4b"], [1.0, "#b40426"]]

SUPPORT_CHOICES = {
    geo.PINNED: "Pinned base ring",
    geo.ROLLER3: "Vertical rollers + 3 lateral restraints",
    geo.ELASTIC: "Elastic springs, 50 kN/mm each way",
}

EXAMPLE_6V = "Example: 6V dome, 30 m span"
EXAMPLE_A = "Verification: Model A against the hand calculation"


# =============================================================================
# analysis (thin wrappers around the geodome API)
# =============================================================================
def solve(model, label, kind, sizing=None, start_mass=None):
    """Run every load combination of ``model`` and time it."""
    t0 = time.perf_counter()
    res = pp.run_analysis(model, n_buckling_modes=6)
    return {"model": model, "res": res, "label": label, "kind": kind,
            "time_s": time.perf_counter() - t0, "sizing": sizing,
            "start_mass": start_mass}


def custom_model(freq, radius, fraction, support, section, dead_kpa,
                 snow_kpa, wind_ms):
    """A geodesic dome built exactly as the desktop Geometry panel builds it."""
    g = geo.generate_dome(int(freq), float(radius), float(fraction))
    m = DomeModel(geom=g, support_type=support)
    m.group_sections = [section] * g.n_strut_types()
    m.loads = LoadSpec(dead_pressure=1e3 * float(dead_kpa),
                       snow_pressure=1e3 * float(snow_kpa),
                       wind_enabled=float(wind_ms) > 0.0,
                       wind_speed=float(wind_ms) if wind_ms > 0 else 45.0)
    return m


def store(run, combo):
    """Keep ``run`` and ask the next "Combination to show" box for ``combo``.

    The box reads the request before it is drawn, so this also works from a
    button placed below the box (Streamlit forbids changing a widget value
    after the widget is drawn in the same run).
    """
    st.session_state["run"] = run
    order = run["res"].order
    st.session_state["pending_combo"] = combo if combo in order else order[0]


# =============================================================================
# formatting helpers
# =============================================================================
def lam_text(value):
    return f"{value:.2f}" if math.isfinite(value) else "none"


def short(name):
    """'LC2  1.2D + 1.6S' -> 'LC2 1.2D + 1.6S' for labels."""
    return " ".join(str(name).split())


def combo_id(name):
    """'LC2  1.2D + 1.6S' -> 'LC2' (the label part of a combination name)."""
    parts = str(name).split()
    return parts[0] if parts else str(name)


def secs(t):
    return f"{t:.2f} s" if t >= 0.1 else f"{t * 1e3:.0f} ms"


def g10(x):
    return f"{x:.10g}"


def extreme(values, mode):
    """``(k, members, value)`` of the largest (mode "max") or smallest value.

    A geodesic dome is symmetric, so several struts often share the extreme
    value to round-off.  ``members`` lists all of them (relative tolerance
    1e-9) and ``k`` is the lowest-numbered one, so the member reported does
    not depend on the round-off of the platform.
    """
    v = np.asarray(values, dtype=float)
    ext = float(v.max() if mode == "max" else v.min())
    tol = 1e-9 * max(abs(ext), 1e-30)
    members = [int(i) for i in np.flatnonzero(np.abs(v - ext) <= tol)]
    return members[0], members, ext


def member_label(text, members):
    extra = len(members) - 1
    return (f"{text}, member {members[0]}"
            + (f" (+{extra} equal)" if extra else ""))


def twins_help(members, what):
    if len(members) < 2:
        return None
    return (f"Members {', '.join(str(k) for k in members)} carry the same "
            f"{what} to round-off (symmetry of the dome). The lowest number "
            f"is shown.")


def to_csv(df):
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    return buf.getvalue().encode("utf-8")


def file_stub(run):
    g = run["model"].geom
    tag = "modelA" if run["kind"] == "modelA" else f"{g.frequency}V"
    return f"geodomelab_{tag}"


# =============================================================================
# 3D view
# =============================================================================
def dome_figure(model, res, combo, height=500):
    """Plain plotly figure dict: struts coloured by axial force [kN]."""
    g = model.geom
    X = np.asarray(g.nodes, dtype=float)
    M = np.asarray(g.members, dtype=int)
    r = res.cases[combo]
    N = np.asarray(r.N, dtype=float) / 1e3
    amax = max(float(np.max(np.abs(N))), 1e-9)

    xs, ys, zs, cs, tx = [], [], [], [], []
    for k, (i, j) in enumerate(M):
        fam = chr(ord("A") + int(g.strut_group[k]))
        kind = "tension" if N[k] > 1e-9 else (
            "compression" if N[k] < -1e-9 else "zero")
        t = (f"Member {k} ({fam}), nodes {i}-{j}<br>"
             f"N = {N[k]:+.2f} kN ({kind})")
        xs += [float(X[i, 0]), float(X[j, 0]), None]
        ys += [float(X[i, 1]), float(X[j, 1]), None]
        zs += [float(X[i, 2]), float(X[j, 2]), None]
        cs += [float(N[k])] * 3
        tx += [t, t, ""]

    members = {
        "type": "scatter3d", "mode": "lines", "x": xs, "y": ys, "z": zs,
        "hovertext": tx, "hoverinfo": "text", "connectgaps": False,
        "line": {"color": cs, "colorscale": COOLWARM, "cmin": -amax,
                 "cmax": amax, "width": 5, "showscale": True,
                 "colorbar": {"title": {"text": "N (kN)", "side": "top"},
                              "orientation": "h", "thickness": 12,
                              "len": 0.8, "x": 0.5, "xanchor": "center",
                              "y": 0.0, "yanchor": "top"}},
        "name": "struts",
    }

    sup = set(int(v) for v in model.support_nodes())
    d = r.displacement_vectors() * 1e3
    free = [k for k in range(g.n_nodes) if k not in sup]
    fixed = sorted(sup)

    def node_trace(ids, size, colour, symbol, name):
        return {
            "type": "scatter3d", "mode": "markers", "name": name,
            "x": [float(X[k, 0]) for k in ids],
            "y": [float(X[k, 1]) for k in ids],
            "z": [float(X[k, 2]) for k in ids],
            "hovertext": [
                f"Node {k}{' (support)' if k in sup else ''}<br>"
                f"x, y, z = {X[k, 0]:.3f}, {X[k, 1]:.3f}, {X[k, 2]:.3f} m<br>"
                f"ux, uy, uz = {d[k, 0]:.3f}, {d[k, 1]:.3f}, {d[k, 2]:.3f} mm"
                for k in ids],
            "hoverinfo": "text",
            "marker": {"size": size, "color": colour, "symbol": symbol},
        }

    data = [members, node_trace(free, 2.5, "#5b6470", "circle", "nodes")]
    if fixed:
        data.append(node_trace(fixed, 4.5, "#1f2a36", "square", "supports"))

    axis = {"showbackground": True, "backgroundcolor": "#f4f6f8",
            "gridcolor": "#d5dbe1", "zerolinecolor": "#b8c0c8",
            "color": "#3a4350"}
    layout = {
        "height": height,
        "margin": {"l": 0, "r": 0, "t": 10, "b": 10},
        "paper_bgcolor": "#ffffff",
        "font": {"family": "Arial, sans-serif", "size": 12,
                 "color": "#28313b"},
        "showlegend": False,
        "hoverlabel": {"bgcolor": "#ffffff"},
        "uirevision": f"{g.n_nodes}-{g.n_members}-{g.radius:.6g}",
        "scene": {
            "aspectmode": "data", "dragmode": "turntable",
            "camera": {"eye": {"x": 1.45, "y": -1.45, "z": 0.95},
                       "up": {"x": 0, "y": 0, "z": 1}},
            "xaxis": {**axis, "title": {"text": "x (m)"}},
            "yaxis": {**axis, "title": {"text": "y (m)"}},
            "zaxis": {**axis, "title": {"text": "z (m)"}},
        },
    }
    return {"data": data, "layout": layout}


PLOT_CONFIG = {"displaylogo": False, "scrollZoom": True, "responsive": True,
               "modeBarButtonsToRemove": ["toImage"]}


def show_figure(fig):
    if HAVE_PLOTLY:
        st.plotly_chart(go.Figure(fig), theme=None, config=PLOT_CONFIG,
                        width="stretch")
        return
    height = int(fig["layout"]["height"])
    html = (
        '<!doctype html><html><head><meta charset="utf-8">'
        '<style>html,body{margin:0;padding:0;background:#fff}</style>'
        f'<script src="{PLOTLY_JS}" charset="utf-8"></script></head><body>'
        f'<div id="g" style="width:100%;height:{height}px"></div><script>'
        f"var f={json.dumps(fig)};"
        # On a narrow screen (a phone held upright) move the camera back so
        # the whole dome fits the width; a wide frame keeps the set view.
        "var el=document.getElementById('g'),e=f.layout.scene.camera.eye,"
        "s=Math.max(1,1.1*el.clientHeight/Math.max(el.clientWidth,1));"
        "f.layout.scene.camera.eye={x:e.x*s,y:e.y*s,z:e.z*s};"
        f"Plotly.newPlot('g',f.data,f.layout,{json.dumps(PLOT_CONFIG)});"
        "</script></body></html>")
    components.html(html, height=height + 6)


# =============================================================================
# tables
# =============================================================================
def members_frame(run, combo):
    res = run["res"]
    strength = res.strength_cases()
    N = res.cases[combo].N / 1e3
    rows = []
    for row in pp.member_table(res):
        k = row["member"]
        gov = (short(strength[int(res.util_case[k])]) if strength else "")
        rows.append({
            "Member": k, "Node i": row["node_i"], "Node j": row["node_j"],
            "Family": row["family"], "Section": row["section"],
            "Length (m)": row["length_m"], "KL/r (-)": row["KL_r"],
            f"N, {short(combo)} (kN)": float(N[k]),
            "N max, strength envelope (kN)": row["N_max_kN"],
            "N min, strength envelope (kN)": row["N_min_kN"],
            "phi Pn tension (kN)": row["phiPn_t_kN"],
            "phi Pn compression (kN)": row["phiPn_c_kN"],
            "Utilisation (-)": row["utilisation"],
            "Governing combination": gov, "Status": row["status"],
        })
    return pd.DataFrame(rows)


def nodes_frame(run, combo):
    res = run["res"]
    R = res.cases[combo].reaction_vectors() / 1e3
    rows = []
    for row in pp.node_table(res, case=combo):
        k = row["node"]
        sup = row["support"] == "yes"
        rows.append({
            "Node": k, "x (m)": row["x_m"], "y (m)": row["y_m"],
            "z (m)": row["z_m"], "Ring": row["ring"],
            "ux (mm)": row["ux_mm"], "uy (mm)": row["uy_mm"],
            "uz (mm)": row["uz_mm"], "|u| (mm)": row["|u|_mm"],
            "Support": "yes" if sup else "",
            "Rx (kN)": float(R[k, 0]) if sup else None,
            "Ry (kN)": float(R[k, 1]) if sup else None,
            "Rz (kN)": float(R[k, 2]) if sup else None,
        })
    return pd.DataFrame(rows)


def number_config(df, digits):
    cfg = {}
    for col, nd in digits.items():
        if col in df.columns:
            cfg[col] = st.column_config.NumberColumn(col, format=f"%.{nd}f")
    return cfg


# =============================================================================
# verification: GeoDomeLab against validation/hand_calc.py (Model A)
# =============================================================================
REL_TOL = 1e-12          # the pass rule of validation/run_validation.py
ABS_TOL = 1e-9           # in the row unit, where the hand answer is zero


@st.cache_data(show_spinner=False)
def model_a_verification():
    """Rows Quantity | Location | Hand | GeoDomeLab | Error (%) for Model A.

    The hand side is validation/hand_calc.py, which never imports geodome.
    The program side is the Model A preset solved by run_analysis, the same
    call the desktop program makes.
    """
    import hand_calc as hc                                  # noqa: WPS433

    t0 = time.perf_counter()
    H = hc.model_a()
    m = bm.validation_model_a()
    res = pp.run_analysis(m, n_buckling_modes=6)
    g = m.geom

    # hand numbering (1-based) -> program numbering, matched by coordinates
    nmap = {}
    for h, xyz in H["nodes"].items():
        dist = np.linalg.norm(g.nodes - xyz, axis=1)
        nmap[h] = int(np.argmin(dist))
    lookup = {frozenset((int(i), int(j))): k
              for k, (i, j) in enumerate(g.members)}
    mmap = {mh: lookup[frozenset((nmap[i], nmap[j]))]
            for mh, (i, j) in H["members"].items()}
    c = nmap[1]

    rows = []

    def add(quantity, unit, location, hand, prog, scale):
        h = float(hand) * scale
        p = float(prog) * scale
        if abs(h) > 1e-9:
            err = abs(p - h) / abs(h) * 100.0
            ok = err <= REL_TOL * 100.0
            err_txt = "0" if err == 0.0 else f"{err:.1e}"
        else:
            err = float("nan")
            ok = abs(p - h) <= ABS_TOL
            diff = abs(p - h)
            err_txt = ("0 (reference 0)" if diff == 0.0
                       else f"reference 0, |diff| {diff:.0e} {unit}")
        rows.append({"Quantity": f"{quantity} ({unit})",
                     "Location": location, "Hand": g10(h),
                     "GeoDomeLab": g10(p), "Error (%)": err_txt,
                     "_err": err, "_pass": ok})

    a1 = H["A1"]
    r1 = res.cases[bm.COMBO_A1]
    crown = f"crown, node {c} (hand 1)"
    for ax, lab in enumerate("xyz"):
        add(f"A1 crown u{lab}", "mm", crown,
            a1["D"][ax], r1.u[3 * c + ax], 1e3)
    for mh in range(1, 6):
        add("A1 axial force N", "kN",
            f"spoke, member {mmap[mh]} (hand {mh})",
            a1["N"][mh], r1.N[mmap[mh]], 1e-3)
    ring = [mmap[mh] for mh in range(6, 11)]
    add("A1 ring struts, largest |N|", "kN",
        f"members {min(ring)} to {max(ring)} (hand 6 to 10)",
        0.0, max(abs(r1.N[k]) for k in ring), 1e-3)
    R1 = r1.reaction_vectors()
    for h in range(2, 7):
        for ax, lab in enumerate("xyz"):
            add(f"A1 reaction R{lab}", "kN",
                f"support, node {nmap[h]} (hand {h})",
                a1["R"][h][ax], R1[nmap[h], ax], 1e-3)
    tot1 = R1.sum(axis=0)
    for ax, lab in ((0, "x"), (2, "z")):
        add(f"A1 sum of reactions R{lab}", "kN", "all 5 supports",
            a1["Rsum"][ax], tot1[ax], 1e-3)

    a2 = H["A2"]
    r2 = res.cases[bm.COMBO_A2]
    add("A2 crown uz", "mm", crown, a2["D"][2],
        r2.u[3 * c + 2], 1e3)
    for mh in range(1, 6):
        add("A2 axial force N", "kN",
            f"spoke, member {mmap[mh]} (hand {mh})",
            a2["N"][mh], r2.N[mmap[mh]], 1e-3)
    R2 = r2.reaction_vectors()
    for h in range(2, 7):
        add("A2 reaction Rz", "kN", f"support, node {nmap[h]} (hand {h})",
            a2["R"][h][2], R2[nmap[h], 2], 1e-3)
    add("A2 total vertical load", "kN", "whole cap", a2["tot_P"],
        r2.F[2::3].sum(), 1e-3)
    add("A2 utilisation, AISC E3", "-",
        f"spoke, member {mmap[1]} (hand 1)",
        H["aisc"]["util"], res.util_env[mmap[1]], 1.0)
    add("A2 load factor λcr", "-", "whole cap",
        H["buckling"]["lam_cr"], res.lambda_cr_by_case[bm.COMBO_A2], 1.0)

    df = pd.DataFrame(rows)
    return {"df": df, "time_s": time.perf_counter() - t0,
            "report": hc.report()}


def show_verification():
    try:
        ver = model_a_verification()
    except Exception as exc:                                # pragma: no cover
        st.error(f"The verification could not run: {exc}")
        return
    df = ver["df"]
    errs = df["_err"].dropna()
    n_pass = int(df["_pass"].sum())
    st.markdown(
        "**Model A**: the 1V pentagonal cap of a 5 m sphere, 6 nodes, "
        "10 CHS 88.9 x 4.0 struts, 5 pinned supports, 3 free DOF. "
        "A1 = crown load (+20, 0, -60) kN. A2 = 1.2D + 1.6S from the area "
        "loads (0.30 kPa cladding, 1.20 kPa snow).")
    a, b, c = st.columns(3)
    a.metric("Quantities compared", f"{len(df)}")
    b.metric("Largest error", f"{errs.max():.1e} %" if len(errs) else "0 %")
    c.metric("Within 1e-10 %", f"{n_pass} of {len(df)}")
    st.caption(
        "Hand = validation/hand_calc.py, written the way the course solves a "
        "truss on paper (direction cosines, element matrices, the 3 x 3 "
        "free-DOF block). It does not import geodome. GeoDomeLab = the "
        "Model A preset solved by the same run_analysis call as the desktop "
        "program.")
    view = df.drop(columns=["_err", "_pass"])
    st.dataframe(view, hide_index=True, width="stretch", height=420)
    st.download_button("Download the comparison (CSV)", to_csv(view),
                       file_name="geodomelab_modelA_verification.csv",
                       mime="text/csv")
    with st.expander("Show the hand calculation printout"):
        st.code(ver["report"], language=None)


# =============================================================================
# page
# =============================================================================
st.set_page_config(page_title="GeoDomeLab", layout="wide",
                   initial_sidebar_state="collapsed")

# Smaller metric numbers, so six key results fit a phone screen.
st.markdown(
    "<style>[data-testid='stMetricValue']{font-size:1.8rem}"
    "@media (max-width:640px){[data-testid='stMetricValue']"
    "{font-size:1.45rem}}</style>", unsafe_allow_html=True)

st.title("GeoDomeLab")
st.markdown(
    f"Single-layer geodesic domes as 3D space trusses, solved by the direct "
    f"stiffness method, with AISC 360-16 member checks and linear buckling. "
    f"Same analysis core as the desktop program, version "
    f"{geodome.__version__}. [Source code on GitHub]({REPO_URL})")
if IN_BROWSER:
    st.caption("Running in your browser (Pyodide). Nothing is sent to a "
               "server.")

col1, col2 = st.columns(2)
if col1.button(EXAMPLE_6V, type="primary", width="stretch",
               help="R = 15 m hemisphere, trial section CHS 114.3 x 5.0, "
                    "dead 0.35 kPa, snow 0.60 kPa, wind 58 m/s, ASCE 7-16 "
                    "LC1 to LC4 and SLS."):
    with st.spinner("Solving the 6V dome..."):
        store(solve(bm.make_preset("case_study_6v"),
                    "6V dome, R = 15 m, 30 m span (trial design)", "example"),
              "LC2  1.2D + 1.6S")
if col2.button(EXAMPLE_A, width="stretch",
               help="1V cap of a 5 m sphere, checked against "
                    "validation/hand_calc.py."):
    with st.spinner("Solving Model A..."):
        store(solve(bm.validation_model_a(),
                    "Model A, 1V cap of a 5 m sphere", "modelA"),
              bm.COMBO_A1)

with st.expander("Custom dome: set the geometry, supports, section and loads"):
    with st.form("custom"):
        a, b, c = st.columns(3)
        freq = a.selectbox("Frequency", [1, 2, 3, 4, 5, 6], index=3,
                           format_func=lambda v: f"{v}V")
        radius = b.number_input("Sphere radius (m)", min_value=2.0,
                                max_value=60.0, value=12.0, step=0.5)
        fraction = c.slider("Dome fraction of the diameter (-)", 0.25, 0.95,
                            0.50, 0.01,
                            help="0.50 is a hemisphere. The cut snaps to the "
                                 "nearest ring of nodes.")
        a, b = st.columns(2)
        support = a.selectbox("Supports", list(SUPPORT_CHOICES),
                              format_func=SUPPORT_CHOICES.get)
        section = b.selectbox("Section, every strut", SECTION_NAMES,
                              index=SECTION_NAMES.index(DEFAULT_SECTION))
        a, b, c = st.columns(3)
        dead = a.number_input("Dead load on the surface (kPa)", 0.0, 10.0,
                              0.30, 0.05, format="%.2f",
                              help="Cladding and services. Self weight is "
                                   "added automatically.")
        snow = b.number_input("Snow load on plan (kPa)", 0.0, 10.0, 1.20,
                              0.05, format="%.2f")
        wind = c.number_input("Wind speed (m/s), 0 = no wind", 0.0, 90.0,
                              45.0, 1.0, format="%.0f")
        size_it = st.checkbox("Auto-size the struts to AISC 360-16 "
                              "(target utilisation 0.95)")
        st.caption("Material ASTM A500 Gr. C. Combinations ASCE 7-16 LC1 to "
                   "LC6 and the service case SLS 1.0D + 1.0S.")
        run_clicked = st.form_submit_button("Run", type="primary")

if run_clicked:
    try:
        with st.spinner("Solving..."):
            m = custom_model(freq, radius, fraction, support, section, dead,
                             snow, wind)
            label = (f"{m.geom.frequency}V dome, R = {m.geom.radius:g} m, "
                     f"{m.geom.base_diameter:.2f} m span")
            sizing, start = None, None
            if size_it:
                start = m.steel_mass()
                sizing = opt.auto_size(m)
                label += " (auto-sized)"
            prev = st.session_state.get("show_combo", "LC2  1.2D + 1.6S")
            store(solve(m, label, "custom", sizing, start), prev)
    except Exception as exc:
        st.error(f"The model could not be solved: {exc}")

run = st.session_state.get("run")

if run is None:
    st.info("Tap an example above to run it, or open Custom dome and press "
            "Run. Each run solves every load combination, checks every strut "
            "to AISC 360-16 and finds the critical buckling load factor. It "
            "takes about a second.")
else:
    model, res = run["model"], run["res"]
    g = model.geom
    st.subheader(run["label"])
    st.caption(
        f"{g.frequency}V Class I, sphere radius {g.radius:.2f} m, span "
        f"{g.base_diameter:.2f} m, rise {g.height:.2f} m. {g.n_nodes} nodes, "
        f"{g.n_members} struts, {g.n_strut_types()} strut families. "
        f"{SUPPORT_CHOICES.get(model.support_type, model.support_type)}. "
        f"{model.material_name}. Sections: "
        f"{', '.join(sorted(set(model.group_sections)))}. "
        f"Solved {len(res.order)} combinations with buckling in "
        f"{secs(run['time_s'])}.")
    for note in g.notes:                  # the snapping notes, not the counts
        if "snapped" in note or "stepped" in note:
            st.caption(note)

    if run.get("sizing") is not None:
        s = run["sizing"]
        st.success(
            f"Auto-sized in {s.iterations} iterations: "
            f"{s.start_mass:,.0f} kg to {s.final_mass:,.0f} kg "
            f"({abs(s.saving) * 100:.1f} % "
            f"{'lighter' if s.saving >= 0 else 'heavier'}).")
        for msg in s.messages:
            st.caption(msg)

    for msg in res.messages:                     # design failures up front
        if msg.startswith(("OVERSTRESSED", "Deflection span/")):
            st.warning(msg)

    if "pending_combo" in st.session_state:
        st.session_state["show_combo"] = st.session_state.pop("pending_combo")
    if st.session_state.get("show_combo") not in res.order:
        st.session_state["show_combo"] = res.order[0]
    combo = st.selectbox("Combination to show", res.order, key="show_combo",
                         format_func=short)
    rc = res.cases[combo]
    N = rc.N / 1e3
    c_node = int(g.crown_node)
    uz = rc.u[3 * c_node + 2] * 1e3

    a, b, c = st.columns(3)
    a.metric(f"Crown displacement uz, node {c_node}", f"{uz:+.3f} mm",
             help=f"Vertical displacement of the crown under {short(combo)}. "
                  "Negative is downward.")
    _, kt, nt = extreme(N, "max")
    if nt > 1e-9:
        b.metric(member_label("Largest tension", kt), f"{nt:+.2f} kN",
                 help=twins_help(kt, "force"))
    else:
        b.metric("Largest tension", "none")
    _, kc, nc = extreme(N, "min")
    if nc < -1e-9:
        c.metric(member_label("Largest compression", kc), f"{nc:+.2f} kN",
                 help=twins_help(kc, "force"))
    else:
        c.metric("Largest compression", "none")

    a, b, c = st.columns(3)
    if res.governing_case:
        ku, mu, uu = extreme(res.util_env, "max")
        gov = res.strength_cases()[int(res.util_case[ku])]
        a.metric(member_label(f"Largest utilisation, {combo_id(gov)}", mu),
                 f"{uu:.3f}",
                 help=f"Envelope of the strength combinations, AISC 360-16 "
                      f"D2 and E3 (E7 for slender walls), governing under "
                      f"{short(gov)}. " + (twins_help(mu, "utilisation")
                                            or ""))
    else:
        a.metric("Largest utilisation", "no strength combination")
    lam = {k: v for k, v in res.lambda_cr_by_case.items()
           if math.isfinite(v)}
    if lam:
        ref = min(lam, key=lam.get)
        b.metric(f"Critical load factor λcr, {combo_id(ref)}",
                 lam_text(lam[ref]),
                 help=f"Lowest positive linear buckling load factor over "
                      f"the strength combinations, under {short(ref)}.")
    else:
        b.metric("Critical load factor λcr", "none")
    c.metric("Steel mass", f"{model.steel_mass():,.0f} kg")

    ok = "OK" if res.deflection_ok else "exceeds the limit"
    st.caption(
        f"Service deflection under {short(res.service_case)}: largest |u| "
        f"{res.max_service_disp * 1e3:.2f} mm at node "
        f"{res.max_service_node}, span/{res.deflection_ratio:,.0f} "
        f"(limit span/{model.deflection_limit:.0f}), {ok}.")

    if run["kind"] != "modelA" and run.get("sizing") is None:
        if st.button("Auto-size the struts (AISC 360-16, target 0.95)"):
            with st.spinner("Auto-sizing..."):
                m2 = copy.deepcopy(model)
                start = m2.steel_mass()
                sizing = opt.auto_size(m2)
                store(solve(m2, run["label"].replace(" (trial design)", "")
                            + " (auto-sized)", run["kind"], sizing, start),
                      combo)
            st.rerun()

    labels = ["3D view", "Members", "Nodes", "Verification"]
    try:
        tabs = st.tabs(labels, default=("Verification"
                                        if run["kind"] == "modelA" else None))
    except TypeError:                                       # older Streamlit
        tabs = st.tabs(labels)

    with tabs[0]:
        show_figure(dome_figure(model, res, combo))
        st.caption(
            f"Axial force under {short(combo)}. Red is tension, blue is "
            f"compression, the scale is symmetric about zero. Drag to turn, "
            f"scroll or pinch to zoom, hover or tap a strut for its force. "
            f"Squares are supports.")

    with tabs[1]:
        dfm = members_frame(run, combo)
        st.dataframe(
            dfm, hide_index=True, width="stretch", height=420,
            column_config=number_config(dfm, {
                "Length (m)": 4, "KL/r (-)": 1,
                f"N, {short(combo)} (kN)": 2,
                "N max, strength envelope (kN)": 2,
                "N min, strength envelope (kN)": 2,
                "phi Pn tension (kN)": 1, "phi Pn compression (kN)": 1,
                "Utilisation (-)": 3}))
        st.caption("Tension is positive. N max and N min are the envelope "
                   "over the strength combinations, as in the desktop "
                   "Members table.")
        st.download_button("Download the member table (CSV)", to_csv(dfm),
                            file_name=f"{file_stub(run)}_members.csv",
                            mime="text/csv")

    with tabs[2]:
        dfn = nodes_frame(run, combo)
        st.dataframe(
            dfn, hide_index=True, width="stretch", height=420,
            column_config=number_config(dfn, {
                "x (m)": 3, "y (m)": 3, "z (m)": 3, "ux (mm)": 3,
                "uy (mm)": 3, "uz (mm)": 3, "|u| (mm)": 3, "Rx (kN)": 2,
                "Ry (kN)": 2, "Rz (kN)": 2}))
        st.caption(f"Displacements and support reactions under "
                   f"{short(combo)}.")
        st.download_button("Download the node table (CSV)", to_csv(dfn),
                           file_name=f"{file_stub(run)}_nodes.csv",
                           mime="text/csv")

    with tabs[3]:
        show_verification()

    if res.messages:
        with st.expander("Program messages"):
            for msg in res.messages:
                st.write(msg)

with st.expander("About this app"):
    st.markdown(
        "Every strut is a two-node pin-ended bar with 3 translations per "
        "node. The global stiffness matrix is assembled from the element "
        "matrices (EA/L) [[cc^T, -cc^T], [-cc^T, cc^T]], the supported DOF are "
        "removed and {P} = [K]{D} is solved for each ASCE 7-16 combination. "
        "Member forces give the AISC 360-16 checks, and the geometric "
        "stiffness of each strength combination gives the linear buckling "
        "load factor λcr from [K + λ Kg]{φ} = 0.\n\n"
        "The desktop program (PyQt6) adds temperature and settlement loads, "
        "corotational arc-length nonlinear analysis with imperfections, "
        "modal analysis, fabrication data and CSV import. "
        f"Install it from the [repository]({REPO_URL}).")
    st.caption(f"GeoDomeLab {geodome.__version__}. 3D view: "
               f"{'plotly' if HAVE_PLOTLY else 'plotly.js (pinned CDN)'}.")
