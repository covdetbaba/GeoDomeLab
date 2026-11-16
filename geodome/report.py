"""
Export: CSV tables, a self-contained HTML calculation report, and PNG figures.

The HTML report embeds its figures as base64 PNGs so that a single file can be
submitted or e-mailed with nothing else attached.

Report sections
---------------
1 Structure, 2 Loading, 3 Analysis results, 4 Auto-sizing history (when a
sizing result is given), 5 Figures (``include_figures``), 6 Solution
verification and 7 Benchmark verification suite (when ``verification`` rows
are given).  These are always written the same way.

``export_html(..., extras=dict)`` adds one section per result of the optional
analyses, after section 7 and in this fixed order:

    key                    section                          built from
    nonlinear_result       8  Geometrically nonlinear       nonlinear.run_nonlinear
    imperfection_result    9  Imperfection sensitivity      nonlinear.imperfection_study
    modal_result           10 Modal analysis                modal.modal_analysis
    fabrication_tables     11 Fabrication schedule          fabrication.fabrication_tables
    study_rows             12 Parameter study               study.parameter_study

A key that is missing, ``None`` or empty adds nothing, so the report without
extras is exactly the 1-7 report.  Unknown keys are ignored.  Each extra
section is built on its own: if one fails, it shows a short note and the rest
of the report is still written.  Every new heading carries an ``id`` (the
second column of ``EXTRA_SECTIONS``) so that it can be linked to.
"""

from __future__ import annotations

import base64
import csv
import datetime as _dt
import html
import io
import math
import os

import numpy as np
from matplotlib.figure import Figure

from . import plotting as plot
from . import postprocess as pp
from . import loads as ld


# ------------------------------------------------------------------- CSV ----
def _write_rows(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def export_csv(res, folder, prefix="dome"):
    """Write the member / node / reaction / schedule tables as CSV files."""
    os.makedirs(folder, exist_ok=True)
    written = []

    p = os.path.join(folder, f"{prefix}_members.csv")
    _write_rows(p, pp.member_table(res))
    written.append(p)

    p = os.path.join(folder, f"{prefix}_nodes.csv")
    _write_rows(p, pp.node_table(res))
    written.append(p)

    p = os.path.join(folder, f"{prefix}_reactions.csv")
    _write_rows(p, pp.reaction_table(res))
    written.append(p)

    p = os.path.join(folder, f"{prefix}_strut_schedule.csv")
    sched = res.model.geom.strut_schedule()
    _write_rows(p, [{"family": s["label"], "length_m": round(s["length"], 6),
                     # empty for a custom truss with no defining sphere
                     "chord_factor": (round(s["chord_factor"], 6)
                                      if math.isfinite(s["chord_factor"])
                                      else ""),
                     "count": s["count"],
                     "section": res.model.group_sections[s["id"]]}
                    for s in sched])
    written.append(p)

    p = os.path.join(folder, f"{prefix}_case_forces.csv")
    rows = []
    for k in range(res.model.geom.n_members):
        row = {"member": k}
        for name in res.order:
            row[f"N_{name}_kN"] = round(float(res.cases[name].N[k] / 1e3), 3)
        rows.append(row)
    _write_rows(p, rows)
    written.append(p)

    p = os.path.join(folder, f"{prefix}_summary.txt")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(pp.summary_text(res))
    written.append(p)
    return written


# ----------------------------------------------------------------- figures --
def _fig3d(draw, *args, size=(7.4, 5.8), dpi=130, colorbar_label=None, **kw):
    fig = Figure(figsize=size, dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111, projection="3d")
    art = draw(ax, *args, **kw)
    if colorbar_label and hasattr(art, "set_array"):
        cb = fig.colorbar(art, ax=ax, shrink=0.60, pad=0.10)
        cb.set_label(colorbar_label, fontsize=9)
        cb.ax.tick_params(labelsize=8)
    fig.tight_layout()
    return fig


def standard_figures(res, model=None):
    """Build the figure set used by the HTML report and by ``export_figures``."""
    model = model or res.model
    gov = res.governing_case or res.order[0]
    figs = {}

    figs["01_model"] = _fig3d(plot.draw_model, model)
    basic = ld.basic_load_vectors(model)
    F_gov = ld.combine(basic, next(c for c in model.combinations if c.name == gov))
    figs["02_loads"] = _fig3d(plot.draw_loads, model, F_gov,
                              title=f"Applied nodal loads -- {gov}")
    figs["03_deflected"] = _fig3d(
        plot.draw_deflected, model, res.cases[res.service_case],
        colorbar_label="nodal displacement [mm]")
    figs["04_axial"] = _fig3d(
        plot.draw_axial, model, res.cases[gov].N,
        title=f"Member axial force -- {gov}",
        colorbar_label="axial force [kN]  (+ tension)")
    figs["05_utilisation"] = _fig3d(plot.draw_utilisation, model, res.util_env)
    if res.buckling is not None and len(res.buckling.load_factors):
        figs["06_buckling"] = _fig3d(plot.draw_buckling_mode, model,
                                     res.buckling, 0,
                                     colorbar_label="normalised mode amplitude")

    fig = Figure(figsize=(7.4, 4.2), dpi=130)
    plot.draw_strut_histogram(fig.add_subplot(111), model)
    fig.tight_layout()
    figs["07_struts"] = fig

    fig = Figure(figsize=(7.4, 4.2), dpi=130)
    plot.draw_utilisation_histogram(fig.add_subplot(111), res)
    fig.tight_layout()
    figs["08_utilisation_hist"] = fig
    return figs


def export_figures(res, folder, prefix="dome", dpi=160):
    os.makedirs(folder, exist_ok=True)
    written = []
    for name, fig in standard_figures(res).items():
        p = os.path.join(folder, f"{prefix}_{name}.png")
        fig.savefig(p, dpi=dpi, facecolor="white")
        written.append(p)
    return written


def _fig_to_b64(fig, dpi=120):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
    return base64.b64encode(buf.getvalue()).decode("ascii")


# ------------------------------------------------------------------ HTML ----
_CSS = """
body{font-family:'Segoe UI',Calibri,Arial,sans-serif;margin:0;background:#eef1f4;color:#1c2733}
.wrap{max-width:1040px;margin:0 auto;background:#fff;padding:38px 48px 60px;
      box-shadow:0 1px 14px rgba(0,0,0,.12)}
h1{font-size:25px;margin:0 0 4px;color:#123a5c;letter-spacing:.2px}
.sub{color:#5b6b7a;font-size:13px;margin-bottom:26px}
h2{font-size:17px;margin:34px 0 10px;color:#123a5c;border-bottom:2px solid #d6dee6;
   padding-bottom:5px}
h3{font-size:14px;margin:20px 0 6px;color:#2b4c6f}
table{border-collapse:collapse;width:100%;font-size:12px;margin:8px 0 16px}
th{background:#2b4c6f;color:#fff;text-align:left;padding:6px 8px;font-weight:600}
td{padding:5px 8px;border-bottom:1px solid #e3e9ee}
tr:nth-child(even) td{background:#f6f9fb}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
pre{background:#f5f7f9;border:1px solid #dde4ea;border-left:4px solid #2b4c6f;
    padding:12px 14px;font-size:11.5px;overflow-x:auto;line-height:1.45}
img{width:100%;border:1px solid #dde4ea;border-radius:3px;margin:6px 0 4px}
.cap{font-size:11.5px;color:#5b6b7a;margin-bottom:20px}
.ok{color:#15803d;font-weight:700}.bad{color:#b4161b;font-weight:700}
.kv{display:grid;grid-template-columns:260px 1fr;gap:3px 14px;font-size:13px;margin:10px 0}
.kv div:nth-child(odd){color:#5b6b7a}
.note{background:#fff7e6;border-left:4px solid #e9a000;padding:9px 13px;
      font-size:12.5px;margin:10px 0}
footer{margin-top:40px;padding-top:14px;border-top:1px solid #d6dee6;
       font-size:11px;color:#7b8a99}
"""


def _table(rows, numeric_cols=(), limit=None, highlight=None):
    if not rows:
        return "<p><i>no data</i></p>"
    shown = rows if limit is None else rows[:limit]
    cols = list(shown[0].keys())
    out = ["<table><thead><tr>"]
    for c in cols:
        cls = ' class="num"' if c in numeric_cols else ""
        out.append(f"<th{cls}>{html.escape(str(c))}</th>")
    out.append("</tr></thead><tbody>")
    for r in shown:
        flag = highlight(r) if highlight else False
        out.append('<tr style="background:#ffecec">' if flag else "<tr>")
        for c in cols:
            v = r.get(c, "")
            cls = ' class="num"' if c in numeric_cols else ""
            if isinstance(v, float):
                v = f"{v:,.3f}"
            out.append(f"<td{cls}>{html.escape(str(v))}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    if limit is not None and len(rows) > limit:
        out.append(f'<p class="cap">showing the first {limit} of '
                   f'{len(rows)} rows -- the full table is in the CSV export</p>')
    return "".join(out)


def export_html(res, path, sizing=None, project_title="Geodesic Dome Analysis",
                author="", include_figures=True, verification=None,
                extras=None):
    """Write a single self-contained HTML calculation report.

    ``extras`` is an optional dict with any of the keys in
    :data:`EXTRA_KEYS` (``nonlinear_result``, ``imperfection_result``,
    ``modal_result``, ``fabrication_tables``, ``study_rows``).  Each present
    value adds one section after section 7; see the module docstring.
    """
    m = res.model
    g = m.geom
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    members = pp.member_table(res)
    over = [r for r in members if r["utilisation"] > 1.0]
    worst = sorted(members, key=lambda r: -r["utilisation"])[:20]

    parts = [f"<!doctype html><html><head><meta charset='utf-8'>"
             f"<title>{html.escape(project_title)}</title>"
             f"<style>{_CSS}</style></head><body><div class='wrap'>"]
    parts.append(f"<h1>{html.escape(project_title)}</h1>")
    parts.append(f"<div class='sub'>GeoDomeLab &ndash; geodesic dome "
                 f"space-truss analyser &middot; generated {now}"
                 + (f" &middot; {html.escape(author)}" if author else "")
                 + "</div>")

    # ---- 1 model ----------------------------------------------------------
    parts.append("<h2>1 &nbsp; Structure</h2><div class='kv'>")
    for line in m.summary_lines():
        k, _, v = line.partition("  ")
        parts.append(f"<div>{html.escape(k.strip())}</div>"
                     f"<div>{html.escape(v.strip())}</div>")
    parts.append("</div>")

    parts.append("<h3>Strut schedule (chord factors)</h3>")
    parts.append(_table(
        [{"family": s["label"],
          "chord factor": (round(s["chord_factor"], 6)
                           if math.isfinite(s["chord_factor"]) else "-"),
          "length [m]": round(s["length"], 4), "count": s["count"],
          "section": m.group_sections[s["id"]]}
         for s in g.strut_schedule()],
        numeric_cols=("chord factor", "length [m]", "count")))

    # ---- 2 loading --------------------------------------------------------
    parts.append("<h2>2 &nbsp; Loading</h2><pre>"
                 + html.escape("\n".join(ld.describe_loads(m))) + "</pre>")
    parts.append(_table([{"combination": c.name, "expression": c.label(),
                          "type": "service" if c.serviceability else "strength"}
                         for c in m.combinations]))

    # ---- 3 results --------------------------------------------------------
    parts.append("<h2>3 &nbsp; Analysis results</h2>")
    rows = []
    for name in res.order:
        r = res.cases[name]
        _, dmax = r.max_displacement()
        rows.append({
            "combination": name,
            "max tension [kN]": round(float(r.N.max() / 1e3), 2),
            "max compression [kN]": round(float(r.N.min() / 1e3), 2),
            "max |u| [mm]": round(dmax * 1e3, 3),
            "max utilisation": round(float(np.max(res.util[name])), 3),
            "vert. reaction [kN]": round(float(r.total_reaction()[2] / 1e3), 2),
        })
    parts.append(_table(rows, numeric_cols=tuple(rows[0].keys())[1:]))

    status = ("<span class='ok'>SATISFIED</span>" if res.governing_util <= 1.0
              else "<span class='bad'>EXCEEDED</span>")
    dstat = ("<span class='ok'>SATISFIED</span>" if res.deflection_ok
             else "<span class='bad'>EXCEEDED</span>")
    parts.append("<div class='kv'>"
                 f"<div>Governing combination</div><div>{html.escape(res.governing_case)}</div>"
                 f"<div>Governing member</div><div>#{res.governing_member}</div>"
                 f"<div>Maximum utilisation</div><div>{res.governing_util:.3f} &nbsp; {status}</div>"
                 f"<div>Members over capacity</div><div>{len(over)} of {g.n_members}</div>"
                 f"<div>Max displacement</div><div>{res.max_service_disp * 1e3:.2f} mm "
                 f"(node {res.max_service_node}, {html.escape(res.service_case)})</div>"
                 f"<div>Deflection ratio</div><div>span/{res.deflection_ratio:.0f} "
                 f"vs limit span/{m.deflection_limit:.0f} &nbsp; {dstat}</div>"
                 "</div>")

    if res.buckling is not None and len(res.buckling.load_factors):
        parts.append("<h3>Linear buckling</h3>")
        parts.append("<p>Critical combination (lowest load factor): "
                     f"{html.escape(res.buckling.reference_case)}</p>")
        parts.append(_table(
            [{"mode": i + 1, "critical load factor": round(float(l), 3)}
             for i, l in enumerate(res.buckling.load_factors)],
            numeric_cols=("critical load factor",)))
        by_case = getattr(res, "lambda_cr_by_case", {}) or {}
        if len(by_case) > 1:
            parts.append(_table(
                [{"combination": n,
                  "first load factor": (round(float(v), 3)
                                        if np.isfinite(v) else "-")}
                 for n, v in by_case.items()],
                numeric_cols=("first load factor",)))
        for w in res.buckling.warnings:
            parts.append(f"<div class='note'>{html.escape(w)}</div>")

    parts.append("<h3>Twenty most highly utilised struts</h3>")
    parts.append(_table(worst, numeric_cols=(
        "length_m", "KL_r", "N_tension_kN", "N_compr_kN", "phiPn_t_kN",
        "phiPn_c_kN", "utilisation"),
        highlight=lambda r: r["utilisation"] > 1.0))

    parts.append("<h3>Support reactions &ndash; "
                 f"{html.escape(res.governing_case or res.order[0])}</h3>")
    parts.append(_table(pp.reaction_table(res), numeric_cols=(
        "x_m", "y_m", "Rx_kN", "Ry_kN", "Rz_kN", "thrust_kN"), limit=40))

    # ---- 4 sizing ---------------------------------------------------------
    if sizing is not None:
        parts.append("<h2>4 &nbsp; Auto-sizing history</h2>")
        parts.append(_table(sizing.table()))
        parts.append(f"<div class='kv'><div>Start mass</div>"
                     f"<div>{sizing.start_mass:,.0f} kg</div>"
                     f"<div>Final mass</div><div>{sizing.final_mass:,.0f} kg</div>"
                     f"<div>Saving</div><div>{sizing.saving * 100:.1f} %</div></div>")
        for msg in sizing.messages:
            parts.append(f"<div class='note'>{html.escape(msg)}</div>")

    # ---- 5 figures --------------------------------------------------------
    if include_figures:
        parts.append("<h2>5 &nbsp; Figures</h2>")
        captions = {
            "01_model": "Generated geodesic geometry, struts coloured by length family.",
            "02_loads": "Nodal load vectors for the governing combination.",
            "03_deflected": "Deflected shape under the serviceability combination.",
            "04_axial": "Member axial forces, governing combination.",
            "05_utilisation": "Demand/capacity envelope over all strength combinations.",
            "06_buckling": "First linear buckling mode.",
            "07_struts": "Fabrication strut schedule.",
            "08_utilisation_hist": "Distribution of the utilisation envelope.",
        }
        for key, fig in sorted(standard_figures(res).items()):
            parts.append(f"<img src='data:image/png;base64,{_fig_to_b64(fig)}'>")
            parts.append(f"<div class='cap'>Figure {key[:2]} &ndash; "
                         f"{html.escape(captions.get(key, ''))}</div>")

    # ---- 6 verification ---------------------------------------------------
    parts.append("<h2>6 &nbsp; Solution verification</h2><pre>")
    vl = []
    for name in res.order:
        r = res.cases[name]
        vl.append(f"{name:<24s} equilibrium residual {r.equilibrium_residual:.3e}   "
                  f"strain-energy mismatch {r.energy_error:.3e}   "
                  f"load/reaction balance {np.max(np.abs(r.reaction_balance)):.3e} N")
    parts.append(html.escape("\n".join(vl)) + "</pre>")
    for msg in res.messages:
        parts.append(f"<div class='note'>{html.escape(msg)}</div>")

    # ---- 7 benchmark suite -------------------------------------------------
    if verification:
        n_pass = sum(1 for r in verification if r["pass"])
        parts.append("<h2>7 &nbsp; Benchmark verification suite</h2>")
        badge = ("<span class='ok'>ALL PASS</span>" if n_pass == len(verification)
                 else f"<span class='bad'>{len(verification) - n_pass} FAILED</span>")
        parts.append(f"<p>{n_pass} of {len(verification)} independent checks "
                     f"passed &nbsp; {badge}. Every reference value below comes "
                     f"from a closed-form solution, a published table or a "
                     f"theoretical identity &ndash; not from this program.</p>")
        parts.append(_table(
            [{"group": r["group"], "check": r["check"],
              "computed": float(r["computed"]), "reference": float(r["reference"]),
              "unit": r["unit"], "rel. error": float(r["rel_error"]),
              "tolerance": float(r["tolerance"]),
              "result": "PASS" if r["pass"] else "FAIL",
              "basis": r["note"]}
             for r in verification],
            numeric_cols=("computed", "reference", "rel. error", "tolerance"),
            highlight=lambda r: r["result"] == "FAIL"))

    # ---- 8 to 12 optional analyses (only when given) ---------------------
    present = _present_extras(extras)
    for key in EXTRA_KEYS:
        if key in present:
            parts.append(_extra_section(key, present[key], res,
                                        include_figures))

    parts.append("<footer>Produced by GeoDomeLab. Linear elastic first-order "
                 "space-truss analysis; member checks to AISC 360 LRFD; "
                 "wind pressures from a simplified domed-roof coefficient "
                 "distribution. Results are for academic use."
                 + _extras_footer(present) + "</footer>")
    parts.append("</div></body></html>")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(parts))
    return path


# ======================================================= optional sections ==
# key -> (heading id, section number, heading, footer phrase)
EXTRA_SECTIONS = {
    "nonlinear_result": ("nonlinear", 8, "Geometrically nonlinear analysis",
                         "geometrically nonlinear path (corotational bars, "
                         "elastic material)"),
    "imperfection_result": ("imperfection", 9, "Imperfection sensitivity",
                            "limit loads of mode-affine imperfect domes"),
    "modal_result": ("modal", 10, "Modal analysis",
                     "undamped natural frequencies and modes"),
    "fabrication_tables": ("fabrication", 11, "Fabrication schedule",
                           "cut list, hub and panel types"),
    "study_rows": ("study", 12, "Parameter study",
                   "parameter study of regenerated domes"),
}
EXTRA_KEYS = tuple(EXTRA_SECTIONS)

_METHOD_NAMES = {
    "arc-length": "Arc-length (Crisfield, cylindrical)",
    "displacement": "Displacement control",
    "load": "Load control",
}


def _is_empty(value):
    if value is None:
        return True
    if isinstance(value, (list, tuple, dict)) and len(value) == 0:
        return True
    return False


def _present_extras(extras):
    """{key: value} for the known keys whose value is given and not empty."""
    if not extras:
        return {}
    try:
        items = dict(extras)
    except (TypeError, ValueError):
        return {}
    return {k: items[k] for k in EXTRA_KEYS
            if k in items and not _is_empty(items[k])}


def _extras_footer(present):
    if not present:
        return ""
    bits = [f"{EXTRA_SECTIONS[k][1]} {EXTRA_SECTIONS[k][3]}"
            for k in EXTRA_KEYS if k in present]
    return (" Extra sections: " + "; ".join(html.escape(b) for b in bits)
            + ".")


def _extra_section(key, value, res, include_figures):
    """Heading plus body of one optional section (never raises)."""
    sid, num, title, _ = EXTRA_SECTIONS[key]
    head = f"<h2 id='{sid}'>{num} &nbsp; {html.escape(title)}</h2>"
    builder = {
        "nonlinear_result": _section_nonlinear,
        "imperfection_result": _section_imperfection,
        "modal_result": _section_modal,
        "fabrication_tables": _section_fabrication,
        "study_rows": _section_study,
    }[key]
    try:
        body = builder(value, res, num, include_figures)
    except Exception as exc:                      # keep the rest of the report
        body = [f"<div class='note'>This section could not be built: "
                f"{html.escape(type(exc).__name__)}: "
                f"{html.escape(str(exc))}</div>"]
    return head + "".join(body)


# ----------------------------------------------------------------- helpers --
def _e(x):
    return html.escape(str(x))


def _fnum(v, fmt="{:.4g}", nan="n/a"):
    """Format a number for a table cell (NaN -> ``nan``, inf -> 'inf')."""
    if v is None:
        return nan
    if isinstance(v, str):
        return v
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(x):
        return nan
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return fmt.format(x)


def _kv(pairs):
    """Key/value grid; keys are escaped, values must already be safe HTML."""
    out = ["<div class='kv'>"]
    for k, v in pairs:
        out.append(f"<div>{_e(k)}</div><div>{v}</div>")
    out.append("</div>")
    return "".join(out)


def _spec_rows(rows, spec):
    """Rows with display headers and formatted cells.

    ``spec`` = [(key, header, fmt)], fmt None (text), "int" or a format
    string.  Returns (rows, numeric headers) for :func:`_table`.
    """
    out = []
    for r in rows:
        d = {}
        for key, head, fmt in spec:
            v = r.get(key, "")
            if fmt is None:
                d[head] = "" if v is None else str(v)
            elif fmt == "int":
                d[head] = _fnum(v, "{:.0f}", nan="-")
            else:
                d[head] = _fnum(v, fmt, nan="-")
        out.append(d)
    numeric = tuple(head for _, head, fmt in spec if fmt is not None)
    return out, numeric


def _img(fig, number, caption):
    return (f"<img src='data:image/png;base64,{_fig_to_b64(fig)}'>"
            f"<div class='cap'>Figure {number} &ndash; {_e(caption)}</div>")


def _notes(lines):
    return [f"<div class='note'>{_e(m)}</div>" for m in (lines or []) if m]


# ------------------------------------------------------- 8 nonlinear path --
def _control_sign(result):
    """-1 when the path moves mainly in the negative control direction."""
    c = np.asarray(result.control, dtype=float)
    c = c[np.isfinite(c)]
    if not len(c):
        return 1.0
    return -1.0 if abs(float(c.min())) > abs(float(c.max())) else 1.0


def _control_axis_label(result, sign):
    label = str(result.control_label or "control DOF")
    if sign > 0:
        return f"{label} [mm]"
    if int(result.control_dof) % 3 == 2:
        return f"Downward displacement, {label} [mm]"
    return f"Displacement, minus {label} [mm]"


def nonlinear_path_figure(result, size=(7.4, 4.6), dpi=130):
    """Load factor against control displacement: corotational, linear and
    P-Delta paths, lambda_cr as a horizontal line, the limit point marked."""
    fig = Figure(figsize=size, dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    sgn = _control_sign(result)
    s = sgn * 1e3
    curves = result.curves()
    cx, clam = (np.asarray(a, dtype=float) for a in curves["corotational"])
    xs = cx * s
    ax.plot(xs, clam, "-o", color="#1f6fb4", lw=1.8, ms=2.6,
            label="Corotational (geometrically nonlinear)")
    if "linear" in curves:
        lx, llam = (np.asarray(a, dtype=float) for a in curves["linear"])
        ax.plot(lx * s, llam, "--", color="0.45", lw=1.2,
                label="Linear (first order)")
    if "p_delta" in curves:
        px, plam = (np.asarray(a, dtype=float) for a in curves["p_delta"])
        ax.plot(px * s, plam, "-.", color="#e9a000", lw=1.3,
                label="P-Delta (linearised second order)")

    fin = np.isfinite(xs) & np.isfinite(clam)
    x_lo = min(0.0, float(xs[fin].min())) if fin.any() else 0.0
    x_hi = max(0.0, float(xs[fin].max())) if fin.any() else 1.0
    y_lo = min(0.0, float(clam[fin].min())) if fin.any() else 0.0
    y_hi = float(clam[fin].max()) if fin.any() else 1.0
    lam_cr = float(result.lambda_cr)
    cr_shown = math.isfinite(lam_cr) and 0.0 < lam_cr <= 8.0 * max(y_hi, 1e-30)
    if cr_shown:
        y_hi = max(y_hi, lam_cr)
        ax.axhline(lam_cr, ls=":", lw=1.4, color="#b4161b",
                   label=f"λ_cr = {lam_cr:.4g} (linear buckling)")
    elif math.isfinite(lam_cr):
        ax.text(0.98, 0.04, f"λ_cr = {lam_cr:.4g} is above the plotted range",
                transform=ax.transAxes, ha="right", fontsize=8, color="#b4161b")

    if result.limit_found and math.isfinite(float(result.lambda_limit)):
        kind = ("bifurcation" if result.limit_type == "bifurcation"
                else "limit point")
        xl = float(result.control_limit) * s
        yl = float(result.lambda_limit)
        ax.plot([xl], [yl], ls="none", marker="*", ms=14, color="#b4161b",
                mec="#1c2733", zorder=6,
                label=f"First load maximum ({kind}), λ = {yl:.4g}")
        ax.annotate(f"{kind}\nλ_limit = {yl:.4g}", (xl, yl),
                    textcoords="offset points", xytext=(10, -26), fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                              alpha=0.9))
    bif = [cp for cp in (result.critical_points or [])
           if cp.get("type") == "bifurcation" and not cp.get("first_peak")]
    if bif:
        ax.plot([float(cp["control"]) * s for cp in bif],
                [float(cp["lambda"]) for cp in bif], ls="none", marker="D",
                ms=6, color="#7b53a6", zorder=6,
                label="Bifurcation on the rising path")

    xpad = 0.05 * max(x_hi - x_lo, 1e-9)
    ypad = 0.08 * max(y_hi - y_lo, 1e-9)
    ax.set_xlim(x_lo - xpad, x_hi + xpad)
    ax.set_ylim(y_lo - (ypad if y_lo < 0 else 0.0), y_hi + ypad)
    ax.axhline(0.0, color="0.75", lw=0.8)
    ax.set_xlabel(_control_axis_label(result, sgn))
    ax.set_ylabel("Load factor λ")
    title = "Load-deflection path"
    if result.combination:
        title += f"  --  {result.combination}"
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best", framealpha=0.9)
    fig.tight_layout()
    return fig


def _section_nonlinear(r, res, num, include_figures):
    out = []
    method = _METHOD_NAMES.get(r.method, str(r.method))
    out.append("<p>Corotational bar elements with engineering strain, full "
               "Newton-Raphson iterations and the path-following method "
               "below. The load factor &lambda; multiplies the whole "
               "combination (the self-straining load T is not part of the "
               "path). The struts stay elastic: member buckling and "
               "yielding are reported as load factors, not modelled.</p>")
    rows = [("Method", _e(method)),
            ("Combination", _e(r.combination or "-")),
            ("Control DOF", _e(r.control_label)),
            ("Path", _e(f"{max(r.n_points - 1, 0)} steps, "
                        f"{int(np.sum(r.iterations))} iterations, "
                        f"{r.elapsed:.2f} s"))]
    if r.limit_found:
        kind = ("bifurcation" if r.limit_type == "bifurcation"
                else "limit point")
        rows.append(("First load maximum λ_limit",
                     _e(f"{r.lambda_limit:.6g} ({kind}) at "
                        f"{r.control_label} = "
                        f"{r.control_limit * 1e3:.3f} mm"
                        + ("" if r.limit_refined else ", not refined"))))
    else:
        rows.append(("First load maximum λ_limit",
                     _e(f"not reached ({r.stop_reason or 'path ended'})")))
    lfc = float(r.lambda_first_critical)
    if math.isfinite(lfc) and r.limit_found and lfc < r.lambda_limit:
        rows.append(("First critical point",
                     _e(f"λ = {lfc:.6g} (bifurcation before the load "
                        f"maximum)")))
    if math.isfinite(float(r.lambda_cr)):
        rows.append(("Linear buckling factor λ_cr", _e(f"{r.lambda_cr:.6g}")))
    if math.isfinite(float(r.limit_ratio)):
        rows.append(("λ_limit / λ_cr", _e(f"{r.limit_ratio:.4f}")))
    if math.isfinite(float(r.lambda_member_euler)):
        rows.append(("First strut at its Euler load",
                     _e(f"λ = {r.lambda_member_euler:.4g} "
                        f"(member {r.member_euler})")))
    if math.isfinite(float(r.lambda_member_design)):
        rows.append(("First strut at its design strength",
                     _e(f"λ = {r.lambda_member_design:.4g} "
                        f"(member {r.member_design})")))
    rows.append(("Largest residual", _e(f"{r.max_residual:.1e}")))
    rows.append(("Path energy mismatch",
                 _e(f"{r.energy_error:.1e} (trapezoidal work vs strain "
                    f"energy)")))
    out.append(_kv(rows))

    if r.critical_points:
        crit = [{"type": cp.get("type", ""),
                 "step": cp.get("step", ""),
                 "load factor": _fnum(cp.get("lambda"), "{:.6g}"),
                 "control displacement [mm]": _fnum(
                     1e3 * float(cp.get("control", float("nan"))), "{:.3f}"),
                 "refined": "yes" if cp.get("refined") else "no"}
                for cp in r.critical_points]
        out.append("<h3>Critical points on the path</h3>")
        out.append(_table(crit, numeric_cols=(
            "step", "load factor", "control displacement [mm]")))
    if include_figures:
        out.append(_img(nonlinear_path_figure(r), f"{num}.1",
                        "Load factor against control displacement. The star "
                        "marks the first load maximum; the dotted line is "
                        "the linear buckling factor."))
    out += _notes(r.messages)
    return out


# ------------------------------------------------------- 9 imperfections --
def _imperfection_rows(value):
    """Rows of ``nonlinear.imperfection_study`` from its list or
    (rows, results) form."""
    if isinstance(value, tuple) and len(value) == 2 and \
            isinstance(value[0], (list, tuple)):
        value = value[0]
    rows = getattr(value, "rows", value)
    rows = list(rows)
    if not rows or not all(isinstance(r, dict) for r in rows):
        raise TypeError("expected the rows of nonlinear.imperfection_study")
    return rows


def imperfection_figure(rows, size=(7.4, 4.2), dpi=130):
    """lambda_limit / lambda_cr against imperfection amplitude / span."""
    rows = _imperfection_rows(rows)
    fig = Figure(figsize=size, dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    a_span = np.array([float(r.get("amplitude_over_span", np.nan))
                       for r in rows])
    if np.all(np.isfinite(a_span)):
        x = a_span
        ax.set_xlabel("Imperfection amplitude / span")
    else:
        x = np.array([1e3 * float(r.get("amplitude_m", np.nan))
                      for r in rows])
        ax.set_xlabel("Imperfection amplitude [mm]")
    ratio = np.array([float(r.get("ratio", np.nan)) for r in rows])
    crit = np.array([float(r.get("ratio_critical", np.nan)) for r in rows])
    if np.any(np.isfinite(ratio)):
        y, yc = ratio, crit
        ax.set_ylabel("λ_limit / λ_cr")
        ax.axhline(1.0, ls=":", color="0.5", lw=1.0)
    else:
        y = np.array([float(r.get("lambda_limit", np.nan)) for r in rows])
        yc = np.array([float(r.get("lambda_critical", np.nan))
                       for r in rows])
        ax.set_ylabel("λ_limit")
    order = np.argsort(x)
    ax.plot(x[order], y[order], "-o", color="#1f6fb4", lw=1.8, ms=5,
            label="First load maximum (governing sign)")
    differ = np.isfinite(yc) & np.isfinite(y) & \
        (np.abs(yc - y) > 1e-9 * np.maximum(np.abs(y), 1e-30))
    if differ.any():
        ax.plot(x[order], yc[order], "--s", color="#d1495b", lw=1.2, ms=4,
                label="First critical point (incl. bifurcation)")
    for xi, yi, r in zip(x, y, rows):
        if not np.isfinite(yi):
            continue
        soa = float(r.get("span_over_amplitude", np.inf))
        tag = "perfect" if not np.isfinite(soa) else f"span/{soa:.0f}"
        ax.annotate(tag, (xi, yi), textcoords="offset points",
                    xytext=(4, 6), fontsize=7.5, color="#2b4c6f")
    top = np.nanmax(np.concatenate([y[np.isfinite(y)], [1.0]]))
    ax.set_ylim(0.0, 1.12 * top)
    ax.set_title("Imperfection sensitivity", fontsize=11, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def _section_imperfection(value, res, num, include_figures):
    rows = _imperfection_rows(value)
    out = ["<p>Each node is moved by a &times; &phi;, where &phi; is a "
           "linear buckling mode scaled so that its largest component is 1, "
           "so a is the largest node offset. Both signs are run and the "
           "lower limit load is kept. The first row is the perfect "
           "structure.</p>"]
    table = []
    for r in rows:
        a = float(r.get("amplitude_m", 0.0))
        soa = float(r.get("span_over_amplitude", np.inf))
        sign = float(r.get("sign", 1.0))
        table.append({
            "amplitude [mm]": _fnum(1e3 * a, "{:.2f}"),
            "span / amplitude": ("perfect" if a == 0.0 or not
                                 math.isfinite(soa) else f"{soa:.0f}"),
            "amplitude / strut": _fnum(r.get("amplitude_over_strut"),
                                       "{:.4f}"),
            "λ_limit": _fnum(r.get("lambda_limit"), "{:.4g}"),
            "λ_cr": _fnum(r.get("lambda_cr"), "{:.4g}"),
            "λ_limit / λ_cr": _fnum(r.get("ratio"), "{:.3f}"),
            "first critical λ": _fnum(r.get("lambda_critical"), "{:.4g}"),
            "sign": "" if a == 0.0 else ("+" if sign > 0 else "−"),
            "λ (+)": _fnum(r.get("lambda_plus"), "{:.4g}"),
            "λ (−)": _fnum(r.get("lambda_minus"), "{:.4g}"),
        })
    out.append(_table(table, numeric_cols=(
        "amplitude [mm]", "span / amplitude", "amplitude / strut", "λ_limit",
        "λ_cr", "λ_limit / λ_cr", "first critical λ", "λ (+)", "λ (−)"),
        highlight=lambda d: d["λ_limit"] in ("n/a", "-")))
    if any(not r.get("found", True) for r in rows):
        out.append("<div class='note'>Rows shown in red reached no load "
                   "maximum within the path limits.</div>")
    if include_figures:
        out.append(_img(imperfection_figure(rows), f"{num}.1",
                        "Limit load of the imperfect dome as a fraction of "
                        "the linear buckling factor."))
    return out


# -------------------------------------------------------------- 10 modal --
def mode_shapes_figure(model, modal, n=3, dpi=130):
    """The first ``n`` mode shapes side by side (3D)."""
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    geom = model.geom
    k_max = int(min(n, modal.n_modes))
    fig = Figure(figsize=(4.0 * max(k_max, 1) + 0.4, 4.4), dpi=dpi)
    fig.patch.set_facecolor("white")
    segs0 = np.stack([geom.nodes[geom.members[:, 0]],
                      geom.nodes[geom.members[:, 1]]], axis=1)
    scale = 0.10 * float(geom.length_scale)
    for k in range(k_max):
        ax = fig.add_subplot(1, k_max, k + 1, projection="3d")
        shape = modal.mode_shape(k)
        defo = geom.nodes + scale * shape
        mag = np.linalg.norm(shape, axis=1)
        f = float(modal.frequencies[k])
        head = (f"Mode {k + 1}   f = {f:.3f} Hz" if math.isfinite(f)
                else f"Mode {k + 1}   unstable")
        plot.prepare_axes(ax, head)
        ax.title.set_fontsize(10)
        ax.set_zlabel("")                  # small panels: keep x, y only
        ax.set_zticklabels([])
        ax.tick_params(labelsize=7)
        ax.add_collection3d(Line3DCollection(segs0, colors="#ccd3da",
                                             linewidths=0.6))
        segs = np.stack([defo[geom.members[:, 0]],
                         defo[geom.members[:, 1]]], axis=1)
        vals = 0.5 * (mag[geom.members[:, 0]] + mag[geom.members[:, 1]])
        lc = Line3DCollection(segs, cmap="viridis", linewidths=1.6,
                              norm=Normalize(0.0, max(vals.max(), 1e-12)))
        lc.set_array(vals)
        ax.add_collection3d(lc)
        b = geom.nodes[geom.base_nodes]
        if len(b):
            ax.scatter(b[:, 0], b[:, 1], b[:, 2], marker="^", s=18,
                       c="#111111", depthshade=False)
        plot.set_equal_aspect(ax, np.vstack([geom.nodes, defo]))
        T = float(modal.periods[k])
        dom = modal.dominant_direction(k)
        info = (f"T = {T:.4f} s" if math.isfinite(T) else "T = inf")
        info += f", mainly {dom}" if dom != "-" else ""
        ax.text2D(0.02, 0.95, info, transform=ax.transAxes, fontsize=7.5,
                  va="top", bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                      ec="0.7", alpha=0.85))
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.04, top=0.92,
                        wspace=0.04)
    return fig


def _section_modal(m, res, num, include_figures):
    out = ["<p>Undamped free vibration, K &phi; = &omega;&sup2; M &phi;, "
           "with mass-normalised modes. The effective modal mass of mode n "
           "in direction d is &Gamma;&sup2; with &Gamma; = &phi;<sup>T</sup> "
           "M &iota;<sub>d</sub>; over all modes it sums to the free mass in "
           "that direction.</p>"]
    tm = np.asarray(m.total_mass, dtype=float)
    fm = np.asarray(m.free_mass, dtype=float)
    pr = np.asarray(m.participating_ratio, dtype=float) * 100.0
    need = m.modes_for_ratio(0.90)
    more = f"more than {m.n_modes}"
    same = bool(len(tm)) and np.allclose(tm, tm[0], rtol=1e-12, atol=0.0)
    rows = [("Mass model", _e(f"{m.mass_type} mass")),
            ("Modes", _e(f"{m.n_modes} of {m.n_dynamic} possible "
                         f"({m.n_free} free DOF)")),
            ("Solver", _e(m.method or "-")),
            ("Total mass", _e(f"{tm[0]:,.1f} kg" if same else
                              " / ".join(f"{v:,.1f}" for v in tm)
                              + " kg (x / y / z)"))]
    if m.mass_parts:
        rows.append(("Mass sources", _e(", ".join(
            f"{k} {v:,.1f} kg" for k, v in m.mass_parts.items()))))
    rows.append(("Free (participating) mass",
                 _e(" / ".join(f"{v:,.1f}" for v in fm) + " kg (x / y / z)")))
    f1 = m.fundamental_frequency
    if math.isfinite(f1) and f1 > 0.0:
        rows.append(("Fundamental frequency",
                     _e(f"{f1:.4f} Hz (T = {1.0 / f1:.4f} s)")))
    rows.append(("Mass captured by these modes",
                 _e(f"x {pr[0]:.1f} %, y {pr[1]:.1f} %, z {pr[2]:.1f} %")))
    rows.append(("Modes for 90 % of the mass", _e(", ".join(
        f"{ax} {need[ax] if need[ax] is not None else more}"
        for ax in ("x", "y", "z")))))
    if m.prestress:
        rows.append(("Pre-stress (geometric stiffness)", _e(m.prestress)))
    rows.append(("Checks", _e(
        f"|Phi'M Phi - I| = {m.orthogonality_error:.1e}, "
        f"|Phi'K Phi - W| = {m.stiffness_error:.1e}, "
        f"residual = {m.residual:.1e}")))
    out.append(_kv(rows))

    table = []
    for r in m.table():
        f = r["frequency_Hz"]
        table.append({
            "mode": r["mode"],
            "f [Hz]": _fnum(f, "{:.4f}", nan="unstable"),
            "T [s]": _fnum(r["period_s"], "{:.4f}", nan="-"),
            "Mx [%]": _fnum(100 * r["ratio_x"], "{:.2f}"),
            "My [%]": _fnum(100 * r["ratio_y"], "{:.2f}"),
            "Mz [%]": _fnum(100 * r["ratio_z"], "{:.2f}"),
            "sum x [%]": _fnum(100 * r["cum_x"], "{:.2f}"),
            "sum y [%]": _fnum(100 * r["cum_y"], "{:.2f}"),
            "sum z [%]": _fnum(100 * r["cum_z"], "{:.2f}"),
            "direction": r["dominant"],
        })
    out.append("<h3>Frequencies and effective modal mass</h3>")
    out.append(_table(table, numeric_cols=(
        "mode", "f [Hz]", "T [s]", "Mx [%]", "My [%]", "Mz [%]",
        "sum x [%]", "sum y [%]", "sum z [%]"), limit=40))

    model = res.model
    if include_figures and m.n_modes:
        if np.asarray(m.modes).shape[0] == 3 * model.geom.n_nodes:
            k = min(3, m.n_modes)
            out.append(_img(mode_shapes_figure(model, m, k), f"{num}.1",
                            f"The first {k} mode shape(s), each scaled so "
                            f"that its largest node movement is 10 % of the "
                            f"model size; colour = relative amplitude."))
        else:
            out.append("<div class='note'>Mode shapes are not drawn: the "
                       "modal result belongs to a different model.</div>")
    out += _notes(m.warnings)
    out += _notes(m.notes)
    return out


# -------------------------------------------------------- 11 fabrication --
_FAB_STRUTS = [
    ("family", "family", None), ("colour", "colour", None),
    ("count", "count", "int"), ("section", "section", None),
    ("length_mm", "length [mm]", "{:.1f}"),
    ("cut_length_mm", "cut length [mm]", "{:.1f}"),
    ("chord_factor", "chord factor", "{:.6f}"),
    ("axial_angle_deg", "axial angle [deg]", "{:.3f}"),
    ("tangent_angle_deg", "tangent angle [deg]", "{:.3f}"),
    ("mass_each_kg", "mass each [kg]", "{:.2f}"),
    ("mass_total_kg", "mass total [kg]", "{:.1f}"),
]
_FAB_HUBS = [
    ("hub_type", "hub type", None), ("location", "location", None),
    ("valence", "valence", "int"), ("struts", "struts", None),
    ("arrangement", "arrangement", None), ("count", "count", "int"),
    ("gaps_deg", "gaps between struts [deg]", None),
    ("example_node", "example node", "int"),
]
_FAB_PANELS = [
    ("panel_type", "panel type", None), ("count", "count", "int"),
    ("side_1_mm", "side 1 [mm]", "{:.1f}"),
    ("side_2_mm", "side 2 [mm]", "{:.1f}"),
    ("side_3_mm", "side 3 [mm]", "{:.1f}"),
    ("angle_1_deg", "angle 1 [deg]", "{:.2f}"),
    ("angle_2_deg", "angle 2 [deg]", "{:.2f}"),
    ("angle_3_deg", "angle 3 [deg]", "{:.2f}"),
    ("area_m2", "area [m2]", "{:.3f}"),
    ("total_area_m2", "total area [m2]", "{:.2f}"),
    ("handedness", "handedness", None),
]
_FAB_DIHEDRALS = [
    ("strut_family", "strut family", None), ("panels", "panels", None),
    ("count", "count", "int"), ("dihedral_deg", "dihedral [deg]", "{:.3f}"),
    ("bend_deg", "bend [deg]", "{:.3f}"),
]
_FAB_PATTERNS = [
    ("section", "section", None), ("pattern", "pattern", None),
    ("bars", "bars", "int"), ("used_mm", "used [mm]", "{:.1f}"),
    ("offcut_mm", "offcut [mm]", "{:.1f}"),
]


def _section_fabrication(tables, res, num, include_figures):
    if not isinstance(tables, dict):
        raise TypeError("expected the dict of fabrication.fabrication_tables")
    out = ["<p>Lengths are hub centre to hub centre; the cut length "
           "subtracts the hub offset at both ends. Chord factor cf = L/R, "
           "axial angle = arccos(L/2R) (strut to hub radius), tangent angle "
           "= 90&deg; minus the axial angle. Angle k of a panel is opposite "
           "side k.</p>"]
    summary = tables.get("summary") or []
    if summary:
        pairs = []
        for r in summary:
            v = r.get("value", "")
            if isinstance(v, (float, np.floating)):
                v = _fnum(v, "{:.6g}")
            unit = str(r.get("unit", "") or "")
            text = f"{v} {unit}" if unit not in ("", "-") else f"{v}"
            if r.get("note"):
                text += f"  ({r['note']})"
            pairs.append((str(r.get("quantity", "")), _e(text)))
        out.append(_kv(pairs))
    for key, spec, title, empty in (
            ("struts", _FAB_STRUTS, "Strut cut list", "No struts."),
            ("hubs", _FAB_HUBS, "Hub types", "No hubs."),
            ("panels", _FAB_PANELS, "Panel types",
             "This model has no panels."),
            ("dihedrals", _FAB_DIHEDRALS, "Dihedral angles between panels",
             None),
            ("stock_patterns", _FAB_PATTERNS, "Cutting patterns", None)):
        rows = tables.get(key) or []
        if not rows and empty is None:
            continue
        out.append(f"<h3>{_e(title)}</h3>")
        if not rows:
            out.append(f"<p><i>{_e(empty)}</i></p>")
            continue
        shown, numeric = _spec_rows(rows, spec)
        out.append(_table(shown, numeric_cols=numeric, limit=60))
    out += _notes([r.get("note", "") for r in tables.get("notes") or []])
    return out


# --------------------------------------------------------------- 12 study --
_STUDY_HEADS = {
    "frequency": "frequency", "radius_m": "R [m]", "fraction": "fraction",
    "nodes": "nodes", "struts": "struts", "families": "families",
    "mass_kg": "mass [kg]", "mass_per_plan_kg_m2": "mass / plan [kg/m2]",
    "max_util": "max util.", "governing_combination": "governing",
    "lambda_cr": "λ_cr", "span_over_deflection": "span / defl.",
    "status": "status", "time_s": "time [s]",
}


def _study_x(rows):
    for key, xk in (("frequency", "frequency"), ("radius_m", "radius_m"),
                    ("fraction_requested", "fraction")):
        if len({r.get(key) for r in rows}) > 1:
            return xk
    return "frequency"


def study_figure(rows):
    """The study charts of ``study.plot_study`` (mass, utilisation,
    lambda_cr, span/deflection against the swept parameter)."""
    from . import study as st
    ys = ["mass_kg", "max_util", "lambda_cr", "span_over_deflection"]
    good = [r for r in rows if r.get("status") != "ERROR"]
    lam = [r.get("lambda_cr") for r in good]
    if not any(v is not None and math.isfinite(float(v)) for v in lam):
        ys.remove("lambda_cr")
    fig, _ = st.plot_study(None, rows, x=_study_x(rows), ys=tuple(ys))
    fig.patch.set_facecolor("white")
    return fig


def _section_study(rows, res, num, include_figures):
    from . import study as st
    source = rows                       # StudyRows keeps .settings/.messages
    rows = list(rows)
    if not all(isinstance(r, dict) for r in rows):
        raise TypeError("expected the rows of study.parameter_study")
    settings = getattr(source, "settings", None) or {}
    out = ["<p>Each row is a regenerated dome analysed with the same "
           "loads, combinations and design checks as the main model. "
           "Rows in red fail a check (status other than OK).</p>"]
    pairs = []
    for key, label in (("frequencies", "Frequencies"), ("radii", "Radii [m]"),
                       ("fractions", "Dome fractions")):
        if settings.get(key) is not None:
            pairs.append((label, _e(", ".join(f"{v:g}" for v in
                                              settings[key]))))
    if settings.get("section"):
        pairs.append(("Start section", _e(settings["section"])))
    if "auto_size" in settings:
        pairs.append(("Auto-sizing", _e("yes" if settings["auto_size"]
                                        else "no")))
    best = st.best_design(rows)
    if best is not None:
        pairs.append(("Lightest design that passes every check", _e(
            f"{int(best['frequency'])}V, R = {float(best['radius_m']):g} m, "
            f"fraction {float(best['fraction']):.3f}, "
            f"{float(best['mass_kg']):,.0f} kg")))
    else:
        pairs.append(("Lightest design that passes every check",
                      "<span class='bad'>none</span>"))
    out.append(_kv(pairs))

    cols = st.DEFAULT_TABLE_COLUMNS
    table = [{_STUDY_HEADS.get(k, k): st._fmt(k, r.get(k)) for k in cols}
             for r in rows]
    numeric = tuple(_STUDY_HEADS.get(k, k) for k in cols if k not in
                    ("governing_combination", "status"))
    out.append(_table(table, numeric_cols=numeric,
                      highlight=lambda d: d.get("status") != "OK"))
    if include_figures and any(r.get("status") != "ERROR" for r in rows):
        out.append(_img(study_figure(rows), f"{num}.1",
                        "Study results against the swept parameter. Filled "
                        "markers pass every check; hollow markers fail "
                        "one."))
    return out + _notes(getattr(source, "messages", None))
