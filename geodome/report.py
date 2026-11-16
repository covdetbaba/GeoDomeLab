"""
Export: CSV tables, a self-contained HTML calculation report, and PNG figures.

The HTML report embeds its figures as base64 PNGs so that a single file can be
submitted or e-mailed with nothing else attached.
"""

from __future__ import annotations

import base64
import csv
import datetime as _dt
import html
import io
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
                     "chord_factor": round(s["chord_factor"], 6),
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
                author="", include_figures=True, verification=None):
    """Write a single self-contained HTML calculation report."""
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
        [{"family": s["label"], "chord factor": round(s["chord_factor"], 6),
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
        parts.append(_table(
            [{"mode": i + 1, "critical load factor": round(float(l), 3)}
             for i, l in enumerate(res.buckling.load_factors)],
            numeric_cols=("critical load factor",)))
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

    parts.append("<footer>Produced by GeoDomeLab. Linear elastic first-order "
                 "space-truss analysis; member checks to AISC 360 LRFD; "
                 "wind pressures from a simplified domed-roof coefficient "
                 "distribution. Results are for academic use.</footer>")
    parts.append("</div></body></html>")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(parts))
    return path
