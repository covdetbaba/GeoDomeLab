"""
Matplotlib 3D renderers.

Every function here draws onto an existing ``Axes3D`` and never touches
``pyplot``, so the same code serves the embedded Qt canvases, the headless
figure export used by the report writer, and Spyder's inline backend without
any backend conflicts.
"""

from __future__ import annotations

import numpy as np
from matplotlib.collections import LineCollection  # noqa: F401  (type hints)
from matplotlib.colors import Normalize, TwoSlopeNorm, ListedColormap
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

FAMILY_COLORS = ["#1f6fb4", "#d1495b", "#2a9d8f", "#e9a000", "#7b53a6",
                 "#4f772d", "#c1666b", "#00798c", "#b56576", "#495057",
                 "#8ac926", "#ff924c"]


# ------------------------------------------------------------------ frame ---
def _segments(nodes, members):
    return np.stack([nodes[members[:, 0]], nodes[members[:, 1]]], axis=1)


def set_equal_aspect(ax, pts, pad=0.05):
    """Give the 3D axes a true 1:1:1 scale without wasting frame height.

    The x and y limits are padded symmetrically, the z limits are clipped just
    below ground level, and the drawing box is then made *proportional to the
    data ranges* rather than cubic.  A hemisphere of radius R therefore fills a
    2R x 2R x R box instead of floating in the middle of a 2R cube, while one
    metre still measures the same along every axis.
    """
    pts = np.asarray(pts, dtype=float)
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    c = 0.5 * (lo + hi)

    rxy = 0.5 * max(hi[0] - lo[0], hi[1] - lo[1]) * (1.0 + pad)
    rxy = max(rxy, 1e-9)
    z_lo = min(lo[2], 0.0) - pad * rxy
    z_hi = hi[2] + pad * rxy
    if z_hi - z_lo < 1e-9:
        z_hi = z_lo + 1.0

    ax.set_xlim(c[0] - rxy, c[0] + rxy)
    ax.set_ylim(c[1] - rxy, c[1] + rxy)
    ax.set_zlim(z_lo, z_hi)
    try:
        ax.set_box_aspect((2.0 * rxy, 2.0 * rxy, z_hi - z_lo))
    except Exception:                                # pragma: no cover
        pass


def prepare_axes(ax, title="", elev=22.0, azim=-58.0, grid=True):
    ax.clear()
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel("x [m]", fontsize=9, labelpad=2)
    ax.set_ylabel("y [m]", fontsize=9, labelpad=2)
    ax.set_zlabel("z [m]", fontsize=9, labelpad=2)
    ax.tick_params(labelsize=8)
    ax.view_init(elev=elev, azim=azim)
    ax.grid(grid)
    try:
        ax.xaxis.pane.set_alpha(0.03)
        ax.yaxis.pane.set_alpha(0.03)
        ax.zaxis.pane.set_alpha(0.06)
    except Exception:                                # pragma: no cover
        pass
    return ax


def _draw_ground(ax, geom):
    """Faint circle marking the base ring footprint."""
    R = geom.base_diameter / 2.0
    t = np.linspace(0, 2 * np.pi, 180)
    ax.plot(R * np.cos(t), R * np.sin(t), np.zeros_like(t),
            color="0.55", lw=0.8, ls="--", zorder=0)


def _draw_supports(ax, geom, size=28):
    b = geom.nodes[geom.base_nodes]
    ax.scatter(b[:, 0], b[:, 1], b[:, 2], marker="^", s=size,
               c="#111111", depthshade=False, zorder=6, label="Support")


# ------------------------------------------------------------- model view ---
def draw_model(ax, model, color_by_family=True, show_nodes=True,
               show_panels=False, show_supports=True, show_labels=False,
               title=None):
    geom = model.geom
    if title is None and getattr(geom, "kind", "geodesic") == "custom":
        title = f"{geom.title or 'Custom truss'}  --  {geom.n_members} struts"
    prepare_axes(ax, title or f"Geodesic dome model  --  {geom.frequency}V, "
                              f"R = {geom.radius:.1f} m")
    segs = _segments(geom.nodes, geom.members)

    if show_panels:
        polys = [geom.nodes[f] for f in geom.faces]
        pc = Poly3DCollection(polys, facecolor="#9ec9e8", alpha=0.18,
                              edgecolor="none")
        ax.add_collection3d(pc)

    if color_by_family:
        cols = [FAMILY_COLORS[g % len(FAMILY_COLORS)] for g in geom.strut_group]
        lc = Line3DCollection(segs, colors=cols, linewidths=1.6)
    else:
        lc = Line3DCollection(segs, colors="#2b4c6f", linewidths=1.4)
    ax.add_collection3d(lc)

    if show_nodes:
        ax.scatter(geom.nodes[:, 0], geom.nodes[:, 1], geom.nodes[:, 2],
                   s=7, c="#1b1b1b", depthshade=False, zorder=5)
    if show_supports:
        _draw_supports(ax, geom)
    if show_labels and geom.n_nodes <= 120:
        for k, p in enumerate(geom.nodes):
            ax.text(p[0], p[1], p[2], str(k), fontsize=6, color="#444444")

    _draw_ground(ax, geom)
    set_equal_aspect(ax, geom.nodes)

    if color_by_family:
        handles = [Line2D([0], [0], color=FAMILY_COLORS[s["id"] % len(FAMILY_COLORS)],
                          lw=2.4,
                          label=f"{s['label']}  {s['length']:.3f} m  "
                                + (f"(CF {s['chord_factor']:.4f})  "
                                   if np.isfinite(s['chord_factor']) else "")
                                + f"x{s['count']}")
                   for s in geom.strut_schedule()]
        ax.legend(handles=handles, loc="upper left", fontsize=7,
                  frameon=True, framealpha=0.9, title="Strut families",
                  title_fontsize=8)
    return ax


# -------------------------------------------------------------- load view ---
def draw_loads(ax, model, F, title="Applied nodal loads", arrow_frac=0.09):
    geom = model.geom
    prepare_axes(ax, title)
    segs = _segments(geom.nodes, geom.members)
    ax.add_collection3d(Line3DCollection(segs, colors="#b9c6d2", linewidths=0.8))

    V = np.asarray(F).reshape(-1, 3)
    mag = np.linalg.norm(V, axis=1)
    mmax = mag.max()
    if mmax > 1e-9:
        scale = arrow_frac * geom.length_scale / mmax
        ax.quiver(geom.nodes[:, 0], geom.nodes[:, 1], geom.nodes[:, 2],
                  V[:, 0] * scale, V[:, 1] * scale, V[:, 2] * scale,
                  color="#c1121f", linewidth=0.9, arrow_length_ratio=0.28)
        ax.text2D(0.02, 0.96,
                  f"max nodal load  {mmax / 1e3:.2f} kN", transform=ax.transAxes,
                  fontsize=8, va="top",
                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                            alpha=0.85))
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, geom.nodes)
    return ax


# --------------------------------------------------------- deflected view ---
def auto_scale(geom, u, target_frac=0.05):
    """Magnification that makes the largest displacement ``target_frac*R``."""
    dmax = np.max(np.linalg.norm(np.asarray(u).reshape(-1, 3), axis=1))
    if dmax < 1e-12:
        return 1.0
    return float(target_frac * geom.length_scale / dmax)


def draw_deflected(ax, model, result, scale=None, show_undeformed=True,
                   title=None):
    geom = model.geom
    if scale is None:
        scale = auto_scale(geom, result.u)
    disp = result.displacement_vectors()
    defo = geom.nodes + scale * disp
    mag = np.linalg.norm(disp, axis=1)

    prepare_axes(ax, title or f"Deflected shape  --  {result.name}")
    if show_undeformed:
        ax.add_collection3d(Line3DCollection(
            _segments(geom.nodes, geom.members),
            colors="#c8cfd6", linewidths=0.7))

    segs = _segments(defo, geom.members)
    vals = 0.5 * (mag[geom.members[:, 0]] + mag[geom.members[:, 1]]) * 1e3
    lc = Line3DCollection(segs, linewidths=1.8, cmap="viridis",
                          norm=Normalize(0.0, max(vals.max(), 1e-9)))
    lc.set_array(vals)
    ax.add_collection3d(lc)
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, np.vstack([geom.nodes, defo]))

    ax.text2D(0.02, 0.96,
              f"magnification x{scale:,.0f}\nmax |u| = {mag.max() * 1e3:.2f} mm",
              transform=ax.transAxes, fontsize=8, va="top",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                        alpha=0.85))
    return lc


# ------------------------------------------------------- axial force view ---
def draw_axial(ax, model, N, title=None, width_scale=True, deformed=None):
    geom = model.geom
    nodes = geom.nodes if deformed is None else deformed
    prepare_axes(ax, title or "Member axial force  (red = tension, "
                              "blue = compression)")
    segs = _segments(nodes, geom.members)
    kN = np.asarray(N) / 1e3
    amax = max(abs(kN).max(), 1e-9)
    norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
    lw = 0.6 + 2.6 * (np.abs(kN) / amax) if width_scale else np.full(len(kN), 1.6)

    lc = Line3DCollection(segs, cmap="coolwarm", norm=norm, linewidths=lw)
    lc.set_array(kN)
    ax.add_collection3d(lc)
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, nodes)

    it = int(np.argmax(kN))
    ic = int(np.argmin(kN))
    ax.text2D(0.02, 0.96,
              f"max tension      {kN[it]:+.1f} kN  (member {it})\n"
              f"max compression  {kN[ic]:+.1f} kN  (member {ic})",
              transform=ax.transAxes, fontsize=8, va="top", family="monospace",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                        alpha=0.85))
    return lc


# ------------------------------------------------------- utilisation view ---
_UTIL_COLORS = ["#2a9d8f", "#8bc34a", "#ffd166", "#f4a261", "#e76f51", "#9d0208"]
_UTIL_BOUNDS = [0.0, 0.25, 0.50, 0.75, 0.90, 1.00, 1e9]


def draw_utilisation(ax, model, util, title=None):
    geom = model.geom
    prepare_axes(ax, title or "Demand / capacity ratio  (AISC LRFD envelope)")
    util = np.asarray(util)
    idx = np.clip(np.searchsorted(_UTIL_BOUNDS, util, side="right") - 1,
                  0, len(_UTIL_COLORS) - 1)
    cols = [_UTIL_COLORS[i] for i in idx]
    lw = np.where(util > 1.0, 3.2, 1.6)
    ax.add_collection3d(Line3DCollection(_segments(geom.nodes, geom.members),
                                         colors=cols, linewidths=lw))
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, geom.nodes)

    labels = ["< 0.25", "0.25-0.50", "0.50-0.75", "0.75-0.90",
              "0.90-1.00", "> 1.00  OVER"]
    handles = [Line2D([0], [0], color=c, lw=3, label=l)
               for c, l in zip(_UTIL_COLORS, labels)]
    ax.legend(handles=handles, loc="upper left", fontsize=7, frameon=True,
              framealpha=0.9, title="Utilisation", title_fontsize=8)
    k = int(np.argmax(util))
    ax.text2D(0.98, 0.96,
              f"max utilisation {util[k]:.3f} at member {k}\n"
              f"{int(np.sum(util > 1.0))} strut(s) over 1.0",
              transform=ax.transAxes, fontsize=8, va="top", ha="right",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                        alpha=0.85))
    return ax


# ---------------------------------------------------------- buckling view ---
def draw_buckling_mode(ax, model, bres, k=0, scale_frac=0.10, title=None):
    geom = model.geom
    if bres is None or not len(bres.load_factors):
        prepare_axes(ax, "Buckling mode -- not available")
        ax.text2D(0.5, 0.5, "No positive buckling mode for this load case",
                  transform=ax.transAxes, ha="center", fontsize=10)
        return ax
    k = int(np.clip(k, 0, len(bres.load_factors) - 1))
    phi = bres.mode(k)
    pmax = np.max(np.linalg.norm(phi, axis=1))
    scale = scale_frac * geom.length_scale / max(pmax, 1e-12)
    defo = geom.nodes + scale * phi
    mag = np.linalg.norm(phi, axis=1)

    prepare_axes(ax, title or f"Buckling mode {k + 1}   "
                              f"lambda_cr = {bres.load_factors[k]:.2f}")
    ax.add_collection3d(Line3DCollection(_segments(geom.nodes, geom.members),
                                         colors="#ccd3da", linewidths=0.7))
    segs = _segments(defo, geom.members)
    vals = 0.5 * (mag[geom.members[:, 0]] + mag[geom.members[:, 1]])
    lc = Line3DCollection(segs, cmap="magma_r", linewidths=1.9,
                          norm=Normalize(0.0, max(vals.max(), 1e-12)))
    lc.set_array(vals)
    ax.add_collection3d(lc)
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, np.vstack([geom.nodes, defo]))
    ax.text2D(0.02, 0.96,
              f"critical combination: {bres.reference_case}\n"
              f"mode shape normalised, shown x{scale:,.1f}",
              transform=ax.transAxes, fontsize=8, va="top",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                        alpha=0.85))
    return lc


# ------------------------------------------------------------- 2D helpers ---
def draw_strut_histogram(ax, model):
    """Bar chart of the strut schedule -- a fabrication deliverable."""
    sched = model.geom.strut_schedule()
    ax.clear()
    labels = [s["label"] for s in sched]
    counts = [s["count"] for s in sched]
    lengths = [s["length"] for s in sched]
    cols = [FAMILY_COLORS[s["id"] % len(FAMILY_COLORS)] for s in sched]
    bars = ax.bar(labels, counts, color=cols, edgecolor="#333333", linewidth=0.6)
    for b, L in zip(bars, lengths):
        ax.annotate(f"{L:.3f} m", (b.get_x() + b.get_width() / 2, b.get_height()),
                    ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Strut family")
    ax.set_ylabel("Number of struts")
    ax.set_title("Strut schedule", fontsize=11, fontweight="bold")
    ax.margins(y=0.18)
    ax.grid(axis="y", alpha=0.3)
    return ax


def draw_utilisation_histogram(ax, res):
    ax.clear()
    u = res.util_env
    ax.hist(u, bins=24, color="#1f6fb4", edgecolor="white")
    ax.axvline(1.0, color="#9d0208", lw=1.6, ls="--", label="capacity")
    ax.set_xlabel("Demand / capacity ratio")
    ax.set_ylabel("Number of struts")
    ax.set_title("Utilisation distribution", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    return ax
