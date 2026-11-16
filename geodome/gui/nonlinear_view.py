"""
Nonlinear analysis view: equilibrium path, path scrubbing and the
imperfection study (the GUI for :mod:`geodome.nonlinear`).

What is plotted
---------------
The load is ``lambda * F_ref`` where ``F_ref`` is the D/S/L/W part of one load
combination (the self-straining load T is not part of the nonlinear path).
Three answers are drawn against the displacement ``u_c`` of the control DOF:

* corotational, large displacement (``nonlinear.run_nonlinear``):
      f(u) = lambda F_ref, traced by arc-length, displacement or load control;
* linearised second order (P-Delta style):
      u = lambda (K + lambda K_g,ref)^-1 F_ref     (singular at lambda_cr);
* first order:
      u = lambda K^-1 F_ref.

``lambda_cr`` (linear buckling) is the dashed horizontal line.  The first load
maximum ``lambda_limit`` is marked; ``lambda_limit / lambda_cr`` is the
knock-down that geometric nonlinearity alone gives.  A slider scrubs along the
path and draws the deformed shape (or the member forces) at that step.

Imperfection study (``nonlinear.imperfection_study``): node coordinates
``X + a phi``, with ``phi`` a linear buckling mode scaled to
max |component| = 1.  Both signs are run and the lower limit load is kept.  The chart is
``lambda_limit / lambda_cr`` against ``a / span``.  The default amplitudes are
span/1000, span/500, span/300 and span/200.

Integration
-----------
``install(win) -> NonlinearExtension``.  It adds, through the public
``MainWindow`` API only:

* Analysis menu: "Nonlinear analysis..." and "Imperfection study...";
* toolbar button "Nonlinear" (after "Auto-size");
* view tab "Nonlinear" (sub-tabs "Load path" and "Imperfection study"); a
  window whose View menu lists the tabs itself (``win.menu_view``) shows it
  there, an older window gets a fixed View menu entry;
* result tabs "Nonlinear" and "Imperfection" (before "Log");
* a Stop button in the status bar, shown while a run is going.

The extension keeps ``nonlinear_result`` (a ``NonlinearResult`` or None),
``imperfection_result`` (the list of row dicts returned by
``imperfection_study``, or None), ``imperfection_runs``
(``{(a, sign): NonlinearResult}``) and ``imperfection_settings`` (dict).
``on_model_changed(model)`` clears them; ``on_results(results)`` does nothing.
Long runs go through :class:`gui.worker.Worker`; ``wait()`` blocks (with a
local event loop) until the current run ends, for scripts and tests.
``install`` registers the extension in ``win.extensions`` when that list
exists, and calling it twice returns the same extension.
"""

from __future__ import annotations

import copy
import csv
import dataclasses
import math
import re
import threading
import traceback

import numpy as np
from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Qt, Signal

from .. import loads as ld
from .. import nonlinear as nl
from .. import plotting as plot
from . import style as st
from .canvas import MplCanvas
from .worker import Worker, plain_message

__all__ = ["install", "NonlinearExtension", "NonlinearDialog",
           "ImperfectionDialog", "default_combination", "parse_amplitudes",
           "nonlinear_job", "imperfection_job", "draw_path",
           "draw_imperfection_ratio", "draw_imperfection_paths",
           "StopRequested"]

METHODS = (("arc-length", "Arc-length (Crisfield)"),
           ("displacement", "Displacement control"),
           ("load", "Load control (Newton-Raphson)"))
AXES = ("ux", "uy", "uz")
DEFAULT_AMPLITUDES_TEXT = "span/1000, span/500, span/300, span/200"
MAX_AMPLITUDES = 12
OBJECT_NAME = "nonlinear_extension"

COL_PATH = st.ACCENT
COL_PDELTA = "#e9a000"
COL_LINEAR = "#6b7785"
COL_CR = st.BAD_RED
COL_EULER = "#7b53a6"
COL_PLUS = "#e9a000"
COL_MINUS = "#2a9d8f"


# ====================================================== core helpers (no Qt) ==
class StopRequested(Exception):
    """Raised inside a running job when the user presses Stop."""


def default_combination(model):
    """The combination ``run_nonlinear`` takes when none is named: the
    strength combination with the largest total vertical load."""
    combos = list(model.combinations)
    if not combos:
        raise ValueError("The model has no load combinations.")
    basic = ld.basic_load_vectors(model)
    strength = [c for c in combos if not c.serviceability] or combos
    return max(strength,
               key=lambda c: abs(float(ld.combine(basic, c)[2::3].sum())))


def find_combination(model, name):
    """Combination by exact name (None = :func:`default_combination`)."""
    if name is None:
        return default_combination(model)
    for c in model.combinations:
        if c.name == name:
            return c
    raise ValueError(f"Combination '{name}' is not in the model.")


def model_span(model):
    """Span used for the imperfection amplitudes [m] (NaN if unknown)."""
    try:
        s = float(model.geom.base_diameter)
    except Exception:
        return float("nan")
    return s if math.isfinite(s) and s > 0.0 else float("nan")


def mode_shape(model, combination, mode=0):
    """Buckling mode ``mode`` (0-based) as an (n, 3) array scaled to
    max |component| = 1, with lambda_cr and F_ref of the combination."""
    F = nl.reference_load(model, combination)
    lam_cr, _, b = nl.linear_buckling_reference(model, F)
    if not len(b.load_factors):
        raise ValueError("This combination does not compress the structure, "
                         "so there is no buckling mode to use as an "
                         "imperfection.")
    if mode >= b.modes.shape[1]:
        raise ValueError(f"Only {b.modes.shape[1]} buckling modes are "
                         f"available. Pick a lower mode number.")
    phi = np.array(b.modes[:, mode], dtype=float)
    phi /= max(float(np.max(np.abs(phi))), 1e-300)
    return phi.reshape(-1, 3), lam_cr, F


def offset_model(model, offset):
    """Shallow copy of ``model`` with node coordinates ``X + offset`` (the
    same construction ``run_nonlinear`` uses for ``node_offset``)."""
    if offset is None:
        return model
    work = copy.copy(model)
    work.geom = dataclasses.replace(
        model.geom, nodes=np.asarray(model.geom.nodes, float) +
        np.asarray(offset, float).reshape(-1, 3))
    return work


_SPAN_TOKEN = re.compile(r"^(?:span|l|d)?/(\d+(?:\.\d*)?|\.\d+)$", re.I)
_MM_TOKEN = re.compile(r"^(\d+(?:\.\d*)?|\.\d+)(?:mm)?$", re.I)


def parse_amplitudes(text, span):
    """Read an amplitude list such as ``"span/1000, span/500, 25"``.

    ``span/n`` (also ``L/n`` or ``/n``) means span / n; a plain number is in
    mm.  Returns ``[(amplitude_m, label), ...]`` sorted and without repeats.
    Raises ValueError with a short message on bad input.
    """
    tokens = [t.replace(" ", "") for t in re.split(r"[,;]", str(text))]
    tokens = [t for t in tokens if t]
    if not tokens:
        raise ValueError("Enter at least one amplitude, for example "
                         "span/1000.")
    out = {}
    for tok in tokens:
        m = _SPAN_TOKEN.match(tok)
        if m:
            n = float(m.group(1))
            if not (math.isfinite(span) and span > 0.0):
                raise ValueError("This model has no span, so span/n cannot "
                                 "be used. Enter the amplitudes in mm.")
            if n <= 0.0:
                raise ValueError(f"'{tok}': the divisor must be positive.")
            a = span / n
            label = f"span/{n:g}"
        else:
            m = _MM_TOKEN.match(tok)
            if not m:
                raise ValueError(f"Cannot read '{tok}'. Use span/1000 or a "
                                 f"value in mm.")
            a = float(m.group(1)) / 1e3
            label = f"{float(m.group(1)):g} mm"
        if not (a > 0.0 and math.isfinite(a)):
            raise ValueError(f"'{tok}': amplitudes must be positive.")
        key = round(a, 12)
        out.setdefault(key, (a, label))
    if len(out) > MAX_AMPLITUDES:
        raise ValueError(f"Use at most {MAX_AMPLITUDES} amplitudes "
                         f"(each one runs twice).")
    return [out[k] for k in sorted(out)]


def _checked_progress(progress, cancel):
    def prog(fraction, message=""):
        if cancel is not None and cancel.is_set():
            raise StopRequested("Stopped by the user.")
        if progress is not None:
            progress(fraction, message)
    return prog


def nonlinear_job(model, settings, cancel=None, progress=None):
    """One nonlinear run (called in the worker thread; no Qt here).

    ``settings`` keys: combination (name or None), method, max_steps,
    stop_after_limit, lambda_max (None = auto), first_increment (None =
    auto), control_node, control_axis, imperfection_mode (1-based),
    imperfection_mm (0 = perfect).  Returns a dict with the result, the
    (possibly imperfect) model to draw, the node offset and the settings.
    """
    prog = _checked_progress(progress, cancel)
    combo = find_combination(model, settings.get("combination"))
    node = settings.get("control_node")
    control = "crown_uz" if node is None else \
        (int(node), str(settings.get("control_axis", "uz")))
    kw = dict(method=settings.get("method", "arc-length"), control=control,
              max_steps=int(settings.get("max_steps", 300)),
              stop_after_limit=settings.get("stop_after_limit", 5))
    if settings.get("first_increment"):
        kw["first_increment"] = float(settings["first_increment"])
    if settings.get("lambda_max"):
        kw["lambda_max"] = float(settings["lambda_max"])
    offset = None
    amp = float(settings.get("imperfection_mm") or 0.0) / 1e3
    if amp != 0.0:
        prog(0.01, "Buckling mode for the imperfection...")
        phi, lam_cr, F = mode_shape(
            model, combo, int(settings.get("imperfection_mode") or 1) - 1)
        offset = amp * phi
        kw.update(node_offset=offset, load_vector=F, lambda_cr=lam_cr)
    prog(0.0, "Nonlinear analysis...")
    res = nl.run_nonlinear(model, combo, progress=prog, **kw)
    return {"result": res, "model": offset_model(model, offset),
            "offset": offset, "settings": dict(settings,
                                               combination=combo.name)}


def imperfection_job(model, settings, cancel=None, progress=None):
    """Imperfection study (worker thread).  ``settings`` keys: combination,
    mode (1-based), amplitudes [(a_m, label)], method, max_steps,
    stop_after_limit."""
    prog = _checked_progress(progress, cancel)
    combo = find_combination(model, settings.get("combination"))
    amps = [float(a) for a, _ in settings["amplitudes"]]
    kw = dict(method=settings.get("method", "arc-length"),
              max_steps=int(settings.get("max_steps", 300)),
              stop_after_limit=settings.get("stop_after_limit", 5),
              control="crown_uz")
    rows, runs = nl.imperfection_study(
        model, combo, amplitudes=amps, mode=int(settings.get("mode", 1)) - 1,
        return_results=True, progress=prog, **kw)
    return {"rows": rows, "runs": runs, "model": model,
            "settings": dict(settings, combination=combo.name),
            "span": model_span(model)}


# ============================================================ drawing ======
def _finite(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _g(x, fmt=".4g"):
    return format(float(x), fmt) if _finite(x) else "-"


def _sign_cell(lam, amplitude):
    """lambda_limit of one imperfection sign: 'failed' when that run found
    no limit point (the perfect row has only one sign: '-')."""
    if _finite(lam):
        return format(float(lam), ".5g")
    return "failed" if amplitude > 0.0 else "-"


def path_sign(res):
    """+1 or -1 so that ``sign * control`` grows with lambda at the start
    (a downward crown deflection is then plotted positive)."""
    lin = getattr(res, "linear", None) or {}
    c = np.asarray(lin.get("control", []), dtype=float)
    if c.size > 1 and c[-1] != c[0]:
        return -1.0 if c[-1] < c[0] else 1.0
    for v in np.asarray(res.control, dtype=float)[1:]:
        if v != 0.0:
            return -1.0 if v < 0.0 else 1.0
    return 1.0


def control_axis_label(res, sign):
    lab = res.control_label or "control DOF"
    if sign < 0:
        side = ", downward positive" if " uz" in f" {lab}" else \
            ", sign reversed"
    else:
        side = ""
    return f"Displacement of {lab}{side} [mm]"


def draw_path(ax, res, step=None, title=None, sign=None):
    """lambda vs control displacement: corotational, P-Delta and linear
    paths, lambda_cr, the first strut at its Euler load, bifurcations and the
    limit point.  Returns ``(marker, sign)``; ``marker`` is the Line2D of
    the current step (None when ``step`` is None)."""
    ax.clear()
    s = path_sign(res) if sign is None else sign
    x = s * np.asarray(res.control, float) * 1e3
    y = np.asarray(res.lam, float)
    curves = res.curves()
    if "linear" in curves and len(curves["linear"][0]) > 1:
        cx, cl = curves["linear"]
        ax.plot(s * cx * 1e3, cl, color=COL_LINEAR, ls=":", lw=1.6,
                label="Linear (first order)")
    if "p_delta" in curves and len(curves["p_delta"][0]) > 1:
        cx, cl = curves["p_delta"]
        ax.plot(s * cx * 1e3, cl, color=COL_PDELTA, ls="--", lw=1.6,
                label="P-Delta (linearised second order)")
    ax.plot(x, y, color=COL_PATH, lw=2.0, marker="o", ms=3.0, zorder=4,
            label="Corotational (large displacement)")

    xmin, xmax = min(0.0, float(x.min())), max(0.0, float(x.max()))
    xr = (xmax - xmin) or 1.0
    ax.set_xlim(xmin - 0.04 * xr, xmax + 0.08 * xr)
    ylo, yhi = min(0.0, float(y.min())), max(0.0, float(y.max()))
    lam_cr = float(res.lambda_cr) if _finite(res.lambda_cr) else float("nan")
    has_cr = _finite(lam_cr) and lam_cr > 0.0
    cr_in = has_cr and lam_cr <= 4.0 * max(yhi, 1e-300)
    if cr_in:
        yhi = max(yhi, lam_cr)
    yr = (yhi - ylo) or 1.0
    ax.set_ylim(ylo - 0.05 * yr, yhi + 0.14 * yr)
    ax.axhline(0.0, color="0.65", lw=0.8, zorder=1)
    if has_cr:
        lab = rf"$\lambda_{{cr}}$ = {lam_cr:.4g} (linear buckling)"
        if cr_in:
            ax.axhline(lam_cr, color=COL_CR, ls="--", lw=1.3, label=lab)
        else:
            ax.plot([], [], color=COL_CR, ls="--", lw=1.3,
                    label=lab + ", above the plot")
    le = getattr(res, "lambda_member_euler", float("nan"))
    if _finite(le) and ylo <= le <= yhi * 1.05:
        ax.axhline(le, color=COL_EULER, ls=":", lw=1.2,
                   label=rf"First strut at its Euler load, "
                         rf"$\lambda$ = {le:.3g}")

    first = True
    for cp in res.critical_points:
        if cp.get("type") == "bifurcation" and not cp.get("first_peak"):
            ax.plot(s * cp["control"] * 1e3, cp["lambda"], "D", color=COL_EULER,
                    ms=7, ls="none", zorder=6,
                    label="Bifurcation (stiffness sign change)" if first
                    else "_nolegend_")
            first = False

    if res.limit_found and _finite(res.lambda_limit):
        xl, yl = s * res.control_limit * 1e3, float(res.lambda_limit)
        kind = "Bifurcation" if res.limit_type == "bifurcation" \
            else "Limit point"
        ax.plot(xl, yl, marker="*", ms=16, color=COL_CR, mec="white",
                mew=0.8, ls="none", zorder=7, label="_nolegend_")
        txt = f"{kind}\n" + rf"$\lambda$ = {yl:.4g}"
        if _finite(res.limit_ratio):
            txt += rf" = {res.limit_ratio:.3f} $\lambda_{{cr}}$"
        _room_above(ax, yl, 54.0)
        right = (xl - xmin) > 0.6 * xr
        label = ax.annotate(
            txt, (xl, yl), xytext=(-14 if right else 14, 12),
            textcoords="offset points", fontsize=8.5,
            ha="right" if right else "left", va="bottom",
            color=COL_CR, fontweight="bold", zorder=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec=COL_CR, lw=0.8, alpha=0.92),
            arrowprops=dict(arrowstyle="-", color=COL_CR, lw=0.8))
    else:
        label, right = None, False

    marker = None
    if step is not None and len(x):
        k = int(np.clip(step, 0, len(x) - 1))
        marker, = ax.plot([x[k]], [y[k]], "o", ms=11, mfc="none",
                          mec="#111111", mew=1.8, zorder=8,
                          label="_nolegend_")
    ax.set_xlabel(control_axis_label(res, s), fontsize=9)
    ax.set_ylabel(r"Load factor $\lambda$", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title(title or f"Equilibrium path  --  {res.combination or ''}",
                 fontsize=10.5, fontweight="bold")
    ax.grid(alpha=0.3)
    legend = ax.legend(fontsize=7.5, loc="best", framealpha=0.92)
    if label is not None:
        _place_label(ax, label, legend, right)
        # kept so that a resize of the canvas can place it again
        ax._limit_label = (label, legend, right, yl)
    else:
        ax._limit_label = None
    return marker, s


def replace_limit_label(ax):
    """Place the limit-point label of :func:`draw_path` again (after the
    canvas changed size: the best spot depends on the pixel size)."""
    info = getattr(ax, "_limit_label", None)
    if not info:
        return
    label, legend, right, yl = info
    _room_above(ax, yl, 54.0)
    _place_label(ax, label, legend, right)


# offsets [points] tried for the limit-point label, best first
_LABEL_SPOTS = [(14, 12, "left", "bottom"), (-14, 12, "right", "bottom"),
                (20, -2, "left", "center"), (-20, -2, "right", "center"),
                (14, -14, "left", "top"), (-14, -14, "right", "top"),
                (14, 42, "left", "bottom"), (-14, 42, "right", "bottom"),
                (14, -46, "left", "top"), (-14, -46, "right", "top"),
                (0, -60, "center", "top"), (0, 36, "center", "bottom")]


def _place_label(ax, label, legend, right=False):
    """Put the limit-point label where it covers least: it must stay
    inside the axes, off the legend and off the plotted lines.

    The legend's "best" place is found first, without the label; each spot
    in :data:`_LABEL_SPOTS` is then scored in screen pixels and the best one
    is used.  The legend keeps "best": at the next draw it also avoids the
    label (matplotlib counts text as an obstacle), so it stays put or moves
    away after a resize.  Without a renderer nothing changes.
    """
    try:
        renderer = ax.figure.canvas.get_renderer()
    except Exception:
        return
    try:
        label.set_visible(False)
        leg_box = None
        if legend is not None:
            leg_box = legend.get_window_extent(renderer)
        label.set_visible(True)
        box_ax = ax.get_window_extent(renderer)
        pts = []
        for ln in ax.get_lines():
            if not ln.get_visible() or ln is label:
                continue
            xy = ln.get_transform().transform(ln.get_xydata())
            if len(xy) > 1:                  # also points along each segment
                t = np.linspace(0.0, 1.0, 9)[:, None, None]
                xy = (xy[:-1][None] * (1 - t) + xy[1:][None] * t
                      ).reshape(-1, 2)
            pts.append(xy)
        pts = np.concatenate(pts) if pts else np.zeros((0, 2))
        pts = pts[np.all(np.isfinite(pts), axis=1)]

        def overlap(a, b):
            w = min(a.x1, b.x1) - max(a.x0, b.x0)
            h = min(a.y1, b.y1) - max(a.y0, b.y0)
            return max(w, 0.0) * max(h, 0.0)

        spots = list(_LABEL_SPOTS)
        if right:                            # keep the old preference first
            spots[0], spots[1] = spots[1], spots[0]
        best = None
        for k, (dx, dy, ha, va) in enumerate(spots):
            label.xyann = (dx, dy)
            label.set_ha(ha)
            label.set_va(va)
            bb = label.get_window_extent(renderer).padded(4.0)
            outside = bb.width * bb.height - overlap(bb, box_ax)
            on_legend = overlap(bb, leg_box) if leg_box is not None else 0.0
            inside = np.count_nonzero(
                (pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1)
                & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)) if len(pts) else 0
            score = 50.0 * (outside + on_legend) + inside + 0.01 * k
            if best is None or score < best[0]:
                best = (score, dx, dy, ha, va)
        _, dx, dy, ha, va = best
        label.xyann = (dx, dy)
        label.set_ha(ha)
        label.set_va(va)
    except Exception:                                     # pragma: no cover
        label.set_visible(True)


def _room_above(ax, y, points):
    """Raise the top of the y axis so that ``points`` (typographic points)
    fit between ``y`` and the top.  The limit-point label sits there, inside
    the axes, where the legend ("best") can see it; above the axes it ran
    into the title and the legend (a limit point at the top of the plot,
    e.g. the star dome)."""
    try:
        h_pt = ax.get_position().height * ax.figure.get_figheight() * 72.0
        y0, y1 = ax.get_ylim()
        f = 1.0 - points / max(h_pt, 1.0)
        if f < 0.4 or y <= y0:
            return
        top = y0 + (y - y0) / f
        if top > y1:
            ax.set_ylim(y0, top)
    except Exception:                                   # pragma: no cover
        pass


def _amp_x(rows):
    """x values for the imperfection chart: a / span, or a [mm]."""
    xs = [r.get("amplitude_over_span", float("nan")) for r in rows]
    if all(_finite(v) for v in xs):
        return np.asarray(xs, float), "Imperfection amplitude / span"
    return (np.asarray([r["amplitude_m"] for r in rows], float) * 1e3,
            "Imperfection amplitude [mm]")


def draw_imperfection_ratio(ax, rows, labels=None, mode=1):
    """lambda_limit / lambda_cr vs a / span (governing sign and each sign)."""
    ax.clear()
    x, xlabel = _amp_x(rows)
    lam_cr = next((r["lambda_cr"] for r in rows if _finite(r["lambda_cr"])),
                  float("nan"))
    y = np.array([r["ratio"] if _finite(r["ratio"]) else np.nan
                  for r in rows], float)

    def rel(v):
        return v / lam_cr if _finite(v) and _finite(lam_cr) else np.nan

    yp = np.array([rel(r.get("lambda_plus")) for r in rows], float)
    ym = np.array([rel(r.get("lambda_minus")) for r in rows], float)
    ax.axhline(1.0, color=COL_CR, ls="--", lw=1.2,
               label=r"Linear buckling ($\lambda_{cr}$)")
    if len(rows) > 1:
        ax.plot(x[1:], yp[1:], "^", ms=7, mfc="none", mec=COL_PLUS, mew=1.4,
                ls="none", label="+ mode shape")
        ax.plot(x[1:], ym[1:], "v", ms=7, mfc="none", mec=COL_MINUS,
                mew=1.4, ls="none", label="- mode shape")
    ax.plot(x, y, "-o", color=COL_PATH, lw=2.0, ms=6, zorder=5,
            label="Governing (lower) sign")
    labels = labels or [""] * len(rows)
    for xi, yi, lab in zip(x, y, labels):
        if lab and np.isfinite(yi):
            ax.annotate(lab, (xi, yi), xytext=(0, 9),
                        textcoords="offset points", ha="center",
                        fontsize=7.5, color="#334155")
    top = np.nanmax(np.concatenate([y, yp, ym, [1.0]]))
    ax.set_ylim(0.0, 1.12 * top if np.isfinite(top) else 1.12)
    xr = (np.nanmax(x) - np.nanmin(x)) or 1.0
    ax.set_xlim(np.nanmin(x) - 0.06 * xr, np.nanmax(x) + 0.08 * xr)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(r"$\lambda_{limit}\;/\;\lambda_{cr}$", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title(f"Imperfection sensitivity  --  buckling mode {mode}",
                 fontsize=10.5, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="best", framealpha=0.92)
    return ax


def draw_imperfection_paths(ax, rows, runs, labels=None):
    """Load paths of the governing sign for every amplitude."""
    ax.clear()
    labels = labels or [f"{r['amplitude_m'] * 1e3:.1f} mm" for r in rows]
    perfect = runs.get((rows[0]["amplitude_m"], 1.0)) if rows else None
    s = path_sign(perfect) if perfect is not None else -1.0
    try:
        from matplotlib import colormaps
        cmap = colormaps["viridis"]
    except Exception:                                   # pragma: no cover
        from matplotlib import cm
        cmap = cm.get_cmap("viridis")
    n = max(len(rows) - 1, 1)
    ytop = 0.0
    for i, r in enumerate(rows):
        a = r["amplitude_m"]
        sg = 1.0 if a == 0.0 else float(r.get("sign", 1.0))
        res = runs.get((a, sg))
        if res is None:
            continue
        col = "#111111" if a == 0.0 else \
            cmap(0.15 + 0.7 * (i - 1) / max(n - 1, 1))
        lab = labels[i] if i < len(labels) else ""
        if a != 0.0:
            lab += " (+)" if sg > 0 else " (-)"
        ax.plot(s * res.control * 1e3, res.lam, color=col, lw=1.6, label=lab)
        if res.limit_found and _finite(res.lambda_limit):
            ax.plot(s * res.control_limit * 1e3, res.lambda_limit, "*",
                    color=col, ms=11, mec="white", mew=0.6, ls="none")
        ytop = max(ytop, float(np.max(res.lam)))
    lam_cr = rows[0]["lambda_cr"] if rows else float("nan")
    if _finite(lam_cr) and 0.0 < lam_cr <= 4.0 * max(ytop, 1e-300):
        ax.axhline(lam_cr, color=COL_CR, ls="--", lw=1.2,
                   label=rf"$\lambda_{{cr}}$ = {lam_cr:.4g}")
    ax.axhline(0.0, color="0.65", lw=0.8)
    if perfect is not None:
        ax.set_xlabel(control_axis_label(perfect, s), fontsize=9)
    ax.set_ylabel(r"Load factor $\lambda$", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_title("Load paths (governing sign)", fontsize=10.5,
                 fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="best", framealpha=0.92)
    return ax


class _Snapshot:
    """Minimal stand-in for a StaticResult so plotting.draw_deflected can
    draw one path point."""

    def __init__(self, u, name):
        self.u = np.asarray(u, float)
        self.name = name

    def displacement_vectors(self):
        return self.u.reshape(-1, 3)


def _segments(nodes, members):
    m = np.asarray(members, int)
    return np.stack([nodes[m[:, 0]], nodes[m[:, 1]]], axis=1)


def _planar_box(ax, nodes):
    """For a truss in the x-z plane: drop the unused y depth from the 3D box
    so the elevation fills the axes.  x and z keep the same scale."""
    xl, zl = ax.get_xlim(), ax.get_zlim()
    xr, zr = xl[1] - xl[0], zl[1] - zl[0]
    yc = float(np.mean(nodes[:, 1]))
    ax.set_ylim(yc - 0.025 * xr, yc + 0.025 * xr)
    ax.set_yticks([])
    ax.set_ylabel("")
    try:
        from matplotlib.ticker import MaxNLocator
        ax.zaxis.set_major_locator(MaxNLocator(4))
        ax.set_box_aspect((xr, 0.05 * xr, zr), zoom=1.25)
        # the z axis title sat on the tick labels of the thin elevation
        ax.set_zlabel(ax.get_zlabel(), labelpad=14)
    except Exception:                                  # pragma: no cover
        pass


def _thin_z_ticks(ax, nodes):
    """Fewer z ticks on a flat structure (shallow dome, planar truss).  At
    true scale its z axis is short and the tick labels ran into each
    other."""
    try:
        from matplotlib.ticker import MaxNLocator
        xr = float(np.ptp(np.asarray(nodes, float)[:, 0])) or 1.0
        z0, z1 = ax.get_zlim()
        ratio = (z1 - z0) / xr
        if ratio < 0.15:
            ax.zaxis.set_major_locator(MaxNLocator(1))
        elif ratio < 0.35:
            ax.zaxis.set_major_locator(MaxNLocator(2))
    except Exception:                                  # pragma: no cover
        pass


# ======================================================== Qt widgets ======
def _exec(dlg):
    fn = getattr(dlg, "exec", None) or getattr(dlg, "exec_")
    return bool(fn())


def _combo_box(names, current=None):
    w = QtWidgets.QComboBox()
    w.addItems(list(names))
    if current in names:
        w.setCurrentIndex(list(names).index(current))
    w.setMinimumWidth(260)
    return w


def _method_box(current="arc-length"):
    w = QtWidgets.QComboBox()
    for key, text in METHODS:
        w.addItem(text, key)
    idx = w.findData(current)
    w.setCurrentIndex(max(idx, 0))
    w.setToolTip("Arc-length passes limit points and snap-through. "
                 "Displacement control passes load maxima. "
                 "Load control stops at the first limit point.")
    return w


def _form():
    f = QtWidgets.QFormLayout()
    f.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
    f.setHorizontalSpacing(10)
    f.setVerticalSpacing(7)
    return f


def _buttons(dlg, ok_text="Run"):
    bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok |
                                    QtWidgets.QDialogButtonBox.Cancel)
    ok = bb.button(QtWidgets.QDialogButtonBox.Ok)
    ok.setText(ok_text)
    ok.setObjectName("primary")
    bb.accepted.connect(dlg.accept)
    bb.rejected.connect(dlg.reject)
    return bb


class NonlinearDialog(QtWidgets.QDialog):
    """Settings for one nonlinear run.  ``values()`` returns the settings
    dict for :func:`nonlinear_job` or raises ValueError with a message."""

    def __init__(self, parent, model, last=None):
        super().__init__(parent)
        self.setWindowTitle("Nonlinear analysis")
        self.model = model
        last = dict(last or {})
        names = [c.name for c in model.combinations]
        try:
            auto = default_combination(model).name
        except Exception:
            auto = names[0] if names else None
        cur = last.get("combination") if last.get("combination") in names \
            else auto
        n = int(model.geom.n_nodes)
        crown = int(model.geom.crown_node)
        span = model_span(model)

        v = QtWidgets.QVBoxLayout(self)
        v.addWidget(st.hint(
            "Traces the load factor λ times the combination through limit "
            "points. Large displacements, elastic bars. Temperature and "
            "settlement (T) are left out."))

        box = QtWidgets.QGroupBox("Load and solver")
        f = _form()
        box.setLayout(f)
        self.combination = _combo_box(names, cur)
        self.combination.setToolTip("λ multiplies the D, S, L and W loads of "
                                    "this combination.")
        f.addRow("Combination", self.combination)
        self.method = _method_box(last.get("method", "arc-length"))
        f.addRow("Method", self.method)
        self.max_steps = st.spin(5, 5000, last.get("max_steps", 300), 10)
        self.max_steps.setToolTip("The run stops after this many steps.")
        f.addRow("Maximum steps", self.max_steps)
        sal = last.get("stop_after_limit", 15)
        self.after_limit = st.spin(0, 2000, 15 if sal is None else sal, 1)
        self.after_limit.setToolTip("Steps traced past the first limit point. "
                                    "Use more to follow the snap-through.")
        f.addRow("Steps past the limit", self.after_limit)

        self.target_auto = QtWidgets.QCheckBox("Auto (2 λcr)")
        self.target = st.dspin(1e-3, 1e6, last.get("lambda_max") or 10.0,
                               1.0, 3)
        self.target.setToolTip("The run stops when λ passes this value.")
        self.target_auto.setChecked(not last.get("lambda_max"))
        f.addRow("Stop at λ", self._pair(self.target_auto, self.target))

        self.inc_auto = QtWidgets.QCheckBox("Auto (λcr / 40)")
        self.increment = st.dspin(1e-6, 1e5, last.get("first_increment") or
                                  0.1, 0.05, 4)
        self.increment.setToolTip("Size of the first step as a load factor. "
                                  "The step then adapts to the iterations.")
        self.inc_auto.setChecked(not last.get("first_increment"))
        f.addRow("First increment Δλ", self._pair(self.inc_auto,
                                                  self.increment))

        node = last.get("control_node", crown)
        node = crown if node is None or not (0 <= int(node) < n) else int(node)
        self.node = st.spin(0, max(n - 1, 0), node, 1)
        self.node.setPrefix("node ")
        self.axis = QtWidgets.QComboBox()
        self.axis.addItems(["uz", "ux", "uy"])
        self.axis.setCurrentText(last.get("control_axis", "uz"))
        crown_btn = QtWidgets.QPushButton("Crown")
        crown_btn.setToolTip(f"Use the crown, node {crown}, direction uz.")
        crown_btn.clicked.connect(lambda: (self.node.setValue(crown),
                                           self.axis.setCurrentText("uz")))
        row = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(row)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(self.node)
        h.addWidget(self.axis)
        h.addWidget(crown_btn)
        h.addStretch(1)
        row.setToolTip("The path is plotted against this displacement. "
                       "Displacement control also steps it.")
        f.addRow("Control DOF", row)
        v.addWidget(box)

        self.imp_box = QtWidgets.QGroupBox("Initial imperfection "
                                           "(buckling mode shape)")
        self.imp_box.setCheckable(True)
        self.imp_box.setChecked(bool(last.get("imperfection_mm")))
        g = _form()
        self.imp_box.setLayout(g)
        self.imp_mode = st.spin(1, 6, last.get("imperfection_mode") or 1, 1)
        self.imp_mode.setToolTip("Linear buckling mode of the chosen "
                                 "combination.")
        g.addRow("Buckling mode", self.imp_mode)
        amp0 = last.get("imperfection_mm") or \
            (round(span, 3) if math.isfinite(span) else 10.0)
        self.imp_amp = st.dspin(-1e5, 1e5, amp0, 1.0, 1, " mm")
        self.imp_amp.setToolTip("Largest node offset. A negative value flips "
                                "the mode.")
        g.addRow("Amplitude", self.imp_amp)
        if math.isfinite(span):
            g.addRow("", st.hint(f"Span {span:.2f} m: span/1000 = "
                                 f"{span:.1f} mm, span/300 = "
                                 f"{span / 0.3:.1f} mm."))
        v.addWidget(self.imp_box)
        v.addWidget(_buttons(self))

    @staticmethod
    def _pair(chk, spin):
        w = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(chk)
        h.addWidget(spin, 1)
        spin.setEnabled(not chk.isChecked())
        chk.toggled.connect(lambda on: spin.setEnabled(not on))
        return w

    def values(self):
        if self.combination.count() == 0:
            raise ValueError("The model has no load combinations.")
        node = int(self.node.value())
        axis = self.axis.currentText()
        dof = 3 * node + AXES.index(axis)
        try:
            fixed, _ = self.model.support_conditions()
            fixed = set(int(d) for d in np.asarray(fixed).ravel())
        except Exception:
            fixed = set()
        if dof in fixed:
            raise ValueError(f"Node {node} {axis} is held by a support. "
                             f"Pick a free node or direction.")
        amp = float(self.imp_amp.value()) if self.imp_box.isChecked() else 0.0
        if self.imp_box.isChecked() and amp == 0.0:
            raise ValueError("Enter a non-zero amplitude, or switch the "
                             "imperfection off.")
        return {
            "combination": self.combination.currentText(),
            "method": self.method.currentData(),
            "max_steps": int(self.max_steps.value()),
            "stop_after_limit": int(self.after_limit.value()),
            "lambda_max": None if self.target_auto.isChecked()
            else float(self.target.value()),
            "first_increment": None if self.inc_auto.isChecked()
            else float(self.increment.value()),
            "control_node": node,
            "control_axis": axis,
            "imperfection_mode": int(self.imp_mode.value())
            if self.imp_box.isChecked() else None,
            "imperfection_mm": amp,
        }

    def accept(self):
        try:
            self.values()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Nonlinear analysis", str(exc))
            return
        super().accept()


class ImperfectionDialog(QtWidgets.QDialog):
    """Settings for the imperfection study."""

    def __init__(self, parent, model, last=None):
        super().__init__(parent)
        self.setWindowTitle("Imperfection study")
        self.model = model
        self.span = model_span(model)
        last = dict(last or {})
        names = [c.name for c in model.combinations]
        try:
            auto = default_combination(model).name
        except Exception:
            auto = names[0] if names else None
        cur = last.get("combination") if last.get("combination") in names \
            else auto

        v = QtWidgets.QVBoxLayout(self)
        v.addWidget(st.hint(
            "Moves the nodes by a buckling mode shape and runs the nonlinear "
            "analysis for each amplitude. Both signs are run and the lower "
            "limit load is kept."))
        box = QtWidgets.QGroupBox("Study")
        f = _form()
        box.setLayout(f)
        self.combination = _combo_box(names, cur)
        f.addRow("Combination", self.combination)
        self.mode = st.spin(1, 6, last.get("mode", 1), 1)
        self.mode.setToolTip("Linear buckling mode used as the imperfection "
                             "shape.")
        f.addRow("Buckling mode", self.mode)
        self.amplitudes = QtWidgets.QLineEdit(
            last.get("amplitudes_text", DEFAULT_AMPLITUDES_TEXT))
        self.amplitudes.setToolTip("Comma separated. span/n means span "
                                   "divided by n. A plain number is in mm.")
        f.addRow("Amplitudes", self.amplitudes)
        self.preview = st.hint("")
        f.addRow("", self.preview)
        self.method = _method_box(last.get("method", "arc-length"))
        f.addRow("Method", self.method)
        self.max_steps = st.spin(5, 5000, last.get("max_steps", 300), 10)
        f.addRow("Maximum steps", self.max_steps)
        self.after_limit = st.spin(0, 2000, last.get("stop_after_limit", 5), 1)
        self.after_limit.setToolTip("Steps traced past each limit point.")
        f.addRow("Steps past the limit", self.after_limit)
        v.addWidget(box)
        v.addWidget(_buttons(self))
        self.amplitudes.textChanged.connect(self._update_preview)
        self._update_preview()

    def _update_preview(self, *_):
        try:
            amps = parse_amplitudes(self.amplitudes.text(), self.span)
        except ValueError as exc:
            self.preview.setText(str(exc))
            return
        mm = ", ".join(f"{a * 1e3:.1f}" for a, _ in amps)
        head = f"Span {self.span:.2f} m: " if math.isfinite(self.span) else ""
        self.preview.setText(f"{head}{mm} mm. {1 + 2 * len(amps)} runs.")

    def values(self):
        if self.combination.count() == 0:
            raise ValueError("The model has no load combinations.")
        amps = parse_amplitudes(self.amplitudes.text(), self.span)
        return {"combination": self.combination.currentText(),
                "mode": int(self.mode.value()),
                "amplitudes": amps,
                "amplitudes_text": self.amplitudes.text(),
                "method": self.method.currentData(),
                "max_steps": int(self.max_steps.value()),
                "stop_after_limit": int(self.after_limit.value())}

    def accept(self):
        try:
            self.values()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Imperfection study", str(exc))
            return
        super().accept()


class _PathCanvas(MplCanvas):
    """Path (2D) on the left, deformed shape (3D) on the right."""

    def __init__(self, parent=None):
        super().__init__(parent, projection=None)
        self.reset_axes()
        self.canvas.mpl_connect("resize_event", self._on_resize)

    def reset_axes(self):
        fig = self.figure
        fig.clear()
        self._colorbar = None
        self.ax_path = fig.add_axes([0.07, 0.11, 0.38, 0.80])
        self.ax_3d = fig.add_axes([0.48, 0.02, 0.44, 0.86], projection="3d")
        self.cax = fig.add_axes([0.935, 0.24, 0.011, 0.50])
        self.cax.set_visible(False)
        self.axes = self.ax_path
        self._fit_path_axes()

    def _on_resize(self, _event=None):
        self._fit_path_axes()
        try:
            replace_limit_label(self.ax_path)
        except Exception:                               # pragma: no cover
            pass

    def _fit_path_axes(self):
        """Keep about 50 px under the path axes for the tick labels and the
        axis title, and 36 px above it for the title, at any canvas height
        (a short view pane used to cut the axis title off)."""
        h = float(self.figure.get_figheight() * self.figure.dpi)
        if h <= 0.0:
            return
        bottom = min(0.30, max(0.11, 50.0 / h))
        top = max(bottom + 0.25, min(0.91, 1.0 - 36.0 / h))
        self.ax_path.set_position([0.07, bottom, 0.38, top - bottom])

    def placeholder(self, text):
        self.reset_axes()
        self.ax_path.set_axis_off()
        self.ax_3d.set_visible(False)
        self.figure.text(0.5, 0.5, text, ha="center", va="center",
                         fontsize=11, color="#64748b")
        self.canvas.draw_idle()


class NonlinearView(QtWidgets.QWidget):
    """The "Nonlinear" view tab: load path page and imperfection page."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        # A scroll area keeps the view (and so the main window) small enough
        # for a 1366 x 768 screen; a larger window shows it without bars.
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setWidget(self.tabs)
        lay.addWidget(scroll)

        # ---- load path page ----------------------------------------------
        page = QtWidgets.QWidget()
        pv = QtWidgets.QVBoxLayout(page)
        pv.setContentsMargins(0, 0, 0, 0)
        pv.setSpacing(0)
        self.path_canvas = _PathCanvas()
        c = self.path_canvas
        self.colour = QtWidgets.QComboBox()
        self.colour.addItem("Displacement", "disp")
        self.colour.addItem("Axial force", "axial")
        self.colour.setToolTip("Colour of the deformed shape.")
        self.scale = st.dspin(1.0, 1e5, 1.0, 1.0, 1)
        self.scale.setPrefix("x ")
        self.scale.setToolTip("Magnification of the displacements. x 1 is "
                              "the true shape.")
        self.auto_btn = QtWidgets.QPushButton("Auto scale")
        self.true_btn = QtWidgets.QPushButton("True scale")
        self.undeformed = QtWidgets.QCheckBox("Show undeformed")
        self.undeformed.setChecked(True)
        for w in (QtWidgets.QLabel("Colour:"), self.colour,
                  QtWidgets.QLabel("Magnification:"), self.scale,
                  self.auto_btn, self.true_btn, self.undeformed):
            c.add_control(w)
        c.add_control_stretch()
        pv.addWidget(c, 1)

        scrub = QtWidgets.QWidget()
        sh = QtWidgets.QHBoxLayout(scrub)
        sh.setContentsMargins(8, 4, 8, 4)
        sh.setSpacing(8)
        sh.addWidget(QtWidgets.QLabel("Step"))
        self.slider = QtWidgets.QSlider(Qt.Horizontal)
        self.slider.setRange(0, 0)
        self.slider.setToolTip("Drag to move along the path.")
        sh.addWidget(self.slider, 1)
        self.step_label = QtWidgets.QLabel("")
        self.step_label.setObjectName("metric")
        self.step_label.setMinimumWidth(330)
        sh.addWidget(self.step_label)
        self.limit_btn = QtWidgets.QPushButton("Limit point")
        self.limit_btn.setToolTip("Jump to the step before the first load "
                                  "maximum.")
        self.play_btn = QtWidgets.QPushButton("Play")
        self.play_btn.setCheckable(True)
        self.play_btn.setToolTip("Run along the path.")
        sh.addWidget(self.limit_btn)
        sh.addWidget(self.play_btn)
        pv.addWidget(scrub)

        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.RichText)
        self.summary.setContentsMargins(10, 4, 10, 6)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        pv.addWidget(self.summary)
        self.tabs.addTab(page, "Load path")

        # ---- imperfection page -------------------------------------------
        page2 = QtWidgets.QWidget()
        iv = QtWidgets.QVBoxLayout(page2)
        iv.setContentsMargins(0, 0, 0, 0)
        iv.setSpacing(0)
        self.imp_canvas = MplCanvas(projection=None)
        self.imp_canvas.figure.clear()
        self.imp_ax_ratio = self.imp_canvas.figure.add_subplot(121)
        self.imp_ax_paths = self.imp_canvas.figure.add_subplot(122)
        self.imp_canvas.axes = self.imp_ax_ratio
        iv.addWidget(self.imp_canvas, 1)
        self.imp_summary = QtWidgets.QLabel("")
        self.imp_summary.setWordWrap(True)
        self.imp_summary.setTextFormat(Qt.RichText)
        self.imp_summary.setContentsMargins(10, 4, 10, 6)
        self.imp_summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        iv.addWidget(self.imp_summary)
        self.tabs.addTab(page2, "Imperfection study")

    def set_controls_enabled(self, on):
        for w in (self.slider, self.limit_btn, self.play_btn, self.colour,
                  self.scale, self.auto_btn, self.true_btn, self.undeformed):
            w.setEnabled(on)

    def imperfection_placeholder(self, text):
        fig = self.imp_canvas.figure
        fig.clear()
        self.imp_ax_ratio = fig.add_subplot(121)
        self.imp_ax_paths = fig.add_subplot(122)
        self.imp_canvas.axes = self.imp_ax_ratio
        for ax in (self.imp_ax_ratio, self.imp_ax_paths):
            ax.set_axis_off()
        fig.text(0.5, 0.5, text, ha="center", va="center", fontsize=11,
                 color="#64748b")
        self.imp_canvas.canvas.draw_idle()

    def current_canvas(self):
        return self.path_canvas if self.tabs.currentIndex() == 0 \
            else self.imp_canvas

    def save(self, path, dpi=200):
        """Save the figure of the page on screen (View > Save as PNG)."""
        self.current_canvas().save(path, dpi=dpi)


def _fill(table, columns, rows, numeric=()):
    """Fill a QTableWidget with pre-formatted strings."""
    table.clearContents()
    table.setColumnCount(len(columns))
    table.setHorizontalHeaderLabels(columns)
    table.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, text in enumerate(row):
            item = QtWidgets.QTableWidgetItem(str(text))
            if j in numeric:
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()


def _find_menu(win, title):
    for act in win.menuBar().actions():
        if act.text().replace("&", "").strip().lower() == title.lower():
            return act.menu()
    return None


def _main_toolbar(win):
    """The window's main toolbar.  The matplotlib navigation bars are
    QToolBars too, so the search must not just take the first one."""
    tb = getattr(win, "toolbar", None)
    if isinstance(tb, QtWidgets.QToolBar):
        return tb
    bars = win.findChildren(QtWidgets.QToolBar)
    anchor = getattr(win, "act_size", None)
    for b in bars:
        if anchor is not None and anchor in b.actions():
            return b
    for b in bars:
        if b.parent() is win:
            return b
    return None


# ========================================================== extension =====
class NonlinearExtension(QtCore.QObject):
    """State and slots of the nonlinear GUI (see the module docstring)."""

    job_finished = Signal(str)          # "done" | "failed" | "cancelled"

    PATH_COLUMNS = ["step", "lambda", "control [mm]", "iterations",
                    "residual", "det sign", "neg. eigenvalues",
                    "min eigenvalue", "stiffness param.",
                    "strain energy [J]"]
    IMP_COLUMNS = ["amplitude", "a [mm]", "a / span", "a / strut",
                   "lambda_limit", "lambda_limit / lambda_cr", "sign",
                   "lambda (+ mode)", "lambda (- mode)",
                   "first critical lambda", "lambda_cr",
                   "hub rise left [mm]"]

    def __init__(self, win):
        super().__init__(win)
        self.setObjectName(OBJECT_NAME)
        self.win = win
        self.nonlinear_result = None
        self.nonlinear_settings = None
        self.nonlinear_model = None          # model drawn (with imperfection)
        self.nonlinear_offset = None
        self.imperfection_result = None      # list of row dicts
        self.imperfection_runs = None        # {(a, sign): NonlinearResult}
        self.imperfection_settings = None
        self.imperfection_labels = None
        self.errors = []
        self.quiet = False                   # tests: no message boxes
        self.last_message = None
        self.last_status = ""
        self._running = None
        self._busy_set = False
        self._cancel = None
        self._last_nl = None
        self._last_imp = None
        self._marker = None
        self._sign = 1.0
        self._reset_view = True
        self._model_changed_during_run = False

        self._shape_timer = QtCore.QTimer(self)
        self._shape_timer.setSingleShot(True)
        self._shape_timer.setInterval(40)
        self._shape_timer.timeout.connect(self._draw_shape_safe)
        self._play_timer = QtCore.QTimer(self)
        self._play_timer.setInterval(140)
        self._play_timer.timeout.connect(self._play_tick)

    # ------------------------------------------------------------ install --
    def _build(self):
        win = self.win
        self.view = NonlinearView()
        win.view_tabs.addTab(self.view, "Nonlinear")
        v = self.view
        v.slider.valueChanged.connect(self._on_slider)
        v.limit_btn.clicked.connect(lambda *_: self.go_to_limit())
        v.play_btn.toggled.connect(self._on_play)
        v.colour.currentIndexChanged.connect(lambda *_: self._redraw_shape())
        v.scale.valueChanged.connect(lambda *_: self._redraw_shape())
        v.undeformed.toggled.connect(lambda *_: self._redraw_shape())
        v.auto_btn.clicked.connect(lambda *_: self._auto_scale())
        v.true_btn.clicked.connect(lambda *_: v.scale.setValue(1.0))

        # result tabs, before "Log" when it exists
        self.result_tab = QtWidgets.QWidget()
        rl = QtWidgets.QVBoxLayout(self.result_tab)
        rl.setContentsMargins(4, 4, 4, 4)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(st.hint("Corotational path, one row per converged "
                              "step. Control displacement in mm, raw sign."))
        top.addStretch(1)
        self.path_csv_btn = QtWidgets.QPushButton("Export path (CSV)...")
        self.path_csv_btn.clicked.connect(lambda *_: self._ask_csv("path"))
        top.addWidget(self.path_csv_btn)
        rl.addLayout(top)
        split = QtWidgets.QSplitter(Qt.Horizontal)
        self.summary_text = QtWidgets.QPlainTextEdit()
        self.summary_text.setReadOnly(True)
        self.path_table = st.table(self.PATH_COLUMNS)
        split.addWidget(self.summary_text)
        split.addWidget(self.path_table)
        split.setSizes([420, 620])
        rl.addWidget(split, 1)

        self.imp_tab = QtWidgets.QWidget()
        il = QtWidgets.QVBoxLayout(self.imp_tab)
        il.setContentsMargins(4, 4, 4, 4)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(st.hint("One row per amplitude. The first row is the "
                              "perfect structure."))
        top.addStretch(1)
        self.imp_csv_btn = QtWidgets.QPushButton("Export table (CSV)...")
        self.imp_csv_btn.clicked.connect(lambda *_: self._ask_csv("imp"))
        top.addWidget(self.imp_csv_btn)
        il.addLayout(top)
        self.imp_table = st.table(self.IMP_COLUMNS)
        il.addWidget(self.imp_table, 1)

        rt = win.result_tabs
        log_idx = -1
        for i in range(rt.count()):
            if rt.tabText(i).replace("&", "").strip() == "Log":
                log_idx = i
        if log_idx >= 0:
            rt.insertTab(log_idx, self.result_tab, "Nonlinear")
            rt.insertTab(log_idx + 1, self.imp_tab, "Imperfection")
        else:
            rt.addTab(self.result_tab, "Nonlinear")
            rt.addTab(self.imp_tab, "Imperfection")

        # actions
        self.act_nonlinear = win._act(
            "&Nonlinear analysis...", lambda *_: self.open_nonlinear_dialog(),
            tip="Large-displacement load path through the limit point")
        self.act_imperfection = win._act(
            "&Imperfection study...",
            lambda *_: self.open_imperfection_dialog(),
            tip="Limit load of domes shaped like a buckling mode")
        self.act_toolbar = win._act(
            "Nonlinear", lambda *_: self.open_nonlinear_dialog(),
            tip="Nonlinear analysis: load path and limit point")
        menu = getattr(win, "menu_analysis", None) or \
            _find_menu(win, "Analysis") or win.menuBar().addMenu("&Analysis")
        menu.addSeparator()
        menu.addAction(self.act_nonlinear)
        menu.addAction(self.act_imperfection)

        tb = _main_toolbar(win)
        if tb is not None:
            acts = tb.actions()
            anchor = getattr(win, "act_size", None)
            before = None
            if anchor in acts:
                i = acts.index(anchor)
                before = acts[i + 1] if i + 1 < len(acts) else None
            if before is not None:
                tb.insertAction(before, self.act_toolbar)
            else:
                tb.addSeparator()
                tb.addAction(self.act_toolbar)

        # A window whose View menu lists the tabs when it opens (menu_view)
        # shows this tab by itself; older windows get a fixed entry.
        vm = None if getattr(win, "menu_view", None) is not None else \
            _find_menu(win, "View")
        if vm is not None:
            act = win._act("Nonlinear", lambda *_: win.view_tabs
                           .setCurrentWidget(self.view))
            seps = [a for a in vm.actions() if a.isSeparator()]
            if seps:
                vm.insertAction(seps[0], act)
            else:
                vm.addAction(act)
            self.act_view = act

        self.stop_btn = QtWidgets.QPushButton("Stop")
        self.stop_btn.setToolTip("Stop the nonlinear run after the current "
                                 "step.")
        self.stop_btn.setStyleSheet("padding: 1px 12px;")
        self.stop_btn.setVisible(False)
        self.stop_btn.clicked.connect(lambda *_: self.stop())
        win.statusBar().addPermanentWidget(self.stop_btn)

        win.installEventFilter(self)
        self._clear_views()

    def eventFilter(self, obj, event):                 # noqa: N802 (Qt API)
        try:
            close = getattr(QtCore.QEvent, "Close", None) or \
                QtCore.QEvent.Type.Close
            if obj is self.win and event.type() == close and \
                    self._cancel is not None:
                self._cancel.set()
        except Exception:
            pass
        return False

    # ---------------------------------------------------------- hooks ----
    def on_model_changed(self, model):
        """A new model: the stored nonlinear results no longer apply."""
        if self._running:
            self._model_changed_during_run = True
            return
        self.nonlinear_result = None
        self.nonlinear_settings = None
        self.nonlinear_model = None
        self.nonlinear_offset = None
        self.imperfection_result = None
        self.imperfection_runs = None
        self.imperfection_settings = None
        self.imperfection_labels = None
        self._last_nl = None
        self._last_imp = None
        self._clear_views()

    def on_results(self, results):
        """Linear results do not change the nonlinear views."""
        return None

    # ---------------------------------------------------------- utils ----
    def _report_error(self, where, exc):
        tb = traceback.format_exc()
        self.errors.append(tb)
        try:
            self.win.log(f"ERROR in the nonlinear view ({where}): {exc}")
            self.win.log(tb)
            self._status(f"Error: {exc}")
        except Exception:
            pass
        self._message("warning", "Nonlinear view", f"{exc}")

    def _message(self, kind, title, text):
        self.last_message = (kind, title, text)
        if self.quiet:
            self.win.log(f"{title}: {text}")
            return
        getattr(QtWidgets.QMessageBox, kind)(self.win, title, text)

    def _status(self, text):
        self.last_status = text
        self.win.status(text)

    def is_running(self):
        return self._running is not None

    def _ready(self):
        if self._running or getattr(self.win, "_busy", False):
            self._status("Busy. Wait for the current job to finish.")
            return False
        return True

    def _model_for_run(self, read_inputs=True):
        """Deep copy of the window's model after reading the panels.  None
        when the window rejected an input (it has shown the message)."""
        if read_inputs:
            try:
                ok = self.win._read_inputs()
            except Exception as exc:
                raise ValueError(f"Could not read the inputs: {exc}") from None
            if ok is False:
                self._status("Check the inputs. Nothing was run.")
                return None
        return copy.deepcopy(self.win.model)

    def _begin(self, kind, message):
        self._running = kind
        self.win._set_busy(True, message)
        self._busy_set = True
        for a in (self.act_nonlinear, self.act_imperfection, self.act_toolbar):
            a.setEnabled(False)
        self.stop_btn.setVisible(True)
        self.stop_btn.setEnabled(True)

    def _end(self):
        self._running = None
        if self._busy_set:
            self._busy_set = False
            self.win._set_busy(False)
        for a in (self.act_nonlinear, self.act_imperfection, self.act_toolbar):
            a.setEnabled(True)
        self.stop_btn.setVisible(False)

    def stop(self):
        """Ask the running job to stop (at its next step or run)."""
        if self._cancel is not None and self._running:
            self._cancel.set()
            self.stop_btn.setEnabled(False)
            self._status("Stopping...")

    def wait(self, timeout_s=600.0):
        """Block with a local event loop until the current run ends.
        Returns the job state ("done", "failed", "cancelled", "idle" or
        "timeout")."""
        if not self._running:
            return "idle"
        loop = QtCore.QEventLoop()
        out = []

        def fin(state):
            out.append(state)
            loop.quit()
        self.job_finished.connect(fin)
        QtCore.QTimer.singleShot(int(timeout_s * 1000), loop.quit)
        try:
            fn = getattr(loop, "exec", None) or getattr(loop, "exec_")
            fn()
        finally:
            try:
                self.job_finished.disconnect(fin)
            except Exception:
                pass
        return out[0] if out else "timeout"

    def _exec_dialog(self, dlg):
        """Run a modal dialog (tests replace this)."""
        return _exec(dlg)

    # -------------------------------------------------------- dialogs ----
    def open_nonlinear_dialog(self):
        try:
            if not self._ready():
                return
            model = self._model_for_run()
            if model is None:
                return
            dlg = NonlinearDialog(self.win, model, self._last_nl)
            if not self._exec_dialog(dlg):
                return
            self.start_nonlinear(dlg.values())
        except ValueError as exc:
            self._message("warning", "Nonlinear analysis", str(exc))
        except Exception as exc:
            self._report_error("nonlinear dialog", exc)

    def open_imperfection_dialog(self):
        try:
            if not self._ready():
                return
            model = self._model_for_run()
            if model is None:
                return
            dlg = ImperfectionDialog(self.win, model, self._last_imp)
            if not self._exec_dialog(dlg):
                return
            self.start_imperfection(dlg.values())
        except ValueError as exc:
            self._message("warning", "Imperfection study", str(exc))
        except Exception as exc:
            self._report_error("imperfection dialog", exc)

    # ----------------------------------------------------------- runs ----
    def start_nonlinear(self, settings, synchronous=False, read_inputs=True):
        """Start a nonlinear run with ``settings`` (see :func:`nonlinear_job`).
        With ``synchronous=True`` it runs in the calling thread (scripts and
        tests).  Returns True when the run was started."""
        if not self._ready():
            return False
        model = self._model_for_run(read_inputs)
        if model is None:
            return False
        self._last_nl = dict(settings)
        self._cancel = threading.Event()
        self._model_changed_during_run = False
        if synchronous:
            self._running = "nonlinear"
            try:
                payload = nonlinear_job(model, settings, self._cancel)
            except Exception as exc:
                self._job_failed(f"{type(exc).__name__}: {exc}",
                                 traceback.format_exc())
                return True
            self._nonlinear_done(payload)
            return True
        self._begin("nonlinear", "Nonlinear analysis...")
        w = Worker(nonlinear_job, model, dict(settings), self._cancel)
        w.progressed.connect(self.win._on_progress)
        w.failed.connect(self._job_failed)
        w.succeeded.connect(self._nonlinear_done)
        self.win._worker = w
        w.start()
        return True

    def start_imperfection(self, settings, synchronous=False,
                           read_inputs=True):
        """Start an imperfection study (see :func:`imperfection_job`)."""
        if not self._ready():
            return False
        model = self._model_for_run(read_inputs)
        if model is None:
            return False
        self._last_imp = dict(settings)
        self._cancel = threading.Event()
        self._model_changed_during_run = False
        if synchronous:
            self._running = "imperfection"
            try:
                payload = imperfection_job(model, settings, self._cancel)
            except Exception as exc:
                self._job_failed(f"{type(exc).__name__}: {exc}",
                                 traceback.format_exc())
                return True
            self._imperfection_done(payload)
            return True
        self._begin("imperfection", "Imperfection study...")
        w = Worker(imperfection_job, model, dict(settings), self._cancel)
        w.progressed.connect(self.win._on_progress)
        w.failed.connect(self._job_failed)
        w.succeeded.connect(self._imperfection_done)
        self.win._worker = w
        w.start()
        return True

    def _job_failed(self, message, tb):
        kind = self._running or "run"
        cancelled = self._cancel is not None and self._cancel.is_set()
        try:
            self._end()
            if self._model_changed_during_run:
                self._model_changed_during_run = False
                self.on_model_changed(self.win.model)
            name = "Imperfection study" if kind == "imperfection" \
                else "Nonlinear analysis"
            if cancelled:
                self.win.log(f"{name} stopped by the user.")
                self._status(f"{name} stopped.")
                return
            text = plain_message(message)
            self.win.log(f"ERROR: {name} failed: {message}")
            self.win.log(tb)
            self._status(f"{name} failed.")
            self._message("critical", f"{name} failed", text)
        except Exception as exc:
            self._report_error("failure handler", exc)
        finally:
            self.job_finished.emit("cancelled" if cancelled else "failed")

    def _discard_if_stale(self):
        """True (and the views cleared) when the model changed mid-run."""
        if not self._model_changed_during_run:
            return False
        self._model_changed_during_run = False
        self.win.log("The model changed during the run, so its result was "
                     "discarded.")
        self._status("Result discarded: the model changed during the run.")
        self.on_model_changed(self.win.model)
        return True

    def _nonlinear_done(self, payload):
        try:
            self._end()
            if self._discard_if_stale():
                return
            res = payload["result"]
            self.nonlinear_result = res
            self.nonlinear_settings = payload["settings"]
            self.nonlinear_model = payload["model"]
            self.nonlinear_offset = payload["offset"]
            self._reset_view = True
            self._show_nonlinear()
            self.win.view_tabs.setCurrentWidget(self.view)
            self.view.tabs.setCurrentIndex(0)
            self.win.result_tabs.setCurrentWidget(self.result_tab)
            s = self.nonlinear_settings
            imp = ""
            if s.get("imperfection_mm"):
                imp = (f", imperfection mode {s.get('imperfection_mode')} "
                       f"{s['imperfection_mm']:g} mm")
            self.win.log(f"Nonlinear analysis ({res.method}, {res.combination}"
                         f"{imp}):")
            for line in res.summary_lines():
                self.win.log("  " + line)
            for msg in res.messages:
                self.win.log("  " + msg)
            if res.limit_found:
                txt = f"λ_limit = {res.lambda_limit:.4g}"
                if _finite(res.limit_ratio):
                    txt += f" = {res.limit_ratio:.3f} λ_cr"
            else:
                txt = f"no limit point ({res.stop_reason})"
            self._status(f"Nonlinear analysis done: {txt}, "
                         f"{res.n_points - 1} steps.")
        except Exception as exc:
            self._report_error("nonlinear result", exc)
        finally:
            self.job_finished.emit("done")

    def _imperfection_done(self, payload):
        try:
            self._end()
            if self._discard_if_stale():
                return
            rows = payload["rows"]
            settings = payload["settings"]
            self.imperfection_result = rows
            self.imperfection_runs = payload["runs"]
            self.imperfection_settings = dict(settings, span=payload["span"])
            labels = ["perfect"]
            amps = settings.get("amplitudes", [])
            if len(amps) == len(rows) - 1:
                labels += [lab for _, lab in amps]
            else:
                labels += [f"{r['amplitude_m'] * 1e3:.1f} mm" for r in rows[1:]]
            self.imperfection_labels = labels
            self._show_imperfection()
            self.win.view_tabs.setCurrentWidget(self.view)
            self.view.tabs.setCurrentIndex(1)
            self.win.result_tabs.setCurrentWidget(self.imp_tab)
            self.win.log(f"Imperfection study ({settings['combination']}, "
                         f"buckling mode {settings.get('mode', 1)}):")
            for lab, r in zip(labels, rows):
                self.win.log(f"  {lab:>12s}  "
                             f"a = {r['amplitude_m'] * 1e3:8.2f} mm"
                             f"  lambda_limit = {_g(r['lambda_limit'])}"
                             f"  ratio = {_g(r['ratio'], '.4f')}")
                if r.get("note"):
                    self.win.log(f"  {lab:>12s}  note: {r['note']}")
            worst = min((r for r in rows if _finite(r["ratio"])),
                        key=lambda r: r["ratio"], default=None)
            txt = (f"lowest λ_limit/λ_cr = {worst['ratio']:.3f}"
                   if worst else "no limit point found")
            n_amp = max(len(rows) - 1, 0)          # rows[0] is the perfect run
            self._status(f"Imperfection study done: {n_amp} "
                         f"amplitude{'s' if n_amp != 1 else ''} plus the "
                         f"perfect dome, {txt}.")
        except Exception as exc:
            self._report_error("imperfection result", exc)
        finally:
            self.job_finished.emit("done")

    # --------------------------------------------------------- drawing ---
    def _clear_views(self):
        self._play_timer.stop()
        v = self.view
        v.play_btn.blockSignals(True)
        v.play_btn.setChecked(False)
        v.play_btn.setText("Play")
        v.play_btn.blockSignals(False)
        v.slider.blockSignals(True)
        v.slider.setRange(0, 0)
        v.slider.setValue(0)
        v.slider.blockSignals(False)
        v.set_controls_enabled(False)
        v.step_label.setText("")
        v.summary.setText("")
        self._marker = None
        v.path_canvas.placeholder(
            "No nonlinear run yet.\n\nAnalysis > Nonlinear analysis... "
            "traces the load path\nthrough the limit point.")
        v.imperfection_placeholder(
            "No imperfection study yet.\n\nAnalysis > Imperfection study... "
            "runs the dome\nwith buckling-mode imperfections.")
        v.imp_summary.setText("")
        self.summary_text.clear()
        for t in (self.path_table, self.imp_table):
            t.clearContents()
            t.setRowCount(0)

    def _show_nonlinear(self):
        res = self.nonlinear_result
        v = self.view
        c = v.path_canvas
        c.reset_axes()
        k0 = res.limit_step if res.limit_found and \
            0 <= res.limit_step < res.n_points else res.n_points - 1
        self._marker, self._sign = draw_path(
            c.ax_path, res, step=k0,
            title=f"Equilibrium path  --  {res.combination.split('  ')[0]}"
                  f" ({res.method})")
        v.set_controls_enabled(res.n_points > 0)
        v.limit_btn.setEnabled(bool(res.limit_found))
        self._auto_scale(redraw=False)
        v.slider.blockSignals(True)
        v.slider.setRange(0, max(res.n_points - 1, 0))
        v.slider.setValue(k0)
        v.slider.blockSignals(False)
        self._update_step_label(k0)
        self._draw_shape()
        v.summary.setText(self._summary_html(res))
        self.summary_text.setPlainText(self._summary_plain(res))
        self._fill_path_table(res)

    def _summary_html(self, res):
        s = self.nonlinear_settings or {}
        parts = []
        if res.limit_found:
            kind = "bifurcation" if res.limit_type == "bifurcation" \
                else "limit point"
            parts.append(f"<b>λ<sub>limit</sub> = {res.lambda_limit:.4g}</b> "
                         f"({kind})")
        else:
            parts.append("<b>No limit point reached</b>")
        parts.append(f"λ<sub>cr</sub> = {_g(res.lambda_cr)}")
        parts.append(f"λ<sub>limit</sub> / λ<sub>cr</sub> = "
                     f"<b>{_g(res.limit_ratio, '.3f')}</b>")
        parts.append(f"{res.n_points - 1} steps, {res.method}, "
                     f"{res.combination.split('  ')[0]}")
        if s.get("imperfection_mm"):
            parts.append(f"imperfection: mode {s.get('imperfection_mode')}, "
                         f"{s['imperfection_mm']:g} mm")
        html = "&nbsp;&nbsp;|&nbsp;&nbsp;".join(parts)
        if res.stop_reason:
            html += (f"<br><span style='color:#64748b'>Stopped: "
                     f"{_html_escape(res.stop_reason)}.</span>")
        notes = [m for m in res.messages
                 if not m.startswith("First load maximum")
                 and not m.startswith("No limit point reached")
                 and "cut-back" not in m]
        for m in notes:
            html += (f"<br><span style='color:{st.WARN}'>Note: "
                     f"{_html_escape(m)}</span>")
        return html

    def _summary_plain(self, res):
        lines = list(res.summary_lines())
        s = self.nonlinear_settings or {}
        if s.get("imperfection_mm"):
            lines.append(f"Imperfection           mode "
                         f"{s.get('imperfection_mode')}, amplitude {s['imperfection_mm']:g} mm "
                         f"(lambda_cr of the perfect structure)")
        lines.append(f"Stop reason            {res.stop_reason}")
        if res.messages:
            lines.append("")
            lines.append("Messages")
            lines += ["  " + m for m in res.messages]
        return "\n".join(lines)

    def _fill_path_table(self, res):
        rows = []
        for k in range(res.n_points):
            rows.append([str(k), f"{res.lam[k]:.6g}",
                         f"{res.control[k] * 1e3:.4f}",
                         str(int(res.iterations[k])),
                         f"{res.residual[k]:.2e}",
                         f"{int(res.det_sign[k]):+d}",
                         str(int(res.n_negative[k])) if res.n_negative[k] >= 0
                         else "-",
                         _g(res.min_eig[k], ".4g"),
                         _g(res.stiffness_parameter[k], ".4f"),
                         _g(res.strain_energy[k], ".5g")])
        _fill(self.path_table, self.PATH_COLUMNS, rows,
              numeric=set(range(len(self.PATH_COLUMNS))))

    def _update_step_label(self, k):
        res = self.nonlinear_result
        if res is None or not res.n_points:
            self.view.step_label.setText("")
            return
        k = int(np.clip(k, 0, res.n_points - 1))
        self.view.step_label.setText(
            f"{k:3d} / {res.n_points - 1}   λ = {res.lam[k]:<9.4g} "
            f"{res.control_label.split(' (')[0]} = "
            f"{res.control[k] * 1e3:.2f} mm")

    def _on_slider(self, k):
        try:
            self._update_step_label(k)
            res = self.nonlinear_result
            if res is not None and self._marker is not None:
                x = self._sign * res.control[k] * 1e3
                self._marker.set_data([x], [res.lam[k]])
                self.view.path_canvas.canvas.draw_idle()
            self._shape_timer.start()
        except Exception as exc:
            self._report_error("slider", exc)

    def show_step(self, k):
        """Move the slider to path point ``k`` and redraw now."""
        self.view.slider.setValue(int(k))
        self._shape_timer.stop()
        self._draw_shape()

    def go_to_limit(self):
        res = self.nonlinear_result
        if res is not None and res.limit_found:
            self.show_step(res.limit_step)

    def _on_play(self, on):
        try:
            v = self.view
            if on and self.nonlinear_result is not None:
                if v.slider.value() >= v.slider.maximum():
                    v.slider.setValue(0)
                v.play_btn.setText("Pause")
                self._play_timer.start()
            else:
                self._play_timer.stop()
                v.play_btn.setText("Play")
        except Exception as exc:
            self._report_error("play", exc)

    def _play_tick(self):
        v = self.view
        if self.nonlinear_result is None or \
                v.slider.value() >= v.slider.maximum():
            v.play_btn.setChecked(False)
            return
        v.slider.setValue(v.slider.value() + 1)

    def _auto_scale(self, redraw=True):
        res = self.nonlinear_result
        model = self.nonlinear_model
        if res is None or model is None or not res.n_points:
            return
        d = np.linalg.norm(np.asarray(res.u).reshape(res.n_points, -1, 3),
                           axis=2)
        dmax = float(np.max(d)) if d.size else 0.0
        s = 1.0 if dmax < 1e-12 else \
            0.10 * float(model.geom.length_scale) / dmax
        s = float(round(min(max(s, 1.0), 1e5)))
        sb = self.view.scale
        sb.blockSignals(True)
        sb.setValue(s)
        sb.blockSignals(False)
        if redraw:
            self._redraw_shape()

    def _redraw_shape(self):
        if self.nonlinear_result is not None:
            self._shape_timer.start()

    def _draw_shape_safe(self):
        try:
            self._draw_shape()
        except Exception as exc:
            self._report_error("deformed shape", exc)

    def _draw_shape(self):
        res = self.nonlinear_result
        model = self.nonlinear_model
        if res is None or model is None or not res.n_points:
            return
        v = self.view
        c = v.path_canvas
        ax = c.ax_3d
        k = int(np.clip(v.slider.value(), 0, res.n_points - 1))
        keep = None if self._reset_view else (ax.elev, ax.azim)
        scale = float(v.scale.value())
        u = np.asarray(res.u[k], float)
        nodes = np.asarray(model.geom.nodes, float)
        title = f"Step {k}:  λ = {res.lam[k]:.4g}"
        if v.colour.currentData() == "axial":
            defo = nodes + scale * u.reshape(-1, 3)
            art = plot.draw_axial(ax, model, res.N[k], title=title,
                                  width_scale=True, deformed=defo)
            if v.undeformed.isChecked():
                from mpl_toolkits.mplot3d.art3d import Line3DCollection
                ax.add_collection3d(Line3DCollection(
                    _segments(nodes, model.geom.members), colors="#c8cfd6",
                    linewidths=0.7, zorder=0))
            ax.text2D(0.02, 0.80, f"magnification x{scale:,.1f}",
                      transform=ax.transAxes, fontsize=8, va="top")
            label = "axial force [kN]  (+ tension)"
        else:
            art = plot.draw_deflected(ax, model, _Snapshot(u, title),
                                      scale=scale,
                                      show_undeformed=v.undeformed.isChecked(),
                                      title=title)
            label = "displacement [mm]"
        ax.title.set_fontsize(10.5)
        planar = np.ptp(nodes[:, 1]) <= 1e-9 * max(float(np.ptp(nodes)),
                                                    1e-300)
        if planar:                       # x-z truss: show it as an elevation
            _planar_box(ax, nodes)
        _thin_z_ticks(ax, nodes)
        if keep is not None:
            ax.view_init(elev=keep[0], azim=keep[1])
        elif planar:
            ax.view_init(elev=0.0, azim=-90.0)
        self._reset_view = False
        c.cax.clear()
        if art is not None and hasattr(art, "get_array") and \
                art.get_array() is not None:
            c.cax.set_visible(True)
            cb = c.figure.colorbar(art, cax=c.cax)
            cb.set_label(label, fontsize=8)
            cb.ax.tick_params(labelsize=7)
        else:
            c.cax.set_visible(False)
        c.canvas.draw_idle()

    def _show_imperfection(self):
        rows = self.imperfection_result
        v = self.view
        fig = v.imp_canvas.figure
        fig.clear()
        v.imp_ax_ratio = fig.add_subplot(121)
        v.imp_ax_paths = fig.add_subplot(122)
        v.imp_canvas.axes = v.imp_ax_ratio
        mode = int(self.imperfection_settings.get("mode", 1))
        draw_imperfection_ratio(v.imp_ax_ratio, rows, self.imperfection_labels,
                                mode=mode)
        draw_imperfection_paths(v.imp_ax_paths, rows, self.imperfection_runs
                                or {}, self.imperfection_labels)
        v.imp_canvas.refresh(tight=True)
        v.imp_summary.setText(self._imperfection_html())
        self._fill_imp_table()

    def _imperfection_html(self):
        rows = self.imperfection_result or []
        s = self.imperfection_settings or {}
        labels = self.imperfection_labels or [""] * len(rows)
        parts = [f"{s.get('combination', '').split('  ')[0]}, buckling mode "
                 f"{s.get('mode', 1)}, {s.get('method', '')}"]
        if rows and _finite(rows[0]["ratio"]):
            parts.append(f"perfect: λ<sub>limit</sub> = "
                         f"{rows[0]['lambda_limit']:.4g} = "
                         f"<b>{rows[0]['ratio']:.3f}</b> λ<sub>cr</sub>")
        ok = [(lab, r) for lab, r in zip(labels, rows)
              if _finite(r["ratio"]) and r["amplitude_m"] > 0.0]
        if ok:
            lab, r = min(ok, key=lambda t: t[1]["ratio"])
            parts.append(f"lowest: {lab} ({r['amplitude_m'] * 1e3:.1f} mm): "
                         f"λ<sub>limit</sub> = {r['lambda_limit']:.4g} = "
                         f"<b>{r['ratio']:.3f}</b> λ<sub>cr</sub>")
        html = "&nbsp;&nbsp;|&nbsp;&nbsp;".join(parts)
        missing = [lab for lab, r in zip(labels, rows) if not r.get("found")]
        half = [lab for lab, r in zip(labels, rows)
                if r.get("found") and r["amplitude_m"] > 0.0 and
                not (_finite(r.get("lambda_plus")) and
                     _finite(r.get("lambda_minus")))]
        if half:
            html += (f"<br><span style='color:{st.WARN}'>Note: one sign "
                     f"found no limit point for {', '.join(half)} (shown as "
                     f"failed in the table), so the other sign is "
                     f"used.</span>")
        if missing:
            html += (f"<br><span style='color:{st.WARN}'>Note: no limit point "
                     f"for {', '.join(missing)}. Raise the step limits.</span>")
        if any(_finite(r.get("ratio_critical")) and _finite(r.get("ratio")) and
               r["ratio_critical"] < r["ratio"] - 1e-9 for r in rows):
            html += (f"<br><span style='color:{st.WARN}'>Note: some runs "
                     f"bifurcate before the load maximum. See the first "
                     f"critical λ in the Imperfection table.</span>")
        return html

    def _fill_imp_table(self):
        rows = self.imperfection_result or []
        labels = self.imperfection_labels or [""] * len(rows)
        out = []
        for lab, r in zip(labels, rows):
            a = r["amplitude_m"]
            sg = r.get("sign", 1.0)
            out.append([lab, f"{a * 1e3:.2f}",
                        _g(r.get("amplitude_over_span"), ".5f"),
                        _g(r.get("amplitude_over_strut"), ".4f"),
                        _g(r["lambda_limit"], ".5g"),
                        _g(r["ratio"], ".4f"),
                        "-" if a == 0.0 else ("+" if sg > 0 else "-") + " mode",
                        _sign_cell(r.get("lambda_plus"), a),
                        _sign_cell(r.get("lambda_minus"), a),
                        _g(r.get("lambda_critical"), ".5g"),
                        _g(r["lambda_cr"], ".5g"),
                        _g(1e3 * r["hub_rise_left_m"], ".1f")
                        if _finite(r.get("hub_rise_left_m")) else "-"])
        _fill(self.imp_table, self.IMP_COLUMNS, out,
              numeric=set(range(1, len(self.IMP_COLUMNS))))

    # ------------------------------------------------------------- CSV ----
    def _ask_csv(self, which):
        try:
            if which == "path" and self.nonlinear_result is None:
                self._message("information", "Export", "Run the nonlinear "
                              "analysis first.")
                return
            if which == "imp" and not self.imperfection_result:
                self._message("information", "Export", "Run the imperfection "
                              "study first.")
                return
            name = "nonlinear_path.csv" if which == "path" else \
                "imperfection_study.csv"
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self.win, "Export CSV", name, "CSV (*.csv)")
            if not path:
                return
            if which == "path":
                self.export_path_csv(path)
            else:
                self.export_imperfection_csv(path)
            self.win.log(f"Wrote {path}")
            self._status("CSV written.")
        except Exception as exc:
            self._report_error("CSV export", exc)

    def export_path_csv(self, path):
        """Write the corotational path (one row per step) to ``path``."""
        res = self.nonlinear_result
        if res is None:
            raise ValueError("No nonlinear result to export.")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["step", "lambda", "control_m", "iterations",
                        "residual", "det_sign", "n_negative", "min_eig",
                        "stiffness_parameter", "strain_energy_J"])
            for k in range(res.n_points):
                w.writerow([k, repr(float(res.lam[k])),
                            repr(float(res.control[k])),
                            int(res.iterations[k]), float(res.residual[k]),
                            int(res.det_sign[k]), int(res.n_negative[k]),
                            float(res.min_eig[k]),
                            float(res.stiffness_parameter[k]),
                            float(res.strain_energy[k])])
        return path

    def export_imperfection_csv(self, path):
        """Write the imperfection-study rows to ``path``."""
        rows = self.imperfection_result
        if not rows:
            raise ValueError("No imperfection study to export.")
        labels = self.imperfection_labels or [""] * len(rows)
        keys = list(rows[0].keys())
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["label"] + keys)
            for lab, r in zip(labels, rows):
                w.writerow([lab] + [r.get(k) for k in keys])
        return path


def _html_escape(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


# ============================================================ install =====
def install(win):
    """Add the nonlinear menu actions, toolbar button, view tab and result
    tabs to ``win`` and return the extension (the existing one if this
    window already has it)."""
    for obj in win.findChildren(QtCore.QObject):
        if obj.objectName() == OBJECT_NAME and \
                hasattr(obj, "open_nonlinear_dialog"):
            return obj
    ext = NonlinearExtension(win)
    ext._build()
    exts = getattr(win, "extensions", None)
    if isinstance(exts, list) and ext not in exts:
        exts.append(ext)
    return ext
