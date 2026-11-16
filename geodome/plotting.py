"""
Matplotlib 3D renderers.

Every function here draws onto an existing ``Axes3D`` and never touches
``pyplot``, so the same code serves the embedded Qt canvases, the headless
figure export used by the report writer, and Spyder's inline backend without
any backend conflicts.
"""

from __future__ import annotations

import re
import weakref

import matplotlib
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, ListedColormap
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import art3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ------------------------------------------------------ fast 3D lines ---
def _mpl_version():
    nums = re.findall(r"\d+", matplotlib.__version__)[:2]
    return tuple(int(n) for n in nums) if len(nums) == 2 else (0, 0)


# Before matplotlib 3.11, Line3DCollection projects its segments one by one
# in Python.  On a 10V dome that loop is most of a redraw (zoom, pan and
# rotate redraw the whole figure), so the collection below projects all
# segments in one numpy step there.  From 3.11 the stock code does the same.
_VECTOR_PROJECTION = _mpl_version() < (3, 11)


def _project_segments(coll):
    """The 2D segments of a Line3DCollection in one numpy step, or None
    when the stock code must do it (ragged, masked, clipped, non-finite)."""
    if getattr(coll, "_axlim_clip", False):
        return None
    segs = coll._segments3d
    if np.ma.isMA(segs) or (isinstance(segs, (list, tuple)) and
                            any(np.ma.isMA(s) for s in segs)):
        return None
    try:
        arr = np.asarray(segs, dtype=float)
    except (TypeError, ValueError):                 # ragged segments
        return None
    if arr.ndim != 3 or arr.shape[-1] != 3 or arr.size == 0 or \
            not np.all(np.isfinite(arr)):
        return None
    n, k, _ = arr.shape
    pts = arr.reshape(-1, 3)
    vec = np.vstack([pts.T, np.ones(len(pts))])
    vecw = np.dot(coll.axes.M, vec)
    w = vecw[3]
    xs, ys, zs = vecw[0] / w, vecw[1] / w, vecw[2] / w
    if not np.all(np.isfinite(zs)):
        return None
    LineCollection.set_segments(
        coll, np.stack([xs, ys], axis=-1).reshape(n, k, 2))
    return min(1e9, float(zs.min()))


class Line3DCollection(art3d.Line3DCollection):
    """matplotlib's Line3DCollection with a vectorised projection before
    matplotlib 3.11.  It draws exactly the same picture."""

    def do_3d_projection(self, *args, **kwargs):
        if _VECTOR_PROJECTION and not args and not kwargs:
            try:
                z = _project_segments(self)
            except Exception:                        # pragma: no cover
                z = None
            if z is not None:
                return z
        return super().do_3d_projection(*args, **kwargs)


def fast_lines(coll):
    """Give a stock Line3DCollection (e.g. made by ``ax.quiver``) the
    vectorised projection of :class:`Line3DCollection`."""
    if type(coll) is art3d.Line3DCollection:
        coll.__class__ = Line3DCollection
    return coll


FAMILY_COLORS = ["#1f6fb4", "#d1495b", "#2a9d8f", "#e9a000", "#7b53a6",
                 "#4f772d", "#c1666b", "#00798c", "#b56576", "#495057",
                 "#8ac926", "#ff924c"]

# Group ids of the collections that hold one segment per member.  The
# interactive viewer (gui/view3d.py) finds them by these ids to hide members.
MEMBERS_GID = "geodome-members"
GHOST_GID = "geodome-ghost"
LEGEND_GID = "geodome-legend"


# ------------------------------------------------------------------ frame ---
def _segments(nodes, members):
    return np.stack([nodes[members[:, 0]], nodes[members[:, 1]]], axis=1)


def family_label(group):
    """Letter of a strut family: A, B, ... Z, then S26, S27, ..."""
    g = int(group)
    return chr(ord("A") + g) if g < 26 else f"S{g}"


# ---------------------------------------------------------- legend column ---
_TITLE_PT = 26.0        # least room kept above the 3D axes for its title [pt]
_TITLE_PAD_PT = 12.0    # title pad plus a margin, added to its height [pt]
_RIGHT_PT = 10.0        # room kept right of the drawing region [pt]
_ZAXIS_PT = 42.0        # room for the z tick labels and the z label [pt]
_LEGEND_PAD_PT = 4.0    # gap between the figure edge and the legend [pt]
_LEGEND_GAP_PT = 8.0    # gap between the legend column and the drawing [pt]
_LEGEND_REF_PT = 8.0    # font size of the trial legends [pt]
_LEGEND_FONTS = (8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0)
_LEGEND_LAYOUTS = ((1, False), (1, True), (2, True), (3, True), (4, True))
_LEGEND_MAX_WIDTH = 0.42    # widest legend column, fraction of the figure
# Measured sizes of trial legends, shared by every layout: a view drawn
# again (a new axes, the same legend) does not measure them again.
_LEGEND_SIZES = {}
_LEGEND_SIZES_MAX = 400


def _figure(artist):
    fig = getattr(artist, "get_figure", None)
    return fig() if callable(fig) else artist.figure


class AxesLayout:
    """Axes locator of a 3D view: legend column on the left, drawing right.

    Matplotlib calls an axes locator at every draw with the renderer, so this
    layout follows the real figure size: on screen after each resize, and in
    a static export at its own size and dpi.  It

    * keeps the legend (when the view has one and it is shown) in a column at
      the left edge of the figure, outside the 3D drawing.  It picks the
      largest font, the long or the short labels and one to four columns so
      that every row fits the figure height.  When even the smallest layout
      is too tall (a custom truss with hundreds of families), it shows the
      rows that fit and ends with a "+ N more" row;
    * gives the 3D axes the rest of the width, keeps room above it for the
      title (measured, so a two-line title fits) and room on the right for
      the z axis labels;
    * stores the drawing region (figure coordinates) in :attr:`region`, a
      Bbox updated in place, so that a viewer can clip the model to it.

    The base rectangle is the axes' subplot cell (so a colour bar made with
    ``figure.colorbar(..., ax=ax)`` keeps its place), or its original
    position for an axes made with ``figure.add_axes``.
    """

    def __init__(self, ax=None):
        self.region = Bbox.from_extents(0.0, 0.0, 1.0, 1.0)
        self.legend_spec = None
        self.legend_visible = True
        self.choice = None          # (columns, short labels, font size)
        self.n_shown = 0            # legend rows shown (without "+ N more")
        self._ax = weakref.ref(ax) if ax is not None else None
        self._fit_key = None
        self._spec_key = None
        self._leg_bb = None
        self._title_state = None

    # ------------------------------------------------------------ legend --
    def set_legend(self, ax, entries, title):
        """Show ``entries`` [(colour, line width, label, short label), ...]
        as the legend of ``ax``, in the left column."""
        spec = {"entries": list(entries), "title": str(title)}
        self.legend_spec = spec
        # the trial sizes are kept per legend content
        self._spec_key = (spec["title"],
                          tuple(tuple(e) for e in spec["entries"]))
        self._fit_key = None
        self._leg_bb = None
        # a first legend at once, so that it exists before the first draw;
        # the draw then fits it to the figure
        self._make(ax, 1, False, 7.0)

    def clear_legend(self):
        self.legend_spec = None
        self._spec_key = None
        self._fit_key = None
        self._leg_bb = None
        self.choice = None
        self.n_shown = 0

    @property
    def has_legend(self):
        return self.legend_spec is not None

    def set_legend_visible(self, ax, visible):
        self.legend_visible = bool(visible)
        leg = ax.get_legend() if ax is not None else None
        if leg is not None and leg.get_gid() == LEGEND_GID:
            leg.set_visible(self.legend_visible)
        self._fit_key = None
        self._leg_bb = None

    def _handles(self, compact, count=None):
        """Legend handles: the first ``count`` entries (all when None) and,
        when some are left out, a last "+ N more" row."""
        entries = self.legend_spec["entries"]
        n = len(entries)
        count = n if count is None else max(1, min(int(count), n))
        out = [Line2D([0], [0], color=c, lw=lw, label=(short if compact
                                                       else full))
               for c, lw, full, short in entries[:count]]
        if count < n:
            out.append(Line2D([0], [0], color="none", lw=0.0,
                              label=f"+ {n - count} more"))
        return out

    def _kwargs(self, ncols, fontsize):
        return dict(loc="upper left", ncols=int(ncols), fontsize=fontsize,
                    title=self.legend_spec["title"],
                    title_fontsize=fontsize + 1.0, frameon=True,
                    framealpha=0.95, borderaxespad=0.0, handlelength=1.6,
                    handletextpad=0.6, columnspacing=1.0, labelspacing=0.45)

    def _anchor(self, fig):
        w = max(fig.bbox.width, 1.0)
        h = max(fig.bbox.height, 1.0)
        pad = _LEGEND_PAD_PT * fig.dpi / 72.0
        return (pad / w, 1.0 - pad / h)

    def _make(self, ax, ncols, compact, fontsize, count=None):
        fig = _figure(ax)
        leg = ax.legend(handles=self._handles(compact, count),
                        bbox_to_anchor=self._anchor(fig),
                        bbox_transform=fig.transFigure,
                        **self._kwargs(ncols, fontsize))
        leg.set_gid(LEGEND_GID)
        leg.set_visible(self.legend_visible)
        self.choice = (int(ncols), bool(compact), float(fontsize))
        n = len(self.legend_spec["entries"])
        self.n_shown = n if count is None else max(1, min(int(count), n))
        self._leg_bb = None
        return leg

    def _size(self, ax, renderer, ncols, compact, count=None):
        """(width, height) [pt] of the legend at the trial font size."""
        key = (self._spec_key, ncols, compact, count,
               round(float(renderer.points_to_pixels(72.0)), 3))
        size = _LEGEND_SIZES.get(key)
        if size is None:
            handles = self._handles(compact, count)
            leg = Legend(ax, handles, [h.get_label() for h in handles],
                         **self._kwargs(ncols, _LEGEND_REF_PT))
            bb = leg.get_window_extent(renderer)
            k = 72.0 / max(renderer.points_to_pixels(72.0), 1e-9)
            size = (bb.width * k, bb.height * k)
            if len(_LEGEND_SIZES) >= _LEGEND_SIZES_MAX:
                _LEGEND_SIZES.clear()
            _LEGEND_SIZES[key] = size
        return size

    def _truncated(self, ax, renderer, avail_h, max_w, n):
        """The layout of a legend that is too tall for the figure even at
        the smallest font: (columns, rows shown).  It takes the most columns
        that fit ``max_w`` and cuts the rows to fit ``avail_h``; a last
        "+ N more" row says how many are left out."""
        fs = _LEGEND_FONTS[-1]
        k = 1.03 * fs / _LEGEND_REF_PT
        cols = sorted((c for c, compact in _LEGEND_LAYOUTS if compact
                       and (c == 1 or (c - 1) * 2 < n + 1)), reverse=True)
        ncols = 1
        for c in cols:
            if self._size(ax, renderer, c, True)[0] * k <= max_w:
                ncols = c
                break
        # the height grows with the rows: estimate from two sizes, then check
        rows_all = -(-n // ncols)
        _, h_all = self._size(ax, renderer, ncols, True)
        _, h_two = self._size(ax, renderer, 1, True, 1)     # one + "more"
        per_row = max((h_all - h_two) / max(rows_all - 2, 1), 1e-6)
        rows = int((avail_h / k - h_two) / per_row) + 2
        count = max(1, min(n - 1, rows * ncols - 1))
        for _ in range(40):
            _, h = self._size(ax, renderer, ncols, True, count)
            if h * k <= avail_h or count <= 1:
                break
            count = max(1, count - ncols)
        return ncols, count

    def _fit(self, ax, renderer, fig):
        """Choose the legend layout for this figure size; returns the
        legend (made again only when the size changed)."""
        leg = ax.get_legend()
        W, H = float(fig.bbox.width), float(fig.bbox.height)
        key = (round(W), round(H), round(float(fig.dpi), 3))
        if key == self._fit_key and leg is not None \
                and leg.get_gid() == LEGEND_GID:
            return leg
        px_per_pt = max(renderer.points_to_pixels(1.0), 1e-9)
        avail_h = (H - 2.0 * _LEGEND_PAD_PT * px_per_pt) / px_per_pt
        max_w = _LEGEND_MAX_WIDTH * W / px_per_pt
        n = len(self.legend_spec["entries"])
        choice = None
        for min_font in (6.0, _LEGEND_FONTS[-1]):
            for ncols, compact in _LEGEND_LAYOUTS:
                if ncols > 1 and (ncols - 1) * 2 >= n + 1:
                    continue            # a column would be (nearly) empty
                w8, h8 = self._size(ax, renderer, ncols, compact)
                for fs in _LEGEND_FONTS:
                    if fs < min_font:
                        break
                    k = 1.03 * fs / _LEGEND_REF_PT     # 3 % safety
                    if h8 * k <= avail_h and w8 * k <= max_w:
                        choice = (ncols, compact, fs)
                        break
                if choice:
                    break
            if choice:
                break
        if choice is not None:
            leg = self._make(ax, *choice)
        else:                       # too many rows: show the rows that fit
            ncols, count = self._truncated(ax, renderer, avail_h, max_w, n)
            leg = self._make(ax, ncols, True, _LEGEND_FONTS[-1], count)
        self._fit_key = key
        return leg

    def _legend_bbox(self, leg, renderer):
        """The legend's window extent, measured once per layout."""
        key = (self._fit_key, id(leg),
               round(float(renderer.points_to_pixels(72.0)), 3))
        c = self._leg_bb
        if c is None or c[0] != key:
            c = (key, leg.get_window_extent(renderer).frozen())
            self._leg_bb = c
        return c[1]

    # ----------------------------------------------------------- locator --
    def get_subplotspec(self):
        """The axes' subplot cell (lets ``tight_layout`` handle the axes)."""
        ax = self._ax() if self._ax is not None else None
        # before matplotlib 3.7 only subplots have get_subplotspec
        get = getattr(ax, "get_subplotspec", None)
        return get() if callable(get) else None

    def _base(self, ax, fig):
        get = getattr(ax, "get_subplotspec", None)
        ss = get() if callable(get) else None
        if ss is not None:
            return ss.get_position(fig)
        return ax.get_position(original=True)

    def _title_room(self, ax, renderer, dpi):
        """Room above the axes for its title [px]: its measured height (any
        number of lines, at the size asked for) plus the pad, and at least
        ``_TITLE_PT``."""
        room = _TITLE_PT * dpi / 72.0
        try:
            t = ax.title
            st = self._title_state
            now = float(t.get_fontsize())
            asked = st[1] if (st is not None and st[0] == t.get_text()
                              and abs(st[2] - now) < 1e-6) else now
            t.set_fontsize(asked)
            h = t.get_window_extent(renderer).height
            t.set_fontsize(now)
            room = max(room, h + _TITLE_PAD_PT * dpi / 72.0)
        except Exception:                            # pragma: no cover
            pass
        return room

    def _fit_title(self, ax, renderer, rect, W, H):
        """Shrink the title (down to 7 pt) when it would run into the legend
        column or off the figure."""
        t = ax.title
        text = t.get_text()
        if not text:
            return
        now = float(t.get_fontsize())
        st = self._title_state
        if st is None or st[0] != text or abs(st[2] - now) > 1e-6:
            req = now                       # the size the caller asked for
        else:
            req = st[1]
        t.set_fontsize(req)
        width = t.get_window_extent(renderer).width
        side = min(rect.width * W, rect.height * H)       # the square axes
        cx = (rect.x0 + 0.5 * rect.width) * W
        avail = 2.0 * min(cx - rect.x0 * W, W - cx) - 6.0
        fs = req
        if width > avail > 0.0 and side > 0.0:
            fs = max(7.0, req * avail / width)
            t.set_fontsize(fs)
        self._title_state = (text, req, fs)

    def __call__(self, ax, renderer):
        fig = _figure(ax)
        W = max(float(fig.bbox.width), 1.0)
        H = max(float(fig.bbox.height), 1.0)
        dpi = float(fig.dpi)
        x0, y0, x1, y1 = self._base(ax, fig).extents
        if ax.get_title():
            y1 = min(y1, 1.0 - self._title_room(ax, renderer, dpi) / H)
        x1 = min(x1, 1.0 - _RIGHT_PT * dpi / 72.0 / W)
        if self.legend_spec is not None and self.legend_visible:
            leg = self._fit(ax, renderer, fig)
            bb = self._legend_bbox(leg, renderer)
            gap = _LEGEND_GAP_PT * dpi / 72.0
            x0 = max(x0, (bb.x1 + gap) / W)
        x1 = max(x1, x0 + 0.05)
        y1 = max(y1, y0 + 0.05)
        self.region.set(Bbox.from_extents(x0, y0, x1, y1))
        # The 3D axes is a square centred in its rectangle, and its z tick
        # labels and z label stand out on its right.  When the square would
        # come closer than _ZAXIS_PT to the right edge, the rectangle is
        # narrowed.  The region stays: a zoomed model still uses the width.
        w, h = (x1 - x0) * W, (y1 - y0) * H
        z = _ZAXIS_PT * dpi / 72.0
        w_ax = min(w, max(w - z, 2.0 * w - 2.0 * z - h))
        w_ax = max(w_ax, min(w, 0.05 * W))
        rect = Bbox.from_extents(x0, y0, x0 + w_ax / W, y1)
        try:
            self._fit_title(ax, renderer, rect, W, H)
        except Exception:                            # pragma: no cover
            pass
        return rect


def axes_layout(ax, create=True):
    """The :class:`AxesLayout` of a 3D axes (made and set when missing).

    Returns None when the axes already has another locator.
    """
    lay = ax.get_axes_locator()
    if isinstance(lay, AxesLayout):
        return lay
    if lay is not None or not create:
        return None
    lay = AxesLayout(ax)
    ax.set_axes_locator(lay)
    return lay


def legend_column(ax, entries, title):
    """Draw ``entries`` as a legend in a column left of the 3D drawing.

    ``entries`` is a list of (colour, line width, label, short label).  The
    short labels are used when the long ones do not fit.
    """
    lay = axes_layout(ax)
    if lay is None:                                  # pragma: no cover
        handles = [Line2D([0], [0], color=c, lw=lw, label=full)
                   for c, lw, full, _ in entries]
        return ax.legend(handles=handles, loc="upper left", fontsize=7,
                         title=title)
    lay.set_legend(ax, entries, title)
    return ax.get_legend()


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
    lay = ax.get_axes_locator()
    if isinstance(lay, AxesLayout):
        lay.clear_legend()                  # each view sets its own legend
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
    lc.set_gid(MEMBERS_GID)
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
        # in a column left of the drawing (AxesLayout), never over the model
        entries = [(FAMILY_COLORS[s["id"] % len(FAMILY_COLORS)], 2.4,
                    f"{s['label']}  {s['length']:.3f} m  "
                    + (f"(CF {s['chord_factor']:.4f})  "
                       if np.isfinite(s['chord_factor']) else "")
                    + f"x{s['count']}",
                    f"{s['label']} {s['length']:.3f} m x{s['count']}")
                   for s in geom.strut_schedule()]
        legend_column(ax, entries, "Strut families")
    return ax


# -------------------------------------------------------------- load view ---
def draw_loads(ax, model, F, title="Applied nodal loads", arrow_frac=0.09):
    geom = model.geom
    prepare_axes(ax, title)
    segs = _segments(geom.nodes, geom.members)
    ax.add_collection3d(Line3DCollection(segs, colors="#b9c6d2", linewidths=0.8,
                                         gid=MEMBERS_GID))

    V = np.asarray(F).reshape(-1, 3)
    mag = np.linalg.norm(V, axis=1)
    mmax = mag.max()
    if mmax > 1e-9:
        scale = arrow_frac * geom.length_scale / mmax
        fast_lines(ax.quiver(
            geom.nodes[:, 0], geom.nodes[:, 1], geom.nodes[:, 2],
            V[:, 0] * scale, V[:, 1] * scale, V[:, 2] * scale,
            color="#c1121f", linewidth=0.9, arrow_length_ratio=0.28))
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


def deflected_nodes(geom, result, scale=None):
    """Drawn node coordinates of the deflected view and the scale used:
    ``(nodes + scale * u, scale)`` (``scale=None``: :func:`auto_scale`)."""
    if scale is None:
        scale = auto_scale(geom, result.u)
    disp = result.displacement_vectors()
    return np.asarray(geom.nodes, float) + float(scale) * disp, float(scale)


def draw_deflected(ax, model, result, scale=None, show_undeformed=True,
                   title=None):
    geom = model.geom
    disp = result.displacement_vectors()
    defo, scale = deflected_nodes(geom, result, scale)
    mag = np.linalg.norm(disp, axis=1)

    prepare_axes(ax, title or f"Deflected shape  --  {result.name}")
    if show_undeformed:
        ax.add_collection3d(Line3DCollection(
            _segments(geom.nodes, geom.members),
            colors="#c8cfd6", linewidths=0.7, gid=GHOST_GID))

    segs = _segments(defo, geom.members)
    vals = 0.5 * (mag[geom.members[:, 0]] + mag[geom.members[:, 1]]) * 1e3
    lc = Line3DCollection(segs, linewidths=1.8, cmap="viridis",
                          norm=Normalize(0.0, max(vals.max(), 1e-9)))
    lc.set_gid(MEMBERS_GID)
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
    lc.set_gid(MEMBERS_GID)
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
                                         colors=cols, linewidths=lw,
                                         gid=MEMBERS_GID))
    _draw_supports(ax, geom)
    _draw_ground(ax, geom)
    set_equal_aspect(ax, geom.nodes)

    labels = ["< 0.25", "0.25-0.50", "0.50-0.75", "0.75-0.90",
              "0.90-1.00", "> 1.00  OVER"]
    legend_column(ax, [(c, 3.0, lab, lab)
                       for c, lab in zip(_UTIL_COLORS, labels)],
                  "Utilisation")
    k = int(np.argmax(util))
    ax.text2D(0.98, 0.96,
              f"max utilisation {util[k]:.3f} at member {k}\n"
              f"{int(np.sum(util > 1.0))} strut(s) over 1.0",
              transform=ax.transAxes, fontsize=8, va="top", ha="right",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                        alpha=0.85))
    return ax


# ---------------------------------------------------------- buckling view ---
def buckling_nodes(geom, bres, k=0, scale_frac=0.10):
    """Drawn node coordinates of buckling mode ``k`` (0-based):
    ``(nodes + scale * phi, scale, phi)``; the largest drawn displacement is
    ``scale_frac * geom.length_scale``."""
    k = int(np.clip(k, 0, len(bres.load_factors) - 1))
    phi = np.asarray(bres.mode(k), float)
    pmax = np.max(np.linalg.norm(phi, axis=1))
    scale = scale_frac * geom.length_scale / max(pmax, 1e-12)
    return np.asarray(geom.nodes, float) + scale * phi, float(scale), phi


def draw_buckling_mode(ax, model, bres, k=0, scale_frac=0.10, title=None):
    geom = model.geom
    if bres is None or not len(bres.load_factors):
        prepare_axes(ax, "Buckling mode -- not available")
        ax.text2D(0.5, 0.5, "No positive buckling mode for this load case",
                  transform=ax.transAxes, ha="center", fontsize=10)
        return ax
    k = int(np.clip(k, 0, len(bres.load_factors) - 1))
    defo, scale, phi = buckling_nodes(geom, bres, k, scale_frac)
    mag = np.linalg.norm(phi, axis=1)

    prepare_axes(ax, title or f"Buckling mode {k + 1}   "
                              f"lambda_cr = {bres.load_factors[k]:.2f}")
    ax.add_collection3d(Line3DCollection(_segments(geom.nodes, geom.members),
                                         colors="#ccd3da", linewidths=0.7,
                                         gid=GHOST_GID))
    segs = _segments(defo, geom.members)
    vals = 0.5 * (mag[geom.members[:, 0]] + mag[geom.members[:, 1]])
    lc = Line3DCollection(segs, cmap="magma_r", linewidths=1.9,
                          norm=Normalize(0.0, max(vals.max(), 1e-12)))
    lc.set_gid(MEMBERS_GID)
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
