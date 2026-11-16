"""Interactive 3D viewer for the GeoDomeLab views.

``install(canvas, ...)`` puts three things on a :class:`~.canvas.MplCanvas`
that holds a 3D axes:

* :class:`ViewController` listens to the matplotlib mouse and key events of
  the canvas.  Left drag rotates, the wheel zooms at the cursor, right or
  middle drag pans, a click selects the member or node under the cursor and
  the keyboard moves the camera.  The full list is :data:`MOUSE_CONTROLS`
  and :data:`KEY_CONTROLS`.
* :class:`ViewBar`, a row of buttons over the drawing: preset views, zoom,
  fit, reset, projection, axes, legend, member and node numbers, a Find box
  and a help button.
* :class:`InfoPanel`, a panel on the right of the drawing with the details
  of the selected member or node.  Its x button closes it.

The controller works on the coordinates actually drawn (the deformed ones
in the deflected and buckling views).  Picking is vectorised: the nodes are
projected with the axes' projection matrix (``ax.get_proj()``) and the
member nearest to the cursor is found from the 2D point-segment distances.

The camera (angles, limits, projection) is kept when a view is drawn again
for the same model, and goes back to the default when the model changes.
The selection is a :class:`Selection` object that several views can share
(the main window shares one between its six 3D views and the result
tables).

Nothing here keeps a window alive: the matplotlib callbacks are bound
methods (matplotlib holds them weakly), the Qt objects are children of the
canvas, and every handler catches its own errors (they are kept in
``ViewController.errors`` and written to stderr).
"""

from __future__ import annotations

import functools
import html
import math
import re
import sys
import traceback
import weakref
import zlib

import numpy as np
from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Qt

from matplotlib.transforms import Bbox, TransformedBbox
from mpl_toolkits.mplot3d import art3d

from .. import plotting as plot
from . import style as st

# ------------------------------------------------------------- settings --
DRAG_PX = 4.0            # a press that moves less than this is a click
MEMBER_PICK_PX = 7.0     # largest distance from the cursor to a member
NODE_PICK_PX = 8.0       # largest distance from the cursor to a node
NODE_SNAP_PX = 5.0       # a node this near the cursor wins over its members
DEPTH_BIAS_PX = 3.0      # nearer items win ties by up to this distance
WHEEL_FACTOR = 1.2       # zoom per wheel notch
WHEEL_MAX_NOTCHES = 4.0  # most notches one wheel event may zoom
KEY_ZOOM = 1.25          # zoom per key press
ROTATE_STEP = 5.0        # degrees per arrow key
ROTATE_STEP_BIG = 15.0   # degrees per arrow key with Shift
PAN_STEP = 0.08          # pan per Ctrl + arrow key, fraction of the view
HOVER_MAX_MEMBERS = 50000    # hover highlight is off above this size
LABEL_CAP = 200          # most member or node numbers shown at once
LABEL_GAP_PX = 18.0      # smallest distance between two numbers on screen
SETTLE_MS = 200          # quiet time before the numbers are placed again

SELECT_COLOUR = "#ff1f8e"
SELECT_EDGE = "#1b1b1b"
HOVER_COLOUR = "#00a6d6"
MEMBER_LABEL_COLOUR = "#123a5c"
NODE_LABEL_COLOUR = "#5b3a00"

# (name, elevation, azimuth) of the preset views, keys 1 to 6
PRESETS = (("Iso", 35.264, -45.0), ("Top", 90.0, -90.0),
           ("Front", 0.0, -90.0), ("Back", 0.0, 90.0),
           ("Left", 0.0, 180.0), ("Right", 0.0, 0.0))

MOUSE_CONTROLS = (
    ("Left drag", "Rotate (orbit) the view."),
    ("Mouse wheel", "Zoom in or out at the cursor."),
    ("Right drag, middle drag, Shift + left drag", "Pan."),
    ("Click", "Select the member under the cursor, or the node when you "
              "click on or next to a node."),
    ("Ctrl + click", "Add to or remove from the selection."),
    ("Click on empty space", "Clear the selection."),
    ("Double-click", "Zoom to the member or node."),
    ("Hover", "Highlight the member or node and show a short tip."),
)

KEY_CONTROLS = (
    ("Arrow keys", "Orbit the camera by 5 degrees (Shift: 15)."),
    ("Ctrl + arrow keys", "Pan."),
    ("Q and E", "Roll the view by 5 degrees (Shift: 15)."),
    ("+ or =, Page Up", "Zoom in."),
    ("-, Page Down", "Zoom out."),
    ("Home or R", "Reset the view."),
    ("F", "Fit the selection, or everything when nothing is selected."),
    ("1 to 6", "Iso, top, front, back, left and right views."),
    ("Esc", "Clear the selection."),
    ("L", "Show or hide the legend."),
    ("M", "Show or hide the member numbers."),
    ("N", "Show or hide the node numbers."),
    ("A", "Show or hide the axes and the grid."),
    ("O or P", "Switch between orthographic and perspective projection."),
    ("H", "Hide the selected members."),
    ("Shift + H", "Show every member again."),
    ("I", "Isolate the selection (hide the other members)."),
    ("Ctrl + C", "Copy the details of the selection."),
)


def help_html():
    """The mouse and keyboard controls as an HTML page (help dialog)."""
    def table(rows):
        return ("<table cellspacing='0' cellpadding='3'>" + "".join(
            f"<tr><td><b>{html.escape(a)}</b></td><td>{html.escape(b)}</td>"
            f"</tr>" for a, b in rows) + "</table>")
    return ("<h3>3D view controls</h3>"
            "<p>Click the drawing first, so that it takes the keys.</p>"
            "<h4>Mouse</h4>" + table(MOUSE_CONTROLS) +
            "<h4>Keyboard</h4>" + table(KEY_CONTROLS) +
            "<h4>View bar</h4>"
            "<p>Iso, Top, Front and Side set the view. + and - zoom, Fit "
            "frames the selection, Reset goes back to the first view. Ortho "
            "switches the projection. Axes, Legend, Member # and Node # show "
            "or hide those parts. Info opens the selection details. Find "
            "selects a member or a node by its number (for example m 12 or "
            "n 5) and frames it.</p>"
            "<p>The matplotlib toolbar under the drawing still saves the "
            "figure. Its own pan and zoom modes work as before. While one "
            "of them is on, the viewer leaves the mouse to it.</p>")


# ======================================================== projection ======
def _figure(artist):
    fn = getattr(artist, "get_figure", None)
    return fn() if callable(fn) else getattr(artist, "figure", None)


def project(ax, pts, M=None):
    """Display coordinates (x, y in pixels) and depth of 3D points.

    Uses the same projection as matplotlib: ``M = ax.get_proj()`` and the
    axes' ``transData``.  A larger depth is further from the viewer.
    """
    P = np.atleast_2d(np.asarray(pts, dtype=float))
    if M is None:
        M = np.asarray(ax.get_proj(), dtype=float)
    h = P @ M[:, :3].T + M[:, 3]
    w = h[:, 3].copy()
    w[np.abs(w) < 1e-300] = 1e-300
    xy = ax.transData.transform(np.column_stack([h[:, 0] / w, h[:, 1] / w]))
    return xy[:, 0], xy[:, 1], h[:, 2] / w


def unproject(ax, x, y, depth, M=None):
    """The 3D point drawn at display point (x, y) with the given depth."""
    if M is None:
        M = np.asarray(ax.get_proj(), dtype=float)
    xd, yd = ax.transData.inverted().transform([[float(x), float(y)]])[0]
    v = np.linalg.solve(M, np.array([xd, yd, float(depth), 1.0]))
    if abs(v[3]) < 1e-300:
        v[3] = 1e-300
    return v[:3] / v[3]


def segment_distance(x, y, ax_, ay_, bx, by):
    """Distance from (x, y) to the segments a-b, and the parameter t of the
    nearest point (arrays)."""
    vx, vy = bx - ax_, by - ay_
    L2 = vx * vx + vy * vy
    with np.errstate(invalid="ignore", divide="ignore"):
        t = np.where(L2 > 1e-12, ((x - ax_) * vx + (y - ay_) * vy) / L2, 0.0)
    t = np.clip(t, 0.0, 1.0)
    cx = ax_ + t * vx
    cy = ay_ + t * vy
    return np.hypot(cx - x, cy - y), t


def geometry_key(nodes, members, token=None):
    """A key that changes when the model changes (camera reset)."""
    n = np.ascontiguousarray(np.asarray(nodes, dtype=float))
    m = np.ascontiguousarray(np.asarray(members, dtype=np.int64))
    return (token, n.shape, m.shape, zlib.crc32(n.tobytes()),
            zlib.crc32(m.tobytes()))


def parse_find(text, default_kind="member"):
    """('member' or 'node', number) from text such as 'm 12', 'node 5', '7'.

    Raises ValueError with a plain message for anything else.
    """
    t = str(text).strip().lower()
    m = re.fullmatch(r"(m|mem|member|n|node)?\s*#?\s*(\d+)", t)
    if not m:
        raise ValueError("Type a number, for example m 12 or n 5.")
    word = m.group(1)
    if word is None:
        kind = default_kind
    else:
        kind = "member" if word.startswith("m") else "node"
    return kind, int(m.group(2))


# ========================================================== selection =====
class Selection(QtCore.QObject):
    """The members and nodes picked in the 3D views.

    One object can be shared by several views and tables; ``changed`` is
    emitted after every change.  ``last`` is the item picked last,
    ('member' or 'node', number), shown in the details panel.
    """

    changed = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.members = set()
        self.nodes = set()
        self.last = None
        self.source = None

    def is_empty(self):
        return not self.members and not self.nodes

    def set(self, members=(), nodes=(), last=None, source=None):
        members = {int(k) for k in members}
        nodes = {int(k) for k in nodes}
        if last is not None and not (
                (last[0] == "member" and last[1] in members)
                or (last[0] == "node" and last[1] in nodes)):
            last = None
        if last is None:
            if members:
                last = ("member", max(members))
            elif nodes:
                last = ("node", max(nodes))
        if (members, nodes, last) == (self.members, self.nodes, self.last):
            return False
        self.members, self.nodes, self.last = members, nodes, last
        self.source = source
        self.changed.emit()
        return True

    def toggle(self, kind, idx, source=None):
        members, nodes = set(self.members), set(self.nodes)
        bag = members if kind == "member" else nodes
        if idx in bag:
            bag.discard(idx)
            last = None
        else:
            bag.add(idx)
            last = (kind, idx)
        return self.set(members, nodes, last, source)

    def clear(self, source=None):
        return self.set((), (), None, source)


# ============================================================== scene ======
class Scene:
    """What a 3D view shows, for picking and for the selection details.

    ``view`` is one of model, loads, deflected, axial, utilisation,
    buckling or generic.  ``nodes`` are the node coordinates as drawn
    (deformed in the deflected and buckling views).  The other arguments
    carry the view's own result: ``case`` (combination name),
    ``case_result`` (its StaticResult), ``loads`` (nodal load vector),
    ``phi`` (buckling mode shape, normalised) and ``mode`` (its number).
    """

    def __init__(self, model=None, view="model", nodes=None, members=None,
                 results=None, case=None, case_result=None, loads=None,
                 phi=None, mode=None, key=None, scale=None):
        self.model = model
        self.view = view
        g = model.geom if model is not None else None
        self.nodes = np.asarray(nodes if nodes is not None else g.nodes,
                                dtype=float)
        self.members = np.asarray(members if members is not None
                                  else g.members, dtype=int).reshape(-1, 2)
        self.results = results
        self.case = case
        self.case_result = case_result
        self.loads = loads
        self.phi = phi
        self.mode = mode
        self.scale = scale
        if key is None:
            src = g if g is not None else self
            key = geometry_key(getattr(src, "nodes", self.nodes),
                               getattr(src, "members", self.members),
                               id(model))
        self.key = key
        self.cache = {}

    @property
    def n_members(self):
        return len(self.members)

    @property
    def n_nodes(self):
        return len(self.nodes)


# ----------------------------------------------------------- details ------
def _fmt(v, digits=3, unit=""):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(f):
        return "n/a"
    return f"{f:,.{digits}f}{(' ' + unit) if unit else ''}"


def _case_names(scene):
    res = scene.results
    if res is None:
        return []
    names = list(res.order)
    if scene.case in names:
        names.remove(scene.case)
        names.insert(0, scene.case)
    return names


def _member_view_rows(scene, k):
    """The view's own value of member k, shown first."""
    rows = []
    g = scene.model.geom
    i, j = (int(v) for v in g.members[k])
    if scene.view == "axial" and scene.case_result is not None:
        N = float(scene.case_result.N[k]) / 1e3
        kind = "tension" if N > 0 else ("compression" if N < 0 else "zero")
        rows.append((f"N ({scene.case})", f"{N:+,.2f} kN  {kind}"))
    elif scene.view == "deflected" and scene.case_result is not None:
        d = scene.case_result.displacement_vectors()
        rows.append((f"|u| at ends ({scene.case})",
                     f"{np.linalg.norm(d[i]) * 1e3:.3f} / "
                     f"{np.linalg.norm(d[j]) * 1e3:.3f} mm"))
    elif scene.view == "utilisation" and scene.results is not None:
        res = scene.results
        u = float(res.util_env[k])
        rows.append(("Utilisation", f"{u:.3f}  "
                                    f"{'OK' if u <= 1.0 else 'OVER'}"))
        gov = _governing_case(res, k)
        if gov:
            rows.append(("Governing combination", gov))
    elif scene.view == "buckling" and scene.phi is not None:
        a = np.linalg.norm(scene.phi, axis=1)
        amax = max(float(a.max()), 1e-300)
        rows.append((f"Mode {scene.mode} amplitude",
                     f"{a[i] / amax:.3f} / {a[j] / amax:.3f} (ends, "
                     f"largest = 1)"))
    return rows


def _governing_case(res, k):
    try:
        strength = res.strength_cases()
        return strength[int(res.util_case[k])] if strength else ""
    except Exception:
        return ""


def _node_loads(scene):
    """Nodal load vector (n, 3) [N] of the view's case, or None."""
    if scene.loads is not None:
        return np.asarray(scene.loads, float).reshape(-1, 3)
    if "loads" in scene.cache:
        return scene.cache["loads"]
    F = None
    model = scene.model
    if scene.case and model is not None:
        try:
            from .. import loads as ld
            combo = next((c for c in model.combinations
                          if c.name == scene.case), None)
            if combo is not None:
                F = np.asarray(ld.combine(ld.basic_load_vectors(model),
                                          combo), float).reshape(-1, 3)
        except Exception:
            F = None
    scene.cache["loads"] = F
    return F


def describe(scene, kind, k):
    """(title, rows) with the details of member or node ``k``.

    ``rows`` is a list of (label, value) strings.  The view's own value
    comes first (axial force, displacement, utilisation or mode amplitude).
    """
    model = scene.model
    if model is None:
        if kind == "member":
            i, j = scene.members[k]
            return f"Member {k}", [("End nodes", f"{i} - {j}")]
        p = scene.nodes[k]
        return f"Node {k}", [("x, y, z", f"{p[0]:.3f}, {p[1]:.3f}, "
                                         f"{p[2]:.3f} m")]
    g = model.geom
    res = scene.results
    if kind == "member":
        i, j = (int(v) for v in g.members[k])
        grp = int(g.strut_group[k])
        L = float(g.strut_lengths[k])
        rows = _member_view_rows(scene, k)
        rows.append(("End nodes", f"{i} - {j}"))
        rows.append(("Family", plot.family_label(grp)))
        rows.append(("Length", f"{L:.4f} m"))
        R = g.radius
        if math.isfinite(R) and R > 0.0:
            rows.append(("Chord factor", f"{L / R:.5f}"))
        try:
            sec = model.section_of_group(grp)
            mat = model.material
            rows.append(("Section", str(model.group_sections[grp])))
            rows.append(("Material", str(model.material_name)))
            rows.append(("Area", f"{sec.A * 1e6:,.1f} mm2"))
            rows.append(("Mass", f"{sec.A * L * mat.rho:,.2f} kg"))
        except Exception:
            pass
        if res is not None:
            try:
                rows.append(("KL/r", _fmt(res.caps["slenderness"][k], 1)))
                if scene.view != "utilisation":
                    u = float(res.util_env[k])
                    rows.append(("Utilisation",
                                 f"{u:.3f}  {'OK' if u <= 1.0 else 'OVER'}"))
                    gov = _governing_case(res, k)
                    if gov:
                        rows.append(("Governing combination", gov))
                rows.append(("phiPn compression",
                             _fmt(res.caps["phiPn_c"][k] / 1e3, 1, "kN")))
                rows.append(("phiPn tension",
                             _fmt(res.caps["phiPn_t"][k] / 1e3, 1, "kN")))
                for name in _case_names(scene):
                    if scene.view == "axial" and name == scene.case:
                        continue
                    rows.append((f"N  {name}",
                                 f"{res.cases[name].N[k] / 1e3:+,.2f} kN"))
            except Exception:
                pass
        return f"Member {k}", rows

    # ---- node
    p = np.asarray(g.nodes[k], float)
    rows = []
    if scene.view in ("deflected", "buckling"):
        q = scene.nodes[k]
        rows_first = []
        if scene.view == "deflected" and scene.case_result is not None:
            d = scene.case_result.displacement_vectors()[k] * 1e3
            rows_first.append((f"u ({scene.case})",
                               f"{d[0]:+.3f}, {d[1]:+.3f}, {d[2]:+.3f} mm"))
            rows_first.append(("|u|", f"{np.linalg.norm(d):.3f} mm"))
        if scene.view == "buckling" and scene.phi is not None:
            a = np.linalg.norm(scene.phi, axis=1)
            rows_first.append((f"Mode {scene.mode} amplitude",
                               f"{a[k] / max(float(a.max()), 1e-300):.3f} "
                               f"(largest = 1)"))
        rows_first.append(("Drawn at", f"{q[0]:.3f}, {q[1]:.3f}, "
                                       f"{q[2]:.3f} m"))
        rows += rows_first
    F = _node_loads(scene) if scene.view in ("loads", "deflected",
                                             "axial") else None
    if F is not None and scene.view == "loads":
        f = F[k] / 1e3
        rows.append((f"Load ({scene.case})", f"{f[0]:+.3f}, {f[1]:+.3f}, "
                                             f"{f[2]:+.3f} kN"))
    rows.append(("x, y, z", f"{p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f} m"))
    try:
        ring = int(g.ring_index[k])
        rows.append(("Ring", str(ring)))
    except Exception:
        pass
    deg = int(np.sum(np.asarray(g.members) == k))
    rows.append(("Members at the node", str(deg)))
    try:
        sup = k in set(np.asarray(model.support_nodes(), int).tolist())
    except Exception:
        sup = False
    rows.append(("Support", str(model.support_type) if sup else "no"))
    if F is not None and scene.view != "loads":
        f = F[k] / 1e3
        if np.any(np.abs(f) > 0.0):
            rows.append((f"Load ({scene.case})",
                         f"{f[0]:+.3f}, {f[1]:+.3f}, {f[2]:+.3f} kN"))
    if res is not None:
        try:
            name = (scene.case if scene.case in res.cases
                    else res.service_case)
            if scene.view != "deflected":
                d = res.cases[name].displacement_vectors()[k] * 1e3
                rows.append((f"u ({name})", f"{d[0]:+.3f}, {d[1]:+.3f}, "
                                            f"{d[2]:+.3f} mm"))
            if sup:
                rname = (scene.case if scene.case in res.cases
                         else (res.governing_case or res.order[0]))
                r = res.cases[rname].reaction_vectors()[k] / 1e3
                rows.append((f"Reaction ({rname})",
                             f"{r[0]:+.2f}, {r[1]:+.2f}, {r[2]:+.2f} kN"))
        except Exception:
            pass
    return f"Node {k}", rows


def tooltip(scene, kind, k):
    """One short line for the hover tip."""
    model = scene.model
    if kind == "member":
        if model is None:
            return f"Member {k}"
        g = model.geom
        head = f"Member {k} ({plot.family_label(g.strut_group[k])})"
        i, j = (int(v) for v in g.members[k])
        if scene.view == "axial" and scene.case_result is not None:
            return f"{head}  N = {scene.case_result.N[k] / 1e3:+,.1f} kN"
        if scene.view == "deflected" and scene.case_result is not None:
            d = scene.case_result.displacement_vectors()
            u = 0.5 * (np.linalg.norm(d[i]) + np.linalg.norm(d[j])) * 1e3
            return f"{head}  |u| = {u:.2f} mm"
        if scene.view == "utilisation" and scene.results is not None:
            return f"{head}  utilisation {scene.results.util_env[k]:.3f}"
        if scene.view == "buckling" and scene.phi is not None:
            a = np.linalg.norm(scene.phi, axis=1)
            amax = max(float(a.max()), 1e-300)
            return f"{head}  amplitude {0.5 * (a[i] + a[j]) / amax:.2f}"
        return f"{head}  L = {g.strut_lengths[k]:.3f} m"
    if model is None:
        return f"Node {k}"
    if scene.view == "loads":
        F = _node_loads(scene)
        if F is not None:
            return f"Node {k}  |F| = {np.linalg.norm(F[k]) / 1e3:.2f} kN"
    if scene.view == "deflected" and scene.case_result is not None:
        d = scene.case_result.displacement_vectors()[k]
        return f"Node {k}  |u| = {np.linalg.norm(d) * 1e3:.2f} mm"
    if scene.view == "buckling" and scene.phi is not None:
        a = np.linalg.norm(scene.phi, axis=1)
        return f"Node {k}  amplitude {a[k] / max(float(a.max()), 1e-300):.2f}"
    return f"Node {k}  z = {model.geom.nodes[k][2]:.2f} m"


# ========================================================= controller =====
def _guard(fn):
    """Run a handler; record an error instead of raising it (a raise inside
    a Qt slot or a matplotlib callback can end the program)."""
    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        try:
            return fn(self, *args, **kwargs)
        except Exception:
            self._report(traceback.format_exc())
            return None
    return wrapper


def _modifiers(event):
    """Set of 'ctrl', 'shift', 'alt' held during a matplotlib event."""
    mods = set()
    key = getattr(event, "key", None) or ""
    for part in str(key).lower().split("+"):
        if part in ("control", "ctrl", "cmd"):
            mods.add("ctrl")
        elif part == "shift":
            mods.add("shift")
        elif part == "alt":
            mods.add("alt")
    for name in getattr(event, "modifiers", None) or ():
        name = str(name).lower()
        mods.add("ctrl" if name in ("control", "ctrl", "cmd") else name)
    gui = getattr(event, "guiEvent", None)
    if gui is not None:
        try:
            qm = gui.modifiers()
            KM = Qt.KeyboardModifier
            if qm & KM.ControlModifier:
                mods.add("ctrl")
            if qm & KM.ShiftModifier:
                mods.add("shift")
            if qm & KM.AltModifier:
                mods.add("alt")
        except Exception:
            pass
    return mods


def _set_lim(fn, lo, hi):
    try:
        fn(float(lo), float(hi), view_margin=0.0)
    except TypeError:                     # matplotlib before 3.9
        fn(float(lo), float(hi))


class ViewController(QtCore.QObject):
    """Mouse, keyboard, picking and camera of one 3D view.

    ``widget`` is the :class:`~.canvas.MplCanvas`.  After each drawing, the
    owner calls :meth:`set_scene` (or :meth:`scene_drawn`) with what was
    drawn.  ``picking=False`` turns selection and hover off (animated
    views).
    """

    camera_changed = QtCore.Signal()
    state_changed = QtCore.Signal()
    info_changed = QtCore.Signal()
    message = QtCore.Signal(str)

    def __init__(self, widget, selection=None, picking=True, hover=True,
                 name=""):
        super().__init__(widget)
        self._widget = weakref.ref(widget)
        self.canvas = widget.canvas
        self.figure = widget.figure
        self.name = name
        self.picking = bool(picking)
        self.hover_enabled = bool(hover) and self.picking
        self.errors = []
        self.ax = None
        self.layout = None
        self.scene = None
        self.home = None
        # camera kept across redraws of the same model
        self._key = None
        self._cam = None
        self._user_view = False
        self._user_limits = False
        # display options, kept across redraws
        self.ortho = False
        self.axes_on = True
        self.legend_on = True
        self.member_labels = False
        self.node_labels = False
        self.hidden = None
        # mouse
        self._press = None
        self._last = None
        self._dragging = False
        self._skip_release = False
        # artists
        self._sel_artists = []
        self._hover_artist = None
        self._label_artists = []
        self._label_keys = None
        self._hover_key = None
        self._hover_note = None
        self._axes_hint = False
        self._bg = None
        self._dirty = False
        self._proj = None
        self._base_lw = {}
        self.selection = selection if selection is not None \
            else Selection(self)
        self.selection.changed.connect(self._on_selection_changed)
        c = self.canvas
        self._cids = [c.mpl_connect(name, getattr(self, fn)) for name, fn in (
            ("button_press_event", "_on_press"),
            ("button_release_event", "_on_release"),
            ("motion_notify_event", "_on_motion"),
            ("scroll_event", "_on_scroll"),
            ("key_press_event", "_on_key"),
            ("draw_event", "_on_draw"),
            ("figure_leave_event", "_on_leave"))]
        self._settle = QtCore.QTimer(self)
        self._settle.setSingleShot(True)
        self._settle.setInterval(SETTLE_MS)
        self._settle.timeout.connect(self._on_settled)
        c.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        c.installEventFilter(self)

    # ------------------------------------------------------------ errors --
    def _report(self, text):
        self.errors.append(text)
        del self.errors[:-20]
        try:
            sys.stderr.write(f"GeoDomeLab 3D view: {text}\n")
        except Exception:                            # pragma: no cover
            pass

    def disconnect_all(self):
        """Disconnect the matplotlib callbacks (the window is closing)."""
        for cid in self._cids:
            try:
                self.canvas.mpl_disconnect(cid)
            except Exception:                        # pragma: no cover
                pass
        self._cids = []

    # ---------------------------------------------------------- the axes --
    def _alive(self):
        ax = self.ax
        return (ax is not None and ax in self.figure.axes
                and ax.get_visible())

    def _live(self):
        """The axes is drawn and holds a model, so the camera may move.  An
        emptied view stays still and keeps the model's camera."""
        return self.scene is not None and self._alive()

    def _attach(self, ax):
        if ax is not self.ax:
            self.ax = ax
            self._hover_artist = None
            self._sel_artists = []
            self._label_artists = []
            self._label_keys = None
            self._bg = None
        try:
            ax.disable_mouse_rotation()          # this controller does it
        except Exception:
            pass
        self.layout = plot.axes_layout(ax)
        self._proj = None

    def set_scene(self, scene, ax=None):
        """Tell the viewer what was just drawn (None: nothing to pick).

        Keeps the camera when the model is the same as before, applies the
        display options and puts the selection highlight on top.
        """
        try:
            w = self._widget()
            ax = ax if ax is not None else getattr(w, "axes", None)
            if ax is None:
                return
            self._attach(ax)
            self.scene = scene
            key = None if scene is None else scene.key
            self._restore_camera(key)
            self._seed_toolbar()
            self._decorate()
        except Exception:
            self._report(traceback.format_exc())

    def _seed_toolbar(self):
        """When a kept camera was put back on a fresh drawing, start the
        matplotlib toolbar's history with the drawing's first view, so
        that the toolbar's Home button still goes back to it."""
        tb = getattr(self.canvas, "toolbar", None)
        if tb is None or self.scene is None or not self.home or \
                not (self._user_view or self._user_limits):
            return
        try:
            stack = getattr(tb, "_nav_stack", None)
            if stack is None or stack() is not None:
                return          # the history is not new (same axes redrawn)
            cur = self.camera()
            self._put_camera(self.home)
            tb.push_current()
            self._put_camera(cur)
            tb.push_current()
        except Exception:
            pass

    def _put_camera(self, cam):
        self._set_angles(cam["elev"], cam["azim"], cam["roll"])
        self._set_limits(
            np.array([cam["xlim"][0], cam["ylim"][0], cam["zlim"][0]]),
            np.array([cam["xlim"][1], cam["ylim"][1], cam["zlim"][1]]))

    def before_clear(self):
        """The canvas is about to clear its figure (a redraw): take over a
        zoom or pan of the matplotlib toolbar that no draw has reported
        yet.  (``Figure.clear`` resets the old axes, so it is read now.)"""
        try:
            if self._alive():
                self._note_camera()
        except Exception:
            self._report(traceback.format_exc())

    @_guard
    def update_scene(self, **fields):
        """Change fields of the current scene without drawing again, for
        example ``results`` after an analysis or when the inputs change.
        The camera and the highlight stay and the details are refreshed."""
        sc = self.scene
        if sc is None:
            return
        for name, value in fields.items():
            setattr(sc, name, value)
        sc.cache.clear()
        self._update_info()

    def scene_drawn(self, ax=None, key=None, nodes=None, members=None):
        """:meth:`set_scene` for a plain drawing (no model details)."""
        scene = None
        if nodes is not None and members is not None:
            scene = Scene(None, view="generic", nodes=nodes, members=members,
                          key=key if key is not None
                          else geometry_key(nodes, members))
        self.set_scene(scene, ax=ax)

    # ------------------------------------------------------------ camera --
    def _roll(self, ax=None):
        r = getattr(ax if ax is not None else self.ax, "roll", 0.0)
        return float(r) if r is not None else 0.0

    def camera(self, ax=None):
        """The current camera: angles, limits and projection."""
        ax = ax if ax is not None else self.ax
        return {"elev": float(ax.elev), "azim": float(ax.azim),
                "roll": self._roll(ax),
                "xlim": tuple(map(float, ax.get_xlim())),
                "ylim": tuple(map(float, ax.get_ylim())),
                "zlim": tuple(map(float, ax.get_zlim())),
                "ortho": bool(self.ortho)}

    @staticmethod
    def _same_angles(a, b):
        return all(abs(float(a[k]) - float(b[k])) < 1e-7
                   for k in ("elev", "azim", "roll"))

    @staticmethod
    def _same_limits(a, b):
        span = max(abs(a["xlim"][1] - a["xlim"][0]), 1e-300)
        return all(np.allclose(a[k], b[k], rtol=0.0, atol=1e-9 * span)
                   for k in ("xlim", "ylim", "zlim"))

    def _note_camera(self, ax=None):
        """Take over a camera change made outside the viewer (the
        matplotlib toolbar's zoom, pan, Home, Back or Forward), so that a
        redraw of the same model keeps it."""
        ax = ax if ax is not None else self.ax
        if ax is None or self._cam is None or self.scene is None:
            return
        cam = self.camera(ax)
        if self._same_angles(cam, self._cam) and \
                self._same_limits(cam, self._cam):
            return
        self._cam = cam
        self._proj = None
        h = self.home
        if h:
            self._user_view = not self._same_angles(cam, h)
            self._user_limits = not self._same_limits(cam, h)

    def _set_angles(self, elev, azim, roll=0.0):
        def wrap(a):
            return (float(a) + 180.0) % 360.0 - 180.0
        e, a, r = wrap(elev), wrap(azim), wrap(roll)
        if abs(float(elev) - 90.0) < 1e-9:
            e = 90.0                            # keep the exact top view
        try:
            self.ax.view_init(elev=e, azim=a, roll=r)
        except TypeError:                        # matplotlib < 3.6
            self.ax.view_init(elev=e, azim=a)

    def _limits(self):
        ax = self.ax
        lo = np.array([ax.get_xlim()[0], ax.get_ylim()[0], ax.get_zlim()[0]])
        hi = np.array([ax.get_xlim()[1], ax.get_ylim()[1], ax.get_zlim()[1]])
        return lo, hi

    def _set_limits(self, lo, hi):
        ax = self.ax
        _set_lim(ax.set_xlim, lo[0], hi[0])
        _set_lim(ax.set_ylim, lo[1], hi[1])
        _set_lim(ax.set_zlim, lo[2], hi[2])

    def _restore_camera(self, key):
        self.home = self.camera()
        if key is None:
            # nothing drawn (results cleared): keep the camera of the last
            # model, so that it comes back when that model is drawn again
            return
        cam = self._cam
        if key == self._key and cam is not None:
            if self._user_view:
                self._set_angles(cam["elev"], cam["azim"], cam["roll"])
            if self._user_limits:
                self._set_limits(
                    np.array([cam["xlim"][0], cam["ylim"][0], cam["zlim"][0]]),
                    np.array([cam["xlim"][1], cam["ylim"][1], cam["zlim"][1]]))
        elif key != self._key:
            self._user_view = False
            self._user_limits = False
            self.hidden = None
            self._hover_key = None
        self._key = key
        self._cam = self.camera()

    def _changed(self, draw=True):
        """The camera moved: remember it, redraw (coalesced)."""
        self._cam = self.camera()
        self._proj = None
        if self._hover_key is not None:
            self._clear_hover(blit=False)
        if draw:
            self.canvas.draw_idle()
        self._settle.start()
        self.camera_changed.emit()

    @_guard
    def _on_settled(self):
        if not self._alive():
            return
        if self.member_labels or self.node_labels:
            # the numbers are 3D texts and follow the camera; place them
            # again (and draw) only when another set of them should show
            self._update_labels(only_if_changed=True)
        tb = getattr(self.canvas, "toolbar", None)
        if tb is not None and self._press is None:
            try:
                tb.push_current()
            except Exception:
                pass

    # public camera actions ------------------------------------------------
    @_guard
    def rotate(self, dazim=0.0, delev=0.0, droll=0.0, draw=True):
        if not self._live():
            return
        self._set_angles(self.ax.elev + delev, self.ax.azim + dazim,
                         self._roll() + droll)
        self._user_view = True
        self._changed(draw)

    @_guard
    def view_preset(self, name):
        """Set one of :data:`PRESETS` by name (Iso, Top, Front, ...)."""
        if not self._live():
            return
        for n, e, a in PRESETS:
            if n.lower() == str(name).lower():
                self._set_angles(e, a, 0.0)
                self._user_view = True
                self._changed()
                return

    def _depth_at(self, x, y):
        """Depth for an unprojection at (x, y): the nearest drawn node
        within 40 px, else the centre of the box."""
        ax = self.ax
        M = np.asarray(ax.get_proj(), float)
        if self.scene is not None and self.scene.n_nodes:
            px, py, pz = self._projected()
            d = np.hypot(px - x, py - y)
            i = int(np.argmin(d))
            if d[i] <= 40.0:
                return float(pz[i]), M
        lo, hi = self._limits()
        _, _, z = project(ax, [0.5 * (lo + hi)], M)
        return float(z[0]), M

    @_guard
    def zoom(self, factor, x=None, y=None, draw=True):
        """Scale the view about the point under (x, y) (default: centre).

        ``factor`` < 1 zooms in.  All three ranges are scaled alike, so the
        1:1:1 scale stays, and the point under the cursor stays put.
        """
        if not self._live():
            return
        ax = self.ax
        if x is None or y is None:
            x, y = ax.bbox.x0 + 0.5 * ax.bbox.width, \
                ax.bbox.y0 + 0.5 * ax.bbox.height
        z, M = self._depth_at(x, y)
        P = unproject(ax, x, y, z, M)
        lo, hi = self._limits()
        span = float(np.max(hi - lo))
        ref = self._home_span()
        f = float(factor)
        if ref > 0.0:
            f = min(max(f, 1e-4 * ref / span), 50.0 * ref / span)
        self._set_limits(P + (lo - P) * f, P + (hi - P) * f)
        self._user_limits = True
        self._changed(draw)
        if self.axes_on and not self._axes_hint and ref > 0.0 and \
                span * f < 0.5 * ref:
            self._axes_hint = True
            self.message.emit("Zoomed in. Press A to hide the axes and the "
                              "grid if they get in the way.")

    def _home_span(self):
        h = self.home
        if not h:
            return 0.0
        return float(max(h["xlim"][1] - h["xlim"][0],
                         h["ylim"][1] - h["ylim"][0],
                         h["zlim"][1] - h["zlim"][0]))

    @_guard
    def pan_pixels(self, dx, dy, draw=True):
        """Move the drawing by (dx, dy) display pixels."""
        if not self._live() or (dx == 0 and dy == 0):
            return
        ax = self.ax
        cx = ax.bbox.x0 + 0.5 * ax.bbox.width
        cy = ax.bbox.y0 + 0.5 * ax.bbox.height
        lo, hi = self._limits()
        M = np.asarray(ax.get_proj(), float)
        _, _, z = project(ax, [0.5 * (lo + hi)], M)
        p0 = unproject(ax, cx, cy, z[0], M)
        p1 = unproject(ax, cx + dx, cy + dy, z[0], M)
        D = p1 - p0
        self._set_limits(lo - D, hi - D)
        self._user_limits = True
        self._changed(draw)

    @_guard
    def reset_view(self):
        """Back to the view the drawing started with."""
        if not self._live() or not self.home:
            return
        h = self.home
        self._set_angles(h["elev"], h["azim"], h["roll"])
        self._set_limits(np.array([h["xlim"][0], h["ylim"][0], h["zlim"][0]]),
                         np.array([h["xlim"][1], h["ylim"][1], h["zlim"][1]]))
        self._changed()
        self._user_view = False
        self._user_limits = False

    def _frame_points(self, pts, fill=0.75):
        """Centre the points and zoom so that they fill ``fill`` of the
        drawing (two passes, for the perspective projection)."""
        ax = self.ax
        pts = np.asarray(pts, float).reshape(-1, 3)
        c3 = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
        lo, hi = self._limits()
        half = 0.5 * (hi - lo)
        self._set_limits(c3 - half, c3 + half)
        ref = self._home_span() or float(np.max(hi - lo))
        for _ in range(2):
            px, py, _ = project(ax, pts)
            w = float(np.ptp(px)) if len(px) > 1 else 0.0
            h = float(np.ptp(py)) if len(py) > 1 else 0.0
            aw = max(ax.bbox.width, 1.0) * fill
            ah = max(ax.bbox.height, 1.0) * fill
            k = max(w / aw, h / ah)
            lo, hi = self._limits()
            span = float(np.max(hi - lo))
            k = min(max(k, 0.02 * ref / span), 50.0 * ref / span)
            if abs(k - 1.0) < 0.02:
                break
            self._set_limits(c3 + (lo - c3) * k, c3 + (hi - c3) * k)
        self._user_limits = True
        self._changed()

    def _item_points(self, items):
        sc = self.scene
        pts = []
        for kind, k in items:
            if kind == "member":
                pts.extend(sc.nodes[sc.members[k]])
            else:
                pts.append(sc.nodes[k])
                # the node and its neighbours, for context
                m = np.nonzero((sc.members == k).any(axis=1))[0]
                for e in m[:8]:
                    pts.extend(sc.nodes[sc.members[e]])
        return np.asarray(pts, float).reshape(-1, 3)

    @_guard
    def fit(self):
        """Frame the selection, or every shown member when nothing is
        selected."""
        if not self._alive() or self.scene is None:
            return
        sel = self.selection
        items = [("member", k) for k in sorted(sel.members)
                 if k < self.scene.n_members]
        items += [("node", k) for k in sorted(sel.nodes)
                  if k < self.scene.n_nodes]
        if items:
            pts = self._item_points(items)
        else:
            mask = self._visible_members()
            used = np.unique(self.scene.members[mask]) if mask.any() \
                else np.arange(self.scene.n_nodes)
            pts = self.scene.nodes[used]
        if len(pts):
            self._frame_points(pts)

    # ------------------------------------------------------ display options --
    def _member_collections(self):
        if not self._alive():
            return []
        return [c for c in self.ax.collections
                if c.get_gid() in (plot.MEMBERS_GID, plot.GHOST_GID)]

    @_guard
    def set_projection(self, ortho):
        self.ortho = bool(ortho)
        if self._alive():
            self.ax.set_proj_type("ortho" if self.ortho else "persp")
            self._changed()
        self.state_changed.emit()

    @_guard
    def set_axes_visible(self, on):
        self.axes_on = bool(on)
        if self._alive():
            if self.axes_on:
                self.ax.set_axis_on()
            else:
                self.ax.set_axis_off()
            self._request_draw()
        self.state_changed.emit()

    @_guard
    def set_legend_visible(self, on):
        self.legend_on = bool(on)
        if self._alive() and self.layout is not None:
            self.layout.set_legend_visible(self.ax, self.legend_on)
            self._request_draw()
        self.state_changed.emit()

    def has_legend(self):
        return bool(self._alive() and self.layout is not None
                    and self.layout.has_legend)

    @_guard
    def set_member_labels(self, on):
        self.member_labels = bool(on)
        self._update_labels()
        self.state_changed.emit()

    @_guard
    def set_node_labels(self, on):
        self.node_labels = bool(on)
        self._update_labels()
        self.state_changed.emit()

    def _visible_members(self):
        sc = self.scene
        n = sc.n_members if sc is not None else 0
        if self.hidden is None or len(self.hidden) != n:
            return np.ones(n, dtype=bool)
        return ~self.hidden

    @_guard
    def hide_selected(self):
        sc = self.scene
        if sc is None or not self.selection.members:
            self.message.emit("Select members first, then press H.")
            return
        hidden = ~self._visible_members()
        for k in self.selection.members:
            if k < sc.n_members:
                hidden[k] = True
        self._set_hidden(hidden)
        self.selection.set((), self.selection.nodes, source=self)

    @_guard
    def isolate(self):
        sc = self.scene
        if sc is None or self.selection.is_empty():
            self.message.emit("Select members or nodes first, then press I.")
            return
        keep = np.zeros(sc.n_members, dtype=bool)
        for k in self.selection.members:
            if k < sc.n_members:
                keep[k] = True
        if self.selection.nodes:
            nodes = np.array(sorted(self.selection.nodes))
            keep |= np.isin(sc.members, nodes).any(axis=1)
        self._set_hidden(~keep)

    @_guard
    def show_all(self):
        self._set_hidden(None)

    def _set_hidden(self, hidden):
        sc = self.scene
        if hidden is not None and (sc is None or not np.any(hidden)):
            hidden = None
        self.hidden = None if hidden is None else np.asarray(hidden, bool)
        self._apply_hidden()
        self._proj = None
        self._update_highlight(draw=False)
        self._update_labels(draw=False)
        self._request_draw()
        n = 0 if self.hidden is None else int(self.hidden.sum())
        self.message.emit(f"{n} member(s) hidden. Shift + H shows them."
                          if n else "Every member is shown.")
        self.state_changed.emit()

    def _apply_hidden(self):
        sc = self.scene
        if sc is None:
            return
        vis = self._visible_members()
        for coll in self._member_collections():
            key = id(coll)
            base = self._base_lw.get(key)
            if base is None or base[0] is not coll:
                lw = np.asarray(coll.get_linewidths(), float).ravel()
                if lw.size == 0:
                    lw = np.array([1.0])
                if lw.size != sc.n_members:
                    lw = np.resize(lw, sc.n_members)
                base = (coll, lw)
                self._base_lw[key] = base
            if len(base[1]) == len(vis):
                coll.set_linewidths(np.where(vis, base[1], 0.0))

    # ---------------------------------------------------------- decorate --
    def _decorate(self):
        """Apply the kept options to a fresh drawing."""
        ax = self.ax
        self._base_lw = {}
        try:
            ax.set_proj_type("ortho" if self.ortho else "persp")
        except Exception:
            pass
        if not self.axes_on:
            ax.set_axis_off()
        if self.layout is not None and self.layout.has_legend:
            self.layout.set_legend_visible(ax, self.legend_on)
        sc = self.scene
        if sc is not None and self.hidden is not None and \
                len(self.hidden) != sc.n_members:
            self.hidden = None
        self._apply_hidden()
        self._apply_clip()
        self._make_hover_artist()
        self._update_highlight(draw=False)
        self._update_labels(draw=False)
        self._update_info()
        self.state_changed.emit()

    def _clip_box(self):
        if self.layout is None:
            return None
        return TransformedBbox(self.layout.region, self.figure.transFigure)

    def _apply_clip(self, artists=None):
        """Clip the model to the drawing region (not to the square axes),
        so that a zoomed model uses the whole width right of the legend."""
        clip = self._clip_box()
        if clip is None or not self._alive():
            return
        arts = artists if artists is not None else \
            list(self.ax.collections) + list(self.ax.lines)
        for a in arts:
            try:
                a.set_clip_box(clip)
            except Exception:
                pass

    def _region_px(self):
        """The drawing region in display pixels."""
        if self.layout is not None:
            r = self.layout.region
        else:
            r = self.ax.get_position(original=True)
        return Bbox(self.figure.transFigure.transform(
            [[r.x0, r.y0], [r.x1, r.y1]]))

    # ------------------------------------------------------------ picking --
    def _projected(self):
        """Projected node coordinates, cached until the view moves."""
        ax = self.ax
        sc = self.scene
        M = np.asarray(ax.get_proj(), float)
        key = (id(ax), id(sc), ax.bbox.bounds)
        c = self._proj
        if c is not None and c[0] == key and np.array_equal(c[1], M):
            return c[2]
        out = project(ax, sc.nodes, M)
        self._proj = (key, M, out)
        return out

    def pick(self, x, y):
        """('member' or 'node', number, distance in px) under display point
        (x, y), or None.

        The nearest member (nearer ones win ties) is the answer, unless a
        node is nearer, or a node lies within a snap radius of the cursor.
        Every member that meets at a node passes through it, so without
        the snap radius a node could only be picked at its exact centre.
        The radius is :data:`NODE_SNAP_PX`, but at most 30 % of the
        node's shortest member on screen, so a short member can still be
        picked between its ends.  A node behind the member under the
        cursor (and not one of its ends) does not snap.
        """
        sc = self.scene
        if sc is None or not self._alive() or x is None or y is None:
            return None
        if not self._region_px().contains(x, y):
            return None
        px, py, pz = self._projected()
        zlo, zhi = float(np.min(pz)), float(np.max(pz))
        zr = max(zhi - zlo, 1e-300)
        best = None
        member_depth = None
        vis = self._visible_members()
        if sc.n_members:
            i, j = sc.members[:, 0], sc.members[:, 1]
            d, t = segment_distance(x, y, px[i], py[i], px[j], py[j])
            depth = (pz[i] + t * (pz[j] - pz[i]) - zlo) / zr
            score = np.where(vis & (d <= MEMBER_PICK_PX),
                             d + DEPTH_BIAS_PX * depth, np.inf)
            k = int(np.argmin(score))
            if np.isfinite(score[k]):
                best = ("member", k, float(d[k]))
                member_depth = float(depth[k])
        dn = np.hypot(px - x, py - y)
        node_ok = dn <= NODE_PICK_PX
        if self.hidden is not None and sc.n_members:
            used = np.zeros(sc.n_nodes, dtype=bool)
            used[np.unique(sc.members[vis])] = True
            node_ok &= used
        if node_ok.any():
            score = np.where(node_ok, dn + DEPTH_BIAS_PX * (pz - zlo) / zr,
                             np.inf)
            k = int(np.argmin(score))
            dk = float(dn[k])
            if best is None or dk <= best[2] + 0.5:
                return ("node", k, dk)
            if dk <= self._snap_radius(k, px, py, vis):
                m = sc.members[best[1]]
                end = k == int(m[0]) or k == int(m[1])
                behind = (pz[k] - zlo) / zr > member_depth + 0.02
                if end or not behind:
                    return ("node", k, dk)
        return best

    def _snap_radius(self, k, px, py, vis):
        """Snap radius of node k [px] (see :meth:`pick`)."""
        sc = self.scene
        at = np.nonzero(vis & (sc.members == k).any(axis=1))[0]
        if not len(at):
            return NODE_SNAP_PX
        m = sc.members[at]
        L = np.hypot(px[m[:, 1]] - px[m[:, 0]], py[m[:, 1]] - py[m[:, 0]])
        return float(min(NODE_SNAP_PX, max(1.5, 0.3 * float(L.min()))))

    # ------------------------------------------------------------- mouse --
    def _toolbar_busy(self):
        tb = getattr(self.canvas, "toolbar", None)
        mode = getattr(tb, "mode", None)
        if mode is None:
            return False
        return bool(str(getattr(mode, "value", mode)))

    def _xy(self, event):
        """The event position in display pixels, as floats.

        matplotlib rounds the position down to a whole pixel.  On a figure
        whose height is not a whole number of pixels that moves a click by
        one pixel, so the exact position is read from the Qt event.
        """
        x, y = getattr(event, "x", None), getattr(event, "y", None)
        if x is None or y is None:
            return None, None
        gui = getattr(event, "guiEvent", None)
        if gui is not None:
            try:
                if hasattr(gui, "position"):          # Qt 6, Qt 5 wheel
                    pos = gui.position()
                elif hasattr(gui, "localPos"):        # Qt 5 mouse
                    pos = gui.localPos()
                else:
                    pos = gui.pos()
                r = float(getattr(self.canvas, "device_pixel_ratio", 1.0)
                          or 1.0)
                h = float(self.figure.bbox.height)
                fx = float(pos.x()) * r
                fy = (h / r - float(pos.y())) * r
                if abs(fx - x) <= 2.0 and abs(fy - y) <= 2.0:
                    return fx, fy
            except Exception:
                pass
        return float(x), float(y)

    def _accept(self, event, xy=None):
        x, y = xy if xy is not None else self._xy(event)
        if not self._alive() or x is None or y is None:
            return False
        if event.inaxes is not None and event.inaxes is not self.ax:
            return False                     # e.g. the colour bar
        return self._region_px().contains(x, y)

    def _quiet_axes(self):
        """Stop matplotlib's own 3D mouse handling of this press.  Before
        matplotlib 3.7, disable_mouse_rotation() left a built-in middle
        button pan on, which panned a second time."""
        try:
            self.ax.button_pressed = None
        except Exception:
            pass

    @_guard
    def _on_press(self, event):
        xy = self._xy(event)
        if self._toolbar_busy() or not self._accept(event, xy):
            self._press = None
            return
        x, y = xy
        button = int(getattr(event, "button", 0) or 0)
        mods = _modifiers(event)
        self._press = (x, y, button, mods)
        self._last = (x, y)
        self._dragging = False
        self._quiet_axes()
        self._clear_hover()
        try:
            self.canvas.setFocus()
        except Exception:
            pass
        if getattr(event, "dblclick", False) and button == 1 and \
                self.picking:
            self._skip_release = True
            hit = self.pick(x, y)
            if hit is not None:
                kind, k, _ = hit
                if "ctrl" not in mods:
                    self.selection.set([k] if kind == "member" else [],
                                       [k] if kind == "node" else [],
                                       (kind, k), source=self)
                self._frame_points(self._item_points([(kind, k)]), fill=0.5)

    @_guard
    def _on_motion(self, event):
        p = self._press
        if p is None:
            if self.hover_enabled:
                self._hover(event)
            return
        x, y = self._xy(event)
        if x is None or y is None or self._toolbar_busy():
            return
        self._quiet_axes()
        if not self._dragging:
            if math.hypot(x - p[0], y - p[1]) < DRAG_PX:
                return
            self._dragging = True
        dx, dy = x - self._last[0], y - self._last[1]
        self._last = (x, y)
        if dx == 0 and dy == 0:
            return
        button, mods = p[2], p[3]
        if button == 1 and "shift" not in mods:
            self._drag_rotate(dx, dy)
        elif button in (1, 2, 3):
            self.pan_pixels(dx, dy)

    def _drag_rotate(self, dx, dy):
        ax = self.ax
        w = max(ax.bbox.width, 1.0)
        h = max(ax.bbox.height, 1.0)
        roll = math.radians(self._roll())
        delev = -(dy / h) * 180.0 * math.cos(roll) + \
            (dx / w) * 180.0 * math.sin(roll)
        dazim = -(dy / h) * 180.0 * math.sin(roll) - \
            (dx / w) * 180.0 * math.cos(roll)
        self.rotate(dazim, delev)

    @_guard
    def _on_release(self, event):
        p, self._press = self._press, None
        if self._skip_release:
            self._skip_release = False
            return
        if p is None:
            return
        if self._dragging:
            self._dragging = False
            self._settle.start()
            return
        if p[2] != 1 or not self.picking or self.scene is None:
            return
        self._click(p[0], p[1], "ctrl" in p[3])

    def _click(self, x, y, ctrl=False):
        hit = self.pick(x, y)
        sel = self.selection
        if hit is None:
            if not ctrl:
                sel.clear(source=self)
            return
        kind, k, _ = hit
        if ctrl:
            sel.toggle(kind, k, source=self)
        else:
            sel.set([k] if kind == "member" else [],
                    [k] if kind == "node" else [], (kind, k), source=self)

    @staticmethod
    def _wheel_notches(event):
        """Wheel notches of a scroll event (fractions for a touchpad).

        The Qt event's angle delta is used (120 per notch).  matplotlib
        passes the pixel delta instead where Qt has one, and a pixel is not
        a notch.  The result is limited to :data:`WHEEL_MAX_NOTCHES`.
        """
        s = None
        gui = getattr(event, "guiEvent", None)
        if gui is not None:
            try:
                a = float(gui.angleDelta().y())
                if a:
                    s = a / 120.0
                else:
                    p = float(gui.pixelDelta().y())
                    s = p / 50.0 if p else None
            except Exception:
                s = None
        if s is None:
            s = float(getattr(event, "step", 0.0) or 0.0)
            if s == 0.0:
                s = 1.0 if getattr(event, "button", "") == "up" else -1.0
        return max(-WHEEL_MAX_NOTCHES, min(WHEEL_MAX_NOTCHES, s))

    @_guard
    def _on_scroll(self, event):
        xy = self._xy(event)
        if not self._accept(event, xy):
            return
        s = self._wheel_notches(event)
        if s:
            self.zoom(WHEEL_FACTOR ** (-s), xy[0], xy[1])

    @_guard
    def _on_leave(self, _event):
        self._clear_hover()

    # ------------------------------------------------------------- keys --
    @_guard
    def _on_key(self, event):
        key = getattr(event, "key", None)
        if not key or not self._alive():
            return
        if key in ("+", "=", "pageup"):
            self.zoom(1.0 / KEY_ZOOM)
            return
        if key in ("-", "pagedown"):
            self.zoom(KEY_ZOOM)
            return
        parts = key.split("+")
        name, mods = parts[-1], set(parts[:-1])
        if name in ("left", "right", "up", "down"):
            sign = 1.0 if name in ("right", "up") else -1.0
            if "ctrl" in mods:
                ax = self.ax
                step = PAN_STEP * min(ax.bbox.width, ax.bbox.height)
                if name in ("left", "right"):
                    self.pan_pixels(sign * step, 0.0)
                else:
                    self.pan_pixels(0.0, sign * step)
                return
            step = ROTATE_STEP_BIG if "shift" in mods else ROTATE_STEP
            if name in ("left", "right"):
                self.rotate(dazim=sign * step)
            else:
                self.rotate(delev=sign * step)
            return
        if mods == {"ctrl"} and name == "c":
            self.copy_selection()
            return
        if mods - {"shift"}:
            return                          # leave other shortcuts alone
        if name in ("home", "r", "R"):
            self.reset_view()
        elif name in ("f", "F"):
            self.fit()
        elif name in ("1", "2", "3", "4", "5", "6"):
            self.view_preset(PRESETS[int(name) - 1][0])
        elif name == "escape":
            self.selection.clear(source=self)
        elif name in ("l", "L"):
            self.set_legend_visible(not self.legend_on)
        elif name in ("m", "M") and self.picking:
            # (not in an animated view: the numbers would stay behind)
            self.set_member_labels(not self.member_labels)
        elif name in ("n", "N") and self.picking:
            self.set_node_labels(not self.node_labels)
        elif name in ("a", "A"):
            self.set_axes_visible(not self.axes_on)
        elif name in ("o", "O", "p", "P"):
            self.set_projection(not self.ortho)
        elif name == "h" and self.picking:
            self.hide_selected()
        elif name == "H" and self.picking:
            self.show_all()
        elif name in ("i", "I") and self.picking:
            self.isolate()
        elif name in ("q", "Q", "e", "E"):
            step = ROTATE_STEP_BIG if name.isupper() else ROTATE_STEP
            self.rotate(droll=step if name in ("e", "E") else -step)

    # ------------------------------------------------------------ hover --
    def _make_hover_artist(self):
        if self._hover_artist is not None:
            self._remove([self._hover_artist])
            self._hover_artist = None
        if not (self.hover_enabled and self.scene is not None):
            self._hover_artist = None
            return
        if self.scene.n_members > HOVER_MAX_MEMBERS:
            self._hover_artist = None
            if self._hover_note != self.scene.key:     # once per model
                self._hover_note = self.scene.key
                self.message.emit(
                    f"Hover is off above {HOVER_MAX_MEMBERS:,} members. "
                    f"Click a member to see its details.")
            return
        line = art3d.Line3D([], [], [], color=HOVER_COLOUR, lw=4.0,
                            alpha=0.9, solid_capstyle="round", zorder=60)
        line.set_animated(True)
        line.set_visible(False)
        self.ax.add_line(line)
        self._apply_clip([line])
        self._hover_artist = line
        self._hover_key = None
        self._bg = None

    def _hover(self, event):
        art = self._hover_artist
        if art is None or self.scene is None or self._toolbar_busy() or \
                not self._alive():
            return
        xy = self._xy(event)
        hit = self.pick(*xy) if self._accept(event, xy) else None
        key = None if hit is None else (hit[0], hit[1])
        if key == self._hover_key:
            return
        self._hover_key = key
        sc = self.scene
        if key is None:
            art.set_visible(False)
            QtWidgets.QToolTip.hideText()
        else:
            kind, k = key
            if kind == "member":
                p = sc.nodes[sc.members[k]]
                art.set_data_3d(p[:, 0], p[:, 1], p[:, 2])
                art.set_linestyle("-")
                art.set_marker("None")
            else:
                p = sc.nodes[k]
                art.set_data_3d([p[0]], [p[1]], [p[2]])
                art.set_linestyle("None")
                art.set_marker("o")
                art.set_markersize(11)
                art.set_markerfacecolor("none")
                art.set_markeredgecolor(HOVER_COLOUR)
                art.set_markeredgewidth(2.5)
            art.set_visible(True)
            self._show_tip(event, tooltip(sc, kind, k))
        self._blit_hover()

    def _show_tip(self, event, text):
        try:
            gui = getattr(event, "guiEvent", None)
            pos = None
            if gui is not None:
                if hasattr(gui, "globalPosition"):
                    pos = gui.globalPosition().toPoint()
                elif hasattr(gui, "globalPos"):
                    pos = gui.globalPos()
            if pos is None:
                r = float(getattr(self.canvas, "device_pixel_ratio", 1.0)
                          or 1.0)
                h = self.figure.bbox.height
                pos = self.canvas.mapToGlobal(QtCore.QPoint(
                    int(event.x / r), int((h - event.y) / r)))
            QtWidgets.QToolTip.showText(pos + QtCore.QPoint(14, 18), text,
                                        self.canvas)
        except Exception:
            pass

    def _clear_hover(self, blit=True):
        if self._hover_key is None:
            return
        self._hover_key = None
        art = self._hover_artist
        try:
            QtWidgets.QToolTip.hideText()
        except Exception:
            pass
        if art is not None:
            art.set_visible(False)
            if blit:
                self._blit_hover()

    def _blit_hover(self):
        art = self._hover_artist
        if art is None or not self._alive():
            return
        c = self.canvas
        if self._bg is not None:
            try:
                c.restore_region(self._bg)
                if art.get_visible():
                    self.ax.draw_artist(art)
                c.blit(self.figure.bbox)
                return
            except Exception:
                self._bg = None
        c.draw_idle()

    @_guard
    def _on_draw(self, event):
        c = self.canvas
        if getattr(event, "canvas", c) is not c:
            return              # a figure export (PDF, SVG) with its canvas
        try:
            if c.is_saving():
                return          # a PNG export: no hover, no background
        except Exception:
            pass
        if self._alive():
            self._note_camera()
        art = self._hover_artist
        if art is None or not self._alive():
            self._bg = None
            return
        try:
            self._bg = self.canvas.copy_from_bbox(self.figure.bbox)
        except Exception:
            self._bg = None
            return
        if art.get_visible():
            self.ax.draw_artist(art)

    # --------------------------------------------------------- selection --
    @_guard
    def _on_selection_changed(self):
        self._update_highlight(draw=False)
        if self.member_labels or self.node_labels:
            self._update_labels(draw=False)
        self._update_info()
        self._request_draw()

    def _request_draw(self):
        try:
            visible = self.canvas.isVisible()
        except RuntimeError:                 # deleted canvas
            return
        if visible:
            self._dirty = False
            self.canvas.draw_idle()
        else:
            self._dirty = True

    def eventFilter(self, obj, event):               # noqa: N802 (Qt name)
        # never raise here (Qt calls it while the window is torn down)
        try:
            if event.type() == QtCore.QEvent.Type.Show and self._dirty:
                self._dirty = False
                self.canvas.draw_idle()
        except Exception:
            pass
        return False

    def _remove(self, artists):
        for a in artists:
            try:
                a.remove()
            except Exception:
                pass

    def _update_highlight(self, draw=True):
        self._remove(self._sel_artists)
        self._sel_artists = []
        sc = self.scene
        if sc is None or not self._alive() or not self.picking:
            return
        vis = self._visible_members()
        mem = [k for k in sorted(self.selection.members)
               if k < sc.n_members and vis[k]]
        nod = [k for k in sorted(self.selection.nodes) if k < sc.n_nodes]
        ax = self.ax
        arts = []
        if mem:
            seg = sc.nodes[sc.members[mem]]                  # (k, 2, 3)
            gap = np.full((len(mem), 1, 3), np.nan)
            xyz = np.concatenate([seg, gap], axis=1).reshape(-1, 3)
            for lw, col, z in ((6.0, SELECT_EDGE, 50), (3.2, SELECT_COLOUR,
                                                           51)):
                ln = art3d.Line3D(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=col,
                                  lw=lw, solid_capstyle="round", zorder=z)
                ax.add_line(ln)
                arts.append(ln)
        if nod:
            p = sc.nodes[nod]
            ln = art3d.Line3D(p[:, 0], p[:, 1], p[:, 2], linestyle="None",
                              marker="o", markersize=11,
                              markerfacecolor=SELECT_COLOUR,
                              markeredgecolor=SELECT_EDGE,
                              markeredgewidth=1.6, zorder=52)
            ax.add_line(ln)
            arts.append(ln)
        self._apply_clip(arts)
        self._sel_artists = arts
        if draw:
            self._request_draw()

    # ------------------------------------------------------------ labels --
    def _label_candidates(self):
        """The numbers that could be shown, best first: arrays (kind, k,
        x, y) with kind 0 for a member and 1 for a node.  Selected items
        come first, then the items nearest to the centre of the view."""
        sc = self.scene
        region = self._region_px()
        cx, cy = 0.5 * (region.x0 + region.x1), 0.5 * (region.y0 + region.y1)
        px, py, _ = self._projected()
        vis = self._visible_members()
        parts = []
        if self.member_labels and sc.n_members:
            i, j = sc.members[:, 0], sc.members[:, 1]
            mx, my = 0.5 * (px[i] + px[j]), 0.5 * (py[i] + py[j])
            ok = vis & (mx >= region.x0) & (mx <= region.x1) & \
                (my >= region.y0) & (my <= region.y1)
            k = np.nonzero(ok)[0]
            sel = np.isin(k, np.fromiter(self.selection.members, int))
            parts.append((np.zeros(len(k), int), k, mx[k], my[k], sel))
        if self.node_labels and sc.n_nodes:
            used = np.ones(sc.n_nodes, dtype=bool)
            if self.hidden is not None and sc.n_members:
                used[:] = False
                used[np.unique(sc.members[vis])] = True
            ok = used & (px >= region.x0) & (px <= region.x1) & \
                (py >= region.y0) & (py <= region.y1)
            k = np.nonzero(ok)[0]
            sel = np.isin(k, np.fromiter(self.selection.nodes, int))
            parts.append((np.ones(len(k), int), k, px[k], py[k], sel))
        if not parts:
            e = np.zeros(0)
            return e.astype(int), e.astype(int), e, e, e.astype(bool)
        kind, k, x, y, sel = (np.concatenate(a) for a in zip(*parts))
        order = np.lexsort((np.hypot(x - cx, y - cy), ~sel))
        return kind[order], k[order], x[order], y[order], sel[order]

    def _update_labels(self, draw=True, only_if_changed=False):
        sc = self.scene
        if sc is None or not self._alive() or not (self.member_labels or
                                                   self.node_labels):
            self._remove(self._label_artists)
            self._label_artists = []
            self._label_keys = None
            if draw:
                self._request_draw()
            return
        ax = self.ax
        kinds, ks, xs, ys, sels = self._label_candidates()
        # keep the numbers apart on screen: one is skipped when one already
        # placed is nearer than LABEL_GAP_PX (grid buckets keep it fast)
        placed, grid = [], {}
        cell = LABEL_GAP_PX
        for kind, k, x, y, sel in zip(kinds.tolist(), ks.tolist(),
                                      xs.tolist(), ys.tolist(),
                                      sels.tolist()):
            gx, gy = int(x // cell), int(y // cell)
            near = False
            if not sel:                          # selected items always show
                for ix in (gx - 1, gx, gx + 1):
                    for iy in (gy - 1, gy, gy + 1):
                        for (qx, qy) in grid.get((ix, iy), ()):
                            if (qx - x) ** 2 + (qy - y) ** 2 < cell * cell:
                                near = True
                                break
                        if near:
                            break
                    if near:
                        break
            if near:
                continue
            grid.setdefault((gx, gy), []).append((x, y))
            placed.append(("member" if kind == 0 else "node", int(k)))
            if len(placed) >= LABEL_CAP:
                break
        keys = tuple(placed)
        if only_if_changed and self._label_keys is not None and \
                keys == self._label_keys:
            return                  # the same numbers: no new drawing
        self._remove(self._label_artists)
        self._label_artists = []
        self._label_keys = keys
        arts = []
        for kind, k in placed:
            if kind == "member":
                p = sc.nodes[sc.members[k]].mean(axis=0)
                t = ax.text(p[0], p[1], p[2], str(k), fontsize=7,
                            color=MEMBER_LABEL_COLOUR, ha="center",
                            va="center", zorder=55,
                            bbox=dict(boxstyle="square,pad=0.08", fc="white",
                                      ec="none", alpha=0.7))
            else:
                p = sc.nodes[k]
                t = ax.text(p[0], p[1], p[2], f" {k}", fontsize=7,
                            color=NODE_LABEL_COLOUR, ha="left",
                            va="bottom", zorder=56)
            t.set_gid("geodome-label")
            t.set_clip_on(True)
            arts.append(t)
        self._apply_clip(arts)
        self._label_artists = arts
        if len(placed) < len(ks):
            self.message.emit(f"{len(placed)} of {len(ks)} numbers "
                              f"shown. Zoom in to see the others.")
        if draw:
            self._request_draw()

    # -------------------------------------------------------------- info --
    def _update_info(self):
        self.info_changed.emit()

    def selection_details(self):
        """(title, rows, summary) of the selection for the info panel."""
        sc = self.scene
        sel = self.selection
        if sc is None or not self._alive():
            return "", [], "Nothing drawn here."
        n_m = len([k for k in sel.members if k < sc.n_members])
        n_n = len([k for k in sel.nodes if k < sc.n_nodes])
        parts = []
        if n_m:
            parts.append(f"{n_m} member{'s' if n_m != 1 else ''}")
        if n_n:
            parts.append(f"{n_n} node{'s' if n_n != 1 else ''}")
        summary = (" and ".join(parts) + " selected.") if parts else \
            "Click a member or a node to see its details."
        last = sel.last
        if last is None:
            return "", [], summary
        kind, k = last
        n = sc.n_members if kind == "member" else sc.n_nodes
        if not 0 <= k < n:
            return "", [], summary
        title, rows = describe(sc, kind, k)
        return title, rows, summary

    def info_html(self):
        title, rows, summary = self.selection_details()
        out = [f"<p style='color:#43566a'>{html.escape(summary)}</p>"]
        if title:
            out.append(f"<p><b style='color:{st.ACCENT_DARK}'>"
                       f"{html.escape(title)}</b></p>")
            out.append("<table cellspacing='0' cellpadding='2'>")
            for i, (a, b) in enumerate(rows):
                bg = "#f3f6f9" if i % 2 else "#ffffff"
                out.append(f"<tr style='background:{bg}'><td>"
                           f"{html.escape(a)}</td><td><b>{html.escape(b)}"
                           f"</b></td></tr>")
            out.append("</table>")
        return "".join(out)

    def info_text(self):
        """The selection details as plain text (Ctrl + C)."""
        title, rows, summary = self.selection_details()
        lines = [summary]
        if title:
            lines.append(title)
            lines += [f"{a}: {b}" for a, b in rows]
        sc = self.scene
        sel = self.selection
        if sc is not None and (len(sel.members) + len(sel.nodes)) > 1:
            if sel.members:
                lines.append("Members: " + ", ".join(
                    str(k) for k in sorted(sel.members)))
            if sel.nodes:
                lines.append("Nodes: " + ", ".join(
                    str(k) for k in sorted(sel.nodes)))
        return "\n".join(lines)

    @_guard
    def copy_selection(self):
        text = self.info_text()
        cb = QtWidgets.QApplication.clipboard()
        if cb is not None:
            cb.setText(text)
        self.message.emit("Selection details copied.")
        return text

    @_guard
    def find(self, text, default_kind="member"):
        """Select the member or node named in ``text`` and frame it.

        Returns True when found.  Emits :attr:`message` otherwise.
        """
        sc = self.scene
        if sc is None or not self._alive():
            self.message.emit("Nothing is drawn in this view yet.")
            return False
        try:
            kind, k = parse_find(text, default_kind)
        except ValueError as exc:
            self.message.emit(str(exc))
            return False
        n = sc.n_members if kind == "member" else sc.n_nodes
        if not 0 <= k < n:
            self.message.emit(f"There is no {kind} {k}. The numbers run from "
                              f"0 to {n - 1}.")
            return False
        if kind == "member" and not self._visible_members()[k]:
            self.show_all()
        self.selection.set([k] if kind == "member" else [],
                           [k] if kind == "node" else [], (kind, k),
                           source=self)
        self._frame_points(self._item_points([(kind, k)]), fill=0.5)
        self.message.emit(f"{kind.capitalize()} {k} selected.")
        return True


# ============================================================ widgets =====
_BAR_CSS = f"""
QToolBar#view3dBar {{
    background: #f7f9fb; border: none; border-bottom: 1px solid #dde4ea;
    spacing: 3px; padding: 2px 4px;
}}
QToolBar#view3dBar QToolButton {{
    background: #ffffff; border: 1px solid #c3cfda; border-radius: 3px;
    padding: 2px 5px; color: #1c2733;
}}
QToolBar#view3dBar QToolButton:hover {{ background: #eef4f9;
                                        border-color: {st.ACCENT}; }}
QToolBar#view3dBar QToolButton:pressed {{ background: #dfe9f1; }}
QToolBar#view3dBar QToolButton:checked {{
    background: #dbe8f3; border-color: {st.ACCENT}; color: {st.ACCENT_DARK};
    font-weight: 600;
}}
QToolBar#view3dBar QToolButton:disabled {{ color: #a7b2bd;
                                           border-color: #dde4ea; }}
QToolBar#view3dBar QLineEdit {{ padding: 1px 4px; min-height: 16px; }}
QToolBar#view3dBar QComboBox {{ padding: 1px 4px; min-height: 16px; }}
QFrame#view3dInfo {{ background: #ffffff; border-left: 1px solid #c8d2db; }}
QFrame#view3dInfo QTextBrowser {{
    border: none; background: #ffffff;
    font-family: "Segoe UI", "Calibri", sans-serif; font-size: 11px;
}}
"""


class ViewBar(QtWidgets.QToolBar):
    """The row of view buttons over a 3D drawing.

    It is a QToolBar, so on a narrow window the buttons that do not fit go
    into the >> menu at its right end instead of widening the window.
    """

    def __init__(self, controller, picking=True, parent=None):
        super().__init__("3D view", parent)
        self.setObjectName("view3dBar")
        self.setStyleSheet(_BAR_CSS)
        self.setMovable(False)
        self.setFloatable(False)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.controller = controller
        self.info_panel = None
        self._help = None
        ctl = controller
        self.preset_actions = {}
        for name, tip in (("Iso", "Isometric view (key 1)"),
                          ("Top", "View from the top (key 2)"),
                          ("Front", "View from the front, along +y (key 3)"),
                          ("Side", "View from the right, along -x (key 6)")):
            a = self._action(name, tip, self._on_preset)
            a.setData("Right" if name == "Side" else name)
            self.preset_actions[name] = a
        self.addSeparator()
        self.zoom_in_act = self._action("+", "Zoom in (wheel up, + or Page "
                                             "Up)", self._zoom_in)
        self.zoom_out_act = self._action("−", "Zoom out (wheel down, - "
                                                   "or Page Down)",
                                         self._zoom_out)
        self.fit_act = self._action("Fit", "Fit the selection, or "
                                           "everything (F)", self._fit)
        self.reset_act = self._action("Reset", "Back to the first view (Home "
                                               "or R)", self._reset)
        self.addSeparator()
        self.ortho_act = self._action("Ortho", "Orthographic projection (O "
                                               "or P). Off: perspective.",
                                      self._on_ortho, True)
        self.axes_act = self._action("Axes", "Axes and grid (A)",
                                     self._on_axes, True)
        self.legend_act = self._action("Legend", "Show the legend (L)",
                                       self._on_legend, True)
        self.member_act = self.node_act = self.info_act = None
        self.find_edit = self.find_kind = None
        if picking:
            self.member_act = self._action("Mem #", "Member numbers (M)",
                                           self._on_member_labels, True)
            self.node_act = self._action("Node #", "Node numbers (N)",
                                         self._on_node_labels, True)
            self.addSeparator()
            self.info_act = self._action("Info", "Details of the selected "
                                                 "member or node",
                                         self._on_info, True)
            self.find_kind = QtWidgets.QComboBox()
            self.find_kind.addItem("m", "member")
            self.find_kind.addItem("n", "node")
            self.find_kind.setToolTip("What a plain number in Find means: "
                                      "m = member, n = node.")
            self.find_edit = QtWidgets.QLineEdit()
            self.find_edit.setPlaceholderText("Find: m 12")
            self.find_edit.setToolTip(
                "Type a member or node number and press Enter: m 12, n 5, "
                "or a plain number (m or n on the left).")
            self.find_edit.setFixedWidth(84)
            self.find_edit.returnPressed.connect(self._on_find)
            self.addWidget(self.find_kind)
            self.addWidget(self.find_edit)
        self.help_act = self._action("?", "Mouse and keyboard controls",
                                     self.show_help)
        self.msg = QtWidgets.QLabel("")
        self.msg.setObjectName("hint")
        self.msg.setContentsMargins(6, 0, 0, 0)
        self.msg.setMinimumWidth(0)
        self.msg.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                               QtWidgets.QSizePolicy.Preferred)
        self.addWidget(self.msg)
        self._msg_timer = QtCore.QTimer(self)
        self._msg_timer.setSingleShot(True)
        self._msg_timer.setInterval(6000)
        self._msg_timer.timeout.connect(self._clear_message)
        ctl.state_changed.connect(self.sync)
        ctl.message.connect(self.show_message)
        self.sync()

    def _action(self, text, tip, slot, checkable=False):
        a = QtWidgets.QAction(text, self)
        a.setToolTip(tip)
        a.setCheckable(checkable)
        if checkable:
            a.toggled.connect(slot)
        else:
            a.triggered.connect(slot)
        self.addAction(a)
        w = self.widgetForAction(a)
        if w is not None:
            w.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        return a

    # ------------------------------------------------------------- slots --
    def _call(self, fn, *args):
        try:
            return fn(*args)
        except Exception:
            self.controller._report(traceback.format_exc())
            return None

    def _on_preset(self, *_):
        a = self.sender()
        name = a.data() if a is not None else "Iso"
        self._call(self.controller.view_preset, name or "Iso")

    def _zoom_in(self, *_):
        self._call(self.controller.zoom, 1.0 / KEY_ZOOM)

    def _zoom_out(self, *_):
        self._call(self.controller.zoom, KEY_ZOOM)

    def _fit(self, *_):
        self._call(self.controller.fit)

    def _reset(self, *_):
        self._call(self.controller.reset_view)

    def _on_ortho(self, on):
        if bool(on) != self.controller.ortho:
            self._call(self.controller.set_projection, bool(on))

    def _on_axes(self, on):
        if bool(on) != self.controller.axes_on:
            self._call(self.controller.set_axes_visible, bool(on))

    def _on_legend(self, on):
        if bool(on) != self.controller.legend_on:
            self._call(self.controller.set_legend_visible, bool(on))

    def _on_member_labels(self, on):
        if bool(on) != self.controller.member_labels:
            self._call(self.controller.set_member_labels, bool(on))

    def _on_node_labels(self, on):
        if bool(on) != self.controller.node_labels:
            self._call(self.controller.set_node_labels, bool(on))

    def _on_info(self, on):
        p = self.info_panel
        if p is not None:
            self._call(p.set_open, bool(on), True)

    def _on_find(self):
        kind = self.find_kind.currentData() if self.find_kind else "member"
        ok = self._call(self.controller.find, self.find_edit.text(), kind)
        if ok:
            self.find_edit.selectAll()

    def show_message(self, text):
        try:
            self.msg.setText(str(text))
            self.msg.setToolTip(str(text))
            self._msg_timer.start()
        except RuntimeError:                         # pragma: no cover
            pass

    def _clear_message(self):
        try:
            self.msg.setText("")
        except RuntimeError:                         # pragma: no cover
            pass

    def sync(self):
        """Show the controller's state on the toggle buttons."""
        ctl = self.controller
        try:
            for a, on in ((self.ortho_act, ctl.ortho),
                          (self.axes_act, ctl.axes_on),
                          (self.legend_act, ctl.legend_on),
                          (self.member_act, ctl.member_labels),
                          (self.node_act, ctl.node_labels)):
                if a is None:
                    continue
                a.blockSignals(True)
                a.setChecked(bool(on))
                a.blockSignals(False)
            self.legend_act.setEnabled(ctl.has_legend())
        except RuntimeError:                         # pragma: no cover
            pass

    def show_help(self, *_):
        """Open the controls list (a window that does not block)."""
        try:
            dlg = self._help
            if dlg is None:
                dlg = QtWidgets.QDialog(self.window())
                dlg.setWindowTitle("3D view controls")
                dlg.resize(560, 620)
                v = QtWidgets.QVBoxLayout(dlg)
                browser = QtWidgets.QTextBrowser()
                browser.setHtml(help_html())
                v.addWidget(browser)
                close = QtWidgets.QPushButton("Close")
                close.clicked.connect(dlg.close)
                v.addWidget(close, 0, Qt.AlignmentFlag.AlignRight)
                dlg.browser = browser
                self._help = dlg
            dlg.show()
            dlg.raise_()
            return dlg
        except Exception:
            self.controller._report(traceback.format_exc())
            return None


class InfoPanel(QtWidgets.QFrame):
    """Details of the selected member or node, on the right of the view.

    It is open from the start, so the drawing does not move when the first
    item is picked (the drawing is limited by its height, so the panel
    hardly ever makes the model smaller).  Its x button or the Info button
    of the bar closes it, and the Info button opens it again.
    """

    WIDTH = 252

    def __init__(self, controller, bar=None, parent=None):
        super().__init__(parent)
        self.setObjectName("view3dInfo")
        self.setStyleSheet(_BAR_CSS)
        self.setFixedWidth(self.WIDTH)
        self.controller = controller
        self.bar = bar
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(8, 6, 6, 6)
        v.setSpacing(4)
        head = QtWidgets.QHBoxLayout()
        head.addWidget(st.heading("Selection"))
        head.addStretch(1)
        self.close_btn = QtWidgets.QToolButton()
        self.close_btn.setText("×")
        self.close_btn.setToolTip("Close the details")
        self.close_btn.setAutoRaise(True)
        self.close_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.close_btn.clicked.connect(self._on_close)
        head.addWidget(self.close_btn)
        v.addLayout(head)
        self.text = QtWidgets.QTextBrowser()
        self.text.setOpenLinks(False)
        v.addWidget(self.text, 1)
        row = QtWidgets.QHBoxLayout()
        self.copy_btn = QtWidgets.QPushButton("Copy")
        self.copy_btn.setToolTip("Copy the details (Ctrl + C on the view)")
        self.copy_btn.clicked.connect(self._on_copy)
        self.clear_btn = QtWidgets.QPushButton("Clear")
        self.clear_btn.setToolTip("Clear the selection (Esc)")
        self.clear_btn.clicked.connect(self._on_clear)
        row.addWidget(self.copy_btn)
        row.addWidget(self.clear_btn)
        row.addStretch(1)
        v.addLayout(row)
        controller.info_changed.connect(self.refresh)

    def is_open(self):
        return not self.isHidden()

    def set_open(self, on, by_user=False):
        self.setVisible(bool(on))
        b = self.bar.info_act if self.bar is not None else None
        if b is not None:
            b.blockSignals(True)
            b.setChecked(bool(on))
            b.blockSignals(False)
        if on:
            self.refresh()

    def _on_close(self, *_):
        try:
            self.set_open(False, by_user=True)
        except Exception:
            self.controller._report(traceback.format_exc())

    def _on_copy(self, *_):
        self.controller.copy_selection()

    def _on_clear(self, *_):
        try:
            self.controller.selection.clear(source=self.controller)
        except Exception:
            self.controller._report(traceback.format_exc())

    def refresh(self):
        try:
            # isHidden, not isVisible: the panel of a view tab that is not
            # on screen is never "visible", but it may be open
            if not self.isHidden():
                self.text.setHtml(self.controller.info_html())
        except Exception:
            self.controller._report(traceback.format_exc())


# ============================================================ install =====
def install(widget, selection=None, picking=True, hover=True, bar=True,
            info=True, name=""):
    """Give an :class:`~.canvas.MplCanvas` the interactive 3D viewer.

    Returns the :class:`ViewController`, also kept as ``widget.view3d``.
    ``picking=False`` makes a view-only controller (no selection, hover or
    details), for animated views.
    """
    ctl = ViewController(widget, selection=selection, picking=picking,
                         hover=hover, name=name)
    widget.view3d = ctl
    ctl.bar = None
    ctl.info_panel = None
    if bar:
        b = ViewBar(ctl, picking=picking)
        widget.insert_bar(b)
        ctl.bar = b
    if info and picking:
        p = InfoPanel(ctl, bar=ctl.bar)
        widget.add_side_panel(p)
        ctl.info_panel = p
        if ctl.bar is not None:
            ctl.bar.info_panel = p
        p.set_open(True)
    return ctl
