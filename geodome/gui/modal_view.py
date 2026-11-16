"""Modal analysis in the GUI: options dialog, "Modes" view and "Modes" table.

``install(win)`` adds

* Analysis > "Modal analysis..." and the toolbar button "Modes".  Both open
  :class:`ModalOptionsDialog` (number of modes, dead load as mass, snow mass
  fraction, lumped or consistent mass) and run
  :func:`geodome.modal.modal_analysis` in the background worker.
* The view tab "Modes": a mode selector, a frequency / period label, the mode
  shape drawn on the undeformed structure, an amplitude spin box and an
  Animate toggle.  The animation draws ``x(t) = x0 + a sin(2 pi t / T_a) phi``
  with a QTimer at about 25 frames per second (T_a = 1.6 s, independent of the
  real period, so slow and fast modes are both easy to watch).  The phase is
  taken from the clock, not from a frame count, so the motion keeps its speed
  even if a frame is dropped.  Each frame redraws only the moving struts over
  a cached image of the rest (matplotlib blitting), so a 6V dome animates at
  the full frame rate.  The animation stops when the tab changes, when the
  model changes and when the window closes.
* The result tab "Modes": mode, f [Hz], T [s], effective modal mass in x, y, z
  [%] and the cumulative sums, with a summary line and a CSV export.

The mode shape is the mass-normalised eigenvector phi_n of
``K phi = omega^2 M phi`` (see :mod:`geodome.modal`), scaled so that the
largest nodal displacement equals the chosen fraction of the structure size.
The extension keeps the last result as ``.modal_result`` (a
:class:`geodome.modal.ModalResult`), the analysed model copy as
``.modal_model`` and the options as ``.modal_options``.
"""

from __future__ import annotations

import copy
import os
import time

import numpy as np
from qtpy import QtCore, QtWidgets
from qtpy.QtCore import QEvent, Qt, QTimer

from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from .. import modal as md
from .. import plotting as plot
from . import ext_util as eu
from . import style as st
from .canvas import Canvas3D
from .worker import Worker

FPS = 25                    # animation frames per second
CYCLE_S = 1.6               # duration of one drawn oscillation [s]
DEFAULT_OPTIONS = {"n_modes": 12, "include_dead": True, "snow_fraction": 0.0,
                   "mass_type": md.LUMPED}

RESULT_COLUMNS = [
    ("mode", "Mode", 0),
    ("f", "f [Hz]", 4),
    ("T", "T [s]", 5),
    ("mx", "Meff x [%]", 2),
    ("my", "Meff y [%]", 2),
    ("mz", "Meff z [%]", 2),
    ("cx", "Cum. x [%]", 2),
    ("cy", "Cum. y [%]", 2),
    ("cz", "Cum. z [%]", 2),
    ("dir", "Direction", None),
]


# ================================================================ job ======
def run_modal_job(model, progress=None, **options):
    """Worker job: modal analysis of ``model`` (a private copy).

    Returns ``(model, ModalResult, options)``.
    """
    if progress is not None:
        progress(0.05, "Modal analysis: building the mass and stiffness "
                       "matrices...")
    res = md.modal_analysis(model, **options)
    if progress is not None:
        progress(1.0, "Modal analysis done.")
    return model, res, dict(options)


# ============================================================== drawing ====
def unit_mode(result, k):
    """(n, 3) mode ``k`` scaled so that the largest nodal displacement is 1."""
    phi = result.mode(k)
    pmax = float(np.max(np.linalg.norm(phi, axis=1))) if len(phi) else 0.0
    return phi / (pmax if pmax > 0.0 else 1.0)


def mode_title(result, k):
    f = float(result.frequencies[k])
    T = float(result.periods[k])
    if not np.isfinite(result.omega[k]):
        return f"Mode {k + 1}   unstable (omega^2 < 0)"
    if f <= 0.0:
        return f"Mode {k + 1}   f = 0 Hz (mechanism)"
    return f"Mode {k + 1}   f = {f:.3f} Hz   T = {T:.4f} s"


def draw_mode_shape(ax, model, result, k=0, amplitude=0.10, phase=1.0,
                    note=""):
    """Draw mode ``k`` of ``result`` on ``model`` into a 3D axes.

    The largest nodal displacement is ``amplitude * geom.length_scale``,
    multiplied by ``phase`` (-1..1, used by the animation).  The axes limits
    cover both extreme positions, so they stay fixed while animating.
    Returns ``(line_collection, base_segments, displacement_segments)``;
    a frame of the animation is ``base + s * displacement``.
    """
    geom = model.geom
    P = np.asarray(geom.nodes, dtype=float)
    M = np.asarray(geom.members, dtype=int)
    shape = unit_mode(result, k)
    d = float(amplitude) * float(geom.length_scale)

    plot.prepare_axes(ax, mode_title(result, k))
    base = np.stack([P[M[:, 0]], P[M[:, 1]]], axis=1)
    disp = np.stack([shape[M[:, 0]], shape[M[:, 1]]], axis=1) * d
    ax.add_collection3d(Line3DCollection(base, colors="#ccd3da",
                                         linewidths=0.7))
    mag = np.linalg.norm(shape, axis=1)
    vals = 0.5 * (mag[M[:, 0]] + mag[M[:, 1]])
    lc = Line3DCollection(base + phase * disp, cmap="viridis",
                          linewidths=1.9,
                          norm=Normalize(0.0, max(float(vals.max()), 1e-12)))
    lc.set_array(vals)
    ax.add_collection3d(lc)

    try:
        sup = np.asarray(model.support_nodes(), dtype=int)
    except Exception:
        sup = np.asarray(geom.base_nodes, dtype=int)
    if len(sup):
        b = P[sup]
        ax.scatter(b[:, 0], b[:, 1], b[:, 2], marker="^", s=28, c="#111111",
                   depthshade=False, zorder=6)
    if not getattr(geom, "is_custom", False):
        R = geom.base_diameter / 2.0
        t = np.linspace(0.0, 2.0 * np.pi, 180)
        ax.plot(R * np.cos(t), R * np.sin(t), np.zeros_like(t), color="0.55",
                lw=0.8, ls="--", zorder=0)
    plot.set_equal_aspect(ax, np.vstack([P, P + d * shape, P - d * shape]))

    r = result.effective_mass_ratio[k] * 100.0
    text = (f"effective mass  x {r[0]:.1f} %, y {r[1]:.1f} %, "
            f"z {r[2]:.1f} %\nlargest displacement drawn = "
            f"{amplitude * 100:.0f} % of the size")
    if note:
        text = note + "\n" + text
    ax.text2D(0.02, 0.96, text, transform=ax.transAxes, fontsize=8,
              va="top", bbox=dict(boxstyle="round,pad=0.3", fc="white",
                                  ec="0.7", alpha=0.85))
    return lc, base, disp


def options_note(options):
    """One line describing the mass options, e.g. 'lumped mass, dead load'."""
    o = dict(DEFAULT_OPTIONS, **(options or {}))
    parts = [f"{o['mass_type']} mass", "struts"]
    if o["include_dead"]:
        parts.append("dead load")
    if o["snow_fraction"] > 0.0:
        parts.append(f"{o['snow_fraction']:g} x snow")
    return parts[0] + ": " + " + ".join(parts[1:])


# =============================================================== dialog ====
class ModalOptionsDialog(QtWidgets.QDialog):
    """Options of one modal analysis run."""

    def __init__(self, parent=None, options=None):
        super().__init__(parent)
        o = dict(DEFAULT_OPTIONS, **(options or {}))
        self.setWindowTitle("Modal analysis")
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(st.heading("Natural frequencies and mode shapes"))
        lay.addWidget(st.hint("Solves K phi = omega^2 M phi for the lowest "
                              "modes. The strut mass is always included."))
        form = QtWidgets.QFormLayout()
        lay.addLayout(form)

        self.n_modes = st.spin(1, 500, o["n_modes"])
        self.n_modes.setToolTip("How many of the lowest modes to compute.")
        form.addRow("Number of modes", self.n_modes)

        self.include_dead = QtWidgets.QCheckBox("Include dead load mass")
        self.include_dead.setChecked(bool(o["include_dead"]))
        self.include_dead.setToolTip(
            "Adds the cladding dead load and the downward D nodal loads as "
            "mass (load / g).")
        form.addRow("", self.include_dead)

        self.snow = st.dspin(0.0, 1.0, o["snow_fraction"], 0.05, 2)
        self.snow.setToolTip("Part of the balanced snow load taken as mass. "
                             "0 = no snow mass, 1 = the full snow load.")
        form.addRow("Snow mass fraction", self.snow)

        self.mass_type = QtWidgets.QComboBox()
        self.mass_type.addItem("Lumped (half of each strut at each end)",
                               md.LUMPED)
        self.mass_type.addItem("Consistent (bar mass matrix)", md.CONSISTENT)
        i = self.mass_type.findData(o["mass_type"])
        self.mass_type.setCurrentIndex(max(i, 0))
        self.mass_type.setToolTip("Lumped mass is the usual choice. Consistent "
                                  "mass gives slightly higher frequencies.")
        form.addRow("Strut mass", self.mass_type)

        bb = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.button(QtWidgets.QDialogButtonBox.Ok).setText("Run")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)

    def options(self):
        return {"n_modes": int(self.n_modes.value()),
                "include_dead": bool(self.include_dead.isChecked()),
                "snow_fraction": float(self.snow.value()),
                "mass_type": str(self.mass_type.currentData())}


# ============================================================ view tab =====
class ModeCanvas(Canvas3D):
    """The "Modes" view: a 3D canvas that reports when it is hidden."""

    hidden = QtCore.Signal()

    def hideEvent(self, event):                      # noqa: N802 (Qt name)
        self.hidden.emit()
        super().hideEvent(event)


# =========================================================== extension =====
class ModalExtension(QtCore.QObject):
    """Modal analysis menu entry, toolbar button, view tab and result tab."""

    def __init__(self, win):
        super().__init__(win)
        self.win = win
        self.modal_result = None
        self.modal_model = None
        self.modal_options = dict(DEFAULT_OPTIONS)
        self.frames = 0                      # animation frames drawn
        self.phase = 1.0                     # s of the drawn frame (-1..1)
        self.segments = None                 # (m, 2, 3) drawn strut ends
        self._generation = 0                 # bumped on every model change
        self._lc = None
        self._base = None
        self._disp = None
        self._t0 = 0.0
        self._bg = None                      # cached static scene (blitting)
        self._blit = True                    # False after a blit failure
        self._build_view()
        self._build_result_tab()
        self.view.canvas.mpl_connect("draw_event", self._on_draw_event)
        self.timer = QTimer(self.view)
        # a coarse timer on Windows rounds 40 ms up to ~47 ms (21 fps)
        self.timer.setTimerType(Qt.PreciseTimer)
        self.timer.setInterval(int(round(1000.0 / FPS)))
        self.timer.timeout.connect(self._advance)
        self._show_placeholder()

    # ------------------------------------------------------------- build --
    def _build_view(self):
        c = self.view = ModeCanvas()
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.setMinimumWidth(190)
        self.mode_combo.setToolTip("Choose the mode to show.")
        self.mode_combo.currentIndexChanged.connect(self._on_mode_selected)
        self.freq_label = QtWidgets.QLabel("")
        self.freq_label.setObjectName("metric")
        self.amplitude = st.dspin(0.01, 0.50, 0.10, 0.01, 2)
        self.amplitude.setPrefix("amplitude ")
        self.amplitude.setToolTip("Largest drawn displacement as a fraction "
                                  "of the structure size.")
        self.amplitude.valueChanged.connect(self._on_amplitude)
        self.animate_btn = QtWidgets.QPushButton("Animate")
        self.animate_btn.setCheckable(True)
        self.animate_btn.setToolTip("Play the mode shape as a motion.")
        self.animate_btn.toggled.connect(self._on_animate_toggled)
        self.run_btn = QtWidgets.QPushButton("Run modal analysis...")
        self.run_btn.setToolTip("Compute the natural frequencies and mode "
                                "shapes.")
        self.run_btn.clicked.connect(self.on_modal_action)
        c.add_control(QtWidgets.QLabel("Mode:"))
        c.add_control(self.mode_combo)
        c.add_control(self.freq_label)
        c.add_control_stretch()
        for w in (self.amplitude, self.animate_btn, self.run_btn):
            c.add_control(w)
        c.hidden.connect(self.stop_animation)

    def _build_result_tab(self):
        w = self.result_widget = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)
        row = QtWidgets.QHBoxLayout()
        self.summary_label = QtWidgets.QLabel("No modal analysis yet.")
        self.summary_label.setObjectName("metric")
        self.summary_label.setWordWrap(True)
        row.addWidget(self.summary_label, 1)
        self.export_btn = QtWidgets.QPushButton("Export CSV...")
        self.export_btn.setToolTip("Save this table as a CSV file.")
        self.export_btn.clicked.connect(lambda: self.export_csv())
        self.export_btn.setEnabled(False)
        row.addWidget(self.export_btn)
        lay.addLayout(row)
        self.table = st.table([c[1] for c in RESULT_COLUMNS])
        self.table.setToolTip("Double-click a mode to show it in the "
                              "Modes view.")
        self.table.cellDoubleClicked.connect(self._on_table_double_click)
        lay.addWidget(self.table, 1)

    def install(self):
        win = self.win
        self.act_modal = win._act("&Modal analysis...", self.on_modal_action,
                                  tip="Natural frequencies and mode shapes")
        menu = (getattr(win, "menu_analysis", None)
                or eu.ensure_menu(win, "Analysis"))
        eu.insert_in_first_group(menu, self.act_modal)
        self.act_modes = win._act("Modes", self.on_modal_action,
                                  tip="Modal analysis: natural frequencies "
                                      "and mode shapes")
        eu.insert_toolbar_action(win, self.act_modes)

        # appended (not inserted) so the index of every existing tab stays
        # the same; the View menu lists the view tabs when it opens
        win.view_tabs.addTab(self.view, "Modes")
        eu.add_result_tab(win, self.result_widget, "Modes")

        win.view_tabs.currentChanged.connect(self._on_view_tab_changed)
        win.installEventFilter(self)
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self.stop_animation)
        eu.register(win, self)
        return self

    # -------------------------------------------------------- extension API --
    def on_model_changed(self, model):
        """A new model: the old modes no longer apply."""
        self._generation += 1
        self.clear()

    def on_results(self, results):
        """Static results do not change the modes."""
        return None

    def clear(self):
        self.stop_animation()
        self.modal_result = None
        self.modal_model = None
        self.mode_combo.blockSignals(True)
        self.mode_combo.clear()
        self.mode_combo.blockSignals(False)
        self.freq_label.setText("")
        self.table.clearContents()
        self.table.setRowCount(0)
        self.summary_label.setText("No modal analysis yet.")
        self.export_btn.setEnabled(False)
        self._show_placeholder()

    # --------------------------------------------------------------- run --
    def on_modal_action(self):
        win = self.win
        if win._busy:
            win.status("Busy: wait for the current run to finish.")
            return
        dlg = ModalOptionsDialog(win, self.modal_options)
        ok = dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()
        if ok:
            self.run_modal(dlg.options())

    def run_modal(self, options=None):
        """Start a modal analysis of the current model in the worker.

        Returns True when the run was started.
        """
        win = self.win
        if win._busy:
            win.status("Busy: wait for the current run to finish.")
            return False
        opts = dict(DEFAULT_OPTIONS, **(options or {}))
        try:
            # returns False (after showing a message) on a bad input
            if win._read_inputs() is False:
                return False
            snapshot = copy.deepcopy(win.model)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Modal analysis",
                                          f"Check the inputs: {exc}")
            return False
        self.modal_options = opts
        self.stop_animation()
        self._job_generation = self._generation
        self.act_modal.setEnabled(False)
        self.act_modes.setEnabled(False)
        self.run_btn.setEnabled(False)
        win._set_busy(True, "Running modal analysis...")
        w = Worker(run_modal_job, snapshot, **opts)
        w.progressed.connect(win._on_progress)
        w.failed.connect(self._job_ended)
        w.failed.connect(win._on_failed)
        w.succeeded.connect(self._modal_done)
        win._worker = w
        w.start()
        return True

    def _job_ended(self, *args):
        for a in (self.act_modal, self.act_modes, self.run_btn):
            a.setEnabled(True)

    def _modal_done(self, payload):
        win = self.win
        win._set_busy(False)
        self._job_ended()
        model, res, opts = payload
        if getattr(self, "_job_generation", None) != self._generation:
            win.log("Modal result discarded: the model changed during the "
                    "run.")
            win.status("Modal result discarded (model changed).")
            return
        if res.n_modes == 0:
            win.log("Modal analysis: no modes. The model has no free degrees "
                    "of freedom with mass.")
            QtWidgets.QMessageBox.information(
                win, "Modal analysis", "No vibration modes: the model has no "
                "free degrees of freedom with mass.")
            return
        self.modal_model = model
        self.modal_result = res
        self.modal_options = opts
        self._fill_result_table()
        self._populate_selector()
        self.draw_mode(0)
        win.view_tabs.setCurrentWidget(self.view)
        win.result_tabs.setCurrentWidget(self.result_widget)
        win.log(res.summary_text())
        f1 = res.fundamental_frequency
        msg = (f"Modal analysis done: {res.n_modes} modes, f1 = {f1:.3f} Hz "
               f"(T1 = {1.0 / f1:.4f} s)." if f1 > 0 else
               f"Modal analysis done: {res.n_modes} modes.")
        if res.warnings:
            msg += f" {len(res.warnings)} warning(s), see the Log."
        win.status(msg)

    # ------------------------------------------------------------- table --
    def table_rows(self):
        res = self.modal_result
        if res is None:
            return []
        rows = []
        for r in res.table():
            k = r["mode"] - 1
            unstable = not np.isfinite(res.omega[k])
            T = r["period_s"]
            rows.append({
                "mode": r["mode"],
                "f": float("nan") if unstable else r["frequency_Hz"],
                "T": T if np.isfinite(T) else float("nan"),
                "mx": 100.0 * r["ratio_x"], "my": 100.0 * r["ratio_y"],
                "mz": 100.0 * r["ratio_z"],
                "cx": 100.0 * r["cum_x"], "cy": 100.0 * r["cum_y"],
                "cz": 100.0 * r["cum_z"],
                "dir": "unstable" if unstable else r["dominant"],
            })
        return rows

    def _fill_result_table(self):
        res = self.modal_result
        eu.fill(self.table, RESULT_COLUMNS, self.table_rows())
        pr = res.participating_ratio * 100.0
        need = res.modes_for_ratio(0.90)
        more = f"> {res.n_modes}"
        tm = res.total_mass
        self.summary_label.setText(
            f"{options_note(self.modal_options)}.  Total mass "
            f"{tm[0]:,.1f} kg.  Captured by {res.n_modes} modes: "
            f"x {pr[0]:.1f} %, y {pr[1]:.1f} %, z {pr[2]:.1f} %.  "
            "Modes for 90 %: " + ", ".join(
                f"{a} {need[a] if need[a] is not None else more}"
                for a in ("x", "y", "z")) + ".")
        self.export_btn.setEnabled(True)

    def _on_table_double_click(self, row, _col):
        if self.modal_result is None or row >= self.mode_combo.count():
            return
        self.mode_combo.setCurrentIndex(row)
        self.win.view_tabs.setCurrentWidget(self.view)

    def export_csv(self, path=None):
        """Write the mode table to CSV (asks for a file when ``path`` is
        None).  Returns the path or None."""
        if self.modal_result is None:
            return None
        if path is None:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self.win, "Export modes", os.path.join(os.getcwd(),
                                                       "modes.csv"),
                "CSV (*.csv)")
            if not path:
                return None
        try:
            md.export_modal_csv(self.modal_result, path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self.win, "Export modes", str(exc))
            return None
        self.win.log(f"Mode table written to {path}")
        self.win.status("Mode table written.")
        return path

    # ------------------------------------------------------------- views --
    def _populate_selector(self):
        res = self.modal_result
        self.mode_combo.blockSignals(True)
        self.mode_combo.clear()
        for k in range(res.n_modes):
            if not np.isfinite(res.omega[k]):
                text = f"{k + 1}:  unstable"
            else:
                text = (f"{k + 1}:  {res.frequencies[k]:.3f} Hz  "
                        f"({res.dominant_direction(k)})")
            self.mode_combo.addItem(text)
        self.mode_combo.setCurrentIndex(0)
        self.mode_combo.blockSignals(False)

    def selected_mode(self):
        return max(self.mode_combo.currentIndex(), 0)

    def _show_placeholder(self):
        ax = self.view.clear()
        plot.prepare_axes(ax, "Mode shapes")
        ax.text2D(0.5, 0.5, "Run Analysis > Modal analysis (toolbar: Modes)\n"
                            "to see the natural frequencies and mode shapes.",
                  transform=ax.transAxes, ha="center", va="center",
                  fontsize=10, color="#43566a")
        self._lc = None
        self.segments = None
        self.view.refresh()

    def _update_label(self, k):
        res = self.modal_result
        if not np.isfinite(res.omega[k]):
            self.freq_label.setText("unstable mode (omega^2 < 0)")
            return
        f = float(res.frequencies[k])
        T = float(res.periods[k])
        Ts = f"{T:.5f} s" if np.isfinite(T) else "-"
        self.freq_label.setText(
            f"f = {f:.4f} Hz   T = {Ts}   omega = {res.omega[k]:.2f} rad/s")

    def draw_mode(self, k=None):
        """Draw mode ``k`` (0-based; default: the selected one)."""
        res = self.modal_result
        if res is None:
            self._show_placeholder()
            return
        k = self.selected_mode() if k is None else int(k)
        k = int(np.clip(k, 0, res.n_modes - 1))
        ax_old = self.view.axes
        view = None
        if self._lc is not None:
            view = (ax_old.elev, ax_old.azim)
        ax = self.view.clear()
        self.phase = self._phase() if self.is_animating() else 1.0
        self._lc, self._base, self._disp = draw_mode_shape(
            ax, self.modal_model, res, k,
            amplitude=self.amplitude.value(), phase=self.phase,
            note=options_note(self.modal_options))
        self.segments = self._base + self.phase * self._disp
        self._bg = None                      # new scene: recapture on draw
        if self.is_animating():
            self._lc.set_animated(True)
        if view is not None:
            ax.view_init(elev=view[0], azim=view[1])
        self.view.add_colorbar(self._lc, "normalised mode amplitude")
        self._update_label(k)
        self.view.refresh()

    def _on_mode_selected(self, _index):
        if self.modal_result is not None:
            self.draw_mode()

    def _on_amplitude(self, _value):
        if self.modal_result is not None:
            self.draw_mode()

    def _on_view_tab_changed(self, _index):
        if self.win.view_tabs.currentWidget() is not self.view:
            self.stop_animation()

    # --------------------------------------------------------- animation --
    def is_animating(self):
        return self.timer.isActive()

    def _phase(self):
        return float(np.sin(2.0 * np.pi * (time.perf_counter() - self._t0)
                            / CYCLE_S))

    def start_animation(self):
        """Start the animation (needs a result). Returns True if running."""
        if self.modal_result is None or self._lc is None:
            self.animate_btn.blockSignals(True)
            self.animate_btn.setChecked(False)
            self.animate_btn.blockSignals(False)
            return False
        if not self.timer.isActive():
            self._t0 = time.perf_counter()
            self._bg = None
            self.timer.start()
            if self._blit:
                self._lc.set_animated(True)
                self.view.canvas.draw_idle()     # caches the static scene
        self.animate_btn.blockSignals(True)
        self.animate_btn.setChecked(True)
        self.animate_btn.setText("Stop")
        self.animate_btn.blockSignals(False)
        return True

    def stop_animation(self):
        """Stop the animation and show the shape at full amplitude."""
        was = self.timer.isActive()
        self.timer.stop()
        self.animate_btn.blockSignals(True)
        self.animate_btn.setChecked(False)
        self.animate_btn.setText("Animate")
        self.animate_btn.blockSignals(False)
        self._bg = None
        if self._lc is not None:
            try:
                self._lc.set_animated(False)
            except Exception:                        # window closing
                pass
            if was:
                self._set_frame(1.0)

    def _set_frame(self, s):
        """Draw the frame x0 + s * a * phi (s in -1..1).

        While animating, only the moving struts are redrawn over the cached
        static scene (matplotlib blitting).  Otherwise, or if blitting fails,
        the canvas is redrawn in full.
        """
        self.phase = float(s)
        self.segments = self._base + self.phase * self._disp
        try:
            self._lc.set_segments(self.segments)
        except Exception:                            # window closing
            return
        if self.is_animating() and self._blit and self._bg is not None:
            try:
                c = self.view.canvas
                c.restore_region(self._bg)
                self._draw_animated()
                c.blit(self.view.figure.bbox)
                return
            except Exception:
                self._stop_blitting()
        try:
            self.view.canvas.draw_idle()
        except Exception:                            # window closing
            pass

    def _draw_animated(self):
        """Project and draw the moving struts on the current canvas."""
        lc = self._lc
        try:
            lc.do_3d_projection()
        except TypeError:                            # matplotlib < 3.5
            lc.do_3d_projection(self.view.canvas.get_renderer())
        self.view.axes.draw_artist(lc)

    def _on_draw_event(self, _event):
        """After a full draw (start, rotation, resize, new mode) cache the
        static scene and put the moving struts on top."""
        if not (self.is_animating() and self._blit and self._lc is not None
                and self._lc.get_animated()):
            return
        try:
            self._bg = self.view.canvas.copy_from_bbox(self.view.figure.bbox)
            self._draw_animated()
        except Exception:
            self._stop_blitting()

    def _stop_blitting(self):
        """Fall back to full redraws for the rest of the session."""
        self._blit = False
        self._bg = None
        try:
            self._lc.set_animated(False)
        except Exception:
            pass

    def _on_animate_toggled(self, checked):
        if checked:
            self.start_animation()
        else:
            self.stop_animation()

    def _advance(self):
        """One animation frame (QTimer slot)."""
        if self._lc is None or self.modal_result is None:
            self.stop_animation()
            return
        self._set_frame(self._phase())
        self.frames += 1

    # ------------------------------------------------------------ events --
    def eventFilter(self, obj, event):               # noqa: N802 (Qt name)
        # Never raise here: Qt still calls this while the window is torn
        # down, when the Python attributes may already be gone, and PyQt6
        # aborts the process on an exception in an event filter.
        try:
            if obj is self.__dict__.get("win") and \
                    event.type() == QEvent.Close:
                self.stop_animation()
        except Exception:                                # pragma: no cover
            pass
        return False


def install(win):
    """Add modal analysis to ``win`` and return the extension object.

    The extension registers itself in ``win.extensions``.  Calling install
    again on the same window returns the extension already there.
    """
    for ext in getattr(win, "extensions", None) or []:
        if isinstance(ext, ModalExtension):
            return ext
    return ModalExtension(win).install()
