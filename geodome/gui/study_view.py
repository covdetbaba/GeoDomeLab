"""Parameter study in the GUI (Analysis > Parameter study...).

The study itself is :func:`geodome.study.parameter_study`: it regenerates
the geodesic dome for every combination of frequency nu, sphere radius R and
dome fraction, and runs the full GeoDomeLab analysis on each one (direct
stiffness solve, AISC 360 member checks, deflection check and, if asked,
linear buckling, optionally after fully stressed auto-sizing).  This module
only adds the front end:

* a dialog for the frequencies (``2, 4, 6, 8``), the radius values and the
  dome fraction values (one value keeps a parameter fixed, several sweep it)
  and the auto-size and buckling options;
* a run in :class:`gui.worker.Worker` with progress and a Cancel button
  (``parameter_study(cancel=...)`` is polled before every study point, so
  Cancel stops after the point that is running);
* a view tab "Study" with four charts against the swept parameter:
  steel mass  M = sum rho A L,  max utilisation  max |N| / (phi Pn),
  buckling factor lambda_cr  (det(K + lambda Kg) = 0)  and
  span / deflection  D_base / max|u|  (with the limit as a dashed line);
* a result tab "Study" with the table, CSV export
  (:func:`geodome.study.study_to_csv`) and "Open selected row as model"
  (:func:`geodome.study.model_from_row`).

``install(win)`` returns the extension object.  It keeps the last study in
``.study_rows`` (a :class:`geodome.study.StudyRows`, or None) and the model
the study started from in ``.study_base_model``.
"""

from __future__ import annotations

import copy
import math
import re

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

from .. import geometry as geo
from .. import study as sy
from . import style as st
from .canvas import Canvas2D
from .worker import Worker

# Colours of the study lines (one per series), starting with the app accent.
SERIES_COLOURS = ("#1f5f8b", "#b45309", "#15803d", "#7c3aed", "#b4161b",
                  "#0e7490", "#4b5563", "#a16207")

# The four charts of the Study tab: (row key, panel title).
CHARTS = (("mass_kg", "Steel mass"),
          ("max_util", "Max utilisation"),
          ("lambda_cr", "Buckling factor lambda_cr"),
          ("span_over_deflection", "Span / deflection"))

# x-axis choices: (label, row key).
X_CHOICES = (("Frequency", "frequency"),
             ("Radius", "radius_m"),
             ("Dome fraction", "fraction"),
             ("Span", "span_m"),
             ("Nodes", "nodes"),
             ("Struts", "struts"))

# Result table: (row key, header).
TABLE_COLUMNS = (("frequency", "Freq [V]"),
                 ("radius_m", "R [m]"),
                 ("fraction", "Fraction"),
                 ("span_m", "Span [m]"),
                 ("nodes", "Nodes"),
                 ("struts", "Struts"),
                 ("families", "Families"),
                 ("sections", "Sections"),
                 ("mass_kg", "Mass [kg]"),
                 ("mass_per_plan_kg_m2", "Mass/plan [kg/m2]"),
                 ("max_util", "Max util"),
                 ("governing_combination", "Governing"),
                 ("lambda_cr", "lambda_cr"),
                 ("span_over_deflection", "Span/defl"),
                 ("status", "Status"),
                 ("time_s", "Time [s]"),
                 ("notes", "Notes"))

FREQ_RANGE = (1, 10)            # same limits as the Geometry tab
RADIUS_RANGE = (1.0, 200.0)     # [m]
MANY_POINTS = 30                # ask before running a study this big

_INT_KEYS = ("frequency", "nodes", "struts", "panels", "families", "supports",
             "members_over", "governing_member", "sizing_iterations")


# ================================================================ parsing ==
def parse_values(text, name, kind=float):
    """A list of numbers from ``"2, 4, 6"`` (commas, semicolons or spaces).

    Raises ``ValueError`` with a short message for bad text.
    """
    parts = [p for p in re.split(r"[,;\s]+", str(text).strip()) if p]
    if not parts:
        raise ValueError(f"{name}: give at least one value.")
    out = []
    for p in parts:
        try:
            v = float(p)
        except ValueError:
            raise ValueError(f"{name}: '{p}' is not a number.") from None
        if not math.isfinite(v):
            raise ValueError(f"{name}: '{p}' is not a finite number.")
        if kind is int:
            if v != round(v):
                raise ValueError(f"{name}: '{p}' is not a whole number.")
            v = int(round(v))
        out.append(v)
    # keep the order, drop repeats
    return list(dict.fromkeys(out))


def check_study_inputs(freq_text, radius_text, fraction_text):
    """Parse and check the three dialog fields.

    Returns ``(frequencies, radii, fractions)``; raises ``ValueError``.
    """
    freqs = parse_values(freq_text, "Frequencies", int)
    lo, hi = FREQ_RANGE
    for f in freqs:
        if not lo <= f <= hi:
            raise ValueError(f"Frequencies: {f} is outside {lo} to {hi}.")
    radii = parse_values(radius_text, "Radius", float)
    rlo, rhi = RADIUS_RANGE
    for r in radii:
        if not rlo <= r <= rhi:
            raise ValueError(f"Radius: {r:g} m is outside {rlo:g} to "
                             f"{rhi:g} m.")
    fractions = parse_values(fraction_text, "Dome fraction", float)
    for f in fractions:
        if not 0.0 < f <= 1.0:
            raise ValueError(f"Dome fraction: {f:g} must be above 0 and at "
                             f"most 1 (0.5 is a hemisphere).")
    return freqs, radii, fractions


def default_x(rows_or_settings):
    """The natural x axis: the first parameter that takes several values."""
    if isinstance(rows_or_settings, dict):
        s = rows_or_settings
        fr, ra, fc = (s.get("frequencies", []), s.get("radii", []),
                      s.get("fractions", []))
    else:
        rows = list(rows_or_settings)
        fr = {r.get("frequency") for r in rows}
        ra = {r.get("radius_m") for r in rows}
        fc = {r.get("fraction_requested") for r in rows}
    if len(fr) > 1:
        return "frequency"
    if len(ra) > 1:
        return "radius_m"
    if len(fc) > 1:
        return "fraction"
    return "frequency"


# =============================================================== format ==
def format_value(key, v):
    """Table text of one study value (the same rules as study.study_table)."""
    if v is None:
        return "-"
    if isinstance(v, str):
        return v
    if isinstance(v, bool):
        return "yes" if v else "no"
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(v):
        return "n/a" if math.isnan(v) else "inf"
    if key in _INT_KEYS:
        return f"{int(v):d}"
    if key in ("fraction", "fraction_requested", "max_util"):
        return f"{v:.3f}"
    if key in ("mass_kg", "span_over_deflection", "deflection_limit"):
        return f"{v:,.0f}"
    if key in ("mass_per_plan_kg_m2", "lambda_cr", "time_s", "radius_m",
               "span_m"):
        return f"{v:,.2f}"
    return f"{v:,.3g}"


def describe_row(row):
    """Short name of a study point, e.g. '6V, R = 15 m, fraction 0.50'."""
    return (f"{int(row['frequency'])}V, R = {float(row['radius_m']):g} m, "
            f"fraction {float(row.get('fraction', float('nan'))):.2f}")


def summary_line(rows):
    """One line under the table: counts and the lightest acceptable design."""
    if rows is None:
        return "No study yet. Use Analysis > Parameter study..."
    n = len(rows)
    if n == 0:
        return "The study has no points (it was stopped before the first one)."
    n_ok = sum(1 for r in rows if r.get("status") == "OK")
    n_err = sum(1 for r in rows if r.get("status") == "ERROR")
    text = f"{n} dome{'s' if n != 1 else ''} analysed, {n_ok} pass every check"
    if n_err:
        text += f", {n_err} failed to run"
    best = sy.best_design(rows)
    if best is not None:
        text += (f". Lightest acceptable design: {describe_row(best)}, "
                 f"{best['mass_kg']:,.0f} kg.")
    else:
        text += ". No design passes every check."
    return text


# ================================================================ charts ==
def _series_keys(rows, x):
    keys = []
    for k in sy.SWEEP_KEYS:
        if k == x or (x == "fraction" and k == "fraction_requested"):
            continue
        if len({r.get(k) for r in rows}) > 1:
            keys.append(k)
    return keys


def _series_label(keys, values):
    parts = []
    for k, v in zip(keys, values):
        if k == "frequency":
            parts.append(f"{int(v)}V")
        elif k == "radius_m":
            parts.append(f"R = {v:g} m")
        else:
            parts.append(f"fraction {v:g}")
    return ", ".join(parts)


def _finite(v):
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def draw_study(fig, rows, x="frequency"):
    """Draw the four study charts on ``fig`` (2 x 2).  Returns the axes.

    Filled markers: the design passes every check; hollow: it fails one.
    Dashed lines: utilisation 1.0, lambda_cr 1.0, the span/deflection limit.
    """
    fig.clear()
    axes = list(fig.subplots(2, 2).ravel())
    good = [r for r in (rows or []) if r.get("status") != "ERROR"]
    xlabel = sy.COLUMN_LABELS.get(x, x)
    if not good:
        msg = ("No study yet.\nUse Analysis > Parameter study..."
               if rows is None else "No study results to plot.")
        for ax, (_, title) in zip(axes, CHARTS):
            ax.set_title(title, fontsize=10, color=st.ACCENT_DARK)
            ax.set_xticks([])
            ax.set_yticks([])
        axes[0].text(0.5, 0.5, msg, transform=axes[0].transAxes,
                     ha="center", va="center", fontsize=10, color="#64748b")
        _tight(fig)
        return axes

    skeys = _series_keys(good, x)
    groups = {}
    for r in good:
        groups.setdefault(tuple(r.get(k) for k in skeys), []).append(r)
    order = sorted(groups)

    for ax, (yk, title) in zip(axes, CHARTS):
        plotted = False
        for gi, gk in enumerate(order):
            grp = [r for r in groups[gk] if _finite(r.get(x))]
            grp.sort(key=lambda r: float(r[x]))
            pts = [(float(r[x]), float(r[yk]), r.get("status") == "OK")
                   for r in grp if _finite(r.get(yk))]
            if not pts:
                continue
            col = SERIES_COLOURS[gi % len(SERIES_COLOURS)]
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            label = _series_label(skeys, gk) if skeys else None
            ax.plot(xs, ys, "-", lw=1.6, color=col, label=label)
            ok = [p for p in pts if p[2]]
            bad = [p for p in pts if not p[2]]
            if ok:
                ax.plot([p[0] for p in ok], [p[1] for p in ok], "o",
                        color=col, ms=6)
            if bad:
                ax.plot([p[0] for p in bad], [p[1] for p in bad], "o",
                        mfc="white", mec=col, mew=1.5, ms=7)
            plotted = True
        if yk in ("max_util", "lambda_cr"):
            ax.axhline(1.0, ls="--", lw=1.0, color="#94a3b8")
        if yk == "span_over_deflection":
            lims = sorted({float(r["deflection_limit"]) for r in good
                           if _finite(r.get("deflection_limit"))})
            for lim in lims:
                ax.axhline(lim, ls="--", lw=1.0, color="#94a3b8",
                           label=f"limit span/{lim:g}")
            if lims:
                handles = [h for h in ax.get_lines()
                           if str(h.get_label()).startswith("limit")]
                ax.legend(handles=handles, fontsize=8, frameon=False,
                          loc="best")
        if not plotted:
            note = ("lambda_cr was not computed"
                    if yk == "lambda_cr" else "no values")
            ax.text(0.5, 0.5, note, transform=ax.transAxes, ha="center",
                    va="center", fontsize=9, color="#64748b")
        ax.set_title(title, fontsize=10, color=st.ACCENT_DARK)
        ax.set_ylabel(sy.COLUMN_LABELS.get(yk, yk), fontsize=9)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, alpha=0.3)
        if x in _INT_KEYS:
            # ticks at the studied values only (no 3V tick for a 2V/4V study)
            xs_all = sorted({float(r[x]) for r in good if _finite(r.get(x))})
            if 0 < len(xs_all) <= 12:
                ax.set_xticks(xs_all)
                ax.set_xticklabels([f"{v:g}" for v in xs_all])
            else:
                from matplotlib.ticker import MaxNLocator
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    if skeys and len(order) > 1:
        axes[0].legend(fontsize=8, frameon=False)
    sized = " (auto-sized)" if any(r.get("auto_sized") for r in good) else ""
    fig.suptitle(f"Parameter study{sized}: filled marker = design OK, "
                 f"hollow = fails a check", fontsize=10)
    _tight(fig, rect=(0, 0, 1, 0.95))
    return axes


def _tight(fig, **kw):
    try:
        fig.tight_layout(h_pad=1.6, w_pad=2.0, **kw)
    except Exception:                                  # pragma: no cover
        pass


# ================================================================ dialog ==
class StudyDialog(QtWidgets.QDialog):
    """Inputs of a parameter study."""

    def __init__(self, model, parent=None, design_note=""):
        super().__init__(parent)
        self.setWindowTitle("Parameter study")
        self.setMinimumWidth(470)
        g = model.geom
        radius = float(g.radius) if _finite(g.radius) and g.radius > 0 \
            else 12.0
        fraction = (float(g.dome_fraction) if _finite(g.dome_fraction)
                    else 0.5)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setSpacing(8)
        lay.addWidget(st.heading("Compare domes of different size and "
                                 "frequency"))
        lay.addWidget(st.hint(
            "Each dome is rebuilt and fully analysed with the current "
            "material, supports, loads and design settings. Give one value "
            "to keep a parameter fixed, or several values separated by "
            "commas to sweep it."))

        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.frequencies = QtWidgets.QLineEdit("2, 4, 6, 8")
        self.frequencies.setToolTip(
            f"Whole numbers from {FREQ_RANGE[0]} to {FREQ_RANGE[1]}, "
            f"for example 2, 4, 6, 8.")
        form.addRow("Frequencies [V]", self.frequencies)
        self.radii = QtWidgets.QLineEdit(f"{radius:g}")
        self.radii.setToolTip(
            "Sphere radius in metres. One value or several, "
            "for example 8, 12, 16.")
        form.addRow("Radius values [m]", self.radii)
        self.fractions = QtWidgets.QLineEdit(f"{fraction:g}")
        self.fractions.setToolTip(
            "Part of the sphere that is kept: 0.5 is a hemisphere. "
            "The cut snaps to the nearest ring of nodes.")
        form.addRow("Dome fraction values", self.fractions)
        lay.addLayout(form)

        self.auto_size = QtWidgets.QCheckBox(
            "Auto-size the strut families of every dome")
        self.auto_size.setToolTip(
            "Fully stressed design for each dome, with the target "
            "utilisation, iterations and deflection option of the Design "
            "tab. Slower, but compares economical designs.")
        lay.addWidget(self.auto_size)
        self.buckling = QtWidgets.QCheckBox(
            "Compute the buckling factor lambda_cr")
        self.buckling.setChecked(True)
        self.buckling.setToolTip("One linear buckling eigenanalysis per dome.")
        lay.addWidget(self.buckling)

        if design_note:
            lay.addWidget(st.hint(design_note))
        self.points_label = QtWidgets.QLabel("")
        self.points_label.setStyleSheet(f"color: {st.ACCENT_DARK};")
        lay.addWidget(self.points_label)
        self.error_label = QtWidgets.QLabel("")
        self.error_label.setWordWrap(True)
        self.error_label.setStyleSheet(f"color: {st.BAD_RED};")
        lay.addWidget(self.error_label)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.run_button = buttons.button(QtWidgets.QDialogButtonBox.Ok)
        self.run_button.setText("Run study")
        self.run_button.setObjectName("primary")
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

        for w in (self.frequencies, self.radii, self.fractions):
            w.textChanged.connect(self._update_count)
        self._update_count()

    # ------------------------------------------------------------------
    def values(self):
        """``dict(frequencies, radii, fractions, auto_size, do_buckling)``.

        Raises ``ValueError`` with a short message for bad input.
        """
        freqs, radii, fractions = check_study_inputs(
            self.frequencies.text(), self.radii.text(), self.fractions.text())
        return {"frequencies": freqs, "radii": radii, "fractions": fractions,
                "auto_size": self.auto_size.isChecked(),
                "do_buckling": self.buckling.isChecked()}

    def _update_count(self):
        try:
            v = self.values()
        except ValueError as exc:
            self.points_label.setText("")
            self.error_label.setText(str(exc))
            self.run_button.setEnabled(False)
            return
        n = len(v["frequencies"]) * len(v["radii"]) * len(v["fractions"])
        self.points_label.setText(f"{n} dome{'s' if n != 1 else ''} to "
                                  f"analyse.")
        self.error_label.setText("")
        self.run_button.setEnabled(True)

    def _accept(self):
        try:
            self.values()
        except ValueError as exc:
            self.error_label.setText(str(exc))
            return
        self.accept()


# ============================================================= extension ==
class StudyExtension(QtCore.QObject):
    """Menu action, view tab, result tab and worker run of the study."""

    def __init__(self, win):
        super().__init__(win)
        self.win = win
        self.study_rows = None          # StudyRows of the last study
        self.study_base_model = None    # deep copy of the model it started from
        self.x_key = "frequency"
        self._cancel_requested = False
        self._running = False
        self._build_action()
        self._build_view_tab()
        self._build_result_tab()
        self._build_cancel_button()
        win.installEventFilter(self)
        self.redraw()

    # ------------------------------------------------------------ build --
    def _build_action(self):
        win = self.win
        self.act_study = win._act(
            "Parameter &study...", self.on_study,
            tip="Analyse a set of domes with other frequencies, radii or "
                "dome fractions and compare them.")
        menu = getattr(win, "menu_analysis", None) or _find_menu(win,
                                                                 "Analysis")
        if menu is None:
            menu = win.menuBar().addMenu("&Analysis")
        menu.addSeparator()
        menu.addAction(self.act_study)

    def _build_view_tab(self):
        self.canvas = Canvas2D()
        self.x_combo = QtWidgets.QComboBox()
        for label, key in X_CHOICES:
            self.x_combo.addItem(label, key)
        self.x_combo.setToolTip("Quantity on the horizontal axis.")
        self.x_combo.currentIndexChanged.connect(self._on_x_changed)
        self.view_new_button = QtWidgets.QPushButton("New study...")
        self.view_new_button.clicked.connect(self.on_study)
        self.canvas.add_control(QtWidgets.QLabel("x axis:"))
        self.canvas.add_control(self.x_combo)
        self.canvas.add_control(self.view_new_button)
        self.canvas.add_control_stretch()
        self.view_tab = self.canvas
        self.win.view_tabs.addTab(self.canvas, "Study")

    def _build_result_tab(self):
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        self.summary_label = QtWidgets.QLabel(summary_line(None))
        self.summary_label.setWordWrap(True)
        lay.addWidget(self.summary_label)
        self.table = st.table([h for _, h in TABLE_COLUMNS])
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SingleSelection)
        self.table.itemSelectionChanged.connect(self._update_buttons)
        lay.addWidget(self.table, 1)
        row = QtWidgets.QHBoxLayout()
        self.csv_button = QtWidgets.QPushButton("Export CSV...")
        self.csv_button.setToolTip("Write every study column to a CSV file.")
        self.csv_button.clicked.connect(lambda: self.export_csv())
        self.open_button = QtWidgets.QPushButton("Open selected row as model")
        self.open_button.setToolTip(
            "Replace the current model with the selected dome, "
            "with its sections.")
        self.open_button.clicked.connect(lambda: self.open_row_as_model())
        self.result_new_button = QtWidgets.QPushButton("New study...")
        self.result_new_button.clicked.connect(self.on_study)
        row.addWidget(self.csv_button)
        row.addWidget(self.open_button)
        row.addStretch(1)
        row.addWidget(self.result_new_button)
        lay.addLayout(row)
        self.result_tab = w
        tabs = self.win.result_tabs
        log = getattr(self.win, "log_view", None)
        idx = tabs.indexOf(log) if log is not None else -1
        if idx >= 0:                    # before Log, so Log stays last
            tabs.insertTab(idx, w, "Study")
        else:
            tabs.addTab(w, "Study")
        self._update_buttons()

    def _build_cancel_button(self):
        self.cancel_button = QtWidgets.QPushButton("Cancel study")
        self.cancel_button.setToolTip("Stop after the dome that is running.")
        self.cancel_button.setStyleSheet("padding: 1px 10px;")
        self.cancel_button.setVisible(False)
        self.cancel_button.clicked.connect(self.cancel)
        self.win.statusBar().addPermanentWidget(self.cancel_button)

    # ------------------------------------------------------- extension API --
    def on_model_changed(self, model):
        """The study does not depend on the model on screen: keep it."""

    def on_results(self, results):
        """Nothing to do: the study runs its own analyses."""

    # --------------------------------------------------------------- run --
    def _design_settings(self):
        dp = getattr(self.win, "design_panel", None)
        out = {}
        for key, attr in (("target_util", "target_util"),
                          ("max_iter", "max_iter")):
            w = getattr(dp, attr, None)
            if w is not None:
                out[key] = w.value()
        w = getattr(dp, "enforce_defl", None)
        if w is not None:
            out["enforce_deflection"] = w.isChecked()
        return out

    def _base_model(self):
        """The current model with every panel value read in (a deep copy).

        None when the window rejected an input (it has shown the message).
        """
        win = self.win
        reader = getattr(win, "_read_inputs", None)
        if reader is not None and reader() is False:
            return None
        return copy.deepcopy(win.model)

    def on_study(self):
        """Analysis > Parameter study...: ask for the inputs, then run."""
        win = self.win
        if win._busy:
            win.status("Busy. Wait for the current run to finish.")
            return
        try:
            base = self._base_model()
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Parameter study", str(exc))
            return
        if base is None:
            return
        problem = self._base_problem(base)
        if problem:
            QtWidgets.QMessageBox.information(win, "Parameter study", problem)
            return
        try:
            _, note = sy.uniform_section(base)
            note = "Sections: " + note + " (unless auto-sized)."
        except Exception:
            note = ""
        dlg = StudyDialog(base, win, design_note=note)
        ok = dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()
        if not ok:
            return
        v = dlg.values()
        n = len(v["frequencies"]) * len(v["radii"]) * len(v["fractions"])
        if n > MANY_POINTS and QtWidgets.QMessageBox.question(
                win, "Parameter study",
                f"This study has {n} domes and may take a few minutes. "
                f"You can cancel it at any time. Continue?") \
                != QtWidgets.QMessageBox.Yes:
            return
        self.run_study(base_model=base, **v)

    @staticmethod
    def _base_problem(model):
        if model.support_type == geo.CUSTOM:
            return ("The study rebuilds geodesic domes, so it cannot use "
                    "per-node (custom) supports. Choose pinned, roller or "
                    "elastic supports on the Geometry tab first.")
        return ""

    def run_study(self, frequencies, radii=None, fractions=None,
                  auto_size=False, do_buckling=True, base_model=None):
        """Start a study in the worker thread (no dialog).

        Returns True when the run started.
        """
        win = self.win
        if win._busy:
            win.status("Busy. Wait for the current run to finish.")
            return False
        if base_model is None:
            base_model = self._base_model()
            if base_model is None:
                return False
        problem = self._base_problem(base_model)
        if problem:
            QtWidgets.QMessageBox.information(win, "Parameter study", problem)
            return False
        self.study_base_model = base_model
        self._cancel_requested = False
        self._running = True
        self._pending_x = default_x({"frequencies": list(frequencies),
                                     "radii": list(radii or []),
                                     "fractions": list(fractions or [])})
        n = (len(frequencies) * max(1, len(radii or [1]))
             * max(1, len(fractions or [1])))
        win._set_busy(True, f"Parameter study: {n} dome"
                            f"{'s' if n != 1 else ''}...")
        self.cancel_button.setText("Cancel study")
        self.cancel_button.setEnabled(True)
        self.cancel_button.setVisible(True)
        win.log(f"Parameter study started: frequencies "
                f"{', '.join(str(f) for f in frequencies)}; radius "
                f"{', '.join(f'{r:g}' for r in (radii or [base_model.geom.radius]))} m; "
                f"fraction "
                f"{', '.join(f'{f:g}' for f in (fractions or [base_model.geom.dome_fraction]))}"
                f"{'; auto-sized' if auto_size else ''}.")
        kw = dict(radii=radii, fractions=fractions, auto_size=auto_size,
                  do_buckling=do_buckling, cancel=self._is_cancelled)
        if auto_size:
            kw.update(self._design_settings())
        worker = Worker(sy.parameter_study, base_model, list(frequencies), **kw)
        worker.progressed.connect(win._on_progress)
        worker.failed.connect(self._study_failed)
        worker.succeeded.connect(self._study_done)
        win._worker = worker
        worker.start()
        return True

    def _is_cancelled(self):
        return self._cancel_requested

    def cancel(self):
        """Stop the study after the dome that is running."""
        if not self._running:
            return
        self._cancel_requested = True
        self.cancel_button.setEnabled(False)
        self.cancel_button.setText("Stopping...")
        self.win.status("Stopping the study after the current dome...")

    def _finish(self):
        self._running = False
        self.cancel_button.setVisible(False)

    def _study_failed(self, message, tb):
        self._finish()
        self.win._on_failed(message, tb)

    def _study_done(self, rows):
        win = self.win
        self._finish()
        win._set_busy(False)
        self.study_rows = rows
        self.x_key = getattr(self, "_pending_x", None) or default_x(rows)
        self.x_combo.blockSignals(True)
        idx = self.x_combo.findData(self.x_key)
        self.x_combo.setCurrentIndex(max(idx, 0))
        self.x_combo.blockSignals(False)
        self.fill_table()
        self.redraw()
        win.view_tabs.setCurrentWidget(self.view_tab)
        win.result_tabs.setCurrentWidget(self.result_tab)
        for msg in getattr(rows, "messages", []):
            win.log("Study: " + msg)
        for r in rows:
            win.log(f"  {describe_row(r)}: {r.get('status')}, "
                    f"{format_value('mass_kg', r.get('mass_kg'))} kg, "
                    f"util {format_value('max_util', r.get('max_util'))}, "
                    f"lambda_cr {format_value('lambda_cr', r.get('lambda_cr'))}")
        stopped = any(m.startswith("Stopped by the user")
                      for m in getattr(rows, "messages", []))
        win.status(("Study stopped. " if stopped else "Study complete. ")
                   + summary_line(rows))

    # ------------------------------------------------------------ views --
    def _on_x_changed(self, _index=None):
        key = self.x_combo.currentData()
        if key:
            self.x_key = key
            self.redraw()

    def redraw(self):
        draw_study(self.canvas.figure, self.study_rows, self.x_key)
        self.canvas.canvas.draw_idle()

    def fill_table(self):
        rows = self.study_rows or []
        t = self.table
        t.setSortingEnabled(False)
        t.clearContents()
        t.setColumnCount(len(TABLE_COLUMNS))
        t.setHorizontalHeaderLabels([h for _, h in TABLE_COLUMNS])
        t.setRowCount(len(rows))
        bad = QColor("#ffe0e0")
        for i, r in enumerate(rows):
            flag = r.get("status") != "OK"
            for j, (key, _) in enumerate(TABLE_COLUMNS):
                v = r.get(key)
                item = QtWidgets.QTableWidgetItem(format_value(key, v))
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if key == "notes" and v:
                    item.setToolTip(str(v))
                if flag:
                    item.setBackground(bad)
                t.setItem(i, j, item)
        t.resizeColumnsToContents()
        self.summary_label.setText(summary_line(self.study_rows))
        self._update_buttons()

    def _update_buttons(self):
        has = bool(self.study_rows)
        self.csv_button.setEnabled(has)
        self.open_button.setEnabled(has and self.selected_index() is not None)

    def selected_index(self):
        rows = self.table.selectionModel().selectedRows() \
            if self.table.selectionModel() is not None else []
        if not rows:
            return None
        i = rows[0].row()
        return i if self.study_rows and 0 <= i < len(self.study_rows) else None

    # ---------------------------------------------------------- actions --
    def export_csv(self, path=None):
        """Write the study table to CSV (asks for the file when no path).

        Returns the path written, or None.
        """
        win = self.win
        if not self.study_rows:
            QtWidgets.QMessageBox.information(
                win, "Parameter study", "Run a parameter study first.")
            return None
        if path is None:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                win, "Export study table", "parameter_study.csv",
                "CSV (*.csv)")
            if not path:
                return None
        try:
            sy.study_to_csv(self.study_rows, path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Export CSV",
                                          f"Could not write the file.\n{exc}")
            return None
        win.log(f"Study table written to {path}")
        win.status("Study table written.")
        return path

    def open_row_as_model(self, index=None, confirm=True):
        """Make the selected study row the current model.

        ``index`` defaults to the selected table row.  Returns the new model
        or None.
        """
        win = self.win
        if win._busy:
            win.status("Busy. Wait for the current run to finish.")
            return None
        if index is None:
            index = self.selected_index()
        if not self.study_rows or index is None:
            QtWidgets.QMessageBox.information(
                win, "Parameter study", "Select a row of the Study table "
                                        "first.")
            return None
        row = self.study_rows[index]
        if row.get("status") == "ERROR":
            QtWidgets.QMessageBox.warning(
                win, "Parameter study",
                "This study point failed to run, so it cannot be opened.")
            return None
        if confirm and QtWidgets.QMessageBox.question(
                win, "Open study row",
                f"Replace the current model with the {describe_row(row)} "
                f"dome? Unsaved changes to the current model are lost.") \
                != QtWidgets.QMessageBox.Yes:
            return None
        base = self.study_base_model or copy.deepcopy(win.model)
        try:
            model = sy.model_from_row(base, row)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Open study row", str(exc))
            return None
        set_window_model(win, model, row)
        win.log(f"Opened study row as the model: {describe_row(row)} "
                f"({row.get('sections', '')}).")
        win.status(f"Model replaced by the {describe_row(row)} dome. "
                   f"Run the analysis (F9).")
        return model

    # ------------------------------------------------------------ events --
    def eventFilter(self, obj, event):
        # Closing the window stops a running study after the current dome.
        # Never raise here: Qt still calls this while the window is torn
        # down, when the Python attributes may already be gone.
        try:
            if self.__dict__.get("_running") and \
                    obj is self.__dict__.get("win") and \
                    event.type() == QtCore.QEvent.Close:
                self._cancel_requested = True
        except Exception:                                # pragma: no cover
            pass
        return False


# =============================================================== helpers ==
def _menu_title(action):
    return action.text().replace("&", "").strip().lower()


def _find_menu(win, title):
    for a in win.menuBar().actions():
        if a.menu() is not None and _menu_title(a) == title.lower():
            return a.menu()
    return None


def set_window_model(win, model, row=None):
    """Replace ``win.model`` and refresh the window.

    Uses ``win.set_model(model)`` when the window has it (it refreshes every
    panel and view and calls ``win.model_changed()``).  Otherwise it does the
    same refresh as File > Open and then calls ``win.model_changed()`` if
    that exists.
    """
    gp = getattr(win, "geometry_panel", None)
    g = model.geom
    if gp is not None and not g.is_custom:
        # the requested fraction rebuilds the same dome on Generate dome
        gp.frequency.setValue(int(g.frequency))
        gp.radius.setValue(float(row["radius_m"] if row else g.radius))
        gp.fraction.setValue(float(row["fraction_requested"] if row
                                   else g.dome_fraction))
        idx = gp.support.findData(model.support_type)
        if idx >= 0:
            gp.support.setCurrentIndex(idx)
    setter = getattr(win, "set_model", None)
    if callable(setter):
        setter(model)
    else:
        win.model = model
        win.results = None
        if hasattr(win, "sizing"):
            win.sizing = None
        for name in ("_clear_results",):
            fn = getattr(win, name, None)
            if callable(fn):
                fn()
        if gp is not None:
            gp.show_model(model)
        p = getattr(win, "section_panel", None)
        if p is not None:
            p.rebuild(model)
        for name in ("load_panel", "design_panel"):
            p = getattr(win, name, None)
            if p is not None:
                p.sync_from_model(model)
        for name in ("redraw_model", "redraw_loads", "model_changed"):
            fn = getattr(win, name, None)
            if callable(fn):
                fn()
    canvas = getattr(win, "canvas_model", None)
    if canvas is not None:
        win.view_tabs.setCurrentWidget(canvas)


def install(win):
    """Add the parameter study to ``win``; returns the extension object
    (the existing one if this window already has it)."""
    ext = getattr(win, "study_extension", None)
    if isinstance(ext, StudyExtension):
        return ext
    ext = StudyExtension(win)
    win.study_extension = ext
    return ext
