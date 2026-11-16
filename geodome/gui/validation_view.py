"""The Validation menu of the GUI.

Four independent checks of the program, each one a click away:

* **Verification suite** (``geodome.verification``, 75 checks): closed-form
  bars and tripods, published chord factors, Euler's formula, equilibrium
  and energy identities.  The window's own ``on_verify`` runs it.
* **Hand-calculation validation**: ``validation/run_validation.py`` (Models
  A and B solved by hand in ``validation/hand_calc.py``) and
  ``validation/run_validation_ext.py`` (cases F1 to F10 of
  ``validation/hand_calc_ext.py``: thermal, settlement, force method,
  springs, von Mises snap-through, natural frequencies, 2V angles).  The
  scripts are loaded by path with importlib (``validation/`` is not a
  package) and only their compare functions are called, so the files in
  ``validation/output`` are not rewritten.  Pass rules:
  run_validation rows pass when |p - h| / |h| <= 1e-12 (1e-10 %), or
  |p| <= 1e-9 in the row unit when the hand value is exactly 0;
  run_validation_ext rows use that script's own rule and tolerance
  (relative 1e-9 for linear and modal quantities, 1e-6 on the nonlinear
  path, exact counts, 5e-6 on published 5-decimal chord factors).
* **OpenSees cross-check**: the quick subset of
  ``validation/opensees_crosscheck.py`` (formulation probes, 1V/2V/4V
  statics and modes, nonlinear paths and exported scripts).  It needs an
  interpreter with openseespy (``geodome.opensees_bridge``); the menu entry
  is disabled with a tooltip when there is none.
* **Export OpenSeesPy script**: a standalone script of the current model
  and combination (``opensees_bridge.write_openseespy_script``).

``install(win)`` returns the extension object.  The last results are kept in
``.hand_result`` and ``.opensees_result`` (dicts, see the run functions).
"""

from __future__ import annotations

import contextlib
import copy
import importlib.util
import os
import sys
import tempfile
import time
from pathlib import Path

from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

import geodome
from . import style as st
from .worker import Worker

# The validation scripts live next to the package in the repository.
VALIDATION_DIR = Path(geodome.__file__).resolve().parents[1] / "validation"

HAND_SCRIPTS = ("run_validation.py", "run_validation_ext.py")
OPENSEES_SCRIPT = "opensees_crosscheck.py"

# The pass rule of run_validation.py (REL_TOL and ABS_TOL there); a test
# checks that these copies stay equal (SPEC: "80/80 within 1e-10 %").
RV_REL_TOL = 1e-12          # |p - h| / |h|
RV_ABS_TOL = 1e-9           # |p| in the row unit when h = 0

MISSING_FOLDER_TIP = ("The validation folder is not here. It comes with the "
                      "full GeoDomeLab repository from GitHub.")

SECTION_TITLES = {"P": "Formulation probes",
                  "S": "Linear static analysis",
                  "M": "Modal analysis",
                  "N": "Geometrically nonlinear analysis",
                  "X": "Exported OpenSeesPy scripts"}

HAND_COLUMNS = ("Case", "Quantity", "Hand", "Program", "Unit",
                "Difference", "Pass", "Tolerance", "Script")
OS_SUMMARY_COLUMNS = ("Section", "Pass", "Worst difference")
OS_COLUMNS = ("Case", "Structure", "Supports", "Combination", "Quantity",
              "Difference", "Pass", "Tolerance", "Note")
MAX_COLUMN_WIDTH = 300      # px; longer text is cut and shown in a tooltip

EXPORT_ANALYSES = (("Linear static", "static"),
                   ("Modal (12 modes, lumped mass)",
                    {"type": "modal", "n_modes": 12, "mass_type": "lumped"}),
                   ("Nonlinear path (crown uz)", "nonlinear"))


# ========================================================= script loading ==
_MODULES = {}


def missing_files(folder=None, names=HAND_SCRIPTS):
    """The script names that are not in the validation folder."""
    folder = Path(folder or VALIDATION_DIR)
    if not folder.is_dir():
        return list(names)
    return [n for n in names if not (folder / n).is_file()]


@contextlib.contextmanager
def _no_backend_switch():
    """run_validation.py calls matplotlib.use("Agg") at import; inside the
    GUI that must not touch the user's backend, so it is a no-op here."""
    try:
        import matplotlib
    except ImportError:                                  # pragma: no cover
        yield
        return
    old = matplotlib.use
    matplotlib.use = lambda *a, **k: None
    try:
        yield
    finally:
        matplotlib.use = old


def load_script(filename, folder=None):
    """Import ``validation/<filename>`` by path (cached per folder)."""
    folder = Path(folder or VALIDATION_DIR).resolve()
    path = folder / filename
    key = str(path)
    if key in _MODULES:
        return _MODULES[key]
    if not path.is_file():
        raise FileNotFoundError(f"{filename} is not in the validation "
                                f"folder.")
    name = "geodome_gui_" + Path(filename).stem
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))     # the scripts import hand_calc
    with _no_backend_switch():
        spec.loader.exec_module(mod)
    _MODULES[key] = mod
    return mod


# ======================================================= hand calculation ==
def _hand_rows_ab(rv):
    """Models A and B: call the compare functions, never main()."""
    rv.ROWS.clear()
    rv.validate_a(rv.hc.model_a())
    rv.validate_b(rv.hc.model_b())
    out = []
    for r in rv.ROWS:
        h, p = float(r["hand"]), float(r["program"])
        if r.get("diff_kind") == "%":
            d = float(r["difference"]) / 100.0
            ok = d <= RV_REL_TOL
            tol = "1e-12 (rel.)"
        else:
            d = float(r["difference"])
            ok = d <= RV_ABS_TOL
            tol = f"1e-9 {r['unit']} (abs.)"
        if "pass" in r:                     # the script's own verdict
            ok = ok and bool(r["pass"])
        out.append({"source": "Models A, B", "case": f"{r['model']}: {r['case']}",
                    "quantity": r["quantity"], "hand": h, "program": p,
                    "unit": r["unit"], "difference": d, "tolerance": tol,
                    "pass": bool(ok)})
    rv.ROWS.clear()
    return out


def _ext_tol_text(kind, tol):
    if kind == "count":
        return "exact"
    if kind == "published":
        return f"{tol:.0e} (abs.)"
    if kind == "identity":
        return f"{tol:.0e} (value)"
    return f"{tol:.0e} (rel.)"


def _hand_rows_ext(ve, progress, cancel):
    """Cases F1 to F10, case by case (run_comparison without the files)."""
    cases = getattr(ve, "CASES", None)
    if cases and hasattr(ve, "Comparison"):
        C = ve.Comparison(ve.HAND)
        for i, fn in enumerate(cases):
            if cancel is not None and cancel():
                return None, {}
            progress(0.15 + 0.8 * i / len(cases),
                     f"Hand-calculation validation: case F{i + 1}...")
            fn(C)
        rows = C.evaluate()
        info = {"missing": sorted(set(ve.HAND) - C.used),
                "n_hand": len(ve.HAND), "n_used": len(C.used)}
    else:                                            # pragma: no cover
        rows, info = ve.run_comparison()
    out = []
    for r in rows:
        out.append({"source": "Cases F1-F10", "case": r["sub"],
                    "quantity": r["quantity"], "hand": float(r["hand"]),
                    "program": float(r["program"]), "unit": r["unit"],
                    "difference": float(r["difference"]),
                    "tolerance": _ext_tol_text(r["kind"], r["tolerance"]),
                    "pass": bool(r["pass"])})
    return out, info


def run_hand_validation(folder=None, progress=None, cancel=None):
    """Compare the program with every hand calculation (no files written).

    Returns ``{"rows", "n", "n_pass", "ab", "ext", "info", "elapsed",
    "cancelled"}`` where ``ab``/``ext`` are ``(n_pass, n)`` per script.
    """
    progress = progress or (lambda f, m="": None)
    t0 = time.perf_counter()
    progress(0.02, "Loading the validation scripts...")
    rv = load_script("run_validation.py", folder)
    ve = load_script("run_validation_ext.py", folder)
    progress(0.08, "Hand-calculation validation: Models A and B...")
    ab = _hand_rows_ab(rv)
    ext, info = _hand_rows_ext(ve, progress, cancel)
    cancelled = ext is None
    rows = ab + (ext or [])
    progress(1.0, "Hand-calculation validation done.")
    return {"rows": rows, "n": len(rows),
            "n_pass": sum(r["pass"] for r in rows),
            "ab": (sum(r["pass"] for r in ab), len(ab)),
            "ext": (sum(r["pass"] for r in ext or []), len(ext or [])),
            "info": info, "elapsed": time.perf_counter() - t0,
            "cancelled": cancelled}


def hand_header(res):
    if res is None:
        return ("No hand-calculation run yet. Use Validation > Run "
                "hand-calculation validation.")
    (pa, na), (pe, ne) = res["ab"], res["ext"]
    text = (f"Hand calculations: {res['n_pass']} of {res['n']} comparisons "
            f"pass.  Models A and B (run_validation.py): {pa} of {na}.  "
            f"Cases F1-F10 (run_validation_ext.py): {pe} of {ne}")
    info = res.get("info") or {}
    if info.get("n_hand"):
        text += (f", {info['n_used']} of {info['n_hand']} hand values "
                 f"compared")
    text += f".  {res['elapsed']:.1f} s."
    if res.get("cancelled"):
        text += "  Stopped by the user before the end."
    return text


# ============================================================== OpenSees ==
def opensees_state():
    """``(available, explanation)`` from ``geodome.opensees_bridge``."""
    try:
        from .. import opensees_bridge as ob
    except Exception as exc:                         # pragma: no cover
        return False, f"The OpenSees bridge did not load: {exc}"
    try:
        status = ob.opensees_status()
    except Exception as exc:                         # pragma: no cover
        return False, f"Could not look for OpenSees: {exc}"
    mode = status.get("mode")
    if mode == "in-process":
        return True, "Uses openseespy in this Python."
    if mode == "subprocess":
        return True, "Uses the separate OpenSees Python found on this computer."
    return False, (f"OpenSees was not found. Install openseespy, or set "
                   f"{ob.ENV_VAR} to a Python that has it.")


def run_opensees_quick(folder=None, progress=None, cancel=None):
    """The quick subset of opensees_crosscheck.py, without its files.

    Returns ``{"rows", "summary", "n", "n_pass", "mode", "elapsed",
    "cancelled"}``.
    """
    from .. import opensees_bridge as ob
    progress = progress or (lambda f, m="": None)
    t0 = time.perf_counter()
    status = ob.opensees_status()
    if status.get("mode") is None:
        raise ob.OpenSeesUnavailable(
            "OpenSees was not found. Install openseespy, or set "
            f"{ob.ENV_VAR} to a Python that has it.")
    progress(0.02, "Loading the OpenSees cross-check...")
    xc = load_script(OPENSEES_SCRIPT, folder)
    cancelled = False
    names = ("section_probes", "section_static", "section_modal",
             "section_nonlinear", "section_scripts")
    if all(hasattr(xc, n) for n in names) and hasattr(xc, "Rows"):
        rows, log = xc.Rows(), []
        steps = (("P", lambda: xc.section_probes(rows, log)),
                 ("S", lambda: xc.section_static(rows, True, log)),
                 ("M", lambda: xc.section_modal(rows, True, log)),
                 ("N", lambda: xc.section_nonlinear(rows, True, log)),
                 ("X", lambda: xc.section_scripts(rows, log)))
        for i, (key, fn) in enumerate(steps):
            if cancel is not None and cancel():
                cancelled = True
                break
            progress(0.05 + 0.9 * i / len(steps),
                     f"OpenSees cross-check: {SECTION_TITLES[key].lower()}...")
            fn()
        raw = rows.rows
    else:                                            # pragma: no cover
        with tempfile.TemporaryDirectory(prefix="geodome_xc_") as tmp:
            rows, _, _ = xc.run(quick=True, out_dir=tmp, verbose=False)
        raw = rows.rows
    tol_text = getattr(xc, "_tol_text", None)
    out = []
    for r in raw:
        tol = (tol_text(r) if tol_text else str(r.get("tolerance")))
        out.append({"section": r["section"], "case": r["case"],
                    "structure": r["structure"], "supports": r["supports"],
                    "combination": r["combination"],
                    "quantity": r["quantity"],
                    "difference": float(r["difference"]),
                    "tolerance": tol, "pass": bool(r["pass"]),
                    "note": r.get("note", ""), "kind": r.get("kind", "max")})
    progress(1.0, "OpenSees cross-check done.")
    return {"rows": out, "summary": opensees_summary(out), "n": len(out),
            "n_pass": sum(r["pass"] for r in out), "mode": status["mode"],
            "elapsed": time.perf_counter() - t0, "cancelled": cancelled}


def opensees_summary(rows):
    """One line per section: cases, comparisons, pass, worst difference."""
    out = []
    for key, title in SECTION_TITLES.items():
        rs = [r for r in rows if r["section"] == key]
        if not rs:
            continue
        worst = None
        for r in rs:
            if r.get("kind", "max") != "max":
                continue
            if worst is None or r["difference"] > worst["difference"]:
                worst = r
        out.append({"section": f"{key}  {title}",
                    "cases": len({r["case"] for r in rs}),
                    "n": len(rs), "n_pass": sum(r["pass"] for r in rs),
                    "worst": worst["difference"] if worst else float("nan"),
                    "worst_quantity": (f"{worst['quantity']} ({worst['case']})"
                                       if worst else "-")})
    return out


def opensees_header(res):
    if res is None:
        return ("No OpenSees cross-check yet. Use Validation > Run OpenSees "
                "cross-check.")
    mode = {"in-process": "openseespy in this Python",
            "subprocess": "a separate OpenSees Python"}.get(res["mode"],
                                                           str(res["mode"]))
    text = (f"OpenSees cross-check (quick subset): {res['n_pass']} of "
            f"{res['n']} comparisons pass.  Run with {mode} in "
            f"{res['elapsed']:.1f} s.  The full matrix: python "
            f"validation/opensees_crosscheck.py")
    if res.get("cancelled"):
        text = "Stopped by the user.  " + text
    return text


# =============================================================== format ==
def _num(v, digits=7):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if v != v:
        return "n/a"
    if v == 0.0:
        return "0"
    return f"{v:.{digits}g}"


def _diff(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if v != v:
        return "n/a"
    if v == 0.0:
        return "0"
    return f"{v:.1e}"


def _fill(table, headers, cells, failed, numeric=()):
    """Fill a QTableWidget from text cells; failed rows are tinted red."""
    table.setSortingEnabled(False)
    table.clearContents()
    table.setColumnCount(len(headers))
    table.setHorizontalHeaderLabels(list(headers))
    table.setRowCount(len(cells))
    bad = QColor("#ffe0e0")
    for i, (row, flag) in enumerate(zip(cells, failed)):
        for j, text in enumerate(row):
            item = QtWidgets.QTableWidgetItem(str(text))
            if j in numeric:
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if len(str(text)) > 40:
                item.setToolTip(str(text))
            if flag:
                item.setBackground(bad)
            table.setItem(i, j, item)
    table.resizeColumnsToContents()
    for j in range(len(headers)):
        if table.columnWidth(j) > MAX_COLUMN_WIDTH:
            table.setColumnWidth(j, MAX_COLUMN_WIDTH)


# ================================================================ dialog ==
class ExportScriptDialog(QtWidgets.QDialog):
    """Combination and analysis type of an exported OpenSeesPy script."""

    def __init__(self, model, current=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export OpenSeesPy script")
        self.setMinimumWidth(420)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(st.hint(
            "Writes a standalone Python script that builds this model in "
            "OpenSees and prints the results. Run it with a Python that has "
            "openseespy."))
        form = QtWidgets.QFormLayout()
        self.combo = QtWidgets.QComboBox()
        for c in model.combinations:
            self.combo.addItem(c.name)
        if current:
            i = self.combo.findText(current)
            if i < 0:                  # "LC2" matches "LC2  1.2D + 1.6S"
                names = [c.name.split()[0] if c.name.split() else ""
                         for c in model.combinations]
                key = str(current).split()[0] if str(current).split() else ""
                i = names.index(key) if key in names else -1
            if i >= 0:
                self.combo.setCurrentIndex(i)
        form.addRow("Load combination", self.combo)
        self.analysis = QtWidgets.QComboBox()
        for label, _ in EXPORT_ANALYSES:
            self.analysis.addItem(label)
        self.analysis.currentIndexChanged.connect(self._update)
        form.addRow("Analysis", self.analysis)
        lay.addLayout(form)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.button(QtWidgets.QDialogButtonBox.Ok).setText("Save script...")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        self._update()

    def _update(self, _i=None):
        # a modal analysis has no load combination
        self.combo.setEnabled(self.analysis.currentIndex() != 1)

    def values(self):
        label, analysis = EXPORT_ANALYSES[self.analysis.currentIndex()]
        combo = None if self.analysis.currentIndex() == 1 \
            else self.combo.currentText()
        return combo, copy.deepcopy(analysis), label


# ============================================================= extension ==
class ValidationExtension(QtCore.QObject):
    """The Validation menu and the Validation result tab."""

    def __init__(self, win, folder=None):
        super().__init__(win)
        self.win = win
        self.folder = Path(folder) if folder else VALIDATION_DIR
        self.hand_result = None
        self.opensees_result = None
        self._cancel_requested = False
        self._running = None
        self._build_menu()
        self._build_result_tab()
        self._build_cancel_button()
        self.refresh_actions()

    # ------------------------------------------------------------ build --
    def _build_menu(self):
        win = self.win
        add_menu = getattr(win, "add_menu", None)
        if callable(add_menu):                 # inserted just before Help
            self.menu = add_menu("&Validation")
        else:
            mb = win.menuBar()
            self.menu = QtWidgets.QMenu("&Validation", win)
            help_action = None
            for a in mb.actions():
                if a.menu() is not None and \
                        a.text().replace("&", "").strip().lower() == "help":
                    help_action = a
                    break
            if help_action is not None:
                mb.insertMenu(help_action, self.menu)
            else:
                mb.addMenu(self.menu)
        self.menu.setToolTipsVisible(True)

        self.act_verify = win._act(
            "Run &verification suite (75 checks)", self.run_verification,
            tip="Closed-form, published and identity checks of the solver.")
        self.act_hand = win._act(
            "Run &hand-calculation validation", self.run_hand_validation,
            tip="Compare the program with the hand calculations of "
                "Models A and B and cases F1-F10.")
        self.act_opensees = win._act(
            "Run &OpenSees cross-check", self.run_opensees_crosscheck,
            tip="Solve a set of models again in OpenSees and compare.")
        # QAction.triggered passes checked=False; it must not become path.
        self.act_export = win._act(
            "&Export OpenSeesPy script...", lambda _=False: self.export_script(),
            tip="Write a standalone OpenSeesPy script of the current model.")
        self.act_folder = win._act(
            "Open validation report (Markdown) &folder", self.open_folder,
            tip="Open validation/output, which holds the Markdown reports.")
        for a in (self.act_verify, self.act_hand, self.act_opensees):
            self.menu.addAction(a)
        self.menu.addSeparator()
        self.menu.addAction(self.act_export)
        self.menu.addSeparator()
        self.menu.addAction(self.act_folder)
        self.menu.aboutToShow.connect(self.refresh_actions)

    def _build_result_tab(self):
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)

        # hand calculations
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        top = QtWidgets.QHBoxLayout()
        self.hand_label = QtWidgets.QLabel(hand_header(None))
        self.hand_label.setWordWrap(True)
        top.addWidget(self.hand_label, 1)
        self.fail_only = QtWidgets.QCheckBox("Show failures only")
        self.fail_only.toggled.connect(self._fill_hand)
        top.addWidget(self.fail_only)
        lay.addLayout(top)
        self.hand_table = st.table(list(HAND_COLUMNS))
        lay.addWidget(self.hand_table, 1)
        self.tabs.addTab(w, "Hand calculations")

        # OpenSees
        w2 = QtWidgets.QWidget()
        lay2 = QtWidgets.QVBoxLayout(w2)
        lay2.setContentsMargins(6, 6, 6, 6)
        lay2.setSpacing(6)
        self.os_label = QtWidgets.QLabel(opensees_header(None))
        self.os_label.setWordWrap(True)
        lay2.addWidget(self.os_label)
        split = QtWidgets.QSplitter(Qt.Horizontal)
        self.os_summary_table = st.table(list(OS_SUMMARY_COLUMNS))
        self.os_summary_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SingleSelection)
        self.os_summary_table.setToolTip(
            "Click a section to show only its rows. Ctrl+click it to show "
            "every row again.")
        self.os_summary_table.itemSelectionChanged.connect(
            self._fill_opensees_rows)
        self.os_table = st.table(list(OS_COLUMNS))
        split.addWidget(self.os_summary_table)
        split.addWidget(self.os_table)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([380, 800])
        lay2.addWidget(split, 1)
        self.tabs.addTab(w2, "OpenSees")

        self.result_tab = self.tabs
        _add_result_tab(self.win, self.tabs, "Validation")

    def _build_cancel_button(self):
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.setToolTip("Stop after the step that is running.")
        self.cancel_button.setStyleSheet("padding: 1px 10px;")
        self.cancel_button.setVisible(False)
        self.cancel_button.clicked.connect(self.cancel)
        self.win.statusBar().addPermanentWidget(self.cancel_button)

    # ------------------------------------------------------- extension API --
    def on_model_changed(self, model):
        """Validation does not depend on the model on screen."""

    def on_results(self, results):
        """Nothing to do."""

    # ------------------------------------------------------------ state --
    def refresh_actions(self):
        """Enable or disable each action, with a tooltip saying why."""
        win = self.win
        has_verify = callable(getattr(win, "on_verify", None))
        self.act_verify.setEnabled(has_verify)
        if not has_verify:                           # pragma: no cover
            self.act_verify.setToolTip("Not available in this window.")

        miss = missing_files(self.folder, HAND_SCRIPTS)
        _set_state(self.act_hand, not miss,
                   "Compare the program with the hand calculations of "
                   "Models A and B and cases F1-F10.",
                   MISSING_FOLDER_TIP if not self.folder.is_dir() else
                   f"Missing in the validation folder: {', '.join(miss)}.")

        if missing_files(self.folder, (OPENSEES_SCRIPT,)):
            ok, why = False, (MISSING_FOLDER_TIP if not self.folder.is_dir()
                              else f"{OPENSEES_SCRIPT} is missing in the "
                                   f"validation folder.")
        else:
            ok, why = opensees_state()
        _set_state(self.act_opensees, ok,
                   "Solve a set of models again in OpenSees and compare. "
                   + why, why)
        self.opensees_reason = why

        out = self.folder / "output"
        _set_state(self.act_folder, self.folder.is_dir(),
                   "Open validation/output, which holds the Markdown "
                   "reports.", MISSING_FOLDER_TIP)
        self.act_export.setEnabled(True)
        state = {"hand": self.act_hand.isEnabled(),
                 "opensees": self.act_opensees.isEnabled(),
                 "folder": self.act_folder.isEnabled(),
                 "output": out.is_dir()}
        if getattr(win, "_busy", False):
            # a run is busy: nothing new may start until it ends
            for a in (self.act_verify, self.act_hand, self.act_opensees,
                      self.act_export):
                a.setEnabled(False)
        return state

    def _busy(self):
        if self.win._busy:
            self.win.status("Busy. Wait for the current run to finish.")
            return True
        return False

    # ------------------------------------------------------------ verify --
    def run_verification(self):
        fn = getattr(self.win, "on_verify", None)
        if callable(fn):
            fn()

    # ------------------------------------------------------ worker runs --
    def _start(self, kind, func, message, done):
        win = self.win
        self._cancel_requested = False
        self._running = kind
        win._set_busy(True, message)
        self.cancel_button.setText("Cancel")
        self.cancel_button.setEnabled(True)
        self.cancel_button.setVisible(True)
        worker = Worker(func, self.folder, cancel=self._is_cancelled)
        worker.progressed.connect(win._on_progress)
        worker.failed.connect(self._failed)
        worker.succeeded.connect(done)
        win._worker = worker
        worker.start()
        return True

    def _is_cancelled(self):
        return self._cancel_requested

    def cancel(self):
        if self._running is None:
            return
        self._cancel_requested = True
        self.cancel_button.setEnabled(False)
        self.cancel_button.setText("Stopping...")
        self.win.status("Stopping after the current step...")

    def _finish(self):
        self._running = None
        self.cancel_button.setVisible(False)

    def _failed(self, message, tb):
        self._finish()
        self.win._on_failed(message, tb)

    def run_hand_validation(self):
        """Validation > Run hand-calculation validation."""
        if self._busy():
            return False
        if not self.refresh_actions()["hand"]:
            QtWidgets.QMessageBox.information(
                self.win, "Validation", self.act_hand.toolTip())
            return False
        return self._start("hand", run_hand_validation,
                           "Running the hand-calculation validation...",
                           self._hand_done)

    def _hand_done(self, res):
        win = self.win
        self._finish()
        win._set_busy(False)
        self.hand_result = res
        self._fill_hand()
        self.tabs.setCurrentIndex(0)
        win.result_tabs.setCurrentWidget(self.result_tab)
        win.log(hand_header(res))
        win.status(f"Hand-calculation validation: {res['n_pass']} of "
                   f"{res['n']} comparisons pass.")
        if res["n_pass"] < res["n"]:
            QtWidgets.QMessageBox.warning(
                win, "Validation",
                f"{res['n'] - res['n_pass']} comparison(s) failed. They are "
                f"shown in red in the Validation tab.")

    def _fill_hand(self, *_):
        res = self.hand_result
        self.hand_label.setText(hand_header(res))
        rows = (res or {}).get("rows", [])
        if self.fail_only.isChecked():
            rows = [r for r in rows if not r["pass"]]
        cells = [(r["case"], r["quantity"], _num(r["hand"]),
                  _num(r["program"]), r["unit"], _diff(r["difference"]),
                  "PASS" if r["pass"] else "FAIL", r["tolerance"],
                  r["source"])
                 for r in rows]
        _fill(self.hand_table, HAND_COLUMNS, cells,
              [not r["pass"] for r in rows], numeric=(2, 3, 5))

    def run_opensees_crosscheck(self):
        """Validation > Run OpenSees cross-check (quick subset)."""
        if self._busy():
            return False
        if not self.refresh_actions()["opensees"]:
            QtWidgets.QMessageBox.information(
                self.win, "OpenSees cross-check", self.opensees_reason)
            return False
        return self._start("opensees", run_opensees_quick,
                           "Running the OpenSees cross-check...",
                           self._opensees_done)

    def _opensees_done(self, res):
        win = self.win
        self._finish()
        win._set_busy(False)
        self.opensees_result = res
        self._fill_opensees()
        self.tabs.setCurrentIndex(1)
        win.result_tabs.setCurrentWidget(self.result_tab)
        win.log(opensees_header(res))
        for r in res["rows"]:
            if not r["pass"]:
                win.log(f"  FAIL {r['case']} {r['combination']} "
                        f"{r['quantity']}: {r['difference']:.3e}")
        win.status(f"OpenSees cross-check: {res['n_pass']} of {res['n']} "
                   f"comparisons pass.")
        if res["n_pass"] < res["n"]:
            QtWidgets.QMessageBox.warning(
                win, "OpenSees cross-check",
                f"{res['n'] - res['n_pass']} comparison(s) failed. They are "
                f"shown in red in the Validation tab.")

    def _fill_opensees(self):
        res = self.opensees_result
        self.os_label.setText(opensees_header(res))
        summ = (res or {}).get("summary", [])
        self.os_summary_table.blockSignals(True)
        _fill(self.os_summary_table, OS_SUMMARY_COLUMNS,
              [(s["section"], f"{s['n_pass']} of {s['n']}",
                _diff(s["worst"])) for s in summ],
              [s["n_pass"] < s["n"] for s in summ], numeric=(1, 2))
        for i, s in enumerate(summ):
            tip = (f"{s['cases']} cases, {s['n']} comparisons. Largest "
                   f"difference: {s['worst_quantity']}.")
            for j in range(len(OS_SUMMARY_COLUMNS)):
                self.os_summary_table.item(i, j).setToolTip(tip)
        self.os_summary_table.blockSignals(False)
        self._fill_opensees_rows()

    def _fill_opensees_rows(self):
        """The detail table: every row, or the selected section only."""
        res = self.opensees_result
        rows = (res or {}).get("rows", [])
        summ = (res or {}).get("summary", [])
        sel = self.os_summary_table.selectionModel()
        picked = sel.selectedRows() if sel is not None else []
        if picked and picked[0].row() < len(summ):
            key = summ[picked[0].row()]["section"].split()[0]
            rows = [r for r in rows if r["section"] == key]
        _fill(self.os_table, OS_COLUMNS,
              [(r["case"], r["structure"], r["supports"], r["combination"],
                r["quantity"], _diff(r["difference"]),
                "PASS" if r["pass"] else "FAIL", r["tolerance"], r["note"])
               for r in rows],
              [not r["pass"] for r in rows], numeric=(5,))

    # ------------------------------------------------------------ export --
    def _current_combination(self):
        win = self.win
        combo = getattr(win, "axial_combo", None)
        if combo is not None and combo.currentText():
            return combo.currentText()
        res = getattr(win, "results", None)
        if res is not None and getattr(res, "governing_case", None):
            return res.governing_case
        return None

    def export_script(self, path=None, combination=None, analysis=None):
        """Validation > Export OpenSeesPy script...

        With ``path`` given no dialog is shown (``combination`` defaults to
        the current one, ``analysis`` to "static").  Returns the path or None.
        """
        win = self.win
        if self._busy():
            return None
        try:
            reader = getattr(win, "_read_inputs", None)
            if reader is not None and reader() is False:
                return None             # the window has shown the message
            model = copy.deepcopy(win.model)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Export OpenSeesPy script",
                                          str(exc))
            return None
        if not model.combinations:
            QtWidgets.QMessageBox.information(
                win, "Export OpenSeesPy script",
                "The model has no load combination.")
            return None
        if path is None:
            dlg = ExportScriptDialog(model, self._current_combination(), win)
            ok = dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()
            if not ok:
                return None
            combination, analysis, _ = dlg.values()
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                win, "Export OpenSeesPy script",
                script_name(model, combination, analysis),
                "Python script (*.py)")
            if not path:
                return None
        else:
            analysis = "static" if analysis is None else analysis
            if combination is None and not _is_modal(analysis):
                combination = (self._current_combination()
                               or model.combinations[0].name)
        if _is_modal(analysis):
            combination = None
        try:
            from .. import opensees_bridge as ob
            ob.write_openseespy_script(model, combination, path, analysis)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                win, "Export OpenSeesPy script",
                f"Could not write the script.\n{type(exc).__name__}: {exc}")
            return None
        win.log(f"OpenSeesPy script written to {path} "
                f"(combination: {combination or 'none, modal'}). Run it with "
                f"a Python that has openseespy: python "
                f"{os.path.basename(path)}")
        win.status("OpenSeesPy script written.")
        return path

    # ------------------------------------------------------------ folder --
    def open_folder(self):
        """Open validation/output (or the validation folder) in the file
        browser."""
        win = self.win
        if not self.folder.is_dir():
            QtWidgets.QMessageBox.information(win, "Validation",
                                              MISSING_FOLDER_TIP)
            return False
        target = self.folder / "output"
        if not target.is_dir():
            target = self.folder
        ok = QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(target)))
        if not ok:
            QtWidgets.QMessageBox.information(
                win, "Validation",
                "Could not open the file browser. The reports are in the "
                "validation/output folder of the repository.")
        else:
            win.status("Opened the validation report folder.")
        return bool(ok)


def _is_modal(analysis):
    return analysis == "modal" or (isinstance(analysis, dict)
                                   and analysis.get("type") == "modal")


def script_name(model, combination, analysis):
    """Suggested file name, e.g. opensees_6V_LC2.py."""
    g = model.geom
    stem = "custom" if g.is_custom else f"{g.frequency}V"
    if _is_modal(analysis):
        return f"opensees_{stem}_modal.py"
    token = "".join(ch if ch.isalnum() else "_"
                    for ch in (str(combination or "").split() or ["load"])[0])
    kind = "_nonlinear" if analysis in ("nonlinear", "path") or (
        isinstance(analysis, dict) and analysis.get("type") == "path") else ""
    return f"opensees_{stem}_{token}{kind}.py"


def _set_state(action, enabled, tip_on, tip_off):
    action.setEnabled(bool(enabled))
    tip = tip_on if enabled else tip_off
    action.setToolTip(tip)
    action.setStatusTip(tip)


def _add_result_tab(win, widget, title):
    """Add a result tab just before the Log tab (Log stays last)."""
    tabs = win.result_tabs
    log = getattr(win, "log_view", None)
    idx = tabs.indexOf(log) if log is not None else -1
    if idx >= 0:
        tabs.insertTab(idx, widget, title)
    else:
        tabs.addTab(widget, title)


def install(win):
    """Add the Validation menu and result tab to ``win``; returns the
    extension object (the existing one if this window already has it)."""
    ext = getattr(win, "validation_extension", None)
    if isinstance(ext, ValidationExtension):
        return ext
    ext = ValidationExtension(win)
    win.validation_extension = ext
    return ext
