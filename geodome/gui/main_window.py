"""The GeoDomeLab main window.

Public hooks for add-on modules (``install(win)``)
--------------------------------------------------
``win.extensions``
    List of extension objects.  Each may define ``on_model_changed(model)``
    and ``on_results(results)``; errors raised there are logged, never
    propagated.
``win.model_changed()``
    Call after any change of model (new, open, example, generate, CSV
    import).  Updates the window title and tells every extension.
``win.set_model(model, title=None, message=None)``
    Replace the model, refresh every panel and view, then call
    ``model_changed()``.
``win.menu_file``, ``win.menu_examples``, ``win.menu_analysis``,
``win.menu_view``, ``win.menu_help``, ``win.toolbar``, ``win.add_menu(title)``
    Menus and the main toolbar.  ``add_menu`` inserts a new top-level menu
    just before Help.
``win.busy_actions``
    Actions and buttons that start a long run.  ``_set_busy`` turns them off
    while a worker runs and on again after it; then it calls
    ``ext.refresh_actions()`` on every extension that has one.

The window loads the add-ons in :data:`EXTENSION_MODULES` order (nonlinear,
modal, fabrication, parameter study, validation) at the end of
``__init__``.  Each extension is kept in ``win.extensions`` and as
``win.<name>_extension``.  ``win.extension_results()`` collects their last
results for the HTML report (``report.export_html(..., extras=...)``).
"""

from __future__ import annotations

import copy
import html
import json
import math
import os
import sys
import traceback
import webbrowser

import numpy as np
from qtpy import QtWidgets, QtGui
from qtpy.QtCore import Qt

from .. import __version__, __author__, __url__
from .. import benchmarks as bm
from .. import csv_import
from .. import geometry as geo
from .. import loads as ld
from .. import plotting as plot
from .. import report as rp
from .. import verification as vf
from ..model import DomeModel, LoadSpec, default_model, \
    default_combinations
from ..optimize import auto_size
from ..postprocess import run_analysis, summary_text, member_table, \
    node_table, reaction_table
from .canvas import Canvas3D, DualCanvas2D
from .panels import (GeometryPanel, SectionPanel, LoadPanel, DesignPanel,
                     InputError)
from .worker import Worker, plain_message
from . import style as st

APP_NAME = "GeoDomeLab"
APP_SUBTITLE = "Geodesic Dome Space-Truss Analyzer"
COURSE_LINE = "CE 7405 Statically Indeterminate Structures, LSU, 2026"
LICENSE_NOTE = ("All rights reserved. Free to use for personal, educational "
                "and research calculations with credit to the author. "
                "See LICENSE.")

# Add-on modules loaded by MainWindow, in this order: (module in geodome.gui,
# window attribute that keeps its extension object).
EXTENSION_MODULES = (
    ("nonlinear_view", "nonlinear_extension"),
    ("modal_view", "modal_extension"),
    ("fabrication_view", "fabrication_extension"),
    ("study_view", "study_extension"),
    ("validation_view", "validation_extension"),
)

CSV_FORMAT_HELP = """
Import truss from CSV -- file format
====================================

File > Import truss from CSV... asks for the node file.  The other files are
read from the same folder.  A common prefix is allowed: picking
tower_nodes.csv also reads tower_members.csv, tower_supports.csv, ...

Separator: comma, semicolon or tab.  A header row is optional.  Blank lines
and lines that start with # are skipped.  Units: m, kN and kN/m.

nodes.csv      (required)   id, x, y, z
    Node id (any whole number) and coordinates in m.  z is up.

members.csv    (required)   i, j[, group]
    Node ids at the two ends of each bar.  The optional group is the strut
    family (give it on every row or on none).  Without groups the families
    are formed by bar length.

supports.csv   (optional)   node, cx, cy, cz
    For each direction: 1 = fixed, 0 = free, any other positive number =
    spring stiffness in kN/m.  The words fixed and free also work.
    Without this file the lowest nodes are pinned.

loads.csv      (optional)   node, Fx, Fy, Fz[, case]
    Nodal force in kN along the global axes.  Case is D, S, L or W
    (default L).

faces.csv      (optional)   i, j, k
    Triangular cladding panels.  Only needed for area loads (dead, snow,
    live, wind).  Without it those loads are zero.

Example (nodes.csv):
    id, x, y, z
    1, 0.0, 0.0, 3.0
    2, -4.0, 0.0, 0.0
    3,  4.0, 0.0, 0.0

Node ids are renumbered 0, 1, 2, ... in id order, because the program
numbers nodes from 0.  The model starts with self weight on, area loads at
zero, wind off, the default section on every family and the ASCE 7-16
strength combinations with live load.  Check the Sections and Loads tabs,
then run the analysis (F9).  A ready example is in examples/csv_truss.
"""


def _check_catalogue_names(model):
    """Raise ValueError when the material or a section is not in the
    program's catalogues (e.g. a file from a newer or edited version)."""
    from ..materials import MATERIALS
    from ..sections import SECTIONS
    if model.material_name not in MATERIALS:
        raise ValueError(f"Unknown material {model.material_name!r}.")
    bad = sorted({s for s in model.group_sections if s not in SECTIONS})
    if bad:
        raise ValueError("Unknown section(s): " + ", ".join(bad) + ".")


# The dome fraction stored in a model file is the snapped one; from the
# 0.25-0.95 input range it lands between 0.237 and 0.970 (1V to 10V).
FILE_FRACTION_RANGE = (0.20, 0.98)


def _model_file_problem(data, freq_range, radius_range):
    """A plain sentence when ``data`` (the parsed JSON of a model file) is
    not a model this window can show, else None.

    A geodesic file must stay inside the Geometry tab's input ranges, so
    that the tab shows the dome that is on screen (and Generate dome builds
    it again).
    """
    if not isinstance(data, dict):
        return "This is not a GeoDomeLab model file."
    gd = data.get("geometry")
    if not isinstance(gd, dict):
        return ('This is not a GeoDomeLab model file: the "geometry" '
                'section is missing.')
    if gd.get("kind", "geodesic") != "geodesic":
        return None
    for key in ("frequency", "radius", "dome_fraction"):
        if key not in gd:
            return (f'This model file is not complete: the geometry has no '
                    f'"{key}".')
    try:
        nu = float(gd["frequency"])
        R = float(gd["radius"])
        f = float(gd["dome_fraction"])
    except (TypeError, ValueError):
        return "The geometry in this file has a value that is not a number."
    lo, hi = freq_range
    if not (math.isfinite(nu) and nu == int(nu) and lo <= nu <= hi):
        return (f"The frequency in this file is {gd['frequency']}. It must "
                f"be a whole number from {lo} to {hi}.")
    lo, hi = radius_range
    if not (math.isfinite(R) and lo <= R <= hi):
        return (f"The sphere radius in this file is {R:g} m. It must be "
                f"between {lo:g} and {hi:g} m.")
    lo, hi = FILE_FRACTION_RANGE
    if not (math.isfinite(f) and lo <= f <= hi):
        return (f"The dome fraction in this file is {f:g}. It must be "
                f"between {lo:g} and {hi:g}.")
    return None


def _model_title(model):
    """A short name for the model: geometry title (+ radius for a dome)."""
    g = model.geom
    if g.is_custom:
        return g.title or "Custom truss"
    title = g.title or f"{g.frequency}V geodesic dome"
    if math.isfinite(g.radius):
        title += f", R = {g.radius:g} m"
    return title


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Light palette first, so every widget (and the matplotlib toolbar
        # icons) is built light even when Windows runs in dark mode.
        st.apply_light_theme()
        # the mouse wheel scrolls panels; it never edits an unfocused box
        st.install_wheel_guard()
        self.setWindowTitle(f"{APP_NAME} {__version__}  --  {APP_SUBTITLE}")
        self.resize(1480, 920)
        self.setStyleSheet(st.STYLESHEET)

        self.model = default_model(4, 12.0, 0.5)
        self.model_title = _model_title(self.model)
        self.results = None
        self.sizing = None
        self.verification_rows = None
        self._worker = None
        self._busy = False
        # add-on modules (see the module docstring)
        self.extensions = []
        # actions and buttons that start a long run; off while a run is busy
        self.busy_actions = []

        self._build_ui()
        self._build_menu()
        self._build_toolbar()

        self.geometry_panel.sync_from_model(self.model)
        self.section_panel.rebuild(self.model)
        self.load_panel.sync_from_model(self.model)
        self.design_panel.sync_from_model(self.model)
        self._fill_case_combos()
        self.redraw_model()
        self._install_extensions()
        self.view_tabs.setCurrentWidget(self.canvas_model)
        self.result_tabs.setCurrentWidget(self.summary_view)
        self.log(f"{APP_NAME} {__version__} ready. "
                 f"A default 4V dome of 12 m radius has been generated.")
        self.status("Ready.")
        self.model_changed()

    # ============================================================ add-ons ==
    def _install_extensions(self):
        """Load the add-on modules in :data:`EXTENSION_MODULES` order.

        Each module's ``install(win)`` adds its menu entries, toolbar button
        and tabs and returns an extension object, kept in
        :attr:`extensions` (once) and as ``win.<name>_extension``.  A module
        that fails to load is logged; the window still opens without it.
        """
        import importlib
        self.extension_errors = []
        for name, attr in EXTENSION_MODULES:
            try:
                mod = importlib.import_module(f".{name}", __package__)
                ext = mod.install(self)
            except Exception as exc:
                self.extension_errors.append(traceback.format_exc())
                self.log(f"ERROR: the add-on {name} could not be loaded: "
                         f"{exc}")
                self.log(traceback.format_exc())
                continue
            if ext is None:
                continue
            if not any(e is ext for e in self.extensions):
                self.extensions.append(ext)
            setattr(self, attr, ext)
        self._arrange_extensions()

    def _arrange_extensions(self):
        """Settle what the add-ons share: menus, toolbar, tabs, busy state."""
        nl = getattr(self, "nonlinear_extension", None)
        md = getattr(self, "modal_extension", None)
        fb = getattr(self, "fabrication_extension", None)
        sy = getattr(self, "study_extension", None)
        va = getattr(self, "validation_extension", None)

        # One "Run verification suite" entry: the Validation menu has it now.
        if va is not None and getattr(va, "act_verify", None) is not None:
            self.menu_analysis.removeAction(self.act_verify_menu)
            self.menu_analysis.removeAction(self._verify_sep)
            va.act_verify.setShortcut("F11")

        # Toolbar: the original buttons first, then Nonlinear and Modes.
        tb = self.toolbar
        extra = [a for a in (getattr(nl, "act_toolbar", None),
                             getattr(md, "act_modes", None)) if a is not None]
        for a in extra:
            tb.removeAction(a)
        if extra:
            tb.addSeparator()
            for a in extra:
                tb.addAction(a)

        # View tabs: Nonlinear, Modes, Study after the original ones.
        order = [getattr(nl, "view", None), getattr(md, "view", None),
                 getattr(sy, "view_tab", None)]
        pos = self.view_tabs.indexOf(self.canvas_charts) + 1
        for w in order:
            i = self.view_tabs.indexOf(w) if w is not None else -1
            if i < 0:
                continue
            if i != pos:
                self.view_tabs.tabBar().moveTab(i, pos)
            pos += 1

        # Result tabs: the new ones before Log, Log last.
        i = self.result_tabs.indexOf(self.log_view)
        last = self.result_tabs.count() - 1
        if 0 <= i != last:
            self.result_tabs.tabBar().moveTab(i, last)

        # Long runs: these are off while a worker runs.
        acts = [self.act_generate_menu, self.act_analyse_menu,
                self.act_size_menu, self.act_verify_menu]
        for ext, names in ((nl, ("act_nonlinear", "act_imperfection",
                                 "act_toolbar")),
                           (md, ("act_modal", "act_modes", "run_btn")),
                           (fb, ()),
                           (sy, ("act_study", "view_new_button",
                                 "result_new_button")),
                           (va, ("act_verify", "act_hand", "act_opensees",
                                 "act_export"))):
            for n in names:
                a = getattr(ext, n, None) if ext is not None else None
                if a is not None and a not in acts:
                    acts.append(a)
        self.busy_actions = acts

    def extension_results(self):
        """The add-ons' last results, keyed as :func:`report.export_html`
        takes them in ``extras`` (only the ones that exist)."""
        out = {}
        for ext in self.extensions:
            for key in rp.EXTRA_KEYS:
                if key in out:
                    continue
                try:
                    value = getattr(ext, key, None)
                except Exception:                    # pragma: no cover
                    value = None
                if value is not None:
                    out[key] = value
        return out

    # =================================================================== UI ==
    def _build_ui(self):
        splitter = QtWidgets.QSplitter(Qt.Horizontal)
        self.setCentralWidget(splitter)

        # ---- left: input dock ---------------------------------------------
        self.input_tabs = QtWidgets.QTabWidget()
        self.geometry_panel = GeometryPanel()
        self.section_panel = SectionPanel()
        self.load_panel = LoadPanel()
        self.design_panel = DesignPanel()
        self.input_tabs.addTab(self.geometry_panel, "Geometry")
        self.input_tabs.addTab(self.section_panel, "Sections")
        self.input_tabs.addTab(self.load_panel, "Loads")
        self.input_tabs.addTab(self.design_panel, "Design")
        self.input_tabs.setMinimumWidth(392)
        self.input_tabs.setMaximumWidth(520)
        splitter.addWidget(self.input_tabs)

        self.geometry_panel.generate_requested.connect(self.on_generate)
        self.geometry_panel.supports_changed.connect(self._on_supports_changed)
        self.section_panel.changed.connect(self.on_section_changed)
        self.load_panel.changed.connect(self._on_loads_edited)
        self.design_panel.analyse_requested.connect(self.on_analyse)
        self.design_panel.size_requested.connect(self.on_auto_size)
        self.design_panel.verify_requested.connect(self.on_verify)

        # ---- right: plots over results ------------------------------------
        right = QtWidgets.QSplitter(Qt.Vertical)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([400, 1080])

        self.view_tabs = QtWidgets.QTabWidget()
        right.addWidget(self.view_tabs)

        # model view
        self.canvas_model = Canvas3D()
        self.chk_panels = QtWidgets.QCheckBox("Cladding panels")
        self.chk_nodes = QtWidgets.QCheckBox("Nodes")
        self.chk_nodes.setChecked(True)
        self.chk_labels = QtWidgets.QCheckBox("Node numbers")
        self.chk_family = QtWidgets.QCheckBox("Colour by strut family")
        self.chk_family.setChecked(True)
        for w in (self.chk_family, self.chk_nodes, self.chk_panels,
                  self.chk_labels):
            self.canvas_model.add_control(w)
            w.toggled.connect(self.redraw_model)
        self.canvas_model.add_control_stretch()
        self.view_tabs.addTab(self.canvas_model, "Model")

        # load view
        self.canvas_loads = Canvas3D()
        self.loadcase_combo = QtWidgets.QComboBox()
        self.loadcase_combo.setMinimumWidth(230)
        self.loadcase_combo.currentIndexChanged.connect(self.redraw_loads)
        self.canvas_loads.add_control(QtWidgets.QLabel("Load case:"))
        self.canvas_loads.add_control(self.loadcase_combo)
        self.canvas_loads.add_control_stretch()
        self.view_tabs.addTab(self.canvas_loads, "Loads")

        # deflected view
        self.canvas_defl = Canvas3D()
        self.defl_combo = QtWidgets.QComboBox()
        self.defl_combo.setMinimumWidth(230)
        self.defl_combo.currentIndexChanged.connect(self.redraw_deflected)
        self.defl_scale = st.dspin(1.0, 100000.0, 200.0, 50.0, 0)
        self.defl_scale.setPrefix("x ")
        self.defl_scale.valueChanged.connect(self.redraw_deflected)
        auto_btn = QtWidgets.QPushButton("Auto scale")
        auto_btn.clicked.connect(self._auto_scale)
        self.chk_undeformed = QtWidgets.QCheckBox("Show undeformed")
        self.chk_undeformed.setChecked(True)
        self.chk_undeformed.toggled.connect(self.redraw_deflected)
        for w in (QtWidgets.QLabel("Case:"), self.defl_combo,
                  QtWidgets.QLabel("Magnification:"), self.defl_scale,
                  auto_btn, self.chk_undeformed):
            self.canvas_defl.add_control(w)
        self.canvas_defl.add_control_stretch()
        self.view_tabs.addTab(self.canvas_defl, "Deflected shape")

        # axial force view
        self.canvas_axial = Canvas3D()
        self.axial_combo = QtWidgets.QComboBox()
        self.axial_combo.setMinimumWidth(230)
        self.axial_combo.currentIndexChanged.connect(self.redraw_axial)
        self.chk_width = QtWidgets.QCheckBox("Line width by force")
        self.chk_width.setChecked(True)
        self.chk_width.toggled.connect(self.redraw_axial)
        for w in (QtWidgets.QLabel("Case:"), self.axial_combo, self.chk_width):
            self.canvas_axial.add_control(w)
        self.canvas_axial.add_control_stretch()
        self.view_tabs.addTab(self.canvas_axial, "Axial force")

        # utilisation view
        self.canvas_util = Canvas3D()
        self.view_tabs.addTab(self.canvas_util, "Utilisation")

        # buckling view
        self.canvas_buck = Canvas3D()
        self.mode_spin = st.spin(1, 20, 1, 1)
        self.mode_spin.setPrefix("mode ")
        self.mode_spin.valueChanged.connect(self.redraw_buckling)
        self.buck_scale = st.dspin(0.01, 0.6, 0.10, 0.01, 2)
        self.buck_scale.setPrefix("amplitude ")
        self.buck_scale.valueChanged.connect(self.redraw_buckling)
        self.buck_label = QtWidgets.QLabel("")
        self.buck_label.setObjectName("metric")
        for w in (self.mode_spin, self.buck_scale, self.buck_label):
            self.canvas_buck.add_control(w)
        self.canvas_buck.add_control_stretch()
        self.view_tabs.addTab(self.canvas_buck, "Buckling mode")

        # charts
        self.canvas_charts = DualCanvas2D()
        self.view_tabs.addTab(self.canvas_charts, "Charts")

        # ---- results ------------------------------------------------------
        self.result_tabs = QtWidgets.QTabWidget()
        right.addWidget(self.result_tabs)
        right.setSizes([620, 300])

        self.summary_view = QtWidgets.QPlainTextEdit()
        self.summary_view.setReadOnly(True)
        self.result_tabs.addTab(self.summary_view, "Summary")

        self.member_table = st.table(["member"])
        self.result_tabs.addTab(self.member_table, "Members")

        self.node_table = st.table(["node"])
        self.result_tabs.addTab(self.node_table, "Nodes")

        self.reaction_table = st.table(["node"])
        self.result_tabs.addTab(self.reaction_table, "Reactions")

        self.sizing_table = st.table(["iteration"])
        self.result_tabs.addTab(self.sizing_table, "Auto-sizing")

        self.verify_table = st.table(["check"])
        self.result_tabs.addTab(self.verify_table, "Verification")

        self.log_view = QtWidgets.QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.result_tabs.addTab(self.log_view, "Log")

        # ---- status bar ---------------------------------------------------
        self.progress = QtWidgets.QProgressBar()
        self.progress.setMaximumWidth(230)
        self.progress.setVisible(False)
        self.status_label = QtWidgets.QLabel("")
        sb = self.statusBar()
        sb.addWidget(self.status_label, 1)
        sb.addPermanentWidget(self.progress)

    # ------------------------------------------------------------- actions --
    def _act(self, text, slot, shortcut=None, tip=None):
        a = QtWidgets.QAction(text, self)
        a.triggered.connect(slot)
        if shortcut:
            a.setShortcut(shortcut)
        if tip:
            a.setToolTip(tip)
            a.setStatusTip(tip)
        return a

    def _build_menu(self):
        mb = self.menuBar()

        m = mb.addMenu("&File")
        self.menu_file = m
        self.act_new = self._act("&New model", self.on_new, "Ctrl+N",
                                 tip="Start again from the default 4V dome")
        self.act_open = self._act("&Open model...", self.on_open, "Ctrl+O",
                                  tip="Open a GeoDomeLab model file (.json)")
        self.act_save = self._act("&Save model...", self.on_save, "Ctrl+S",
                                  tip="Save the model with every input "
                                      "(.json)")
        self.act_import = self._act("&Import truss from CSV...",
                                    self.on_import_csv,
                                    tip="Read nodes, members, supports and "
                                        "loads from CSV files")
        # off while a run is busy (the run uses the model)
        self.file_actions = (self.act_new, self.act_open, self.act_save,
                             self.act_import)
        for a in self.file_actions:
            m.addAction(a)
        m.addSeparator()
        m.addAction(self._act("Export results to &CSV...", self.on_export_csv))
        m.addAction(self._act("Export &figures (PNG)...", self.on_export_figs))
        m.addAction(self._act("Export &HTML report...", self.on_export_html,
                              "Ctrl+R"))
        m.addAction(self._act("Export &everything to a folder...",
                              self.on_export_all))
        m.addSeparator()
        m.addAction(self._act("E&xit", self.close, "Ctrl+Q"))

        m = mb.addMenu("&Examples")
        self.menu_examples = m
        m.setToolTipsVisible(True)
        self.example_actions = {}
        for key in bm.preset_keys():
            title, desc = bm.preset_info(key)
            a = self._act(title, lambda _=False, k=key: self.load_preset(k))
            a.setStatusTip(desc)
            a.setToolTip(f"<p>{html.escape(desc)}</p>")
            a.setData(key)
            m.addAction(a)
            self.example_actions[key] = a

        m = mb.addMenu("&Analysis")
        self.menu_analysis = m
        self.act_generate_menu = self._act("&Generate dome", self.on_generate,
                                           "F5")
        self.act_analyse_menu = self._act("&Run analysis", self.on_analyse,
                                          "F9")
        self.act_size_menu = self._act("&Auto-size strut families",
                                       self.on_auto_size, "F10")
        for a in (self.act_generate_menu, self.act_analyse_menu,
                  self.act_size_menu):
            m.addAction(a)
        # Moved to the Validation menu (with F11) when that add-on loads.
        self._verify_sep = m.addSeparator()
        self.act_verify_menu = self._act("Run &verification suite",
                                         self.on_verify, "F11")
        m.addAction(self.act_verify_menu)

        m = mb.addMenu("&View")
        self.menu_view = m
        # Rebuilt each time it opens, so view tabs added later by add-on
        # modules are listed too.
        m.aboutToShow.connect(self._fill_view_menu)
        m.triggered.connect(self._on_view_menu)
        self._fill_view_menu()

        m = mb.addMenu("&Help")
        self.menu_help = m
        m.addAction(self._act("&Theory and method", self.on_theory))
        m.addAction(self._act("&CSV import format", self.on_csv_help))
        m.addAction(self._act("&About", self.on_about))

    def _fill_view_menu(self):
        """List the view tabs.  The actions belong to the menu (clear()
        deletes them) and carry the tab index as data; one slot,
        :meth:`_on_view_menu`, handles them all.  (Lambda slots that hold
        ``self`` kept every closed window alive, and collecting such a
        window later could crash the process.)"""
        m = self.menu_view
        m.clear()
        self._view_actions = []
        for i in range(self.view_tabs.count()):
            a = QtWidgets.QAction(self.view_tabs.tabText(i), m)
            a.setData(i)
            m.addAction(a)
            self._view_actions.append(a)
        m.addSeparator()
        a = QtWidgets.QAction("Save current view as PNG...", m)
        a.setData(-1)
        m.addAction(a)
        self._view_actions.append(a)

    def _on_view_menu(self, action):
        k = action.data()
        if k == -1:
            self.on_save_view()
        elif isinstance(k, int) and 0 <= k < self.view_tabs.count():
            self.view_tabs.setCurrentIndex(k)

    def add_menu(self, title):
        """Add a top-level menu just before Help and return it."""
        menu = QtWidgets.QMenu(title, self)
        self.menuBar().insertMenu(self.menu_help.menuAction(), menu)
        return menu

    def _build_toolbar(self):
        tb = QtWidgets.QToolBar("Main")
        tb.setToolButtonStyle(Qt.ToolButtonTextOnly)
        tb.setMovable(False)
        self.addToolBar(tb)
        self.toolbar = tb
        self.act_generate = self._act("Generate dome", self.on_generate,
                                      tip="Rebuild the geometry (F5)")
        self.act_analyse = self._act("Run analysis", self.on_analyse,
                                     tip="Solve every load combination (F9)")
        self.act_size = self._act("Auto-size", self.on_auto_size,
                                  tip="Fully stressed design (F10)")
        self.act_verify = self._act("Verify", self.on_verify,
                                    tip="Run the benchmark suite (F11)")
        self.act_report = self._act("Export report", self.on_export_html,
                                    tip="Self-contained HTML calculation "
                                        "report (Ctrl+R)")
        for a in (self.act_generate, self.act_analyse, self.act_size):
            tb.addAction(a)
        tb.addSeparator()
        tb.addAction(self.act_verify)
        tb.addSeparator()
        tb.addAction(self.act_report)

    # =============================================================== utils ==
    def log(self, text):
        self.log_view.appendPlainText(str(text))

    def status(self, text):
        self.status_label.setText(str(text))

    def _set_busy(self, busy, message=""):
        self._busy = busy
        for w in (self.input_tabs, self.act_generate, self.act_analyse,
                  self.act_size, self.act_verify, self.act_report,
                  *getattr(self, "file_actions", ())):
            w.setEnabled(not busy)
        # the add-ons' long runs (Nonlinear, Modes, study, validation, ...)
        for w in getattr(self, "busy_actions", ()):
            try:
                w.setEnabled(not busy)
            except RuntimeError:                     # deleted C++ object
                pass
        # a model must not be swapped while a run uses it
        for m in (getattr(self, "menu_examples", None),):
            if m is not None:
                m.setEnabled(not busy)
        if not busy:
            # an add-on may keep an action off for its own reason (for
            # example no OpenSees found); let it set its state again
            for ext in list(getattr(self, "extensions", ())):
                fn = getattr(ext, "refresh_actions", None)
                if callable(fn):
                    try:
                        fn()
                    except Exception as exc:
                        self.log(f"ERROR in {type(ext).__name__}."
                                 f"refresh_actions: {exc}")
        self.progress.setVisible(busy)
        if busy:
            self.progress.setValue(0)
            self.status(message or "Working...")
            QtWidgets.QApplication.setOverrideCursor(
                QtGui.QCursor(Qt.WaitCursor))
        else:
            QtWidgets.QApplication.restoreOverrideCursor()

    def _on_progress(self, fraction, message):
        self.progress.setValue(int(round(100 * fraction)))
        if message:
            self.status(message)

    def _on_failed(self, message, tb):
        self._set_busy(False)
        self.status("Failed.")
        self.log("ERROR: " + message)
        self.log(tb)
        QtWidgets.QMessageBox.critical(self, "Analysis failed",
                                       plain_message(message))

    def _read_inputs(self):
        """Pull every panel's values into the model.

        Returns True when every input is valid.  On a bad input it shows a
        message, logs it and returns False.  It never raises.
        """
        try:
            self.geometry_panel.read_into(self.model)
            self.section_panel.read_into(self.model)
            self.design_panel.read_into(self.model)
            self.load_panel.read_into(self.model)
        except InputError as exc:
            self._input_error(str(exc))
            return False
        except Exception as exc:                     # never crash on input
            self._input_error(f"The inputs could not be read: {exc}")
            self.log(traceback.format_exc())
            return False
        return True

    def _input_error(self, message):
        self.log("INPUT: " + message.replace("\n\n", "\n"))
        self.status("Check the inputs. Nothing was run.")
        QtWidgets.QMessageBox.warning(self, "Check the inputs", message)

    def _busy_message(self):
        """True (with a status message) while a run is busy.  The exports
        that draw figures must wait: they switch the busy state off when
        they end."""
        if self._busy:
            self.status("Busy. Wait for the current run to finish.")
            return True
        return False

    def _require_results(self):
        if self.results is None:
            QtWidgets.QMessageBox.information(
                self, APP_NAME, "Run the analysis first (F9).")
            return False
        return True

    def _inputs_changed(self):
        """The model was edited in place (sections, material, supports or
        auto-sizing): drop every result that belongs to the old inputs,
        here and in the add-ons (modal, nonlinear, ...)."""
        had = self.results is not None
        self.results = None
        self.sizing = None
        self._clear_results()
        self.model_changed()
        if had:
            self.status("Inputs changed. Press F9 to analyse again.")
            self.log("Inputs changed: the results were cleared. Press F9 "
                     "to analyse again.")

    def _on_loads_edited(self):
        """A Loads tab input changed.  The model takes it at the next run;
        until then the results on screen belong to the old loads."""
        if self.results is not None and not self._busy:
            self.status("Inputs changed. Press F9 to analyse again.")

    # ======================================================= model change ==
    def model_changed(self):
        """Public hook: call after any change of model.

        Updates the window title and calls ``on_model_changed(model)`` on
        every extension in :attr:`extensions`.  An error in an extension is
        logged and does not stop the others.
        """
        self.setWindowTitle(f"{APP_NAME} {__version__}  --  "
                            f"{self.model_title}")
        for ext in list(self.extensions):
            fn = getattr(ext, "on_model_changed", None)
            if fn is None:
                continue
            try:
                fn(self.model)
            except Exception as exc:
                self.log(f"ERROR in {type(ext).__name__}.on_model_changed: "
                         f"{exc}")
                self.log(traceback.format_exc())

    def _notify_results(self):
        for ext in list(self.extensions):
            fn = getattr(ext, "on_results", None)
            if fn is None:
                continue
            try:
                fn(self.results)
            except Exception as exc:
                self.log(f"ERROR in {type(ext).__name__}.on_results: {exc}")
                self.log(traceback.format_exc())

    def set_model(self, model, title=None, message=None, sync_inputs=True):
        """Replace the model, refresh every panel and view, then call
        :meth:`model_changed`.

        ``title`` is shown in the window title (default: the geometry title).
        ``sync_inputs=False`` keeps what the user typed in the Loads and
        Design tabs and the geodesic inputs (used by Generate dome, which has
        already read them into the model).
        """
        self.model = model
        self.model_title = title or _model_title(model)
        self.results = None
        self.sizing = None
        self.geometry_panel.sync_from_model(model, set_inputs=sync_inputs)
        self.section_panel.rebuild(model)
        if sync_inputs:
            self.load_panel.sync_from_model(model)
            self.design_panel.sync_from_model(model)
        else:
            self.load_panel.show_model_info(model)
        self._clear_results()
        self._fill_case_combos()
        self.redraw_model()
        self.redraw_loads()
        try:
            self.canvas_charts.ax_top.clear()
            plot.draw_strut_histogram(self.canvas_charts.ax_top, model)
            self.canvas_charts.refresh()
        except Exception as exc:                     # pragma: no cover
            self.log(f"Strut chart not drawn: {exc}")
        if message:
            self.log(message)
        self.model_changed()

    def _fill_case_combos(self):
        """List the model's combinations in the Loads view (before a run)."""
        names = [c.name for c in self.model.combinations]
        self.loadcase_combo.blockSignals(True)
        self.loadcase_combo.clear()
        self.loadcase_combo.addItems(names)
        self.loadcase_combo.blockSignals(False)
        for combo in (self.defl_combo, self.axial_combo):
            combo.blockSignals(True)
            combo.clear()
            combo.blockSignals(False)

    # ============================================================= actions ==
    def on_new(self):
        if self._busy_message():
            return
        self.geometry_panel.frequency.setValue(4)
        self.geometry_panel.radius.setValue(12.0)
        self.geometry_panel.fraction.setValue(0.5)
        self.set_model(default_model(4, 12.0, 0.5),
                       message="New model created.")
        self.view_tabs.setCurrentWidget(self.canvas_model)
        self.status("New model.")

    def load_preset(self, key):
        """Load the example ``key`` of :data:`geodome.benchmarks.PRESETS`."""
        if self._busy:
            return False
        try:
            title, desc = bm.preset_info(key)
            model = bm.make_preset(key)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self, "Examples", f"The example could not be built:\n{exc}")
            return False
        self.set_model(model, title=title,
                       message=f"Example loaded: {title}.\n  {desc}")
        for note in model.geom.notes:
            self.log("  " + note)
        self.view_tabs.setCurrentWidget(self.canvas_model)
        self.status(f"Example loaded: {title}.  Press F9 to analyse.")
        return True

    def _clear_results(self):
        self.summary_view.clear()
        for t in (self.member_table, self.node_table, self.reaction_table,
                  self.sizing_table):
            t.clearContents()
            t.setRowCount(0)
        for c in (self.canvas_defl, self.canvas_axial, self.canvas_util,
                  self.canvas_buck):
            c.clear()
            c.refresh()
        self.canvas_charts.ax_top.clear()
        self.canvas_charts.ax_bottom.clear()
        self.canvas_charts.refresh()

    def on_generate(self):
        if self._busy:
            return
        freq, radius, frac = self.geometry_panel.params()
        try:
            g = geo.generate_dome(freq, radius, frac)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Geometry", str(exc))
            return
        old = self.model
        from_custom = old.geom.is_custom
        # Per-node data belongs to the old node numbering.  Keep it only when
        # the numbering cannot have changed (same geodesic breakdown).
        same_numbering = (not old.geom.is_custom
                          and old.geom.frequency == g.frequency
                          and old.geom.n_nodes == g.n_nodes)
        # An example truss has its own loads (often none but a nodal load
        # that does not exist on the dome): a dome gets the default loads.
        loads = LoadSpec() if from_custom else old.loads.copy()
        combinations = (default_combinations() if from_custom
                        else old.combinations)
        if not same_numbering:
            had = (self.load_panel.settle_table.rowCount()
                   + self.load_panel.nodal_table.rowCount())
            self.load_panel.clear_node_tables()
            loads.settlements = []
            loads.nodal_loads = []
            if had:
                self.log("Settlements and nodal loads were cleared: the node "
                         "numbers changed.")
        if self.geometry_panel.support.currentData() == geo.CUSTOM:
            self.geometry_panel.set_support_type(geo.PINNED)
            self.log("Custom supports belong to the old model.  The new dome "
                     "uses the pinned base ring.")
        keep = list(old.group_sections)
        self.model = DomeModel(geom=g, material_name=old.material_name,
                               group_sections=keep,
                               support_type=geo.PINNED
                               if old.support_type == geo.CUSTOM
                               else old.support_type,
                               loads=loads,
                               combinations=combinations,
                               K_effective=old.K_effective,
                               deflection_limit=old.deflection_limit,
                               phi_c=old.phi_c, phi_t=old.phi_t,
                               design_thickness_factor=getattr(
                                   old, "design_thickness_factor", 1.0),
                               spring_stiffness=old.spring_stiffness)
        # Rebuild the section rows before reading the panels back: the old
        # rows belong to the previous geometry and may be fewer than the new
        # strut families (e.g. 4V has 6, 6V has 9).
        self.section_panel.rebuild(self.model)
        if from_custom:
            # show the default loads first: _read_inputs reads the tab back
            self.load_panel.sync_from_model(self.model)
            self.log("The example's loads and combinations were replaced by "
                     "the default dome loads.")
        self._read_inputs()
        self.set_model(self.model, sync_inputs=False)
        self.view_tabs.setCurrentWidget(self.canvas_model)
        for note in g.notes:
            self.log(note)
        self.status(f"Generated {g.frequency}V dome: {g.n_nodes} nodes, "
                    f"{g.n_members} struts.")

    def _on_supports_changed(self):
        """Keep the model summary in step with the support inputs."""
        if self._busy:
            return
        try:
            self.geometry_panel.read_into(self.model)
        except Exception:
            return          # e.g. Custom with no data: reported at Analyse
        self.geometry_panel.show_model(self.model)
        self.load_panel.show_model_info(self.model)
        self._inputs_changed()

    def on_section_changed(self):
        self.section_panel.read_into(self.model)
        self.section_panel.update_mass(self.model)
        self.geometry_panel.show_model(self.model)
        if self.chk_family.isChecked():
            self.redraw_model()
        self._inputs_changed()

    def _model_state(self):
        try:
            return json.dumps(self.model.to_dict(), sort_keys=True)
        except Exception:                            # pragma: no cover
            return None

    def on_analyse(self):
        if self._busy:
            return
        before = self._model_state()
        if not self._read_inputs():
            return
        if before is None or self._model_state() != before:
            # Loads or design inputs changed: the add-ons' results (modal,
            # nonlinear, ...) belong to the old inputs.
            self.model_changed()
        self._set_busy(True, "Running analysis...")
        # The results keep their own copy of the model they were computed
        # for, so later edits never pair new inputs with old results.
        self._worker = Worker(
            run_analysis, copy.deepcopy(self.model),
            n_buckling_modes=self.design_panel.n_modes.value(),
            do_buckling=self.design_panel.do_buckling.isChecked())
        self._worker.progressed.connect(self._on_progress)
        self._worker.failed.connect(self._on_failed)
        self._worker.succeeded.connect(self._analysis_done)
        self._worker.start()

    def _analysis_done(self, res):
        self._set_busy(False)
        self.results = res
        self.summary_view.setPlainText(summary_text(res))

        rows = member_table(res)
        st.fill_table(self.member_table, rows,
                      highlight=lambda r: r["utilisation"] > 1.0)
        st.fill_table(self.node_table, node_table(res))
        st.fill_table(self.reaction_table, reaction_table(res))

        for combo in (self.loadcase_combo, self.defl_combo, self.axial_combo):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(res.order)
            combo.blockSignals(False)
        gov = res.governing_case or res.order[0]
        self.axial_combo.setCurrentText(gov)
        self.loadcase_combo.setCurrentText(gov)
        self.defl_combo.setCurrentText(res.service_case)

        nmodes = (len(res.buckling.load_factors)
                  if res.buckling is not None else 0)
        self.mode_spin.setMaximum(max(1, nmodes))
        self.mode_spin.setEnabled(nmodes > 0)

        self._auto_scale()
        self.redraw_loads()
        self.redraw_axial()
        self.redraw_utilisation()
        self.redraw_buckling()
        self.redraw_charts()
        self.view_tabs.setCurrentWidget(self.canvas_util)
        self.result_tabs.setCurrentIndex(0)

        for msg in res.messages:
            self.log(msg)
        flag = "OK" if res.governing_util <= 1.0 else "OVERSTRESSED"
        strength = res.strength_cases()
        loaded = [n for n in strength
                  if np.any(res.cases[n].F != 0.0)
                  or res.cases[n].eps0 is not None
                  or res.cases[n].prescribed]
        if not strength:
            note = ("No strength combination: untick Service on at least "
                    "one combination.")
            self.log(note)
            self.status(f"Analysis complete. {note}")
        elif not loaded:
            note = "No load is applied. Check the Loads tab."
            self.log(note)
            self.status(f"Analysis complete. {note}")
        else:
            self.status(f"Analysis complete. Max utilisation "
                        f"{res.governing_util:.3f} ({flag}) under "
                        f"{res.governing_case}.")
        self._notify_results()

    def on_auto_size(self):
        if self._busy:
            return
        if not self._read_inputs():
            return
        self._set_busy(True, "Auto-sizing...")
        self._worker = Worker(
            auto_size, self.model,
            target_util=self.design_panel.target_util.value(),
            max_iter=self.design_panel.max_iter.value(),
            enforce_deflection=self.design_panel.enforce_defl.isChecked(),
            slenderness_limit=self.design_panel.slender_limit.value())
        self._worker.progressed.connect(self._on_progress)
        self._worker.failed.connect(self._on_failed)
        self._worker.succeeded.connect(self._sizing_done)
        self._worker.start()

    def _sizing_done(self, sizing):
        self._set_busy(False)
        # auto_size changed the sections in place: the old results (and the
        # add-ons' modal and nonlinear results) belong to the old sections
        self._inputs_changed()
        self.sizing = sizing
        st.fill_table(self.sizing_table, sizing.table())
        self.section_panel.sync_from_model(self.model)
        self.geometry_panel.show_model(self.model)
        self.result_tabs.setCurrentWidget(self.sizing_table)
        direction = "lighter" if sizing.saving >= 0 else "heavier"
        self.log(f"Auto-sizing {'converged' if sizing.converged else 'stopped'} "
                 f"after {sizing.iterations} iterations: "
                 f"{sizing.start_mass:,.0f} kg -> {sizing.final_mass:,.0f} kg "
                 f"({abs(sizing.saving) * 100:.1f} % {direction}).")
        for m in sizing.messages:
            self.log(m)
        self.status(f"Auto-sizing done: {sizing.final_mass:,.0f} kg "
                    f"({abs(sizing.saving) * 100:.1f} % {direction}). "
                    f"Re-analysing...")
        self.on_analyse()

    def on_verify(self):
        if self._busy:
            return
        self._set_busy(True, "Running verification suite...")
        self._worker = Worker(vf.run_all)
        self._worker.progressed.connect(self._on_progress)
        self._worker.failed.connect(self._on_failed)
        self._worker.succeeded.connect(self._verify_done)
        self._worker.start()

    def _verify_done(self, rows):
        self._set_busy(False)
        self.verification_rows = rows
        table_rows = [{
            "group": r["group"],
            "check": r["check"],
            "computed": r["computed"],
            "reference": r["reference"],
            "unit": r["unit"],
            "rel. error": r["rel_error"],
            "tolerance": r["tolerance"],
            "result": "PASS" if r["pass"] else "FAIL",
            "basis": r["note"],
        } for r in rows]
        st.fill_table(self.verify_table, table_rows,
                      highlight=lambda r: r["result"] == "FAIL",
                      numeric_digits=6)
        self.result_tabs.setCurrentWidget(self.verify_table)
        n_pass = sum(1 for r in rows if r["pass"])
        self.log(vf.format_report(rows))
        self.status(f"Verification: {n_pass} of {len(rows)} checks passed.")
        if n_pass < len(rows):
            QtWidgets.QMessageBox.warning(
                self, "Verification",
                f"{len(rows) - n_pass} check(s) failed. See the "
                f"Verification tab.")

    # =============================================================== plots ==
    def redraw_model(self):
        c = self.canvas_model
        ax = c.clear()
        plot.draw_model(ax, self.model,
                        color_by_family=self.chk_family.isChecked(),
                        show_nodes=self.chk_nodes.isChecked(),
                        show_panels=self.chk_panels.isChecked(),
                        show_labels=self.chk_labels.isChecked())
        c.refresh()

    def redraw_loads(self):
        c = self.canvas_loads
        ax = c.clear()
        try:
            basic = ld.basic_load_vectors(self.model)
            name = self.loadcase_combo.currentText()
            combo = None
            for cc in self.model.combinations:
                if cc.name == name:
                    combo = cc
                    break
            if combo is None:
                combo = self.model.combinations[0]
            F = ld.combine(basic, combo)
            title = f"Applied nodal loads -- {combo.name}"
            if abs(getattr(combo, "T", 0.0)) > 0.0 and \
                    self.model.loads.has_self_straining:
                title += "\n(+ temperature / settlement, not drawn)"
            plot.draw_loads(ax, self.model, F, title=title)
        except Exception as exc:
            ax = c.clear()
            ax.text2D(0.5, 0.5, f"Loads cannot be drawn:\n{exc}",
                      transform=ax.transAxes, ha="center", va="center",
                      fontsize=9, color="#9d0208")
            self.log(f"Loads view: {exc}")
        c.refresh()

    def _auto_scale(self):
        if self.results is None:
            return
        name = self.defl_combo.currentText() or self.results.service_case
        r = self.results.cases.get(name)
        if r is None:
            return
        s = plot.auto_scale(self.model.geom, r.u)
        self.defl_scale.blockSignals(True)
        self.defl_scale.setValue(min(max(s, 1.0), 100000.0))
        self.defl_scale.blockSignals(False)
        self.redraw_deflected()

    def redraw_deflected(self):
        if self.results is None:
            return
        name = self.defl_combo.currentText()
        r = self.results.cases.get(name)
        if r is None:
            return
        c = self.canvas_defl
        ax = c.clear()
        art = plot.draw_deflected(ax, self.model, r,
                                  scale=self.defl_scale.value(),
                                  show_undeformed=self.chk_undeformed.isChecked())
        c.add_colorbar(art, "nodal displacement [mm]")
        c.refresh()

    def redraw_axial(self):
        if self.results is None:
            return
        name = self.axial_combo.currentText()
        r = self.results.cases.get(name)
        if r is None:
            return
        c = self.canvas_axial
        ax = c.clear()
        art = plot.draw_axial(ax, self.model, r.N,
                              title=f"Member axial force -- {name}",
                              width_scale=self.chk_width.isChecked())
        c.add_colorbar(art, "axial force [kN]   (+ tension)")
        c.refresh()

    def redraw_utilisation(self):
        if self.results is None:
            return
        c = self.canvas_util
        ax = c.clear()
        plot.draw_utilisation(ax, self.model, self.results.util_env)
        c.refresh()

    def redraw_buckling(self):
        if self.results is None:
            return
        c = self.canvas_buck
        ax = c.clear()
        b = self.results.buckling
        art = plot.draw_buckling_mode(ax, self.model, b,
                                      k=self.mode_spin.value() - 1,
                                      scale_frac=self.buck_scale.value())
        if b is not None and len(b.load_factors):
            c.add_colorbar(art, "normalised mode amplitude")
            self.buck_label.setText(
                "lambda_cr = " + ", ".join(f"{x:.2f}"
                                           for x in b.load_factors[:6]))
        else:
            self.buck_label.setText("")
        c.refresh()

    def redraw_charts(self):
        cc = self.canvas_charts
        plot.draw_strut_histogram(cc.ax_top, self.model)
        if self.results is not None:
            plot.draw_utilisation_histogram(cc.ax_bottom, self.results)
        cc.refresh()

    # ============================================================= file I/O ==
    def _default_file_name(self):
        g = self.model.geom
        if g.is_custom:
            stem = "".join(ch if ch.isalnum() else "_"
                           for ch in (g.title or "truss").lower())
            stem = "_".join(p for p in stem.split("_") if p)[:40] or "truss"
            return f"{stem}.json"
        return "dome_model.json"

    def on_save(self):
        if self._busy_message():
            return
        if not self._read_inputs():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save model", self._default_file_name(),
            "GeoDomeLab model (*.json)")
        if not path:
            return
        if not path.lower().endswith(".json"):
            path += ".json"
        try:
            self.model.save(path)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(
                self, "Save failed", f"The model could not be saved:\n{exc}")
            return
        self.log(f"Model saved to {os.path.basename(path)}")
        self.status("Model saved.")

    def on_open(self):
        if self._busy_message():
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open model", "",
            "GeoDomeLab model (*.json);;All files (*)")
        if not path:
            return
        gp = self.geometry_panel
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            problem = _model_file_problem(
                data, (gp.frequency.minimum(), gp.frequency.maximum()),
                (gp.radius.minimum(), gp.radius.maximum()))
            if problem:
                raise ValueError(problem)
            model = DomeModel.from_dict(data)
            _check_catalogue_names(model)
        except json.JSONDecodeError:
            text = "This file is not a valid model file (JSON)."
        except UnicodeDecodeError:
            text = "This file is not a text model file (JSON)."
        except KeyError as exc:
            text = (f"This model file is not complete: the "
                    f"\"{exc.args[0] if exc.args else exc}\" entry is "
                    f"missing.")
        except OSError as exc:
            text = f"The file could not be read: {exc.strerror or exc}."
        except Exception as exc:
            text = plain_message(f"{type(exc).__name__}: {exc}")
        else:
            text = None
        if text is not None:
            self.log(f"Open failed: {os.path.basename(path)}: {text}")
            QtWidgets.QMessageBox.critical(
                self, "Open failed",
                f"{os.path.basename(path)} could not be opened:\n{text}")
            return
        # The geometry is rebuilt from the file as it is: a custom truss is
        # never regenerated as a dome.
        name = os.path.basename(path)
        self.set_model(model, title=f"{_model_title(model)}  [{name}]",
                       message=f"Model loaded from {name}")
        for note in model.geom.notes:
            if model.geom.is_custom:
                self.log("  " + note)
        self.view_tabs.setCurrentWidget(self.canvas_model)
        self.status("Model loaded.")

    def on_import_csv(self):
        if self._busy_message():
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Import truss from CSV: choose the node file (nodes.csv)",
            "", "Node file (*nodes.csv);;CSV files (*.csv);;All files (*)")
        if not path:
            return
        try:
            model, notes = csv_import.import_truss_csv(path)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self, "Import failed",
                f"The truss could not be imported:\n\n{exc}\n\n"
                "See Help > CSV import format.")
            self.log(f"CSV import failed: {exc}")
            return
        g = model.geom
        self.set_model(model, message=f"Imported {g.title}: {g.n_nodes} "
                                      f"nodes, {g.n_members} members.")
        for note in notes:
            self.log("  " + note)
        self.view_tabs.setCurrentWidget(self.canvas_model)
        self.status(f"Imported {g.n_nodes} nodes and {g.n_members} members. "
                    f"Check the Sections and Loads tabs, then press F9.")

    def _project_title(self):
        """Title of the HTML report (geometry title for a custom truss)."""
        g = self.model.geom
        if g.is_custom:
            return f"Space Truss Analysis -- {g.title or 'Custom truss'}"
        return (f"Geodesic Dome Analysis -- {g.frequency}V, "
                f"R = {g.radius:.1f} m")

    def _ask_folder(self, title):
        return QtWidgets.QFileDialog.getExistingDirectory(
            self, title, os.getcwd())

    def on_export_csv(self):
        if not self._require_results():
            return
        folder = self._ask_folder("Choose a folder for the CSV tables")
        if not folder:
            return
        try:
            files = rp.export_csv(self.results, folder)
        except Exception as exc:
            self._export_failed("The CSV tables", exc)
            return
        self.log("Wrote:\n  " + "\n  ".join(files))
        self.status(f"{len(files)} CSV/text files written.")

    def _export_failed(self, what, exc):
        self.log(f"ERROR: {what} could not be written: {exc}")
        self.log(traceback.format_exc())
        self.status("Export failed.")
        QtWidgets.QMessageBox.critical(
            self, "Export failed", f"{what} could not be written:\n{exc}")

    def on_export_figs(self):
        if self._busy_message():
            return
        if not self._require_results():
            return
        folder = self._ask_folder("Choose a folder for the PNG figures")
        if not folder:
            return
        self._set_busy(True, "Rendering figures...")
        try:
            files = rp.export_figures(self.results, folder)
        except Exception as exc:
            self._set_busy(False)
            self._export_failed("The figures", exc)
            return
        finally:
            if self._busy:
                self._set_busy(False)
        self.log("Wrote:\n  " + "\n  ".join(files))
        self.status(f"{len(files)} figures written.")

    def on_export_html(self):
        if self._busy_message():
            return
        if not self._require_results():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export HTML report", "dome_report.html", "HTML (*.html)")
        if not path:
            return
        self._set_busy(True, "Building report...")
        try:
            rp.export_html(self.results, path, sizing=self.sizing,
                           verification=self.verification_rows,
                           project_title=self._project_title(),
                           extras=self.extension_results())
        except Exception as exc:
            self._set_busy(False)
            self._export_failed("The HTML report", exc)
            return
        finally:
            if self._busy:
                self._set_busy(False)
        self.log(f"Report written to {path}")
        self.status("Report written.")
        if QtWidgets.QMessageBox.question(
                self, "Report", "Report written. Open it in your browser?") \
                == QtWidgets.QMessageBox.Yes:
            webbrowser.open("file:///" + os.path.abspath(path).replace("\\", "/"))

    def on_export_all(self):
        if self._busy_message():
            return
        if not self._require_results():
            return
        folder = self._ask_folder("Choose an output folder")
        if not folder:
            return
        self._set_busy(True, "Exporting everything...")
        try:
            files = rp.export_csv(self.results, folder)
            files += rp.export_figures(self.results, folder)
            report = os.path.join(folder, "dome_report.html")
            rp.export_html(self.results, report, sizing=self.sizing,
                           verification=self.verification_rows,
                           project_title=self._project_title(),
                           extras=self.extension_results())
            files.append(report)
            # the model the results were computed for
            self.results.model.save(os.path.join(folder, "dome_model.json"))
            files.append(os.path.join(folder, "dome_model.json"))
        except Exception as exc:
            self._set_busy(False)
            self._export_failed("The output files", exc)
            return
        finally:
            if self._busy:
                self._set_busy(False)
        self.log("Wrote:\n  " + "\n  ".join(files))
        self.status(f"{len(files)} files written to {folder}.")

    def on_save_view(self):
        w = self.view_tabs.currentWidget()
        if not hasattr(w, "save"):
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save view", "dome_view.png", "PNG (*.png)")
        if path:
            try:
                w.save(path)
            except Exception as exc:
                self._export_failed("The view", exc)
                return
            self.log(f"View saved to {path}")

    # ================================================================ help ==
    def about_text(self):
        """The About box text (HTML)."""
        return (
            f"<h3>{APP_NAME} {html.escape(__version__)}</h3>"
            f"<p>{APP_SUBTITLE}</p>"
            "<p>Class I icosahedral geodesic domes and general space trusses. "
            "Direct stiffness analysis with temperature, settlement and "
            "spring supports, AISC 360-16 LRFD member checks, linear "
            "buckling, geometrically nonlinear load paths, modal analysis, "
            "fabrication schedules, parameter studies and fully stressed "
            "auto-sizing.</p>"
            f"<p>Author: {html.escape(__author__)}<br>"
            f"{COURSE_LINE}<br>"
            f"<a href=\"{html.escape(__url__)}\">"
            f"{html.escape(__url__)}</a></p>"
            f"<p>{LICENSE_NOTE}</p>"
            "<p>Built on NumPy, SciPy, matplotlib and Qt (via qtpy).</p>")

    def on_about(self):
        QtWidgets.QMessageBox.about(self, f"About {APP_NAME}",
                                    self.about_text())

    def _text_dialog(self, title, text, size=(760, 680)):
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle(title)
        dlg.resize(*size)
        v = QtWidgets.QVBoxLayout(dlg)
        te = QtWidgets.QPlainTextEdit(text.strip())
        te.setReadOnly(True)
        v.addWidget(te)
        b = QtWidgets.QPushButton("Close")
        b.clicked.connect(dlg.accept)
        v.addWidget(b, 0, Qt.AlignRight)
        dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()

    def on_csv_help(self):
        self._text_dialog("CSV import format", CSV_FORMAT_HELP, (700, 640))

    def theory_dialog(self):
        """The Help > Theory window: :mod:`geodome.gui.theory` in a read-only
        text browser.  The contents links jump inside the page; web links
        open in the browser."""
        from .theory import theory_html
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("Theory and method")
        dlg.resize(820, 760)
        v = QtWidgets.QVBoxLayout(dlg)
        browser = QtWidgets.QTextBrowser()
        browser.setOpenExternalLinks(True)
        browser.setHtml(theory_html())
        v.addWidget(browser)
        b = QtWidgets.QPushButton("Close")
        b.clicked.connect(dlg.accept)
        v.addWidget(b, 0, Qt.AlignRight)
        dlg.browser = browser
        return dlg

    def on_theory(self):
        dlg = self.theory_dialog()
        dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()

    # ---------------------------------------------------------------- close --
    def closeEvent(self, event):
        if self._worker is not None and self._worker.isRunning():
            self._worker.wait(3000)
        QtWidgets.QApplication.restoreOverrideCursor()
        event.accept()


# ================================================================== launch ==
_WINDOW = None          # module-level reference keeps the window alive


def _safe_excepthook(etype, value, tb):
    """Report an unexpected error in the window instead of closing it."""
    text = "".join(traceback.format_exception(etype, value, tb))
    sys.stderr.write(text)
    win = _WINDOW
    if win is None:
        return
    try:
        win.log("UNEXPECTED ERROR (the program keeps running):\n" + text)
        win.status(f"Unexpected error: {etype.__name__}: {value}")
        if win._busy and (win._worker is None or not win._worker.isRunning()):
            win._set_busy(False)
        QtWidgets.QMessageBox.warning(
            win, "Unexpected error",
            f"{etype.__name__}: {value}\n\nThe details are in the Log tab.")
    except Exception:                                # pragma: no cover
        pass


def launch(show=True):
    """Create (or reuse) the QApplication and show the main window.

    Safe to call from a plain ``python`` run, from Spyder's IPython console
    with the Qt event loop already running, and from an interactive session.
    """
    app = QtWidgets.QApplication.instance()
    created = app is None
    if created:
        app = QtWidgets.QApplication(sys.argv)
        app.setApplicationName(APP_NAME)
    try:
        app.setStyle("Fusion")
    except Exception:
        pass

    st.apply_light_theme(app)

    # A plain ``python run_geodome.py`` run: PyQt aborts the whole program
    # on an exception inside a slot.  Log it and keep the window open
    # instead.  (Spyder and test runners install their own hook; leave it.)
    if sys.excepthook is sys.__excepthook__:
        sys.excepthook = _safe_excepthook

    global _WINDOW
    _WINDOW = MainWindow()
    if show:
        _WINDOW.show()
        _WINDOW.raise_()
        _WINDOW.activateWindow()

    if created:
        exec_fn = getattr(app, "exec", None) or getattr(app, "exec_")
        exec_fn()
    return _WINDOW
