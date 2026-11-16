"""The GeoDomeLab main window."""

from __future__ import annotations

import os
import sys
import traceback
import webbrowser

import numpy as np
from qtpy import QtWidgets, QtGui
from qtpy.QtCore import Qt

from .. import __version__
from .. import geometry as geo
from .. import loads as ld
from .. import plotting as plot
from .. import report as rp
from .. import verification as vf
from ..model import DomeModel, default_model
from ..optimize import auto_size
from ..postprocess import run_analysis, summary_text, member_table, \
    node_table, reaction_table
from .canvas import Canvas3D, DualCanvas2D
from .panels import GeometryPanel, SectionPanel, LoadPanel, DesignPanel
from .worker import Worker
from . import style as st

APP_NAME = "GeoDomeLab"
APP_SUBTITLE = "Geodesic Dome Space-Truss Analyzer"


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"{APP_NAME} {__version__}  --  {APP_SUBTITLE}")
        self.resize(1480, 920)
        self.setStyleSheet(st.STYLESHEET)

        self.model = default_model(4, 12.0, 0.5)
        self.results = None
        self.sizing = None
        self.verification_rows = None
        self._worker = None
        self._busy = False

        self._build_ui()
        self._build_menu()
        self._build_toolbar()

        self.geometry_panel.show_model(self.model)
        self.section_panel.rebuild(self.model)
        self.load_panel.sync_from_model(self.model)
        self.design_panel.sync_from_model(self.model)
        self.redraw_model()
        self.log(f"{APP_NAME} {__version__} ready. "
                 f"A default 4V dome of 12 m radius has been generated.")
        self.status("Ready.")

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
        self.section_panel.changed.connect(self.on_section_changed)
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
        m.addAction(self._act("&New model", self.on_new, "Ctrl+N"))
        m.addAction(self._act("&Open model...", self.on_open, "Ctrl+O"))
        m.addAction(self._act("&Save model...", self.on_save, "Ctrl+S"))
        m.addSeparator()
        m.addAction(self._act("Export results to &CSV...", self.on_export_csv))
        m.addAction(self._act("Export &figures (PNG)...", self.on_export_figs))
        m.addAction(self._act("Export &HTML report...", self.on_export_html,
                              "Ctrl+R"))
        m.addAction(self._act("Export &everything to a folder...",
                              self.on_export_all))
        m.addSeparator()
        m.addAction(self._act("E&xit", self.close, "Ctrl+Q"))

        m = mb.addMenu("&Analysis")
        m.addAction(self._act("&Generate dome", self.on_generate, "F5"))
        m.addAction(self._act("&Run analysis", self.on_analyse, "F9"))
        m.addAction(self._act("&Auto-size strut families", self.on_auto_size,
                              "F10"))
        m.addSeparator()
        m.addAction(self._act("Run &verification suite", self.on_verify, "F11"))

        m = mb.addMenu("&View")
        for i in range(self.view_tabs.count()):
            name = self.view_tabs.tabText(i)
            m.addAction(self._act(name,
                                  lambda _=False, k=i:
                                  self.view_tabs.setCurrentIndex(k)))
        m.addSeparator()
        m.addAction(self._act("Save current view as PNG...", self.on_save_view))

        m = mb.addMenu("&Help")
        m.addAction(self._act("&Theory and method", self.on_theory))
        m.addAction(self._act("&About", self.on_about))

    def _build_toolbar(self):
        tb = QtWidgets.QToolBar("Main")
        tb.setToolButtonStyle(Qt.ToolButtonTextOnly)
        tb.setMovable(False)
        self.addToolBar(tb)
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
                  self.act_size, self.act_verify, self.act_report):
            w.setEnabled(not busy)
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
        QtWidgets.QMessageBox.critical(self, "Analysis failed", message)

    def _read_inputs(self):
        """Pull every panel's values into the model."""
        self.geometry_panel.read_into(self.model)
        self.section_panel.read_into(self.model)
        self.load_panel.read_into(self.model)
        self.design_panel.read_into(self.model)

    def _require_results(self):
        if self.results is None:
            QtWidgets.QMessageBox.information(
                self, APP_NAME, "Run the analysis first (F9).")
            return False
        return True

    # ============================================================= actions ==
    def on_new(self):
        self.model = default_model(4, 12.0, 0.5)
        self.results = None
        self.sizing = None
        self.geometry_panel.frequency.setValue(4)
        self.geometry_panel.radius.setValue(12.0)
        self.geometry_panel.fraction.setValue(0.5)
        self.geometry_panel.show_model(self.model)
        self.section_panel.rebuild(self.model)
        self.load_panel.sync_from_model(self.model)
        self.design_panel.sync_from_model(self.model)
        self._clear_results()
        self.redraw_model()
        self.log("New model created.")

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
        keep = list(self.model.group_sections)
        self.model = DomeModel(geom=g, material_name=self.model.material_name,
                               group_sections=keep,
                               support_type=self.model.support_type,
                               loads=self.model.loads,
                               combinations=self.model.combinations)
        # Rebuild the section rows before reading the panels back: the old
        # rows belong to the previous geometry and may be fewer than the new
        # strut families (e.g. 4V has 6, 6V has 9).
        self.section_panel.rebuild(self.model)
        self._read_inputs()
        self.results = None
        self.sizing = None
        self.geometry_panel.show_model(self.model)
        self._clear_results()
        self.redraw_model()
        self.redraw_loads()
        self.canvas_charts.ax_top.clear()
        plot.draw_strut_histogram(self.canvas_charts.ax_top, self.model)
        self.canvas_charts.refresh()
        self.view_tabs.setCurrentWidget(self.canvas_model)
        for note in g.notes:
            self.log(note)
        self.status(f"Generated {g.frequency}V dome: {g.n_nodes} nodes, "
                    f"{g.n_members} struts.")

    def on_section_changed(self):
        self.section_panel.read_into(self.model)
        self.section_panel.update_mass(self.model)
        self.geometry_panel.show_model(self.model)
        if self.chk_family.isChecked():
            self.redraw_model()

    def on_analyse(self):
        if self._busy:
            return
        self._read_inputs()
        self._set_busy(True, "Running analysis...")
        self._worker = Worker(
            run_analysis, self.model,
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
        self.status(f"Analysis complete. Max utilisation "
                    f"{res.governing_util:.3f} ({flag}) under "
                    f"{res.governing_case}.")

    def on_auto_size(self):
        if self._busy:
            return
        self._read_inputs()
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
        plot.draw_loads(ax, self.model, F,
                        title=f"Applied nodal loads -- {combo.name}")
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
    def on_save(self):
        self._read_inputs()
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save model", "dome_model.json", "JSON (*.json)")
        if not path:
            return
        self.model.save(path)
        self.log(f"Model saved to {path}")
        self.status("Model saved.")

    def on_open(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open model", "", "JSON (*.json)")
        if not path:
            return
        try:
            self.model = DomeModel.load(path)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Open failed", str(exc))
            return
        g = self.model.geom
        self.geometry_panel.frequency.setValue(g.frequency)
        self.geometry_panel.radius.setValue(g.radius)
        self.geometry_panel.fraction.setValue(g.dome_fraction)
        idx = self.geometry_panel.support.findData(self.model.support_type)
        if idx >= 0:
            self.geometry_panel.support.setCurrentIndex(idx)
        self.geometry_panel.show_model(self.model)
        self.section_panel.rebuild(self.model)
        self.load_panel.sync_from_model(self.model)
        self.design_panel.sync_from_model(self.model)
        self.results = None
        self._clear_results()
        self.redraw_model()
        self.log(f"Model loaded from {path}")
        self.status("Model loaded.")

    def _ask_folder(self, title):
        return QtWidgets.QFileDialog.getExistingDirectory(
            self, title, os.getcwd())

    def on_export_csv(self):
        if not self._require_results():
            return
        folder = self._ask_folder("Choose a folder for the CSV tables")
        if not folder:
            return
        files = rp.export_csv(self.results, folder)
        self.log("Wrote:\n  " + "\n  ".join(files))
        self.status(f"{len(files)} CSV/text files written.")

    def on_export_figs(self):
        if not self._require_results():
            return
        folder = self._ask_folder("Choose a folder for the PNG figures")
        if not folder:
            return
        self._set_busy(True, "Rendering figures...")
        try:
            files = rp.export_figures(self.results, folder)
        finally:
            self._set_busy(False)
        self.log("Wrote:\n  " + "\n  ".join(files))
        self.status(f"{len(files)} figures written.")

    def on_export_html(self):
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
                           project_title=f"Geodesic Dome Analysis -- "
                                         f"{self.model.geom.frequency}V, "
                                         f"R = {self.model.geom.radius:.1f} m")
        finally:
            self._set_busy(False)
        self.log(f"Report written to {path}")
        self.status("Report written.")
        if QtWidgets.QMessageBox.question(
                self, "Report", "Report written. Open it in your browser?") \
                == QtWidgets.QMessageBox.Yes:
            webbrowser.open("file:///" + os.path.abspath(path).replace("\\", "/"))

    def on_export_all(self):
        if not self._require_results():
            return
        folder = self._ask_folder("Choose an output folder")
        if not folder:
            return
        self._set_busy(True, "Exporting everything...")
        try:
            files = rp.export_csv(self.results, folder)
            files += rp.export_figures(self.results, folder)
            html = os.path.join(folder, "dome_report.html")
            rp.export_html(self.results, html, sizing=self.sizing,
                           verification=self.verification_rows)
            files.append(html)
            self.model.save(os.path.join(folder, "dome_model.json"))
            files.append(os.path.join(folder, "dome_model.json"))
        finally:
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
            w.save(path)
            self.log(f"View saved to {path}")

    # ================================================================ help ==
    def on_about(self):
        QtWidgets.QMessageBox.about(
            self, f"About {APP_NAME}",
            f"<h3>{APP_NAME} {__version__}</h3>"
            f"<p>{APP_SUBTITLE}</p>"
            "<p>Class I icosahedral geodesic dome generation, three-dimensional "
            "space-truss analysis by the direct stiffness method, AISC 360 "
            "LRFD member checks, linear buckling eigenanalysis and "
            "fully-stressed auto-sizing.</p>"
            "<p>Built on NumPy, SciPy, matplotlib and Qt (via qtpy).</p>")

    def on_theory(self):
        text = """
GeoDomeLab -- method summary
============================

1  GEOMETRY
   Class I ("alternate") icosahedral breakdown.  The 20 faces of a regular
   icosahedron are each divided into v^2 triangles by barycentric
   interpolation and every new point is projected radially onto the sphere.
   A complete sphere therefore has 10v^2 + 2 nodes, 30v^2 struts and 20v^2
   panels, which the verification suite checks together with Euler's formula
   V - E + F = 2.  The dome is truncated on a horizontal plane snapped to an
   existing ring of nodes.

2  ELEMENT
   Two-node pin-ended bar, three translational DOF per node.  With direction
   cosines {c},

              EA | [c c^T]  -[c c^T] |
       [k] =  -- |                   |
              L  |-[c c^T]   [c c^T] |

3  SOLUTION
   Sparse assembly, partition into free (f) and supported (s) DOF,

       [Kff]{uf} = {Ff}          reactions  {R} = [Ksf]{uf} - {Fs}
       N = (EA/L) {c}^T ({u2} - {u1})          tension positive

4  LOADS
   Dead on the faceted surface area, snow and live on the projected area,
   both lumped to the nodes through one third of each adjacent panel.  Wind
   uses qz = 0.613 Kz Kzt Kd V^2 with a meridional Cp distribution applied
   normal to the surface.  ASCE 7 LRFD combinations are enveloped.

5  MEMBER CHECKS  (AISC 360)
       tension      phi_t Pn = 0.90 Fy Ag
       compression  Fe  = pi^2 E / (KL/r)^2
                    Fcr = 0.658^(Fy/Fe) Fy    if KL/r <= 4.71 sqrt(E/Fy)
                    Fcr = 0.877 Fe            otherwise
                    phi_c Pn = 0.90 Fcr Ag
   Round HSS with D/t > 0.11E/Fy get the E7 slender-element factor Q.

6  GLOBAL BUCKLING
   Geometric stiffness of a bar carrying axial force N:

              N  | [I - c c^T]  -[I - c c^T] |
       [kg] = -  |                           |
              L  |-[I - c c^T]   [I - c c^T] |

   Eigenproblem [K + lambda Kg]{phi} = 0, recast as (-Kg){phi} =
   (1/lambda) K {phi} so the largest eigenvalue of a symmetric pair with a
   positive-definite right-hand side gives the smallest critical load factor.

7  AUTO-SIZING
   Fully stressed design on the strut families: analyse, pick the lightest
   catalogue section covering the envelope demands at the target utilisation
   subject to KL/r <= 200, repeat; step every family up if the deflection
   limit is not met.

8  VERIFICATION
   Closed-form bar and tripod solutions, published chord-factor tables,
   Euler's polyhedron formula, method-of-sections free bodies, membrane shell
   theory (n_phi = pR/2), five-fold symmetry, hoop reversal, strain-energy and
   equilibrium identities, and the linearity of the buckling eigenproblem.
   Run it from the Analysis menu (F11).
"""
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("Theory and method")
        dlg.resize(760, 680)
        v = QtWidgets.QVBoxLayout(dlg)
        te = QtWidgets.QPlainTextEdit(text.strip())
        te.setReadOnly(True)
        v.addWidget(te)
        b = QtWidgets.QPushButton("Close")
        b.clicked.connect(dlg.accept)
        v.addWidget(b, 0, Qt.AlignRight)
        dlg.exec() if hasattr(dlg, "exec") else dlg.exec_()

    # ---------------------------------------------------------------- close --
    def closeEvent(self, event):
        if self._worker is not None and self._worker.isRunning():
            self._worker.wait(3000)
        QtWidgets.QApplication.restoreOverrideCursor()
        event.accept()


# ================================================================== launch ==
_WINDOW = None          # module-level reference keeps the window alive


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
