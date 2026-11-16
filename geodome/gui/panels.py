"""The four input panels of the control dock."""

from __future__ import annotations

from qtpy import QtWidgets
from qtpy.QtCore import Qt, Signal

from .. import geometry as geo
from ..materials import MATERIALS, DEFAULT_MATERIAL
from ..model import LoadCombination, default_combinations
from ..sections import SECTION_NAMES, SECTIONS
from . import style as st


# ============================================================== geometry =====
class GeometryPanel(QtWidgets.QWidget):
    """Frequency, radius, truncation and support condition."""

    generate_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)

        box = QtWidgets.QGroupBox("Geodesic breakdown")
        form = QtWidgets.QFormLayout(box)
        form.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)

        self.frequency = st.spin(1, 10, 4, 1, "V")
        self.frequency.setToolTip(
            "Class I icosahedral frequency.  A full sphere has 10v^2 + 2 nodes "
            "and 30v^2 struts, so the model grows quickly: 2V is a hand-check "
            "model, 4V a realistic small dome, 8V a large one.")
        form.addRow("Frequency", self.frequency)

        self.radius = st.dspin(1.0, 200.0, 12.0, 0.5, 2, " m")
        form.addRow("Sphere radius", self.radius)

        frac_row = QtWidgets.QWidget()
        fr = QtWidgets.QHBoxLayout(frac_row)
        fr.setContentsMargins(0, 0, 0, 0)
        fr.setSpacing(6)
        self.fraction = st.dspin(0.25, 0.95, 0.50, 0.01, 2, "")
        self.fraction.setMaximumWidth(78)
        self.fraction_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.fraction_slider.setRange(25, 95)
        self.fraction_slider.setValue(50)
        fr.addWidget(self.fraction)
        fr.addWidget(self.fraction_slider, 1)
        form.addRow("Dome fraction", frac_row)
        self.fraction.valueChanged.connect(
            lambda v: self.fraction_slider.setValue(int(round(v * 100))))
        self.fraction_slider.valueChanged.connect(
            lambda v: self.fraction.setValue(v / 100.0))

        self.support = QtWidgets.QComboBox()
        for key, text in geo.SUPPORT_TYPES.items():
            self.support.addItem(text, key)
        form.addRow("Supports", self.support)
        lay.addWidget(box)

        lay.addWidget(st.hint(
            "0.50 is a hemisphere.  The requested cut is snapped onto the "
            "nearest ring of real nodes, so the base ring is always a closed "
            "polygon; odd frequencies have no equator ring and will snap "
            "noticeably."))

        self.generate_btn = QtWidgets.QPushButton("Generate dome")
        self.generate_btn.setObjectName("primary")
        self.generate_btn.clicked.connect(self.generate_requested.emit)
        lay.addWidget(self.generate_btn)

        info = QtWidgets.QGroupBox("Generated geometry")
        iv = QtWidgets.QVBoxLayout(info)
        self.info = QtWidgets.QPlainTextEdit()
        self.info.setReadOnly(True)
        self.info.setMaximumHeight(188)
        iv.addWidget(self.info)
        self.schedule = st.table(["Fam", "Length [m]", "Chord factor", "Count"])
        self.schedule.setMaximumHeight(188)
        iv.addWidget(self.schedule)
        lay.addWidget(info, 1)

    # ------------------------------------------------------------- binding --
    def read_into(self, model):
        model.support_type = self.support.currentData()

    def params(self):
        return (self.frequency.value(), self.radius.value(),
                self.fraction.value())

    def show_model(self, model):
        self.info.setPlainText("\n".join(model.summary_lines()))
        rows = [{"Fam": s["label"],
                 "Length [m]": round(s["length"], 4),
                 "Chord factor": round(s["chord_factor"], 6),
                 "Count": s["count"]}
                for s in model.geom.strut_schedule()]
        st.fill_table(self.schedule, rows, numeric_digits=6)


# =============================================================== sections ====
class SectionPanel(QtWidgets.QWidget):
    """Material choice and a section per strut family."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._combos = []
        self._block = False

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)

        box = QtWidgets.QGroupBox("Material")
        form = QtWidgets.QFormLayout(box)
        self.material = QtWidgets.QComboBox()
        self.material.addItems(list(MATERIALS.keys()))
        self.material.setCurrentText(DEFAULT_MATERIAL)
        self.material.currentTextChanged.connect(self._emit)
        form.addRow("Grade", self.material)
        self.material_info = QtWidgets.QLabel("")
        self.material_info.setObjectName("metric")
        self.material_info.setWordWrap(True)
        form.addRow("", self.material_info)
        self.material.currentTextChanged.connect(self._update_material_info)
        lay.addWidget(box)

        box2 = QtWidgets.QGroupBox("Section per strut family")
        v2 = QtWidgets.QVBoxLayout(box2)
        apply_row = QtWidgets.QHBoxLayout()
        self.apply_combo = QtWidgets.QComboBox()
        self.apply_combo.addItems(SECTION_NAMES)
        self.apply_combo.setCurrentText("CHS 88.9 x 4.0")
        btn = QtWidgets.QPushButton("Apply to all")
        btn.clicked.connect(self._apply_all)
        apply_row.addWidget(self.apply_combo, 1)
        apply_row.addWidget(btn)
        v2.addLayout(apply_row)

        self.family_table = QtWidgets.QTableWidget(0, 4)
        self.family_table.setHorizontalHeaderLabels(
            ["Family", "Length [m]", "Count", "Section"])
        self.family_table.verticalHeader().setVisible(False)
        self.family_table.setAlternatingRowColors(True)
        self.family_table.horizontalHeader().setSectionResizeMode(
            3, QtWidgets.QHeaderView.Stretch)
        v2.addWidget(self.family_table, 1)
        lay.addWidget(box2, 1)

        self.mass_label = QtWidgets.QLabel("")
        self.mass_label.setObjectName("metric")
        lay.addWidget(self.mass_label)

        lay.addWidget(st.hint(
            "Geodesic domes are fabricated from a small set of repeated strut "
            "lengths, so sections are assigned by family rather than by "
            "individual member.  That is also what the auto-sizing routine "
            "optimises."))
        self._update_material_info(self.material.currentText())

    # ------------------------------------------------------------ internal --
    def _emit(self, *_):
        if not self._block:
            self.changed.emit()

    def _update_material_info(self, name):
        self.material_info.setText(MATERIALS[name].describe())

    def _apply_all(self):
        name = self.apply_combo.currentText()
        self._block = True
        for c in self._combos:
            c.setCurrentText(name)
        self._block = False
        self.changed.emit()

    # ------------------------------------------------------------- binding --
    def rebuild(self, model):
        """Re-create one row per strut family after the geometry changes."""
        self._block = True
        sched = model.geom.strut_schedule()
        self._combos = []
        self.family_table.setRowCount(len(sched))
        for i, s in enumerate(sched):
            self.family_table.setItem(
                i, 0, QtWidgets.QTableWidgetItem(s["label"]))
            it = QtWidgets.QTableWidgetItem(f"{s['length']:.4f}")
            it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.family_table.setItem(i, 1, it)
            it = QtWidgets.QTableWidgetItem(str(s["count"]))
            it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.family_table.setItem(i, 2, it)

            combo = QtWidgets.QComboBox()
            combo.addItems(SECTION_NAMES)
            combo.setCurrentText(model.group_sections[s["id"]])
            combo.currentTextChanged.connect(self._emit)
            self.family_table.setCellWidget(i, 3, combo)
            self._combos.append(combo)
        self.family_table.resizeColumnsToContents()
        self.family_table.horizontalHeader().setSectionResizeMode(
            3, QtWidgets.QHeaderView.Stretch)
        self.material.setCurrentText(model.material_name)
        self._block = False
        self.update_mass(model)

    def sync_from_model(self, model):
        """Push section names from the model into the combo boxes."""
        self._block = True
        for i, c in enumerate(self._combos):
            if i < len(model.group_sections):
                c.setCurrentText(model.group_sections[i])
        self.material.setCurrentText(model.material_name)
        self._block = False
        self.update_mass(model)

    def read_into(self, model):
        model.material_name = self.material.currentText()
        if self._combos:
            model.group_sections = [c.currentText() for c in self._combos]

    def update_mass(self, model):
        self.mass_label.setText(
            f"Structural mass  {model.steel_mass():,.0f} kg     "
            f"({model.mass_per_plan_area():.1f} kg per m2 of plan area)")


# ================================================================== loads ====
class LoadPanel(QtWidgets.QWidget):
    """Gravity, snow and wind input plus the load combination table."""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        inner = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(inner)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)
        scroll.setWidget(inner)
        outer.addWidget(scroll)

        # ---- gravity -------------------------------------------------------
        box = QtWidgets.QGroupBox("Gravity loads")
        form = QtWidgets.QFormLayout(box)
        self.self_weight = QtWidgets.QCheckBox("Include strut self weight")
        self.self_weight.setChecked(True)
        form.addRow("", self.self_weight)
        self.dead = st.dspin(0.0, 20.0, 0.30, 0.05, 3, " kPa")
        self.dead.setToolTip("Cladding, glazing, services.  Applied on the "
                             "faceted surface area.")
        form.addRow("Superimposed dead", self.dead)
        self.snow = st.dspin(0.0, 20.0, 1.20, 0.05, 3, " kPa")
        self.snow.setToolTip("Applied on the horizontally projected area.")
        form.addRow("Snow", self.snow)
        self.snow_unbalanced = QtWidgets.QCheckBox(
            "Unbalanced snow (1.5 / 0.5 split)")
        form.addRow("", self.snow_unbalanced)
        self.live = st.dspin(0.0, 20.0, 0.0, 0.05, 3, " kPa")
        form.addRow("Live (projected)", self.live)
        self.crown_load = st.dspin(0.0, 5000.0, 0.0, 1.0, 2, " kN")
        self.crown_load.setToolTip("Downward point load at the crown node, "
                                   "e.g. a hung service or a test load.")
        form.addRow("Crown point load", self.crown_load)
        lay.addWidget(box)

        # ---- wind ----------------------------------------------------------
        self.wind_box = QtWidgets.QGroupBox("Wind")
        self.wind_box.setCheckable(True)
        self.wind_box.setChecked(True)
        wf = QtWidgets.QFormLayout(self.wind_box)
        self.wind_speed = st.dspin(5.0, 110.0, 45.0, 1.0, 1, " m/s")
        self.wind_speed.setToolTip("3-second gust basic wind speed.")
        wf.addRow("Basic wind speed V", self.wind_speed)
        self.wind_Kz = st.dspin(0.3, 2.0, 0.90, 0.05, 2)
        wf.addRow("Exposure Kz", self.wind_Kz)
        self.wind_Kzt = st.dspin(0.5, 3.0, 1.00, 0.05, 2)
        wf.addRow("Topographic Kzt", self.wind_Kzt)
        self.wind_Kd = st.dspin(0.5, 1.0, 0.95, 0.05, 2)
        wf.addRow("Directionality Kd", self.wind_Kd)
        self.wind_G = st.dspin(0.5, 1.5, 0.85, 0.05, 2)
        wf.addRow("Gust factor G", self.wind_G)
        self.wind_az = st.dspin(0.0, 360.0, 0.0, 15.0, 1, " deg")
        wf.addRow("Wind direction", self.wind_az)
        self.qz_label = QtWidgets.QLabel("")
        self.qz_label.setObjectName("metric")
        wf.addRow("", self.qz_label)
        lay.addWidget(self.wind_box)
        lay.addWidget(st.hint(
            "qz = 0.613 Kz Kzt Kd V^2.  The pressure coefficient varies along "
            "the meridian from +0.80 at the windward springing through -1.00 "
            "at the crown to -0.40 leeward, applied normal to the surface."))

        # ---- combinations --------------------------------------------------
        box3 = QtWidgets.QGroupBox("Load combinations")
        v3 = QtWidgets.QVBoxLayout(box3)
        self.combos = QtWidgets.QTableWidget(0, 6)
        self.combos.setHorizontalHeaderLabels(
            ["Name", "D", "S", "L", "W", "Service"])
        self.combos.verticalHeader().setVisible(False)
        self.combos.setAlternatingRowColors(True)
        self.combos.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.Stretch)
        self.combos.setMinimumHeight(165)
        v3.addWidget(self.combos)
        row = QtWidgets.QHBoxLayout()
        for text, slot in (("Add", self._add_row),
                           ("Remove", self._remove_row),
                           ("Reset to ASCE 7", self._reset_rows)):
            b = QtWidgets.QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        row.addStretch(1)
        v3.addLayout(row)
        lay.addWidget(box3)
        lay.addWidget(st.hint(
            "Factors are editable.  Tick 'Service' to use a row for the "
            "deflection check instead of the strength envelope."))
        lay.addStretch(1)

        for w in (self.wind_speed, self.wind_Kz, self.wind_Kzt, self.wind_Kd):
            w.valueChanged.connect(self._update_qz)
        self.wind_box.toggled.connect(self._update_qz)
        self._reset_rows()
        self._update_qz()

        for w in (self.dead, self.snow, self.live, self.crown_load,
                  self.wind_speed, self.wind_az):
            w.valueChanged.connect(self.changed.emit)
        self.self_weight.toggled.connect(self.changed.emit)
        self.snow_unbalanced.toggled.connect(self.changed.emit)
        self.wind_box.toggled.connect(self.changed.emit)

    # ------------------------------------------------------------ internal --
    def _update_qz(self, *_):
        q = (0.613 * self.wind_Kz.value() * self.wind_Kzt.value()
             * self.wind_Kd.value() * self.wind_speed.value() ** 2)
        self.qz_label.setText(f"qz = {q / 1e3:.3f} kPa    "
                              f"qz*G = {q * self.wind_G.value() / 1e3:.3f} kPa")

    def _set_row(self, i, combo):
        self.combos.setItem(i, 0, QtWidgets.QTableWidgetItem(combo.name))
        for j, key in enumerate(("D", "S", "L", "W"), start=1):
            it = QtWidgets.QTableWidgetItem(f"{getattr(combo, key):g}")
            it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.combos.setItem(i, j, it)
        chk = QtWidgets.QTableWidgetItem()
        chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
        chk.setCheckState(Qt.Checked if combo.serviceability else Qt.Unchecked)
        self.combos.setItem(i, 5, chk)

    def _add_row(self):
        i = self.combos.rowCount()
        self.combos.insertRow(i)
        self._set_row(i, LoadCombination(f"LC{i + 1}", D=1.0))

    def _remove_row(self):
        i = self.combos.currentRow()
        if i >= 0 and self.combos.rowCount() > 1:
            self.combos.removeRow(i)

    def _reset_rows(self):
        rows = default_combinations()
        self.combos.setRowCount(len(rows))
        for i, c in enumerate(rows):
            self._set_row(i, c)
        self.combos.resizeColumnsToContents()
        self.combos.horizontalHeader().setSectionResizeMode(
            0, QtWidgets.QHeaderView.Stretch)

    # ------------------------------------------------------------- binding --
    def read_into(self, model):
        ls = model.loads
        ls.include_self_weight = self.self_weight.isChecked()
        ls.dead_pressure = self.dead.value() * 1e3
        ls.snow_pressure = self.snow.value() * 1e3
        ls.snow_unbalanced = self.snow_unbalanced.isChecked()
        ls.live_pressure = self.live.value() * 1e3
        ls.crown_point_load = self.crown_load.value() * 1e3
        ls.wind_enabled = self.wind_box.isChecked()
        ls.wind_speed = self.wind_speed.value()
        ls.wind_Kz = self.wind_Kz.value()
        ls.wind_Kzt = self.wind_Kzt.value()
        ls.wind_Kd = self.wind_Kd.value()
        ls.wind_G = self.wind_G.value()
        ls.wind_azimuth = self.wind_az.value()

        # The table has no T column yet: keep the self-straining factor T of
        # the combination already in the model at the same row and name, so a
        # loaded v2 model does not silently lose its temperature cases.
        old = list(getattr(model, "combinations", []) or [])
        combos = []
        for i in range(self.combos.rowCount()):
            def val(j):
                it = self.combos.item(i, j)
                try:
                    return float(it.text()) if it and it.text().strip() else 0.0
                except ValueError:
                    return 0.0
            name_item = self.combos.item(i, 0)
            name = (name_item.text().strip() if name_item else "") or f"LC{i + 1}"
            chk = self.combos.item(i, 5)
            T = (float(getattr(old[i], "T", 0.0))
                 if i < len(old) and old[i].name == name else 0.0)
            combos.append(LoadCombination(
                name, D=val(1), S=val(2), L=val(3), W=val(4),
                serviceability=bool(chk and chk.checkState() == Qt.Checked),
                T=T))
        if combos:
            model.combinations = combos

    def sync_from_model(self, model):
        ls = model.loads
        self.self_weight.setChecked(ls.include_self_weight)
        self.dead.setValue(ls.dead_pressure / 1e3)
        self.snow.setValue(ls.snow_pressure / 1e3)
        self.snow_unbalanced.setChecked(ls.snow_unbalanced)
        self.live.setValue(ls.live_pressure / 1e3)
        self.crown_load.setValue(ls.crown_point_load / 1e3)
        self.wind_box.setChecked(ls.wind_enabled)
        self.wind_speed.setValue(ls.wind_speed)
        self.wind_Kz.setValue(ls.wind_Kz)
        self.wind_Kzt.setValue(ls.wind_Kzt)
        self.wind_Kd.setValue(ls.wind_Kd)
        self.wind_G.setValue(ls.wind_G)
        self.wind_az.setValue(ls.wind_azimuth)
        self.combos.setRowCount(len(model.combinations))
        for i, c in enumerate(model.combinations):
            self._set_row(i, c)
        self._update_qz()


# ================================================================= design ====
class DesignPanel(QtWidgets.QWidget):
    """Code parameters, analysis options and the action buttons."""

    analyse_requested = Signal()
    size_requested = Signal()
    verify_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)

        box = QtWidgets.QGroupBox("Design parameters (AISC 360 LRFD)")
        form = QtWidgets.QFormLayout(box)
        self.K = st.dspin(0.5, 2.0, 1.0, 0.05, 2)
        self.K.setToolTip("Effective length factor.  K = 1.0 suits the pinned "
                          "node connections normally used in a single-layer "
                          "geodesic dome.")
        form.addRow("Effective length K", self.K)
        self.defl_limit = st.spin(60, 1000, 240, 10)
        self.defl_limit.setPrefix("span / ")
        form.addRow("Deflection limit", self.defl_limit)
        self.phi_c = st.dspin(0.5, 1.0, 0.90, 0.05, 2)
        form.addRow("phi compression", self.phi_c)
        self.phi_t = st.dspin(0.5, 1.0, 0.90, 0.05, 2)
        form.addRow("phi tension", self.phi_t)
        lay.addWidget(box)

        box2 = QtWidgets.QGroupBox("Stability analysis")
        f2 = QtWidgets.QFormLayout(box2)
        self.do_buckling = QtWidgets.QCheckBox(
            "Run linear buckling eigenanalysis")
        self.do_buckling.setChecked(True)
        f2.addRow("", self.do_buckling)
        self.n_modes = st.spin(1, 20, 6, 1)
        f2.addRow("Modes to extract", self.n_modes)
        lay.addWidget(box2)
        lay.addWidget(st.hint(
            "Solves [K + lambda Kg]{phi} = 0.  Single-layer domes are "
            "imperfection sensitive: the real limit load can be well below "
            "lambda_cr; quantify it with the nonlinear imperfection study."))

        box3 = QtWidgets.QGroupBox("Auto-sizing (fully stressed design)")
        f3 = QtWidgets.QFormLayout(box3)
        self.target_util = st.dspin(0.30, 1.00, 0.95, 0.05, 2)
        f3.addRow("Target utilisation", self.target_util)
        self.max_iter = st.spin(2, 40, 12, 1)
        f3.addRow("Max iterations", self.max_iter)
        self.slender_limit = st.dspin(60.0, 400.0, 200.0, 10.0, 0)
        self.slender_limit.setPrefix("KL/r <= ")
        f3.addRow("Slenderness cap", self.slender_limit)
        self.enforce_defl = QtWidgets.QCheckBox(
            "Enforce the deflection limit as well")
        self.enforce_defl.setChecked(True)
        f3.addRow("", self.enforce_defl)
        lay.addWidget(box3)

        self.analyse_btn = QtWidgets.QPushButton("Run analysis")
        self.analyse_btn.setObjectName("primary")
        self.analyse_btn.clicked.connect(self.analyse_requested.emit)
        lay.addWidget(self.analyse_btn)

        self.size_btn = QtWidgets.QPushButton("Auto-size the strut families")
        self.size_btn.clicked.connect(self.size_requested.emit)
        lay.addWidget(self.size_btn)

        self.verify_btn = QtWidgets.QPushButton("Run verification suite")
        self.verify_btn.clicked.connect(self.verify_requested.emit)
        lay.addWidget(self.verify_btn)
        lay.addStretch(1)

    def read_into(self, model):
        model.K_effective = self.K.value()
        model.deflection_limit = float(self.defl_limit.value())
        model.phi_c = self.phi_c.value()
        model.phi_t = self.phi_t.value()

    def sync_from_model(self, model):
        self.K.setValue(model.K_effective)
        self.defl_limit.setValue(int(model.deflection_limit))
        self.phi_c.setValue(model.phi_c)
        self.phi_t.setValue(model.phi_t)
