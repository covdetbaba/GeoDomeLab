"""The four input panels of the control dock.

Every panel has ``sync_from_model(model)`` (model -> widgets) and
``read_into(model)`` (widgets -> model).  ``read_into`` raises
:class:`InputError` with a plain message when an input is not valid; the
main window shows that message and does not run anything.
"""

from __future__ import annotations

import html
import math

from qtpy import QtWidgets
from qtpy.QtCore import Qt, Signal

from .. import geometry as geo
from ..materials import MATERIALS, DEFAULT_MATERIAL
from ..model import LoadCombination, default_combinations
from ..sections import SECTION_NAMES, SECTIONS
from . import style as st

KN = 1.0e3          # kN -> N, kN/m -> N/m
MM = 1.0e-3         # mm -> m
LOAD_CASES = ("D", "S", "L", "W")


class InputError(ValueError):
    """A user input is not valid.  The message is shown as it is."""


def _scroll_wrap(widget, inner_layout_owner):
    """Put ``inner_layout_owner`` in a frameless scroll area inside widget."""
    outer = QtWidgets.QVBoxLayout(widget)
    outer.setContentsMargins(0, 0, 0, 0)
    scroll = QtWidgets.QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
    scroll.setWidget(inner_layout_owner)
    outer.addWidget(scroll)
    return scroll


def _spin_value(spin, old, scale=1.0):
    """``spin.value() * scale``, or ``old`` when the spin box only shows
    ``old`` rounded (so a value from a model file reads back exactly)."""
    v = spin.value() * scale
    try:
        old = float(old)
    except (TypeError, ValueError):
        return v
    tol = 0.51 * 10.0 ** (-spin.decimals()) * abs(scale)
    return old if math.isfinite(old) and abs(v - old) <= tol else v


def _support_text(c):
    """One custom-support entry as text: fixed / free / spring in kN/m."""
    if isinstance(c, bool) or c is None:
        return "fixed" if c else "free"
    if isinstance(c, str):
        return c
    k = float(c)
    if math.isinf(k):
        return "fixed"
    if k <= 0.0:
        return "free"
    return f"{k / KN:,.6g} kN/m"


def custom_support_counts(custom_supports):
    """(nodes, fixed DOFs, spring DOFs, free DOFs) of a custom support dict."""
    nf = ns = nfree = 0
    for conds in (custom_supports or {}).values():
        for c in conds:
            t = _support_text(c)
            if t == "fixed":
                nf += 1
            elif t == "free":
                nfree += 1
            else:
                ns += 1
    return len(custom_supports or {}), nf, ns, nfree


# ============================================================== geometry =====
class GeometryPanel(QtWidgets.QWidget):
    """Frequency, radius, truncation and support condition."""

    generate_requested = Signal()
    # the support type or a spring stiffness was changed by the user
    supports_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._syncing = False
        inner = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(inner)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)
        _scroll_wrap(self, inner)
        self._custom = False

        # ---- custom-model notice (hidden for a geodesic dome) -------------
        self.custom_note = st.note("")
        self.custom_note.setTextFormat(Qt.RichText)
        self.custom_note.setVisible(False)
        lay.addWidget(self.custom_note)

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
            self.support.setItemData(self.support.count() - 1, text,
                                     Qt.ToolTipRole)
        # The labels are long: let the box shrink to the panel width and
        # give the drop-down list room for the full text.
        self.support.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.support.setMinimumContentsLength(16)
        self.support.view().setMinimumWidth(
            self.support.view().sizeHintForColumn(0) + 28)
        self.support.setToolTip(
            "How the structure is held.  Elastic springs model a flexible "
            "foundation.  Custom per-node supports come from a model file "
            "or a CSV import.")
        form.addRow("Supports", self.support)
        self._form = form
        lay.addWidget(box)

        lay.addWidget(st.hint(
            "0.50 is a hemisphere.  The requested cut is snapped onto the "
            "nearest ring of real nodes, so the base ring is always a closed "
            "polygon; odd frequencies have no equator ring and will snap "
            "noticeably."))

        # ---- elastic supports (shown only when selected) ------------------
        self.spring_box = QtWidgets.QGroupBox("Elastic supports")
        sf = QtWidgets.QFormLayout(self.spring_box)
        sf.setLabelAlignment(Qt.AlignRight | Qt.AlignVCenter)
        k0 = [k / KN for k in geo.DEFAULT_SPRING_STIFFNESS]
        self.spring_kx = st.dspin(0.0, 1.0e9, k0[0], 1000.0, 1, " kN/m")
        self.spring_ky = st.dspin(0.0, 1.0e9, k0[1], 1000.0, 1, " kN/m")
        self.spring_kz = st.dspin(0.0, 1.0e9, k0[2], 1000.0, 1, " kN/m")
        for w, axis in ((self.spring_kx, "x"), (self.spring_ky, "y"),
                        (self.spring_kz, "z")):
            w.setToolTip(f"Spring stiffness in {axis} at every base node, "
                         f"in kN/m.  0 means free in {axis}.")
        sf.addRow("kx", self.spring_kx)
        sf.addRow("ky", self.spring_ky)
        sf.addRow("kz", self.spring_kz)
        sf.addRow(st.hint("One spring per direction at every base node.  "
                          "0 = free in that direction."))
        self.spring_box.setVisible(False)
        lay.addWidget(self.spring_box)

        # ---- custom supports (read-only) -----------------------------------
        self.custom_box = QtWidgets.QGroupBox("Custom supports (from the model)")
        cv = QtWidgets.QVBoxLayout(self.custom_box)
        self.custom_summary = QtWidgets.QLabel("")
        self.custom_summary.setWordWrap(True)
        cv.addWidget(self.custom_summary)
        self.custom_table = st.table(["Node", "x", "y", "z"])
        self.custom_table.setMaximumHeight(150)
        cv.addWidget(self.custom_table)
        self.custom_box.setVisible(False)
        lay.addWidget(self.custom_box)
        self.support.currentIndexChanged.connect(self._support_changed)
        self.support.currentIndexChanged.connect(self._emit_supports)
        for w in (self.spring_kx, self.spring_ky, self.spring_kz):
            w.valueChanged.connect(self._emit_supports)

        self.generate_btn = QtWidgets.QPushButton("Generate dome")
        self.generate_btn.setObjectName("primary")
        self.generate_btn.setToolTip("Build a geodesic dome from the inputs "
                                     "above (F5).")
        self.generate_btn.clicked.connect(
            lambda *_: self.generate_requested.emit())
        lay.addWidget(self.generate_btn)

        info = QtWidgets.QGroupBox("Generated geometry")
        iv = QtWidgets.QVBoxLayout(info)
        self.info = QtWidgets.QPlainTextEdit()
        self.info.setReadOnly(True)
        self.info.setMaximumHeight(188)
        iv.addWidget(self.info)
        self.schedule = st.table(["Fam", "Length [m]", "Chord factor", "Count"])
        # the table takes the free height of the box (no blank bands)
        self.schedule.setMinimumHeight(120)
        iv.addWidget(self.schedule, 1)
        lay.addWidget(info, 1)
        self._info_box = info

    # ------------------------------------------------------------ internal --
    def _support_changed(self, *_):
        key = self.support.currentData()
        self.spring_box.setVisible(key == geo.ELASTIC)
        self.custom_box.setVisible(key == geo.CUSTOM)

    def _emit_supports(self, *_):
        if not self._syncing:
            self.supports_changed.emit()

    def set_support_type(self, key):
        """Select a support type in the combo (no-op for an unknown key)."""
        idx = self.support.findData(key)
        if idx >= 0:
            self.support.setCurrentIndex(idx)
        self._support_changed()

    def _show_custom_supports(self, model):
        cs = model.custom_supports or {}
        n, nf, ns, nfree = custom_support_counts(cs)
        if n:
            self.custom_summary.setText(
                f"{n} supported node(s): {nf} fixed, {ns} spring and "
                f"{nfree} free direction(s).  Read only: edit the model file "
                f"or the CSV to change them.")
        else:
            self.custom_summary.setText(
                "This model has no per-node supports.  They come from a "
                "model file or a CSV import.  Choose another support type.")
        rows = [{"Node": int(nd),
                 "x": _support_text(c[0]), "y": _support_text(c[1]),
                 "z": _support_text(c[2])}
                for nd, c in sorted(cs.items())]
        st.fill_table(self.custom_table, rows)
        if not rows:
            self.custom_table.setColumnCount(4)
            self.custom_table.setHorizontalHeaderLabels(["Node", "x", "y", "z"])

    def _set_custom_mode(self, model):
        g = model.geom
        self._custom = bool(g.is_custom)
        for w in (self.frequency, self.radius, self.fraction,
                  self.fraction_slider):
            w.setEnabled(not self._custom)
        if self._custom:
            title = html.escape(g.title or "Custom truss")
            self.custom_note.setText(
                f"<b>{title}</b><br>"
                f"{g.n_nodes} nodes, {g.n_members} members, "
                f"{g.n_faces} panels.<br>"
                "Custom truss: geometry comes from the model file; use "
                "File &gt; New or Generate dome to return to a geodesic dome.")
            self._info_box.setTitle("Model geometry")
        else:
            self._info_box.setTitle("Generated geometry")
        self.custom_note.setVisible(self._custom)

    # ------------------------------------------------------------- binding --
    @property
    def is_custom(self):
        return self._custom

    def read_into(self, model):
        key = self.support.currentData()
        if key == geo.CUSTOM and not model.custom_supports:
            raise InputError(
                "Custom supports need per-node data from a model file or a "
                "CSV import.  This model has none.  Choose another support "
                "type in the Geometry tab.")
        model.support_type = key
        old = tuple(model.spring_stiffness)
        model.spring_stiffness = tuple(
            _spin_value(w, k, KN) for w, k in
            zip((self.spring_kx, self.spring_ky, self.spring_kz), old))

    def params(self):
        return (self.frequency.value(), self.radius.value(),
                self.fraction.value())

    def sync_from_model(self, model, set_inputs=True):
        """Show a model: geodesic inputs, custom notice, supports, summary.

        ``set_inputs=False`` keeps the frequency / radius / fraction the user
        typed (used right after Generate dome).
        """
        g = model.geom
        self._syncing = True
        try:
            if set_inputs and not g.is_custom:
                self.frequency.setValue(int(g.frequency))
                if math.isfinite(g.radius):
                    self.radius.setValue(float(g.radius))
                if math.isfinite(g.dome_fraction):
                    self.fraction.setValue(float(g.dome_fraction))
            self._set_custom_mode(model)
            for w, k in zip((self.spring_kx, self.spring_ky, self.spring_kz),
                            model.spring_stiffness):
                w.setValue(float(k) / KN)
            self.set_support_type(model.support_type)
        finally:
            self._syncing = False
        self.show_model(model)

    def show_model(self, model):
        try:
            lines = model.summary_lines()
        except Exception as exc:                     # never break the panel
            lines = [f"Summary not available: {exc}"]
        self.info.setPlainText("\n".join(lines))
        rows = [{"Fam": s["label"],
                 "Length [m]": round(s["length"], 4),
                 "Chord factor": float(s["chord_factor"]),
                 "Count": s["count"]}
                for s in model.geom.strut_schedule()]
        st.fill_table(self.schedule, rows, numeric_digits=6)
        self._show_custom_supports(model)


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
        if name in MATERIALS:
            mat = MATERIALS[name]
            text = mat.describe()
            if not getattr(mat, "aisc_checks_apply", True):
                text += ("\nThe member checks follow AISC 360, a steel code. "
                         "They do not apply to this material; use the "
                         "results for forces and displacements only.")
            self.material_info.setText(text)

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
        per_area = model.mass_per_plan_area()
        if math.isfinite(per_area):
            extra = f"({per_area:.1f} kg per m2 of plan area)"
        else:
            extra = "(no panels, so no plan area)"
        self.mass_label.setText(
            f"Structural mass  {model.steel_mass():,.0f} kg     {extra}")


# ================================================================== loads ====
def _cell_text(table, row, col):
    w = table.cellWidget(row, col)
    if isinstance(w, QtWidgets.QComboBox):
        return w.currentText()
    it = table.item(row, col)
    return it.text().strip() if it is not None else ""


def _fmt(value):
    return f"{float(value):.6g}"


def _num_item(value, integer=False, si=None):
    """A right-aligned number cell showing ``value`` (display units).

    ``si`` is the exact model value behind it (default ``value``).  It is
    kept in the item, so a value loaded from a model reads back unchanged
    unless the user edits the cell (the cell shows six significant digits).
    """
    text = str(int(value)) if integer else _fmt(value)
    it = QtWidgets.QTableWidgetItem(text)
    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    if not integer:
        it.setData(Qt.UserRole, [text, float(value if si is None else si)])
    return it


def _cell_exact(table, row, col):
    """The stored model value of a cell if its text was not edited."""
    it = table.item(row, col)
    if it is None:
        return None
    data = it.data(Qt.UserRole)
    try:
        text, si = data
        si = float(si)
    except (TypeError, ValueError):
        return None
    return si if it.text().strip() == str(text) else None


class LoadPanel(QtWidgets.QWidget):
    """Gravity, snow, wind, temperature, settlements, nodal loads and the
    load combination table."""

    changed = Signal()

    SETTLE_COLS = ["Node", "ux [mm]", "uy [mm]", "uz [mm]"]
    NODAL_COLS = ["Node", "Fx [kN]", "Fy [kN]", "Fz [kN]", "Case"]
    COMBO_COLS = ["Name", "D", "S", "L", "W", "T", "Service"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self._n_nodes = 0
        self._crown = 0
        self._support_nodes = []

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

        # ---- no-panel notice ------------------------------------------------
        self.no_faces_note = st.note(
            "This model has no panels, so area loads (dead, snow, live and "
            "wind) are zero.  Use self weight, the crown load and nodal "
            "loads.")
        self.no_faces_note.setVisible(False)
        lay.addWidget(self.no_faces_note)

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
        self.snow_unbalanced.setToolTip(
            "1.5 x snow on the half x >= 0 and 0.5 x snow on x < 0.  "
            "A simple test pattern, not the ASCE 7 curved-roof rule.")
        form.addRow("", self.snow_unbalanced)
        self.live = st.dspin(0.0, 20.0, 0.0, 0.05, 3, " kPa")
        self.live.setToolTip(
            "Applied on the horizontally projected area.  It enters only "
            "the combinations with an L factor (LC6 in the default set).")
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
        self.wind_az.setToolTip(
            "Direction the wind blows towards, in degrees from +x "
            "(counter-clockwise).  0 means wind from -x towards +x.")
        wf.addRow("Wind direction", self.wind_az)
        self.qz_label = QtWidgets.QLabel("")
        self.qz_label.setObjectName("metric")
        wf.addRow("", self.qz_label)
        lay.addWidget(self.wind_box)
        lay.addWidget(st.hint(
            "qz = 0.613 Kz Kzt Kd V^2.  The pressure coefficient varies along "
            "the meridian from +0.80 at the windward springing through -1.00 "
            "at the crown to -0.40 leeward, applied normal to the surface."))

        # ---- self-straining load T -------------------------------------------
        tbox = QtWidgets.QGroupBox("Temperature and settlement (load T)")
        tv = QtWidgets.QVBoxLayout(tbox)
        tf = QtWidgets.QFormLayout()
        self.delta_T = st.dspin(-200.0, 200.0, 0.0, 5.0, 2, " °C")
        self.delta_T.setToolTip(
            "Uniform temperature change of every strut.  + is heating.  "
            "It acts through the T factor of each combination.")
        tf.addRow("Temperature change ΔT", self.delta_T)
        tv.addLayout(tf)
        tv.addWidget(QtWidgets.QLabel("Support settlements"))
        self.settle_table = st.input_table(self.SETTLE_COLS)
        self.settle_table.setToolTip(
            "Prescribed movement of a supported node, in mm.  Only supported "
            "directions can move.  Negative uz is downward.")
        tv.addWidget(self.settle_table)
        tv.addLayout(self._add_remove_row(self._add_settlement,
                                          self.settle_table))
        lay.addWidget(tbox)
        lay.addWidget(st.hint(
            "Temperature and settlement form the self-straining load T.  "
            "Its factor is the T column of the combination table."))

        # ---- nodal loads -------------------------------------------------------
        nbox = QtWidgets.QGroupBox("Nodal loads")
        nv = QtWidgets.QVBoxLayout(nbox)
        self.nodal_table = st.input_table(self.NODAL_COLS)
        self.nodal_table.setToolTip(
            "Point loads in kN along the global axes (z is up).  Each load "
            "is added to the load case chosen in the Case column.")
        nv.addWidget(self.nodal_table)
        nv.addLayout(self._add_remove_row(self._add_nodal_load,
                                          self.nodal_table))
        lay.addWidget(nbox)
        self.node_hint = st.hint("")
        lay.addWidget(self.node_hint)

        # ---- combinations --------------------------------------------------
        box3 = QtWidgets.QGroupBox("Load combinations")
        v3 = QtWidgets.QVBoxLayout(box3)
        self.combos = QtWidgets.QTableWidget(0, len(self.COMBO_COLS))
        self.combos.setHorizontalHeaderLabels(self.COMBO_COLS)
        self.combos.horizontalHeaderItem(5).setToolTip(
            "Factor on the self-straining load T (temperature + settlement).")
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
            "Factors are editable.  T is the factor on temperature and "
            "settlement.  Tick 'Service' to use a row for the deflection "
            "check instead of the strength envelope."))
        lay.addStretch(1)

        for w in (self.wind_speed, self.wind_Kz, self.wind_Kzt, self.wind_Kd):
            w.valueChanged.connect(self._update_qz)
        self.wind_box.toggled.connect(self._update_qz)
        self._reset_rows()
        self._update_qz()

        for w in (self.dead, self.snow, self.live, self.crown_load,
                  self.wind_speed, self.wind_az, self.delta_T):
            # the lambdas drop the value: PySide6 refuses to pass it on to
            # a signal without arguments
            w.valueChanged.connect(lambda *_: self.changed.emit())
        for w in (self.self_weight, self.snow_unbalanced, self.wind_box):
            w.toggled.connect(lambda *_: self.changed.emit())
        self.settle_table.itemChanged.connect(lambda *_: self.changed.emit())
        self.nodal_table.itemChanged.connect(lambda *_: self.changed.emit())

    # ------------------------------------------------------------ internal --
    def _add_remove_row(self, add_slot, table):
        row = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("Add")
        add.clicked.connect(add_slot)
        rem = QtWidgets.QPushButton("Remove")
        rem.clicked.connect(lambda _=False, t=table: self._remove_from(t))
        row.addWidget(add)
        row.addWidget(rem)
        row.addStretch(1)
        return row

    @staticmethod
    def _remove_from(table):
        i = table.currentRow()
        if i < 0:
            i = table.rowCount() - 1
        if i >= 0:
            table.removeRow(i)

    def _update_qz(self, *_):
        q = (0.613 * self.wind_Kz.value() * self.wind_Kzt.value()
             * self.wind_Kd.value() * self.wind_speed.value() ** 2)
        self.qz_label.setText(f"qz = {q / 1e3:.3f} kPa    "
                              f"qz*G = {q * self.wind_G.value() / 1e3:.3f} kPa")

    def _set_row(self, i, combo):
        name = QtWidgets.QTableWidgetItem(combo.name)
        name.setToolTip(combo.name)
        self.combos.setItem(i, 0, name)
        for j, key in enumerate(("D", "S", "L", "W", "T"), start=1):
            self.combos.setItem(
                i, j, _num_item(float(getattr(combo, key, 0.0) or 0.0)))
        chk = QtWidgets.QTableWidgetItem()
        chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
        chk.setCheckState(Qt.Checked if combo.serviceability else Qt.Unchecked)
        self.combos.setItem(i, 6, chk)

    def _add_row(self):
        i = self.combos.rowCount()
        self.combos.insertRow(i)
        self._set_row(i, LoadCombination(f"LC{i + 1}", D=1.0))

    def _remove_row(self):
        i = self.combos.currentRow()
        if i >= 0 and self.combos.rowCount() > 1:
            self.combos.removeRow(i)

    def _fit_columns(self):
        """Name takes the spare width; factors and Service stay readable."""
        hh = self.combos.horizontalHeader()
        hh.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        for j in range(1, 6):
            hh.setSectionResizeMode(j, QtWidgets.QHeaderView.Fixed)
            self.combos.setColumnWidth(j, 33)
        hh.setSectionResizeMode(6, QtWidgets.QHeaderView.ResizeToContents)

    def _reset_rows(self):
        rows = default_combinations()
        self.combos.setRowCount(len(rows))
        for i, c in enumerate(rows):
            self._set_row(i, c)
        self._fit_columns()

    # ---- settlement / nodal-load rows ---------------------------------------
    def add_settlement(self, node, ux=0.0, uy=0.0, uz=0.0, si=None):
        """Append a settlement row (node number, movements in mm).

        ``si`` = (ux, uy, uz) in m: the exact model values behind the row.
        """
        t = self.settle_table
        i = t.rowCount()
        t.insertRow(i)
        t.setItem(i, 0, _num_item(node, integer=True))
        exact = si or [float(v) * MM for v in (ux, uy, uz)]
        for j, (v, e) in enumerate(zip((ux, uy, uz), exact), start=1):
            t.setItem(i, j, _num_item(float(v), si=e))
        return i

    def add_nodal_load(self, node, Fx=0.0, Fy=0.0, Fz=0.0, case="L",
                       si=None):
        """Append a nodal-load row (node number, forces in kN, case).

        ``si`` = (Fx, Fy, Fz) in N: the exact model values behind the row.
        """
        t = self.nodal_table
        i = t.rowCount()
        t.insertRow(i)
        t.setItem(i, 0, _num_item(node, integer=True))
        exact = si or [float(v) * KN for v in (Fx, Fy, Fz)]
        for j, (v, e) in enumerate(zip((Fx, Fy, Fz), exact), start=1):
            t.setItem(i, j, _num_item(float(v), si=e))
        cb = QtWidgets.QComboBox()
        cb.setObjectName("cellcombo")
        cb.addItems(list(LOAD_CASES))
        case = str(case or "L").strip().upper()
        cb.setCurrentText(case if case in LOAD_CASES else "L")
        cb.setToolTip("D dead, S snow, L live, W wind")
        cb.currentTextChanged.connect(lambda *_: self.changed.emit())
        t.setCellWidget(i, 4, cb)
        return i

    def _add_settlement(self):
        used = set()
        for i in range(self.settle_table.rowCount()):
            try:
                used.add(int(float(_cell_text(self.settle_table, i, 0))))
            except ValueError:
                pass
        free = [n for n in self._support_nodes if n not in used]
        node = free[0] if free else (self._support_nodes[0]
                                     if self._support_nodes else 0)
        self.add_settlement(node)
        self.settle_table.setCurrentCell(self.settle_table.rowCount() - 1, 3)

    def _add_nodal_load(self):
        self.add_nodal_load(self._crown)
        self.nodal_table.setCurrentCell(self.nodal_table.rowCount() - 1, 3)

    def clear_node_tables(self):
        """Empty the settlement and nodal-load tables (node numbers changed)."""
        self.settle_table.setRowCount(0)
        self.nodal_table.setRowCount(0)

    def show_model_info(self, model):
        """Update the notices that depend on the geometry (not the inputs)."""
        g = model.geom
        self._n_nodes = int(g.n_nodes)
        self._crown = int(g.crown_node)
        try:
            self._support_nodes = [int(n) for n in model.support_nodes()]
        except Exception:
            self._support_nodes = [int(n) for n in g.base_nodes]
        self.no_faces_note.setVisible(g.n_faces == 0)
        self.node_hint.setText(
            f"Nodes are numbered 0 to {g.n_nodes - 1}; the crown is node "
            f"{g.crown_node}.  Tick 'Node numbers' in the Model view to see "
            f"them.")

    # ---- parsing ------------------------------------------------------------
    @staticmethod
    def _float(table, i, j, what, errors, scale=1.0):
        """Cell value in model units (display value x ``scale``)."""
        exact = _cell_exact(table, i, j)
        if exact is not None:
            return exact
        text = _cell_text(table, i, j)
        if text == "":
            return 0.0
        try:
            v = float(text)
        except ValueError:
            errors.append(f"{what}, row {i + 1}: '{text}' is not a number.")
            return 0.0
        if not math.isfinite(v):
            errors.append(f"{what}, row {i + 1}: '{text}' is not a finite "
                          f"number.")
            return 0.0
        return v * scale

    def _node(self, table, i, what, n, errors):
        text = _cell_text(table, i, 0)
        if text == "":
            errors.append(f"{what}, row {i + 1}: the node number is missing.")
            return None
        try:
            v = float(text)
            if v != int(v):
                raise ValueError
            node = int(v)
        except ValueError:
            errors.append(f"{what}, row {i + 1}: node '{text}' is not a "
                          f"whole number.")
            return None
        if node < 0 or node >= n:
            errors.append(f"{what}, row {i + 1}: node {node} does not exist. "
                          f"Nodes are 0 to {n - 1}.")
            return None
        return node

    def _read_settlements(self, model, errors):
        n = model.geom.n_nodes
        out = []
        rows = []
        for i in range(self.settle_table.rowCount()):
            node = self._node(self.settle_table, i, "Settlements", n, errors)
            vals = [self._float(self.settle_table, i, j, "Settlements", errors,
                                scale=MM) for j in (1, 2, 3)]
            if node is None:
                continue
            rows.append((i, node, vals))
            out.append({"node": node, "ux": vals[0], "uy": vals[1],
                        "uz": vals[2]})
        if not rows:
            return out
        try:
            fixed, springs = model.support_conditions()
        except ValueError as exc:
            errors.append(f"Settlements cannot be checked: {exc}")
            return out
        held = set(int(d) for d in fixed) | set(int(d) for d in springs)
        supported = sorted({d // 3 for d in held})
        for i, node, vals in rows:
            if node not in supported:
                shown = ", ".join(str(s) for s in supported[:12])
                more = " ..." if len(supported) > 12 else ""
                errors.append(f"Settlements, row {i + 1}: node {node} is not "
                              f"a supported node.  Supported nodes: "
                              f"{shown}{more}.")
                continue
            for a, v in enumerate(vals):
                if v != 0.0 and 3 * node + a not in held:
                    errors.append(f"Settlements, row {i + 1}: node {node} is "
                                  f"free in {'xyz'[a]}, so u{'xyz'[a]} must "
                                  f"be 0.")
        return out

    def _read_nodal_loads(self, model, errors):
        n = model.geom.n_nodes
        out = []
        for i in range(self.nodal_table.rowCount()):
            node = self._node(self.nodal_table, i, "Nodal loads", n, errors)
            vals = [self._float(self.nodal_table, i, j, "Nodal loads", errors,
                                scale=KN) for j in (1, 2, 3)]
            case = (_cell_text(self.nodal_table, i, 4) or "L").upper()
            if case not in LOAD_CASES:
                errors.append(f"Nodal loads, row {i + 1}: case '{case}' must "
                              f"be D, S, L or W.")
                continue
            if node is None:
                continue
            out.append({"node": node, "Fx": vals[0], "Fy": vals[1],
                        "Fz": vals[2], "case": case})
        return out

    def _read_combinations(self, errors):
        combos = []
        names = {}
        for i in range(self.combos.rowCount()):
            name_item = self.combos.item(i, 0)
            name = (name_item.text().strip() if name_item else "") \
                or f"LC{i + 1}"
            f = [self._float(self.combos, i, j, "Load combinations", errors)
                 for j in range(1, 6)]
            if name in names:
                errors.append(f"Load combinations, rows {names[name] + 1} and "
                              f"{i + 1}: two combinations are named "
                              f"'{name}'.  Give each a different name.")
            names.setdefault(name, i)
            chk = self.combos.item(i, 6)
            combos.append(LoadCombination(
                name, D=f[0], S=f[1], L=f[2], W=f[3],
                serviceability=bool(chk and chk.checkState() == Qt.Checked),
                T=f[4]))
        if not combos:
            errors.append("Load combinations: add at least one row.")
        elif all(c.serviceability for c in combos):
            errors.append("Load combinations: untick Service on at least one "
                          "row.  The strength check needs one.")
        return combos

    # ------------------------------------------------------------- binding --
    def read_into(self, model):
        """Copy the inputs into ``model``; raise InputError if any is bad.

        Nothing is written to the model when an input is bad.
        """
        errors = []
        settlements = self._read_settlements(model, errors)
        nodal = self._read_nodal_loads(model, errors)
        combos = self._read_combinations(errors)
        if errors:
            raise InputError("Please check the Loads tab:\n\n- "
                             + "\n- ".join(errors))

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
        ls.temperature_change = _spin_value(self.delta_T,
                                            ls.temperature_change)
        ls.settlements = settlements
        ls.nodal_loads = nodal
        model.combinations = combos
        self.show_model_info(model)

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
        self.delta_T.setValue(float(getattr(ls, "temperature_change", 0.0)
                                    or 0.0))
        self.clear_node_tables()
        for s in getattr(ls, "settlements", None) or []:
            si = [float(s.get(k, 0.0) or 0.0) for k in ("ux", "uy", "uz")]
            self.add_settlement(int(s["node"]), *[v / MM for v in si], si=si)
        for p in getattr(ls, "nodal_loads", None) or []:
            si = [float(p.get(k, 0.0) or 0.0) for k in ("Fx", "Fy", "Fz")]
            self.add_nodal_load(int(p["node"]), *[v / KN for v in si],
                                case=p.get("case", "L"), si=si)
        self.combos.setRowCount(len(model.combinations))
        for i, c in enumerate(model.combinations):
            self._set_row(i, c)
        self._update_qz()
        self.show_model_info(model)


# ================================================================= design ====
class DesignPanel(QtWidgets.QWidget):
    """Code parameters, analysis options and the action buttons."""

    analyse_requested = Signal()
    size_requested = Signal()
    verify_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        # in a scroll area, so the window fits a 1366 x 768 screen
        inner = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(inner)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(9)
        _scroll_wrap(self, inner)

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
        self.thickness_factor = st.dspin(0.50, 1.00, 1.00, 0.01, 2)
        self.thickness_factor.setToolTip(
            "Design wall thickness = factor x nominal wall thickness.  It is "
            "used in the strength checks only.  AISC 360-16 B4.2 uses 0.93 "
            "for ERW HSS such as ASTM A500.  1.00 uses the nominal wall.")
        form.addRow("Design wall thickness factor", self.thickness_factor)
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
        self.analyse_btn.clicked.connect(
            lambda *_: self.analyse_requested.emit())
        lay.addWidget(self.analyse_btn)

        self.size_btn = QtWidgets.QPushButton("Auto-size the strut families")
        self.size_btn.clicked.connect(lambda *_: self.size_requested.emit())
        lay.addWidget(self.size_btn)

        self.verify_btn = QtWidgets.QPushButton("Run verification suite")
        self.verify_btn.clicked.connect(
            lambda *_: self.verify_requested.emit())
        lay.addWidget(self.verify_btn)
        lay.addStretch(1)

    def read_into(self, model):
        model.K_effective = self.K.value()
        model.deflection_limit = float(self.defl_limit.value())
        model.phi_c = self.phi_c.value()
        model.phi_t = self.phi_t.value()
        model.design_thickness_factor = _spin_value(
            self.thickness_factor,
            getattr(model, "design_thickness_factor", 1.0))

    def sync_from_model(self, model):
        self.K.setValue(model.K_effective)
        self.defl_limit.setValue(int(round(model.deflection_limit)))
        self.phi_c.setValue(model.phi_c)
        self.phi_t.setValue(model.phi_t)
        self.thickness_factor.setValue(
            float(getattr(model, "design_thickness_factor", 1.0) or 1.0))
