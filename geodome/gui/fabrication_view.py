"""Fabrication schedules in the GUI, and the CAD (DXF) export.

``install(win)`` adds

* the result tab "Fabrication" with the sub-tabs Summary, Struts (cut list),
  Hubs, Panels and Dihedrals, from :func:`geodome.fabrication.
  fabrication_tables`;
* a hub offset input [mm]: the distance from a hub centre to the end of a
  tube, so that ``cut length = L - 2 x hub offset`` (L = hub-to-hub length);
* "Export cut list (CSV)...", which writes every table with
  :func:`geodome.fabrication.export_fabrication_csv`;
* "File > Export geometry to CAD (DXF)...", which writes the geometry with
  :func:`geodome.cad.export_dxf` (struts on STRUT_<family> layers, panels as
  3DFACE, hubs as POINT; when analysis results exist the user can choose to
  colour the struts by utilisation).

The tables are rebuilt when the model changes, after each analysis (the
auto-sizer may have changed the sections, and so the tube masses), when a
section or the hub offset changes, and on demand (Refresh).  The extension
keeps the last tables as ``.fabrication_tables`` (None when the hub offset is
too large for the shortest strut).  NaN values, for example the chord factor
of a custom truss that is not spherical, show as "-".
"""

from __future__ import annotations

import os
import re

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import QTimer

from .. import cad
from .. import fabrication as fab
from . import ext_util as eu
from . import style as st

SUBTABS = ("summary", "struts", "hubs", "panels", "dihedrals")
SUBTAB_TITLES = {"summary": "Summary", "struts": "Struts", "hubs": "Hubs",
                 "panels": "Panels", "dihedrals": "Dihedrals"}

# (key in the fabrication table, header, decimals or None for text)
COLUMNS = {
    "summary": [("quantity", "Quantity", None), ("value", "Value", 3),
                ("unit", "Unit", None), ("note", "Note", None)],
    "struts": [("family", "Family", None), ("colour", "Colour", None),
               ("count", "Count", 0), ("section", "Section", None),
               ("length_mm", "Length [mm]", 1),
               ("cut_length_mm", "Cut length [mm]", 1),
               ("chord_factor", "Chord factor", 5),
               ("axial_angle_deg", "Axial angle [deg]", 2),
               ("tangent_angle_deg", "Tangent angle [deg]", 2),
               ("mass_each_kg", "Mass each [kg]", 2),
               ("mass_total_kg", "Mass total [kg]", 1)],
    "hubs": [("hub_type", "Hub", None), ("location", "Location", None),
             ("valence", "Struts at hub", 0), ("struts", "Strut mix", None),
             ("arrangement", "Order around hub", None),
             ("count", "Count", 0), ("gaps_deg", "Gaps [deg]", None),
             ("example_node", "Example node", 0)],
    "panels": [("panel_type", "Panel", None), ("count", "Count", 0),
               ("side_1_mm", "Side 1 [mm]", 1),
               ("side_2_mm", "Side 2 [mm]", 1),
               ("side_3_mm", "Side 3 [mm]", 1),
               ("angle_1_deg", "Angle 1 [deg]", 2),
               ("angle_2_deg", "Angle 2 [deg]", 2),
               ("angle_3_deg", "Angle 3 [deg]", 2),
               ("area_m2", "Area [m2]", 3),
               ("total_area_m2", "Total area [m2]", 2),
               ("handedness", "Handedness", None)],
    "dihedrals": [("strut_family", "Strut family", None),
                  ("panels", "Panels", None), ("count", "Count", 0),
                  ("dihedral_deg", "Dihedral [deg]", 2),
                  ("bend_deg", "Bend [deg]", 2)],
}

EMPTY_TEXT = {
    "panels": "No panels: this truss has no faces.",
    "dihedrals": "No dihedral angles: no strut is shared by two panels.",
}

_UNIT_DECIMALS = {"mm": 1, "m": 3, "m2": 3, "kg": 1, "deg": 4, "-": 0}

DXF_FILTER_FAMILY = "DXF, struts coloured by family (*.dxf)"
DXF_FILTER_UTIL = "DXF, struts coloured by utilisation (*.dxf)"


def summary_rows(rows):
    """Summary rows with the value formatted by its unit (as text)."""
    out = []
    for r in rows:
        v = r.get("value", "")
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            dec = _UNIT_DECIMALS.get(r.get("unit", ""), 3)
            if isinstance(v, int):
                text = eu.fmt(v)
            else:
                # the closed-form checks are ~1e-14: show them, not 0.0000
                text = eu.fmt(float(v), dec, small_sci=True)
        else:
            text = str(v)
        out.append(dict(r, value=text))
    return out


class FabricationExtension(QtCore.QObject):
    """The Fabrication result tab and the DXF export."""

    def __init__(self, win):
        super().__init__(win)
        self.win = win
        self.fabrication_tables = None
        self.error = ""
        self._build()
        self._debounce = QTimer(self.widget)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(250)
        self._debounce.timeout.connect(self.refresh)

    # ------------------------------------------------------------- build --
    def _build(self):
        w = self.widget = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Hub offset:"))
        self.hub_offset = st.dspin(0.0, 1000.0, 0.0, 5.0, 1, " mm")
        self.hub_offset.setToolTip(
            "Distance from the hub centre to the end of the tube. "
            "Cut length = strut length - 2 x hub offset.")
        self.hub_offset.valueChanged.connect(
            lambda _v: self._debounce.start())
        row.addWidget(self.hub_offset)
        self.refresh_btn = QtWidgets.QPushButton("Refresh")
        self.refresh_btn.setToolTip("Rebuild the schedules from the current "
                                    "model.")
        self.refresh_btn.clicked.connect(lambda: self.refresh())
        row.addWidget(self.refresh_btn)
        row.addStretch(1)
        self.export_btn = QtWidgets.QPushButton("Export cut list (CSV)...")
        self.export_btn.setToolTip("Save every fabrication table as a CSV "
                                   "file in a folder.")
        self.export_btn.clicked.connect(lambda: self.export_csv())
        row.addWidget(self.export_btn)
        lay.addLayout(row)

        self.info = st.hint("")
        lay.addWidget(self.info)

        self.subtabs = QtWidgets.QTabWidget()
        self.tables = {}
        for key in SUBTABS:
            t = st.table([c[1] for c in COLUMNS[key]])
            self.tables[key] = t
            self.subtabs.addTab(t, SUBTAB_TITLES[key])
        self.subtabs.setCurrentWidget(self.tables["struts"])
        lay.addWidget(self.subtabs, 1)

    def install(self):
        win = self.win
        eu.add_result_tab(win, self.widget, "Fabrication")
        self.act_dxf = win._act(
            "Export geometry to CAD (&DXF)...",
            lambda _=False: self.export_dxf(),
            tip="Struts, panels and hubs as a DXF file for CAD (units: m)")
        menu = getattr(win, "menu_file", None) or eu.ensure_menu(win, "File")
        eu.insert_before_exit(menu, self.act_dxf)
        sp = getattr(win, "section_panel", None)
        if sp is not None and hasattr(sp, "changed"):
            sp.changed.connect(lambda *a: self._debounce.start())
        eu.register(win, self)
        self.refresh()
        return self

    # -------------------------------------------------------- extension API --
    def on_model_changed(self, model):
        self.refresh(model)

    def on_results(self, results):
        # Auto-sizing changes the sections and then re-analyses: refresh the
        # tube masses.
        self.refresh()

    # ----------------------------------------------------------- refresh --
    def hub_offset_m(self):
        return float(self.hub_offset.value()) / 1000.0

    def refresh(self, model=None):
        """Rebuild the tables for ``model`` (default: the window's model).

        Returns the tables, or None when they cannot be built (the message is
        shown above the tables).
        """
        self._debounce.stop()
        model = self.win.model if model is None else model
        try:
            tables = fab.fabrication_tables(model,
                                            hub_offset=self.hub_offset_m())
        except Exception as exc:
            self.fabrication_tables = None
            self.error = str(exc)
            for key in SUBTABS:
                eu.fill(self.tables[key], COLUMNS[key], [])
            self.info.setText("Cannot build the schedules: " + self.error)
            self.info.setStyleSheet(f"color: {st.BAD_RED};")
            return None
        self.fabrication_tables = tables
        self.error = ""
        for key in SUBTABS:
            rows = tables.get(key) or []
            if key == "summary":
                rows = summary_rows(rows)
            eu.fill(self.tables[key], COLUMNS[key], rows,
                    empty_text=EMPTY_TEXT.get(key, "No rows."))
        self.info.setStyleSheet("")
        self.info.setText(self._info_text(model, tables))
        return tables

    def _info_text(self, model, tables):
        g = model.geom
        name = (g.title or "Custom truss") if g.is_custom else \
            f"{g.frequency}V dome"
        n_cut = len(tables.get("struts") or [])
        n_hub = len(tables.get("hubs") or [])
        n_pan = len(tables.get("panels") or [])
        mass = sum(r.get("mass_total_kg", 0.0) or 0.0
                   for r in tables.get("struts") or []
                   if not eu.is_missing(r.get("mass_total_kg")))
        text = (f"{name}: {g.n_members} struts in {n_cut} cut lengths, "
                f"{g.n_nodes} hubs of {n_hub} types, "
                + (f"{g.n_faces} panels of {n_pan} types. " if n_pan else
                   "no panels. ")
                + f"Tube mass {mass:,.0f} kg (hubs not included).")
        notes = [r["note"] for r in tables.get("notes") or []]
        if notes:
            text += "  Note: " + "  ".join(notes[:2])
        return text

    # ------------------------------------------------------------ export --
    def _prefix(self):
        g = self.win.model.geom
        if g.is_custom:
            p = re.sub(r"[^A-Za-z0-9]+", "_", g.title or "truss").strip("_")
            return (p or "truss").lower()[:40]
        return f"dome_{g.frequency}V"

    def export_csv(self, folder=None):
        """Write the cut list and every other table as CSV files.

        Asks for a folder when ``folder`` is None.  Returns the written paths.
        """
        win = self.win
        tables = self.fabrication_tables or self.refresh()
        if tables is None:
            QtWidgets.QMessageBox.warning(
                win, "Export cut list",
                "The schedules could not be built: " + self.error)
            return []
        if folder is None:
            folder = QtWidgets.QFileDialog.getExistingDirectory(
                win, "Choose a folder for the cut list", os.getcwd())
            if not folder:
                return []
        try:
            files = fab.export_fabrication_csv(tables, folder,
                                               prefix=self._prefix())
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Export cut list", str(exc))
            return []
        win.log("Fabrication tables written:\n  " + "\n  ".join(files))
        win.status(f"{len(files)} fabrication CSV files written.")
        return files

    def _usable_results(self):
        """The window's results if they belong to the current model."""
        res = getattr(self.win, "results", None)
        if res is None or getattr(res, "util_env", None) is None:
            return None
        if getattr(res, "model", None) is not self.win.model:
            return None
        if len(res.util_env) != self.win.model.geom.n_members:
            return None
        return res

    def export_dxf(self, path=None, by_utilisation=False):
        """Write the geometry to a DXF file.

        Asks for the file (and the colouring) when ``path`` is None.  Returns
        the path or None.
        """
        win = self.win
        results = self._usable_results()
        if path is None:
            filters = [DXF_FILTER_FAMILY]
            if results is not None:
                filters.append(DXF_FILTER_UTIL)
            name = self._prefix() + ".dxf"
            path, chosen = QtWidgets.QFileDialog.getSaveFileName(
                win, "Export geometry to CAD", os.path.join(os.getcwd(), name),
                ";;".join(filters))
            if not path:
                return None
            by_utilisation = chosen == DXF_FILTER_UTIL
        if not path.lower().endswith(".dxf"):
            path += ".dxf"
        use = results if (by_utilisation and results is not None) else None
        if by_utilisation and results is None:
            win.log("DXF: no analysis results for this model, so the struts "
                    "are coloured by family.")
        try:
            cad.export_dxf(win.model, path, results=use)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(win, "Export to CAD", str(exc))
            return None
        g = win.model.geom
        how = "utilisation" if use is not None else "strut family"
        win.log(f"DXF written to {path}: {g.n_members} struts, {g.n_faces} "
                f"panels, {g.n_nodes} hubs; units m; struts coloured by "
                f"{how}.")
        win.status("Geometry exported to DXF.")
        return path


def install(win):
    """Add the Fabrication tab and the DXF export to ``win``.

    The extension registers itself in ``win.extensions``.  Calling install
    again on the same window returns the extension already there.
    """
    for ext in getattr(win, "extensions", None) or []:
        if isinstance(ext, FabricationExtension):
            return ext
    return FabricationExtension(win).install()
