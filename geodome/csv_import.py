"""Read a general pin-jointed space truss from plain CSV files.

This is the reader behind the GUI's *File > Import truss from CSV...*.  It
turns a few small tables, as a spreadsheet would save them, into a
:class:`~geodome.model.DomeModel` with a custom geometry
(:func:`geodome.geometry.custom_geometry`).  No Qt and no matplotlib.

File format
-----------
One table per file.  Comma, semicolon or tab separated.  An optional header
row is skipped, as are blank lines and lines that start with ``#``.  Units
are those of the GUI: metres, kN and kN/m.

``nodes.csv`` (required)
    ``id, x, y, z`` -- node id (any integer) and coordinates [m], z up.
``members.csv`` (required)
    ``i, j[, group]`` -- node ids of each bar and, optionally, its strut
    family (give a group on every row or on none; without groups the
    families are formed by bar length).
``supports.csv`` (optional)
    ``node, cx, cy, cz`` -- per direction: ``1`` = fixed, ``0`` = free, any
    other positive number = spring stiffness [kN/m]; the words ``fixed`` and
    ``free`` work too.  Without this file the lowest nodes are pinned.
``loads.csv`` (optional)
    ``node, Fx, Fy, Fz[, case]`` -- nodal force [kN] in global axes and its
    load case ``D``, ``S``, ``L`` or ``W`` (default ``L``).
``faces.csv`` (optional)
    ``i, j, k`` -- triangular cladding panels.  Only needed for area loads
    (dead, snow, live, wind); without it they are zero.

The files sit in one folder.  A common prefix is allowed: picking
``tower_nodes.csv`` reads ``tower_members.csv``, ``tower_supports.csv``, ...

Node ids are renumbered 0, 1, ..., n-1 in ascending id order, because the
program numbers nodes from 0 (a note says so when the ids change).

The imported model starts with self weight on, area loads at zero, wind off,
the default section on every strut family and the ASCE 7-16 Sec. 2.3.1
strength combinations with live load (:func:`import_combinations`).
"""

from __future__ import annotations

import math
import os

from . import geometry as geo
from .model import DomeModel, LoadCombination, LoadSpec

KN = 1.0e3                        # kN -> N  (and kN/m -> N/m)
CASES = ("D", "S", "L", "W")
TABLES = ("nodes", "members", "supports", "loads", "faces")


class CSVImportError(ValueError):
    """Raised with a plain message that names the file and line."""


# ------------------------------------------------------------- low level --
def _split(line, delim):
    if delim is None:
        return line.split()
    return [c.strip() for c in line.split(delim)]


def _detect_delimiter(line):
    for d in (";", "\t", ","):
        if d in line:
            return d
    return None                   # whitespace


def read_rows(path):
    """Return ``[(line_number, [cells])]`` of a CSV file, header skipped.

    The delimiter is detected from the first data line.  With a semicolon
    delimiter a decimal comma is accepted (``1,5`` reads 1.5).
    """
    name = os.path.basename(path)
    try:
        with open(path, "r", encoding="utf-8-sig") as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        raise CSVImportError(f"Cannot read {name}: {exc.strerror or exc}.") \
            from None
    rows = []
    delim = None
    first = True
    for k, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if first:
            delim = _detect_delimiter(line)
        cells = _split(line, delim)
        if delim == ";":
            cells = [c.replace(",", ".") for c in cells]
        while cells and cells[-1] == "":
            cells.pop()
        if first:
            first = False
            if cells and not _is_number(cells[0]):
                continue          # header row
        rows.append((k, cells))
    return rows


def _is_number(text):
    try:
        float(text)
    except (TypeError, ValueError):
        return False
    return True


def _num(text, name, line, what):
    try:
        v = float(text)
    except (TypeError, ValueError):
        raise CSVImportError(f"{name}, line {line}: {what} must be a number, "
                             f"not {text!r}.") from None
    if not math.isfinite(v):
        raise CSVImportError(f"{name}, line {line}: {what} is not finite.")
    return v


def _int(text, name, line, what):
    v = _num(text, name, line, what)
    if v != round(v):
        raise CSVImportError(f"{name}, line {line}: {what} must be a whole "
                             f"number, not {text!r}.")
    return int(round(v))


def _need(cells, n, name, line, columns):
    if len(cells) < n:
        raise CSVImportError(f"{name}, line {line}: expected {columns}, "
                             f"found {len(cells)} value(s).")


def _support_entry(text, name, line, axis):
    """1 / fixed -> True, 0 / free / blank -> False, k > 0 -> spring [N/m]."""
    t = (text or "").strip().lower()
    if t in ("", "0", "free", "false", "no"):
        return False
    if t in ("fixed", "true", "yes", "pinned", "rigid"):
        return True
    k = _num(text, name, line, f"c{axis}")
    if k < 0.0:
        raise CSVImportError(f"{name}, line {line}: c{axis} = {text} is "
                             f"negative. Use 1 (fixed), 0 (free) or a "
                             f"positive spring stiffness in kN/m.")
    if k == 0.0:
        return False
    if k == 1.0:
        return True
    return k * KN


# ---------------------------------------------------------------- files --
def find_files(path):
    """Locate the table files from ``nodes.csv`` (or its folder).

    Returns ``{"nodes": path, "members": path or None, ...}``.  ``path`` may
    be the folder, ``.../nodes.csv`` or ``.../<prefix>nodes.csv``.
    """
    if os.path.isdir(path):
        folder, prefix = path, ""
    else:
        folder, base = os.path.split(path)
        stem = base[:-4] if base.lower().endswith(".csv") else base
        if stem.lower().endswith("nodes"):
            prefix = stem[:-5]
        else:
            raise CSVImportError(f"Pick the node file (nodes.csv), not "
                                 f"{base}.")
    listing = {f.lower(): f for f in os.listdir(folder or ".")}
    out = {}
    for t in TABLES:
        want = f"{prefix}{t}.csv".lower()
        out[t] = (os.path.join(folder, listing[want]) if want in listing
                  else None)
    if out["nodes"] is None:
        raise CSVImportError(f"No {prefix}nodes.csv in the chosen folder.")
    if out["members"] is None:
        raise CSVImportError(f"No {prefix}members.csv next to "
                             f"{prefix}nodes.csv.")
    return out


def import_combinations():
    """ASCE 7-16 Sec. 2.3.1 strength combinations with live load, + service.

    Basic combinations 1 to 5 with roof live/rain taken as snow S.  The
    "L or 0.5W" term of combination 3 gives two rows: LC3 with L and LC6
    with 0.5W.  ``SLS`` is unfactored and used for the deflection check.
    """
    return [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6L + 0.5S", D=1.2, L=1.6, S=0.5),
        LoadCombination("LC3  1.2D + 1.6S + 1.0L", D=1.2, S=1.6, L=1.0),
        LoadCombination("LC4  1.2D + 1.0W + 1.0L + 0.5S",
                        D=1.2, W=1.0, L=1.0, S=0.5),
        LoadCombination("LC5  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("LC6  1.2D + 1.6S + 0.5W", D=1.2, S=1.6, W=0.5),
        LoadCombination("SLS  1.0D + 1.0L + 1.0S", D=1.0, L=1.0, S=1.0,
                        serviceability=True),
    ]


# ---------------------------------------------------------------- reader --
def read_truss_csv(nodes, members, supports=None, loads=None, faces=None,
                   title="Imported truss"):
    """Build a :class:`DomeModel` from the CSV files (paths; None = absent).

    Returns ``(model, notes)``; ``notes`` are short lines for the log.
    Raises :class:`CSVImportError` with the file and line of the problem.
    """
    notes = []
    # ---- nodes ------------------------------------------------------------
    name = os.path.basename(nodes)
    ids, coords = [], []
    seen = {}
    for line, c in read_rows(nodes):
        _need(c, 4, name, line, "id, x, y, z")
        nid = _int(c[0], name, line, "id")
        if nid in seen:
            raise CSVImportError(f"{name}, line {line}: node id {nid} is used "
                                 f"twice (first on line {seen[nid]}).")
        seen[nid] = line
        ids.append(nid)
        coords.append([_num(c[1], name, line, "x"),
                       _num(c[2], name, line, "y"),
                       _num(c[3], name, line, "z")])
    if len(ids) < 2:
        raise CSVImportError(f"{name}: at least two nodes are needed.")
    order = sorted(range(len(ids)), key=lambda k: ids[k])
    index = {ids[k]: i for i, k in enumerate(order)}
    xyz = [coords[k] for k in order]
    if sorted(ids) != list(range(len(ids))):
        notes.append(f"Node ids {min(ids)}..{max(ids)} were renumbered "
                     f"0..{len(ids) - 1} in id order (the program numbers "
                     f"nodes from 0).")

    def node_of(text, fname, line, what="node"):
        nid = _int(text, fname, line, what)
        if nid not in index:
            raise CSVImportError(f"{fname}, line {line}: {what} {nid} is not "
                                 f"in {name}.")
        return index[nid]

    # ---- members ----------------------------------------------------------
    mname = os.path.basename(members)
    bars, groups = [], []
    for line, c in read_rows(members):
        _need(c, 2, mname, line, "i, j[, group]")
        i = node_of(c[0], mname, line, "node i")
        j = node_of(c[1], mname, line, "node j")
        bars.append([i, j])
        groups.append(_int(c[2], mname, line, "group")
                      if len(c) > 2 and c[2] != "" else None)
    if not bars:
        raise CSVImportError(f"{mname}: no members found.")
    have = [g is not None for g in groups]
    if all(have):
        grp = groups
    elif not any(have):
        grp = None
    else:
        raise CSVImportError(f"{mname}: give a group on every member or on "
                             f"none.")

    # ---- faces ------------------------------------------------------------
    tri = None
    if faces:
        fname = os.path.basename(faces)
        tri = []
        for line, c in read_rows(faces):
            _need(c, 3, fname, line, "i, j, k")
            tri.append([node_of(c[a], fname, line, "node") for a in range(3)])

    # ---- supports ---------------------------------------------------------
    custom = {}
    if supports:
        sname = os.path.basename(supports)
        for line, c in read_rows(supports):
            _need(c, 4, sname, line, "node, cx, cy, cz")
            nd = node_of(c[0], sname, line)
            if nd in custom:
                raise CSVImportError(f"{sname}, line {line}: node {c[0]} is "
                                     f"listed twice.")
            custom[nd] = [_support_entry(c[1 + a], sname, line, "xyz"[a])
                          for a in range(3)]
        if not custom:
            raise CSVImportError(f"{sname}: no supports found.")

    # ---- loads ------------------------------------------------------------
    nodal = []
    if loads:
        lname = os.path.basename(loads)
        for line, c in read_rows(loads):
            _need(c, 4, lname, line, "node, Fx, Fy, Fz[, case]")
            nd = node_of(c[0], lname, line)
            case = (c[4].strip().upper() if len(c) > 4 and c[4].strip()
                    else "L")
            if case not in CASES:
                raise CSVImportError(f"{lname}, line {line}: case {c[4]!r} "
                                     f"must be D, S, L or W.")
            nodal.append({"node": nd,
                          "Fx": _num(c[1], lname, line, "Fx") * KN,
                          "Fy": _num(c[2], lname, line, "Fy") * KN,
                          "Fz": _num(c[3], lname, line, "Fz") * KN,
                          "case": case})

    # ---- model ------------------------------------------------------------
    try:
        g = geo.custom_geometry(xyz, bars, faces=tri,
                                base_nodes=sorted(custom) if custom else None,
                                groups=grp, title=title)
    except ValueError as exc:
        raise CSVImportError(f"The truss is not valid: {exc}") from None
    notes.extend(g.notes)
    model = DomeModel(
        geom=g,
        support_type=geo.CUSTOM if custom else geo.PINNED,
        loads=LoadSpec(include_self_weight=True, dead_pressure=0.0,
                       snow_pressure=0.0, live_pressure=0.0,
                       wind_enabled=False, nodal_loads=nodal),
        combinations=import_combinations(),
        custom_supports=custom)
    if not custom:
        notes.append(f"No supports file: the {len(g.base_nodes)} lowest "
                     f"node(s) are pinned.")
    try:
        fixed, springs = model.support_conditions()
    except ValueError as exc:
        raise CSVImportError(str(exc)) from None
    if len(fixed) + len(springs) < 3:
        notes.append("Warning: fewer than three restrained directions; the "
                     "truss cannot be stable.")
    notes.append(f"{len(nodal)} nodal load(s) read.")
    return model, notes


def import_truss_csv(path, title=None):
    """Find the files from ``nodes.csv`` (or its folder) and read them.

    The title defaults to the folder name (and the file prefix, if any).
    Returns ``(model, notes)``.
    """
    files = find_files(path)
    if title is None:
        folder = os.path.basename(os.path.dirname(files["nodes"])) or "truss"
        stem = os.path.basename(files["nodes"])[:-len("nodes.csv")]
        stem = stem.rstrip("_- .")
        title = f"Imported truss ({stem or folder})"
    return read_truss_csv(files["nodes"], files["members"],
                          files["supports"], files["loads"], files["faces"],
                          title=title)
