"""
Tests for geodome/cad.py, the ASCII DXF R12 export.

The file is read back with a small DXF reader written in this test (no
third-party package), and checked against the model:

* the file is pure ASCII, a clean list of (group code, value) pairs, with the
  sections HEADER, TABLES, BLOCKS, ENTITIES and a final EOF;
* $ACADVER is AC1009 (R12) and $INSUNITS is 6 (metres);
* the LAYER table holds "0", STRUT_<family> for every family, PANELS and HUBS,
  each with a valid colour, and every entity sits on a layer of the table;
* one LINE per strut (on its family layer), one 3DFACE per panel (a triangle
  repeats its third corner), one POINT per node;
* every coordinate reads back to within 1e-9 m;
* a custom truss without faces writes no 3DFACE;
* with analysis results every LINE carries its utilisation colour.

If ezdxf happens to be installed, the file is also opened with it.

Run with pytest::

    pytest -q tests/test_cad.py

or standalone::

    python tests/test_cad.py
"""

from __future__ import annotations

import atexit
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import benchmarks as bm                                # noqa: E402
from geodome import cad                                             # noqa: E402
from geodome import fabrication as fab                              # noqa: E402
from geodome.model import default_model                             # noqa: E402

TOL = 1e-9          # metres

_TMP = tempfile.mkdtemp(prefix="geodome_cad_test_")
atexit.register(shutil.rmtree, _TMP, ignore_errors=True)


# ============================================================ DXF reader ====
def read_dxf(path):
    """Parse an ASCII DXF file into a small dict.

    Returns ``{"sections": [names], "header": {var: [(code, value)]},
    "layers": {name: {"colour": int, "linetype": str, "flags": int}},
    "layer_count": int (group 70 of the LAYER table), "entities": [dict]}``.
    Each entity dict maps "type" to the entity name and every group code to
    its value (codes are unique inside the entities written by cad.py).
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    text = raw.decode("ascii")                        # raises if not ASCII
    lines = text.splitlines()
    assert len(lines) % 2 == 0, "odd number of lines: broken group pairs"
    pairs = []
    for i in range(0, len(lines), 2):
        code = int(lines[i].strip())                  # raises if not a code
        pairs.append((code, lines[i + 1].strip()))
    assert pairs[-1] == (0, "EOF"), "file does not end with EOF"

    out = {"sections": [], "header": {}, "layers": {}, "layer_count": None,
           "entities": [], "pairs": pairs, "crlf": b"\r\n" in raw}
    i = 0
    while i < len(pairs):
        code, val = pairs[i]
        if (code, val) == (0, "EOF"):
            assert i == len(pairs) - 1, "data after EOF"
            break
        assert (code, val) == (0, "SECTION"), f"expected SECTION at pair {i}"
        assert pairs[i + 1][0] == 2
        name = pairs[i + 1][1]
        out["sections"].append(name)
        i += 2
        body = []
        while pairs[i] != (0, "ENDSEC"):
            body.append(pairs[i])
            i += 1
        i += 1
        if name == "HEADER":
            var = None
            for c, v in body:
                if c == 9:
                    var = v
                    out["header"][var] = []
                else:
                    out["header"][var].append((c, v))
        elif name == "TABLES":
            _read_tables(body, out)
        elif name == "ENTITIES":
            ent = None
            for c, v in body:
                if c == 0:
                    ent = {"type": v}
                    out["entities"].append(ent)
                else:
                    assert c not in ent, f"repeated code {c} in {ent['type']}"
                    ent[c] = v
    return out


def _read_tables(body, out):
    table = None
    rec = None
    for c, v in body:
        if c == 0 and v == "TABLE":
            table, rec = "?", None
        elif c == 0 and v == "ENDTAB":
            table, rec = None, None
        elif c == 0:
            rec = {"type": v}
            if v == "LAYER":
                rec["flags"] = None
        elif rec is None and c == 2:
            table = v
        elif rec is None and c == 70 and table == "LAYER":
            out["layer_count"] = int(v)
        elif rec is not None and rec["type"] == "LAYER":
            if c == 2:
                rec["name"] = v
                out["layers"][v] = rec
            elif c == 62:
                rec["colour"] = int(v)
            elif c == 6:
                rec["linetype"] = v
            elif c == 70:
                rec["flags"] = int(v)


def _xyz(ent, base):
    return np.array([float(ent[10 + base]), float(ent[20 + base]),
                     float(ent[30 + base])])


def _by_type(dxf, kind):
    return [e for e in dxf["entities"] if e["type"] == kind]


def _write(model, results=None, name="model.dxf"):
    path = os.path.join(_TMP, name)
    assert cad.export_dxf(model, path, results=results) == path
    return path


def _check_file_against_model(path, model):
    g = model.geom if hasattr(model, "geom") else model
    dxf = read_dxf(path)
    assert dxf["sections"] == ["HEADER", "TABLES", "BLOCKS", "ENTITIES"]
    assert dxf["crlf"]

    # header
    assert dxf["header"]["$ACADVER"] == [(1, "AC1009")]
    assert dxf["header"]["$INSUNITS"] == [(70, "6")]
    ext_lo = np.array([float(v) for _, v in dxf["header"]["$EXTMIN"]])
    ext_hi = np.array([float(v) for _, v in dxf["header"]["$EXTMAX"]])
    assert np.allclose(ext_lo, g.nodes.min(axis=0), atol=TOL, rtol=0)
    assert np.allclose(ext_hi, g.nodes.max(axis=0), atol=TOL, rtol=0)

    # layer table
    nfam = g.n_strut_types()
    expect = ["0"] + [f"STRUT_{s['label']}" for s in g.strut_schedule()] \
        + ["PANELS", "HUBS"]
    assert list(dxf["layers"]) == expect, list(dxf["layers"])
    assert dxf["layer_count"] == len(expect)
    for name, rec in dxf["layers"].items():
        assert 1 <= rec["colour"] <= 255, (name, rec)
        assert rec["linetype"] == "CONTINUOUS"
        assert rec["flags"] == 0
    # the strut families have different colours (up to 12 families)
    fam_cols = [dxf["layers"][f"STRUT_{s['label']}"]["colour"]
                for s in g.strut_schedule()]
    assert len(set(fam_cols[:12])) == min(nfam, 12)
    for e in dxf["entities"]:
        assert e[8] in dxf["layers"], e

    # entity counts
    lines = _by_type(dxf, "LINE")
    faces = _by_type(dxf, "3DFACE")
    points = _by_type(dxf, "POINT")
    assert len(lines) == g.n_members
    assert len(faces) == g.n_faces
    assert len(points) == g.n_nodes
    assert len(dxf["entities"]) == g.n_members + g.n_faces + g.n_nodes

    # per-layer strut counts = strut schedule counts
    for s in g.strut_schedule():
        n = sum(1 for e in lines if e[8] == f"STRUT_{s['label']}")
        assert n == s["count"], (s["label"], n, s["count"])

    # coordinates round trip, entity k = member / face / node k
    worst = 0.0
    for k, (i, j) in enumerate(g.members):
        e = lines[k]
        assert e[8] == f"STRUT_{fab.family_labels(nfam)[g.strut_group[k]]}"
        worst = max(worst, np.max(np.abs(_xyz(e, 0) - g.nodes[i])),
                    np.max(np.abs(_xyz(e, 1) - g.nodes[j])))
    for k, (a, b, c) in enumerate(g.faces):
        e = faces[k]
        assert e[8] == "PANELS"
        corners = [_xyz(e, q) for q in range(4)]
        assert np.array_equal(corners[2], corners[3]), "triangle: 4th = 3rd"
        for q, node in enumerate((a, b, c)):
            worst = max(worst, np.max(np.abs(corners[q] - g.nodes[node])))
    for k, p in enumerate(g.nodes):
        e = points[k]
        assert e[8] == "HUBS"
        worst = max(worst, np.max(np.abs(_xyz(e, 0) - p)))
    assert worst < TOL, worst
    return dxf


# ================================================================= tests ====
def test_2v_dome_counts_layers_and_coordinates():
    """2V dome: 65 LINE, 40 3DFACE, 26 POINT, two strut layers."""
    m = default_model(2, 5.0, 0.5)
    dxf = _check_file_against_model(_write(m), m)
    assert len(_by_type(dxf, "LINE")) == 65
    assert len(_by_type(dxf, "3DFACE")) == 40
    assert len(_by_type(dxf, "POINT")) == 26
    assert dxf["layers"]["STRUT_A"]["colour"] == cad.FAMILY_ACI[0]   # red
    assert dxf["layers"]["STRUT_B"]["colour"] == cad.FAMILY_ACI[1]   # blue
    # no colour override without results
    assert all(62 not in e for e in dxf["entities"])


def test_6v_case_study_round_trip():
    """The 6V case study (9 families, 555 struts) reads back exactly."""
    m = bm.make_preset("case_study_6v")
    dxf = _check_file_against_model(_write(m), m)
    assert len(dxf["layers"]) == 1 + 9 + 2


def test_large_coordinates_round_trip():
    """A 40 m 10V dome (27 families) still round-trips to 1e-9 m."""
    m = default_model(10, 40.0, 0.5)
    _check_file_against_model(_write(m), m)


def test_custom_truss_without_faces():
    """Custom trusses have no panels: no 3DFACE, NaN radius is harmless."""
    for key in ("three_bar", "von_mises"):
        m = bm.make_preset(key)
        assert m.geom.is_custom and m.geom.n_faces == 0
        dxf = _check_file_against_model(_write(m, name=f"{key}.dxf"), m)
        assert _by_type(dxf, "3DFACE") == []
        assert len(_by_type(dxf, "LINE")) == m.geom.n_members
        assert len(_by_type(dxf, "POINT")) == m.geom.n_nodes


def test_custom_truss_with_faces():
    """The star dome is custom but has 12 faces: 12 3DFACE entities."""
    m = bm.make_preset("star_dome")
    dxf = _check_file_against_model(_write(m), m)
    assert len(_by_type(dxf, "3DFACE")) == 12


def test_bare_geometry_is_accepted():
    m = default_model(3, 6.0, 0.5)
    _check_file_against_model(_write(m.geom), m.geom)


def test_utilisation_colours():
    """With results every LINE gets the ACI colour of its utilisation band."""
    from geodome.postprocess import run_analysis
    m = bm.make_preset("validation_b")
    res = run_analysis(m, do_buckling=False)
    dxf = _check_file_against_model(_write(m, results=res), m)
    lines = _by_type(dxf, "LINE")
    expect = cad.utilisation_aci(res.util_env)
    got = np.array([int(e[62]) for e in lines])
    assert np.array_equal(got, expect)
    # points and faces keep the layer colour
    assert all(62 not in e for e in dxf["entities"] if e["type"] != "LINE")
    # a plain array works too, and a wrong length is refused clearly
    path = _write(m, results=np.full(m.geom.n_members, 1.2), name="u.dxf")
    lines = _by_type(read_dxf(path), "LINE")
    assert {int(e[62]) for e in lines} == {cad.UTIL_ACI[-1]}
    try:
        cad.dxf_text(m, results=np.zeros(3))
    except ValueError as exc:
        assert "Run the analysis again" in str(exc)
    else:
        raise AssertionError("wrong-length utilisation was accepted")


def test_utilisation_bands():
    u = [0.0, 0.25, 0.2501, 0.5, 0.74, 0.9, 0.95, 1.0, 1.0001, 3.0, -0.1]
    got = list(cad.utilisation_aci(u))
    b = cad.UTIL_ACI
    assert got == [b[0], b[0], b[1], b[1], b[2], b[3], b[4], b[4], b[5],
                   b[5], b[0]], got


def test_number_format_has_no_exponent():
    """Tiny and negative-zero values are written in fixed notation."""
    assert cad._num(1e-17) == "0.0"
    assert cad._num(-0.0) == "0.0"
    assert cad._num(-1e-13) == "0.0"
    assert cad._num(2.5) == "2.5"
    assert cad._num(-12.000000000001) == "-12.000000000001"
    assert cad._num(7654321.25) == "7654321.25"
    for bad in (float("nan"), float("inf")):
        try:
            cad._num(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("non-finite value accepted")


def test_layer_names_beyond_z():
    """Families 26 and up use the S26 labels of geometry.strut_schedule."""
    assert cad.strut_layer(0) == "STRUT_A"
    assert cad.strut_layer(25) == "STRUT_Z"
    assert cad.strut_layer(26) == "STRUT_S26"


def test_ezdxf_reads_the_file_if_installed():
    """Optional: open the file with ezdxf when it is installed."""
    try:
        import ezdxf                                   # noqa: F401
    except Exception:
        print("        (ezdxf not installed: skipped)")
        return
    m = default_model(2, 5.0, 0.5)
    doc = ezdxf.readfile(_write(m))
    msp = doc.modelspace()
    assert len(msp.query("LINE")) == m.geom.n_members
    assert len(msp.query("3DFACE")) == m.geom.n_faces
    assert len(msp.query("POINT")) == m.geom.n_nodes
    assert doc.header.get("$INSUNITS") == 6
    for name in ("STRUT_A", "STRUT_B", "PANELS", "HUBS"):
        assert doc.layers.has_entry(name)


def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = 0
    for name, fn in fns:
        try:
            fn()
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        else:
            npass += 1
            print(f"  pass  {name}")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
