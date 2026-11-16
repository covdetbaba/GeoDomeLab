"""
Tests for the extended validation (validation/run_validation_ext.py).

The comparison runs every F1-F10 hand calculation of
validation/hand_calc_ext.py through GeoDomeLab.  These tests require that

* the hand files are independent of the program (no geodome import, checked
  on their source code with the ast module);
* every hand reference value is compared, and every comparison passes at its
  stated tolerance (linear 1e-9, modal 1e-9, nonlinear 1e-6, counts exact,
  published chord factors 5e-6, identities 1e-10);
* the program identities (equilibrium, energy, load balance, modal
  orthonormality) hold for every solve used;
* the von Mises path really goes through both limit points and the
  snap-through target;
* the comparison is not vacuous: a wrong value or a wrong material is caught;
* the CSV / markdown outputs are complete and read back exactly.

Run with pytest or standalone::

    python tests/test_validation_ext.py
"""

from __future__ import annotations

import ast
import csv
import math
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAL = os.path.join(ROOT, "validation")
for p in (VAL, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

import run_validation_ext as ve                                    # noqa: E402

_CACHE = {}


def _rows():
    if "rows" not in _CACHE:
        _CACHE["rows"], _CACHE["info"] = ve.run_comparison()
    return _CACHE["rows"], _CACHE["info"]


def _imports(path):
    tree = ast.parse(open(path, encoding="utf-8").read(), filename=path)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            names.add(node.module or "")
    return names


# ============================================================ independence ==
def test_hand_files_do_not_import_geodome():
    for name in ("hand_calc.py", "hand_calc_ext.py"):
        mods = _imports(os.path.join(VAL, name))
        assert not any(m.split(".")[0] == "geodome" for m in mods), (name,
                                                                    mods)
    # and the hand modules hold no reference to a geodome object
    for mod in (ve.hc, ve.hx):
        for v in vars(mod).values():
            assert not getattr(v, "__module__", "").startswith("geodome")
            if hasattr(v, "__name__") and type(v).__name__ == "module":
                assert not v.__name__.startswith("geodome")


def test_hand_internal_checks_all_pass():
    chk = ve.hx.checks()
    assert len(chk) >= 90
    bad = [c for c in chk if not c[2] <= c[3]]
    assert not bad, bad


# ============================================================== comparison ==
def test_every_hand_value_is_compared():
    rows, info = _rows()
    assert info["missing"] == [], info["missing"]
    assert info["n_used"] == info["n_hand"] == len(ve.hx.results())
    keys = {r["key"] for r in rows}
    for case in ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
                 "F10"):
        assert any(k.startswith(case + ".") for k in keys), case
    # every row is filed under its own case (F10 is not F1)
    assert {r["case"] for r in rows} == {f"F{k}" for k in range(1, 11)}


def test_all_comparisons_pass():
    rows, _ = _rows()
    bad = [(r["sub"], r["quantity"], r["hand"], r["program"], r["difference"])
           for r in rows if not r["pass"]]
    assert not bad, bad
    assert len(rows) > 300


def test_tolerances_by_kind():
    rows, _ = _rows()
    limits = {"linear": 1e-9, "modal": 1e-9, "nonlinear": 1e-6,
              "count": 0.0, "published": 5e-6, "identity": 1e-10}
    for r in rows:
        assert r["tolerance"] == limits[r["kind"]], r
        assert math.isfinite(r["difference"]), r
        assert r["difference"] <= limits[r["kind"]], r
    kinds = {r["kind"] for r in rows}
    assert kinds == set(limits), kinds


def test_nonzero_references_use_plain_relative_difference():
    rows, _ = _rows()
    for r in rows:
        if r["kind"] in ("linear", "modal", "nonlinear") and r["hand"] != 0:
            assert r["scale"] == abs(r["hand"])
            assert abs(r["difference"] - r["rel_diff"]) <= \
                1e-15 * max(1.0, r["rel_diff"])


def test_program_identities_for_every_solve():
    rows, _ = _rows()
    ident = [r for r in rows if r["kind"] == "identity"]
    subs = {r["sub"] for r in ident}
    want = {"F1a", "F1b", "F2a", "F2b", "F3a", "F3b", "F4 load",
            "F4 temperature", "F4 settlement", "F4 combined", "F5a", "F5b",
            "F5c", "F6 linear", "F7 sw", "F7 dead", "F9 wind Model B"}
    assert subs == want, sorted(subs ^ want)
    assert all(r["program"] <= 1e-10 for r in ident)


# ==================================================== specific cases ======
def test_wind_case_f9_matches_independent_numbers():
    """Review V02: the wind field on Model B (default wind data) against
    numbers worked out independently of both files (qz = 0.613 Kz Kzt Kd
    V^2, Cp interpolated along the meridian, F = -qz G Cp A_trib n)."""
    h = ve.HAND
    assert abs(h["F9.qz"]["value"] - 1061.332875) < 1e-9
    assert abs(h["F9.Fz_1"]["value"] - 17994.811) < 1e-3       # crown
    for key, want in (("F9.Fx_2", 7863.239), ("F9.Fz_2", 3931.619),
                      ("F9.Fx_3", 4441.627), ("F9.Fy_3", 13669.924),
                      ("F9.Fz_3", 7186.704), ("F9.Fx_4", -2134.404),
                      ("F9.Fy_4", 1550.736), ("F9.Fz_4", 1319.135)):
        assert abs(h[key]["value"] - want) < 1e-3, (key, h[key]["value"])
    assert abs(h["F9.Cp_1"]["value"] + 1.0) < 1e-15
    rows, _ = _rows()
    f9 = [r for r in rows if r["case"] == "F9" and r["kind"] == "linear"]
    assert len(f9) == 49 and all(r["pass"] for r in f9)


def test_aisc_case_f10_matches_independent_numbers():
    """Review V06: inelastic E3, B4.2 design wall and E7 effective area
    (AISC 360-16), worked out independently: 239.578, 223.996 and
    1112.935 kN; zero above the 0.45 E/Fy scope limit."""
    h = ve.HAND
    for key, want in (("F10.phiPn_inel", 239578.0), ("F10.phiPn_b42",
                                                      223996.0),
                      ("F10.phiPn_e7", 1112935.0)):
        assert abs(h[key]["value"] - want) < 1.0, (key, h[key]["value"])
    assert abs(h["F10.Q_e7"]["value"] - 0.93871) < 1e-5
    assert h["F10.phiPn_scope"]["value"] == 0.0
    rows, _ = _rows()
    f10 = [r for r in rows if r["case"] == "F10"]
    assert len(f10) == 14 and all(r["pass"] for r in f10)
    assert max(r["difference"] for r in f10) <= 1e-12

def test_von_mises_path_covers_both_limits_and_snap():
    rows, _ = _rows()
    hand = ve.HAND
    path = [r for r in rows if r["key"] == "path"]
    assert len(path) == 2 and all(r["pass"] for r in path)
    # the traced path reached the inverted branch beyond the snap target
    assert path[1]["program"] == 1.0
    lim = [r for r in rows if r["key"] == "F6.P_limit"]
    assert lim and abs(lim[0]["program"] - hand["F6.P_limit"]["value"]) \
        <= 1e-6 * hand["F6.P_limit"]["value"]
    w2 = [r for r in rows if r["key"] == "F6.w_second_limit"]
    assert w2 and w2[0]["difference"] < 1e-9
    snap = [r for r in rows if r["key"] == "F6.w_snap"]
    assert snap and snap[0]["program"] > 2.0 * ve.hx.VM_H


def test_self_straining_cases_use_the_program_inputs():
    """The heated / settled cases really carry their load: a zero program
    answer would not match the non-zero hand forces."""
    rows, _ = _rows()
    for key in ("F1.N", "F2.N_ring", "F3.vert_N1", "F3.spread_N6",
                "F4.N2_temperature", "F4.N2_settlement", "F5.thermal_N"):
        sel = [r for r in rows if r["key"] == key]
        assert sel, key
        for r in sel:
            assert abs(r["hand"]) > 1e3 and r["pass"], r


def test_full_sphere_counts_and_kenner():
    rows, _ = _rows()
    want = {"F8.nodes": 42, "F8.struts": 120, "F8.panels": 80,
            "F8.hubs_valence5": 12, "F8.hubs_valence6": 30}
    for k, v in want.items():
        sel = [r for r in rows if r["key"] == k]
        assert sel and all(r["program"] == v for r in sel), k
    ken = [r for r in rows if r["kind"] == "published"]
    assert len(ken) == 4 and all(r["pass"] for r in ken)


# ================================================== not vacuous ===========
def test_comparison_catches_a_wrong_value():
    C = ve.Comparison(ve.HAND)
    h = ve.HAND["F1.N"]["value"]
    C.add("F1a", "F1.N", "right", h * (1 + 1e-12), "test")
    C.add("F1a", "F1.N", "wrong by 1e-8", h * (1 + 1e-8), "test")
    C.add("F1a", "F1.R1_x", "zero-free check", 0.0, "test")
    C.add("F6 nonlinear", "F6.P_limit", "nonlinear off by 2e-6",
          ve.HAND["F6.P_limit"]["value"] * (1 + 2e-6), "test",
          kind="nonlinear")
    rows = C.evaluate()
    assert [r["pass"] for r in rows] == [True, False, False, False]


def test_comparison_catches_a_wrong_material():
    """F1 run with EN S355 (alpha = 12e-6, E = 210 GPa) must fail."""
    d = ve.DEFS["F1a"]
    m = ve.custom_model(d["nodes"], d["members"], [ve.SEC_MAIN],
                        {0: ve.PIN, 1: ve.PIN},
                        loads=ve._loads(temperature_change=d["dT"]),
                        material="EN S355")
    r = ve.analyse(m)
    C = ve.Comparison(ve.HAND)
    C.add("F1a", "F1.N", "N with the wrong material", r.N[0], "test")
    assert not C.evaluate()[0]["pass"]


# =================================================== outputs ==============
def test_outputs_complete_and_exact():
    rows, info = _rows()
    with tempfile.TemporaryDirectory() as tmp:
        pc = ve.write_csv(rows, os.path.join(tmp, "c.csv"))
        pm = ve.write_markdown(rows, info, os.path.join(tmp, "c.md"))
        with open(pc, encoding="utf-8") as fh:
            back = list(csv.DictReader(fh))
        assert len(back) == len(rows)
        for a, b in zip(rows, back):
            assert float(b["program"]) == a["program"]
            assert float(b["hand"]) == a["hand"]
            assert b["pass"] == ("yes" if a["pass"] else "NO")
        text = open(pm, encoding="utf-8").read()
        assert f"{len(rows)} of {len(rows)} comparisons pass" in text
        assert "## Failures" in text and "None." in text
        for case in ve.CASE_TITLES:
            assert f"| {case} |" in text


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
    import matplotlib
    matplotlib.use("Agg")
    raise SystemExit(_run_standalone())
