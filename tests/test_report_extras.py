"""
Tests for the optional sections of the HTML report (``report.export_html``
with ``extras=``) and for the Theory help page (``geodome/gui/theory.py``).

Every extra is built with the real core functions on small models:

* nonlinear path      ``nonlinear.run_nonlinear`` on the von Mises preset and
                      on a 2V dome (R = 6 m, hemisphere);
* imperfection study  ``nonlinear.imperfection_study`` on the 2V dome;
* modal analysis      ``modal.modal_analysis`` on validation Model A and the
                      2V dome;
* fabrication         ``fabrication.fabrication_tables`` on the 2V dome and
                      on the von Mises truss (no panels, not spherical);
* parameter study     ``study.parameter_study`` over 2V and 3V.

Checks: the HTML is well formed (``html.parser``: every tag closed in order,
no duplicate ids), each new section appears exactly when its extra is given
and in the fixed order 8-12, sections 1-7 are unchanged byte for byte (after
masking the time stamp), the key numbers in each section equal the result
objects, bad or empty extras never break the report, and the 2V fabrication
section shows the closed-form chord factors 2 sin(acos(1/sqrt5)/4) and
2 sin(18 deg) and axial angles acos(cf/2).

Run with pytest::

    pytest -q tests/test_report_extras.py

or standalone, without pytest installed::

    python tests/test_report_extras.py
"""

from __future__ import annotations

import importlib.util
import math
import os
import re
import sys
import tempfile
from html.parser import HTMLParser

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

from geodome import benchmarks as bm                                # noqa: E402
from geodome import fabrication as fb                               # noqa: E402
from geodome import modal as md                                     # noqa: E402
from geodome import nonlinear as nl                                 # noqa: E402
from geodome import report as rp                                    # noqa: E402
from geodome import study as st                                     # noqa: E402
from geodome.model import default_model                             # noqa: E402
from geodome.postprocess import run_analysis                        # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OLD_H2 = ["1   Structure", "2   Loading",
          "3   Analysis results", "5   Figures",
          "6   Solution verification"]
NEW_IDS = ["nonlinear", "imperfection", "modal", "fabrication", "study"]
FOOTER = ("<footer>Produced by GeoDomeLab. Linear elastic first-order "
          "space-truss analysis; member checks to AISC 360 LRFD; "
          "wind pressures from a simplified domed-roof coefficient "
          "distribution. Results are for academic use.")


# ================================================================ helpers ==
class _Checker(HTMLParser):
    """Stack-based well-formedness check plus a list of the headings."""

    VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input",
            "link", "meta", "param", "source", "track", "wbr"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack, self.errors, self.ids = [], [], []
        self.headings = []            # [(tag, id or None, text)]
        self._head = None

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if d.get("id") is not None:
            self.ids.append(d["id"])
        if tag in ("h2", "h3"):
            self._head = [tag, d.get("id"), ""]
        if tag not in self.VOID:
            self.stack.append(tag)

    def handle_startendtag(self, tag, attrs):
        if tag not in self.VOID:
            self.errors.append(f"self-closed <{tag}/>")

    def handle_endtag(self, tag):
        if tag in self.VOID:
            self.errors.append(f"end tag for void <{tag}>")
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append(f"</{tag}> closes <{self.stack[-1] if self.stack else None}>")
            if tag in self.stack:                # resynchronise
                while self.stack and self.stack.pop() != tag:
                    pass
            return
        self.stack.pop()
        if tag in ("h2", "h3") and self._head is not None:
            self.headings.append(tuple(self._head))
            self._head = None

    def handle_data(self, data):
        if self._head is not None:
            self._head[2] += data


def _parse(text):
    c = _Checker()
    c.feed(text)
    c.close()
    return c


def _assert_well_formed(text):
    assert text.startswith("<!doctype html><html>"), text[:40]
    assert text.endswith("</div></body></html>")
    c = _parse(text)
    assert not c.errors, c.errors[:5]
    assert not c.stack, f"unclosed tags {c.stack}"
    assert len(c.ids) == len(set(c.ids)), f"duplicate ids {c.ids}"
    return c


def _mask(text):
    return re.sub(r"generated \d{4}-\d{2}-\d{2} \d{2}:\d{2}", "generated T",
                  text)


def _export(res, extras="__absent__", **kw):
    kw.setdefault("include_figures", False)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "r.html")
        if extras == "__absent__":
            rp.export_html(res, p, **kw)
        else:
            rp.export_html(res, p, extras=extras, **kw)
        with open(p, encoding="utf-8") as fh:
            return _mask(fh.read())


def _section(text, sid):
    """HTML of the section whose heading has id ``sid`` (to the next h2 or
    the footer).  The whole page is checked for well-formedness first."""
    _assert_well_formed(text)
    i = text.index(f"<h2 id='{sid}'>")
    j = text.find("<h2", i + 4)
    if j < 0:
        j = text.index("<footer>")
    return text[i:j]


def _old_part(text):
    """Everything before the first optional section (or the footer)."""
    i = text.find("<h2 id=")
    return text[:i] if i >= 0 else text[:text.index("<footer>")]


_FX = {}


def _fx():
    """Models, results and extras, built once with the real core code."""
    if _FX:
        return _FX
    m2 = default_model(2, 6.0, 0.5)
    vm = bm.make_preset("von_mises")
    va = bm.make_preset("validation_a")
    _FX.update(
        m2=m2, res2=run_analysis(m2), vm=vm, res_vm=run_analysis(vm),
        va=va, res_va=run_analysis(va),
        nl2=nl.run_nonlinear(m2), nl_vm=nl.run_nonlinear(vm),
        imp=nl.imperfection_study(m2, return_results=True),
        modal_va=md.modal_analysis(va), modal2=md.modal_analysis(m2),
        fab2=fb.fabrication_tables(m2, hub_offset=0.05, stock_length=6.0),
        fab_vm=fb.fabrication_tables(vm),
        study=st.parameter_study(m2, [2, 3]),
    )
    _FX["extras2"] = dict(
        nonlinear_result=_FX["nl2"], imperfection_result=_FX["imp"][0],
        modal_result=_FX["modal2"], fabrication_tables=_FX["fab2"],
        study_rows=_FX["study"])
    return _FX


# ============================================================ structure ===
def test_no_extras_keeps_the_old_report():
    fx = _fx()
    base = _export(fx["res2"], include_figures=True)
    c = _assert_well_formed(base)
    h2 = [(i, t) for tag, i, t in c.headings if tag == "h2"]
    assert h2 == [(None, t) for t in OLD_H2], h2
    assert c.ids == []                               # no ids before or now
    assert base[base.index("<footer>"):] == FOOTER + "</footer></div></body></html>"
    for key in ("extras=None", "{}", "all None", "all empty", "unknown"):
        extras = {"extras=None": None, "{}": {},
                  "all None": {k: None for k in rp.EXTRA_KEYS},
                  "all empty": {"imperfection_result": [],
                                "fabrication_tables": {}, "study_rows": []},
                  "unknown": {"foo": fx["nl2"], "nonlinear": fx["nl2"]}}[key]
        assert _export(fx["res2"], extras, include_figures=True) == base, key


def test_each_extra_alone():
    fx = _fx()
    base = _export(fx["res2"])
    # sections 1-7 with figures on, to compare against the figure-on exports
    old = _export(fx["res2"], include_figures=True)
    n_img = {"nonlinear_result": 1, "imperfection_result": 1,
             "modal_result": 1, "fabrication_tables": 0, "study_rows": 1}
    for key, sid in zip(rp.EXTRA_KEYS, NEW_IDS):
        text = _export(fx["res2"], {key: fx["extras2"][key]},
                       include_figures=True)
        c = _assert_well_formed(text)
        new = [(i, t) for tag, i, t in c.headings if tag == "h2" and i]
        num = rp.EXTRA_SECTIONS[key][1]
        assert new == [(sid, f"{num}   {rp.EXTRA_SECTIONS[key][2]}")], new
        assert c.ids == [sid]
        sec = _section(text, sid)
        assert "could not be built" not in sec, sec[:400]
        assert sec.count("<img ") == n_img[key], (key, sec.count("<img "))
        assert sec.count("Figure ") == n_img[key]
        assert _old_part(text) == old[:old.index("<footer>")], key
        foot = text[text.index("<footer>"):]
        assert foot.startswith(FOOTER + " Extra sections: "
                               f"{num} "), foot[:300]
        assert foot.endswith(".</footer></div></body></html>")
        # without figures the section has no image at all
        lean = _export(fx["res2"], {key: fx["extras2"][key]})
        assert "<img " not in _section(lean, sid)
        assert _old_part(lean) == base[:base.index("<footer>")]


def test_all_extras_in_fixed_order():
    fx = _fx()
    # pass them in reverse order: the report order must not change
    extras = dict(reversed(list(fx["extras2"].items())))
    text = _export(fx["res2"], extras, include_figures=True)
    c = _assert_well_formed(text)
    h2 = [(i, t) for tag, i, t in c.headings if tag == "h2"]
    assert [t for i, t in h2 if i is None] == OLD_H2
    assert [i for i, t in h2 if i] == NEW_IDS
    nums = [int(t.split()[0]) for i, t in h2 if i]
    assert nums == [8, 9, 10, 11, 12]
    for sid in NEW_IDS:
        assert text.count(f"<h2 id='{sid}'>") == 1
        assert "could not be built" not in _section(text, sid)
    assert text.count("<img ") == 8 + 4             # 8 old + 4 new figures
    assert "Extra sections: 8 " in text and "; 12 " in text


def test_extras_with_old_optional_sections():
    """Sizing (4) and benchmark (7) sections stay where they were."""
    fx = _fx()
    from geodome.optimize import auto_size
    m = default_model(2, 6.0, 0.5)
    sizing = auto_size(m, target_util=0.95, max_iter=6)
    res = run_analysis(m)
    ver = [{"group": "g", "check": "c", "computed": 1.0, "reference": 1.0,
            "unit": "-", "rel_error": 0.0, "tolerance": 1e-9, "pass": True,
            "note": "identity"}]
    old = _export(res, sizing=sizing, verification=ver)
    new = _export(res, {"modal_result": md.modal_analysis(m)},
                  sizing=sizing, verification=ver)
    _assert_well_formed(new)
    assert _old_part(new) == old[:old.index("<footer>")]
    h2 = [t for tag, i, t in _parse(new).headings if tag == "h2"]
    assert [t.split()[0] for t in h2] == ["1", "2", "3", "4", "6", "7", "10"]
    assert fx is not None


# ============================================================== content ===
def _plain(html_text):
    text = re.sub(r"<[^>]+>", " ", html_text)
    import html as _html
    return re.sub(r"\s+", " ", _html.unescape(text))


def test_nonlinear_section_numbers_and_figure():
    fx = _fx()
    r = fx["nl_vm"]
    assert r.limit_found and math.isfinite(r.lambda_cr)
    # the result itself: von Mises closed form (benchmarks reference)
    a, h, EA, P_ref = bm.von_mises_parameters(fx["vm"])
    lim = bm.von_mises_limit_point(a, h, EA)
    assert abs(r.lambda_limit * P_ref - lim["P"]) < 1e-6 * lim["P"]

    text = _export(fx["res_vm"], {"nonlinear_result": r})
    sec = _plain(_section(text, "nonlinear"))
    assert f"{r.lambda_limit:.6g} (limit point)" in sec, sec
    assert f"Linear buckling factor λ_cr {r.lambda_cr:.6g}" in sec
    assert f"λ_limit / λ_cr {r.limit_ratio:.4f}" in sec
    assert "Arc-length (Crisfield, cylindrical)" in sec
    assert r.control_label in sec
    for m in r.messages:
        assert _plain(m) in sec

    fig = rp.nonlinear_path_figure(r)
    ax = fig.axes[0]
    sgn = -1.0                                   # the apex moves down
    star = [ln for ln in ax.lines if ln.get_marker() == "*"]
    assert len(star) == 1
    x, y = star[0].get_xdata()[0], star[0].get_ydata()[0]
    assert abs(x - sgn * r.control_limit * 1e3) < 1e-9
    assert abs(y - r.lambda_limit) < 1e-12
    hlines = [ln for ln in ax.lines if len(ln.get_ydata()) == 2
              and ln.get_ydata()[0] == ln.get_ydata()[1]
              and abs(ln.get_ydata()[0] - r.lambda_cr) < 1e-12]
    assert len(hlines) == 1                      # lambda_cr line
    labels = [ln.get_label() for ln in ax.lines]
    assert any(lb.startswith("Corotational") for lb in labels)
    assert any(lb.startswith("Linear") for lb in labels)
    assert any(lb.startswith("P-Delta") for lb in labels)
    assert ax.get_ylim()[1] >= r.lambda_cr       # lambda_cr is in view
    corot = next(ln for ln in ax.lines if ln.get_label().startswith("Corot"))
    assert np.allclose(corot.get_xdata(), sgn * r.control * 1e3)
    assert np.allclose(corot.get_ydata(), r.lam)


def test_nonlinear_bifurcation_rows_on_the_dome():
    fx = _fx()
    r = fx["nl2"]
    sec = _section(_export(fx["res2"], {"nonlinear_result": r}), "nonlinear")
    rows = re.findall(r"<tr>", sec)
    # header row + one row per critical point in the table
    assert len(rows) == 1 + len(r.critical_points)
    if r.lambda_first_critical < r.lambda_limit:
        assert "bifurcation before the load maximum" in _plain(sec)


def test_imperfection_section_list_and_tuple_forms():
    fx = _fx()
    rows, results = fx["imp"]
    a = _section(_export(fx["res2"], {"imperfection_result": rows}),
                 "imperfection")
    b = _section(_export(fx["res2"], {"imperfection_result": (rows, results)}),
                 "imperfection")
    assert a == b
    assert "could not be built" not in a
    assert a.count("<tr") == 1 + len(rows)
    sec = _plain(a)
    assert "perfect" in sec
    for r in rows:
        assert f"{r['lambda_limit']:.4g}" in sec
        if r["amplitude_m"] > 0:
            assert f"{r['span_over_amplitude']:.0f}" in sec
            assert f"{r['ratio']:.3f}" in sec
    # the knock-down falls with amplitude (2V dome, both signs run)
    ratios = [r["ratio"] for r in rows]
    assert all(x > y for x, y in zip(ratios, ratios[1:])), ratios
    fig = rp.imperfection_figure(rows)
    main = next(ln for ln in fig.axes[0].lines
                if ln.get_label().startswith("First load maximum"))
    assert np.allclose(sorted(main.get_ydata()), sorted(ratios))
    assert np.allclose(sorted(main.get_xdata()),
                       sorted(r["amplitude_over_span"] for r in rows))


def test_modal_section_table_and_mode_shapes():
    fx = _fx()
    m = fx["modal_va"]
    text = _export(fx["res_va"], {"modal_result": m}, include_figures=True)
    sec = _section(text, "modal")
    assert "could not be built" not in sec
    plain = _plain(sec)
    assert sec.count("<tr") == 1 + m.n_modes
    for f in m.frequencies:
        assert f"{f:.4f}" in plain
    assert f"{m.fundamental_frequency:.4f} Hz" in plain
    tm = float(np.asarray(m.total_mass)[0])
    assert f"{tm:,.1f} kg" in plain
    assert sec.count("<img ") == 1
    fig = rp.mode_shapes_figure(fx["va"], m, 3)
    assert len(fig.axes) == min(3, m.n_modes)
    for k, ax in enumerate(fig.axes):
        assert ax.get_title().startswith(f"Mode {k + 1}")
    # a modal result from another model: table yes, shapes no, no crash
    other = _section(_export(fx["res2"], {"modal_result": m},
                             include_figures=True), "modal")
    assert "belongs to a different model" in other
    assert "<img " not in other
    # 2V dome: twelve modes, the first three drawn
    big = _section(_export(fx["res2"], {"modal_result": fx["modal2"]},
                           include_figures=True), "modal")
    assert big.count("<tr") == 1 + fx["modal2"].n_modes
    assert "The first 3 mode shape(s)" in big


def test_fabrication_section_2v_closed_forms():
    fx = _fx()
    t = fx["fab2"]
    sec = _section(_export(fx["res2"], {"fabrication_tables": t}),
                   "fabrication")
    assert "could not be built" not in sec
    titles = re.findall(r"<h3>([^<]*)</h3>", sec)
    assert titles == ["Strut cut list", "Hub types", "Panel types",
                      "Dihedral angles between panels", "Cutting patterns"]
    tables = re.findall(r"<table>.*?</table>", sec, re.S)
    assert [tb.count("<tr") - 1 for tb in tables] == [
        len(t["struts"]), len(t["hubs"]), len(t["panels"]),
        len(t["dihedrals"]), len(t["stock_patterns"])]
    plain = _plain(sec)
    cf_a = 2.0 * math.sin(math.acos(1.0 / math.sqrt(5.0)) / 4.0)
    cf_b = 2.0 * math.sin(math.radians(18.0))
    for cf in (cf_a, cf_b):
        assert f"{cf:.6f}" in plain                     # 0.546533, 0.618034
        assert f"{math.degrees(math.acos(cf / 2.0)):.3f}" in plain
    assert "74.141" in plain and "72.000" in plain
    # cut length = L - 2 x 50 mm hub offset
    for r in t["struts"]:
        assert f"{r['length_mm'] - 100.0:.1f}" in plain
    for r in t["hubs"]:
        assert r["arrangement"] in plain
    for r in t["panels"]:
        assert r["panel_type"] in plain
    # von Mises truss: no panels, not spherical, still a clean section
    vm = _section(_export(fx["res_vm"], {"fabrication_tables": fx["fab_vm"]}),
                  "fabrication")
    assert "could not be built" not in vm
    assert "This model has no panels." in vm
    assert "Dihedral angles" not in vm and "Cutting patterns" not in vm
    for n in fx["fab_vm"]["notes"]:
        assert _plain(n["note"])[:40] in _plain(vm)


def test_study_section_table_and_charts():
    fx = _fx()
    rows = fx["study"]
    text = _export(fx["res2"], {"study_rows": rows}, include_figures=True)
    sec = _section(text, "study")
    assert "could not be built" not in sec
    table = re.findall(r"<table>.*?</table>", sec, re.S)
    assert len(table) == 1 and table[0].count("<tr") == 1 + len(rows)
    plain = _plain(sec)
    for r in rows:
        assert st._fmt("mass_kg", r["mass_kg"]) in plain
        assert st._fmt("lambda_cr", r["lambda_cr"]) in plain
    best = st.best_design(rows)
    assert f"{int(best['frequency'])}V, R = " in plain
    assert "Frequencies 2, 3" in plain
    for msg in rows.messages:
        assert _plain(msg)[:40] in plain
    assert sec.count("<img ") == 1
    fig = rp.study_figure(rows)
    assert len(fig.axes) == 4
    # a plain list (no .settings / .messages) also works
    lean = _section(_export(fx["res2"], {"study_rows": list(rows)}), "study")
    assert "could not be built" not in lean
    assert lean.count("<tr") == 1 + len(rows)


def test_bad_extra_gives_a_note_not_a_crash():
    fx = _fx()
    text = _export(fx["res2"], {"nonlinear_result": object(),
                                "fabrication_tables": {"struts": "x"},
                                "modal_result": fx["modal2"]},
                   include_figures=True)
    _assert_well_formed(text)
    assert "could not be built" in _section(text, "nonlinear")
    assert "could not be built" in _section(text, "fabrication")
    assert "could not be built" not in _section(text, "modal")
    assert text.endswith("</footer></div></body></html>")


def test_report_module_stays_headless_safe():
    """report.py may import matplotlib's Figure (it always did) but never
    pyplot or Qt."""
    with open(os.path.join(ROOT, "geodome", "report.py"), encoding="utf-8") as fh:
        src = fh.read()
    assert "pyplot" not in src
    assert "qtpy" not in src and "PyQt" not in src
    assert not re.search(r"[A-Za-z]:\\Users", src)


# ================================================================= theory ==
def _theory_module():
    """Load gui/theory.py by path (it imports nothing), so the test does not
    need Qt."""
    path = os.path.join(ROOT, "geodome", "gui", "theory.py")
    spec = importlib.util.spec_from_file_location("_theory_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, path


def test_custom_truss_report_and_csv_have_no_nan():
    """Review GUI-11: a truss without a sphere and without panels shows no
    raw 'nan' (plan mass, chord factors) and its summary is headed SPACE
    TRUSS."""
    m = bm.make_preset("three_bar")
    res = run_analysis(m, do_buckling=False)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "r.html")
        rp.export_html(res, p, include_figures=False)
        with open(p, encoding="utf-8") as fh:
            page = fh.read()
        text = _plain(page)
        assert not re.search(r"\bnan\b", text, re.I), \
            re.findall(r".{30}\bnan\b.{30}", text, re.I)[:3]
        assert "no panels, so no plan area" in text
        files = rp.export_csv(res, d)
        sched = [f for f in files if f.endswith("_strut_schedule.csv")][0]
        with open(sched, encoding="utf-8") as fh:
            assert "nan" not in fh.read().lower()
        summ = [f for f in files if f.endswith("_summary.txt")][0]
        with open(summ, encoding="utf-8") as fh:
            body = fh.read()
        assert "SPACE TRUSS -- ANALYSIS SUMMARY" in body
        assert not re.search(r"\bnan\b", body, re.I)


def test_theory_page_is_well_formed_and_complete():
    mod, path = _theory_module()
    page = mod.theory_html()
    assert isinstance(mod.THEORY_HTML, str) and page == mod.THEORY_HTML
    c = _parse(page)
    assert not c.errors, c.errors[:5]
    assert not c.stack, c.stack
    assert len(c.ids) == len(set(c.ids))
    text = _plain(page)
    for needle in ("direct stiffness", "ASCE 7-16", "self-straining",
                   "E3", "E7", "D2", "K + λ Kg", "corotational",
                   "Newton-Raphson", "displacement control",
                   "arc-length", "effective modal mass", "chord factor",
                   "axial angle", "dihedral", "OpenSees",
                   "hand calculation", "References"):
        assert needle.lower() in text.lower(), needle
    for ref in ("McGuire", "Gallagher", "Ziemian", "2000", "Crisfield",
                "1991", "1981", "Chopra", "2023", "Kenner", "1976",
                "AISC 360-16", "ASCE/SEI 7-16", "McKenna", "2010"):
        assert ref in text, ref
    # every anchor in the contents list points at an id on the page
    for target in re.findall(r"href='#([^']+)'", page):
        assert target in c.ids, target
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    import ast
    assert not [n for n in ast.walk(ast.parse(src))
                if isinstance(n, (ast.Import, ast.ImportFrom))], \
        "theory.py must not import anything"
    # the contents list names every section heading
    assert len(re.findall(r"href='#", page)) == \
        len([h for h in c.headings if h[0] == "h2"])
    assert not re.search(r"[A-Za-z]:\\Users", src)
    assert "@" not in src                         # no e-mail addresses


# ============================================================== driver ====
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
    sys.stdout.reconfigure(encoding="utf-8")
    raise SystemExit(_run_standalone())
