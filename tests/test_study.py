"""
Tests for the parameter study (geodome/study.py).

The study only orchestrates the existing analysis, so the tests check it in
three independent ways:

1. closed forms: Class I hemisphere counts (nodes, struts, panels, Euler),
   the 2V chord factors / strut counts and the resulting steel mass and plan
   area, computed from trigonometry rather than from the program;
2. exact dimensional-analysis scaling laws (u ~ R^3, lambda_cr ~ R^-2, ...)
   between two radii, which no single-model test can see;
3. plumbing identities: every row equals a direct ``run_analysis`` /
   ``auto_size`` of the same dome, the solver identities hold in every row,
   the base model is not mutated, and CSV / plots / tables read back what the
   rows hold.

Run with pytest, or standalone::

    python tests/test_study.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import geometry as geo                                 # noqa: E402
from geodome import sections as sx                                  # noqa: E402
from geodome import study as st                                     # noqa: E402
from geodome.model import DomeModel, LoadCombination, \
    default_model                                                   # noqa: E402
from geodome.optimize import auto_size                              # noqa: E402
from geodome.postprocess import run_analysis                        # noqa: E402


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _pressure_only_model(nu=4, R=10.0, section="CHS 88.9 x 4.0"):
    """Pressure loads only (no self weight, no point load), one ULS combo."""
    m = DomeModel(geom=geo.generate_dome(nu, R, 0.5))
    m.group_sections = [section] * m.geom.n_strut_types()
    m.loads.include_self_weight = False
    m.loads.crown_point_load = 0.0
    m.combinations = [
        LoadCombination("ULS 1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("SLS 1.0D + 1.0S", D=1.0, S=1.0, serviceability=True),
    ]
    return m


# ============================================================ closed forms ==
def test_hemisphere_counts_closed_form():
    """V = (10v^2+5v+2)/2, E = (30v^2+5v)/2, F = 10v^2, V - E + F = 1."""
    base = default_model(4, 15.0, 0.5)
    rows = st.parameter_study(base, [2, 4, 6, 8], radii=[15.0],
                              fractions=[0.5], do_buckling=False)
    assert len(rows) == 4
    for r in rows:
        nu = r["frequency"]
        assert r["nodes"] == (10 * nu ** 2 + 5 * nu + 2) // 2
        assert r["struts"] == (30 * nu ** 2 + 5 * nu) // 2
        assert r["panels"] == 10 * nu ** 2
        assert r["nodes"] - r["struts"] + r["panels"] == 1      # disc
        assert r["supports"] == 5 * nu                          # equator ring
        assert abs(r["fraction"] - 0.5) < 1e-12
        assert abs(r["span_m"] - 30.0) < 1e-9                   # 2R
        assert abs(r["rise_m"] - 15.0) < 1e-9                   # R
    # published strut-length counts: 2V = 2, 4V = 6 (test_geodome checks the
    # 4V chord factors against the published table)
    assert rows[0]["families"] == 2 and rows[1]["families"] == 6


def test_2v_mass_and_plan_area_closed_form():
    """2V hemisphere: 30 struts of 2 sin(atan(2)/4) R, 35 of 2 sin(18 deg) R.

    Mass = rho A R (30 cf_A + 35 cf_B); plan area = area of the regular
    10-gon inscribed in the base circle = (10/2) R^2 sin(36 deg).
    """
    R = 15.0
    sec = "CHS 114.3 x 6.3"
    base = default_model(4, R, 0.5)
    rows = st.parameter_study(base, [2], radii=[R], fractions=[0.5],
                              section=sec, do_buckling=False)
    r = rows[0]
    cf_a = 2.0 * math.sin(math.atan(2.0) / 4.0)          # 0.546533
    cf_b = 2.0 * math.sin(math.radians(18.0))            # 0.618034
    assert abs(cf_a - 0.546533) < 1e-6 and abs(cf_b - 0.618034) < 1e-6
    A = sx.SECTIONS[sec].A
    rho = base.material.rho
    mass = rho * A * R * (30 * cf_a + 35 * cf_b)
    plan = 0.5 * 10 * R ** 2 * math.sin(math.radians(36.0))
    assert _rel(r["mass_kg"], mass) < 1e-12, (r["mass_kg"], mass)
    assert _rel(r["mass_per_plan_kg_m2"], mass / plan) < 1e-12
    assert _rel(r["total_strut_length_m"], R * (30 * cf_a + 35 * cf_b)) < 1e-12
    assert _rel(r["strut_min_m"], cf_a * R) < 1e-12
    assert _rel(r["strut_max_m"], cf_b * R) < 1e-12
    assert r["sections"] == f"all: {sec}"
    # the per-family counts, through the rebuilt model of the row
    sched = st.model_from_row(base, r).geom.strut_schedule()
    assert [s["count"] for s in sched] == [30, 35]


def test_dimensional_scaling_laws():
    """Fixed section, pressure loads only: R -> 2R gives exact ratios.

    K ~ EA/L ~ 1/R, F ~ p R^2  =>  u ~ R^3, N ~ R^2, Kg ~ N/L ~ R,
    lambda_cr ~ K/Kg ~ R^-2, span/deflection ~ R^-2, mass ~ R,
    mass / plan area ~ 1/R.
    """
    base = _pressure_only_model()
    rows = st.parameter_study(base, [4], radii=[10.0, 20.0], fractions=[0.5])
    a, b = rows
    assert a["radius_m"] == 10.0 and b["radius_m"] == 20.0
    for k in ("nodes", "struts", "families", "governing_combination",
              "buckling_case"):
        assert a[k] == b[k], k
    checks = {
        "mass_kg": 2.0,
        "mass_per_plan_kg_m2": 0.5,
        "total_strut_length_m": 2.0,
        "strut_max_m": 2.0,
        "max_service_disp_mm": 8.0,
        "crown_disp_mm": 8.0,
        "span_over_deflection": 0.25,
        "lambda_cr": 0.25,
    }
    for k, ratio in checks.items():
        got = b[k] / a[k]
        assert _rel(got, ratio) < 1e-9, (k, got, ratio)


# ======================================================= plumbing checks ====
def test_row_equals_direct_analysis():
    base = default_model(4, 12.0, 0.5)
    before = list(base.group_sections)
    rows = st.parameter_study(base, [4], radii=[12.0], fractions=[0.5])
    r = rows[0]
    res = run_analysis(base, n_buckling_modes=1)
    assert r["governing_combination"] == res.governing_case
    assert r["governing_member"] == res.governing_member
    assert _rel(r["max_util"], res.governing_util) < 1e-12
    assert _rel(r["lambda_cr"], res.buckling.load_factors[0]) < 1e-12
    assert _rel(r["mass_kg"], base.steel_mass()) < 1e-12
    assert _rel(r["span_over_deflection"], res.deflection_ratio) < 1e-12
    assert _rel(r["crown_disp_mm"], res.crown_disp * 1e3) < 1e-12
    # the study reports the lowest load factor over the strength
    # combinations, which need not be the governing combination
    assert r["buckling_case"] == res.lambda_cr_case
    assert _rel(r["lambda_cr"], min(v for v in res.lambda_cr_by_case.values()
                                    if v == v)) < 1e-12
    assert base.group_sections == before
    # the solver identities hold
    assert r["equilibrium_residual"] < 1e-9
    assert r["energy_error"] < 1e-9


def test_study_model_copies_settings():
    base = default_model(4, 12.0, 0.5)
    base.material_name = "EN S355"
    base.support_type = geo.ELASTIC
    base.spring_stiffness = (1e7, 2e7, 3e7)
    base.K_effective = 0.9
    base.deflection_limit = 300.0
    base.phi_c, base.phi_t = 0.85, 0.88
    base.loads.snow_pressure = 900.0
    base.loads.temperature_change = 25.0
    base.combinations.append(LoadCombination("T  1.0D + 1.0T", D=1.0, T=1.0))
    m = st.study_model(base, 6, 9.0, 0.5)
    assert m.geom.frequency == 6 and m.geom.radius == 9.0
    assert m.material_name == "EN S355"
    assert m.support_type == geo.ELASTIC
    assert m.spring_stiffness == (1e7, 2e7, 3e7)
    assert (m.K_effective, m.deflection_limit, m.phi_c, m.phi_t) == \
        (0.9, 300.0, 0.85, 0.88)
    assert m.loads.snow_pressure == 900.0
    assert m.loads.temperature_change == 25.0
    assert [c.name for c in m.combinations] == [c.name for c in
                                                base.combinations]
    assert m.combinations[-1].T == 1.0
    assert len(set(m.group_sections)) == 1
    assert len(m.group_sections) == m.geom.n_strut_types()
    # deep copies: editing the study model leaves the base alone
    m.loads.snow_pressure = 1.0
    m.combinations[0].D = 9.9
    assert base.loads.snow_pressure == 900.0
    assert base.combinations[0].D != 9.9


def test_uniform_section_rule():
    base = default_model(4, 12.0, 0.5)
    assert st.uniform_section(base)[0] == base.group_sections[0]
    base.group_sections[0] = "CHS 60.3 x 3.2"
    base.group_sections[1] = "CHS 139.7 x 6.3"
    sec, why = st.uniform_section(base)
    assert sec == "CHS 139.7 x 6.3" and "heaviest" in why
    assert st.uniform_section(base, "CHS 48.3 x 3.2")[0] == "CHS 48.3 x 3.2"
    try:
        st.uniform_section(base, "no such tube")
    except ValueError:
        pass
    else:
        raise AssertionError("unknown section accepted")


def _uniform_scan(model, target, check_deflection):
    """Brute force: the lightest catalogue section that, put on every
    family, passes the final check (utilisation <= target, KL/r <= 200 and,
    optionally, the deflection limit).  Uses run_analysis only for the
    member forces; the selection logic of optimize is not used."""
    for name in sx.SECTION_NAMES:
        m = DomeModel.from_dict(model.to_dict())      # a fresh copy
        m.group_sections = [name] * m.geom.n_strut_types()
        res = run_analysis(m, do_buckling=False)
        if res.governing_util > target + 1e-12:
            continue
        if np.max(res.caps["slenderness"]) > 200.0:
            continue
        if check_deflection and not res.deflection_ok:
            continue
        return name
    return None


def test_auto_size_uses_the_design_wall_thickness():
    """Review F3 / V03: the selection must use the B4.2 design wall
    (factor * t) like the final check, so the result meets the target."""
    from geodome import benchmarks as bm
    from geodome import buckling as bk
    from geodome.optimize import _select_section
    # unit level: a 3 m strut with 100 kN compression, factor 0.93
    m = default_model(4, 12.0, 0.5)
    m.design_thickness_factor = 0.93
    name = _select_section(m, 3.0, 0.0, -100e3, 1.0)
    i = sx.section_index(name)
    cap = lambda k: bk.compression_capacity(                      # noqa: E731
        bk.design_section(sx.CATALOGUE[k], 0.93), m.material, 3.0)[0]
    assert cap(i) >= 100e3
    assert i == 0 or cap(i - 1) < 100e3
    # model level: the case study with factor 0.93 and target 1.00
    cs = bm.make_preset("case_study_6v")
    cs.design_thickness_factor = 0.93
    s = auto_size(cs, target_util=1.0)
    fin = run_analysis(cs, do_buckling=False)
    assert s.converged
    assert fin.governing_util <= 1.0 + 1e-12, fin.governing_util
    assert not any("heaviest" in t for t in s.messages), s.messages
    # one-family determinate Model A: exactly the brute-force answer
    a = bm.make_preset("validation_a")
    assert a.geom.n_strut_types() == 1
    a.design_thickness_factor = 0.93
    want = _uniform_scan(a, 1.0, False)
    s = auto_size(a, target_util=1.0, enforce_deflection=False)
    assert a.group_sections == [want], (a.group_sections, want)


def test_auto_size_deflection_step_up_is_kept():
    """Review F4: a deflection-governed design must not oscillate.  With
    span/8000 the 4V dome needs CHS 114.3 x 6.3 on every family (the
    lightest uniform section that passes both checks)."""
    m = default_model(4, 12.0, 0.5)
    m.deflection_limit = 8000.0
    base = default_model(4, 12.0, 0.5)
    base.deflection_limit = 8000.0
    want = _uniform_scan(base, 0.95, True)
    s = auto_size(m, max_iter=12)
    fin = run_analysis(m, do_buckling=False)
    assert s.converged, s.messages
    assert fin.deflection_ok
    assert set(m.group_sections) == {want}, (m.group_sections, want)
    assert want == "CHS 114.3 x 6.3"
    assert any("deflection limit governs" in t for t in s.messages)
    # after the first deflection step no family ever goes down again
    idx = [[sx.section_index(n) for n in h.sections] for h in s.history]
    k0 = next(k for k, h in enumerate(s.history)
              if "deflection" in h.note)
    for a, b in zip(idx[k0 + 1:], idx[k0 + 2:]):
        assert all(y >= x for x, y in zip(a, b)), (a, b)


def test_auto_size_needs_a_strength_combination():
    """Review GUI-08: with every combination a service case there is
    nothing to size for; auto_size must refuse, not pick the heaviest
    tube."""
    m = default_model(3, 8.0, 0.5)
    for c in m.combinations:
        c.serviceability = True
    try:
        auto_size(m)
    except ValueError as exc:
        assert "Service" in str(exc)
    else:
        raise AssertionError("auto_size ran without a strength combination")
    res = run_analysis(m, do_buckling=False)
    assert res.governing_case == ""
    assert any("No strength combination" in t for t in res.messages)


def test_auto_size_rows_equal_direct_sizing():
    base = default_model(4, 10.0, 0.5)
    before = list(base.group_sections)
    rows = st.parameter_study(base, [3, 4], radii=[10.0], fractions=[0.5],
                              auto_size=True)
    assert base.group_sections == before                # base untouched
    for r in rows:
        m = st.study_model(base, r["frequency"], 10.0, 0.5)
        s = auto_size(m)
        res = run_analysis(m, n_buckling_modes=1)
        assert r["auto_sized"] is True
        assert r["group_sections"] == m.group_sections
        assert r["sizing_converged"] == s.converged
        assert r["sizing_iterations"] == s.iterations
        assert _rel(r["mass_kg"], m.steel_mass()) < 1e-12
        assert _rel(r["max_util"], res.governing_util) < 1e-12
        assert _rel(r["lambda_cr"], res.buckling.load_factors[0]) < 1e-12
        if r["sizing_converged"]:
            # fully-stressed design: every family at or below the target
            assert r["max_util"] <= 0.95 + 1e-12
        # the rebuilt row model reproduces the row
        m2 = st.model_from_row(base, r)
        assert m2.group_sections == r["group_sections"]
        assert _rel(m2.steel_mass(), r["mass_kg"]) < 1e-12


def test_snapped_fractions_reuse_the_geometry():
    base = default_model(4, 10.0, 0.5)
    rows = st.parameter_study(base, [4], radii=[10.0],
                              fractions=[0.5, 0.52])
    a, b = rows
    assert a["fraction_requested"] == 0.5 and b["fraction_requested"] == 0.52
    assert a["fraction"] == b["fraction"] == 0.5
    assert "reused" in b["notes"]
    for k in ("nodes", "struts", "mass_kg", "max_util", "lambda_cr",
              "span_over_deflection"):
        assert a[k] == b[k], k


def test_node_specific_inputs_are_dropped_with_a_message():
    base = default_model(4, 10.0, 0.5)
    base.loads.settlements = [{"node": int(base.geom.base_nodes[0]),
                               "uz": -0.01}]
    base.loads.nodal_loads = [{"node": 0, "Fz": -5e3, "case": "L"}]
    base.loads.temperature_change = 30.0
    base.combinations = [
        LoadCombination("ULS 1.2D + 1.6S + 1.0T", D=1.2, S=1.6, T=1.0),
        LoadCombination("SLS 1.0D + 1.0T", D=1.0, T=1.0, serviceability=True),
    ]
    rows = st.parameter_study(base, [4], radii=[10.0], fractions=[0.5])
    text = " ".join(rows.messages)
    assert "settlements" in text and "nodal loads" in text
    m = st.study_model(base, 4, 10.0, 0.5)
    assert m.loads.settlements == [] and m.loads.nodal_loads == []
    assert m.loads.temperature_change == 30.0           # uniform: carried
    assert len(base.loads.settlements) == 1             # base untouched
    res = run_analysis(m, n_buckling_modes=1)
    r = rows[0]
    assert _rel(r["crown_disp_mm"], res.crown_disp * 1e3) < 1e-12
    # the temperature really reaches the study point
    base.loads.temperature_change = 0.0
    r0 = st.parameter_study(base, [4], radii=[10.0], fractions=[0.5],
                            do_buckling=False)[0]
    assert abs(r0["crown_disp_mm"] - r["crown_disp_mm"]) > 1.0


def test_custom_supports_and_bad_inputs_raise():
    base = default_model(4, 10.0, 0.5)
    bad = [
        dict(frequencies=[]),
        dict(frequencies=[2.5]),
        dict(frequencies=[0]),
        dict(frequencies=[4], radii=[-1.0]),
        dict(frequencies=[4], radii=[float("nan")]),
        dict(frequencies=[4], fractions=[0.0]),
        dict(frequencies=[4], fractions=[1.2]),
        dict(frequencies=[4], section="CHS 1 x 1"),
        dict(frequencies=[4], auto_size=True, target_util=1.5),
    ]
    for kw in bad:
        try:
            st.parameter_study(base, **kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {kw}")
    base.support_type = geo.CUSTOM
    base.custom_supports = {0: [True, True, True]}
    try:
        st.parameter_study(base, [4])
    except ValueError as exc:
        assert "CUSTOM" in str(exc)
    else:
        raise AssertionError("CUSTOM supports accepted")


def test_scalar_inputs_and_defaults():
    base = default_model(3, 8.0, 0.5)
    rows = st.parameter_study(base, 3, do_buckling=False)
    assert len(rows) == 1
    r = rows[0]
    assert r["radius_m"] == 8.0
    assert r["fraction_requested"] == base.geom.dome_fraction
    assert r["nodes"] == base.geom.n_nodes
    assert math.isnan(r["lambda_cr"]) and r["buckling_case"] == ""


def test_failed_point_gives_error_row_and_study_continues():
    base = default_model(4, 10.0, 0.5)
    base.combinations = []                    # run_analysis refuses this
    rows = st.parameter_study(base, [2, 3], radii=[10.0], fractions=[0.5])
    assert len(rows) == 2
    for r in rows:
        assert r["status"] == "ERROR"
        assert "No load combinations" in r["notes"]
        assert math.isnan(r["mass_kg"])
    assert sum("failed" in m for m in rows.messages) == 2
    # the outputs cope with error rows
    with tempfile.TemporaryDirectory() as d:
        st.study_to_csv(rows, os.path.join(d, "err.csv"))
    st.study_table(rows)
    assert st.best_design(rows) is None


def test_progress_and_cancel():
    base = default_model(4, 10.0, 0.5)
    seen = []
    st.parameter_study(base, [2, 3, 4], radii=[10.0], fractions=[0.5],
                       auto_size=True, progress=lambda f, m: seen.append(f))
    assert seen and all(0.0 <= f <= 1.0 for f in seen)
    assert all(b >= a - 1e-12 for a, b in zip(seen, seen[1:]))   # monotonic
    assert seen[-1] == 1.0

    calls = {"n": 0}

    def cancel():
        calls["n"] += 1
        return calls["n"] > 1                 # allow one point, then stop
    rows = st.parameter_study(base, [2, 3, 4], radii=[10.0], fractions=[0.5],
                              do_buckling=False, cancel=cancel)
    assert len(rows) == 1 and rows[0]["frequency"] == 2
    assert any("Stopped" in m for m in rows.messages)


def test_status_and_best_design():
    base = default_model(4, 15.0, 0.5)
    rows = st.parameter_study(base, [2, 4, 6, 8], radii=[15.0],
                              fractions=[0.5])
    for r in rows:
        ok = (r["max_util"] <= 1.0
              and r["span_over_deflection"] >= r["deflection_limit"]
              and not r["lambda_cr"] < 1.0)
        assert (r["status"] == "OK") == ok, r["status"]
    best = st.best_design(rows)
    good = [r for r in rows if r["status"] == "OK"]
    assert best is not None and best["mass_kg"] == min(r["mass_kg"]
                                                       for r in good)
    # an impossible stability demand leaves nothing
    assert st.best_design(rows, min_lambda_cr=1e9) is None


# =========================================================== outputs ========
def test_csv_round_trip():
    base = default_model(4, 10.0, 0.5)
    rows = st.parameter_study(base, [2, 3], radii=[10.0, 12.0],
                              fractions=[0.5], auto_size=True)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "study.csv")
        st.study_to_csv(rows, path)
        with open(path, newline="", encoding="utf-8") as fh:
            rd = csv.DictReader(fh)
            assert rd.fieldnames == st.COLUMN_KEYS
            back = list(rd)
    assert len(back) == len(rows)
    for r, b in zip(rows, back):
        for k in st.COLUMN_KEYS:
            v = r[k]
            if isinstance(v, bool):
                assert b[k] == ("yes" if v else "no"), k
            elif isinstance(v, float):
                if math.isnan(v):
                    assert b[k] == "nan", k
                else:
                    assert float(b[k]) == v, k          # repr: exact
            elif isinstance(v, int):
                assert int(b[k]) == v, k
            elif v is None:
                assert b[k] == "", k
            else:
                assert b[k] == str(v), k


def test_plot_study_draws_every_series():
    base = default_model(4, 10.0, 0.5)
    rows = st.parameter_study(base, [2, 3, 4], radii=[10.0, 15.0],
                              fractions=[0.5])
    ys = ["mass_kg", "max_util", "lambda_cr", "span_over_deflection"]
    fig, axes = st.plot_study(None, rows, x="frequency", ys=ys)
    assert len(axes) == len(ys)
    for ax, yk in zip(axes, ys):
        lines = [ln for ln in ax.get_lines() if ln.get_linestyle() == "-"]
        assert len(lines) == 2                         # one per radius
        for ln, R in zip(lines, (10.0, 15.0)):
            want = [r[yk] for r in rows if r["radius_m"] == R]
            assert np.allclose(ln.get_ydata(), want)
            assert np.allclose(ln.get_xdata(), [2, 3, 4])
    assert axes[0].get_legend() is not None
    # radius on the x axis: one line per frequency
    fig2, ax2 = st.plot_study(None, rows, x="radius_m", ys="mass_kg")
    solid = [ln for ln in ax2[0].get_lines() if ln.get_linestyle() == "-"]
    assert len(solid) == 3
    # an existing figure is reused
    fig3, _ = st.plot_study(fig, rows, ys=["mass_kg"])
    assert fig3 is fig and len(fig.axes) == 1


def test_study_table_text():
    base = default_model(4, 10.0, 0.5)
    rows = st.parameter_study(base, [2, 4], radii=[10.0], fractions=[0.5],
                              do_buckling=False)
    lines = st.study_table(rows)
    assert len(lines) == 3 + len(rows)
    assert all(ln.isascii() for ln in lines)
    assert "lambda_cr" in lines[0] and "n/a" in lines[3]   # no buckling run


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
