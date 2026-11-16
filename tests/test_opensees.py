"""
Tests for the OpenSees bridge (spec task E, ``geodome/opensees_bridge.py``
and ``geodome/_opensees_runner.py``).

Without OpenSees (always run):

* the spec is plain JSON and carries exactly the model's inputs (load
  vector of the combination, alpha*dT*T, T*settlement, supports, springs,
  masses);
* the runner imports nothing but the standard library and openseespy;
* the OpenSees command list for every support ingredient (fix, sp on a
  settling restraint, zeroLength springs to ground nodes, a moving spring
  ground, InitStrainMaterial(-eps0), mass);
* piecewise displacement control segments of a path that snaps back;
* MAC / subspace MAC identities; interpreter discovery; a clear error when
  no OpenSees interpreter exists; the exported script parses.

With OpenSees (skipped cleanly when it is not available):

* corotTruss strain = engineering strain (vs Green and log, decided
  numerically); InitStrainMaterial sign on a heated bar vs closed form;
* linear static: every support type, every combination, temperature,
  settlement, custom springs -- GeoDomeLab = OpenSees to 1e-10;
* a 1 ppm change of E in the OpenSees input is detected (sensitivity);
* modal: all modes, lumped and consistent mass (frequencies, MAC);
* nonlinear: von Mises truss path = OpenSees = closed form; star dome
  through both snap-backs (piecewise displacement control);
* exported standalone scripts run and reproduce GeoDomeLab;
* in-process and subprocess execution give the same numbers.

Run with pytest::

    pytest -q tests/test_opensees.py

or standalone, without pytest installed::

    python tests/test_opensees.py
"""

from __future__ import annotations

import ast
import json
import math
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import _opensees_runner as runner                     # noqa: E402
from geodome import benchmarks as bm                                # noqa: E402
from geodome import geometry as geo                                 # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import modal as md                                     # noqa: E402
from geodome import nonlinear as nl                                 # noqa: E402
from geodome import opensees_bridge as ob                           # noqa: E402
from geodome import solver as slv                                   # noqa: E402
from geodome.model import LoadCombination, default_model            # noqa: E402


def _need_opensees():
    if not ob.opensees_available():
        raise unittest.SkipTest("OpenSees (openseespy) is not available; "
                                f"set {ob.ENV_VAR} to an interpreter that "
                                f"has it.")


def _t_model():
    """3V dome on elastic supports with temperature and settlement."""
    m = default_model(3, 10.0, 0.5)
    m.support_type = geo.ELASTIC
    b = m.geom.base_nodes
    m.loads.temperature_change = 20.0
    m.loads.settlements = [{"node": int(b[0]), "uz": -0.004}]
    m.combinations.append(LoadCombination("T1  1.2D + 1.6S + 1.5T", D=1.2,
                                          S=1.6, T=1.5))
    return m


def _custom_model():
    """3V dome on mixed custom supports, cooling and settling supports."""
    m = default_model(3, 10.0, 0.5)
    cs = {}
    for k, nd in enumerate(m.geom.base_nodes):
        cs[int(nd)] = ([True, True, True] if k % 3 == 0 else
                       [2e7, 3e7, True] if k % 3 == 1 else
                       [False, 1e8, 5e7])
    m.support_type = geo.CUSTOM
    m.custom_supports = cs
    m.loads.temperature_change = -25.0
    nd1 = [nd for nd, v in cs.items() if v[1] == 1e8][0]
    m.loads.settlements = [{"node": int(m.geom.base_nodes[0]), "uz": -0.008},
                           {"node": nd1, "uy": 0.004, "uz": 0.002}]
    m.combinations.append(LoadCombination("T1  1.2D + 1.6S + 1.0T", D=1.2,
                                          S=1.6, T=1.0))
    return m


# ===================================================== no OpenSees needed ===
def test_spec_is_json_and_carries_the_inputs():
    m = _t_model()
    c = ob.find_combination(m, "T1")
    spec = ob.to_opensees_spec(m, "T1")
    back = json.loads(json.dumps(spec))
    assert back == spec
    g = m.geom
    assert len(spec["nodes"]) == g.n_nodes
    assert len(spec["members"]) == g.n_members
    F = ld.combine(ld.basic_load_vectors(m), c)
    assert np.array_equal(np.asarray(spec["load"]), F)
    eps0 = np.asarray(spec["eps0"])
    assert np.allclose(eps0, 1.5 * m.material.alpha * 20.0, rtol=0,
                       atol=1e-18)
    dof = 3 * int(g.base_nodes[0]) + 2
    assert spec["prescribed"] == [[dof, 1.5 * -0.004]]
    fixed, springs = m.support_conditions()
    assert spec["fixed"] == [int(d) for d in fixed]
    assert dict((d, k) for d, k in spec["springs"]) == \
        {int(d): float(k) for d, k in springs.items()}
    assert np.allclose(spec["A"], m.member_areas())
    assert np.allclose(spec["E"], m.material.E)
    mm = md.lumped_masses(m)
    assert np.allclose(spec["line_mass"], m.material.rho * m.member_areas())
    assert np.allclose(spec["nodal_mass"], mm.non_structural)


def test_spec_without_combination_has_no_load():
    m = _t_model()
    spec = ob.to_opensees_spec(m, None)
    assert spec["combination"] is None
    assert not any(spec["load"])
    assert spec["eps0"] is None and spec["prescribed"] == []


def test_find_combination_variants():
    m = _t_model()
    c = m.combinations[1]
    assert ob.find_combination(m, c.name) is c
    assert ob.find_combination(m, "LC2") is c
    assert ob.find_combination(m, "  " + c.name.strip() + " ") is c
    assert ob.find_combination(m, c) is c
    try:
        ob.find_combination(m, "LC99")
    except KeyError:
        pass
    else:
        raise AssertionError("unknown combination accepted")


def test_runner_is_self_contained():
    """The runner must run in a bare OpenSees interpreter (no geodome)."""
    with open(runner.__file__, "r", encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            mods.add((node.module or "").split(".")[0])
    assert mods <= {"__future__", "json", "math", "sys", "time",
                    "openseespy"}, mods


def test_command_list_for_every_support_ingredient():
    m = _custom_model()
    spec = ob.to_opensees_spec(m, "T1", include_mass=True)
    cmds = runner.model_commands(spec, "Truss", mass_type="lumped")
    loads = runner.load_commands(spec)
    names = [c[0] for c in cmds]
    n = m.geom.n_nodes
    fixes = {a[0]: a[1:] for nm, a in cmds if nm == "fix"}
    base0 = int(m.geom.base_nodes[0])
    # node base0 is fully restrained but settles in z: z is not fixed, an
    # sp carries the settlement
    assert fixes[base0 + 1] == (1, 1, 0)
    sps = [a for nm, a in loads if nm == "sp"]
    assert (base0 + 1, 3, -0.008) in sps
    # zeroLength springs: one ground node per sprung node, fixed except a
    # moving ground end
    zl = [a for nm, a in cmds if nm == "element" and a[0] == "zeroLength"]
    sprung = sorted({d // 3 for d, _k in spec["springs"]
                     if d not in set(spec["fixed"])})
    assert len(zl) == len(sprung)
    nd1 = [nd for nd, v in m.custom_supports.items() if v[1] == 1e8][0]
    gtag = [a[2] for a in zl if a[3] == nd1 + 1][0]
    assert gtag > n
    assert fixes[gtag] == (1, 0, 0)            # y and z ground ends move
    assert (gtag, 2, 0.004) in sps and (gtag, 3, 0.002) in sps
    # thermal strain: InitStrainMaterial wraps Elastic with -eps0
    ism = [a for nm, a in cmds if nm == "uniaxialMaterial"
           and a[0] == "InitStrainMaterial"]
    assert len(ism) == 1
    assert ism[0][3] == -m.material.alpha * -25.0
    # lumped strut mass via -rho and nodal masses
    el = [a for nm, a in cmds if nm == "element" and a[0] == "Truss"]
    assert len(el) == m.geom.n_members
    assert all("-rho" in a and "-cMass" not in a for a in el)
    assert names.count("mass") == int(np.count_nonzero(spec["nodal_mass"]))
    cmds_c = runner.model_commands(spec, "Truss", mass_type="consistent")
    assert all("-cMass" in a for nm, a in cmds_c
               if nm == "element" and a[0] == "Truss")
    # a prescribed value on a free DOF is refused
    bad = dict(spec, prescribed=[[0, 0.01]])
    try:
        runner.model_commands(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("settlement of a free DOF accepted")


def test_path_segments_follow_a_snap_back():
    # dof 0 rises then falls (snap-back), dof 1 rises throughout
    t = np.linspace(0.0, 1.0, 11)
    U = np.zeros((11, 3))
    U[:, 0] = np.where(t <= 0.5, t, 1.0 - t)
    U[:, 1] = 2.0 * t
    segs = ob.path_segments(U, free=[0, 1], preferred=0)
    idx = [i for s in segs for i in s["index"]]
    assert idx == list(range(1, 11))
    for s in segs:
        assert s["targets"] == [float(U[i, s["control_dof"]])
                                for i in s["index"]]
    # dof 1 is monotone over the whole path, so it wins over the preferred
    assert len(segs) == 1 and segs[0]["control_dof"] == 1
    # with only dof 0 allowed the path splits at the reversal
    segs = ob.path_segments(U, free=[0], preferred=0)
    assert [len(s["index"]) for s in segs] == [5, 5]
    # repeated states are skipped
    U2 = np.vstack([U[:3], U[2:3], U[3:]])
    segs = ob.path_segments(U2, free=[0, 1], preferred=0)
    assert sum(len(s["index"]) for s in segs) == 10


def test_monotone_targets_stop_at_reversal():
    c = [0.0, -1.0, -2.0, -2.0, -3.0, -2.5, -4.0]
    assert ob.monotone_targets(c) == [1, 2, 4]


def test_mac_identities():
    rng = np.random.default_rng(0)
    a = rng.normal(size=12)
    assert abs(ob.mac(a, -3.0 * a) - 1.0) < 1e-15
    b = rng.normal(size=12)
    b -= (a @ b) / (a @ a) * a
    assert abs(ob.mac(a, b)) < 1e-28
    B = np.column_stack([a, b])
    th = 0.7
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th),
                                                     math.cos(th)]])
    assert abs(ob.subspace_mac(a, B @ R) - 1.0) < 1e-14
    assert ob.frequency_clusters([1.0, 1.0 + 1e-12, 2.0, 3.0, 3.0]) == \
        [[0, 1], [2], [3, 4]]


def test_interpreter_discovery_and_missing_interpreter():
    old = os.environ.get(ob.ENV_VAR)
    try:
        os.environ[ob.ENV_VAR] = os.path.join("no", "such", "python.exe")
        cands = ob.opensees_python_candidates()
        assert cands[0] == os.environ[ob.ENV_VAR]
        assert any(c.endswith(os.path.join("envs", "opensees",
                                           ob._python_name()))
                   for c in cands)
    finally:
        if old is None:
            os.environ.pop(ob.ENV_VAR, None)
        else:
            os.environ[ob.ENV_VAR] = old
    m = default_model(1, 5.0, 0.5)
    try:
        ob.run_opensees(ob.to_opensees_spec(m, "LC1"), "static",
                        python=os.path.join("no", "such", "python.exe"),
                        in_process=False)
    except ob.OpenSeesUnavailable:
        pass
    else:
        raise AssertionError("missing interpreter not reported")


def test_exported_script_parses():
    m = _t_model()
    with tempfile.TemporaryDirectory() as tmp:
        for an in ("static", {"type": "modal", "n_modes": 5}, "nonlinear"):
            p = ob.write_openseespy_script(
                m, "T1" if an == "static" else "LC2",
                os.path.join(tmp, "s.py"), an)
            with open(p, "r", encoding="utf-8") as fh:
                src = fh.read()
            ast.parse(src)
            assert "import openseespy.opensees as ops" in src
            if an == "static":
                assert "InitStrainMaterial" in src and "zeroLength" in src
                assert "ops.sp(" in src
            if an == "nonlinear":
                assert "corotTruss" in src and "DisplacementControl" in src


# ==================================================== OpenSees required ===
def test_corottruss_uses_engineering_strain():
    _need_opensees()
    s = ob.strain_measure_check()
    assert s["ok"], s
    assert s["measure"] == "engineering"
    assert s["rel_engineering"] < 1e-12
    assert s["rel_green"] > 0.05 and s["rel_log"] > 0.05
    assert s["equilibrium"] < 1e-12


def test_initstrainmaterial_sign_on_heated_bar():
    _need_opensees()
    t = ob.thermal_sign_check()
    assert t["ok"], t
    assert t["rel_u"] < 1e-12 and t["rel_N"] < 1e-12
    assert t["N1_os"] < 0.0                  # heating -> compression


def test_static_every_support_type_and_combination():
    _need_opensees()
    jobs, idx = [], []
    for st in (geo.PINNED, geo.ROLLER3, geo.ELASTIC):
        m = default_model(2, 12.0, 0.5)
        m.support_type = st
        for c in m.combinations:
            jobs.append(ob.static_job(m, c.name))
            idx.append((m, c.name))
    res = ob.run_opensees_batch(jobs)
    for (m, c), r in zip(idx, res):
        x = ob.compare_static(m, c, r)
        assert x["ok"], x
        for k in ("du", "dN", "dR", "dspring", "os_residual"):
            assert x[k] < 1e-10, (m.support_type, c, k, x[k])


def test_static_self_straining_and_custom_supports():
    _need_opensees()
    for m, combos in ((_custom_model(), ["LC2", "T1"]),
                      (_t_model(), ["T1"]),
                      (bm.make_preset("three_bar"),
                       [c.name for c in bm.make_preset("three_bar")
                        .combinations])):
        for c in combos:
            x = ob.crosscheck_static(m, c)
            assert x["ok"], x
            for k in ("du", "dN", "dR", "dspring"):
                assert x[k] < 1e-10, (c, k, x[k])
            assert x["gdl"].energy_error < 1e-10


def test_three_bar_opensees_matches_force_method():
    """OpenSees against GeoDomeLab's independent force-method solution."""
    _need_opensees()
    m = bm.make_preset("three_bar")
    c = ob.find_combination(m, "P + T")
    r = ob.run_opensees(*ob.static_job(m, c.name))
    nodal = m.loads.nodal_loads[0]
    ls = m.loads
    alpha = m.material.alpha
    fm = bm.three_bar_force_method(
        m, load=(nodal["Fx"], nodal["Fz"]),
        dT=(ls.temperature_change,) * 3,
        settlement={2: (0.0, ls.settlements[0]["uz"])})
    N = np.asarray(fm["N"])
    assert np.max(np.abs(r["N"] - N)) / np.max(np.abs(N)) < 1e-10
    assert alpha > 0.0


def test_sensitivity_one_ppm_is_detected():
    _need_opensees()
    m = default_model(2, 12.0, 0.5)
    spec, an = ob.static_job(m, "LC2")
    spec["E"] = [e * (1.0 + 1e-6) for e in spec["E"]]
    x = ob.compare_static(m, "LC2", ob.run_opensees(spec, an))
    assert 0.5e-6 < x["du"] < 2e-6, x["du"]


def test_modal_all_modes_lumped_and_consistent():
    _need_opensees()
    m = default_model(2, 12.0, 0.5)
    for mt in (md.LUMPED, md.CONSISTENT):
        x = ob.crosscheck_modal(m, "all", mt)
        assert x["ok"], x
        assert x["n_compared"] == 3 * m.geom.n_nodes - len(m.fixed_dofs())
        assert x["df"] < 1e-9, x["df"]
        assert 1.0 - x["mac_min"] < 1e-9 and 1.0 - x["submac_min"] < 1e-9
        assert x["n_clustered"] > 0               # symmetric pairs exist
    m.support_type = geo.ELASTIC
    x = ob.crosscheck_modal(m, 12)
    assert x["df"] < 1e-9 and 1.0 - x["mac_min"] < 1e-9


def test_von_mises_path_matches_opensees_and_closed_form():
    _need_opensees()
    m = bm.make_preset("von_mises")
    a, h, EA, P = bm.von_mises_parameters(m)
    x = ob.crosscheck_path(m, None, gdl_kwargs=dict(
        control_max=2.4 * h, lambda_max=1e9, lambda_min=-1e9,
        stop_after_limit=1000, max_steps=800))
    assert x["ok"], x["message"]
    assert x["n_matched"] == x["n_targets"] > 10
    assert x["dlam"] < 1e-8 and x["du"] < 1e-8
    assert x["res_os_in_gdl"] < 1e-9
    w = -np.asarray(x["os"]["control"])
    assert w.max() > 2.0 * h                      # past both limit points
    Pex = bm.von_mises_load(w, a, h, EA)
    Pos = np.asarray(x["os"]["lam"]) * P
    assert np.max(np.abs(Pos - Pex)) / np.max(np.abs(Pex)) < 1e-9
    lp = bm.von_mises_limit_point(a, h, EA)
    assert abs(x["lambda_limit_os"] * P - lp["P"]) / lp["P"] < 1e-9


def test_star_dome_through_both_snap_backs():
    _need_opensees()
    m = bm.make_preset("star_dome")
    x = ob.crosscheck_path(m, None, gdl_kwargs=dict(
        control_max=0.17, lambda_max=1e9, lambda_min=-1e9,
        stop_after_limit=1000, max_steps=800))
    assert x["ok"], x["message"]
    assert x["n_segments"] > 1                    # control had to switch
    assert x["dlam"] < 1e-8 and x["du"] < 1e-8
    assert x["res_os_in_gdl"] < 1e-9
    assert x["dlimit"] < 1e-10
    p_ea = x["lambda_limit_os"] * bm.STAR_LOAD_RATIO
    assert abs(p_ea - bm.STAR_P_OVER_EA_LIMIT) <= 5e-10


def test_exported_scripts_run_and_match():
    _need_opensees()
    m = _t_model()
    c = ob.find_combination(m, "T1")
    with tempfile.TemporaryDirectory() as tmp:
        r = ob.run_script(ob.write_openseespy_script(
            m, "T1", os.path.join(tmp, "s.py")))
        g = slv.solve_model(m, ld.combine(ld.basic_load_vectors(m), c),
                            T_factor=c.T)
        assert ob.rel_inf(g.u, r["u"]) < 1e-10
        assert ob.rel_inf(g.N, r["N"]) < 1e-10
        assert ob.rel_inf(g.reactions, r["reactions"]) < 1e-10
        r = ob.run_script(ob.write_openseespy_script(
            m, None, os.path.join(tmp, "m.py"),
            {"type": "modal", "n_modes": 8}))
        f = md.modal_analysis(m, n_modes=8).frequencies
        assert np.max(np.abs(np.asarray(r["frequencies"]) - f) / f) < 1e-9
        r = ob.run_script(ob.write_openseespy_script(
            m, "LC2", os.path.join(tmp, "n.py"),
            {"type": "path", "control_max": 0.3, "n_steps": 10}))
        system = nl.model_system(m)
        F = nl.reference_load(m, ob.find_combination(m, "LC2"))
        assert len(r["lam"]) == 10
        assert max(ob.gdl_residual(system, F, u, lam)
                   for u, lam in zip(r["u"], r["lam"])) < 1e-9


def test_inprocess_and_subprocess_agree():
    _need_opensees()
    m = _custom_model()
    spec, an = ob.static_job(m, "T1")
    results = []
    py = ob.find_opensees_python()
    if py:
        results.append(ob.run_opensees(spec, an, python=py,
                                       in_process=False))
    if ob.inprocess_available():
        results.append(ob.run_opensees(spec, an, in_process=True))
    if len(results) < 2:
        raise unittest.SkipTest("only one execution mode is available")
    a, b = results
    assert np.array_equal(a["u"], b["u"]) and np.array_equal(a["N"], b["N"])


# ================================================================ driver ===
def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = nskip = 0
    for name, fn in fns:
        try:
            fn()
        except unittest.SkipTest as exc:
            nskip += 1
            print(f"  skip  {name}: {exc}")
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        else:
            npass += 1
            print(f"  pass  {name}")
    print(f"\n{npass} of {len(fns)} tests passed"
          f"{f', {nskip} skipped' if nskip else ''}.")
    return 0 if npass + nskip == len(fns) else 1


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    raise SystemExit(_run_standalone())
