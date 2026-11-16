"""
Unit tests for GeoDomeLab.

Run with pytest::

    pytest -q

or standalone, without pytest installed::

    python tests/test_geodome.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import geometry as geo                                 # noqa: E402
from geodome import buckling as bk                                  # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import sections as sx                                  # noqa: E402
from geodome import verification as vf                              # noqa: E402
from geodome.materials import MATERIALS                             # noqa: E402
from geodome.model import DomeModel, LoadCombination, LoadSpec, \
    default_model                                                   # noqa: E402
from geodome.optimize import auto_size                              # noqa: E402
from geodome.postprocess import run_analysis, member_table, \
    node_table, reaction_table, summary_text                        # noqa: E402
from geodome.solver import solve_static, direction_cosines, \
    assemble_stiffness, free_dofs                                   # noqa: E402


# =========================================================== geometry =======
def test_icosahedron_is_regular():
    v, f = icosa = geo.icosahedron()
    v, f = icosa
    assert v.shape == (12, 3)
    assert f.shape == (20, 3)
    assert np.allclose(np.linalg.norm(v, axis=1), 1.0)
    edges = geo.edges_from_faces(f)
    assert len(edges) == 30
    L = np.linalg.norm(v[edges[:, 1]] - v[edges[:, 0]], axis=1)
    assert np.allclose(L, L[0])                      # all edges equal


def test_class1_topology():
    v, f = geo.icosahedron()
    for nu in (1, 2, 3, 4, 5, 6):
        pts, tris = (v, f) if nu == 1 else geo.subdivide_class1(v, f, nu)
        e = geo.edges_from_faces(tris)
        assert len(pts) == 10 * nu ** 2 + 2
        assert len(tris) == 20 * nu ** 2
        assert len(e) == 30 * nu ** 2
        assert len(pts) - len(e) + len(tris) == 2       # Euler
        assert np.allclose(np.linalg.norm(pts, axis=1), 1.0)


def test_chord_factors_match_published_tables():
    g = geo.generate_dome(2, 1.0, 0.5)
    cf = [s["chord_factor"] for s in g.strut_schedule()]
    assert abs(cf[0] - 0.546533) < 2e-5
    assert abs(cf[1] - 0.618034) < 2e-5              # 1/golden ratio
    g = geo.generate_dome(4, 1.0, 0.5)
    cf = [s["chord_factor"] for s in g.strut_schedule()]
    for got, ref in zip(cf, [0.253185, 0.294531, 0.295242,
                             0.298588, 0.312869, 0.324920]):
        assert abs(got - ref) < 2e-5


def test_hemisphere_dimensions():
    R = 9.0
    g = geo.generate_dome(6, R, 0.5)
    assert abs(g.height - R) < 1e-9
    assert abs(g.base_diameter - 2 * R) < 1e-9
    assert abs(g.nodes[:, 2].min()) < 1e-9           # base sits on z = 0
    assert g.surface_area < 2 * math.pi * R ** 2     # faceted < smooth
    assert g.surface_area > 0.97 * 2 * math.pi * R ** 2


def test_every_truncation_is_stable():
    """No dome in the usable range may contain an unsupported mechanism."""
    for nu in (3, 4, 6):
        for frac in (0.25, 0.35, 0.5, 0.65, 0.8):
            g = geo.generate_dome(nu, 10.0, frac)
            deg = np.zeros(g.n_nodes, dtype=int)
            for i, j in g.members:
                deg[i] += 1
                deg[j] += 1
            free = np.setdiff1d(np.arange(g.n_nodes), g.base_nodes)
            assert len(g.base_nodes) >= 3
            # Every unsupported node needs at least three struts.
            assert deg[free].min() >= 3, (nu, frac, deg[free].min())


def test_boundary_nodes_are_the_free_edge():
    g = geo.generate_dome(4, 10.0, 0.5)
    bnd = geo.boundary_nodes(g.faces)
    assert set(bnd.tolist()) == set(g.base_nodes.tolist())
    assert np.allclose(g.nodes[bnd][:, 2], 0.0)      # hemisphere: level ring


def test_support_dofs():
    g = geo.generate_dome(4, 10.0, 0.5)
    pinned = geo.support_dofs(g, geo.PINNED)
    assert len(pinned) == 3 * len(g.base_nodes)
    roller = geo.support_dofs(g, geo.ROLLER3)
    assert len(roller) == len(g.base_nodes) + 6


# ============================================================= sections =====
def test_section_properties():
    s = sx.get_section("CHS 114.3 x 5.0")
    D, t = 0.1143, 0.005
    d = D - 2 * t
    assert abs(s.A - math.pi / 4 * (D ** 2 - d ** 2)) < 1e-15
    assert abs(s.I - math.pi / 64 * (D ** 4 - d ** 4)) < 1e-18
    assert abs(s.r - math.sqrt(s.I / s.A)) < 1e-15


def test_catalogue_is_sorted_by_area():
    areas = [s.A for s in sx.CATALOGUE]
    assert areas == sorted(areas)


# =============================================================== solver =====
def test_single_bar_closed_form():
    E, A, L, P = 200e9, 2e-3, 3.0, 50e3
    nodes = np.array([[0, 0, 0], [L, 0, 0]], float)
    r = solve_static(nodes, np.array([[0, 1]]), np.array([E * A]),
                     np.array([0, 0, 0, P, 0, 0.0]), [0, 1, 2, 4, 5])
    assert abs(r.u[3] - P * L / (E * A)) < 1e-18
    assert abs(r.N[0] - P) < 1e-6


def test_tripod_closed_form():
    E, A, a, h, P = 200e9, 2e-3, 2.0, 4.0, 80e3
    ang = np.deg2rad([0.0, 120.0, 240.0])
    nodes = np.vstack([np.column_stack([a * np.cos(ang), a * np.sin(ang),
                                        np.zeros(3)]), [0, 0, h]])
    F = np.zeros(12)
    F[11] = -P
    r = solve_static(nodes, np.array([[0, 3], [1, 3], [2, 3]]),
                     np.full(3, E * A), F, list(range(9)))
    L = math.hypot(a, h)
    assert abs(r.N[0] - (-P * L / (3 * h))) < 1e-6
    assert abs(r.u[11] - (-P * L ** 3 / (3 * E * A * h ** 2))) < 1e-15
    assert abs(r.total_reaction()[2] - P) < 1e-6


def test_stiffness_matrix_is_symmetric_positive_semidefinite():
    g = geo.generate_dome(3, 8.0, 0.5)
    EA = np.full(g.n_members, 200e9 * 1e-3)
    K = assemble_stiffness(g.nodes, g.members, EA)
    Kd = K.toarray() if hasattr(K, "toarray") else K
    assert np.allclose(Kd, Kd.T)
    w = np.linalg.eigvalsh(Kd)
    assert w.min() > -1e-8 * abs(w).max()            # no negative stiffness


def test_closed_sphere_is_maxwell_isostatic():
    """A complete Class I sphere satisfies m = 3n - 6 exactly.

    n = 10v^2 + 2 and m = 30v^2, so 3n - 6 = 30v^2 = m for every frequency:
    the closed geodesic sphere is statically determinate in Maxwell's count.
    A truncated dome is not -- it is short of struts and relies on its
    supports, which is exactly what the next test confirms.
    """
    v, f = geo.icosahedron()
    for nu in (1, 2, 3, 4, 6):
        pts, tris = (v, f) if nu == 1 else geo.subdivide_class1(v, f, nu)
        m = len(geo.edges_from_faces(tris))
        n = len(pts)
        assert m == 3 * n - 6


def test_supported_dome_has_no_mechanism():
    """The free partition of the stiffness matrix must be positive definite.

    An open triangulated cap floating in space *is* a mechanism (Maxwell's
    count falls short once the closed sphere is truncated), so this is the
    check that matters: once the free edge is supported, the remaining system
    has no zero-energy mode and is not merely well conditioned by luck.
    """
    for nu in (2, 3, 4, 6):
        for support in (geo.PINNED, geo.ROLLER3):
            g = geo.generate_dome(nu, 10.0, 0.5)
            EA = np.full(g.n_members, 200e9 * 1e-3)
            K = assemble_stiffness(g.nodes, g.members, EA)
            Kd = K.toarray() if hasattr(K, "toarray") else K
            fd = free_dofs(3 * g.n_nodes, geo.support_dofs(g, support))
            w = np.linalg.eigvalsh(Kd[np.ix_(fd, fd)])
            assert w.min() > 1e-10 * w.max(), (nu, support, w.min() / w.max())


def test_equilibrium_and_energy_identities():
    m = default_model(4, 10.0, 0.5)
    res = run_analysis(m, do_buckling=False)
    for r in res.cases.values():
        assert r.equilibrium_residual < 1e-9
        assert r.energy_error < 1e-8
        assert np.max(np.abs(r.reaction_balance)) < 1e-5


def test_direction_cosines_are_unit():
    g = geo.generate_dome(4, 10.0, 0.5)
    L, C = direction_cosines(g.nodes, g.members)
    assert np.allclose(np.linalg.norm(C, axis=1), 1.0)
    assert np.all(L > 0)


# ================================================================ loads =====
def test_load_resultants():
    m = default_model(4, 10.0, 0.5)
    m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                       snow_pressure=1000.0, wind_enabled=False)
    S = ld.snow_load(m)
    fx, fy, fz = ld.load_totals(S)
    assert abs(fx) < 1e-6 and abs(fy) < 1e-6
    # Snow acts on the projected area, so the resultant is p times the plan area.
    assert abs(fz + 1000.0 * m.geom.plan_area) < 1e-6


def test_self_weight_equals_mass_times_g():
    m = default_model(4, 10.0, 0.5)
    _, _, fz = ld.load_totals(ld.self_weight(m))
    assert abs(-fz - m.steel_mass() * 9.80665) < 1e-6


def test_wind_is_mostly_uplift_on_a_dome():
    m = default_model(4, 10.0, 0.5)
    m.loads.wind_enabled = True
    fx, fy, fz = ld.load_totals(ld.wind_load(m))
    assert fz > 0.0                      # net uplift
    assert fx > 0.0                      # net drag downwind
    Cp, theta = ld.wind_pressure_coefficients(m)
    assert Cp.max() <= 0.80 + 1e-9 and Cp.min() >= -1.00 - 1e-9


def test_combination_is_linear():
    m = default_model(3, 8.0, 0.5)
    basic = ld.basic_load_vectors(m)
    c = LoadCombination("t", D=1.2, S=1.6)
    F = ld.combine(basic, c)
    assert np.allclose(F, 1.2 * basic["D"] + 1.6 * basic["S"])


# ======================================================== code checks =======
def test_aisc_branches():
    mat = MATERIALS["ASTM A500 Gr. C (round)"]
    s = sx.get_section("CHS 114.3 x 6.3")
    _, slend, Fe, Fcr, _ = bk.compression_capacity(s, mat, 12.0)
    assert slend > 4.71 * math.sqrt(mat.E / mat.Fy)
    assert abs(Fcr - 0.877 * Fe) < 1e-6
    _, slend, Fe, Fcr, _ = bk.compression_capacity(s, mat, 2.0)
    assert slend < 4.71 * math.sqrt(mat.E / mat.Fy)
    assert abs(Fcr - (0.658 ** (mat.Fy / Fe)) * mat.Fy) < 1e-6
    assert Fcr < mat.Fy


def test_capacity_decreases_with_length():
    mat = MATERIALS["EN S355"]
    s = sx.get_section("CHS 88.9 x 4.0")
    caps = [bk.compression_capacity(s, mat, L)[0] for L in (1, 2, 3, 4, 5)]
    assert all(caps[i] > caps[i + 1] for i in range(len(caps) - 1))


def test_buckling_scaling_is_linear():
    m = default_model(3, 8.0, 0.5)
    res = run_analysis(m, do_buckling=False)
    N = res.cases[res.order[0]].N
    b1 = bk.linear_buckling(m, N, 3)
    b2 = bk.linear_buckling(m, 2.0 * N, 3)
    assert abs(b2.load_factors[0] - b1.load_factors[0] / 2.0) < 1e-8


def test_buckling_modes_are_zero_at_supports():
    m = default_model(4, 10.0, 0.5)
    res = run_analysis(m, n_buckling_modes=3)
    fixed = m.fixed_dofs()
    for k in range(len(res.buckling.load_factors)):
        assert np.allclose(res.buckling.modes[fixed, k], 0.0)


# ============================================================== model =======
def test_model_round_trip(tmp_path=None):
    import tempfile
    m = default_model(4, 11.0, 0.5)
    m.group_sections = ["CHS 76.1 x 4.0"] * m.geom.n_strut_types()
    m.material_name = "EN S355"
    m.loads.snow_pressure = 1800.0
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.json")
        m.save(p)
        m2 = DomeModel.load(p)
    assert m2.geom.frequency == m.geom.frequency
    assert m2.geom.n_nodes == m.geom.n_nodes
    assert m2.group_sections == m.group_sections
    assert m2.material_name == "EN S355"
    assert abs(m2.loads.snow_pressure - 1800.0) < 1e-9
    assert abs(m2.steel_mass() - m.steel_mass()) < 1e-6


def test_group_sections_resized_with_geometry():
    m = DomeModel(geom=geo.generate_dome(2, 10.0, 0.5),
                  group_sections=["CHS 60.3 x 4.0"])
    assert len(m.group_sections) == m.geom.n_strut_types()


# ========================================================= postprocess ======
def test_analysis_tables_and_summary():
    m = default_model(4, 10.0, 0.5)
    res = run_analysis(m, n_buckling_modes=3)
    assert len(member_table(res)) == m.geom.n_members
    assert len(node_table(res)) == m.geom.n_nodes
    assert len(reaction_table(res)) == len(m.geom.base_nodes) + 1
    assert "ANALYSIS SUMMARY" in summary_text(res)
    assert 0.0 <= res.governing_util
    assert res.governing_case in res.order
    # The envelope must bound every individual case.
    for name in res.strength_cases():
        assert np.all(res.util[name] <= res.util_env + 1e-12)


def test_vertical_reaction_balances_applied_load():
    m = default_model(4, 10.0, 0.5)
    res = run_analysis(m, do_buckling=False)
    for r in res.cases.values():
        applied = r.F.reshape(-1, 3).sum(axis=0)
        reac = r.total_reaction()
        assert np.max(np.abs(applied + reac)) < 1e-5


# =========================================================== optimize =======
def test_auto_size_reduces_mass_and_respects_limits():
    m = default_model(4, 12.0, 0.5)
    m.group_sections = ["CHS 168.3 x 8.0"] * m.geom.n_strut_types()
    start = m.steel_mass()
    s = auto_size(m, target_util=0.95, max_iter=10)
    assert m.steel_mass() < start
    assert abs(s.final_mass - m.steel_mass()) < 1e-6
    res = run_analysis(m, do_buckling=False)
    assert res.governing_util <= 1.0
    caps = res.caps
    assert caps["slenderness"].max() <= 200.0 + 1e-6


# ======================================================== verification ======
def test_full_verification_suite_passes():
    rows = vf.run_all()
    failed = [r["check"] for r in rows if not r["pass"]]
    assert not failed, f"failed checks: {failed}"


# =============================================================== driver =====
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
