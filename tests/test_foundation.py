"""
Tests for the GeoDomeLab 2.0 data-model foundation (spec task S1).

Covered here, each against a closed-form or hand-derived answer:

* thermal load (initial strain alpha*dT) on fixed-fixed and two-segment bars,
  and on validation Model A (1V cap) through the full model pipeline;
* bar + spring in series, spring forces in the reactions, spring energy;
* support settlement of a statically determinate two-bar truss (rigid-body
  rotation, N = 0) and of an indeterminate two-bar line (N from compatibility);
* superposition of load + temperature + settlement, and the generalised
  equilibrium / energy identities in every case;
* custom (general) truss geometry, JSON version 2 round trip, version-1 load;
* ELASTIC and CUSTOM supports, including the rigid limits that must equal the
  PINNED model exactly, and springs in the buckling eigenproblem;
* nodal loads, the T factor of combinations, materials' alpha;
* the GUI still opens, generates, and shows a custom model.

Run with pytest::

    pytest -q tests/test_foundation.py

or standalone, without pytest installed::

    python tests/test_foundation.py
"""

from __future__ import annotations

import json
import math
import os
import sys
import tempfile

import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

from geodome import geometry as geo                                 # noqa: E402
from geodome import buckling as bk                                  # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import solver as slv                                   # noqa: E402
from geodome.materials import MATERIALS, Material                   # noqa: E402
from geodome.model import DomeModel, LoadCombination, LoadSpec, \
    default_model, default_combinations                             # noqa: E402
from geodome.postprocess import run_analysis, member_table, \
    node_table, reaction_table, summary_text                        # noqa: E402
from geodome.solver import solve_static, solve_model               # noqa: E402

E = 200e9
A = 1.0e-3
EA = E * A
ALPHA = 11.7e-6


def _check_identities(r, tol=1e-10):
    """Equilibrium, reaction balance and energy identity of one solve."""
    assert r.equilibrium_residual < tol, r.equilibrium_residual
    scale = max(np.abs(r.F).max(), np.abs(r.reactions).max(), 1.0)
    assert np.max(np.abs(r.reaction_balance)) < 1e-10 * scale, \
        r.reaction_balance
    assert r.energy_error < tol, r.energy_error


def _fix_all_but(n_nodes, free):
    """Restrained DOFs = everything except the listed free DOFs."""
    return sorted(set(range(3 * n_nodes)) - set(free))


# ========================================================== materials =======
def test_material_alpha_values():
    assert MATERIALS["ASTM A500 Gr. C (round)"].alpha == 11.7e-6
    assert MATERIALS["ASTM A500 Gr. B (round)"].alpha == 11.7e-6
    assert MATERIALS["ASTM A53 Gr. B pipe"].alpha == 11.7e-6
    assert MATERIALS["EN S235"].alpha == 12e-6
    assert MATERIALS["EN S355"].alpha == 12e-6
    assert MATERIALS["Aluminium 6061-T6"].alpha == 23e-6
    # AISC: 6.5e-6 per deg F == 11.7e-6 per deg C
    from geodome.units import ALPHA_PER_F_TO_PER_C
    assert abs(6.5e-6 * ALPHA_PER_F_TO_PER_C - 11.7e-6) < 1e-18
    # review V04: the database values against the specified values
    # (ASTM A500 Gr. B / C round 317 / 345 MPa, A53 Gr. B 240 MPa, E 200 GPa;
    # EN 10025-2 S235 / S355 and EN 1993-1-1 E = 210 GPa; 6061-T6 E 69 GPa)
    spec = {"ASTM A500 Gr. B (round)": (200e9, 317e6),
            "ASTM A500 Gr. C (round)": (200e9, 345e6),
            "ASTM A53 Gr. B pipe": (200e9, 240e6),
            "EN S235": (210e9, 235e6), "EN S355": (210e9, 355e6),
            "Aluminium 6061-T6": (69e9, 240e6)}
    assert set(spec) == set(MATERIALS)
    for name, (E_, Fy_) in spec.items():
        assert MATERIALS[name].E == E_ and MATERIALS[name].Fy == Fy_, name
        assert MATERIALS[name].rho in (7850.0, 2700.0)
    # the old five-argument constructor still works
    m = Material("X", 200e9, 300e6, 400e6, 7850.0)
    assert m.alpha == 11.7e-6
    assert "alpha" in m.describe()


# ============================================================ thermal =======
def test_heated_fixed_fixed_bar():
    """N = -EA alpha dT, no displacement, reactions +-EA alpha dT."""
    dT = 40.0
    nodes = np.array([[0.0, 0, 0], [3.0, 0, 0], [6.0, 0, 0]])
    mem = np.array([[0, 1], [1, 2]])
    fixed = _fix_all_but(3, [3])                       # middle node, x only
    r = solve_static(nodes, mem, np.full(2, EA), np.zeros(9), fixed,
                     eps0=ALPHA * dT)
    N_exact = -EA * ALPHA * dT
    assert np.allclose(r.N, N_exact, rtol=1e-12, atol=0)
    assert abs(r.u[3]) < 1e-15
    # reaction = force of the wall on the bar: R_start = -N, R_end = +N
    # (a compressed bar pushes both walls outwards, the walls push back)
    assert abs(r.reactions[0] + N_exact) < 1e-9 * abs(N_exact)
    assert abs(r.reactions[6] - N_exact) < 1e-9 * abs(N_exact)
    _check_identities(r)
    # both sides of the energy identity vanish here (N e = 0) ...
    assert abs(r.strain_energy_ext) < 1e-12 * r.energy_scale
    # ... while the stored energy 1/2 N^2 L/EA is not zero
    assert r.energy_scale > 0.0


def test_heated_two_segment_bar_closed_form():
    """Unequal segments between rigid walls: compatibility e1 + e2 = 0.

        N  = -alpha dT (L1 + L2) / (L1/EA1 + L2/EA2)
        u1 = L1 (N/EA1 + alpha dT)
    """
    dT, L1, L2, EA1, EA2 = 35.0, 2.0, 5.0, 3.0e8, 1.2e8
    nodes = np.array([[0.0, 0, 0], [L1, 0, 0], [L1 + L2, 0, 0]])
    mem = np.array([[0, 1], [1, 2]])
    r = solve_static(nodes, mem, np.array([EA1, EA2]), np.zeros(9),
                     _fix_all_but(3, [3]), eps0=ALPHA * dT)
    N = -ALPHA * dT * (L1 + L2) / (L1 / EA1 + L2 / EA2)
    u1 = L1 * (N / EA1 + ALPHA * dT)
    assert np.allclose(r.N, N, rtol=1e-12)
    assert abs(r.u[3] - u1) < 1e-12 * abs(u1)
    _check_identities(r)


def test_free_bar_expands_without_force():
    dT, L = -25.0, 4.0                                 # cooling
    nodes = np.array([[0.0, 0, 0], [0, 0, L]])
    r = solve_static(nodes, np.array([[0, 1]]), np.array([EA]), np.zeros(6),
                     _fix_all_but(2, [5]), eps0=ALPHA * dT)
    assert abs(r.N[0]) < 1e-6
    assert abs(r.u[5] - ALPHA * dT * L) < 1e-15
    _check_identities(r)


def test_initial_strain_load_is_self_equilibrated():
    m = default_model(3, 8.0, 0.5)
    F0 = slv.initial_strain_load(m.geom.nodes, m.geom.members, m.EA(), 1e-4)
    assert np.allclose(F0.reshape(-1, 3).sum(axis=0), 0.0, atol=1e-6)


def _model_A(dT=0.0, settlements=None):
    """Validation Model A: 1V cap cut at the upper ring (crown + 5 pins)."""
    frac = 0.5 * (1.0 - 1.0 / math.sqrt(5.0))
    g = geo.generate_dome(1, 6.0, frac)
    m = DomeModel(geom=g, group_sections=["CHS 88.9 x 4.0"])
    m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                       snow_pressure=0.0, live_pressure=0.0,
                       wind_enabled=False, temperature_change=dT,
                       settlements=list(settlements or []))
    m.combinations = [LoadCombination("T only", T=1.0)]
    return m


def test_model_A_uniform_temperature():
    """Crown free expansion; ring bars between pins carry -EA alpha dT.

    The five crown bars (length l, rise h of the crown above the ring) are
    free to lengthen, so N = 0 and the crown rises by alpha dT l^2 / h.  The
    ring bars join two pinned nodes, so N = -EA alpha dT exactly.
    """
    dT = 30.0
    m = _model_A(dT)
    g = m.geom
    assert g.n_nodes == 6 and g.n_members == 10
    res = run_analysis(m, do_buckling=False)
    r = res.cases["T only"]
    c = g.crown_node
    alpha = m.material.alpha
    EAm = m.EA()
    crown_bars = np.array([c in mm for mm in g.members])
    l = float(g.member_lengths()[crown_bars][0])
    h = float(g.nodes[c, 2] - g.nodes[g.base_nodes, 2].mean())
    assert np.allclose(r.N[crown_bars], 0.0, atol=1e-6)
    assert np.allclose(r.N[~crown_bars], -EAm[~crown_bars] * alpha * dT,
                       rtol=1e-12)
    assert abs(r.u[3 * c + 2] - alpha * dT * l ** 2 / h) < 1e-12 * l
    assert np.allclose(r.u[3 * c:3 * c + 2], 0.0, atol=1e-15)
    _check_identities(r)
    # T factor scales linearly; a combination without T sees nothing
    m.combinations = [LoadCombination("half", T=0.5),
                      LoadCombination("none", D=1.0)]
    res2 = run_analysis(m, do_buckling=False)
    assert np.allclose(res2.cases["half"].N, 0.5 * r.N, rtol=1e-12, atol=1e-6)
    assert np.allclose(res2.cases["none"].N, 0.0)


# ============================================================ springs =======
def test_bar_and_spring_in_series():
    """Spring support (k) -> bar (EA, L) -> load P.

        u0 = P/k,  u1 = P (1/k + L/EA),  N = P,  spring force -P
        strain energy = P^2 L / (2 EA) + P^2 / (2 k)
    """
    P, k, L = 12.0e3, 4.0e6, 2.5
    nodes = np.array([[0.0, 0, 0], [L, 0, 0]])
    F = np.zeros(6)
    F[3] = P
    fixed = [1, 2, 4, 5]                               # x free at both nodes
    r = solve_static(nodes, np.array([[0, 1]]), np.array([EA]), F, fixed,
                     springs={0: k})
    assert abs(r.u[0] - P / k) < 1e-12 * P / k
    assert abs(r.u[3] - P * (1 / k + L / EA)) < 1e-12 * r.u[3]
    assert abs(r.N[0] - P) < 1e-9 * P
    assert abs(r.spring_forces[0] + P) < 1e-9 * P
    assert abs(r.reactions[0] + P) < 1e-9 * P          # spring force reported
    assert abs(r.strain_energy_springs - P ** 2 / (2 * k)) < 1e-9 * P ** 2 / k
    _check_identities(r)
    # total stored energy = work of the load (Clapeyron)
    U = 0.5 * P ** 2 * L / EA + 0.5 * P ** 2 / k
    assert abs(0.5 * P * r.u[3] - U) < 1e-12 * U


def test_heated_bar_against_spring():
    """Wall -> heated bar -> spring.

        EA (u1/L - alpha dT) = -k u1  ->  u1 = alpha dT EA / (EA/L + k),
        N = -k u1
    """
    dT, k, L = 50.0, 2.0e7, 3.0
    nodes = np.array([[0.0, 0, 0], [L, 0, 0]])
    r = solve_static(nodes, np.array([[0, 1]]), np.array([EA]), np.zeros(6),
                     [0, 1, 2, 4, 5], springs={3: k}, eps0=ALPHA * dT)
    u1 = ALPHA * dT * EA / (EA / L + k)
    assert abs(r.u[3] - u1) < 1e-12 * u1
    assert abs(r.N[0] + k * u1) < 1e-9 * k * u1
    _check_identities(r)


def test_spring_ground_movement():
    """A settling spring support: ground end moves by delta.

    Wall -> bar -> node -> spring -> ground (moved by delta), no load:
    u = k delta / (EA/L + k), N = EA u / L.
    """
    k, L, delta = 3.0e7, 2.0, -0.01
    nodes = np.array([[0.0, 0, 0], [L, 0, 0]])
    r = solve_static(nodes, np.array([[0, 1]]), np.array([EA]), np.zeros(6),
                     [0, 1, 2, 4, 5], springs={3: k}, prescribed={3: delta})
    u = k * delta / (EA / L + k)
    assert abs(r.u[3] - u) < 1e-12 * abs(u)
    assert abs(r.N[0] - EA * u / L) < 1e-9 * abs(EA * u / L)
    _check_identities(r)


# ========================================================= settlement =======
def test_settlement_two_bar_truss_determinate():
    """Two-bar truss A-C-B (span 2a, rise h), out-of-plane y restrained.

    Pin B settles by delta: the truss is statically determinate, so it turns
    as a rigid body about A with omega = delta / 2a:
        N = 0,   u_C = (delta h / 2a, 0, -delta / 2)
    """
    a, h, delta = 3.0, 2.0, 0.012
    nodes = np.array([[0.0, 0, 0], [2 * a, 0, 0], [a, 0, h]])
    mem = np.array([[0, 2], [1, 2]])
    fixed = [0, 1, 2, 3, 4, 5, 7]                      # C free in x and z
    r = solve_static(nodes, mem, np.full(2, EA), np.zeros(9), fixed,
                     prescribed={5: -delta})
    assert np.allclose(r.N, 0.0, atol=1e-6)
    assert abs(r.u[6] - delta * h / (2 * a)) < 1e-14
    assert abs(r.u[8] + delta / 2) < 1e-14
    assert abs(r.u[5] + delta) == 0.0                  # prescribed exactly
    _check_identities(r)


def test_settlement_and_temperature_indeterminate_line():
    """Two bars in line between walls; wall B moves by delta, bars heated.

    Compatibility: e1 + e2 = delta, with e_i = L_i (N/EA_i + alpha dT):
        N = (delta - alpha dT (L1 + L2)) / (L1/EA1 + L2/EA2)
    """
    L1, L2, EA1, EA2, delta, dT = 2.0, 3.0, 2e8, 1e8, 0.002, 20.0
    nodes = np.array([[0.0, 0, 0], [L1, 0, 0], [L1 + L2, 0, 0]])
    mem = np.array([[0, 1], [1, 2]])
    fixed = _fix_all_but(3, [3])
    flex = L1 / EA1 + L2 / EA2
    for d, t in ((delta, 0.0), (0.0, dT), (delta, dT)):
        r = solve_static(nodes, mem, np.array([EA1, EA2]), np.zeros(9), fixed,
                         prescribed={6: d}, eps0=ALPHA * t)
        N = (d - ALPHA * t * (L1 + L2)) / flex
        assert np.allclose(r.N, N, rtol=1e-12), (d, t, r.N, N)
        assert abs(r.reactions[6] - N) < 1e-9 * abs(N)     # R_end = +N
        _check_identities(r)


def test_superposition_load_temperature_settlement():
    """The response is linear: u(F, eps0, us) = u(F) + u(eps0) + u(us)."""
    m = default_model(3, 8.0, 0.5)
    g = m.geom
    nodes, mem, EAm = g.nodes, g.members, m.EA()
    fixed = m.fixed_dofs()
    F = ld.combine(ld.basic_load_vectors(m), LoadCombination("x", D=1.2, S=1.6))
    b0 = int(g.base_nodes[0])
    pres = {3 * b0 + 2: -0.01, 3 * int(g.base_nodes[2]): 0.004}
    e0 = ALPHA * 25.0
    r_all = solve_static(nodes, mem, EAm, F, fixed, prescribed=pres, eps0=e0)
    r_F = solve_static(nodes, mem, EAm, F, fixed)
    r_T = solve_static(nodes, mem, EAm, 0 * F, fixed, eps0=e0)
    r_S = solve_static(nodes, mem, EAm, 0 * F, fixed, prescribed=pres)
    for x in ("u", "N", "reactions"):
        tot = getattr(r_F, x) + getattr(r_T, x) + getattr(r_S, x)
        scale = np.abs(tot).max()
        assert np.max(np.abs(getattr(r_all, x) - tot)) < 1e-9 * scale, x
    for r in (r_all, r_T, r_S):
        _check_identities(r)


def test_settlement_through_model_and_free_dof_error():
    """LoadSpec.settlements x T factor; a free direction raises clearly."""
    m = _model_A(0.0, [{"node": 0, "uz": 0.0}])
    b = int(m.geom.base_nodes[0])
    m.loads.settlements = [{"node": b, "uz": -0.008}]
    r = solve_model(m, np.zeros(3 * m.geom.n_nodes), T_factor=1.0)
    assert r.u[3 * b + 2] == -0.008
    _check_identities(r)
    r2 = solve_model(m, np.zeros(3 * m.geom.n_nodes), T_factor=2.0)
    assert np.allclose(r2.N, 2.0 * r.N, rtol=1e-12, atol=1e-6)
    # crown is free: a "settlement" there is an input error
    m.loads.settlements = [{"node": m.geom.crown_node, "uz": -0.01}]
    try:
        solve_model(m, np.zeros(3 * m.geom.n_nodes), T_factor=1.0)
    except ValueError as exc:
        assert "neither restrained nor a spring" in str(exc)
    else:
        raise AssertionError("settlement of a free DOF was accepted")


# =========================================================== supports =======
def _same_result(r1, r2, rtol=1e-12):
    for x in ("u", "N", "reactions"):
        a, b = getattr(r1, x), getattr(r2, x)
        scale = max(np.abs(a).max(), 1e-30)
        assert np.max(np.abs(a - b)) <= rtol * scale, x


def test_elastic_rigid_limit_equals_pinned():
    m = default_model(4, 10.0, 0.5)
    res_p = run_analysis(m, n_buckling_modes=2)
    m.support_type = geo.ELASTIC
    m.spring_stiffness = (math.inf, math.inf, math.inf)
    fixed, springs = m.support_conditions()
    assert not springs
    assert np.array_equal(fixed, geo.support_dofs(m.geom, geo.PINNED))
    res_e = run_analysis(m, n_buckling_modes=2)
    for n in res_p.order:
        _same_result(res_p.cases[n], res_e.cases[n], rtol=0.0)
    assert np.array_equal(res_p.buckling.load_factors,
                          res_e.buckling.load_factors)


def test_elastic_supports_converge_to_pinned_and_balance():
    m = default_model(4, 10.0, 0.5)
    F = ld.combine(ld.basic_load_vectors(m), m.combinations[1])
    r_p = solve_model(m, F)
    m.support_type = geo.ELASTIC
    prev = None
    for k in (1e7, 1e9, 1e11):
        m.spring_stiffness = (k, k, k)
        r = solve_model(m, F)
        _check_identities(r, tol=1e-9)
        fixed, springs = m.support_conditions()
        assert len(fixed) == 0 and len(springs) == 3 * len(m.geom.base_nodes)
        sd = np.array(sorted(springs))
        assert np.allclose(r.spring_forces[sd], -k * r.u[sd], rtol=1e-12)
        assert np.allclose(r.reactions[sd], r.spring_forces[sd])
        err = np.max(np.abs(r.N - r_p.N)) / np.max(np.abs(r_p.N))
        if prev is not None:
            assert err < prev                           # monotone approach
        prev = err
    assert prev < 1e-3
    # softer supports -> larger crown deflection
    m.spring_stiffness = (1e6, 1e6, 1e6)
    r_soft = solve_model(m, F)
    c = m.geom.crown_node
    assert abs(r_soft.u[3 * c + 2]) > abs(r_p.u[3 * c + 2])


def test_custom_supports_match_pinned_and_mix_springs():
    m = default_model(3, 8.0, 0.5)
    F = ld.combine(ld.basic_load_vectors(m), m.combinations[1])
    r_p = solve_model(m, F)
    m.support_type = geo.CUSTOM
    m.custom_supports = {int(b): [True, True, True] for b in m.geom.base_nodes}
    assert np.array_equal(m.fixed_dofs(), geo.support_dofs(m.geom, geo.PINNED))
    _same_result(r_p, solve_model(m, F), rtol=0.0)
    # mix: one node a vertical spring, another free horizontally
    b = [int(x) for x in m.geom.base_nodes]
    m.custom_supports[b[0]] = [True, True, 2.0e7]
    m.custom_supports[b[1]] = [False, False, True]
    fixed, springs = m.support_conditions()
    assert springs == {3 * b[0] + 2: 2.0e7}
    assert 3 * b[1] not in fixed and 3 * b[1] + 2 in fixed
    r = solve_model(m, F)
    _check_identities(r)
    assert abs(r.reactions[3 * b[1]]) == 0.0
    assert abs(r.reactions[3 * b[0] + 2] + 2.0e7 * r.u[3 * b[0] + 2]) < 1e-6
    assert set(m.support_nodes().tolist()) == set(b)
    rows = reaction_table(run_analysis(m, do_buckling=False))
    assert len(rows) == len(b) + 1
    # an empty table is an error with a helpful message
    m.custom_supports = {}
    try:
        m.support_conditions()
    except ValueError as exc:
        assert "custom_supports" in str(exc)
    else:
        raise AssertionError("empty CUSTOM table accepted")


def _expect_mechanism(fn, text):
    try:
        fn()
    except np.linalg.LinAlgError as exc:
        assert text in str(exc), str(exc)
        assert "Unstable structure" in str(exc), str(exc)
    else:
        raise AssertionError("a mechanism was solved without an error")


def test_mechanisms_are_refused():
    """Review F1 / GUI-02: supports that let the structure move as a rigid
    body (non-zero diagonal, singular K) must raise, not return garbage."""
    # Elastic supports with kx = ky = 0: slides in x, y and turns about z.
    m = default_model(4, 12.0, 0.5)
    m.support_type = geo.ELASTIC
    m.spring_stiffness = (0.0, 0.0, 5e7)
    fixed, springs = m.support_conditions()
    assert slv.unrestrained_rigid_modes(
        m.geom.nodes, list(fixed) + list(springs)) == 3
    _expect_mechanism(lambda: run_analysis(m, do_buckling=False),
                      "3 rigid-body motion")
    # every spring zero: nothing is supported at all
    m.spring_stiffness = (0.0, 0.0, 0.0)
    _expect_mechanism(lambda: run_analysis(m, do_buckling=False),
                      "6 rigid-body motion")
    # the star dome on vertical rollers (only uz held)
    from geodome import benchmarks as bm
    s = bm.star_dome()
    s.support_type = geo.CUSTOM
    s.custom_supports = {n: [False, False, True] for n in s.geom.base_nodes}
    _expect_mechanism(lambda: run_analysis(s, do_buckling=False),
                      "rigid-body")
    # the stable supports still solve exactly
    for st_ in (geo.PINNED, geo.ROLLER3):
        m = default_model(3, 8.0, 0.5)
        m.support_type = st_
        fixed, springs = m.support_conditions()
        assert slv.unrestrained_rigid_modes(m.geom.nodes, fixed) == 0
        res = run_analysis(m, do_buckling=False)
        assert max(r.equilibrium_residual for r in res.cases.values()) < 1e-10
    # a collinear bar fixed at all DOFs but its axis is legal (rank R = 5)
    nodes = np.array([[0.0, 0, 0], [2.0, 0, 0]])
    r = solve_static(nodes, np.array([[0, 1]]), np.array([EA]),
                     np.array([0, 0, 0, 1e3, 0, 0.0]), [0, 1, 2, 4, 5])
    assert abs(r.N[0] - 1e3) < 1e-9


def test_internal_mechanism_is_refused():
    """A planar four-bar parallelogram (no diagonal) is held against rigid
    motion but sways under a lateral load.  The singular solve must raise,
    with the frame on the axes (exactly singular) and turned by 17 deg
    (singular up to round-off: caught by the residual check)."""
    base = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0],
                     [3.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    mem = np.array([[0, 3], [1, 2], [2, 3]])
    for ang in (0.0, 17.0):
        c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        nodes = base @ rot.T
        fixed = [0, 1, 2, 3, 4, 5, 8, 11]          # 0, 1 pinned; plane z
        assert slv.unrestrained_rigid_modes(nodes, fixed) == 0
        F = np.zeros(12)
        F[6:8] = 1e4 * rot[:2, 0]                   # push node 2 sideways
        _expect_mechanism(
            lambda: solve_static(nodes, mem, np.full(3, EA), F, fixed),
            "mechanism")
        # a diagonal makes it a stable truss again
        r = solve_static(nodes, np.vstack([mem, [0, 2]]), np.full(4, EA), F,
                         fixed)
        _check_identities(r)


def test_default_combinations_follow_asce7_2_3_1():
    """Review F5 / V01: the default strength set has ASCE 7-16 Sec. 2.3.1
    combination 3 with 0.5W and combination 2 with the live load L."""
    want = [("LC1  1.4D", dict(D=1.4)),
            ("LC2  1.2D + 1.6S", dict(D=1.2, S=1.6)),
            ("LC3  1.2D + 1.0W + 0.5S", dict(D=1.2, W=1.0, S=0.5)),
            ("LC4  0.9D + 1.0W", dict(D=0.9, W=1.0)),
            ("LC5  1.2D + 1.6S + 0.5W", dict(D=1.2, S=1.6, W=0.5)),
            ("LC6  1.2D + 1.6L + 0.5S", dict(D=1.2, L=1.6, S=0.5)),
            ("SLS  1.0D + 1.0S", dict(D=1.0, S=1.0))]
    got = default_combinations()
    assert [c.name for c in got] == [n for n, _ in want]
    for c, (_, f) in zip(got, want):
        full = {k: f.get(k, 0.0) for k in ("D", "S", "L", "W", "T")}
        assert c.factors() == full, (c.name, c.factors())
        assert c.serviceability == c.name.startswith("SLS")
    from geodome import csv_import
    names = [c.name for c in csv_import.import_combinations()]
    assert any("0.5W" in n for n in names), names


def test_unused_basic_load_is_reported():
    """Review F2: a non-zero load case that no strength combination uses
    must be named in the run messages, not dropped silently."""
    m = default_model(3, 8.0, 0.5)
    m.loads.live_pressure = 2000.0
    m.loads.crown_point_load = 50e3
    res = run_analysis(m, do_buckling=False)
    assert not any("no strength combination has" in t for t in res.messages)
    # the live load acts in LC6: 1.2D + 1.6L + 0.5S
    F = ld.combine(ld.basic_load_vectors(m), m.combinations[5])
    assert abs(ld.load_totals(F)[2] - ld.load_totals(
        1.2 * ld.dead_load(m) + 1.6 * ld.live_load(m)
        + 0.5 * ld.snow_load(m))[2]) < 1e-6
    # combinations without L: the live load is reported as unused
    m.combinations = [c for c in m.combinations if c.L == 0.0]
    res = run_analysis(m, do_buckling=False)
    assert any("live load" in t and "no strength combination has an L"
               in t for t in res.messages), res.messages
    # an L nodal load alone is reported too
    m.loads.live_pressure = 0.0
    m.loads.crown_point_load = 0.0
    m.loads.nodal_loads = [{"node": 0, "Fz": -300e3}]       # case L
    res = run_analysis(m, do_buckling=False)
    assert any("live load" in t for t in res.messages)
    # wind used only in a service row is reported as well
    m.loads.nodal_loads = []
    m.combinations = [LoadCombination("U", D=1.2, S=1.6),
                      LoadCombination("SV", D=1.0, W=1.0,
                                      serviceability=True)]
    res = run_analysis(m, do_buckling=False)
    assert any("wind load" in t for t in res.messages)
    assert not any("snow load" in t for t in res.messages)


def test_aluminium_member_checks_are_flagged():
    """Review V07: AISC 360 is a steel code; aluminium runs must say so."""
    m = default_model(3, 8.0, 0.5)
    res = run_analysis(m, do_buckling=False)
    assert not any("not valid for" in t for t in res.messages)
    assert MATERIALS[m.material_name].aisc_checks_apply
    m.material_name = "Aluminium 6061-T6"
    assert not m.material.aisc_checks_apply
    res = run_analysis(m, do_buckling=False)
    assert any("Member checks use AISC 360" in t and "Aluminium" in t
               for t in res.messages), res.messages
    assert "not valid for Aluminium 6061-T6" in summary_text(res)


def test_wind_azimuth_is_the_direction_the_wind_blows_towards():
    """Review F7 / V13: azimuth 0 = wind from -x, so the windward base node
    (theta = 0, the largest Cp) is at x = -R; azimuth 90 puts it at y = -R."""
    R = 10.0
    m = default_model(4, R, 0.5)
    base = m.geom.base_nodes
    for az, where in ((0.0, (-R, 0.0)), (90.0, (0.0, -R))):
        m.loads.wind_azimuth = az
        Cp, theta = ld.wind_pressure_coefficients(m)
        k = base[int(np.argmax(Cp[base]))]
        assert np.allclose(m.geom.nodes[k, :2], where, atol=1e-9), \
            (az, m.geom.nodes[k])
        assert abs(theta[k]) < 1e-6 and abs(Cp[k] - 0.80) < 1e-12
        fx, fy, _ = ld.load_totals(ld.wind_load(m))
        drag = fx * math.cos(math.radians(az)) + fy * math.sin(math.radians(az))
        assert drag > 0.0                    # drag along the wind direction


def test_unbalanced_snow_pattern_on_model_B():
    """Review V08: unbalanced snow is 1.5 p on x >= 0 and 0.5 p on x < 0,
    on projected tributary areas computed here from the panels."""
    from geodome import benchmarks as bm
    m = bm.make_preset("validation_b")
    p = 1000.0
    m.loads.snow_pressure = p
    m.loads.snow_unbalanced = True
    g = m.geom
    # one third of the plan area of each adjacent upward-facing panel
    # (outward = away from the sphere centre; Model B is a super-hemisphere
    # and 5 of its band panels face down)
    centre = np.array([0.0, 0.0, g.nodes[:, 2].max() - g.radius])
    Ap = np.zeros(g.n_nodes)
    n_down = 0
    for a, b, c in g.faces:
        cr = np.cross(g.nodes[b] - g.nodes[a], g.nodes[c] - g.nodes[a])
        cen = (g.nodes[a] + g.nodes[b] + g.nodes[c]) / 3.0
        if cr @ (cen - centre) < 0.0:
            cr = -cr
        n_down += cr[2] < 0.0
        for k in (a, b, c):
            Ap[k] += max(cr[2], 0.0) / 6.0
    assert n_down == 5
    fac = np.where(g.nodes[:, 0] >= 0.0, 1.5, 0.5)
    F = ld.snow_load(m)
    assert np.allclose(F[2::3], -p * Ap * fac, rtol=1e-12, atol=1e-9)
    assert np.all(F[0::3] == 0.0) and np.all(F[1::3] == 0.0)
    assert fac[g.crown_node] == 1.5                 # crown (x = 0) is heavy
    tot = -p * (1.5 * Ap[g.nodes[:, 0] >= 0].sum()
                + 0.5 * Ap[g.nodes[:, 0] < 0].sum())
    assert abs(ld.load_totals(F)[2] - tot) < 1e-9 * abs(tot)


def test_wind_vertex_lumping_error_is_as_documented():
    """Review V09: with Cp = 1 everywhere the vertex-normal lumping gives
    sum(A_i n_i) = plan area x (0.949, 0.986, 0.994, 0.996) at 2V-8V,
    while the faceted surface gives exactly the plan area (documented in
    loads.py and the theory text)."""
    doc = {2: 0.94876, 4: 0.98616, 6: 0.99376, 8: 0.99647}
    for nu, ratio in doc.items():
        g = geo.generate_dome(nu, 10.0, 0.5)
        lump = (g.trib_area[:, None] * g.normals).sum(axis=0)
        exact = np.zeros(3)
        c0 = np.array([0.0, 0.0, 0.0])
        for a, b, c in g.faces:
            cr = np.cross(g.nodes[b] - g.nodes[a], g.nodes[c] - g.nodes[a])
            cen = (g.nodes[a] + g.nodes[b] + g.nodes[c]) / 3.0
            if cr @ (cen - c0) < 0:
                cr = -cr
            exact += 0.5 * cr
        assert abs(exact[2] - g.plan_area) < 1e-9 * g.plan_area
        assert abs(lump[0]) < 1e-9 and abs(lump[1]) < 1e-9
        assert abs(lump[2] / g.plan_area - ratio) < 2e-5, (nu, lump[2] /
                                                           g.plan_area)


def test_super_hemisphere_snow_only_on_upward_panels():
    """Review F6: panels below the equator face down and carry no snow; the
    plan area is the upward-facing projection (hemispheres unchanged)."""
    R, p = 10.0, 1000.0
    half = geo.generate_dome(4, R, 0.5)
    _, old = geo._tributary_areas(half.nodes, half.faces)
    assert np.array_equal(old, half.trib_area_proj)  # |n_z| = n_z up there
    for frac in (0.626, 0.724, 0.925):
        g = geo.generate_dome(4, R, frac)
        assert g.dome_fraction > 0.55
        # the 4V mesh has an equator ring, so the upward panels are exactly
        # the hemisphere's panels
        assert abs(g.plan_area - half.plan_area) < 1e-9 * half.plan_area
        assert g.plan_area <= math.pi * R ** 2
        m = DomeModel(geom=g)
        m.loads.snow_pressure = p
        assert abs(ld.load_totals(ld.snow_load(m))[2]
                   + p * g.plan_area) < 1e-9 * p * g.plan_area
        zc = g.nodes[:, 2].max() - R                # sphere centre height
        low = g.nodes[:, 2] < zc - 1e-6
        assert np.any(low)
        assert np.all(g.trib_area_proj[low] == 0.0)  # below the equator


def test_custom_truss_text_has_no_nan():
    """Review GUI-11: a truss without panels prints 'no plan area', and
    its summary is headed SPACE TRUSS."""
    from geodome import benchmarks as bm
    m = bm.make_preset("three_bar")
    lines = "\n".join(m.summary_lines())
    assert "nan" not in lines.lower(), lines
    assert "no panels, so no plan area" in lines
    res = run_analysis(m, do_buckling=False)
    txt = summary_text(res)
    assert "SPACE TRUSS -- ANALYSIS SUMMARY" in txt
    assert "GEODESIC DOME" not in txt
    assert "kg/m2 of plan" in "\n".join(default_model(3, 8.0, 0.5)
                                        .summary_lines())


def test_support_types_listed():
    assert set(geo.SUPPORT_TYPES) == {geo.PINNED, geo.ROLLER3, geo.ELASTIC,
                                      geo.CUSTOM}
    m = default_model(2, 6.0, 0.5)
    m.support_type = geo.ELASTIC
    assert m.spring_stiffness == (5e7, 5e7, 5e7)
    assert "kN/mm" in "\n".join(m.summary_lines())


def test_buckling_with_springs():
    """Springs enter K: rigid limit = pinned, finite springs -> lower lambda.

    For the same reference forces, finite springs only remove stiffness from
    the pinned problem (a Rayleigh-quotient argument), so lambda_cr cannot
    rise.
    """
    m = default_model(4, 10.0, 0.5)
    N = run_analysis(m, do_buckling=False).cases["LC2  1.2D + 1.6S"].N
    lam_p = bk.linear_buckling(m, N, 3).load_factors
    m.support_type = geo.ELASTIC
    m.spring_stiffness = (5e6, 5e6, 5e6)
    b = bk.linear_buckling(m, N, 3)
    assert b.load_factors[0] <= lam_p[0] * (1 + 1e-9)
    # spring DOFs are free, so the mode may move them
    sd = sorted(m.support_conditions()[1])
    assert np.max(np.abs(b.modes[sd, 0])) > 0.0
    m.spring_stiffness = (1e14, 1e14, 1e14)
    lam_stiff = bk.linear_buckling(m, N, 3).load_factors
    assert abs(lam_stiff[0] - lam_p[0]) < 1e-3 * lam_p[0]


def test_buckling_warning_has_no_unsourced_knockdown():
    m = default_model(4, 10.0, 0.5)
    N = run_analysis(m, do_buckling=False).cases["LC2  1.2D + 1.6S"].N
    lam = bk.linear_buckling(m, N, 1).load_factors[0]
    b = bk.linear_buckling(m, N * lam / 1.5, 1)    # scale so lambda = 1.5
    assert abs(b.load_factors[0] - 1.5) < 1e-6
    text = " ".join(b.warnings)
    assert "imperfection" in text and "0.3" not in text
    assert "knock-down" not in text


# ===================================================== loads & combos =======
def test_combination_T_and_labels():
    c = LoadCombination("x", D=1.2, T=1.0)
    assert list(c.factors()) == ["D", "S", "L", "W", "T"]
    assert c.label() == "1.2D + 1T"
    # positional calls of the old signature still work
    c2 = LoadCombination("y", 1.2, 1.6, 0.0, 0.0, True)
    assert c2.serviceability and c2.T == 0.0
    assert [x.name for x in default_combinations()][0] == "LC1  1.4D"
    assert all(x.T == 0.0 for x in default_combinations())
    m = default_model(2, 6.0, 0.5)
    basic = ld.basic_load_vectors(m)
    assert "T" not in basic
    F = ld.combine(basic, LoadCombination("z", D=1.0, T=5.0))
    assert np.array_equal(F, ld.combine(basic, LoadCombination("z", D=1.0)))


def test_nodal_loads_go_to_their_case():
    m = default_model(3, 8.0, 0.5)
    base = ld.basic_load_vectors(m)
    c = m.geom.crown_node
    m.loads.nodal_loads = [
        {"node": c, "Fz": -10e3, "case": "L"},
        {"node": c, "Fx": 2e3, "Fy": -1e3, "case": "W"},
        {"node": 0, "Fz": -500.0, "case": "D"},
        {"node": 1, "Fz": -1.0},                    # default case = L
    ]
    new = ld.basic_load_vectors(m)
    dL = new["L"] - base["L"]
    assert dL[3 * c + 2] == -10e3 and dL[3 * 1 + 2] == -1.0
    assert np.count_nonzero(dL) == 2
    dW = new["W"] - base["W"]
    assert abs(dW[3 * c] - 2e3) < 1e-9 and abs(dW[3 * c + 1] + 1e3) < 1e-9
    assert (new["D"] - base["D"])[2] == -500.0
    assert np.array_equal(new["S"], base["S"])
    m.loads.nodal_loads = [{"node": 0, "Fz": -1.0, "case": "T"}]
    try:
        ld.basic_load_vectors(m)
    except ValueError as exc:
        assert "case" in str(exc)
    else:
        raise AssertionError("nodal load in case T accepted")


def test_describe_loads_mentions_T():
    m = default_model(2, 6.0, 0.5)
    assert any("Self-straining T" in s and "none" in s
               for s in ld.describe_loads(m))
    m.loads.temperature_change = 30.0
    m.loads.settlements = [{"node": int(m.geom.base_nodes[0]), "uz": -0.005}]
    text = "\n".join(ld.describe_loads(m))
    assert "dT = +30.0 C" in text and "5.00 mm" in text


def test_loadspec_copy_is_deep():
    ls = LoadSpec(settlements=[{"node": 1, "uz": -0.01}],
                  nodal_loads=[{"node": 2, "Fz": -5.0, "case": "L"}])
    cp = ls.copy()
    cp.settlements[0]["uz"] = 99.0
    cp.nodal_loads.append({"node": 3})
    assert ls.settlements[0]["uz"] == -0.01 and len(ls.nodal_loads) == 1
    assert LoadSpec().settlements is not LoadSpec().settlements


def test_run_analysis_messages_for_T():
    m = _model_A(20.0)
    m.combinations = [LoadCombination("D", D=1.0)]
    res = run_analysis(m, do_buckling=False)
    assert any("no combination has a T factor" in s for s in res.messages)
    m.combinations = [LoadCombination("DT", D=1.0, T=1.0)]
    res = run_analysis(m, do_buckling=False)
    assert any("applied in 1 combination" in s for s in res.messages)


# ===================================================== custom geometry ======
def _star_dome():
    """24-member shallow star dome geometry (classic dimensions, mm -> m)."""
    nodes = [[0.0, 0.0, 0.0821]]
    for k in range(6):                              # inner hexagon
        a = math.radians(60.0 * k)
        nodes.append([0.25 * math.cos(a), 0.25 * math.sin(a), 0.0621])
    for k in range(6):                              # outer supports
        a = math.radians(60.0 * k + 30.0)
        nodes.append([0.4330 * math.cos(a), 0.4330 * math.sin(a), 0.0])
    members = []
    for k in range(6):
        members.append([0, 1 + k])
        members.append([1 + k, 1 + (k + 1) % 6])
        members.append([1 + k, 7 + k])
        members.append([1 + (k + 1) % 6, 7 + k])
    faces = []
    for k in range(6):
        faces.append([0, 1 + k, 1 + (k + 1) % 6])
        faces.append([1 + k, 7 + k, 1 + (k + 1) % 6])
    return np.array(nodes), np.array(members), np.array(faces)


def test_custom_geometry_attributes():
    nodes, members, faces = _star_dome()
    g = geo.custom_geometry(nodes, members, faces, title="Star dome")
    assert g.kind == "custom" and g.is_custom and g.title == "Star dome"
    assert g.frequency == 0 and math.isnan(g.dome_fraction)
    assert g.n_nodes == 13 and g.n_members == 24 and g.n_faces == 12
    assert g.crown_node == 0
    assert sorted(g.base_nodes.tolist()) == list(range(7, 13))
    assert np.isfinite(g.radius) and g.radius > 0.5     # shallow -> big R
    # tributary areas sum to the faceted area; projected to the plan area
    area = sum(0.5 * np.linalg.norm(np.cross(nodes[b] - nodes[a],
                                             nodes[c] - nodes[a]))
               for a, b, c in faces)
    assert abs(g.surface_area - area) < 1e-12
    assert g.plan_area < g.surface_area
    # outward normals point up at the crown
    assert g.normals[0, 2] > 0.99
    # families by length: radial, hexagon (both 0.25) and the two diagonals
    assert g.n_strut_types() == len(np.unique(np.round(g.strut_lengths, 6)))
    cf = [s["chord_factor"] for s in g.strut_schedule()]
    assert all(np.isfinite(cf))
    assert abs(g.base_diameter - 2 * 0.4330) < 1e-12


def test_custom_planar_truss_nan_radius_and_consumers():
    """A planar truss has no sphere: NaN radius, every consumer survives."""
    nodes = [[0, 0, 0], [4, 0, 0], [8, 0, 0], [2, 0, 3], [6, 0, 3]]
    members = [[0, 1], [1, 2], [0, 3], [1, 3], [1, 4], [2, 4], [3, 4]]
    g = geo.custom_geometry(nodes, members, title="Plane truss")
    assert math.isnan(g.radius)
    assert all(math.isnan(s["chord_factor"]) for s in g.strut_schedule())
    m = DomeModel(geom=g, support_type=geo.CUSTOM)
    m.custom_supports = {0: [True, True, True], 2: [False, True, True],
                         1: [False, True, False], 3: [False, True, False],
                         4: [False, True, False]}
    m.loads = LoadSpec(include_self_weight=True, dead_pressure=0.0,
                       snow_pressure=0.0, wind_enabled=False,
                       nodal_loads=[{"node": 3, "Fz": -50e3, "case": "L"}])
    m.combinations = [LoadCombination("1.2D+1.6L", D=1.2, L=1.6),
                      LoadCombination("SLS", D=1.0, L=1.0,
                                      serviceability=True)]
    text = "\n".join(m.summary_lines())
    assert "n/a" in text and "nan" not in text        # plan mass is NaN
    assert "no panels, so no plan area" in text
    res = run_analysis(m, n_buckling_modes=2)
    r = res.cases["1.2D+1.6L"]
    _check_identities(r)
    assert len(member_table(res)) == 7 and len(node_table(res)) == 5
    assert "ANALYSIS SUMMARY" in summary_text(res)
    # report export does not crash on NaN radius / zero plan area
    from geodome import report as rp
    with tempfile.TemporaryDirectory() as d:
        rp.export_csv(res, d)
        rp.export_html(res, os.path.join(d, "r.html"))
        assert os.path.getsize(os.path.join(d, "r.html")) > 1000


def test_custom_geometry_validation_errors():
    bad = [
        dict(nodes=[[0, 0, 0], [1, 0, 0]], members=[[0, 0]]),
        dict(nodes=[[0, 0, 0], [1, 0, 0]], members=[[0, 1], [1, 0]]),
        dict(nodes=[[0, 0, 0], [1, 0, 0]], members=[[0, 2]]),
        dict(nodes=[[0, 0, 0], [0, 0, 0]], members=[[0, 1]]),
        dict(nodes=[[0, 0, 0], [1, 0, 0]], members=[[0, 1]], groups=[0, 1]),
    ]
    for kw in bad:
        try:
            geo.custom_geometry(**kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {kw}")


def test_best_fit_sphere_exact():
    rng = np.random.default_rng(1)
    v = rng.normal(size=(40, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    c0, R0 = np.array([1.0, -2.0, 0.5]), 7.25
    c, R, rms = geo.best_fit_sphere(c0 + R0 * v)
    assert abs(R - R0) < 1e-10 and np.allclose(c, c0, atol=1e-10)
    assert rms < 1e-10
    # a generated dome's nodes give back its radius
    g = geo.generate_dome(4, 11.0, 0.5)
    _, R, _ = geo.best_fit_sphere(g.nodes)
    assert abs(R - 11.0) < 1e-9


# ================================================================ JSON ======
def test_custom_model_json_round_trip():
    nodes, members, faces = _star_dome()
    g = geo.custom_geometry(nodes, members, faces, title="Star dome")
    m = DomeModel(geom=g, material_name="EN S355",
                  group_sections=["CHS 60.3 x 4.0"] * g.n_strut_types(),
                  support_type=geo.CUSTOM)
    m.custom_supports = {int(b): [True, True, 3.0e7] for b in g.base_nodes}
    m.custom_supports[7] = [True, math.inf, False]
    m.spring_stiffness = (1e6, 2e6, 3e6)
    m.loads = LoadSpec(temperature_change=-15.0,
                       settlements=[{"node": 8, "ux": 0.001}],
                       nodal_loads=[{"node": 0, "Fz": -2e3, "case": "L"}])
    m.combinations = [LoadCombination("DLT", D=1.2, L=1.6, T=1.2),
                      LoadCombination("SLS", D=1.0, L=1.0,
                                      serviceability=True)]
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "star.json")
        m.save(p)
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)                    # strict JSON
        m2 = DomeModel.load(p)
    assert raw["version"] == 2 and raw["geometry"]["kind"] == "custom"
    assert raw["custom_supports"]["7"] == [True, True, False]   # inf -> fixed
    g2 = m2.geom
    assert g2.kind == "custom" and g2.title == "Star dome"
    assert np.array_equal(g2.nodes, g.nodes)
    assert np.array_equal(g2.members, g.members)
    assert np.array_equal(g2.faces, g.faces)
    assert np.array_equal(g2.base_nodes, g.base_nodes)
    assert np.array_equal(g2.strut_group, g.strut_group)
    assert m2.spring_stiffness == (1e6, 2e6, 3e6)
    assert m2.custom_supports[8] == [True, True, 3.0e7]
    assert m2.loads.temperature_change == -15.0
    assert m2.loads.settlements == [{"node": 8, "ux": 0.001}]
    assert m2.loads.nodal_loads == m.loads.nodal_loads
    assert m2.combinations[0].T == 1.2
    r1 = run_analysis(m, do_buckling=False)
    r2 = run_analysis(m2, do_buckling=False)
    for n in r1.order:
        _same_result(r1.cases[n], r2.cases[n], rtol=0.0)
        _check_identities(r1.cases[n])


def test_geodesic_model_json_v2_round_trip():
    m = default_model(4, 11.0, 0.5)
    m.support_type = geo.ELASTIC
    m.spring_stiffness = (4e7, 4e7, 9e7)
    m.loads.temperature_change = 25.0
    m.combinations.append(LoadCombination("T", D=1.0, T=1.0))
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.json")
        m.save(p)
        m2 = DomeModel.load(p)
    assert m2.to_dict()["geometry"]["kind"] == "geodesic"
    assert m2.support_type == geo.ELASTIC
    assert m2.spring_stiffness == (4e7, 4e7, 9e7)
    assert m2.combinations[-1].T == 1.0
    _same_result(run_analysis(m, do_buckling=False).cases["T"],
                 run_analysis(m2, do_buckling=False).cases["T"], rtol=0.0)


V1_MODEL = {
    "format": "GeoDomeLab model",
    "version": 1,
    "geometry": {"frequency": 3, "radius": 9.0, "dome_fraction": 0.5},
    "material": "EN S355",
    "group_sections": ["CHS 76.1 x 4.0"],
    "support_type": "roller3",
    "loads": {"include_self_weight": True, "dead_pressure": 350.0,
              "snow_pressure": 1500.0, "snow_unbalanced": False,
              "live_pressure": 0.0, "wind_enabled": True, "wind_speed": 40.0,
              "wind_Kz": 0.9, "wind_Kzt": 1.0, "wind_Kd": 0.95,
              "wind_G": 0.85, "wind_azimuth": 30.0, "crown_point_load": 0.0},
    "combinations": [{"name": "LC2  1.2D + 1.6S", "D": 1.2, "S": 1.6,
                      "L": 0.0, "W": 0.0, "serviceability": False},
                     {"name": "SLS", "D": 1.0, "S": 1.0, "L": 0.0, "W": 0.0,
                      "serviceability": True}],
    "design": {"K_effective": 1.0, "deflection_limit": 240.0,
               "phi_c": 0.9, "phi_t": 0.9},
}


def test_version1_file_loads_unchanged():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "v1.json")
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(V1_MODEL, fh)
        m = DomeModel.load(p)
    assert m.geom.kind == "geodesic" and m.geom.frequency == 3
    assert m.support_type == geo.ROLLER3
    assert m.loads.temperature_change == 0.0
    assert m.loads.settlements == [] and m.loads.nodal_loads == []
    assert all(c.T == 0.0 for c in m.combinations)
    assert m.spring_stiffness == geo.DEFAULT_SPRING_STIFFNESS
    assert m.custom_supports == {}
    # identical to the same model built directly in code
    ref = DomeModel(geom=geo.generate_dome(3, 9.0, 0.5),
                    material_name="EN S355",
                    group_sections=["CHS 76.1 x 4.0"],
                    support_type=geo.ROLLER3,
                    loads=LoadSpec(dead_pressure=350.0, snow_pressure=1500.0,
                                   wind_speed=40.0, wind_azimuth=30.0),
                    combinations=[LoadCombination("LC2  1.2D + 1.6S",
                                                  D=1.2, S=1.6),
                                  LoadCombination("SLS", D=1.0, S=1.0,
                                                  serviceability=True)])
    ra, rb = run_analysis(m, do_buckling=False), run_analysis(ref,
                                                               do_buckling=False)
    for n in ra.order:
        _same_result(ra.cases[n], rb.cases[n], rtol=0.0)


def test_unknown_load_keys_are_ignored():
    d = json.loads(json.dumps(V1_MODEL))
    d["loads"]["some_future_field"] = 3
    d["combinations"][0]["future_factor"] = 1.0
    m = DomeModel.from_dict(d)
    assert m.loads.dead_pressure == 350.0


# ================================================================= GUI ======
def test_gui_opens_generates_and_shows_custom_model():
    import traceback
    from qtpy import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from geodome.gui.main_window import MainWindow
    errors = []
    old = sys.excepthook
    sys.excepthook = lambda t, v, tb: errors.append(
        "".join(traceback.format_exception(t, v, tb)))
    try:
        win = MainWindow()
        sup = win.geometry_panel.support
        keys = [sup.itemData(i) for i in range(sup.count())]
        assert keys == [geo.PINNED, geo.ROLLER3, geo.ELASTIC, geo.CUSTOM]
        win.geometry_panel.frequency.setValue(3)
        win.geometry_panel.generate_btn.click()
        app.processEvents()
        assert not errors, errors[-1]
        assert win.model.geom.frequency == 3
        # a v2 combination with T survives a read of the GUI panels
        win.model.combinations = [LoadCombination("DT", D=1.0, T=1.0)]
        win.load_panel.sync_from_model(win.model)
        win._read_inputs()
        assert win.model.combinations[0].T == 1.0
        # open a custom truss through File > Open (dialog patched)
        nodes, members, faces = _star_dome()
        m = DomeModel(geom=geo.custom_geometry(nodes, members, faces,
                                               title="Star dome"))
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "star.json")
            m.save(p)
            orig = QtWidgets.QFileDialog.getOpenFileName
            QtWidgets.QFileDialog.getOpenFileName = \
                staticmethod(lambda *a, **k: (p, "JSON (*.json)"))
            try:
                win.on_open()
            finally:
                QtWidgets.QFileDialog.getOpenFileName = orig
        app.processEvents()
        assert not errors, errors[-1]
        assert win.model.geom.kind == "custom"
        assert win.geometry_panel.schedule.rowCount() == \
            win.model.geom.n_strut_types()
        win.model.loads = LoadSpec(include_self_weight=True,
                                   wind_enabled=False)
        win._analysis_done(run_analysis(win.model, n_buckling_modes=2))
        app.processEvents()
        assert not errors, errors[-1]
        assert win.member_table.rowCount() == win.model.geom.n_members
        win.close()
    finally:
        sys.excepthook = old


# ============================================ critical buckling combination ==
def test_lambda_cr_is_the_minimum_over_strength_combinations():
    """run_analysis buckles every strength combination and reports the
    lowest load factor, which need not be under the combination that
    governs the member checks (default 4V dome: LC5 governs the members,
    LC2 has the lowest load factor)."""
    from geodome.buckling import linear_buckling
    from geodome.model import default_model
    from geodome.postprocess import run_analysis, summary_text

    m = default_model(4, 15.0, 0.5)
    res = run_analysis(m, n_buckling_modes=3)
    strength = [c.name for c in m.combinations if not c.serviceability]
    by_case = res.lambda_cr_by_case
    assert list(by_case) == strength
    finite = {k: v for k, v in by_case.items() if np.isfinite(v)}
    low = min(finite, key=finite.get)
    assert res.lambda_cr_case == low
    assert res.lambda_cr == finite[low]
    assert len(res.buckling.load_factors) == 3
    # each entry equals a direct eigenanalysis of that combination
    for name in strength:
        b = linear_buckling(m, res.cases[name].N, 1)
        ref = float(b.load_factors[0]) if len(b.load_factors) else np.nan
        assert (np.isnan(ref) and np.isnan(by_case[name])) or \
            abs(by_case[name] - ref) <= 1e-12 * abs(ref)
    assert res.lambda_cr_case != res.governing_case
    assert any(msg.startswith("The lowest critical load factor")
               for msg in res.messages)
    assert f"critical combination: {low}" in summary_text(res)


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
    raise SystemExit(_run_standalone())
