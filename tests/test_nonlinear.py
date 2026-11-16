"""
Tests for geodome.nonlinear (geometrically nonlinear truss analysis).

Every reference value is closed form or derived independently in this file:

* corotational element: finite-difference tangent, rigid-body objectivity,
  mirror (reflection) states, engineering strain on a straight bar;
* von Mises two-bar truss: exact engineering-strain path
  P(w) = 2 EA y (1/l - 1/L), limit point l* = (b^2 L)^(1/3), linearised
  lambda_cr = 2 EA h^3 / (L b^2) and the closed-form P-Delta path;
* spring-braced column: exact bifurcation load P_b = kL / (1 + kL/EA);
* 24-member star dome (Hangai & Kawamata geometry, docs/VALIDATION.md):
  first limit P/EA = 3.15655e-4 at w = 0.76844 cm, lambda_cr / EA = 1.52097e-3,
  exact zero-load mirror state at w = 4 cm, second peak 8.8654e-3;
* small-load agreement with the linear solver; 6V hemispherical dome under
  1.2D + 1.6S: arc-length and displacement control agree on the limit point.

Run with pytest or standalone::

    python tests/test_nonlinear.py
"""

from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from geodome import geometry as geo                                 # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import nonlinear as nl                                 # noqa: E402
from geodome import solver as slv                                   # noqa: E402
from geodome.model import DomeModel, LoadCombination, LoadSpec, \
    default_model                                                   # noqa: E402


# ============================================================== helpers ====
B, H, EA_VM = 1.0, 0.1, 1.0                    # von Mises truss (closed form derived here)


def _vm_system(EA=EA_VM, b=B, h=H):
    X = [[-b, 0, 0], [0, 0, h], [b, 0, 0]]
    M = [[0, 1], [1, 2]]
    fixed = [0, 1, 2, 4, 6, 7, 8]               # apex free in x and z only
    F = np.zeros(9)
    F[5] = -1.0                                 # unit downward apex load
    return nl.TrussSystem(X, M, [EA, EA], fixed), F


def _vm_closed(w, EA=EA_VM, b=B, h=H):
    """Independent closed form: P(w) = 2 EA y (1/l - 1/L)."""
    L = math.hypot(b, h)
    y = h - np.asarray(w, dtype=float)
    l = np.sqrt(b * b + y * y)
    return 2.0 * EA * y * (1.0 / l - 1.0 / L)


def _vm_limit(EA=EA_VM, b=B, h=H):
    L = math.hypot(b, h)
    ls = (b * b * L) ** (1.0 / 3.0)
    ys = math.sqrt(ls * ls - b * b)
    return h - ys, 2.0 * EA * ys * (1.0 / ls - 1.0 / L)


def _star_geometry():
    """24-member star dome, metres (docs/VALIDATION.md, exact hexagon)."""
    P = [[0.0, 0.0, 8.216]]
    P += [[25 * math.cos(math.radians(60 * k)),
           25 * math.sin(math.radians(60 * k)), 6.216] for k in range(6)]
    c30 = 50 * math.cos(math.radians(30))
    P += [[c30, 25, 0], [0, 50, 0], [-c30, 25, 0], [-c30, -25, 0],
          [0, -50, 0], [c30, -25, 0]]
    M = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (2, 3), (3, 4),
         (4, 5), (5, 6), (6, 7), (7, 2), (2, 8), (3, 8), (3, 9), (4, 9),
         (4, 10), (5, 10), (5, 11), (6, 11), (6, 12), (7, 12), (7, 13),
         (2, 13)]
    return np.array(P) / 100.0, np.array(M) - 1


def _star_model():
    X, M = _star_geometry()
    g = geo.custom_geometry(X, M, title="Star dome")
    m = DomeModel(geom=g, support_type=geo.PINNED)
    m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                       snow_pressure=0.0, wind_enabled=False,
                       nodal_loads=[{"node": 0, "Fz": -1.0, "case": "L"}])
    m.combinations = [LoadCombination("P", L=1.0)]
    return m


def _vm_model():
    g = geo.custom_geometry([[-B, 0, 0], [0, 0, H], [B, 0, 0]],
                            [[0, 1], [1, 2]], title="von Mises truss")
    m = DomeModel(geom=g, support_type=geo.CUSTOM,
                  custom_supports={0: [True, True, True],
                                   2: [True, True, True],
                                   1: [False, True, False]})
    m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                       snow_pressure=0.0, wind_enabled=False,
                       nodal_loads=[{"node": 1, "Fz": -1.0, "case": "L"}])
    m.combinations = [LoadCombination("P", L=1.0)]
    return m


def _worked_model():
    sys.path.insert(0, os.path.join(ROOT, "examples"))
    from worked_example import build_model
    return build_model()


# ===================================================== element level ======
def test_tangent_matches_finite_difference():
    rng = np.random.default_rng(7)
    X, M = _star_geometry()
    springs = {3 * 8 + 0: 4.0e5, 3 * 9 + 2: 7.0e5}   # springs at supports
    fixed = [d for d in range(3 * 7, 3 * 13) if d not in springs]
    EA = rng.uniform(0.5e6, 2.0e6, len(M))
    s = nl.TrussSystem(X, M, EA, fixed, springs=springs, eps0=1e-4)
    u = np.zeros(s.ndof)
    u[s.free] = 0.01 * rng.standard_normal(s.nf)      # large-ish motion
    K = s.tangent_full(u)
    h = 1e-7
    Kfd = np.zeros_like(K)
    for d in range(s.ndof):
        up, um = u.copy(), u.copy()
        up[d] += h
        um[d] -= h
        Kfd[:, d] = (s.state(up).f_full - s.state(um).f_full) / (2 * h)
    st = s.state(u)
    assert np.max(np.abs(st.N)) > 1e3                  # geometric part active
    err = np.max(np.abs(K - Kfd)) / np.max(np.abs(K))
    assert err < 1e-7, err
    assert np.allclose(K, K.T, atol=1e-9 * np.max(np.abs(K)))
    # free-free sparse tangent equals the free partition of the full one
    Kf = s.tangent(st).toarray()
    assert np.allclose(Kf, K[np.ix_(s.free, s.free)], rtol=0, atol=1e-9 *
                       np.max(np.abs(K)))


def test_tangent_at_rest_is_linear_stiffness():
    m = default_model(3, 8.0, 0.5)
    s = nl.model_system(m)
    K0 = s.tangent_full(np.zeros(s.ndof))
    K = slv.assemble_stiffness(m.geom.nodes, m.geom.members, m.EA()).toarray()
    assert np.allclose(K0, K, rtol=0, atol=1e-6 * np.abs(K).max())


def test_rigid_body_motion_is_stress_free():
    """Corotational objectivity: a large rigid rotation + translation
    produces no strain and no internal force."""
    X, M = _star_geometry()
    s = nl.TrussSystem(X, M, 1e6, [])
    a = math.radians(70.0)
    R = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0],
                  [0, 0, 1]]) @ np.array([[1, 0, 0],
                                          [0, math.cos(0.4), -math.sin(0.4)],
                                          [0, math.sin(0.4), math.cos(0.4)]])
    Xn = X @ R.T + np.array([0.3, -1.2, 0.5])
    st = s.state((Xn - X).ravel())
    assert np.max(np.abs(st.eps)) < 1e-14
    assert np.max(np.abs(st.f_full)) < 1e-6           # EA = 1e6 -> 1e-8 rel


def test_mirror_states_are_stress_free():
    """Star dome: crown reflected about the ring plane (w = 4 cm) and the
    whole dome reflected about the base (w = 16.432 cm) leave every bar at
    its original length."""
    X, M = _star_geometry()
    fixed = list(range(21, 39))
    s = nl.TrussSystem(X, M, 1e6, fixed)
    u = np.zeros(s.ndof)
    u[2] = -0.04
    assert np.max(np.abs(s.state(u).N)) < 1e-6
    u = np.zeros(s.ndof)
    u[2::3][:7] = -2.0 * X[:7, 2]
    st = s.state(u)
    assert np.max(np.abs(st.N)) < 1e-6
    assert np.max(np.abs(st.f_free)) < 1e-6


def test_engineering_strain_and_spring_are_exact():
    """Vertical bar on a vertical spring, crushed to 30 % strain: the
    engineering-strain corotational bar is exactly linear in w,
    lambda = (EA/L + k) w (a Green-strain bar would not be)."""
    L, EA, k = 2.0, 3.0e6, 5.0e5
    X = [[0, 0, 0], [0, 0, L]]
    s = nl.TrussSystem(X, [[0, 1]], EA, [0, 1, 2, 3, 4], springs={5: k})
    F = np.zeros(6)
    F[5] = -1.0
    r = nl.trace_path(s, F, "displacement", control_dof=5,
                      first_increment=0.02 * (EA / L + k) * L,
                      control_max=0.3 * L, stop_after_limit=None)
    w = -r.control
    assert w.max() > 0.3 * L
    assert np.max(np.abs(r.lam - (EA / L + k) * w)) < 1e-9 * r.lam.max()
    assert np.allclose(r.N[:, 0], -EA * w / L, rtol=0, atol=1e-6)


def test_initial_strain_state():
    """eps0 on a bar restrained by a spring: the lambda = 0 state solves
    EA(u/L - eps0) + k u = 0, i.e. u = EA eps0 / (EA/L + k)."""
    L, EA, k, e0 = 1.5, 2.0e6, 8.0e5, 1.2e-3
    s = nl.TrussSystem([[0, 0, 0], [L, 0, 0]], [[0, 1]], EA,
                       [0, 1, 2, 4, 5], springs={3: k}, eps0=e0)
    F = np.zeros(6)
    F[3] = 1.0
    r = nl.trace_path(s, F, "load", control_dof=3, first_increment=100.0,
                      max_steps=2)
    u0 = EA * e0 / (EA / L + k)
    assert abs(r.control[0] - u0) < 1e-12 * u0


# ================================================ von Mises truss ========
def test_von_mises_path_matches_closed_form_through_snap():
    s, F = _vm_system()
    ws, Ps = _vm_limit()
    r = nl.trace_path(s, F, "arc-length", control_dof=5,
                      first_increment=Ps / 20, stop_after_limit=None,
                      control_max=2.5 * H, lambda_max=20 * Ps, max_steps=500)
    w = -r.control
    err = np.max(np.abs(r.lam - _vm_closed(w))) / Ps
    assert err < 1e-9, err
    assert r.max_residual < 1e-8
    # traced through the unstable branch (lambda < 0 between w = h and 2h)
    assert r.lam.min() < -0.9 * Ps
    # ... and onto the stable inverted branch above the snap load
    assert w.max() > 2.0 * H and r.lam[-1] > Ps
    # stiffness sign: negative between the two limit points, positive outside
    k_un = np.where((w > ws + 1e-3) & (w < 2 * H - ws - 1e-3))[0]
    assert np.all(r.det_sign[k_un] == -1)
    assert r.det_sign[-1] == 1 and r.det_sign[1] == 1


def test_von_mises_limit_point_all_methods():
    s, F = _vm_system()
    ws, Ps = _vm_limit()
    for method in ("arc-length", "displacement"):
        r = nl.trace_path(s, F, method, control_dof=5,
                          first_increment=Ps / 15, stop_after_limit=3)
        assert r.limit_found and r.limit_refined and r.limit_type == "limit"
        assert abs(r.lambda_limit - Ps) / Ps < 1e-6, (method, r.lambda_limit)
        assert abs(-r.control_limit - ws) / ws < 1e-4
        assert r.limit_load_projection > 0.99           # 1-DOF: phi // F
    r = nl.trace_path(s, F, "load", control_dof=5, first_increment=Ps / 15)
    assert r.limit_found and "load control cannot pass" in r.stop_reason
    assert r.lambda_limit <= Ps * (1 + 1e-12)
    assert (Ps - r.lambda_limit) / Ps < 1e-4


def test_von_mises_model_lambda_cr_and_p_delta():
    m = _vm_model()
    EA = float(m.EA()[0])
    L = math.hypot(B, H)
    ws, Ps = _vm_limit(EA)
    r = nl.run_nonlinear(m, "P", control=(1, "uz"), stop_after_limit=3)
    lam_cr = 2 * EA * H ** 3 / (L * B ** 2)
    assert abs(r.lambda_cr - lam_cr) / lam_cr < 1e-10
    assert abs(r.lambda_limit - Ps) / Ps < 1e-6
    assert abs(r.limit_ratio - Ps / lam_cr) < 1e-6
    # linear path: initial stiffness 2 EA h^2 / L^3
    kz = 2 * EA * H ** 2 / L ** 3
    lin = r.linear
    assert np.allclose(-lin["control"], lin["lam"] / kz, rtol=1e-10)
    # P-Delta closed form: u = lambda / (kz - lambda b^2 / (h L^2))
    pd = r.p_delta
    assert len(pd["lam"]) > 10 and pd["lam"].max() < lam_cr
    ref = pd["lam"] / (kz - pd["lam"] * B ** 2 / (H * L ** 2))
    assert np.allclose(-pd["control"], ref, rtol=1e-9, atol=0)
    # dedicated P-Delta sweep up to just below lambda_cr, then singular
    s = nl.model_system(m)
    lams = np.linspace(0, 1.2 * lam_cr, 121)
    sw = nl.pdelta_path(s, r.F_ref, lams, 3 * 1 + 2)
    assert sw["lam"].max() < lam_cr <= sw["lam"].max() + lams[1] + 1e-12
    curves = r.curves()
    assert set(curves) == {"corotational", "linear", "p_delta"}


def test_energy_identity_converges_quadratically():
    s, F = _vm_system()
    _, Ps = _vm_limit()
    errs = []
    for inc in (Ps / 10, Ps / 20, Ps / 40):
        r = nl.trace_path(s, F, "displacement", control_dof=5,
                          first_increment=inc, desired_iterations=100,
                          max_arc_factor=1.0, stop_after_limit=None,
                          control_max=1.5 * H)
        errs.append(r.energy_error)
    assert errs[2] < errs[1] < errs[0]
    assert 3.0 < errs[0] / errs[1] < 5.0 and 3.0 < errs[1] / errs[2] < 5.0
    assert errs[2] < 1e-3


# ===================================================== bifurcation =======
def _tent(b=1.0, h=4.0, EA=1e7, P=1e6):
    """Three-bar tent: apex (0, 0, h) on three pinned legs at radius b,
    120 deg apart; load P down at the apex.  Only the apex is free."""
    X = [[0.0, 0.0, h]] + [[b * math.cos(math.radians(a)),
                            b * math.sin(math.radians(a)), 0.0]
                           for a in (0.0, 120.0, 240.0)]
    M = [[1, 0], [2, 0], [3, 0]]
    F = np.zeros(12)
    F[2] = -P
    return nl.TrussSystem(X, M, [EA] * 3, list(range(3, 12))), F


def _tent_closed(b=1.0, h=4.0, EA=1e7, P=1e6):
    """Independent closed form of the tent on its symmetric path.

    Apex drop w, y = h - w, l = sqrt(b^2 + y^2), L = sqrt(b^2 + h^2),
    N = EA (l - L)/L, vertical equilibrium lambda = -3 N y / (l P).
    Lateral stiffness (double, x and y):
        K_xx = (EA/L)(3/2)(b/l)^2 + (N/l)(3 - (3/2)(b/l)^2)
    Vertical stiffness K_zz = 3 [(EA/L)(y/l)^2 + (N/l)(b/l)^2].
    Returns (lambda at K_xx = 0, lambda at K_zz = 0)."""
    from scipy.optimize import brentq
    L = math.hypot(b, h)

    def parts(w):
        y = h - w
        l = math.hypot(b, y)
        N = EA * (l - L) / L
        lam = -3.0 * N * y / (l * P)
        kxx = EA / L * 1.5 * (b / l) ** 2 + N / l * (3.0 - 1.5 * (b / l) ** 2)
        kzz = 3.0 * (EA / L * (y / l) ** 2 + N / l * (b / l) ** 2)
        return lam, kxx, kzz

    w_b = brentq(lambda w: parts(w)[1], 1e-9, 0.5 * h, xtol=1e-15, rtol=1e-15)
    w_l = brentq(lambda w: parts(w)[2], 1e-9, h - 1e-9, xtol=1e-15,
                 rtol=1e-15)
    return parts(w_b)[0], parts(w_l)[0]


def test_tent_double_bifurcation_bracket():
    """Review F8: next to the double bifurcation of the three-bar tent the
    bisection must not reject converged stable points (closed form
    lambda_b = 0.9104292470; the old bracket missed it by 2.7e-7)."""
    lam_b, lam_l = _tent_closed()
    assert abs(lam_b - 0.9104292470) < 1e-9
    s, F = _tent()
    r = nl.trace_path(s, F, "arc-length", control_dof=2,
                      first_increment=0.05, lambda_max=2.0 * lam_l,
                      stop_after_limit=3)
    bif = [cp for cp in r.critical_points if cp["type"] == "bifurcation"]
    assert bif and bif[0]["refined"], r.critical_points
    lo, hi = bif[0]["lambda_bracket"]
    assert hi - lo <= 1e-9 * lam_b
    assert abs(bif[0]["lambda"] - lam_b) < 2e-9 * lam_b, \
        (bif[0]["lambda"], lam_b)
    assert abs(r.lambda_first_critical - lam_b) < 2e-9 * lam_b
    assert not any("det K_T only" in t for t in r.messages)


def test_tent_determinant_sign_path_warns():
    """Review F9: above eig_limit the stability check is the determinant
    sign, which cannot see a double crossing; the run must say so."""
    lam_b, lam_l = _tent_closed()
    s, F = _tent()
    r = nl.trace_path(s, F, "arc-length", control_dof=2,
                      first_increment=0.05, lambda_max=2.0 * lam_l,
                      stop_after_limit=3, eig_limit=0)
    assert any("det K_T only" in t for t in r.messages), r.messages
    # the double bifurcation is missed: the limit point comes first
    assert r.limit_found and r.limit_type == "limit"
    assert abs(r.lambda_limit - lam_l) < 1e-7 * lam_l, (r.lambda_limit, lam_l)
    assert r.lambda_first_critical > 10.0 * lam_b


def test_first_arc_step_never_goes_to_a_negative_load():
    """Review GUI-06: a first arc-length step that is far too long may
    converge on the branch below lambda = 0; it must be cut back, and a
    path that dips below zero before its maximum is not a limit point."""
    # The case study at span/200 with the dialog defaults: the + run used
    # to start at lambda = -12.3 and report a spurious 57.1 (> lambda_cr).
    from geodome import benchmarks as bm
    m = bm.make_preset("case_study_6v")
    span = m.geom.base_diameter
    rows, runs = nl.imperfection_study(m, amplitudes=[span / 200.0],
                                       return_results=True)
    perfect, row = rows
    for (a, sg), r in runs.items():
        k = r.limit_step if r.limit_found else r.n_points - 1
        assert np.all(r.lam[:k + 1] >= 0.0), (a, sg, r.lam[:4])
        if r.limit_found:
            assert r.lambda_limit < r.lambda_cr, (a, sg, r.lambda_limit)
    assert row["found"]
    assert row["lambda_limit"] < perfect["lambda_limit"]
    # 150 mm along the crown dimple leaves the crown hub almost flat with
    # its neighbours: the first limit is close to zero load (the 0.288 found
    # before the capped first step was the later maximum of the inverted
    # crown, see test_imperfect_crown_snap_near_zero_load_is_not_jumped_over)
    assert row["ratio"] < 1e-3, row["ratio"]
    assert row["sign"] == 1.0
    # 150 mm leaves almost none of the crown's 0.198 m rise: the row says so
    assert row["hub_rise_left_m"] < 0.01 and "99 %" in row["note"]


def test_spring_braced_column_bifurcation():
    """Bar on a fixed pin, top braced by a lateral spring k: the straight
    path loses stability at P_b = kL / (1 + kL/EA) (lateral stiffness
    k + N/l = 0 with N = -P and l = L(1 - P/EA))."""
    L, EA, k = 3.0, 1.0e7, 2.0e5
    s = nl.TrussSystem([[0, 0, 0], [0, 0, L]], [[0, 1]], EA, [0, 1, 2, 4],
                       springs={3: k})
    F = np.zeros(6)
    F[5] = -1.0
    Pb = k * L / (1 + k * L / EA)
    r = nl.trace_path(s, F, "arc-length", control_dof=5,
                      first_increment=Pb / 7, lambda_max=2.0 * Pb)
    bif = [cp for cp in r.critical_points if cp["type"] == "bifurcation"]
    assert len(bif) == 1 and bif[0]["refined"]
    assert abs(bif[0]["lambda"] - Pb) / Pb < 1e-8, bif[0]["lambda"] / Pb
    assert bif[0]["load_projection"] < 1e-10           # mode orthogonal to F
    assert not r.limit_found
    assert abs(r.lambda_first_critical - Pb) / Pb < 1e-8


# ======================================================== star dome ======
def test_star_dome_limit_point_and_linear_buckling():
    m = _star_model()
    EA = float(m.EA()[0])
    for method in ("arc-length", "displacement"):
        r = nl.run_nonlinear(m, "P", method=method, stop_after_limit=3,
                             first_increment=4e-5 * EA)
        assert r.limit_found and r.limit_type == "limit" and r.limit_refined
        assert abs(r.lambda_limit / EA - 3.15655e-4) < 1e-9, \
            r.lambda_limit / EA
        assert abs(-r.control_limit * 100 - 0.76844) < 2e-4
        assert abs(r.lambda_cr / EA - 1.52097e-3) < 1e-8
        assert abs(r.limit_ratio - 0.2075) < 1e-3


def test_star_dome_symmetric_path_landmarks():
    """Displacement control on the crown follows the symmetric path: zero
    load in the exact mirror state w = 4 cm, second maximum 8.8654e-3 EA
    near 10.54 cm, failure where the crown displacement snaps back
    (near 12.96 cm, docs/VALIDATION.md)."""
    m = _star_model()
    EA = float(m.EA()[0])
    r = nl.run_nonlinear(m, "P", method="displacement", stop_after_limit=None,
                         first_increment=2e-5 * EA, control_max=0.2,
                         lambda_max=0.05 * EA, max_steps=3000, compare=False,
                         member_checks=False)
    w = -r.control
    lam = r.lam / EA
    k = int(np.searchsorted(w, 0.04))
    sl = slice(k - 2, k + 2)                       # local cubic through 4 points
    lam4 = np.polyval(np.polyfit(w[sl] - 0.04, lam[sl], 3), 0.0)
    assert abs(lam4) < 1e-6 * lam.max(), lam4
    i2 = int(np.argmax(lam))
    assert abs(lam[i2] - 8.8654e-3) / 8.8654e-3 < 1e-3
    assert abs(w[i2] * 100 - 10.537) < 0.15
    assert 0.1290 < w.max() < 0.1300 and r.stop_reason.startswith("step")


# ================================================ consistency checks ====
def test_small_load_agrees_with_linear_solver():
    for support in (geo.PINNED, geo.ELASTIC):
        m = default_model(3, 8.0, 0.5)
        m.support_type = support
        m.spring_stiffness = (2e7, 3e7, 4e7)
        combo = m.combinations[1]
        F = ld.combine(ld.basic_load_vectors(m), combo)
        lin = slv.solve_model(m, F)
        s = nl.model_system(m)
        eps = 1e-6
        r = nl.trace_path(s, F, "load", control_dof=3 * m.geom.crown_node + 2,
                          first_increment=eps, max_steps=1)
        u_nl = r.u[1] / eps
        rel = np.linalg.norm(u_nl - lin.u) / np.linalg.norm(lin.u)
        assert rel < 1e-5, (support, rel)
        assert np.allclose(r.N[1] / eps, lin.N, rtol=0,
                           atol=1e-4 * np.abs(lin.N).max())


def test_self_straining_excluded_with_message():
    m = _vm_model()
    m.loads.temperature_change = 30.0
    m.combinations = [LoadCombination("PT", L=1.0, T=1.0)]
    r = nl.run_nonlinear(m, "PT", control=(1, "uz"), stop_after_limit=2,
                         compare=False)
    assert any("not included" in t for t in r.messages)
    _, Ps = _vm_limit(float(m.EA()[0]))
    assert abs(r.lambda_limit - Ps) / Ps < 1e-6


def test_result_is_json_serialisable():
    m = _vm_model()
    r = nl.run_nonlinear(m, "P", control=(1, "uz"), stop_after_limit=2)
    d = r.to_dict()
    txt = json.dumps(d, allow_nan=False)
    back = json.loads(txt)
    assert back["lambda_limit"] == r.lambda_limit
    assert "u" not in back
    assert len(json.loads(json.dumps(r.to_dict(snapshots=True),
                                     allow_nan=False))["u"]) == r.n_points
    assert any(line.startswith("First load maximum")
               for line in r.summary_lines())


# =================================================== 6V case-study dome ==
def test_6v_dome_limit_point_arc_length_vs_displacement():
    """The 30 m 6V hemisphere of the worked example under 1.2D + 1.6S is
    traced to its first limit point in well under 30 s; two independent
    path-following methods agree on it, and the tangent is singular there."""
    m = _worked_model()
    t0 = time.perf_counter()
    ra = nl.run_nonlinear(m, "LC2")
    ta = time.perf_counter() - t0
    assert ta < (120.0 if os.environ.get("CI") else 30.0), ta
    assert ra.limit_found and ra.limit_refined and ra.limit_type == "limit"
    assert ra.max_residual < 1e-8
    assert 0.0 < ra.limit_ratio < 1.0
    rd = nl.run_nonlinear(m, "LC2", method="displacement", compare=False,
                          member_checks=False)
    assert abs(ra.lambda_limit - rd.lambda_limit) / ra.lambda_limit < 1e-8
    # stability count goes 0 -> 1 across the limit (between the samples
    # around the peak), and K_T(u_limit) is singular
    k = ra.limit_step
    assert ra.n_negative[k - 1] == 0 and ra.n_negative[k + 1] == 1
    s = nl.model_system(m)
    ev = np.linalg.eigvalsh(s.tangent(s.state(ra.u_limit)).toarray())
    assert abs(ev[0]) < 1e-4 * ev[1]
    # first strut reaches its Euler load well before the global limit
    assert ra.lambda_member_euler < ra.lambda_limit
    # regression value (cross-checked with OpenSees corotTruss, 3e-12)
    assert abs(ra.lambda_limit - 34.2369577) / 34.2369577 < 1e-6


def test_imperfection_study_6v_dome():
    m = _worked_model()
    span = m.geom.base_diameter
    rows = nl.imperfection_study(m, "LC2",
                                 amplitudes=[span / 1000, span / 500])
    assert len(rows) == 3 and rows[0]["amplitude_m"] == 0.0
    assert all(r["found"] for r in rows)
    lam = [r["lambda_limit"] for r in rows]
    assert lam[0] > lam[1] > lam[2] > 0.0            # imperfection sensitive
    for r in rows:
        assert abs(r["ratio"] - r["lambda_limit"] / r["lambda_cr"]) < 1e-12
        assert r["lambda_critical"] <= r["lambda_limit"] + 1e-9
    assert abs(rows[1]["amplitude_over_span"] - 1e-3) < 1e-12
    assert rows[1]["lambda_limit"] == min(rows[1]["lambda_plus"],
                                          rows[1]["lambda_minus"])
    # the crown is the hub of the dimple; span/1000 and span/500 take 20 %
    # and 40 % of its 0.198 m rise, below the 50 % note threshold
    c = m.geom.crown_node
    assert all(r["hub"] == c for r in rows)
    assert all(abs(r["hub_rise_m"] - 0.1982) < 5e-4 for r in rows)
    left = [r["hub_rise_left_m"] for r in rows]
    assert left[0] == rows[0]["hub_rise_m"] and left[0] > left[1] > left[2]
    assert all(r["note"] == "" for r in rows)


def test_local_rise_and_hub_note():
    """local_rise of the 6V crown equals the closed geometry (crown height
    above the plane of its five neighbours, which is horizontal by
    symmetry); the note appears only when more than half of the rise is
    gone, and says 'inverted' when the hub passes its neighbours."""
    import types
    from geodome import benchmarks as bm
    m = bm.make_preset("case_study_6v")
    g = m.geom
    c = g.crown_node
    nb = sorted({int(j) if int(i) == c else int(i) for i, j in g.members
                 if c in (int(i), int(j))})
    assert len(nb) == 5
    h, n = nl.local_rise(g, c)
    assert abs(h - (g.nodes[c, 2] - np.mean(g.nodes[nb, 2]))) < 1e-12
    assert abs(h - 0.1982) < 5e-4 and n[2] > 0.999999
    assert nl._hub_note(c, h, 0.6 * h) == ""
    note = nl._hub_note(c, h, 0.3 * h)
    assert "70 %" in note and "inverted" not in note
    off = np.zeros_like(g.nodes)
    off[c, 2] = -1.5 * h                       # push the crown through
    h2, _ = nl.local_rise(types.SimpleNamespace(nodes=g.nodes + off,
                                                members=g.members), c)
    assert abs(h2 + 0.5 * h) < 1e-12
    assert "inverted" in nl._hub_note(c, h, h2)
    assert np.isnan(nl.local_rise(types.SimpleNamespace(
        nodes=np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0]]),
        members=np.array([[0, 1], [0, 2]])), 0)[0])


def test_buckling_mode_sign_is_platform_independent():
    """Eigenvector signs differ between LAPACK builds (Windows CI with
    OpenBLAS gave the opposite crown dimple to MKL).  linear_buckling scales
    every mode to a largest component of +1, and linear_buckling_reference
    turns each mode to do positive work with the load, so a positive
    imperfection always points the way the load pushes."""
    from geodome import benchmarks as bm
    from geodome import buckling as bk
    m = bm.make_preset("case_study_6v")
    F = nl.reference_load(m, nl._find_combination(m, "LC2"))
    r = slv.solve_model(m, F, name="ref", T_factor=0.0)
    b = bk.linear_buckling(m, r.N, n_modes=4)
    for k in range(b.modes.shape[1]):
        v = b.modes[:, k]                  # (a symmetric mode may also
        assert abs(np.max(v) - 1.0) < 1e-12           # reach -1 elsewhere)
        assert abs(np.max(np.abs(v)) - 1.0) < 1e-12
    lam, phi, b2 = nl.linear_buckling_reference(m, F)
    for k in range(b2.modes.shape[1]):
        assert float(b2.modes[:, k] @ F) >= -1e-9 * np.linalg.norm(F)
    c = m.geom.crown_node
    assert phi[3 * c + 2] == -1.0          # the crown dimple points down
    # linear_buckling alone points the crown up (largest component +1);
    # the reference turns it to point with the load
    assert b.modes[3 * c + 2, 0] == 1.0
    assert np.allclose(phi, -b.modes[:, 0], atol=1e-9)


def test_imperfect_crown_snap_near_zero_load_is_not_jumped_over():
    """The 6V case-study trial dome, mode-1 imperfection (a crown dimple).
    At span/200 = 150 mm the crown hub is almost flat with its five
    neighbours (local rise 0.198 m), so its first limit load is close to
    zero.  A nominal first arc-length step (about 0.3 m) jumped over that
    snap onto the inverted-crown branch and reported its later maximum
    (5.85 on the sized dome).  The default capped first step must find the
    same first limit as an explicit 1 mm step."""
    from geodome import benchmarks as bm
    m = bm.make_preset("case_study_6v")
    combo = "LC2"
    F = nl.reference_load(m, nl._find_combination(m, combo))
    lam_cr, _, b = nl.linear_buckling_reference(m, F)
    phi = b.modes[:, 0] / np.max(np.abs(b.modes[:, 0]))
    off = (m.geom.base_diameter / 200.0 * phi).reshape(-1, 3)
    kw = dict(node_offset=off, load_vector=F, lambda_cr=lam_cr,
              compare=False, member_checks=False)
    r_def = nl.run_nonlinear(m, combo, **kw)
    r_small = nl.run_nonlinear(m, combo, arc_length=1e-3, max_steps=2000,
                               **kw)
    assert r_def.limit_found and r_small.limit_found
    assert r_def.lambda_limit < 1e-3 * lam_cr, r_def.lambda_limit
    assert abs(r_def.lambda_limit - r_small.lambda_limit) <= \
        1e-6 * lam_cr
    # the first step is the capped one, and the run still finishes quickly
    assert np.linalg.norm(r_def.u[1] - r_def.u[0]) <= \
        2e-5 * m.geom.length_scale * (1 + 1e-9)
    assert len(r_def.lam) < 300
    # a perfect run keeps the nominal first step (no change for it)
    r0 = nl.run_nonlinear(m, combo, compare=False, member_checks=False,
                          lambda_cr=lam_cr, load_vector=F)
    assert np.linalg.norm(r0.u[1] - r0.u[0]) > 1e-2


# =============================================================== driver =====
def _run_standalone():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    npass = 0
    for name, fn in fns:
        t0 = time.perf_counter()
        try:
            fn()
        except Exception as exc:
            print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        else:
            npass += 1
            print(f"  pass  {name}  ({time.perf_counter() - t0:.2f} s)")
    print(f"\n{npass} of {len(fns)} tests passed.")
    return 0 if npass == len(fns) else 1


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    raise SystemExit(_run_standalone())
