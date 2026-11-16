"""
Tests for the modal analysis module (spec task B, ``geodome/modal.py``).

Every check is against a closed-form result or an identity:

* single bar + end mass, omega^2 = EA / (L m), lumped and consistent bar
  mass, and an inclined bar whose two transverse modes are mechanisms;
* two-DOF spring-mass chain: closed-form eigenvalues, mode shapes,
  participation factors and effective masses;
* exact static condensation of massless nodes (springs in series);
* fixed-free rod: the exact discrete spectra of the lumped and consistent
  mass models, the bounds and O(h^2) convergence to the continuous rod, and
  the continuum effective modal mass 8 / ((2j - 1)^2 pi^2);
* validation Model A (1V cap): crown frequencies from the hand formulas;
* mass = dead load / g; M- and K-orthogonality against an independently
  assembled K and M and an independent numpy eigensolve;
* sum of the effective masses over all modes = free mass in x, y and z for
  every support type, lumped and consistent;
* modal superposition with all modes reproduces the static solution;
* sparse shift-invert solver = dense solver; repeated-frequency pairs are
  aligned with x and y without changing their total effective mass;
* pre-stress: beaded string omega^2 = 2 N / (a m), and the lowest frequency
  vanishes at lambda_cr of the linear buckling analysis (buckling mode);
* the 6V 15 m case-study dome, input validation, CSV export, headless import.

Run with pytest::

    pytest -q tests/test_modal.py

or standalone, without pytest installed::

    python tests/test_modal.py
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geodome import geometry as geo                                 # noqa: E402
from geodome import buckling as bk                                  # noqa: E402
from geodome import loads as ld                                     # noqa: E402
from geodome import modal as md                                     # noqa: E402
from geodome.model import DomeModel, LoadCombination, LoadSpec      # noqa: E402
from geodome.sections import SECTIONS                               # noqa: E402
from geodome.solver import solve_model                              # noqa: E402
from geodome.units import G_ACCEL                                   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _check_identities(r, tol=1e-10):
    assert r.orthogonality_error < tol, r.orthogonality_error
    assert r.stiffness_error < tol, r.stiffness_error
    assert r.residual < 1e-12, r.residual


def _line_x(n_el, h, free_axis=0):
    """Nodes 0..n_el on the x axis; node 0 fixed, others free in one axis."""
    nodes = np.array([[k * h, 0.0, 0.0] for k in range(n_el + 1)])
    members = np.array([[k, k + 1] for k in range(n_el)])
    fixed = [0, 1, 2] + [3 * k + a for k in range(1, n_el + 1)
                         for a in range(3) if a != free_axis]
    return nodes, members, fixed


def _dome_model(freq, R=10.0, frac=0.5, section="CHS 88.9 x 4.0",
                support=geo.PINNED, dead=300.0, snow=1200.0):
    g = geo.generate_dome(freq, R, frac)
    m = DomeModel(geom=g, group_sections=[section] * g.n_strut_types(),
                  support_type=support)
    m.loads = LoadSpec(include_self_weight=True, dead_pressure=dead,
                       snow_pressure=snow, wind_enabled=False)
    return m


# ===================================================== single bar ===========
def test_single_bar_end_mass_closed_form():
    """omega^2 = EA / (L m): one axial DOF, massless bar."""
    EA, L, m = 2.0e8, 2.5, 400.0
    nodes, members, fixed = _line_x(1, L)
    r = md.solve_modes(nodes, members, [EA], [0.0, m], fixed)
    assert r.n_modes == 1 and r.n_free == 1
    assert _rel(r.omega2[0], EA / (L * m)) < 1e-13
    assert _rel(r.frequencies[0], math.sqrt(EA / (L * m)) / (2 * math.pi)) \
        < 1e-13
    assert _rel(r.periods[0], 2 * math.pi / math.sqrt(EA / (L * m))) < 1e-13
    # mass-normalised: phi = 1/sqrt(m), Gamma_x = sqrt(m), M*_x = m
    assert _rel(r.modes[3, 0], 1.0 / math.sqrt(m)) < 1e-13
    assert _rel(r.participation[0, 0], math.sqrt(m)) < 1e-13
    assert _rel(r.effective_mass[0, 0], m) < 1e-13
    assert _rel(r.effective_mass_ratio[0, 0], 1.0) < 1e-13
    assert np.allclose(r.free_mass, [m, 0.0, 0.0])
    assert np.allclose(r.total_mass, [m, m, m])
    assert not r.warnings
    _check_identities(r)


def test_inclined_bar_axial_mode_and_two_mechanisms():
    """A free node on one bar: one axial mode, two zero-frequency modes."""
    EA, m = 3.0e8, 250.0
    c = np.array([1.0, 2.0, 2.0]) / 3.0              # L = 3 m
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 2.0]])
    r = md.solve_modes(nodes, [[0, 1]], [EA], [0.0, m], [0, 1, 2],
                       n_modes=None)
    assert r.n_modes == 3
    assert r.omega2[0] == r.omega2[0] and r.frequencies[0] == 0.0
    assert r.frequencies[1] == 0.0 and math.isinf(r.periods[1])
    assert _rel(r.omega2[2], EA / (3.0 * m)) < 1e-13
    phi = r.mode(2)[1]
    assert abs(abs(phi @ c) * math.sqrt(m) - 1.0) < 1e-12   # along the bar
    assert _rel(r.effective_mass[2, 0], m * c[0] ** 2) < 1e-12
    assert _rel(r.effective_mass[2, 2], m * c[2] ** 2) < 1e-12
    assert np.allclose(r.effective_mass.sum(axis=0), m, rtol=1e-12)
    assert any("mechanism" in w for w in r.warnings), r.warnings
    _check_identities(r)


def test_bar_with_self_mass_through_the_model():
    """Custom geometry + custom supports + nodal mass, lumped and consistent.

    Lumped:     omega^2 = EA / (L (m_end + rho A L / 2))
    Consistent: omega^2 = EA / (L (m_end + rho A L / 3))   (2/6 at the end)
    """
    L, m_end = 3.0, 120.0
    g = geo.custom_geometry([[0.0, 0, 0], [L, 0, 0]], [[0, 1]],
                            base_nodes=[0], title="bar")
    model = DomeModel(geom=g, group_sections=["CHS 88.9 x 4.0"],
                      support_type=geo.CUSTOM,
                      custom_supports={0: [True, True, True],
                                       1: [False, True, True]})
    A = SECTIONS["CHS 88.9 x 4.0"].A
    E, rho = model.material.E, model.material.rho
    for mtype, frac in (("lumped", 0.5), ("consistent", 1.0 / 3.0)):
        r = md.modal_analysis(model, nodal_masses={1: m_end},
                              mass_type=mtype)
        exact = E * A / (L * (m_end + frac * rho * A * L))
        assert r.n_modes == 1
        assert _rel(r.omega2[0], exact) < 1e-12, (mtype, r.omega2[0], exact)
        _check_identities(r)
    # lumped totals: struts rho A L + the end mass
    assert _rel(r.total_mass[0], rho * A * L + m_end) < 1e-12
    assert _rel(r.mass_parts["added"], m_end) < 1e-12


# ======================================================== two DOF ===========
def test_two_dof_chain_closed_form():
    """Fixed - k1 - m1 - k2 - m2 (x only): eigenvalues, shapes, Gamma, M*."""
    EA1, a, EA2, b = 2.0e8, 2.0, 1.0e8, 1.5
    m1, m2 = 300.0, 150.0
    k1, k2 = EA1 / a, EA2 / b
    nodes = np.array([[0.0, 0, 0], [a, 0, 0], [a + b, 0, 0]])
    members = [[0, 1], [1, 2]]
    fixed = [0, 1, 2, 4, 5, 7, 8]
    r = md.solve_modes(nodes, members, [EA1, EA2], [0.0, m1, m2], fixed,
                       n_modes=None)
    # det(K - lam M) = m1 m2 lam^2 - (m1 k2 + m2 (k1 + k2)) lam + k1 k2 = 0
    B = m1 * k2 + m2 * (k1 + k2)
    disc = math.sqrt(B * B - 4.0 * m1 * m2 * k1 * k2)
    lam = [(B - disc) / (2 * m1 * m2), (B + disc) / (2 * m1 * m2)]
    for n in range(2):
        assert _rel(r.omega2[n], lam[n]) < 1e-12
        ratio = (k1 + k2 - lam[n] * m1) / k2          # u2 / u1
        u1 = 1.0 / math.sqrt(m1 + m2 * ratio ** 2)
        u2 = ratio * u1
        sgn = np.sign(r.modes[3, n]) * np.sign(u1)
        assert _rel(sgn * r.modes[3, n], u1) < 1e-12
        assert _rel(sgn * r.modes[6, n], u2) < 1e-12
        Gam = m1 * u1 + m2 * u2
        assert _rel(abs(r.participation[n, 0]), abs(Gam)) < 1e-12
        assert _rel(r.effective_mass[n, 0], Gam ** 2) < 1e-12
    assert _rel(r.effective_mass[:, 0].sum(), m1 + m2) < 1e-13
    assert r.participation[0, 0] > 0.0                # sign convention
    _check_identities(r)


def test_massless_node_is_condensed_exactly():
    """Fixed - k1 - (no mass) - k2 - m: omega^2 = k1 k2 / ((k1 + k2) m)."""
    EA1, a, EA2, b, m = 2.0e8, 2.0, 1.0e8, 1.5, 500.0
    k1, k2 = EA1 / a, EA2 / b
    nodes = np.array([[0.0, 0, 0], [a, 0, 0], [a + b, 0, 0]])
    r = md.solve_modes(nodes, [[0, 1], [1, 2]], [EA1, EA2], [0.0, 0.0, m],
                       [0, 1, 2, 4, 5, 7, 8])
    assert r.n_free == 2 and r.n_dynamic == 1
    assert _rel(r.omega2[0], k1 * k2 / ((k1 + k2) * m)) < 1e-13
    # the massless node follows the static relation u1 = k2/(k1+k2) u2
    assert _rel(r.modes[3, 0], k2 / (k1 + k2) * r.modes[6, 0]) < 1e-13
    assert "condensed" in r.method
    _check_identities(r)


# ============================================================ rod ===========
def _rod(N, Lr=10.0, EA=2.1e8, mu=7.85, mtype="lumped", **kw):
    nodes, members, fixed = _line_x(N, Lr / N)
    return md.solve_modes(nodes, members, np.full(N, EA), None, fixed,
                          line_mass=mu, mass_type=mtype, **kw)


def test_rod_lumped_and_consistent_exact_discrete_spectra():
    """Fixed-free rod of N bars (h = L/N, c^2 = EA/mu), all N modes:

    lumped      omega_j = (2c/h) sin((2j-1) pi / (4N))
    consistent  omega_j^2 = (6c^2/h^2)(1 - cos t)/(2 + cos t),
                t = (2j-1) pi / (2N)
    (the free end is the mirror plane of a fixed-fixed rod of 2N bars).
    """
    Lr, EA, mu = 10.0, 2.1e8, 7.85
    c = math.sqrt(EA / mu)
    for N in (1, 2, 5, 12):
        h = Lr / N
        j = np.arange(1, N + 1)
        r = _rod(N, Lr, EA, mu, "lumped", n_modes=None)
        exact = (2 * c / h) * np.sin((2 * j - 1) * math.pi / (4 * N))
        assert np.allclose(r.omega, exact, rtol=1e-12, atol=0), N
        _check_identities(r)
        r = _rod(N, Lr, EA, mu, "consistent", n_modes=None)
        t = (2 * j - 1) * math.pi / (2 * N)
        exact2 = 6 * c * c / h ** 2 * (1 - np.cos(t)) / (2 + np.cos(t))
        assert np.allclose(r.omega2, exact2, rtol=1e-12, atol=0), N
        _check_identities(r)


def test_rod_converges_to_continuum_with_bounds():
    """omega_1 = pi c / (2 L): lumped below, consistent above, error ~ h^2;
    first-mode effective mass -> 8/pi^2 of the rod mass."""
    Lr, EA, mu = 10.0, 2.1e8, 7.85
    w_exact = math.pi * math.sqrt(EA / mu) / (2 * Lr)
    err = {}
    for mtype in ("lumped", "consistent"):
        e = []
        for N in (8, 16, 32):
            w = _rod(N, Lr, EA, mu, mtype, n_modes=3).omega[0]
            if mtype == "lumped":
                assert w < w_exact
            else:
                assert w > w_exact
            e.append(abs(w - w_exact) / w_exact)
        err[mtype] = e
        for k in range(2):
            assert 3.8 < e[k] / e[k + 1] < 4.2, (mtype, e)   # O(h^2)
    # continuum effective mass of mode j: 8 / ((2j - 1)^2 pi^2) of mu L
    r = _rod(200, Lr, EA, mu, "lumped", n_modes=3)
    for jj in range(3):
        cont = 8.0 / ((2 * jj + 1) ** 2 * math.pi ** 2)
        assert abs(r.effective_mass[jj, 0] / (mu * Lr) - cont) < 2e-3
    # the sparse shift-invert path gives the same numbers
    rs = _rod(200, Lr, EA, mu, "lumped", n_modes=3, dense_limit=0)
    assert "sparse" in rs.method
    assert np.allclose(rs.omega2, r.omega2, rtol=1e-10, atol=0)
    _check_identities(rs)


def test_chain_with_massless_nodes_dense_and_sparse():
    """2N equal bars, mass m on every second node (m/2 at the free end):
    a lumped chain of N masses on springs k/2, so
    omega_j = 2 sqrt(k/(2m)) sin((2j-1) pi / (4N)).  Both solvers."""
    N, h, EA, m = 40, 0.5, 1.0e8, 80.0
    k = EA / h
    nodes, members, fixed = _line_x(2 * N, h)
    masses = np.zeros(2 * N + 1)
    masses[2::2] = m
    masses[-1] = 0.5 * m
    j = np.arange(1, 7)
    exact = 2 * math.sqrt(k / (2 * m)) * np.sin((2 * j - 1) * math.pi
                                                / (4 * N))
    for dl in (2400, 0):
        r = md.solve_modes(nodes, members, np.full(2 * N, EA), masses,
                           fixed, n_modes=6, dense_limit=dl)
        assert r.n_dynamic == N and r.n_free == 2 * N
        assert np.allclose(r.omega, exact, rtol=1e-11, atol=0), r.method
        _check_identities(r)


# ======================================================== Model A ===========
def test_model_A_crown_frequencies_hand_formulas():
    """Validation Model A (1V cap, crown free): from hand formulas,

        k = EA/L, s = a/L, t = H/L,  m = 5 rho A L / 2 + q (5/3)(sqrt3/4 L^2)/g
        omega_z^2 = 5 k t^2 / m,   omega_x^2 = omega_y^2 = (5/2) k s^2 / m
    """
    R, q = 5.0, 300.0
    frac = 0.5 * (1.0 - 1.0 / math.sqrt(5.0))
    g = geo.generate_dome(1, R, frac)
    model = DomeModel(geom=g, group_sections=["CHS 88.9 x 4.0"])
    model.loads = LoadSpec(include_self_weight=False, dead_pressure=q,
                           snow_pressure=0.0, wind_enabled=False)
    sec = SECTIONS["CHS 88.9 x 4.0"]
    E, rho = model.material.E, model.material.rho
    a = 2 * R / math.sqrt(5)
    H = R - R / math.sqrt(5)
    L = R * math.sqrt(2 - 2 / math.sqrt(5))
    k = E * sec.A / L
    s, t = a / L, H / L
    for include_dead in (False, True):
        m = 5 * rho * sec.A * L / 2
        if include_dead:
            m += q * 5.0 / 3.0 * (math.sqrt(3) / 4 * L * L) / G_ACCEL
        r = md.modal_analysis(model, include_dead=include_dead)
        assert r.n_modes == 3 and r.n_free == 3
        wz2, wx2 = 5 * k * t * t / m, 2.5 * k * s * s / m
        assert wz2 < wx2                       # vertical mode is the lowest
        assert _rel(r.omega2[0], wz2) < 1e-11
        assert _rel(r.omega2[1], wx2) < 1e-11
        assert _rel(r.omega2[2], wx2) < 1e-11
        # every mode is a pure crown translation carrying all of the mass
        assert np.allclose(r.free_mass, m, rtol=1e-12)
        assert np.allclose(np.diag(r.effective_mass_ratio[:, [2, 0, 1]]),
                           1.0, rtol=0, atol=1e-12)
        crown = g.crown_node
        assert _rel(r.mode(0)[crown, 2], 1 / math.sqrt(m)) < 1e-12
        assert _rel(r.mode(1)[crown, 0], 1 / math.sqrt(m)) < 1e-12  # pure x
        assert _rel(r.mode(2)[crown, 1], 1 / math.sqrt(m)) < 1e-12  # pure y
        assert [r.dominant_direction(i) for i in range(3)] == ["z", "x", "y"]
        _check_identities(r)


# ================================================================ mass ======
def test_mass_equals_dead_load_over_g():
    model = _dome_model(3, R=8.0, dead=400.0, snow=1000.0)
    crown = model.geom.crown_node
    model.loads.nodal_loads = [
        {"node": crown, "Fz": -5000.0, "case": "D"},
        {"node": 7, "Fx": 100.0, "Fz": -2000.0, "case": "D"},
        {"node": 3, "Fz": -1000.0, "case": "L"},          # live: no mass
    ]
    mm = md.lumped_masses(model)
    D = ld.basic_load_vectors(model)["D"]
    assert _rel(mm.total * G_ACCEL, -D[2::3].sum()) < 1e-12
    assert np.allclose(mm.nodal * G_ACCEL, -D[2::3], rtol=1e-12, atol=1e-9)
    assert _rel(mm.totals()["struts"], model.steel_mass()) < 1e-12
    assert _rel(mm.totals()["cladding"],
                400.0 * model.geom.surface_area / G_ACCEL) < 1e-12
    assert _rel(mm.totals()["dead nodal loads"], 7000.0 / G_ACCEL) < 1e-12
    assert any("Horizontal" in n for n in mm.notes)
    parts = sum(v for v in mm.parts.values())
    assert np.allclose(parts, mm.nodal, rtol=1e-14, atol=0)
    # snow fraction on the plan area; strut mass even without self weight
    mm2 = md.lumped_masses(model, include_dead=False, snow_fraction=0.2,
                           nodal_masses={crown: 50.0})
    assert _rel(mm2.totals()["snow"],
                0.2 * 1000.0 * model.geom.plan_area / G_ACCEL) < 1e-12
    model.loads.include_self_weight = False
    mm3 = md.lumped_masses(model, include_dead=False)
    assert _rel(mm3.total, model.steel_mass()) < 1e-12
    assert _rel(mm2.total, model.steel_mass() + mm2.totals()["snow"]
                + 50.0) < 1e-12
    # the modal result carries the same totals
    r = md.modal_analysis(model, n_modes=3, nodal_masses=[{"node": crown,
                                                           "mass": 50.0}])
    assert _rel(r.total_mass[2], mm.total + 50.0) < 1e-12
    assert _rel(r.mass_parts["dead nodal loads"], 7000.0 / G_ACCEL) < 1e-12


def _independent_KM(model, mm):
    """Dense K (with springs) and lumped M by plain loops, free DOFs only."""
    g = model.geom
    n = g.n_nodes
    K = np.zeros((3 * n, 3 * n))
    EA = model.EA()
    for e, (i, j) in enumerate(g.members):
        v = g.nodes[j] - g.nodes[i]
        L = math.sqrt(v @ v)
        c = v / L
        kk = EA[e] / L * np.outer(c, c)
        for (p, q_), sgn in (((i, i), 1), ((j, j), 1), ((i, j), -1),
                             ((j, i), -1)):
            K[3 * p:3 * p + 3, 3 * q_:3 * q_ + 3] += sgn * kk
    fixed, springs = model.support_conditions()
    for d, kv in springs.items():
        K[d, d] += kv
    Mdiag = np.repeat(mm.nodal, 3)
    fd = np.setdiff1d(np.arange(3 * n), fixed)
    return K[np.ix_(fd, fd)], Mdiag[fd], fd


def test_orthogonality_against_independent_assembly():
    """Phi'M Phi = I, Phi'K Phi = diag(omega^2) with K, M built by loops;
    omega^2 = numpy eigvalsh of M^-1/2 K M^-1/2 (independent solver)."""
    for support in (geo.PINNED, geo.ELASTIC):
        model = _dome_model(3, support=support)
        model.spring_stiffness = (2e7, 3e7, 8e7)
        r = md.modal_analysis(model, n_modes=15)
        mm = md.lumped_masses(model)
        K, Md, fd = _independent_KM(model, mm)
        assert np.array_equal(fd, r.free_dofs)
        P = r.modes[fd]
        assert np.max(np.abs(P.T @ (Md[:, None] * P) - np.eye(15))) < 1e-12
        H = P.T @ K @ P
        assert np.max(np.abs(H - np.diag(r.omega2))) < 1e-10 * r.omega2[-1]
        Dm = 1.0 / np.sqrt(Md)
        lam = np.linalg.eigvalsh(Dm[:, None] * K * Dm[None, :])[:15]
        assert np.allclose(r.omega2, lam, rtol=1e-9, atol=0), support
        _check_identities(r)


def test_effective_mass_sums_to_free_mass_every_support():
    """sum over ALL modes of M*_nd = iota_d' M iota_d (= free mass)."""
    custom = None
    for support in (geo.PINNED, geo.ROLLER3, geo.ELASTIC, geo.CUSTOM):
        model = _dome_model(2, R=6.0, support=support)
        if support == geo.CUSTOM:
            b = [int(x) for x in model.geom.base_nodes]
            custom = {nd: ([True, True, True] if k % 2 == 0
                           else [4e7, False, True]) for k, nd in enumerate(b)}
            model.custom_supports = custom
        for mtype in ("lumped", "consistent"):
            r = md.modal_analysis(model, n_modes=None, mass_type=mtype)
            assert r.n_modes == r.n_dynamic == r.n_free
            assert np.allclose(r.effective_mass.sum(axis=0), r.free_mass,
                               rtol=1e-10, atol=0), (support, mtype)
            assert np.allclose(r.participating_ratio, 1.0, rtol=0,
                               atol=1e-10)
            assert np.allclose(r.cumulative_ratio[-1], 1.0, atol=1e-10)
            _check_identities(r)
        # lumped free mass = total minus the restrained directions' masses
        r = md.modal_analysis(model, n_modes=None)
        fixed, _ = model.support_conditions()
        node_m = md.lumped_masses(model).nodal
        for a in range(3):
            held = sum(node_m[d // 3] for d in fixed if d % 3 == a)
            assert _rel(r.free_mass[a], node_m.sum() - held) < 1e-12
    assert custom is not None


def test_modal_superposition_reproduces_static_solution():
    """u = sum_n phi_n phi_n' F / omega_n^2 with all modes = K^-1 F."""
    for freq, support in ((2, geo.ELASTIC), (6, geo.PINNED)):
        model = _dome_model(freq, R=15.0, support=support,
                            section="CHS 114.3 x 5.0", dead=350.0)
        r = md.modal_analysis(model, n_modes=None)
        F = ld.basic_load_vectors(model)["D"] + ld.basic_load_vectors(
            model)["S"]
        st = solve_model(model, F)
        u_modal = r.modes @ ((r.modes.T @ F) / r.omega2)
        scale = np.max(np.abs(st.u))
        assert np.max(np.abs(u_modal - st.u)) < 1e-9 * scale, freq


def test_elastic_supports_converge_to_pinned():
    """Springs soften the dome; the error vanishes like 1/k as k -> inf."""
    pin = md.modal_analysis(_dome_model(3), n_modes=8)
    errs = []
    for k in (1e8, 1e10, 1e13):
        model = _dome_model(3, support=geo.ELASTIC)
        model.spring_stiffness = (k, k, k)
        r = md.modal_analysis(model, n_modes=8)
        assert np.all(r.omega2 <= pin.omega2 * (1 + 1e-12))   # softer
        errs.append(np.max(np.abs(r.omega2 / pin.omega2 - 1.0)))
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 2e-6
    assert 800.0 < errs[1] / errs[2] < 1200.0           # first order in 1/k


# ================================================== solvers / symmetry ======
def test_sparse_solver_matches_dense():
    model = _dome_model(4, R=12.0, section="CHS 114.3 x 5.0")
    d = md.modal_analysis(model, n_modes=14)
    s = md.modal_analysis(model, n_modes=14, dense_limit=0)
    assert "dense" in d.method and "sparse" in s.method
    assert np.allclose(s.omega2, d.omega2, rtol=1e-10, atol=0)
    _check_identities(s)
    # each repeated-frequency cluster spans the same subspace
    Mf = np.repeat(md.lumped_masses(model).nodal, 3)
    for a, b in md._clusters(d.omega2, d.omega2[-1]):
        if b > 13:
            continue                          # cluster may be cut at the end
        sv = np.linalg.svd(d.modes[:, a:b].T @ (Mf[:, None]
                                                * s.modes[:, a:b]),
                           compute_uv=False)
        assert np.allclose(sv, 1.0, atol=1e-8), (a, b, sv)
    assert np.allclose(s.effective_mass, d.effective_mass, rtol=1e-7,
                       atol=1e-7 * d.free_mass.max())


def test_repeated_sway_pair_aligned_with_x_and_y():
    sys.path.insert(0, os.path.join(ROOT, "examples"))
    from worked_example import build_model
    model = build_model()
    r = md.modal_analysis(model, n_modes=12)
    raw = md.modal_analysis(model, n_modes=12, align_degenerate=False)
    assert np.allclose(r.omega2, raw.omega2, rtol=1e-10)
    # modes 1 and 2 are the fundamental sway pair (equal frequency)
    assert _rel(r.omega2[1], r.omega2[0]) < 1e-8
    R = r.effective_mass_ratio
    assert R[0, 0] > 0.5 and R[0, 1] < 1e-12 and R[0, 2] < 1e-12
    assert R[1, 1] > 0.5 and R[1, 0] < 1e-12 and R[1, 2] < 1e-12
    assert _rel(R[0, 0], R[1, 1]) < 1e-8                 # five-fold symmetry
    # the rotation does not change the pair's total effective mass
    assert np.allclose(raw.effective_mass[:2].sum(axis=0),
                       r.effective_mass[:2].sum(axis=0), rtol=1e-10,
                       atol=1e-9 * r.free_mass.max())
    assert np.allclose(raw.effective_mass.sum(axis=0),
                       r.effective_mass.sum(axis=0), rtol=1e-9)
    _check_identities(r)


# ===================================================== pre-stress ==========
def test_prestressed_beaded_string():
    """Mass m between two bars of length a under tension N (node free in x
    and y): omega_y^2 = 2N/(a m) (geometric stiffness), omega_x^2 = 2EA/(a m);
    compression gives omega_y^2 < 0, reported as unstable."""
    EA, a, m, N = 1.0e7, 1.2, 30.0, 5.0e3
    nodes = np.array([[0.0, 0, 0], [a, 0, 0], [2 * a, 0, 0]])
    members = [[0, 1], [1, 2]]
    fixed = [0, 1, 2, 5, 6, 7, 8]
    r = md.solve_modes(nodes, members, [EA, EA], [0, m, 0], fixed,
                       n_modes=None, N_prestress=[N, N])
    assert _rel(r.omega2[0], 2 * N / (a * m)) < 1e-12
    assert _rel(r.omega2[1], 2 * EA / (a * m)) < 1e-12
    assert r.dominant_direction(0) == "y" and r.dominant_direction(1) == "x"
    _check_identities(r)
    c = md.solve_modes(nodes, members, [EA, EA], [0, m, 0], fixed,
                       n_modes=None, N_prestress=[-N, -N])
    assert _rel(c.omega2[0], -2 * N / (a * m)) < 1e-12
    assert c.unstable[0] and math.isnan(c.frequencies[0])
    assert not c.unstable[1]
    assert any("unstable" in w for w in c.warnings)


def test_prestress_at_lambda_cr_gives_zero_frequency():
    """det(K + lambda_cr Kg) = 0: the lowest frequency vanishes at the
    linear buckling load and its mode is the buckling mode."""
    model = _dome_model(3, R=10.0, section="CHS 60.3 x 3.2", snow=1500.0)
    combo = LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6)
    model.combinations = [LoadCombination("LC1  1.4D", D=1.4), combo]
    F = ld.combine(ld.basic_load_vectors(model), combo)
    st = solve_model(model, F)
    buck = bk.linear_buckling(model, st.N, n_modes=4)
    lam = buck.critical_factor
    assert lam > 1.0
    w0 = md.modal_analysis(model, n_modes=4)
    half = md.modal_analysis(model, n_modes=4, prestress="LC2",
                             prestress_factor=0.5 * lam)
    crit = md.modal_analysis(model, n_modes=4, prestress=combo.name,
                             prestress_factor=lam)
    assert "LC2" in crit.prestress
    # compression softens: 0 < omega_1^2(lambda_cr/2) < omega_1^2(0)
    assert 0.0 < half.omega2[0] < w0.omega2[0]
    assert abs(crit.omega2[0]) < 1e-7 * w0.omega2[0], crit.omega2[0]
    # the zero-frequency mode lies in the span of the critical buckling modes
    same = np.abs(buck.load_factors / lam - 1.0) < 1e-6
    B = buck.modes[:, same]
    phi = crit.modes[:, 0]
    coef, *_ = np.linalg.lstsq(B, phi, rcond=None)
    assert np.linalg.norm(B @ coef - phi) < 1e-5 * np.linalg.norm(phi)
    # forces given as an array give the same result
    arr = md.modal_analysis(model, n_modes=4, prestress=st.N,
                            prestress_factor=0.5 * lam)
    assert np.allclose(arr.omega2, half.omega2, rtol=1e-12)


# ============================================================ case study =====
def test_case_study_dome():
    sys.path.insert(0, os.path.join(ROOT, "examples"))
    from worked_example import build_model
    model = build_model()
    r = md.modal_analysis(model, n_modes=30)
    _check_identities(r)
    assert r.n_free == 3 * (model.geom.n_nodes - len(model.geom.base_nodes))
    assert 1.0 < r.fundamental_frequency < 100.0
    assert np.all(np.diff(r.omega2) >= -1e-9 * r.omega2[-1])
    assert _rel(r.total_mass[0], md.lumped_masses(model).total) < 1e-14
    # the fundamental sway mode carries most of the horizontal mass, but a
    # lattice shell spreads the rest over many high (membrane) modes
    assert r.effective_mass_ratio[0, 0] > 0.8
    assert r.modes_for_ratio(0.90)["x"] is None        # not within 30 modes
    rall = md.modal_analysis(model, n_modes=None)
    need = rall.modes_for_ratio(0.90)
    assert need["x"] is not None and abs(need["x"] - need["y"]) <= 1
    assert np.allclose(rall.participating_ratio, 1.0, atol=1e-10)
    assert np.allclose(rall.omega2[:30], r.omega2, rtol=1e-10)
    # consistent mass: higher frequencies, same total mass
    rc = md.modal_analysis(model, n_modes=6, mass_type="consistent")
    assert np.all(rc.omega2 > r.omega2[:6])
    _check_identities(rc)
    txt = r.summary_text()
    assert "Modal analysis" in txt and "Checks" in txt


# ======================================================= inputs / io ========
def test_input_validation():
    model = _dome_model(2, R=6.0)
    bad = [
        dict(nodal_masses={0: -1.0}),
        dict(nodal_masses={999: 1.0}),
        dict(nodal_masses=[1.0, 2.0]),
        dict(mass_type="diagonal"),
        dict(n_modes=0),
        dict(n_modes="some"),
        dict(snow_fraction=1.5),
        dict(prestress="LC9 does not exist"),
        dict(prestress=np.zeros(3)),
    ]
    for kw in bad:
        try:
            md.modal_analysis(model, **kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f"no ValueError for {kw}")
    r = md.modal_analysis(model, n_modes="all")
    assert r.n_modes == r.n_dynamic
    # a model with no free mass at all
    nodes = np.array([[0.0, 0, 0], [1.0, 0, 0]])
    try:
        md.solve_modes(nodes, [[0, 1]], [1e6], None, [0, 1, 2, 4, 5])
    except ValueError:
        pass
    else:
        raise AssertionError("massless model accepted")


def test_mechanism_on_both_solvers():
    """A bar chain along x with every node free in y has no transverse
    stiffness: zero frequencies, flagged, on the dense and on the sparse
    path (sigma = 0 fails on the singular K and a negative shift is used).
    The degenerate y cluster puts its y participation into its first mode."""
    N, h = 30, 0.5
    nodes, members, _ = _line_x(N, h)
    fixed = [0, 1, 2] + [3 * k + 2 for k in range(1, N + 1)]
    out = []
    for dl in (2400, 0):
        r = md.solve_modes(nodes, members, np.full(N, 1e8), None, fixed,
                           line_mass=5.0, n_modes=6, dense_limit=dl)
        assert np.all(r.frequencies == 0.0)
        assert any("mechanism" in w for w in r.warnings)
        assert r.effective_mass_ratio[0, 1] > 0.0
        assert np.all(r.effective_mass_ratio[1:, 1] < 1e-20)
        _check_identities(r)
        out.append(r)
    assert "sigma = -" in out[1].method
    # the axial modes are found when more modes are asked for
    r = md.solve_modes(nodes, members, np.full(N, 1e8), None, fixed,
                       line_mass=5.0, n_modes=None)
    pos = r.omega2[r.omega2 > 0]
    c = math.sqrt(1e8 / 5.0)
    assert _rel(math.sqrt(pos[0]), (2 * c / h) * math.sin(math.pi / (4 * N))) \
        < 1e-11


def test_frequency_load_curve():
    model = _dome_model(3, R=10.0, section="CHS 60.3 x 3.2", snow=1500.0)
    combo = LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6)
    model.combinations = [combo]
    F = ld.combine(ld.basic_load_vectors(model), combo)
    lam = bk.linear_buckling(model, solve_model(model, F).N,
                             n_modes=2).critical_factor
    cur = md.frequency_load_curve(model, "LC2", [0.0, 0.5 * lam, lam],
                                  n_modes=2)
    assert cur["omega2"].shape == (3, 2)
    w0 = md.modal_analysis(model, n_modes=2).omega2
    assert np.allclose(cur["omega2"][0], w0, rtol=1e-12)
    one = md.modal_analysis(model, n_modes=2, prestress=combo,
                            prestress_factor=0.5 * lam)
    assert np.allclose(cur["omega2"][1], one.omega2, rtol=1e-12)
    assert abs(cur["omega2"][2, 0]) < 1e-7 * w0[0]
    assert cur["combination"] == combo.name


def test_disconnected_node_is_left_out():
    nodes = np.array([[0.0, 0, 0], [2.0, 0, 0], [5.0, 5.0, 5.0]])
    r = md.solve_modes(nodes, [[0, 1]], [1e8], [0, 10.0, 0.0],
                       [0, 1, 2, 4, 5])
    assert r.n_free == 1
    assert _rel(r.omega2[0], 1e8 / (2.0 * 10.0)) < 1e-13
    assert any("neither stiffness nor mass" in w for w in r.warnings)


def test_csv_export_and_table():
    model = _dome_model(2, R=6.0)
    r = md.modal_analysis(model, n_modes=5)
    rows = r.table()
    assert len(rows) == 5 and rows[0]["mode"] == 1
    sh = r.mode_shape(0)
    assert abs(np.max(np.abs(sh)) - 1.0) < 1e-15
    with tempfile.TemporaryDirectory() as tmp:
        path = md.export_modal_csv(r, os.path.join(tmp, "modes.csv"))
        with open(path, newline="", encoding="utf-8") as fh:
            got = list(csv.DictReader(fh))
    assert len(got) == 5
    assert abs(float(got[0]["frequency_Hz"]) - r.frequencies[0]) \
        < 1e-12 * r.frequencies[0]


def test_core_import_is_headless():
    code = ("import sys; import geodome.modal; "
            "bad = [m for m in ('matplotlib', 'qtpy', 'PyQt6', 'PyQt5', "
            "'PySide6') if m in sys.modules]; print(bad); "
            "sys.exit(1 if bad else 0)")
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT,
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stdout + out.stderr


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
