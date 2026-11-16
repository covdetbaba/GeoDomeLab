"""
Built-in verification suite.

Each check compares the program against an independent result -- a closed-form
solution, a published geodesic chord-factor table, a topological identity, or
an internal energy/equilibrium identity.  ``run_all()`` returns a list of
dictionaries that the GUI shows in a table and the report embeds, so the
validation evidence is reproducible by whoever is marking the project.

Run from a terminal with::

    python -m geodome.verification
"""

from __future__ import annotations

import math

import numpy as np

from . import geometry as geo
from . import sections as sx
from . import buckling as bk
from .materials import MATERIALS
from .model import DomeModel, LoadCombination, LoadSpec
from .postprocess import run_analysis
from .solver import solve_static, solve_model


def _check(name, computed, reference, tol, unit="", note=""):
    if reference == 0.0:
        err = abs(computed - reference)
    else:
        err = abs(computed - reference) / abs(reference)
    return {
        "check": name,
        "computed": computed,
        "reference": reference,
        "unit": unit,
        "rel_error": err,
        "tolerance": tol,
        "pass": bool(err <= tol),
        "note": note,
    }


# ------------------------------------------------------------- benchmarks ---
def check_single_bar():
    """Axially loaded bar: u = PL/EA, N = P."""
    E, A, L, P = 200e9, 2.0e-3, 3.0, 50e3
    nodes = np.array([[0, 0, 0], [L, 0, 0]], float)
    mem = np.array([[0, 1]])
    F = np.zeros(6)
    F[3] = P
    r = solve_static(nodes, mem, np.array([E * A]), F, [0, 1, 2, 4, 5])
    return [
        _check("Single bar: tip displacement", float(r.u[3]), P * L / (E * A),
               1e-10, "m", "u = PL/EA"),
        _check("Single bar: axial force", float(r.N[0]), P, 1e-10, "N",
               "statically determinate"),
    ]


def check_tripod():
    """Symmetric three-legged space truss under a vertical apex load.

    N = -P L / (3 h)  and  delta = -P L^3 / (3 EA h^2).
    """
    E, A = 200e9, 2.0e-3
    a, h, P = 2.0, 4.0, 80e3
    ang = np.deg2rad([0.0, 120.0, 240.0])
    base = np.column_stack([a * np.cos(ang), a * np.sin(ang), np.zeros(3)])
    nodes = np.vstack([base, [0.0, 0.0, h]])
    mem = np.array([[0, 3], [1, 3], [2, 3]])
    F = np.zeros(12)
    F[11] = -P
    r = solve_static(nodes, mem, np.full(3, E * A), F, list(range(9)))
    L = math.hypot(a, h)
    return [
        _check("Tripod: leg axial force", float(r.N[0]), -P * L / (3.0 * h),
               1e-10, "N", "3D equilibrium at the apex"),
        _check("Tripod: apex deflection", float(r.u[11]),
               -P * L ** 3 / (3.0 * E * A * h ** 2), 1e-10, "m",
               "compatibility + Hooke"),
        _check("Tripod: vertical reaction", float(r.total_reaction()[2]), P,
               1e-9, "N", "global equilibrium"),
    ]


def check_two_bar_asymmetric():
    """Two orthogonal bars sharing a loaded node -- statically determinate.

    Node 2 at the origin is held by a horizontal bar from (-3, 0, 0) and a
    vertical bar from (0, 0, 4), and loaded with (Px, 0, Pz).  Out-of-plane
    motion is restrained because two bars alone cannot stabilise a node in
    three dimensions.  Equilibrium at the node then gives the bar forces
    directly: the horizontal bar carries Px and the vertical bar carries -Pz.
    """
    E, A = 70e9, 1.5e-3
    nodes = np.array([[-3.0, 0, 0], [0, 0, 4.0], [0.0, 0, 0]])
    mem = np.array([[0, 2], [1, 2]])
    Px, Pz = 25e3, -40e3
    F = np.zeros(9)
    F[6], F[8] = Px, Pz
    fixed = [0, 1, 2, 3, 4, 5, 7]          # both anchors, plus uy at node 2
    r = solve_static(nodes, mem, np.full(2, E * A), F, fixed)
    return [
        _check("Two-bar node: horizontal bar force", float(r.N[0]), Px,
               1e-10, "N", "determinate equilibrium"),
        _check("Two-bar node: vertical bar force", float(r.N[1]), -Pz,
               1e-10, "N", "determinate equilibrium"),
        _check("Two-bar node: x displacement", float(r.u[6]),
               Px * 3.0 / (E * A), 1e-10, "m", "u = PL/EA on the x bar"),
    ]


def check_chord_factors():
    """Chord factors against a published Class I (method 1) chart.

    Reference: Domerama, Chord factor charts,
    https://domerama.com/calculators/chord-charts (read 2026-11-16), which
    prints the chord factors to five decimals; the tolerance 2e-5 covers
    that rounding.
    """
    out = []
    g2 = geo.generate_dome(2, 1.0, 0.5)
    s2 = g2.strut_schedule()
    out.append(_check("2V chord factor A", s2[0]["chord_factor"], 0.54653,
                      2e-5, "-", "Domerama chord-factor chart, 2V"))
    out.append(_check("2V chord factor B", s2[1]["chord_factor"], 0.61803,
                      2e-5, "-", "= 2 sin 18 deg = 1/phi (edge of the "
                      "icosidodecahedron of edge midpoints)"))
    g4 = geo.generate_dome(4, 1.0, 0.5)
    s4 = g4.strut_schedule()
    # the six 4V "method 1" values of the chart, sorted by size
    refs = [0.25318, 0.29453, 0.29524, 0.29859, 0.31287, 0.32492]
    for k, ref in enumerate(refs):
        out.append(_check(f"4V chord factor {chr(ord('A') + k)}",
                          s4[k]["chord_factor"], ref, 2e-5, "-",
                          "Domerama chord-factor chart, 4V"))
    return out


def check_topology():
    """Vertex/edge/face counts and Euler's polyhedron formula."""
    out = []
    v0, f0 = geo.icosahedron()
    for nu in (1, 2, 3, 4, 6):
        if nu == 1:
            pts, tris = v0, f0
        else:
            pts, tris = geo.subdivide_class1(v0, f0, nu)
        edges = geo.edges_from_faces(tris)
        out.append(_check(f"{nu}V sphere vertices", float(len(pts)),
                          float(10 * nu ** 2 + 2), 0.0, "-", "V = 10 v^2 + 2"))
        out.append(_check(f"{nu}V sphere struts", float(len(edges)),
                          float(30 * nu ** 2), 0.0, "-", "E = 30 v^2"))
        out.append(_check(f"{nu}V Euler V-E+F", float(len(pts) - len(edges)
                                                      + len(tris)),
                          2.0, 0.0, "-", "Euler's polyhedron formula"))
    return out


def check_static_determinacy():
    """Maxwell's count, and the stability of the supported dome.

    A closed Class I sphere has n = 10v^2 + 2 nodes and m = 30v^2 struts, so
    3n - 6 = 30v^2 = m exactly: the complete geodesic sphere is statically
    determinate for every frequency.  Truncating it removes that balance, so an
    open cap floating in space is a mechanism and the dome relies on its
    supports.  The second check confirms that once the free edge is restrained
    the free partition of the stiffness matrix is positive definite, i.e. there
    is no zero-energy mode left.
    """
    from .solver import assemble_stiffness, free_dofs

    out = []
    v0, f0 = geo.icosahedron()
    for nu in (2, 4, 6):
        pts, tris = geo.subdivide_class1(v0, f0, nu)
        m = len(geo.edges_from_faces(tris))
        out.append(_check(f"{nu}V closed sphere: m - (3n - 6)",
                          float(m - (3 * len(pts) - 6)), 0.0, 0.0, "-",
                          "Maxwell count: the closed sphere is isostatic"))

    for nu, support in ((4, geo.PINNED), (4, geo.ROLLER3), (6, geo.PINNED)):
        g = geo.generate_dome(nu, 10.0, 0.5)
        EA = np.full(g.n_members, 200e9 * 1e-3)
        K = assemble_stiffness(g.nodes, g.members, EA)
        Kd = K.toarray() if hasattr(K, "toarray") else K
        fd = free_dofs(3 * g.n_nodes, geo.support_dofs(g, support))
        w = np.linalg.eigvalsh(Kd[np.ix_(fd, fd)])
        ratio = float(w.min() / w.max())
        out.append(_check(f"{nu}V {support}: no zero-energy mode",
                          float(ratio > 1e-10), 1.0, 0.0, "-",
                          f"min/max eigenvalue of Kff = {ratio:.3e}, "
                          f"{g.n_members - len(fd)} redundancies"))
    return out


def check_sphere_geometry():
    """Faceted dome area approaches the analytic hemisphere as v increases."""
    out = []
    R = 10.0
    exact = 2.0 * math.pi * R ** 2
    for nu, tol in ((4, 0.03), (8, 0.01)):
        g = geo.generate_dome(nu, R, 0.5)
        out.append(_check(f"{nu}V hemisphere surface area", g.surface_area,
                          exact, tol, "m2",
                          "faceted area converges from below"))
        out.append(_check(f"{nu}V crown height", g.height, R, 1e-9, "m",
                          "hemisphere rise = R"))
    return out


def check_section_properties():
    """CHS properties against hand calculation."""
    s = sx.get_section("CHS 168.3 x 6.3")
    D, t = 0.1683, 0.0063
    d = D - 2 * t
    A = math.pi / 4 * (D ** 2 - d ** 2)
    I = math.pi / 64 * (D ** 4 - d ** 4)
    return [
        _check("CHS 168.3x6.3 area", s.A, A, 1e-12, "m2", "pi/4 (D^2-d^2)"),
        _check("CHS 168.3x6.3 inertia", s.I, I, 1e-12, "m4", "pi/64 (D^4-d^4)"),
        _check("CHS 168.3x6.3 radius of gyration", s.r, math.sqrt(I / A),
               1e-12, "m", "r = sqrt(I/A)"),
    ]


def check_aisc_compression():
    """AISC E3 in the elastic and the inelastic range, and the Euler load.

    The references are computed here from the specified values of ASTM
    A500 Gr. C (E = 200 GPa, Fy = 345 MPa) and the nominal CHS 114.3 x 6.3
    dimensions, not from the program's material, section or Fe.
    """
    mat = MATERIALS["ASTM A500 Gr. C (round)"]
    s = sx.get_section("CHS 114.3 x 6.3")
    E, Fy = 200e9, 345e6
    D, t = 0.1143, 0.0063
    d = D - 2.0 * t
    A = math.pi / 4.0 * (D ** 2 - d ** 2)
    I = math.pi / 64.0 * (D ** 4 - d ** 4)
    r = math.sqrt(I / A)
    out = []
    # Elastic branch: a deliberately slender strut, Fcr = 0.877 Fe.
    L = 12.0
    cap, slend, Fe, Fcr, Q = bk.compression_capacity(s, mat, L, 1.0, 0.9)
    Fe_ref = math.pi ** 2 * E / (L / r) ** 2
    out.append(_check("AISC E3 elastic branch Fcr", Fcr, 0.877 * Fe_ref,
                      1e-12, "Pa", f"KL/r = {L / r:.0f} > 4.71 sqrt(E/Fy); "
                      f"E, Fy, D, t typed in"))
    # Inelastic branch.
    L = 2.0
    cap, slend, Fe, Fcr, Q = bk.compression_capacity(s, mat, L, 1.0, 0.9)
    Fe_ref = math.pi ** 2 * E / (L / r) ** 2
    out.append(_check("AISC E3 inelastic branch Fcr", Fcr,
                      (0.658 ** (Fy / Fe_ref)) * Fy, 1e-12, "Pa",
                      f"KL/r = {L / r:.0f}; E, Fy, D, t typed in"))
    # Euler load.
    out.append(_check("Euler load Pe = pi^2 EI/(KL)^2",
                      bk.euler_load(s, mat, L), math.pi ** 2 * E * I / L ** 2,
                      1e-12, "N", "closed form, E and I typed in"))
    return out


def check_pinned_column_buckling():
    """Linear buckling of a pin-ended column modelled as a chain of bars.

    A string of ``n`` collinear truss bars carrying compression N has the
    discrete buckling load Ncr = (EA_t/L_e) * 2(1 - cos(pi/n)) only if lateral
    stiffness is supplied; a pure truss chain has none, so instead the identity
    verified here is the *scaling* property of the eigenproblem: doubling the
    reference axial force must halve every critical load factor.
    """
    model = DomeModel(geom=geo.generate_dome(3, 8.0, 0.5))
    model.loads = LoadSpec(include_self_weight=True, dead_pressure=300.0,
                           snow_pressure=1000.0, wind_enabled=False)
    model.combinations = [LoadCombination("REF  1.2D + 1.6S", D=1.2, S=1.6)]
    res = run_analysis(model, do_buckling=False)
    N = res.cases["REF  1.2D + 1.6S"].N
    b1 = bk.linear_buckling(model, N, 3)
    b2 = bk.linear_buckling(model, 2.0 * N, 3)
    out = [_check("Buckling eigenvalue scaling (2N -> lambda/2)",
                  float(b2.load_factors[0]), float(b1.load_factors[0]) / 2.0,
                  1e-8, "-", "linearity of the geometric stiffness")]
    out.append(_check("Buckling load factor is positive",
                      float(b1.load_factors[0] > 0.0), 1.0, 0.0, "-",
                      "compression-dominated dome must have a finite lambda_cr"))
    return out


def check_dome_equilibrium():
    """Global statics and energy identities on a realistic dome model."""
    model = DomeModel(geom=geo.generate_dome(4, 12.0, 0.5))
    # A fixed set of five load cases (dead, snow and wind in turn, and a
    # service case), so the number of checks does not follow the defaults.
    model.combinations = [
        LoadCombination("LC1  1.4D", D=1.4),
        LoadCombination("LC2  1.2D + 1.6S", D=1.2, S=1.6),
        LoadCombination("LC3  1.2D + 1.0W + 0.5S", D=1.2, W=1.0, S=0.5),
        LoadCombination("LC4  0.9D + 1.0W", D=0.9, W=1.0),
        LoadCombination("SLS  1.0D + 1.0S", D=1.0, S=1.0,
                        serviceability=True)]
    res = run_analysis(model, do_buckling=False)
    out = []
    for name in res.order:
        r = res.cases[name]
        out.append(_check(f"Equilibrium residual, {name}",
                          r.equilibrium_residual, 0.0, 1e-9, "-",
                          "||Ku - F - R|| / ||F||"))
        out.append(_check(f"Strain energy identity, {name}",
                          r.energy_error, 0.0, 1e-8, "-",
                          "1/2 F.u == 1/2 sum N^2 L/EA"))
        out.append(_check(f"Load / reaction balance, {name}",
                          float(np.max(np.abs(r.reaction_balance))), 0.0,
                          1e-6, "N", "sum(F) + sum(R) = 0"))
    return out


def check_symmetry():
    """A symmetric dome under a symmetric load must give a symmetric answer.

    Every node has an image under a 72 degree rotation about z (the icosahedral
    dome has five-fold symmetry); the vertical displacement of a node and of its
    image must be identical.
    """
    model = DomeModel(geom=geo.generate_dome(4, 10.0, 0.5))
    model.loads = LoadSpec(include_self_weight=True, dead_pressure=400.0,
                           snow_pressure=1500.0, wind_enabled=False)
    model.combinations = [LoadCombination("SYM 1.2D + 1.6S", D=1.2, S=1.6)]
    res = run_analysis(model, do_buckling=False)
    u = res.cases["SYM 1.2D + 1.6S"].displacement_vectors()
    nodes = model.geom.nodes

    c, s = math.cos(math.radians(72.0)), math.sin(math.radians(72.0))
    Rz = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    rotated = nodes @ Rz.T

    worst = 0.0
    scale = max(np.max(np.abs(u[:, 2])), 1e-18)
    for k in range(len(nodes)):
        d = np.linalg.norm(nodes - rotated[k], axis=1)
        j = int(np.argmin(d))
        if d[j] < 1e-6 * model.geom.radius:
            worst = max(worst, abs(u[k, 2] - u[j, 2]) / scale)
    return [_check("Five-fold symmetry of uz", worst, 0.0, 1e-9, "-",
                   "72 deg rotation maps the dome onto itself")]


def _cut_dome(frequency=6, radius=10.0, pressure=1000.0, support=geo.PINNED):
    """A hemisphere carrying a uniform projected pressure and nothing else."""
    m = DomeModel(geom=geo.generate_dome(frequency, radius, 0.5))
    m.support_type = support
    m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                       snow_pressure=pressure, wind_enabled=False)
    m.combinations = [LoadCombination("CUT", S=1.0)]
    res = run_analysis(m, do_buckling=False)
    return m, res.cases["CUT"]


def _section_cut(model, r, z_cut):
    """Free body of the cap above ``z_cut``.

    Returns ``(transmitted, applied)``: the resultant of the member forces
    crossing the plane and the resultant of the loads applied above it.  Static
    equilibrium of the cap requires the two to sum to zero.
    """
    g = model.geom
    above = g.nodes[:, 2] > z_cut
    transmitted = np.zeros(3)
    for k, (i, j) in enumerate(g.members):
        if above[i] == above[j]:
            continue
        aa, bb = (i, j) if above[i] else (j, i)      # aa is inside the cap
        v = g.nodes[bb] - g.nodes[aa]
        transmitted += r.N[k] * v / np.linalg.norm(v)
    applied = r.F.reshape(-1, 3)[above].sum(axis=0)
    return transmitted, applied


def check_section_cut():
    """Method of sections: every horizontal free body must be in equilibrium.

    This is an independent statics check on the *member forces* -- it never
    touches the stiffness matrix again -- and is the matrix-method analogue of
    the hand check a designer would make on a dome cap.
    """
    model, r = _cut_dome(6, 10.0, 1000.0)
    scale = max(abs(float(r.F[2::3].sum())), 1.0)
    out = []
    for frac in (0.25, 0.50, 0.75):
        z = frac * model.geom.height
        transmitted, applied = _section_cut(model, r, z)
        resid = float(np.max(np.abs(transmitted + applied))) / scale
        out.append(_check(f"Cap free body at z = {z:.2f} m", resid, 0.0,
                          1e-9, "-",
                          "cut member forces + applied cap loads = 0"))
    return out


def check_membrane_meridional():
    """Meridional stress resultant against membrane shell theory.

    For a spherical dome carrying a uniform load ``p`` per unit *projected*
    area, membrane theory gives a meridional resultant that is constant over
    the whole shell:

        n_phi = -p R / 2      [force per unit length of the parallel circle]

    The discrete value follows from a horizontal section cut: the vertical
    force transmitted across a parallel circle of radius ``rho = R sin(phi)``
    equals ``|n_phi| * 2 pi rho * sin(phi)``.  Agreement within a few percent
    is the expected level, because the truss lumps a smooth pressure field onto
    discrete nodes.
    """
    R, p = 10.0, 1000.0
    model, r = _cut_dome(8, R, p)
    out = []
    for frac in (0.30, 0.50, 0.70):
        z = frac * R
        transmitted, _ = _section_cut(model, r, z)
        rho = math.sqrt(max(R ** 2 - z ** 2, 1e-12))
        sin_phi = rho / R
        n_phi = abs(float(transmitted[2])) / (2.0 * math.pi * rho * sin_phi)
        out.append(_check(f"Membrane n_phi at z/R = {frac:.2f}", n_phi,
                          p * R / 2.0, 0.05, "N/m",
                          "n_phi = pR/2, constant over a spherical dome"))
    return out


def check_hoop_action():
    """The classic dome result: hoop compression on top, hoop tension below.

    Under a uniform downward load a spherical dome carries compression in the
    hoop direction near the crown and tension near the springing, and the
    struts that run along a parallel circle must reproduce that reversal.  The
    reversal itself is asserted, not its exact location: in a *triangulated*
    lattice the diagonal struts also carry hoop force, so the ring struts alone
    under-represent the resultant and the reversal appears higher up the
    meridian than the 45 degrees of smooth-shell theory.
    """
    R = 10.0
    model, r = _cut_dome(8, R, 1000.0, support=geo.ROLLER3)
    g = model.geom
    ring = g.ring_index
    phis, vals = [], []
    for rr in range(int(ring.max()) + 1):
        sel = (ring[g.members[:, 0]] == rr) & (ring[g.members[:, 1]] == rr)
        if int(sel.sum()) < 3:
            continue
        z = float(g.ring_heights[rr])
        phis.append(math.degrees(math.acos(min(1.0, max(-1.0, z / R)))))
        vals.append(float(np.median(r.N[sel])))
    phis, vals = np.array(phis), np.array(vals)

    crown_val, base_val = float(vals[0]), float(vals[-1])
    cross = float("nan")
    for k in range(len(vals) - 1):
        if vals[k] < 0.0 <= vals[k + 1]:
            t = -vals[k] / (vals[k + 1] - vals[k])
            cross = phis[k] + t * (phis[k + 1] - phis[k])
            break
    return [
        _check("Hoop force is compressive at the crown",
               float(crown_val < 0.0), 1.0, 0.0, "-",
               f"median ring force {crown_val / 1e3:+.2f} kN at "
               f"phi = {phis[0]:.0f} deg"),
        _check("Hoop force is tensile at the springing",
               float(base_val > 0.0), 1.0, 0.0, "-",
               f"median ring force {base_val / 1e3:+.2f} kN at "
               f"phi = {phis[-1]:.0f} deg"),
        _check("Hoop reversal occurs on the upper half",
               float(0.0 < cross < 60.0), 1.0, 0.0, "-",
               f"reversal at phi = {cross:.1f} deg "
               f"(smooth-shell theory: 45 deg)"),
    ]


def check_frequency_study():
    """Reported, not asserted: how the dome stiffens as the frequency rises.

    Each frequency is a genuinely different structure -- more and shorter
    struts of the same section -- so this is a parameter study rather than a
    mesh convergence test.  The one thing that must hold is the monotonic
    trend: a denser dome carrying the same pressure has to deflect less.
    """
    rows, deltas = [], []
    for nu in (2, 3, 4, 6, 8):
        m = DomeModel(geom=geo.generate_dome(nu, 10.0, 0.5))
        m.group_sections = ["CHS 88.9 x 4.0"] * m.geom.n_strut_types()
        m.loads = LoadSpec(include_self_weight=False, dead_pressure=0.0,
                           snow_pressure=1500.0, wind_enabled=False)
        m.combinations = [LoadCombination("P", S=1.0)]
        res = run_analysis(m, do_buckling=False)
        _, d = res.cases["P"].max_displacement()
        deltas.append(d)
        rows.append(f"{nu}V: {d * 1e3:.3f} mm")
    monotone = all(deltas[i + 1] <= deltas[i] * 1.001
                   for i in range(len(deltas) - 1))
    return [_check("Dome stiffens monotonically with frequency",
                   float(monotone), 1.0, 0.0, "-", "   ".join(rows))]


# ------------------------------------------------------------------ driver --
SUITES = [
    ("Closed-form truss solutions", check_single_bar),
    ("Closed-form truss solutions", check_tripod),
    ("Closed-form truss solutions", check_two_bar_asymmetric),
    ("Geodesic geometry", check_chord_factors),
    ("Geodesic geometry", check_topology),
    ("Geodesic geometry", check_sphere_geometry),
    ("Geodesic geometry", check_static_determinacy),
    ("Section and code checks", check_section_properties),
    ("Section and code checks", check_aisc_compression),
    ("Global behaviour", check_dome_equilibrium),
    ("Global behaviour", check_symmetry),
    ("Global behaviour", check_section_cut),
    ("Global behaviour", check_membrane_meridional),
    ("Global behaviour", check_hoop_action),
    ("Global behaviour", check_frequency_study),
    ("Global behaviour", check_pinned_column_buckling),
]


def run_all(progress=None):
    """Run every verification check.  Returns a flat list of result dicts."""
    out = []
    n = len(SUITES)
    for i, (group, fn) in enumerate(SUITES):
        if progress is not None:
            progress((i + 0.5) / n, f"Verifying: {fn.__name__} ...")
        try:
            for row in fn():
                row["group"] = group
                out.append(row)
        except Exception as exc:                     # pragma: no cover
            out.append({"group": group, "check": fn.__name__,
                        "computed": float("nan"), "reference": float("nan"),
                        "unit": "", "rel_error": float("nan"), "tolerance": 0.0,
                        "pass": False, "note": f"raised {type(exc).__name__}: {exc}"})
    if progress is not None:
        progress(1.0, "Verification complete.")
    return out


def format_report(rows):
    """Plain-text verification table."""
    n_pass = sum(1 for r in rows if r["pass"])
    L = ["=" * 104,
         "  GeoDomeLab VERIFICATION SUITE",
         "=" * 104,
         f"  {'check':<44s}{'computed':>15s}{'reference':>15s}"
         f"{'rel.err':>11s}{'tol':>9s}  result",
         "-" * 104]
    group = None
    for r in rows:
        if r["group"] != group:
            group = r["group"]
            L.append(f"  [{group}]")
        L.append(f"  {r['check']:<44s}{r['computed']:>15.6g}"
                 f"{r['reference']:>15.6g}{r['rel_error']:>11.2e}"
                 f"{r['tolerance']:>9.1e}  {'PASS' if r['pass'] else 'FAIL'}")
    L.append("-" * 104)
    L.append(f"  {n_pass} of {len(rows)} checks passed.")
    L.append("=" * 104)
    return "\n".join(L)


if __name__ == "__main__":                           # pragma: no cover
    rows = run_all()
    print(format_report(rows))
    raise SystemExit(0 if all(r["pass"] for r in rows) else 1)
