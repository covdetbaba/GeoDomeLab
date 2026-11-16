"""
Member strength checks (AISC 360 LRFD) and global linear buckling analysis.

Member checks
-------------
Tension, yielding on the gross section (Chapter D):

    phi_t Pn = 0.90 * Fy * Ag

Compression, flexural buckling of a doubly symmetric member (Chapter E3):

    Fe   = pi^2 E / (K L / r)^2
    Fcr  = 0.658^(Fy/Fe) * Fy            if  K L / r <= 4.71 sqrt(E/Fy)
    Fcr  = 0.877 Fe                      otherwise
    phi_c Pn = 0.90 * Fcr * Ag

Round HSS with D/t > 0.11 E/Fy have slender walls (AISC 360-16 Table B4.1a,
case 9).  For them AISC 360-16 Section E7 replaces the gross area by the
effective area (Eq. E7-7):

    Ae = Ag * [ 0.038 E / (Fy D/t) + 2/3 ]      0.11 E/Fy < D/t < 0.45 E/Fy
    phi_c Pn = 0.90 * Fcr * Ae                  (Fcr from E3 with Fy)

Above D/t = 0.45 E/Fy the specification does not apply and the capacity is
reported as zero.  The ratio Ae/Ag is returned in the ``Q`` slot for the
diagnostics tables (it is 1 for a non-slender wall).  Every CHS in the
built-in catalogue has D/t <= 27.4, below the slender limit of all built-in
steels, so this clause only acts for user-defined thin sections.

``model.design_thickness_factor`` (default 1.0) scales the wall thickness
used in the strength checks.  AISC 360-16 Section B4.2 requires 0.93 t_nom
for ERW HSS such as ASTM A500; the analysis stiffness always uses the nominal
section.

Global buckling
---------------
A single-layer dome fails by *system* instability -- snap-through of a node or
of a whole cap -- long before the individual struts reach their Euler loads.
That mode is captured by the linearised (eigenvalue) buckling problem

    [ K + lambda Kg(N) ] {phi} = {0}

where ``Kg`` is the geometric stiffness assembled from the axial forces of the
reference load case.  Rearranged into a symmetric generalised problem with a
positive-definite right-hand side,

    (-Kg) {phi} = (1/lambda) K {phi}

so that the *largest* eigenvalue ``mu`` of the pair gives the *smallest*
critical load factor ``lambda_cr = 1/mu``.  ``lambda_cr`` is the multiple of
the applied load at which the dome becomes unstable.

Supports come from ``model.support_conditions()``: rigid restraints are
removed from the problem and grounded springs (ELASTIC / CUSTOM supports) add
their stiffness ``k`` to the diagonal of ``K`` (they carry no axial force, so
they add nothing to ``Kg``).

``lambda_cr`` is the *bifurcation* load of the perfect structure.  Single-
layer domes are imperfection sensitive: the limit load of the real, slightly
imperfect dome can be well below ``lambda_cr``, and by how much depends on the
geometry and on the imperfection.  GeoDomeLab quantifies it with the
geometrically nonlinear imperfection study (``nonlinear.imperfection_study``)
instead of applying a generic knock-down factor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

try:
    import scipy.linalg as sla
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False

from .solver import assemble_stiffness, assemble_geometric_stiffness, \
    direction_cosines, free_dofs


# ------------------------------------------------------- member capacities --
def tension_capacity(section, material, phi_t=0.90):
    """phi_t * Pn for yielding on the gross section [N]."""
    return phi_t * material.Fy * section.A


def local_buckling_Q(section, material):
    """Effective-area ratio Ae/Ag of a round HSS (AISC 360-16 Eq. E7-7).

    1.0 for a non-slender wall, NaN above the scope limit 0.45 E/Fy.
    """
    E, Fy = material.E, material.Fy
    Dt = section.D_over_t
    if Dt <= 0.11 * E / Fy:
        return 1.0
    if Dt > 0.45 * E / Fy:
        return float("nan")          # outside the scope of the specification
    return 0.038 * E / (Fy * Dt) + 2.0 / 3.0


def compression_capacity(section, material, L, K=1.0, phi_c=0.90):
    """phi_c * Pn for flexural buckling [N], plus diagnostics.

    AISC 360-16 E3 (flexural buckling, Fcr from Fy) with the E7 effective
    area for slender round HSS.  Returns ``(phiPn, slenderness, Fe, Fcr, Q)``
    where ``Q = Ae/Ag``.
    """
    E, Fy = material.E, material.Fy
    slend = K * L / section.r
    slend = max(slend, 1e-9)
    Fe = math.pi ** 2 * E / slend ** 2
    Q = local_buckling_Q(section, material)
    if not math.isfinite(Q):
        Q = 0.0
    if slend <= 4.71 * math.sqrt(E / Fy):
        Fcr = (0.658 ** (Fy / Fe)) * Fy
    else:
        Fcr = 0.877 * Fe
    return phi_c * Fcr * Q * section.A, slend, Fe, Fcr, Q


def design_section(section, factor=1.0):
    """Section with the design wall thickness ``factor * t`` (AISC B4.2)."""
    if factor == 1.0:
        return section
    from .sections import Section
    return Section(section.name, section.D, section.t * float(factor))


def euler_load(section, material, L, K=1.0):
    """Elastic Euler critical load of one strut [N]."""
    return math.pi ** 2 * material.E * section.I / (K * L) ** 2


def member_capacities(model):
    """Vectorised capacities for every member of a model.

    Returns a dict of arrays: ``phiPn_t``, ``phiPn_c``, ``slenderness``,
    ``Fe``, ``Fcr``, ``Pe_euler``, ``Q``.
    """
    geom = model.geom
    mat = model.material
    secs = model.member_sections()
    L = geom.member_lengths()
    K = model.K_effective

    n = len(secs)
    out = {k: np.zeros(n) for k in
           ("phiPn_t", "phiPn_c", "slenderness", "Fe", "Fcr", "Pe_euler", "Q")}
    tf = float(getattr(model, "design_thickness_factor", 1.0))
    for m, s in enumerate(secs):
        s = design_section(s, tf)
        out["phiPn_t"][m] = tension_capacity(s, mat, model.phi_t)
        cap, sl, Fe, Fcr, Q = compression_capacity(s, mat, L[m], K, model.phi_c)
        out["phiPn_c"][m] = cap
        out["slenderness"][m] = sl
        out["Fe"][m] = Fe
        out["Fcr"][m] = Fcr
        out["Q"][m] = Q
        out["Pe_euler"][m] = euler_load(s, mat, L[m], K)
    return out


def utilisation(model, N, caps=None):
    """Demand/capacity ratio for every member (tension and compression)."""
    if caps is None:
        caps = member_capacities(model)
    cap = np.where(N >= 0.0, caps["phiPn_t"], caps["phiPn_c"])
    return np.abs(N) / np.maximum(cap, 1e-9), caps


# ------------------------------------------------------- global buckling ----
@dataclass
class BucklingResult:
    load_factors: np.ndarray                      # (k,) ascending, positive
    modes: np.ndarray                             # (3n, k) full-DOF vectors
    reference_case: str = ""
    info: str = ""
    warnings: list = field(default_factory=list)

    @property
    def critical_factor(self):
        return float(self.load_factors[0]) if len(self.load_factors) else float("nan")

    def mode(self, k):
        return self.modes[:, k].reshape(-1, 3)


def linear_buckling(model, N_reference, n_modes=6, dense_limit=2400):
    """Linearised buckling analysis about the stress state ``N_reference``.

    ``N_reference`` is the member axial force array from a static solve of the
    reference load combination.  The returned load factors multiply *that*
    combination.
    """
    geom = model.geom
    nodes, members = geom.nodes, geom.members
    ndof = 3 * len(nodes)
    L, C = direction_cosines(nodes, members)

    K = assemble_stiffness(nodes, members, model.EA(), L, C)
    Kg = assemble_geometric_stiffness(nodes, members, np.asarray(N_reference), L, C)

    fixed, springs = model.support_conditions()
    fd = free_dofs(ndof, fixed)
    warnings = []
    if springs:
        fixed_set = set(np.asarray(fixed, dtype=int).tolist())
        sd = [int(d) for d in springs if int(d) not in fixed_set]
        kv = [float(springs[d]) for d in sd]
        if _HAVE_SCIPY and sp.issparse(K):
            K = (K + sp.coo_matrix((kv, (sd, sd)), shape=K.shape)).tocsr()
        else:                                        # pragma: no cover
            K = np.array(K, dtype=float)
            K[sd, sd] += kv

    if _HAVE_SCIPY and sp.issparse(K):
        Kff = K[fd][:, fd]
        Kgff = Kg[fd][:, fd]
    else:                                            # pragma: no cover
        Kff = K[np.ix_(fd, fd)]
        Kgff = Kg[np.ix_(fd, fd)]

    nf = len(fd)
    # The dense solver returns every eigenpair; only ARPACK needs k < nf - 1.
    cap = nf if nf <= dense_limit else nf - 2
    n_modes = int(max(1, min(n_modes, cap)))

    if nf <= dense_limit:
        A = (-Kgff).toarray() if _HAVE_SCIPY and sp.issparse(Kgff) else -np.asarray(Kgff)
        B = Kff.toarray() if _HAVE_SCIPY and sp.issparse(Kff) else np.asarray(Kff)
        A = 0.5 * (A + A.T)
        B = 0.5 * (B + B.T)
        if _HAVE_SCIPY:
            mu, V = sla.eigh(A, B)
        else:                                        # pragma: no cover
            mu, V = np.linalg.eigh(np.linalg.solve(B, A))
        info = f"Dense symmetric generalised eigensolve, {nf} free DOF."
    else:
        # Lanczos on the pair (-Kg, K) in regular generalised mode (no
        # shift); the algebraically largest mu = 1/lambda come first.
        mu, V = spla.eigsh(-Kgff.tocsc(), k=min(n_modes + 4, nf - 2),
                           M=Kff.tocsc(), which="LA",
                           maxiter=5000, tol=1e-8)
        info = f"Sparse Lanczos eigensolve (ARPACK), {nf} free DOF."

    order = np.argsort(mu)[::-1]                     # largest mu first
    mu, V = mu[order], V[:, order]

    keep = mu > 1e-12
    if not np.any(keep):
        warnings.append("No positive buckling load factor was found: under this "
                        "load case the dome is predominantly in tension.")
        return BucklingResult(np.array([]), np.zeros((ndof, 0)),
                              info=info, warnings=warnings)

    mu, V = mu[keep], V[:, keep]
    lam = 1.0 / mu
    k = min(n_modes, len(lam))
    lam, V = lam[:k], V[:, :k]

    # Scale each mode so that its largest component is +1.  The sign of an
    # eigenvector is arbitrary (it differs between LAPACK builds); this
    # makes the modes the same on every platform.
    modes = np.zeros((ndof, k))
    for i in range(k):
        v = V[:, i]
        j = int(np.argmax(np.abs(v)))
        piv = v[j] if v[j] != 0.0 else 1.0
        modes[fd, i] = v / piv

    if lam[0] < 1.0:
        warnings.append(f"Critical load factor {lam[0]:.3f} < 1.0: the dome is "
                        f"unstable under the applied load.")
    elif lam[0] < 2.0:
        warnings.append(f"Critical load factor {lam[0]:.3f} is low. Single-layer "
                        f"domes are imperfection sensitive, so the real limit "
                        f"load can be well below this bifurcation value; "
                        f"quantify it with the nonlinear imperfection study.")

    return BucklingResult(lam, modes, info=info, warnings=warnings)
