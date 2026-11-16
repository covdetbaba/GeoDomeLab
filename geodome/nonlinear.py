"""
Geometrically nonlinear analysis of pin-jointed space trusses.

Large displacements, small strains, linear elastic bars.  This module traces
the equilibrium path ``lambda -> u`` of a truss under a proportionally
increasing reference load ``lambda * F_ref`` through limit points and
snap-through, and compares it with the first-order (linear) and the
linearised second-order (P-Delta) answers.

Corotational bar element (Crisfield 1991, Vol. 1, Ch. 3; Yaw 2011)
------------------------------------------------------------------
For bar ``i -> j`` with original vector ``d0 = X_j - X_i`` (length ``L``)
and current vector ``d = d0 + u_j - u_i`` (length ``l``):

    c   = d / l                               current unit vector
    e   = l - L = (l^2 - L^2) / (l + L)       (cancellation-free form)
    eps = e / L                               engineering strain
    N   = EA (eps - eps0)                     axial force, tension +
    f   = N [ -c ; +c ]                       internal (resisting) forces

Differentiating ``f`` exactly (``dl = c . dd`` and ``dc = (I - c c^T) dd / l``)
gives the consistent tangent

    k3  = (EA/L) c c^T + (N/l) (I - c c^T)
    k_T = [ k3  -k3 ; -k3  k3 ]               (6 x 6)

i.e. a material part and the geometric (initial-stress) part.  It is the
form used by OpenSees ``CorotTruss`` (engineering strain on the original
length, same tangent); the finite-difference test in
``tests/test_nonlinear.py`` checks it to round-off.  At ``u = 0`` (N = 0)
``k_T`` is the linear bar stiffness of ``solver.element_stiffness``.

Grounded springs ``{dof: k}`` (ELASTIC / CUSTOM supports) stay linear and
fixed in direction: they add ``k u`` to the internal forces and ``k`` to the
tangent diagonal.  Restrained DOFs are removed.

Equilibrium and the solution procedures
---------------------------------------
Out-of-balance force on the free DOFs ``r(u, lambda) = lambda F_ref - f(u)``.
Every procedure uses full Newton-Raphson (tangent re-formed every iteration)
with ``K_T dv_r = r`` and ``K_T dv_t = F_ref`` (McGuire, Gallagher & Ziemian
2000, Ch. 12; Bathe 2014; Crisfield 1991, Ch. 9).

* **Load control** -- ``lambda`` fixed in a step, ``du = dv_r``.  It cannot
  pass a limit point; the step is cut back until it brackets the limit load.
* **Displacement control** on one DOF ``c`` -- ``du_c`` fixed in a step:
  predictor ``dlambda = du_c / v_t,c``, corrector
  ``dlambda = -v_r,c / v_t,c``, ``du = v_r + dlambda v_t``.  Passes load
  limit points; fails at a snap-back in ``u_c``.
* **Cylindrical arc-length** (Crisfield 1981; Crisfield 1991, Ch. 9; the
  idea goes back to Riks 1979) -- the step increment ``(Du, Dlambda)`` obeys

      Du^T Du = Dl^2

  Predictor ``Dlambda = s Dl / |v_t|`` with the sign ``s`` of
  ``Du_prev^T v_t`` (the previous converged increment), so the path keeps
  going forward through limit points.  Corrector: with
  ``Du_new = Du + v_r + dlambda v_t`` the constraint is the quadratic

      a dlambda^2 + b dlambda + c = 0,   a = v_t.v_t,
      b = 2 v_t.(Du + v_r),   c = |Du + v_r|^2 - Dl^2

  and the root giving the smallest angle between ``Du`` and ``Du_new`` is
  kept (Crisfield 1981; Hellweg & Crisfield 1998 discuss its limits at sharp
  snap-backs).  Complex roots or no convergence -> the step is cut back
  (``Dl / 2``).  After convergence ``Dl <- Dl sqrt(n_desired / n_iter)``
  (clamped to [0.25, 2] per step).

Convergence of a step needs both

    |r| <= tol * max(|lambda F_ref|, lambda_ref |F_ref|, |f_full|)
    |du_last| <= tol_disp * |Du|

where ``lambda_ref`` is the largest |lambda| reached so far and ``f_full``
the internal force vector including the support DOFs.  The two extra terms
in the force scale only matter where the path crosses ``lambda = 0`` (the
von Mises truss at ``w = h`` and ``w = 2h``) and keep the test meaningful
there; elsewhere the criterion is ``|r| / |lambda F_ref| <= tol``.

Stability indicators (every converged point)
--------------------------------------------
* ``det_sign`` -- sign of det K_T from the sparse LU factors (pivots and
  row/column permutations).  It flips each time an odd number of
  eigenvalues crosses zero; a pair of equal eigenvalues crossing together
  (a double bifurcation, common in symmetric domes) leaves it unchanged.
* ``n_negative``, ``min_eig`` -- dense eigenvalues of K_T for systems up to
  ``eig_limit`` free DOFs (-1 / NaN above it).  The negative count also sees
  a double crossing.  Above ``eig_limit`` only the determinant sign is
  available, and the run says so in its messages.
* ``stiffness_parameter`` S = F^T K_0^-1 F / F^T K_T^-1 F, 1 at the start.
  S -> 0 at a **limit point** but stays finite at a **bifurcation**.  When
  the critical mode does little work on the load, S only collapses very
  close to the limit, so the classification uses the load projection of the
  critical eigenvector phi (smallest |eigenvalue| of K_T, inverse
  iteration) at the refined critical point:

      rho = |phi^T F_ref| / (|phi| |F_ref|)

  rho is finite at a limit point and at round-off level (~1e-14) when the
  critical mode is orthogonal to the load (bifurcation of a symmetric
  structure).  rho >= 1e-6 -> limit point.  A stability change while
  ``lambda`` keeps rising is reported as a bifurcation.

The first load maximum is refined after the path has passed it: ``lambda``
is maximised (Brent's method) over a one-parameter family of fully converged
equilibrium points (cylindrical distance from the last point before the peak,
or the control displacement), so the reported ``lambda_limit`` is an exact
equilibrium point, not a sample.  When that family is not smooth (a sharp
limit or a branch switch at a bifurcation) the stability change is bracketed
instead by load-controlled bisection along the rising primary path (the
bracket closes to a relative width of 1e-10, but next to a bifurcation the
nearly singular tangent limits the accuracy of the bracketed load to about
1e-9; a point is accepted when its force residual has converged, and a load
step without equilibrium counts as beyond the critical point); the same
bisection refines bifurcations found earlier on the rising path.  ``lambda_limit`` is the maximum load of the
traced path: a limit point, or the bifurcation load when the path turns down
onto a secondary branch there (``limit_type`` says which).

Comparison paths
----------------
* Linear (first order):           u(lambda) = lambda K^-1 F_ref
* Linearised second order
  (P-Delta style):                u(lambda) = lambda (K + lambda K_g,ref)^-1 F_ref
  with ``K_g,ref`` the geometric stiffness of the linear member forces
  under ``F_ref`` (``solver.assemble_geometric_stiffness``).  It becomes
  singular at the linear buckling factor ``lambda_cr`` of
  ``buckling.linear_buckling`` (McGuire et al. 2000, Ch. 9), so it is
  traced up to just below ``lambda_cr``.

Energy check: for this conservative system the work of the external load
along the traced path (trapezoidal rule, O(h^2)) equals the strain energy
``sum EA (l - L)^2 / (2L) + sum k u^2 / 2``; ``energy_error`` reports the
mismatch at the last point.

Scope and limits
----------------
* Self-straining loads T (temperature, settlement) are NOT part of the
  nonlinear path: ``run_nonlinear`` uses the D/S/L/W part of the combination
  only and says so in ``messages``.  The element routines accept a constant
  initial strain ``eps0``; a path started from such a state first solves the
  ``lambda = 0`` equilibrium.
* Elastic material, no member (Euler) buckling or yielding inside the truss
  model.  ``run_nonlinear`` therefore also reports the load factors at which
  the first strut reaches its Euler load and its AISC design strength along
  the path: the global limit load governs only when it is below those.
* The load is proportional (``lambda`` multiplies the whole combination,
  dead load included) and conservative (fixed direction).

References
----------
Crisfield, M. A. (1981). A fast incremental/iterative solution procedure that
handles "snap-through". Computers & Structures 13(1-3), 55-62.
Crisfield, M. A. (1991). Non-linear Finite Element Analysis of Solids and
Structures, Vol. 1. Wiley (Ch. 3 truss strain measures, Ch. 9 arc-length).
Riks, E. (1979). An incremental approach to the solution of snapping and
buckling problems. Int. J. Solids Struct. 15(7), 529-551.
Hellweg, H.-B., Crisfield, M. A. (1998). A new arc-length method for handling
sharp snap-backs. Computers & Structures 66(5), 704-709.
McGuire, W., Gallagher, R. H., Ziemian, R. D. (2000). Matrix Structural
Analysis, 2nd ed. Wiley (Ch. 9 geometric stiffness, Ch. 12 solution of
nonlinear equilibrium equations).
Bathe, K. J. (2014). Finite Element Procedures, 2nd ed.
Yaw, L. L. (2011). 3D Co-rotational Truss Formulation. Walla Walla Univ.
Pecknold, D. A., Ghaboussi, J., Healey, T. J. (1985). Snap-through and
bifurcation in a simple structure. J. Eng. Mech. 111(7), 909-922
(von Mises truss closed form).
Bruno, L., Sassone, M., Venuti, F. (2016). Effects of the equivalent
geometric nodal imperfections on the stability of single layer grid shells.
Eng. Struct. 112, 184-199 (mode-affine imperfections, L/1000 .. L/300).
"""

from __future__ import annotations

import copy
import dataclasses
import math
import time
from dataclasses import dataclass, field

import numpy as np

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False

from . import solver as slv

_AXIS = ("ux", "uy", "uz")


# ============================================================ the system ===
@dataclass
class TrussState:
    """Corotational state of every bar for one displacement vector."""

    u: np.ndarray            # (3n,) displacements (zeros at restrained DOFs)
    l: np.ndarray            # (m,) current lengths
    c: np.ndarray            # (m, 3) current unit vectors
    eps: np.ndarray          # (m,) engineering strain (l - L) / L
    N: np.ndarray            # (m,) axial force, tension +
    f_full: np.ndarray       # (3n,) internal forces incl. springs, all DOFs
    f_free: np.ndarray       # (nf,) internal forces on the free DOFs


class TrussSystem:
    """Pin-jointed truss with restrained DOFs and grounded linear springs.

    Parameters
    ----------
    nodes : (n, 3) original coordinates [m]
    members : (m, 2) node indices
    EA : (m,) axial rigidity [N]
    fixed : restrained DOF numbers (3*node + axis)
    springs : {dof: k} grounded springs [N/m] (ignored on restrained DOFs)
    eps0 : scalar or (m,) constant initial strain (default none)
    """

    def __init__(self, nodes, members, EA, fixed, springs=None, eps0=None):
        self.X = np.array(nodes, dtype=float)
        self.members = np.array(members, dtype=int)
        m = len(self.members)
        self.EA = np.broadcast_to(np.asarray(EA, dtype=float), (m,)).copy()
        self.n_nodes = len(self.X)
        self.ndof = 3 * self.n_nodes
        self.d0 = self.X[self.members[:, 1]] - self.X[self.members[:, 0]]
        self.L = np.linalg.norm(self.d0, axis=1)
        if np.any(self.L <= 0.0):
            raise ValueError(f"Member {int(np.argmin(self.L))} has zero length.")
        self.eps0 = None
        if eps0 is not None:
            e0 = np.broadcast_to(np.asarray(eps0, dtype=float), (m,)).copy()
            if np.any(e0 != 0.0):
                self.eps0 = e0

        self.fixed = np.unique(np.asarray(fixed, dtype=int))
        self.free = slv.free_dofs(self.ndof, self.fixed)
        self.nf = len(self.free)
        if self.nf == 0:
            raise ValueError("Every degree of freedom is restrained.")
        self.fmap = -np.ones(self.ndof, dtype=int)
        self.fmap[self.free] = np.arange(self.nf)

        fixed_set = set(self.fixed.tolist())
        spr = {}
        for dof, k in (springs or {}).items():
            dof, k = int(dof), float(k)
            if dof in fixed_set:
                continue
            if not (np.isfinite(k) and k > 0.0):
                raise ValueError(f"Spring at DOF {dof} must have k > 0.")
            spr[dof] = spr.get(dof, 0.0) + k
        self.spr_dofs = np.asarray(sorted(spr), dtype=int)
        self.spr_k = np.asarray([spr[d] for d in self.spr_dofs], dtype=float)
        self.spr_free = self.fmap[self.spr_dofs]

        # sparse pattern of the free-free tangent
        dofs = slv.element_dofs(self.members)               # (m, 6)
        self.dofs = dofs
        rows = np.repeat(dofs, 6, axis=1).ravel()
        cols = np.tile(dofs, (1, 6)).ravel()
        fr, fc = self.fmap[rows], self.fmap[cols]
        self._keep = (fr >= 0) & (fc >= 0)
        self._rows = fr[self._keep]
        self._cols = fc[self._keep]

    # ------------------------------------------------------------ helpers --
    def expand(self, uf):
        """Free-DOF vector -> full (3n,) vector with zeros at restraints."""
        u = np.zeros(self.ndof)
        u[self.free] = uf
        return u

    def length_scale(self):
        return float(max(np.ptp(self.X, axis=0).max(), self.L.max()))

    # -------------------------------------------------------------- state --
    def state(self, u):
        """Corotational member state and internal forces for displacement u."""
        u = np.asarray(u, dtype=float)
        mem = self.members
        du = u.reshape(-1, 3)
        delta = du[mem[:, 1]] - du[mem[:, 0]]
        d = self.d0 + delta
        l = np.linalg.norm(d, axis=1)
        if np.any(l <= 1e-12 * self.L):
            raise FloatingPointError("A bar has collapsed to zero length.")
        # l^2 - L^2 = 2 d0.delta + delta.delta  (no cancellation for small u)
        dl2 = 2.0 * np.einsum("ij,ij->i", self.d0, delta) + \
            np.einsum("ij,ij->i", delta, delta)
        e = dl2 / (l + self.L)
        eps = e / self.L
        N = self.EA * (eps - self.eps0) if self.eps0 is not None \
            else self.EA * eps
        c = d / l[:, None]
        fv = N[:, None] * c                                   # (m, 3)
        f = np.bincount(self.dofs[:, 3:6].ravel(), weights=fv.ravel(),
                        minlength=self.ndof)
        f -= np.bincount(self.dofs[:, 0:3].ravel(), weights=fv.ravel(),
                         minlength=self.ndof)
        if len(self.spr_dofs):
            f[self.spr_dofs] += self.spr_k * u[self.spr_dofs]
        return TrussState(u=u, l=l, c=c, eps=eps, N=N, f_full=f,
                          f_free=f[self.free])

    def element_tangents(self, st):
        """(m, 6, 6) consistent element tangents k_T for state ``st``."""
        c = st.c
        cc = c[:, :, None] * c[:, None, :]
        I3 = np.eye(3)[None, :, :]
        k3 = (self.EA / self.L)[:, None, None] * cc + \
            (st.N / st.l)[:, None, None] * (I3 - cc)
        ke = np.empty((len(c), 6, 6))
        ke[:, 0:3, 0:3] = k3
        ke[:, 0:3, 3:6] = -k3
        ke[:, 3:6, 0:3] = -k3
        ke[:, 3:6, 3:6] = k3
        return ke

    def tangent(self, st):
        """Free-free tangent stiffness K_T (sparse CSC, or dense ndarray)."""
        vals = self.element_tangents(st).reshape(-1)[self._keep]
        if _HAVE_SCIPY:
            K = sp.coo_matrix((vals, (self._rows, self._cols)),
                              shape=(self.nf, self.nf)).tocsc()
            if len(self.spr_dofs):
                K = K + sp.coo_matrix(
                    (self.spr_k, (self.spr_free, self.spr_free)),
                    shape=(self.nf, self.nf)).tocsc()
            return K.tocsc()
        K = np.zeros((self.nf, self.nf))                    # pragma: no cover
        np.add.at(K, (self._rows, self._cols), vals)
        K[self.spr_free, self.spr_free] += self.spr_k
        return K

    def tangent_full(self, u):
        """Full (3n x 3n) dense tangent incl. springs (for checks/tests)."""
        st = self.state(u)
        K = np.zeros((self.ndof, self.ndof))
        ke = self.element_tangents(st)
        for m, d in enumerate(self.dofs):
            K[np.ix_(d, d)] += ke[m]
        K[self.spr_dofs, self.spr_dofs] += self.spr_k
        return K

    def strain_energy(self, st):
        """sum EA L (eps - eps0)^2 / 2 + sum k u^2 / 2   [J]."""
        e = st.eps - self.eps0 if self.eps0 is not None else st.eps
        U = 0.5 * float(np.sum(self.EA * self.L * e ** 2))
        if len(self.spr_dofs):
            U += 0.5 * float(np.sum(self.spr_k * st.u[self.spr_dofs] ** 2))
        return U


# ===================================================== linear algebra ====
def _perm_sign(p):
    """Sign of a permutation given as an index array."""
    p = np.asarray(p, dtype=int)
    seen = np.zeros(len(p), dtype=bool)
    sign = 1
    for i in range(len(p)):
        if seen[i]:
            continue
        j, n = i, 0
        while not seen[j]:
            seen[j] = True
            j = p[j]
            n += 1
        if n % 2 == 0:
            sign = -sign
    return sign


class _Factor:
    """LU factorisation of K_T with solve() and the sign of det K_T."""

    def __init__(self, K):
        self.sparse = _HAVE_SCIPY and sp.issparse(K)
        if self.sparse:
            self.lu = spla.splu(K.tocsc())      # raises RuntimeError if singular
        else:                                            # pragma: no cover
            self.K = np.asarray(K)

    def solve(self, b):
        x = self.lu.solve(b) if self.sparse else np.linalg.solve(self.K, b)
        if not np.all(np.isfinite(x)):
            raise FloatingPointError("non-finite solution")
        return x

    def det_sign(self):
        if not self.sparse:                              # pragma: no cover
            return int(np.sign(np.linalg.slogdet(self.K)[0]))
        d = self.lu.U.diagonal()
        if np.any(d == 0.0):
            return 0
        s = -1 if (np.sum(d < 0.0) % 2) else 1
        return s * _perm_sign(self.lu.perm_r) * _perm_sign(self.lu.perm_c)


def _factor(K):
    try:
        return _Factor(K)
    except (RuntimeError, np.linalg.LinAlgError) as exc:
        raise FloatingPointError(f"singular tangent ({exc})") from None


# ============================================================ results =====
@dataclass
class NonlinearResult:
    """Equilibrium path and diagnostics (plain numpy arrays, plot-ready).

    Arrays indexed by path point k = 0 .. n_points-1 (k = 0 is the unloaded
    state).  ``control`` is the raw displacement of the control DOF [m]
    (negative = downward for ``uz``).
    """

    method: str
    control_dof: int
    control_label: str
    lam: np.ndarray
    control: np.ndarray
    u: np.ndarray                    # (k, 3n) displacement snapshots
    N: np.ndarray                    # (k, m) member forces
    iterations: np.ndarray
    residual: np.ndarray
    det_sign: np.ndarray
    n_negative: np.ndarray
    min_eig: np.ndarray
    stiffness_parameter: np.ndarray
    strain_energy: np.ndarray
    external_work: np.ndarray
    F_ref: np.ndarray
    combination: str = ""

    limit_found: bool = False
    limit_refined: bool = False
    limit_type: str = ""             # "limit" | "bifurcation" | ""
    limit_stiffness_parameter: float = float("nan")
    limit_load_projection: float = float("nan")   # |phi.F|/|F| at the peak
    lambda_limit: float = float("nan")
    control_limit: float = float("nan")
    u_limit: np.ndarray = None
    N_limit: np.ndarray = None
    limit_step: int = -1
    critical_points: list = field(default_factory=list)

    lambda_cr: float = float("nan")
    buckling_mode: np.ndarray = None
    linear: dict = field(default_factory=dict)      # {"lam", "control"}
    p_delta: dict = field(default_factory=dict)     # {"lam", "control"}

    lambda_member_euler: float = float("nan")
    member_euler: int = -1
    lambda_member_design: float = float("nan")
    member_design: int = -1

    stop_reason: str = ""
    messages: list = field(default_factory=list)
    elapsed: float = 0.0

    # ------------------------------------------------------------ queries --
    @property
    def n_points(self):
        return len(self.lam)

    @property
    def limit_ratio(self):
        """lambda_limit / lambda_cr (NaN when either is unknown)."""
        if self.limit_found and np.isfinite(self.lambda_cr) and \
                self.lambda_cr > 0.0:
            return self.lambda_limit / self.lambda_cr
        return float("nan")

    @property
    def energy_error(self):
        if len(self.strain_energy) < 2:
            return 0.0
        scale = max(float(np.max(np.abs(self.strain_energy))), 1e-300)
        return abs(float(self.external_work[-1] - self.strain_energy[-1])) \
            / scale

    @property
    def lambda_first_critical(self):
        """Lowest critical load factor on the rising path: the first load
        maximum or an earlier bifurcation (NaN when neither was found)."""
        vals = [cp["lambda"] for cp in self.critical_points
                if cp["type"] == "bifurcation" and not cp.get("first_peak")]
        if self.limit_found:
            vals.append(self.lambda_limit)
        return float(min(vals)) if vals else float("nan")

    @property
    def max_residual(self):
        return float(np.max(self.residual)) if len(self.residual) else 0.0

    def curves(self):
        """{name: (control [m], lambda)} for the comparison plot."""
        out = {"corotational": (self.control.copy(), self.lam.copy())}
        if self.linear:
            out["linear"] = (np.asarray(self.linear["control"]),
                             np.asarray(self.linear["lam"]))
        if self.p_delta:
            out["p_delta"] = (np.asarray(self.p_delta["control"]),
                              np.asarray(self.p_delta["lam"]))
        return out

    def summary_lines(self):
        L = [f"Method                 {self.method}",
             f"Combination            {self.combination or '-'}",
             f"Control DOF            {self.control_label}",
             f"Path points            {self.n_points - 1} steps, "
             f"{int(np.sum(self.iterations))} iterations, "
             f"{self.elapsed:.2f} s",
             f"Max residual           {self.max_residual:.1e}",
             f"Path energy mismatch   {self.energy_error:.1e} "
             f"(trapezoidal work vs strain energy)"]
        if self.limit_found:
            kind = "bifurcation" if self.limit_type == "bifurcation" \
                else "limit point"
            L.append(f"First load maximum     lambda = {self.lambda_limit:.6g} "
                     f"({kind}) at {self.control_label} = "
                     f"{self.control_limit * 1e3:.3f} mm"
                     + ("" if self.limit_refined else " (not refined)"))
            if np.isfinite(self.limit_load_projection):
                L.append(f"Critical mode          load projection "
                         f"|phi.F|/|F| = {self.limit_load_projection:.2e}")
        else:
            L.append("First limit point      not reached "
                     f"({self.stop_reason})")
        for cp in self.critical_points:
            if cp["type"] == "bifurcation" and not cp.get("first_peak"):
                L.append(f"Bifurcation            lambda = {cp['lambda']:.6g} "
                         f"(stiffness sign change on the rising path)")
        if np.isfinite(self.lambda_cr):
            L.append(f"Linear buckling        lambda_cr = {self.lambda_cr:.6g}")
            if self.limit_found:
                L.append(f"lambda_limit / lambda_cr  {self.limit_ratio:.4f}")
        if np.isfinite(self.lambda_member_euler):
            L.append(f"First strut at Euler   lambda = "
                     f"{self.lambda_member_euler:.4g} (member "
                     f"{self.member_euler})")
        if np.isfinite(self.lambda_member_design):
            L.append(f"First strut at phiPn   lambda = "
                     f"{self.lambda_member_design:.4g} (member "
                     f"{self.member_design})")
        return L

    def to_dict(self, snapshots=False):
        """JSON-friendly dict (lists; NaN -> None)."""
        def cl(x):
            if x is None:
                return None
            if isinstance(x, np.ndarray):
                return [cl(v) for v in x.tolist()]
            if isinstance(x, (list, tuple)):
                return [cl(v) for v in x]
            if isinstance(x, dict):
                return {str(k): cl(v) for k, v in x.items()}
            if isinstance(x, (np.integer,)):
                return int(x)
            if isinstance(x, (float, np.floating)):
                x = float(x)
                return x if math.isfinite(x) else None
            if isinstance(x, np.bool_):
                return bool(x)
            return x
        d = {}
        for f in dataclasses.fields(self):
            if f.name in ("u", "N") and not snapshots:
                continue
            d[f.name] = cl(getattr(self, f.name))
        d["limit_ratio"] = cl(self.limit_ratio)
        d["energy_error"] = cl(self.energy_error)
        return d


# ===================================================== path tracing =======
_METHODS = {"arc-length": "arc-length", "arclength": "arc-length",
            "arc_length": "arc-length", "crisfield": "arc-length",
            "displacement": "displacement", "displacement-control":
            "displacement", "displacement_control": "displacement",
            "load": "load", "load-control": "load", "load_control": "load",
            "newton": "load", "newton-raphson": "load"}


class _Fail(Exception):
    pass


class _Tracer:
    """Incremental-iterative path following on a TrussSystem."""

    def __init__(self, system, F_ref, control_dof, tol, tol_disp, max_iter,
                 eig_limit):
        self.s = system
        self.F = np.asarray(F_ref, dtype=float)
        self.Ff = self.F[system.free]
        self.nF = float(np.linalg.norm(self.Ff))
        if self.nF <= 0.0:
            raise ValueError("The reference load is zero on every free DOF.")
        self.c = int(control_dof)
        self.ic = int(system.fmap[self.c])
        if self.ic < 0:
            raise ValueError(f"Control DOF {self.c} is restrained.")
        self.tol = float(tol)
        self.tol_disp = float(tol_disp)
        self.max_iter = int(max_iter)
        self.eig_limit = int(eig_limit)
        self.lam_ref = 0.0
        self.Lc = system.length_scale()
        self.du_ref = 0.0          # typical step norm (set by trace_path)

    # ------------------------------------------------------- convergence --
    def _rnorm(self, st, lam):
        r = lam * self.Ff - st.f_free
        scale = max(abs(lam) * self.nF, self.lam_ref * self.nF,
                    float(np.linalg.norm(st.f_full)), 1e-300)
        return r, float(np.linalg.norm(r)) / scale

    def _converged(self, rn, corr, Dnorm, tol, force_only=False):
        if not np.isfinite(rn):
            raise _Fail("non-finite residual")
        if rn > tol:
            return False
        if corr is None:                 # predictor only: accept if exact
            return rn <= 1e-3 * tol
        if force_only:
            # Next to a critical point the tangent is nearly singular and
            # turns round-off into corrections that never pass the
            # displacement test, although equilibrium holds to round-off.
            return True
        return corr <= self.tol_disp * max(Dnorm, 1e-3 * self.du_ref,
                                           1e-14 * self.Lc)

    # -------------------------------------------------------- correctors --
    def correct(self, kind, u0, lam0, Du, DL, param, tol=None,
                force_only=False):
        """Iterate from predictor (Du, DL) to an equilibrium point.

        kind 'arc': |Du| = param; 'disp': Du[ic] = param; 'load': DL fixed.
        ``force_only`` accepts a point once the force residual has converged
        (used by the critical-point bisection).
        Returns (u_full, lam, state, iterations, rnorm).
        """
        tol = self.tol if tol is None else tol
        corr = None
        Du = np.array(Du, dtype=float)
        DL = float(DL)
        if kind == "disp":
            Du[self.ic] = param
        for it in range(self.max_iter + 1):
            try:
                st = self.s.state(u0 + self.s.expand(Du))
            except FloatingPointError as exc:
                raise _Fail(str(exc)) from None
            lam = lam0 + DL
            r, rn = self._rnorm(st, lam)
            if self._converged(rn, corr, float(np.linalg.norm(Du)), tol,
                               force_only):
                return st.u, lam, st, it, rn
            if it == self.max_iter:
                break
            try:
                fac = _factor(self.s.tangent(st))
                vr = fac.solve(r)
                vt = fac.solve(self.Ff)
            except FloatingPointError as exc:
                raise _Fail(str(exc)) from None
            if kind == "load":
                dl = 0.0
            elif kind == "disp":
                if abs(vt[self.ic]) < 1e-300:
                    raise _Fail("control DOF insensitive to the load")
                dl = -vr[self.ic] / vt[self.ic]
            else:
                w = Du + vr
                a = float(vt @ vt)
                b = 2.0 * float(vt @ w)
                cq = float(w @ w) - param ** 2
                disc = b * b - 4.0 * a * cq
                if a <= 0.0 or disc < 0.0:
                    raise _Fail("complex arc-length roots")
                sq = math.sqrt(disc)
                q = -0.5 * (b + math.copysign(sq, b))
                roots = [q / a] + ([cq / q] if q != 0.0 else [])
                best, best_cos = None, -np.inf
                nD = float(np.linalg.norm(Du))
                for root in roots:
                    Dn = w + root * vt
                    cs = float(Du @ Dn) / max(nD * float(np.linalg.norm(Dn)),
                                              1e-300)
                    if cs > best_cos:
                        best, best_cos = root, cs
                dl = best
            step = vr + dl * vt
            Du = Du + step
            DL += dl
            corr = float(np.linalg.norm(step))
        raise _Fail(f"no convergence in {self.max_iter} iterations "
                    f"(residual {rn:.1e})")

    # ----------------------------------------------------- point analysis --
    def analyse_point(self, u):
        """Tangent factor, det sign, eigen indicator, v_t at a converged u."""
        st = self.s.state(u)
        K = self.s.tangent(st)
        fac = _factor(K)
        v = fac.solve(self.Ff)
        neg, mine = -1, float("nan")
        if self.s.nf <= self.eig_limit:
            Kd = K.toarray() if (_HAVE_SCIPY and sp.issparse(K)) else K
            ev = np.linalg.eigvalsh(0.5 * (Kd + Kd.T))
            mine = float(ev[0])
            neg = int(np.sum(ev < -1e-12 * max(abs(ev[-1]), 1e-300)))
        return st, fac, v, fac.det_sign(), neg, mine


def _interp_predictor(ts, us, lams, t):
    """Quadratic (Lagrange) interpolation of path samples in parameter t."""
    ts = np.asarray(ts, dtype=float)
    k = len(ts)
    w = np.ones(k)
    for i in range(k):
        for j in range(k):
            if i != j:
                den = ts[i] - ts[j]
                w[i] *= (t - ts[j]) / den if den != 0.0 else 0.0
    u = sum(w[i] * us[i] for i in range(k))
    lam = float(sum(w[i] * lams[i] for i in range(k)))
    return u, lam


def _golden_max(f, a, b, xtol, maxit=200):
    """Maximise a unimodal f on [a, b]; Brent via SciPy, golden otherwise."""
    try:
        from scipy.optimize import minimize_scalar
        res = minimize_scalar(lambda x: -f(x), bounds=(a, b),
                              method="bounded",
                              options={"xatol": xtol, "maxiter": maxit})
        return float(res.x)
    except ImportError:                                  # pragma: no cover
        g = (math.sqrt(5.0) - 1.0) / 2.0
        c, d = b - g * (b - a), a + g * (b - a)
        fc, fd = f(c), f(d)
        for _ in range(maxit):
            if abs(b - a) < xtol:
                break
            if fc > fd:
                b, d, fd = d, c, fc
                c = b - g * (b - a)
                fc = f(c)
            else:
                a, c, fc = c, d, fd
                d = a + g * (b - a)
                fd = f(d)
        return 0.5 * (a + b)


def _resolve_method(method):
    key = str(method).strip().lower()
    if key not in _METHODS:
        raise ValueError(f"Unknown method {method!r}: use 'arc-length', "
                         f"'displacement' or 'load'.")
    return _METHODS[key]


def trace_path(system, F_ref, method="arc-length", control_dof=None, *,
               first_increment=0.1, arc_length=None, max_arc_factor=5.0,
               first_arc_length=None,
               lambda_max=np.inf, lambda_min=-np.inf, control_max=np.inf,
               stop_after_limit=5, max_steps=300, tol=1e-8, tol_disp=1e-6,
               max_iter=30, desired_iterations=5, max_cutbacks=12,
               refine_limit=True, refine_bifurcation=True, eig_limit=600,
               max_turn=0.5, turn_cutbacks=4, progress=None, label=""):
    """Trace the equilibrium path of ``system`` under ``lambda * F_ref``.

    Parameters
    ----------
    method : 'arc-length' (Crisfield cylindrical), 'displacement' (control
        DOF), or 'load' (Newton-Raphson, load control).
    control_dof : global DOF recorded as the path abscissa (and controlled
        by the displacement method).  Default: the free DOF with the largest
        linear displacement.
    first_increment : load-factor size of the first step; it fixes the arc
        length ``Dl = first_increment * |K0^-1 F_ref|`` (arc-length), the
        control increment (displacement) or the load step (load).
    arc_length : explicit initial arc length [m] (overrides first_increment).
    first_arc_length : arc-length method only: cap [m] on the FIRST step
        alone.  The step then grows (by the iteration-count rule) up to the
        usual maximum, ``max_arc_factor`` times the nominal arc length.  A
        small first step resolves a limit point that lies very close to the
        start of the path (an imperfect node that is almost flat with its
        neighbours), which a nominal step would jump over.
    max_arc_factor : arc length (or control/load step) is kept within
        [initial * 1e-6, initial * max_arc_factor].
    lambda_max, lambda_min, control_max : stop when lambda leaves
        [lambda_min, lambda_max] or |u_control| exceeds control_max.
    stop_after_limit : number of steps traced after the first limit point
        before stopping (None: keep going).
    tol, tol_disp : force and displacement-increment tolerances (see the
        module docstring).
    eig_limit : dense eigenvalue indicator up to this many free DOFs.
    max_turn, turn_cutbacks : step-acceptance check against jumping over a
        critical point or onto a remote branch: a converged increment ``Du``
        with ``|Du - Du_pred| > max_turn |Du_pred|`` is rejected and the step
        halved, up to ``turn_cutbacks`` times (after that it is accepted, so
        a genuine corner such as a branch switch cannot stall the trace).
        Load control always uses ``max_turn = 1`` and no acceptance limit.
    progress : optional callable(fraction, message).

    Returns a :class:`NonlinearResult` (without lambda_cr / comparison paths).
    """
    t0 = time.perf_counter()
    method = _resolve_method(method)
    s = system
    F = np.asarray(F_ref, dtype=float)
    if F.shape != (s.ndof,):
        raise ValueError(f"F_ref must have length {s.ndof}.")

    # initial tangent and linear response
    st0 = s.state(np.zeros(s.ndof))
    fac0 = _factor(s.tangent(st0))
    v0 = fac0.solve(F[s.free])
    if control_dof is None:
        control_dof = int(s.free[int(np.argmax(np.abs(v0)))])
    tr = _Tracer(s, F, control_dof, tol, tol_disp, max_iter, eig_limit)
    msgs = []
    if s.nf > int(eig_limit):
        msgs.append(
            f"This model has {s.nf} free DOF, more than eig_limit = "
            f"{int(eig_limit)}, so stability is checked by the sign of "
            f"det K_T only.  A pair of equal eigenvalues crossing zero "
            f"together (a double bifurcation, common in symmetric domes) "
            f"leaves that sign unchanged and is missed.  Compare "
            f"lambda_limit with lambda_cr.")

    u = np.zeros(s.ndof)
    lam = 0.0
    if s.eps0 is not None:                     # self-strained initial state
        try:
            u, lam, _, _, _ = tr.correct("load", u, 0.0, np.zeros(s.nf),
                                         0.0, None)
        except _Fail as exc:
            raise RuntimeError(f"Initial self-strained state did not "
                               f"converge: {exc}") from None
        msgs.append("Initial-strain state solved at lambda = 0.")
    st, fac, v, dsg, neg, mine = tr.analyse_point(u)
    FtK0F = float(tr.Ff @ v)

    P = {k: [] for k in ("lam", "control", "u", "N", "it", "res", "det",
                         "neg", "mine", "sp", "U", "W")}

    def record(st, lam, it, res, dsg, neg, mine, v):
        P["lam"].append(lam)
        P["control"].append(float(st.u[tr.c]))
        P["u"].append(st.u.copy())
        P["N"].append(st.N.copy())
        P["it"].append(it)
        P["res"].append(res)
        P["det"].append(dsg)
        P["neg"].append(neg)
        P["mine"].append(mine)
        ftv = float(tr.Ff @ v)
        P["sp"].append(FtK0F / ftv if ftv != 0.0 else float("inf"))
        P["U"].append(s.strain_energy(st))
        if P["W"]:
            du = st.u - P["u"][-2]
            P["W"].append(P["W"][-1] + 0.5 * (P["lam"][-2] + lam) *
                          float(F @ du))
        else:
            P["W"].append(P["U"][-1])

    record(st, lam, 0, 0.0, dsg, neg, mine, v)

    # step size
    nv0 = float(np.linalg.norm(v0))
    if method == "arc-length":
        size0 = float(arc_length) if arc_length else \
            abs(first_increment) * nv0
    elif method == "displacement":
        size0 = abs(first_increment * v0[tr.ic])
        if size0 <= 0.0:
            raise ValueError("The control DOF does not move under F_ref; "
                             "choose another control DOF.")
    else:
        size0 = abs(first_increment)
    size, size_min, size_max = size0, size0 * 1e-6, size0 * max_arc_factor
    if method == "arc-length" and first_arc_length and \
            0.0 < float(first_arc_length) < size0:
        size = float(first_arc_length)
        size_min = size * 1e-6
    tr.du_ref = abs(first_increment) * nv0 if method != "arc-length" \
        else size0
    dsign = 1.0
    if method == "displacement":
        dsign = math.copysign(1.0, v0[tr.ic] * first_increment)
    Du_prev = None
    DL_prev = None
    limit_k = None
    steps_after = 0
    stop = "maximum number of steps reached"
    cut_total = 0

    for step in range(1, max_steps + 1):
        cut = 0
        while True:
            try:
                if method == "arc-length":
                    sgn = 1.0
                    if Du_prev is not None:
                        dot = float(Du_prev @ v)
                        sgn = math.copysign(1.0, dot) if dot != 0.0 else \
                            math.copysign(1.0, DL_prev)
                    DL = sgn * size / max(float(np.linalg.norm(v)), 1e-300)
                    Du = DL * v
                    un, lamn, stn, it, rn = tr.correct("arc", u, lam, Du, DL,
                                                       size)
                    if Du_prev is None and not lamn > lam:
                        # A first step that is too long can converge on
                        # the branch below lambda = 0 (the path then runs
                        # backwards); cut it back.
                        raise _Fail("the first step went to a lower load")
                elif method == "displacement":
                    dc = dsign * size
                    if abs(v[tr.ic]) < 1e-300:
                        raise _Fail("control DOF insensitive to the load")
                    DL = dc / v[tr.ic]
                    Du = DL * v
                    un, lamn, stn, it, rn = tr.correct("disp", u, lam, Du, DL,
                                                       dc)
                else:
                    DL = size
                    Du = DL * v
                    un, lamn, stn, it, rn = tr.correct("load", u, lam, Du, DL,
                                                       None)
                    # Newton under load control can converge to a remote
                    # equilibrium (the far side of a snap-through) without
                    # ever crossing the limit point.  Accept a step only
                    # when its correction is smaller than its predictor.
                    if np.linalg.norm((un - u)[s.free] - Du) > \
                            np.linalg.norm(Du):
                        raise _Fail("load step jumped to a remote "
                                    "equilibrium (snap-through)")
                if method != "load" and cut < turn_cutbacks:
                    dev = np.linalg.norm((un - u)[s.free] - Du)
                    if dev > max_turn * np.linalg.norm(Du):
                        raise _Fail("converged increment turned away from "
                                    "the predictor")
                stn, facn, vn, dsgn, negn, minen = tr.analyse_point(un)
                if method == "load" and dsgn != dsg:
                    raise _Fail("stiffness sign changed in a load step")
                break
            except (_Fail, FloatingPointError) as exc:
                cut += 1
                cut_total += 1
                size *= 0.5
                if cut > max_cutbacks or size < size_min:
                    stop = (f"step {step} failed after {cut - 1} "
                            f"cut-back(s), step size {size * 2:.3g} ({exc})")
                    if method == "load" and limit_k is None:
                        stop = ("load control cannot pass the limit point; "
                                f"step {step} failed at lambda = {lam:.6g}")
                    size = None
                    break
        if size is None:
            break

        DL_step = lamn - lam
        Du_prev = (un - u)[s.free]
        DL_prev = DL_step if DL_step != 0.0 else (DL_prev or 1.0)
        u, lam, v, dsg = un, lamn, vn, dsgn
        tr.lam_ref = max(tr.lam_ref, abs(lam))
        record(stn, lam, it, rn, dsgn, negn, minen, vn)

        # adapt the step
        if cut == 0:
            fct = math.sqrt(desired_iterations / max(it, 1))
            size = min(max(size * min(max(fct, 0.25), 2.0), size_min),
                       size_max)

        # first limit point: lambda stops rising
        k = len(P["lam"]) - 1
        if limit_k is None and k >= 2 and P["lam"][k] < P["lam"][k - 1] \
                and P["lam"][k - 1] > 0.0 and \
                P["lam"][k - 1] >= P["lam"][k - 2]:
            limit_k = k - 1
        elif limit_k is not None:
            steps_after += 1

        if progress is not None:
            progress(min(0.95, step / max_steps),
                     f"{label}step {step}: lambda = {lam:.4g}")

        if lam > lambda_max:
            stop = f"lambda exceeded lambda_max = {lambda_max:g}"
            break
        if lam < lambda_min:
            stop = f"lambda fell below lambda_min = {lambda_min:g}"
            break
        if abs(P["control"][-1]) > control_max:
            stop = f"|control displacement| exceeded {control_max:g} m"
            break
        if limit_k is not None and stop_after_limit is not None and \
                steps_after >= stop_after_limit:
            stop = (f"first limit point passed ({stop_after_limit} steps "
                    f"traced beyond it)")
            break
    else:
        stop = f"maximum number of steps ({max_steps}) reached"

    arr = {k: np.asarray(v_) for k, v_ in P.items()}
    res = NonlinearResult(
        method=method, control_dof=tr.c, control_label=_dof_label(tr.c),
        lam=arr["lam"].astype(float), control=arr["control"].astype(float),
        u=np.vstack(P["u"]), N=np.vstack(P["N"]),
        iterations=arr["it"].astype(int), residual=arr["res"].astype(float),
        det_sign=arr["det"].astype(int), n_negative=arr["neg"].astype(int),
        min_eig=arr["mine"].astype(float),
        stiffness_parameter=arr["sp"].astype(float),
        strain_energy=arr["U"].astype(float),
        external_work=arr["W"].astype(float), F_ref=F.copy(),
        stop_reason=stop, messages=msgs)
    if cut_total:
        msgs.append(f"{cut_total} step cut-back(s) during the trace.")

    # ---- critical points ---------------------------------------------------
    tr.FtK0F = FtK0F
    lamv = res.lam
    key = _crit_key(res.det_sign, res.n_negative)
    if limit_k is not None:
        res.limit_found = True
        res.limit_step = limit_k
        pk = _analyse_peak(tr, method, res, limit_k, key, refine_limit, msgs)
        res.limit_type = pk["type"]
        res.limit_refined = pk["refined"]
        res.lambda_limit = float(pk["lambda"])
        res.u_limit = np.array(pk["u"], dtype=float)
        res.N_limit = s.state(res.u_limit).N.copy()
        res.control_limit = float(res.u_limit[tr.c])
        res.limit_stiffness_parameter = float(pk["sp"])
        res.limit_load_projection = float(pk["rho"])
        cp = {"type": pk["type"], "step": limit_k,
              "lambda": res.lambda_limit, "control": res.control_limit,
              "refined": pk["refined"], "stiffness_parameter": pk["sp"],
              "load_projection": pk["rho"], "first_peak": True}
        if "lambda_bracket" in pk:
            cp["lambda_bracket"] = pk["lambda_bracket"]
        res.critical_points.append(cp)
        if pk["type"] == "bifurcation":
            msgs.append(
                f"The load maximum of the traced path is a bifurcation "
                f"(critical-mode load projection {pk['rho']:.1e}): the "
                f"critical mode is orthogonal to the load and the path "
                f"continues on a falling secondary branch.  lambda_limit is "
                f"the bifurcation load of the perfect structure.")
    elif method == "load" and stop.startswith("load control"):
        # The cut-backs bracket the limit: last converged lambda.
        res.limit_found = True
        res.limit_type = "limit"
        res.limit_step = res.n_points - 1
        res.lambda_limit = float(lamv[-1])
        res.control_limit = float(res.control[-1])
        res.u_limit = res.u[-1].copy()
        res.N_limit = res.N[-1].copy()
        res.limit_stiffness_parameter = float(res.stiffness_parameter[-1])
        msgs.append("Load control: the limit load is bracketed by the last "
                    "converged step (within the smallest load step); use "
                    "arc-length to pass it.")
        res.critical_points.append({
            "type": "limit", "step": res.limit_step,
            "lambda": res.lambda_limit, "control": res.control_limit,
            "refined": False, "first_peak": True,
            "stiffness_parameter": res.limit_stiffness_parameter})

    # stiffness changes on the rising path before the first peak
    k_stop = limit_k if limit_k is not None else res.n_points
    for k in range(1, res.n_points):
        if key[k] == key[k - 1]:
            continue
        if limit_k is not None and k in (limit_k, limit_k + 1):
            continue                     # belongs to the first peak
        if lamv[k] > lamv[k - 1] and k < k_stop:
            cp = {"type": "bifurcation", "step": k,
                  "lambda": float(0.5 * (lamv[k - 1] + lamv[k])),
                  "control": float(0.5 * (res.control[k - 1] +
                                          res.control[k])),
                  "refined": False, "first_peak": False}
            if refine_bifurcation:
                try:
                    b = _bisect_critical(tr, res, k, key)
                    cp.update({"lambda": b["lambda"], "refined": True,
                               "control": float(b["u"][tr.c]),
                               "lambda_bracket": b["lambda_bracket"],
                               "stiffness_parameter": b["sp"],
                               "load_projection": b["rho"]})
                except Exception as exc:                  # pragma: no cover
                    msgs.append(f"Bifurcation near step {k} not refined "
                                f"({exc}).")
            res.critical_points.append(cp)
    res.critical_points.sort(key=lambda d: d["step"])
    res.elapsed = time.perf_counter() - t0
    return res


def _dof_label(dof):
    return f"node {dof // 3} {_AXIS[dof % 3]}"


def _point_solver(tr, method, res, base):
    """Equilibrium point at parameter t measured from path point ``base``.

    Arc-length: t = cylindrical distance |u - u_base| (free DOFs).
    Displacement control: t = u_c - u_c,base.
    Predictor: quadratic interpolation of up to three path samples.
    """
    s = tr.s
    ub = res.u[base]
    lb = float(res.lam[base])
    idx = [i for i in (base, base + 1, base + 2) if i < res.n_points]
    if method == "arc-length":
        ts = [float(np.linalg.norm((res.u[i] - ub)[s.free])) for i in idx]
    else:
        ts = [float(res.u[i][tr.c] - ub[tr.c]) for i in idx]
    us = [res.u[i][s.free] - ub[s.free] for i in idx]
    ls = [float(res.lam[i]) - lb for i in idx]
    tol = min(tr.tol, 1e-11)

    def solve(t):
        Du, DL = _interp_predictor(ts, us, ls, t)
        if method == "arc-length":
            return tr.correct("arc", ub, lb, Du, DL, t, tol=tol)
        return tr.correct("disp", ub, lb, Du, DL, t, tol=tol)
    return solve, ts


def _refine_limit(tr, method, res, k):
    """Maximise lambda over converged points between samples k-1 and k+1."""
    base = k - 1
    solve, ts = _point_solver(tr, method, res, base)
    a, b = 0.0, ts[-1]
    cache = {}

    def lam_of(t):
        if t not in cache:
            u, lam, _, _, _ = solve(t)
            cache[t] = (u, lam)
        return cache[t][1]

    lo, hi = (a, b) if b > a else (b, a)
    lam_of(_golden_max(lam_of, lo, hi, xtol=1e-10 * abs(b - a)))
    # Keep the best evaluated point (Brent may stop next to it).
    tb = max(cache, key=lambda t: cache[t][1])
    return cache[tb][0].copy(), float(cache[tb][1])


_RHO_LIMIT = 1e-6     # load projection of the critical mode: limit vs bifurcation


def _crit_key(det_sign, n_negative):
    """Per-point stability key: the negative-eigenvalue count when it is
    known (it also catches a double eigenvalue crossing, which leaves the
    determinant sign unchanged), else the determinant sign."""
    neg = np.asarray(n_negative)
    if len(neg) and np.all(neg >= 0):
        return neg.astype(int)
    return np.asarray(det_sign).astype(int)


def _point_key(neg, dsg):
    return int(neg) if neg >= 0 else int(dsg)


def _sp_of(tr, v):
    ftv = float(tr.Ff @ v)
    return tr.FtK0F / ftv if ftv != 0.0 else 0.0


def critical_mode(tr, u, n_iter=60):
    """Eigenpair of K_T(u) with the smallest |eigenvalue| (inverse iteration
    with the LU factor) and its load projection rho = |phi.F| / (|phi| |F|).

    rho is O(1e-14) when the critical mode is orthogonal to the load
    (bifurcation of a symmetric structure) and finite at a limit point.
    """
    st = tr.s.state(u)
    K = tr.s.tangent(st)
    fac = _factor(K)
    rng = np.random.default_rng(12345)
    x = rng.standard_normal(tr.s.nf)
    x /= np.linalg.norm(x)
    mu_old = None
    for _ in range(n_iter):
        y = fac.solve(x)
        x = y / np.linalg.norm(y)
        mu = float(x @ (K @ x))
        if mu_old is not None and abs(mu - mu_old) <= 1e-12 * max(abs(mu), 1e-300):
            break
        mu_old = mu
    rho = abs(float(x @ tr.Ff)) / tr.nF
    return mu, x, rho


def _analyse_peak(tr, method, res, k, key, refine, msgs):
    """Classify and refine the first maximum of lambda (sample k).

    1. Maximise lambda over a smooth one-parameter family of converged
       points (Brent).  A smooth maximum whose critical mode does work on
       the load (rho >= _RHO_LIMIT) is a limit point.
    2. Otherwise (or if 1 fails) bracket the stability change along the
       rising primary path by load-controlled bisection; the load
       projection rho of the critical mode there decides limit point
       (rho finite) or bifurcation (rho ~ round-off).
    """
    out = {"type": "limit", "u": res.u[k].copy(), "lambda": float(res.lam[k]),
           "refined": False, "sp": float(res.stiffness_parameter[k]),
           "rho": float("nan")}
    if not refine or method == "load":
        return out
    why = ""
    try:
        ul, ll = _refine_limit(tr, method, res, k)
        _, _, v, _, _, _ = tr.analyse_point(ul)
        _, _, rho = critical_mode(tr, ul)
        if ll >= out["lambda"] - 1e-12 * max(abs(ll), 1e-300) and \
                rho >= _RHO_LIMIT:
            return {"type": "limit", "u": ul, "lambda": float(ll),
                    "refined": True, "sp": float(_sp_of(tr, v)),
                    "rho": float(rho)}
        why = f"critical-mode load projection {rho:.1e}"
    except Exception as exc:
        why = str(exc)
    j = next((i for i in (k, k + 1) if i < res.n_points and
              key[i] != key[i - 1]), None)
    if j is None:
        msgs.append(f"First peak not refined ({why}); the largest path "
                    f"sample is reported.")
        return out
    try:
        b = _bisect_critical(tr, res, j, key)
    except Exception as exc:
        msgs.append(f"First peak not refined ({why}; {exc}); the largest "
                    f"path sample is reported.")
        return out
    if b["lambda"] < out["lambda"] - 1e-9 * abs(out["lambda"]):
        msgs.append("Bracketed critical point lies below the largest path "
                    "sample; the sample is reported.")
        return out
    typ = "limit" if b["rho"] >= _RHO_LIMIT else "bifurcation"
    return {"type": typ, "u": b["u"], "lambda": b["lambda"],
            "refined": True, "sp": b["sp"], "rho": b["rho"],
            "lambda_bracket": b["lambda_bracket"]}


def _bisect_critical(tr, res, j, key, rtol=1e-10, max_bisect=80):
    """Load-controlled bisection for the stability change at sample j.

    Starts at the last point with the old stability key (j-1), climbs the
    primary path with load-controlled Newton steps (doubling) until the key
    changes or Newton fails (beyond a limit load), then bisects the load
    factor.  A point counts as stable when its force residual converges
    (``force_only``: the displacement-increment test is not used, because
    the nearly singular tangent next to the critical point amplifies
    round-off) and its stability key is still the old one.  The bracket
    closes to ``rtol``; the bracketed load is accurate to about 1e-9 next
    to a double bifurcation.  Returns {"lambda", "lambda_bracket", "u"
    (last stable point), "sp" (stiffness parameter there)}.
    """
    key0 = int(key[j - 1])
    u_lo = res.u[j - 1].copy()
    l_lo = float(res.lam[j - 1])
    _, _, v_lo, _, _, _ = tr.analyse_point(u_lo)
    scale = max(abs(l_lo), abs(float(res.lam[j])), 1e-300)
    d = max(abs(float(res.lam[j]) - l_lo), 1e-6 * scale)
    tol = min(tr.tol, 1e-11)

    def attempt(lt):
        Du = (lt - l_lo) * v_lo
        try:
            u, _, _, _, _ = tr.correct("load", u_lo, l_lo, Du, lt - l_lo,
                                       None, tol=tol, force_only=True)
        except _Fail:
            return None
        if np.linalg.norm((u - u_lo)[tr.s.free] - Du) > np.linalg.norm(Du):
            return None                                   # jumped
        _, _, v, dsg, neg, _ = tr.analyse_point(u)
        if _point_key(neg, dsg) != key0:
            return None
        return u, v

    l_hi = None
    for _ in range(60):
        r = attempt(l_lo + d)
        if r is None:
            l_hi = l_lo + d
            break
        u_lo, v_lo = r
        l_lo = l_lo + d
        d *= 2.0
    if l_hi is None:
        raise RuntimeError("no stability change found above the sample")
    for _ in range(max_bisect):
        if l_hi - l_lo <= rtol * scale:
            break
        lt = 0.5 * (l_lo + l_hi)
        r = attempt(lt)
        if r is None:
            l_hi = lt
        else:
            u_lo, v_lo = r
            l_lo = lt
    _, _, rho = critical_mode(tr, u_lo)
    return {"lambda": float(0.5 * (l_lo + l_hi)),
            "lambda_bracket": (float(l_lo), float(l_hi)),
            "u": u_lo, "sp": float(_sp_of(tr, v_lo)), "rho": float(rho)}


# ===================================================== comparison paths ===
def linear_path(system, F_ref, lams, control_dof):
    """First-order path: u(lambda) = lambda K^-1 F_ref (control DOF only)."""
    st0 = system.state(np.zeros(system.ndof))
    v = _factor(system.tangent(st0)).solve(
        np.asarray(F_ref, dtype=float)[system.free])
    vc = float(system.expand(v)[int(control_dof)])
    lams = np.asarray(lams, dtype=float)
    return {"lam": lams.copy(), "control": lams * vc}


def pdelta_path(system, F_ref, lams, control_dof):
    """Linearised second-order path u = lambda (K + lambda Kg_ref)^-1 F_ref.

    ``Kg_ref`` is the geometric stiffness of the linear member forces under
    ``F_ref`` on the undeformed geometry.  Points where ``K + lambda Kg`` has
    lost positive definiteness (lambda >= lambda_cr) are dropped.
    """
    s = system
    F = np.asarray(F_ref, dtype=float)
    st0 = s.state(np.zeros(s.ndof))
    K0 = s.tangent(st0)
    fac = _factor(K0)
    u1 = s.expand(fac.solve(F[s.free]))
    N1 = _linear_forces(s, u1)
    Kg = slv.assemble_geometric_stiffness(s.X, s.members, N1, s.L,
                                          s.d0 / s.L[:, None])
    if _HAVE_SCIPY and sp.issparse(Kg):
        Kg = Kg.tocsr()[s.free][:, s.free].tocsc()
    else:                                                # pragma: no cover
        Kg = np.asarray(Kg)[np.ix_(s.free, s.free)]
    sgn0 = fac.det_sign()
    out_l, out_c = [], []
    for lam in np.asarray(lams, dtype=float):
        try:
            f = _factor(K0 + lam * Kg)
            if f.det_sign() != sgn0:
                break
            uf = lam * f.solve(F[s.free])
        except FloatingPointError:
            break
        out_l.append(lam)
        out_c.append(float(s.expand(uf)[int(control_dof)]))
    return {"lam": np.asarray(out_l), "control": np.asarray(out_c)}


def _linear_forces(s, u):
    """Small-displacement member forces N = EA/L c0.(u_j - u_i)."""
    du = u.reshape(-1, 3)
    delta = du[s.members[:, 1]] - du[s.members[:, 0]]
    c0 = s.d0 / s.L[:, None]
    return s.EA / s.L * np.einsum("ij,ij->i", c0, delta)


# ======================================================= model level ======
def _find_combination(model, combination_name):
    from .model import LoadCombination
    if isinstance(combination_name, LoadCombination):
        return combination_name
    combos = list(model.combinations)
    if not combos:
        raise ValueError("The model has no load combinations.")
    if combination_name is None:
        from . import loads as ld
        basic = ld.basic_load_vectors(model)
        strength = [c for c in combos if not c.serviceability] or combos
        return max(strength,
                   key=lambda c: abs(float(ld.combine(basic, c)[2::3].sum())))
    for c in combos:
        if c.name == combination_name:
            return c
    names = [c.name.strip() for c in combos]
    key = str(combination_name).strip()
    for c, nm in zip(combos, names):
        if nm == key or nm.split()[0] == key:
            return c
    raise KeyError(f"No combination named {combination_name!r} "
                   f"(have: {', '.join(c.name for c in combos)}).")


def resolve_control(model, control):
    """Control DOF from 'crown_uz' / 'crown_ux' / int DOF / (node, axis)."""
    if control is None:
        control = "crown_uz"
    if isinstance(control, (int, np.integer)):
        return int(control)
    if isinstance(control, str):
        key = control.strip().lower()
        if key.startswith("crown_"):
            ax = _AXIS.index(key.split("_", 1)[1])
            return 3 * int(model.geom.crown_node) + ax
        raise ValueError(f"Unknown control {control!r}.")
    node, ax = control
    if isinstance(ax, str):
        ax = _AXIS.index(ax.lower())
    return 3 * int(node) + int(ax)


def model_system(model, node_offset=None):
    """TrussSystem of a DomeModel (supports and springs from the model)."""
    fixed, springs = model.support_conditions()
    X = np.array(model.geom.nodes, dtype=float)
    if node_offset is not None:
        X = X + np.asarray(node_offset, dtype=float).reshape(-1, 3)
    return TrussSystem(X, model.geom.members, model.EA(), fixed, springs)


def reference_load(model, combination):
    """D/S/L/W load vector of a combination (T excluded)."""
    from . import loads as ld
    return ld.combine(ld.basic_load_vectors(model), combination)


def linear_buckling_reference(model, F_ref):
    """(lambda_cr, mode) of ``F_ref`` from buckling.linear_buckling.

    The sign of an eigenvector is arbitrary and differs between LAPACK
    builds.  Each mode is turned so that it does positive work with the
    reference load (phi . F_ref > 0): a positive imperfection amplitude is
    then always in the direction the load pushes (for a dome under gravity,
    the crown dimple points down), on every platform.  A mode orthogonal to
    the load keeps the sign of buckling.linear_buckling (largest component
    positive).
    """
    from . import buckling as bk
    r = slv.solve_model(model, F_ref, name="ref", T_factor=0.0)
    b = bk.linear_buckling(model, r.N, n_modes=6)
    F = np.asarray(F_ref, dtype=float)
    nF = float(np.linalg.norm(F))
    for k in range(b.modes.shape[1]):
        v = b.modes[:, k]
        w = float(v @ F)
        if abs(w) > 1e-9 * float(np.linalg.norm(v)) * nF and w < 0.0:
            b.modes[:, k] = -v
    if len(b.load_factors):
        return float(b.load_factors[0]), b.modes[:, 0].copy(), b
    return float("nan"), None, b


def run_nonlinear(model, combination_name=None, method="arc-length",
                  control="crown_uz", *, first_increment=None,
                  lambda_max=None, control_max=None, stop_after_limit=5,
                  max_steps=300, tol=1e-8, tol_disp=1e-6, max_iter=30,
                  desired_iterations=5, refine_limit=True,
                  refine_bifurcation=True, compare=True, buckling=True,
                  lambda_cr=None, n_compare=60, eig_limit=600,
                  node_offset=None, load_vector=None, member_checks=True,
                  progress=None, **trace_options):
    """Geometrically nonlinear analysis of a :class:`DomeModel`.

    ``lambda`` multiplies the D/S/L/W part of the combination (self-straining
    T is excluded, see the module docstring).  Supports and springs come
    from ``model.support_conditions()``.

    Parameters (beyond :func:`trace_path`)
    --------------------------------------
    combination_name : name (or leading token such as 'LC2') of a model
        combination, a LoadCombination, or None = the strength combination
        with the largest vertical load.
    control : 'crown_uz' (default), 'crown_ux', an int DOF or (node, axis).
    first_increment : default lambda_cr / 40 (0.05 when lambda_cr is unknown).
    lambda_max : default 2 lambda_cr (20 when unknown).
    control_max : default half the structure's length scale.
    buckling : run buckling.linear_buckling for lambda_cr (or pass
        ``lambda_cr`` to reuse a known value).
    compare : also compute the linear and P-Delta paths.
    node_offset : (n, 3) coordinate perturbation (imperfections).
    load_vector : override F_ref (used by the imperfection study so that the
        imperfect structure carries exactly the perfect structure's loads).
    member_checks : report the first strut reaching its Euler load and its
        AISC design strength along the path.
    trace_options : any other :func:`trace_path` keyword (arc_length,
        max_arc_factor, max_turn, lambda_min, ...).
    """
    t0 = time.perf_counter()
    msgs = []
    combo = _find_combination(model, combination_name)
    F = reference_load(model, combo) if load_vector is None else \
        np.asarray(load_vector, dtype=float)
    if abs(getattr(combo, "T", 0.0)) > 1e-12 and model.loads.has_self_straining:
        msgs.append(f"The self-straining load T (factor {combo.T:g}) of "
                    f"{combo.name} is not included in the nonlinear path.")

    work = model
    arc = _METHODS.get(str(method).lower(), "") == "arc-length"
    if node_offset is not None and arc and \
            "arc_length" not in trace_options and \
            "first_arc_length" not in trace_options:
        # An imperfection can put a node almost flat with its neighbours;
        # its local snap then lies within a millimetre or so of the start
        # of the path.  Start with a step of 2e-5 times the length scale
        # (0.3 mm for a 15 m radius); the step grows to the nominal size.
        trace_options["first_arc_length"] = 2e-5 * float(
            model.geom.length_scale)
    if node_offset is not None:
        work = copy.copy(model)
        g = dataclasses.replace(
            model.geom, nodes=np.asarray(model.geom.nodes, float) +
            np.asarray(node_offset, float).reshape(-1, 3))
        work.geom = g
    system = model_system(work)
    cdof = resolve_control(model, control)

    mode = None
    lam_cr = float("nan") if lambda_cr is None else float(lambda_cr)
    if lambda_cr is None and buckling:
        if progress:
            progress(0.02, "Linear buckling reference...")
        try:
            lam_cr, mode, _ = linear_buckling_reference(work, F)
        except Exception as exc:
            msgs.append(f"Linear buckling reference failed: {exc}")
    finite_cr = np.isfinite(lam_cr) and lam_cr > 0.0
    if first_increment is None:
        first_increment = lam_cr / 40.0 if finite_cr else 0.05
    if lambda_max is None:
        lambda_max = 2.0 * lam_cr if finite_cr else 20.0
    if control_max is None:
        control_max = 0.5 * float(work.geom.length_scale)

    res = trace_path(system, F, method=method, control_dof=cdof,
                     first_increment=first_increment, lambda_max=lambda_max,
                     control_max=control_max,
                     stop_after_limit=stop_after_limit, max_steps=max_steps,
                     tol=tol, tol_disp=tol_disp, max_iter=max_iter,
                     desired_iterations=desired_iterations,
                     refine_limit=refine_limit,
                     refine_bifurcation=refine_bifurcation,
                     eig_limit=eig_limit, progress=progress,
                     **trace_options)
    res.combination = combo.name
    res.lambda_cr = lam_cr
    res.buckling_mode = mode
    res.control_label = f"node {cdof // 3} {_AXIS[cdof % 3]}" + \
        (" (crown)" if cdof // 3 == int(model.geom.crown_node) else "")
    res.messages = msgs + res.messages

    if compare:
        top = float(np.max(np.abs(res.lam))) if len(res.lam) else 1.0
        lam_end = 1.25 * top
        if finite_cr:
            lam_end = min(lam_end, 0.98 * lam_cr)
        lams = np.linspace(0.0, lam_end, int(max(n_compare, 2)))
        res.linear = linear_path(system, F, lams, cdof)
        res.p_delta = pdelta_path(system, F, lams, cdof)

    if member_checks:
        try:
            _member_checks(work, res)
        except Exception as exc:                          # pragma: no cover
            msgs.append(f"Member checks skipped ({exc}).")

    if res.limit_found and res.limit_step is not None and \
            np.any(res.lam[:int(res.limit_step) + 1] <
                   -1e-9 * max(abs(res.lambda_limit), 1.0)):
        # The path ran backwards (to negative loads) before its first
        # maximum: that maximum is not the limit load of the primary path.
        res.limit_found = False
        res.stop_reason = ("the path went to a negative load before its "
                           "first maximum")
        res.messages.append(
            "The path went to a negative load before its first maximum, so "
            "that maximum is not used as a limit point.  Try a smaller first "
            "increment.")
    if not res.limit_found:
        res.messages.append(f"No limit point reached: {res.stop_reason}.")
    else:
        kind = ("bifurcation" if res.limit_type == "bifurcation"
                else "limit point")
        txt = (f"First load maximum ({kind}) lambda = "
               f"{res.lambda_limit:.4g}")
        if finite_cr:
            txt += (f" = {res.limit_ratio:.3f} x lambda_cr "
                    f"({lam_cr:.4g}, linear buckling)")
        res.messages.append(txt + ".")
    bif = [cp for cp in res.critical_points if cp["type"] == "bifurcation"
           and not cp.get("first_peak")]
    if bif:
        res.messages.append(
            f"Stiffness sign change on the rising path (bifurcation) at "
            f"lambda = {bif[0]['lambda']:.4g}: the perfect structure can "
            f"buckle into another mode there; run imperfection_study.")
    if np.isfinite(res.lambda_member_euler) and (
            not res.limit_found or res.lambda_member_euler < res.lambda_limit):
        res.messages.append(
            f"Strut {res.member_euler} reaches its Euler load at lambda = "
            f"{res.lambda_member_euler:.4g}, before the global limit: member "
            f"buckling governs (the truss model does not include it).")
    res.elapsed = time.perf_counter() - t0
    return res


def _first_crossing(lams, ratios):
    """First lambda (linear interpolation) where max ratio reaches 1."""
    r = np.max(ratios, axis=1)
    for k in range(1, len(r)):
        if r[k] >= 1.0 > r[k - 1]:
            f = (1.0 - r[k - 1]) / (r[k] - r[k - 1])
            return float(lams[k - 1] + f * (lams[k] - lams[k - 1])), \
                int(np.argmax(ratios[k]))
    if len(r) and r[0] >= 1.0:
        return float(lams[0]), int(np.argmax(ratios[0]))
    return float("nan"), -1


def _member_checks(model, res):
    from . import buckling as bk
    caps = bk.member_capacities(model)
    k_end = res.limit_step if res.limit_found else res.n_points - 1
    lam = res.lam[:k_end + 1]
    N = res.N[:k_end + 1]
    comp = np.maximum(-N, 0.0)
    res.lambda_member_euler, res.member_euler = _first_crossing(
        lam, comp / np.maximum(caps["Pe_euler"], 1e-300))
    util = np.where(N >= 0.0, N / np.maximum(caps["phiPn_t"], 1e-300),
                    comp / np.maximum(caps["phiPn_c"], 1e-300))
    res.lambda_member_design, res.member_design = _first_crossing(lam, util)


# ===================================================== imperfections ======
def imperfection_study(model, combination_name=None, amplitudes=None, mode=0,
                       signs=(1.0, -1.0), return_results=False, progress=None,
                       **kwargs):
    """Limit load of mode-affine imperfect domes.

    Node coordinates are perturbed by ``a * phi`` where ``phi`` is linear
    buckling mode ``mode`` of the combination, scaled to max |component| = 1
    (so ``a`` is the largest coordinate offset).  Both signs are run and the
    lower limit load is kept (Bruno et al. 2016; EN 1993-1-6 uses the
    eigenmode-affine pattern).  The imperfect structure carries the perfect
    structure's load vector.

    amplitudes : list of amplitudes [m]; default span/1000, span/500, span/300
        (the range used for single-layer grid shells, Bruno et al. 2016).
        When the buckling mode is local (one node, as at the 5-valent crown
        of a geodesic dome) the amplitude is better judged against the strut
        length, so ``amplitude_over_strut`` is reported too.
    kwargs : passed to :func:`run_nonlinear`.

    Returns a list of row dicts: amplitude_m, amplitude_over_span,
    span_over_amplitude, amplitude_over_strut, lambda_limit (first load
    maximum of the governing sign), lambda_cr, ratio (lambda_limit /
    lambda_cr), lambda_critical and ratio_critical (first critical point of
    that run: the load maximum or an earlier bifurcation into another mode),
    sign (+1/-1, the sign with the lower lambda_limit), lambda_plus,
    lambda_minus, found.  The first row is the perfect structure (a = 0).
    With ``return_results=True`` returns ``(rows, results)`` where results
    maps (amplitude, sign) -> NonlinearResult.
    """
    combo = _find_combination(model, combination_name)
    F = reference_load(model, combo)
    lam_cr, _, b = linear_buckling_reference(model, F)
    if not len(b.load_factors):
        raise ValueError("No positive linear buckling factor: the "
                         "combination does not compress the structure.")
    mode = int(mode)
    if mode >= b.modes.shape[1]:
        raise ValueError(f"Only {b.modes.shape[1]} buckling modes available.")
    phi = b.modes[:, mode].copy()
    phi /= max(float(np.max(np.abs(phi))), 1e-300)
    span = float(model.geom.base_diameter)
    strut = float(np.mean(model.geom.member_lengths()))
    if amplitudes is None:
        amplitudes = [span / d for d in (1000.0, 500.0, 300.0)]
    amps = [0.0] + [float(a) for a in amplitudes if float(a) != 0.0]

    rows, results = [], {}
    kw = dict(kwargs)
    kw.setdefault("compare", False)
    kw.setdefault("refine_bifurcation", True)
    kw.setdefault("member_checks", False)
    n_runs = 1 + 2 * (len(amps) - 1)
    run = 0
    for a in amps:
        lam_s, crit_s = {}, {}
        for sgn in ((1.0,) if a == 0.0 else signs):
            if progress:
                progress(run / n_runs, f"Imperfection {a * 1e3:.1f} mm, "
                                       f"sign {sgn:+.0f}")
            run += 1
            r = run_nonlinear(model, combo, node_offset=(sgn * a * phi)
                              .reshape(-1, 3) if a else None,
                              load_vector=F, lambda_cr=lam_cr, **kw)
            results[(a, sgn)] = r
            lam_s[sgn] = r.lambda_limit if r.limit_found else float("nan")
            crit_s[sgn] = r.lambda_first_critical
        finite = {k: v for k, v in lam_s.items() if np.isfinite(v)}
        if finite:
            sg = min(finite, key=finite.get)
            lam_lim = finite[sg]
        else:
            sg, lam_lim = 1.0, float("nan")
        lam_c = crit_s.get(sg, float("nan"))
        rows.append({
            "amplitude_m": a,
            "amplitude_over_span": a / span if span > 0 else float("nan"),
            "span_over_amplitude": span / a if a > 0 else float("inf"),
            "amplitude_over_strut": a / strut if strut > 0 else float("nan"),
            "lambda_limit": lam_lim,
            "lambda_cr": lam_cr,
            "ratio": lam_lim / lam_cr if np.isfinite(lam_lim) else float("nan"),
            "lambda_critical": lam_c,
            "ratio_critical": lam_c / lam_cr if np.isfinite(lam_c)
            else float("nan"),
            "sign": sg,
            "lambda_plus": lam_s.get(1.0, float("nan")),
            "lambda_minus": lam_s.get(-1.0, float("nan")),
            "found": bool(np.isfinite(lam_lim)),
        })
    if progress:
        progress(1.0, "Imperfection study complete.")
    return (rows, results) if return_results else rows


# ===================================================== closed forms =======
def von_mises_load(b, h, EA, w):
    """Exact engineering-strain path of the symmetric two-bar truss.

    Half span ``b``, rise ``h``, apex deflection ``w`` (down +):
    ``P(w) = 2 EA y (1/l - 1/L)`` with ``y = h - w``, ``l = sqrt(b^2 + y^2)``,
    ``L = sqrt(b^2 + h^2)`` (closed form derived here for engineering
    strain; the problem studied by Pecknold et al. 1985).
    """
    w = np.asarray(w, dtype=float)
    L = math.hypot(b, h)
    y = h - w
    l = np.sqrt(b * b + y * y)
    return 2.0 * EA * y * (1.0 / l - 1.0 / L)


def von_mises_limit(b, h, EA):
    """(w*, P*) of the first limit point: l* = (b^2 L)^(1/3)."""
    L = math.hypot(b, h)
    ls = (b * b * L) ** (1.0 / 3.0)
    ys = math.sqrt(max(ls * ls - b * b, 0.0))
    return h - ys, 2.0 * EA * ys * (1.0 / ls - 1.0 / L)
