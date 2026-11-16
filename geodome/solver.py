"""
Three-dimensional space-truss direct stiffness solver.

Element formulation
-------------------
Each strut is a two-node pin-ended bar with three translational degrees of
freedom per node.  With direction cosines

    {c} = (x2 - x1)/L ,  (y2 - y1)/L ,  (z2 - z1)/L

the 6x6 element stiffness matrix in global coordinates is

          EA  |  [c c^T]   -[c c^T] |
    k  =  --  |                     |
          L   | -[c c^T]    [c c^T] |

which is the familiar ``k = T^T k_local T`` with the 1x6 transformation
``T = [-c^T  c^T]`` and ``k_local = EA/L``.

Solution
--------
The global matrix is assembled in sparse COO form, partitioned into free (f)
and supported (s) degrees of freedom, and the free partition solved:

    [Kff Kfs] {uf}   {Ff}
    [Ksf Kss] {us} = {Fs}      with  us = 0

    Kff uf = Ff            ->  uf
    R = Ksf uf - Fs            ->  support reactions

Member axial force follows from the element compatibility relation

    N = (EA/L) * {c}^T ({u2} - {u1})

with tension positive.  Three independent verification quantities are returned
with every solution: the global equilibrium residual, the reaction/applied
load balance, and the strain-energy identity

    1/2 {F}^T {u}  ==  1/2 sum_m  N_m^2 L_m / (E A)_m
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False


# ---------------------------------------------------------------- results ---
@dataclass
class StaticResult:
    """Everything produced by one static solve."""

    name: str
    u: np.ndarray                 # (3n,) nodal displacements         [m]
    N: np.ndarray                 # (m,)  member axial forces, + = tension [N]
    stress: np.ndarray            # (m,)  axial stress                [Pa]
    reactions: np.ndarray         # (3n,) reaction vector (0 at free dofs) [N]
    F: np.ndarray                 # (3n,) applied load vector         [N]
    lengths: np.ndarray           # (m,)  member lengths              [m]

    # verification / diagnostics
    equilibrium_residual: float = 0.0
    reaction_balance: np.ndarray = field(default_factory=lambda: np.zeros(3))
    strain_energy_ext: float = 0.0
    strain_energy_int: float = 0.0
    solve_info: str = ""

    # ------------------------------------------------------------ queries --
    def displacement_vectors(self):
        return self.u.reshape(-1, 3)

    def displacement_magnitude(self):
        return np.linalg.norm(self.displacement_vectors(), axis=1)

    def max_displacement(self):
        d = self.displacement_magnitude()
        i = int(np.argmax(d))
        return i, float(d[i])

    def reaction_vectors(self):
        return self.reactions.reshape(-1, 3)

    def total_reaction(self):
        R = self.reaction_vectors()
        return R.sum(axis=0)

    def max_tension(self):
        i = int(np.argmax(self.N))
        return i, float(self.N[i])

    def max_compression(self):
        i = int(np.argmin(self.N))
        return i, float(self.N[i])

    @property
    def energy_error(self):
        denom = max(abs(self.strain_energy_ext), 1e-30)
        return abs(self.strain_energy_ext - self.strain_energy_int) / denom


# --------------------------------------------------------------- geometry ---
def direction_cosines(nodes, members):
    """Return ``(L, C)`` with lengths (m,) and unit direction cosines (m, 3)."""
    p = nodes[members[:, 0]]
    q = nodes[members[:, 1]]
    v = q - p
    L = np.linalg.norm(v, axis=1)
    if np.any(L <= 0.0):
        bad = int(np.argmin(L))
        raise ValueError(f"Member {bad} has zero length.")
    return L, v / L[:, None]


def element_stiffness(EA, L, c):
    """6x6 global element stiffness matrix for one bar."""
    cc = np.outer(c, c)
    k = np.empty((6, 6))
    k[0:3, 0:3] = cc
    k[0:3, 3:6] = -cc
    k[3:6, 0:3] = -cc
    k[3:6, 3:6] = cc
    return (EA / L) * k


def element_dofs(members):
    """(m, 6) array of global DOF numbers for every member."""
    i = members[:, 0]
    j = members[:, 1]
    return np.column_stack([3 * i, 3 * i + 1, 3 * i + 2,
                            3 * j, 3 * j + 1, 3 * j + 2])


# --------------------------------------------------------------- assembly ---
def assemble_stiffness(nodes, members, EA, L=None, C=None):
    """Assemble the global stiffness matrix (sparse CSR if SciPy is present)."""
    if L is None or C is None:
        L, C = direction_cosines(nodes, members)
    ndof = 3 * len(nodes)
    dofs = element_dofs(members)

    nm = len(members)
    rows = np.empty(36 * nm, dtype=np.int64)
    cols = np.empty(36 * nm, dtype=np.int64)
    vals = np.empty(36 * nm, dtype=float)

    for m in range(nm):
        ke = element_stiffness(EA[m], L[m], C[m])
        d = dofs[m]
        s = 36 * m
        rows[s:s + 36] = np.repeat(d, 6)
        cols[s:s + 36] = np.tile(d, 6)
        vals[s:s + 36] = ke.ravel()

    if _HAVE_SCIPY:
        K = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    else:                                            # pragma: no cover
        K = np.zeros((ndof, ndof))
        np.add.at(K, (rows, cols), vals)
    return K


def assemble_geometric_stiffness(nodes, members, N, L=None, C=None):
    """Global geometric (initial stress) stiffness matrix.

    For a bar carrying axial force ``N`` the consistent geometric stiffness is

              N  |  [I - c c^T]   -[I - c c^T] |
        kg =  -  |                             |
              L  | -[I - c c^T]    [I - c c^T] |

    which captures the transverse restoring/destabilising effect of the axial
    force and is what makes the linearised buckling eigenproblem possible.
    """
    if L is None or C is None:
        L, C = direction_cosines(nodes, members)
    ndof = 3 * len(nodes)
    dofs = element_dofs(members)
    I3 = np.eye(3)

    nm = len(members)
    rows = np.empty(36 * nm, dtype=np.int64)
    cols = np.empty(36 * nm, dtype=np.int64)
    vals = np.empty(36 * nm, dtype=float)

    for m in range(nm):
        c = C[m]
        M = I3 - np.outer(c, c)
        kg = np.empty((6, 6))
        kg[0:3, 0:3] = M
        kg[0:3, 3:6] = -M
        kg[3:6, 0:3] = -M
        kg[3:6, 3:6] = M
        kg *= (N[m] / L[m])
        d = dofs[m]
        s = 36 * m
        rows[s:s + 36] = np.repeat(d, 6)
        cols[s:s + 36] = np.tile(d, 6)
        vals[s:s + 36] = kg.ravel()

    if _HAVE_SCIPY:
        return sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
    Kg = np.zeros((ndof, ndof))                      # pragma: no cover
    np.add.at(Kg, (rows, cols), vals)
    return Kg


# ------------------------------------------------------------------ solve ---
def free_dofs(ndof, fixed):
    mask = np.ones(ndof, dtype=bool)
    mask[np.asarray(fixed, dtype=int)] = False
    return np.where(mask)[0]


def solve_static(nodes, members, EA, F, fixed, name="LC"):
    """Solve one static load case and return a :class:`StaticResult`."""
    nodes = np.asarray(nodes, dtype=float)
    members = np.asarray(members, dtype=int)
    EA = np.asarray(EA, dtype=float)
    F = np.asarray(F, dtype=float)

    ndof = 3 * len(nodes)
    L, C = direction_cosines(nodes, members)
    K = assemble_stiffness(nodes, members, EA, L, C)

    fd = free_dofs(ndof, fixed)
    sd = np.asarray(sorted(set(range(ndof)) - set(fd.tolist())), dtype=int)

    info = []
    u = np.zeros(ndof)
    if len(fd) == 0:
        info.append("Every degree of freedom is restrained.")
    else:
        if _HAVE_SCIPY:
            Kff = K[fd][:, fd].tocsc()
            # Guard against an unstable mechanism: a zero diagonal means a node
            # has no stiffness in that direction at all.
            diag = Kff.diagonal()
            if np.any(np.abs(diag) < 1e-12 * max(np.abs(diag).max(), 1.0)):
                bad = fd[np.where(np.abs(diag) < 1e-12 * np.abs(diag).max())[0]]
                raise np.linalg.LinAlgError(
                    f"Unstable structure: zero stiffness at DOF(s) "
                    f"{bad[:8].tolist()} (node {bad[0] // 3}).")
            try:
                uf = spla.spsolve(Kff, F[fd])
                if not np.all(np.isfinite(uf)):
                    raise RuntimeError("non-finite solution")
                info.append("Sparse LU (SuperLU).")
            except Exception as exc:                 # pragma: no cover
                info.append(f"Sparse solve fell back to dense ({exc}).")
                uf = np.linalg.solve(Kff.toarray(), F[fd])
        else:                                        # pragma: no cover
            uf = np.linalg.solve(K[np.ix_(fd, fd)], F[fd])
            info.append("Dense LAPACK solve.")
        u[fd] = uf

    # ---- member forces ----------------------------------------------------
    d = element_dofs(members)
    du = u[d[:, 3:6]] - u[d[:, 0:3]]
    elong = np.einsum("ij,ij->i", C, du)
    N = EA / L * elong

    # ---- reactions --------------------------------------------------------
    Ku = K @ u if _HAVE_SCIPY else K.dot(u)
    reactions = np.zeros(ndof)
    reactions[sd] = Ku[sd] - F[sd]

    # ---- verification quantities -----------------------------------------
    residual = Ku - F - reactions
    rnorm = float(np.linalg.norm(residual) / max(np.linalg.norm(F), 1e-12))

    applied = F.reshape(-1, 3).sum(axis=0)
    reac = reactions.reshape(-1, 3).sum(axis=0)
    balance = applied + reac

    Ue = 0.5 * float(F @ u)
    Ui = 0.5 * float(np.sum(N ** 2 * L / EA))

    return StaticResult(
        name=name, u=u, N=N, stress=np.zeros_like(N), reactions=reactions,
        F=F, lengths=L, equilibrium_residual=rnorm, reaction_balance=balance,
        strain_energy_ext=Ue, strain_energy_int=Ui, solve_info=" ".join(info))


def solve_model(model, F, name="LC"):
    """Solve ``model`` for a given global load vector."""
    res = solve_static(model.geom.nodes, model.geom.members, model.EA(),
                       F, model.fixed_dofs(), name=name)
    res.stress = res.N / model.member_areas()
    return res
