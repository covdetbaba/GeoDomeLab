"""
OpenSees bridge: export a GeoDomeLab model to OpenSees and cross-check it.

OpenSees (McKenna, Scott & Fenves 2010; Python interface OpenSeesPy, Zhu,
McKenna & Scott 2018) is an independent, widely used open-source finite
element framework.  Solving the same structure in OpenSees and in
GeoDomeLab, with nothing shared except the input data, is a direct check of
GeoDomeLab's stiffness assembly, support handling, self-straining loads,
eigen solution and geometrically nonlinear path tracing.

What is shared and what is independent
--------------------------------------
:func:`to_opensees_spec` turns a :class:`DomeModel` and one load combination
into a plain, JSON-able dict: node coordinates, connectivity, E and A of
every strut, restrained DOFs, spring stiffnesses, the factored nodal load
vector (D/S/L/W), the free strains ``eps0 = T alpha dT`` and the prescribed
settlements ``T delta``, and the masses (``rho A`` of every strut and the
non-structural nodal masses of ``modal.lumped_masses``).  These are the
*inputs*, so the cross-check cannot see an error in them.  Checked by hand
elsewhere: the dead and snow tributary loads and the crown load on Models A
and B (``validation/run_validation.py``), the wind field on Model B (case F9
of ``validation/run_validation_ext.py``), and E, alpha and rho through the
typed-in constants of the hand files.  Not checked by a hand calculation:
the unbalanced snow pattern and the default combination factors (unit tests
only), and the 7-point Cp table itself (an illustrative simplification, not
ASCE 7-16 Fig. 27.3-2).  Everything after the inputs is computed
independently by OpenSees:

* stiffness: ``element('Truss')`` (linear) or ``element('corotTruss')``
  (geometrically nonlinear), ``uniaxialMaterial('Elastic')``;
* thermal strain: ``uniaxialMaterial('InitStrainMaterial', ..., -eps0)``,
  i.e. N = EA (eps - eps0) -- the sign is confirmed on a heated fixed bar
  (:func:`thermal_sign_check`);
* settlement: ``sp`` inside the load pattern on an unfixed DOF;
* elastic supports: ``zeroLength`` springs to a fixed ground node (its
  reaction = the spring force on the structure);
* strut mass: ``'-rho', rho*A`` (lumped rho A L/2 per end, or consistent
  with ``'-cMass', 1``), non-structural mass with ``mass``;
* solutions: ``analyze`` (Linear algorithm, LAPACK or UMFPACK), ``eigen``
  (ARPACK or LAPACK generalised), ``DisplacementControl`` and ``ArcLength``
  Newton path following.

Strain measure of ``corotTruss``
--------------------------------
GeoDomeLab's corotational bar uses engineering strain on the original
length, eps = (l - L)/L.  :func:`strain_measure_check` determines what
OpenSees' ``corotTruss`` uses *numerically*: one bar is given a large rigid
rotation (transverse end displacement 0.8 L) and a 10 % stretch; the
OpenSees axial force is compared with

    engineering      N = EA (l - L)/L
    Green-Lagrange   N = EA (l^2 - L^2)/(2 L^2)       (2nd Piola-Kirchhoff)
    logarithmic      N = EA ln(l/L)

At a 10 % stretch these differ by about 5 %, so the test is decisive (the
installed OpenSeesPy 3.8.0 matches engineering strain to round-off, as the
source file CorotTruss.cpp says).  Paths of the two programs can therefore
be compared at any strain, not only in the small-strain regime.

Comparison measures
-------------------
Pointwise relative differences are meaningless for quantities that pass
through zero (the y-reaction of a symmetric dome, a lambda path crossing
lambda = 0), so differences are normalised by the largest magnitude of the
quantity in the case (normwise, infinity norm):

    d_u = max_i |u_GDL,i - u_OS,i| / max_i |u_GDL,i|

and likewise for member forces N (scale: the larger of max |N| and the fully
restrained thermal force max EA |eps0|), reactions R (scale: the larger of
max |R| and max |F|), frequencies (pointwise, f > 0), the load factor along
a path (scale: max |lambda| on the path).  Mode shapes are compared with the
modal assurance criterion

    MAC(a, b) = (a.b)^2 / ((a.a)(b.b))

and, inside a cluster of repeated frequencies (a symmetric dome has many
exactly repeated pairs whose basis is arbitrary), with the subspace MAC
||Q^T a||^2 / ||a||^2 (Q = orthonormal basis of the other program's modes of
that cluster).

Running OpenSees
----------------
:func:`run_opensees` uses ``openseespy`` in-process when it imports in the
current interpreter, otherwise it runs ``geodome/_opensees_runner.py`` in an
OpenSees interpreter as a subprocess (JSON in and out).  That interpreter is
the environment variable ``GEODOME_OPENSEES_PYTHON`` if set, else the first
existing of ``<conda base>/envs/opensees/python(.exe)`` (and the same under
``~/anaconda3`` and ``~/miniconda3``).  :func:`write_openseespy_script`
writes a standalone, readable OpenSeesPy script of any model that anyone
can run with ``pip install openseespy``.

References
----------
F. McKenna, M. H. Scott, G. L. Fenves (2010). Nonlinear finite-element
analysis software architecture using object composition. J. Computing in
Civil Engineering 24(1), 95-107. DOI 10.1061/(ASCE)CP.1943-5487.0000002.
M. Zhu, F. McKenna, M. H. Scott (2018). OpenSeesPy: Python library for the
OpenSees finite element framework. SoftwareX 7, 6-11.
DOI 10.1016/j.softx.2017.10.009.
OpenSees source: SRC/element/truss/CorotTruss.cpp (M. H. Scott),
SRC/material/uniaxial/InitStrainMaterial.cpp (F. McKenna),
SRC/analysis/integrator/ArcLength.cpp, DisplacementControl.cpp.
M. A. Crisfield (1991). Non-linear Finite Element Analysis of Solids and
Structures, Vol. 1. Wiley, Ch. 3 (truss strain measures).
R. J. Allemang (2003). The modal assurance criterion -- twenty years of use
and abuse. Sound and Vibration 37(8), 14-23.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

from . import _opensees_runner as runner

ENV_VAR = "GEODOME_OPENSEES_PYTHON"
RUNNER_PATH = os.path.abspath(runner.__file__)
_AXIS = ("ux", "uy", "uz")


class OpenSeesUnavailable(RuntimeError):
    """Neither in-process openseespy nor an OpenSees interpreter was found."""


# ================================================== finding OpenSees ======
_INPROCESS = None          # cached: True / False


def inprocess_available():
    """True when ``openseespy`` imports in this interpreter (cached)."""
    global _INPROCESS
    if _INPROCESS is None:
        try:
            import openseespy.opensees  # noqa: F401
            _INPROCESS = True
        except Exception:                             # ImportError, RuntimeError
            _INPROCESS = False
    return _INPROCESS


def _python_name():
    return "python.exe" if os.name == "nt" else os.path.join("bin", "python")


def opensees_python_candidates():
    """Interpreters that may hold openseespy, in the order they are tried."""
    out = []
    env = os.environ.get(ENV_VAR, "").strip()
    if env:
        out.append(env)
    exe = _python_name()
    bases = []
    for p in (os.environ.get("CONDA_EXE"),):
        if p:                       # <base>/Scripts/conda.exe or <base>/bin/conda
            bases.append(os.path.dirname(os.path.dirname(p)))
    for p in (sys.base_prefix, sys.prefix):
        bases.append(p)
        # inside an env: <base>/envs/<name>
        parent = os.path.dirname(p)
        if os.path.basename(parent).lower() == "envs":
            bases.append(os.path.dirname(parent))
    home = os.path.expanduser("~")
    bases += [os.path.join(home, "anaconda3"), os.path.join(home, "miniconda3"),
              os.path.join(home, "miniforge3")]
    for b in bases:
        cand = os.path.join(b, "envs", "opensees", exe)
        if cand not in out:
            out.append(cand)
    return out


def find_opensees_python():
    """Path of an existing OpenSees interpreter, or None."""
    for c in opensees_python_candidates():
        if c and os.path.isfile(c):
            return c
    return None


def opensees_status():
    """``{"mode": "in-process" | "subprocess" | None, "python": path}``."""
    if inprocess_available():
        return {"mode": "in-process", "python": sys.executable}
    py = find_opensees_python()
    if py:
        return {"mode": "subprocess", "python": py}
    return {"mode": None, "python": None}


def opensees_available():
    return opensees_status()["mode"] is not None


# ============================================================ the spec ====
def find_combination(model, combination_name):
    """A model combination by exact name, stripped name, leading token or
    unique leading part of the name.

    ``None`` returns ``None`` (no load: a modal-only spec).  A
    :class:`LoadCombination` is returned unchanged.
    """
    from .model import LoadCombination
    if combination_name is None:
        return None
    if isinstance(combination_name, LoadCombination):
        return combination_name
    key = str(combination_name).strip()
    for c in model.combinations:
        if c.name == combination_name or c.name.strip() == key:
            return c
    for c in model.combinations:
        if c.name.split() and c.name.split()[0] == key:
            return c
    pref = [c for c in model.combinations if c.name.strip().startswith(key)]
    if len(pref) == 1:                    # a unique leading part of the name
        return pref[0]
    raise KeyError(f"No combination named {combination_name!r} (have: "
                   f"{', '.join(c.name for c in model.combinations)}).")


def to_opensees_spec(model, combination_name=None, *, include_mass=True,
                     include_dead=True, nodal_masses=None, snow_fraction=0.0):
    """Plain JSON-able description of ``model`` under one combination.

    Keys: ``nodes`` (n x 3, m), ``members`` (m x 2, 0-based), ``E`` [Pa] and
    ``A`` [m^2] per member, ``fixed`` (restrained global DOFs, DOF =
    3 node + axis), ``springs`` ([dof, k N/m]), ``load`` (3n factored D/S/L/W
    vector [N]), ``eps0`` (per-member free strain T alpha dT, or None),
    ``prescribed`` ([dof, m] = T x settlement), ``line_mass`` (rho A [kg/m])
    and ``nodal_mass`` (non-structural lumped mass [kg] per node, from
    ``modal.lumped_masses``), plus descriptive metadata.  With
    ``combination_name=None`` there is no load and no T.
    """
    from . import loads as ld
    from .solver import model_self_straining
    g = model.geom
    n = g.n_nodes
    combo = find_combination(model, combination_name)
    if combo is None:
        F = np.zeros(3 * n)
        eps0, presc = None, None
    else:
        F = ld.combine(ld.basic_load_vectors(model), combo)
        eps0, presc = model_self_straining(model, float(combo.T))
    fixed, springs = model.support_conditions()
    E = np.full(g.n_members, float(model.material.E))
    A = np.asarray(model.member_areas(), dtype=float)
    spec = {
        "format": "GeoDomeLab OpenSees spec", "version": 1,
        "title": (g.title if getattr(g, "is_custom", False) else
                  f"{g.frequency}V geodesic dome, R = {g.radius:g} m, "
                  f"fraction {g.dome_fraction:g}"),
        "combination": combo.name if combo is not None else None,
        "factors": combo.factors() if combo is not None else None,
        "units": "N, m, kg, s",
        "support_type": model.support_type,
        "crown_node": int(g.crown_node),
        "nodes": [[float(v) for v in row] for row in g.nodes],
        "members": [[int(v) for v in row] for row in g.members],
        "E": E.tolist(),
        "A": A.tolist(),
        "fixed": [int(d) for d in np.asarray(fixed, dtype=int)],
        "springs": [[int(d), float(k)] for d, k in sorted(springs.items())],
        "load": [float(v) for v in F],
        "eps0": None if eps0 is None else [float(v) for v in eps0],
        "prescribed": [[int(d), float(v)] for d, v in
                       sorted((presc or {}).items())],
    }
    if include_mass:
        from . import modal as md
        mm = md.lumped_masses(model, include_dead=include_dead,
                              nodal_masses=nodal_masses,
                              snow_fraction=snow_fraction)
        spec["line_mass"] = [float(v) for v in mm.line_mass]
        spec["nodal_mass"] = [float(v) for v in mm.non_structural]
    return spec


# ===================================================== running OpenSees ===
_ARRAY_KEYS = ("u", "N", "reactions", "spring_forces", "lam", "control",
               "omega2", "frequencies", "iterations")


def _to_numpy(res):
    out = dict(res)
    for k in _ARRAY_KEYS:
        if k in out and isinstance(out[k], list):
            out[k] = np.asarray(out[k], dtype=float)
    if out.get("type") == "modal" and isinstance(out.get("modes"), list):
        M = np.asarray(out["modes"], dtype=float)
        out["modes"] = M.T if M.ndim == 2 else M.reshape(-1, 0)
    return out


def run_opensees_batch(jobs, *, python=None, in_process=None, timeout=1800):
    """Run several ``(spec, analysis)`` jobs; one subprocess for all.

    Returns one result dict per job (numpy arrays for vectors).  Raises
    :class:`OpenSeesUnavailable` when OpenSees cannot be run at all.
    """
    jobs = [{"spec": s, "analysis": a} for s, a in jobs]
    if not jobs:
        return []
    use_inproc = inprocess_available() if in_process is None else in_process
    if use_inproc and python is None:
        if not inprocess_available():
            raise OpenSeesUnavailable("openseespy does not import here.")
        res = runner.run_jobs(jobs)
        return [_to_numpy(r) for r in res]
    py = python or find_opensees_python()
    if not py:
        raise OpenSeesUnavailable(
            "No OpenSees interpreter found: install openseespy here "
            "(pip install openseespy) or set the environment variable "
            f"{ENV_VAR} to a Python that has it.")
    if not os.path.isfile(py):
        raise OpenSeesUnavailable(f"OpenSees interpreter {py!r} does not "
                                  f"exist.")
    with tempfile.TemporaryDirectory(prefix="geodome_ops_") as tmp:
        fin = os.path.join(tmp, "job.json")
        fout = os.path.join(tmp, "result.json")
        with open(fin, "w", encoding="utf-8") as fh:
            json.dump({"jobs": jobs}, fh)
        proc = subprocess.run([py, RUNNER_PATH, fin, fout],
                              capture_output=True, text=True,
                              timeout=timeout)
        if proc.returncode != 0 or not os.path.isfile(fout):
            tail = (proc.stderr or proc.stdout or "")[-2000:]
            raise RuntimeError(f"OpenSees runner failed (exit "
                               f"{proc.returncode}) with {py}:\n{tail}")
        with open(fout, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    return [_to_numpy(r) for r in data["results"]]


def run_opensees(spec, analysis="static", **kwargs):
    """Run one analysis of ``spec`` in OpenSees; returns a result dict.

    ``analysis``: "static", or a dict ``{"type": "modal", "n_modes": 12,
    "mass_type": "lumped"}``, ``{"type": "path", "control_dof": d,
    "targets": [...]}``, ``{"type": "arclength", "arc_length": s, ...}``,
    ``{"type": "strain_probe"}`` or ``{"type": "thermal_probe"}`` (see
    ``_opensees_runner``).  Keyword arguments go to
    :func:`run_opensees_batch` (``python``, ``in_process``, ``timeout``).
    In-process when openseespy imports here, else a subprocess.
    """
    return run_opensees_batch([(spec, analysis)], **kwargs)[0]


# ======================================================== comparisons =====
def _scale(*arrays, floor=1e-300):
    s = 0.0
    for a in arrays:
        a = np.asarray(a, dtype=float)
        if a.size:
            s = max(s, float(np.max(np.abs(a))))
    return max(s, floor)


def rel_inf(a, b, scale=None):
    """max |a - b| / scale (default scale = max |a|)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0:
        return 0.0
    s = _scale(a) if scale is None else max(float(scale), 1e-300)
    return float(np.max(np.abs(a - b))) / s


def mac(a, b):
    """Modal assurance criterion of two vectors."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = float(a @ a) * float(b @ b)
    return float((a @ b) ** 2 / den) if den > 0.0 else float("nan")


def subspace_mac(a, B):
    """||Q^T a||^2 / ||a||^2 with Q an orthonormal basis of the columns of B."""
    a = np.asarray(a, dtype=float)
    B = np.asarray(B, dtype=float).reshape(len(a), -1)
    Q, _ = np.linalg.qr(B)
    na = float(a @ a)
    return float(np.sum((Q.T @ a) ** 2) / na) if na > 0.0 else float("nan")


def frequency_clusters(freqs, rtol=1e-7):
    """Groups of indices whose frequencies agree within ``rtol``."""
    f = np.asarray(freqs, dtype=float)
    groups, cur = [], [0] if len(f) else []
    for i in range(1, len(f)):
        if abs(f[i] - f[cur[0]]) <= rtol * max(abs(f[cur[0]]), 1e-300):
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    if cur:
        groups.append(cur)
    return groups


# ------------------------------------------------------------- static ----
def static_job(model, combination_name):
    return (to_opensees_spec(model, combination_name, include_mass=False),
            {"type": "static"})


def compare_static(model, combination_name, os_result, gdl_result=None):
    """GeoDomeLab vs OpenSees for one linear static combination.

    Returns a dict: ``du``, ``dN``, ``dR``, ``dspring`` (normwise relative
    differences, see the module doc), the scales, the largest values, the
    OpenSees free-DOF residual and the GeoDomeLab identities.
    """
    from .solver import solve_model
    from . import loads as ld
    combo = find_combination(model, combination_name)
    if gdl_result is None:
        F = ld.combine(ld.basic_load_vectors(model), combo)
        gdl_result = solve_model(model, F, name=combo.name,
                                 T_factor=float(combo.T))
    g = gdl_result
    if not os_result.get("ok"):
        return {"ok": False, "error": os_result.get("error") or
                "OpenSees analysis failed"}
    uo = np.asarray(os_result["u"])
    No = np.asarray(os_result["N"])
    Ro = np.asarray(os_result["reactions"])
    So = np.asarray(os_result["spring_forces"])
    EA = model.EA()
    sN = _scale(g.N, (EA * np.abs(g.eps0)) if g.eps0 is not None else [])
    sR = _scale(g.reactions, g.F)
    sF = _scale(g.F, g.reactions)
    return {
        "ok": True,
        "du": rel_inf(g.u, uo),
        "dN": rel_inf(g.N, No, sN),
        "dR": rel_inf(g.reactions, Ro, sR),
        "dspring": rel_inf(g.spring_forces, So, sR)
        if np.any(g.spring_forces != 0.0) else 0.0,
        "u_max": _scale(g.u, floor=0.0),
        "N_max": _scale(g.N, floor=0.0),
        "R_max": _scale(g.reactions, floor=0.0),
        "os_residual": float(os_result.get("free_dof_residual", 0.0)) / sF,
        "gdl_residual": float(g.equilibrium_residual),
        "gdl_energy_error": float(g.energy_error),
        "gdl_balance": float(np.max(np.abs(g.reaction_balance))) / sF,
        "system": os_result.get("system", ""),
        "gdl": g, "os": os_result,
    }


def crosscheck_static(model, combination_name, **run_kwargs):
    spec, an = static_job(model, combination_name)
    return compare_static(model, combination_name,
                          run_opensees(spec, an, **run_kwargs))


# -------------------------------------------------------------- modal ----
def modal_job(model, n_modes=12, mass_type="lumped", include_dead=True,
              extra=4):
    """OpenSees modal job; ``extra`` more modes than compared (clusters)."""
    spec = to_opensees_spec(model, None, include_mass=True,
                            include_dead=include_dead)
    k = None if n_modes in (None, "all") else int(n_modes) + int(extra)
    return spec, {"type": "modal", "n_modes": k if k is not None else "all",
                  "mass_type": mass_type}


def compare_modal(model, os_result, n_modes=12, mass_type="lumped",
                  include_dead=True, gdl_result=None, cluster_rtol=1e-7):
    """GeoDomeLab vs OpenSees natural frequencies and mode shapes.

    Returns ``df`` (max pointwise relative frequency difference over the
    first ``n_modes``), ``mac_min`` (smallest MAC of non-repeated modes),
    ``submac_min`` (smallest subspace MAC inside repeated clusters), the
    frequencies of both programs and the number of modes compared.
    """
    from . import modal as md
    if not os_result.get("ok"):
        return {"ok": False, "error": os_result.get("error") or
                "OpenSees eigen analysis failed"}
    fo = np.sort(np.asarray(os_result["frequencies"]))
    order = np.argsort(np.asarray(os_result["frequencies"]))
    Vo = np.asarray(os_result["modes"])[:, order] \
        if np.asarray(os_result["modes"]).size else None
    k_all = len(fo)
    if gdl_result is None:
        gdl_result = md.modal_analysis(model, n_modes=k_all,
                                       include_dead=include_dead,
                                       mass_type=mass_type,
                                       align_degenerate=False)
    fg = np.asarray(gdl_result.frequencies)
    Vg = np.asarray(gdl_result.modes)
    k_all = min(k_all, len(fg))
    k = k_all if n_modes in (None, "all") else min(int(n_modes), k_all)
    df_each = np.abs(fg[:k] - fo[:k]) / np.maximum(np.abs(fg[:k]), 1e-300)
    mac_min, sub_min, n_single, n_clustered = 1.0, 1.0, 0, 0
    if Vo is not None:
        for grp in frequency_clusters(fg[:k_all], cluster_rtol):
            grp_in = [i for i in grp if i < k]
            if not grp_in:
                continue
            if len(grp) == 1:
                i = grp[0]
                mac_min = min(mac_min, mac(Vg[:, i], Vo[:, i]))
                n_single += 1
            elif grp[-1] < k_all:
                for i in grp_in:
                    sub_min = min(sub_min, subspace_mac(Vg[:, i], Vo[:, grp]))
                    n_clustered += 1
    return {"ok": True, "df": float(np.max(df_each)) if k else 0.0,
            "df_each": df_each, "mac_min": mac_min, "submac_min": sub_min,
            "n_compared": k, "n_single": n_single,
            "n_clustered": n_clustered, "f_gdl": fg[:k], "f_os": fo[:k],
            "f1": float(fg[0]) if k else float("nan"),
            "solver": os_result.get("solver", ""),
            "gdl": gdl_result, "os": os_result}


def crosscheck_modal(model, n_modes=12, mass_type="lumped", **run_kwargs):
    spec, an = modal_job(model, n_modes, mass_type)
    return compare_modal(model, run_opensees(spec, an, **run_kwargs),
                         n_modes=n_modes, mass_type=mass_type)


# ------------------------------------------------------- nonlinear path ---
def monotone_targets(control, min_step=0.0):
    """Leading strictly monotone part of a control-displacement history.

    Returns the indices (into ``control``) of the points after the start
    that keep moving in the initial direction; stops at the first reversal
    (a snap-back in the control DOF, where displacement control fails).
    """
    c = np.asarray(control, dtype=float)
    idx = []
    if len(c) < 2:
        return idx
    sgn = 0.0
    last = c[0]
    for i in range(1, len(c)):
        d = c[i] - last
        if sgn == 0.0:
            if abs(d) <= min_step:
                continue
            sgn = math.copysign(1.0, d)
        if d * sgn <= min_step:
            if d * sgn < 0.0:
                break
            continue
        idx.append(i)
        last = c[i]
    return idx


def path_segments(U, free, preferred=None, min_frac=0.05):
    """Split a traced path into segments of piecewise displacement control.

    ``U`` (k, 3n): the displacement states of the path (row 0 = start).
    Greedy: from the last matched point, every free DOF whose increment
    keeps one sign and is at least ``min_frac`` of the largest component of
    the step (so that it is a well-conditioned control) is a candidate; the
    one that stays monotone for the most points wins (``preferred``, the
    user's control DOF, wins ties).  Repeated states (zero step) are
    skipped.  Returns ``[{"control_dof", "targets", "index"}]`` where
    ``index`` are the rows of ``U`` reached by that segment.
    """
    U = np.asarray(U, dtype=float)
    free = np.asarray(free, dtype=int)
    k = len(U)
    keep = [0]
    scale = max(float(np.max(np.abs(U))) if U.size else 0.0, 1e-300)
    for i in range(1, k):
        if np.max(np.abs(U[i] - U[keep[-1]])) > 1e-12 * scale:
            keep.append(i)
    rows = np.asarray(keep)
    D = np.diff(U[rows][:, free], axis=0)              # (p-1, nf)
    big = np.max(np.abs(D), axis=1) if len(D) else np.zeros(0)
    segs = []
    s = 0
    pref = None
    if preferred is not None and int(preferred) in set(free.tolist()):
        pref = int(np.where(free == int(preferred))[0][0])
    while s < len(D):
        sgn = np.sign(D[s])
        ok = (D[s:] * sgn[None, :]) >= min_frac * big[s:, None]
        ok &= sgn[None, :] != 0.0
        run = np.where(ok.all(axis=0), ok.shape[0], np.argmin(ok, axis=0))
        best = int(np.argmax(run))
        if pref is not None and run[pref] >= run[best]:
            best = pref
        r = int(run[best])
        if r == 0:                                    # pragma: no cover
            best = int(np.argmax(np.abs(D[s])))
            r = 1
        idx = rows[s + 1:s + 1 + r]
        d = int(free[best])
        segs.append({"control_dof": d,
                     "targets": [float(U[i, d]) for i in idx],
                     "index": [int(i) for i in idx]})
        s += r
    return segs


def gdl_residual(system, F_ref, u, lam):
    """GeoDomeLab equilibrium residual of an OpenSees state (u, lambda).

    ``|lambda F_ref - f(u)|_free / max(|F_ref|, |lambda F_ref|)`` with the
    corotational internal forces of ``nonlinear.TrussSystem``.  It is at
    round-off level only if both programs use the same strain measure and
    the same kinematics.
    """
    st = system.state(np.asarray(u, dtype=float))
    Ff = np.asarray(F_ref, dtype=float)[system.free]
    r = lam * Ff - st.f_free
    return float(np.linalg.norm(r) / max(np.linalg.norm(Ff),
                                         abs(lam) * np.linalg.norm(Ff),
                                         1e-300))


def _free_dofs(model):
    from .nonlinear import model_system
    return model_system(model).free


def path_jobs(model, gdl_result, *, arclength=None, tol=1e-10,
              include_limit=True, min_frac=0.05):
    """OpenSees jobs that follow a GeoDomeLab nonlinear result.

    The load is the D/S/L/W reference load of ``gdl_result.combination``
    (T excluded, as in ``nonlinear.run_nonlinear``).

    1. ``path``: piecewise DisplacementControl (:func:`path_segments`)
       through every point of the GeoDomeLab path, with the refined first
       limit point inserted in order, so that OpenSees lands on the control
       value of every GeoDomeLab state exactly -- through limit points and
       snap-backs.
    2. optionally ``arclength`` (dict of runner options): an independent
       OpenSees ArcLength trace.

    Returns ``(jobs, segments, rows)``; ``rows[j]`` is the GeoDomeLab state
    of the j-th OpenSees point (a path index, or "limit").
    """
    spec = to_opensees_spec(model, gdl_result.combination,
                            include_mass=False)
    spec["eps0"] = None                  # T is not part of the nonlinear path
    spec["prescribed"] = []
    r = gdl_result
    U = np.asarray(r.u, dtype=float)
    labels = list(range(len(U)))
    if include_limit and r.limit_found and r.u_limit is not None:
        # the limit state lies between path points; put it between the two
        # neighbours that bracket its control displacement
        ul = np.asarray(r.u_limit, dtype=float)
        c = U[:, r.control_dof]
        cl = ul[r.control_dof]
        pos = len(U)
        for i in range(1, len(U)):
            if (c[i - 1] - cl) * (c[i] - cl) <= 0.0:
                pos = i
                break
        U = np.insert(U, pos, ul, axis=0)
        labels.insert(pos, "limit")
    segs = path_segments(U, _free_dofs(model), preferred=r.control_dof,
                         min_frac=min_frac)
    rows = [labels[i] for sgm in segs for i in sgm["index"]]
    opts = {"type": "path", "tol": tol, "max_iter": 60,
            "segments": [{"control_dof": sgm["control_dof"],
                          "targets": sgm["targets"]} for sgm in segs]}
    jobs = [(spec, opts)]
    if arclength:
        a = {"type": "arclength", "control_dof": int(r.control_dof),
             "tol": tol, "max_iter": 60}
        a.update(arclength)
        jobs.append((spec, a))
    return jobs, segs, rows


def compare_path(model, gdl_result, os_path, rows, os_arc=None):
    """Compare GeoDomeLab's path with OpenSees at matched states.

    ``rows`` comes from :func:`path_jobs`.  Returns ``dlam`` (max
    |lambda_GDL - lambda_OS| / max |lambda| on the path), ``dlam_point``
    (pointwise relative where |lambda| > 5 % of that maximum), ``du`` (max
    over points of |u_GDL - u_OS|_inf / max |u_GDL|), ``dlimit`` (relative
    difference of lambda at the refined limit state), the GeoDomeLab
    residual of the OpenSees states (``res_os_in_gdl``), the number of
    control segments, and the residual of the arc-length points when
    ``os_arc`` is given.
    """
    from .nonlinear import model_system, reference_load
    from . import nonlinear as nl
    r = gdl_result
    combo = nl._find_combination(model, r.combination)
    F = reference_load(model, combo)
    system = model_system(model)
    lo = np.asarray(os_path.get("lam", []), dtype=float)
    uo = os_path.get("u", [])
    uo = np.asarray(uo, dtype=float) if len(uo) else np.zeros((0, len(F)))
    lam_scale = _scale(r.lam)
    u_scale = _scale(r.u)
    dl = dlp = du = ctrl_err = 0.0
    lim_os = float("nan")
    matched = 0
    pairs = []
    dofs = list(os_path.get("control_dofs", []))
    for j in range(len(lo)):
        lab = rows[j]
        if lab == "limit":
            ug, lg = np.asarray(r.u_limit, dtype=float), float(r.lambda_limit)
            lim_os = float(lo[j])
        else:
            ug, lg = np.asarray(r.u[lab], dtype=float), float(r.lam[lab])
        matched += 1
        if dofs:
            d = int(dofs[j])
            ctrl_err = max(ctrl_err, abs(uo[j][d] - ug[d]) / u_scale)
        dl = max(dl, abs(lg - lo[j]) / lam_scale)
        if abs(lg) > 0.05 * lam_scale:
            dlp = max(dlp, abs(lg - lo[j]) / abs(lg))
        du = max(du, float(np.max(np.abs(ug - uo[j]))) / u_scale)
        pairs.append((lab, float(ug[r.control_dof]), lg, float(lo[j])))
    res_os = 0.0
    for j in range(len(lo)):
        res_os = max(res_os, gdl_residual(system, F, uo[j], lo[j]))
    seg = os_path.get("segment", [])
    out = {
        "ok": bool(os_path.get("ok")) and matched == len(rows),
        "message": os_path.get("message", "") or os_path.get("error", ""),
        "dlam": dl, "dlam_point": dlp, "du": du, "control_error": ctrl_err,
        "n_targets": len(rows), "n_matched": matched,
        "n_segments": len(set(int(v) for v in seg)) if len(seg) else 0,
        "control_dofs": sorted(set(int(d) for d in dofs)),
        "lambda_scale": lam_scale, "pairs": pairs,
        "lambda_limit_gdl": float(r.lambda_limit) if r.limit_found
        else float("nan"),
        "control_limit_gdl": float(r.control_limit) if r.limit_found
        else float("nan"),
        "lambda_limit_os": lim_os,
        "dlimit": (abs(r.lambda_limit - lim_os) / abs(r.lambda_limit)
                   if r.limit_found and np.isfinite(lim_os)
                   else float("nan")),
        "res_os_in_gdl": res_os,
        "control_range": (float(np.min(r.control)), float(np.max(r.control)))
        if len(r.control) else (0.0, 0.0),
        "lambda_range": (float(np.min(r.lam)), float(np.max(r.lam)))
        if len(r.lam) else (0.0, 0.0),
    }
    if os_arc is not None:
        la = np.asarray(os_arc.get("lam", []), dtype=float)
        ua = os_arc.get("u", [])
        ra = 0.0
        for j in range(len(la)):
            if len(ua):
                ra = max(ra, gdl_residual(system, F, ua[j], la[j]))
        ca = np.asarray(os_arc.get("control", []), dtype=float)
        out.update({"arc_ok": bool(os_arc.get("ok")),
                    "arc_points": int(len(la)),
                    "arc_message": os_arc.get("message", ""),
                    "res_arc_in_gdl": ra, "arc_lam": la, "arc_control": ca,
                    "arc_lam_max": float(np.max(la)) if len(la)
                    else float("nan")})
    return out


def crosscheck_path(model, combination_name=None, *, gdl_kwargs=None,
                    arclength=None, **run_kwargs):
    """Run GeoDomeLab's nonlinear path and follow it in OpenSees.

    ``gdl_kwargs`` go to ``nonlinear.run_nonlinear`` (method, control,
    control_max, lambda_max, ...).  Returns the :func:`compare_path` dict
    plus the two results (``gdl``, ``os``).
    """
    from .nonlinear import run_nonlinear
    kw = dict(compare=False, member_checks=False)
    kw.update(gdl_kwargs or {})
    r = run_nonlinear(model, combination_name, **kw)
    jobs, _segs, rows = path_jobs(model, r, arclength=arclength)
    res = run_opensees_batch(jobs, **run_kwargs)
    out = compare_path(model, r, res[0], rows,
                       res[1] if len(res) > 1 else None)
    out["gdl"] = r
    out["os"] = res[0]
    return out


# -------------------------------------------------- formulation probes ----
def strain_measure_check(rot=0.8, load=0.1, **run_kwargs):
    """Identify the strain measure of OpenSees ``corotTruss`` numerically.

    Returns ``{"N", "l_over_L", "engineering", "green", "log",
    "rel_engineering", "rel_green", "rel_log", "measure"}``; ``measure`` is
    the prediction that matches to 1e-9 (or "unknown").
    """
    r = run_opensees(None, {"type": "strain_probe", "rot": rot,
                            "load": load}, **run_kwargs)
    if not r.get("ok"):
        return {"ok": False, "error": r.get("error", "probe failed")}
    L, EA, l, N = r["L"], r["EA"], r["l"], r["N"]
    pred = {"engineering": EA * (l - L) / L,
            "green": EA * (l * l - L * L) / (2.0 * L * L),
            "log": EA * math.log(l / L)}
    rel = {k: abs(N - v) / abs(v) for k, v in pred.items()}
    best = min(rel, key=rel.get)
    eq = abs(N * r["x"] / l - r["Px"]) / r["Px"]
    return {"ok": True, "N": N, "l_over_L": l / L,
            "rotation_deg": math.degrees(math.atan2(r["y"], r["x"])),
            **pred, **{f"rel_{k}": v for k, v in rel.items()},
            "measure": best if rel[best] < 1e-9 else "unknown",
            "equilibrium": eq, "opensees_version":
            r.get("opensees_version", "")}


def thermal_sign_check(eps0=11.7e-6 * 30.0, ratio=2.0, **run_kwargs):
    """Heated bar between a wall and a second bar (InitStrainMaterial sign).

    Closed form (bar 1 heated, bar 2 with ``ratio`` times the area, both
    ends fixed): u2 = eps0 L / (1 + r), N1 = N2 = -EA eps0 r / (1 + r).
    Returns the OpenSees values, the closed form and the relative errors.
    """
    r = run_opensees(None, {"type": "thermal_probe", "eps0": eps0,
                            "ratio": ratio}, **run_kwargs)
    if not r.get("ok"):
        return {"ok": False, "error": r.get("error", "probe failed")}
    L, EA, rr = r["L"], r["EA"], r["ratio"]
    u2 = eps0 * L / (1.0 + rr)
    N = -EA * eps0 * rr / (1.0 + rr)
    return {"ok": True, "u2_os": r["u2"], "u2_exact": u2,
            "N1_os": r["N1"], "N2_os": r["N2"], "N_exact": N,
            "R1_os": r["R1"], "R3_os": r["R3"],
            "rel_u": abs(r["u2"] - u2) / abs(u2),
            "rel_N": max(abs(r["N1"] - N), abs(r["N2"] - N)) / abs(N)}


# ===================================================== script export ======
_SCRIPT_HEAD = '''"""
OpenSeesPy model exported by GeoDomeLab {version}.

Structure:    {title}
Combination:  {combo}
Analysis:     {analysis}
Units:        N, m, kg, s  (DOF 1 = x, 2 = y, 3 = z; tension positive)

Run it with any Python that has OpenSeesPy (pip install openseespy):

    python {fname}                  print the results
    python {fname} results.json     also write them as JSON

Node tag = GeoDomeLab node + 1, element tag = GeoDomeLab member + 1.
Ground nodes of spring supports (if any) have tags above the last node and
are linked to the structure by zeroLength springs.  Thermal strain eps0 is
modelled with InitStrainMaterial(-eps0) so that N = EA (eps - eps0);
settlements are sp constraints inside the load pattern.
"""

import json
import math
import sys

import openseespy.opensees as ops

N_NODES = {n}
N_MEMBERS = {m}
'''


def _script_analysis_block(spec, analysis):
    a = runner._normalise(analysis)
    kind = a.get("type", "static")
    n = len(spec["nodes"])
    lines = []
    if kind == "static":
        _f, _s, springs, ground = runner._support_layout(spec)
        fixed = sorted(int(d) for d in spec.get("fixed", []))
        spring_dofs = [(node, ax) for node, lst in sorted(springs.items())
                       for ax, _k in lst]
        lines += [
            "# ---- linear static analysis to load factor 1 ----",
            "ops.constraints('Transformation')",
            "ops.numberer('RCM')",
            f"ops.system({'FullGeneral' if 3 * n <= 900 else 'UmfPack'!r})",
            "ops.test('NormUnbalance', 1e-6, 10)",
            "ops.algorithm('Linear')",
            "ops.integrator('LoadControl', 1.0)",
            "ops.analysis('Static')",
            "assert ops.analyze(1) == 0, 'analysis failed'",
            "ops.reactions()",
            "",
            f"FIXED_DOFS = {fixed!r}          # restrained GeoDomeLab DOFs",
            f"SPRINGS = {[(nd, ax, ground[nd]) for nd, ax in spring_dofs]!r}"
            "   # (node, axis, ground tag)",
            "u = [ops.nodeDisp(i + 1, a + 1) for i in range(N_NODES) "
            "for a in range(3)]",
            "N = [ops.eleResponse(k + 1, 'axialForce')[0] "
            "for k in range(N_MEMBERS)]",
            "R = [0.0] * (3 * N_NODES)",
            "for d in FIXED_DOFS:",
            "    R[d] = ops.nodeReaction(d // 3 + 1, d % 3 + 1)",
            "for node, ax, gtag in SPRINGS:",
            "    R[3 * node + ax] += ops.nodeReaction(gtag, ax + 1)",
            "RESULTS = {'type': 'static', 'u': u, 'N': N, 'reactions': R}",
            "",
            "print('node      ux [mm]       uy [mm]       uz [mm]')",
            "for i in range(N_NODES):",
            "    print(f'{i:4d} {u[3*i]*1e3:13.6f} {u[3*i+1]*1e3:13.6f} "
            "{u[3*i+2]*1e3:13.6f}')",
            "print('member    N [kN] (tension +)')",
            "for k in range(N_MEMBERS):",
            "    print(f'{k:4d} {N[k]/1e3:13.6f}')",
            "print('reactions [kN]: ' + ', '.join(f'DOF {d}: {R[d]/1e3:.6f}' "
            "for d in range(3 * N_NODES) if R[d] != 0.0))",
        ]
    elif kind == "modal":
        k = a.get("n_modes", 12)
        nfree = 3 * n - len(set(spec.get("fixed", [])))
        solver = a.get("solver")
        if k in (None, "all") or int(k) >= nfree - 1:
            k, solver = nfree, solver or "-fullGenLapack"
        solver = solver or "-genBandArpack"
        lines += [
            "# ---- natural frequencies (eigen returns omega^2) ----",
            f"omega2 = ops.eigen({solver!r}, {int(k)})",
            "freq = [math.sqrt(max(w2, 0.0)) / (2.0 * math.pi) "
            "for w2 in omega2]",
            "modes = [[ops.nodeEigenvector(i + 1, j + 1, a + 1) "
            "for i in range(N_NODES) for a in range(3)] "
            "for j in range(len(omega2))]",
            "RESULTS = {'type': 'modal', 'omega2': list(omega2), "
            "'frequencies': freq, 'modes': modes}",
            "print('mode   f [Hz]          T [s]')",
            "for j, f in enumerate(freq):",
            "    print(f'{j + 1:4d} {f:14.8f} {1.0 / f if f > 0 else 0.0:14.8f}')",
        ]
    elif kind == "path":
        cdof = int(a["control_dof"])
        targets = a.get("targets")
        if not targets:
            ns = int(a.get("n_steps", 50))
            cmax = float(a.get("control_max", 0.1))
            targets = [-cmax * (i + 1) / ns for i in range(ns)]
        Fn = math.sqrt(sum(float(v) ** 2 for v in spec["load"]))
        tol = float(a.get("tol", 1e-9)) * max(Fn, 1e-30)
        lines += [
            "# ---- geometrically nonlinear path: corotTruss + "
            "DisplacementControl ----",
            f"CONTROL_NODE, CONTROL_DOF = {cdof // 3 + 1}, {cdof % 3 + 1}"
            f"   # GeoDomeLab DOF {cdof} = node {cdof // 3} "
            f"{_AXIS[cdof % 3]}",
            f"TARGETS = {[float(t) for t in targets]!r}",
            "ops.constraints('Transformation')",
            "ops.numberer('RCM')",
            f"ops.system({'FullGeneral' if 3 * n <= 900 else 'UmfPack'!r})",
            f"ops.test('NormUnbalance', {tol!r}, 60)",
            "ops.algorithm('Newton')",
            "ops.integrator('DisplacementControl', CONTROL_NODE, "
            "CONTROL_DOF, 0.0)",
            "ops.analysis('Static')",
            "lam, ctrl, U = [], [], []",
            "current = 0.0",
            "for w in TARGETS:",
            "    nsub = 1",
            "    while nsub > 0:                  # sub-step on failure",
            "        ops.integrator('DisplacementControl', CONTROL_NODE, "
            "CONTROL_DOF, (w - current) / nsub)",
            "        if ops.analyze(1) == 0:",
            "            current = ops.nodeDisp(CONTROL_NODE, CONTROL_DOF)",
            "            nsub -= 1",
            "        else:",
            "            nsub *= 2",
            "            if nsub > 256:",
            "                raise SystemExit(f'no convergence before {w}')",
            "    lam.append(ops.getLoadFactor(1))",
            "    ctrl.append(current)",
            "    U.append([ops.nodeDisp(i + 1, a + 1) for i in range(N_NODES) "
            "for a in range(3)])",
            "RESULTS = {'type': 'path', 'lam': lam, 'control': ctrl, 'u': U}",
            "print('control [mm]     load factor')",
            "for c, l in zip(ctrl, lam):",
            "    print(f'{c * 1e3:12.6f} {l:16.10f}')",
        ]
    else:
        raise ValueError(f"write_openseespy_script: analysis {kind!r} is not "
                         f"supported (static, modal, path).")
    lines += ["", "if len(sys.argv) > 1:",
              "    with open(sys.argv[1], 'w') as fh:",
              "        json.dump(RESULTS, fh)", ""]
    return lines


def write_openseespy_script(model, combination_name, path, analysis="static"):
    """Write a standalone OpenSeesPy script of ``model`` and return ``path``.

    ``analysis``: "static" (linear static of the combination, with T),
    "modal" or ``{"type": "modal", "n_modes": k, "mass_type": ...}``, or
    "nonlinear" / ``{"type": "path", "control_dof": d, "targets": [...]}``
    (corotTruss with DisplacementControl; by default the crown uz in 50
    steps to 10 % of ``geom.length_scale``, the sphere radius of a
    geodesic dome).  The script prints the results and, given a
    file name argument, writes them as JSON.
    """
    from . import __version__
    a = runner._normalise(analysis)
    if a.get("type") == "nonlinear":
        a["type"] = "path"
    kind = a.get("type", "static")
    if kind == "path":
        from .nonlinear import resolve_control
        a.setdefault("control_dof",
                     resolve_control(model, a.pop("control", "crown_uz")))
        a.setdefault("control_max", 0.1 * float(model.geom.length_scale))
    mass_type = a.get("mass_type", "lumped") if kind == "modal" else None
    spec = to_opensees_spec(model, combination_name if kind != "modal"
                            else None, include_mass=(kind == "modal"))
    if kind == "path":
        spec["eps0"], spec["prescribed"] = None, []
    element = "corotTruss" if kind == "path" else "Truss"
    cmds = runner.model_commands(spec, element, mass_type=mass_type)
    if kind != "modal":
        cmds += runner.load_commands(spec, with_sp=(kind == "static"))
    fname = os.path.basename(path)
    head = _SCRIPT_HEAD.format(
        version=__version__, title=spec["title"],
        combo=spec["combination"] or "none (modal analysis)",
        analysis={"static": "linear static (Truss)",
                  "modal": f"modal ({mass_type} mass)",
                  "path": "geometrically nonlinear (corotTruss)"}[kind],
        fname=fname, n=len(spec["nodes"]), m=len(spec["members"]))
    body = ["", "# ---- model ----"]
    body += [runner.format_command(nm, args) for nm, args in cmds]
    body += [""] + _script_analysis_block(spec, a)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(head + "\n".join(body))
    return path


def run_script(path, python=None, timeout=900):
    """Run an exported script with an OpenSees interpreter; returns RESULTS."""
    py = python or (sys.executable if inprocess_available()
                    else find_opensees_python())
    if not py:
        raise OpenSeesUnavailable("No OpenSees interpreter found.")
    with tempfile.TemporaryDirectory(prefix="geodome_ops_") as tmp:
        out = os.path.join(tmp, "results.json")
        proc = subprocess.run([py, path, out], capture_output=True,
                              text=True, timeout=timeout)
        if proc.returncode != 0 or not os.path.isfile(out):
            raise RuntimeError(f"Script failed (exit {proc.returncode}):\n"
                               f"{(proc.stderr or proc.stdout)[-2000:]}")
        with open(out, "r", encoding="utf-8") as fh:
            return json.load(fh)


__all__ = [
    "ENV_VAR", "OpenSeesUnavailable", "inprocess_available",
    "opensees_python_candidates", "find_opensees_python", "opensees_status",
    "opensees_available", "find_combination", "to_opensees_spec",
    "run_opensees", "run_opensees_batch", "rel_inf", "mac", "subspace_mac",
    "frequency_clusters", "static_job", "compare_static",
    "crosscheck_static", "modal_job", "compare_modal", "crosscheck_modal",
    "monotone_targets", "path_segments", "gdl_residual", "path_jobs",
    "compare_path", "crosscheck_path",
    "strain_measure_check", "thermal_sign_check", "write_openseespy_script",
    "run_script",
]
