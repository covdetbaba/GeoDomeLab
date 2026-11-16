"""
OpenSees (OpenSeesPy) runner for GeoDomeLab cross-checks.

This file is deliberately self-contained: it imports only the Python
standard library and ``openseespy``, never ``geodome`` or numpy.  It can
therefore run

* in-process, when ``openseespy`` imports in the GeoDomeLab interpreter, or
* in a separate OpenSees interpreter as a script::

      python _opensees_runner.py job.json result.json

  where ``job.json`` holds ``{"jobs": [{"spec": {...}, "analysis": {...}},
  ...]}`` and ``result.json`` receives ``{"results": [...]}``.

``geodome.opensees_bridge`` writes the spec (a plain dict, see
``to_opensees_spec``) and calls this runner.  The OpenSees model is built
from a list of ``(command, args)`` tuples (:func:`model_commands`), so the
same list can be executed here or printed as a standalone OpenSeesPy script
(``opensees_bridge.write_openseespy_script``).

OpenSees modelling of each GeoDomeLab ingredient
------------------------------------------------
* Node ``i`` -> tag ``i + 1``; member ``m`` -> element tag ``m + 1``.
  ``model('basic', '-ndm', 3, '-ndf', 3)``.
* Strut: ``element('Truss', ...)`` (small displacements, linear and modal
  analyses) or ``element('corotTruss', ...)`` (geometrically nonlinear),
  with ``uniaxialMaterial('Elastic', E)``.  The corotTruss strain is
  engineering strain (l - L)/L on the original length (CorotTruss.cpp,
  M. H. Scott); the bridge confirms it numerically (``strain_probe``).
* Initial (thermal) strain eps0: ``uniaxialMaterial('InitStrainMaterial',
  tag, elastic, -eps0)``.  InitStrainMaterial evaluates the wrapped material
  at ``strain + epsInit``, so ``epsInit = -eps0`` gives
  N = EA (eps - eps0) (InitStrainMaterial.cpp, F. McKenna); the sign is
  checked on a heated fixed bar.  It is present from the first step and is
  not scaled by the load factor; the linear static solve is one step to
  load factor 1, so this is immaterial there.
* Rigid restraint: ``fix``.  Prescribed movement of a restrained DOF
  (settlement): that DOF is not fixed; ``sp(node, dof, value)`` inside the
  load pattern imposes ``value`` x load factor.
* Grounded spring k on DOF d of node i: a ground node at the same point,
  fully fixed, and ``element('zeroLength', ground, node, '-mat', ...,
  '-dir', ...)`` with an Elastic material of stiffness k.  A movement of the
  spring's ground end is an ``sp`` on the ground node.  The spring force on
  the structure equals the reaction at the ground node.
* Mass: strut mass ``'-rho', rho*A`` (mass per unit length; lumped
  rho A L / 2 per end and direction, or consistent with ``'-cMass', 1``) and
  non-structural nodal masses with ``mass(node, m, m, m)``.

Analyses (``analysis["type"]``)
-------------------------------
``static``      linear static solve to load factor 1 (``algorithm Linear``)
                -> u, N, reactions, spring forces, residual at free DOFs.
``modal``       ``eigen`` (ARPACK ``-genBandArpack`` or ``-fullGenLapack``)
                -> omega^2, frequencies, mode shapes.
``path``        corotTruss, ``DisplacementControl`` to a list of target
                control displacements (each target reached exactly, with
                sub-stepping on failure) -> lambda, u, N at every target.
``arclength``   corotTruss, ``ArcLength`` (alpha = 0: cylindrical) or
                ``MinUnbalDispNorm`` stepping with step halving on failure
                -> the traced (lambda, u) points.
``strain_probe`` one corotTruss bar with a large prescribed rotation and an
                axial load -> current length and axial force, to identify the
                strain measure.

References: F. McKenna, M. H. Scott, G. L. Fenves (2010), J. Comput. Civ.
Eng. 24(1), 95-107; M. Zhu, F. McKenna, M. H. Scott (2018), SoftwareX 7,
6-11; OpenSees source files CorotTruss.cpp, InitStrainMaterial.cpp,
ArcLength.cpp, DisplacementControl.cpp.
"""

from __future__ import annotations

import json
import math
import sys
import time

AXES = ("ux", "uy", "uz")
PATTERN = 1
TSERIES = 1


def _ops():
    import openseespy.opensees as ops                # noqa: WPS433
    return ops


def opensees_version():
    try:
        v = _ops().version()
        return str(v) if v is not None else "unknown"
    except Exception:                                 # pragma: no cover
        return "unknown"


# ======================================================= model commands ===
def _support_layout(spec):
    """Fixed flags, sp list, springs and ground nodes from a spec.

    Returns ``(fix_flags {node: [0/1]*3}, sps [(tag, dof1, value)],
    springs {node: [(axis, k)]}, ground {node: ground_tag})``.
    """
    n = len(spec["nodes"])
    fixed = set(int(d) for d in spec.get("fixed", []))
    presc = {}
    for d, v in spec.get("prescribed", []) or []:
        d = int(d)
        presc[d] = presc.get(d, 0.0) + float(v)
    springs = {}
    for d, k in spec.get("springs", []) or []:
        d = int(d)
        if d in fixed:
            continue                     # the rigid restraint governs
        springs.setdefault(d // 3, []).append((d % 3, float(k)))

    fix_flags, sps = {}, []
    for d in sorted(fixed):
        node = d // 3
        flags = fix_flags.setdefault(node, [0, 0, 0])
        if presc.get(d, 0.0) != 0.0:
            sps.append((node + 1, d % 3 + 1, presc[d]))
        else:
            flags[d % 3] = 1
    ground = {}
    tag = n + 1
    for node in sorted(springs):
        ground[node] = tag
        tag += 1
    for node, lst in springs.items():
        gflags = [1, 1, 1]
        for ax, _k in lst:
            dl = presc.get(3 * node + ax, 0.0)
            if dl != 0.0:
                gflags[ax] = 0
                sps.append((ground[node], ax + 1, dl))
        fix_flags[("g", node)] = gflags
    for d, v in presc.items():
        if v != 0.0 and d not in fixed and (d // 3 not in springs or
                                             d % 3 not in
                                             [a for a, _ in springs[d // 3]]):
            raise ValueError(f"Prescribed value at DOF {d} which is neither "
                             f"restrained nor a spring support.")
    return fix_flags, sps, springs, ground


def model_commands(spec, element="Truss", mass_type=None):
    """List of ``(command, args)`` that builds the OpenSees model of ``spec``.

    ``element``: "Truss" (small displacement) or "corotTruss".
    ``mass_type``: None (no mass), "lumped" or "consistent".
    """
    nodes = spec["nodes"]
    members = spec["members"]
    n, m = len(nodes), len(members)
    E = spec["E"]
    A = spec["A"]
    eps0 = spec.get("eps0")
    cmds = [("wipe", ()), ("model", ("basic", "-ndm", 3, "-ndf", 3))]
    for i, (x, y, z) in enumerate(nodes):
        cmds.append(("node", (i + 1, float(x), float(y), float(z))))
    fix_flags, _sps, springs, ground = _support_layout(spec)
    for node, gtag in sorted(ground.items()):
        x, y, z = nodes[node]
        cmds.append(("node", (gtag, float(x), float(y), float(z))))
    for key in sorted(fix_flags, key=lambda k: (isinstance(k, tuple),
                                                k if not isinstance(k, tuple)
                                                else k[1])):
        flags = fix_flags[key]
        if not any(flags):
            continue
        tag = ground[key[1]] if isinstance(key, tuple) else key + 1
        cmds.append(("fix", (tag, *flags)))

    # materials: one Elastic per distinct E, one InitStrain per (E, eps0)
    mat_of = {}
    next_mat = 1
    el_tag = {}
    for k in range(m):
        e = float(E[k])
        if e not in el_tag:
            el_tag[e] = next_mat
            cmds.append(("uniaxialMaterial", ("Elastic", next_mat, e)))
            next_mat += 1
    for k in range(m):
        e = float(E[k])
        e0 = float(eps0[k]) if eps0 is not None else 0.0
        key = (e, e0)
        if key in mat_of:
            continue
        if e0 == 0.0:
            mat_of[key] = el_tag[e]
        else:
            mat_of[key] = next_mat
            cmds.append(("uniaxialMaterial",
                         ("InitStrainMaterial", next_mat, el_tag[e], -e0)))
            next_mat += 1

    rho = spec.get("line_mass")
    for k, (i, j) in enumerate(members):
        e0 = float(eps0[k]) if eps0 is not None else 0.0
        args = [element, k + 1, int(i) + 1, int(j) + 1, float(A[k]),
                mat_of[(float(E[k]), e0)]]
        if mass_type and rho is not None and float(rho[k]) > 0.0:
            args += ["-rho", float(rho[k])]
            if mass_type == "consistent":
                args += ["-cMass", 1]
        cmds.append(("element", tuple(args)))

    # grounded springs
    etag = m + 1
    for node in sorted(springs):
        mats, dirs = [], []
        for ax, kk in springs[node]:
            cmds.append(("uniaxialMaterial", ("Elastic", next_mat, kk)))
            mats.append(next_mat)
            dirs.append(ax + 1)
            next_mat += 1
        cmds.append(("element", ("zeroLength", etag, ground[node], node + 1,
                                 "-mat", *mats, "-dir", *dirs)))
        etag += 1

    if mass_type:
        nm = spec.get("nodal_mass") or [0.0] * n
        for i, mi in enumerate(nm):
            if float(mi) > 0.0:
                cmds.append(("mass", (i + 1, float(mi), float(mi),
                                      float(mi))))
    return cmds


def load_commands(spec, scale=1.0, with_sp=True):
    """Time series, load pattern, nodal loads and sp settlements."""
    cmds = [("timeSeries", ("Linear", TSERIES)),
            ("pattern", ("Plain", PATTERN, TSERIES))]
    F = spec.get("load") or []
    for i in range(len(spec["nodes"])):
        f = [float(F[3 * i + a]) * scale if F else 0.0 for a in range(3)]
        if any(v != 0.0 for v in f):
            cmds.append(("load", (i + 1, *f)))
    if with_sp:
        _f, sps, _s, _g = _support_layout(spec)
        for tag, dof, val in sps:
            cmds.append(("sp", (tag, dof, float(val))))
    return cmds


def apply_commands(ops, cmds):
    for name, args in cmds:
        getattr(ops, name)(*args)


def format_command(name, args):
    """One command as an OpenSeesPy source line."""
    return f"ops.{name}({', '.join(repr(a) for a in args)})"


# ============================================================== helpers ===
def _system(ops, ndof, kind=None):
    """Linear system: dense LAPACK for small models, UMFPACK otherwise."""
    if kind is None:
        kind = "FullGeneral" if ndof <= 900 else "UmfPack"
    try:
        ops.system(kind)
    except Exception:                                 # pragma: no cover
        kind = "BandGeneral"
        ops.system(kind)
    return kind


def _disp(ops, n):
    u = []
    for i in range(n):
        u.extend(float(v) for v in ops.nodeDisp(i + 1))
    return u


def _axial(ops, m):
    out = []
    for k in range(m):
        r = ops.eleResponse(k + 1, "axialForce")
        out.append(float(r[0] if isinstance(r, (list, tuple)) else r))
    return out


def _reactions(ops, spec):
    """Reactions at restrained DOFs and spring forces (ground reactions)."""
    n = len(spec["nodes"])
    ops.reactions()
    R = [0.0] * (3 * n)
    spring_f = [0.0] * (3 * n)
    resid = 0.0
    fixed = set(int(d) for d in spec.get("fixed", []))
    _f, _sps, springs, ground = _support_layout(spec)
    spring_dofs = set()
    for node, lst in springs.items():
        g = ops.nodeReaction(ground[node])
        for ax, _k in lst:
            d = 3 * node + ax
            spring_f[d] = float(g[ax])
            spring_dofs.add(d)
    for i in range(n):
        r = ops.nodeReaction(i + 1)
        for a in range(3):
            d = 3 * i + a
            if d in fixed:
                R[d] = float(r[a])
            else:
                resid = max(resid, abs(float(r[a])))
    for d in spring_dofs:
        R[d] += spring_f[d]
    return R, spring_f, resid


# ============================================================= analyses ===
def run_static(spec, opts=None):
    """Linear static analysis to load factor 1."""
    opts = opts or {}
    ops = _ops()
    t0 = time.perf_counter()
    n, m = len(spec["nodes"]), len(spec["members"])
    apply_commands(ops, model_commands(spec, "Truss"))
    apply_commands(ops, load_commands(spec))
    ops.constraints(opts.get("constraints", "Transformation"))
    ops.numberer("RCM")
    sysname = _system(ops, 3 * n, opts.get("system"))
    ops.test("NormUnbalance", 1e-6, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    u = _disp(ops, n)
    N = _axial(ops, m)
    R, spring_f, resid = _reactions(ops, spec)
    out = {"type": "static", "ok": ok == 0, "u": u, "N": N,
           "reactions": R, "spring_forces": spring_f,
           "free_dof_residual": resid,
           "load_factor": float(ops.getLoadFactor(PATTERN)),
           "system": sysname, "elapsed": time.perf_counter() - t0}
    ops.wipe()
    return out


def run_modal(spec, opts=None):
    """Eigenvalues omega^2 and mode shapes."""
    opts = opts or {}
    ops = _ops()
    t0 = time.perf_counter()
    n = len(spec["nodes"])
    mass_type = opts.get("mass_type", "lumped")
    apply_commands(ops, model_commands(spec, "Truss", mass_type=mass_type))
    nfree = 3 * n - len(set(int(d) for d in spec.get("fixed", [])))
    k = opts.get("n_modes", 12)
    solver = opts.get("solver")
    if k in (None, "all") or int(k) >= nfree - 1:
        k = nfree
        solver = solver or "-fullGenLapack"
    k = int(k)
    solver = solver or "-genBandArpack"
    lam = ops.eigen(solver, k)
    lam = [float(v) for v in lam]
    modes = []
    if opts.get("modes", True):
        for j in range(len(lam)):
            vec = []
            for i in range(n):
                vec.extend(float(v) for v in ops.nodeEigenvector(i + 1, j + 1))
            modes.append(vec)
    freq = [math.sqrt(max(v, 0.0)) / (2.0 * math.pi) for v in lam]
    out = {"type": "modal", "ok": len(lam) == k, "omega2": lam,
           "frequencies": freq, "modes": modes, "solver": solver,
           "mass_type": mass_type, "elapsed": time.perf_counter() - t0}
    ops.wipe()
    return out


def _nl_setup(ops, spec, opts):
    n = len(spec["nodes"])
    apply_commands(ops, model_commands(spec, "corotTruss"))
    apply_commands(ops, load_commands(spec, with_sp=False))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    _system(ops, 3 * n, opts.get("system"))
    Fn = math.sqrt(sum(float(v) ** 2 for v in spec["load"]))
    tol = float(opts.get("tol", 1e-9)) * max(Fn, 1e-30)
    ops.test("NormUnbalance", tol, int(opts.get("max_iter", 50)))
    ops.algorithm("Newton")
    return n, len(spec["members"])


def run_path(spec, opts):
    """corotTruss + DisplacementControl to every target control value.

    ``opts``: either ``control_dof`` (0-based global DOF) and ``targets``
    (control displacements [m], in order), or ``segments`` = a list of
    ``{"control_dof", "targets"}`` run one after the other from the state
    the previous segment ended in (piecewise displacement control: a path
    that snaps back in one DOF is followed with another DOF that is
    monotone there).  ``tol`` (unbalance tolerance relative to |F_ref|,
    default 1e-9), ``max_iter``, ``max_split`` (sub-step halvings).
    """
    ops = _ops()
    t0 = time.perf_counter()
    n, m = _nl_setup(ops, spec, opts)
    segments = opts.get("segments") or [
        {"control_dof": opts["control_dof"], "targets": opts["targets"]}]
    max_split = int(opts.get("max_split", 8))
    keep_u = bool(opts.get("keep_u", True))
    lam, ctrl, U, NN, iters, seg_of, dofs = [], [], [], [], [], [], []
    message = ""
    n_targets = sum(len(s["targets"]) for s in segments)
    first = segments[0]["control_dof"]
    ops.integrator("DisplacementControl", int(first) // 3 + 1,
                   int(first) % 3 + 1, 0.0)
    ops.analysis("Static")
    for si, seg in enumerate(segments):
        cdof = int(seg["control_dof"])
        node, dof = cdof // 3 + 1, cdof % 3 + 1
        current = float(ops.nodeDisp(node, dof))
        for w in [float(v) for v in seg["targets"]]:
            # ``nsub`` equal sub-steps remain between the current state and
            # w.  A failed step is reverted by OpenSees (StaticAnalysis
            # reverts the domain to the last commit); the rest is split
            # finer.
            nsub, its = 1, 0
            while nsub > 0:
                step = (w - current) / nsub
                ops.integrator("DisplacementControl", node, dof, step)
                if ops.analyze(1) == 0:
                    current = float(ops.nodeDisp(node, dof))
                    its += int(ops.testIter())
                    nsub -= 1
                else:
                    nsub *= 2
                    if nsub > 2 ** max_split:
                        break
            if nsub > 0:
                message = (f"DisplacementControl failed before the target "
                           f"{w:g} m of DOF {cdof} (reached {current:g} m).")
                break
            lam.append(float(ops.getLoadFactor(PATTERN)))
            ctrl.append(current)
            if keep_u:
                U.append(_disp(ops, n))
                NN.append(_axial(ops, m))
            iters.append(its)
            seg_of.append(si)
            dofs.append(cdof)
        if message:
            break
    out = {"type": "path", "ok": len(lam) == n_targets, "lam": lam,
           "control": ctrl, "control_dofs": dofs, "segment": seg_of,
           "u": U, "N": NN, "iterations": iters, "message": message,
           "elapsed": time.perf_counter() - t0}
    ops.wipe()
    return out


def run_arclength(spec, opts):
    """corotTruss traced by ArcLength (or MinUnbalDispNorm) stepping.

    ``opts``: ``arc_length`` [m] (initial step), ``alpha`` (load scaling in
    the constraint, 0 = cylindrical), ``integrator`` ("ArcLength" |
    "MinUnbalDispNorm"), ``n_steps``, ``control_dof``, ``control_max``
    (stop when |control| exceeds it), ``lambda_min`` (stop below),
    ``min_arc`` (smallest step after halving), ``grow`` (step growth after a
    success, default 1.0 = back to the initial size only).
    """
    ops = _ops()
    t0 = time.perf_counter()
    n, m = _nl_setup(ops, spec, opts)
    s0 = float(opts["arc_length"])
    alpha = float(opts.get("alpha", 0.0))
    integ = opts.get("integrator", "ArcLength")
    n_steps = int(opts.get("n_steps", 400))
    cdof = opts.get("control_dof")
    cmax = float(opts.get("control_max", float("inf")))
    lmin = float(opts.get("lambda_min", -float("inf")))
    lmax = float(opts.get("lambda_max", float("inf")))
    min_arc = float(opts.get("min_arc", s0 / 2 ** 10))
    keep_u = bool(opts.get("keep_u", True))
    lam, ctrl, U, iters = [0.0], [0.0], [[0.0] * (3 * n)] if keep_u else [], [0]
    s = s0
    message = "step limit reached"

    def set_integrator(step):
        if integ == "MinUnbalDispNorm":
            ops.integrator("MinUnbalDispNorm", step)
        else:
            ops.integrator("ArcLength", step, alpha)
        ops.analysis("Static")

    set_integrator(s)
    for _k in range(n_steps):
        if ops.analyze(1) != 0:
            if s / 2.0 < min_arc:
                message = f"no convergence with step {s:g}"
                break
            s /= 2.0
            set_integrator(s)
            continue
        lam.append(float(ops.getLoadFactor(PATTERN)))
        c = float(ops.nodeDisp(int(cdof) // 3 + 1, int(cdof) % 3 + 1)) \
            if cdof is not None else 0.0
        ctrl.append(c)
        if keep_u:
            U.append(_disp(ops, n))
        iters.append(int(ops.testIter()))
        if s < s0:
            s = min(s0, 2.0 * s)
            set_integrator(s)
        if abs(c) > cmax:
            message = "control limit reached"
            break
        if lam[-1] < lmin or lam[-1] > lmax:
            message = "load-factor limit reached"
            break
    out = {"type": "arclength", "ok": len(lam) > 1, "lam": lam,
           "control": ctrl, "u": U, "iterations": iters,
           "message": message, "integrator": integ,
           "elapsed": time.perf_counter() - t0}
    ops.wipe()
    return out


def run_strain_probe(opts):
    """Identify the corotTruss strain measure on one bar.

    Bar from (0,0,0) (pinned) to (L,0,0); its end gets a prescribed
    transverse displacement ``uy = rot * L`` (sp, large rigid rotation) and
    an axial-direction load ``Px = load * EA``; uz is fixed, ux free.  Both
    ramp with the load factor to 1 in ``steps`` Newton steps.  Returns the
    current length ``l`` and the element axial force ``N`` so the caller can
    compare N with EA (l - L)/L (engineering), EA (l^2 - L^2)/(2 L^2)
    (Green-Lagrange, 2nd Piola-Kirchhoff) and EA ln(l/L) (logarithmic).
    """
    ops = _ops()
    L = float(opts.get("L", 2.0))
    E = float(opts.get("E", 200e9))
    A = float(opts.get("A", 1e-3))
    rot = float(opts.get("rot", 0.8))
    load = float(opts.get("load", 0.1))
    steps = int(opts.get("steps", 50))
    element = opts.get("element", "corotTruss")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element(element, 1, 1, 2, A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 2, rot * L)
    ops.load(2, load * E * A, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-9 * E * A, 100)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    ok = ops.analyze(steps)
    x = L + float(ops.nodeDisp(2, 1))
    y = float(ops.nodeDisp(2, 2))
    r = ops.eleResponse(1, "axialForce")
    N = float(r[0] if isinstance(r, (list, tuple)) else r)
    out = {"type": "strain_probe", "ok": ok == 0, "L": L, "EA": E * A,
           "x": x, "y": y, "l": math.hypot(x, y), "N": N,
           "Px": load * E * A, "element": element,
           "lambda": float(ops.getLoadFactor(1))}
    ops.wipe()
    return out


def run_thermal_probe(opts):
    """Heated fixed-fixed bar: InitStrainMaterial(-eps0) must give
    N = -EA eps0 (compression for heating) and zero displacement.

    Two collinear bars 1-2-3, both ends fixed, node 2 free in x only; the
    free strain ``eps0`` is on bar 1 only.  Closed form: u2 = eps0 L A1 /
    (A1 + A2) ... with A2 = r A1:  u2 = eps0 L / (1 + r),
    N1 = N2 = -EA1 eps0 r / (1 + r).
    """
    ops = _ops()
    L = float(opts.get("L", 2.0))
    E = float(opts.get("E", 200e9))
    A = float(opts.get("A", 1e-3))
    r = float(opts.get("ratio", 2.0))
    eps0 = float(opts.get("eps0", 11.7e-6 * 30.0))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, L, 0.0, 0.0)
    ops.node(3, 2 * L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(3, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.uniaxialMaterial("InitStrainMaterial", 2, 1, -eps0)
    ops.element("Truss", 1, 1, 2, A, 2)
    ops.element("Truss", 2, 2, 3, r * A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-6, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    ops.reactions()
    out = {"type": "thermal_probe", "ok": ok == 0, "L": L, "EA": E * A,
           "ratio": r, "eps0": eps0, "u2": float(ops.nodeDisp(2, 1)),
           "N1": _axial(ops, 2)[0], "N2": _axial(ops, 2)[1],
           "R1": float(ops.nodeReaction(1, 1)),
           "R3": float(ops.nodeReaction(3, 1))}
    ops.wipe()
    return out


RUNNERS = {"static": run_static, "modal": run_modal, "path": run_path,
           "arclength": run_arclength}
PROBES = {"strain_probe": run_strain_probe,
          "thermal_probe": run_thermal_probe}


def _normalise(analysis):
    if isinstance(analysis, str):
        return {"type": analysis}
    return dict(analysis)


def run(spec, analysis="static"):
    """Run one analysis; returns a JSON-able dict (``ok``, ``error``...)."""
    a = _normalise(analysis)
    kind = a.get("type", "static")
    try:
        if kind in PROBES:
            out = PROBES[kind](a)
        elif kind in RUNNERS:
            out = RUNNERS[kind](spec, a)
        else:
            raise ValueError(f"Unknown analysis type {kind!r}.")
        out["error"] = ""
    except Exception as exc:
        try:
            _ops().wipe()
        except Exception:                              # pragma: no cover
            pass
        out = {"type": kind, "ok": False,
               "error": f"{type(exc).__name__}: {exc}"}
    out["opensees_version"] = opensees_version()
    return out


def run_jobs(jobs):
    return [run(j.get("spec"), j.get("analysis", "static")) for j in jobs]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("usage: python _opensees_runner.py job.json result.json")
        return 2
    with open(argv[0], "r", encoding="utf-8") as fh:
        job = json.load(fh)
    jobs = job["jobs"] if "jobs" in job else [job]
    results = run_jobs(jobs)
    with open(argv[1], "w", encoding="utf-8") as fh:
        json.dump({"results": results,
                   "python": sys.version.split()[0],
                   "opensees_version": opensees_version()}, fh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
