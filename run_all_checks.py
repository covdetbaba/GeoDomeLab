#!/usr/bin/env python
"""
Run every check of GeoDomeLab with one command.

    python run_all_checks.py

It runs, each in its own Python process:

1. every ``tests/test_*.py`` file on its own (no pytest needed);
2. ``run_geodome.py --verify``, the built-in verification suite;
3. ``validation/run_validation.py``, Models A and B against the hand
   calculation in ``validation/hand_calc.py``;
4. ``validation/hand_calc_ext.py``, the extended hand calculations F1-F10
   (this file must not import the program);
5. ``validation/run_validation_ext.py``, the program against F1-F10;
6. ``validation/opensees_crosscheck.py --quick``, the cross-check against
   OpenSees, only when OpenSees is found (``geodome.opensees_bridge``).
   The quick run writes its report to a temporary folder, so the full report
   in ``validation/output`` is kept.

Steps 3 to 5 rewrite their own reports in ``validation/output``, as they do
when you run them by hand.

The children run with Qt offscreen (``QT_QPA_PLATFORM=offscreen``) and the
Agg matplotlib backend, so no window opens.  The script works from any
folder.  It prints one line per step, then a table with the counts, then
PASS or FAIL.  The exit code is 0 when everything passes (skips are allowed)
and 1 otherwise.

Options
-------
    --jobs N           run the test files N at a time (default 1)
    --only TEXT        run only the steps whose name contains TEXT
    --skip-opensees    do not run the OpenSees cross-check
    --opensees-full    run the full OpenSees matrix (rewrites its report in
                       validation/output; takes a few minutes)
    --timeout S        time limit per step in seconds (default 1800)
    --verbose          print the full output of every step
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VALIDATION = ROOT / "validation"
TESTS = ROOT / "tests"

# The pass rule of run_validation.py (REL_TOL, ABS_TOL there) and of the
# GUI's Validation menu: |p - h| / |h| <= 1e-12 (that is 1e-10 %), or
# |p| <= 1e-9 in the row unit when the hand answer is exactly zero.  A test
# checks that the copies stay equal.
RV_REL_TOL_PERCENT = 1e-10
RV_ABS_TOL = 1e-9


# ------------------------------------------------------------ environment --
def child_env():
    """Environment for every child: headless Qt and matplotlib, UTF-8."""
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["MPLBACKEND"] = "Agg"
    env["PYTHONIOENCODING"] = "utf-8"
    if os.name == "nt" and not env.get("QT_QPA_FONTDIR"):
        fonts = os.path.join(env.get("WINDIR", r"C:\Windows"), "Fonts")
        if os.path.isdir(fonts):
            env["QT_QPA_FONTDIR"] = fonts.replace("\\", "/")
    root = str(ROOT)
    old = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = root + (os.pathsep + old if old else "")
    return env


def run_child(args, env, timeout):
    """Run ``python *args`` in the repository root.

    Returns ``(returncode or None, output text, seconds)``; ``None`` means
    the step timed out or could not start.
    """
    t0 = time.perf_counter()
    try:
        cp = subprocess.run([sys.executable, *args], cwd=str(ROOT), env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace",
                            timeout=timeout)
        return cp.returncode, cp.stdout or "", time.perf_counter() - t0
    except subprocess.TimeoutExpired as exc:
        out = exc.stdout or ""
        if isinstance(out, bytes):
            out = out.decode("utf-8", "replace")
        return None, out + f"\n[timed out after {timeout} s]", \
            time.perf_counter() - t0
    except OSError as exc:
        return None, f"[could not start: {exc}]", time.perf_counter() - t0


# ----------------------------------------------------------------- result --
class Result:
    """Outcome of one step."""

    def __init__(self, name, group, status, passed=None, total=None,
                 skipped=0, seconds=0.0, note="", output=""):
        self.name = name
        self.group = group
        self.status = status            # PASS, FAIL, SKIP or ERROR
        self.passed = passed
        self.total = total
        self.skipped = skipped
        self.seconds = seconds
        self.note = note
        self.output = output

    @property
    def ok(self):
        return self.status in ("PASS", "SKIP")


def _last(pattern, text):
    m = None
    for m in re.finditer(pattern, text):
        pass
    return m


def _status(code, passed, total, extra_ok=True):
    if code is None:
        return "ERROR"
    if passed is None or total is None:
        return "ERROR"
    return "PASS" if (code == 0 and passed == total and extra_ok) else "FAIL"


# ------------------------------------------------------------------ steps --
def step_test(path, env, timeout):
    rel = path.relative_to(ROOT).as_posix()
    code, out, sec = run_child([rel], env, timeout)
    m = _last(r"(\d+) of (\d+) tests passed(?:, (\d+) skipped)?", out)
    if m:
        passed, total = int(m.group(1)), int(m.group(2))
        skipped = int(m.group(3) or 0)
        # a skipped test counts as neither passed nor failed
        ok_count = passed + skipped == total
        status = "ERROR" if code is None else (
            "PASS" if code == 0 and ok_count else "FAIL")
        return Result(rel, "tests", status, passed, total - skipped, skipped,
                      sec, "", out)
    return Result(rel, "tests", "ERROR" if code else "FAIL", None, None, 0,
                  sec, "no test summary line in the output", out)


def step_verify(env, timeout):
    code, out, sec = run_child(["run_geodome.py", "--verify"], env, timeout)
    m = _last(r"(\d+) of (\d+) checks passed", out)
    p, t = (int(m.group(1)), int(m.group(2))) if m else (None, None)
    return Result("run_geodome.py --verify", "verification",
                  _status(code, p, t), p, t, 0, sec,
                  "" if m else "no summary line", out)


def step_validation_ab(env, timeout):
    code, out, sec = run_child([str(Path("validation", "run_validation.py"))],
                               env, timeout)
    m = _last(r"(\d+) quantities compared\.\s+Largest relative difference "
              r"([0-9.eE+-]+) %\.\s+Largest absolute value where the hand "
              r"answer is zero ([0-9.eE+-]+)", out)
    note = ""
    passed = total = None
    path = VALIDATION / "output" / "comparison.csv"
    if code == 0 and path.is_file():
        try:
            with open(path, newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
            total = len(rows)
            passed = 0
            for r in rows:
                d = float(r["difference"])
                if r.get("diff_kind") == "%":
                    passed += d <= RV_REL_TOL_PERCENT
                else:
                    passed += d <= RV_ABS_TOL
        except Exception as exc:                    # pragma: no cover
            note = f"comparison.csv not read: {exc}"
    if m:
        note = (f"worst {float(m.group(2)):.1e} %"
                + (f"; {note}" if note else ""))
        if total is not None and int(m.group(1)) != total:
            note += f"; printed {m.group(1)} rows, csv {total}"
    return Result("validation/run_validation.py", "hand",
                  _status(code, passed, total, m is not None), passed, total,
                  0, sec, note or "no summary line", out)


def step_hand_ext(env, timeout):
    code, out, sec = run_child([str(Path("validation", "hand_calc_ext.py"))],
                               env, timeout)
    m = _last(r"(\d+) reference values, (\d+) of (\d+) internal checks", out)
    imp = _last(r"geodome imported: (True|False)", out)
    independent = imp is not None and imp.group(1) == "False"
    p, t = (int(m.group(2)), int(m.group(3))) if m else (None, None)
    note = (f"{m.group(1)} reference values" if m else "no summary line")
    if not independent:
        note += "; the hand file must not import geodome"
    return Result("validation/hand_calc_ext.py", "hand",
                  _status(code, p, t, independent), p, t, 0, sec, note, out)


def step_validation_ext(env, timeout):
    code, out, sec = run_child(
        [str(Path("validation", "run_validation_ext.py"))], env, timeout)
    m = _last(r"(\d+) of (\d+) comparisons pass; (\d+) of (\d+) hand "
              r"reference values compared", out)
    p, t = (int(m.group(1)), int(m.group(2))) if m else (None, None)
    note = (f"{m.group(3)} of {m.group(4)} hand values used" if m
            else "no summary line")
    all_used = m is not None and m.group(3) == m.group(4)
    return Result("validation/run_validation_ext.py", "hand",
                  _status(code, p, t, all_used), p, t, 0, sec, note, out)


def opensees_mode(env):
    """'in-process', 'subprocess' or None, as geodome.opensees_bridge sees
    it in a child process (the same interpreter as the checks)."""
    code, out, _ = run_child(
        ["-c", "from geodome import opensees_bridge as ob; "
               "print('MODE=' + str(ob.opensees_status()['mode']))"],
        env, 300)
    m = _last(r"MODE=(\S+)", out)
    if code != 0 or not m or m.group(1) == "None":
        return None
    return m.group(1)


def step_opensees(env, timeout, mode, full=False):
    name = "validation/opensees_crosscheck.py" + ("" if full else " --quick")
    if mode is None:
        return Result(name, "opensees", "SKIP", note="OpenSees not found "
                      "(pip install openseespy, or set "
                      "GEODOME_OPENSEES_PYTHON)")
    tmp = None
    args = [str(Path("validation", "opensees_crosscheck.py"))]
    if not full:
        tmp = tempfile.mkdtemp(prefix="geodome_os_")
        args += ["--quick", "--out", tmp]
    try:
        code, out, sec = run_child(args, env, timeout)
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    if code == 2:
        return Result(name, "opensees", "SKIP", seconds=sec,
                      note="OpenSees unavailable: " + out.strip()[-120:],
                      output=out)
    m = _last(r"(\d+) of (\d+) comparisons pass", out)
    p, t = (int(m.group(1)), int(m.group(2))) if m else (None, None)
    return Result(name, "opensees", _status(code, p, t), p, t, 0, sec,
                  f"OpenSees {mode}" if m else "no summary line", out)


# ----------------------------------------------------------------- report --
GROUP_TITLES = {
    "tests": "Unit and GUI tests",
    "verification": "Built-in verification suite",
    "hand": "Hand-calculation validation",
    "opensees": "OpenSees cross-check",
}


def _n(v):
    return "-" if v is None else str(v)


def print_table(results):
    w = max([len(r.name) for r in results] + [30])
    head = (f"  {'Check':<{w}}  {'Result':<6}  {'Passed':>6}  {'Total':>6}"
            f"  {'Skipped':>7}  {'Time [s]':>8}  Note")
    print(head)
    print("  " + "-" * (len(head) + 18))
    for r in results:
        print(f"  {r.name:<{w}}  {r.status:<6}  {_n(r.passed):>6}  "
              f"{_n(r.total):>6}  {r.skipped or '':>7}  {r.seconds:>8.1f}"
              f"  {r.note}")
    print()
    print("  Summary by kind of check")
    for g, title in GROUP_TITLES.items():
        rs = [r for r in results if r.group == g]
        if not rs:
            continue
        p = sum(r.passed or 0 for r in rs)
        t = sum(r.total or 0 for r in rs)
        s = sum(r.skipped or 0 for r in rs)
        bad = [r for r in rs if not r.ok]
        skipped_steps = [r for r in rs if r.status == "SKIP"]
        state = "FAIL" if bad else ("SKIP" if len(skipped_steps) == len(rs)
                                    else "PASS")
        extra = f", {s} skipped" if s else ""
        files = (f" in {len(rs)} file{'s' if len(rs) != 1 else ''}"
                 if g == "tests" else "")
        print(f"    {title:<30} {state:<5} {p} of {t} pass{extra}{files}")


def main(argv=None):
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:                              # pragma: no cover
            pass
    ap = argparse.ArgumentParser(
        description="Run every GeoDomeLab test, verification and validation "
                    "check, and print one table.")
    ap.add_argument("--jobs", type=int, default=1,
                    help="run the test files N at a time (default 1)")
    ap.add_argument("--only", default="",
                    help="run only the steps whose name contains this text")
    ap.add_argument("--skip-opensees", action="store_true",
                    help="do not run the OpenSees cross-check")
    ap.add_argument("--opensees-full", action="store_true",
                    help="run the full OpenSees matrix (rewrites its report)")
    ap.add_argument("--timeout", type=float, default=1800.0,
                    help="time limit per step in seconds (default 1800)")
    ap.add_argument("--verbose", action="store_true",
                    help="print the full output of every step")
    a = ap.parse_args(argv)

    env = child_env()
    t_start = time.perf_counter()
    sel = a.only.lower()

    def wanted(name):
        return not sel or sel in name.lower()

    try:
        import geodome                                 # only for the title
        version = geodome.__version__
    except Exception:                                  # pragma: no cover
        sys.path.insert(0, str(ROOT))
        try:
            import geodome
            version = geodome.__version__
        except Exception:
            version = "?"
    print(f"GeoDomeLab {version}: all checks")
    print(f"Python {sys.version.split()[0]} on {sys.platform}; "
          f"Qt offscreen, matplotlib Agg; repository {ROOT.name}/")

    tests = sorted(TESTS.glob("test_*.py"))
    jobs = []                       # (name, callable)
    for p in tests:
        name = p.relative_to(ROOT).as_posix()
        if wanted(name):
            jobs.append((name, lambda p=p: step_test(p, env, a.timeout)))
    later = [("run_geodome.py --verify",
              lambda: step_verify(env, a.timeout)),
             ("validation/run_validation.py",
              lambda: step_validation_ab(env, a.timeout)),
             ("validation/hand_calc_ext.py",
              lambda: step_hand_ext(env, a.timeout)),
             ("validation/run_validation_ext.py",
              lambda: step_validation_ext(env, a.timeout))]
    later = [(n, f) for n, f in later if wanted(n)]
    os_name = "validation/opensees_crosscheck.py"
    run_os = wanted(os_name) and not a.skip_opensees
    n_steps = len(jobs) + len(later) + (1 if (run_os or a.skip_opensees)
                                        and wanted(os_name) else 0)
    if n_steps == 0:
        print(f"No step matches --only {a.only!r}.")
        return 1
    print(f"{n_steps} step{'s' if n_steps != 1 else ''}.\n")

    results = []
    counter = [0]

    def report(r):
        counter[0] += 1
        if r.total is not None:
            counts = f"{r.passed} of {r.total}"
            if r.skipped:
                counts += f", {r.skipped} skipped"
        else:
            counts = r.note
        print(f"[{counter[0]:>2}/{n_steps}] {r.name:<44} {r.status:<5} "
              f"({counts}, {r.seconds:.1f} s)", flush=True)
        if a.verbose and r.output:
            print(r.output.rstrip())
            print()

    # test files (optionally several at a time), then the validation steps
    if a.jobs > 1 and len(jobs) > 1:
        with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
            futs = {ex.submit(fn): name for name, fn in jobs}
            done = {}
            for f in cf.as_completed(futs):
                r = f.result()
                done[futs[f]] = r
                report(r)
        results += [done[name] for name, _ in jobs]
    else:
        for name, fn in jobs:
            r = fn()
            report(r)
            results.append(r)
    for name, fn in later:
        r = fn()
        report(r)
        results.append(r)
    if wanted(os_name):
        if a.skip_opensees:
            r = Result(os_name, "opensees", "SKIP",
                       note="skipped (--skip-opensees)")
        else:
            mode = opensees_mode(env)
            r = step_opensees(env, a.timeout, mode, full=a.opensees_full)
        report(r)
        results.append(r)

    total_s = time.perf_counter() - t_start
    print()
    print_table(results)
    bad = [r for r in results if not r.ok]
    for r in bad:
        print()
        print(f"---- output of {r.name} (last 40 lines) " + "-" * 20)
        lines = r.output.rstrip().splitlines()
        print("\n".join(lines[-40:]) if lines else "(no output)")
    print()
    n_pass = sum(r.passed or 0 for r in results)
    n_total = sum(r.total or 0 for r in results)
    n_skip = sum(r.skipped or 0 for r in results)
    skipped_steps = sum(r.status == "SKIP" for r in results)
    verdict = "PASS" if not bad else "FAIL"
    print(f"OVERALL: {verdict}.  {n_pass} of {n_total} checks pass in "
          f"{len(results)} step{'s' if len(results) != 1 else ''}"
          + (f"; {n_skip} test{'s' if n_skip != 1 else ''} skipped"
             if n_skip else "")
          + (f"; {skipped_steps} step(s) skipped" if skipped_steps else "")
          + (f"; {len(bad)} step(s) failed" if bad else "")
          + f".  {total_s:.0f} s.")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
