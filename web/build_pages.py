"""
Build the static web site that runs web/streamlit_app.py in the browser.

The site uses stlite (Streamlit compiled for Pyodide, the WebAssembly build of
Python), so it needs no server and no account: GitHub Pages serves plain
files and the visitor's browser runs the analysis.  The site holds

    index.html              mounts the app with a pinned @stlite/browser
    web/streamlit_app.py    the app, unchanged
    geodome/*.py            only the core modules the app imports (found by
                            reading the import statements, so the list
                            follows the code)
    validation/hand_calc.py the independent hand calculation (Model A)
    LICENSE                 the GeoDomeLab licence (attribution)
    .nojekyll               stops GitHub Pages from hiding __init__.py
    demo/                   with --video: a page that plays the 3.5-minute
                            demo video (demo/index.html, the mp4, a poster)

Python packages: numpy and scipy come from the Pyodide distribution.  The app
draws its 3D view with plotly.js from a pinned CDN file when Python plotly is
missing, so plotly (a 9.7 MB wheel) is not installed in the browser.

Usage (from the repository root):

    python web/build_pages.py --out <folder>            # build
    python web/build_pages.py --out <folder> --check    # build + self-test
    python web/build_pages.py --out <folder> --video GeoDomeLab_demo.mp4 \
        --poster poster.jpg                              # + demo video page

The video is not kept in the repository (it is 31 MB).  It is attached to the
GitHub release and copied into the site when --video is given.

--check imports the copied modules from the site folder only, in an isolated
interpreter, and runs both examples and the hand calculation, so a module
missing from the site is caught here and not in the visitor's browser.

Preview locally:  python -m http.server 8600 --directory <folder>
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import html
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

STLITE_VERSION = "1.12.1"         # pinned: @stlite/browser on npm, 2026-11-05
STLITE = f"https://cdn.jsdelivr.net/npm/@stlite/browser@{STLITE_VERSION}/build"
REQUIREMENTS = ["numpy", "scipy"]  # both are prebuilt Pyodide packages
APP = "web/streamlit_app.py"
EXTRA_FILES = ["validation/hand_calc.py"]
STATIC_FILES = ["LICENSE"]         # copied, not mounted into Python
REPO_URL = "https://github.com/covdetbaba/GeoDomeLab"
RELEASE_URL = REPO_URL + "/releases/tag/v2.0.0"
DEMO_DIR = "demo"
VIDEO_NAME = "GeoDomeLab_demo.mp4"
POSTER_NAME = "poster.jpg"
# start of each part of the demo video, in seconds
CHAPTERS = [
    (0, "Title"),
    (5, "Browser version: the 6V example, no installation"),
    (47, "Desktop program: the input"),
    (77, "Run, deflected shape, axial force, pick a strut"),
    (112, "Auto-size the struts"),
    (133, "Nonlinear path to the limit point"),
    (172, "Imperfection study (shown faster)"),
    (199, "End: try it yourself"),
]


# =============================================================================
# which geodome modules does the app need?
# =============================================================================
def _geodome_imports_of_app(path):
    """Top-level geodome modules imported by the app file."""
    tree = ast.parse(open(path, encoding="utf-8").read(), path)
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "geodome":
                mods.update(a.name for a in node.names)
            elif node.module.startswith("geodome."):
                mods.add(node.module.split(".")[1])
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith("geodome."):
                    mods.add(a.name.split(".")[1])
    return mods


def _relative_imports(path):
    """Sibling modules named by ``from . import x`` / ``from .x import y``
    anywhere in the file (function-level imports included)."""
    tree = ast.parse(open(path, encoding="utf-8").read(), path)
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 1:
            if node.module:
                out.add(node.module.split(".")[0])
            else:
                out.update(a.name for a in node.names)
    return out


def geodome_modules(app_path=os.path.join(ROOT, APP)):
    """Sorted list of geodome/<name>.py files the app needs (closure)."""
    pkg = os.path.join(ROOT, "geodome")
    todo = [m for m in _geodome_imports_of_app(app_path)
            if os.path.isfile(os.path.join(pkg, m + ".py"))]
    seen = {"__init__"}
    while todo:
        m = todo.pop()
        if m in seen:
            continue
        if m == "gui":
            raise RuntimeError("The web app must not import geodome.gui.")
        seen.add(m)
        for dep in _relative_imports(os.path.join(pkg, m + ".py")):
            if os.path.isfile(os.path.join(pkg, dep + ".py")):
                todo.append(dep)
    return sorted(f"geodome/{m}.py" for m in seen)


# =============================================================================
# index.html
# =============================================================================
def _short_hash(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()[:10]


INDEX = """\
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
<title>GeoDomeLab in the browser</title>
<meta name="description" content="GeoDomeLab {version}: geodesic dome space-truss analysis by the direct stiffness method, running in your browser.">
<link rel="stylesheet" href="{stlite}/stlite.css">
<style>
  html, body {{ margin: 0; background: #ffffff; }}
  #gdl-note {{
    position: fixed; left: 50%; bottom: 16px; transform: translateX(-50%);
    z-index: 1000; box-sizing: border-box; width: max-content;
    max-width: calc(100% - 32px); padding: 12px 16px; border-radius: 8px;
    background: #1f2a36; color: #ffffff;
    font: 15px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif;
    box-shadow: 0 2px 12px rgba(0, 0, 0, 0.25);
  }}
  #gdl-note a {{ color: #9fd0ff; }}
  #gdl-time {{ font-variant-numeric: tabular-nums; }}
</style>
</head>
<body>
<div id="gdl-note" role="status">
  <strong>GeoDomeLab {version} is loading<span id="gdl-time"></span>.</strong>
  The first load takes about 20 to 40 s while Python, NumPy and SciPy
  download into your browser. The analysis then runs on your device.
  <a href="{repo}">Source on GitHub</a> ·
  <a href="demo/">Demo video (3.5 min)</a>
</div>
<div id="root"></div>
<script type="module">
import {{ mount }} from "{stlite}/stlite.js";

const note = document.getElementById("gdl-note");
const clock = document.getElementById("gdl-time");
mount({options}, document.getElementById("root"));

// Count the seconds while loading, and hide the note once the app has drawn
// its title.  The time is kept in window.__gdlReadyMs (milliseconds from
// navigation start).
const timer = setInterval(() => {{
  if (document.querySelector("#root h1")) {{
    clearInterval(timer);
    window.__gdlReadyMs = Math.round(performance.now());
    console.log("GeoDomeLab ready after " + window.__gdlReadyMs + " ms");
    note.style.display = "none";
  }} else {{
    clock.textContent = " (" + Math.round(performance.now() / 1000) + " s)";
  }}
}}, 250);
setTimeout(() => {{
  if (!window.__gdlReadyMs) {{
    note.innerHTML = "Still loading. A slow connection can take a few " +
      "minutes. If nothing appears, reload the page or use the desktop " +
      'program from <a href="{repo}">GitHub</a>.';
  }}
}}, 120000);
</script>
</body>
</html>
"""


DEMO = """\
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GeoDomeLab demo video</title>
<meta name="description" content="A 3.5-minute recording of GeoDomeLab {version}: the browser version and the desktop program, from input to output.">
<style>
  :root {{ --navy: #13243a; --gold: #f0a24b; --ink: #1c2733; --muted: #5b6b7c;
           --line: #d9e0e8; --bg: #f4f6f9; }}
  html, body {{ margin: 0; background: var(--bg); color: var(--ink);
    font: 16px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, Arial, sans-serif; }}
  header {{ background: var(--navy); color: #fff; padding: 20px 16px 16px;
    border-bottom: 4px solid var(--gold); }}
  header h1 {{ margin: 0; font-size: 1.6rem; }}
  header p {{ margin: 4px 0 0; color: #dbe4f1; }}
  main {{ max-width: 1100px; margin: 0 auto; padding: 16px; }}
  video {{ width: 100%; height: auto; background: #000; border-radius: 6px;
    display: block; }}
  .links {{ margin: 12px 0 4px; display: flex; flex-wrap: wrap; gap: 8px; }}
  .links a {{ padding: 8px 14px; border-radius: 6px; text-decoration: none;
    background: #fff; border: 1px solid var(--line); color: var(--navy); }}
  .links a.main {{ background: #e5484d; border-color: #e5484d; color: #fff; }}
  h2 {{ font-size: 1.1rem; margin: 20px 0 6px; }}
  ol {{ margin: 0; padding-left: 0; list-style: none; }}
  li button {{ width: 100%; text-align: left; border: 0; background: none;
    padding: 6px 4px; font: inherit; color: var(--ink); cursor: pointer;
    border-bottom: 1px solid var(--line); }}
  li button:hover {{ background: #fff; }}
  li .t {{ display: inline-block; min-width: 3.2em; color: var(--muted);
    font-variant-numeric: tabular-nums; }}
  footer {{ color: var(--muted); font-size: 0.9rem; margin: 20px 0; }}
</style>
</head>
<body>
<header>
  <h1>GeoDomeLab {version}: demo video</h1>
  <p>Geodesic dome space-truss analysis by the direct stiffness method.
  The recording runs one complete example from input to output.</p>
</header>
<main>
  <video id="v" controls playsinline preload="metadata" poster="{poster}">
    <source src="{video}" type="video/mp4">
    Your browser cannot play this video. <a href="{video}">Download it</a>.
  </video>
  <div class="links">
    <a class="main" href="../">Try it in your browser</a>
    <a href="{repo}">Source code on GitHub</a>
    <a href="{release}">Release v{version}</a>
    <a href="{video}" download>Download the video (mp4)</a>
  </div>
  <h2>Parts of the video</h2>
  <ol>
{chapters}
  </ol>
  <footer>3 min 25 s, 1920 &times; 1080, no sound (captions in the picture).
  Recorded on 22 November 2026 from the real program: the browser version at
  covdetbaba.github.io/GeoDomeLab and the desktop window of release
  v{version}.</footer>
</main>
<script>
  const v = document.getElementById("v");
  for (const b of document.querySelectorAll("button[data-t]")) {{
    b.addEventListener("click", () => {{
      v.currentTime = Number(b.dataset.t);
      v.play();
      v.scrollIntoView({{behavior: "smooth", block: "nearest"}});
    }});
  }}
</script>
</body>
</html>
"""


def demo_html(version):
    rows = []
    for t, name in CHAPTERS:
        stamp = f"{t // 60}:{t % 60:02d}"
        rows.append(f'    <li><button data-t="{t}"><span class="t">{stamp}'
                    f'</span>{html.escape(name)}</button></li>')
    return DEMO.format(version=html.escape(version), poster=POSTER_NAME,
                       video=VIDEO_NAME, repo=REPO_URL, release=RELEASE_URL,
                       chapters="\n".join(rows))


def index_html(files, version):
    options = {
        "entrypoint": APP,
        "requirements": REQUIREMENTS,
        "files": {f: {"url": f"{f}?v={h}"} for f, h in files.items()},
        "streamlitConfig": {"client.toolbarMode": "viewer"},
        # The note above shows the progress.  stlite's own progress toasts
        # stay on screen when the page loses focus and cover a phone screen.
        "disableProgressToasts": True,
    }
    opts = json.dumps(options, indent=2)
    return INDEX.format(stlite=STLITE, options=opts, repo=REPO_URL,
                        version=html.escape(version))


# =============================================================================
# build and check
# =============================================================================
def _version():
    tree = ast.parse(open(os.path.join(ROOT, "geodome", "__init__.py"),
                          encoding="utf-8").read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", "") == "__version__" for t in node.targets):
            return node.value.value
    return "?"


def build(out, video=None, poster=None):
    mods = geodome_modules()
    mounted = [APP] + mods + EXTRA_FILES
    if os.path.isdir(out):
        for name in os.listdir(out):          # keep a .git folder if present
            if name == ".git":
                continue
            p = os.path.join(out, name)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    os.makedirs(out, exist_ok=True)

    hashes = {}
    for rel in mounted + STATIC_FILES:
        src = os.path.join(ROOT, rel)
        dst = os.path.join(out, rel)
        os.makedirs(os.path.dirname(dst) or out, exist_ok=True)
        with open(src, "rb") as fh:                 # text files, LF endings
            data = fh.read().replace(b"\r\n", b"\n")
        with open(dst, "wb") as fh:
            fh.write(data)
        if rel in mounted:
            hashes[rel] = _short_hash(dst)

    with open(os.path.join(out, "index.html"), "w", encoding="utf-8",
              newline="\n") as fh:
        fh.write(index_html(hashes, _version()))
    open(os.path.join(out, ".nojekyll"), "w").close()

    if video:                               # the demo video page
        demo = os.path.join(out, DEMO_DIR)
        os.makedirs(demo, exist_ok=True)
        shutil.copyfile(video, os.path.join(demo, VIDEO_NAME))
        if poster:
            shutil.copyfile(poster, os.path.join(demo, POSTER_NAME))
        with open(os.path.join(demo, "index.html"), "w", encoding="utf-8",
                  newline="\n") as fh:
            fh.write(demo_html(_version()))
        for name in sorted(os.listdir(demo)):
            print(f"  {DEMO_DIR}/{name:27s} "
                  f"{os.path.getsize(os.path.join(demo, name)):11,d} bytes")

    total = 0
    print(f"Site written to {out}")
    for rel in ["index.html"] + mounted + STATIC_FILES:
        size = os.path.getsize(os.path.join(out, rel))
        total += size
        print(f"  {rel:32s} {size:9,d} bytes")
    print(f"  {len(mounted)} Python files mounted, {total:,d} bytes in all, "
          f"stlite {STLITE_VERSION}, requirements {REQUIREMENTS}")
    return mounted


CHECK = r"""
import os, sys, time
site = sys.argv[1]
sys.path.insert(0, site)
sys.path.insert(0, os.path.join(site, "validation"))
import numpy as np
from geodome import benchmarks as bm, postprocess as pp, optimize as opt
import copy, hand_calc
for name, mod in list(sys.modules.items()):
    f = getattr(mod, "__file__", None) or ""
    if name == "geodome" or name.startswith("geodome."):
        assert os.path.abspath(f).startswith(os.path.abspath(site)), (name, f)
t = time.perf_counter()
m = bm.make_preset("case_study_6v")
r = pp.run_analysis(m)
N = r.cases["LC2  1.2D + 1.6S"].N / 1e3
print("6V trial: %d nodes, %d struts, %.0f kg, util %.3f (%s), "
      "lambda_cr %.2f, LC2 N %+.2f (m%d) / %+.2f (m%d) kN"
      % (m.geom.n_nodes, m.geom.n_members, m.steel_mass(), r.governing_util,
         r.governing_case.split()[0], min(r.lambda_cr_by_case.values()),
         N.max(), N.argmax(), N.min(), N.argmin()))
m2 = copy.deepcopy(m); s = opt.auto_size(m2); r2 = pp.run_analysis(m2)
print("6V sized: %.0f kg, util %.3f (%s), lambda_cr %.2f"
      % (m2.steel_mass(), r2.governing_util, r2.governing_case.split()[0],
         min(r2.lambda_cr_by_case.values())))
H = hand_calc.model_a()
ma = bm.validation_model_a(); ra = pp.run_analysis(ma)
c = ma.geom.crown_node
uz = ra.cases[bm.COMBO_A1].u[3 * c + 2]
err = abs(uz - H["A1"]["D"][2]) / abs(H["A1"]["D"][2]) * 100
print("Model A: crown uz %.6f mm, hand %.6f mm, error %.1e %%"
      % (uz * 1e3, H["A1"]["D"][2] * 1e3, err))
assert err < 1e-10
print("geodome modules loaded:", sorted(n for n in sys.modules
                                         if n.startswith("geodome.")))
print("check time %.2f s" % (time.perf_counter() - t))
"""


def check(out):
    """Import and run the core from the site folder only (isolated mode)."""
    # -I: ignore PYTHONPATH and the user site, -B: no __pycache__ in the site
    cmd = [sys.executable, "-I", "-B", "-c", CHECK, os.path.abspath(out)]
    proc = subprocess.run(cmd, cwd=out, capture_output=True, text=True)
    for folder, dirs, _ in os.walk(out):
        if "__pycache__" in dirs:
            shutil.rmtree(os.path.join(folder, "__pycache__"))
    print(proc.stdout)
    if proc.returncode != 0:
        print(proc.stderr)
        raise SystemExit("Site check FAILED")
    print("Site check passed.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--check", action="store_true",
                    help="run the core from the built site as a self-test")
    ap.add_argument("--video", help="demo video (mp4) for the demo/ page")
    ap.add_argument("--poster", help="poster image (jpg) for the video")
    a = ap.parse_args(argv)
    build(a.out, video=a.video, poster=a.poster)
    if a.check:
        check(a.out)


if __name__ == "__main__":
    main()
