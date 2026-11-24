# GeoDomeLab in the browser

A small Streamlit front end for GeoDomeLab. It imports the same `geodome`
package as the desktop program, so the numbers are the desktop numbers. There
is no copy of the analysis code here and no Qt.

**Open it with no installation:** <https://covdetbaba.github.io/GeoDomeLab/>

The first visit downloads about 45 MB (Python, NumPy, SciPy, pandas and
Streamlit for Pyodide). That took 10 s on a fast university connection, so
allow up to a minute on a phone or on conference Wi-Fi. A second visit loads from
the browser cache in about 7 s. A run of the 6V example then takes about one
second on a laptop, and nothing is sent to a server.

## What the app does

* **Example: 6V dome, 30 m span.** The case study of the report: R = 15 m
  hemisphere, 196 nodes, 555 struts, trial section CHS 114.3 x 5.0, ASCE 7-16
  LC1 to LC4 and SLS. One more click auto-sizes the struts (22,445 kg to
  7,504 kg, utilisation 0.808 under LC4, critical load factor 20.67 under LC2).
* **Verification: Model A against the hand calculation.** The 1V cap of
  `validation/hand_calc.py` (which does not import geodome) against the
  program: crown displacements, member forces, reactions, utilisation and the
  buckling load factor, with the error in per cent.
* **Custom dome.** Frequency 1V to 6V, radius, dome fraction, supports, one
  section for every strut, dead, snow and wind loads, optional auto-sizing.

Outputs: crown displacement (mm), the largest tension and compression (kN)
with member numbers, the largest utilisation and its combination, the
critical load factor, a 3D view coloured by axial force (turn and zoom with
the mouse or by touch), member and node tables with units, and CSV downloads.

When several struts share the extreme value because of the symmetry of the
dome, the app names the lowest member number and says how many struts share
it (the desktop program reports whichever one round-off puts first).

## Run it on your own computer

From the repository root:

```
pip install -r web/requirements.txt
streamlit run web/streamlit_app.py
```

The browser opens at <http://localhost:8501>.

## Deploy on Streamlit Community Cloud (three steps)

1. Sign in at <https://share.streamlit.io> with your own GitHub account and
   allow it to read the repository.
2. Choose **Create app**, then **Deploy a public app from GitHub**, and enter
   repository `covdetbaba/GeoDomeLab`, branch `main`, main file path
   `web/streamlit_app.py`. You can pick an address such as
   `geodomelab.streamlit.app` under **App URL**.
3. Press **Deploy**. The first build takes a few minutes. Cloud reads
   `web/requirements.txt` (next to the main file), not the desktop
   requirements in the root folder.

A Community Cloud app goes to sleep after some days without visitors. Open
it once shortly before a presentation so the audience does not wait for it to
wake up.

## Rebuild the GitHub Pages version

`build_pages.py` writes a static site that runs the app with
[stlite](https://github.com/whitphx/stlite) (Streamlit on Pyodide). It copies
only the geodome modules the app imports (found from the import statements),
`validation/hand_calc.py` and an `index.html` that mounts the app with a
pinned `@stlite/browser` version.

```
python web/build_pages.py --out ../geodomelab-site --check
python -m http.server 8600 --directory ../geodomelab-site
```

`--check` runs both examples from the built folder alone, so a missing module
shows up before the site is published. To publish, commit the folder to the
`gh-pages` branch. GitHub Pages serves that branch at
<https://covdetbaba.github.io/GeoDomeLab/>.

The demo video page <https://covdetbaba.github.io/GeoDomeLab/demo/> is added
with `--video GeoDomeLab_demo.mp4 --poster poster.jpg`. The video (31 MB) is
not kept in the repository. Download it from the
[v2.0.0 release](https://github.com/covdetbaba/GeoDomeLab/releases/tag/v2.0.0)
to rebuild the page.

## Files

| File | Purpose |
| --- | --- |
| `streamlit_app.py` | The app (Streamlit, imports `geodome`) |
| `requirements.txt` | Packages for Streamlit Community Cloud or a local run |
| `build_pages.py` | Builds the no-install GitHub Pages site (stlite) |
