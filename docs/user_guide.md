# GeoDomeLab — User Guide

---

## 1. Starting the program

**Spyder:** open `run_geodome.py` and press **F5**.
If the window freezes, set `Tools → Preferences → IPython console → Graphics →
Backend` to **Qt5** (or *Automatic*) and restart the kernel.

**Terminal:** `python run_geodome.py`

The program opens with a 4V dome of 12 m radius already generated, so there is
always something on screen to work with.

---

## 2. The window

```
┌───────────────┬──────────────────────────────────────────┐
│  INPUT DOCK   │           3D / chart views               │
│               │  Model · Loads · Deflected · Axial       │
│  Geometry     │  Utilisation · Buckling · Charts         │
│  Sections     ├──────────────────────────────────────────┤
│  Loads        │           RESULTS                        │
│  Design       │  Summary · Members · Nodes · Reactions   │
│               │  Auto-sizing · Verification · Log        │
└───────────────┴──────────────────────────────────────────┘
```

---

## 3. A typical session

### Step 1 — Geometry

| Field | Meaning |
|---|---|
| **Frequency** | Class I subdivision ν. 2V is a hand-check model, 4V a realistic small dome, 6V–8V a large one. A full sphere has 10ν² + 2 nodes and 30ν² struts. |
| **Sphere radius** | Radius of the generating sphere, in metres. |
| **Dome fraction** | Fraction of the sphere *diameter* kept, from the crown down. 0.50 is a hemisphere. |
| **Supports** | *Pinned base ring* fixes all three translations at every edge node. *Vertical rollers + three lateral restraints* lets the dome spread, which is the condition membrane shell theory assumes. |

Press **Generate dome** (or F5). The panel then shows the model summary and the
**strut schedule**: one row per length family, with its chord factor and count.
That table is the fabrication drawing — a 4V dome needs six strut lengths, not
250 different ones.

> **About the dome fraction.** A geodesic mesh can only be cut along its own
> strut lines, so the requested fraction is snapped to the nearest available
> cut. Even frequencies have an exact equator, so 0.50 always gives a flat
> springing line; other cuts may give a slightly stepped base, and the program
> says so in the Log. Supports are always placed on the true free edge of the
> mesh, whatever the cut.

### Step 2 — Sections

Choose a material grade, then a circular hollow section **per strut family**.
Use *Apply to all* to set a trial section everywhere. The structural mass updates
live at the bottom of the panel.

Sections are assigned by family rather than by member because that is how domes
are fabricated, and it is also what keeps auto-sizing to six to nine design
variables instead of hundreds.

### Step 3 — Loads

| Group | Notes |
|---|---|
| **Self weight** | Member mass, lumped half to each end node. |
| **Superimposed dead** | kPa on the *faceted surface* area — cladding, glazing, services. |
| **Snow** | kPa on the *horizontally projected* area, as codes define it. Tick *Unbalanced* for a 1.5 / 0.5 split across the dome. |
| **Live** | kPa on the projected area. |
| **Crown point load** | kN downward at the apex node. |
| **Wind** | Basic 3-second gust speed and the ASCE coefficients. The *q*z readout updates as you type. Untick the group header to switch wind off. |

**Load combinations** are fully editable. Type a name and the factors on D, S, L
and W. Tick **Service** to use a row for the deflection check instead of the
strength envelope. *Reset to ASCE 7* restores the defaults.

### Step 4 — Design

Effective length factor *K*, deflection limit, resistance factors, number of
buckling modes, and the auto-sizing settings.

Then:

* **Run analysis** (F9) — solves every combination, envelopes the strength cases,
  checks every strut, and runs the buckling eigenanalysis.
* **Auto-size** (F10) — fully stressed design on the strut families, then
  re-analyses automatically.
* **Verify** (F11) — runs the 75-check benchmark suite.

---

## 4. Reading the views

| View | What to look for |
|---|---|
| **Model** | Struts coloured by length family; the legend doubles as the strut schedule. Black triangles are supports. |
| **Loads** | Nodal load vectors for the selected combination. On a wind case the crown arrows point *outwards* — that is suction, and it is correct. |
| **Deflected shape** | Magnified displacements over a ghost of the undeformed shape. *Auto scale* picks a sensible magnification. |
| **Axial force** | Blue = compression, red = tension; line width scales with magnitude. A dome shows meridional compression everywhere and hoop tension in the lower rings. |
| **Utilisation** | Demand/capacity envelope over all strength combinations. Anything above 1.0 is drawn thick dark red. |
| **Buckling mode** | The critical mode shape. Step through modes with the spinner; near-identical λ values in pairs are the five-fold symmetry of the dome, not a bug. |
| **Charts** | Strut schedule bar chart and the utilisation histogram. |

Every view has the standard matplotlib toolbar — rotate by dragging, zoom, and
save the current figure.

---

## 5. Reading the results

**Summary** is the full calculation narrative: geometry, load resultants, a row
per combination, the governing member and combination, the serviceability check,
the buckling load factors, and any warnings.

**Members** lists every strut with its family, section, length, *KL*/*r*, the
envelope tension and compression, both capacities, the utilisation and an
OK/OVER flag. Overstressed rows are shaded red.

**Nodes** gives coordinates, ring number and the three displacement components.
**Reactions** gives the support forces and the horizontal thrust per node, with a
total row.

**Verification** shows each check with its computed value, its independent
reference, the relative error and the tolerance.

---

## 6. Exporting

| Menu item | Produces |
|---|---|
| Export results to CSV | members, nodes, reactions, strut schedule, per-case forces, summary text |
| Export figures (PNG) | the eight standard figures at 160 dpi |
| **Export HTML report** (Ctrl+R) | one self-contained file with all tables and figures embedded |
| Export everything to a folder | all of the above plus the JSON model |
| Save / Open model | JSON, round-trips the complete definition |

The HTML report needs no other files, so it can be submitted or e-mailed on its
own.

---

## 7. Running without the GUI

```bash
python run_geodome.py --verify               # verification suite
python run_geodome.py --demo                 # the worked case study
python run_geodome.py --batch 6 15 0.5 out/  # frequency, radius, fraction, folder
```

Or script it directly:

```python
from geodome.model import default_model
from geodome.postprocess import run_analysis, summary_text
from geodome.optimize import auto_size
from geodome import report as rp

model = default_model(frequency=6, radius=15.0, dome_fraction=0.5)
model.loads.snow_pressure = 900.0          # Pa — the API is strict SI
auto_size(model, target_util=0.95)
res = run_analysis(model)
print(summary_text(res))
rp.export_html(res, "report.html")
```

> **Units:** the Python API is strict SI — metres, newtons, pascals. Only the
> GUI uses kPa and kN.

---

## 8. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| Window opens then Spyder hangs | Inline graphics backend. Set the backend to Qt5 and restart the kernel. |
| `ModuleNotFoundError: qtpy` | `pip install qtpy` (Spyder normally provides it). |
| "Unstable structure: zero stiffness at DOF…" | A node has no stiffness in some direction. Regenerate the geometry; if it persists at a particular dome fraction, report it — the truncation logic is meant to prevent exactly this. |
| Buckling says "no positive load factor" | The reference case is predominantly tension, e.g. a pure uplift combination. Nothing is wrong; there is no compressive instability to find. |
| Auto-sizing stops well below the target utilisation | The catalogue is discrete and slender struts are buckling-governed. The Log explains which family and which next-smaller section failed. |
| Everything is over 1.0 after auto-sizing | The heaviest catalogue section is not enough. Raise the frequency (shorter struts buckle at higher loads) rather than reaching for bigger tubes. |
| Figures look cramped | Drag the splitter between the views and the result tables. |
