# GeoDomeLab 2.0 user guide

This guide covers every menu, panel and tab of the GeoDomeLab window, the two
input file formats, and the analyses. The methods behind the numbers are in
[THEORY.md](THEORY.md). How the program is checked is in
[VALIDATION.md](VALIDATION.md).

Contents

1. [Install and start](#1-install-and-start)
2. [Units](#2-units)
3. [The window](#3-the-window)
4. [Your first analysis, step by step](#4-your-first-analysis-step-by-step)
5. [Input tabs](#5-input-tabs)
6. [View tabs](#6-view-tabs)
7. [Result tabs](#7-result-tabs)
8. [Menus and toolbar](#8-menus-and-toolbar)
9. [Nonlinear analysis and the imperfection study](#9-nonlinear-analysis-and-the-imperfection-study)
10. [Modal analysis](#10-modal-analysis)
11. [Fabrication schedules](#11-fabrication-schedules)
12. [Parameter study](#12-parameter-study)
13. [Validation menu](#13-validation-menu)
14. [Export: report, CSV, figures, CAD and OpenSees](#14-export-report-csv-figures-cad-and-opensees)
15. [Model files: JSON format version 2](#15-model-files-json-format-version-2)
16. [Importing a truss from CSV files](#16-importing-a-truss-from-csv-files)
17. [Running without the window](#17-running-without-the-window)
18. [Troubleshooting](#18-troubleshooting)

---

## 1. Install and start

You need Python 3.10 or newer.

```bash
git clone https://github.com/covdetbaba/GeoDomeLab.git
cd GeoDomeLab
pip install -r requirements.txt
python run_geodome.py
```

Other ways to start it:

| How | Command |
|---|---|
| Installed command | `pip install -e .` once, then `geodomelab` |
| As a module | `python -m geodome` |
| Conda | `conda env create -f environment.yml`, `conda activate geodomelab`, then `python run_geodome.py` |
| Spyder | open `run_geodome.py` and press F5 |

In Spyder, set Tools > Preferences > IPython console > Graphics > Backend to
Qt (or Automatic) if the window freezes, then restart the kernel. With the
Inline backend the program still works, but Spyder waits until you close the
window.

The window opens with a 4V dome of 12 m radius. It is already generated, so
you can run an analysis at once.

The OpenSees cross-check is optional. It needs `pip install openseespy`, or a
separate Python with openseespy named in the environment variable
`GEODOME_OPENSEES_PYTHON`.

## 2. Units

The window uses engineering units. Inside the program, in model files and in
the Python API everything is SI.

| Quantity | In the window | In model files and the API |
|---|---|---|
| Length, coordinates | m | m |
| Displacement | mm | m |
| Settlement | mm | m |
| Force, point and nodal loads | kN | N |
| Pressure (dead, snow, live) | kPa | Pa |
| Spring stiffness | kN/m | N/m |
| Stress, E, Fy | MPa (results), GPa (E) | Pa |
| Temperature change | °C | °C |
| Wind speed | m/s | m/s |
| Mass | kg | kg |
| Frequency, period | Hz, s | Hz, s |

Tension is positive. The z axis points up. Nodes are numbered from 0. The
degree of freedom (DOF) number of node n is 3n (x), 3n + 1 (y) and 3n + 2 (z).

The CSV import uses m, kN and kN/m.

## 3. The window

```
+------------------+----------------------------------------------------+
| Input tabs       | View tabs                                          |
|                  | Model, Loads, Deflected shape, Axial force,        |
|  Geometry        | Utilisation, Buckling mode, Charts, Nonlinear,     |
|  Sections        | Modes, Study                                       |
|  Loads           +----------------------------------------------------+
|  Design          | Result tabs                                        |
|                  | Summary, Members, Nodes, Reactions, Auto-sizing,   |
|                  | Verification, Nonlinear, Imperfection, Modes,      |
|                  | Fabrication, Study, Validation, Log                |
+------------------+----------------------------------------------------+
| status bar: messages, progress bar, Stop / Cancel during long runs    |
+-----------------------------------------------------------------------+
```

Menus: File, Examples, Analysis, View, Validation, Help.
Toolbar: Generate dome, Run analysis, Auto-size, Verify, Export report,
Nonlinear, Modes.

Drag the splitters to give the plots or the tables more room. Every plot has
the matplotlib toolbar: drag to rotate a 3D view, zoom, pan and save the
figure.

Long runs (analysis, auto-sizing, nonlinear, modal, study, validation) run in
the background. The window stays responsive, and the inputs and the File
actions New, Open, Save and Import are locked until the run ends.

Results always belong to the inputs that made them. A change of section,
material or supports, or an auto-size, clears the results, including the
modal and nonlinear ones: press F9 again. A change in the Loads or Design
tab is used at the next run, and the status bar says so. When the supports
cannot hold the structure (a mechanism), the analysis stops with a message
that says what is free to move.

The window can shrink to about 1090 x 610 pixels, so it fits a 1366 x 768
screen; the input tabs and the Nonlinear view scroll when they need more
room.

## 4. Your first analysis, step by step

This example analyses a 6V hemisphere of 15 m radius (30 m span).

1. **Geometry tab.** Set Frequency to 6V, Sphere radius to 15 m and Dome
   fraction to 0.50. Keep Supports at "Pinned base ring". Press
   **Generate dome** (F5). The table under it lists the 9 strut families,
   their lengths, chord factors and counts.
2. **Sections tab.** Keep the grade ASTM A500 Gr. C (round). Choose
   CHS 114.3 x 5.0 in the "Apply to all" box and press **Apply to all**.
   The structural mass at the bottom updates at once.
3. **Loads tab.** Set Superimposed dead to 0.35 kPa and Snow to 0.60 kPa.
   Keep wind on, set the basic wind speed to 58 m/s and Exposure Kz to 0.95.
   Keep the default combinations (LC1 to LC6 and SLS).
4. **Design tab.** Keep K = 1.0, span/240, φ = 0.90 and the buckling
   analysis on.
5. Press **Run analysis** (F9). The status bar shows the progress. When it
   ends, the Summary tab reports the governing member, the governing
   combination, the largest utilisation, the deflection check and the
   buckling load factor λ_cr.
6. Look at the **Axial force** view (meridians in compression, lower hoops in
   tension), the **Utilisation** view and the **Buckling mode** view.
7. Press **Auto-size** (F10) to find lighter sections. The Auto-sizing tab
   shows each iteration and the mass.
8. Press **Export report** (Ctrl+R) to write one HTML file with every table
   and figure.

Examples > Case study - 6V dome, 30 m span loads this same model in one
click. It keeps its own combination list (LC1 to LC4 and SLS, the list of the
report), so its trial utilisation is 0.094 where the default list, with
1.2D + 1.6S + 0.5W, gives 0.103.

## 5. Input tabs

### 5.1 Geometry

| Input | Meaning |
|---|---|
| Frequency | Class I icosahedral frequency ν, from 1V to 10V. A full sphere has 10ν² + 2 nodes and 30ν² struts, so the model grows fast. 2V is a hand-check size, 4V a small dome, 6V to 8V a large one. |
| Sphere radius | Radius of the sphere the nodes lie on, 1 to 200 m. |
| Dome fraction | Part of the sphere's diameter that is kept, measured down from the crown. 0.50 is a hemisphere. Range 0.25 to 0.95. |
| Supports | How the base is held (see below). |
| Generate dome | Builds the geodesic dome from these inputs (F5). |

The requested cut is snapped to the nearest ring of real nodes, so the base
ring is always a closed polygon. Even frequencies have an equator ring. Odd
frequencies have none, and a cut near 0.50 snaps noticeably. The Log says
which fraction was used.

Support types:

| Support | Meaning |
|---|---|
| Pinned base ring | ux, uy and uz are fixed at every base node. |
| Vertical rollers + three lateral restraints | uz is fixed at every base node, and three nodes at about 120° hold the dome sideways. The dome can spread, as membrane theory assumes. |
| Elastic springs | A spring kx, ky, kz [kN/m] at every base node. 0 means free in that direction. The kx, ky, kz boxes appear only for this choice. |
| Custom per-node supports | Fixed, free or spring in each direction of each listed node. They come from a model file or a CSV import and are shown read-only. |

The **Generated geometry** box lists the size of the model: frequency,
radius, fraction, base diameter, crown height, nodes, struts, panels, strut
families, supported nodes, material, surface and plan areas. The table under
it lists each strut family with its length, chord factor and count.

**Custom trusses.** A model that is not a geodesic dome (an example, a model
file or a CSV import) is a custom truss. A note at the top shows its name and
size. The geodesic inputs are greyed out. Chord factors show "-" when the
nodes do not lie on a sphere. Generate dome or File > New returns to a
geodesic dome.

### 5.2 Sections

| Input | Meaning |
|---|---|
| Grade | Material: ASTM A500 Gr. B or C (round HSS), ASTM A53 Gr. B pipe, EN S235, EN S355 or aluminium 6061-T6. The line below shows E, Fy, density and the thermal coefficient α. The member checks follow AISC 360, a steel code: for aluminium they do not apply, and the window and every run say so. |
| Apply to all | Sets the chosen section on every strut family. |
| Section per strut family | One circular hollow section (CHS) per family. The catalogue runs from CHS 26.9 x 2.3 to CHS 323.9 x 12.5. |

The structural mass and the mass per m² of plan area update as you change a
section. Sections are set per family, not per member, because domes are
built that way. It also keeps auto-sizing to a few design variables.

### 5.3 Loads

| Group | Inputs |
|---|---|
| Gravity loads | Include strut self weight (half of each strut's weight to each end node). Superimposed dead [kPa] on the faceted surface area. Snow [kPa] on the horizontally projected area of the panels that face up, with an optional unbalanced split (1.5 on the half x ≥ 0, 0.5 on the other half; a simple test pattern, not the ASCE 7 curved-roof rule). Live [kPa] on the projected area; it acts only in combinations with an L factor (LC6 in the default list). Crown point load [kN], downward at the crown node, in case L. |
| Wind | Basic wind speed V [m/s] (3-second gust), exposure Kz, topographic Kzt, directionality Kd, gust factor G and the wind direction [°]: the direction the wind blows towards, measured from +x (0 = wind from −x). The line below shows qz and qz·G. Untick the group title to switch wind off. |
| Temperature and settlement (load T) | Temperature change ΔT [°C] of every strut, + is heating. A settlement table: node, ux, uy, uz [mm]. Use Add and Remove. Only supported nodes, and only their supported directions, can move. |
| Nodal loads | Node, Fx, Fy, Fz [kN] along the global axes (z up) and the load case D, S, L or W that the load belongs to. |
| Load combinations | Name and the factors on D, S, L, W and T. Tick Service to use a row for the deflection check instead of the strength envelope; at least one row must stay a strength row. Add, Remove, and Reset to ASCE 7 (restores LC1 1.4D, LC2 1.2D + 1.6S, LC3 1.2D + 1.0W + 0.5S, LC4 0.9D + 1.0W, LC5 1.2D + 1.6S + 0.5W, LC6 1.2D + 1.6L + 0.5S and SLS 1.0D + 1.0S, from ASCE 7-16 Sec. 2.3.1). With a live load, ASCE 7 also adds 1.0L to combinations 3 and 4: add those rows. The Log names any load case that no strength combination uses. |

Temperature and settlement together are the self-straining load T. They act
only through the T factor of a combination. Use a factor of at least 1.0 on T
(ASCE 7-16 Sec. 2.3.4).

A model without panels (for example a CSV import without faces.csv) has no
area, so dead, snow, live and wind pressures give no load. A note says so.
Use self weight and nodal loads for such a model.

Tick "Node numbers" in the Model view to see the node numbers. When you run,
the program checks every node number, every settlement direction and every
combination name. A bad entry gives a message that names the table and the
row, and nothing is run.

### 5.4 Design

| Group | Inputs |
|---|---|
| Design parameters (AISC 360 LRFD) | Effective length factor K (1.0 suits pin-ended struts). Deflection limit, as span/n (default span/240). Resistance factors φ for compression and tension (0.90). Design wall thickness factor, 1.00 by default. AISC 360-16 B4.2 uses 0.93 for ERW HSS such as ASTM A500. It affects the strength checks only. |
| Stability analysis | Run the linear buckling eigenanalysis, and the number of modes to extract (1 to 20). |
| Auto-sizing (fully stressed design) | Target utilisation (0.95), maximum iterations (12), slenderness cap KL/r ≤ 200, and whether to enforce the deflection limit too. |
| Buttons | Run analysis (F9), Auto-size the strut families (F10), Run verification suite (F11). |

## 6. View tabs

| View | What it shows |
|---|---|
| Model | The struts coloured by strut family, nodes, supports (black triangles). Tick Cladding panels, Nodes or Node numbers. |
| Loads | The nodal load arrows of the chosen combination. On a wind case the arrows near the crown point outward. That is suction, and it is correct. |
| Deflected shape | The magnified displaced shape over the undeformed one, coloured by displacement. Choose the case, the magnification or press Auto scale. |
| Axial force | Blue is compression, red is tension. Line width can follow the force. Choose the case. |
| Utilisation | Demand over capacity of every strut, enveloped over the strength combinations. Struts above 1.0 are drawn thick and dark red. |
| Buckling mode | The buckling mode shape. Step through the modes and set the drawn amplitude. Pairs of nearly equal λ come from the symmetry of the dome. |
| Charts | The strut schedule bar chart and the utilisation histogram. |
| Nonlinear | The load path and the deformed shape at any step (Section 9). Sub-tabs Load path and Imperfection study. |
| Modes | The mode shapes and their animation (Section 10). |
| Study | Four charts of a parameter study (Section 12). |

View > Save current view as PNG saves the figure of the current view tab.

## 7. Result tabs

| Tab | Content |
|---|---|
| Summary | The calculation narrative: geometry, load totals, one row per combination, the governing member and combination, the deflection check, λ_cr and any warnings. |
| Members | Every strut: nodes, family, section, length, KL/r, envelope maximum and minimum axial force over the strength combinations (N max and N min; a strut always in compression shows two negative values), both design strengths, utilisation and OK or OVER. Rows over 1.0 are shaded red. |
| Nodes | Coordinates, ring number and the three displacements. |
| Reactions | Support forces per node (spring forces included) with a total row. |
| Auto-sizing | One row per iteration of the fully stressed design: sections and mass. |
| Verification | The 75 built-in checks: computed value, reference, relative error, tolerance and PASS or FAIL. |
| Nonlinear | One row per converged step of the load path, and a text summary. Export path (CSV). |
| Imperfection | One row per imperfection amplitude. The first row is the perfect structure. Export table (CSV). |
| Modes | Mode, f, T, effective mass in x, y, z and the running sums. Double-click a mode to show it. Export CSV. |
| Fabrication | Hub offset, Refresh, Export cut list (CSV), and the sub-tabs Summary, Struts, Hubs, Panels and Dihedrals (Section 11). |
| Study | The parameter study table. Export CSV, Open selected row as model, New study. |
| Validation | Sub-tabs Hand calculations and OpenSees (Section 13). |
| Log | Every message of the session, with the details of any error. |

## 8. Menus and toolbar

### File

| Entry | Shortcut | Action |
|---|---|---|
| New model | Ctrl+N | Back to the default 4V dome of 12 m radius. |
| Open model... | Ctrl+O | Opens a GeoDomeLab model file (.json, Section 15). |
| Save model... | Ctrl+S | Saves every input as a .json model file. |
| Import truss from CSV... | | Reads a truss from CSV files (Section 16). |
| Export results to CSV... | | Members, nodes, reactions, strut schedule, forces per case and the summary. |
| Export figures (PNG)... | | The eight standard figures at 160 dpi. |
| Export HTML report... | Ctrl+R | One self-contained HTML calculation report. |
| Export everything to a folder... | | CSV, PNG, the HTML report and the model file. |
| Export geometry to CAD (DXF)... | | Struts, panels and hubs as an ASCII DXF file in metres. |
| Exit | Ctrl+Q | Closes the program. |

### Examples

Each entry loads a ready model with known answers. Hover over an entry to
read what it checks.

| Example | What it is |
|---|---|
| Validation Model A - 1V pentagonal cap | 6 nodes, 10 struts, 3 free DOF. Solved by hand in `validation/hand_calc.py`. |
| Validation Model B - 1V dome on its lower ring | 11 nodes, 25 struts, 18 free DOF. Solved by hand with five-fold symmetry. |
| von Mises two-bar truss (snap-through) | Half span 2.5 m, rise 0.125 m. Run the nonlinear analysis to see the snap-through. |
| 24-member shallow star dome (snap-through) | The classic benchmark of Greco & Venturini (2006) and Iu & Bradford (2012). The load factor reads 10⁴ P/EA. |
| Planar three-bar truss - temperature and settlement | An indeterminate truss with ΔT = +30 °C and a 5 mm settlement, checked by the force method. |
| Case study - 6V dome, 30 m span | The worked example: 6V hemisphere, R = 15 m, snow and wind. |

Loading an example replaces the model, refreshes every panel and view and
shows the model name in the window title.

### Analysis

| Entry | Shortcut | Action |
|---|---|---|
| Generate dome | F5 | Builds a geodesic dome from the Geometry tab. |
| Run analysis | F9 | Solves every combination, checks every strut and runs the buckling analysis. |
| Auto-size strut families | F10 | Fully stressed design, then a new analysis. |
| Modal analysis... | | Natural frequencies and mode shapes (Section 10). |
| Nonlinear analysis... | | Large-displacement load path through the limit point (Section 9). |
| Imperfection study... | | Limit loads of domes shaped like a buckling mode (Section 9). |
| Parameter study... | | Analyses a set of domes and compares them (Section 12). |

### View

One entry per view tab, and Save current view as PNG.

### Validation

| Entry | Shortcut | Action |
|---|---|---|
| Run verification suite (75 checks) | F11 | Closed-form, published and identity checks. Result tab Verification. |
| Run hand-calculation validation | | Models A and B and cases F1 to F10 against the hand calculations. |
| Run OpenSees cross-check | | Solves a set of models again in OpenSees and compares. Greyed out when no OpenSees is found. |
| Export OpenSeesPy script... | | A standalone OpenSeesPy script of the current model. |
| Open validation report (Markdown) folder | | Opens `validation/output`. |

### Help

| Entry | Content |
|---|---|
| Theory and method | A short page of the methods, equations and references. |
| CSV import format | The CSV format of Section 16. |
| About | Version, author, repository and license. |

### Toolbar

Generate dome (F5), Run analysis (F9), Auto-size (F10), Verify (F11), Export
report (Ctrl+R), Nonlinear and Modes. While a nonlinear run, a study or a
validation run is busy, a Stop or Cancel button appears in the status bar.

## 9. Nonlinear analysis and the imperfection study

**Analysis > Nonlinear analysis...** (or the Nonlinear button) traces the
load factor λ times one combination through the limit point. The bars are
corotational: displacements can be large, strains stay small and the
material stays elastic. Temperature and settlement (T) are left out of this
analysis.

Dialog inputs:

| Input | Meaning |
|---|---|
| Combination | λ multiplies the D, S, L and W loads of this combination. The default is the strength combination with the largest vertical load. |
| Method | Arc-length passes limit points and snap-through. Displacement control passes load maxima but stops at a snap-back. Load control stops at the first limit point. |
| Maximum steps | The run stops after this many steps (default 300). |
| Steps past the limit | Steps traced after the first limit point (default 15). Use more to follow the snap-through. |
| Stop at λ | Auto uses 2 λ_cr. |
| First increment Δλ | Auto uses λ_cr / 40. The step then adapts to the number of iterations. |
| Control DOF | The node and direction the path is plotted against (crown uz by default). Displacement control also steps it. It must be a free DOF. |
| Initial imperfection | Optional. Moves the nodes by a buckling mode shape scaled to the amplitude [mm]. A negative amplitude flips the mode. The dialog shows span/1000 and span/300 for scale. |

The **Nonlinear** view tab, sub-tab **Load path**, shows:

- the corotational path (blue), the P-Delta path (dashed) and the linear path
  (dotted), as λ against the control displacement,
- λ_cr as a dashed horizontal line, when it fits on the plot,
- the first limit point, marked and labelled with λ_limit and λ_limit/λ_cr,
- the deformed shape at the chosen step, coloured by displacement or axial
  force.

Drag the **Step** slider to move along the path. **Limit point** jumps to the
step before the first load maximum. **Play** runs along the path.
**Magnification**, **Auto scale** and **True scale** set the drawn
displacement. The line under the plots lists λ_limit, λ_cr, the ratio, the
number of steps and why the run stopped. The **Nonlinear** result tab lists
every step (λ, control displacement, iterations, residual, sign of det K_T,
number of negative eigenvalues, smallest eigenvalue, stiffness parameter,
strain energy).

**Analysis > Imperfection study...** moves the nodes by a buckling mode shape
and runs the nonlinear analysis for each amplitude. Both signs are run and
the lower limit load is kept.

| Input | Meaning |
|---|---|
| Combination | As above. |
| Buckling mode | The linear buckling mode used as the shape (1 to 6). |
| Amplitudes | Comma separated. `span/n` means the span divided by n. A plain number is in mm. Default span/1000, span/500, span/300, span/200. |
| Method, Maximum steps, Steps past the limit | As above (default 5 steps past each limit). |

The sub-tab **Imperfection study** plots λ_limit/λ_cr against amplitude/span.
The **Imperfection** result tab lists each amplitude, both signs, the
first critical λ and the local rise left at the hub where the mode is
largest (usually the crown). A positive amplitude always points the way the
load pushes. When an amplitude removes more than half of that hub's rise
above its neighbours, the Log says so: the limit load then measures a local
snap of that hub, so judge the amplitude against the hub's rise or the strut
length rather than the span (in the 6V case study the crown stands 0.198 m
above its neighbours, and span/300 = 100 mm removes 66 % of it).

Try it on Examples > 24-member shallow star dome: arc-length finds the limit
point at 10⁴ P/EA = 3.157, only 0.21 of λ_cr.

## 10. Modal analysis

**Analysis > Modal analysis...** (or the Modes button) solves
K φ = ω² M φ for the lowest modes.

| Input | Meaning |
|---|---|
| Number of modes | How many of the lowest modes to compute (default 12). |
| Include dead load mass | Adds the cladding dead load and the downward D nodal loads as mass (load / g). The strut mass is always included. |
| Snow mass fraction | Part of the balanced snow load taken as mass, 0 to 1. |
| Strut mass | Lumped (half of each strut at each end, the usual choice) or consistent (bar mass matrix, slightly higher frequencies). |

The **Modes** view shows the chosen mode on the undeformed structure, with
f, T and ω. **Animate** plays the mode as a motion. It stops when you change
the tab or the model. The **Modes** result tab lists f, T and the effective
modal mass in x, y and z with the running sums. The line above it says how
many modes reach 90 % of the mass. Double-click a row to show that mode.

## 11. Fabrication schedules

The **Fabrication** result tab is filled for every model and updates after
each change and each analysis.

| Sub-tab | Content |
|---|---|
| Summary | The counts of struts, cut lengths, hubs and panels, the tube mass, and notes. |
| Struts | The cut list: family, colour code, count, section, length, cut length, chord factor, axial angle, tangent angle, mass each and in total. |
| Hubs | Hub types by location (interior, edge, base), number of struts, strut mix (for example 2A+4B), the order of the families around the hub, count, the angular gaps and an example node. |
| Panels | Triangle types by the families of their sides, the three sides and angles, area, total area, and handedness for scalene types. |
| Dihedrals | The fold angle between the two panels on each strut family, and the bend angle. |

**Hub offset** [mm] is the distance from the hub centre to the end of the
tube. Cut length = strut length − 2 × hub offset. **Export cut list (CSV)...**
writes every table as a CSV file into a folder. Hub and connector mass is not
included in the tube mass.

## 12. Parameter study

**Analysis > Parameter study...** rebuilds the dome for every combination of
the values you give and runs the full analysis on each one. It uses the
current material, supports, loads and design settings.

| Input | Meaning |
|---|---|
| Frequencies [V] | Whole numbers, for example `2, 4, 6, 8`. |
| Radius values [m] | One value keeps the radius fixed. Several values sweep it. |
| Dome fraction values | One value or several, as above. 0.5 is a hemisphere. |
| Auto-size the strut families of every dome | Fully stressed design for each dome, with the Design tab settings. Slower, but compares economical designs. |
| Compute the buckling factor λ_cr | One linear buckling analysis per dome. |

Without auto-sizing every family of a dome gets one section, the heaviest one
of the current model. Settlements, nodal loads and custom supports are tied
to node numbers, so they do not carry over to a rebuilt dome. The Log says
so.

The **Study** view tab draws four charts against the parameter you choose in
"x axis": steel mass, largest utilisation, λ_cr and span/deflection. Filled
markers pass every check. Hollow markers fail one. The **Study** result tab
lists every dome. **Open selected row as model** loads that dome with its
sections. **Cancel study** in the status bar stops after the dome that is
running.

## 13. Validation menu

- **Run verification suite (75 checks)** runs the built-in checks
  (`geodome/verification.py`) and fills the Verification tab.
- **Run hand-calculation validation** compares the program with
  `validation/run_validation.py` (Models A and B, 80 quantities) and
  `validation/run_validation_ext.py` (cases F1 to F10, 422 comparisons). The
  table shows the hand value, the program value, the difference, the
  tolerance and PASS or FAIL. Tick "Show failures only" to filter. This run
  does not rewrite the files in `validation/output`.
- **Run OpenSees cross-check** runs the quick set of
  `validation/opensees_crosscheck.py`. It needs openseespy in this Python or
  in the Python named by `GEODOME_OPENSEES_PYTHON`. The OpenSees sub-tab shows
  one summary row per section and the comparisons below it. Click a section
  to show only its rows. Ctrl+click it to show every row again.
- **Export OpenSeesPy script...** writes a standalone script of the current
  model: linear static, modal (12 modes, lumped mass) or the nonlinear path
  of the crown uz. Run it with any Python that has openseespy.
- **Open validation report (Markdown) folder** opens `validation/output`,
  where the full reports live.

These entries need the `validation/` folder of the repository. When it is
missing they are greyed out with a tooltip.

## 14. Export: report, CSV, figures, CAD and OpenSees

| Output | How | Content |
|---|---|---|
| HTML report | File > Export HTML report (Ctrl+R) | One self-contained file: structure, loading, results, auto-sizing history, figures, solution checks and the verification suite. When you have run them, it also has sections for the nonlinear analysis, the imperfection study, the modal analysis, the fabrication schedule and the parameter study. |
| CSV tables | File > Export results to CSV | `dome_members.csv`, `dome_nodes.csv`, `dome_reactions.csv`, `dome_strut_schedule.csv`, `dome_case_forces.csv`, `dome_summary.txt`. |
| Figures | File > Export figures (PNG) | Model, loads, deflected shape, axial force, utilisation, buckling mode, strut chart, utilisation histogram. |
| Everything | File > Export everything to a folder | All of the above and `dome_model.json`. |
| CAD | File > Export geometry to CAD (DXF) | ASCII DXF (R12) in metres: one LINE per strut on layer `STRUT_<family>`, one 3DFACE per panel on layer `PANELS`, one POINT per node on layer `HUBS`. After an analysis, pick the file type "struts coloured by utilisation" in the save dialog to colour the struts by utilisation. |
| Fabrication | Fabrication tab > Export cut list (CSV) | One CSV per fabrication table. |
| Study, modes, path | Export buttons on those tabs | CSV. |
| OpenSees | Validation > Export OpenSeesPy script | A standalone Python script. |

The HTML report needs no other file, so you can send it on its own.

## 15. Model files: JSON format version 2

File > Save model writes a JSON file. File > Open reads it back exactly.
Version 1 files (GeoDomeLab 1.x) still open. All values are SI: m, N, Pa,
N/m and °C.

A geodesic dome (two of the seven default combinations are shown):

```json
{
  "format": "GeoDomeLab model",
  "version": 2,
  "geometry": {"kind": "geodesic", "frequency": 4, "radius": 12.0,
               "dome_fraction": 0.5},
  "material": "ASTM A500 Gr. C (round)",
  "group_sections": ["CHS 88.9 x 4.0", "CHS 88.9 x 4.0", "CHS 88.9 x 4.0",
                     "CHS 88.9 x 4.0", "CHS 88.9 x 4.0", "CHS 88.9 x 4.0"],
  "support_type": "pinned",
  "spring_stiffness": [50000000.0, 50000000.0, 50000000.0],
  "custom_supports": {},
  "loads": {
    "include_self_weight": true,
    "dead_pressure": 300.0,
    "snow_pressure": 1200.0,
    "snow_unbalanced": false,
    "live_pressure": 0.0,
    "wind_enabled": true,
    "wind_speed": 45.0, "wind_Kz": 0.9, "wind_Kzt": 1.0, "wind_Kd": 0.95,
    "wind_G": 0.85, "wind_azimuth": 0.0,
    "crown_point_load": 0.0,
    "temperature_change": 0.0,
    "settlements": [],
    "nodal_loads": []
  },
  "combinations": [
    {"name": "LC2  1.2D + 1.6S", "D": 1.2, "S": 1.6, "L": 0.0, "W": 0.0,
     "serviceability": false, "T": 0.0},
    {"name": "SLS  1.0D + 1.0S", "D": 1.0, "S": 1.0, "L": 0.0, "W": 0.0,
     "serviceability": true, "T": 0.0}
  ],
  "design": {"K_effective": 1.0, "deflection_limit": 240.0, "phi_c": 0.9,
             "phi_t": 0.9, "design_thickness_factor": 1.0}
}
```

A custom truss replaces the `geometry` block:

```json
"geometry": {
  "kind": "custom",
  "title": "Planar three-bar indeterminate truss",
  "nodes": [[-4.0, 0.0, 0.0], [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 3.0]],
  "members": [[0, 3], [1, 3], [2, 3]],
  "faces": [],
  "base_nodes": [0, 1, 2],
  "groups": [0, 1, 2]
}
```

Keys:

| Key | Meaning |
|---|---|
| `geometry.kind` | `"geodesic"` (frequency, radius [m], dome_fraction) or `"custom"`. |
| `nodes` | Coordinates [m], one `[x, y, z]` per node, numbered from 0. |
| `members` | Node pairs `[i, j]`, 0-based. |
| `faces` | Triangles `[i, j, k]` for area loads. May be empty. |
| `base_nodes` | Nodes supported by the pinned, roller or elastic support types. |
| `groups` | Strut family of each member (0, 1, 2, ...). |
| `material` | A grade name from the Sections tab. |
| `group_sections` | One catalogue section per strut family, in family order. |
| `support_type` | `"pinned"`, `"roller3"`, `"elastic"` or `"custom"`. |
| `spring_stiffness` | `[kx, ky, kz]` in N/m, for `"elastic"`. |
| `custom_supports` | For `"custom"`: `{"node": [cx, cy, cz]}` with `true` = fixed, `false` = free, a positive number = spring in N/m. The node keys are strings. |
| `loads` | Pressures in Pa, crown point load in N, wind speed in m/s, azimuth in degrees, ΔT in °C. |
| `settlements` | `[{"node": 2, "ux": 0.0, "uy": 0.0, "uz": -0.005}]` in m. |
| `nodal_loads` | `[{"node": 3, "Fx": 15000.0, "Fy": 0.0, "Fz": -60000.0, "case": "L"}]` in N. |
| `combinations` | Name, factors on D, S, L, W and T, and `serviceability`. |
| `design` | K, deflection limit (span/n), φ_c, φ_t and the design wall thickness factor. |

Ready examples are in `examples/models/`. The Examples menu builds the same
models in code, so it does not need these files. In Python,
`geodome.benchmarks.load_preset_file(key)` reads them from the `examples/models`
folder next to the `geodome` package (in a clone or an editable install), or
from `examples/models` under the current folder.

## 16. Importing a truss from CSV files

File > Import truss from CSV... reads any pin-jointed space truss from a few
small tables. Pick the nodes file. The other files are read from the same
folder. A common prefix is allowed: picking `tower_nodes.csv` also reads
`tower_members.csv`, `tower_supports.csv` and so on.

Separator: comma, semicolon or tab. A header row is optional. Blank lines and
lines that start with `#` are skipped. With a semicolon separator a decimal
comma also works. Units: m, kN and kN/m.

| File | Columns | Notes |
|---|---|---|
| `nodes.csv` (required) | `id, x, y, z` | Any whole-number ids. Coordinates in m, z up. |
| `members.csv` (required) | `i, j[, group]` | Node ids at the two ends. The optional group is the strut family. Give it on every row or on none. Without groups the families are formed by length. |
| `supports.csv` (optional) | `node, cx, cy, cz` | Per direction: 1 = fixed, 0 = free, another positive number = spring in kN/m. The words `fixed` and `free` also work. Without this file the lowest nodes are pinned. |
| `loads.csv` (optional) | `node, Fx, Fy, Fz[, case]` | Force in kN along the global axes. Case D, S, L or W (default L). |
| `faces.csv` (optional) | `i, j, k` | Triangular panels. Needed only for area loads. |

Example, a pyramid roof (the full set is in `examples/csv_truss/`):

`nodes.csv`
```
id, x, y, z
1, -3.0, -3.0, 0.0
2,  3.0, -3.0, 0.0
3,  3.0,  3.0, 0.0
4, -3.0,  3.0, 0.0
5,  0.0,  0.0, 2.5
```

`members.csv`
```
i, j, group
1, 5, 1
2, 5, 1
3, 5, 1
4, 5, 1
1, 2, 2
2, 3, 2
3, 4, 2
4, 1, 2
```

`supports.csv` (node 4 sits on a 20 000 kN/m spring in z)
```
node, cx, cy, cz
1, 1, 1, 1
2, 0, 1, 1
3, 0, 0, 1
4, 1, 0, 20000
```

`loads.csv`
```
node, Fx, Fy, Fz, case
5, 0, 0, -20, D
5, 0, 0, -30, S
5, 5, 0, 0, W
```

`faces.csv` (the four sloping roof panels, for area loads)
```
i, j, k
1, 2, 5
2, 3, 5
3, 4, 5
4, 1, 5
```

The ids are renumbered 0, 1, 2, ... in id order, because the program numbers
nodes from 0. The Log says so when the ids change. The imported model starts
with self weight on, area loads at zero, wind off, the default section on
every family and the ASCE 7-16 strength combinations with live load. Check
the Sections and Loads tabs, then run the analysis (F9). A file error names
the file and the line.

## 17. Running without the window

```bash
python run_geodome.py --verify               # the verification suite
python run_geodome.py --demo                 # the worked case study
python run_geodome.py --batch 6 15 0.5 out   # frequency, radius, fraction, folder
python run_all_checks.py                     # every test and validation
```

`geodomelab` and `python -m geodome` take the same options.

From Python (SI units):

```python
from geodome.model import default_model
from geodome.optimize import auto_size
from geodome.postprocess import run_analysis, summary_text
from geodome import report as rp

model = default_model(6, 15.0, 0.5)        # frequency, radius [m], fraction
model.loads.snow_pressure = 600.0          # Pa
auto_size(model, target_util=0.95)
res = run_analysis(model)
print(summary_text(res))
rp.export_html(res, "report.html")
```

The other modules work the same way, for example
`geodome.nonlinear.run_nonlinear(model)`,
`geodome.modal.modal_analysis(model)` and
`geodome.fabrication.fabrication_tables(model)`. Their docstrings give the
options and the equations.

## 18. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| The window opens, then Spyder hangs | The Inline graphics backend. Set it to Qt (or Automatic) and restart the kernel. |
| `ModuleNotFoundError: qtpy` or "No Qt bindings could be found" | Run `pip install -r requirements.txt` (it installs qtpy and PyQt6). |
| Windows: pip stops with "No such file or directory" for a file deep inside `PyQt6\Qt6\qml` | The path is longer than 260 characters. Create the environment in a folder with a shorter path, or turn on Windows long path support. |
| Linux: "Could not load the Qt platform plugin" | Install the system Qt libraries, for example `sudo apt-get install libegl1 libgl1 libxkbcommon0 libfontconfig1 libdbus-1-3`. On a server without a display, set `QT_QPA_PLATFORM=offscreen`. |
| "Check the inputs" message | A node number, settlement, nodal load or combination is wrong. The message names the table and the row. Nothing was run. |
| "Unstable structure" or a zero-stiffness message | A node or direction has no stiffness, for example a custom truss with too few supports or a node with one bar. Add supports or members. |
| Buckling says there is no positive load factor | The reference combination is mostly tension, for example pure wind uplift. There is no compressive instability to find. |
| Auto-sizing stops below the target utilisation | The catalogue is discrete and slender struts are governed by buckling. The Log names the family. |
| Every strut is over 1.0 after auto-sizing | The heaviest section is not enough. A higher frequency gives shorter struts that buckle at higher loads. |
| The nonlinear run stops before the limit point | Raise Maximum steps or Stop at λ, or use arc-length. Displacement control stops where the control displacement snaps back. |
| Run OpenSees cross-check is greyed out | No openseespy was found. Run `pip install openseespy`, or set `GEODOME_OPENSEES_PYTHON` to a Python that has it. |
| The Validation entries are greyed out | The `validation/` folder is missing. Run the program from a clone of the repository. |
| `--demo` says the worked example was not found | `python -m geodome --demo` needs the `examples/` folder of the clone. Use an editable install (`pip install -e .`). |
| Plots look cramped | Drag the splitters between the inputs, the views and the tables. |
