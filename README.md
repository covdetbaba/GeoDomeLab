# GeoDomeLab — Geodesic Dome Space-Truss Analyzer

A matrix-structural-analysis tool for single-layer Class I (icosahedral) geodesic
domes, written in Python with a Qt graphical interface.

The user defines the dome (frequency, radius, truncation, supports), the members
(material and a circular hollow section per strut family) and the loading (dead,
snow, live, wind, point loads and LRFD combinations). The program generates the
geometry, solves it by the **direct stiffness method** for three-dimensional
space trusses, checks every strut against **AISC 360 LRFD**, runs a **linear
buckling eigenanalysis**, optionally **auto-sizes** the strut families to minimum
weight, and exports a complete calculation report.

---

## 1. Quick start

### In Spyder

1. Open `run_geodome.py`.
2. Press **F5**.

If the window freezes, set the console's graphics backend to Qt:
`Tools → Preferences → IPython console → Graphics → Backend: Qt5` (or
`Automatic`), then restart the kernel. With the *Inline* backend the program
still runs correctly, but Spyder blocks until the window is closed.

### From a terminal

```bash
python run_geodome.py                      # start the GUI
python run_geodome.py --verify             # run the 75-check verification suite
python run_geodome.py --demo               # run the worked case study, headless
python run_geodome.py --batch 6 15 0.5 out # analyse and export without a GUI
python tests/test_geodome.py               # run the unit tests (no pytest needed)
pytest -q                                  # ... or with pytest
```

## 2. Requirements

| Package    | Minimum | Notes                                        |
|------------|---------|----------------------------------------------|
| Python     | 3.9     | developed on 3.13                            |
| numpy      | 1.20    |                                              |
| scipy      | 1.7     | sparse solver and eigensolver                |
| matplotlib | 3.5     | `backend_qtagg`, `mpl_toolkits.mplot3d`      |
| qtpy       | 2.0     | binding-independent Qt layer                 |
| PyQt5 **or** PyQt6 **or** PySide6 | — | any one of them |

Everything above ships with a standard Anaconda/Spyder installation. `qtpy`
selects whichever binding is present, so the same source runs unmodified on
PyQt5, PyQt6 and PySide6.

```bash
pip install -r requirements.txt
```

## 3. What the program does

| Stage | Method |
|-------|--------|
| Geometry | Class I icosahedral breakdown: each of the 20 faces split into ν² triangles by barycentric interpolation, projected onto the sphere, duplicate vertices merged with a KD-tree |
| Truncation | Cut along an existing ring of struts; supports placed on the **topological free edge** of the retained cap |
| Analysis | 3D space truss, 3 DOF/node, sparse assembly, partitioned solve, reactions by back-substitution |
| Loads | Dead on the faceted surface, snow/live on the projected area, wind normal to the surface with a meridional *C*ₚ distribution; ASCE 7 LRFD combinations enveloped |
| Member design | AISC 360 Chapter D (tension) and E3/E7 (compression, round HSS) |
| Stability | Linear buckling, `[K + λ Kg]{φ} = 0`, recast as a symmetric generalised eigenproblem |
| Optimisation | Fully stressed design on the strut families, subject to utilisation, *KL/r* ≤ 200 and a deflection limit |
| Output | Interactive 3D views, sortable tables, CSV, PNG, JSON model, self-contained HTML report |

## 4. Repository layout

```
GeoDomeLab/
├── run_geodome.py            launcher (GUI / --verify / --demo / --batch)
├── requirements.txt
├── README.md
├── geodome/
│   ├── geometry.py           Class I generation, truncation, chord factors
│   ├── sections.py           CHS catalogue and section properties
│   ├── materials.py          steel and aluminium grades
│   ├── model.py              DomeModel container, JSON save/load
│   ├── loads.py              dead / snow / live / wind, LRFD combination
│   ├── solver.py             direct stiffness solver, geometric stiffness
│   ├── buckling.py           AISC member checks, eigenvalue buckling
│   ├── optimize.py           fully stressed auto-sizing
│   ├── postprocess.py        multi-case envelope, summary, result tables
│   ├── plotting.py           matplotlib 3D renderers
│   ├── report.py             CSV / PNG / HTML export
│   ├── verification.py       75-check benchmark suite
│   └── gui/                  Qt front end (canvas, panels, worker, window)
├── examples/worked_example.py   the case study used in the report and demo
├── tests/test_geodome.py        30 unit tests
├── docs/                        proposal, final report, presentation, user guide
└── output/                      generated results (safe to delete)
```

## 5. Units

Everything inside the program is **SI**: metres, newtons, pascals, kilograms.
The interface converts at its boundary only — pressures are entered in kPa,
point loads in kN, and results are displayed in mm, kN and MPa.

## 6. Verification

`python run_geodome.py --verify` runs 75 independent checks:

* closed-form bar, tripod and two-bar solutions
* published Class I chord-factor tables (2V and 4V)
* Euler's polyhedron formula and the vertex/strut counts at five frequencies
* Maxwell's count — the closed geodesic sphere is exactly isostatic
* positive definiteness of the supported stiffness matrix (no mechanisms)
* method-of-sections free bodies at three cut heights
* membrane shell theory, *n*φ = *pR*/2, at three heights
* five-fold rotational symmetry of the displacement field
* equilibrium, load/reaction balance and strain-energy identities for every case
* linearity of the buckling eigenproblem

All 75 pass. The full table is reproduced in `docs/final_report.md` and can be
exported from the GUI.

### Hand-calculation validation

`validation/` holds a check against hand calculations of two small geodesic
models, built by the program's own generator at 1V:

* **Model A**, 1V pentagonal cap: 6 nodes, 10 struts, 3 free DOF, solved by
  the direct stiffness method with every matrix written out.
* **Model B**, 1V dome on its lower ring: 11 nodes, 25 struts, 18 free DOF,
  solved by virtual work with five-fold symmetry (3 unknowns).

```bash
python validation/hand_calc.py               # the hand calculation (no geodome import)
python validation/run_validation.py          # program vs hand, 80 quantities
python validation/open_validation_in_gui.py B1   # open a validation model in the GUI
```

All 80 compared quantities agree; the largest difference is 1.2e-13 %.
The comparison table, sketches, plots and GUI-ready model files are written to
`validation/output/`. Nothing in `geodome/` was changed for this check.

**GUI fix, 27 October 2026.** *Generate dome* failed silently when the new
dome had more strut families than the one on screen (for example 1V or 4V
followed by 6V): the section rows of the old dome were read back before they
were rebuilt. `on_generate` in `geodome/gui/main_window.py` now rebuilds them
first. Only the order of two lines changed, and only in the GUI; the solver,
geometry and design code are untouched, and the 30 tests, 75 checks and 80
hand-calc comparisons give the same results. The original file is kept in
`_revert/main_window_2026-10-01.py`. `tests/test_gui_generate.py` clicks
Generate dome through the real window for every frequency from 1V to 10V,
up and down; it fails on the original file and passes on the fixed one:

```bash
python tests/test_gui_generate.py
```

### Proposal documents

The proposal deliverables sit in this folder and are rebuilt by scripts in
`docs/build/` (Word and PowerPoint must be installed; the PDFs are exported
by Word and PowerPoint themselves):

| File | Rebuild with |
|---|---|
| `GeoDomeLab_Proposal.docx` / `.pdf` (CE 7405 written-proposal template) | `python docs/build/proposal_figures.py` then `python docs/build/build_proposal.py` |
| `GeoDomeLab_Validation_Report.docx` / `.pdf` | `python docs/build/build_report.py` |
| `GeoDomeLab_Proposal_Presentation.pptx` (speech in the speaker notes) | `python docs/build/build_presentation.py` |
| `GeoDomeLab_Presentation_Speech.docx` / `.pdf` | `python docs/build/build_speech.py` |

The written proposal replaces the older `docs/proposal.md`. One reference list,
`docs/build/references.py`, feeds the proposal, the report and the slides;
every entry was checked against library or publisher records. The speech text
lives in `docs/build/speech.py`; edit it there and rebuild both the
presentation and the speech so they stay identical. The student name, LSU ID
and dates are set at the top of `build_proposal.py` (`STUDENT`, `LSU_ID`,
`PRES_DAY`, `SUBMITTED`), `build_report.py` (`AUTHOR`) and `speech.py`
(`PRESENTER`). The proposal copies the template's look exactly (Times New
Roman 11 pt, LSU purple tables, the template footer); the 10-26 version is
kept in `_revert/`.

## 7. AI assistance

This project was developed with the assistance of Anthropic's Claude, in line
with the course AI policy. See `docs/ai_acknowledgment.md` for the full
disclosure of what was AI-assisted and what was independently verified.
