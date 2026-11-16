# Supplementary runs for the final report (Chapter 6)

These five scripts reproduce the supplementary numbers that Chapter 6 of the
GeoDomeLab 2.0 final report quotes. Each one uses the `geodome` package of
this repository and writes its results to a JSON file next to itself (or to
the path given as the first argument). The JSON files are git-ignored.

| Script | What it computes |
|---|---|
| `b42_wall.py` | Effect of the AISC 360-16 B4.2 design wall (0.93 t) on the sized 6V case study, with the same sections and after re-sizing |
| `freq_imperfection.py` | Perfect and imperfect limit loads of the auto-sized 4V, 6V and 8V hemispheres (R = 15 m, LC2) with span-based and hub-scale (10 mm) amplitudes, and the member-strength load factor of each |
| `hub_tolerance.py` | Imperfection sensitivity of the sized 6V dome at hub-scale amplitudes (3 to 60 mm), mode-affine and with the crown hub alone |
| `mode_shapes.py` | Where the first LC2 buckling mode of the sized 4V, 6V and 8V domes sits, and how much of the crown rise each amplitude removes |
| `ring_and_rise.py` | Crown rise above its five neighbours for 2V to 8V hemispheres against span/1000, and the base-ring force under +30 °C against the horizontal support stiffness |

Run from the repository root, for example:

```bash
python examples/report_supplementary/hub_tolerance.py
```
