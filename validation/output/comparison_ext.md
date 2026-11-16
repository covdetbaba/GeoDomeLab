# GeoDomeLab extended validation: hand calculations F1-F10 vs the program

422 of 422 comparisons pass. 290 of 290 hand reference values are compared (every one of them). The hand file's own 120 internal checks: 120 within tolerance.

Produced by `validation/run_validation_ext.py`. The reference values come from `validation/hand_calc_ext.py`, which does not import the geodome package; they were computed in this run before geodome was imported (geodome loaded before the hand values: False). The program values come from the normal GeoDomeLab API (column `source`).

## Comparison rule

For a non-zero reference d = |program - hand| / |hand|, the plain relative difference. Where the exact answer is zero (a stress-free spoke, a zero reaction, a symmetric displacement) a relative difference is not defined, so d = |program| / s with s the largest |hand value| of the same unit in the same sub-case (column `scale` of the CSV). A row passes when d is within its tolerance:

- linear statics, linear buckling, geometry: d <= 1e-9
- natural frequencies and modal quantities: d <= 1e-9
- nonlinear path (limit point, snap-through target): d <= 1e-6
- counts: exact
- published 5-decimal chord factors: |p - h| <= 5e-6
- program identity (residual, energy mismatch, load balance): value <= 1e-10

## Summary by case

| Case | Structure | Rows | Pass | Worst d (linear) | Worst d (modal) | Worst d (nonlinear) | Worst identity |
|---|---|---|---|---|---|---|---|
| F1 | Heated fixed-fixed bars | 21 | 21 | 1.6e-15 | - | - | 7.4e-17 |
| F2 | Model A under a temperature change | 52 | 52 | 6.6e-16 | - | - | 2.2e-16 |
| F3 | Model A with a support movement | 66 | 66 | 1.2e-15 | - | - | 1.6e-16 |
| F4 | Three-bar truss by the force method | 78 | 78 | 7.2e-16 | - | - | 3.0e-16 |
| F5 | Bar and elastic spring support | 21 | 21 | 2.2e-16 | - | - | 1.7e-16 |
| F6 | von Mises truss snap-through | 32 | 32 | 1.2e-13 | - | 4.9e-10 | 0.0e+00 |
| F7 | Model A natural frequencies | 34 | 34 | - | 7.0e-16 | - | 4.4e-16 |
| F8 | 2V geodesic geometry and fabrication angles | 52 | 52 | 3.2e-15 | - | - | - |
| F9 | Wind on Model B | 52 | 52 | 2.6e-15 | - | - | 2.8e-16 |
| F10 | AISC 360-16 compression, B4.2 and E7 | 14 | 14 | 0.0e+00 | - | - | - |

## How each case is run in GeoDomeLab

- **F1** custom 2-node and 3-node trusses (`custom_geometry`, CUSTOM supports) with `LoadSpec.temperature_change = 30`, combination T = 1.0, `run_analysis`.
- **F2** `benchmarks.validation_model_a()` with a uniform temperature change through `run_analysis`; the one-spoke case uses `solve_static` with a per-member initial strain (the GUI load input heats every strut uniformly).
- **F3** Model A with `LoadSpec.settlements` (vertical settlement of hand node 2, radial spread of hand node 3), `run_analysis`.
- **F4** `benchmarks.three_bar_truss` (support C settles 5 mm): load case and settlement case through `run_analysis`, bar 2 heated alone and the combined case through `solve_static`. The flexibility coefficient f11 is obtained independently of the hand force method by a stiffness solve of the primary structure. The program's own force-method routine (`benchmarks.three_bar_force_method`, used by the built-in preset check) is compared too.
- **F5** CUSTOM supports with a spring stiffness entry (5e7 N/m) in series and in parallel with the bar, and a heated bar on a spring.
- **F6** `benchmarks.von_mises_truss()`: linear solve and `buckling.linear_buckling`; `nonlinear.run_nonlinear` (Crisfield arc length) traced through the limit point, the unstable branch and onto the inverted stable branch; every traced point is compared with the closed form P(w). The sample loads P(w) are the program's corotational internal force at the symmetric state; both limit points are also located as zeros of the program's consistent tangent; the snap-through target is the program's equilibrium at lambda* found by Newton iteration on its residual. `linear_buckling` returns n_free - 2 modes, one for this 2-DOF truss, so the sway root uses the same program matrices K and K_g.
- **F7** `modal.modal_analysis` on Model A, lumped mass, with and without the cladding dead load as mass.
- **F8** a full 2V Class I sphere built with the program's own `icosahedron`, `subdivide_class1` and `custom_geometry`, and the 2V hemisphere from `generate_dome`, both through `fabrication.fabrication_tables`; chord factors also against the printed values of the Domerama chord-factor chart.
- **F9** `benchmarks.validation_model_b()` with the default `LoadSpec` (wind on): `loads.velocity_pressure`, `wind_pressure_coefficients` and `wind_load` at every node, the resultant, and a `run_analysis` of 1.0W for the identities. This checks the implementation of the wind field; the 7-point Cp table itself is an illustrative simplification, not ASCE 7-16 Fig. 27.3-2.
- **F10** `buckling.compression_capacity` (inelastic E3 branch, E7 slender wall with a custom `Section`, the 0.45 E/Fy scope limit) and `buckling.member_capacities` of a one-bar model with `design_thickness_factor = 0.93` (AISC 360-16 B4.2).

## Failures

None.

## Full table

| Sub-case | Quantity | Hand key | Hand | GeoDomeLab | Unit | d | Tol | Pass | Program source |
|---|---|---|---|---|---|---|---|---|---|
| F1a | bar length L | F1.L | 3 | 3 | m | 0.0e+00 | 1e-09 | yes | custom_geometry |
| F1a | axial force N | F1.N | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | axial stress N/A | F1.sigma | -70200000 | -70200000 | Pa | 2.1e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Rx at node 1 | F1.R1_x | 49930.21169 | 49930.21169 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Rx at node 2 | F1.R2_x | -49930.21169 | -49930.21169 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Ry at node 1 | F1.R1_y | 24965.10584 | 24965.10584 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Ry at node 2 | F1.R2_y | -24965.10584 | -24965.10584 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Rz at node 1 | F1.R1_z | 49930.21169 | 49930.21169 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | reaction Rz at node 2 | F1.R2_z | -49930.21169 | -49930.21169 | N | 2.9e-16 | 1e-09 | yes | custom_geometry + CUSTOM pins, LoadSpec.temperature_change, run_analysis (T = 1.0) |
| F1a | free thermal elongation eps0 L | F1.free_elongation | 0.001053 | 0.001053 | m | 0.0e+00 | 1e-09 | yes | StaticResult.eps0 x L |
| F1a | temperature rise at the Euler load | F1.dT_euler | 84.63689564 | 84.63689564 | degC | 0.0e+00 | 1e-09 | yes | buckling.euler_load / (EA alpha) |
| F1a | F1a: equilibrium residual | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F1a | F1a: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F1a | F1a: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F1b | force in segment 1 | F1.stepped_N | -69599.82832 | -69599.82832 | N | 0.0e+00 | 1e-09 | yes | two sections, run_analysis (T = 1.0) |
| F1b | force in segment 2 | F1.stepped_N | -69599.82832 | -69599.82832 | N | 0.0e+00 | 1e-09 | yes | two sections, run_analysis (T = 1.0) |
| F1b | x displacement of the middle node | F1.stepped_u2 | 4.9635058e-05 | 4.9635058e-05 | m | 1.6e-15 | 1e-09 | yes | two sections, run_analysis (T = 1.0) |
| F1b | x reaction at node 1 | F1.stepped_R1x | 69599.82832 | 69599.82832 | N | 0.0e+00 | 1e-09 | yes | two sections, run_analysis (T = 1.0) |
| F1b | F1b: equilibrium residual | identity | 0 | 7.356496861e-17 | - | 7.4e-17 | 1e-10 | yes | solver.StaticResult |
| F1b | F1b: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F1b | F1b: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F2a | crown uz | F2.crown_uz | 0.00351 | 0.00351 | m | 3.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | crown ux | F2.crown_ux | 0 | -2.915534066e-19 | m | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | crown uy (symmetry, 0) | F2.crown_ux | 0 | 9.909771544e-20 | m | 2.8e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N spoke 1 | F2.N_spoke | 0 | -1.455191523e-11 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N spoke 2 | F2.N_spoke | 0 | 1.455191523e-11 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N spoke 3 | F2.N_spoke | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N spoke 4 | F2.N_spoke | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N spoke 5 | F2.N_spoke | 0 | 1.455191523e-11 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N ring member 6 | F2.N_ring | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N ring member 7 | F2.N_ring | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N ring member 8 | F2.N_ring | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N ring member 9 | F2.N_ring | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | N ring member 10 | F2.N_ring | -74895.31753 | -74895.31753 | N | 1.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | radial reaction at hand node 2 | F2.R_radial | -88044.72622 | -88044.72622 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | vertical reaction at hand node 2 | F2.R_z | 0 | 7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rx at hand node 2 | F2.R2_x | -88044.72622 | -88044.72622 | N | 3.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Ry at hand node 2 | F2.R2_y | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rz at hand node 2 | F2.R2_z | 0 | 7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | radial reaction at hand node 3 | F2.R_radial | -88044.72622 | -88044.72622 | N | 6.6e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | vertical reaction at hand node 3 | F2.R_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rx at hand node 3 | F2.R3_x | -27207.31667 | -27207.31667 | N | 1.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Ry at hand node 3 | F2.R3_y | -83735.5106 | -83735.5106 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rz at hand node 3 | F2.R3_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | radial reaction at hand node 4 | F2.R_radial | -88044.72622 | -88044.72622 | N | 3.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | vertical reaction at hand node 4 | F2.R_z | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rx at hand node 4 | F2.R4_x | 71229.67978 | 71229.67978 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Ry at hand node 4 | F2.R4_y | -51751.39162 | -51751.39162 | N | 5.6e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rz at hand node 4 | F2.R4_z | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | radial reaction at hand node 5 | F2.R_radial | -88044.72622 | -88044.72622 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | vertical reaction at hand node 5 | F2.R_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rx at hand node 5 | F2.R5_x | 71229.67978 | 71229.67978 | N | 2.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Ry at hand node 5 | F2.R5_y | 51751.39162 | 51751.39162 | N | 4.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rz at hand node 5 | F2.R5_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | radial reaction at hand node 6 | F2.R_radial | -88044.72622 | -88044.72622 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | vertical reaction at hand node 6 | F2.R_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rx at hand node 6 | F2.R6_x | -27207.31667 | -27207.31667 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Ry at hand node 6 | F2.R6_y | 83735.5106 | 83735.5106 | N | 3.5e-16 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | Rz at hand node 6 | F2.R6_z | 0 | -7.275957614e-12 | N | 8.3e-17 | 1e-09 | yes | validation_model_a, LoadSpec.temperature_change, run_analysis |
| F2a | F2a: equilibrium residual | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F2a | F2a: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F2a | F2a: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 1.652786697e-16 | - | 1.7e-16 | 1e-10 | yes | solver.StaticResult |
| F2b | crown ux | F2.oneSpoke_crown_ux | -0.0008677197202 | -0.0008677197202 | m | 1.2e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | crown uy | F2.oneSpoke_crown_uy | 0 | -2.201319473e-20 | m | 2.5e-17 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | crown uz | F2.oneSpoke_crown_uz | 0.000702 | 0.000702 | m | 6.2e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | N spoke 1 | F2.oneSpoke_N1 | -29958.12701 | -29958.12701 | N | 2.4e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | N spoke 2 | F2.oneSpoke_N2 | 24236.63387 | 24236.63387 | N | 0.0e+00 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | N spoke 3 | F2.oneSpoke_N3 | -9257.570367 | -9257.570367 | N | 2.0e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | N spoke 4 | F2.oneSpoke_N4 | -9257.570367 | -9257.570367 | N | 2.0e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | N spoke 5 | F2.oneSpoke_N5 | 24236.63387 | 24236.63387 | N | 1.5e-16 | 1e-09 | yes | solver.solve_static, eps0 on spoke 1 only |
| F2b | F2b: equilibrium residual | identity | 0 | 1.906648933e-33 | - | 1.9e-33 | 1e-10 | yes | solver.StaticResult |
| F2b | F2b: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F2b | F2b: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 2.172303471e-16 | - | 2.2e-16 | 1e-10 | yes | solver.StaticResult |
| F3a | crown ux | F3.vert_crown_ux | 0.001236067977 | 0.001236067977 | m | 1.8e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | crown uy | F3.vert_crown_uy | 0 | 3.135782725e-20 | m | 2.5e-17 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | crown uz | F3.vert_crown_uz | -0.001 | -0.001 | m | 2.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 1 | F3.vert_N1 | 42675.39461 | 42675.39461 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 2 | F3.vert_N2 | -34525.11948 | -34525.11948 | N | 4.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 3 | F3.vert_N3 | 13187.42218 | 13187.42218 | N | 5.5e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 4 | F3.vert_N4 | 13187.42218 | 13187.42218 | N | 5.5e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 5 | F3.vert_N5 | -34525.11948 | -34525.11948 | N | 2.1e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 6 | F3.vert_N6 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 7 | F3.vert_N7 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 8 | F3.vert_N8 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 9 | F3.vert_N9 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | N member 10 | F3.vert_N10 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rx at hand node 2 | F3.vert_R2_x | 36301.85892 | 36301.85892 | N | 2.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Ry at hand node 2 | F3.vert_R2_y | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rz at hand node 2 | F3.vert_R2_z | -22435.78267 | -22435.78267 | N | 1.6e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rx at hand node 3 | F3.vert_R3_x | -9075.46473 | -9075.46473 | N | 4.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Ry at hand node 3 | F3.vert_R3_y | -27931.40839 | -27931.40839 | N | 3.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rz at hand node 3 | F3.vert_R3_z | 18150.92946 | 18150.92946 | N | 6.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rx at hand node 4 | F3.vert_R4_x | -9075.46473 | -9075.46473 | N | 2.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Ry at hand node 4 | F3.vert_R4_y | 6593.711088 | 6593.711088 | N | 8.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rz at hand node 4 | F3.vert_R4_z | -6933.038126 | -6933.038126 | N | 5.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rx at hand node 5 | F3.vert_R5_x | -9075.46473 | -9075.46473 | N | 8.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Ry at hand node 5 | F3.vert_R5_y | -6593.711088 | -6593.711088 | N | 1.2e-15 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rz at hand node 5 | F3.vert_R5_z | -6933.038126 | -6933.038126 | N | 9.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rx at hand node 6 | F3.vert_R6_x | -9075.46473 | -9075.46473 | N | 4.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Ry at hand node 6 | F3.vert_R6_y | 27931.40839 | 27931.40839 | N | 1.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | Rz at hand node 6 | F3.vert_R6_z | 18150.92946 | 18150.92946 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3a | strain energy sum N^2 L/(2EA) | F3.vert_energy | 56.08945667 | 56.08945667 | J | 2.5e-16 | 1e-09 | yes | StaticResult.strain_energy_int |
| F3a | work 1/2 R . us (Clapeyron) | F3.vert_energy | 56.08945667 | 56.08945667 | J | 0.0e+00 | 1e-09 | yes | StaticResult.strain_energy_ext |
| F3a | F3a: equilibrium residual | identity | 0 | 1.078307462e-16 | - | 1.1e-16 | 1e-10 | yes | solver.StaticResult |
| F3a | F3a: energy identity mismatch | identity | 0 | 1.013442138e-16 | - | 1.0e-16 | 1e-10 | yes | solver.StaticResult |
| F3a | F3a: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 9.719732111e-17 | - | 9.7e-17 | 1e-10 | yes | solver.StaticResult |
| F3b | crown ux | F3.spread_crown_ux | 0.0006180339887 | 0.0006180339887 | m | 5.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | crown uy | F3.spread_crown_uy | 0.001902113033 | 0.001902113033 | m | 2.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | crown uz | F3.spread_crown_uz | -0.001618033989 | -0.001618033989 | m | 1.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 1 | F3.spread_N1 | -55862.81678 | -55862.81678 | N | 3.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 2 | F3.spread_N2 | 69050.23896 | 69050.23896 | N | 8.4e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 3 | F3.spread_N3 | -55862.81678 | -55862.81678 | N | 1.3e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 4 | F3.spread_N4 | 21337.6973 | 21337.6973 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 5 | F3.spread_N5 | 21337.6973 | 21337.6973 | N | 1.7e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 6 | F3.spread_N6 | 119281.3541 | 119281.3541 | N | 2.4e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 7 | F3.spread_N7 | 119281.3541 | 119281.3541 | N | 2.4e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 8 | F3.spread_N8 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 9 | F3.spread_N9 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | N member 10 | F3.spread_N10 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rx at hand node 2 | F3.spread_R2_x | 22592.07058 | 22592.07058 | N | 3.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Ry at hand node 2 | F3.spread_R2_y | -96500.64261 | -96500.64261 | N | 3.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rz at hand node 2 | F3.spread_R2_z | 29368.82079 | 29368.82079 | N | 1.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rx at hand node 3 | F3.spread_R3_x | 61482.41775 | 61482.41775 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Ry at hand node 3 | F3.spread_R3_y | 189223.4249 | 189223.4249 | N | 3.1e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rz at hand node 3 | F3.spread_R3_z | -36301.85892 | -36301.85892 | N | 6.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rx at hand node 4 | F3.spread_R4_x | -74999.0236 | -74999.0236 | N | 0.0e+00 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Ry at hand node 4 | F3.spread_R4_y | -64791.37393 | -64791.37393 | N | 1.1e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rz at hand node 4 | F3.spread_R4_z | 29368.82079 | 29368.82079 | N | 2.5e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rx at hand node 5 | F3.spread_R5_x | -14684.4104 | -14684.4104 | N | 5.0e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Ry at hand node 5 | F3.spread_R5_y | -10668.84865 | -10668.84865 | N | 3.4e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rz at hand node 5 | F3.spread_R5_z | -11217.89133 | -11217.89133 | N | 4.9e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rx at hand node 6 | F3.spread_R6_x | 5608.945667 | 5608.945667 | N | 8.1e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Ry at hand node 6 | F3.spread_R6_y | -17262.55974 | -17262.55974 | N | 4.2e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | Rz at hand node 6 | F3.spread_R6_z | -11217.89133 | -11217.89133 | N | 6.5e-16 | 1e-09 | yes | validation_model_a, LoadSpec.settlements, run_analysis |
| F3b | strain energy sum N^2 L/(2EA) | F3.spread_energy | 497.4032081 | 497.4032081 | J | 4.6e-16 | 1e-09 | yes | StaticResult.strain_energy_int |
| F3b | work 1/2 R . us (Clapeyron) | F3.spread_energy | 497.4032081 | 497.4032081 | J | 2.3e-16 | 1e-09 | yes | StaticResult.strain_energy_ext |
| F3b | F3b: equilibrium residual | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F3b | F3b: energy identity mismatch | identity | 0 | 1.584111618e-16 | - | 1.6e-16 | 1e-10 | yes | solver.StaticResult |
| F3b | F3b: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 8.592786905e-17 | - | 8.6e-17 | 1e-10 | yes | solver.StaticResult |
| F4 load | N bar 1 | F4.N1_load | -6956.712954 | -6956.712954 | N | 3.9e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | N bar 2 | F4.N2_load | -35260.60186 | -35260.60186 | N | 2.1e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | N bar 3 | F4.N3_load | -29083.82568 | -29083.82568 | N | 2.5e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | ux of joint D | F4.ux_load | 0.0002340532636 | 0.0002340532636 | m | 3.5e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | uz of joint D | F4.uz_load | -0.0005837620405 | -0.0005837620405 | m | 5.6e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rx at support A | F4.RAx_load | 5565.370364 | 5565.370364 | N | 6.5e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rz at support A | F4.RAz_load | 4174.027773 | 4174.027773 | N | 4.4e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rx at support B | F4.RBx_load | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rz at support B | F4.RBz_load | 35260.60186 | 35260.60186 | N | 2.1e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rx at support C | F4.RCx_load | -20565.37036 | -20565.37036 | N | 1.8e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | Rz at support C | F4.RCz_load | 20565.37036 | 20565.37036 | N | 1.8e-16 | 1e-09 | yes | three_bar_truss, run_analysis 'P' |
| F4 load | F4 load: equilibrium residual | identity | 0 | 1.176452618e-16 | - | 1.2e-16 | 1e-10 | yes | solver.StaticResult |
| F4 load | F4 load: energy identity mismatch | identity | 0 | 1.843816489e-16 | - | 1.8e-16 | 1e-10 | yes | solver.StaticResult |
| F4 load | F4 load: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 2.37121606e-16 | - | 2.4e-16 | 1e-10 | yes | solver.StaticResult |
| F4 settlement | N bar 1 | F4.N1_settlement | 49180.88585 | 49180.88585 | N | 4.4e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | N bar 2 | F4.N2_settlement | -68853.2402 | -68853.2402 | N | 0.0e+00 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | N bar 3 | F4.N3_settlement | 55641.82063 | 55641.82063 | N | 2.6e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | ux of joint D | F4.ux_settlement | 0.002295483985 | 0.002295483985 | m | 1.9e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | uz of joint D | F4.uz_settlement | -0.001139909867 | -0.001139909867 | m | 1.9e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rx at support A | F4.RAx_settlement | -39344.70868 | -39344.70868 | N | 3.7e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rz at support A | F4.RAz_settlement | -29508.53151 | -29508.53151 | N | 4.9e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rx at support B | F4.RBx_settlement | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rz at support B | F4.RBz_settlement | 68853.2402 | 68853.2402 | N | 0.0e+00 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rx at support C | F4.RCx_settlement | 39344.70868 | 39344.70868 | N | 3.7e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | Rz at support C | F4.RCz_settlement | -39344.70868 | -39344.70868 | N | 3.7e-16 | 1e-09 | yes | three_bar_truss (dT = 0), run_analysis 'T' |
| F4 settlement | F4 settlement: equilibrium residual | identity | 0 | 7.68932273e-17 | - | 7.7e-17 | 1e-10 | yes | solver.StaticResult |
| F4 settlement | F4 settlement: energy identity mismatch | identity | 0 | 1.808376589e-16 | - | 1.8e-16 | 1e-10 | yes | solver.StaticResult |
| F4 settlement | F4 settlement: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 2.962343021e-16 | - | 3.0e-16 | 1e-10 | yes | solver.StaticResult |
| F4 temperature | N bar 1 | F4.N1_temperature | 18125.61548 | 18125.61548 | N | 2.0e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | N bar 2 | F4.N2_temperature | -25375.86167 | -25375.86167 | N | 1.4e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | N bar 3 | F4.N3_temperature | 20506.79299 | 20506.79299 | N | 3.5e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | ux of joint D | F4.ux_temperature | 5.625062267e-05 | 5.625062267e-05 | m | 7.2e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | uz of joint D | F4.uz_temperature | 0.0006328862184 | 0.0006328862184 | m | 3.4e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rx at support A | F4.RAx_temperature | -14500.49239 | -14500.49239 | N | 3.8e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rz at support A | F4.RAz_temperature | -10875.36929 | -10875.36929 | N | 1.7e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rx at support B | F4.RBx_temperature | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rz at support B | F4.RBz_temperature | 25375.86167 | 25375.86167 | N | 1.4e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rx at support C | F4.RCx_temperature | 14500.49239 | 14500.49239 | N | 2.5e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | Rz at support C | F4.RCz_temperature | -14500.49239 | -14500.49239 | N | 2.5e-16 | 1e-09 | yes | solver.solve_static, eps0 on bar 2 only |
| F4 temperature | F4 temperature: equilibrium residual | identity | 0 | 8.088964356e-17 | - | 8.1e-17 | 1e-10 | yes | solver.StaticResult |
| F4 temperature | F4 temperature: energy identity mismatch | identity | 0 | 1.060912562e-16 | - | 1.1e-16 | 1e-10 | yes | solver.StaticResult |
| F4 temperature | F4 temperature: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 1.541183952e-16 | - | 1.5e-16 | 1e-10 | yes | solver.StaticResult |
| F4 combined | N bar 1 | F4.N1_combined | 60349.78838 | 60349.78838 | N | 0.0e+00 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | N bar 2 | F4.N2_combined | -129489.7037 | -129489.7037 | N | 1.1e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | N bar 3 | F4.N3_combined | 47064.78794 | 47064.78794 | N | 3.1e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | ux of joint D | F4.ux_combined | 0.002585787871 | 0.002585787871 | m | 0.0e+00 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | uz of joint D | F4.uz_combined | -0.00109078569 | -0.00109078569 | m | 0.0e+00 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rx at support A | F4.RAx_combined | -48279.83071 | -48279.83071 | N | 1.5e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rz at support A | F4.RAz_combined | -36209.87303 | -36209.87303 | N | 2.0e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rx at support B | F4.RBx_combined | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rz at support B | F4.RBz_combined | 129489.7037 | 129489.7037 | N | 1.1e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rx at support C | F4.RCx_combined | 33279.83071 | 33279.83071 | N | 2.2e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | Rz at support C | F4.RCz_combined | -33279.83071 | -33279.83071 | N | 2.2e-16 | 1e-09 | yes | solver.solve_static, load + eps0 (bar 2) + settlement |
| F4 combined | F4 combined: equilibrium residual | identity | 0 | 5.786805085e-17 | - | 5.8e-17 | 1e-10 | yes | solver.StaticResult |
| F4 combined | F4 combined: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F4 combined | F4 combined: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 1.081640729e-16 | - | 1.1e-16 | 1e-10 | yes | solver.StaticResult |
| F4 flexibility | f11 = L2/EA2 + c2 . K13^-1 c2 (primary structure) | F4.f11 | 4.149612784e-08 | 4.149612784e-08 | m/N | 1.6e-16 | 1e-09 | yes | solver.solve_static on bars 1 and 3 |
| F4 force-method load | N bar 1 | F4.N1_load | -6956.712954 | -6956.712954 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method load | N bar 2 | F4.N2_load | -35260.60186 | -35260.60186 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method load | N bar 3 | F4.N3_load | -29083.82568 | -29083.82568 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method load | ux of joint D | F4.ux_load | 0.0002340532636 | 0.0002340532636 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method load | uz of joint D | F4.uz_load | -0.0005837620405 | -0.0005837620405 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method temperature | N bar 1 | F4.N1_temperature | 18125.61548 | 18125.61548 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method temperature | N bar 2 | F4.N2_temperature | -25375.86167 | -25375.86167 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method temperature | N bar 3 | F4.N3_temperature | 20506.79299 | 20506.79299 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method temperature | ux of joint D | F4.ux_temperature | 5.625062267e-05 | 5.625062267e-05 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method temperature | uz of joint D | F4.uz_temperature | 0.0006328862184 | 0.0006328862184 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method settlement | N bar 1 | F4.N1_settlement | 49180.88585 | 49180.88585 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method settlement | N bar 2 | F4.N2_settlement | -68853.2402 | -68853.2402 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method settlement | N bar 3 | F4.N3_settlement | 55641.82063 | 55641.82063 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method settlement | ux of joint D | F4.ux_settlement | 0.002295483985 | 0.002295483985 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method settlement | uz of joint D | F4.uz_settlement | -0.001139909867 | -0.001139909867 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method combined | N bar 1 | F4.N1_combined | 60349.78838 | 60349.78838 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method combined | N bar 2 | F4.N2_combined | -129489.7037 | -129489.7037 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method combined | N bar 3 | F4.N3_combined | 47064.78794 | 47064.78794 | N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method combined | ux of joint D | F4.ux_combined | 0.002585787871 | 0.002585787871 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method combined | uz of joint D | F4.uz_combined | -0.00109078569 | -0.00109078569 | m | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F4 force-method load | f11 | F4.f11 | 4.149612784e-08 | 4.149612784e-08 | m/N | 0.0e+00 | 1e-09 | yes | benchmarks.three_bar_force_method |
| F5a | u1 (spring node) | F5.series_u1 | 0.001 | 0.001 | m | 2.2e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5a | u2 (loaded node) | F5.series_u2 | 0.001702981197 | 0.001702981197 | m | 1.3e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5a | bar force | F5.series_N | 50000 | 50000 | N | 1.5e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5a | spring reaction | F5.series_spring_reaction | -50000 | -50000 | N | 1.5e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5a | StaticResult.spring_forces | F5.series_spring_reaction | -50000 | -50000 | N | 1.5e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5a | F5a: equilibrium residual | identity | 0 | 1.455191523e-16 | - | 1.5e-16 | 1e-10 | yes | solver.StaticResult |
| F5a | F5a: energy identity mismatch | identity | 0 | 1.668938534e-16 | - | 1.7e-16 | 1e-10 | yes | solver.StaticResult |
| F5a | F5a: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 1.455191523e-16 | - | 1.5e-16 | 1e-10 | yes | solver.StaticResult |
| F5b | u2 | F5.parallel_u2 | 0.0004127944563 | 0.0004127944563 | m | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5b | bar force | F5.parallel_N | 29360.27719 | 29360.27719 | N | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5b | spring reaction | F5.parallel_spring_reaction | -20639.72281 | -20639.72281 | N | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5b | fixed reaction at node 1 | F5.parallel_R1x | -29360.27719 | -29360.27719 | N | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5b | F5b: equilibrium residual | identity | 0 | 7.275957614e-17 | - | 7.3e-17 | 1e-10 | yes | solver.StaticResult |
| F5b | F5b: energy identity mismatch | identity | 0 | 1.721299124e-16 | - | 1.7e-16 | 1e-10 | yes | solver.StaticResult |
| F5b | F5b: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F5c | u2 | F5.thermal_u2 | 0.0006183274376 | 0.0006183274376 | m | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5c | bar force | F5.thermal_N | -30916.37188 | -30916.37188 | N | 1.2e-16 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5c | spring reaction | F5.thermal_spring_reaction | -30916.37188 | -30916.37188 | N | 0.0e+00 | 1e-09 | yes | CUSTOM supports with a spring entry, run_analysis |
| F5c | F5c: equilibrium residual | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F5c | F5c: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F5c | F5c: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F6 linear | bar length | F6.L | 2.503123049 | 2.503123049 | m | 0.0e+00 | 1e-09 | yes | custom_geometry |
| F6 linear | axial rigidity EA | F6.EA | 213376973 | 213376973 | N | 0.0e+00 | 1e-09 | yes | DomeModel.EA |
| F6 linear | reference load (downward) | F6.P_ref | 10000 | 10000 | N | 0.0e+00 | 1e-09 | yes | LoadSpec.nodal_loads -> loads.combine |
| F6 linear | linear bar force under P_ref | F6.N_ref | -100124.922 | -100124.922 | N | 1.5e-16 | 1e-09 | yes | run_analysis |
| F6 linear | initial stiffness P_ref / w_lin | F6.K0 | 425158.6052 | 425158.6052 | N/m | 1.4e-16 | 1e-09 | yes | run_analysis |
| F6 linear | F6 linear: equilibrium residual | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F6 linear | F6 linear: energy identity mismatch | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F6 linear | F6 linear: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 0 | - | 0.0e+00 | 1e-10 | yes | solver.StaticResult |
| F6 linear | linear buckling factor | F6.lambda_cr | 5.327768772 | 5.327768772 | - | 1.7e-16 | 1e-09 | yes | buckling.linear_buckling |
| F6 linear | lambda_cr P_ref | F6.P_cr_linear | 53277.68772 | 53277.68772 | N | 2.7e-16 | 1e-09 | yes | buckling.linear_buckling |
| F6 linear | lowest root of det(K + lam Kg) | F6.lambda_cr | 5.327768772 | 5.327768772 | - | 1.7e-16 | 1e-09 | yes | solver.assemble_stiffness / assemble_geometric_stiffness |
| F6 linear | sway root of det(K + lam Kg) | F6.lambda_sway | 852443.0035 | 852443.0035 | - | 1.2e-13 | 1e-09 | yes | solver.assemble_stiffness / assemble_geometric_stiffness |
| F6 nonlinear | limit load lambda_limit P_ref | F6.P_limit | 10240.4996 | 10240.4996 | N | 3.6e-16 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | apex deflection at the limit | F6.w_limit | 0.05286125123 | 0.05286125126 | m | 4.9e-10 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | bar length at the limit | F6.l_limit | 2.501040583 | 2.501040583 | m | 3.0e-13 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | bar force at the limit | F6.N_limit | -177518.3623 | -177518.3624 | N | 3.6e-10 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | strain at the limit | F6.strain_limit | -0.0008319471394 | -0.0008319471397 | - | 3.6e-10 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | lambda_cr / lambda_limit | F6.ratio_linear_to_limit | 5.202645361 | 5.202645361 | - | 6.8e-16 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | lambda_cr reported by run_nonlinear | F6.lambda_cr | 5.327768772 | 5.327768772 | - | 1.7e-16 | 1e-09 | yes | run_nonlinear (linear_buckling) |
| F6 nonlinear | max \|P_program - P_exact(w)\| / P* over all 63 path points (w up to 2.84 h, 25 on the unstable branch) | path | 0 | 2.842032284e-15 | - | 2.8e-15 | 1e-06 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 nonlinear | path passes the snap-through target (max w / w_snap >= 1) | path | 1 | 1 | - | 0.0e+00 | 0e+00 | yes | run_nonlinear arc-length (62 steps, 0.09 s) |
| F6 tangent | first limit point w* | F6.w_limit | 0.05286125123 | 0.05286125123 | m | 1.6e-14 | 1e-06 | yes | zero of the program's tangent K_T,zz (TrussSystem.tangent) |
| F6 tangent | second limit point h + y* | F6.w_second_limit | 0.1971387488 | 0.1971387488 | m | 4.4e-15 | 1e-06 | yes | zero of the program's tangent K_T,zz (TrussSystem.tangent) |
| F6 path | P at w = 0.1 h | F6.P_at_w0.1h | 4545.497953 | 4545.497953 | N | 0.0e+00 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 0.25 h | F6.P_at_w0.25h | 8726.21398 | 8726.21398 | N | 4.2e-16 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 0.5 h | F6.P_at_w0.5h | 9978.654483 | 9978.654483 | N | 1.8e-16 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 1 h | F6.P_at_w1h | 0 | 0 | N | 0.0e+00 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 1.5 h | F6.P_at_w1.5h | -9978.654483 | -9978.654483 | N | 1.8e-16 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 2 h | F6.P_at_w2h | 0 | 0 | N | 0.0e+00 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 path | P at w = 2.5 h | F6.P_at_w2.5h | 49707.09372 | 49707.09372 | N | 0.0e+00 | 1e-06 | yes | TrussSystem.state (corotational element) |
| F6 snap | apex deflection after snap-through | F6.w_snap | 0.2693676147 | 0.2693676147 | m | 2.1e-16 | 1e-06 | yes | Newton on the program residual at lambda* from the traced path |
| F6 snap | bar force after snap-through | F6.N_snap | 88814.58642 | 88814.58642 | N | 2.3e-15 | 1e-06 | yes | Newton on the program residual at lambda* from the traced path |
| F7 sw | lumped crown mass | F7.m_crown_sw | 110.0755587 | 110.0755587 | kg | 3.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | free mass in x (= crown mass) | F7.m_crown_sw | 110.0755587 | 110.0755587 | kg | 3.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | f of the x mode | F7.f_xy_sw | 129.9834146 | 129.9834146 | Hz | 4.4e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | omega of the x mode | F7.omega_xy_sw | 816.7098805 | 816.7098805 | rad/s | 5.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | f of the y mode | F7.f_xy_sw | 129.9834146 | 129.9834146 | Hz | 6.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | omega of the y mode | F7.omega_xy_sw | 816.7098805 | 816.7098805 | rad/s | 7.0e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | f of the z mode | F7.f_z_sw | 113.6096702 | 113.6096702 | Hz | 6.3e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | omega of the z mode | F7.omega_z_sw | 713.8306102 | 713.8306102 | rad/s | 6.4e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | period of the z mode | F7.T_z_sw | 0.008802067629 | 0.008802067629 | s | 5.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | max \|phi\| of the x mode | F7.phi_sw | 0.09531352924 | 0.09531352924 | 1/sqrt(kg) | 2.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | effective mass of the x mode in x | F7.meff_x_sw | 110.0755587 | 110.0755587 | kg | 1.3e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | max \|phi\| of the y mode | F7.phi_sw | 0.09531352924 | 0.09531352924 | 1/sqrt(kg) | 2.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | effective mass of the y mode in y | F7.meff_y_sw | 110.0755587 | 110.0755587 | kg | 1.3e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | max \|phi\| of the z mode | F7.phi_sw | 0.09531352924 | 0.09531352924 | 1/sqrt(kg) | 1.5e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | effective mass of the z mode in z | F7.meff_z_sw | 110.0755587 | 110.0755587 | kg | 5.2e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=False) |
| F7 sw | M-orthonormality error max\|Phi^T M Phi - I\| | identity | 0 | 2.220446049e-16 | - | 2.2e-16 | 1e-10 | yes | modal.ModalResult |
| F7 sw | eigen residual (normwise backward error) | identity | 0 | 2.000302751e-16 | - | 2.0e-16 | 1e-10 | yes | modal.ModalResult |
| F7 dead | lumped crown mass | F7.m_crown_dead | 720.2827512 | 720.2827512 | kg | 1.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | free mass in x (= crown mass) | F7.m_crown_dead | 720.2827512 | 720.2827512 | kg | 1.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | f of the x mode | F7.f_xy_dead | 50.81382449 | 50.81382449 | Hz | 1.4e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | omega of the x mode | F7.omega_xy_dead | 319.2726754 | 319.2726754 | rad/s | 1.8e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | f of the y mode | F7.f_xy_dead | 50.81382449 | 50.81382449 | Hz | 2.8e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | omega of the y mode | F7.omega_xy_dead | 319.2726754 | 319.2726754 | rad/s | 1.8e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | f of the z mode | F7.f_z_dead | 44.41291113 | 44.41291113 | Hz | 4.8e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | omega of the z mode | F7.omega_z_dead | 279.0545507 | 279.0545507 | rad/s | 4.1e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | period of the z mode | F7.T_z_dead | 0.02251597507 | 0.02251597507 | s | 4.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | max \|phi\| of the x mode | F7.phi_dead | 0.03726048406 | 0.03726048406 | 1/sqrt(kg) | 0.0e+00 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | effective mass of the x mode in x | F7.meff_x_dead | 720.2827512 | 720.2827512 | kg | 6.3e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | max \|phi\| of the y mode | F7.phi_dead | 0.03726048406 | 0.03726048406 | 1/sqrt(kg) | 1.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | effective mass of the y mode in y | F7.meff_y_dead | 720.2827512 | 720.2827512 | kg | 1.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | max \|phi\| of the z mode | F7.phi_dead | 0.03726048406 | 0.03726048406 | 1/sqrt(kg) | 1.9e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | effective mass of the z mode in z | F7.meff_z_dead | 720.2827512 | 720.2827512 | kg | 1.6e-16 | 1e-09 | yes | modal.modal_analysis (lumped, include_dead=True) |
| F7 dead | M-orthonormality error max\|Phi^T M Phi - I\| | identity | 0 | 4.440892099e-16 | - | 4.4e-16 | 1e-10 | yes | modal.ModalResult |
| F7 dead | eigen residual (normwise backward error) | identity | 0 | 2.582002532e-16 | - | 2.6e-16 | 1e-10 | yes | modal.ModalResult |
| F8 sphere | nodes | F8.nodes | 42 | 42 | - | 0.0e+00 | 0e+00 | yes | icosahedron + subdivide_class1 + custom_geometry |
| F8 sphere | struts | F8.struts | 120 | 120 | - | 0.0e+00 | 0e+00 | yes | edges_from_faces |
| F8 sphere | panels | F8.panels | 80 | 80 | - | 0.0e+00 | 0e+00 | yes | subdivide_class1 |
| F8 sphere | chord factor A | F8.CF_A | 0.5465330578 | 0.5465330578 | - | 6.1e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | strut A length (R = 5 m) | F8.L_A_R5 | 2.732665289 | 2.732665289 | m | 1.6e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | axial angle A | F8.axial_A_deg | 74.14126279 | 74.14126279 | deg | 1.9e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | tangent (bevel) angle A | F8.bevel_A_deg | 15.85873721 | 15.85873721 | deg | 9.0e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | dihedral angle across A | F8.dihedral_A_deg | 157.5410761 | 157.5410761 | deg | 0.0e+00 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | bend angle across A | F8.bend_A_deg | 22.45892392 | 22.45892392 | deg | 0.0e+00 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | chord factor A vs Domerama chart | F8.CF_A | 0.54653 | 0.5465330578 | - | 3.1e-06 | 5e-06 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | chord factor B | F8.CF_B | 0.6180339887 | 0.6180339887 | - | 9.0e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | strut B length (R = 5 m) | F8.L_B_R5 | 3.090169944 | 3.090169944 | m | 5.7e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | axial angle B | F8.axial_B_deg | 72 | 72 | deg | 2.0e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | tangent (bevel) angle B | F8.bevel_B_deg | 18 | 18 | deg | 7.9e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | dihedral angle across B | F8.dihedral_B_deg | 161.9708979 | 161.9708979 | deg | 1.8e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | bend angle across B | F8.bend_B_deg | 18.02910211 | 18.02910211 | deg | 1.6e-15 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | chord factor B vs Domerama chart | F8.CF_B | 0.61803 | 0.6180339887 | - | 4.0e-06 | 5e-06 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | AAB apex angle (opposite B) | F8.face_AAB_apex_deg | 68.86197642 | 68.86197642 | deg | 6.2e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | AAB base angle 1 | F8.face_AAB_base_deg | 55.56901179 | 55.56901179 | deg | 1.3e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | AAB base angle 2 | F8.face_AAB_base_deg | 55.56901179 | 55.56901179 | deg | 1.3e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | AAB area / R^2 | F8.area_AAB_R1 | 0.139300145 | 0.139300145 | - | 6.0e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | BBB area / R^2 | F8.area_BBB_R1 | 0.1653961346 | 0.1653961346 | - | 1.0e-15 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | total panel area / R^2 | F8.area_total_R1 | 11.66593139 | 11.66593139 | - | 6.1e-16 | 1e-09 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | A struts | F8.struts_A | 60 | 60 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | B struts | F8.struts_B | 60 | 60 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | AAB panels | F8.panels_AAB | 60 | 60 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | BBB panels | F8.panels_BBB | 20 | 20 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | 5-valent hubs | F8.hubs_valence5 | 12 | 12 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 sphere | 6-valent hubs | F8.hubs_valence6 | 30 | 30 | - | 0.0e+00 | 0e+00 | yes | fabrication_tables (full 2V sphere) |
| F8 dome | chord factor A | F8.CF_A | 0.5465330578 | 0.5465330578 | - | 2.0e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | strut A length (R = 5 m) | F8.L_A_R5 | 2.732665289 | 2.732665289 | m | 1.6e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | axial angle A | F8.axial_A_deg | 74.14126279 | 74.14126279 | deg | 1.9e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | tangent (bevel) angle A | F8.bevel_A_deg | 15.85873721 | 15.85873721 | deg | 9.0e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | dihedral angle across A | F8.dihedral_A_deg | 157.5410761 | 157.5410761 | deg | 3.6e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | bend angle across A | F8.bend_A_deg | 22.45892392 | 22.45892392 | deg | 2.5e-15 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | chord factor A vs Domerama chart | F8.CF_A | 0.54653 | 0.5465330578 | - | 3.1e-06 | 5e-06 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | chord factor B | F8.CF_B | 0.6180339887 | 0.6180339887 | - | 5.4e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | strut B length (R = 5 m) | F8.L_B_R5 | 3.090169944 | 3.090169944 | m | 5.7e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | axial angle B | F8.axial_B_deg | 72 | 72 | deg | 3.9e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | tangent (bevel) angle B | F8.bevel_B_deg | 18 | 18 | deg | 1.6e-15 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | dihedral angle across B | F8.dihedral_B_deg | 161.9708979 | 161.9708979 | deg | 3.5e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | bend angle across B | F8.bend_B_deg | 18.02910211 | 18.02910211 | deg | 3.2e-15 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | chord factor B vs Domerama chart | F8.CF_B | 0.61803 | 0.6180339887 | - | 4.0e-06 | 5e-06 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | AAB apex angle (opposite B) | F8.face_AAB_apex_deg | 68.86197642 | 68.86197642 | deg | 2.1e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | AAB base angle 1 | F8.face_AAB_base_deg | 55.56901179 | 55.56901179 | deg | 2.6e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | AAB base angle 2 | F8.face_AAB_base_deg | 55.56901179 | 55.56901179 | deg | 2.6e-16 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | AAB area / R^2 | F8.area_AAB_R1 | 0.139300145 | 0.139300145 | - | 0.0e+00 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 dome | BBB area / R^2 | F8.area_BBB_R1 | 0.1653961346 | 0.1653961346 | - | 1.0e-15 | 1e-09 | yes | generate_dome(2V, 0.5) -> fabrication_tables |
| F8 sphere | axial angle A, program formula column | F8.axial_A_deg | 74.14126279 | 74.14126279 | deg | 1.9e-16 | 1e-09 | yes | fabrication_tables axial_formula_deg |
| F8 sphere | axial angle B, program formula column | F8.axial_B_deg | 72 | 72 | deg | 2.0e-16 | 1e-09 | yes | fabrication_tables axial_formula_deg |
| F8 sphere | dihedral A, program formula column | F8.dihedral_A_deg | 157.5410761 | 157.5410761 | deg | 1.8e-16 | 1e-09 | yes | fabrication_tables dihedral_formula_deg |
| F8 sphere | dihedral B, program formula column | F8.dihedral_B_deg | 161.9708979 | 161.9708979 | deg | 1.8e-16 | 1e-09 | yes | fabrication_tables dihedral_formula_deg |
| F9 wind Model B | velocity pressure qz | F9.qz | 1061.332875 | 1061.332875 | Pa | 0.0e+00 | 1e-09 | yes | loads.velocity_pressure (LoadSpec defaults) |
| F9 wind Model B | qz G | F9.qG | 902.1329437 | 902.1329437 | Pa | 0.0e+00 | 1e-09 | yes | loads.velocity_pressure x LoadSpec.wind_G |
| F9 wind Model B | Cp at hand node 1 | F9.Cp_1 | -1 | -1 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 1 | F9.Fx_1 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 1 | F9.Fy_1 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 1 | F9.Fz_1 | 17994.81084 | 17994.81084 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 2 | F9.Cp_2 | -0.4885501706 | -0.4885501706 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 2 | F9.Fx_2 | 7863.238501 | 7863.238501 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 2 | F9.Fy_2 | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 2 | F9.Fz_2 | 3931.61925 | 3931.61925 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 3 | F9.Cp_3 | -0.8930329524 | -0.8930329524 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 3 | F9.Fx_3 | 4441.627443 | 4441.627443 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 3 | F9.Fy_3 | 13669.92366 | 13669.92366 | N | 2.7e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 3 | F9.Fz_3 | 7186.704168 | 7186.704168 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 4 | F9.Cp_4 | -0.1639180562 | -0.1639180562 | - | 1.0e-15 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 4 | F9.Fx_4 | -2134.404458 | -2134.404458 | N | 1.1e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 4 | F9.Fy_4 | 1550.735611 | 1550.735611 | N | 1.2e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 4 | F9.Fz_4 | 1319.134501 | 1319.134501 | N | 1.0e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 5 | F9.Cp_5 | -0.1639180562 | -0.1639180562 | - | 2.4e-15 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 5 | F9.Fx_5 | -2134.404458 | -2134.404458 | N | 2.3e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 5 | F9.Fy_5 | -1550.735611 | -1550.735611 | N | 2.3e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 5 | F9.Fz_5 | 1319.134501 | 1319.134501 | N | 2.6e-15 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 6 | F9.Cp_6 | -0.8930329524 | -0.8930329524 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 6 | F9.Fx_6 | 4441.627443 | 4441.627443 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 6 | F9.Fy_6 | -13669.92366 | -13669.92366 | N | 1.3e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 6 | F9.Fz_6 | 7186.704168 | 7186.704168 | N | 2.5e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 7 | F9.Cp_7 | -0.6364692711 | -0.6364692711 | - | 1.7e-16 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 7 | F9.Fx_7 | 4972.543773 | 4972.543773 | N | 1.8e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 7 | F9.Fy_7 | 3612.764523 | 3612.764523 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 7 | F9.Fz_7 | -3073.201062 | -3073.201062 | N | 3.0e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 8 | F9.Cp_8 | -0.7860659049 | -0.7860659049 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 8 | F9.Fx_8 | -2345.766896 | -2345.766896 | N | 1.9e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 8 | F9.Fy_8 | 7219.528158 | 7219.528158 | N | 2.5e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 8 | F9.Fz_8 | -3795.530567 | -3795.530567 | N | 1.2e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 9 | F9.Cp_9 | 0.2686989765 | 0.2686989765 | - | 0.0e+00 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 9 | F9.Fx_9 | 2594.833772 | 2594.833772 | N | 1.8e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 9 | F9.Fy_9 | 0 | -3.177754874e-13 | N | 1.2e-17 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 9 | F9.Fz_9 | 1297.416886 | 1297.416886 | N | 1.8e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 10 | F9.Cp_10 | -0.7860659049 | -0.7860659049 | - | 2.8e-16 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 10 | F9.Fx_10 | -2345.766896 | -2345.766896 | N | 1.9e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 10 | F9.Fy_10 | -7219.528158 | -7219.528158 | N | 1.3e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 10 | F9.Fz_10 | -3795.530567 | -3795.530567 | N | 3.6e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Cp at hand node 11 | F9.Cp_11 | -0.6364692711 | -0.6364692711 | - | 1.7e-16 | 1e-09 | yes | loads.wind_pressure_coefficients |
| F9 wind Model B | Fx at hand node 11 | F9.Fx_11 | 4972.543773 | 4972.543773 | N | 0.0e+00 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fy at hand node 11 | F9.Fy_11 | -3612.764523 | -3612.764523 | N | 2.5e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | Fz at hand node 11 | F9.Fz_11 | -3073.201062 | -3073.201062 | N | 3.0e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | wind resultant x | F9.Rx | 20326.072 | 20326.072 | N | 5.4e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | wind resultant y | F9.Ry | 0 | -3.637978807e-12 | N | 1.4e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | wind resultant z | F9.Rz | 26498.06105 | 26498.06105 | N | 1.4e-16 | 1e-09 | yes | loads.wind_load (LoadSpec defaults) |
| F9 wind Model B | 1.0W solve: equilibrium residual | identity | 0 | 2.832544991e-16 | - | 2.8e-16 | 1e-10 | yes | solver.StaticResult |
| F9 wind Model B | 1.0W solve: energy identity mismatch | identity | 0 | 1.962862267e-16 | - | 2.0e-16 | 1e-10 | yes | solver.StaticResult |
| F9 wind Model B | 1.0W solve: \|sum F + sum R\| / max(\|F\|, \|R\|) | identity | 0 | 2.777773583e-16 | - | 2.8e-16 | 1e-10 | yes | solver.StaticResult |
| F10 AISC | phi Pn, inelastic | F10.phiPn_inel | 239578.0961 | 239578.0961 | N | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity |
| F10 AISC | Fcr, inelastic | F10.Fcr_inel | 249509476.9 | 249509476.9 | Pa | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity |
| F10 AISC | Ae/Ag, inelastic | F10.Q_inel | 1 | 1 | - | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity |
| F10 AISC | KL/r, inelastic | F10.KLr_inel | 66.5557865 | 66.5557865 | - | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity |
| F10 AISC | phi Pn, B4.2 wall | F10.phiPn_b42 | 223995.9127 | 223995.9127 | N | 0.0e+00 | 1e-09 | yes | buckling.member_capacities, design_thickness_factor 0.93 |
| F10 AISC | Fcr, B4.2 wall | F10.Fcr_b42 | 250015616.3 | 250015616.3 | Pa | 0.0e+00 | 1e-09 | yes | buckling.member_capacities, design_thickness_factor 0.93 |
| F10 AISC | Ae/Ag, B4.2 wall | F10.Q_b42 | 1 | 1 | - | 0.0e+00 | 1e-09 | yes | buckling.member_capacities, design_thickness_factor 0.93 |
| F10 AISC | KL/r, B4.2 wall | F10.KLr_b42 | 66.34735148 | 66.34735148 | - | 0.0e+00 | 1e-09 | yes | buckling.member_capacities, design_thickness_factor 0.93 |
| F10 AISC | phi_t Pn, B4.2 wall | F10.phiPn_t_b42 | 309095.0519 | 309095.0519 | N | 0.0e+00 | 1e-09 | yes | buckling.member_capacities, design_thickness_factor 0.93 |
| F10 AISC | phi Pn, E7 slender wall | F10.phiPn_e7 | 1112934.907 | 1112934.907 | N | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity (custom Section) |
| F10 AISC | Fcr, E7 slender wall | F10.Fcr_e7 | 327695234.8 | 327695234.8 | Pa | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity (custom Section) |
| F10 AISC | Ae/Ag, E7 (Eq. E7-7) | F10.Q_e7 | 0.9387134158 | 0.9387134158 | - | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity (custom Section) |
| F10 AISC | KL/r, E7 slender wall | F10.KLr_e7 | 26.52271999 | 26.52271999 | - | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity (custom Section) |
| F10 AISC | phi Pn, D/t above 0.45 E/Fy | F10.phiPn_scope | 0 | 0 | N | 0.0e+00 | 1e-09 | yes | buckling.compression_capacity (custom Section) |

Run time 0.1 s.
