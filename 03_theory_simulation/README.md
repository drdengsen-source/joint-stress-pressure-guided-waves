# Theory and simulation (Sections 2 and 3.2, Appendix A, Supplementary S1–S3)

Two independent routes to the delay sensitivities, plus the transient EMAT model:

| Folder | Method | Software |
|---|---|---|
| `1_theory/` (theory) | prestressed fluid-filled-pipe SAFE model, first-order perturbation (Eqs. (4)–(7)), closed-form limits (Appendix A) | MATLAB only |
| `2_simulation/` T4_* | Floquet–Bloch unit-cell FEM under load (independent validation, Supplementary S1.4) | COMSOL 6.4 + LiveLink for MATLAB |
| `2_simulation/` E3_* | 2D-axisymmetric transient EMAT model, two probes 100 mm apart (Section 3.2, Supplementary S3.2) | COMSOL 6.4 (model file); MATLAB for post-processing |

Keep the folder structure: the scripts find each other through relative paths (`2_simulation` → `../1_theory`).

## 1_theory — SAFE model (MATLAB, no toolboxes)

```
cd 1_theory
T1_thin_shell        % closed-form limits                       (~1 s)
T3_sensitivity       % dispersion and sensitivities, 0–800 kHz -> T3_sweep.mat (a few minutes)
T3_check_perturb     % first-order perturbation vs direct re-solution
```

| File | Content | Manuscript |
|---|---|---|
| `params_pipe.m` | tube, oil and TOE parameters; oil c_f(p, T), B/A, β_f, C_p | Table 1, Appendix A.2 |
| `murnaghan_A.m`, `tensor_to_D5.m` | tangent stiffness A(σ⁰) with Murnaghan constants, reduced to the axisymmetric 5×5 form | Eq. (4) |
| `safe_pipe_setup.m`, `safe_pipe_solve.m` | SAFE discretization of the liquid-filled pipe (quadratic elements, wall 8 / liquid 40), modes, power flow, c_g, liquid energy fraction η | Eq. (5) |
| `safe_perturb.m` | first-order wavenumber perturbation | Eqs. (6), (7), Appendix A.1 |
| `s0_reduced_modulus.m`, `lamb_S0_phase.m` | S0 plate-wave limit | Appendix A.3 |
| `T1_thin_shell.m` | K_z, K_θ, K_r, K_p^ae; S11 = 1.16, wall-only S12 = −0.74 ns/MPa; Korteweg Ω = 0.320, c_K = 1262 m/s, S22 = −176 ns/MPa before the hoop-stiffness term of the wall (−168 ns/MPa with it, Supplementary S1.3), S21 ≈ −0.24 ns/MPa; flat-plate slope −0.64 | Appendix A.2–A.3, dashed lines in Fig. 2, flat-plate model in Fig. 10(d) |
| `T3_sensitivity.m` | dispersion curves coloured by η, sensitivities S_αj(f), clean windows, operating points (520 kHz, η_w = 0.052), temperature coefficients, condition number | Fig. 2, Eq. (8), Section 2.3 |
| `T3_check_perturb.m` | perturbation vs direct re-solution (≤ 1.2 %), Korteweg check (1261 vs 1262 m/s) | Appendix A.1, Table S3 |
| `T3_sweep.mat` | output of `T3_sensitivity` (reference) | Fig. 2 |
| `branches_operating.mat` | tracked wall-borne mode and upper / working / lowest liquid-borne branches near the operating frequencies (f, k, c_p, c_g, η, phase- and group-timed sensitivities) | Fig. 2(d), (e); arrival times in Figs. 6 and 8 |

## 2_simulation — Floquet–Bloch FEM (T4)

```
T4_all        % all load cases -> T4_runs/T4_<case>.csv, T4_all.mat (≈ 30 min)
T4_lin        % load linearity (5/10/20 MPa, 50/100 MPa) -> T4_runs/T4_lin.mat
T4_report     % comparison with SAFE (T3_sweep.mat) and the closed forms
```

| Load case (`T4_runs/T4_<case>.csv`) | Definition | Manuscript |
|---|---|---|
| L0 | unloaded | open circles in Fig. 2(a) |
| L1, L2, L2f, L3, L4 | σ_z,ext = 100 MPa; p = 20 MPa closed end; p = 20 MPa with σ_z = 0; both loads; ΔT = +10 K | FEM checks in Fig. 2(b), (c), (e); Table S3 |
| L2p, L2pf | as L2 / L2f with liquid properties frozen at p = 0 | wall-acoustoelastic vs liquid-loading split, Fig. 2(c) |
| L0w1, L1w1 | wall without third-order constants (W1) | Supplementary S1.4 |
| L0cf1400, L0cf1500 | oil sound speed 1400 / 1500 m/s | clean-window shift with c_f, Section 4.1 |

CSV columns: `f_Hz, k_ref_1perm, cp_nat, cg_nat, eta_f, load_id` (natural velocities, wavenumber in material coordinates).
`T4_lin.mat`: secant and tangent pressure coefficients (Table S3) and the axial-stress curvature used in Fig. 9(b).

## 2_simulation — transient EMAT model (E3)

`E3_base.mph` is the complete COMSOL 6.4 model (2D axisymmetric; static magnetic field `mf`, transient
magnetic field `mf2`, solid mechanics in the wall, pressure acoustics in the oil; two probes with the D1T stack and
both coil sets at z_T = 250 mm and z_R = 350 mm; low-reflecting ends). It contains the mesh and the solved
static field (Study 1). Each transient run uses the same model with the global parameters below
(`E3_runs/run_parameters.csv`); set them and compute Study 2.

| Run | f0 | Nc | txw | txl | t_end | frad | Purpose |
|---|---|---|---|---|---|---|---|
| A | 520 kHz | 5 | 1 | 0 | 130 µs | 1 | wall-mode coil set |
| AZ | 520 kHz | 5 | 1 | 0 | 130 µs | 0 | as A, radial Lorentz force off |
| B | 300 kHz | 5 | 0 | 1 | 200 µs | 1 | liquid-mode set, initial operating point |
| BP (B′) | 280 kHz | 5 | 0 | 1 | 200 µs | 1 | liquid-mode set, final operating point |
| C | 325 kHz | 2 | 0 | 1 | 200 µs | 1 | broadband burst, transmission efficiency |

Exported results (output step 0.2 µs):

| File | Export definition | Manuscript |
|---|---|---|
| `E3_runs/static_B.csv` | Study 1, `mf.Br`, `mf.Bz` on r = 6.4995 mm, z = z_T ± 32 mm, step 0.5 mm. RMS deviation from the analytical field (`02_transducer_design/data/E4_paperfig_data.csv`) 3.2 %, peak 0.73 T | Fig. 4(b), Table S6 |
| `E3_runs/run{A,AZ,B,BP,C}_global.csv` | global variables `Ibrst` (transmit current), `V_rx_w`, `V_rx_l` (receive-coil voltages) vs t | Fig. 6(a), (c); Table S6 |
| `E3_runs/runA_line_uz.csv`, `runBP_line_ur.csv` | `w` (run A) and `u` (run BP) on r = 6.4995 mm, z = 275–435 mm, step 1 mm, vs t | f–k maps in Fig. 6(b), (d) |
| `E3_runs/runA_lineneg_uz.csv`, `runBP_lineneg_ur.csv` | same on the opposite side, z = 225–65 mm | ±z symmetry, Table S6 |
| `E3_extra_pp.m` | MATLAB post-processing of the CSVs: phase-timed c_g (residual +20 ns), ±z symmetry (0.9998 / 0.9999), branch shares at 280 kHz (0.08 / 0.52 / 0.40) | Table S6 |
| `plot_fig6.py` | Python: receive-coil voltages, f–k maps with the SAFE branches, f–k purity of run A (99.9 %) | Fig. 6 |

```
cd 2_simulation
python plot_fig6.py        # writes fig6_reproduced.png / .pdf (numpy, scipy, matplotlib)
```

Fig. 5 shows snapshots of `w` (wall) and `p` (oil) from runs A (t = 20 µs) and BP (t = 70 µs).

If `mphstart` is not on the MATLAB path, edit the LiveLink path in `T4_all.m` and `T4_lin.m` (`D:\comsol64\Multiphysics\mli`).
