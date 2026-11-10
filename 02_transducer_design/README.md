# Dual-polarization EMAT — transducer design (Section 3.1, Fig. 4, Supplementary S3.1)

Analytical model of the magnet stack and the two coil sets: static field of radially magnetized rings,
eddy current induced by each coil set, Lorentz force densities and their spatial spectra.

## Quick start (no MATLAB needed)

```
python plot_fig4cd.py     # Fig. 4(c), (d) from data/E4_paperfig_data.csv; prints the design values
```
Requires Python 3 with numpy and matplotlib. Output: `fig4cd_reproduced.png / .pdf`.

## Full regeneration (MATLAB R2021a or later; runs in about 1–2 min)

```
cd matlab
E1_force_field      % static field of all magnet stacks + force densities   -> E1_fields.mat
E2_coil_spectrum    % coil-width sweep and final coil design, spatial spectra -> E2_coils.mat
E4_paperfig_data    % exports the D1T curves                                -> E4_paperfig_data.csv
```
The regenerated `E4_paperfig_data.csv` is identical to `data/E4_paperfig_data.csv`.

## Files and where they appear in the manuscript

| File | What it computes | Manuscript |
|---|---|---|
| `matlab/loopB.m` | field of a circular current loop (elliptic integrals) | basis of all fields |
| `matlab/ringB.m` | field of a uniformly magnetized ring (radial or axial) via equivalent surface currents | Fig. 4(b) analytical counterpart |
| `matlab/coilK.m` | eddy-current sheet induced in the pipe wall by a coil element (image approximation) | Fig. 4(c) |
| `matlab/params_pipe.m` | pipe and liquid parameters (Table 1) | — |
| `matlab/E1_force_field.m` | B_r, B_z on the outer wall and force densities for the candidate stacks (thin rings, same / alternating polarity, axially magnetized rings, thick three-ring D1T stack) | Supplementary S3.1 and Table S5 (stack comparison) |
| `matlab/E2_coil_spectrum.m` | element-width sweep and spatial spectra of the final coil design for the thin two-ring and thick three-ring stacks | Supplementary Fig. S1, Table S5 (10 dB / 12 dB gains); choice of the λ_w/2 element width |
| `matlab/E4_paperfig_data.m` | exports z, B_r, B_z and the four force densities of the D1T stack | Fig. 4(c), (d) |
| `data/E4_paperfig_data.csv` | the exported curves (3201 points, −32…32 mm) | Fig. 4(c), (d) |
| `plot_fig4cd.py` | force densities and spatial spectra from the CSV | Fig. 4(c), (d) |
| `D1T_design_spec.md` | magnet, coil, bobbin, electronics and acquisition settings | Sections 3.1, 4.1; Supplementary S3.3 |
| `cad/probe_A_assembly.step` | complete probe assembly (half-shells, 4 × M3 screws and nuts, three ring magnets, two coil sets, bobbin, pipe section); dimensions in mm | Fig. 4(a) |
| `cad/half_shell_upper.stl`, `cad/half_shell_lower_socket.stl` | the two printed half-shells; the ribs in the pockets set the 10.4 mm magnet pitch (assembly coordinates) | Fig. 4(a) |
| `cad/bobbin.stl` | coil bobbin, inner Ø13.2 / outer Ø14.0 / length 52 | Fig. 4(a) |

The numerical static field of Fig. 4(b) and its 3.2 % RMS comparison with this analytical field are part of the
transient FEM model (folder `03`, `E3_static.m`, `E3_runs/static_B.csv`).

## Model assumptions

Uniform magnetization (B_r = 1.35 T), relative permeability 1 everywhere; eddy current from the image
approximation of a plane conductor (curvature error O(h/r_o) ≈ 8 %). Fig. 4(b) uses the full numerical model.
