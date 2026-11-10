# Joint guided-wave measurement of axial stress and internal pressure in liquid-filled pipes

Data, models and scripts accompanying the manuscript

> *Joint Non-Contact Guided-Wave Measurement of Axial Stress and Internal Pressure in Liquid-Filled Pipes*,
> submitted to *Nondestructive Testing and Evaluation* (Taylor & Francis).

The method measures the axial stress and the internal pressure of a liquid-filled pipe jointly with a
wall-borne and a liquid-borne guided wave, both excited by a single dual-polarization EMAT.
This repository contains representative experimental records (Fig. 8), the transducer design
files (Fig. 4), the prestressed fluid–structure-coupled waveguide model (Section 2, Appendix A), and
the finite-element models used to check the model and the transducer (Fig. 2, Fig. 6).

## Repository structure

```
01_fig8_signals/           three measured oscilloscope records and the script that plots Fig. 8
02_transducer_design/      analytical EMAT design, design specification and CAD files (Fig. 4)
03_theory_simulation/
  ├── 1_theory/            SAFE waveguide model and first-order sensitivities (Fig. 2, Eqs. (4)–(8), Appendix A)
  └── 2_simulation/        Floquet–Bloch FEM validation (T4_*) and transient EMAT model (E3_*, Fig. 6)
```

Each folder has its own README with a file-by-file description and the related manuscript items.

## Contents and related manuscript items

| Manuscript item | Folder | Entry point | Software |
|---|---|---|---|
| Fig. 8 (typical signals, spectrograms, predicted arrivals) | `01_fig8_signals` | `python plot_fig8.py` | Python 3 |
| Fig. 4(c), (d), design values of the coil sets | `02_transducer_design` | `python plot_fig4cd.py` | Python 3 |
| Magnet and coil design, Supplementary S3.1 | `02_transducer_design/matlab` | `E1_force_field`, `E2_coil_spectrum`, `E4_paperfig_data` | MATLAB |
| Fig. 4(a) probe assembly | `02_transducer_design/cad` | STEP / STL files | any CAD |
| Fig. 2, operating points, Eq. (8), Appendix A | `03_theory_simulation/1_theory` | `T1_thin_shell`, `T3_sensitivity`, `T3_check_perturb` | MATLAB |
| FEM points in Fig. 2, FEM curve in Fig. 9(b) | `03_theory_simulation/2_simulation` | `T4_all`, `T4_report` (results included in `T4_runs/`) | COMSOL + LiveLink for MATLAB |
| Fig. 6 (transient EMAT simulation) | `03_theory_simulation/2_simulation` | `python plot_fig6.py` (results included in `E3_runs/`), model `E3_base.mph` | Python 3; COMSOL 6.4 to rerun |

The Python scripts need only `numpy`, `scipy` and `matplotlib`:

```
pip install numpy scipy matplotlib
python 01_fig8_signals/plot_fig8.py
python 02_transducer_design/plot_fig4cd.py
python 03_theory_simulation/2_simulation/plot_fig6.py
```

The MATLAB scripts need no toolboxes (R2021a or later). Keep the folder structure unchanged,
because the scripts locate each other through relative paths. Comments in the MATLAB files are
partly in Chinese; the READMEs and the Python scripts are in English.

## Citation

If you use these files, please cite the article (reference to be added after publication).

## License

Code (`*.m`, `*.py`) is released under the MIT License (see `LICENSE`).
Data, model and CAD files are released under the Creative Commons Attribution 4.0 International
License (CC BY 4.0).
