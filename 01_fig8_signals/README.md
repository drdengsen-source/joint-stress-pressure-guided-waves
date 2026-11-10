# Fig. 8 — typical received signals (data and reproduction script)

This folder reproduces Fig. 8 of the manuscript: (a) wall-borne record at zero load, (b) liquid-borne
records at 0 and 18.35 MPa, (c), (d) STFT spectrograms, with the time gates and predicted arrival times.

```
python plot_fig8.py            # writes fig8_reproduced.png / .pdf
```
Requires Python 3 with numpy, scipy and matplotlib. `fig8_reproduced.png` is the output of this script.

## Files

| File | Content | Used for |
|---|---|---|
| `data/W520_zero-load.csv` | Record 55: wall-mode coil set driven (5 cycles, 520 kHz), σ_z,ext = 0, p = 0 | Fig. 8(a), (c) |
| `data/L280_zero-load.csv` | Record 56: liquid-mode coil set driven (5 cycles, 280 kHz), σ_z,ext = 0, p = 0 | Fig. 8(b) "0 MPa", (d) |
| `data/L280_sz200MPa_p18.35MPa.csv` | Record 204: liquid-mode coil set driven, σ_z,ext = 200 MPa, p = 18.35 MPa | Fig. 8(b) "18.35 MPa" |
| `data/branch_arrivals.csv` | f0, phase and group velocity, liquid energy fraction of the wall-borne mode and the three liquid-borne branches at the operating frequencies (SAFE model, Section 2; extracted from `03_theory_simulation/1_theory/branches_operating.mat`) | predicted arrival lines |
| `data/records_used.csv` | Load state, wall temperature and instrument settings of the three records | provenance |
| `plot_fig8.py` | Reproduction script | — |
| `readMSO24csv.m` | MATLAB reader for the oscilloscope CSV format | optional |

Oscilloscope CSV (Tektronix MSO24): 13 metadata lines, then `TIME, CH1, CH2, CH3`; 25 000 samples,
8 ns interval, −10 to 190 µs, high-resolution mode, 512 averages.
CH1 = transmit monitor, **CH2 = wall-mode receive coil (used in Fig. 8; both families are received
on this coil, Section 4.1)**, CH3 = liquid-mode receive coil.

## Processing (as in the manuscript)

- Time gates: 15–40 µs (wall-borne), 60–86 µs (liquid-borne), Section 4.2.
- Envelope: Hilbert envelope after a 4th-order Butterworth band-pass of 0.4–1.6 f0.
- Spectrograms: Hann-window STFT (8 µs for 520 kHz, 14 µs for 280 kHz), normalised to the maximum, 6 dB contours down to −30 dB.
- Predicted arrival of the burst centre: t = (N/2)/f0 + L/c_g with N = 5 cycles, L = 100 mm and the
  group velocities in `branch_arrivals.csv`: wall-borne 24.7 µs; lowest 67.0 µs, working 76.4 µs, upper 92.4 µs.
