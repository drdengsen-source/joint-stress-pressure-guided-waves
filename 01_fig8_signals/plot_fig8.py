"""Reproduce Fig. 8 of the manuscript (typical received signals and STFT spectrograms).

Usage:  python plot_fig8.py [data_dir] [out_prefix]
        data_dir defaults to ./data next to this script; out_prefix to ./fig8_reproduced

Inputs (in data_dir):
  W520_zero-load.csv              wall-mode coil set driven, 5 cycles at 520 kHz, sigma_z,ext = 0, p = 0
  L280_zero-load.csv              liquid-mode coil set driven, 5 cycles at 280 kHz, sigma_z,ext = 0, p = 0
  L280_sz200MPa_p18.35MPa.csv     liquid-mode coil set driven, sigma_z,ext = 200 MPa, p = 18.35 MPa
  branch_arrivals.csv             f0 and group velocity of each branch (from the SAFE model)
Oscilloscope CSV columns: TIME (s), CH1 transmit monitor (V), CH2 wall-mode receive coil (V),
CH3 liquid-mode receive coil (V). Fig. 8 uses CH2 (both families are received by the wall-mode coil).
Requires numpy, scipy, matplotlib.
"""
import csv
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.signal import butter, sosfiltfilt, hilbert, stft

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "fig8_reproduced")

L_GAUGE = 0.100      # m, centre-to-centre probe distance (Section 4.1)
N_CYC = 5            # cycles per burst
GATE_W = (15, 40)    # us, wall-borne time gate (Section 4.2)
GATE_L = (60, 86)    # us, liquid-borne time gate


def read_mso(path):
    """Read a Tektronix MSO24 CSV (metadata lines followed by a TIME,CH1,... table)."""
    with open(path, encoding="utf-8-sig") as fh:
        lines = fh.readlines()
    i0 = next(i for i, ln in enumerate(lines) if ln.startswith("TIME"))
    data = np.genfromtxt(lines[i0 + 1:], delimiter=",")
    cols = lines[i0].strip().split(",")
    return {c: data[:, j] for j, c in enumerate(cols)}


def envelope(x, fs, f0):
    """Hilbert envelope after the 0.4-1.6 f0 band-pass used for timing."""
    sos = butter(4, [0.4 * f0, 1.6 * f0], btype="band", fs=fs, output="sos")
    return np.abs(hilbert(sosfiltfilt(sos, x)))


def spectrogram_db(x, fs, win_us):
    nper = int(round(win_us * 1e-6 * fs))
    f, t, Z = stft(x, fs=fs, window="hann", nperseg=nper, noverlap=nper - nper // 16,
                   nfft=8 * nper, boundary="zeros", padded=True)
    S = 20 * np.log10(np.abs(Z) / np.abs(Z).max() + 1e-12)
    return f, t, S


def arrival_us(f0, cg):
    """Predicted arrival of the burst centre: t0 + L/cg with t0 = (N/2)/f0."""
    return (N_CYC / 2 / f0 + L_GAUGE / cg) * 1e6


# ---------------------------------------------------------------- data
w0 = read_mso(os.path.join(DATA, "W520_zero-load.csv"))
l0 = read_mso(os.path.join(DATA, "L280_zero-load.csv"))
l1 = read_mso(os.path.join(DATA, "L280_sz200MPa_p18.35MPa.csv"))
with open(os.path.join(DATA, "branch_arrivals.csv"), encoding="utf-8") as fh:
    br = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
cg = {r["branch"]: float(r["cg0_m_per_s"]) for r in br}
f0 = {r["branch"]: float(r["f0_Hz"]) for r in br}
t_wall = arrival_us(f0["wall"], cg["wall"])
t_low = arrival_us(f0["lowest"], cg["lowest"])
t_work = arrival_us(f0["working"], cg["working"])
t_up = arrival_us(f0["upper"], cg["upper"])
print(f"predicted arrivals (us): wall {t_wall:.2f}, lowest {t_low:.2f}, working {t_work:.2f}, upper {t_up:.2f}")

tw = w0["TIME"] * 1e6
tl = l0["TIME"] * 1e6
fs = 1 / (w0["TIME"][1] - w0["TIME"][0])
vw = w0["CH2"] * 1e3            # mV
vl0 = l0["CH2"] * 1e3
vl1 = l1["CH2"] * 1e3

# ---------------------------------------------------------------- style
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                     "font.size": 9, "axes.linewidth": 0.8})
ORANGE, BLUE, PINK = "#D55E00", "#0072B2", "#CC79A7"
cmap = ListedColormap(["#0b1d5c", "#1f3f7a", "#4a6a8e", "#7d968d", "#bccb9c", "#f7f3d5"])
cmap.set_under("#0b1d5c")
levels = [-30, -24, -18, -12, -6, 0]

fig, ax = plt.subplots(2, 2, figsize=(7.48, 3.45), sharex="col",
                       gridspec_kw=dict(hspace=0.05, wspace=0.2, left=0.075, right=0.985,
                                        bottom=0.12, top=0.95))

# (a) wall-borne record at zero load
a = ax[0, 0]
pre = tw < GATE_W[0]
a.plot(tw[pre], vw[pre], color=ORANGE, lw=0.6, alpha=0.45)
a.plot(tw[~pre], vw[~pre], color=ORANGE, lw=0.6)
a.plot(tw, envelope(w0["CH2"], fs, 520e3) * 1e3, color="k", lw=0.9)
a.axvspan(*GATE_W, color="0.88", zorder=0)
for g in GATE_W:
    a.axvline(g, color="0.6", lw=0.8)
a.axvline(t_wall, color="k", ls="--", lw=0.8)
a.set(xlim=(0, 120), ylim=(-620, 620), yticks=range(-600, 601, 200), ylabel="wall-mode coil (mV)")
a.text(43, 520, f"gate {GATE_W[0]}–{GATE_W[1]} μs")
a.text(0.98, 0.17, "wall set, 5 cycles at 520 kHz", transform=a.transAxes, ha="right")
a.text(0.98, 0.06, f"dashed: predicted arrival {t_wall:.1f} μs", transform=a.transAxes, ha="right")

# (b) liquid-borne records at 0 and 18.35 MPa
b = ax[0, 1]
b.plot(tl, vl0, color=BLUE, lw=0.6)
b.plot(tl, vl1, color=PINK, lw=0.6)
b.plot(tl, envelope(l0["CH2"], fs, 280e3) * 1e3, color="k", lw=0.9)
b.axvspan(*GATE_L, color="0.88", zorder=0)
for g in GATE_L:
    b.axvline(g, color="0.6", lw=0.8)
for tt, lab, yy in [(t_low, "lowest", 190), (t_work, "working", 150)]:
    b.axvline(tt, color="k", ls=":", lw=0.8)
    b.text(tt, yy, lab, ha="center", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
b.set(xlim=(0, 190), ylim=(-210, 210), ylabel="wall-mode coil (mV)")
b.text(0.14, 0.92, "0 MPa", color=BLUE, transform=b.transAxes)
b.text(0.14, 0.80, "18.35 MPa", color=PINK, transform=b.transAxes)
b.text(GATE_L[1] + 2, 180, f"gate {GATE_L[0]}–{GATE_L[1]} μs")
b.text(0.02, 0.05, "liquid set, 5 cycles at 280 kHz", transform=b.transAxes,
       bbox=dict(fc="w", ec="none", pad=0.5))

# (c), (d) STFT spectrograms of the zero-load records
for axx, rec, win, flim, gate in [(ax[1, 0], w0, 8.0, (150, 1050), GATE_W),
                                  (ax[1, 1], l0, 14.0, (80, 560), GATE_L)]:
    f, t, S = spectrogram_db(rec["CH2"], fs, win)
    t = t * 1e6 + rec["TIME"][0] * 1e6
    m = (f >= flim[0] * 1e3) & (f <= flim[1] * 1e3)
    cs = axx.contourf(t, f[m] / 1e3, S[m], levels=levels, cmap=cmap, extend="min")
    axx.contour(t, f[m] / 1e3, S[m], levels=levels[1:-1], colors="w", linewidths=0.3,
                linestyles="--", alpha=0.6)
    for g in gate:
        axx.axvline(g, color="w", ls=":", lw=1.0)
    axx.set(ylim=flim, xlabel=r"$t$ (μs)", ylabel=r"$f$ (kHz)")
    axx.set_facecolor("#0b1d5c")
cax = ax[1, 0].inset_axes([0.62, 0.24, 0.25, 0.06])
cb = fig.colorbar(cs, cax=cax, orientation="horizontal", ticks=[-30, -15, 0])
cb.ax.tick_params(labelsize=8, colors="w")
cax.text(1.08, 0.5, "dB", color="w", fontsize=8, va="center", transform=cax.transAxes)

for axx, lab in zip(ax.flat, "abcd"):
    axx.text(-0.005, 1.0, f"({lab})", transform=axx.transAxes, ha="right", va="bottom" if lab in "ab" else "top",
             fontweight="bold", color="k") if lab in "ab" else axx.text(0.01, 0.97, f"({lab})", transform=axx.transAxes,
             va="top", fontweight="bold", color="w")

fig.savefig(OUT + ".png", dpi=600)
fig.savefig(OUT + ".pdf")
print("saved", OUT + ".png", OUT + ".pdf")
