"""Reproduce Fig. 6 of the manuscript from the transient EMAT simulation (runs A and B').

Usage:  python plot_fig6.py [uz|ur]   (component for the run B' f-k map, default from FK_COMP_BP)
Inputs (relative to this script):
  E3_runs/runA_global.csv, E3_runs/runBP_global.csv   t, I_tx, V_rx_wall, V_rx_liq
  E3_runs/runA_line_uz.csv, E3_runs/runBP_line_ur.csv outer-wall displacement on r = 6.4995 mm, z = 275-435 mm
  ../1_theory/T3_sweep.mat                               SAFE dispersion (overlay)
  ../1_theory/branches_operating.mat is not needed; arrival times use the group velocities below.
Requires numpy, scipy, matplotlib.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "E3_runs")
FK_COMP_BP = sys.argv[1] if len(sys.argv) > 1 else "ur"
L = 0.100                                   # m, transmitter-receiver distance
CG_280 = {"lowest": 1722.75, "working": 1482.86, "upper": 1197.36}   # m/s, SAFE, 280 kHz


def read_global(name):
    d = np.genfromtxt(os.path.join(RUNS, name), delimiter=",", comments="%")
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def read_line(name):
    with open(os.path.join(RUNS, name), encoding="utf-8") as fh:
        hdr = [next(fh) for _ in range(2)][1]
    t = np.array([float(s.split("t=")[1]) for s in hdr.split(",") if "t=" in s])
    d = np.genfromtxt(os.path.join(RUNS, name), delimiter=",", comments="%")
    return d[:, 1], t, d[:, 2:]


def fk_map(U, z, t):
    """2D FFT with a Hann window in z and 5 % cosine tapers in t (as E3_extra_pp.m)."""
    nz, nt = U.shape
    wz = 0.5 * (1 - np.cos(2 * np.pi * np.arange(nz) / (nz - 1)))[:, None]
    wt = np.ones(nt)
    nr = max(2, round(0.05 * nt))
    ramp = 0.5 * (1 - np.cos(np.pi * np.arange(nr) / (nr - 1)))
    wt[:nr] = ramp
    wt[-nr:] = ramp[::-1]
    Uw = (U - U.mean(axis=1, keepdims=True)) * wz * wt[None, :]
    nk = int(2 ** np.ceil(np.log2(4 * nz)))
    nf = int(2 ** np.ceil(np.log2(2 * nt)))
    F = np.fft.fft2(Uw, (nk, nf))
    dz, dt = z[1] - z[0], t[1] - t[0]
    kax = np.fft.fftfreq(nk, dz)
    fax = np.arange(nf) / (nf * dt)
    fsel = fax <= 1 / (2 * dt)
    A = np.abs(F[:, fsel]).T
    k = -2 * np.pi * kax
    o = np.argsort(k)
    return fax[fsel], k[o], A[:, o]


def env(x):
    return np.abs(hilbert(x))


# ------------------------------------------------------------------ data
tA, IA, VwA, VlA = read_global("runA_global.csv")
tB, IB, VwB, VlB = read_global("runBP_global.csv")
zA, tlA, UA = read_line("runA_line_uz.csv")
zB, tlB, UB = read_line(f"runBP_line_{FK_COMP_BP}.csv")
fA, kA, FA = fk_map(UA, zA, tlA)
fB, kB, FB = fk_map(UB, zB, tlB)

T3 = loadmat(os.path.join(HERE, "..", "1_theory", "T3_sweep.mat"), squeeze_me=True)
sf, scp, seta = T3["allF"], T3["allCP"], T3["allETA"]
sk = 2 * np.pi * sf / scp                    # rad/m

# purity of run A (Table S6): within 470-570 kHz, fraction of f-k energy with |k - 2 pi f / 5403| <= 15 %
fsel = (fA >= 470e3) & (fA <= 570e3)
kw = 2 * np.pi * fA[fsel] / 5403
inband = np.abs(kA[None, :] - kw[:, None]) <= 0.15 * kw[:, None]
E = FA[fsel] ** 2
print("run A: energy fraction within c_p = 5403 m/s +-15 %% (470-570 kHz): %.1f %%" % (100 * E[inband].sum() / E.sum()))

# predicted arrivals of the burst centre at 280 kHz
t0B = 2.5 / 280e3
arr = {b: (t0B + L / c) * 1e6 for b, c in CG_280.items()}
print("run B' predicted arrivals (us):", {b: round(v, 1) for b, v in arr.items()})

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9})
ORANGE, BLUE, GREEN = "#D55E00", "#0072B2", "#009E73"
fig, ax = plt.subplots(2, 2, figsize=(7.48, 3.9), gridspec_kw=dict(hspace=0.42, wspace=0.18,
                                                                   left=0.07, right=0.985, top=0.95, bottom=0.11))


def tpanel(a, t, I, Vmain, Vcross, col, gate, xmax, lab_main, lab_cross):
    t = t * 1e6
    n = np.max(np.abs(Vmain))
    a.axvspan(*gate, color="0.88", zorder=0)
    a.plot(t, I / np.max(np.abs(I)), ":", color="0.45", lw=0.9)
    a.plot(t, Vcross / n, "-.", color=GREEN, lw=0.7)
    a.plot(t, Vmain / n, color=col, lw=0.7)
    a.plot(t, env(Vmain) / n, color="k", lw=1.0)
    a.set(xlim=(0, xmax), ylim=(-1.3, 1.3), xlabel=r"$t$ (μs)", ylabel="normalised")
    a.text(0.98, 0.93, lab_main + ", envelope in black", color=col, transform=a.transAxes, ha="right", fontsize=8)
    a.text(0.98, 0.05, r"$I_\mathrm{tx}$ (dotted)   " + lab_cross + " (dash-dot)", color="0.4",
           transform=a.transAxes, ha="right", fontsize=8)
    a.text(gate[0] + 1, 1.1, f"gate {gate[0]}–{gate[1]} μs", fontsize=8)


tpanel(ax[0, 0], tA, IA, VwA, VlA, ORANGE, (15, 40), 130, "wall-mode coil", "liquid-mode coil")
tpanel(ax[0, 1], tB, IB, VlB, VwB, BLUE, (60, 86), 200, "liquid-mode coil", "wall-mode coil")
for b, ls in [("lowest", ":"), ("working", "--"), ("upper", ":")]:
    ax[0, 1].axvline(arr[b], color="k", ls=ls, lw=0.8)
    ax[0, 1].text(arr[b], -1.25, f"{b} {arr[b]:.1f}", rotation=90, fontsize=7, va="bottom", ha="right")


def fkpanel(a, f, k, A, flim, klim):
    S = 20 * np.log10(A / A.max() + 1e-12)
    m = (f >= flim[0] * 1e3) & (f <= flim[1] * 1e3)
    S = np.ma.masked_less(S[m], -30)
    a.pcolormesh(k / 1e3, f[m] / 1e3, S, cmap="jet", vmin=-30, vmax=0, shading="auto", rasterized=True)
    liq = seta >= 0.5
    a.plot(sk[liq] / 1e3, sf[liq] / 1e3, "o", mfc="none", mec="0.25", ms=2.5, mew=0.5)
    a.plot(sk[~liq] / 1e3, sf[~liq] / 1e3, "s", color="0.15", ms=2.2)
    a.set(xlim=klim, ylim=flim, xlabel=r"$k$ (rad/mm)", ylabel=r"$f$ (kHz)")


fkpanel(ax[1, 0], fA, kA, FA, (300, 1000), (0, 1.8))
fkpanel(ax[1, 1], fB, kB, FB, (50, 460), (0, 1.6))
for a, lab in zip(ax.flat, "acbd"):
    a.text(-0.005, 1.0, f"({lab})", transform=a.transAxes, ha="right", va="bottom", fontweight="bold")

out = os.path.join(HERE, "fig6_reproduced")
fig.savefig(out + ".png", dpi=600)
fig.savefig(out + ".pdf")
print("saved", out + ".png")
