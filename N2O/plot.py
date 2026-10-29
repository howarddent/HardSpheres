# This file is part of HardSpheres.
# Copyright (C) 2026 Howard Dent
# SPDX-License-Identifier: GPL-3.0-or-later
#
# HardSpheres is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version. It is distributed WITHOUT ANY WARRANTY; see
# the file LICENSE at the top of the repository for details.

"""P-V isotherms of N2O from the stepped-LJ sphere model, mapped by
corresponding states with the model's own critical point.

Data in data/ (or the folder given as the first argument):
  iso2.csv        1M-sphere fixed-volume isotherm points
  coex.csv        200k coexistence runs at rho 0.31: P (the vapour pressure below Tc)
  histall_T.csv   their whole-box cube-density histograms: the coexisting densities
"""
import csv, glob, math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

SRC = os.path.dirname(os.path.abspath(__file__))
HERE = sys.argv[1] if len(sys.argv) > 1 else f"{SRC}/data"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{SRC}/n2o_isotherms.png"

TC_N2O, PC_N2O, M_N2O = 309.52, 72.447, 44.013     # K, bar, g/mol (NIST)
R = 0.083144626                                     # L bar / (mol K)
BETA = 0.325
TEMPS = [(0.9134, "0 °C"), (0.9803, "20 °C"), (1.0350, "36.4 °C  (critical)"), (1.2478, "100 °C")]
SUBCRIT = (0.9134, 0.9803)


def read(path):
    with open(path) as f:
        return [[float(x) for x in r] for r in csv.reader(f) if r and not r[0].startswith("#")]


# below rho 0.1 the first sweep (the walls hardly matter there); the rest from
# the rerun, at the interior density
iso3 = read(f"{HERE}/iso3.csv")
iso = [r for r in read(f"{HERE}/iso2.csv") if r[1] < 0.1 - 1e-6]
iso += [[r[0], r[11]] + r[2:11] for r in iso3]
coexP = {round(r[0], 4): r[1] for r in read(f"{HERE}/coex.csv")}


# coexisting densities: the two peaks of each whole-box histogram (smoothed
# over 4 bins: cube counts alias every other bin), split at rho 0.29
def peaks(T):
    h = np.array(read(f"{HERE}/histall_{T:.4f}.csv"))
    rho, fr = h[:, 0], np.convolve(h[:, 1], np.ones(4) / 4, mode="same")
    lo, hi = rho < 0.29, rho > 0.29
    return rho[lo][np.argmax(fr[lo])], rho[hi][np.argmax(fr[hi])]


COEX_T = [0.9134, 0.96, 0.9803, 0.99, 1.02]
coex = {T: peaks(T) for T in COEX_T}

# scaling fit: (rho_l - rho_v)^(1/beta) = A (Tc - T); diameter rho_c + D (Tc - T)
t = np.array(COEX_T)
y = np.array([(coex[T][1] - coex[T][0]) ** (1 / BETA) for T in COEX_T])
# Tc pinned at the critical isotherm's temperature (the 1M isotherm is flat
# there); the 200k histograms alone put it ~0.007 lower (finite-size blur)
TC = 1.035
A = np.sum(y * (TC - t)) / np.sum((TC - t) ** 2)
dm = np.array([(coex[T][0] + coex[T][1]) / 2 for T in COEX_T])
D, RHOC = np.polyfit(TC - t, dm, 1)
# vapour pressure: ln P = a - b/T through the coexistence runs; Pc at Tc
pts = [(T, coexP[round(T, 4)]) for T in COEX_T]
b_, a_ = np.polyfit([1 / T for T, _ in pts], [math.log(p) for _, p in pts], 1)
psat = lambda T: math.exp(a_ + b_ / T)
PC = psat(TC)

KP = PC_N2O / PC                         # bar per model pressure unit
KV = R * TC_N2O / TC / KP                # L/mol per model volume per sphere (keeps Z)
V = lambda rho: KV / np.asarray(rho)

# the model's dome
Td = np.concatenate([np.linspace(0.9134, TC - 0.01, 50), TC - np.geomspace(0.01, 1e-6, 30)])
half = (A * (TC - Td)) ** BETA / 2
mid = RHOC + D * (TC - Td)
Pd = KP * np.array([psat(T) for T in Td])
dome_V = np.concatenate([V(mid + half), [V(RHOC)], V(mid - half)[::-1]])
dome_P = np.concatenate([Pd, [PC_N2O], Pd[::-1]])

# N2O's own saturation curve (NIST), from the app's data file
nist = [tuple(float(x) for x in m.groups()) for m in re.finditer(
    r"\(([\d.]+), ([\d.]+), ([\d.]+), ([\d.]+)\)",
    open(f"{SRC}/n2odata.inc").read())]
nist_V = [M_N2O / r[2] for r in nist] + [M_N2O / r[3] for r in nist[::-1]]
nist_P = [r[1] for r in nist] + [r[1] for r in nist[::-1]]

# --- the plot -------------------------------------------------------------
INK, INK2, MUTED, GRID, SURF, DOME = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb", "#efeeea"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
PMAX, VMIN, VMAX = 150, 0.035, 7.0

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
fig, ax = plt.subplots(figsize=(11.5, 7.5), dpi=160)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
ax.set_xscale("log"); ax.set_xlim(VMIN, VMAX); ax.set_ylim(0, PMAX)
ax.grid(True, which="major", color=GRID, lw=0.8)
ax.grid(True, which="minor", axis="x", color=GRID, lw=0.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

ax.fill(dome_V, dome_P, color=DOME, zorder=0, lw=0)
ax.plot(dome_V, dome_P, color=INK2, lw=1.5, zorder=2, label="Liquid–vapour dome (model)")
ax.plot(nist_V, nist_P, color=MUTED, lw=1.2, ls=(0, (4, 3)), zorder=2,
        label="N₂O's actual dome (NIST), for comparison")

report = []
# ideal gas P = RT/V at each isotherm's real temperature (the volume scale keeps Z,
# so this is also the model's ideal gas)
Vi = np.geomspace(VMIN, VMAX, 200)
for i, ((T, _), col) in enumerate(zip(TEMPS, SERIES)):
    ax.plot(Vi, R * T * TC_N2O / TC / Vi, color=col, lw=1.1, ls=(0, (1.5, 2.5)), alpha=0.9, zorder=3,
            label="Ideal gas, P = RT/V (dotted, same colours)" if i == 0 else None)

for (T, name), col in zip(TEMPS, SERIES):
    p = sorted((r for r in iso if abs(r[0] - T) < 1e-4), key=lambda r: r[1])
    rho = np.array([r[1] for r in p]); P = KP * np.array([r[2] for r in p])
    if T in SUBCRIT:
        ps = KP * coexP[round(T, 4)]
        rv, rl = coex[T]
        gas = rho < 0.3
        liq = (rho > 0.3) & (P >= ps)
        meta = (rho > 0.3) & (P < ps)
        xs = np.concatenate([V(rho[liq])[::-1][::-1], [V(rl), V(rv)], V(rho[gas])[::-1]])
        ys = np.concatenate([P[liq], [ps, ps], P[gas][::-1]])
        order = np.argsort(xs)
        ax.plot(xs[order], ys[order], color=col, lw=2, zorder=4)
        ax.plot(V(rho[gas | liq]), P[gas | liq], "o", ms=5.5, color=col, mec=SURF, mew=1.2, zorder=5,
                label=f"{name}  (T* = {T:.3f})")
        if meta.any():
            ax.plot(V(rho[meta]), P[meta], "o", ms=5.5, mfc=SURF, mec=col, mew=1.4, zorder=5)
        report.append(f"{name}: vapour pressure {ps:.1f} bar; coexisting {V(rl):.4f} and {V(rv):.3f} L/mol"
                      f"{'; ' + str(meta.sum()) + ' metastable liquid point(s)' if meta.any() else ''}")
        xr, yr = math.sqrt(V(rl) * V(rv)) * 1.0, ps
    else:
        ax.plot(V(rho), P, color=col, lw=2, zorder=4)
        ax.plot(V(rho), P, "o", ms=5.5, color=col, mec=SURF, mew=1.2, zorder=5,
                label=f"{name}  (T* = {T:.3f})")
        k = np.argmin(abs(V(rho) - (0.13 if T > 1.1 else 0.3)))
        xr, yr = V(rho[k]), P[k]
    ax.annotate(name.split("  ")[0], (xr, yr), xytext=(4, 8), textcoords="offset points", va="bottom",
                ha="left", color=INK, fontsize=10, fontweight="bold", zorder=6, annotation_clip=False)

ax.plot([V(RHOC)], [PC_N2O], "o", ms=7, color=INK, zorder=7)
ax.annotate("critical point", (V(RHOC), PC_N2O), xytext=(0, -12), textcoords="offset points",
            ha="center", va="top", color=INK, fontsize=10, zorder=7)

kw = dict(color=INK, fontsize=13, fontweight="bold", ha="center", va="center", zorder=6)
ax.text(0.0425, 70, "LIQUID", rotation=90, **kw)
ax.text(0.125, 42, "LIQUID + VAPOUR", **kw)
ax.text(0.85, 10, "VAPOUR", **kw)
ax.text(0.9, 125, "GAS", **kw)
ax.text(0.9, 116, "above the critical temperature", color=INK2, fontsize=9.5, ha="center", va="center", zorder=6)

ax.set_xlabel("Molar volume (L/mol, log scale)")
ax.set_ylabel("Pressure (bar)")
ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1, 2, 5])
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
ax.set_title("Isotherms of nitrous oxide, from a 1,000,000-sphere simulation",
             loc="left", color=INK, fontsize=14, fontweight="bold", pad=24)
ax.text(0, 1.018, "Stepped Lennard-Jones spheres mapped onto N₂O by corresponding states "
        f"(model Tc* = {TC:.3f}, Pc* = {PC:.4f}).",
        transform=ax.transAxes, color=INK2, fontsize=9.5)
leg = ax.legend(loc="upper right", frameon=True, fontsize=9)
leg.get_frame().set_facecolor(SURF); leg.get_frame().set_edgecolor(GRID)
fig.tight_layout()
fig.savefig(OUT, facecolor=SURF)

print(f"Tc* {TC:.4f}, rho_c* {RHOC:.4f}, A {A:.3f}, D {D:.3f}; ln P = {a_:.3f} {b_:+.3f}/T, Pc* {PC:.5f}")
print(f"KP {KP:.1f} bar, KV {KV:.5f} L/mol; Vc {V(RHOC):.4f} L/mol (NIST 0.0974), Zc {PC / (RHOC * TC):.3f} (N2O 0.274)")
for T in COEX_T:
    print(f"  coexistence T* {T}: vapour {coex[T][0]:.3f}, liquid {coex[T][1]:.3f}, P {coexP[round(T, 4)]:.5f}")
for line in report:
    print(" ", line)
for r in iso:
    bl = r[3:7]
    if abs(bl[3] - bl[0]) / max(abs(r[2]), 1e-9) > 0.03:
        print(f"  drift: T* {r[0]} rho {r[1]}: blocks {bl}")
print("wrote", OUT)
