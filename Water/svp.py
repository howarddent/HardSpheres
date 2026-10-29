# This file is part of HardSpheres.
# Copyright (C) 2026 Howard Dent
# SPDX-License-Identifier: GPL-3.0-or-later
#
# HardSpheres is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version. It is distributed WITHOUT ANY WARRANTY; see
# the file LICENSE at the top of the repository for details.

"""Saturated vapour pressure of water against temperature, from the NIST
Chemistry WebBook (IAPWS-95 via REFPROP), triple point to critical point."""
import csv, os, numpy as np, matplotlib
matplotlib.use("Agg")
HERE = os.path.dirname(os.path.abspath(__file__))
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

rows = list(csv.reader(open(os.path.join(HERE, "water_sat_T.tsv")), delimiter="\t"))[1:]
T = np.array([float(r[0]) for r in rows]) - 273.15
P = np.array([float(r[1]) for r in rows]) * 100    # bar -> kPa

INK, INK2, MUTED, GRID, SURF, BLUE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), dpi=160)
fig.patch.set_facecolor(SURF)
pts = [(0.01, 0.6117, "triple point\n0.01 °C, 0.611 kPa", (10, 12)),
       (100.0, 101.418, "100 °C: 101.4 kPa\n(boils at 1 atm, 99.97 °C)", (10, 8)),
       (373.946, 22064, "critical point\n373.9 °C, 22 064 kPa", (-10, -6))]
for ax, log in zip(axes, (False, True)):
    ax.set_facecolor(SURF)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.plot(T, P, color=BLUE, lw=2.2, zorder=3)
    ax.set_xlim(0, 390); ax.set_xlabel("Temperature (°C)")
    if log:
        ax.set_yscale("log"); ax.set_ylim(0.3, 60000)
        ax.grid(True, which="minor", axis="y", color=GRID, lw=0.4)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
        ax.set_title("Log pressure: the whole range", loc="left", color=INK, fontsize=11.5, fontweight="bold")
    else:
        ax.set_ylim(0, 24000)
        ax.set_title("Linear pressure", loc="left", color=INK, fontsize=11.5, fontweight="bold")
    ax.set_ylabel("Saturated vapour pressure (kPa)")
    for t, p, txt, off in pts:
        if not log and p < 5 and t < 50: continue
        ax.plot([t], [p], "o", ms=7, color=BLUE, mec=SURF, mew=1.5, zorder=4)
        ha = "right" if off[0] < 0 else "left"
        ax.annotate(txt, (t, p), xytext=off, textcoords="offset points", ha=ha,
                    va="top" if off[1] < 0 else "bottom", fontsize=9, color=INK)
fig.suptitle("Saturated vapour pressure of water", x=0.012, ha="left", y=0.99,
             fontsize=14, fontweight="bold", color=INK)
fig.text(0.012, 0.925, "NIST Chemistry WebBook, IAPWS-95 formulation: 1 K steps from the triple "
         "point to the critical point.", color=INK2, fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(HERE, "water_svp.png"), facecolor=SURF)
for t in (0.01, 20, 25, 37, 50, 100, 150, 200, 300, 373.946):
    print(f"{t:7.2f} °C  {np.exp(np.interp(t, T, np.log(P))):.5g} kPa")
