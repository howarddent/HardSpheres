# This file is part of HardSpheres.
# Copyright (C) 2026 Howard Dent
# SPDX-License-Identifier: GPL-3.0-or-later
#
# HardSpheres is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version. It is distributed WITHOUT ANY WARRANTY; see
# the file LICENSE at the top of the repository for details.

"""Saturated vapour pressure of water, 0-100 °C, linear axes (NIST, IAPWS-95)."""
import csv, os, numpy as np, matplotlib
matplotlib.use("Agg")
HERE = os.path.dirname(os.path.abspath(__file__))
import matplotlib.pyplot as plt

rows = list(csv.reader(open(os.path.join(HERE, "water_sat_T.tsv")), delimiter="\t"))[1:]
T = np.array([float(r[0]) for r in rows]) - 273.15
P = np.array([float(r[1]) for r in rows]) * 100    # bar -> kPa
keep = T <= 100.5
T, P = T[keep], P[keep]
psat = lambda t: float(np.exp(np.interp(t, T, np.log(P))))
# fill the last half-kelvin so the curve ends exactly at 100 °C
p100 = psat(100.0)
below = T < 100
T, P = np.append(T[below], 100.0), np.append(P[below], p100)

INK, INK2, MUTED, GRID, SURF, BLUE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(True, color=GRID, lw=0.8)
ax.plot(T, P, color=BLUE, lw=2.2, zorder=3)
ax.set_xlim(0, 100); ax.set_ylim(0, 110)
ax.set_xticks(range(0, 101, 10))
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Saturated vapour pressure (kPa)")
ax.axhline(101.325, color=MUTED, lw=1, ls=(0, (4, 3)), zorder=2)
ax.text(2, 101.325, "1 standard atmosphere, 101.3 kPa", va="bottom", ha="left", fontsize=9, color=INK2)
for t, txt, off in [(20, None, (-8, 8)), (37.4, "body temperature", (-8, 8)), (50, None, (8, -6)),
                    (80, None, (-8, 8)), (100, None, (-8, -14))]:
    p = psat(t)
    ax.plot([t], [p], "o", ms=6.5, color=BLUE, mec=SURF, mew=1.5, zorder=4)
    label = f"{t:g} °C: {p:.3g} kPa" if p < 100 else f"{t:g} °C: {p:.1f} kPa"
    if txt: label = f"{txt}\n{label}"
    ax.annotate(label, (t, p), xytext=off, textcoords="offset points", ha="right" if off[0] < 0 else "left",
                va="bottom" if off[1] > 0 else "top", fontsize=9, color=INK)
fig.text(0.012, 0.975, "Saturated vapour pressure of water, 0–100 °C", fontsize=14, fontweight="bold",
         color=INK, va="top")
fig.text(0.012, 0.93, "NIST Chemistry WebBook, IAPWS-95 formulation, 1 K steps.", color=INK2,
         fontsize=9.5, va="top")
fig.tight_layout(rect=(0, 0, 1, 0.91))
fig.savefig(os.path.join(HERE, "water_svp_0-100.png"), facecolor=SURF)
