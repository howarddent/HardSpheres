"""The stepped Lennard-Jones pair potential (Sticky/ustickysim.pas) in N2O units.
sigma and epsilon come from the same corresponding-states mapping as the P-V plot:
eps/k = Tc / Tc*, sigma^3 = KV / N_A."""
import os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TC, TCS = 309.52, 1.035           # K; model Tc*
KV = 0.030154                     # L/mol per sigma^3 (plot.py)
NA, KB = 6.02214076e23, 1.380649e-23
EPSK = TC / TCS                   # K
EPS = EPSK * KB * NA / 1000       # kJ/mol
SIG = (KV * 1e-3 / NA) ** (1 / 3) * 1e9   # nm

edges = np.array([1.0, 1.06, 1.12, 1.2, 1.3, 1.4, 1.5, 1.65, 1.8, 2.0])
F = lambda r: 4 * (-r ** -9 / 9 + r ** -3 / 3)
U = (F(edges[1:]) - F(edges[:-1])) / ((edges[1:] ** 3 - edges[:-1] ** 3) / 3)

INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
BLUE, CORE = "#2a78d6", "#efeeea"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
fig, ax = plt.subplots(figsize=(10, 6.2), dpi=160)
fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(True, color=GRID, lw=0.8)

XMAX, YMIN, YMAX = 2.4 * SIG, -3.0, 2.0
ax.axvspan(0.75 * SIG, SIG, color=CORE, lw=0, zorder=0)
ax.axhline(0, color=MUTED, lw=0.8)

r = np.linspace(0.97, 2.4, 600)
ax.plot(r * SIG, EPS * 4 * (r ** -12 - r ** -6), color=MUTED, lw=1.4, ls=(0, (4, 3)), zorder=2,
        label="Full Lennard-Jones, for comparison")
# staircase: wall at sigma, the shells, zero beyond the cutoff
xs = [SIG, SIG] + sum([[a * SIG, b * SIG] for a, b in zip(edges[:-1], edges[1:])], []) + [2 * SIG, XMAX]
ys = [YMAX + 1, EPS * U[0]] + sum([[EPS * u, EPS * u] for u in U], [])[1:] + [0, 0]
ax.plot(xs[:len(ys)], ys, color=BLUE, lw=2.2, drawstyle="default", zorder=4,
        label="Simulated potential: hard core + 9 flat shells")
# step verticals
xx, yy = [SIG, SIG], [YMAX + 1, EPS * U[0]]
for k, u in enumerate(U):
    xx += [edges[k] * SIG, edges[k + 1] * SIG]; yy += [EPS * u, EPS * u]
xx += [2 * SIG, XMAX]; yy += [0, 0]
ax.lines[-1].remove()
ax.plot(xx, yy, color=BLUE, lw=2.2, zorder=4, label="Simulated potential: hard core + 9 flat shells")

kmin = np.argmin(U)
ax.annotate(f"well floor {EPS * U[kmin]:.2f} kJ/mol\n({-U[kmin]:.2f} ε)",
            ((edges[kmin] + edges[kmin + 1]) / 2 * SIG, EPS * U[kmin]), xytext=(30, -38),
            textcoords="offset points", fontsize=9.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax.text(0.875 * SIG, -2.2, "hard core\n(overlap\nforbidden)", ha="center", va="center", color=INK2, fontsize=9.5)
ax.annotate(f"σ = {SIG:.3f} nm", (SIG, 1.2), xytext=(8, 0), textcoords="offset points", va="center",
            color=INK, fontsize=10, fontweight="bold")
ax.annotate(f"cutoff 2σ = {2 * SIG:.3f} nm", (2 * SIG, 0), xytext=(4, 14), textcoords="offset points",
            ha="left", color=INK, fontsize=9.5)
ax.annotate(f"bond edge 1.5σ", (1.5 * SIG, EPS * U[6]), xytext=(6, -22), textcoords="offset points",
            ha="left", color=INK2, fontsize=9)
ax.axvline(1.5 * SIG, color=MUTED, lw=0.7, ls=":", zorder=1)

ax.set_xlim(0.75 * SIG, XMAX); ax.set_ylim(YMIN, YMAX)
ax.set_xlabel("Distance between sphere centres (nm)")
ax.set_ylabel("Pair potential energy (kJ/mol)")
top = ax.secondary_xaxis("top", functions=(lambda x: x / SIG, lambda s: s * SIG))
top.set_xlabel("r / σ", color=INK2); top.tick_params(colors=INK2)
right = ax.secondary_yaxis("right", functions=(lambda u: u / EPS, lambda e: e * EPS))
right.set_ylabel("U / ε", color=INK2); right.tick_params(colors=INK2)
for a in (top, right): a.spines[:].set_color(MUTED)
fig.text(0.012, 0.975, "Pair potential of the N₂O model spheres", color=INK, fontsize=14,
         fontweight="bold", va="top")
fig.text(0.012, 0.935, f"Stepped Lennard-Jones, scaled to N₂O by its critical point: ε/k = {EPSK:.0f} K "
        f"(ε = {EPS:.2f} kJ/mol), σ = {SIG * 10:.2f} Å", color=INK2, fontsize=9.5, va="top")
leg = ax.legend(loc="lower right", fontsize=9, frameon=True)
leg.get_frame().set_facecolor(SURF); leg.get_frame().set_edgecolor(GRID)
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "n2o_potential.png"), facecolor=SURF)
print(f"eps/k {EPSK:.1f} K, eps {EPS:.3f} kJ/mol, sigma {SIG:.4f} nm")
for a, b, u in zip(edges[:-1], edges[1:], U):
    print(f"  {a:.2f}-{b:.2f} sigma ({a*SIG:.3f}-{b*SIG:.3f} nm): U {u:+.4f} eps = {u*EPS:+.3f} kJ/mol")
