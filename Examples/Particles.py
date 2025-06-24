import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc("text", usetex=True)
rc("font", **{"family": "serif", "serif": ["DejaVu Serif Display"]})
plt.rcParams.update({"font.size": 20})

colors = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

try:
    folder = Path(sys.argv[1])
except IndexError:
    folder = Path("Output")

lowfrac = np.genfromtxt(folder / "lowfrac.dat")
lowfrac_cool = np.genfromtxt(folder / "lowfrac_cool.dat")
highfrac = np.genfromtxt(folder / "highfrac.dat")
highfrac_cool = np.genfromtxt(folder / "highfrac_cool.dat")
kdist = np.genfromtxt(folder / "kdist.dat")
kdist_cool = np.genfromtxt(folder / "kdist_cool.dat")
bkndist = np.genfromtxt(folder / "bkndist.dat")
bkndist_cool = np.genfromtxt(folder / "bkndist_cool.dat")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6))
ax1.loglog(
    lowfrac.T[0],
    lowfrac.T[2] * lowfrac.T[0],
    color=colors[0],
    linewidth=2.5,
    label="Mixed $f_{\\rm pl}=0.1$",
)
ax1.loglog(
    lowfrac_cool.T[0],
    lowfrac_cool.T[2] * lowfrac_cool.T[0],
    color=colors[0],
    linewidth=2.5,
    linestyle="dashed",
)
ax1.loglog(
    highfrac.T[0],
    0.03 * highfrac.T[2] * highfrac.T[0],
    color=colors[1],
    linewidth=2.5,
    label="Mixed $f_{\\rm pl}=0.9$",
)
ax1.loglog(
    highfrac_cool.T[0],
    0.03 * highfrac_cool.T[2] * highfrac_cool.T[0],
    color=colors[1],
    linewidth=2.5,
    linestyle="dashed",
)
ax1.loglog(
    kdist.T[0],
    0.001 * kdist.T[2] * kdist.T[0],
    color=colors[2],
    linewidth=2.5,
    label="$\\kappa$",
)
ax1.loglog(
    kdist_cool.T[0],
    0.001 * kdist_cool.T[2] * kdist_cool.T[0],
    color=colors[2],
    linewidth=2.5,
    linestyle="dashed",
)
ax1.loglog(
    bkndist.T[0],
    0.0001 * bkndist.T[2] * bkndist.T[0],
    color=colors[3],
    linewidth=2.5,
    label="Broken PL",
)
ax1.loglog(
    bkndist_cool.T[0],
    0.0001 * bkndist_cool.T[2] * bkndist_cool.T[0],
    color=colors[3],
    linewidth=2.5,
    linestyle="dashed",
)
ax1.set_ylim([1e-8, 2.0])
ax1.set_xlabel("Particle momentum", fontsize=24)
ax1.set_ylabel("Number density (Arb. units)", fontsize=24)

ax2.loglog(
    lowfrac.T[1],
    lowfrac.T[3] * lowfrac.T[1],
    color=colors[0],
    linewidth=2.5,
    label="$f_{\\rm nth}=0.1$",
)
ax2.loglog(
    lowfrac_cool.T[1],
    lowfrac_cool.T[3] * lowfrac_cool.T[1],
    color=colors[0],
    linewidth=2.5,
    linestyle="dashed",
)
ax2.loglog(
    highfrac.T[1],
    highfrac.T[3] * highfrac.T[1],
    color=colors[1],
    linewidth=2.5,
    label="$f_{\\rm nth}=0.9$",
)
ax2.loglog(
    highfrac_cool.T[1],
    highfrac_cool.T[3] * highfrac_cool.T[1],
    color=colors[1],
    linewidth=2.5,
    linestyle="dashed",
)
ax2.loglog(
    kdist.T[1],
    kdist.T[3] * kdist.T[1],
    color=colors[2],
    linewidth=2.5,
    label="$\\kappa$",
)
ax2.loglog(
    kdist_cool.T[1],
    kdist_cool.T[3] * kdist_cool.T[1],
    color=colors[2],
    linewidth=2.5,
    linestyle="dashed",
)
ax2.loglog(
    bkndist.T[1],
    bkndist.T[3] * bkndist.T[1],
    color=colors[3],
    linewidth=2.5,
    label="Broken PL",
)
ax2.loglog(
    bkndist_cool.T[1],
    bkndist_cool.T[3] * bkndist_cool.T[1],
    color=colors[3],
    linewidth=2.5,
    linestyle="dashed",
)
ax2.set_ylim([4e-5, 2.0])
ax2.set_xlabel("Particle Lorenz factor", fontsize=24)
ax2.yaxis.tick_right()
ax2.legend(loc="lower left", fontsize=24)
ax2.yaxis.set_label_position("right")

fig.tight_layout()
fig.subplots_adjust(wspace=0)

plt.savefig(folder / "particle_examples.pdf")
# plt.show()
