from pathlib import Path
import sys

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

mjy = 1.0e-26
kpc = 3.0e21

try:
    folder = Path(sys.argv[1])
except IndexError:
    folder = Path("Output")

particles = np.genfromtxt(folder / "Singlezone_Particles.dat")
syn = np.genfromtxt(folder / "Singlezone_Syn.dat")
ssc = np.genfromtxt(folder / "Singlezone_SSC.dat")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6))

ax1.loglog(syn.T[0], syn.T[0] * syn.T[1], linewidth=2.5, color=colors[0])
ax1.loglog(ssc.T[0], ssc.T[0] * ssc.T[1], linewidth=2.5, color=colors[1])
ax1.set_ylim(1e37, 5e41)
ax1.set_xlim(1e8, 2e28)
ax1.set_xlabel(r"Frequency ($\rm{Hz}$)", fontsize=24)
ax1.set_ylabel(r"Luminosity ($\rm{erg\,s^{-1}}$)", fontsize=24)

ax2.loglog(
    particles.T[1], particles.T[3] * particles.T[1], color=colors[2], linewidth=2.5
)
ax2.set_ylabel(r"Number density ($\rm{\#/cm^{-3}}$)", fontsize=24)
ax2.set_xlabel(r"Particle Lorenz factor", fontsize=24)

fig.tight_layout()

plt.savefig(folder / "singlezone_examples.pdf")
# plt.show()
