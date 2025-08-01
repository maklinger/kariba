import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import math

from matplotlib import rc

rc("text", usetex=True)
rc("font", **{"family": "serif", "serif": ["DejaVu Serif Display"]})
plt.rcParams.update({"font.size": 20})

plot_colors = [
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

pars = np.genfromtxt("Input/ip.dat", comments="//")
pars[26] = 1

kevconv = 1.0
# kevconv = 2.41*10**17

# styles = [(0,(1,6)),(0,(2,5)),(0,(3,4)),(0,(4,3)),(0,(5,2)),(0,(6,1))]

styles = [
    (0, (1, 5)),
    (0, (1, 1)),
    (0, (1, 1)),
    (0, (5, 10)),
    (0, (5, 5)),
    (0, (5, 1)),
    (0, (3, 10, 1, 10)),
    (0, (3, 5, 1, 5)),
    (0, (3, 1, 1, 1)),
    (0, (3, 5, 1, 5, 1, 5)),
    (0, (3, 10, 1, 10, 1, 10)),
    (0, (3, 1, 1, 1, 1, 1)),
]

# fluxconv = 1.
fluxconv = 4.0 * math.pi * (pars[2] * 3.0 * 10**21) ** 2

if pars[0] >= 1.0e3:
    ulim_f = 5.0e25
    blim_f = 1e9
    ulim_f2 = 1.0e13
    blim_f2 = 0.5e9
    ulim_fl = 8.0e-10
    blim_fl = 3.0e-14
    ulim_fd = 5.0e2
    blim_fd = 0.3
else:
    ulim_f = 0.9e21
    blim_f = 0.5e9
    ulim_f2 = 1.0e13
    blim_f2 = 0.5e9
    ulim_fl = 1e-7
    blim_fl = 1.0e-16
    ulim_fd = 5.0e3
    blim_fd = 0.2e-3

zmin = 2.0
zmax = pars[8]

zheight = np.zeros(5)

zheight = np.zeros(5)  # the array that we use to print the labels in colorbar
for n in range(0, len(zheight)):
    zheight[n] = np.log10(
        int((pow(10.0, math.log10(zmax / zmin) * (n * 0.25) + math.log10(zmin))))
    )

zheight = np.around(zheight, decimals=2)

z = pars[3]
plotcheck = pars[26]

if plotcheck >= 1:
    Presyn = np.genfromtxt("Output/Presyn.dat")
    Postsyn = np.genfromtxt("Output/Postsyn.dat")
    Precom = np.genfromtxt("Output/Precom.dat")
    Postcom = np.genfromtxt("Output/Postcom.dat")
    Disk = np.genfromtxt("Output/Disk.dat")
    BB = np.genfromtxt("Output/BB.dat")
    Total = np.genfromtxt("Output/Total.dat")
else:
    exit()

mjy = 1.0e-26

i = 0
j = 0
nzones = 100

colors = pl.cm.magma(np.linspace(0.0, 0.9, nzones))

size_cyclo_arr = np.zeros(nzones)
size_com_arr = np.zeros(nzones)

if plotcheck >= 2:
    Cyclosyn_zones = np.genfromtxt("Output/Cyclosyn_zones.dat")
    Compton_zones = np.genfromtxt("Output/Compton_zones.dat")
    Numdens = np.genfromtxt("Output/Numdens.dat")
    size3 = len(Numdens.T[0]) // nzones
    i = 0
    j = 0
    k = 0
    for j in range(len(Cyclosyn_zones.T[0]) - 1):
        if Cyclosyn_zones.T[0][j + 1] > Cyclosyn_zones.T[0][j]:
            k = k + 1
        else:
            size_cyclo_arr[i] = k + 1
            i = i + 1
            k = 0
    i = 0
    j = 0
    k = 0
    for j in range(len(Compton_zones.T[0]) - 1):
        if Compton_zones.T[0][j + 1] > Compton_zones.T[0][j]:
            k = k + 1
        else:
            size_com_arr[i] = k + 1
            i = i + 1
            k = 0
    garr = np.zeros(size3)
    ngarr = np.zeros(size3)
    parr = np.zeros(size3)
    nparr = np.zeros(size3)

totindex1 = 0
totindex2 = 0
fig1, (ax1) = plt.subplots(1, 1, figsize=(7.5, 6))
for i in range(nzones - 1):
    nu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))
    lnu_cyclosyn = np.zeros(int(size_cyclo_arr[i]))
    nu_compton = np.zeros(int(size_com_arr[i]))
    lnu_compton = np.zeros(int(size_com_arr[i]))
    for j in range(int(size_cyclo_arr[i])):
        nu_cyclosyn[j] = Cyclosyn_zones.T[0][totindex1 + j] / kevconv
        lnu_cyclosyn[j] = (
            Cyclosyn_zones.T[1][totindex1 + j] * nu_cyclosyn[j] * mjy * kevconv
        )
    for k in range(int(size_com_arr[i])):
        nu_compton[k] = Compton_zones.T[0][totindex2 + k] / kevconv
        lnu_compton[k] = (
            Compton_zones.T[1][totindex2 + k] * nu_compton[k] * mjy * kevconv
        )
    totindex1 = totindex1 + int(size_cyclo_arr[i])
    totindex2 = totindex2 + int(size_com_arr[i])
    if (i % 15) == 0 and i != 0:
        number = int(i / 15 - 1)
        ax1.plot(
            nu_cyclosyn,
            lnu_cyclosyn * fluxconv,
            linewidth=2.0,
            color=colors[i],
            zorder=nzones - i,
            linestyle=styles[number],
        )  # ,linestyle='dashed')
        ax1.plot(
            nu_compton,
            lnu_compton * fluxconv,
            linewidth=2.0,
            color=colors[i],
            zorder=nzones - i,
            linestyle=styles[number],
        )  # ,linestyle='dashed')
# ax1.plot(Postsyn.T[0]/kevconv,Postsyn.T[1]*Postsyn.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[0],zorder=nzones+1)
# ax1.plot(Postcom.T[0]/kevconv,Postcom.T[1]*Postcom.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[1],zorder=nzones+1)
# ax1.plot(Disk.T[0]/kevconv,Disk.T[1]*Disk.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[2],zorder=nzones+1)
# ax1.plot(BB.T[0]/kevconv,BB.T[1]*BB.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[3],zorder=nzones+1)
# ax1.plot(Presyn.T[0]/kevconv,Presyn.T[1]*Presyn.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[4],zorder=nzones+1)
# ax1.plot(Precom.T[0]/kevconv,Precom.T[1]*Precom.T[0]*mjy*fluxconv,linewidth=2.5,color=plot_colors[5],zorder=nzones+1)
ax1.plot(
    Total.T[0] / kevconv,
    Total.T[1] * Total.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color="black",
    zorder=nzones + 1,
)
ax1.set_ylim([0.001 * blim_fl * fluxconv, 0.1 * ulim_fl * fluxconv])
ax1.set_xlim([blim_f / kevconv, ulim_f / kevconv])
ax1.set_xscale("log", base=10)
ax1.set_yscale("log", base=10)
if kevconv != 1:
    ax1.set_xlabel("Energy (kev)", fontsize=18)
else:
    ax1.set_xlabel("Frequency (Hz)", fontsize=18)
if fluxconv != 1:
    ax1.set_ylabel("Luminosity (erg/s)", fontsize=18)
else:
    ax1.set_ylabel("Flux (erg/s/cm^2)", fontsize=18)

fig2, (ax2) = plt.subplots(1, 1, figsize=(7.5, 6))
# ax2.plot(Total.T[0]/kevconv,Total.T[1]*Total.T[0]*mjy*fluxconv,linewidth=1.5,color='black',zorder=nzones+1,label='Total')
ax2.plot(
    Postsyn.T[0] / kevconv,
    Postsyn.T[1] * Postsyn.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[0],
    zorder=nzones + 1,
    linestyle=(0, (5, 1)),
    label="Syn, $z\\geq z_{\\rm diss}$",
)
ax2.plot(
    Postcom.T[0] / kevconv,
    Postcom.T[1] * Postcom.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[1],
    zorder=nzones + 1,
    linestyle=(0, (5, 3)),
    label="IC, $z\\geq z_{\\rm diss}$",
)
ax2.plot(
    Disk.T[0] / kevconv,
    Disk.T[1] * Disk.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[2],
    zorder=nzones + 1,
    linestyle=(0, (3, 1, 1, 1)),
    label="Disk",
)
ax2.plot(
    BB.T[0] / kevconv,
    BB.T[1] * BB.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[3],
    zorder=nzones + 1,
    linestyle=(0, (3, 3, 1, 3)),
    label="Torus",
)
ax2.plot(
    Presyn.T[0] / kevconv,
    Presyn.T[1] * Presyn.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[4],
    zorder=nzones + 1,
    linestyle=(0, (3, 2, 1, 2, 1, 2)),
    label="Syn, $z<z_{\\rm diss}$",
)
ax2.plot(
    Precom.T[0] / kevconv,
    Precom.T[1] * Precom.T[0] * mjy * fluxconv,
    linewidth=2.5,
    color=plot_colors[5],
    zorder=nzones + 1,
    linestyle=(0, (3, 1, 1, 1, 1, 1)),
    label="IC, $z<z_{\\rm diss}$",
)
ax2.set_ylim([blim_fl * fluxconv, ulim_fl * fluxconv])
ax2.set_xlim([blim_f / kevconv, ulim_f / kevconv])
ax2.legend(loc="upper left", fontsize=18, ncol=2)
ax2.set_xscale("log", base=10)
ax2.set_yscale("log", base=10)
if kevconv != 1:
    ax2.set_xlabel("Energy (kev)", fontsize=18)
else:
    ax2.set_xlabel("Frequency (Hz)", fontsize=18)
if fluxconv != 1:
    ax2.set_ylabel("Luminosity (erg/s)", fontsize=18)
else:
    ax2.set_ylabel("Flux (erg/s/cm^2)", fontsize=18)

i = 0
j = 0
if plotcheck >= 2:
    plt.figure(2)
    fig4, (ax4) = plt.subplots(1, 1, figsize=(9.0, 6))
    for j in range(nzones - 1):
        for i in range(size3):
            parr[i] = Numdens.T[0][i + size3 * j]
            nparr[i] = Numdens.T[2][i + size3 * j]
        if (j % 3) == 0:  # this plots one in five zones for clarity
            ax4.plot(parr, nparr * parr, linewidth=3.5, color=colors[j], zorder=150 - j)
    ax4.set_xlabel("Momentum (erg*s/cm)", fontsize=18)
    ax4.set_ylabel("Number density ($\\#/\\rm{cm}^3$)", fontsize=18)
    ax4.tick_params(axis="both", which="major", labelsize=20)
    ax4.set_xscale("log", base=10)
    ax4.set_yscale("log", base=10)

    fig4.subplots_adjust(wspace=0)
    # We create the colorbar:
    sm2 = plt.cm.ScalarMappable(cmap="magma")
    sm2.set_array(colors)
    bar2 = fig4.colorbar(sm2, pad=0.10, ax=[ax4], ticks=[0.01, 0.25, 0.5, 0.75, 0.99])
    bar2.set_label(r"$log_{10}(z)\,(r_g$)", fontsize=18)
    bar2.ax.set_yticklabels(zheight, fontsize=15)

    i = 0
    j = 0
    plt.figure(3)
    fig5, (ax5) = plt.subplots(1, 1, figsize=(7.5, 6))
    for j in range(nzones - 1):
        for i in range(size3):
            garr[i] = Numdens.T[1][i + size3 * j]
            ngarr[i] = Numdens.T[3][i + size3 * j]
        if (j % 3) == 0:  # this plots one in five zones for clarity
            ax5.plot(garr, ngarr * garr, linewidth=3.5, color=colors[j], zorder=150 - j)
    ax5.set_xlabel("Electron $\\gamma$", fontsize=18)
    ax5.set_ylabel("Number density ($\\#/cm^3$)", fontsize=18)
    ax5.set_xscale("log", base=10)
    ax5.set_yscale("log", base=10)
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")

plt.tight_layout()
plt.savefig("separate.pdf")
# plt.show()
