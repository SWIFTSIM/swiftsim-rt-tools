#!/usr/bin/env python3

import swiftsimio
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys
import stromgren_plotting_tools as spt
import unyt
from scipy import stats

#  mpl.use('Gtk')


# Plot parameters
#  params = {
#      "axes.labelsize": 14,
#      "axes.titlesize": 14,
#      "font.size": 14,
#      "legend.fontsize": 14,
#      "xtick.labelsize": 12,
#      "ytick.labelsize": 12,
#      "xtick.direction": "in",
#      "ytick.direction": "in",
#      "xtick.top": True,
#      "ytick.right": True,
#      "xtick.major.width": 1.5,
#      "ytick.major.width": 1.5,
#      "axes.linewidth": 1.5,
#      "text.usetex": True,
#      "figure.figsize": (5, 4),
#      "figure.subplot.left": 0.045,
#      "figure.subplot.right": 0.99,
#      "figure.subplot.bottom": 0.05,
#      "figure.subplot.top": 0.99,
#      "figure.subplot.wspace": 0.15,
#      "figure.subplot.hspace": 0.12,
#      "lines.markersize": 1,
#      "lines.linewidth": 2.0,
#  }
#  mpl.rcParams.update(params)


scatterplot_kwargs = {
    "alpha": 0.6,
    "s": 2,
    "marker": "o",
    "linewidth": 0.0,
}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
# WARNING: The reference solution is comparable with snapshot_500 only
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1


plot_SPHM1RT = False
plot_Helium = False

snapshot_base = "output"
if plot_SPHM1RT:
    snapshot_base = "SPHM1RT_output"


def plot_solution(filename):

    # Read in data first
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / (0.5 * boxsize[0])

    # Get mass fractions
    imf = spt.get_imf(scheme, data)
    xH = imf.HI + imf.HII
    xHI = imf.HI / xH
    xHII = imf.HII / xH

    # get temperature
    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    T = spt.gas_temperature(data.gas.internal_energies, mu, meta.gas_gamma)

    # get profiles
    # max r should be sqrt(3) * boxlen
    #  nbins = 100
    #  r_bin_edges = np.linspace(0., 2., nbins + 1)
    #  r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
    #  xHI_binned, _, _ = stats.binned_statistic(r, xHI, statistic="mean", bins=r_bin_edges, range=(0., 1.))
    #  xHII_binned, _, _ = stats.binned_statistic(r, xHII, statistic="mean", bins=r_bin_edges, range=(0., 1.))
    #  T_binned, _, _ = stats.binned_statistic(r, T, statistic="mean", bins=r_bin_edges, range=(0., 1.))

    fig = plt.figure(figsize=(10, 5.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.scatter(r, imf.HI, **scatterplot_kwargs, label="HI")
    ax1.scatter(r, imf.HII, **scatterplot_kwargs, label="HII")
    if plot_Helium:
        ax1.scatter(r, imf.HeI, **scatterplot_kwargs, label="HeI")
        ax1.scatter(r, imf.HeII, **scatterplot_kwargs, label="HeII")
        ax1.scatter(r, imf.HeIII, **scatterplot_kwargs, label="HeIII")
    #  ax1.semilogy(r_bin_centers, xHI_binned, label=r"$x_{HI}$")
    #  ax1.semilogy(r_bin_centers, xHII_binned, label=r"$x_{HII}$")
    ax1.set_ylim(5e-7, 1.5)


    ax2.scatter(r, T, **scatterplot_kwargs)
    #  ax2.semilogy(r_bin_centers, T_binned, label=r"gas temperature")


    ax1.set_xlabel("r / (L / 2)")
    ax1.set_ylabel("Hydrogen Fractions")
    ax1.set_xlim(0., 1.4)
    ax1.set_yscale("log")
    ax1.legend()

    ax2.set_xlabel("r / (L / 2)")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_xlim(0., 1.)
    ax2.set_yscale("log")
    #  ax2.legend()

    #  plt.xlim([0, boxsize[0] / 2.0])
    plt.tight_layout()
    figname = filename[:-5]
    figname += ".png"

    if plot_SPHM1RT:
        figname = filename[8:-5]
        figname += "-SPHM1RT.png"

    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_solution(f)
