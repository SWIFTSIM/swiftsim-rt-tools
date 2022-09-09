#!/usr/bin/env python3

# ---------------------------------------------------
# Plots temperature and mass fractions without any
# references.
# ---------------------------------------------------

import swiftsimio
import matplotlib as mpl

mpl.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import sys
import stromgren_plotting_tools as spt
import unyt
from scipy import stats

# plot individual particle data?
plot_particles = True


#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "axes.linewidth": 1.5,
    #  "text.usetex": True,
    "figure.subplot.left": 0.045,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)


scatterplot_kwargs = {"alpha": 0.1, "s": 1, "marker": ".", "linewidth": 0.0}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1

snapshot_base = "output"


def plot_solution(filename):

    # Read in data first
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    # This is the original test setup
    boxsize_ref = 6.6 * unyt.kpc

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    ids = data.gas.particle_ids
    mask = ids < 1000000000
    xpart = xpart[mask]
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / boxsize_ref

    # Get mass fractions
    imf = spt.get_imf(scheme, data)
    imf.HI = imf.HI[mask]
    imf.HII = imf.HII[mask]
    imf.HeI = imf.HeI[mask]
    imf.HeII = imf.HeII[mask]
    imf.HeIII = imf.HeIII[mask]
    xH = imf.HI + imf.HII
    xHI = imf.HI / xH
    xHII = imf.HII / xH

    # get profiles
    # max r should be sqrt(3) * boxlen
    nbins = 100
    r_bin_edges = np.linspace(0.0, 1.1, nbins + 1)
    r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
    xHI_binned, _, _ = stats.binned_statistic(
        r, xHI, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHI_std, _, _ = stats.binned_statistic(
        r, xHI, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHII_binned, _, _ = stats.binned_statistic(
        r, xHII, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHII_std, _, _ = stats.binned_statistic(
        r, xHII, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )

    fig = plt.figure(figsize=(5.5, 5.5))
    ax1 = fig.add_subplot(111)

    if plot_particles:
        ax1.scatter(r, xHI, **scatterplot_kwargs, zorder=0)
        ax1.scatter(r, xHII, **scatterplot_kwargs, zorder=0)
        # black background lines
        ax1.semilogy(
            r_bin_centers,
            xHI_binned,
            c="k",
            linewidth=params["lines.linewidth"] * 2.0,
            zorder=1,
        )
        ax1.semilogy(
            r_bin_centers,
            xHII_binned,
            c="k",
            linewidth=params["lines.linewidth"] * 2.0,
            zorder=1,
        )
    ax1.semilogy(r_bin_centers, xHI_binned, label=r"$x_{\mathrm{HI}}$", zorder=2)
    ax1.semilogy(r_bin_centers, xHII_binned, label=r"$x_{\mathrm{HII}}$", zorder=2)
    #  ax1.errorbar(r_bin_centers, xHI_binned, yerr=xHI_std, label=r"GEARRT $x_{\mathrm{HI}}$", capsize=2,  zorder=20)
    #  ax1.errorbar(r_bin_centers, xHII_binned, yerr=xHII_std, label=r"GEARRT $x_{\mathrm{HII}}$", capsize=2, zorder=20)

    for ax in fig.axes:
        ax.set_xlabel("r / L")
        # note: L is the box size of the original test, not the actual run
        # with GEARRT
        ax.set_xlim(0.0, 1.01)
        ax.set_yscale("log")

    ax1.legend(loc="lower right")
    ax1.set_ylabel("Hydrogen Fractions")
    ax1.set_ylim(5e-5, 1.2)

    fig.suptitle("Iliev+06 Test 1, $t$ = {0:.0f}".format(meta.time.to("Myr")))
    plt.tight_layout()
    figname = filename[:-5]
    figname += "-NoRef.png"
    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_solution(f)
