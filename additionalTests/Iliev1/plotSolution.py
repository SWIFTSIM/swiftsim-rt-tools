#!/usr/bin/env python3

# --------------------------------------------
# Plot temperature and neutral fractions.
# --------------------------------------------

import sys

import matplotlib as mpl

mpl.use("Agg")  # use this in non-interactive environments
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy import stats

import stromgren_plotting_tools as spt

# set here which reference you want to use
# (use string!)

#  ref = "10Myr"
#  ref = "30Myr"
#  ref = "100Myr"
#  ref = "200Myr"
ref = "500Myr"

# plot individual particle data?
plot_particles = False
# if True, use only a few references, and add labels.
# Otherwise use all references, and make them grey
label_refs = False


#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
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
    "text.usetex": True,
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


scatterplot_kwargs = {"alpha": 0.1, "s": 2, "marker": ".", "linewidth": 0.0}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
# WARNING: Make sure to set the correct reference time above
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
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    # This is the original test setup
    boxsize_ref = 6.6 * unyt.kpc

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / boxsize_ref

    # Get mass fractions
    imf = spt.get_imf(scheme, data)
    xH = imf.HI + imf.HII
    xHI = imf.HI / xH
    xHII = imf.HII / xH

    # get profiles
    # max r should be sqrt(3) * boxlen
    nbins = 100
    r_bin_edges = np.linspace(0.0, 1.4, nbins + 1)
    r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
    xHI_binned, _, _ = stats.binned_statistic(
        r, xHI, statistic="mean", bins=r_bin_edges, range=(0.0, 1.4)
    )
    xHI_std, _, _ = stats.binned_statistic(
        r, xHI, statistic="std", bins=r_bin_edges, range=(0.0, 1.4)
    )
    xHII_binned, _, _ = stats.binned_statistic(
        r, xHII, statistic="mean", bins=r_bin_edges, range=(0.0, 1.4)
    )
    xHII_std, _, _ = stats.binned_statistic(
        r, xHII, statistic="std", bins=r_bin_edges, range=(0.0, 1.4)
    )

    fig = plt.figure(figsize=(5.5, 5.5), dpi=200)
    ax1 = fig.add_subplot(111)

    # First the references
    if label_refs:
        codes = ["C2Ray", "RSPH", "OTVET"]
        for c, code in enumerate(codes):
            fname_xHI = "reference/" + code + "_" + ref + "_profile_xHI.dat"

            xHI_ref, xHI_std_ref = np.loadtxt(fname_xHI, delimiter=",", unpack=True)
            r_ref = np.linspace(0.0, 1.0, xHI_ref.shape[0])
            # shift to bin center
            r_ref += (r_ref[1] - r_ref[0]) * 0.5

            ax1.errorbar(
                r_ref,
                xHI_ref,
                yerr=xHI_std_ref,
                label=code,
                alpha=0.6,
                capsize=2,
                zorder=1 + c,
            )

    else:
        codes = [
            "ART",
            "C2Ray",
            "Crash",
            "FFTE",
            "Flash",
            "IFT",
            "OTVET",
            "RSPH",
            "Zeus",
        ]
        for c, code in enumerate(codes):
            fname_xHI = "reference/" + code + "_" + ref + "_profile_xHI.dat"

            xHI_ref, xHI_std_ref = np.loadtxt(fname_xHI, delimiter=",", unpack=True)
            r_ref = np.linspace(0.0, 1.0, xHI_ref.shape[0])
            # shift to bin center
            r_ref += (r_ref[1] - r_ref[0]) * 0.5

            ax1.errorbar(
                r_ref,
                xHI_ref,
                yerr=xHI_std_ref,
                c="grey",
                alpha=0.6,
                capsize=2,
                zorder=1 + c,
            )

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

    ax1.errorbar(
        r_bin_centers,
        xHI_binned,
        yerr=xHI_std,
        label=r"GEARRT $x_{\mathrm{HI}}$",
        capsize=2,
        alpha=0.6,
        zorder=20,
    )
    ax1.errorbar(
        r_bin_centers,
        xHII_binned,
        yerr=xHII_std,
        label=r"GEARRT $x_{\mathrm{HII}}$",
        capsize=2,
        alpha=0.6,
        zorder=20,
    )

    ax1.set_xlabel("r / L")
    # note: L is the box size of the original test, not the actual run
    # with GEARRT
    ax1.set_xlim(0.0, 1.0)
    ax1.set_yscale("log")
    ax1.legend()
    ax1.grid()

    ax1.set_ylabel("Hydrogen Fractions [1]")

    fig.suptitle("Iliev+06 Test 1, $t$ = {0:.0f}".format(meta.time.to("Myr")))
    plt.tight_layout()
    figname = filename[:-5]
    figname += ".png"
    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_solution(f)
