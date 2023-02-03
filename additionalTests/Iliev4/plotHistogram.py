#!/usr/bin/env python3

# ----------------------------------------------------
# Plot histogram of neutral hydrogen and temperature
# in entire box.
#
# Usage: ./plotHistogram.py <output_nr>
# ----------------------------------------------------

import sys

import matplotlib as mpl

#  mpl.use("Agg")
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
from swiftsimio.visualisation.slice import slice_gas

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "output"

# parameters for imshow plots
imshow_kwargs = {
    "origin": "lower",
    #  "cmap": rainbow4
    "cmap": "cividis",
}

# which reference to plot
# must be a string!
ref = "0.05Myr"
#  ref = "0.1Myr"
#  ref = "0.2Myr"
#  ref = "0.3Myr"
#  ref = "0.4Myr"

logTmin = 1.9
logTmax = 5.
logxHImin = -2.5
logxHImax = 0.
nbins = 100

#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
    #  "legend.fontsize": 14,
    #  "xtick.labelsize": 12,
    #  "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    #  "xtick.major.width": 1.5,
    #  "ytick.major.width": 1.5,
    "axes.linewidth": 1.5,
    "text.usetex": True,
    "figure.subplot.left": 0.10,
    "figure.subplot.right": 0.98,
    "figure.subplot.bottom": 0.125,
    "figure.subplot.top": 0.85,
    "figure.subplot.wspace": 0.18,
    "figure.subplot.hspace": 0.06,
    #  "lines.markersize": 1,
    "lines.linewidth": 3.0,
    "mpl_toolkits.legacy_colorbar": False,
}
mpl.rcParams.update(params)


# -----------------------------------------------------------------------


def get_ref_data(code, quantity):

    datafile = "reference/" + code + "_" + ref + "_histogram_" + quantity + ".dat"
    print("reading", datafile)
    centers, counts = np.loadtxt(datafile, delimiter=",", unpack=True)
    return centers, counts


def plot_result(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    # Grab your own data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    gamma = meta.gas_gamma

    # remove boundary particles
    real_parts = data.gas.particle_ids <= 1000000000

    imf = spt.get_imf(scheme, data)

    logxHI = np.log10(imf.HI[real_parts])
    logxHII = np.log10(imf.HII[real_parts])

    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    T = spt.gas_temperature(data.gas.internal_energies, mu, gamma)

    logT = np.log10(T[real_parts])

    histlogxHI, edges_xHI = np.histogram(logxHI, bins=nbins, range=(logxHImin, logxHImax))
    histlogxHII, edges_xHII = np.histogram(logxHII, bins=nbins, range=(logxHImin, logxHImax))
    histlogT, edges_T = np.histogram(logT, bins=nbins, range=(logTmin, logTmax))

    # normalize
    histlogxHI = histlogxHI / np.sum(histlogxHI)
    histlogxHII = histlogxHII / np.sum(histlogxHII)
    histlogT = histlogT / np.sum(histlogT)

    centers_xHI = 0.5*(edges_xHI[:-1] + edges_xHI[1:])
    centers_xHII = 0.5*(edges_xHII[:-1] + edges_xHII[1:])
    centers_T = 0.5*(edges_T[:-1] + edges_T[1:])



    fig = plt.figure(figsize=(10, 5), dpi=200)
    figname = filename[:-5] + "-Histogram.png"

    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    #  ax1.semilogy(centers_xHII, histlogxHII, label="GEARRT", zorder=20)
    ax1.semilogy(centers_xHI, histlogxHI, label="GEARRT", zorder=20)
    ax2.semilogy(centers_T, histlogT, label="GEARRT", zorder=20)

    references = ["C2Ray", "Crash", "FFTE"]
    for ref in references:

        label = None
        col = "grey"
        alpha = 0.6
        zorder = 1

        if ref == references[-1]:
            label = "references"
        if ref == "C2Ray":
            label = "C2Ray"
            col = "C3"
            alpha = 0.8
            zorder = 10

        centers, counts = get_ref_data(ref, "logXHI")
        counts = counts / np.sum(counts)
        ax1.semilogy(centers, counts, label=label, alpha=alpha, zorder=zorder, c=col)

        centers, counts = get_ref_data(ref, "logT")
        counts = counts / np.sum(counts)
        ax2.semilogy(centers, counts, label=label, alpha=alpha, zorder=zorder, c=col)


    #  ax1.set_title("Ionized Hydrogen Mass Fraction")
    #  ax1.set_xlabel("log(x_{\mathrm{HII}})")
    #  ax1.set_ylabel(r"$N(x_{\mathrm{HII}}) / N_{\mathrm{tot}}$")
    ax1.legend()
    ax1.set_title("Neutral Hydrogen Mass Fraction")
    ax1.set_xlabel(r"$\log(x_{\mathrm{HI}})$")
    ax1.set_ylabel(r"$N(x_{\mathrm{HI}}) / N_{\mathrm{tot}}$")

    ax2.set_title("Temperature")
    ax2.set_xlabel(r"$\log(T)$")
    ax2.set_ylabel(r"$N(T) / N_{\mathrm{tot}}$")


    title = "Iliev+06 Test 4"
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_result(snap)
