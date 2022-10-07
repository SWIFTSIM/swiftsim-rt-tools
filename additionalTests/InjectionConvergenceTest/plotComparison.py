#!/usr/bin/env python3

# -------------------------------------------------------------
# Plot temperature and neutral fractions of all resolutions.
# Plot reference solutions in background if selected to do so.
#
# This script is intended to overplot several convergence tests
# over each other. Say you run a flux injection scheme test
# with individual directories "Flux1", "Flux2", "Flux3", each
# of them having resolutions 16, 32, 64, 128.
# -------------------------------------------------------------

import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy import stats

import stromgren_plotting_tools as spt

# which resolutions to use
resolutions = [128, 64, 32, 16]

# which subdirectories to look for
subdirs = ["results_fullFlux", "results_FinterHalf", "results_noF"]
subdir_labels = ["full flux injection", "modified flux injection", "no flux injection"]

# plot vertical lines of star's smoothing length?
# NOTE: this actually plots compact support radius of star, i.e.
# the outer edge into which the radiation gets injected.
plot_smoothing_length = True
# you might need to adapt the kernel gamma value. This one
# is valid for 3D cubic spline kernel.
kernel_gamma = 1.825742

# shall I plot the reference solutions?
plot_refs = False

# set here which reference you want to use
# (use string!)
ref = "10Myr"

# if True, use only a few references, and add labels.
# Otherwise, use all references, and make them grey
label_refs = False

# plot individual particle data?
plot_particles = False

#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "legend.fontsize": 11,
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
    "figure.figsize": (5, 4),
    "figure.subplot.left": 0.045,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 1,
    "lines.linewidth": 1.0,
}
mpl.rcParams.update(params)

scatterplot_kwargs = {"alpha": 0.1, "s": 2, "marker": ".", "linewidth": 0.0}
refplotkwargs = {"alpha": 0.6, "capsize": 2}
if not label_refs:
    refplotkwargs["color"] = "grey"

# Read in cmdline arg: Are we plotting only one snapshot, or all?
# WARNING: Make sure to set the correct reference time above
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1

snapshot_base = "output"


def get_snapshot_number_list(snapshot_basename=snapshot_base):
    """
    Find the snapshot numbers that are to be plotted 
    and return the numbers as list
    """

    snapnumbers = []

    base = snapshot_basename + "-" + str(resolutions[0])

    dirlist = os.listdir(subdirs[0])
    for f in dirlist:
        if f.startswith(base) and f.endswith("hdf5"):
            snapnr = f[len(base) + 1 : -5]
            snapnumbers.append(snapnr)

    snapnumbers = sorted(snapnumbers)

    return snapnumbers


def plot_solution(snapnr):
    """
    Plot the solution for the given snapnr.
    snapnr is expected to be 4 digit zero padded string,
    e.g. 0003
    """

    fig = plt.figure(figsize=(10, 5.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Are we plotting references in the background?
    # ------------------------------------------------
    if plot_refs:
        # First the references
        if label_refs:
            codes = ["C2Ray", "RSPH", "OTVET"]
            for c, code in enumerate(codes):
                fname_xHI = "reference/" + code + "_" + ref + "_profile_xHI.dat"
                fname_T = "reference/" + code + "_" + ref + "_profile_T.dat"

                xHI_ref, xHI_std_ref = np.loadtxt(fname_xHI, delimiter=",", unpack=True)
                T_ref, T_std_ref = np.loadtxt(fname_T, delimiter=",", unpack=True)
                r_ref = np.linspace(0.0, 1.0, T_ref.shape[0])
                # shift to bin center
                r_ref += (r_ref[1] - r_ref[0]) * 0.5

                ax1.errorbar(
                    r_ref,
                    xHI_ref,
                    yerr=xHI_std_ref,
                    label=code,
                    zorder=1 + c,
                    **refplotkwargs,
                )
                ax2.errorbar(
                    r_ref,
                    T_ref,
                    yerr=T_std_ref,
                    label=code,
                    zorder=1 + c,
                    **refplotkwargs,
                )

        else:
            codes = ["ART", "Crash", "IFT", "RSPH", "C2Ray", "FFTE", "OTVET", "Zeus"]
            for c, code in enumerate(codes):
                fname_xHI = "reference/" + code + "_" + ref + "_profile_xHI.dat"
                fname_T = "reference/" + code + "_" + ref + "_profile_T.dat"

                xHI_ref, xHI_std_ref = np.loadtxt(fname_xHI, delimiter=",", unpack=True)
                T_ref, T_std_ref = np.loadtxt(fname_T, delimiter=",", unpack=True)
                r_ref = np.linspace(0.0, 1.0, T_ref.shape[0])
                # shift to bin center
                r_ref += (r_ref[1] - r_ref[0]) * 0.5

                ax1.errorbar(
                    r_ref, xHI_ref, yerr=xHI_std_ref, zorder=1 + c, **refplotkwargs
                )
                ax2.errorbar(
                    r_ref, T_ref, yerr=T_std_ref, zorder=1 + c, **refplotkwargs
                )

    # Plot actual results
    # -------------------------------

    # gather min/max values for axis limits
    xHI_min = 2.0
    xHI_max = 0.0
    T_min = 1e38
    T_max = 0.0

    linestyles = [":", "--", "-"]

    for s, subdir in enumerate(subdirs):
        for i, res in enumerate(resolutions):

            res_str = str(res)
            filename = snapshot_base + "-" + str(res_str) + "_" + snapnr + ".hdf5"
            filename = os.path.join(subdir, filename)

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
            #  hstar = data.stars.smoothing_lengths[0]
            hstar = data.stars.smoothing_lengths[0] * kernel_gamma
            hstar_L = hstar / boxsize_ref

            # Get mass fractions
            imf = spt.get_imf(scheme, data)
            xH = imf.HI + imf.HII
            xHI = imf.HI / xH
            xHII = imf.HII / xH

            # get temperature
            mu = spt.mean_molecular_weight(
                imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII
            )
            T = spt.gas_temperature(data.gas.internal_energies, mu, meta.gas_gamma)

            xHI_min = min(xHI_min, xHI.min())
            xHI_max = max(xHI_max, xHI.max())
            T_min = min(T_min, T.min())
            T_max = max(T_max, T.max())

            # get profiles
            # max r should be sqrt(3) * boxlen
            nbins = int(res * 3)
            outer_edge = 1.01
            #  outer_edge = 1.4
            r_bin_edges = np.linspace(0.0, outer_edge, nbins + 1)
            r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
            xHI_binned, _, _ = stats.binned_statistic(
                r, xHI, statistic="mean", bins=r_bin_edges, range=(0.0, outer_edge)
            )
            #  xHI_std, _, _ = stats.binned_statistic(
            #      r, xHI, statistic="std", bins=r_bin_edges, range=(0.0, outer_edge)
            #  )
            #  xHII_binned, _, _ = stats.binned_statistic(
            #      r, xHII, statistic="mean", bins=r_bin_edges, range=(0.0, outer_edge)
            #  )
            #  xHII_std, _, _ = stats.binned_statistic(
            #      r, xHII, statistic="std", bins=r_bin_edges, range=(0.0, outer_edge)
            #  )
            T_binned, _, _ = stats.binned_statistic(
                r, T, statistic="mean", bins=r_bin_edges, range=(0.0, outer_edge)
            )
            #  T_std, _, _ = stats.binned_statistic(
            #      r, T, statistic="std", bins=r_bin_edges, range=(0.0, outer_edge)
            #  )

            if plot_smoothing_length:
                ax1.plot(
                    [hstar_L, hstar_L],
                    [1e-12, 5],
                    ls="-.",
                    zorder=0,
                    c="C" + str(i),
                    lw=1,
                )
                ax2.plot(
                    [hstar_L, hstar_L],
                    [1, 1e12],
                    ls="-.",
                    zorder=0,
                    c="C" + str(i),
                    lw=1,
                )

            if plot_particles:
                ax1.scatter(r, xHI, **scatterplot_kwargs, zorder=1, c="C" + str(i))
                ax1.scatter(r, xHII, **scatterplot_kwargs, zorder=1, c="C" + str(i))
                # black background lines
                ax1.semilogy(
                    r_bin_centers,
                    xHI_binned,
                    c="k",
                    linewidth=params["lines.linewidth"] * 2.0,
                    zorder=1,
                )
                # ax1.semilogy(
                #     r_bin_centers,
                #     xHII_binned,
                #     c="k",
                #     linewidth=params["lines.linewidth"] * 2.0,
                #     zorder=1,
                # )

            label_left = None
            if s == 0:
                label_left = str(res) + r"$^3$"
            label_right = None
            if i == 0:
                label_right = subdir_labels[s]
            ax1.semilogy(
                r_bin_centers,
                xHI_binned,
                ls=linestyles[s],
                label=label_left,
                zorder=10,
                c="C" + str(i),
            )

            #  ax1.semilogy(r_bin_centers, xHII_binned, label=r"GEARRT $x_{\mathrm{HII}}$", zorder=10)
            #  ax1.errorbar(
            #      r_bin_centers,
            #      xHI_binned,
            #      yerr=xHI_std,
            #      label=str(res) + r"$^3$",
            #      capsize=2,
            #      alpha=0.6,
            #      zorder=10,
            #      c="C" + str(i),
            #  )
            #  ax1.errorbar(
            #      r_bin_centers,
            #      xHII_binned,
            #      yerr=xHII_std,
            #      label=str(res)+r"$^3$",
            #      capsize=2,
            #      alpha=0.6,
            #      zorder=20, c="C"+str(i),
            #  )

            if plot_particles:
                ax2.scatter(r, T, **scatterplot_kwargs, zorder=0)
                # black background lines
                ax2.semilogy(
                    r_bin_centers,
                    T_binned,
                    c="k",
                    linewidth=params["lines.linewidth"] * 2.0,
                    zorder=1,
                )
            #  ax2.errorbar(
            #      r_bin_centers,
            #      T_binned,
            #      yerr=T_std,
            #      label=str(res) + r"$^3$",
            #      capsize=2,
            #      zorder=10,
            #      c="C" + str(i),
            #  )
            ax2.semilogy(
                r_bin_centers,
                T_binned,
                ls=linestyles[s],
                label=label_right,
                zorder=10,
                c="C" + str(i),
            )

    # Cosmetics and save fig
    # ------------------------

    ax1.set_ylim(xHI_min * 0.7, xHI_max * 2)
    ax2.set_ylim(T_min * 0.7, T_max * 2)

    for ax in fig.axes:
        ax.set_xlabel("r / L")
        # note: L is the box size of the original test, not the actual run
        # with GEARRT
        ax.set_xlim(0.0, 1.0)
        ax.set_yscale("log")
        ax.legend()

    ax1.set_ylabel("Neutral Hydrogen Fractions [1]")
    ax2.set_ylabel("Temperature [K]")

    fig.suptitle("Iliev+06 Test 2, $t$ = {0:.0f}".format(meta.time.to("Myr")))
    plt.tight_layout()
    figname = "comparison_" + snapnr
    if plot_refs:
        figname += "-ref"
    if plot_particles:
        figname += "-particles"
    figname += ".png"
    plt.savefig(figname, dpi=200)
    print("saved", figname)
    plt.close()

    return


if __name__ == "__main__":

    if plot_all:
        snapnumbers = get_snapshot_number_list(snapshot_base)
    else:
        snapnumbers = [str(snapnr).zfill(4)]

    for n in snapnumbers:
        plot_solution(n)
