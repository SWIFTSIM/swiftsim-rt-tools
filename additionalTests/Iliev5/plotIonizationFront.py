#!/usr/bin/env python3

import os

import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy import stats

import stromgren_plotting_tools as spt

# --------------------------------------------------------
# This script is intended to reproduce Figure 18
# of Iliev et al. 2009
# (https://ui.adsabs.harvard.edu/abs/2009MNRAS.400.1283I)
# --------------------------------------------------------

# in some cases, the final snapshot is very close in time
# to the second to last snapshot, and due to roundoff errors,
# we get much too high ionization front velocities. So you might
# want to skip plotting the last velocity estimate.
skip_last = True

# if True, use only a few references, and add labels.
# Otherwise, use all references, and make them grey
label_refs = False
#  label_refs = True

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
    "figure.subplot.left": 0.145,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.075,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.0,
    "figure.dpi": 200,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)

snapshot_base = "output"


# In case you want to make sanity checks:
# t_rec = 1./ (alpha_B(T) n_H), Iliev+09 eq 9
t_rec = 122.4 * unyt.Myr

# r_S = (3 * Ndot_gamma / (4*pi*alpha_B(T) * n_H**2)) ** (1/3), Iliev+09 eq 8
r_S = 5.4 * unyt.kpc


def plot_ionization_fronts_from_snapshots(snaplist, fig):
    """
    Use snapshot information to estimate where the ionization front
    it and what its velocity is

    snaplist: list of snapshot names to read in
    fig: pyplot.figure object to plot into
    """

    snaptimes = []

    r_ifront = []
    v_ifront = []

    time = 0 * unyt.Myr
    rI = 0 * unyt.kpc
    for filename in snaplist:
        time_prev = time
        rI_prev = rI
        data = swiftsimio.load(filename)
        meta = data.metadata
        time = meta.time
        boxsize = meta.boxsize
        if time == 0.0 * unyt.Myr:
            continue

        scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
        xstar = data.stars.coordinates
        xpart = data.gas.coordinates
        ids = data.gas.particle_ids
        # cut out boundary particles
        mask = ids < 1000000000
        xpart = xpart[mask]
        dxp = xpart - xstar
        r = np.sqrt(np.sum(dxp ** 2, axis=1))

        # Get mass fractions
        imf = spt.get_imf(scheme, data)
        imf.HI = imf.HI[mask]
        imf.HII = imf.HII[mask]
        imf.HeI = imf.HeI[mask]
        imf.HeII = imf.HeII[mask]
        imf.HeIII = imf.HeIII[mask]
        xH = imf.HI + imf.HII
        xHI = imf.HI / xH
        # xHII = imf.HII / xH

        # get profiles
        r = r.to(boxsize.units)
        #  nbins = 10 # use only few bins for MiniStromgren example
        nbins = 100
        r_bin_edges = np.linspace(0.0, 0.5 * boxsize[0].v, nbins + 1) * boxsize.units
        r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
        xHI_binned, _, _ = stats.binned_statistic(
            r, xHI, statistic="mean", bins=r_bin_edges, range=(0.0, 0.5 * boxsize[0])
        )

        rI = None
        for i in range(nbins):
            if xHI_binned[i] > 0.5:
                if i == 0:
                    rI = r_bin_edges[0]
                else:
                    # interpolate radius
                    dx = xHI_binned[i] - xHI_binned[i - 1]
                    x1 = xHI_binned[i - 1]
                    dr = r_bin_centers[i] - r_bin_centers[i - 1]
                    a = dx / dr
                    b = x1 - dx / dr * r_bin_centers[i - 1]

                    rI = (0.5 - b) / a
                break

        if rI is None:
            print("No rI_current found for snap", filename, "time", time)
            quit()

        vI = (rI - rI_prev) / (time - time_prev)

        r_ifront.append(rI)
        v_ifront.append(vI)

        snaptimes.append(time)

    times_array = unyt.unyt_array(snaptimes)
    times_trec = times_array / t_rec

    r_ifront = unyt.unyt_array(r_ifront)
    v_ifront = unyt.unyt_array(v_ifront)

    timesV_trec = times_trec
    if skip_last:
        timesV_trec = times_trec[:-1]
        v_ifront = v_ifront[:-1]

    axes = fig.axes
    ax1 = axes[0]
    ax2 = axes[1]

    ax1.plot(
        times_trec,
        r_ifront.to("kpc") / r_S.to("kpc"),
        label="GEARRT",
        alpha=0.7,
        zorder=30,
    )
    ax2.plot(
        timesV_trec,
        v_ifront.to("kpc/kyr") / (r_S.to("kpc") / t_rec.to("kyr")).to("kpc/kyr"),
        label="GEARRT",
        alpha=0.7,
        zorder=30,
    )

    return


def plot_ionization_fronts_from_log(fig):
    """
    Read in data collected each step by swift and plot it.

    fig: pyplot.figure object to plot into
    """

    datafile = "ionization_front_radius.txt"

    if not os.path.exists(datafile):
        print("can't find file", datafile)
        print("skipping plotting of stats from each step")
        return

    time, radius = np.loadtxt(datafile, unpack=True)
    length_units = None
    time_units = None
    with open(datafile) as f:
        firstline = f.readline()
        firstline = firstline.strip()
        first_split = firstline.split(" ")
        length_units = float(first_split[-1])

        secondline = f.readline()
        secondline = secondline.strip()
        second_split = secondline.split(" ")
        time_units = float(second_split[-1])

    if length_units is None:
        print("failed to get length units?")
        quit()
    if time_units is None:
        print("failed to get length units?")
        quit()

    length_units = length_units * unyt.cm
    length_units = length_units.to("kpc")
    time_units = time_units * unyt.s
    time_units = time_units.to("kyr")

    times = unyt.unyt_array(time, time_units)
    rI = unyt.unyt_array(radius, length_units)

    axes = fig.axes
    ax1 = axes[0]
    ax2 = axes[1]

    # get velocity
    dr = rI[1:] - rI[:-1]
    dt = times[1:] - times[:-1]
    vI_full = dr / dt
    stepsize = 32
    n = int(vI_full.shape[0] / stepsize + 0.5)
    vI = np.zeros(n)
    timesV = np.zeros(n)
    for i in range(n):
        v_av = np.sum(vI_full[i * stepsize : (i + 1) * stepsize]) / stepsize
        vI[i] = v_av
        timesV[i] = times[int((i + 0.5) * stepsize)]

    vI = vI * vI_full.units
    timesV = timesV * times.units

    ax1.plot(times / t_rec, rI.to("kpc") / r_S.to("kpc"), label="GEARRT", zorder=31)
    ax2.plot(
        timesV / t_rec,
        vI.to("kpc/kyr") / (r_S.to("kpc") / t_rec.to("kyr")).to("kpc/kyr"),
        label="GEARRT",
        zorder=31,
    )

    return


def plot_ionization_fronts_reference(fig):
    """
    Read in and plot reference data

    fig: pyplot.figure object to plot into
    """

    references = [
        "Enzo_Ifront.dat",
        "Flash_Ifront.dat",
        "Hart_Ifront.dat",
        "Licorice_Ifront.dat",
        "C2Ray_Ifront.dat",
        "RH1D_Ifront.dat",
        "RSPH_Ifront.dat",
        "Zeus-MP_Ifront.dat",
    ]

    for ref in references:
        fname = "reference/" + ref

        times, rI, vI = np.loadtxt(fname, unpack=True)

        times = times * unyt.Myr
        rI = rI * unyt.kpc
        vI = vI * unyt.km / unyt.s

        axes = fig.axes
        ax1 = axes[0]
        ax2 = axes[1]

        if label_refs:
            col = None
            label = ref[:-11]
            alpha = 0.8
        else:
            col = "grey"
            alpha = 0.6
            if ref == references[-1]:
                label = "reference"
            else:
                label = None

        ax1.plot(
            times / t_rec,
            rI.to("kpc") / r_S.to("kpc"),
            label=label,
            c=col,
            alpha=alpha,
            zorder=1,
        )

        mask = vI > 0.0

        ax2.plot(
            times[mask] / t_rec,
            vI[mask].to("kpc/kyr") / (r_S.to("kpc") / t_rec).to("kpc/kyr"),
            label=label,
            c=col,
            alpha=alpha,
            zorder=1,
        )


if __name__ == "__main__":

    fig = plt.figure(figsize=(5.5, 8.5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    #  snaplist = spt.get_snapshot_list(snapshot_base, True, -1)
    #  plot_ionization_fronts_from_snapshots(snaplist, fig)
    plot_ionization_fronts_from_log(fig)
    plot_ionization_fronts_reference(fig)

    ax1.set_ylim(0.0, 1.75)

    ax1.set_ylabel(r"$r_I/r_{S,0}$")
    ax2.set_ylabel(r"$v_I/(r_{S,0}/t_{\mathrm{rec},0})$")
    ax2.set_yscale("log")
    figname = "ionization_fronts.png"

    ax1.grid()
    ax2.grid()

    ax2.set_xlabel(r"$t/t_{\mathrm{rec}}$")
    ax1.legend()

    #  plt.show()
    plt.savefig(figname)
