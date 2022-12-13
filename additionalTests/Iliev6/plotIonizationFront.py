#!/usr/bin/env python3

# --------------------------------------------------------
# This script is intended to reproduce Figure 19
# of Iliev et al. 2009
# (https://ui.adsabs.harvard.edu/abs/2009MNRAS.400.1283I)
# --------------------------------------------------------

import os

import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import unyt
from matplotlib import pyplot as plt

# in some cases, the final snapshot is very close in time
# to the second to last snapshot, and due to roundoff errors,
# we get much too high ionization front velocities. So you might
# want to skip plotting the last velocity estimate.
skip_last = True

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
    time_units = time_units.to("Myr")

    times = unyt.unyt_array(time, time_units)
    rI = unyt.unyt_array(radius, length_units)

    axes = fig.axes
    ax1 = axes[0]
    ax2 = axes[1]

    # get velocity
    dr = rI[1:] - rI[:-1]
    dt = times[1:] - times[:-1]
    vI_full = dr / dt
    stepsize = 1024
    n = int(vI_full.shape[0] / stepsize + 0.5)
    vI = np.zeros(n)
    timesV = np.zeros(n)
    for i in range(n):
        v_av = np.sum(vI_full[i * stepsize : (i + 1) * stepsize]) / stepsize
        vI[i] = v_av
        timesV[i] = times[int((i + 0.5) * stepsize)]

    vI = vI * vI_full.units
    timesV = timesV * times.units

    ax1.plot(times.to("Myr"), rI.to("kpc"), label="GEARRT", zorder=31)
    ax2.plot(timesV.to("Myr"), vI.to("km/s"), label="GEARRT", zorder=31)

    return


def plot_ionization_fronts_reference(fig):
    """
    Read in and plot reference data

    fig: pyplot.figure object to plot into
    """

    references = [
        "Flash_Ifront.dat",
        "Hart_Ifront.dat",
        "Licorice_Ifront.dat",
        "C2Ray+Capreole_Ifront.dat",
        "C2Ray+TVD_Ifront.dat",
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
            times.to("Myr"), rI.to("kpc"), label=label, c=col, alpha=alpha, zorder=1
        )

        mask = vI > 0.0

        ax2.plot(
            times[mask].to("Myr"),
            vI[mask].to("km/s"),
            label=label,
            c=col,
            alpha=alpha,
            zorder=1,
        )


if __name__ == "__main__":

    fig = plt.figure(figsize=(5.5, 8.5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    plot_ionization_fronts_reference(fig)
    plot_ionization_fronts_from_log(fig)

    ax1.set_ylim(0.0, 1.0)

    ax1.set_ylabel(r"$r_I$ [kpc]")
    ax2.set_ylabel(r"$v_I$ [km/s]")
    ax2.set_yscale("log")
    ax2.set_ylim([1, 1e3])
    figname = "ionization_fronts.png"

    ax1.grid()
    ax2.grid()

    ax2.set_xlabel(r"$t$ [Myr]")
    ax2.set_xlim([-1, 30])
    ax1.legend()

    plt.savefig(figname)
