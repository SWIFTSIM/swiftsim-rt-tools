#!/usr/bin/env python3

import swiftsimio
import matplotlib as mpl

#  mpl.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import sys, os
import stromgren_plotting_tools as spt
import unyt
from scipy import stats


#-------------------------------------------------
# This script is intended to reproduce Figure 3
# of Iliev et al 2006 
# (https://arxiv.org/pdf/astro-ph/0508416.pdf)
#-------------------------------------------------

# in some cases, the final snapshot is very close in time 
# to the second to last snapshot, and due to roundoff errors,
# we get much too high ionization front velocities. So you might 
# want to skip plotting the last velocity estimate.
skip_last = True

# plot actual values for rI and vI instead of ratios with
# analytical solutions
plot_actual_values = False

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
    "figure.subplot.left": 0.145,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.075,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.,
    "figure.dpi": 200,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)

snapshot_base = "output"



# In case you want to make sanity checks:
# t_rec = 1./ (alpha_B(T) n_H), Iliev+06 eq 9
t_rec = 122.4 * unyt.Myr

# r_S = (3 * Ndot_gamma / (4*pi*alpha_B(T) * n_H**2)) ** (1/3), Iliev+06 eq 7
r_S = 5.4 * unyt.kpc

def r_I(t):
    """
    Expected analytical I-front radius assuming T=1e4K
    (Iliev+06 eq 5)
    """

    return r_S * (1 - np.exp(-t/t_rec))**(1./3.)

def v_I(t):
    """
    Expected analytical I-front velocity assuming T=1e4K
    (Iliev+06 eq 6)
    """
    expterm = np.exp(-t/t_rec)

    return r_S / (3. * t_rec) * expterm / (1 - expterm)**(2./3.)


def plot_ionization_fronts_from_snapshots(snaplist, fig):
    """
    Use snapshot information to estimate where the ionization front
    it and what its velocity is

    snaplist: list of snapshot names to read in
    fig: pyplot.figure object to plot into
    """

    snaptimes = []
    r_ifront_ana = []
    v_ifront_ana = []

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
        if (time == 0. * unyt.Myr): continue

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
        xHII = imf.HII / xH

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
                    dx = xHI_binned[i] - xHI_binned[i-1]
                    x1 = xHI_binned[i-1]
                    dr = r_bin_centers[i] - r_bin_centers[i-1] 
                    a = dx / dr
                    b =  x1 - dx/dr * r_bin_centers[i-1]

                    rI = (0.5 - b ) / a
                break

        if rI is None:
            print("No rI_current found for snap", filename, "time", time)
            quit()

        vI = (rI - rI_prev) / (time - time_prev)

        r_ifront.append(rI)
        v_ifront.append(vI)

        r_ifront_ana.append(r_I(time))
        v_ifront_ana.append(v_I(time))

        snaptimes.append(time)


    times_array = unyt.unyt_array(snaptimes)
    times_trec = times_array / t_rec

    r_ifront = unyt.unyt_array(r_ifront)
    v_ifront = unyt.unyt_array(v_ifront)
    r_ifront_ana = unyt.unyt_array(r_ifront_ana)
    v_ifront_ana = unyt.unyt_array(v_ifront_ana)
    
    timesV_trec = times_trec
    if skip_last:
        timesV_trec = times_trec[:-1]
        v_ifront = v_ifront[:-1]
        v_ifront_ana = v_ifront_ana[:-1]

    axes = fig.axes
    ax1 = axes[0]
    ax2 = axes[1]

    if plot_actual_values:
        ax1.plot(times_trec, r_ifront_ana.to("kpc"), label="analytical")
        ax1.plot(times_trec, r_ifront.to("kpc"), label="GEARRT")
        ax2.plot(timesV_trec, v_ifront_ana.to("kpc/kyr"), label="analytical")
        ax2.plot(timesV_trec, v_ifront.to("kpc/kyr"), label="GEARRT")
    
    else:
        ax1.plot(times_trec, r_ifront / r_ifront_ana, label="GEARRT")
        ax2.plot(timesV_trec, v_ifront / v_ifront_ana, label="GEARRT")

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

    dr = rI[1:] - rI[:-1]
    dt = times[1:] - times[:-1]
    vI = dr/dt
    timesV = times[1:]

    if plot_actual_values:
        ax1.plot(times/t_rec, rI.to("kpc"), ls="--", label="GEARRT each step")
        ax2.plot(timesV/t_rec, vI.to("kpc/kyr"), ls="--", label="GEARRT each step")

    else:
        r_ifront_ana = r_I(times)
        v_ifront_ana = v_I(timesV)
        ax1.plot(times/t_rec, rI.to("kpc")/r_ifront_ana.to("kpc"), ls="--", label="GEARRT each step")
        ax2.plot(timesV/t_rec, vI.to("kpc/kyr")/v_ifront_ana.to("kpc/kyr"), ls="--", label="GEARRT each step")

    return

   


if __name__ == "__main__":

    fig = plt.figure(figsize=(5.5, 7.5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    snaplist = spt.get_snapshot_list(snapshot_base, True, -1)
    plot_ionization_fronts_from_snapshots(snaplist, fig)
    plot_ionization_fronts_from_log(fig)

    if plot_actual_values:
        ax1.set_ylabel("$r_I$ [kpc]")
        ax2.set_ylabel("$v_I$ [kpc/kyr]")
    else:
        ax1.set_ylabel("$r_I/r_{analyt,0}$")
        ax2.set_ylabel("$v_I/v_{analyt,0}$")

    #  for ax in fig.axes:
    ax2.set_xlabel("$t/t_{rec}$")
    ax2.legend()

    #  plt.show()
    plt.savefig("ionization_fronts.png")
