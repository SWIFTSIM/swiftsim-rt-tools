#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Plot the speedup w.r.t. no subcycling for different runs.
#
# This script assumes that you have all your timings in multiple subdirs
# provided by 'subdirs' below. They follow the following naming scheme:
#   timefile_base-$NSUBCYCLES
#
# Alternately, you can provide time step files. They follow the following
# naming scheme:
#   timesteps-$NSUBCYCLES.txt
#
# Specify below whether you're working with SLURM timefiles or SWIFT
# timestep files.
# ----------------------------------------------------------------------------

import os
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from swift_extract_subcycling_timing import get_times

# Parameters users should/may tweak

subdirs = ["subcycling-stromgren", "subcycling-stromgren-grav"]
subdir_labels = ["RT + hydro", "RT + hydro + gravity"]

timefile_base = "timing-128-MFHHe"
from_timesteps = True
#  from_timesteps = False

# if reading in timing data from timesteps files, filter out the
# steps where restarts or snapshots are being written
filter_dump_steps = True
#  filter_dump_steps = False


lineplot_kwargs = {
        "linewidth": 2,
        "alpha": 0.3
        }

# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True


def plot_speedup(times, nsubcycles, label, fig=None):
    """
    Plot the individual timings.
    fig: pyplot figure object. If None, a new one will be generated.

    """




if __name__ == "__main__":

    fig = plt.figure(figsize=(5,5), dpi=200)
    ax = fig.add_subplot(111)

    ymin = 1e30

    for s, subdir in enumerate(subdirs):
        times, nsubcycles = get_times(subdir, timefile_base, from_timesteps, filter_dump_steps)

        times = np.array(times)
        times = times / times[0]
        ymin_new = max(0., times.min() * 0.5)
        ymin = min(ymin, ymin_new)

        ax.semilogx(nsubcycles, times, label=subdir_labels[s])
        ax.scatter(nsubcycles, times)

        # Plot vertical lines
        for i, n in enumerate(nsubcycles):
            ax.plot([n, n], [1e-9, times[i]], c="gray", linewidth=1, alpha=0.25, zorder=-3)


    ax.set_xlabel("Number of subcycles")
    ax.set_ylabel("Time to Solution relative to no subcycling")
    ax.set_ylim(ymin, 1.1)
    ax.set_xticks(nsubcycles)
    ax.set_xticklabels(nsubcycles)
    ax.legend()

    # save plot.
    plt.tight_layout()
    if from_timesteps:
        if filter_dump_steps:
            figname = "compare-"+subdir+"-timesteps.png"
        else:
            figname = "compare-"+subdir+"-timesteps-unfiltered.png"
    else:
        figname = "compare-"+subdir+"-timefiles.png"
    plt.savefig(figname)

    print("Finished", figname)
