#!/usr/bin/env python3

# -----------------------------------------
# Plot the speedup w.r.t. no subcycling.
# -----------------------------------------

import os
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt

# Parameters users should/may tweak
directory_base = "stromgren_test"
# Which snapshot are we plotting from each directory?
timefile_base = "timing"

# use gravity runs?
#  with_gravity = False
with_gravity = True

# sim name in filenames
#  sim = "MF"
sim = "MFHHe"

# resolution
#  res = "64"
res = "128"

lineplot_kwargs = {"linewidth": 2, "alpha": 0.3}

# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True


def get_timing_files(directory_base):
    """
    Get list of all timings we're going to plot

    Returns:
    filelist:   list of filenames to be read in
    nsubcycles: list of number of subcycles used to obtain snapshot
    """

    filelist = []
    nsubcycles = []

    ls = os.listdir()
    for entry in ls:
        if entry.startswith(directory_base) and os.path.isdir(entry):
            print("found directory", entry)

            # extract nsubcycles from directory name
            # Add +1 for underscore after base
            nstr = entry[len(directory_base) + 1 :]
            if with_gravity:
                tfilename = "-".join((timefile_base, res, sim, nstr, "grav"))
            else:
                tfilename = "-".join((timefile_base, res, sim, nstr))
            filepath = os.path.join(entry, tfilename)
            if os.path.exists(filepath):
                filelist.append(filepath)
                nsubcycles.append(int(nstr))
            else:
                print("Didn't find", filepath, "; skipping")

    if len(filelist) == 0:
        print("didn't find any time files")
        exit()

    z = sorted(zip(nsubcycles, filelist))
    filelist = [f for n, f in z]
    nsubcycles = [n for n, f in z]

    return filelist, nsubcycles


def extract_timing(file):
    """
    Extract the timing from a file.
    """
    time = None
    with open(file) as f:
        firstline = f.readline()
        secline = f.readline()
        timingline = f.readline()
        timingline = timingline.strip()
        tsplit = timingline.split(" ")
        time = tsplit[-1].strip()

    # time is in format HH:MM:SS
    #  time = np.loadtxt(file, skiprows=2, usecols=[2], dtype=str)
    seconds = sum(x * int(t) for x, t in zip([3600, 60, 1], time.split(":")))
    return seconds


def plot_speedup(timefiles, nsubcycles):
    """
    Plot the individual timings.
    """

    times = []
    for f in timefiles:
        t = extract_timing(f)
        times.append(t)

    for i in range(len(timefiles)):
        print(i, nsubcycles[i], times[i], timefiles[i])

    fig = plt.figure(figsize=(5, 5), dpi=200)
    ax = fig.add_subplot(111)

    times = np.array(times)
    times = times / times[0]
    ymin = max(0.0, times.min() * 0.5)

    ax.semilogx(nsubcycles, times)
    ax.scatter(nsubcycles, times)

    # Plot vertical lines
    for i, n in enumerate(nsubcycles):
        ax.plot([n, n], [ymin, times[i]], c="gray", linewidth=1, alpha=0.5, zorder=-3)
    ax.set_xlabel("Number of subcycles")
    ax.set_ylabel("Time to Solution relative to no subcycling")
    ax.set_ylim(ymin, 1.1)
    ax.set_xticks(nsubcycles)
    ax.set_xticklabels(nsubcycles)

    #  plt.show()
    plt.tight_layout()
    if with_gravity:
        figname = "subcycles-speedup-stromgren-gravity.png"
    else:
        figname = "subcycles-speedup-stromgren.png"
    plt.savefig(figname)


if __name__ == "__main__":

    timefiles, nsubcycles = get_timing_files(directory_base)
    plot_speedup(timefiles, nsubcycles)
