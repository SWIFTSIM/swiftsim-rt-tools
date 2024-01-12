#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Plot the scaling results.
#
# This script assumes that you have all your timings in a single subdir
# provided by 'subdir' below. They follow the following naming scheme:
#
#   timesteps-$REPLICATION.txt
#
# This script assumes that extra measurement data is available. This
# extra data comes from a modified version of SWIFT. See the README
# for more details.
# ----------------------------------------------------------------------------

import os
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from swift_extract_scaling_timing import get_times_with_extra_data

# Parameters users should/may tweak
subdir = "results"

# How many cores do we use for a scale?
base_cores = 128

lineplot_kwargs = {"linewidth": 2, "alpha": 0.3}

# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True


def plot_scaling(data, replications):
    """
    Plot the individual timings.
    """

    gs = gridspec.GridSpec(
        2,
        1,
        hspace=0.0,
        height_ratios=[
            1.0,
            3.0,
        ],
        left=0.1,
        right=0.99,
    )
    fig = plt.figure(figsize=(8, 4.5), dpi=200)

    ax2 = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharex=ax2)

    times = np.array([d.time for d in data])

    deadtime = np.array([d.deadtime for d in data])
    #  deadtime = deadtime / times
    #  deadtime = deadtime / deadtime[0]
    deadtime = deadtime / times[0]

    drifts = np.array([d.drifts for d in data])
    #  drifts = drifts / times
    #  drifts = drifts / drifts[0]
    drifts = drifts / times[0]

    sorts = np.array([d.sorts for d in data])
    #  sorts = sorts / times
    #  sorts = sorts / sorts[0]
    sorts = sorts / times[0]

    hydro = np.array([d.hydro for d in data])
    #  hydro = hydro / times
    hydro = hydro / hydro[0]

    hydro_ghost = np.array([d.hydro_ghost for d in data])
    #  hydro_ghost = hydro_ghost / times
    hydro_ghost = hydro_ghost / hydro_ghost[0]

    hydro_density = np.array([d.hydro_density for d in data])
    #  hydro_density = hydro_density / times
    hydro_density = hydro_density / hydro_density[0]

    rt_propagation = np.array([d.rt_propagation for d in data])
    #  rt_propagation = rt_propagation / times
    #  rt_propagation = rt_propagation / rt_propagation[0]
    rt_propagation = rt_propagation / times[0]

    rt_tchem = np.array([d.rt_tchem for d in data])
    #  rt_tchem = rt_tchem / times
    #  rt_tchem = rt_tchem / rt_propagation[0]
    rt_tchem = rt_tchem / times[0]

    gravity = np.array([d.gravity for d in data])
    #  gravity = gravity / times
    #  gravity = gravity / gravity[0]
    gravity = gravity / times[0]

    feedback = np.array([d.feedback for d in data])
    #  feedback = feedback/times
    #  feedback = feedback / feedback[0]
    feedback = feedback / times[0]

    integration = np.array([d.integration for d in data])
    #  integration = integration / times
    #  integration = integration / integration[0]
    integration = integration / times[0]

    mpi = np.array([d.mpi for d in data])
    #  mpi = mpi / times
    #  mpi = mpi / mpi[0]
    mpi = mpi / times[0]

    others = np.array([d.others for d in data])
    #  others = others / times
    #  others = others / others[0]
    others = others / times[0]

    tree = np.array([d.tree for d in data])
    #  tree = tree / times
    tree = tree / times[0]

    rt_full = np.array([d.rt_full for d in data])
    #  rt_full = rt_full / times
    rt_full = rt_full / times[0]

    hydro_full = np.array([d.hydro_full for d in data])
    #  hydro_full = hydro_full / times
    hydro_full = hydro_full / times[0]

    times = times / times[0]

    #  print("WARNING: Scaling down times to particle number so plots will look similar to what they should.")
    #  print("WARNING: Disable this once you're actually plotting results.")
    #  times /= 8.**(np.array(replications) - 1)

    #  ymin = max(0., times.min() * 0.5)
    ymin = 1.0e-4
    ymax = 1.25 * times.max()

    fulltimelinestyle = {
        "label": "Full Time",
        "c": "k",
        "alpha": 1.0,
    }
    fulltimelinescatterstyle = {
        "c": "k",
        "alpha": 1.0,
    }

    linestyle = {
        "alpha": 0.6,
        "lw": 2.0,
    }

    scatterstyle = {
        "alpha": 0.6,
        "zorder": 3,
        "s": 12,
    }

    cores = [base_cores * (rep) ** 3 for rep in replications]

    ax.semilogx(cores, times, **fulltimelinestyle)
    ax.scatter(cores, times, **fulltimelinescatterstyle)

    #  ax.semilogx(cores, deadtime, label="Deadtime", **linestyle)
    #  ax.scatter(cores, deadtime, **scatterstyle)

    everything_else = drifts + sorts + feedback + integration + others + mpi

    #  ax.semilogx(cores, drifts, label="Drifts", **linestyle)
    #  ax.scatter(cores, drifts, **scatterstyle)

    #  ax.semilogx(cores, sorts, label="Sorts", **linestyle)
    #  ax.scatter(cores, sorts, **scatterstyle)

    #  ax.semilogx(cores, hydro_density, label="Hydro Density", ls="--", **linestyle)
    #  ax.scatter(cores, hydro_density, **scatterstyle)

    #  ax.semilogx(cores, hydro_ghost, label="Hydro Ghost", ls="--", **linestyle)
    #  ax.scatter(cores, hydro_ghost, **scatterstyle)

    #  ax.semilogx(cores, rt_propagation, label="RT Propagation", ls="--", **linestyle)
    #  ax.scatter(cores, rt_propagation, **scatterstyle)

    #  ax.semilogx(cores, rt_tchem, label="RT Thermochemistry", ls="--", **linestyle)
    #  ax.scatter(cores, rt_tchem, **scatterstyle)

    ax.semilogx(cores, gravity, label="Gravity", **linestyle)
    ax.scatter(cores, gravity, **scatterstyle)

    #  ax.semilogx(cores, feedback, label="Feedback", **linestyle)
    #  ax.scatter(cores, feedback, **scatterstyle)

    #  ax.semilogx(cores, integration, label="Integration", **linestyle)
    #  ax.scatter(cores, integration, **scatterstyle)

    #  ax.semilogx(cores, mpi, label="MPI", **linestyle)
    #  ax.scatter(cores, mpi, **scatterstyle)

    #  ax.semilogx(cores, others, label="Others", **linestyle)
    #  ax.scatter(cores, others, **scatterstyle)

    ax.semilogx(cores, tree, label="Tree", **linestyle)
    ax.scatter(cores, tree, **scatterstyle)

    ax.semilogx(cores, rt_full, label="RT", **linestyle)
    ax.scatter(cores, rt_full, **scatterstyle)

    ax.semilogx(cores, hydro_full, label="Hydro", **linestyle)
    ax.scatter(cores, hydro_full, **scatterstyle)

    ax.semilogx(cores, everything_else, label="Everything Else", **linestyle)
    ax.scatter(cores, everything_else, **scatterstyle)

    # Plot vertical lines
    for i, c in enumerate(cores):
        ax.plot([c, c], [ymin, times[i]], c="gray", linewidth=1, alpha=0.5, zorder=-3)
    # Plot horizontal line
    ax.plot(
        [cores[0], cores[-1]],
        [times[0], times[0]],
        linewidth=1,
        c="gray",
        zorder=-2,
        linestyle=":",
    )

    # ----------------------------------
    # Upper Plot
    # ----------------------------------
    ax2.semilogx(cores, times, **fulltimelinestyle)
    ax2.scatter(cores, times, **fulltimelinescatterstyle)
    plt.sca(ax2)
    plt.tick_params(axis="x", which="both", labelbottom=False, direction="in")
    ax2.set_ylim([0.9, 2.1])
    for i, c in enumerate(cores):
        ax2.plot([c, c], [ymin, times[i]], c="gray", linewidth=1, alpha=0.5, zorder=-3)
    # Plot horizontal line
    ax2.plot(
        [cores[0], cores[-1]],
        [times[0], times[0]],
        linewidth=1,
        c="gray",
        zorder=-2,
        linestyle=":",
    )

    plt.sca(ax)

    ax.set_xlabel("Cores [-]")
    ax.set_ylabel("Time to Solution relative to 128 cores")
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(cores)
    ax.set_xticklabels(["{0:.0f}".format(c) for c in cores])
    ax.legend()

    ax.set_yscale("log")

    #  plt.show()
    #  plt.tight_layout(hspace=0.)
    figname = subdir + ".png"
    plt.savefig(figname)

    print("Finished", figname)
    return


if __name__ == "__main__":

    times, replications = get_times_with_extra_data(subdir, filter_dump_steps=True)
    plot_scaling(times, replications)
