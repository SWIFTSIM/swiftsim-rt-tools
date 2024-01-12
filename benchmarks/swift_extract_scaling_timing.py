#!/usr/bin/env python3

# -------------------------------------------------------------
# Functions related to extracting timing for SWIFT RT scaling
# benchmarks
# -------------------------------------------------------------

import os
import numpy as np


class ExtraData(object):
    def __init__(self):
        """
        Container for timing data with
        extra output written.
        """

        self.time = None
        self.deadtime = None

        self.drifts = None
        self.sorts = None
        self.hydro = None
        self.hydro_ghost = None
        self.hydro_density = None
        self.rt_propagation = None
        self.rt_tchem = None
        self.gravity = None
        self.feedback = None
        self.integration = None
        self.mpi = None
        self.others = None
        self.tree = None

        self.rt_full = None
        self.hydro_full = None

        return


def extract_extra_timing_from_timesteps_file(file, filter_dump_steps):
    """
    Read in the extra timestep data from timesteps file.

    file: str
        file to read in

    filter_dump_steps:
        if True, filter out steps that write stuff.

    returns:
        timing in seconds.
    """

    if not os.path.exists(file):
        raise FileNotFoundError(f"no file {file} found.")

    data = np.loadtxt(
        file,
        dtype=float,
        usecols=[12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
    )

    wallclock = data[:, 0]
    props = data[:, 1]
    deadtime = data[:, 2]
    drifts = data[:, 3]
    sorts = data[:, 4]
    hydro = data[:, 5]
    ghost = data[:, 6]
    density = data[:, 7]
    grav = data[:, 8]
    feedback = data[:, 9]
    integration = data[:, 10]
    mpi = data[:, 11]
    rt_propagation = data[:, 12]
    rt_tchem = data[:, 13]
    others = data[:, 14]
    tree = data[:, 15]

    props = props.astype(int)

    use_timing_mask = np.ones(wallclock.shape[0], dtype=bool)

    if filter_dump_steps:

        is_snapshot = 0b00010000
        is_restart = 0b00100000

        for skip in [is_snapshot, is_restart]:
            mask = np.logical_not(props & skip)
            use_timing_mask = np.logical_and(mask, use_timing_mask)

        steps = np.array(range(0, props.shape[0]))
        skipped = np.logical_not(use_timing_mask)

        np.count_nonzero(skipped)
        print(f"File {file}: Skipping:")
        print("    Steps:", steps[skipped])
        print("    Props:", props[skipped])

    wallclock_tot = np.sum(wallclock[use_timing_mask]) * 1e-3
    deadtime_tot = np.sum(deadtime[use_timing_mask]) * 1e-3
    drifts_tot = np.sum(drifts[use_timing_mask]) * 1e-3
    sorts_tot = np.sum(sorts[use_timing_mask]) * 1e-3
    hydro_tot = np.sum(hydro[use_timing_mask]) * 1e-3
    ghost_tot = np.sum(ghost[use_timing_mask]) * 1e-3
    density_tot = np.sum(density[use_timing_mask]) * 1e-3
    grav_tot = np.sum(grav[use_timing_mask]) * 1e-3
    feedback_tot = np.sum(feedback[use_timing_mask]) * 1e-3
    integration_tot = np.sum(integration[use_timing_mask]) * 1e-3
    mpi_tot = np.sum(mpi[use_timing_mask]) * 1e-3
    rt_propagation_tot = np.sum(rt_propagation[use_timing_mask]) * 1e-3
    rt_tchem_tot = np.sum(rt_tchem[use_timing_mask]) * 1e-3
    others_tot = np.sum(others[use_timing_mask]) * 1e-3
    tree_tot = np.sum(tree[use_timing_mask]) * 1e-3

    rt_full_tot = rt_propagation_tot + rt_tchem_tot
    hydro_full_tot = hydro_tot + ghost_tot + density_tot

    data = ExtraData()
    data.time = wallclock_tot
    data.deadtime = deadtime_tot
    data.drifts = drifts_tot
    data.sorts = sorts_tot
    data.hydro = hydro_tot
    data.hydro_ghost = ghost_tot
    data.hydro_density = density_tot
    data.rt_propagation = rt_propagation_tot
    data.rt_tchem = rt_tchem_tot
    data.gravity = grav_tot
    data.feedback = feedback_tot
    data.integration = integration_tot
    data.mpi = mpi_tot
    data.others = others_tot
    data.tree = tree_tot

    data.rt_full = rt_full_tot
    data.hydro_full = hydro_full_tot

    return data


def extract_timing_from_timesteps_file(file, filter_dump_steps):
    """
    Read in the timestep data from timesteps file.

    file: str
        file to read in

    filter_dump_steps:
        if True, filter out steps that write stuff.

    returns:
        timing in seconds.
    """

    if not os.path.exists(file):
        raise FileNotFoundError(f"no file {file} found.")

    wallclock, props = np.loadtxt(file, dtype=float, usecols=[12, 13], unpack=True)

    props = props.astype(int)

    use_timing_mask = np.ones(wallclock.shape[0], dtype=bool)

    if filter_dump_steps:

        is_snapshot = 0b00010000
        is_restart = 0b00100000

        for skip in [is_snapshot, is_restart]:
            mask = np.logical_not(props & skip)
            use_timing_mask = np.logical_and(mask, use_timing_mask)

        steps = np.array(range(0, props.shape[0]))
        skipped = np.logical_not(use_timing_mask)

        np.count_nonzero(skipped)
        print(f"File {file}: Skipping:")
        print("    Steps:", steps[skipped])
        print("    Props:", props[skipped])

    seconds = np.sum(wallclock[use_timing_mask]) * 1e-3

    return seconds


def get_timestep_files_from_single_subdir(subdir):
    """
    Get list of all timestep files we're going to plot.

    This function assumes that the timesteps files are
    named using the following format:

        timesteps-$REPLICATIONS.txt


    Parameters
    ----------

    subdir: str
        which subdir to search for timing files


    Returns
    -------

    filelist: [str]
        list of filenames to be read in

    replications: [int]
        list of number of replications used to obtain snapshot
    """

    filelist = []
    replications = []

    ls = os.listdir(subdir)
    for entry in ls:

        if entry.startswith("timesteps"):

            # extract replications from directory name
            # Add +1 for underscore after base
            repstring = entry[len("timesteps-") : -4]
            tfilename = "timesteps-" + repstring + ".txt"

            if tfilename != entry:
                print(
                    f"Something wrong with format of entry '{entry}', I expected '{tfilename}'"
                )

            filepath = os.path.join(subdir, tfilename)
            if os.path.exists(filepath):
                filelist.append(filepath)
                replications.append(int(repstring))
            else:
                print("Didn't find", filepath, "; skipping")

    if len(filelist) == 0:
        raise FileNotFoundError("didn't find any timestep files")

    z = sorted(zip(replications, filelist))
    filelist = [f for r, f in z]
    replications = [r for r, f in z]

    return filelist, replications


def get_times(subdir, filter_dump_steps=True):
    """
    Get list of all times and numbers of sub-cycles used we're going to plot.

    Parameters
    ----------

    subdir: str
        which subdir to search for timing files

    filter_dump_steps: bool
        if True, filter out (= don't include) the steps where restarts
        or snapshots are being written


    Returns
    -------

    times: [float]
        list of times for different runs, in seconds

    replications: [int]
        list of number of replications of the entire box
        used to obtain snapshot.
    """

    timefiles, replications = get_timestep_files_from_single_subdir(subdir)

    times = []
    for f in timefiles:
        t = extract_timing_from_timesteps_file(f, filter_dump_steps)
        times.append(t)

    return times, replications


def get_times_with_extra_data(subdir, filter_dump_steps=True):
    """
    Get list of all times and numbers of sub-cycles used we're going to plot.

    Parameters
    ----------

    subdir: str
        which subdir to search for timing files

    filter_dump_steps: bool
        if True, filter out (= don't include) the steps where restarts
        or snapshots are being written


    Returns
    -------

    times: [ExtraData]
        list of timining for different runs, in seconds.

    replications: [int]
        list of number of replications of the entire box
        used to obtain snapshot.
    """

    timefiles, replications = get_timestep_files_from_single_subdir(subdir)

    times = []
    for f in timefiles:
        t = extract_extra_timing_from_timesteps_file(f, filter_dump_steps)
        times.append(t)

    return times, replications
