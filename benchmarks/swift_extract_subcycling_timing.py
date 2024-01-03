#!/usr/bin/env python3

# -------------------------------------------------------------
# Functions related to extracting timing for SWIFT sub-cycling
# benchmarks
# -------------------------------------------------------------

import os
import numpy as np


def extract_timing_from_timingfile(file):
    """
    Function to read in data from SLURM timing files.
    Expects output written using the following command
    (called at the end of your batch script):
    sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,Elapsed | tee timingfile

    returns:
        timing in seconds.
    """

    if not os.path.exists(file):
        raise FileNotFoundError(f"no file {file} found.")

    time = None
    with open(file) as f:
        firstline = f.readline()
        secline = f.readline()
        timingline = f.readline()
        timingline = timingline.strip()
        tsplit = timingline.split(" ")
        time = tsplit[-1].strip()

    days = 0
    if "-" in time:
        # remove days part from time string
        # output may be in format DD-HH:MM:SS
        days, rest = time.split("-")
        days = int(days)
        time = rest

    # time is in format HH:MM:SS
    seconds = sum(x * int(t) for x, t in zip([3600, 60, 1], time.split(":")))
    seconds += days * 3600 * 24

    return seconds


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


def get_timingfiles_from_single_subdir(subdir, timefile_base):
    """
    Get list of all timing files we're going to plot.
    This function expects that all the timefiles have been
    named with the following format:

        timefile_base-$NSUBCYCLES

    and generated using the following SLURM command:

        $ sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,Elapsed


    Parameters
    ----------

    subdir: str
        which subdir to search for timing files

    timefile_base: str
        base filename for SLURM timing files.


    Returns
    -------

    filelist: [str]
        list of filenames to be read in

    nsubcycles: [int]
        list of number of subcycles used to obtain snapshot
    """

    filelist = []
    nsubcycles = []

    ls = os.listdir(subdir)
    for entry in ls:

        if entry.startswith(timefile_base):

            # extract nsubcycles from directory name
            # Add +1 for underscore after base
            nstr = entry[len(timefile_base) + 1 :]
            tfilename = "-".join((timefile_base, nstr))

            if tfilename != entry:
                print(
                    f"Something wrong with format of entry '{entry}', I expected '{tfilename}'"
                )

            filepath = os.path.join(subdir, tfilename)
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


def get_timestep_files_from_single_subdir(subdir):
    """
    Get list of all timestep files we're going to plot.

    This function assumes that the timesteps files are
    named using the following format:

        timesteps-$NSUBCYCLES.txt


    Parameters
    ----------

    subdir: str
        which subdir to search for timing files


    Returns
    -------

    filelist: [str]
        list of filenames to be read in

    nsubcycles: [int]
        list of number of subcycles used to obtain snapshot
    """

    filelist = []
    nsubcycles = []

    ls = os.listdir(subdir)
    for entry in ls:

        if entry.startswith("timesteps"):

            # extract nsubcycles from directory name
            # Add +1 for underscore after base
            nstr = entry[len("timesteps") + 1 : -4]
            tfilename = "timesteps-" + nstr + ".txt"

            if tfilename != entry:
                print(
                    f"Something wrong with format of entry '{entry}', I expected '{tfilename}'"
                )

            filepath = os.path.join(subdir, tfilename)
            if os.path.exists(filepath):
                filelist.append(filepath)
                nsubcycles.append(int(nstr))
            else:
                print("Didn't find", filepath, "; skipping")

    if len(filelist) == 0:
        raise FileNotFoundError("didn't find any timestep files")

    z = sorted(zip(nsubcycles, filelist))
    filelist = [f for n, f in z]
    nsubcycles = [n for n, f in z]

    return filelist, nsubcycles


def get_times(subdir, timefile_base, from_timesteps, filter_dump_steps=True):
    """
    Get list of all times and numbers of sub-cycles used we're going to plot.

    Parameters
    ----------

    subdir: str
        which subdir to search for timing files

    timefile_base: str
        base filename for SLURM timing files.

    from_timesteps: bool
        if True, use timesteps.txt files. Otherwise, use SLURM timing files.

    filter_dump_steps: bool
        if True and from_timesteps==True, when reading in timing data from
        timesteps files, filter out (= don't include) the steps where restarts
        or snapshots are being written


    Returns
    -------

    times: [float]
        list of times for different runs, in seconds

    nsubcycles: [int]
        list of number of subcycles used to obtain snapshot
    """

    if from_timesteps:
        timefiles, nsubcycles = get_timestep_files_from_single_subdir(subdir)
    else:
        timefiles, nsubcycles = get_timingfiles_from_single_subdir(
            subdir, timefile_base
        )

    times = []
    for f in timefiles:
        if from_timesteps:
            t = extract_timing_from_timesteps_file(f, filter_dump_steps)
        else:
            t = extract_timing_from_timingfile(f)
        times.append(t)

    return times, nsubcycles


def get_times_from_timesteps(subdir, filter_dump_steps=True):
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

    nsubcycles: [int]
        list of number of subcycles used to obtain snapshot
    """

    return get_times(
        subdir,
        timefile_base="",
        from_timesteps=True,
        filter_dump_steps=filter_dump_steps,
    )
