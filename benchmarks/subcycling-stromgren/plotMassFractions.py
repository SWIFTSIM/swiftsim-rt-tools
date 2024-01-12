#!/usr/bin/env python3

# --------------------------------------------
# Plot mass fractions of different subcycle
# numbers on top of each other
# --------------------------------------------

import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import swiftsimio
import os
import unyt

# You might wanna tweak these parameters

# directory basename for each subcycle counter run
directory_base = "stromgren_test"

# snapshot basename
snapshot_basename = "output_HHe"

# wich snapshot to plot?
snapnr = 7

# plot Helium too?
plot_Helium = True

# which subcycles to plot?
nsubcycles = [1, 2, 4, 8, 16, 32, 64, 128]
nsubcycles = [1, 16, 64, 128]
#  nsubcycles = [2, 4]

plotalpha = 0.4


# ------------------------------------------------------


linestyles = [
    "-",
    "--",
    "-.",
    ":",
]  # , "dashdotdotted", "loosely dashdotted", "loosely dashdotdotted", "densely dashdotdotted"]


def get_output_files(directory_base):
    """
    Get list of all snapshots we're going to plot,
    as well as the number of subcycles it was used to produce

    Returns:
    filelist:   list of filenames to be read in
    """

    filelist = []

    ls = os.listdir()

    outputfilename = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"

    # only use number of subcycles specified by user
    for n in nsubcycles:
        dirbase = directory_base + "_" + str(n)

        #  for entry in ls:
        #      if entry.startswith(directory_base) and os.path.isdir(entry):
        #          print("found directory", entry)

        # extract nsubcycles from directory name
        # Add +1 for underscore after base
        #  nstr = entry[len(directory_base)+1:]
        filepath = os.path.join(dirbase, outputfilename)
        if os.path.exists(filepath):
            filelist.append(filepath)
        else:
            print("Didn't find", filepath, "; skipping")

    if len(filelist) == 0:
        print("didn't find any time files")
        exit()

    return filelist


def mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp):
    """
    Determines the mean molecular weight for given
    mass fractions of
        hydrogen:   XH0
        H+:         XHp
        He:         XHe0
        He+:        XHep
        He++:       XHepp

    returns:
        mu: mean molecular weight [in atomic mass units]
        NOTE: to get the actual mean mass, you still need
        to multiply it by m_u, as is tradition in the formulae
    """

    # 1/mu = sum_j X_j / A_j * (1 + E_j)
    # A_H    = 1, E_H    = 0
    # A_Hp   = 1, E_Hp   = 1
    # A_He   = 4, E_He   = 0
    # A_Hep  = 4, E_Hep  = 1
    # A_Hepp = 4, E_Hepp = 2
    one_over_mu = XH0 + 2 * XHp + 0.25 * XHe0 + 0.5 * XHep + 0.75 * XHepp

    return 1.0 / one_over_mu


def gas_temperature(u, mu, gamma):
    """
    Compute the gas temperature given the specific internal
    energy u and the mean molecular weight mu
    """

    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    T = u * (gamma - 1) * mu * unyt.atomic_mass_unit / unyt.boltzmann_constant

    return T.to("K")


def read_data(filename):
    """
    Extract data from swift hdf5 file

    return:
    r: distance from star
    imf: ion mass fractions
    T: gas temperature
    """

    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / (0.5 * boxsize[0])

    # Get mass fractions
    imf = data.gas.ion_mass_fractions

    # get temperature
    mu = mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    T = gas_temperature(data.gas.internal_energies, mu, meta.gas_gamma)

    return r, imf, T


class imf_container:
    """
    Container for ion mass fractions.
    """

    def __init__(self, HI, HII, HeI, HeII, HeIII):
        self.HI = HI
        self.HII = HII
        self.HeI = HeI
        self.HeII = HeII
        self.HeIII = HeIII
        return


def get_bin_averages(r, imf, T):
    """
    Get the averages in radial bins

    returns:
    r_bins: center of radial bins
    imf_av: averaged ion mass fractions
    T_av: averaged gas temperature
    """

    # get profiles
    nbins = 100
    xHI_binned, bin_edges, _ = stats.binned_statistic(
        r, imf.HI, statistic="mean", bins=nbins
    )
    xHII_binned, bin_edges, _ = stats.binned_statistic(
        r, imf.HII, statistic="mean", bins=nbins
    )
    xHeI_binned, bin_edges, _ = stats.binned_statistic(
        r, imf.HeI, statistic="mean", bins=nbins
    )
    xHeII_binned, bin_edges, _ = stats.binned_statistic(
        r, imf.HeII, statistic="mean", bins=nbins
    )
    xHeIII_binned, bin_edges, _ = stats.binned_statistic(
        r, imf.HeIII, statistic="mean", bins=nbins
    )
    T_av, bin_edges, _ = stats.binned_statistic(r, T, statistic="mean", bins=nbins)

    imf_av = imf_container(
        xHI_binned, xHII_binned, xHeI_binned, xHeII_binned, xHeIII_binned
    )

    r_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return r_bins, imf_av, T_av


def plot_results(outputfiles):
    """
    Extract data and create plots
    """

    radii = []
    imfs = []
    Ts = []

    radii_bins = []
    imfs_average = []
    T_average = []

    for f in outputfiles:
        r, imf, T = read_data(f)
        r_bins, imf_av, T_av = get_bin_averages(r, imf, T)
        radii_bins.append(r_bins)
        imfs_average.append(imf_av)
        T_average.append(T_av)

    fig = plt.figure(figsize=(10, 5.5))

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    for i, n in enumerate(nsubcycles):
        if i == 0:
            labels = ["HI", "HII", "HeI", "HeII", "HeIII"]
        else:
            labels = [None, None, None, None, None]
        r = radii_bins[i]
        imf = imfs_average[i]
        T = T_average[i]
        ls = linestyles[i]

        ax1.semilogy(r, imf.HI, c="C0", label=labels[0], ls=ls, alpha=plotalpha)
        ax1.semilogy(r, imf.HII, c="C1", label=labels[1], ls=ls, alpha=plotalpha)
        if plot_Helium:
            ax1.semilogy(r, imf.HeI, c="C2", label=labels[2], ls=ls, alpha=plotalpha)
            ax1.semilogy(r, imf.HeII, c="C3", label=labels[3], ls=ls, alpha=plotalpha)
            ax1.semilogy(r, imf.HeIII, c="C4", label=labels[4], ls=ls, alpha=plotalpha)

        ax2.semilogy(
            r, T, c="C0", label="nsub = {0:d}".format(n), ls=ls, alpha=plotalpha
        )

    for ax in fig.axes:
        ax.legend()
        ax.set_xlabel("r / (L/2)")
        ax.set_xlim(0, 1.1)

    ax1.set_ylabel("Mass Fractions")
    ax1.set_ylabel("Gas Temperature")

    plt.savefig(
        snapshot_basename + "_" + str(snapnr).zfill(4) + "-nsubcycles-gravity.png",
        dpi=200,
    )


if __name__ == "__main__":

    outputfiles = get_output_files(directory_base)
    plot_results(outputfiles)
