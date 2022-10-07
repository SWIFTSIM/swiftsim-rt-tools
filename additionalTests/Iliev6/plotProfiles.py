#!/usr/bin/env python3

# ---------------------------------------------------
# Plot profiles of mass fractions, number density,
# pressure, temperature, and mach number
# ---------------------------------------------------

import swiftsimio
import matplotlib as mpl

#  mpl.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import sys
import stromgren_plotting_tools as spt
import unyt
from scipy import stats
import h5py

# plot references?
plot_refs = True
# if plotting references, label them or make them all grey?
label_refs = False

# chose which reference to plot
ref = "1Myr"
#  ref = "3Myr"
#  ref = "10Myr"
#  ref = "25Myr"
#  ref = "75Myr"


snapshot_base = "output"


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
    "figure.figsize": (5, 4),
    "figure.subplot.left": 0.045,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.93,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.22,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)


errorbarkwargs = {"capsize": 2}


def read_reference(code, ref):
    """
    Read in the reference file

    """

    filename = "reference/" + code + "_" + ref + "_profiles.dat"

    data = np.loadtxt(filename)
    xHI = data[:, 0]
    xHI_std = data[:, 1]
    xHII = data[:, 2]
    xHII_std = data[:, 3]
    n = data[:, 4]
    n_std = data[:, 5]
    T = data[:, 6]
    T_std = data[:, 7]
    P = data[:, 8]
    P_std = data[:, 9]
    mach = data[:, 10]
    mach_std = data[:, 11]

    return T, T_std, P, P_std, xHI, xHI_std, xHII, xHII_std, n, n_std, mach, mach_std


def plot_solution(filename):

    # Read in data first
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    # This is the original test setup
    boxsize_ref = 0.8 * unyt.kpc

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    ids = data.gas.particle_ids
    mask = ids > 1000000000
    xpart = xpart[mask]
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r / boxsize_ref

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

    P = data.gas.pressures[mask].to("g/cm/s**2")

    vels = data.gas.velocities[mask]
    vnorm = np.sqrt(np.sum(vels ** 2, axis=1))
    cs = spt.get_soundspeed_from_density_pressure(data)
    # use formula cs = sqrt(p/rho) for *isothermal* sound speed
    cs = cs / np.sqrt(meta.gas_gamma)
    cs = cs[mask]
    mach = vnorm / cs

    number_density = data.gas.densities[mask]
    number_density = number_density.to("kg/cm**3") / unyt.proton_mass

    # get temperature
    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    T = spt.gas_temperature(data.gas.internal_energies[mask], mu, meta.gas_gamma)

    # get profiles
    # max r should be sqrt(3) * boxlen
    nbins = 100
    r_bin_edges = np.linspace(0.0, 1.1, nbins + 1)
    r_bin_centers = 0.5 * (r_bin_edges[:-1] + r_bin_edges[1:])
    xHI_binned, _, _ = stats.binned_statistic(
        r, xHI, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHI_std, _, _ = stats.binned_statistic(
        r, xHI, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHII_binned, _, _ = stats.binned_statistic(
        r, xHII, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    xHII_std, _, _ = stats.binned_statistic(
        r, xHII, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    T_binned, _, _ = stats.binned_statistic(
        r, T, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    T_std, _, _ = stats.binned_statistic(
        r, T, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    P_binned, _, _ = stats.binned_statistic(
        r, P, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    P_std, _, _ = stats.binned_statistic(
        r, P, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    n_binned, _, _ = stats.binned_statistic(
        r, number_density, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    n_std, _, _ = stats.binned_statistic(
        r, number_density, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )
    mach_binned, _, _ = stats.binned_statistic(
        r, mach, statistic="mean", bins=r_bin_edges, range=(0.0, 1.1)
    )
    mach_std, _, _ = stats.binned_statistic(
        r, mach, statistic="std", bins=r_bin_edges, range=(0.0, 1.1)
    )

    fig = plt.figure(figsize=(18, 11))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    ax1.set_title("Neutral Hydrogen Mass Fraction [1]")
    ax2.set_title("Ionized Hydrogen Mass Fraction [1]")
    ax3.set_title(r"Hydrogen Number Density [cm$^{-3}$]")
    ax4.set_title(r"Temperature [K]")
    ax5.set_title(r"Pressure [g cm$^{-1}$ s$^{-2}$]")
    ax6.set_title(r"Mach Number [1]")

    if plot_refs:

        codes = [
            "C2Ray+Capreole",
            "C2Ray+TVD",
            "Flash",
            "HART",
            "Licorice",
            "RH1D",
            "RSPH",
            "Zeus-MP",
        ]

        for code in codes:

            Tref, Tref_std, Pref, Pref_std, xHIref, xHIref_std, xHIIref, xHIIref_std, nref, nref_std, machref, machref_std = read_reference(
                code, ref
            )

            rref = np.linspace(0.0, 1.0, Tref.shape[0])
            dx = 0.5 * (rref[1] - rref[0])
            rref += dx

            if label_refs:
                label = code
                ax1.semilogy(rref, xHIref, label=label, alpha=0.6)
                ax2.semilogy(rref, xHIIref, label=label, alpha=0.6)
                ax3.semilogy(rref, nref, label=label, alpha=0.6)
                ax4.semilogy(rref, Tref, label=label, alpha=0.6)
                ax5.semilogy(rref, Pref, label=label, alpha=0.6)
                ax6.semilogy(rref, machref, label=label, alpha=0.6)

            else:
                if code == codes[-1]:
                    label = "reference"
                else:
                    label = None
                ax1.errorbar(
                    rref,
                    xHIref,
                    yerr=xHIref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )
                ax2.errorbar(
                    rref,
                    xHIIref,
                    yerr=xHIIref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )
                ax3.errorbar(
                    rref,
                    nref,
                    yerr=nref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )
                ax4.errorbar(
                    rref,
                    Tref,
                    yerr=Tref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )
                ax5.errorbar(
                    rref,
                    Pref,
                    yerr=Pref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )
                ax6.errorbar(
                    rref,
                    machref,
                    yerr=machref_std,
                    label=label,
                    **errorbarkwargs,
                    alpha=0.4,
                    c="grey",
                )

    label = r"GEARRT"
    if label_refs:
        ax1.semilogy(r_bin_centers, xHI_binned, label=label)
        ax2.semilogy(r_bin_centers, xHII_binned, label=label)
        ax3.semilogy(r_bin_centers, n_binned, label=label)
        ax4.semilogy(r_bin_centers, T_binned, label=label)
        ax5.semilogy(r_bin_centers, P_binned, label=label)
        ax6.semilogy(r_bin_centers, mach_binned, label=label)

    else:
        ax1.errorbar(
            r_bin_centers, xHI_binned, yerr=xHI_std, label=label, **errorbarkwargs
        )
        ax2.errorbar(
            r_bin_centers, xHII_binned, yerr=xHII_std, label=label, **errorbarkwargs
        )
        ax3.errorbar(r_bin_centers, n_binned, yerr=n_std, label=label, **errorbarkwargs)
        ax4.errorbar(r_bin_centers, T_binned, yerr=T_std, label=label, **errorbarkwargs)
        ax5.errorbar(r_bin_centers, P_binned, yerr=P_std, label=label, **errorbarkwargs)
        ax6.errorbar(
            r_bin_centers, mach_binned, yerr=mach_std, label=label, **errorbarkwargs
        )

    for ax in fig.axes:
        ax.set_xlabel("r / L")
        # note: L is the box size of the original test, not the actual run
        # with GEARRT
        ax.set_xlim(0.0, 1.01)
        ax.grid()
    ax1.set_yscale("log")
    ax1.set_ylim(5e-6, 1.2)

    ax2.set_yscale("log")
    ax2.set_ylim(5e-6, 1.2)

    ax3.set_yscale("log")
    ax3.set_ylim(1.0e-3, 6.0)

    ax4.set_yscale("log")
    ax4.set_ylim(80, 5e4)

    ax5.set_yscale("log")
    ax5.set_ylim(1e-16, 5e-11)

    ax6.set_yscale("log")
    ax6.set_ylim(1e-3, 5.0)
    ax1.legend()

    fig.suptitle("Iliev+09 Test 6, $t$ = {0:.0f}".format(meta.time.to("Myr")))
    #  plt.tight_layout()
    figname = filename[:-5]
    figname += "-Profiles.png"
    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_solution(snap)
