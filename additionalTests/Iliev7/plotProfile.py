#!/usr/bin/env python3

# ---------------------------------------------------
# Plot profiles of mass fractions, number density,
# pressure, temperature, and mach number
# ---------------------------------------------------

import sys

import h5py
import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from swiftsimio.visualisation.slice import slice_gas

import stromgren_plotting_tools as spt

# plot references?
plot_refs = True
# if plotting references, label them or make them all grey?
label_refs = False

# chose which reference to plot
#  ref = "1Myr"
#  ref = "3Myr"
#  ref = "10Myr"
#  ref = "25Myr"
ref = "50Myr"


snapshot_base = "output"

# parameters for swiftsimio slices
slice_kwargs = {"resolution": 512, "parallel": True}

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


def read_reference(code, ref):
    """
    Read in the reference file

    """

    filename = "reference/" + code + ".hdf5"
    f = h5py.File(filename, "r")
    profiles = f["profiles"]
    profiles_age = profiles[ref]

    P = unyt.unyt_array(profiles_age["P"], profiles_age["P"].attrs["unyts"])
    T = unyt.unyt_array(profiles_age["T"], profiles_age["T"].attrs["unyts"])
    xHI = unyt.unyt_array(profiles_age["xHI"], profiles_age["xHI"].attrs["unyts"])
    xHII = unyt.unyt_array(profiles_age["xHII"], profiles_age["xHII"].attrs["unyts"])
    n = unyt.unyt_array(profiles_age["n"], profiles_age["n"].attrs["unyts"])
    mach = unyt.unyt_array(profiles_age["mach"], profiles_age["mach"].attrs["unyts"])

    f.close()

    return T, P, xHI, xHII, n, mach


def plot_solution(filename):

    # Read in data first
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    ntot = data.gas.masses.shape[0]
    # assume base number of parts is a power of 2
    # you got extra parts as boundaries
    npart_goal = int(ntot ** (1.0 / 3.0))
    npart = 1
    while npart < npart_goal:
        npart *= 2
    npart /= 2
    print("Found resolution", npart)

    mass_map = slice_gas(
        data, project="masses", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    gamma = meta.gas_gamma

    imf = spt.get_imf(scheme, data)

    data.gas.mXHI = imf.HI * data.gas.masses.to("M_Sun")
    data.gas.mXHII = imf.HII * data.gas.masses.to("M_Sun")
    data.gas.mm = data.gas.masses.to("M_Sun") ** 2

    vels = data.gas.velocities
    vnorm = np.sqrt(np.sum(vels ** 2, axis=1))
    cs = spt.get_soundspeed_from_internal_energy(data)
    # use formula cs = sqrt(p/rho) for *isothermal* sound speed
    cs = cs / np.sqrt(meta.gas_gamma)
    mach = vnorm / cs
    data.gas.mmach = mach * data.gas.masses.to("M_Sun")

    data.gas.mP = data.gas.pressures * data.gas.masses.to("M_Sun")

    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    data.gas.mT = spt.gas_temperature(
        data.gas.internal_energies, mu, gamma
    ) * data.gas.masses.to("M_Sun")

    mass_weighted_HI_map = slice_gas(
        data, project="mXHI", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_HII_map = slice_gas(
        data, project="mXHII", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_mach_map = slice_gas(
        data, project="mmach", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_temperature_map = slice_gas(
        data, project="mT", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_pressure_map = slice_gas(
        data, project="mP", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    HI_map = mass_weighted_HI_map / mass_map

    HII_map = mass_weighted_HII_map / mass_map

    # use 1e30 proton masses here to avoid overflows
    number_density_map = mass_map.to("M_Sun/pc**3") / (
        1e30 * unyt.proton_mass.to("M_Sun")
    )
    number_density_map = number_density_map.to("1e20*pc**-3")
    number_density_map = number_density_map.to("cm**-3")
    number_density_map = number_density_map * 1e30  # fix the 1e30 proton masses used
    number_density_map = number_density_map.to("cm**(-3)")

    pressure_map = mass_weighted_pressure_map / mass_map
    pressure_map = pressure_map.to("g/cm/s**2")

    temperature_map = mass_weighted_temperature_map / mass_map

    mach_map = mass_weighted_mach_map / mass_map

    n = HI_map.shape[0]
    x = np.linspace(0.5 / n, (n - 0.5), n) / n * meta.boxsize[0]

    L_test = 6.6 * unyt.kpc
    shift = (meta.boxsize[0] - L_test) * 0.5
    shiftint = int(shift / meta.boxsize[0] * n)
    x = x[shiftint:-shiftint] - shift

    xHI_profile = HI_map.T[int(n / 2), shiftint:-shiftint]
    xHII_profile = HII_map.T[int(n / 2), shiftint:-shiftint]
    n_profile = number_density_map.T[int(n / 2), shiftint:-shiftint]
    T_profile = temperature_map.T[int(n / 2), shiftint:-shiftint]
    P_profile = pressure_map.T[int(n / 2), shiftint:-shiftint]
    mach_profile = mach_map.T[int(n / 2), shiftint:-shiftint]

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

        codes = ["C2Ray+Capreole", "Coral", "Flash", "Licorice", "RSPH", "Zeus-MP"]

        for code in codes:

            Tref, Pref, xHIref, xHIIref, nref, machref = read_reference(code, ref)

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
                ax6.plot(rref, machref, label=label, alpha=0.6)

            else:
                if code == codes[-1]:
                    label = "reference"
                else:
                    label = None
                ax1.semilogy(rref, xHIref, label=label, alpha=0.4, c="grey")
                ax2.semilogy(rref, xHIIref, label=label, alpha=0.4, c="grey")
                ax3.semilogy(rref, nref, label=label, alpha=0.4, c="grey")
                ax4.semilogy(rref, Tref, label=label, alpha=0.4, c="grey")
                ax5.semilogy(rref, Pref, label=label, alpha=0.4, c="grey")
                ax6.plot(rref, machref, label=label, alpha=0.4, c="grey")

    label = r"GEARRT"
    ax1.semilogy(x / L_test, xHI_profile, label=label)
    ax2.semilogy(x / L_test, xHII_profile, label=label)
    ax3.semilogy(x / L_test, n_profile, label=label)
    ax4.semilogy(x / L_test, T_profile, label=label)
    ax5.semilogy(x / L_test, P_profile, label=label)
    ax6.plot(x / L_test, mach_profile, label=label)

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
    ax3.set_ylim(1.0e-4, 1.0)

    ax4.set_yscale("log")
    ax4.set_ylim(40, 5e4)

    ax5.set_yscale("log")
    ax5.set_ylim(1e-16, 1e-12)

    #  ax6.set_yscale("log")
    ax6.set_ylim(0.0, 4.0)
    ax1.legend()

    fig.suptitle("Iliev+09 Test 7, $t$ = {0:.0f}".format(meta.time.to("Myr")))
    #  plt.tight_layout()
    figname = filename[:-5]
    figname += "-Profiles.png"
    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_solution(snap)
