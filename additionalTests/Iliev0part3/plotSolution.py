#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

# -------------------------------------------
# Plot the gas temperature, mean molecular
# weight, and mass fractions
# -------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import unyt
import swiftsimio
import os
import copy
import stromgren_plotting_tools as spt

# arguments for plots of results
plotkwargs = {"alpha": 0.5}

# basename
snapshot_basename = "output"

# label individual references, or make them grey?
label_refs = False

# Plot parameters
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
    "text.usetex": True,
    "font.family": "serif",
    "figure.figsize": (10, 4),
    "figure.subplot.left": 0.045,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.0,
    "lines.linewidth": 1.2,
}
mpl.rcParams.update(params)

if label_refs:
    ref_kwargs = {"alpha": 0.8}
else:
    ref_kwargs = {"alpha": 0.6, "c": "grey"}


# -----------------------------------------------------------------------

energy_units = unyt.Msun * unyt.kpc ** 2 / unyt.kyr ** 2
flux_units = unyt.erg / (unyt.cm ** 2 * unyt.s)
mass_units = unyt.Msun
time_units = unyt.yr


def get_full_snapshot_data(snaplist):
    """
    Extract the relevant data from the list of snapshots.

    Returns:
        numpy arrays of:
            time
            temperatures 
            mean molecular weights
            mass fractions
    """

    nsnaps = len(snaplist)
    firstdata = swiftsimio.load(snaplist[0])
    ngroups = int(firstdata.metadata.subgrid_scheme["PhotonGroupNumber"])

    times = np.zeros(nsnaps) * time_units
    temperatures = np.zeros(nsnaps) * unyt.K
    volumes = np.zeros(nsnaps) * unyt.kpc ** 3
    mean_molecular_weights = np.zeros(nsnaps)
    mass_fractions = np.zeros((nsnaps, 5))
    internal_energies = np.zeros(nsnaps) * energy_units
    photon_energies = np.zeros((ngroups, nsnaps)) * energy_units

    for i, snap in enumerate(snaplist):

        data = swiftsimio.load(snap)
        gamma = data.gas.metadata.gas_gamma[0]
        time = data.metadata.time
        scheme = str(data.metadata.subgrid_scheme["RT Scheme"].decode("utf-8"))
        gas = data.gas

        u = gas.internal_energies[:].to(energy_units / mass_units)
        masses = gas.masses[:].to(mass_units)
        volumes_parts = gas.masses[:] / gas.densities[:]
        volumes_parts = volumes_parts.to("kpc**3")
        imf = spt.get_imf(scheme, data)
        mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
        T = spt.gas_temperature(u, mu, gamma).to("K")
        um = u.to(energy_units / mass_units) * masses
        um.to(energy_units)

        times[i] = time.to(time_units)
        temperatures[i] = np.mean(T)
        mean_molecular_weights[i] = np.mean(mu)
        internal_energies[i] = np.mean(um)
        volumes[i] = np.mean(volumes_parts)

        mass_fractions[i, 0] = np.mean(imf.HI)
        mass_fractions[i, 1] = np.mean(imf.HII)
        mass_fractions[i, 2] = np.mean(imf.HeI)
        mass_fractions[i, 3] = np.mean(imf.HeII)
        mass_fractions[i, 4] = np.mean(imf.HeIII)

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            en = en[:].to(energy_units)
            photon_energies[g, i] = en.sum() / en.shape[0]

    c_reduced = data.metadata.reduced_lightspeed

    return (
        times,
        temperatures,
        mean_molecular_weights,
        mass_fractions,
        internal_energies,
        photon_energies,
        volumes,
        c_reduced,
    )


def get_snapshot_data(snaplist):
    """
    Extract the relevant data from the list of snapshots.

    Returns:
        numpy arrays of:
            time
            temperatures 
    """

    nsnaps = len(snaplist)

    times = np.zeros(nsnaps) * time_units
    temperatures = np.zeros(nsnaps) * unyt.K
    xHI = np.zeros(nsnaps)

    for i, snap in enumerate(snaplist):

        data = swiftsimio.load(snap)
        gamma = data.gas.metadata.gas_gamma[0]
        time = data.metadata.time
        gas = data.gas
        scheme = str(data.metadata.subgrid_scheme["RT Scheme"].decode("utf-8"))

        u = gas.internal_energies[:].to(energy_units / mass_units)
        imf = spt.get_imf(scheme, data)
        mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
        T = spt.gas_temperature(u, mu, gamma).to("K")

        times[i] = time.to(time_units)
        temperatures[i] = np.mean(T)

        xHI[i] = np.mean(imf.HI / (imf.HI + imf.HII))

    return times, xHI, temperatures


def get_grackle_reference():
    """
    Read in the temperatures from the reference file
    """

    #  resultfile = "references/reference-grackle-standalone-caseArecombination.dat"
    resultfile = "references/reference-grackle-standalone-caseBrecombination.dat"

    # Read in units.
    # This should still work, but shouldn't be
    # necessary in the current state of this example.
    #  f = open(resultfile, "r")
    #  firstline = f.readline()
    #  massline = f.readline()
    #  lengthline = f.readline()
    #  velline = f.readline()
    #  f.close()
    #  units = []
    #  for l in [massline, lengthline, velline]:
    #      before, after = l.split("used:")
    #      val, unit = after.split("[")
    #      val = val.strip()
    #      units.append(float(val))

    #  mass_units = units[0]
    #  length_units = units[1]
    #  velocity_units = units[2]
    #  time_units = velocity_units / length_units
    #  density_units = mass_units / length_units ** 3

    # Read in all other data
    data = np.loadtxt(resultfile)

    Time = data[:, 1]
    #  Time_Myr = Time * 1e-6
    #  dt = data[:, 2]
    Temperature = data[:, 3]
    #  mu = data[:, 4]
    #  tot_density = data[:, 5]  # mass density
    HI_density = data[:, 6]
    HII_density = data[:, 7]
    #  HeI_density = data[:, 8]
    #  HeII_density = data[:, 9]
    #  HeIII_density = data[:, 10]
    #  e_density = data[:, 11]  # number density

    xHI = HI_density / (HI_density + HII_density)

    return Time, xHI, Temperature


if __name__ == "__main__":

    # ------------------
    # Plot figures
    # ------------------

    snaplist = spt.get_snapshot_list(snapshot_basename)
    #  t, T, mu, mass_fraction, u, photon_energies, volumes, c_reduced = get_full_snapshot_data(
    #      snaplist
    #  )
    t, xHI, T = get_snapshot_data(snaplist)

    fig = plt.figure(figsize=(7, 5), dpi=300)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    reflist = ["FFTE", "RSPH", "FLASH", "ART", "CRASH", "C2RAY"]
    for ref in reflist:
        t_ref, xHI_ref, T_ref = np.loadtxt("references/" + ref + ".dat", unpack=True)
        if label_refs:
            ax1.loglog(t_ref, xHI_ref, label=ref, **ref_kwargs)
            ax2.loglog(t_ref, T_ref, label=ref, **ref_kwargs)
        else:
            label = None
            if ref == reflist[-1]:
                label = "references"
            ax1.loglog(t_ref, xHI_ref, label=label, **ref_kwargs)
            ax2.loglog(t_ref, T_ref, label=label, **ref_kwargs)

    if label_refs:
        grackle_label = "grackle reference"
    else:
        grackle_label = None
    t_ref, xHI_ref, T_ref = get_grackle_reference()
    ax1.loglog(t_ref, xHI_ref, label=grackle_label, **ref_kwargs)
    ax2.loglog(t_ref, T_ref, label=grackle_label, **ref_kwargs)

    ax1.loglog(t, xHI, label="GEARRT")
    ax2.loglog(t, T, label="GEARRT")

    # Cosmetics
    # ------------

    ax1.set_ylabel(r"neutral fraction [1]", usetex=True)
    ax2.set_ylabel(r"gas temperature [K]", usetex=True)
    ax1.legend()
    ax1.grid()
    ax2.grid()
    ax1.tick_params(labelbottom=False)
    #  ax1.set_xlim(0.5 * t.min(), 1.2 * t.max())
    ax1.set_xlim(1.0e-6, 6.0e6)

    ax2.set_xlabel(r"Time [$" + time_units.latex_representation() + "$]")

    plt.tight_layout(h_pad=0.0)
    plt.savefig("ilievTest0part3.png")
    #  plt.show()
