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

import copy
import os

import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt

# arguments for plots of results
plotkwargs = {"alpha": 0.5}
# arguments for legends
legendprops = {"size": 8}
# snapshot basenames
snapshot_base = "output"
# skip the zeroth snapshot?
skip_zeroth_snapshot = False


# -----------------------------------------------------------------------

energy_units = unyt.Msun * unyt.kpc ** 2 / unyt.kyr ** 2
flux_units = unyt.erg / (unyt.cm ** 2 * unyt.s)
mass_units = unyt.Msun
time_units = unyt.yr


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

    return T


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snap_list = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            if skip_zeroth_snapshot and f == snapshot_base + "_0000.hdf5":
                # skip snapshot zero.
                print("skipping", f)
                continue
            snap_list.append(f)

    if len(snap_list) == 0:
        raise FileNotFoundError(
            "Didn't find any snapshots with basename '" + snapshot_basename + "'"
        )

    snap_list = sorted(snap_list)

    return snap_list


def get_ion_mass_fractions(swiftsimio_loaded_data):
    """
    Returns the ion mass fractions according to
    the used scheme.

    swiftsimio_loaded_data: the swiftsimio.load() object
    """

    data = swiftsimio_loaded_data
    meta = data.metadata
    gas = data.gas
    with_rt = True
    scheme = ""
    try:
        scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    except KeyError:
        # allow to read in solutions with only cooling, without RT
        with_rt = False

    if with_rt:
        if scheme.startswith("GEAR M1closure"):
            imf = data.gas.ion_mass_fractions
        elif scheme.startswith("SPH M1closure"):
            # atomic mass
            mamu = {
                "e": 0.0,
                "HI": 1.0,
                "HII": 1.0,
                "HeI": 4.0,
                "HeII": 4.0,
                "HeIII": 4.0,
            }
            mass_function_hydrogen = data.gas.rt_element_mass_fractions.hydrogen
            imf = copy.deepcopy(data.gas.rt_species_abundances)
            named_columns = data.gas.rt_species_abundances.named_columns
            for column in named_columns:
                # abundance is in n_X/n_H unit. We convert it to mass fraction by multipling mass fraction of H
                mass_function = (
                    getattr(data.gas.rt_species_abundances, column)
                    * mass_function_hydrogen
                    * mamu[column]
                )
                setattr(imf, column, mass_function)
        else:
            print("Error: Unknown scheme", scheme)
            quit(1)
    else:
        # try to find solutions for cooling only runs
        imf = {
            "HI": gas.hi[:],
            "HII": gas.hii[:],
            "HeI": gas.he_i[:],
            "HeII": gas.he_ii[:],
            "HeIII": gas.he_iii[:],
        }

    return imf


def get_snapshot_data(snaplist):
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
        gas = data.gas

        u = gas.internal_energies[:].to(energy_units / mass_units)
        masses = gas.masses[:].to(mass_units)
        volumes_parts = gas.masses[:] / gas.densities[:]
        volumes_parts = volumes_parts.to("kpc**3")
        imf = get_ion_mass_fractions(data)
        mu = mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
        T = gas_temperature(u, mu, gamma).to("K")
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


def get_reference():
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
    #  HI_density = data[:, 6]
    #  HII_density = data[:, 7]
    #  HeI_density = data[:, 8]
    #  HeII_density = data[:, 9]
    #  HeIII_density = data[:, 10]
    #  e_density = data[:, 11]  # number density

    return Time, Temperature


if __name__ == "__main__":

    # ------------------
    # Plot figures
    # ------------------

    snaplist = get_snapshot_list(snapshot_base)
    t, T, mu, mass_fraction, u, photon_energies, volumes, c_reduced = get_snapshot_data(
        snaplist
    )
    ngroups = photon_energies.shape[0]

    t_ref, T_ref = get_reference()

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    ax1.loglog(t, T, label="SWIFT results")
    ax1.loglog(t_ref, T_ref, label="grackle reference", linestyle="--")
    ax1.set_ylabel("gas temperature [K]")
    ax1.legend(prop=legendprops)
    ax1.grid()
    ax1.set_xlim(0.9 * t.min(), 1.01 * t.max())

    ax2.plot(t, mu, label="obtained results")
    ax2.set_ylabel("mean molecular weight")
    ax2.legend(prop=legendprops)
    ax2.grid()

    total_mass_fraction = np.sum(mass_fraction, axis=1)
    ax3.plot(t, total_mass_fraction, "k", label="total", ls="-")

    ax3.loglog(t, mass_fraction[:, 0], label="HI", ls=":", **plotkwargs, zorder=1)
    ax3.loglog(t, mass_fraction[:, 1], label="HII", ls="-.", **plotkwargs, zorder=1)
    ax3.loglog(t, mass_fraction[:, 2], label="HeI", ls=":", **plotkwargs, zorder=1)
    ax3.loglog(t, mass_fraction[:, 3], label="HeII", ls="-.", **plotkwargs, zorder=1)
    ax3.loglog(t, mass_fraction[:, 4], label="HeIII", ls="--", **plotkwargs, zorder=1)
    ax3.legend(loc="upper right", prop=legendprops)
    ax3.set_ylabel("gas mass fractions [1]")
    ax3.grid()

    # for GEAR-RT, the flux corresponds to the energy density flux,
    # i.e. c * E_{snapshot}/V_{snapshot} = F_{snapshot}
    # Try to reconstruct the fluxes that you've injected in the test to verify
    # that you injected the right amount.
    # I expect:
    #  Bin   0:  3.288e+15 -  5.945e+15 [Hz]  Luminosity/cm^2 = 1.350e+01 [erg/s/cm^2]
    #  Bin   1:  5.945e+15 -  1.316e+16 [Hz]  Luminosity/cm^2 = 2.779e+01 [erg/s/cm^2]
    #  Bin   2:  1.316e+16 -  5.879e+17 [Hz]  Luminosity/cm^2 = 6.152e+00 [erg/s/cm^2]

    Ec = photon_energies * c_reduced / volumes
    Ec.convert_to_units(flux_units)
    tot_flux = np.sum(Ec, axis=0)
    ax4.semilogx(
        t, tot_flux, label="total radiation flux", color="k", ls="--", **plotkwargs
    )
    for g in range(ngroups):
        ax4.plot(t, Ec[g, :], label=f"radiation flux group {g+1}", **plotkwargs)
    ax4.set_ylabel(
        r"total radiation flux $E \times \tilde{c}$ [$"
        + Ec.units.latex_representation()
        + "$]",
        usetex=True,
    )
    ax4.legend(prop=legendprops)
    ax4.grid()

    for ax in fig.axes:
        ax.set_xlabel("Time [$" + time_units.latex_representation() + "$]")

    plt.tight_layout()
    plt.savefig("ilievTest0part3-Full.png")
