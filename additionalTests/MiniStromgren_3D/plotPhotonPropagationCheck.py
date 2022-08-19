#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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

# Script adapted from one provided in StromgrenSphere_3D/ example in SWIFT repository

# ----------------------------------------------------------------------
# plots
#   - radiation energies of particles as function of radius
#   - total energy in radial bins
# and compare with expected propagation speed solution.
# Usage:
#   give snapshot number as cmdline arg to plot
#   single snapshot, otherwise this script plots
#   all snapshots available in the workdir.
#   Make sure to select the photon group to plot that
#   doesn't interact with gas to check the *propagation*
#   correctly.
# ----------------------------------------------------------------------

import gc
import os
import sys
import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "propagation_output"

# additional anisotropy estimate plot?
plot_anisotropy_estimate = False

# which photon group to use.
# NOTE: array index, not group number (which starts at 1 for GEAR)
group_index = 0

scatterplot_kwargs = {
    "alpha": 0.1,
    "s": 1,
    "marker": ".",
    "linewidth": 0.0,
    "facecolor": "blue",
}

lineplot_kwargs = {"linewidth": 2}

# -----------------------------------------------------------------------

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True


def analytical_integrated_energy_solution(L, time, r, rmax):
    """
    Compute analytical solution for the sum of the energy
    in bins for given injection rate <L> at time <time> 
    at bin edges <r> and maximal radius <rmax>
    """

    r_center = 0.5 * (r[:-1] + r[1:])
    r0 = r[0]
    Etot = L * time

    if rmax == 0:
        return r_center, np.zeros(r.shape[0] - 1) * Etot.units

    E = np.zeros(r.shape[0] - 1) * Etot.units
    mask = r_center <= rmax
    E[mask] = Etot / (rmax - r0) * (r[1:] - r[:-1])[mask]

    return r_center, E


def analytical_energy_solution(L, time, r, rmax):
    """
    Compute analytical solution for the energy distribution
    for given injection rate <L> at time <time> at radii <r>
    """

    r_center = 0.5 * (r[:-1] + r[1:])
    r0 = r[0]
    Etot = L * time

    if rmax == 0:
        return r_center, np.zeros(r.shape[0] - 1) * Etot.units

    E_fraction_bin = np.zeros(r.shape[0] - 1) * Etot.units
    mask = r_center <= rmax
    dr = r[1:] ** 2 - r[:-1] ** 2
    E_fraction_bin[mask] = 1.0 / (rmax ** 2 - r0 ** 2) * dr[mask]
    bin_surface = dr
    total_weight = Etot / np.sum(E_fraction_bin / bin_surface)
    E = E_fraction_bin / bin_surface * total_weight

    return r_center, E


def analytical_flux_magnitude_solution(L, time, r, rmax, scheme):
    """
    For radiation that doesn't interact with the gas, the
    flux should correspond to the free streaming (optically
    thin) limit. So compute and return that.
    """
    r, E = analytical_energy_solution(L, time, r, rmax)
    if scheme.startswith("GEAR M1closure"):
        F = unyt.c.to(r.units / time.units) * E / r.units ** 3
    elif scheme.startswith("SPH M1closure"):
        F = unyt.c.to(r.units / time.units) * E

    return r, F


def x2(x, a, b):
    return a * x ** 2 + b


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def plot_photons(filename, emin, emax, fmin, fmax):
    """
    Create the actual plot.

    filename: file to work with
    emin: list of minimal nonzero energy of all snapshots
    emax: list of maximal energy of all snapshots
    fmin: list of minimal flux magnitude of all snapshots
    fmax: list of maximal flux magnitude of all snapshots
    """

    print("working on", filename)

    length_units = unyt.pc
    energy_units = 1e30 * unyt.erg
    time_units = unyt.kyr

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    boxsize = meta.boxsize.to(length_units)
    edgelen = min(boxsize[0], boxsize[1])
    edgelen = edgelen.to(length_units)

    xstar = data.stars.coordinates.to(length_units)
    xpart = data.gas.coordinates.to(length_units)
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))
    r = r.to_physical()
    r = r.to(length_units)

    time = meta.time.to(time_units)
    r_expect = meta.time * meta.reduced_lightspeed
    r_expect = r_expect.to(length_units)

    L = None

    use_const_emission_rates = False
    if scheme.startswith("GEAR M1closure"):
        use_const_emission_rates = bool(
            meta.parameters["GEARRT:use_const_emission_rates"]
        )
    elif scheme.startswith("SPH M1closure"):
        use_const_emission_rates = bool(
            meta.parameters["SPHM1RT:use_const_emission_rates"]
        )
    else:
        print("Error: Unknown RT scheme " + scheme)
        exit()

    if use_const_emission_rates:
        # read emission rate parameter as string
        if scheme.startswith("GEAR M1closure"):
            const_emission_rates = (
                spt.trim_paramstr(
                    meta.parameters["GEARRT:star_emission_rates_LSol"].decode("utf-8")
                )
                * unyt.L_Sun
            )
            L = const_emission_rates[group_index]
        elif scheme.startswith("SPH M1closure"):
            units = data.units
            unit_l_in_cgs = units.length.in_cgs()
            unit_v_in_cgs = (units.length / units.time).in_cgs()
            unit_m_in_cgs = units.mass.in_cgs()
            const_emission_rates = (
                spt.trim_paramstr(
                    meta.parameters["SPHM1RT:star_emission_rates"].decode("utf-8")
                )
                * unit_m_in_cgs
                * unit_v_in_cgs ** 3
                / unit_l_in_cgs
            )
            L = const_emission_rates[group_index]
        else:
            print("Error: Unknown RT scheme " + scheme)
            exit()

    L = L.to(energy_units / time_units)

    ncols = 2
    fig = plt.figure(figsize=(5 * ncols, 5.5), dpi=200)

    nbins = 20
    # all these three arrays have unyts
    r_bin_edges = np.linspace(0.5 * edgelen * 1e-3, 0.507 * edgelen, nbins + 1)
    r_bin_centres = 0.5 * (r_bin_edges[1:] + r_bin_edges[:-1])
    r_analytical_bin_edges = np.linspace(
        0.5 * edgelen * 1e-6, 0.507 * edgelen, nbins + 1
    )

    # --------------------------
    # Read in and process data
    # --------------------------

    energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
    energies = energies.to(energy_units)

    particle_count, _ = np.histogram(
        r,
        bins=r_analytical_bin_edges,
        range=(r_analytical_bin_edges[0], r_analytical_bin_edges[-1]),
    )
    L = L.to(energies.units / time.units)

    xlabel_units_str = boxsize.units.latex_representation()
    energy_units_str = energies.units.latex_representation()

    # ------------------------
    # Plot photon energies
    # ------------------------
    ax1 = fig.add_subplot(1, ncols, 1)
    ax1.set_title("Particle Radiation Energies")
    ax1.set_ylabel("Photon Energy [$" + energy_units_str + "$]")

    # don't expect more than float precision
    emin_to_use = max(emin, 1e-5 * emax)

    if use_const_emission_rates:
        # plot entire expected solution
        rA, EA = analytical_energy_solution(L, time, r_analytical_bin_edges, r_expect)
        EA = EA.to(energy_units)
        rA = rA.to(length_units)

        mask = particle_count > 0
        if mask.any():
            EA = EA[mask].to(energy_units)
            rA = rA[mask].to(length_units)
            pcount = particle_count[mask]

            # the particle bin counts will introduce noise.
            # So use a linear fit for the plot. I assume here
            # that the particle number per bin increases
            # proprtional to r^2, which should roughly be the
            # case for the underlying glass particle distribution.
            popt, _ = curve_fit(x2, rA, pcount)
            # I've got no idea why this plot ignores the units and just
            # uses erg. Just make sure you transformed everything explicitly
            # manually before.
            ax1.plot(
                rA.v,
                EA.v / x2(rA.v, popt[0], popt[1]),
                **lineplot_kwargs,
                linestyle="--",
                c="red",
                label="Analytical Solution",
            )

    else:
        # just plot where photon front should be
        ax1.plot(
            [r_expect, r_expect],
            [emin_to_use, emax * 1.1],
            label="expected photon front",
            color="red",
        )

    ax1.scatter(r, energies, **scatterplot_kwargs)
    energies_binned, _, _ = stats.binned_statistic(
        r,
        energies,
        statistic="mean",
        bins=r_bin_edges,
        range=(r_bin_edges[0], r_bin_edges[-1]),
    )
    ax1.plot(
        r_bin_centres, energies_binned, **lineplot_kwargs, label="Mean Radiation Energy"
    )

    ax1.set_ylim(emin_to_use.to(energy_units).v, emax.to(energy_units).v * 1.1)

    # ------------------------------
    # Plot binned photon energies
    # ------------------------------
    ax2 = fig.add_subplot(1, ncols, 2)
    ax2.set_title("Total Radiation Energy in radial bins")
    ax2.set_ylabel("Total Photon Energy [$" + energy_units_str + "$]")

    energies_summed_bin, _, _ = stats.binned_statistic(
        r,
        energies,
        statistic="sum",
        bins=r_bin_edges,
        range=(r_bin_edges[0], r_bin_edges[-1]),
    )
    ax2.plot(
        r_bin_centres,
        energies_summed_bin,
        **lineplot_kwargs,
        label="Total Energy in Bin",
    )
    current_ylims = ax2.get_ylim()
    ax2.set_ylim(emin_to_use, current_ylims[1])

    if use_const_emission_rates:
        # plot entire expected solution
        # Note: you need to use the same bins as for the actual results
        rA, EA = analytical_integrated_energy_solution(L, time, r_bin_edges, r_expect)

        ax2.plot(
            rA,
            EA.to(energies.units),
            **lineplot_kwargs,
            linestyle="--",
            c="red",
            label="Analytical Solution",
        )
    else:
        # just plot where photon front should be
        ax2.plot(
            [r_expect, r_expect],
            ax2.get_ylim(r),
            label="Expected Photon Front",
            color="red",
        )

    # -------------------------------------------
    # Cosmetics that all axes have in common
    # -------------------------------------------
    for ax in fig.axes:
        ax.set_xlabel("r [$" + xlabel_units_str + "$]")
        ax.set_yscale("log")
        ax.set_xlim(0.0, 0.501 * edgelen)
        ax.legend(fontsize="x-small")

    # Add title
    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2f}".format(meta.z)
    title += ", $t$ = {0:.1f}".format(meta.time.to(time_units))
    fig.suptitle(title)

    plt.tight_layout()
    figname = filename[:-5]
    figname += "-PhotonPropagation.png"
    plt.savefig(figname)
    plt.close()
    gc.collect()

    return


def get_plot_boundaries(filenames):
    """
    Get minimal and maximal nonzero photon energy values
    """

    data = swiftsimio.load(filenames[0])
    energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
    emaxguess = energies.max()

    emin = emaxguess
    emax = 0.0
    fmagmin = 1e30
    fmagmax = -10.0

    for f in filenames:
        data = swiftsimio.load(f)

        energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
        mask = energies > 0.0

        if mask.any():

            nonzero_energies = energies[mask]
            this_emin = nonzero_energies.min()
            emin = min(this_emin, emin)

            this_emax = energies.max()
            emax = max(emax, this_emax)

        fx = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "X")
        fy = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "Y")
        fmag = np.sqrt(fx ** 2 + fy ** 2)

        fmagmin = min(fmagmin, fmag.min())
        fmagmax = max(fmagmax, fmag.max())

    return emin, emax, fmagmin, fmagmax


if __name__ == "__main__":

    print(
        "REMINDER: Make sure you selected the correct photon group",
        "to plot, which is hardcoded in this script.",
    )
    snaplist = get_snapshot_list(snapshot_base)
    emin, emax, fmagmin, fmagmax = get_plot_boundaries(snaplist)
    for f in snaplist:
        plot_photons(f, emin, emax, fmagmin, fmagmax)
