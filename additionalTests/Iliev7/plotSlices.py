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


# ----------------------------------------------------
# Plot slices of hydrogen mass fractions, gas number
# densities, temperature, pressure, and mach number
# in combination with references
# ----------------------------------------------------

import gc
import sys

import h5py
import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
from swiftsimio.visualisation.slice import slice_gas

from rainbow4_colormap import rainbow4
import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# chose which reference to plot
#  ref = "1Myr"
#  ref = "3Myr"
#  ref = "10Myr"
ref = "25Myr"
#  ref = "50Myr"


# snapshot basename
snapshot_base = "output"

# parameters for imshow plots
imshow_kwargs = {"origin": "lower"}

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
    "figure.subplot.left": 0.05,
    "figure.subplot.right": 0.90,
    "figure.subplot.bottom": 0.035,
    "figure.subplot.top": 0.93,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.20,
    "figure.dpi": 200,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
    "mpl_toolkits.legacy_colorbar": False,
}
mpl.rcParams.update(params)


# -----------------------------------------------------------------------


def read_reference(code, ref, quantity):
    """
    Read in the reference file

    code: which code to read in
    ref: which refernce (age) to read in
    quantity: which quantity to read in
    """

    filename = "reference/" + code + ".hdf5"
    f = h5py.File(filename, "r")
    slices = f["slices"]
    slice_age = slices[ref]

    quant = unyt.unyt_array(slice_age[quantity], slice_age[quantity].attrs["unyts"])
    f.close()

    return quant


def plot_result(
    swift_slice, simulation_time, filename, quantity, quantity_label, norm_limits
):
    """
    Create and save the plot
    swift_slice:  slice image data of swift output
    simulation_time: unyt quantity of output time
    filname: file name that was rad in
    quantity: which quantity is being plotted
    quantity_label: label of the quantity for the title
    norm_limits: tuple with min/max val for colorbar
    """

    fig = plt.figure(figsize=(14, 10), dpi=200)
    figname = filename[:-5] + "-Slice-" + quantity + ".png"

    # grab all slice data into a single list
    slices = [swift_slice]
    slice_names = ["GEARRT"]

    # You need to coordinate nr codes + SWIFT = ncols * nrows
    # otherwise, not all reference codes will be plotted.
    codes = ["C2Ray+Capreole", "Coral", "Flash", "Licorice", "RSPH", "Zeus-MP"]

    for code in codes:
        refslice = read_reference(code, ref, quantity)
        slices.append(refslice)
        slice_names.append(code)

    ncols = 3
    nrows = 2
    axrows = [[] for _ in range(nrows)]
    for r in range(nrows):

        # set up every column
        axcols = ImageGrid(
            fig,
            (nrows, 1, r + 1),
            nrows_ncols=(1, ncols),
            axes_pad=0.1,
            share_all=True,
            label_mode="L",
            cbar_mode="edge",
            cbar_location="right",
            cbar_size="7%",
            cbar_pad="2%",
        )

        # and store it
        axrows[r] = axcols

    minnorm = norm_limits[0]
    maxnorm = norm_limits[1]

    ind = 0
    for r in range(nrows):
        axcols = axrows[r]

        for c in range(ncols):
            ax = axcols[c]

            data = slices[ind]
            im = ax.imshow(
                data.T,
                **imshow_kwargs,
                norm=LogNorm(vmin=minnorm, vmax=maxnorm),
                #  cmap=rainbow4,
                cmap="cividis",
            )
            ax.set_title(slice_names[ind])

            ind += 1
            ax.set_ylabel("[kpc]")
            ax.set_xlabel("[kpc]")

        axcols.cbar_axes[0].colorbar(im)

    title = r"Iliev+09 Test 7"
    title += r", $t$ = {0:.0f}".format(simulation_time.to("Myr"))
    title += r", " + quantity_label
    fig.suptitle(title)

    #  plt.tight_layout()
    plt.savefig(figname)
    print("saved figure", figname)
    plt.close()
    gc.collect()
    return


def get_swift_output(filename):
    """
    Read in the swift data, and make slices.
    """

    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    simulation_time = meta.time

    ntot = data.gas.masses.shape[0]
    # assume base number of parts is a power of 2
    # you got extra parts as boundaries
    npart_goal = int(ntot ** (1.0 / 3.0))
    npart = 1
    while npart < npart_goal:
        npart *= 2
    npart /= 2
    print("Found resolution", npart)

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0 * meta.boxsize[0].v,
        (1.0 - 8 / (npart + 8)) * meta.boxsize[0].to("kpc").v,
        0.0 * meta.boxsize[1].v,
        (1.0 - 8 / (npart + 8)) * meta.boxsize[0].to("kpc").v,
    ]

    # cut off boundary particles
    cutoff = int(4 / (npart + 8) * slice_kwargs["resolution"])

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
    cs = spt.get_soundspeed_from_density_pressure(data)
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
    HI_map = HI_map[cutoff:-cutoff, cutoff:-cutoff]

    HII_map = mass_weighted_HII_map / mass_map
    HII_map = HII_map[cutoff:-cutoff, cutoff:-cutoff]

    # use 1e30 proton masses here to avoid overflows
    number_density_map = mass_map.to("M_Sun/pc**3") / (
        1e30 * unyt.proton_mass.to("M_Sun")
    )
    number_density_map = number_density_map.to("1e20*pc**-3")
    number_density_map = number_density_map.to("cm**-3")
    number_density_map = number_density_map * 1e30  # fix the 1e30 proton masses used

    number_density_map = number_density_map[cutoff:-cutoff, cutoff:-cutoff]
    number_density_map = number_density_map.to("cm**(-3)")

    pressure_map = mass_weighted_pressure_map / mass_map
    pressure_map = pressure_map[cutoff:-cutoff, cutoff:-cutoff]
    pressure_map = pressure_map.to("g/cm/s**2")

    temperature_map = mass_weighted_temperature_map / mass_map
    temperature_map = temperature_map[cutoff:-cutoff, cutoff:-cutoff]

    mach_map = mass_weighted_mach_map / mass_map
    mach_map = mach_map[cutoff:-cutoff, cutoff:-cutoff]

    slices = {
        "xHI": HI_map,
        "xHII": HII_map,
        "n": number_density_map,
        "T": temperature_map,
        "P": pressure_map,
        "mach": mach_map,
    }

    return slices, simulation_time


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"

    # mach and n might need some rework.
    quantities = [
        "xHI",
        "xHII",
        #  "n",
        "T",
        "P",
        #  "mach",
    ]
    quantity_labels = {
        "xHI": "Neutral Hydrogen Mass Fraction [1]",
        "xHII": "Ionized Hydrogen Mass Fraction [1]",
        "n": r"Hydrogen Number Density [cm$^{-3}$]",
        "T": r"Temperature [K]",
        "P": r"Pressure [g cm$^{-1}$ s$^{-2}$]",
        "mach": r"Mach Number [1]",
    }
    norm_limits = {
        "xHI": (1e-6, 1.2),
        "xHII": (1e-6, 1.2),
        "n": (1.0e-4, 1.0),
        "T": (40, 1e5),
        "P": (1e-16, 1e-12),
        "mach": (1e-5, 4.0),
    }

    swift_slices, simtime = get_swift_output(snap)

    for q in quantities:
        plot_result(
            swift_slices[q], simtime, snap, q, quantity_labels[q], norm_limits[q]
        )
