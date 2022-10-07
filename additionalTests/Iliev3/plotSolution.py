#!/usr/bin/env python3

# ----------------------------------------------------
# Plot slices of the neutral hydrogen number density,
# and temperature for SWIFT and reference solutions
#
# Usage: ./plotSolution.py <output_nr>
# ----------------------------------------------------

import sys

import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
from swiftsimio.visualisation.slice import slice_gas

import stromgren_plotting_tools as spt

# Parameters users should/may tweak


# which reference to use?
#  ref = "1Myr"
#  ref = "2Myr"
#  ref = "3Myr"
#  ref = "5Myr"
ref = "15Myr"

resolution = 128
#  resolution = 64

# snapshot basename
snapshot_base = "output"

# parameters for imshow plots
imshow_kwargs = {
    "origin": "lower",
    #  "cmap": rainbow4
    "cmap": "cividis",
}

# parameters for swiftsimio slices
slice_kwargs = {"resolution": 512, "parallel": True}


Tmin = 30
Tmax = 1e5
xHImin = 1e-6
xHImax = 1.2

#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
    #  "legend.fontsize": 14,
    #  "xtick.labelsize": 12,
    #  "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    #  "xtick.major.width": 1.5,
    #  "ytick.major.width": 1.5,
    "axes.linewidth": 1.5,
    "text.usetex": True,
    "figure.subplot.left": 0.015,
    "figure.subplot.right": 0.98,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.98,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.06,
    #  "lines.markersize": 1,
    #  "lines.linewidth": 2.0,
    "mpl_toolkits.legacy_colorbar": False,
}
mpl.rcParams.update(params)


# -----------------------------------------------------------------------


def get_ref_data(code, quantity):

    datafile = "reference/" + code + "_" + ref + "_slice_" + quantity + "_z=64.dat"
    print("reading", datafile)
    res = np.loadtxt(datafile, delimiter=",")
    return res


def plot_result(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    # Grab your own data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    # cut off boundary particles
    cutoff = int(4 / (resolution + 8) * slice_kwargs["resolution"])

    mass_map = slice_gas(
        data, project="masses", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    gamma = meta.gas_gamma

    imf = spt.get_imf(scheme, data)

    data.gas.mXHI = imf.HI * data.gas.masses.to("M_Sun")

    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    data.gas.mT = spt.gas_temperature(
        data.gas.internal_energies, mu, gamma
    ) * data.gas.masses.to("M_Sun")

    mass_weighted_HI_map = slice_gas(
        data, project="mXHI", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_temperature_map = slice_gas(
        data, project="mT", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    mass_map = mass_map[cutoff:-cutoff, cutoff:-cutoff]
    mass_weighted_temperature_map = mass_weighted_temperature_map[
        cutoff:-cutoff, cutoff:-cutoff
    ]
    mass_weighted_HI_map = mass_weighted_HI_map[cutoff:-cutoff, cutoff:-cutoff]

    HI_map = mass_weighted_HI_map / mass_map
    temperature_map = mass_weighted_temperature_map / mass_map

    references = ["C2Ray", "Coral", "Crash", "FFTE", "Flash", "IFT", "RSPH"]
    nrows = 2
    ncols = len(references) + 1

    fig = plt.figure(figsize=(4.5 * ncols, nrows * 5.6), dpi=200)
    figname = filename[:-5] + ".png"

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0 * meta.boxsize[0].v,
        (1.0 - 8 / (resolution + 8)) * meta.boxsize[0].to("kpc").v,
        0.0 * meta.boxsize[1].v,
        (1.0 - 8 / (resolution + 8)) * meta.boxsize[0].to("kpc").v,
    ]

    axrows = [[] for _ in range(nrows)]
    # loop over each row (T, xHI)
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

    # Plot Neutral Fractions
    axcols = axrows[0]
    for col, ax in enumerate(axcols):

        if col == 0:
            # plot your own data
            im = ax.imshow(
                HI_map.T, **imshow_kwargs, norm=LogNorm(vmin=xHImin, vmax=xHImax)
            )
            ax.set_title("Neutral Hydrogen fraction GEARRT")

        else:
            # plot reference data
            code = references[col - 1]
            refdata = get_ref_data(code, "xHI")
            im = ax.imshow(
                refdata.T, **imshow_kwargs, norm=LogNorm(vmin=xHImin, vmax=xHImax)
            )
            ax.set_title("Neutral Hydrogen fraction " + code)

        # Add colorbar to every row
        axcols.cbar_axes[0].colorbar(im)

        ax.set_xlabel("[kpc]")
        ax.set_ylabel("[kpc]")

    # Plot Temperature
    axcols = axrows[1]
    for col, ax in enumerate(axcols):

        if col == 0:
            # plot your own data
            im = ax.imshow(
                temperature_map.T, **imshow_kwargs, norm=LogNorm(vmin=Tmin, vmax=Tmax)
            )
            ax.set_title("Temperature GEARRT")

        else:
            # plot reference data
            code = references[col - 1]
            refdata = get_ref_data(code, "T")
            im = ax.imshow(
                refdata.T, **imshow_kwargs, norm=LogNorm(vmin=Tmin, vmax=Tmax)
            )
            ax.set_title("Temperature " + code)

        # Add colorbar to every row
        axcols.cbar_axes[0].colorbar(im)

        ax.set_xlabel("[kpc]")
        ax.set_ylabel("[kpc]")

    title = "Iliev+06 Test 3"
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.1f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_result(snap)
