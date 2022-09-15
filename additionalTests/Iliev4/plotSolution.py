#!/usr/bin/env python3

# ----------------------------------------------------
# Plot slices of the neutral hydrogen number density,
# and temperature for SWIFT and reference solutions
#
# Usage: ./plotSolution.py <output_nr>
# ----------------------------------------------------

import sys
import swiftsimio
import unyt
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size, ImageGrid
from matplotlib.colors import LogNorm
from swiftsimio.visualisation.slice import slice_gas
from rainbow4_colormap import rainbow4

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "output"

# parameters for imshow plots
imshow_kwargs = {
    "origin": "lower",
    "cmap": rainbow4
    #  "cmap":"cividis",
}

# parameters for swiftsimio slices
slice_kwargs = {"resolution": 512, "parallel": True}


# which reference to plot
# must be a string!
#  ref = "0.05Myr"
#  ref = "0.1Myr"
#  ref = "0.2Myr"
#  ref = "0.3Myr"
ref = "0.4Myr"

Tmin = 100
Tmax = 1e5
xHImin = 1e-6
xHImax = 1.2

# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
snapnr = -1
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True
mpl.rcParams["mpl_toolkits.legacy_colorbar"] = False


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
    cutoff = int(4 / (128 + 8) * slice_kwargs["resolution"])

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

    references = ["C2Ray", "Crash", "FFTE"]
    nrows = 2
    ncols = len(references) + 1

    fig = plt.figure(figsize=(5.0 * ncols, nrows * 5.6), dpi=200)
    figname = filename[:-5] + ".png"

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0 * meta.boxsize[0].v,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
        0.0 * meta.boxsize[1].v,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
    ]

    axrows = [[] for r in range(nrows)]
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

    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":

    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)

    for f in snaplist:
        plot_result(f)