#!/usr/bin/env python3

# ----------------------------------------------------
# Plot contours of ionization front (defined by ion
# mass fraction = 0.5)
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from swiftsimio.visualisation.slice import slice_gas

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

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


# which reference to plot
# must be a string!
ref = "0.05Myr"
#  ref = "0.1Myr"
#  ref = "0.2Myr"
#  ref = "0.3Myr"
#  ref = "0.4Myr"

Tmin = 100
Tmax = 1e5
xHImin = 1e-6
xHImax = 1.1

#  Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "font.family": "serif",
    "legend.fontsize": 7,
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
    "figure.subplot.left": 0.05,
    "figure.subplot.right": 0.97,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.93,
    "figure.subplot.wspace": 0.20,
    "figure.subplot.hspace": 0.06,
    #  "lines.markersize": 1,
    #  "lines.linewidth": 2.0,
    "mpl_toolkits.legacy_colorbar": False,
}
mpl.rcParams.update(params)


# -----------------------------------------------------------------------


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


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

    global imshow_kwargs
    extent = [
        0.0,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
        0.0,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
    ]
    imshow_kwargs["extent"] = extent

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

    fig = plt.figure(figsize=(11, 5.5), dpi=200)

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Plot background data first

    im1 = ax1.imshow(HI_map.T, **imshow_kwargs, norm=LogNorm(vmin=xHImin, vmax=xHImax))
    ax1.set_title("Neutral Hydrogen Fraction")
    ax1.set_xlabel("[kpc]")
    ax1.set_ylabel("[kpc]")
    set_colorbar(ax1, im1)

    im2 = ax2.imshow(
        temperature_map.T, **imshow_kwargs, norm=LogNorm(vmin=Tmin, vmax=Tmax)
    )

    ax2.set_title("Temperature")
    ax2.set_xlabel("[kpc]")
    ax2.set_ylabel("[kpc]")
    set_colorbar(ax2, im2)

    # Get SWIFT contours

    edge = np.linspace(extent[0], extent[1], HI_map.shape[0])
    x, y = np.meshgrid(edge, edge)
    levels_xHI = [0.01, 0.5]
    level_labels_xHI = ["0.01", "0.5"]
    levels_T = [1e4, 4e4]
    level_labels_T = [r"$10^4$K", r"$4 \times 10^4$K"]
    linestyles = ["dashed", "solid"]
    cnt1 = ax1.contour(
        x, y, HI_map.T, colors="C0", levels=levels_xHI, linestyles=linestyles
    )
    for i, l in enumerate(level_labels_xHI):
        cnt1.collections[i].set_label("GEARRT xHI=" + l)

    cnt2 = ax2.contour(
        x, y, temperature_map.T, colors="C0", levels=levels_T, linestyles=linestyles
    )
    for i, l in enumerate(level_labels_T):
        cnt2.collections[i].set_label("GEARRT T=" + l)

    # plot references now.
    #  references = ["C2Ray", "Crash", "FFTE"]
    references = ["C2Ray"]

    refedge = np.linspace(extent[0], extent[1], 128)
    xref, yref = np.meshgrid(refedge, refedge)
    for c, code in enumerate(references):
        xHI_ref = get_ref_data(code, "xHI")
        cnt1 = ax1.contour(
            xref,
            yref,
            xHI_ref.T,
            colors="C" + str(c + 1),
            levels=levels_xHI,
            linestyles=linestyles,
        )
        for i, l in enumerate(level_labels_xHI):
            cnt1.collections[i].set_label("{0} xHI=".format(code) + l)

        T_ref = get_ref_data(code, "T")
        cnt2 = ax2.contour(
            xref,
            yref,
            T_ref.T,
            colors="C" + str(c + 1),
            levels=levels_T,
            linestyles=linestyles,
        )
        for i, l in enumerate(level_labels_T):
            cnt2.collections[i].set_label("{0} T=".format(code) + l)

    ax1.legend(loc="lower left", ncol=2)
    ax2.legend(loc="lower left", ncol=2)

    title = "Iliev+06 Test 4"
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    figname = filename[:-5] + "-Contours.png"
    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_result(snap)
