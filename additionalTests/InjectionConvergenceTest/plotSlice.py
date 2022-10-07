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


# ---------------------------------------------------------------------
# Plot slices of the neutral hydrogen number density and temperature
# ---------------------------------------------------------------------

import gc
import sys

import matplotlib as mpl

mpl.use("Agg")
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from swiftsimio.visualisation.slice import slice_gas

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "output-128"

# parameters for imshow plots
imshow_kwargs = {"origin": "lower"}

# parameters for swiftsimio slices
slice_kwargs = {"resolution": 1024, "parallel": True}

# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
snapnr = -1
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_result(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0,
        meta.boxsize[0].to("kpc").v,
        0.0,
        meta.boxsize[1].to("kpc").v,
    ]

    mass_map = slice_gas(
        data, project="masses", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    gamma = meta.gas_gamma

    imf = spt.get_imf(scheme, data)
    data.gas.mXHI = imf.HI * data.gas.masses
    data.gas.mXHII = imf.HII * data.gas.masses

    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    data.gas.mT = (
        spt.gas_temperature(data.gas.internal_energies, mu, gamma) * data.gas.masses
    )
    data.gas.mT = data.gas.mT.to(data.gas.masses.units * unyt.K)

    mass_weighted_HI_map = slice_gas(
        data, project="mXHI", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_HII_map = slice_gas(
        data, project="mXHII", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )
    mass_weighted_temperature_map = slice_gas(
        data, project="mT", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    HI_map = mass_weighted_HI_map / mass_map
    HII_map = mass_weighted_HII_map / mass_map
    temperature_map = mass_weighted_temperature_map / mass_map
    temperature_map = temperature_map.to("K")

    fig = plt.figure(figsize=(10, 5), dpi=200)
    figname = filename[:-5] + "-slice.png"

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    im1 = ax1.imshow(
        HII_map.T,
        #  HI_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=1.0e-5, vmax=1.0),
        cmap="cividis",
    )
    set_colorbar(ax1, im1)
    #  ax1.set_title("Neutral Hydrogen Mass Fraction [1]")
    ax1.set_title("Ionized Hydrogen Mass Fraction [1]")

    im2 = ax2.imshow(
        temperature_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=1e2, vmax=5e4),
        cmap="inferno",
    )
    set_colorbar(ax2, im2)
    ax2.set_title(r"Temperature [K]")

    #  for ax in fig.axes:
    #      ax.set_xlabel("[kpc]")
    #      ax.set_ylabel("[kpc]")

    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += r", $z$ = {0:.2e}".format(meta.z)
    title += r", $t$ = {0:.1f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()
    return


if __name__ == "__main__":

    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)

    for f in snaplist:
        plot_result(f)
