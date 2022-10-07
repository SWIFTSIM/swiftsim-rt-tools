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
# Plot slices of the neutral hydrogen number density,
# and temperature
# ----------------------------------------------------

import sys

import matplotlib as mpl

mpl.use("Agg")

import swiftsimio
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from swiftsimio.visualisation.slice import slice_gas

from rainbow4_colormap import rainbow4

# Parameters users should/may tweak

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

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0 * meta.boxsize[0].v,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
        0.0 * meta.boxsize[1].v,
        (1.0 - 8 / (128 + 8)) * meta.boxsize[0].to("kpc").v,
    ]

    # cut off boundary particles
    cutoff = int(4 / (128 + 8) * slice_kwargs["resolution"])

    mass_map = slice_gas(
        data, project="masses", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    data.gas.m2 = data.gas.masses * data.gas.masses
    mass_weighted_m_map = slice_gas(
        data, project="m2", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    data.gas.mu = data.gas.internal_energy * data.gas.masses
    mass_weighted_u_map = slice_gas(
        data, project="mu", z_slice=0.5 * meta.boxsize[2], **slice_kwargs
    )

    m_img = mass_weighted_m_map / mass_map
    m_img = m_img[cutoff:-cutoff, cutoff:-cutoff]
    u_img = mass_weighted_u_map / mass_map
    u_img = u_img[cutoff:-cutoff, cutoff:-cutoff]

    fig = plt.figure(figsize=(12, 6), dpi=200)
    figname = filename[:-5] + ".png"

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    im1 = ax1.imshow(
        m_img.T,
        **imshow_kwargs,
        norm=LogNorm(),
        cmap=rainbow4,
        #  cmap = "cividis"
    )
    set_colorbar(ax1, im1)
    ax1.set_title("mass")

    im2 = ax2.imshow(
        u_img.T,
        **imshow_kwargs,
        norm=LogNorm(),
        cmap=rainbow4,
        #  cmap = "inferno"
    )
    set_colorbar(ax2, im2)
    ax2.set_title(r"internal energy")

    ax1.set_xlabel("[kpc]")
    ax1.set_ylabel("[kpc]")
    ax2.set_xlabel("[kpc]")
    ax2.set_ylabel("[kpc]")

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":

    filename = "ilievTest3.hdf5"
    plot_result(filename)
