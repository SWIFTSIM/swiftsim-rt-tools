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


# ---------------------------------------------------------------
# Plot slices of neutral hydrogen mass fraction, hydrogen ion
# mass fraction, gas number density, temperature, pressure, and
# mach number without reference solutions.
# ---------------------------------------------------------------

import sys
import swiftsimio
import gc
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from swiftsimio.visualisation.slice import slice_gas
import numpy as np

from rainbow4_colormap import rainbow4
import stromgren_plotting_tools as spt

# Parameters users should/may tweak

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
    "figure.subplot.left": 0.145,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.075,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.0,
    "figure.dpi": 200,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)


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
    # todo: get actual mach nr

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
    mach_map = mass_weighted_mach_map[cutoff:-cutoff, cutoff:-cutoff]

    fig = plt.figure(figsize=(18, 12), dpi=200)
    figname = filename[:-5] + "-NoRef.png"

    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    im1 = ax1.imshow(
        HI_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=1.0e-7, vmax=1.2),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax1, im1)
    ax1.set_title("Neutral Hydrogen Mass Fraction [1]")

    im2 = ax2.imshow(
        HII_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=1.0e-7, vmax=1.2),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax2, im2)
    ax2.set_title("Ionized Hydrogen Mass Fraction [1]")

    im3 = ax3.imshow(
        number_density_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=2.5e-4, vmax=2.5e-3),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax3, im3)
    ax3.set_title(r"Hydrogen Number Density [cm$^{-3}$]")

    im4 = ax4.imshow(
        temperature_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=80, vmax=1e5),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax4, im4)
    ax4.set_title(r"Temperature [K]")

    im5 = ax5.imshow(
        pressure_map.T,
        **imshow_kwargs,
        norm=LogNorm(vmin=1e-17, vmax=1e-14),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax5, im5)
    ax5.set_title(r"Pressure [g cm$^{-1}$ s$^{-2}$]")

    im6 = ax6.imshow(
        mach_map.T,
        **imshow_kwargs,
        norm=LogNorm(),
        #  norm=LogNorm(vmin=1e-3, vmax=1e1),
        #  cmap=rainbow4,
        cmap="cividis",
    )
    set_colorbar(ax6, im6)
    ax6.set_title(r"Mach Number [1]")

    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.set_xlabel("[kpc]")
        ax.set_ylabel("[kpc]")

    title = "Iliev+09 Test 6"
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.0f}".format(meta.time.to("Myr"))
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
