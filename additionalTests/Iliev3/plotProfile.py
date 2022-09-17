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
# Plot xHI and T profiles along axis of symmetry
# Usage: ./plotProfile.py <snapshot number>
# ----------------------------------------------------

import sys
import swiftsimio
import gc
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from swiftsimio.visualisation.slice import slice_gas
import numpy as np

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# plot references?
plot_refs = True

# label references, or make them all grey?
label_refs = False

# which reference to use?
#  ref = "1Myr"
#  ref = "2Myr"
#  ref = "3Myr"
#  ref = "5Myr"
ref = "15Myr"

# snapshot basename
snapshot_base = "output"

# parameters for swiftsimio slices
slice_kwargs = {"resolution": 64, "parallel": True}
#  slice_kwargs = {"resolution": 128, "parallel": True}

# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?

mpl.rcParams["text.usetex"] = True


def plot_result(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

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

    HI_map = mass_weighted_HI_map / mass_map

    temperature_map = mass_weighted_temperature_map / mass_map

    n = HI_map.shape[0]
    x = np.linspace(0.5 / n, (n - 0.5), n) / n * meta.boxsize[0]

    L_test = 6.6 * unyt.kpc
    shift = (meta.boxsize[0] - L_test) * 0.5
    shiftint = int(shift / meta.boxsize[0] * n)
    x = x[shiftint:-shiftint] - shift

    HI_profile = HI_map.T[int(n / 2), shiftint:-shiftint]
    T_profile = temperature_map.T[int(n / 2), shiftint:-shiftint]

    fig = plt.figure(figsize=(12, 6), dpi=200)
    if plot_refs:
        figname = filename[:-5] + "-Profiles.png"
    else:
        figname = filename[:-5] + "-ProfilesNoRef.png"

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # first the references
    if plot_refs:
        codes = ["C2Ray", "Coral", "Crash", "FFTE", "Flash", "IFT", "RSPH"]
        for code in codes:
            Trefname = "reference/" + code + "_" + ref + "_profile_" + "T.dat"
            xHIrefname = "reference/" + code + "_" + ref + "_profile_" + "xHI.dat"

            Tref = np.loadtxt(Trefname)
            xHIref = np.loadtxt(xHIrefname)
            nref = Tref.shape[0]

            xref = np.linspace(0.0, 1.0, nref)

            if label_refs:
                ax1.plot(xref, xHIref, label=code)
                ax2.plot(xref, Tref, label=code)
            else:
                ax1.plot(xref, xHIref, c="grey", alpha=0.6)
                ax2.plot(xref, Tref, c="grey", alpha=0.6)

    ax1.semilogy(x / L_test, HI_profile, label="GEARRT")
    ax1.set_title("Neutral Hydrogen Mass Fraction [1]")
    ax1.legend()

    ax2.semilogy(x / L_test, T_profile, label="GEARRT")
    ax2.set_title(r"Temperature [K]")

    for ax in fig.axes:
        ax.set_xlabel("x/L")

    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2f}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()
    return


if __name__ == "__main__":

    snapnr = int(sys.argv[1])
    snap = snapshot_base + "_" + str(snapnr).zfill(4) + ".hdf5"
    plot_result(snap)
